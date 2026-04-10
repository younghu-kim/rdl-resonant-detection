#!/usr/bin/env python3
"""
=============================================================================
Supervised-Ψ Phase Experiment — rem:psi_defect 의 진단 검정 (2026-04-08)
=============================================================================
phase_L4_experiment.py 의 단일 변경 검정:

  (a) Ψ_target 을 IdealPhaseTargetGenerator (||X||-종속 합성 곡선) 에서
      arg ξ(½+it) (= xi_phase) 로 교체.
  (b) cos² PQO 를 학습 손실에서 제거 (λ_pqo = 0). PQO 는 측정 도구로만.
  (c) 동일 seed/hidden/epochs/dataset 으로 L=4 negative 베이스라인과 비교.

귀무가설 H0: Ψ-decoupling 진단이 틀렸다 — supervised Ψ + λ_pqo=0 으로
             교체해도 영점/비영점 phase MSE 분리는 일어나지 않는다.

대립가설 H1: 진단이 옳다 — supervised Ψ 가 ξ-영점 정보를 학습 신호로
             주입하면, arg(Z) 가 영점 근처에서 ±π 점프를 추적하기 시작하고
             영점/비영점 phase MSE 가 통계적으로 분리된다.

측정:
  1. arg(Z) vs arg(ξ) 의 phase MSE — 영점 vs 비영점 (분리되면 신호)
  2. cos² L=4 cyclotomic k 분포 — 학습 손실 밖에서 측정만
  3. 영점 근처 ±π 점프 검출 (sign(Re Z) 변화)
"""
import os
import sys
import json
import time
import copy
import math
import numpy as np
import torch
from scipy.stats import chi2_contingency, ttest_ind

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_ROOT)

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.pipeline.xi_feature_dataset import XiFeatureDataset

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs", "analysis")
OVERNIGHT_DIR = os.path.join(PROJECT_ROOT, "outputs", "overnight")
JSON_OUT = os.path.join(OUTPUT_DIR, "phase_supervised_experiment.json")
TXT_OUT = os.path.join(OUTPUT_DIR, "phase_supervised_experiment.txt")

# 실험 구성 — phase_L4_experiment.py 와 동일 (Ψ + λ_pqo 만 변경)
L_PQO = 4
HIDDEN = 32
SEEDS = [42, 123, 777]
EPOCHS = 200
LR = 0.001
IN_FEATURES = 64
BATCH_SIZE = 32
T_MIN, T_MAX = 100.0, 200.0
N_NONZERO = 200
NONZERO_GAP = 0.5

ATTRACTORS = np.array([(2 * k + 1) * math.pi / L_PQO for k in range(L_PQO)])


class XiFeatureDatasetSupervised(XiFeatureDataset):
    """
    부모 클래스가 (features, [Re ξ, Im ξ]) 를 반환하는 것에 더해,
    arg ξ(½+it) (= xi_phase) 를 supervisory phase 로 함께 반환.

    이 phase 가 학습 루프에서 모델의 broken Ψ_target 을 대체한다.
    """
    def __getitem__(self, idx):
        features = self._build_features(self.t[idx])
        xi_target = torch.stack([self.xi_real[idx], self.xi_imag[idx]])
        # arg ξ(½+it) — 임계선 위에서 ξ 는 실수이므로 위상은 {0, π}
        # 영점 교차에서 π 점프가 일어남 (eq:phase_crit, unified_en.tex)
        psi_super = self.xi_phase[idx].unsqueeze(0)  # shape [1]
        return features, xi_target, psi_super


def load_zeros():
    with open(os.path.join(OVERNIGHT_DIR, "result_t100-200.json")) as f:
        return json.load(f)["zeros_list"]


def sample_nonzero_t(zeros, n=N_NONZERO, seed=0):
    rng = np.random.default_rng(seed)
    candidates = rng.uniform(T_MIN, T_MAX, size=n * 5)
    z = np.array(zeros)
    keep = []
    for t in candidates:
        if np.min(np.abs(z - t)) >= NONZERO_GAP:
            keep.append(t)
            if len(keep) >= n:
                break
    return sorted(keep)


def cyclotomic_label(phi_array):
    phi_mod = np.mod(phi_array, 2 * math.pi)
    diffs = np.abs(phi_mod[:, None] - ATTRACTORS[None, :])
    diffs = np.minimum(diffs, 2 * math.pi - diffs)
    return np.argmin(diffs, axis=1)


def phase_distance(phi, target):
    """원형 거리: min(|d|, 2π-|d|), d = phi - target"""
    d = np.mod(phi - target + math.pi, 2 * math.pi) - math.pi
    return np.abs(d)


def train_model(train_loader, val_loader, hidden, seed, device, label):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = MasterResonantNetwork(
        in_features=IN_FEATURES, hidden_features=hidden, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    # 핵심 변경: λ_pqo = 0 — PQO 는 학습 손실에서 완전히 제거
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.0,
        pqo_L=L_PQO, pqo_mode="cos2",
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    best_val = float("inf")
    best_state = None
    t0 = time.time()

    for epoch in range(1, EPOCHS + 1):
        model.train()
        for X_batch, _, psi_super in train_loader:
            X_in = X_batch.to(device, dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            psi_super_b = psi_super.to(device, dtype=PrecisionManager.REAL_DTYPE)

            optimizer.zero_grad(set_to_none=True)
            outputs = model(X_in)

            # 핵심 변경: 모델이 산출한 broken Ψ_target 을
            # arg ξ(½+it) supervisory phase 로 덮어쓴다.
            outputs["Psi_target"] = psi_super_b

            total_loss, _ = loss_fn(**outputs)
            if torch.isnan(total_loss) or torch.isinf(total_loss):
                continue
            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()

        model.eval()
        val_acc = 0.0
        n_val = 0
        with torch.enable_grad():
            for X_batch, _, psi_super in val_loader:
                X_in = X_batch.to(device, dtype=PrecisionManager.REAL_DTYPE)
                X_in.requires_grad_(True)
                psi_super_b = psi_super.to(device, dtype=PrecisionManager.REAL_DTYPE)
                out = model(X_in)
                out["Psi_target"] = psi_super_b
                vl, _ = loss_fn(**out)
                if not (torch.isnan(vl) or torch.isinf(vl)):
                    val_acc += vl.item()
                    n_val += 1
        val_loss = val_acc / max(1, n_val)
        if val_loss < best_val:
            best_val = val_loss
            best_state = copy.deepcopy(model.state_dict())
        if epoch % 50 == 0 or epoch == EPOCHS:
            print(f"    {label} ep {epoch:03d}: val={val_loss:.4f}")

    if best_state:
        model.load_state_dict(best_state)
    return model, best_val, time.time() - t0


def extract_phase_and_target(model, dataset, t_values, device):
    """모델 출력 phase + 동일 t 의 ground-truth arg ξ 를 함께 반환"""
    feats = dataset.get_features_at_t(t_values).to(device)
    model.eval()
    with torch.enable_grad():
        X = feats.clone().requires_grad_(True)
        out = model(X)
    phi_pred = torch.angle(out["Z_out"]).mean(dim=-1).detach().cpu().numpy()

    # ground-truth phase: dataset 캐시에서 보간 없이 가장 가까운 인덱스 사용
    t_arr = dataset.t.cpu().numpy()
    phi_true = []
    for t in t_values:
        idx = int(np.argmin(np.abs(t_arr - t)))
        phi_true.append(float(dataset.xi_phase[idx]))
    return phi_pred, np.array(phi_true)


def main():
    PrecisionManager.setup_precision()
    device = torch.device("cpu")

    print("=" * 70)
    print(f"  Supervised-Ψ 검정 실험 — Ψ := arg ξ(½+it), λ_pqo=0")
    print(f"  hidden={HIDDEN}, seeds={SEEDS}, epochs={EPOCHS}, t∈[{T_MIN},{T_MAX}]")
    print("=" * 70)

    zeros = load_zeros()
    nonzeros = sample_nonzero_t(zeros)
    print(f"  영점 {len(zeros)}개, 비영점 {len(nonzeros)}개")

    cache_data = torch.load(
        os.path.join(OVERNIGHT_DIR, "xi_cache_t100-200_n500.pt"),
        weights_only=True
    )
    dataset = XiFeatureDatasetSupervised(cache_data, in_features=IN_FEATURES)

    val_size = int(len(dataset) * 0.2)
    train_size = len(dataset) - val_size
    train_ds, val_ds = torch.utils.data.random_split(
        dataset, [train_size, val_size],
        generator=torch.Generator().manual_seed(0)
    )
    train_loader = torch.utils.data.DataLoader(
        train_ds, batch_size=BATCH_SIZE, shuffle=True, drop_last=True)
    val_loader = torch.utils.data.DataLoader(
        val_ds, batch_size=BATCH_SIZE, shuffle=False, drop_last=False)

    results = {}
    for seed in SEEDS:
        label = f"sup,h={HIDDEN},s={seed}"
        print(f"\n  [훈련] {label}")
        model, best_val, train_time = train_model(
            train_loader, val_loader, HIDDEN, seed, device, label
        )
        print(f"    완료: best_val={best_val:.4f}, time={train_time:.1f}s")

        phi_z_pred, phi_z_true = extract_phase_and_target(model, dataset, zeros, device)
        phi_nz_pred, phi_nz_true = extract_phase_and_target(model, dataset, nonzeros, device)

        # phase MSE: 모델 예측 vs ground-truth arg ξ
        mse_z = phase_distance(phi_z_pred, phi_z_true) ** 2
        mse_nz = phase_distance(phi_nz_pred, phi_nz_true) ** 2

        # cyclotomic k 분포 (측정만, 학습에는 안 씀)
        k_z = cyclotomic_label(phi_z_pred)
        k_nz = cyclotomic_label(phi_nz_pred)

        results[seed] = {
            "phi_z_pred": phi_z_pred.tolist(),
            "phi_z_true": phi_z_true.tolist(),
            "phi_nz_pred": phi_nz_pred.tolist(),
            "phi_nz_true": phi_nz_true.tolist(),
            "mse_z": mse_z.tolist(),
            "mse_nz": mse_nz.tolist(),
            "k_z": k_z.tolist(),
            "k_nz": k_nz.tolist(),
            "best_val_loss": float(best_val),
            "train_time": float(train_time),
        }

    out_data = {
        "config": {
            "L_pqo": L_PQO, "lambda_pqo": 0.0, "psi_source": "arg xi(1/2+it)",
            "hidden": HIDDEN, "seeds": SEEDS, "epochs": EPOCHS,
            "t_min": T_MIN, "t_max": T_MAX,
            "n_zeros": len(zeros), "n_nonzeros": len(nonzeros),
            "nonzero_gap": NONZERO_GAP, "attractors": ATTRACTORS.tolist(),
        },
        "zeros": zeros, "nonzeros": nonzeros,
        "results": {str(k): v for k, v in results.items()},
    }
    with open(JSON_OUT, "w") as f:
        json.dump(out_data, f, indent=2, ensure_ascii=False)

    # ===== 통계 분석 =====
    lines = []
    lines.append("=" * 70)
    lines.append(f"  Supervised-Ψ 실험 결과 — Ψ := arg ξ(½+it), λ_pqo=0")
    lines.append(f"  hidden={HIDDEN}, seeds={SEEDS}, epochs={EPOCHS}")
    lines.append(f"  영점 {len(zeros)}개, 비영점 {len(nonzeros)}개 (gap≥{NONZERO_GAP})")
    lines.append("=" * 70)

    # 1. phase MSE 영점 vs 비영점
    lines.append("")
    lines.append("=" * 70)
    lines.append("  1. Phase MSE: arg(Z) vs arg(ξ)  — 영점 vs 비영점 (핵심)")
    lines.append("=" * 70)
    lines.append(f"  {'seed':>6s} {'영점 MSE':>14s} {'비영점 MSE':>14s} "
                 f"{'비율(z/nz)':>12s} {'t-test p':>12s}")
    for seed in SEEDS:
        mz = np.array(results[seed]["mse_z"])
        mn = np.array(results[seed]["mse_nz"])
        ratio = mz.mean() / max(mn.mean(), 1e-12)
        try:
            t_stat, p_val = ttest_ind(mz, mn, equal_var=False)
        except Exception:
            p_val = float('nan')
        lines.append(f"  {seed:>6d} {mz.mean():>14.6f} {mn.mean():>14.6f} "
                     f"{ratio:>12.4f} {p_val:>12.4e}")

    # 2. cyclotomic k 분포 (순수 측정)
    lines.append("")
    lines.append("=" * 70)
    lines.append("  2. Cyclotomic k(t) 측정 (학습 손실 외부)")
    lines.append("=" * 70)
    for seed in SEEDS:
        kz = np.array(results[seed]["k_z"])
        knz = np.array(results[seed]["k_nz"])
        z_counts = np.bincount(kz, minlength=L_PQO)
        nz_counts = np.bincount(knz, minlength=L_PQO)
        lines.append(f"  s={seed} 영점:   " + " ".join(f"k{k}={c:>3d}" for k, c in enumerate(z_counts)))
        lines.append(f"  s={seed} 비영점: " + " ".join(f"k{k}={c:>3d}" for k, c in enumerate(nz_counts)))
        table = np.array([z_counts, nz_counts])
        if (table == 0).any():
            table = table + 1
        try:
            chi2, p, dof, _ = chi2_contingency(table)
            lines.append(f"     χ²={chi2:.4f}, dof={dof}, p={p:.4e}")
        except Exception as e:
            lines.append(f"     χ² 검정 실패: {e}")

    # 3. 시드 풀 통합
    lines.append("")
    lines.append("=" * 70)
    lines.append("  3. 시드 풀 통합 phase MSE")
    lines.append("=" * 70)
    all_mz = np.concatenate([np.array(results[s]["mse_z"]) for s in SEEDS])
    all_mn = np.concatenate([np.array(results[s]["mse_nz"]) for s in SEEDS])
    pool_ratio = all_mz.mean() / max(all_mn.mean(), 1e-12)
    try:
        pool_t, pool_p = ttest_ind(all_mz, all_mn, equal_var=False)
    except Exception:
        pool_p = float('nan')
    lines.append(f"  영점 MSE 평균   = {all_mz.mean():.6f} (n={len(all_mz)})")
    lines.append(f"  비영점 MSE 평균 = {all_mn.mean():.6f} (n={len(all_mn)})")
    lines.append(f"  비율 z/nz       = {pool_ratio:.4f}")
    lines.append(f"  Welch t-test p  = {pool_p:.4e}")

    # 4. 판정
    lines.append("")
    lines.append("=" * 70)
    lines.append("  4. 판정")
    lines.append("=" * 70)
    if pool_p < 0.01 and pool_ratio < 0.5:
        verdict = (
            "→ 강한 양성: supervised Ψ 가 영점 정보를 학습에 주입하는 데 성공.\n"
            "  Ψ-decoupling 진단 (rem:psi_defect) 옳음. cos²+L=2/4 음성 결과는\n"
            "  Ψ 결함의 결과였음이 확인됨."
        )
    elif pool_p < 0.05 and pool_ratio < 0.8:
        verdict = (
            "→ 약한 양성: supervised Ψ 로 분리는 일어나지만 미약함.\n"
            "  진단 방향은 옳으나 추가 손실 항/구조 조정 필요."
        )
    else:
        verdict = (
            "→ 음성: supervised Ψ 로도 분리 안 됨.\n"
            "  진단이 더 깊다 — Ψ 교체만으로 부족하거나 모델 자체의 표현력이\n"
            "  arg(Z)→±π 점프를 학습할 수 없음."
        )
    lines.append(f"  {verdict}")

    output = "\n".join(lines)
    print("\n" + output)
    with open(TXT_OUT, "w") as f:
        f.write(output)
    print(f"\n저장: {JSON_OUT}")
    print(f"저장: {TXT_OUT}")


if __name__ == "__main__":
    main()
