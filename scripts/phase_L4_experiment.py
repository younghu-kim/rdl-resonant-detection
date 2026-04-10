#!/usr/bin/env python3
"""
=============================================================================
L=4 Cyclotomic PQO 실험 — 어제 발견의 진짜 검정
=============================================================================
가설: cos² PQO를 L=2 → L=4로 일반화하면, 끌개가 4개 (π/4, 3π/4, 5π/4, 7π/4)
가 되어 영점이 어느 시클로토믹 섹터 ω_k에 떨어지는지 라벨링 가능.

귀무가설 H0: 영점과 비영점의 k(t) 분포가 동일 (=L=4도 분해 못 함)
대립가설 H1: 영점은 특정 k 섹터에 편중 (= 진짜 신호)

label k(t) = arg min_k |arg Z(t) - (2k+1)π/4|, k ∈ {0,1,2,3}
대응 단위근:
  k=0: ω_0 = e^{iπ/4}     (1+i)/√2
  k=1: ω_1 = e^{i3π/4}    (-1+i)/√2
  k=2: ω_2 = e^{i5π/4}    (-1-i)/√2
  k=3: ω_3 = e^{i7π/4}    (1-i)/√2

설정: hidden=32, seeds=[42,123,777], EPOCHS=200 (이전 실험과 동일)
"""
import os
import sys
import json
import time
import copy
import math
import numpy as np
import torch
from scipy.stats import chi2_contingency, entropy

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_ROOT)

from gdl.rdl.constants import R_CONST, PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.pipeline.xi_feature_dataset import XiFeatureDataset

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs", "analysis")
OVERNIGHT_DIR = os.path.join(PROJECT_ROOT, "outputs", "overnight")
JSON_OUT = os.path.join(OUTPUT_DIR, "phase_L4_experiment.json")
TXT_OUT = os.path.join(OUTPUT_DIR, "phase_L4_experiment.txt")

L_PQO = 4   # ← 핵심 변경: L=2 → L=4
HIDDEN = 32
SEEDS = [42, 123, 777]
EPOCHS = 200
LR = 0.001
IN_FEATURES = 64
BATCH_SIZE = 32
T_MIN, T_MAX = 100.0, 200.0
N_NONZERO = 200
NONZERO_GAP = 0.5

# L=4 cyclotomic 끌개 위치
ATTRACTORS = np.array([(2 * k + 1) * math.pi / L_PQO for k in range(L_PQO)])
# = [π/4, 3π/4, 5π/4, 7π/4]


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
    """φ → k ∈ {0,1,2,3}, 가장 가까운 ω_k 섹터."""
    phi_mod = np.mod(phi_array, 2 * math.pi)
    diffs = np.abs(phi_mod[:, None] - ATTRACTORS[None, :])
    diffs = np.minimum(diffs, 2 * math.pi - diffs)  # 원형 거리
    return np.argmin(diffs, axis=1)


def train_model(train_loader, val_loader, hidden, seed, device, label):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = MasterResonantNetwork(
        in_features=IN_FEATURES, hidden_features=hidden, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.5,
        pqo_L=L_PQO,         # ← 핵심: L=4
        pqo_mode="cos2",
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    best_val = float("inf")
    best_state = None
    t0 = time.time()
    for epoch in range(1, EPOCHS + 1):
        model.train()
        for X_batch, _ in train_loader:
            X_in = X_batch.to(device, dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            optimizer.zero_grad(set_to_none=True)
            outputs = model(X_in)
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
            for X_batch, _ in val_loader:
                X_in = X_batch.to(device, dtype=PrecisionManager.REAL_DTYPE)
                X_in.requires_grad_(True)
                out = model(X_in)
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


def extract_phase(model, dataset, t_values, device):
    feats = dataset.get_features_at_t(t_values).to(device)
    model.eval()
    with torch.enable_grad():
        X = feats.clone().requires_grad_(True)
        out = model(X)
    return torch.angle(out["Z_out"]).mean(dim=-1).detach().cpu().numpy()


def main():
    PrecisionManager.setup_precision()
    device = torch.device("cpu")

    print("=" * 70)
    print(f"  L={L_PQO} Cyclotomic PQO 실험")
    print(f"  끌개: π/{L_PQO}, 3π/{L_PQO}, ..., (2L-1)π/{L_PQO}")
    print("=" * 70)

    zeros = load_zeros()
    nonzeros = sample_nonzero_t(zeros)
    print(f"  영점 {len(zeros)}개, 비영점 {len(nonzeros)}개")

    cache_data = torch.load(
        os.path.join(OVERNIGHT_DIR, "xi_cache_t100-200_n500.pt"),
        weights_only=True
    )
    dataset = XiFeatureDataset(cache_data, in_features=IN_FEATURES)

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
        label = f"L=4,h={HIDDEN},s={seed}"
        print(f"\n  [훈련] {label}")
        model, best_val, train_time = train_model(
            train_loader, val_loader, HIDDEN, seed, device, label
        )
        print(f"    완료: best_val={best_val:.4f}, time={train_time:.1f}s")

        phi_zero = extract_phase(model, dataset, zeros, device)
        phi_nonzero = extract_phase(model, dataset, nonzeros, device)
        k_zero = cyclotomic_label(phi_zero)
        k_nonzero = cyclotomic_label(phi_nonzero)

        results[seed] = {
            "phi_zero": phi_zero.tolist(),
            "phi_nonzero": phi_nonzero.tolist(),
            "k_zero": k_zero.tolist(),
            "k_nonzero": k_nonzero.tolist(),
            "best_val_loss": float(best_val),
            "train_time": float(train_time),
        }

    out_data = {
        "config": {
            "L_pqo": L_PQO, "hidden": HIDDEN, "seeds": SEEDS, "epochs": EPOCHS,
            "t_min": T_MIN, "t_max": T_MAX,
            "n_zeros": len(zeros), "n_nonzeros": len(nonzeros),
            "nonzero_gap": NONZERO_GAP,
            "attractors": ATTRACTORS.tolist(),
        },
        "zeros": zeros,
        "nonzeros": nonzeros,
        "results": {str(k): v for k, v in results.items()},
    }
    with open(JSON_OUT, "w") as f:
        json.dump(out_data, f, indent=2, ensure_ascii=False)

    # ===== 통계 분석 =====
    lines = []
    lines.append("=" * 70)
    lines.append(f"  L={L_PQO} Cyclotomic PQO 실험 결과")
    lines.append(f"  hidden={HIDDEN}, seeds={SEEDS}, epochs={EPOCHS}")
    lines.append(f"  영점 {len(zeros)}개, 비영점 {len(nonzeros)}개 (gap≥{NONZERO_GAP})")
    lines.append(f"  L={L_PQO} 끌개: {[f'{a:.4f}' for a in ATTRACTORS]}")
    lines.append("=" * 70)

    # 1. 시드별 φ 통계
    lines.append("")
    lines.append("=" * 70)
    lines.append("  1. 시드별 φ 통계 (영점 vs 비영점)")
    lines.append("=" * 70)
    lines.append(f"  {'seed':>6s} {'영점 φ평균':>12s} {'영점 φstd':>12s} "
                 f"{'비영점 φ평균':>14s} {'비영점 φstd':>14s}")
    for seed in SEEDS:
        z_phi = np.array(results[seed]["phi_zero"])
        nz_phi = np.array(results[seed]["phi_nonzero"])
        lines.append(f"  {seed:>6d} {z_phi.mean():>12.6f} {z_phi.std():>12.6f} "
                     f"{nz_phi.mean():>14.6f} {nz_phi.std():>14.6f}")

    # 2. 시드별 k(t) 분포 — 핵심
    lines.append("")
    lines.append("=" * 70)
    lines.append("  2. 시드별 k(t) 시클로토믹 섹터 분포 (핵심)")
    lines.append("=" * 70)
    lines.append(f"  {'seed':>6s} {'set':>10s} "
                 + " ".join(f"k={k}".rjust(8) for k in range(L_PQO)) + "  엔트로피")
    for seed in SEEDS:
        kz = np.array(results[seed]["k_zero"])
        knz = np.array(results[seed]["k_nonzero"])
        z_counts = np.bincount(kz, minlength=L_PQO)
        nz_counts = np.bincount(knz, minlength=L_PQO)
        z_probs = z_counts / z_counts.sum()
        nz_probs = nz_counts / nz_counts.sum()
        ent_z = entropy(z_probs + 1e-12)
        ent_nz = entropy(nz_probs + 1e-12)
        lines.append(f"  {seed:>6d} {'영점':>10s} "
                     + " ".join(f"{c:>8d}" for c in z_counts) + f"  {ent_z:.4f}")
        lines.append(f"  {seed:>6d} {'비영점':>10s} "
                     + " ".join(f"{c:>8d}" for c in nz_counts) + f"  {ent_nz:.4f}")

    lines.append("")
    lines.append(f"  최대 엔트로피 (균등 분포 H0) = log({L_PQO}) = {math.log(L_PQO):.4f}")
    lines.append("  → 영점 엔트로피가 최대보다 유의하게 작으면 = 섹터 편중 (신호)")

    # 3. χ² 검정 (영점 vs 비영점 분포 차이)
    lines.append("")
    lines.append("=" * 70)
    lines.append("  3. χ² 검정: 영점 vs 비영점 k 분포 차이")
    lines.append("=" * 70)
    for seed in SEEDS:
        kz = np.array(results[seed]["k_zero"])
        knz = np.array(results[seed]["k_nonzero"])
        z_counts = np.bincount(kz, minlength=L_PQO)
        nz_counts = np.bincount(knz, minlength=L_PQO)
        # 2×L 분할표
        table = np.array([z_counts, nz_counts])
        # 0 셀이 있으면 검정 불가능 → 1로 패딩
        if (table == 0).any():
            table = table + 1
        try:
            chi2, p, dof, _ = chi2_contingency(table)
            lines.append(f"  s={seed}: χ²={chi2:.4f}, dof={dof}, p={p:.4e}")
        except Exception as e:
            lines.append(f"  s={seed}: 검정 실패 — {e}")

    # 4. 시드 합산 (전체 풀)
    lines.append("")
    lines.append("=" * 70)
    lines.append("  4. 시드 풀 전체 분포")
    lines.append("=" * 70)
    all_kz = np.concatenate([np.array(results[s]["k_zero"]) for s in SEEDS])
    all_knz = np.concatenate([np.array(results[s]["k_nonzero"]) for s in SEEDS])
    z_total = np.bincount(all_kz, minlength=L_PQO)
    nz_total = np.bincount(all_knz, minlength=L_PQO)
    lines.append(f"  영점 합계 (n={len(all_kz)}):")
    for k in range(L_PQO):
        bar = "#" * int(40 * z_total[k] / max(1, z_total.max()))
        pct = 100 * z_total[k] / z_total.sum()
        lines.append(f"    k={k} (φ={ATTRACTORS[k]:.4f}): {z_total[k]:>4d} "
                     f"({pct:>5.1f}%) {bar}")
    lines.append(f"  비영점 합계 (n={len(all_knz)}):")
    for k in range(L_PQO):
        bar = "#" * int(40 * nz_total[k] / max(1, nz_total.max()))
        pct = 100 * nz_total[k] / nz_total.sum()
        lines.append(f"    k={k} (φ={ATTRACTORS[k]:.4f}): {nz_total[k]:>4d} "
                     f"({pct:>5.1f}%) {bar}")

    table_total = np.array([z_total, nz_total])
    if (table_total == 0).any():
        table_total = table_total + 1
    chi2, p, dof, _ = chi2_contingency(table_total)
    lines.append("")
    lines.append(f"  전체 χ² 검정: χ²={chi2:.4f}, dof={dof}, p={p:.4e}")

    # 5. 끌개 잠금 정도 (|φ - 가장 가까운 ω_k|)
    lines.append("")
    lines.append("=" * 70)
    lines.append("  5. 잠금 정도: |φ - nearest attractor|")
    lines.append("=" * 70)
    for seed in SEEDS:
        z_phi = np.array(results[seed]["phi_zero"])
        nz_phi = np.array(results[seed]["phi_nonzero"])
        z_phi_mod = np.mod(z_phi, 2 * math.pi)
        nz_phi_mod = np.mod(nz_phi, 2 * math.pi)
        z_dist = np.min(np.abs(z_phi_mod[:, None] - ATTRACTORS[None, :]), axis=1)
        nz_dist = np.min(np.abs(nz_phi_mod[:, None] - ATTRACTORS[None, :]), axis=1)
        lines.append(f"  s={seed}: 영점 잔차 평균={z_dist.mean():.6f}, "
                     f"비영점={nz_dist.mean():.6f}")

    # 6. 판정
    lines.append("")
    lines.append("=" * 70)
    lines.append("  6. 판정")
    lines.append("=" * 70)
    z_ent = entropy((z_total + 1e-12) / z_total.sum())
    max_ent = math.log(L_PQO)
    ent_ratio = z_ent / max_ent
    lines.append(f"  영점 분포 엔트로피 / 최대 = {z_ent:.4f}/{max_ent:.4f} "
                 f"= {ent_ratio:.4f}")
    if ent_ratio < 0.7 and p < 0.05:
        verdict = "→ 강한 섹터 편중 + 영점/비영점 분포 차이 (신호 가능)"
    elif p < 0.05:
        verdict = "→ 영점/비영점 분포 차이 있음 (약한 신호)"
    elif ent_ratio < 0.7:
        verdict = "→ 분포 편중 있으나 비영점도 같은 편중 (모델 아티팩트)"
    else:
        verdict = "→ 균등 분포, 신호 없음 또는 L=4도 분해 부족"
    lines.append(f"  {verdict}")

    output = "\n".join(lines)
    print("\n" + output)
    with open(TXT_OUT, "w") as f:
        f.write(output)
    print(f"\n저장: {JSON_OUT}")
    print(f"저장: {TXT_OUT}")


if __name__ == "__main__":
    main()
