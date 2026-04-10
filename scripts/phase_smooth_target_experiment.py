#!/usr/bin/env python3
"""
=============================================================================
Phase Smooth-Target Experiment — Task #17, rem:isolation 후속 (2026-04-09)
=============================================================================
isolation 실험 (phase_isolation_experiment.py) 은 λ_res=λ_curv=λ_pqo=0 +
circular mean 으로 supervised chain 건전성을 확인했으나, 영점 근방의 ±π
이산 점프를 부드러운 R-값 네트워크가 표현할 수 없다는 **출력 표현 결함**을
드러냈다 (영점 MSE ≈ 2.17 ≈ (π/2)², 비영점 MSE ≈ 0.12).

rem:isolation 의 제안 해법 중 첫 번째를 시험한다:

  Ψ_smooth(t) = (π/2) · (1 - σ_sign(t))

  σ_sign(t) := ∏_{k=1}^{K} tanh((t - t_k) / σ)

여기서 {t_k} 는 t∈[100,200] 구간의 리만 영점 50개이며, σ 는 전이 폭(t-공간
단위). 이 함수는 arg ξ(½+it) 의 {0, π} 사각파와 모든 영점에서 부호가
일치하지만, 영점 근방에서 σ-폭의 부드러운 S-커브로 전이한다. 따라서
매끄러운 네트워크도 이 목표를 표현 가능하다.

가설:
  H_smooth: isolation 실험의 영점 MSE ≈ (π/2)² 는 순수 표현력 병목이었다.
            목표를 σ-부드럽게 바꾸면 영점 MSE 가 σ→0 극한에서 (π/2)² 로
            회귀하고, σ 가 합리적인 범위(≈0.2)일 때는 비영점 MSE 수준으로
            떨어진다.

스크립트는 σ ∈ {0.1, 0.2, 0.3} 에 대해 isolation 실험과 동일한 구조/시드로
재학습하고, 세 가지 모두에서 영점/비영점 phase MSE 비율과 Welch p 를
기록한다.
"""
import os
import sys
import json
import time
import copy
import math
import numpy as np
import torch
from scipy.stats import ttest_ind

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_ROOT)

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.pipeline.xi_feature_dataset import XiFeatureDataset
from scripts.phase_supervised_experiment import (
    load_zeros, sample_nonzero_t, phase_distance,
    L_PQO, IN_FEATURES, BATCH_SIZE, T_MIN, T_MAX,
)
from scripts.phase_isolation_experiment import circular_mean

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs", "analysis")
OVERNIGHT_DIR = os.path.join(PROJECT_ROOT, "outputs", "overnight")
JSON_OUT = os.path.join(OUTPUT_DIR, "phase_smooth_target_experiment.json")
TXT_OUT = os.path.join(OUTPUT_DIR, "phase_smooth_target_experiment.txt")

HIDDEN = 32
SEEDS = [42, 123, 777]
EPOCHS = 200
LR = 0.001
SIGMAS = [0.1, 0.2, 0.3]  # 전이 폭 sweep (t-공간)


def build_smooth_psi(t_values: np.ndarray, zeros: list, sigma: float,
                    xi_real_ref: np.ndarray) -> np.ndarray:
    """
    Ψ_smooth(t) = (π/2)(1 - s(t)),
        s(t) = global_sign · ∏_k tanh((t - z_k)/σ)

    global_sign 은 xi_real_ref 의 부호와 일치시키기 위해 필요할 경우 뒤집는다
    (∏ tanh 의 절대 부호는 영점 개수의 패리티와 t=0 쪽 부호에 의존).
    """
    z = np.asarray(zeros, dtype=np.float64)
    # shape [T, K]
    arg = (t_values[:, None] - z[None, :]) / sigma
    th = np.tanh(arg)
    s = np.prod(th, axis=1)
    # xi_real_ref 와 부호 정합
    sign_xr = np.sign(xi_real_ref)
    agree = (np.sign(s) == sign_xr).mean()
    if agree < 0.5:
        s = -s
    psi = (math.pi / 2.0) * (1.0 - s)
    return psi  # [0, π] 범위


class XiFeatureDatasetSmooth(XiFeatureDataset):
    """
    부모 클래스에 더해 미리 계산된 smooth Ψ 를 반환한다.

    psi_smooth_full: shape [N_t], t 인덱스별 부드러운 목표.
    """
    def __init__(self, cache_data, in_features, psi_smooth_full):
        super().__init__(cache_data, in_features=in_features)
        self.psi_smooth_full = torch.as_tensor(
            psi_smooth_full, dtype=PrecisionManager.REAL_DTYPE
        )

    def __getitem__(self, idx):
        features = self._build_features(self.t[idx])
        xi_target = torch.stack([self.xi_real[idx], self.xi_imag[idx]])
        psi_super = self.psi_smooth_full[idx].unsqueeze(0)
        return features, xi_target, psi_super

    def get_features_at_t(self, t_values):
        """isolation 실험과 호환: XiFeatureDataset 의 기존 메서드 사용."""
        return super().get_features_at_t(t_values)


def train_smooth(train_loader, val_loader, hidden, seed, device, label):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = MasterResonantNetwork(
        in_features=IN_FEATURES, hidden_features=hidden, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    loss_fn = TotalResonanceLoss(
        lambda_res=0.0, lambda_curv=0.0, lambda_tgt=1.0, lambda_pqo=0.0,
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
            print(f"    {label} ep {epoch:03d}: val={val_loss:.6f}", flush=True)

    if best_state:
        model.load_state_dict(best_state)
    return model, best_val, time.time() - t0


def extract_phases_with_target(model, dataset, t_values, device,
                                psi_smooth_eval):
    """
    모델 출력의 circular-mean phase 를 추출하고, 동일 t 위치의 smooth 목표값을
    반환한다. isolation 실험의 extract_phases 와 동일한 구조지만 목표를
    {0, π} 대신 부드러운 값으로 사용한다.
    """
    feats = dataset.get_features_at_t(t_values).to(device)
    model.eval()
    with torch.enable_grad():
        X = feats.clone().requires_grad_(True)
        out = model(X)
    Z = out["Z_out"].detach().cpu().numpy()
    angles = np.angle(Z)
    while angles.ndim > 2:
        angles = angles.reshape(angles.shape[0], -1)
    if angles.ndim == 1:
        phi_circ = angles
    else:
        phi_circ = circular_mean(angles, axis=-1)

    t_arr = dataset.t.cpu().numpy()
    phi_true = []
    for t in t_values:
        idx = int(np.argmin(np.abs(t_arr - t)))
        phi_true.append(float(psi_smooth_eval[idx]))
    return phi_circ, np.array(phi_true)


def main():
    PrecisionManager.setup_precision()
    device = torch.device("cpu")

    print("=" * 70, flush=True)
    print("  Smooth-Target 실험 — Task #17", flush=True)
    print(f"  σ sweep={SIGMAS}, seeds={SEEDS}, epochs={EPOCHS}", flush=True)
    print("=" * 70, flush=True)

    zeros = load_zeros()
    nonzeros = sample_nonzero_t(zeros)
    print(f"  영점 {len(zeros)}개, 비영점 {len(nonzeros)}개", flush=True)

    cache_data = torch.load(
        os.path.join(OVERNIGHT_DIR, "xi_cache_t100-200_n500.pt"),
        weights_only=True
    )
    t_all = cache_data["t"].numpy()
    xi_real_all = cache_data["xi_real"].numpy()

    results = {}
    for sigma in SIGMAS:
        print(f"\n{'='*70}\n  σ = {sigma}\n{'='*70}", flush=True)
        psi_smooth = build_smooth_psi(t_all, zeros, sigma, xi_real_all)
        print(f"  smooth target range: [{psi_smooth.min():.3f}, {psi_smooth.max():.3f}]",
              flush=True)

        dataset = XiFeatureDatasetSmooth(
            cache_data, in_features=IN_FEATURES,
            psi_smooth_full=psi_smooth,
        )
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

        sigma_key = f"sigma_{sigma}"
        results[sigma_key] = {}
        for seed in SEEDS:
            label = f"smooth,σ={sigma},s={seed}"
            print(f"\n  [훈련] {label}", flush=True)
            model, best_val, train_time = train_smooth(
                train_loader, val_loader, HIDDEN, seed, device, label
            )
            print(f"    완료: best_val={best_val:.6f}, time={train_time:.1f}s",
                  flush=True)

            phi_z, phi_z_true = extract_phases_with_target(
                model, dataset, zeros, device, psi_smooth
            )
            phi_nz, phi_nz_true = extract_phases_with_target(
                model, dataset, nonzeros, device, psi_smooth
            )

            mse_z = phase_distance(phi_z, phi_z_true) ** 2
            mse_nz = phase_distance(phi_nz, phi_nz_true) ** 2

            # isolation 실험과의 비교를 위해 arg ξ ({0,π}) 대비 MSE 도 기록
            xi_phase_ref = np.where(xi_real_all >= 0, 0.0, math.pi)
            phi_z_discrete = []
            for t in zeros:
                idx = int(np.argmin(np.abs(t_all - t)))
                phi_z_discrete.append(xi_phase_ref[idx])
            phi_nz_discrete = []
            for t in nonzeros:
                idx = int(np.argmin(np.abs(t_all - t)))
                phi_nz_discrete.append(xi_phase_ref[idx])
            mse_z_disc = phase_distance(phi_z, np.array(phi_z_discrete)) ** 2
            mse_nz_disc = phase_distance(phi_nz, np.array(phi_nz_discrete)) ** 2

            results[sigma_key][seed] = {
                "mse_z_smooth": mse_z.tolist(),
                "mse_nz_smooth": mse_nz.tolist(),
                "mse_z_discrete": mse_z_disc.tolist(),
                "mse_nz_discrete": mse_nz_disc.tolist(),
                "best_val_loss": float(best_val),
                "train_time": float(train_time),
            }

    out_data = {
        "config": {
            "task": "Task #17 — smooth zero-crossing target",
            "psi_formula": "Ψ_smooth = (π/2)(1 - sign · ∏ tanh((t-z_k)/σ))",
            "sigmas": SIGMAS,
            "lambda_res": 0.0, "lambda_curv": 0.0,
            "lambda_tgt": 1.0, "lambda_pqo": 0.0,
            "hidden": HIDDEN, "seeds": SEEDS, "epochs": EPOCHS,
            "t_min": T_MIN, "t_max": T_MAX,
            "n_zeros": len(zeros), "n_nonzeros": len(nonzeros),
        },
        "zeros": zeros, "nonzeros": nonzeros,
        "results": {
            sk: {str(k): v for k, v in seed_res.items()}
            for sk, seed_res in results.items()
        },
    }
    with open(JSON_OUT, "w") as f:
        json.dump(out_data, f, indent=2, ensure_ascii=False)

    # ===== 통계 리포트 =====
    lines = []
    lines.append("=" * 70)
    lines.append("  Smooth-Target 실험 결과 (Task #17)")
    lines.append(f"  σ ∈ {SIGMAS}, seeds={SEEDS}, epochs={EPOCHS}")
    lines.append("  Ψ_smooth = (π/2)(1 - sign · ∏_k tanh((t-z_k)/σ))")
    lines.append("=" * 70)

    for sigma in SIGMAS:
        sk = f"sigma_{sigma}"
        lines.append("")
        lines.append("=" * 70)
        lines.append(f"  σ = {sigma}")
        lines.append("=" * 70)
        for metric_tag, key_z, key_nz in [
            ("smooth target", "mse_z_smooth", "mse_nz_smooth"),
            ("discrete {0,π} target (isolation 비교용)",
             "mse_z_discrete", "mse_nz_discrete"),
        ]:
            lines.append(f"  --- {metric_tag} ---")
            lines.append(f"  {'seed':>6s} {'영점 MSE':>14s} {'비영점 MSE':>14s} "
                         f"{'비율':>10s} {'p':>12s}")
            for seed in SEEDS:
                mz = np.array(results[sk][seed][key_z])
                mn = np.array(results[sk][seed][key_nz])
                ratio = mz.mean() / max(mn.mean(), 1e-12)
                try:
                    _, p_val = ttest_ind(mz, mn, equal_var=False)
                except Exception:
                    p_val = float("nan")
                lines.append(f"  {seed:>6d} {mz.mean():>14.6f} {mn.mean():>14.6f} "
                             f"{ratio:>10.4f} {p_val:>12.4e}")
            all_mz = np.concatenate([np.array(results[sk][s][key_z]) for s in SEEDS])
            all_mn = np.concatenate([np.array(results[sk][s][key_nz]) for s in SEEDS])
            ratio = all_mz.mean() / max(all_mn.mean(), 1e-12)
            try:
                _, p = ttest_ind(all_mz, all_mn, equal_var=False)
            except Exception:
                p = float("nan")
            lines.append(f"  pool: 영점={all_mz.mean():.6f} 비영점={all_mn.mean():.6f} "
                         f"비율={ratio:.4f} p={p:.4e}")

    lines.append("")
    lines.append("=" * 70)
    lines.append("  판정 및 해석")
    lines.append("=" * 70)
    iso_ref = {"z": 2.1708, "nz": 0.1195}
    lines.append(f"  isolation baseline (arg ξ {{0,π}} target, 참고):")
    lines.append(f"    영점 MSE={iso_ref['z']:.4f}, 비영점 MSE={iso_ref['nz']:.4f}, "
                 f"비율={iso_ref['z']/iso_ref['nz']:.2f}")
    for sigma in SIGMAS:
        sk = f"sigma_{sigma}"
        all_mz = np.concatenate([np.array(results[sk][s]["mse_z_smooth"])
                                  for s in SEEDS])
        all_mn = np.concatenate([np.array(results[sk][s]["mse_nz_smooth"])
                                  for s in SEEDS])
        all_mz_d = np.concatenate([np.array(results[sk][s]["mse_z_discrete"])
                                    for s in SEEDS])
        lines.append(f"  σ={sigma}: 영점 MSE(smooth)={all_mz.mean():.4f}, "
                     f"영점 MSE(discrete)={all_mz_d.mean():.4f}, "
                     f"비영점 MSE(smooth)={all_mn.mean():.4f}")
        if all_mz.mean() < 2 * all_mn.mean():
            lines.append(f"    → σ={sigma} 에서 영점 MSE 가 비영점 수준으로 하락. "
                         f"표현력 병목 해소.")
        else:
            lines.append(f"    → σ={sigma} 에서도 영점/비영점 gap 잔존.")

    output = "\n".join(lines)
    print("\n" + output, flush=True)
    with open(TXT_OUT, "w") as f:
        f.write(output)
    print(f"\n저장: {JSON_OUT}", flush=True)
    print(f"저장: {TXT_OUT}", flush=True)


if __name__ == "__main__":
    main()
