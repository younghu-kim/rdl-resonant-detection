#!/usr/bin/env python3
"""
=============================================================================
Phase Isolation Experiment — supervised Ψ + ONLY target loss (2026-04-08)
=============================================================================
phase_supervised_experiment.py 가 음성으로 끝남.
판정 직후의 두 가지 의심:

  (H_a) 손실 항 간섭: λ_res=1, λ_curv=0.1 가 target term 을 압도 → 학습이
        arg ξ 를 따라가지 못함.
  (H_b) 측정 결함: extract_phase_and_target 에서 mean(angle(Z)) 를 쓰는데
        이는 ±π wrap 신호를 평균 0 으로 망가뜨림. circular mean 필요.

이 스크립트는 두 가지를 동시에 차단:
  • λ_res = λ_curv = λ_pqo = 0, λ_tgt = 1.0 (또는 더 큰 값으로 sweep)
  • 측정에서 circular mean 사용: atan2(Σ sin θ, Σ cos θ)
  • Ψ_target 은 여전히 arg ξ(½+it)
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
from scripts.phase_supervised_experiment import (
    XiFeatureDatasetSupervised, load_zeros, sample_nonzero_t, cyclotomic_label,
    phase_distance, L_PQO, IN_FEATURES, BATCH_SIZE, T_MIN, T_MAX, ATTRACTORS,
)

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs", "analysis")
OVERNIGHT_DIR = os.path.join(PROJECT_ROOT, "outputs", "overnight")
JSON_OUT = os.path.join(OUTPUT_DIR, "phase_isolation_experiment.json")
TXT_OUT = os.path.join(OUTPUT_DIR, "phase_isolation_experiment.txt")

HIDDEN = 32
SEEDS = [42, 123, 777]
EPOCHS = 200
LR = 0.001


def circular_mean(angles, axis=-1):
    """θ_mean = atan2(mean sin θ, mean cos θ). ±π wrap 안전."""
    s = np.sin(angles).mean(axis=axis)
    c = np.cos(angles).mean(axis=axis)
    return np.arctan2(s, c)


def train_isolated(train_loader, val_loader, hidden, seed, device, label):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = MasterResonantNetwork(
        in_features=IN_FEATURES, hidden_features=hidden, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    # 모든 다른 항 차단 — 오직 target phase matching 만 학습
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


def extract_phases(model, dataset, t_values, device):
    """원시 phase (channel 별) + circular mean + arithmetic mean 모두 반환."""
    feats = dataset.get_features_at_t(t_values).to(device)
    model.eval()
    with torch.enable_grad():
        X = feats.clone().requires_grad_(True)
        out = model(X)
    Z = out["Z_out"].detach().cpu().numpy()  # 복소
    angles = np.angle(Z)  # shape [B, ...]
    # 마지막 축을 채널로 보고 두 가지 mean 계산
    while angles.ndim > 2:
        angles = angles.reshape(angles.shape[0], -1)
    if angles.ndim == 1:
        phi_circ = angles
        phi_arith = angles
    else:
        phi_circ = circular_mean(angles, axis=-1)
        phi_arith = angles.mean(axis=-1)

    t_arr = dataset.t.cpu().numpy()
    phi_true = []
    for t in t_values:
        idx = int(np.argmin(np.abs(t_arr - t)))
        phi_true.append(float(dataset.xi_phase[idx]))
    return phi_circ, phi_arith, np.array(phi_true)


def main():
    PrecisionManager.setup_precision()
    device = torch.device("cpu")

    print("=" * 70, flush=True)
    print("  Isolation 실험 — λ_res=λ_curv=λ_pqo=0, λ_tgt=1.0", flush=True)
    print(f"  hidden={HIDDEN}, seeds={SEEDS}, epochs={EPOCHS}, t∈[{T_MIN},{T_MAX}]", flush=True)
    print("=" * 70, flush=True)

    zeros = load_zeros()
    nonzeros = sample_nonzero_t(zeros)
    print(f"  영점 {len(zeros)}개, 비영점 {len(nonzeros)}개", flush=True)

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
        label = f"iso,h={HIDDEN},s={seed}"
        print(f"\n  [훈련] {label}", flush=True)
        model, best_val, train_time = train_isolated(
            train_loader, val_loader, HIDDEN, seed, device, label
        )
        print(f"    완료: best_val={best_val:.6f}, time={train_time:.1f}s", flush=True)

        phi_z_c, phi_z_a, phi_z_true = extract_phases(model, dataset, zeros, device)
        phi_nz_c, phi_nz_a, phi_nz_true = extract_phases(model, dataset, nonzeros, device)

        mse_z_c = phase_distance(phi_z_c, phi_z_true) ** 2
        mse_nz_c = phase_distance(phi_nz_c, phi_nz_true) ** 2
        mse_z_a = phase_distance(phi_z_a, phi_z_true) ** 2
        mse_nz_a = phase_distance(phi_nz_a, phi_nz_true) ** 2

        k_z = cyclotomic_label(phi_z_c)
        k_nz = cyclotomic_label(phi_nz_c)

        results[seed] = {
            "phi_z_circ": phi_z_c.tolist(),
            "phi_z_true": phi_z_true.tolist(),
            "phi_nz_circ": phi_nz_c.tolist(),
            "phi_nz_true": phi_nz_true.tolist(),
            "mse_z_circ": mse_z_c.tolist(),
            "mse_nz_circ": mse_nz_c.tolist(),
            "mse_z_arith": mse_z_a.tolist(),
            "mse_nz_arith": mse_nz_a.tolist(),
            "k_z": k_z.tolist(),
            "k_nz": k_nz.tolist(),
            "best_val_loss": float(best_val),
            "train_time": float(train_time),
        }

    out_data = {
        "config": {
            "lambda_res": 0.0, "lambda_curv": 0.0,
            "lambda_tgt": 1.0, "lambda_pqo": 0.0,
            "psi_source": "arg xi(1/2+it)",
            "hidden": HIDDEN, "seeds": SEEDS, "epochs": EPOCHS,
            "t_min": T_MIN, "t_max": T_MAX,
            "n_zeros": len(zeros), "n_nonzeros": len(nonzeros),
        },
        "zeros": zeros, "nonzeros": nonzeros,
        "results": {str(k): v for k, v in results.items()},
    }
    with open(JSON_OUT, "w") as f:
        json.dump(out_data, f, indent=2, ensure_ascii=False)

    # ===== 통계 =====
    lines = []
    lines.append("=" * 70)
    lines.append("  Isolation 실험 결과 — only λ_tgt, supervised Ψ")
    lines.append(f"  hidden={HIDDEN}, seeds={SEEDS}, epochs={EPOCHS}")
    lines.append("=" * 70)

    for tag, key_z, key_nz in [
        ("circular mean", "mse_z_circ", "mse_nz_circ"),
        ("arithmetic mean", "mse_z_arith", "mse_nz_arith"),
    ]:
        lines.append("")
        lines.append("=" * 70)
        lines.append(f"  Phase MSE — {tag}")
        lines.append("=" * 70)
        lines.append(f"  {'seed':>6s} {'영점 MSE':>14s} {'비영점 MSE':>14s} "
                     f"{'비율':>10s} {'p':>12s}")
        for seed in SEEDS:
            mz = np.array(results[seed][key_z])
            mn = np.array(results[seed][key_nz])
            ratio = mz.mean() / max(mn.mean(), 1e-12)
            try:
                _, p_val = ttest_ind(mz, mn, equal_var=False)
            except Exception:
                p_val = float('nan')
            lines.append(f"  {seed:>6d} {mz.mean():>14.6f} {mn.mean():>14.6f} "
                         f"{ratio:>10.4f} {p_val:>12.4e}")
        all_mz = np.concatenate([np.array(results[s][key_z]) for s in SEEDS])
        all_mn = np.concatenate([np.array(results[s][key_nz]) for s in SEEDS])
        ratio = all_mz.mean() / max(all_mn.mean(), 1e-12)
        try:
            _, p = ttest_ind(all_mz, all_mn, equal_var=False)
        except Exception:
            p = float('nan')
        lines.append(f"  pool: 영점={all_mz.mean():.6f} 비영점={all_mn.mean():.6f} "
                     f"비율={ratio:.4f} p={p:.4e}")

    lines.append("")
    lines.append("=" * 70)
    lines.append("  Cyclotomic k(t) 측정 (circular mean)")
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

    # 판정
    lines.append("")
    lines.append("=" * 70)
    lines.append("  판정 & 진단")
    lines.append("=" * 70)
    all_mz_c = np.concatenate([np.array(results[s]["mse_z_circ"]) for s in SEEDS])
    all_mn_c = np.concatenate([np.array(results[s]["mse_nz_circ"]) for s in SEEDS])
    all_mz_a = np.concatenate([np.array(results[s]["mse_z_arith"]) for s in SEEDS])
    all_mn_a = np.concatenate([np.array(results[s]["mse_nz_arith"]) for s in SEEDS])
    pool_mean_c = (all_mz_c.mean() + all_mn_c.mean()) / 2
    pool_mean_a = (all_mz_a.mean() + all_mn_a.mean()) / 2
    lines.append(f"  pooled mean MSE (circular)  = {pool_mean_c:.6f}")
    lines.append(f"  pooled mean MSE (arithmetic)= {pool_mean_a:.6f}")
    if pool_mean_c < 0.5:
        lines.append("  → 모델이 ψ=arg ξ 를 잘 학습. 이전 음성 결과는 손실 항 간섭 (H_a 옳음).")
    elif pool_mean_c < pool_mean_a * 0.5:
        lines.append("  → circular mean 으로만 신호 회복. 측정 결함 (H_b 옳음).")
    else:
        lines.append("  → 두 가지 모두 차단해도 학습 안 됨. 모델 표현력/구조 결함 의심.")

    output = "\n".join(lines)
    print("\n" + output, flush=True)
    with open(TXT_OUT, "w") as f:
        f.write(output)
    print(f"\n저장: {JSON_OUT}", flush=True)
    print(f"저장: {TXT_OUT}", flush=True)


if __name__ == "__main__":
    main()
