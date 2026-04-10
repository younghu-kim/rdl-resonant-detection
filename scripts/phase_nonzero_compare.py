#!/usr/bin/env python3
"""
=============================================================================
φ 비영점 비교 실험 — φ → π/2 가 영점 특이적 신호인가?
=============================================================================
귀무가설 H0: φ는 모든 t에서 π/2 근처 (모델이 π/2 상수 출력)
대립가설 H1: φ는 영점에서만 π/2로 수렴, 비영점에서는 다른 분포

설정:
  - hidden=32, 3 seeds (기존 f2_model_independence와 동일)
  - 구간: t∈[100,200]
  - 영점: 캐시된 50개
  - 비영점: 균일 200점, 어떤 영점에서도 ≥0.5 떨어진 곳

출력:
  - JSON: 영점/비영점 위상 배열
  - TXT: 평균/분산/KS 검정/Wilcoxon
"""
import os
import sys
import json
import time
import copy
import math
import numpy as np
import torch
from scipy.stats import ks_2samp, mannwhitneyu

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_ROOT)

from gdl.rdl.constants import R_CONST, PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.pipeline.xi_feature_dataset import XiFeatureDataset

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs", "analysis")
OVERNIGHT_DIR = os.path.join(PROJECT_ROOT, "outputs", "overnight")
JSON_OUT = os.path.join(OUTPUT_DIR, "phase_nonzero_compare.json")
TXT_OUT = os.path.join(OUTPUT_DIR, "phase_nonzero_compare.txt")

PI_HALF = math.pi / 2.0
HIDDEN = 32
SEEDS = [42, 123, 777]
EPOCHS = 200
LR = 0.001
IN_FEATURES = 64
BATCH_SIZE = 32
T_MIN, T_MAX = 100.0, 200.0
N_NONZERO = 200
NONZERO_GAP = 0.5  # 영점에서 최소 떨어져야 할 거리


def load_zeros():
    path = os.path.join(OVERNIGHT_DIR, "result_t100-200.json")
    with open(path) as f:
        data = json.load(f)
    return data["zeros_list"]


def sample_nonzero_t(zeros, n=N_NONZERO, seed=0):
    """t∈[T_MIN, T_MAX]에서 균일 후보 생성 후 영점 ±NONZERO_GAP 제외."""
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


def train_model(train_loader, val_loader, hidden, seed, device, label):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = MasterResonantNetwork(
        in_features=IN_FEATURES, hidden_features=hidden, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.5,
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
    Z_out = out["Z_out"]
    phase = torch.angle(Z_out).mean(dim=-1).detach().cpu().numpy()
    return phase


def main():
    PrecisionManager.setup_precision()
    device = torch.device("cpu")

    print("=" * 70)
    print("  φ 비영점 비교 실험")
    print("=" * 70)

    zeros = load_zeros()
    nonzeros = sample_nonzero_t(zeros)
    print(f"  영점 {len(zeros)}개, 비영점 {len(nonzeros)}개")

    # 데이터셋
    cache_path = os.path.join(OVERNIGHT_DIR, "xi_cache_t100-200_n500.pt")
    cache_data = torch.load(cache_path, weights_only=True)
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
        label = f"h={HIDDEN},s={seed}"
        print(f"\n  [훈련] {label}")
        model, best_val, train_time = train_model(
            train_loader, val_loader, HIDDEN, seed, device, label
        )
        print(f"    완료: best_val={best_val:.4f}, time={train_time:.1f}s")

        phi_zero = extract_phase(model, dataset, zeros, device)
        phi_nonzero = extract_phase(model, dataset, nonzeros, device)

        results[seed] = {
            "phi_zero": phi_zero.tolist(),
            "phi_nonzero": phi_nonzero.tolist(),
            "best_val_loss": float(best_val),
            "train_time": float(train_time),
        }

    # 저장 (즉시)
    out_data = {
        "config": {
            "hidden": HIDDEN, "seeds": SEEDS, "epochs": EPOCHS,
            "t_min": T_MIN, "t_max": T_MAX,
            "n_zeros": len(zeros), "n_nonzeros": len(nonzeros),
            "nonzero_gap": NONZERO_GAP,
        },
        "zeros": zeros,
        "nonzeros": nonzeros,
        "results": {str(k): v for k, v in results.items()},
    }
    with open(JSON_OUT, "w") as f:
        json.dump(out_data, f, indent=2, ensure_ascii=False)

    # 통계 분석
    lines = []
    lines.append("=" * 70)
    lines.append("  φ 비영점 비교 — 결과")
    lines.append(f"  hidden={HIDDEN}, seeds={SEEDS}, epochs={EPOCHS}")
    lines.append(f"  영점 {len(zeros)}개, 비영점 {len(nonzeros)}개 (gap≥{NONZERO_GAP})")
    lines.append(f"  π/2 = {PI_HALF:.10f}")
    lines.append("=" * 70)

    lines.append("")
    lines.append("=" * 70)
    lines.append("  1. 시드별 영점 vs 비영점 φ 통계")
    lines.append("=" * 70)
    lines.append(f"  {'seed':>6s} {'영점 φ평균':>14s} {'비영점 φ평균':>14s} "
                 f"{'영점 |φ-π/2|':>15s} {'비영점 |φ-π/2|':>16s}")
    for seed in SEEDS:
        z_phi = np.array(results[seed]["phi_zero"])
        nz_phi = np.array(results[seed]["phi_nonzero"])
        z_dev = np.abs(z_phi - PI_HALF).mean()
        nz_dev = np.abs(nz_phi - PI_HALF).mean()
        lines.append(f"  {seed:>6d} {z_phi.mean():>14.6f} {nz_phi.mean():>14.6f} "
                     f"{z_dev:>15.6f} {nz_dev:>16.6f}")

    lines.append("")
    lines.append("=" * 70)
    lines.append("  2. 분포 검정 (KS, Mann-Whitney U)")
    lines.append("=" * 70)
    for seed in SEEDS:
        z_phi = np.array(results[seed]["phi_zero"])
        nz_phi = np.array(results[seed]["phi_nonzero"])
        ks_stat, ks_p = ks_2samp(z_phi, nz_phi)
        mw_stat, mw_p = mannwhitneyu(z_phi, nz_phi, alternative="two-sided")
        lines.append(f"  s={seed}: KS={ks_stat:.4f} (p={ks_p:.2e}), "
                     f"MWU p={mw_p:.2e}")

    lines.append("")
    lines.append("=" * 70)
    lines.append("  3. 비영점 φ 분포 (시드 풀)")
    lines.append("=" * 70)
    all_nz = np.concatenate([np.array(results[s]["phi_nonzero"]) for s in SEEDS])
    all_z = np.concatenate([np.array(results[s]["phi_zero"]) for s in SEEDS])
    lines.append(f"  전체 비영점: 평균={all_nz.mean():.6f}, std={all_nz.std():.6f}")
    lines.append(f"  전체 영점:   평균={all_z.mean():.6f}, std={all_z.std():.6f}")
    lines.append(f"  비영점 |φ-π/2| 평균: {np.abs(all_nz - PI_HALF).mean():.6f}")
    lines.append(f"  영점 |φ-π/2| 평균:   {np.abs(all_z - PI_HALF).mean():.6f}")
    ratio = np.abs(all_nz - PI_HALF).mean() / max(1e-12, np.abs(all_z - PI_HALF).mean())
    lines.append(f"  비/영 비율: {ratio:.2f}배")

    lines.append("")
    lines.append("  히스토그램 (φ - π/2, 비영점):")
    bins = np.linspace(-math.pi, math.pi, 21)
    hist, _ = np.histogram(all_nz - PI_HALF, bins=bins)
    for i in range(len(hist)):
        lo, hi = bins[i], bins[i + 1]
        bar = "#" * int(40 * hist[i] / max(1, hist.max()))
        lines.append(f"  [{lo:+.3f},{hi:+.3f}) {hist[i]:>4d} {bar}")

    lines.append("")
    lines.append("=" * 70)
    lines.append("  4. 판정")
    lines.append("=" * 70)
    z_dev_mean = np.abs(all_z - PI_HALF).mean()
    nz_dev_mean = np.abs(all_nz - PI_HALF).mean()
    if nz_dev_mean > 5 * z_dev_mean:
        verdict = "→ φ→π/2는 영점 특이적 신호 (대립가설 채택)"
    elif nz_dev_mean > 2 * z_dev_mean:
        verdict = "→ 영점에서 더 강하지만 비영점도 부분 수렴"
    else:
        verdict = "→ φ는 모든 t에서 π/2 근처 (귀무가설, 신호 사소)"
    lines.append(f"  영점 잔차={z_dev_mean:.6f}, 비영점 잔차={nz_dev_mean:.6f}")
    lines.append(f"  {verdict}")

    output = "\n".join(lines)
    print("\n" + output)
    with open(TXT_OUT, "w") as f:
        f.write(output)
    print(f"\n저장: {JSON_OUT}")
    print(f"저장: {TXT_OUT}")


if __name__ == "__main__":
    main()
