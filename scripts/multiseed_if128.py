#!/usr/bin/env python3
"""
multiseed_if128.py — in_features=128 multi-seed 검증.

단일 seed=42 결과 (basis_bandwidth_scaling):
  cluster ratio 0.90×, t=173.41 |F₂|=0.00014

5-seed 앙상블로 재현성 검증. ep=100, hidden=64.
"""
import copy
import json
import os
import sys
import time
import math
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "5")
os.environ.setdefault("MKL_NUM_THREADS", "5")
import numpy as np
import torch
torch.set_num_threads(5)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.pipeline.xi_feature_dataset import get_xi_feature_dataloaders

OUTPUT_DIR = PROJECT_ROOT / "outputs" / "overnight"
ANALYSIS = PROJECT_ROOT / "outputs" / "analysis"
OUT_JSON = ANALYSIS / "multiseed_if128.json"
OUT_TXT = ANALYSIS / "multiseed_if128.txt"

T_MIN, T_MAX = 100, 200
IN_FEATURES = 128
HIDDEN = 64
EPOCHS = 100
NUM_POINTS = 500
SEEDS = [42, 7, 123, 2026, 314]
CLUSTER_RANGE = (165.0, 173.5)


def train_one(train_loader, val_loader, seed, label):
    torch.manual_seed(seed); np.random.seed(seed)
    model = MasterResonantNetwork(
        in_features=IN_FEATURES, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1,
        lambda_tgt=1.0, lambda_pqo=0.5, pqo_mode="cos2",
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    best_val = float("inf"); best_state = None; best_ep = 0
    t0 = time.time()

    for epoch in range(1, EPOCHS + 1):
        model.train()
        for X_batch, _ in train_loader:
            X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
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
        va = 0.0; nv = 0
        with torch.enable_grad():
            for X_batch, _ in val_loader:
                X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
                X_in.requires_grad_(True)
                o = model(X_in)
                vl, _ = loss_fn(**o)
                if not (torch.isnan(vl) or torch.isinf(vl)):
                    va += vl.item(); nv += 1
        vl_m = va / max(1, nv)
        if vl_m < best_val:
            best_val = vl_m; best_state = copy.deepcopy(model.state_dict()); best_ep = epoch
        if epoch % 25 == 0 or epoch == EPOCHS:
            print(f"    {label} ep{epoch:03d}: val={vl_m:.5f} best={best_val:.5f}@{best_ep}", flush=True)

    if best_state:
        model.load_state_dict(best_state)
    return model, best_val, best_ep, time.time() - t0


def eval_F2(model, zero_features):
    model.eval()
    with torch.enable_grad():
        X = zero_features.clone().requires_grad_(True)
        out = model(X)
    phi = out["phi"].detach()
    psi = out["psi"].detach()
    L_G = out["L_G"].detach()
    phi_real = phi.to(dtype=PrecisionManager.REAL_DTYPE)
    rot = torch.complex(torch.cos(phi_real), -torch.sin(phi_real))
    psi_c = psi.to(dtype=PrecisionManager.COMPLEX_DTYPE)
    return (rot * (L_G - psi_c)).imag.mean(dim=-1).cpu().numpy()


def main():
    PrecisionManager.setup_precision()
    lines = []
    def log(s=""):
        print(s, flush=True); lines.append(s)

    log("=" * 72)
    log("  multiseed_if128 — in_features=128 5-seed 재현성 검증")
    log("=" * 72)
    log(f"  t∈[{T_MIN},{T_MAX}], if={IN_FEATURES}, hidden={HIDDEN}, ep={EPOCHS}, n={NUM_POINTS}")
    log(f"  seeds: {SEEDS}")

    with open(OUTPUT_DIR / "result_t100-200.json") as f:
        rd = json.load(f)
    zeros_list = rd["zeros_list"]
    zeros_arr = np.array(zeros_list)
    cluster_mask = (zeros_arr >= CLUSTER_RANGE[0]) & (zeros_arr <= CLUSTER_RANGE[1])
    n_cluster = int(cluster_mask.sum())
    log(f"  영점 수: {len(zeros_list)}  군집: {n_cluster}")

    train_loader, val_loader, dataset = get_xi_feature_dataloaders(
        in_features=IN_FEATURES, batch_size=32,
        t_min=T_MIN, t_max=T_MAX, num_points=NUM_POINTS,
        cache_dir=str(OUTPUT_DIR), zeros_list=zeros_list,
    )
    zero_features = dataset.get_features_at_t(zeros_list)

    per_seed = {}
    t_start = time.time()

    for seed in SEEDS:
        log(f"\n  === seed={seed} ===")
        model, bv, be, dt = train_one(train_loader, val_loader, seed, f"s={seed}")
        f2 = eval_F2(model, zero_features)
        ab = np.abs(f2)
        gl = float(ab.mean())
        cl = float(ab[cluster_mask].mean())
        ot = float(ab[~cluster_mask].mean())
        ratio = cl / ot if ot > 0 else None
        idx_173 = int(np.argmin(np.abs(zeros_arr - 173.4115)))
        f2_173 = float(abs(f2[idx_173]))
        log(f"    best_val={bv:.5f}@ep{be}  time={dt:.1f}s")
        log(f"    |F₂| global={gl:.4f} cluster={cl:.4f} other={ot:.4f} ratio={ratio:.2f}×")
        log(f"    t=173.41: |F₂|={f2_173:.5f}")
        per_seed[seed] = {
            "seed": seed, "best_val": bv, "best_epoch": be,
            "train_time": dt, "F2": f2.tolist(),
            "global": gl, "cluster": cl, "other": ot, "ratio": ratio,
            "F2_at_173": f2_173,
        }

    total_time = time.time() - t_start
    log(f"\n  총 학습 시간: {total_time:.1f}s")

    # ==================== 통계 ====================
    log("\n" + "=" * 72)
    log("  앙상블 통계")
    log("=" * 72)

    ratios = [per_seed[s]["ratio"] for s in SEEDS]
    globals_ = [per_seed[s]["global"] for s in SEEDS]
    clusters = [per_seed[s]["cluster"] for s in SEEDS]
    others = [per_seed[s]["other"] for s in SEEDS]
    f2_173s = [per_seed[s]["F2_at_173"] for s in SEEDS]

    log(f"\n  {'metric':>15s} {'mean':>10s} {'std':>10s} {'min':>10s} {'max':>10s}")
    for name, arr in [("ratio", ratios), ("global", globals_),
                      ("cluster", clusters), ("other", others),
                      ("|F₂|@173", f2_173s)]:
        a = np.array(arr)
        log(f"  {name:>15s} {a.mean():>10.4f} {a.std():>10.4f} {a.min():>10.4f} {a.max():>10.4f}")

    # 핵심 질문: ratio < 1 이 재현되는가?
    n_resolved = sum(1 for r in ratios if r < 1.2)
    log(f"\n  ratio < 1.2 (군집 해소): {n_resolved}/{len(SEEDS)} seeds")
    log(f"  ratio 평균: {np.mean(ratios):.3f} ± {np.std(ratios):.3f}")

    # if=64 baseline과 비교
    try:
        with open(ANALYSIS / "basis_bandwidth_scaling.json") as f:
            prev = json.load(f)
        f2_64 = np.abs(np.array(prev["series"]["64"]["F2"]))
        cl_64 = float(f2_64[cluster_mask].mean())
        ot_64 = float(f2_64[~cluster_mask].mean())
        ratio_64 = cl_64 / ot_64
        log(f"\n  if=64 baseline (seed=42): ratio={ratio_64:.2f}× cluster={cl_64:.4f}")
        log(f"  if=128 앙상블:            ratio={np.mean(ratios):.2f}× cluster={np.mean(clusters):.4f}")
        improvement = (1 - np.mean(clusters) / cl_64) * 100
        log(f"  군집 |F₂| 개선: {improvement:+.1f}%")
    except Exception as e:
        log(f"  (if=64 baseline 로드 실패: {e})")

    # 영점별 앙상블 평균 |F₂|
    all_f2 = np.array([np.abs(np.array(per_seed[s]["F2"])) for s in SEEDS])
    mean_f2 = all_f2.mean(axis=0)
    std_f2 = all_f2.std(axis=0)

    log("\n  군집 영점별 앙상블 |F₂| (mean ± std):")
    cluster_idx = np.where(cluster_mask)[0]
    log(f"    {'t':>10s} {'mean':>10s} {'std':>10s} {'cv':>8s}")
    for idx in cluster_idx:
        cv = std_f2[idx] / mean_f2[idx] if mean_f2[idx] > 0 else 0
        log(f"    {zeros_list[idx]:>10.4f} {mean_f2[idx]:>10.5f} {std_f2[idx]:>10.5f} {cv:>7.1%}")

    # 판정
    log("\n" + "=" * 72)
    if n_resolved >= 4:
        log("  판정: ✓ in_features=128 군집 해소 재현 확인 (≥4/5 seeds)")
    elif n_resolved >= 3:
        log("  판정: △ 대체로 재현 (3/5 seeds)")
    else:
        log("  판정: ✗ 재현 실패 — single-seed artifact 가능성")
    log("=" * 72)

    out = {
        "in_features": IN_FEATURES, "hidden": HIDDEN, "epochs": EPOCHS,
        "seeds": SEEDS, "cluster_range": list(CLUSTER_RANGE),
        "per_seed": {str(s): per_seed[s] for s in SEEDS},
        "ensemble_mean_F2": mean_f2.tolist(),
        "ensemble_std_F2": std_f2.tolist(),
        "zeros_list": zeros_list,
        "summary": {
            "ratio_mean": float(np.mean(ratios)),
            "ratio_std": float(np.std(ratios)),
            "cluster_mean": float(np.mean(clusters)),
            "n_resolved": n_resolved,
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")


if __name__ == "__main__":
    main()
