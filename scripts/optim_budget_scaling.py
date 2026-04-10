#!/usr/bin/env python3
"""
optim_budget_scaling.py — 3단계. H_optimization 분리.

B2 실험 결과: in_features=256 은 params 많아서 ep=100 으로 수렴 못함
(best_val=0.23, ep100 발산 11.9). 128 도 ep100 에서 발산 조짐 (0.155).

가설: 군집 어려움은 "모델 크기 vs 학습 예산" 상호작용.
  → ep↑ 로 128, 256 이 수렴하면 군집이 풀림 → 본질적 어려움 아님
  → 수렴해도 여전히 군집 남음 → 구조적 병목

실험:
  in_features ∈ {128, 256}, epochs=300, seed=42, hidden=64, n=500
  best-val 기준 checkpoint 저장, 매 50ep 마다 |F₂| 중간 관측
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
OUT_JSON = ANALYSIS / "optim_budget_scaling.json"
OUT_TXT = ANALYSIS / "optim_budget_scaling.txt"

T_MIN, T_MAX = 100, 200
HIDDEN = 64
EPOCHS = 300
SEED = 42
NUM_POINTS = 500
IN_FEATURES_LIST = [128, 256]
HARD_TARGETS = [173.4115, 107.1686]
CLUSTER_RANGE = (165.0, 173.5)
CHECKPOINT_EVERY = 50


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


def train_long(train_loader, val_loader, zero_features, in_features, seed,
               zeros_arr, cluster_mask, label):
    torch.manual_seed(seed); np.random.seed(seed)
    model = MasterResonantNetwork(
        in_features=in_features, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    n_params = sum(p.numel() for p in model.parameters())
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1,
        lambda_tgt=1.0, lambda_pqo=0.5, pqo_mode="cos2",
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    best_val = float("inf"); best_state = None
    best_epoch = 0
    t0 = time.time()
    checkpoints = []  # (epoch, val, |F2|_global, |F2|_cluster, |F2|_other)

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
            best_val = vl_m
            best_state = copy.deepcopy(model.state_dict())
            best_epoch = epoch

        if epoch % CHECKPOINT_EVERY == 0 or epoch == EPOCHS:
            # best-val 스냅샷 기준 F2 측정
            snap = copy.deepcopy(model.state_dict())
            model.load_state_dict(best_state)
            f2 = eval_F2(model, zero_features)
            ab = np.abs(f2)
            g = float(ab.mean()); cl = float(ab[cluster_mask].mean()); ot = float(ab[~cluster_mask].mean())
            checkpoints.append({
                "epoch": epoch, "best_val_so_far": best_val, "best_epoch": best_epoch,
                "F2_global": g, "F2_cluster": cl, "F2_other": ot, "ratio": cl/ot if ot > 0 else None,
                "F2_at_173": float(abs(f2[int(np.argmin(np.abs(zeros_arr - 173.4115)))])),
            })
            print(f"    {label} ep{epoch:03d}: val={vl_m:.5f} best={best_val:.5f}@{best_epoch}  "
                  f"|F₂|g={g:.4f} c={cl:.4f} o={ot:.4f}", flush=True)
            model.load_state_dict(snap)

    if best_state:
        model.load_state_dict(best_state)
    return model, best_val, best_epoch, time.time() - t0, n_params, checkpoints


def main():
    PrecisionManager.setup_precision()
    lines = []
    def log(s=""):
        print(s, flush=True); lines.append(s)

    log("=" * 72)
    log("  optim_budget_scaling — 3단계: ep=300 확장")
    log("=" * 72)
    log(f"  t∈[{T_MIN},{T_MAX}], hidden={HIDDEN}, seed={SEED}, ep={EPOCHS}, n={NUM_POINTS}")
    log(f"  in_features: {IN_FEATURES_LIST}   (ep=100 baseline 은 basis_bandwidth_scaling 에서 가져옴)")

    with open(OUTPUT_DIR / "result_t100-200.json") as f:
        rd = json.load(f)
    zeros_list = rd["zeros_list"]
    zeros_arr = np.array(zeros_list)
    cluster_mask = (zeros_arr >= CLUSTER_RANGE[0]) & (zeros_arr <= CLUSTER_RANGE[1])
    log(f"  영점 수: {len(zeros_list)}  군집: {int(cluster_mask.sum())}")

    # ep=100 결과 로드
    try:
        with open(ANALYSIS / "basis_bandwidth_scaling.json") as f:
            prev = json.load(f)
        log("\n  ep=100 baseline (basis_bandwidth_scaling):")
        log(f"    {'in_feat':>8s} {'global':>10s} {'cluster':>10s} {'other':>10s} {'ratio':>8s} {'val_best':>10s}")
        for ifx in IN_FEATURES_LIST:
            s = prev["series"][str(ifx)]
            arr = np.abs(np.array(s["F2"]))
            cl = arr[cluster_mask].mean(); ot = arr[~cluster_mask].mean()
            log(f"    {ifx:>8d} {arr.mean():>10.4f} {cl:>10.4f} {ot:>10.4f} {cl/ot:>7.2f}× "
                f"{s.get('best_val', 0):>10.4f}")
    except Exception as e:
        log(f"  (baseline 로드 실패: {e})")

    series = {}
    t_start = time.time()
    for ifx in IN_FEATURES_LIST:
        log(f"\n  === in_features={ifx}  ep={EPOCHS} ===")
        train_loader, val_loader, dataset = get_xi_feature_dataloaders(
            in_features=ifx, batch_size=32,
            t_min=T_MIN, t_max=T_MAX, num_points=NUM_POINTS,
            cache_dir=str(OUTPUT_DIR), zeros_list=zeros_list,
        )
        zero_features = dataset.get_features_at_t(zeros_list)
        model, bv, be, dt, np_, chks = train_long(
            train_loader, val_loader, zero_features, ifx, SEED,
            zeros_arr, cluster_mask, f"if={ifx}"
        )
        f2 = eval_F2(model, zero_features)
        log(f"    params={np_:,}  best_val={bv:.5f}@ep{be}  time={dt:.1f}s")
        ab = np.abs(f2)
        log(f"    최종 |F₂| global={ab.mean():.4f} cluster={ab[cluster_mask].mean():.4f} "
            f"other={ab[~cluster_mask].mean():.4f}")
        log(f"    t=173.4115: |F₂|={float(abs(f2[int(np.argmin(np.abs(zeros_arr-173.4115)))])):.5f}")
        series[ifx] = {
            "in_features": ifx, "F2": f2.tolist(),
            "n_params": np_, "best_val": bv, "best_epoch": be,
            "train_time_sec": dt, "checkpoints": chks,
        }

    log(f"\n  총 학습 시간: {time.time() - t_start:.1f}s")

    # ------------------ 판정 ------------------
    log("\n" + "=" * 72)
    log("  수렴 궤적 & 판정")
    log("=" * 72)
    for ifx in IN_FEATURES_LIST:
        log(f"\n  in_features={ifx} checkpoint 궤적:")
        log(f"    {'ep':>5s} {'best_val':>10s} {'|F₂|_g':>9s} {'|F₂|_c':>9s} {'ratio':>7s} {'|F₂|@173':>10s}")
        for c in series[ifx]["checkpoints"]:
            log(f"    {c['epoch']:>5d} {c['best_val_so_far']:>10.5f} "
                f"{c['F2_global']:>9.4f} {c['F2_cluster']:>9.4f} "
                f"{(c['ratio'] or 0):>6.2f}× {c['F2_at_173']:>10.5f}")

    # baseline ep=100 vs 최종 ep=300 비교
    log("\n  ep=100 → ep=300 개선:")
    try:
        prev_series = prev["series"]
        for ifx in IN_FEATURES_LIST:
            arr100 = np.abs(np.array(prev_series[str(ifx)]["F2"]))
            arr300 = np.abs(np.array(series[ifx]["F2"]))
            g1 = arr100.mean(); g3 = arr300.mean()
            c1 = arr100[cluster_mask].mean(); c3 = arr300[cluster_mask].mean()
            log(f"    if={ifx}: global {g1:.4f} → {g3:.4f} ({(1-g3/g1)*100:+.1f}%)  "
                f"cluster {c1:.4f} → {c3:.4f} ({(1-c3/c1)*100:+.1f}%)")
    except Exception:
        pass

    out = {
        "epochs": EPOCHS, "in_features_list": IN_FEATURES_LIST,
        "cluster_range": list(CLUSTER_RANGE),
        "series": {str(k): v for k, v in series.items()},
        "zeros_list": zeros_list,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")


if __name__ == "__main__":
    main()
