#!/usr/bin/env python3
"""
capacity_scaling.py — N2b. hidden 크기 스케일링.

가설: 난이도는 시드 다양성이 아닌 모델 함수족의 대역폭 한계.
방법: hidden ∈ {32, 128} 로 학습 (기존 hidden=64 결과와 합쳐 3 포인트).
지표: t=173.41, 107.17 에서 |F₂| 가 hidden 증가와 함께 감소하는가?

단일 seed=42, epochs=100. 기존 ensemble_hard_zeros.json 의 seed 42 결과 재사용.
"""
import copy
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "5")
os.environ.setdefault("MKL_NUM_THREADS", "5")
import numpy as np
import torch
torch.set_num_threads(5)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from gdl.rdl.constants import R_CONST, PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.pipeline.xi_feature_dataset import get_xi_feature_dataloaders

OUTPUT_DIR = PROJECT_ROOT / "outputs" / "overnight"
ANALYSIS = PROJECT_ROOT / "outputs" / "analysis"
OUT_JSON = ANALYSIS / "capacity_scaling.json"
OUT_TXT = ANALYSIS / "capacity_scaling.txt"

T_MIN, T_MAX = 100, 200
EPOCHS = 100
SEED = 42
IN_FEATURES = 64
HIDDENS = [32, 128]      # 64 는 ensemble_hard_zeros 에서 가져옴
HARD_TARGETS = [173.4115, 107.1686]


def train_one(train_loader, val_loader, hidden, seed, label):
    torch.manual_seed(seed); np.random.seed(seed)
    model = MasterResonantNetwork(
        in_features=IN_FEATURES, hidden_features=hidden, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    n_params = sum(p.numel() for p in model.parameters())
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1,
        lambda_tgt=1.0, lambda_pqo=0.5, pqo_mode="cos2",
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    best_val = float("inf"); best_state = None
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
            best_val = vl_m
            best_state = copy.deepcopy(model.state_dict())
        if epoch % 25 == 0 or epoch == EPOCHS:
            print(f"    {label} ep{epoch:03d}: val={vl_m:.5f}", flush=True)

    if best_state:
        model.load_state_dict(best_state)
    return model, best_val, time.time() - t0, n_params


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
    log("  capacity_scaling — N2b  (threads=5)")
    log("=" * 72)
    log(f"  t∈[{T_MIN},{T_MAX}], seed={SEED}, ep={EPOCHS}, hiddens={HIDDENS}")

    with open(OUTPUT_DIR / "result_t100-200.json") as f:
        rd = json.load(f)
    zeros_list = rd["zeros_list"]
    log(f"  영점 수: {len(zeros_list)}")

    train_loader, val_loader, dataset = get_xi_feature_dataloaders(
        in_features=IN_FEATURES, batch_size=32,
        t_min=T_MIN, t_max=T_MAX, num_points=500,
        cache_dir=str(OUTPUT_DIR), zeros_list=zeros_list,
    )
    zero_features = dataset.get_features_at_t(zeros_list)

    # 기존 hidden=64 결과 로드
    series = {}
    try:
        with open(ANALYSIS / "ensemble_hard_zeros.json") as f:
            ens = json.load(f)
        f2_64 = np.array(ens["per_seed_F2"]["42"])
        series[64] = {
            "hidden": 64, "seed": 42, "F2": f2_64.tolist(),
            "n_params": None, "best_val": None, "source": "ensemble_hard_zeros.json",
        }
        log(f"  기존 hidden=64 seed=42 로드: |F₂| mean={np.abs(f2_64).mean():.4f}")
    except Exception as e:
        log(f"  (hidden=64 로드 실패: {e})")

    t_start = time.time()
    for h in HIDDENS:
        log(f"\n  === hidden={h} ===")
        model, bv, dt, np_ = train_one(train_loader, val_loader, h, SEED, f"h={h}")
        f2 = eval_F2(model, zero_features)
        log(f"    params={np_:,}  best_val={bv:.5f}  time={dt:.1f}s")
        log(f"    |F₂| mean={np.abs(f2).mean():.4f}  max={np.abs(f2).max():.4f}")
        series[h] = {
            "hidden": h, "seed": SEED, "F2": f2.tolist(),
            "n_params": np_, "best_val": bv, "train_time_sec": dt,
        }

    log(f"\n  총 학습 시간: {time.time() - t_start:.1f}s")

    # ------------------ 분석 ------------------
    log("\n" + "=" * 72)
    log("  난이도 focus — hidden 스케일링")
    log("=" * 72)

    hs = sorted(series.keys())
    for t_h in HARD_TARGETS:
        idx = int(np.argmin(np.abs(np.array(zeros_list) - t_h)))
        log(f"\n  t={zeros_list[idx]:.4f}:")
        log(f"    {'hidden':>8s} {'F₂':>12s} {'|F₂|':>12s}")
        for h in hs:
            f2v = series[h]["F2"][idx]
            log(f"    {h:>8d} {f2v:>+12.5f} {abs(f2v):>12.5f}")

    # 전역 요약
    log("\n  전역 |F₂| mean/max:")
    log(f"    {'hidden':>8s} {'mean':>12s} {'max':>12s} {'params':>10s}")
    for h in hs:
        arr = np.abs(np.array(series[h]["F2"]))
        np_s = f"{series[h]['n_params']:,}" if series[h]['n_params'] else "—"
        log(f"    {h:>8d} {arr.mean():>12.4f} {arr.max():>12.4f} {np_s:>10s}")

    out = {"hiddens": hs, "series": series,
           "zeros_list": zeros_list, "hard_targets": HARD_TARGETS}
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")


if __name__ == "__main__":
    main()
