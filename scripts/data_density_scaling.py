#!/usr/bin/env python3
"""
data_density_scaling.py — N2c. 훈련 데이터 밀도 스케일링.

가설 구분:
  (A) 모델 함수족 대역폭 한계 → 샘플 밀도 늘려도 개선 없음
  (B) 샘플 부족으로 국소 구조 학습 못함 → 밀도 늘리면 개선

hidden=64 고정, num_points ∈ {500, 1500} 로 학습 비교.
단일 seed=42, epochs=100.
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
OUT_JSON = ANALYSIS / "data_density_scaling.json"
OUT_TXT = ANALYSIS / "data_density_scaling.txt"

T_MIN, T_MAX = 100, 200
HIDDEN = 64
EPOCHS = 100
SEED = 42
IN_FEATURES = 64
NUM_POINTS_LIST = [1500]   # 기존 500 은 ensemble_hard_zeros 에서 로드
HARD_TARGETS = [173.4115, 107.1686]


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
    return model, best_val, time.time() - t0


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
    log("  data_density_scaling — N2c  (threads=5)")
    log("=" * 72)
    log(f"  hidden={HIDDEN}, seed={SEED}, ep={EPOCHS}")
    log(f"  num_points: [500 (기존 재사용), 1500]")

    with open(OUTPUT_DIR / "result_t100-200.json") as f:
        rd = json.load(f)
    zeros_list = rd["zeros_list"]

    series = {}
    # 기존 500 로드
    try:
        with open(ANALYSIS / "ensemble_hard_zeros.json") as f:
            ens = json.load(f)
        f2_500 = np.array(ens["per_seed_F2"]["42"])
        series[500] = {
            "num_points": 500, "F2": f2_500.tolist(),
            "source": "ensemble_hard_zeros.json",
        }
        log(f"  기존 num_points=500 seed=42 로드: |F₂| mean={np.abs(f2_500).mean():.4f}")
    except Exception as e:
        log(f"  (500 로드 실패: {e})")

    t_start = time.time()
    for npts in NUM_POINTS_LIST:
        log(f"\n  === num_points={npts} ===")
        # 별도 cache 파일 이름으로 빠져나가도록 cache_dir 은 공유, 하지만 파일은
        # 자동으로 n{npts} 접미 붙어 저장됨 (xi_feature_dataset 내부)
        train_loader, val_loader, dataset = get_xi_feature_dataloaders(
            in_features=IN_FEATURES, batch_size=32,
            t_min=T_MIN, t_max=T_MAX, num_points=npts,
            cache_dir=str(OUTPUT_DIR), zeros_list=zeros_list,
        )
        zero_features = dataset.get_features_at_t(zeros_list)
        model, bv, dt = train_one(train_loader, val_loader, SEED, f"n={npts}")
        f2 = eval_F2(model, zero_features)
        log(f"    best_val={bv:.5f}  time={dt:.1f}s")
        log(f"    |F₂| mean={np.abs(f2).mean():.4f}  max={np.abs(f2).max():.4f}")
        series[npts] = {
            "num_points": npts, "F2": f2.tolist(),
            "best_val": bv, "train_time_sec": dt,
        }

    log(f"\n  총 학습 시간: {time.time() - t_start:.1f}s")

    # ------------------ 분석 ------------------
    log("\n  난이도 focus:")
    nps = sorted(series.keys())
    for t_h in HARD_TARGETS:
        idx = int(np.argmin(np.abs(np.array(zeros_list) - t_h)))
        log(f"\n  t={zeros_list[idx]:.4f}:")
        log(f"    {'num_points':>12s} {'F₂':>12s} {'|F₂|':>12s}")
        for n in nps:
            f2v = series[n]["F2"][idx]
            log(f"    {n:>12d} {f2v:>+12.5f} {abs(f2v):>12.5f}")

    log(f"\n  전역 |F₂| mean/max:")
    log(f"    {'num_points':>12s} {'mean':>12s} {'max':>12s}")
    for n in nps:
        arr = np.abs(np.array(series[n]["F2"]))
        log(f"    {n:>12d} {arr.mean():>12.4f} {arr.max():>12.4f}")

    out = {"num_points_list": nps, "series": series,
           "zeros_list": zeros_list, "hard_targets": HARD_TARGETS}
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")


if __name__ == "__main__":
    main()
