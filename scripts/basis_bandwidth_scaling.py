#!/usr/bin/env python3
"""
basis_bandwidth_scaling.py — B2 검증: in_features(=푸리에 기저 주파수 수) 스케일링.

가설 H3' (basis bandwidth):
  모델 대역폭 천장 ω_max = 2π · num_freqs / t_range
  군집 t∈[165,173] 의 0.8 간격 pair 를 풀려면 ω_req ≈ 3.93 rad/unit t
  num_freqs = (in_features-1)//2, t_range=100 에서
    in_features= 32 → num_freqs= 15 → ω_max=0.942  → 재앙
    in_features= 64 → num_freqs= 31 → ω_max=1.948  → 군집 불가 (기존)
    in_features=128 → num_freqs= 63 → ω_max=3.958  → 0.8 pair 경계 통과
    in_features=256 → num_freqs=127 → ω_max=7.980  → 여유

예측:
  (a) 64→128 에서 t=169 근처 |F₂| 급감
  (b) 128→256 에서 포화
  (c) 32 는 전체 재앙

hidden=64, seed=42, ep=100, n=500 고정.
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
OUT_JSON = ANALYSIS / "basis_bandwidth_scaling.json"
OUT_TXT = ANALYSIS / "basis_bandwidth_scaling.txt"

T_MIN, T_MAX = 100, 200
HIDDEN = 64
EPOCHS = 100
SEED = 42
NUM_POINTS = 500
IN_FEATURES_LIST = [32, 128, 256]   # 64 는 ensemble_hard_zeros 에서 재사용
HARD_TARGETS = [173.4115, 107.1686]
CLUSTER_RANGE = (165.0, 173.5)


def train_one(train_loader, val_loader, in_features, seed, label):
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


def omega_max(in_features, t_range):
    num_freqs = (in_features - 1) // 2
    return 2.0 * math.pi * num_freqs / t_range, num_freqs


def main():
    PrecisionManager.setup_precision()
    lines = []
    def log(s=""):
        print(s, flush=True); lines.append(s)

    log("=" * 72)
    log("  basis_bandwidth_scaling — B2 검증 (in_features 스케일링)")
    log("=" * 72)
    log(f"  t∈[{T_MIN},{T_MAX}], hidden={HIDDEN}, seed={SEED}, ep={EPOCHS}, n={NUM_POINTS}")
    log(f"  in_features 목록: {[32, 64, 128, 256]}  (64 는 재사용)")

    t_range = float(T_MAX - T_MIN)
    log("\n  이론 대역폭 (t_range=100):")
    log(f"    {'in_feat':>8s} {'num_freqs':>10s} {'ω_max':>10s} {'min_period':>12s}")
    for ifx in [32, 64, 128, 256]:
        om, nf = omega_max(ifx, t_range)
        log(f"    {ifx:>8d} {nf:>10d} {om:>10.4f} {2*math.pi/om:>12.4f}")
    log("\n  군집 임계: 0.8 간격 pair → 요구 ω ≈ 3.927 rad/unit t")
    log("    in_features=64  ω_max=1.948 → 불가 (2.0× 부족)")
    log("    in_features=128 ω_max=3.958 → 경계 통과")
    log("    in_features=256 ω_max=7.980 → 여유")

    with open(OUTPUT_DIR / "result_t100-200.json") as f:
        rd = json.load(f)
    zeros_list = rd["zeros_list"]
    log(f"\n  영점 수: {len(zeros_list)}")

    zeros_arr = np.array(zeros_list)
    cluster_mask = (zeros_arr >= CLUSTER_RANGE[0]) & (zeros_arr <= CLUSTER_RANGE[1])

    # 기존 hidden=64, in_features=64 seed=42 결과 로드
    series = {}
    try:
        with open(ANALYSIS / "ensemble_hard_zeros.json") as f:
            ens = json.load(f)
        f2_64 = np.array(ens["per_seed_F2"]["42"])
        series[64] = {
            "in_features": 64, "F2": f2_64.tolist(),
            "source": "ensemble_hard_zeros.json",
        }
        log(f"  기존 in_features=64 로드: |F₂| mean={np.abs(f2_64).mean():.4f}")
    except Exception as e:
        log(f"  (in_features=64 로드 실패: {e})")

    t_start = time.time()
    for ifx in IN_FEATURES_LIST:
        log(f"\n  === in_features={ifx} ===")
        train_loader, val_loader, dataset = get_xi_feature_dataloaders(
            in_features=ifx, batch_size=32,
            t_min=T_MIN, t_max=T_MAX, num_points=NUM_POINTS,
            cache_dir=str(OUTPUT_DIR), zeros_list=zeros_list,
        )
        zero_features = dataset.get_features_at_t(zeros_list)
        model, bv, dt, np_ = train_one(train_loader, val_loader, ifx, SEED, f"if={ifx}")
        f2 = eval_F2(model, zero_features)
        log(f"    params={np_:,}  best_val={bv:.5f}  time={dt:.1f}s")
        log(f"    전역 |F₂| mean={np.abs(f2).mean():.4f}  max={np.abs(f2).max():.4f}")
        cl = np.abs(f2[cluster_mask]).mean()
        ot = np.abs(f2[~cluster_mask]).mean()
        log(f"    군집 |F₂| mean={cl:.4f}  대조군 mean={ot:.4f}  비율={cl/ot:.2f}×")
        series[ifx] = {
            "in_features": ifx, "F2": f2.tolist(),
            "n_params": np_, "best_val": bv, "train_time_sec": dt,
        }

    log(f"\n  총 학습 시간: {time.time() - t_start:.1f}s")

    # ------------------ 분석 ------------------
    log("\n" + "=" * 72)
    log("  결과 요약")
    log("=" * 72)

    ifs = sorted(series.keys())

    # focus 영점
    for t_h in HARD_TARGETS:
        idx = int(np.argmin(np.abs(zeros_arr - t_h)))
        log(f"\n  t={zeros_list[idx]:.4f}:")
        log(f"    {'in_feat':>8s} {'F₂':>12s} {'|F₂|':>12s}")
        for ifx in ifs:
            f2v = series[ifx]["F2"][idx]
            log(f"    {ifx:>8d} {f2v:>+12.5f} {abs(f2v):>12.5f}")

    # 군집 상세 — 0.8 gap pair 특히
    log("\n  군집 영점별 |F₂| (0.8 gap pair 주목: t=169.1, 169.9 근처):")
    cluster_idx = np.where(cluster_mask)[0]
    header = "    " + f"{'t':>10s}" + "".join(f"{f'if={ifx}':>10s}" for ifx in ifs)
    log(header)
    for idx in cluster_idx:
        row = f"    {zeros_list[idx]:>10.4f}"
        for ifx in ifs:
            row += f"{abs(series[ifx]['F2'][idx]):>10.4f}"
        log(row)

    # 전역/군집 요약
    log("\n  전역 & 군집 |F₂| mean:")
    log(f"    {'in_feat':>8s} {'ω_max':>8s} {'global':>10s} {'cluster':>10s} {'other':>10s} {'ratio':>8s}")
    for ifx in ifs:
        om, _ = omega_max(ifx, t_range)
        arr = np.abs(np.array(series[ifx]["F2"]))
        cl = arr[cluster_mask].mean()
        ot = arr[~cluster_mask].mean()
        log(f"    {ifx:>8d} {om:>8.3f} {arr.mean():>10.4f} {cl:>10.4f} {ot:>10.4f} {cl/ot:>7.2f}×")

    # 예측 검증
    log("\n  B2 예측 검증:")
    g64 = np.abs(np.array(series[64]["F2"]))[cluster_mask].mean()
    g128 = np.abs(np.array(series[128]["F2"]))[cluster_mask].mean() if 128 in series else None
    g256 = np.abs(np.array(series[256]["F2"]))[cluster_mask].mean() if 256 in series else None
    g32 = np.abs(np.array(series[32]["F2"]))[cluster_mask].mean() if 32 in series else None

    if g32 is not None:
        log(f"    (a) 32 재앙? cluster={g32:.4f}  (64 의 {g32/g64:.2f}×)")
    if g128 is not None:
        log(f"    (b) 64→128 군집 급감? {g64:.4f} → {g128:.4f}  ({(1-g128/g64)*100:+.1f}%)")
    if g256 is not None and g128 is not None:
        log(f"    (c) 128→256 포화? {g128:.4f} → {g256:.4f}  ({(1-g256/g128)*100:+.1f}%)")

    out = {
        "in_features_list": ifs,
        "series": series,
        "zeros_list": zeros_list,
        "cluster_range": list(CLUSTER_RANGE),
        "hard_targets": HARD_TARGETS,
        "theory": {
            str(ifx): {
                "num_freqs": (ifx - 1) // 2,
                "omega_max": omega_max(ifx, t_range)[0],
                "min_period": 2 * math.pi / omega_max(ifx, t_range)[0],
            } for ifx in [32, 64, 128, 256]
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")


if __name__ == "__main__":
    main()
