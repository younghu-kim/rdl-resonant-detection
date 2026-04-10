#!/usr/bin/env python3
"""
ensemble_hard_zeros.py — 난이도 영점 앙상블 실험 (재정의된 N2).

N4 결과: Lehmer 압축/gap 요인은 전역 난이도를 설명 못함. 잔존 가설:
  "어려운 영점 (t=173.41, 107.17) = 손실 지형의 다중 국소 최소"
  → 시드마다 다른 국소 표현 → cstd 에서 관찰된 분산 (1.62 vs 평균 0.26)

검증 방법:
  1. 동일 아키텍처 / 동일 데이터로 K 개 시드 독립 학습
  2. 각 모델의 per-zero |F₂| 측정
  3. 시드 분산 σ_F2(t) vs |F₂|(t) 상관
  4. 앙상블 (feature-space 평균 Z_out) 의 F₂ vs 개별 F₂ 비교
  5. 어려운 영점에서 앙상블 이득이 쉬운 영점보다 큰가?

만약 다중 국소 최소 가설이 맞다면:
  - σ_F2(173.41) ≫ σ_F2(easy)
  - 앙상블 F₂(173.41) ≪ mean(개별 F₂(173.41))
  - 쉬운 영점은 앙상블 이득 미미
"""
import copy
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from gdl.rdl.constants import R_CONST, PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.pipeline.xi_feature_dataset import get_xi_feature_dataloaders

OUTPUT_DIR = PROJECT_ROOT / "outputs" / "overnight"
ANALYSIS = PROJECT_ROOT / "outputs" / "analysis"
OUT_JSON = ANALYSIS / "ensemble_hard_zeros.json"
OUT_TXT = ANALYSIS / "ensemble_hard_zeros.txt"

T_MIN, T_MAX = 100, 200
HIDDEN = 64
EPOCHS = 100       # 절반 (overnight 는 200) — 앙상블 효과 입증에 충분
IN_FEATURES = 64
SEEDS = [42, 123, 777, 1337, 2024]
HARD_TARGETS = [173.4115, 107.1686]


def train_one(train_loader, val_loader, seed, label):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = MasterResonantNetwork(
        in_features=IN_FEATURES, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1,
        lambda_tgt=1.0, lambda_pqo=0.5,
        pqo_mode="cos2",
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    best_val = float("inf")
    best_state = None
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


def eval_at_zeros(model, zero_features):
    """한 모델로 영점 위치의 F₂, Z_out, phi, L_G, psi 를 뽑아 온다."""
    model.eval()
    with torch.enable_grad():
        X = zero_features.clone().requires_grad_(True)
        out = model(X)

    Z_out = out["Z_out"].detach()
    phi = out["phi"].detach()
    psi = out["psi"].detach()
    L_G = out["L_G"].detach()

    phi_real = phi.to(dtype=PrecisionManager.REAL_DTYPE)
    rotation = torch.complex(torch.cos(phi_real), -torch.sin(phi_real))
    psi_c = psi.to(dtype=PrecisionManager.COMPLEX_DTYPE)
    residual = rotation * (L_G - psi_c)
    F2 = residual.imag.mean(dim=-1)
    return {
        "F2": F2.cpu().numpy(),           # [N_zero]
        "Z_out": Z_out.cpu().numpy(),     # [N_zero, ...]
        "phi": phi.cpu().numpy(),         # [N_zero, ...]
        "psi": psi.cpu().numpy(),
        "L_G_real": L_G.real.cpu().numpy(),
        "L_G_imag": L_G.imag.cpu().numpy(),
    }


def main():
    PrecisionManager.setup_precision()
    device = torch.device("cpu")

    lines = []
    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    log("=" * 72)
    log("  ensemble_hard_zeros — N2 재정의: 다중 국소 최소 가설 검증")
    log("=" * 72)
    log(f"  t∈[{T_MIN},{T_MAX}], hidden={HIDDEN}, epochs={EPOCHS}, seeds={SEEDS}")

    # 영점 목록: overnight result 로부터 재사용 (mpmath 재계산 회피)
    with open(OUTPUT_DIR / "result_t100-200.json") as f:
        rd = json.load(f)
    zeros_list = rd["zeros_list"]
    log(f"  영점 수: {len(zeros_list)}")

    # 데이터 로더
    train_loader, val_loader, dataset = get_xi_feature_dataloaders(
        in_features=IN_FEATURES, batch_size=32,
        t_min=T_MIN, t_max=T_MAX, num_points=500,
        cache_dir=str(OUTPUT_DIR),
        zeros_list=zeros_list,
    )
    zero_features = dataset.get_features_at_t(zeros_list).to(device)

    # K 개 시드 학습
    per_seed = {}
    t_start = time.time()
    for i, seed in enumerate(SEEDS):
        label = f"s={seed}"
        log(f"\n  [{i+1}/{len(SEEDS)}] 학습 시작: {label}")
        model, bv, dt = train_one(train_loader, val_loader, seed, label)
        log(f"    best val={bv:.5f}  ({dt:.1f}s)")
        ev = eval_at_zeros(model, zero_features)
        per_seed[seed] = ev
        log(f"    |F₂| mean={np.abs(ev['F2']).mean():.4f}  "
            f"max={np.abs(ev['F2']).max():.4f}")

    total_t = time.time() - t_start
    log(f"\n  총 학습 시간: {total_t:.1f}s")

    # ------------------------------------------------------------------
    # 시드 간 분석
    # ------------------------------------------------------------------
    log("\n" + "=" * 72)
    log("  분석 1 — 시드 간 |F₂| 분산")
    log("=" * 72)

    # [K, N_zero] 배열
    F2_stack = np.stack([per_seed[s]["F2"] for s in SEEDS], axis=0)  # [K, N]
    absF2 = np.abs(F2_stack)
    per_zero_mean = absF2.mean(axis=0)    # mean over seeds (arithmetic)
    per_zero_std = absF2.std(axis=0)

    # 시드 간 F2 부호/값 범위
    F2_std_signed = F2_stack.std(axis=0)

    # 앙상블: Z_out 공간에서 복소수 평균 → F₂ 재계산
    Z_stack = np.stack([per_seed[s]["Z_out"] for s in SEEDS], axis=0)  # [K, N, ...]
    # 각 모델의 L_G, psi, phi 는 다르므로 앙상블은 residual 공간에서 평균이 자연스러움
    # residual_k = rot_k * (L_G_k - psi_k), F2_k = mean(Im(residual_k))
    # 앙상블 F2 = mean_k F2_k (단순)
    F2_ensemble_simple = F2_stack.mean(axis=0)  # [N]
    abs_ens_simple = np.abs(F2_ensemble_simple)

    # 영점 인덱스 찾기
    hard_idx = {t: int(np.argmin(np.abs(np.array(zeros_list) - t)))
                for t in HARD_TARGETS}

    log(f"\n  영점별 요약 ({len(zeros_list)}개):")
    log(f"    |F₂| 시드 평균: mean={per_zero_mean.mean():.4f}  "
        f"max={per_zero_mean.max():.4f}  at t={zeros_list[int(per_zero_mean.argmax())]:.4f}")
    log(f"    |F₂| 시드 표준편차: mean={per_zero_std.mean():.4f}  "
        f"max={per_zero_std.max():.4f}  at t={zeros_list[int(per_zero_std.argmax())]:.4f}")
    log(f"    앙상블 (mean F₂) 의 |F₂|: mean={abs_ens_simple.mean():.4f}  "
        f"max={abs_ens_simple.max():.4f}")

    # 개별-평균 vs 앙상블 비율 (앙상블 이득)
    with np.errstate(divide='ignore', invalid='ignore'):
        gain = np.where(per_zero_mean > 1e-9,
                        abs_ens_simple / per_zero_mean, 1.0)

    log("\n  난이도 focus 영점:")
    log(f"  {'t':>10s} {'|F₂| 평균':>12s} {'|F₂| std':>12s} "
        f"{'앙상블':>10s} {'이득':>8s} {'seeds':>40s}")
    for t_h, idx in hard_idx.items():
        seeds_fmt = " ".join(f"{F2_stack[k, idx]:+.3f}" for k in range(len(SEEDS)))
        log(f"  {zeros_list[idx]:>10.4f} {per_zero_mean[idx]:>12.4f} "
            f"{per_zero_std[idx]:>12.4f} {abs_ens_simple[idx]:>10.4f} "
            f"{gain[idx]:>7.2f}× {seeds_fmt:>40s}")

    # 쉬운 영점 대조군 (per_zero_mean 하위 5)
    easy_idx_sorted = np.argsort(per_zero_mean)[:5]
    log("\n  쉬운 영점 대조군 (|F₂| 최소 5):")
    log(f"  {'t':>10s} {'|F₂| 평균':>12s} {'|F₂| std':>12s} "
        f"{'앙상블':>10s} {'이득':>8s}")
    for idx in easy_idx_sorted:
        log(f"  {zeros_list[idx]:>10.4f} {per_zero_mean[idx]:>12.4f} "
            f"{per_zero_std[idx]:>12.4f} {abs_ens_simple[idx]:>10.4f} "
            f"{gain[idx]:>7.2f}×")

    # 어려운 영점 상위 5 (per_zero_mean 기준)
    hard_idx_sorted = np.argsort(-per_zero_mean)[:5]
    log("\n  어려운 영점 (|F₂| 평균 상위 5):")
    log(f"  {'t':>10s} {'|F₂| 평균':>12s} {'|F₂| std':>12s} "
        f"{'앙상블':>10s} {'이득':>8s}")
    for idx in hard_idx_sorted:
        log(f"  {zeros_list[idx]:>10.4f} {per_zero_mean[idx]:>12.4f} "
            f"{per_zero_std[idx]:>12.4f} {abs_ens_simple[idx]:>10.4f} "
            f"{gain[idx]:>7.2f}×")

    # 분산 vs 난이도 상관
    corr_mean_std = float(np.corrcoef(per_zero_mean, per_zero_std)[0, 1])
    log(f"\n  corr(|F₂| 평균, |F₂| std) = {corr_mean_std:+.4f}")
    log(f"    (>0 이면 어려운 영점 = 시드 분산 큼 → 다중 국소 최소 가설 지지)")

    # 전역 앙상블 효과
    log(f"\n  전역 앙상블 이득:")
    log(f"    평균 이득 (쉬운): {gain[easy_idx_sorted].mean():.3f}×")
    log(f"    평균 이득 (어려운 top5): {gain[hard_idx_sorted].mean():.3f}×")

    # ------------------------------------------------------------------
    # 저장
    # ------------------------------------------------------------------
    out = {
        "config": {
            "t_range": [T_MIN, T_MAX], "hidden": HIDDEN,
            "epochs": EPOCHS, "seeds": SEEDS,
            "hard_targets": HARD_TARGETS,
        },
        "zeros_list": zeros_list,
        "per_seed_F2": {str(s): per_seed[s]["F2"].tolist() for s in SEEDS},
        "per_zero_absF2_mean": per_zero_mean.tolist(),
        "per_zero_absF2_std": per_zero_std.tolist(),
        "ensemble_F2": F2_ensemble_simple.tolist(),
        "ensemble_absF2": abs_ens_simple.tolist(),
        "ensemble_gain": gain.tolist(),
        "corr_mean_std": corr_mean_std,
        "hard_targets_detail": {
            str(t): {
                "idx": idx, "t_actual": zeros_list[idx],
                "per_seed_F2": [float(F2_stack[k, idx]) for k in range(len(SEEDS))],
                "abs_mean": float(per_zero_mean[idx]),
                "abs_std": float(per_zero_std[idx]),
                "ensemble_absF2": float(abs_ens_simple[idx]),
                "gain": float(gain[idx]),
            } for t, idx in hard_idx.items()
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")
    log(f"저장: {OUT_TXT}")


if __name__ == "__main__":
    main()
