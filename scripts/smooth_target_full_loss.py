#!/usr/bin/env python3
"""
=============================================================================
Smooth Target + Full Loss Pipeline 실험 (2026-04-10)
=============================================================================
배경:
  - smooth_zero_crossing.py (isolation): smooth 타겟 xi_real/max|xi_real| 으로
    영점/비영점 MSE 갭 18x -> 0.03x 달성 (강한 양성).
  - basis_bandwidth_scaling.py: 표준 TotalResonanceLoss (4항 전체)로 50/50 검출.

질문:
  isolation 에서 확인된 smooth 타겟의 이점이 FULL pipeline (4항 손실 전체 활성)
  에서도 F2 기반 영점 검출 성능을 개선하는가?

설계:
  (A) Baseline: 표준 TotalResonanceLoss (lambda_res=1, lambda_curv=0.1,
      lambda_tgt=1, lambda_pqo=0.5, pqo_mode="cos2") — 기존 파이프라인.
  (B) Smooth: 동일 TotalResonanceLoss 에서 lambda_tgt=0 (angular target 비활성)
      + 커스텀 MSE(phase_out/pi vs xi_real/max|xi_real|) 항 추가 (lambda_smooth=1.0).

접근법 (b): 기존 손실 구조 유지 + 타겟 채널만 교체.
  - L_res, L_curv, L_pqo 는 그대로 유지
  - L_tgt (atan2 angular loss) 를 0으로 끄고
  - MSE(normalized_phase_out, smooth_target) 을 대체 항으로 추가

메트릭:
  - |F2| global mean, cluster mean (t in [165, 173.5]), other mean, ratio
  - 영점 검출 수 (|F2| < threshold)
  - baseline vs smooth 비교
"""
import os
import sys
import json
import time
import copy
import math
import numpy as np
import torch
import torch.nn as nn

os.environ.setdefault("OMP_NUM_THREADS", "5")
os.environ.setdefault("MKL_NUM_THREADS", "5")
torch.set_num_threads(5)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_ROOT)

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.pipeline.xi_feature_dataset import (
    XiFeatureDataset,
    get_xi_feature_dataloaders,
)

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs", "overnight")
ANALYSIS_DIR = os.path.join(PROJECT_ROOT, "outputs", "analysis")
JSON_OUT = os.path.join(ANALYSIS_DIR, "smooth_target_full_loss.json")
TXT_OUT = os.path.join(ANALYSIS_DIR, "smooth_target_full_loss.txt")

# ── 실험 상수 ──
IN_FEATURES = 128
HIDDEN = 64
SEEDS = [42, 7, 123]
EPOCHS = 100
LR = 0.001
BATCH_SIZE = 32
T_MIN, T_MAX = 100.0, 200.0
NUM_POINTS = 500
CLUSTER_RANGE = (165.0, 173.5)
F2_THRESHOLD = 0.01


# ─────────────────────────────────────────────────────────────────────
# Smooth 타겟 데이터셋: XiFeatureDataset + smooth target 추가 반환
# ─────────────────────────────────────────────────────────────────────
class XiSmoothTargetDataset(XiFeatureDataset):
    """
    표준 XiFeatureDataset 을 상속하되, __getitem__ 에서
    (features, xi_target, smooth_target) 3-튜플을 반환.
    smooth_target = xi_real / max|xi_real|
    """

    def __init__(self, cache_data, in_features=128):
        super().__init__(cache_data, in_features=in_features)
        # smooth 타겟 사전 계산: xi_real / max|xi_real|
        max_abs = torch.max(torch.abs(self.xi_real)).item()
        max_abs = max(max_abs, 1e-12)
        self.smooth_target = self.xi_real / max_abs

    def __getitem__(self, idx):
        features = self._build_features(self.t[idx])
        xi_target = torch.stack([self.xi_real[idx], self.xi_imag[idx]])
        smooth_tgt = self.smooth_target[idx].unsqueeze(0)  # [1]
        return features, xi_target, smooth_tgt


# ─────────────────────────────────────────────────────────────────────
# F2 평가 (basis_bandwidth_scaling.py 와 동일)
# ─────────────────────────────────────────────────────────────────────
def eval_F2(model, zero_features):
    """모델 출력에서 F2 = Im{ e^{-i*phi} (L_G - psi) } 계산"""
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


# ─────────────────────────────────────────────────────────────────────
# 훈련: Baseline (표준 TotalResonanceLoss, 4항 전체)
# ─────────────────────────────────────────────────────────────────────
def train_baseline(train_loader, val_loader, seed, label):
    """표준 TotalResonanceLoss 로 학습 (lambda_tgt=1.0 활성)"""
    torch.manual_seed(seed)
    np.random.seed(seed)

    model = MasterResonantNetwork(
        in_features=IN_FEATURES, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1,
        lambda_tgt=1.0, lambda_pqo=0.5, pqo_mode="cos2",
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
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

        # 검증
        model.eval()
        va = 0.0
        nv = 0
        with torch.enable_grad():
            for X_batch, _ in val_loader:
                X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
                X_in.requires_grad_(True)
                o = model(X_in)
                vl, _ = loss_fn(**o)
                if not (torch.isnan(vl) or torch.isinf(vl)):
                    va += vl.item()
                    nv += 1
        vl_m = va / max(1, nv)
        if vl_m < best_val:
            best_val = vl_m
            best_state = copy.deepcopy(model.state_dict())
        if epoch % 25 == 0 or epoch == EPOCHS:
            print(f"    {label} ep{epoch:03d}: val={vl_m:.5f}", flush=True)

    if best_state:
        model.load_state_dict(best_state)
    return model, best_val, time.time() - t0


# ─────────────────────────────────────────────────────────────────────
# 훈련: Smooth (TotalResonanceLoss lambda_tgt=0 + MSE smooth target)
# ─────────────────────────────────────────────────────────────────────
def train_smooth(train_loader, val_loader, seed, label, lambda_smooth=1.0):
    """
    TotalResonanceLoss (lambda_tgt=0) + 커스텀 MSE smooth 타겟.
    train_loader 는 XiSmoothTargetDataset 기반 (3-튜플 반환).
    """
    torch.manual_seed(seed)
    np.random.seed(seed)

    model = MasterResonantNetwork(
        in_features=IN_FEATURES, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    # lambda_tgt=0: angular target loss 비활성
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1,
        lambda_tgt=0.0, lambda_pqo=0.5, pqo_mode="cos2",
    )
    mse_fn = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    best_val = float("inf")
    best_state = None
    t0 = time.time()

    for epoch in range(1, EPOCHS + 1):
        model.train()
        for X_batch, _, smooth_tgt in train_loader:
            X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            tgt = smooth_tgt.to(dtype=PrecisionManager.REAL_DTYPE)  # [B, 1]

            optimizer.zero_grad(set_to_none=True)
            outputs = model(X_in)

            # (1) TotalResonanceLoss (lambda_tgt=0 이므로 L_res + L_curv + L_pqo 만)
            total_loss, loss_dict = loss_fn(**outputs)

            # (2) Smooth MSE: phase_out 채널 평균 / pi -> [-1,1] vs smooth target
            phase_out = outputs["phase_out"]  # [B, out_features]
            pred = phase_out.mean(dim=-1, keepdim=True) / math.pi  # [B, 1]
            smooth_loss = mse_fn(pred, tgt)

            # (3) 합산
            combined_loss = total_loss + lambda_smooth * smooth_loss

            if torch.isnan(combined_loss) or torch.isinf(combined_loss):
                continue
            combined_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()

        # 검증
        model.eval()
        va = 0.0
        nv = 0
        with torch.enable_grad():
            for X_batch, _, smooth_tgt in val_loader:
                X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
                X_in.requires_grad_(True)
                tgt = smooth_tgt.to(dtype=PrecisionManager.REAL_DTYPE)
                o = model(X_in)
                vl, _ = loss_fn(**o)
                pred = o["phase_out"].mean(dim=-1, keepdim=True) / math.pi
                sl = mse_fn(pred, tgt)
                combined = vl + lambda_smooth * sl
                if not (torch.isnan(combined) or torch.isinf(combined)):
                    va += combined.item()
                    nv += 1
        vl_m = va / max(1, nv)
        if vl_m < best_val:
            best_val = vl_m
            best_state = copy.deepcopy(model.state_dict())
        if epoch % 25 == 0 or epoch == EPOCHS:
            print(f"    {label} ep{epoch:03d}: val={vl_m:.5f}", flush=True)

    if best_state:
        model.load_state_dict(best_state)
    return model, best_val, time.time() - t0


# ─────────────────────────────────────────────────────────────────────
# F2 통계 분석
# ─────────────────────────────────────────────────────────────────────
def compute_f2_stats(f2_arr, zeros_arr, cluster_range, threshold):
    """F2 배열에서 global/cluster/other 평균 및 검출 수 계산"""
    abs_f2 = np.abs(f2_arr)
    cluster_mask = (zeros_arr >= cluster_range[0]) & (zeros_arr <= cluster_range[1])

    global_mean = float(abs_f2.mean())
    cluster_mean = float(abs_f2[cluster_mask].mean()) if cluster_mask.any() else float("nan")
    other_mean = float(abs_f2[~cluster_mask].mean()) if (~cluster_mask).any() else float("nan")
    ratio = cluster_mean / max(other_mean, 1e-12)
    detected = int((abs_f2 < threshold).sum())

    return {
        "global_mean": global_mean,
        "cluster_mean": cluster_mean,
        "other_mean": other_mean,
        "ratio": ratio,
        "detected": detected,
        "total": len(f2_arr),
        "threshold": threshold,
    }


# ─────────────────────────────────────────────────────────────────────
# 메인
# ─────────────────────────────────────────────────────────────────────
def main():
    PrecisionManager.setup_precision()
    device = torch.device("cpu")
    os.makedirs(ANALYSIS_DIR, exist_ok=True)

    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    log("=" * 72)
    log("  Smooth Target + Full Loss Pipeline 실험")
    log("=" * 72)
    log(f"  in_features={IN_FEATURES}, hidden={HIDDEN}, seeds={SEEDS}, ep={EPOCHS}")
    log(f"  t in [{T_MIN}, {T_MAX}], batch={BATCH_SIZE}, lr={LR}")
    log(f"  cluster range: {CLUSTER_RANGE}")
    log(f"  F2 threshold: {F2_THRESHOLD}")
    log("")
    log("  Baseline: TotalResonanceLoss(res=1, curv=0.1, tgt=1, pqo=0.5)")
    log("  Smooth:   TotalResonanceLoss(res=1, curv=0.1, tgt=0, pqo=0.5)")
    log("            + MSE(phase_out/pi, xi_real/max|xi_real|) * lambda_smooth=1.0")
    log("=" * 72)

    # ── 영점 목록 로드 ──
    result_path = os.path.join(OUTPUT_DIR, "result_t100-200.json")
    with open(result_path) as f:
        rd = json.load(f)
    zeros_list = rd["zeros_list"]
    zeros_arr = np.array(zeros_list)
    log(f"\n  영점 수: {len(zeros_list)}")

    # ── 캐시 로드 ──
    cache_path = os.path.join(OUTPUT_DIR, "xi_cache_t100-200_n500.pt")
    cache_data = torch.load(cache_path, weights_only=True)

    # ── 데이터셋 준비 ──
    # Baseline 용: 표준 XiFeatureDataset (2-튜플)
    # Smooth 용: XiSmoothTargetDataset (3-튜플)
    baseline_dataset = XiFeatureDataset(cache_data, in_features=IN_FEATURES)
    smooth_dataset = XiSmoothTargetDataset(cache_data, in_features=IN_FEATURES)

    # smooth 타겟 통계
    stgt = smooth_dataset.smooth_target.numpy()
    log(f"  smooth 타겟 범위: [{stgt.min():.4f}, {stgt.max():.4f}], "
        f"mean={stgt.mean():.4f}, std={stgt.std():.4f}")

    # 영점 근처 smooth 타겟 확인
    t_arr = smooth_dataset.t.cpu().numpy()
    zero_smooth_vals = []
    for z in zeros_list:
        idx = int(np.argmin(np.abs(t_arr - z)))
        zero_smooth_vals.append(float(smooth_dataset.smooth_target[idx]))
    log(f"  영점 근처 smooth 타겟 |평균|: {np.abs(zero_smooth_vals).mean():.6f} (0에 가까워야 함)")

    # 영점 특징 벡터 (F2 평가용)
    zero_features = baseline_dataset.get_features_at_t(zeros_list)

    # ── 시드별 실험 ──
    all_results = {"baseline": {}, "smooth": {}}
    t_total_start = time.time()

    for seed in SEEDS:
        log(f"\n{'─' * 72}")
        log(f"  seed={seed}")
        log(f"{'─' * 72}")

        # ── Baseline 훈련 ──
        log(f"\n  [Baseline] 훈련 시작...")
        # 데이터 분할 (seed 에 따라 일관적)
        val_size = int(len(baseline_dataset) * 0.2)
        train_size = len(baseline_dataset) - val_size
        bl_train, bl_val = torch.utils.data.random_split(
            baseline_dataset, [train_size, val_size],
            generator=torch.Generator().manual_seed(seed),
        )
        bl_train_loader = torch.utils.data.DataLoader(
            bl_train, batch_size=BATCH_SIZE, shuffle=True, drop_last=True,
        )
        bl_val_loader = torch.utils.data.DataLoader(
            bl_val, batch_size=BATCH_SIZE, shuffle=False, drop_last=False,
        )

        bl_model, bl_val_loss, bl_time = train_baseline(
            bl_train_loader, bl_val_loader, seed, f"BL,s={seed}",
        )
        bl_f2 = eval_F2(bl_model, zero_features)
        bl_stats = compute_f2_stats(bl_f2, zeros_arr, CLUSTER_RANGE, F2_THRESHOLD)

        log(f"    Baseline 완료: val={bl_val_loss:.5f}, time={bl_time:.1f}s")
        log(f"    |F2| global={bl_stats['global_mean']:.4f}, "
            f"cluster={bl_stats['cluster_mean']:.4f}, "
            f"other={bl_stats['other_mean']:.4f}, "
            f"ratio={bl_stats['ratio']:.4f}")
        log(f"    검출: {bl_stats['detected']}/{bl_stats['total']} "
            f"(threshold={F2_THRESHOLD})")

        all_results["baseline"][str(seed)] = {
            "F2": bl_f2.tolist(),
            "stats": bl_stats,
            "best_val_loss": float(bl_val_loss),
            "train_time": float(bl_time),
        }

        # ── Smooth 훈련 ──
        log(f"\n  [Smooth] 훈련 시작...")
        sm_train, sm_val = torch.utils.data.random_split(
            smooth_dataset, [train_size, val_size],
            generator=torch.Generator().manual_seed(seed),
        )
        sm_train_loader = torch.utils.data.DataLoader(
            sm_train, batch_size=BATCH_SIZE, shuffle=True, drop_last=True,
        )
        sm_val_loader = torch.utils.data.DataLoader(
            sm_val, batch_size=BATCH_SIZE, shuffle=False, drop_last=False,
        )

        sm_model, sm_val_loss, sm_time = train_smooth(
            sm_train_loader, sm_val_loader, seed, f"SM,s={seed}",
        )
        sm_f2 = eval_F2(sm_model, zero_features)
        sm_stats = compute_f2_stats(sm_f2, zeros_arr, CLUSTER_RANGE, F2_THRESHOLD)

        log(f"    Smooth 완료: val={sm_val_loss:.5f}, time={sm_time:.1f}s")
        log(f"    |F2| global={sm_stats['global_mean']:.4f}, "
            f"cluster={sm_stats['cluster_mean']:.4f}, "
            f"other={sm_stats['other_mean']:.4f}, "
            f"ratio={sm_stats['ratio']:.4f}")
        log(f"    검출: {sm_stats['detected']}/{sm_stats['total']} "
            f"(threshold={F2_THRESHOLD})")

        all_results["smooth"][str(seed)] = {
            "F2": sm_f2.tolist(),
            "stats": sm_stats,
            "best_val_loss": float(sm_val_loss),
            "train_time": float(sm_time),
        }

    total_time = time.time() - t_total_start
    log(f"\n  총 실행 시간: {total_time:.1f}s")

    # ─────────────────────────────────────────────────────────────────
    # 시드 풀 통합 분석
    # ─────────────────────────────────────────────────────────────────
    log("\n" + "=" * 72)
    log("  시드 풀 통합 결과")
    log("=" * 72)

    summary = {}
    for variant in ["baseline", "smooth"]:
        # 시드별 F2 를 모두 수집하여 평균
        all_f2_per_zero = []  # shape: [n_seeds, n_zeros]
        all_detected = 0
        all_total = 0
        for seed_key, res in all_results[variant].items():
            all_f2_per_zero.append(np.array(res["F2"]))
            all_detected += res["stats"]["detected"]
            all_total += res["stats"]["total"]

        # 시드 평균 F2
        f2_stack = np.stack(all_f2_per_zero, axis=0)  # [n_seeds, n_zeros]
        f2_mean_per_zero = np.abs(f2_stack).mean(axis=0)  # [n_zeros]
        cluster_mask = (zeros_arr >= CLUSTER_RANGE[0]) & (zeros_arr <= CLUSTER_RANGE[1])

        pool_global = float(f2_mean_per_zero.mean())
        pool_cluster = float(f2_mean_per_zero[cluster_mask].mean()) if cluster_mask.any() else float("nan")
        pool_other = float(f2_mean_per_zero[~cluster_mask].mean()) if (~cluster_mask).any() else float("nan")
        pool_ratio = pool_cluster / max(pool_other, 1e-12)

        summary[variant] = {
            "global_mean": pool_global,
            "cluster_mean": pool_cluster,
            "other_mean": pool_other,
            "ratio": pool_ratio,
            "detected_total": all_detected,
            "evaluated_total": all_total,
        }

        log(f"\n  [{variant.upper()}] 시드 풀 통합:")
        log(f"    |F2| global={pool_global:.4f}, cluster={pool_cluster:.4f}, "
            f"other={pool_other:.4f}")
        log(f"    cluster/other 비율: {pool_ratio:.4f}")
        log(f"    총 검출: {all_detected}/{all_total} "
            f"(시드x영점, threshold={F2_THRESHOLD})")

    # ─────────────────────────────────────────────────────────────────
    # 비교 테이블
    # ─────────────────────────────────────────────────────────────────
    log("\n" + "=" * 72)
    log("  비교 요약: Baseline vs Smooth")
    log("=" * 72)
    log(f"  {'variant':>12s} {'global':>10s} {'cluster':>10s} {'other':>10s} "
        f"{'ratio':>8s} {'detected':>10s}")

    for variant in ["baseline", "smooth"]:
        s = summary[variant]
        det_str = f"{s['detected_total']}/{s['evaluated_total']}"
        log(f"  {variant:>12s} {s['global_mean']:>10.4f} {s['cluster_mean']:>10.4f} "
            f"{s['other_mean']:>10.4f} {s['ratio']:>8.4f} {det_str:>10s}")

    # 개선 분석
    bl_s = summary["baseline"]
    sm_s = summary["smooth"]
    global_delta = (1.0 - sm_s["global_mean"] / max(bl_s["global_mean"], 1e-12)) * 100
    cluster_delta = (1.0 - sm_s["cluster_mean"] / max(bl_s["cluster_mean"], 1e-12)) * 100
    ratio_delta = sm_s["ratio"] - bl_s["ratio"]

    log(f"\n  global |F2| 변화: {global_delta:+.1f}% ({'개선' if global_delta > 0 else '악화'})")
    log(f"  cluster |F2| 변화: {cluster_delta:+.1f}% ({'개선' if cluster_delta > 0 else '악화'})")
    log(f"  ratio 변화: {bl_s['ratio']:.4f} -> {sm_s['ratio']:.4f} ({ratio_delta:+.4f})")
    log(f"  검출 변화: {bl_s['detected_total']} -> {sm_s['detected_total']}")

    # ── 판정 ──
    log("\n" + "=" * 72)
    log("  판정")
    log("=" * 72)

    # 기준: smooth 의 global |F2| 가 baseline 보다 낮으면 개선
    # 또는 검출 수가 더 많으면 개선
    if sm_s["global_mean"] < bl_s["global_mean"] * 0.8:
        verdict = (
            f"강한 양성: smooth 타겟이 full pipeline 에서도 |F2| 를 유의미하게 감소시킴.\n"
            f"    global |F2|: {bl_s['global_mean']:.4f} -> {sm_s['global_mean']:.4f} ({global_delta:+.1f}%)\n"
            f"    smooth 타겟 기반 full pipeline 이 기존 angular target 을 대체 가능."
        )
    elif sm_s["global_mean"] < bl_s["global_mean"] * 0.95:
        verdict = (
            f"약한 양성: smooth 타겟이 소폭 개선을 보임.\n"
            f"    global |F2|: {bl_s['global_mean']:.4f} -> {sm_s['global_mean']:.4f} ({global_delta:+.1f}%)\n"
            f"    추가 하이퍼파라미터 튜닝으로 개선 여지 있음."
        )
    elif sm_s["global_mean"] < bl_s["global_mean"] * 1.05:
        verdict = (
            f"중립: smooth 타겟과 baseline 이 유사한 성능.\n"
            f"    global |F2|: {bl_s['global_mean']:.4f} -> {sm_s['global_mean']:.4f} ({global_delta:+.1f}%)\n"
            f"    smooth 타겟이 angular 을 대체해도 성능 저하 없음 (대체 가능)."
        )
    else:
        verdict = (
            f"음성: smooth 타겟이 full pipeline 에서는 성능을 악화시킴.\n"
            f"    global |F2|: {bl_s['global_mean']:.4f} -> {sm_s['global_mean']:.4f} ({global_delta:+.1f}%)\n"
            f"    angular target loss (L_tgt) 가 L_res/L_pqo 와 시너지를 이루는 것으로 판단.\n"
            f"    smooth 타겟은 isolation 에서만 유효."
        )

    log(f"  {verdict}")

    # 다음 단계
    log("")
    log("-" * 72)
    log("  다음 단계:")
    log("  1. 강한 양성: lambda_smooth 스케일링 탐색 (0.5, 1.0, 2.0)")
    log("  2. 약한 양성/중립: L_tgt 와 smooth MSE 동시 활성 (하이브리드 손실)")
    log("  3. 음성: angular target 유지, smooth 타겟은 보조 진단 지표로만 활용")
    log("  4. 모든 경우: smooth 타겟의 cluster 내 상세 영점별 |F2| 패턴 분석")
    log("-" * 72)

    # ─────────────────────────────────────────────────────────────────
    # 시드별 상세 (군집 영점별 |F2|)
    # ─────────────────────────────────────────────────────────────────
    cluster_idx = np.where(
        (zeros_arr >= CLUSTER_RANGE[0]) & (zeros_arr <= CLUSTER_RANGE[1])
    )[0]

    if len(cluster_idx) > 0:
        log("\n" + "=" * 72)
        log("  군집 영점별 |F2| 상세 (시드 평균)")
        log("=" * 72)
        header = f"    {'t':>10s} {'baseline':>10s} {'smooth':>10s} {'delta%':>8s}"
        log(header)

        for ci in cluster_idx:
            bl_vals = [abs(all_results["baseline"][str(s)]["F2"][ci]) for s in SEEDS]
            sm_vals = [abs(all_results["smooth"][str(s)]["F2"][ci]) for s in SEEDS]
            bl_m = np.mean(bl_vals)
            sm_m = np.mean(sm_vals)
            delta = (1.0 - sm_m / max(bl_m, 1e-12)) * 100
            log(f"    {zeros_list[ci]:>10.4f} {bl_m:>10.4f} {sm_m:>10.4f} {delta:>+7.1f}%")

    # ─────────────────────────────────────────────────────────────────
    # JSON 저장
    # ─────────────────────────────────────────────────────────────────
    out_data = {
        "config": {
            "in_features": IN_FEATURES,
            "hidden": HIDDEN,
            "seeds": SEEDS,
            "epochs": EPOCHS,
            "lr": LR,
            "batch_size": BATCH_SIZE,
            "t_min": T_MIN,
            "t_max": T_MAX,
            "num_points": NUM_POINTS,
            "cluster_range": list(CLUSTER_RANGE),
            "f2_threshold": F2_THRESHOLD,
            "baseline_config": {
                "lambda_res": 1.0, "lambda_curv": 0.1,
                "lambda_tgt": 1.0, "lambda_pqo": 0.5,
                "pqo_mode": "cos2",
            },
            "smooth_config": {
                "lambda_res": 1.0, "lambda_curv": 0.1,
                "lambda_tgt": 0.0, "lambda_pqo": 0.5,
                "lambda_smooth": 1.0, "pqo_mode": "cos2",
                "smooth_target": "xi_real / max|xi_real|",
            },
        },
        "zeros_list": zeros_list,
        "results": all_results,
        "summary": summary,
        "verdict": verdict,
    }

    with open(JSON_OUT, "w") as f:
        json.dump(out_data, f, indent=2, ensure_ascii=False)

    with open(TXT_OUT, "w") as f:
        f.write("\n".join(lines))

    log(f"\n저장: {JSON_OUT}")
    log(f"저장: {TXT_OUT}")


if __name__ == "__main__":
    main()
