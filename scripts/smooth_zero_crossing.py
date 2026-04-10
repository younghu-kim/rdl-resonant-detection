#!/usr/bin/env python3
"""
=============================================================================
Smooth Zero-Crossing Experiment (2026-04-10)
=============================================================================
phase_isolation_experiment.py 결과:
  - 비영점 MSE = 0.12 (학습 성공)
  - 영점 MSE   = 2.17 ≈ (π/2)² (학습 실패)
  - 비율       = 18×

진단: arg ξ(½+it) 는 임계선 위에서 {0, π} 계단 함수이므로, 영점에서
±π 불연속 점프가 발생. 연속 네트워크로는 표현 불가.

가설: ±π 점프가 진짜 병목이라면, 동일한 부호 정보를 담되 불연속이 없는
SMOOTH 타겟으로 교체하면 영점/비영점 MSE 갭이 극적으로 줄어들어야 한다.

두 가지 후보 타겟:
  (A) sign(ξ_real) · exp(-ξ_real²/σ²)  — 매끄러운 bump, 영점에서 0 교차
  (B) ξ_real / max|ξ_real|               — 정규화된 Z-유사 함수 (본질적 연속)

설정: isolation (λ_res=λ_curv=λ_pqo=0), 순수 MSE 회귀, circular mean 불필요
(타겟이 각도가 아니라 실수값이므로).

비교 기준선: isolation 실험의 arg ξ 타겟 결과 (phase_isolation_experiment.json)
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
from scipy.stats import ttest_ind

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_ROOT)

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.pipeline.xi_feature_dataset import XiFeatureDataset

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs", "analysis")
OVERNIGHT_DIR = os.path.join(PROJECT_ROOT, "outputs", "overnight")
JSON_OUT = os.path.join(OUTPUT_DIR, "smooth_zero_crossing.json")
TXT_OUT = os.path.join(OUTPUT_DIR, "smooth_zero_crossing.txt")

# 실험 구성 — bandwidth-clean 설정
IN_FEATURES = 128
HIDDEN = 64
SEEDS = [42, 7, 123]
EPOCHS = 100
LR = 0.001
BATCH_SIZE = 32
T_MIN, T_MAX = 100.0, 200.0
N_NONZERO = 200
NONZERO_GAP = 0.5


# ─────────────────────────────────────────────────────────────────────
# 데이터셋: smooth 타겟을 반환하는 래퍼
# ─────────────────────────────────────────────────────────────────────
class XiSmoothTargetDataset(XiFeatureDataset):
    """
    부모 클래스의 푸리에 기저 특징을 재사용하되,
    타겟으로 smooth 실수값(bump 또는 normalized xi_real)을 반환.

    target_mode:
      "bump"       → sign(ξ_real) · exp(-ξ_real²/σ²), σ = median(|ξ_real|)
      "normalized" → ξ_real / max(|ξ_real|)
    """
    def __init__(self, cache_data, in_features=128, target_mode="normalized"):
        super().__init__(cache_data, in_features=in_features)
        self.target_mode = target_mode

        # smooth 타겟 사전 계산
        xr = self.xi_real  # 텐서
        if target_mode == "bump":
            sigma = torch.median(torch.abs(xr)).item()
            sigma = max(sigma, 1e-12)
            self.smooth_target = torch.sign(xr) * torch.exp(-xr ** 2 / sigma ** 2)
        elif target_mode == "normalized":
            max_abs = torch.max(torch.abs(xr)).item()
            max_abs = max(max_abs, 1e-12)
            self.smooth_target = xr / max_abs
        else:
            raise ValueError(f"알 수 없는 target_mode: {target_mode}")

    def __getitem__(self, idx):
        features = self._build_features(self.t[idx])
        xi_target = torch.stack([self.xi_real[idx], self.xi_imag[idx]])
        smooth_tgt = self.smooth_target[idx].unsqueeze(0)  # shape [1]
        return features, xi_target, smooth_tgt


# ─────────────────────────────────────────────────────────────────────
# 유틸리티
# ─────────────────────────────────────────────────────────────────────
def load_zeros():
    """result_t100-200.json 에서 영점 목록 로드"""
    with open(os.path.join(OVERNIGHT_DIR, "result_t100-200.json")) as f:
        return json.load(f)["zeros_list"]


def sample_nonzero_t(zeros, n=N_NONZERO, seed=0):
    """영점에서 충분히 떨어진 비영점 t 샘플링"""
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


def get_smooth_target_at_t(dataset, t_values):
    """임의의 t 값들에 대해 smooth 타겟값 반환 (가장 가까운 인덱스 사용)"""
    t_arr = dataset.t.cpu().numpy()
    targets = []
    for t in t_values:
        idx = int(np.argmin(np.abs(t_arr - t)))
        targets.append(float(dataset.smooth_target[idx]))
    return np.array(targets)


# ─────────────────────────────────────────────────────────────────────
# 훈련 루프: 순수 MSE 회귀
# ─────────────────────────────────────────────────────────────────────
def train_smooth(train_loader, val_loader, hidden, seed, device, label):
    """
    MasterResonantNetwork 의 phase_out (마지막 채널 평균) 을 smooth 타겟에
    MSE 회귀. TotalResonanceLoss 대신 순수 MSE 사용.
    """
    torch.manual_seed(seed)
    np.random.seed(seed)

    model = MasterResonantNetwork(
        in_features=IN_FEATURES, hidden_features=hidden, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    mse_loss = nn.MSELoss()
    best_val = float("inf")
    best_state = None
    t0 = time.time()

    for epoch in range(1, EPOCHS + 1):
        model.train()
        for X_batch, _, smooth_tgt in train_loader:
            X_in = X_batch.to(device, dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            tgt = smooth_tgt.to(device, dtype=PrecisionManager.REAL_DTYPE)  # [B, 1]

            optimizer.zero_grad(set_to_none=True)
            outputs = model(X_in)

            # 모델 출력에서 실수 스칼라 추출:
            # phase_out = angle(Z_out), shape [B, out_features]
            # 채널 평균 → [B, 1] 로 매칭
            phase_out = outputs["phase_out"]  # [B, out_features]
            # phase_out 을 [-pi, pi] → [-1, 1] 로 정규화하여 타겟 범위와 매칭
            pred = phase_out.mean(dim=-1, keepdim=True) / math.pi  # [B, 1]

            loss = mse_loss(pred, tgt)
            if torch.isnan(loss) or torch.isinf(loss):
                continue
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()

        # 검증
        model.eval()
        val_acc = 0.0
        n_val = 0
        with torch.enable_grad():
            for X_batch, _, smooth_tgt in val_loader:
                X_in = X_batch.to(device, dtype=PrecisionManager.REAL_DTYPE)
                X_in.requires_grad_(True)
                tgt = smooth_tgt.to(device, dtype=PrecisionManager.REAL_DTYPE)
                out = model(X_in)
                pred = out["phase_out"].mean(dim=-1, keepdim=True) / math.pi
                vl = mse_loss(pred, tgt)
                if not (torch.isnan(vl) or torch.isinf(vl)):
                    val_acc += vl.item()
                    n_val += 1
        val_loss = val_acc / max(1, n_val)
        if val_loss < best_val:
            best_val = val_loss
            best_state = copy.deepcopy(model.state_dict())
        if epoch % 25 == 0 or epoch == EPOCHS:
            print(f"    {label} ep {epoch:03d}: val_mse={val_loss:.6f}", flush=True)

    if best_state:
        model.load_state_dict(best_state)
    return model, best_val, time.time() - t0


# ─────────────────────────────────────────────────────────────────────
# 추론: 영점/비영점 MSE 측정
# ─────────────────────────────────────────────────────────────────────
def evaluate_at_t(model, dataset, t_values, device):
    """
    t_values 에서 모델 예측 vs smooth 타겟 MSE 벡터 반환.
    """
    feats = dataset.get_features_at_t(t_values).to(device)
    targets = get_smooth_target_at_t(dataset, t_values)

    model.eval()
    with torch.enable_grad():
        X = feats.clone().requires_grad_(True)
        out = model(X)
    pred = out["phase_out"].mean(dim=-1).detach().cpu().numpy() / math.pi  # [B]

    # 포인트별 MSE
    mse_per_point = (pred - targets) ** 2
    return pred, targets, mse_per_point


# ─────────────────────────────────────────────────────────────────────
# 기준선 로드: isolation 실험 결과
# ─────────────────────────────────────────────────────────────────────
def load_isolation_baseline():
    """phase_isolation_experiment.json 에서 기준선 MSE 추출"""
    baseline_path = os.path.join(OUTPUT_DIR, "phase_isolation_experiment.json")
    if not os.path.exists(baseline_path):
        print(f"  [경고] 기준선 파일 없음: {baseline_path}", flush=True)
        return None
    with open(baseline_path) as f:
        data = json.load(f)
    # 전체 시드 풀 MSE
    all_mz = []
    all_mn = []
    for seed_key, res in data["results"].items():
        all_mz.extend(res["mse_z_circ"])
        all_mn.extend(res["mse_nz_circ"])
    return {
        "zero_mse_mean": float(np.mean(all_mz)),
        "nonzero_mse_mean": float(np.mean(all_mn)),
        "ratio": float(np.mean(all_mz)) / max(float(np.mean(all_mn)), 1e-12),
        "seeds": data["config"]["seeds"],
        "hidden": data["config"]["hidden"],
        "epochs": data["config"]["epochs"],
    }


# ─────────────────────────────────────────────────────────────────────
# 메인
# ─────────────────────────────────────────────────────────────────────
def main():
    PrecisionManager.setup_precision()
    device = torch.device("cpu")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 70, flush=True)
    print("  Smooth Zero-Crossing 실험", flush=True)
    print(f"  가설: ±pi 불연속이 병목이면, smooth 타겟으로 MSE 갭이 닫힌다", flush=True)
    print(f"  hidden={HIDDEN}, in_features={IN_FEATURES}, seeds={SEEDS}, ep={EPOCHS}", flush=True)
    print(f"  t in [{T_MIN}, {T_MAX}]", flush=True)
    print("=" * 70, flush=True)

    # 영점/비영점 준비
    zeros = load_zeros()
    nonzeros = sample_nonzero_t(zeros)
    print(f"  영점 {len(zeros)}개, 비영점 {len(nonzeros)}개", flush=True)

    # 기준선 로드
    baseline = load_isolation_baseline()
    if baseline:
        print(f"  기준선 (arg xi, isolation): "
              f"영점={baseline['zero_mse_mean']:.4f}, "
              f"비영점={baseline['nonzero_mse_mean']:.4f}, "
              f"비율={baseline['ratio']:.2f}x", flush=True)

    # 캐시 로드
    cache_data = torch.load(
        os.path.join(OVERNIGHT_DIR, "xi_cache_t100-200_n500.pt"),
        weights_only=True,
    )

    # ──── 두 타겟 모드 순회 ────
    TARGET_MODES = {
        "bump": "sign(xi_real) * exp(-xi_real^2/sigma^2)",
        "normalized": "xi_real / max|xi_real|",
    }

    all_results = {}
    for mode, description in TARGET_MODES.items():
        print(f"\n{'=' * 70}", flush=True)
        print(f"  타겟: {mode} — {description}", flush=True)
        print(f"{'=' * 70}", flush=True)

        dataset = XiSmoothTargetDataset(cache_data, in_features=IN_FEATURES, target_mode=mode)

        # 타겟 통계 출력
        tgt = dataset.smooth_target.numpy()
        print(f"  타겟 범위: [{tgt.min():.4f}, {tgt.max():.4f}], "
              f"mean={tgt.mean():.4f}, std={tgt.std():.4f}", flush=True)
        zero_tgts = get_smooth_target_at_t(dataset, zeros)
        print(f"  영점 근처 타겟 절대값 평균: {np.abs(zero_tgts).mean():.6f} "
              f"(0에 가까워야 함)", flush=True)

        # 데이터 분할
        val_size = int(len(dataset) * 0.2)
        train_size = len(dataset) - val_size
        train_ds, val_ds = torch.utils.data.random_split(
            dataset, [train_size, val_size],
            generator=torch.Generator().manual_seed(0),
        )
        train_loader = torch.utils.data.DataLoader(
            train_ds, batch_size=BATCH_SIZE, shuffle=True, drop_last=True,
        )
        val_loader = torch.utils.data.DataLoader(
            val_ds, batch_size=BATCH_SIZE, shuffle=False, drop_last=False,
        )

        mode_results = {}
        for seed in SEEDS:
            label = f"{mode},h={HIDDEN},s={seed}"
            print(f"\n  [훈련] {label}", flush=True)

            model, best_val, train_time = train_smooth(
                train_loader, val_loader, HIDDEN, seed, device, label,
            )
            print(f"    완료: best_val={best_val:.6f}, time={train_time:.1f}s", flush=True)

            # 영점/비영점 평가
            pred_z, tgt_z, mse_z = evaluate_at_t(model, dataset, zeros, device)
            pred_nz, tgt_nz, mse_nz = evaluate_at_t(model, dataset, nonzeros, device)

            mode_results[seed] = {
                "pred_z": pred_z.tolist(),
                "tgt_z": tgt_z.tolist(),
                "mse_z": mse_z.tolist(),
                "pred_nz": pred_nz.tolist(),
                "tgt_nz": tgt_nz.tolist(),
                "mse_nz": mse_nz.tolist(),
                "best_val_loss": float(best_val),
                "train_time": float(train_time),
            }

            print(f"    영점 MSE={mse_z.mean():.6f}, "
                  f"비영점 MSE={mse_nz.mean():.6f}, "
                  f"비율={mse_z.mean() / max(mse_nz.mean(), 1e-12):.4f}", flush=True)

        all_results[mode] = mode_results

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
            "n_zeros": len(zeros),
            "n_nonzeros": len(nonzeros),
            "nonzero_gap": NONZERO_GAP,
            "target_modes": list(TARGET_MODES.keys()),
            "hypothesis": "±pi 불연속이 병목 → smooth 타겟으로 MSE 갭 닫힘",
        },
        "baseline": baseline,
        "zeros": zeros,
        "nonzeros": nonzeros,
        "results": {
            mode: {str(k): v for k, v in mode_res.items()}
            for mode, mode_res in all_results.items()
        },
    }
    with open(JSON_OUT, "w") as f:
        json.dump(out_data, f, indent=2, ensure_ascii=False)

    # ─────────────────────────────────────────────────────────────────
    # 통계 분석 & 텍스트 리포트
    # ─────────────────────────────────────────────────────────────────
    lines = []
    lines.append("=" * 70)
    lines.append("  Smooth Zero-Crossing 실험 결과")
    lines.append(f"  in_features={IN_FEATURES}, hidden={HIDDEN}, seeds={SEEDS}, ep={EPOCHS}")
    lines.append(f"  t in [{T_MIN}, {T_MAX}], 영점 {len(zeros)}개, 비영점 {len(nonzeros)}개")
    lines.append("=" * 70)

    # 기준선
    if baseline:
        lines.append("")
        lines.append("-" * 70)
        lines.append("  기준선: arg xi (isolation 실험, phase_isolation_experiment)")
        lines.append(f"    영점 MSE   = {baseline['zero_mse_mean']:.6f}")
        lines.append(f"    비영점 MSE = {baseline['nonzero_mse_mean']:.6f}")
        lines.append(f"    비율 z/nz  = {baseline['ratio']:.4f}")
        lines.append(f"    (h={baseline['hidden']}, seeds={baseline['seeds']}, ep={baseline['epochs']})")
        lines.append("-" * 70)

    # 각 타겟 모드별 결과
    summary = {}
    for mode in TARGET_MODES:
        mode_res = all_results[mode]
        lines.append("")
        lines.append("=" * 70)
        lines.append(f"  타겟: {mode} — {TARGET_MODES[mode]}")
        lines.append("=" * 70)
        lines.append(f"  {'seed':>6s} {'영점 MSE':>14s} {'비영점 MSE':>14s} "
                     f"{'비율(z/nz)':>12s} {'t-test p':>12s}")

        for seed in SEEDS:
            mz = np.array(mode_res[seed]["mse_z"])
            mn = np.array(mode_res[seed]["mse_nz"])
            ratio = mz.mean() / max(mn.mean(), 1e-12)
            try:
                _, p_val = ttest_ind(mz, mn, equal_var=False)
            except Exception:
                p_val = float("nan")
            lines.append(f"  {seed:>6d} {mz.mean():>14.6f} {mn.mean():>14.6f} "
                         f"{ratio:>12.4f} {p_val:>12.4e}")

        # 시드 풀 통합
        all_mz = np.concatenate([np.array(mode_res[s]["mse_z"]) for s in SEEDS])
        all_mn = np.concatenate([np.array(mode_res[s]["mse_nz"]) for s in SEEDS])
        pool_ratio = all_mz.mean() / max(all_mn.mean(), 1e-12)
        try:
            _, pool_p = ttest_ind(all_mz, all_mn, equal_var=False)
        except Exception:
            pool_p = float("nan")

        lines.append(f"  풀: 영점={all_mz.mean():.6f}, 비영점={all_mn.mean():.6f}, "
                     f"비율={pool_ratio:.4f}, p={pool_p:.4e}")

        summary[mode] = {
            "zero_mse": float(all_mz.mean()),
            "nonzero_mse": float(all_mn.mean()),
            "ratio": float(pool_ratio),
            "p": float(pool_p),
        }

    # ──── 비교 테이블 ────
    lines.append("")
    lines.append("=" * 70)
    lines.append("  비교 요약: 기준선 vs smooth 타겟")
    lines.append("=" * 70)
    lines.append(f"  {'타겟':>14s} {'영점 MSE':>14s} {'비영점 MSE':>14s} "
                 f"{'비율(z/nz)':>12s} {'갭 축소?':>10s}")

    if baseline:
        lines.append(f"  {'arg xi':>14s} {baseline['zero_mse_mean']:>14.6f} "
                     f"{baseline['nonzero_mse_mean']:>14.6f} "
                     f"{baseline['ratio']:>12.4f} {'(기준선)':>10s}")
    for mode in TARGET_MODES:
        s = summary[mode]
        if baseline:
            gap_reduction = (1.0 - s["ratio"] / baseline["ratio"]) * 100
            gap_str = f"{gap_reduction:+.1f}%"
        else:
            gap_str = "N/A"
        lines.append(f"  {mode:>14s} {s['zero_mse']:>14.6f} "
                     f"{s['nonzero_mse']:>14.6f} "
                     f"{s['ratio']:>12.4f} {gap_str:>10s}")

    # ──── 판정 ────
    lines.append("")
    lines.append("=" * 70)
    lines.append("  판정")
    lines.append("=" * 70)

    # 판정 기준: baseline 비율 18x → smooth 비율 < 3x 이면 강한 양성
    best_mode = min(summary, key=lambda m: summary[m]["ratio"])
    best_ratio = summary[best_mode]["ratio"]

    if baseline:
        baseline_ratio = baseline["ratio"]
        if best_ratio < 2.0:
            verdict = (
                f"강한 양성: {best_mode} 타겟의 비율 {best_ratio:.2f}x "
                f"(기준선 {baseline_ratio:.1f}x 대비 극적 감소).\n"
                f"  ±pi 불연속이 핵심 병목이었음을 확인. smooth 타겟으로\n"
                f"  영점/비영점 MSE 갭이 사실상 소멸."
            )
        elif best_ratio < 5.0:
            verdict = (
                f"양성: {best_mode} 타겟의 비율 {best_ratio:.2f}x "
                f"(기준선 {baseline_ratio:.1f}x 대비 유의미한 감소).\n"
                f"  ±pi 불연속이 주요 병목이나, 잔존 갭은 다른 요인 기여."
            )
        elif best_ratio < baseline_ratio * 0.5:
            verdict = (
                f"약한 양성: {best_mode} 타겟의 비율 {best_ratio:.2f}x "
                f"(기준선 {baseline_ratio:.1f}x 대비 절반 이하).\n"
                f"  불연속이 부분 기여. 모델 표현력 한계도 존재."
            )
        else:
            verdict = (
                f"음성: smooth 타겟으로도 갭이 줄지 않음 "
                f"(best ratio {best_ratio:.2f}x vs baseline {baseline_ratio:.1f}x).\n"
                f"  ±pi 불연속이 유일 병목이 아님. 모델 구조/특징 설계 근본 문제."
            )
    else:
        if best_ratio < 2.0:
            verdict = f"양성 시사: {best_mode} 타겟 비율 {best_ratio:.2f}x (기준선 없이 절대 판단)."
        else:
            verdict = f"불확실: {best_mode} 타겟 비율 {best_ratio:.2f}x. 기준선 비교 필요."

    lines.append(f"  {verdict}")

    # 다음 단계 제안
    lines.append("")
    lines.append("-" * 70)
    lines.append("  다음 단계 제안:")
    lines.append("  1. 양성 시: smooth 타겟 기반 새 손실 함수 설계 (연속 부호 학습)")
    lines.append("  2. 양성 시: sign-aware readout — tanh 출력 + BCE 로 부호 분류")
    lines.append("  3. 음성 시: 모델 아키텍처 자체의 표현력 한계 조사")
    lines.append("-" * 70)

    output = "\n".join(lines)
    print("\n" + output, flush=True)
    with open(TXT_OUT, "w") as f:
        f.write(output)

    print(f"\n저장: {JSON_OUT}", flush=True)
    print(f"저장: {TXT_OUT}", flush=True)


if __name__ == "__main__":
    main()
