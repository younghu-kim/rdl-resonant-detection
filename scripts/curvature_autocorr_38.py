#!/usr/bin/env python3
"""
실험 #38: κ 자기상관 함수 C(τ)와 GUE 두점 상관의 관계

설계:
  - σ₀ = 1/2 + 0.03 (δ=0.03), t ∈ [1000, 5000]
  - 균등 격자 Δt = 0.5, N = 8000점
  - C(τ) = Corr(κ(t), κ(t+τ)) 자기상관 계산
  - τ_unf = τ × d̄(t),  d̄ = log(t/2π)/(2π) 평균 영점 밀도
  - GUE 예측: C_pred(τ_unf) ∝ −(sin πτ_unf / πτ_unf)²
  - 비교: 첫 번째 영교차 위치 τ_unf ≈ 1.0 ± 0.2 여부, 감쇠 형태 R²

주의:
  - ★ κ: 해석적 G(s)+ζ'/ζ (h=1e-20 수치미분 금지)
  - ★ δ=0.03 오프셋 필수
  - ★ τ_max/total = 100/4000 = 0.025 < 0.25 ✓
  - ★ 5σ 아웃라이어 클리핑 필수
  - ★ except: pass 금지
  - ★ np.trapezoid (trapz 금지)

성공 기준:
  - 양성: 첫 번째 영교차 τ_unf ≈ 1.0 ± 0.2
  - 강한 양성: + 형태 피팅 R² > 0.80
  - 음성: 영교차 없음 또는 τ_unf > 2.0 또는 비단조 감쇠
"""

import os
import sys
import time
import numpy as np
import mpmath
from scipy import stats as scipy_stats
from scipy.optimize import curve_fit

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_PATH = os.path.join(BASE_DIR, 'results', 'curvature_autocorr_38.txt')
CACHE_PATH = os.path.join(BASE_DIR, 'outputs', 'cache', 'autocorr_kappa_38.npy')

DELTA = 0.03
T_MIN = 1000.0
T_MAX = 5000.0
DT = 0.5
K_MAX = 200   # 최대 lag (τ_max = K_MAX * DT = 100)
SIGMA_CLIP = 5.0  # 아웃라이어 클리핑

t0_global = time.time()
lines = []

def log(msg=""):
    s = str(msg)
    print(s, flush=True)
    lines.append(s)

def elapsed():
    return time.time() - t0_global


# ═══════════════════════════════════════════════════════════════════════════
# 해석적 κ 계산 (bundle_utils h=1e-20 금지 — 직접 구현)
# ═══════════════════════════════════════════════════════════════════════════

def get_dps(t):
    if t > 10000: return 100
    if t > 5000:  return 80
    if t > 1000:  return 60
    return max(50, int(30 + t / 200))


def G_component(s):
    """G(s) = 1/s + 1/(s−1) − log(π)/2 + ψ(s/2)/2 (해석적)"""
    return (mpmath.mpf(1) / s
            + mpmath.mpf(1) / (s - 1)
            - mpmath.log(mpmath.pi) / 2
            + mpmath.digamma(s / 2) / 2)


def zeta_log_deriv(s):
    """ζ'/ζ(s) — h=1e-6 중앙 차분 (mpmath.diff 금지, h=1e-20 금지)"""
    z = mpmath.zeta(s)
    if abs(z) < mpmath.mpf(10) ** (-mpmath.mp.dps + 10):
        return mpmath.mpc(1e8, 0)
    h = mpmath.mpf('1e-6')
    zd = (mpmath.zeta(s + h) - mpmath.zeta(s - h)) / (2 * h)
    return zd / z


def compute_kappa_analytic(t, delta=DELTA):
    """해석적 공식으로 κ = |ξ'/ξ(1/2+δ+it)|² 계산"""
    dps = get_dps(t)
    with mpmath.workdps(dps):
        s = mpmath.mpc(0.5 + delta, t)
        xi_ld = G_component(s) + zeta_log_deriv(s)
        kappa = float(abs(xi_ld) ** 2)
    return kappa


# ═══════════════════════════════════════════════════════════════════════════
# GUE 두점 상관 예측
# ═══════════════════════════════════════════════════════════════════════════

def gue_connected_pair(r):
    """
    GUE 연결 쌍 상관: b₂^c(r) = −(sin πr / πr)²
    r=0일 때: 극한값 −1 (sinc(0) = 1)
    κ 자기상관 예측: C_pred(τ_unf) ∝ b₂^c(τ_unf)
    """
    r = np.asarray(r, dtype=float)
    result = np.zeros_like(r)
    nonzero = r != 0
    x = np.pi * r[nonzero]
    result[nonzero] = -(np.sin(x) / x) ** 2
    result[~nonzero] = -1.0  # 극한값
    return result


def mean_zero_density(t):
    """평균 영점 밀도 d̄(t) = log(t/2π) / (2π) [Riemann-Siegel]"""
    return np.log(t / (2.0 * np.pi)) / (2.0 * np.pi)


# ═══════════════════════════════════════════════════════════════════════════
# 메인 실행
# ═══════════════════════════════════════════════════════════════════════════

def main():
    log("=" * 70)
    log("실험 #38: κ 자기상관 함수 C(τ) vs GUE 두점 상관")
    log("=" * 70)
    log(f"σ₀ = 1/2 + {DELTA} = {0.5+DELTA:.2f}")
    log(f"t ∈ [{T_MIN:.0f}, {T_MAX:.0f}], Δt = {DT}, K_max = {K_MAX}")
    log(f"τ_max = {K_MAX * DT:.1f}, total_range = {T_MAX - T_MIN:.0f}")
    log(f"τ_max/total = {K_MAX * DT / (T_MAX - T_MIN):.4f} (< 0.25 ✓)")
    log()

    # ─────────────────────────────────────────────────────────────────────
    # 파트 A: κ 계산 (캐시 우선)
    # ─────────────────────────────────────────────────────────────────────
    log("─" * 50)
    log("파트 A: κ 계산 (N=8000, t∈[1000,5000])")
    log("─" * 50)

    t_grid = np.arange(T_MIN, T_MAX + DT / 2, DT)
    N = len(t_grid)
    log(f"  격자점 수: N = {N}")
    log(f"  t 범위: [{t_grid[0]:.1f}, {t_grid[-1]:.1f}]")

    if os.path.exists(CACHE_PATH):
        log(f"  캐시 발견: {CACHE_PATH}")
        kappa_arr = np.load(CACHE_PATH)
        if len(kappa_arr) == N:
            log(f"  ✓ 캐시 로드 완료 (N={len(kappa_arr)})")
        else:
            log(f"  ⚠️ 캐시 크기 불일치 ({len(kappa_arr)} ≠ {N}) — 재계산")
            kappa_arr = None
    else:
        kappa_arr = None

    if kappa_arr is None:
        log(f"  κ 계산 시작 (예상 ~{int(N*0.1/60+0.5)}분)...")
        kappa_arr = np.zeros(N)
        fail_count = 0
        t_start_calc = time.time()
        for i, t in enumerate(t_grid):
            try:
                kappa_arr[i] = compute_kappa_analytic(t, DELTA)
                if not np.isfinite(kappa_arr[i]) or kappa_arr[i] <= 0:
                    log(f"  WARNING: t={t:.1f} κ 비정상 ({kappa_arr[i]:.4g}), 재설정=NaN")
                    kappa_arr[i] = np.nan
                    fail_count += 1
            except Exception as e:
                log(f"  WARNING: t={t:.1f} 계산 실패: {e}")
                kappa_arr[i] = np.nan
                fail_count += 1

            if (i + 1) % 500 == 0:
                elapsed_sec = time.time() - t_start_calc
                rate = (i + 1) / elapsed_sec
                eta = (N - i - 1) / rate
                log(f"  {i+1}/{N} ({100*(i+1)/N:.1f}%) | {elapsed_sec:.0f}s 경과 | ETA {eta:.0f}s")

            # 중간 저장 (500점마다)
            if (i + 1) % 500 == 0:
                np.save(CACHE_PATH + '.tmp', kappa_arr)

        if fail_count > N // 2:
            log(f"  ⚠️ 실패 {fail_count}/{N} — 절반 이상 실패. 중단.")
            sys.exit(1)

        log(f"  κ 계산 완료: {time.time()-t_start_calc:.0f}s, 실패={fail_count}")
        np.save(CACHE_PATH, kappa_arr)
        log(f"  캐시 저장: {CACHE_PATH}")

    # NaN 처리
    valid_mask = np.isfinite(kappa_arr) & (kappa_arr > 0)
    log(f"  유효점: {valid_mask.sum()}/{N} ({100*valid_mask.sum()/N:.1f}%)")

    if valid_mask.sum() < N // 2:
        log("  ⚠️ 유효점 50% 미만 — 중단")
        sys.exit(1)

    # ─────────────────────────────────────────────────────────────────────
    # 파트 B: 5σ 아웃라이어 클리핑
    # ─────────────────────────────────────────────────────────────────────
    log()
    log("─" * 50)
    log("파트 B: 5σ 아웃라이어 클리핑 (κ 발산점 제거)")
    log("─" * 50)

    kappa_valid = kappa_arr.copy()
    kappa_valid[~valid_mask] = np.nan

    # 유효값에서 5σ 클리핑
    finite_vals = kappa_valid[valid_mask]
    kappa_mean = np.nanmean(finite_vals)
    kappa_std = np.nanstd(finite_vals)
    threshold_hi = kappa_mean + SIGMA_CLIP * kappa_std
    threshold_lo = kappa_mean - SIGMA_CLIP * kappa_std

    log(f"  클리핑 전: mean={kappa_mean:.2f}, std={kappa_std:.2f}")
    log(f"  임계값: [{threshold_lo:.2f}, {threshold_hi:.2f}]")

    outlier_hi = valid_mask & (kappa_arr > threshold_hi)
    outlier_lo = valid_mask & (kappa_arr < threshold_lo) & (kappa_arr > 0)
    outlier_mask = outlier_hi | outlier_lo

    log(f"  아웃라이어: 상위 {outlier_hi.sum()}개, 하위 {outlier_lo.sum()}개")
    log(f"  총 클리핑: {outlier_mask.sum()}개 ({100*outlier_mask.sum()/N:.2f}%)")

    # 클리핑: 아웃라이어를 경계값으로 교체
    kappa_clipped = kappa_arr.copy()
    kappa_clipped[outlier_hi] = threshold_hi
    kappa_clipped[outlier_lo] = threshold_lo
    kappa_clipped[~valid_mask] = kappa_mean  # NaN → 평균 대체 (자기상관 연속성 위해)

    kappa_mean_c = np.mean(kappa_clipped)
    kappa_std_c = np.std(kappa_clipped)
    log(f"  클리핑 후: mean={kappa_mean_c:.2f}, std={kappa_std_c:.2f}")
    log(f"  κ 기준값 1/δ² = {1/DELTA**2:.1f} (참고)")

    # ─────────────────────────────────────────────────────────────────────
    # 파트 C: 자기상관 함수 계산 (FFT 방법)
    # ─────────────────────────────────────────────────────────────────────
    log()
    log("─" * 50)
    log("파트 C: 자기상관 함수 C(k) 계산 (FFT)")
    log("─" * 50)

    # 평균 제거
    kappa_detrended = kappa_clipped - kappa_mean_c

    # FFT 기반 자기상관 (정확하고 빠름)
    n = len(kappa_detrended)
    # zero-padding으로 circular 효과 방지
    fft_k = np.fft.rfft(kappa_detrended, n=2 * n)
    acf_raw = np.fft.irfft(fft_k * np.conj(fft_k))[:n]
    # 정규화: k번째 lag는 (N-k)개 쌍으로 추정
    norm = np.arange(n, 0, -1, dtype=float)
    acf_biased = acf_raw / n           # biased estimator
    acf_unbiased = acf_raw / norm      # unbiased estimator (큰 lag에서 노이즈 ↑)

    # k=0..K_MAX만 사용
    lag_indices = np.arange(0, K_MAX + 1)
    acf_vals = acf_biased[lag_indices]

    # 정규화: C(0) = 1
    C = acf_vals / acf_vals[0]

    log(f"  총 lag 수: {len(lag_indices)} (k=0..{K_MAX})")
    log(f"  C(0) = {C[0]:.6f} (= 1.0 by design)")
    log(f"  C(1) = {C[1]:.6f}")
    log(f"  C(10) = {C[10]:.6f}")
    log(f"  C(50) = {C[50]:.6f}")
    log(f"  C(200) = {C[200]:.6f}")

    # ─────────────────────────────────────────────────────────────────────
    # 파트 D: τ → τ_unf 변환 (unfolded 좌표)
    # ─────────────────────────────────────────────────────────────────────
    log()
    log("─" * 50)
    log("파트 D: τ → τ_unf 변환")
    log("─" * 50)

    # 평균 영점 밀도 d̄: t∈[1000,5000] 격자에서 평균
    d_bar_grid = mean_zero_density(t_grid)
    d_bar_mean = np.mean(d_bar_grid)
    d_bar_mid = mean_zero_density(0.5 * (T_MIN + T_MAX))  # 중점에서 계산

    log(f"  d̄_mean = {d_bar_mean:.6f} (격자 평균)")
    log(f"  d̄_mid  = {d_bar_mid:.6f} (t=3000 중점)")
    log(f"  평균 영점 간격 1/d̄ = {1/d_bar_mean:.4f}")

    # τ_unf = k * Δt * d̄_mean
    tau_unf = lag_indices * DT * d_bar_mean

    log(f"  τ_unf(k=1)  = {tau_unf[1]:.4f}")
    log(f"  τ_unf(k=10) = {tau_unf[10]:.4f}")
    log(f"  τ_unf(k={K_MAX}) = {tau_unf[K_MAX]:.4f}")
    log(f"  (GUE 첫 영교차 예측: τ_unf ≈ 1.0)")

    # ─────────────────────────────────────────────────────────────────────
    # 파트 E: GUE 두점 상관 비교
    # ─────────────────────────────────────────────────────────────────────
    log()
    log("─" * 50)
    log("파트 E: GUE 두점 상관 예측과 비교")
    log("─" * 50)

    # GUE 연결 쌍 상관: b₂^c(r) = −(sin πr / πr)²
    # τ_unf > 0에서의 예측 (k>0)
    tau_fit = tau_unf[1:]   # k=1..K_MAX (k=0 제외: C(0)=1은 별도)
    C_emp_fit = C[1:]

    gue_pred_raw = gue_connected_pair(tau_fit)  # shape 동일, 값 ∈ [-1, 0]

    # 1) 첫 번째 영교차 탐색
    log("  [E1] 첫 번째 영교차 탐색")
    zero_cross_k = None
    zero_cross_tau_unf = None
    for k in range(0, K_MAX):
        if C[k] >= 0 and C[k + 1] < 0:
            # 선형 보간
            frac = C[k] / (C[k] - C[k + 1])
            zero_cross_k = k + frac
            zero_cross_tau_unf = zero_cross_k * DT * d_bar_mean
            break
        elif C[k] <= 0 and C[k + 1] > 0:
            # 반대 방향 교차
            frac = -C[k] / (C[k + 1] - C[k])
            zero_cross_k = k + frac
            zero_cross_tau_unf = zero_cross_k * DT * d_bar_mean
            break

    if zero_cross_tau_unf is not None:
        log(f"  첫 번째 영교차: k ≈ {zero_cross_k:.1f}, τ_unf ≈ {zero_cross_tau_unf:.4f}")
        gue_zero = 1.0  # GUE 예측 τ_unf = 1.0
        dev = abs(zero_cross_tau_unf - gue_zero)
        log(f"  GUE 예측값:      τ_unf = {gue_zero:.1f}")
        log(f"  편차:            |{zero_cross_tau_unf:.4f} - {gue_zero:.1f}| = {dev:.4f}")
        if dev <= 0.2:
            log(f"  ✅ 영교차 일치! (편차 ≤ 0.2)")
        elif dev <= 0.5:
            log(f"  ⚠️ 영교차 근접 (편차 0.2 < {dev:.4f} ≤ 0.5)")
        else:
            log(f"  ❌ 영교차 불일치 (편차 {dev:.4f} > 0.5)")
    else:
        log(f"  ❌ k=1..{K_MAX} 범위에서 영교차 없음")
        zero_cross_tau_unf = None

    # 2) 형태 비교: 선형 회귀 (C_emp vs C_pred_shape)
    log()
    log("  [E2] 형태 비교: C_emp(τ_unf) ~ A × b₂^c(τ_unf)")
    # 단순 선형 회귀: C_emp ≈ A × gue_pred_raw
    # (offset=0으로 설정, 형태만 비교)
    # 단, C(0)=1 ↔ b²^c(0)=-1이므로 k=1..K_MAX만 사용

    # R² 계산: C_emp vs A × gue_pred_raw (단일 파라미터 최소제곱)
    mask_fit = tau_fit <= 20  # τ_unf ≤ 20 범위에서 피팅 (통계 신뢰 구간 내)
    if mask_fit.sum() < 5:
        mask_fit = np.ones(len(tau_fit), dtype=bool)

    C_fit = C_emp_fit[mask_fit]
    G_fit = gue_pred_raw[mask_fit]

    # 1파라미터 최소제곱: A = (C_fit · G_fit) / (G_fit · G_fit)
    A_best = np.dot(C_fit, G_fit) / np.dot(G_fit, G_fit)
    C_pred_fit = A_best * G_fit

    SS_res = np.sum((C_fit - C_pred_fit) ** 2)
    SS_tot = np.sum((C_fit - np.mean(C_fit)) ** 2)
    R2_shape = 1.0 - SS_res / SS_tot if SS_tot > 0 else 0.0

    log(f"  피팅 범위: τ_unf ≤ 20 (N_fit = {mask_fit.sum()}점)")
    log(f"  최적 진폭: A = {A_best:.4f}")
    log(f"  R² (형태) = {R2_shape:.4f}")

    if R2_shape > 0.80:
        log(f"  ✅ R² > 0.80 → 강한 형태 일치")
    elif R2_shape > 0.50:
        log(f"  ⚠️ 0.50 < R² = {R2_shape:.4f} ≤ 0.80 → 중간 형태 일치")
    else:
        log(f"  ❌ R² = {R2_shape:.4f} ≤ 0.50 → 형태 불일치")

    # 3) 감쇠 형태: log|C(τ)| vs τ_unf 기울기 (단조 감쇠 확인)
    log()
    log("  [E3] 감쇠 형태 분석")
    k_decay = lag_indices[1:50]
    C_decay = C[1:50]
    C_abs = np.abs(C_decay)

    # log|C| vs τ_unf 선형 회귀 (첫 영교차 이전 구간)
    k_before_zero = lag_indices[1:] if zero_cross_k is None else lag_indices[1:int(zero_cross_k)]
    if len(k_before_zero) > 3:
        C_before = C[1:len(k_before_zero)+1]
        tau_before = k_before_zero * DT * d_bar_mean
        # 단조 감쇠 여부
        is_monotone = np.all(np.diff(C_before) <= 0)
        log(f"  첫 영교차 이전 단조 감소: {'✅ 예' if is_monotone else '❌ 아니오'}")

        # 감쇠율 추정: log|C| ~ -B * τ_unf
        log_C = np.log(np.abs(C_before) + 1e-10)
        if len(tau_before) > 3:
            slope, intercept, r_val, p_val, _ = scipy_stats.linregress(tau_before, log_C)
            log(f"  감쇠율 추정: B = {-slope:.4f} (τ_decay = {-1/slope:.4f})")
            log(f"  log|C|-τ R² = {r_val**2:.4f}")
    else:
        log("  (영교차 이전 구간 너무 짧음)")
        is_monotone = None

    # 4) 상세 자기상관 테이블
    log()
    log("  [E4] 자기상관 테이블 (k=0, 1, 2, ... 20, 50, 100, 200)")
    log(f"  {'k':>5}  {'τ':>7}  {'τ_unf':>7}  {'C(k)':>10}  {'b₂^c':>10}  판정")
    sample_ks = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 50, 100, 150, 200]
    for k in sample_ks:
        if k > K_MAX:
            continue
        tau_k = k * DT
        tau_unf_k = tau_k * d_bar_mean
        c_k = C[k]
        if k == 0:
            b2c_k = -1.0
        else:
            b2c_k = gue_connected_pair(np.array([tau_unf_k]))[0]
        flag = ""
        if k > 0 and abs(tau_unf_k - 1.0) < 0.1:
            flag = " ← GUE 예측 영교차"
        log(f"  {k:>5}  {tau_k:>7.2f}  {tau_unf_k:>7.4f}  {c_k:>+10.5f}  {b2c_k:>+10.5f}{flag}")

    # ─────────────────────────────────────────────────────────────────────
    # 파트 F: 2차 분석 — 단기 구조 + 분산 스펙트럼
    # ─────────────────────────────────────────────────────────────────────
    log()
    log("─" * 50)
    log("파트 F: 추가 분석")
    log("─" * 50)

    # F1: C(τ)의 파워 스펙트럼 (주기성 확인)
    freqs = np.fft.rfftfreq(len(C), d=DT * d_bar_mean)
    power = np.abs(np.fft.rfft(C)) ** 2
    # 주요 주파수 (f=1 ↔ τ_unf=1 주기성 여부)
    freq_resolution = freqs[1] if len(freqs) > 1 else 0
    log(f"  주파수 해상도: Δf = {freq_resolution:.6f}")

    # f ≈ 1 근처 (τ_unf ≈ 1 주기) 파워
    idx_f1 = np.argmin(np.abs(freqs - 1.0))
    if 0 < idx_f1 < len(power):
        power_f1 = power[idx_f1]
        power_total = np.sum(power[1:])
        log(f"  f≈1.0 (τ_unf≈1 주기) 파워: {power_f1:.2f} ({100*power_f1/power_total:.1f}%)")
    log(f"  저주파(f<0.1) 파워 비율: {100*np.sum(power[freqs<0.1])/np.sum(power[1:]):.1f}%")

    # F2: C(τ) 특성 지표
    log()
    # 양의 영역 면적 (C>0인 구간)
    pos_area = np.trapezoid(np.maximum(C[1:], 0), tau_unf[1:])
    neg_area = np.trapezoid(np.minimum(C[1:], 0), tau_unf[1:])
    log(f"  C>0 면적 (τ_unf): {pos_area:.6f}")
    log(f"  C<0 면적 (τ_unf): {neg_area:.6f}")
    log(f"  총 적분: {pos_area + neg_area:.6f}")
    log(f"  |음수/양수| 비율: {abs(neg_area)/max(pos_area, 1e-10):.4f}")

    # F3: 정규화된 C_pred와의 상관 (다른 구간별 R²)
    log()
    log("  R² 피팅 결과 (τ_unf 범위별):")
    ranges = [(0, 5), (0, 10), (0, 20), (0, 50), (0, K_MAX * DT * d_bar_mean)]
    for r_min, r_max in ranges:
        mask_r = (tau_fit >= r_min) & (tau_fit <= r_max)
        if mask_r.sum() < 5:
            continue
        C_r = C_emp_fit[mask_r]
        G_r = gue_pred_raw[mask_r]
        A_r = np.dot(C_r, G_r) / max(np.dot(G_r, G_r), 1e-20)
        SS_res_r = np.sum((C_r - A_r * G_r) ** 2)
        SS_tot_r = np.sum((C_r - np.mean(C_r)) ** 2)
        R2_r = 1.0 - SS_res_r / SS_tot_r if SS_tot_r > 0 else 0.0
        log(f"  τ_unf ∈ [{r_min:.0f}, {r_max:.1f}]: R² = {R2_r:.4f}, A = {A_r:.4f}")

    # ─────────────────────────────────────────────────────────────────────
    # 파트 G: 종합 판정
    # ─────────────────────────────────────────────────────────────────────
    log()
    log("═" * 70)
    log("파트 G: 종합 판정")
    log("═" * 70)

    log()
    log("성공 기준 체크:")
    criteria = {}

    # 기준 1: 첫 번째 영교차 τ_unf ≈ 1.0 ± 0.2
    if zero_cross_tau_unf is not None:
        dev1 = abs(zero_cross_tau_unf - 1.0)
        criteria['zero_cross'] = dev1 <= 0.2
        criteria['zero_cross_near'] = dev1 <= 0.5
        log(f"  C1. 첫 영교차 τ_unf = {zero_cross_tau_unf:.4f}  (GUE=1.0, 편차={dev1:.4f})")
        log(f"      → {'✅ 통과' if criteria['zero_cross'] else ('⚠️ 근접' if criteria['zero_cross_near'] else '❌ 실패')}")
    else:
        criteria['zero_cross'] = False
        criteria['zero_cross_near'] = False
        log(f"  C1. 첫 영교차: 없음 → ❌ 실패")

    # 기준 2: 형태 피팅 R² > 0.80
    criteria['shape_strong'] = R2_shape > 0.80
    criteria['shape_weak'] = R2_shape > 0.50
    log(f"  C2. 형태 R² = {R2_shape:.4f}")
    log(f"      → {'✅ 강한 일치(>0.80)' if criteria['shape_strong'] else ('⚠️ 중간 일치(0.50-0.80)' if criteria['shape_weak'] else '❌ 일치 없음')}")

    # 기준 3: 단조 감쇠 (영교차 이전)
    criteria['monotone'] = is_monotone if is_monotone is not None else False
    log(f"  C3. 단조 감쇠: {'✅' if criteria['monotone'] else '❌'}")

    log()
    # 최종 판정
    if criteria['zero_cross'] and criteria['shape_strong']:
        verdict = "★★ 강한 양성"
        detail = "영교차 일치(±0.2) + R²>0.80"
    elif criteria['zero_cross'] and criteria['shape_weak']:
        verdict = "★ 양성"
        detail = "영교차 일치(±0.2) + 형태 중간 일치"
    elif criteria['zero_cross']:
        verdict = "★ 양성 (형태 불일치)"
        detail = "영교차 일치(±0.2)이나 R² 낮음"
    elif criteria['zero_cross_near']:
        verdict = "⚠️ 중립"
        detail = f"영교차 근접(τ_unf={zero_cross_tau_unf:.3f}±0.5) 이나 편차 큼"
    else:
        verdict = "❌ 음성"
        if zero_cross_tau_unf is None:
            detail = "영교차 없음"
        else:
            detail = f"영교차 τ_unf={zero_cross_tau_unf:.3f} (GUE=1.0, 편차>{0.5})"

    log(f"  ╔══════════════════════════════════╗")
    log(f"  ║  최종 판정: {verdict:<20}║")
    log(f"  ║  근거: {detail:<26}║")
    log(f"  ╚══════════════════════════════════╝")

    log()
    zc_str = f"{zero_cross_tau_unf:.4f}" if zero_cross_tau_unf is not None else "N/A"
    log(f"  첫 영교차 τ_unf = {zc_str}")
    log(f"  GUE 예측 τ_unf  = 1.0000")
    log(f"  형태 R² = {R2_shape:.4f}")
    log(f"  진폭 A = {A_best:.4f}")
    log(f"  총 소요: {elapsed():.1f}초 ({elapsed()/60:.1f}분)")

    # 결과 저장
    os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
    with open(RESULTS_PATH, 'w') as f:
        f.write('\n'.join(lines))
    log()
    log(f"결과 저장: {RESULTS_PATH}")

    # 자기상관 배열도 저장 (추후 분석용)
    acf_save = np.column_stack([lag_indices, lag_indices * DT, tau_unf, C])
    np.save(CACHE_PATH.replace('.npy', '_acf.npy'), acf_save)
    log(f"자기상관 배열 저장: {CACHE_PATH.replace('.npy', '_acf.npy')}")


if __name__ == '__main__':
    main()
