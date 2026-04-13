#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lehmer 쌍 곡률 분석 — 근접 영점에서의 κ 거동
==============================================
사이클 #22 수학자 지시:
  - 다발 프레임워크의 "정량적 예측" 검증
  - 근접 영점 → 위상 체류시간 증가 → κ 프로파일 변화
  - GUE 간격 분포와 κ/|F₂| 수렴 속도의 상관

실험:
  1. t∈[0,200] 영점 ~80개 수집, 정규화 간격 계산
  2. 최근접 쌍 5개 vs 최격리 5개: κ 프로파일 비교
  3. |F₂| 수렴 속도: Hadamard 부분합 수렴 곡선
  4. 상관 분석: g_norm vs κ 특성값 (Spearman)
"""

import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import mpmath
from scipy import stats

from bundle_utils import (
    xi_func, connection_zeta, curvature_zeta,
    find_zeros_zeta, _near_zero,
)

mpmath.mp.dps = 100  # 고정밀도 (t~200 대비)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(SCRIPT_DIR, '..', 'results', 'lehmer_curvature_analysis.txt')


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 유틸리티
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

class Tee:
    """stdout + 파일 동시 출력"""
    def __init__(self, *files):
        self.files = files

    def write(self, data):
        for f in self.files:
            f.write(data)
            f.flush()

    def flush(self):
        for f in self.files:
            f.flush()


def log(msg=""):
    print(msg, flush=True)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1단계: 영점 수집 + 정규화 간격
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def collect_zeros_and_spacings(t_max=200.0):
    """영점 수집 + 정규화 간격 g_norm 계산"""
    log("=" * 70)
    log("1단계: 영점 수집 및 정규화 간격")
    log("=" * 70)

    zeros = find_zeros_zeta(0, t_max)
    N = len(zeros)
    log(f"  t in [0, {t_max}] 영점 {N}개 수집")

    if N < 10:
        log("  ⚠️ 영점 부족 — 탐색 로직 점검 필요")
        return None, None, None, None, None, None

    # 연속 영점 간격 Δt_n = t_{n+1} - t_n
    spacings = np.diff(zeros)

    # 정규화: g_norm = Δt / <Δt>_local, 여기서 <Δt> = 2π / log(t/(2π))
    t_mid = (zeros[:-1] + zeros[1:]) / 2.0
    # t < 2π ≈ 6.28 이면 log(t/2π) < 0 → 무효
    valid_mask = t_mid > 2 * np.pi
    mean_spacing = np.full_like(spacings, np.nan)
    mean_spacing[valid_mask] = 2.0 * np.pi / np.log(t_mid[valid_mask] / (2.0 * np.pi))

    g_norm = np.full_like(spacings, np.nan)
    g_norm[valid_mask] = spacings[valid_mask] / mean_spacing[valid_mask]

    valid_idx = np.where(valid_mask)[0]
    g_valid = g_norm[valid_idx]
    log(f"  유효 간격 수: {len(valid_idx)}")
    log(f"  g_norm 통계: min={np.nanmin(g_valid):.4f}, max={np.nanmax(g_valid):.4f}, "
        f"mean={np.nanmean(g_valid):.4f}, std={np.nanstd(g_valid):.4f}")

    # 가장 근접한 쌍 5개 / 최격리 5개
    sorted_by_g = valid_idx[np.argsort(g_valid)]
    closest_5 = sorted_by_g[:5]
    isolated_5 = sorted_by_g[-5:]

    log(f"\n  ★ 최근접 쌍 5개 (Lehmer-like):")
    for idx in closest_5:
        log(f"    n={idx:3d}: t={zeros[idx]:.6f}, t+1={zeros[idx+1]:.6f}, "
            f"Δt={spacings[idx]:.6f}, g_norm={g_norm[idx]:.4f}")

    log(f"\n  ★ 최격리 영점 5개:")
    for idx in isolated_5:
        log(f"    n={idx:3d}: t={zeros[idx]:.6f}, t+1={zeros[idx+1]:.6f}, "
            f"Δt={spacings[idx]:.6f}, g_norm={g_norm[idx]:.4f}")

    return zeros, spacings, g_norm, closest_5, isolated_5, valid_idx


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# κ 프로파일 계산 (영점 중심 ±half_width, 영점 위 δ 오프셋)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_kappa_profile(t_center, half_width=2.0, n_points=200, delta=0.05,
                          avoid_zeros=None):
    """
    t_center 중심 ±half_width에서 κ(0.5 + it) 프로파일 계산.
    avoid_zeros: 영점 좌표 리스트 — 이 점 근처(δ 이내)에서 오프셋 적용.
    """
    ts = np.linspace(t_center - half_width, t_center + half_width, n_points)
    kappas = np.zeros(n_points)

    for i, t in enumerate(ts):
        t_eval = t
        # 알려진 영점과의 거리 < δ 이면 오프셋
        if avoid_zeros is not None:
            for tz in avoid_zeros:
                if abs(t - tz) < delta:
                    t_eval = tz + delta
                    break

        try:
            s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_eval))
            xi_val = xi_func(s)
            if _near_zero(xi_val):
                kappas[i] = 1e8  # 영점 위 — 매우 큰 값으로 표시
                continue
            L = connection_zeta(s)
            kappa = float(abs(L)**2)
            if not np.isfinite(kappa):
                kappas[i] = 0.0
            else:
                kappas[i] = kappa
        except Exception as e:
            log(f"    WARNING: κ({t_eval:.4f}) 실패: {e}")
            kappas[i] = 0.0

    return ts, kappas


def extract_profile_chars(ts, kappas):
    """κ 프로파일에서 특성 추출: 피크높이, FWHM, 피크수"""
    # 유효 값만 사용 (1e8 cap 값 제외)
    valid = kappas < 1e7
    if np.sum(valid) < 3:
        return {'peak_height': 0, 'n_peaks': 0, 'fwhm': 0, 'mean_kappa': 0}

    kappas_clean = kappas.copy()
    kappas_clean[~valid] = 0.0
    median_k = np.median(kappas_clean[valid])
    if median_k <= 0:
        median_k = 1.0

    # 피크 찾기 (배경 대비 3배 이상)
    threshold = max(median_k * 3.0, 10.0)
    peaks = []
    for i in range(1, len(kappas_clean) - 1):
        if kappas_clean[i] > kappas_clean[i-1] and kappas_clean[i] > kappas_clean[i+1]:
            if kappas_clean[i] > threshold:
                peaks.append(i)

    peak_height = float(np.max(kappas_clean)) if len(kappas_clean) > 0 else 0.0
    n_peaks = len(peaks)

    # FWHM 계산 (가장 높은 피크)
    fwhm = 0.0
    if n_peaks > 0:
        best_peak = peaks[np.argmax([kappas_clean[p] for p in peaks])]
        half_max = kappas_clean[best_peak] / 2.0
        left = best_peak
        while left > 0 and kappas_clean[left] > half_max:
            left -= 1
        right = best_peak
        while right < len(kappas_clean) - 1 and kappas_clean[right] > half_max:
            right += 1
        fwhm = float(ts[min(right, len(ts)-1)] - ts[max(left, 0)])

    return {
        'peak_height': peak_height,
        'n_peaks': n_peaks,
        'fwhm': fwhm,
        'mean_kappa': float(np.mean(kappas_clean[valid])),
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2단계: 근접쌍 vs 격리영점 κ 프로파일 비교
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def experiment_2(zeros, spacings, closest_5, isolated_5):
    log("\n" + "=" * 70)
    log("2단계: κ 프로파일 비교 (근접 5쌍 vs 격리 5영점)")
    log("=" * 70)

    close_results = []
    iso_results = []

    log("\n--- 근접 쌍 κ 프로파일 ---")
    for idx in closest_5:
        t1, t2 = zeros[idx], zeros[idx + 1]
        gap = t2 - t1
        t_mid = (t1 + t2) / 2.0
        # 창 크기: 간격의 10배 또는 최소 2.0
        half_w = max(2.0, gap * 10.0)
        log(f"\n  쌍 n={idx}: t₁={t1:.6f}, t₂={t2:.6f}, Δt={gap:.6f}")

        ts, kappas = compute_kappa_profile(
            t_mid, half_width=half_w, n_points=200, delta=0.03,
            avoid_zeros=[t1, t2]
        )
        chars = extract_profile_chars(ts, kappas)
        close_results.append({
            'idx': idx, 't1': t1, 't2': t2, 'gap': gap,
            'chars': chars,
        })
        log(f"    피크수={chars['n_peaks']}, 피크높이={chars['peak_height']:.2f}, "
            f"FWHM={chars['fwhm']:.4f}, 평균κ={chars['mean_kappa']:.2f}")

    log("\n--- 격리 영점 κ 프로파일 ---")
    for idx in isolated_5:
        t1 = zeros[idx]
        log(f"\n  영점 n={idx}: t={t1:.6f}, Δt(좌)={spacings[idx]:.6f}")

        ts, kappas = compute_kappa_profile(
            t1, half_width=2.0, n_points=200, delta=0.03,
            avoid_zeros=[t1]
        )
        chars = extract_profile_chars(ts, kappas)
        iso_results.append({
            'idx': idx, 't1': t1, 'gap': spacings[idx],
            'chars': chars,
        })
        log(f"    피크수={chars['n_peaks']}, 피크높이={chars['peak_height']:.2f}, "
            f"FWHM={chars['fwhm']:.4f}, 평균κ={chars['mean_kappa']:.2f}")

    # 비교 요약
    log("\n--- 프로파일 비교 요약 ---")
    cH = [r['chars']['peak_height'] for r in close_results]
    iH = [r['chars']['peak_height'] for r in iso_results]
    cF = [r['chars']['fwhm'] for r in close_results]
    iF = [r['chars']['fwhm'] for r in iso_results]
    cN = [r['chars']['n_peaks'] for r in close_results]
    iN = [r['chars']['n_peaks'] for r in iso_results]

    log(f"  근접쌍:  피크높이 {np.mean(cH):.1f} ± {np.std(cH):.1f}, "
        f"FWHM {np.mean(cF):.4f} ± {np.std(cF):.4f}, "
        f"피크수 {np.mean(cN):.1f} ± {np.std(cN):.1f}")
    log(f"  격리:    피크높이 {np.mean(iH):.1f} ± {np.std(iH):.1f}, "
        f"FWHM {np.mean(iF):.4f} ± {np.std(iF):.4f}, "
        f"피크수 {np.mean(iN):.1f} ± {np.std(iN):.1f}")

    # Mann-Whitney U (표본 작지만 시도)
    if len(cH) >= 3 and len(iH) >= 3:
        try:
            _, p_h = stats.mannwhitneyu(cH, iH, alternative='two-sided')
            _, p_f = stats.mannwhitneyu(cF, iF, alternative='two-sided')
            log(f"  피크높이 Mann-Whitney p = {p_h:.4f}")
            log(f"  FWHM Mann-Whitney p = {p_f:.4f}")
        except Exception as e:
            log(f"  Mann-Whitney 실패: {e}")

    return close_results, iso_results


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3단계: |F₂| 수렴 속도 — Hadamard 부분합
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def experiment_3(zeros, closest_5, isolated_5):
    """
    |F₂| 정의: ξ'/ξ(s) = B + Σ_ρ [1/(s-ρ) + 1/ρ]  (Hadamard 전개)
    L_K(s) = Σ_{k=1}^{K} [1/(s-ρ_k) + 1/(s-ρ̄_k)]  (부분합)
    |F₂|(t,K) = |ξ'/ξ(s) - B - L_K(s)|  (잔차)
    B는 보정 상수로, 큰 K에서 경험적으로 추정.
    """
    log("\n" + "=" * 70)
    log("3단계: |F₂| 수렴 속도 (Hadamard 부분합)")
    log("=" * 70)

    K_values = [16, 32, 64, 128, 256]

    # Hadamard 부분합용 영점 수집 (넉넉히)
    all_zeros = find_zeros_zeta(0, 600)
    n_all = len(all_zeros)
    log(f"  Hadamard 계산용 영점 {n_all}개 수집 (t<600)")

    if n_all < 256:
        log(f"  ⚠️ 영점 부족 ({n_all} < 256)")

    def hadamard_partial(s, K):
        """첫 K개 영점 쌍의 Hadamard 부분합"""
        L_K = mpmath.mpc(0, 0)
        K_use = min(K, n_all)
        for k in range(K_use):
            rho = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(all_zeros[k]))
            rho_bar = mpmath.mpf('0.5') - 1j * mpmath.mpf(str(all_zeros[k]))
            d1 = s - rho
            d2 = s - rho_bar
            # 영점에 너무 가까우면 건너뜀
            if abs(d1) < mpmath.mpf('0.01') or abs(d2) < mpmath.mpf('0.01'):
                continue
            L_K += 1 / d1 + 1 / d2
        return L_K

    # B 추정: 안전한 점 s=0.5+200i에서 큰 K로 보정
    log("\n  B (보정 상수) 추정...")
    s_calib = mpmath.mpf('0.5') + 1j * mpmath.mpf('200')
    L_exact_calib = connection_zeta(s_calib)
    L_K_calib = hadamard_partial(s_calib, min(n_all, 300))
    B_est = L_exact_calib - L_K_calib
    log(f"  B ≈ {float(mpmath.re(B_est)):.6f} + {float(mpmath.im(B_est)):.6f}i")

    # 테스트 영점들
    test_indices = list(closest_5) + list(isolated_5)
    test_labels = ['close'] * len(closest_5) + ['isolated'] * len(isolated_5)
    delta = 0.05  # 영점 오프셋

    results = {}
    close_slopes = []
    iso_slopes = []

    for idx, label in zip(test_indices, test_labels):
        t0 = zeros[idx]
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t0 + delta))
        L_exact = connection_zeta(s)
        key = f"{label}_n={idx}_t={t0:.2f}"
        log(f"\n  {key}")

        residuals = []
        for K in K_values:
            try:
                L_K = hadamard_partial(s, K)
                r = float(abs(L_exact - B_est - L_K))
                residuals.append(r)
                log(f"    K={K:3d}: |F₂| = {r:.6e}")
            except Exception as e:
                log(f"    K={K:3d}: 실패 — {e}")
                residuals.append(np.nan)

        # 수렴 기울기 (log-log)
        r_arr = np.array(residuals)
        valid = ~np.isnan(r_arr) & (r_arr > 0)
        slope = np.nan
        if np.sum(valid) >= 3:
            log_k = np.log10(np.array(K_values, dtype=float)[valid])
            log_r = np.log10(r_arr[valid])
            res = stats.linregress(log_k, log_r)
            slope = res.slope
            log(f"    수렴 기울기: {slope:.3f} (R²={res.rvalue**2:.3f})")
            if label == 'close':
                close_slopes.append(slope)
            else:
                iso_slopes.append(slope)

        results[key] = {
            't': t0, 'label': label, 'residuals': residuals, 'slope': slope,
        }

    # 비교
    log("\n--- |F₂| 수렴 속도 비교 ---")
    if close_slopes and iso_slopes:
        log(f"  근접쌍 평균 기울기: {np.mean(close_slopes):.3f} ± {np.std(close_slopes):.3f}")
        log(f"  격리   평균 기울기: {np.mean(iso_slopes):.3f} ± {np.std(iso_slopes):.3f}")
        log(f"  (더 큰 음수 = 더 빠른 수렴)")
        if len(close_slopes) >= 3 and len(iso_slopes) >= 3:
            try:
                _, p_slope = stats.mannwhitneyu(close_slopes, iso_slopes, alternative='two-sided')
                log(f"  기울기 Mann-Whitney p = {p_slope:.4f}")
            except Exception as e:
                log(f"  기울기 비교 실패: {e}")
    else:
        log("  데이터 부족으로 비교 불가")

    return results


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4단계: 전체 영점 상관 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def experiment_4(zeros, g_norm, valid_idx):
    """80개 영점에서 g_norm vs κ 특성 Spearman 상관"""
    log("\n" + "=" * 70)
    log("4단계: 상관 분석 (g_norm vs κ 특성, 전체 유효 영점)")
    log("=" * 70)

    g_list = []
    peak_h_list = []
    fwhm_list = []
    delta = 0.03

    n_total = len(valid_idx)
    for count, idx in enumerate(valid_idx):
        t0 = zeros[idx]
        g = g_norm[idx]
        if np.isnan(g):
            continue

        # 효율: 영점 근방 좁은 창에서 소수 점만 계산
        # ±0.5 범위, 30점 → 영점당 ~30 × 3 = 90 mpmath 호출
        try:
            ts_local, kappas_local = compute_kappa_profile(
                t0, half_width=0.5, n_points=30, delta=delta,
                avoid_zeros=[t0]
            )
            chars = extract_profile_chars(ts_local, kappas_local)
            g_list.append(g)
            peak_h_list.append(chars['peak_height'])
            fwhm_list.append(chars['fwhm'])

            if (count + 1) % 10 == 0 or count == 0:
                log(f"  [{count+1}/{n_total}] t={t0:.2f}, g={g:.4f}, "
                    f"피크높이={chars['peak_height']:.1f}, FWHM={chars['fwhm']:.4f}")
        except Exception as e:
            log(f"  [{count+1}/{n_total}] t={t0:.2f} 실패: {e}")

    g_arr = np.array(g_list)
    ph_arr = np.array(peak_h_list)
    fw_arr = np.array(fwhm_list)
    log(f"\n  유효 데이터: {len(g_arr)}개 / 전체 {n_total}개")

    if len(g_arr) < 5:
        log("  ⚠️ 데이터 부족으로 상관 분석 불가")
        return None

    # Spearman 상관
    rho_h, p_h = stats.spearmanr(g_arr, ph_arr)
    rho_f, p_f = stats.spearmanr(g_arr, fw_arr)
    log(f"\n  g_norm vs 피크높이: Spearman ρ = {rho_h:.4f}, p = {p_h:.4e}")
    log(f"  g_norm vs FWHM:     Spearman ρ = {rho_f:.4f}, p = {p_f:.4e}")

    # 추가: g_norm vs 1/피크높이 (반비례 가설)
    inv_ph = 1.0 / (ph_arr + 1.0)  # +1 안전장치
    rho_inv, p_inv = stats.spearmanr(g_arr, inv_ph)
    log(f"  g_norm vs 1/(피크높이+1): Spearman ρ = {rho_inv:.4f}, p = {p_inv:.4e}")

    # 사분위별 통계
    q25 = np.percentile(g_arr, 25)
    q75 = np.percentile(g_arr, 75)
    close_mask = g_arr <= q25
    iso_mask = g_arr >= q75
    if np.sum(close_mask) >= 3 and np.sum(iso_mask) >= 3:
        log(f"\n  하위 25% (g≤{q25:.3f}, n={np.sum(close_mask)}):")
        log(f"    피크높이: {np.mean(ph_arr[close_mask]):.1f} ± {np.std(ph_arr[close_mask]):.1f}")
        log(f"    FWHM:     {np.mean(fw_arr[close_mask]):.4f} ± {np.std(fw_arr[close_mask]):.4f}")
        log(f"  상위 25% (g≥{q75:.3f}, n={np.sum(iso_mask)}):")
        log(f"    피크높이: {np.mean(ph_arr[iso_mask]):.1f} ± {np.std(ph_arr[iso_mask]):.1f}")
        log(f"    FWHM:     {np.mean(fw_arr[iso_mask]):.4f} ± {np.std(fw_arr[iso_mask]):.4f}")

    return {
        'g_norms': g_arr, 'peak_heights': ph_arr, 'fwhms': fw_arr,
        'rho_height': rho_h, 'p_height': p_h,
        'rho_fwhm': rho_f, 'p_fwhm': p_f,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()
    log(f"Lehmer 쌍 곡률 분석 시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"mpmath.mp.dps = {mpmath.mp.dps}")
    log(f"설정: 프로파일 200점, 상관분석 30점/영점, δ=0.03-0.05")
    log("")

    # === 1단계 ===
    result = collect_zeros_and_spacings(200.0)
    if result[0] is None:
        log("영점 수집 실패 — 중단")
        return
    zeros, spacings, g_norm, closest_5, isolated_5, valid_idx = result

    # === 2단계 ===
    t2_start = time.time()
    close_res, iso_res = experiment_2(zeros, spacings, closest_5, isolated_5)
    log(f"\n  [2단계 소요: {time.time()-t2_start:.0f}초]")

    # === 3단계 ===
    t3_start = time.time()
    f2_results = experiment_3(zeros, closest_5, isolated_5)
    log(f"\n  [3단계 소요: {time.time()-t3_start:.0f}초]")

    # === 4단계 ===
    t4_start = time.time()
    corr = experiment_4(zeros, g_norm, valid_idx)
    log(f"\n  [4단계 소요: {time.time()-t4_start:.0f}초]")

    # === 최종 판정 ===
    elapsed = time.time() - t_start
    log("\n" + "=" * 70)
    log("최종 판정")
    log("=" * 70)

    verdicts = []

    # 기준 1: 프로파일 정성적 차이
    cN = [r['chars']['n_peaks'] for r in close_res]
    iN = [r['chars']['n_peaks'] for r in iso_res]
    cF = [r['chars']['fwhm'] for r in close_res]
    iF = [r['chars']['fwhm'] for r in iso_res]
    cH = [r['chars']['peak_height'] for r in close_res]
    iH = [r['chars']['peak_height'] for r in iso_res]

    profile_diff = False
    reasons = []
    if np.mean(cN) > np.mean(iN) + 0.3:
        profile_diff = True
        reasons.append(f"피크수 근접={np.mean(cN):.1f} > 격리={np.mean(iN):.1f}")
    if np.mean(cF) > np.mean(iF) * 1.3:
        profile_diff = True
        reasons.append(f"FWHM 근접={np.mean(cF):.4f} > 격리×1.3={np.mean(iF)*1.3:.4f}")
    if np.mean(cH) > np.mean(iH) * 1.5:
        profile_diff = True
        reasons.append(f"피크높이 근접={np.mean(cH):.1f} > 격리×1.5={np.mean(iH)*1.5:.1f}")

    if profile_diff:
        verdicts.append(f"양성: 프로파일 정성적 차이 [{'; '.join(reasons)}]")
    else:
        verdicts.append(f"음성: 프로파일 유사 (근접/격리 차이 미약)")

    # 기준 2: Spearman |ρ| > 0.5
    if corr is not None:
        rh = corr['rho_height']
        rf = corr['rho_fwhm']
        if abs(rh) > 0.5 or abs(rf) > 0.5:
            verdicts.append(f"양성: Spearman 상관 유의 (ρ_h={rh:.3f}, ρ_f={rf:.3f})")
        elif abs(rh) > 0.3 or abs(rf) > 0.3:
            verdicts.append(f"중립: Spearman 약한 상관 (ρ_h={rh:.3f}, ρ_f={rf:.3f})")
        else:
            verdicts.append(f"음성: Spearman 상관 약함 (ρ_h={rh:.3f}, ρ_f={rf:.3f})")
    else:
        verdicts.append("미결정: 상관 분석 데이터 부족")

    # 기준 3: |F₂| 수렴 속도 차이
    close_slopes = [d['slope'] for d in f2_results.values()
                    if d['label'] == 'close' and not np.isnan(d['slope'])]
    iso_slopes = [d['slope'] for d in f2_results.values()
                  if d['label'] == 'isolated' and not np.isnan(d['slope'])]
    if close_slopes and iso_slopes:
        cs_mean = np.mean(close_slopes)
        is_mean = np.mean(iso_slopes)
        # 근접쌍이 더 느린 수렴 = 기울기가 덜 음수 (더 큰 값)
        if cs_mean > is_mean + 0.2:
            verdicts.append(f"양성: |F₂| 수렴 근접쌍 느림 ({cs_mean:.3f} vs {is_mean:.3f})")
        elif abs(cs_mean - is_mean) < 0.2:
            verdicts.append(f"중립: |F₂| 수렴 속도 유사 ({cs_mean:.3f} vs {is_mean:.3f})")
        else:
            verdicts.append(f"음성: |F₂| 수렴 근접쌍이 오히려 빠름 ({cs_mean:.3f} vs {is_mean:.3f})")
    else:
        verdicts.append("미결정: |F₂| 기울기 데이터 부족")

    log("\n  개별 판정:")
    for v in verdicts:
        log(f"    {v}")

    pos = sum(1 for v in verdicts if v.startswith("양성"))
    neg = sum(1 for v in verdicts if v.startswith("음성"))
    neu = sum(1 for v in verdicts if v.startswith("중립"))
    total = len(verdicts)

    if pos >= 2:
        final = f"양성 ({pos}/{total} 기준 충족)"
    elif pos == 0 and neg >= 2:
        final = f"음성 ({neg}/{total} 기준 불충족)"
    else:
        final = f"중립 (양성 {pos}, 음성 {neg}, 중립 {neu} / {total})"

    log(f"\n  ★ 최종: {final}")
    log(f"\n  총 소요 시간: {elapsed:.1f}초 ({elapsed/60:.1f}분)")
    log(f"  완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == '__main__':
    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    f = open(OUT, 'w', encoding='utf-8')
    sys.stdout = Tee(sys.__stdout__, f)

    try:
        main()
    except Exception as e:
        import traceback
        log(f"\n!!! 치명적 오류: {e}")
        traceback.print_exc()
    finally:
        f.close()
        sys.stdout = sys.__stdout__
        print(f"\n결과 저장: {OUT}")
