#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lehmer κ 고해상도 + δ-강건성 검증
====================================
사이클 #23 수학자 지시:
  - 4단계 FWHM 해상도 수정: half_width=0.3, n_points=200 (step=0.003)
  - δ∈{0.01, 0.03, 0.05} 강건성: 각 δ에서 g_norm vs κ 피크높이 Spearman ρ
  - g_norm vs FWHM 상관 측정 (고해상도)
  - 결과 파일: results/lehmer_highres_robustness.txt

사이클 #22에서 확인된 한계:
  - 4단계 FWHM=0 (step=0.033 > 피크폭 ~0.06) → 고해상도 재실행 필요
  - δ-강건성 미검증 → δ∈{0.01, 0.03, 0.05} 각각 ρ 측정
  - 결과 #16 "조건부 양성" → 두 기준 충족 시 "양성"으로 격상 가능
"""

import sys
import os
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import mpmath
from scipy import stats

from bundle_utils import (
    xi_func, connection_zeta,
    find_zeros_zeta, _near_zero,
)

mpmath.mp.dps = 100  # 고정밀도 (δ=0.01 근접 측정 대비)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(SCRIPT_DIR, '..', 'results', 'lehmer_highres_robustness.txt')


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
# 단계 1: 영점 수집 + 정규화 간격
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def collect_zeros_and_spacings(t_max=200.0):
    """영점 수집 + 정규화 간격 g_norm 계산 (사이클 #22와 동일)"""
    log("=" * 70)
    log("단계 1: 영점 수집 및 정규화 간격")
    log("=" * 70)

    zeros = find_zeros_zeta(0, t_max)
    N = len(zeros)
    log(f"  t in [0, {t_max}] 영점 {N}개 수집")

    if N < 10:
        log("  ⚠️ 영점 0개 또는 부족 — 탐색 로직 점검 필요")
        return None

    # 연속 영점 간격 Δt_n = t_{n+1} - t_n
    spacings = np.diff(zeros)

    # 정규화: g_norm = Δt / <Δt>_local, 여기서 <Δt> = 2π / log(t/(2π))
    t_mid = (zeros[:-1] + zeros[1:]) / 2.0
    valid_mask = t_mid > 2 * np.pi
    mean_spacing = np.full_like(spacings, np.nan)
    mean_spacing[valid_mask] = 2.0 * np.pi / np.log(t_mid[valid_mask] / (2.0 * np.pi))

    g_norm = np.full_like(spacings, np.nan)
    g_norm[valid_mask] = spacings[valid_mask] / mean_spacing[valid_mask]

    valid_idx = np.where(valid_mask)[0]
    g_valid = g_norm[valid_idx]
    log(f"  유효 간격 수: {len(valid_idx)}")
    log(f"  g_norm: min={np.nanmin(g_valid):.4f}, max={np.nanmax(g_valid):.4f}, "
        f"mean={np.nanmean(g_valid):.4f}, std={np.nanstd(g_valid):.4f}")

    return zeros, spacings, g_norm, valid_idx


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# κ 측정 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_kappa_at(t_eval):
    """κ(0.5 + it_eval) 단일점 측정. 실패 시 None 반환."""
    try:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_eval))
        xi_val = xi_func(s)
        if _near_zero(xi_val):
            return None  # 영점 위
        L = connection_zeta(s)
        kappa = float(abs(L) ** 2)
        if not np.isfinite(kappa):
            return None
        return kappa
    except Exception as e:
        log(f"    WARNING: κ({t_eval:.8f}) 실패: {e}")
        return None


def compute_kappa_profile_hr(t_center, zeros_nearby, half_width=0.3, n_points=200, avoid_delta=0.005):
    """
    고해상도 κ 프로파일 계산.
    t_center ± half_width, n_points=200 (step≈0.003)
    zeros_nearby: 이 범위 안의 모든 영점 리스트 (오프셋 회피용)
    avoid_delta: 영점으로부터 이 거리 이내면 zeros_nearby의 점을 + avoid_delta로 이동
    """
    ts = np.linspace(t_center - half_width, t_center + half_width, n_points)
    kappas = np.zeros(n_points)

    for i, t in enumerate(ts):
        t_eval = float(t)
        # 영점 근처 → 가장 가까운 영점 + avoid_delta로 이동
        for tz in zeros_nearby:
            if abs(t_eval - tz) < avoid_delta:
                t_eval = tz + avoid_delta
                break

        kappa = compute_kappa_at(t_eval)
        if kappa is None:
            kappas[i] = 0.0
        else:
            kappas[i] = min(kappa, 5e6)  # 스파이크 cap (1e8 대신 더 낮게)

    return ts, kappas


def extract_fwhm(ts, kappas):
    """
    κ 프로파일에서 FWHM 추출.
    전역 최대값 기반: 평탄 구간(plateau)에도 대응.
    반환: (peak_height, fwhm, n_peaks)
    """
    # 유효값 필터
    valid = np.isfinite(kappas) & (kappas > 0)
    if np.sum(valid) < 5:
        return 0.0, 0.0, 0

    kappas_c = kappas.copy()
    kappas_c[~valid] = 0.0

    peak_height = float(np.max(kappas_c))
    if peak_height <= 0:
        return 0.0, 0.0, 0

    step = float(ts[1] - ts[0])

    # FWHM: 전역 최대값의 절반 수준을 넘는 연속 구간의 총 폭
    half_max = peak_height / 2.0
    above = kappas_c >= half_max

    # 가장 넓은 연속 구간 찾기 (전역 최대값 포함)
    best_max_idx = int(np.argmax(kappas_c))
    # 최대값 위치에서 좌우로 half_max 기준선 탐색
    left = best_max_idx
    while left > 0 and kappas_c[left - 1] >= half_max:
        left -= 1
    right = best_max_idx
    while right < len(kappas_c) - 1 and kappas_c[right + 1] >= half_max:
        right += 1

    fwhm = float((right - left)) * step

    # 피크 수: 배경 대비 5배 이상인 로컬 최대값 수
    background = float(np.median(kappas_c[valid]))
    threshold = max(background * 5.0, peak_height * 0.1)
    n_peaks = 0
    for i in range(1, len(kappas_c) - 1):
        if (kappas_c[i] >= kappas_c[i - 1] and
                kappas_c[i] > kappas_c[i + 1] and
                kappas_c[i] > threshold):
            n_peaks += 1
        elif (kappas_c[i] > kappas_c[i - 1] and
              kappas_c[i] >= kappas_c[i + 1] and
              kappas_c[i] > threshold):
            n_peaks += 1

    return peak_height, fwhm, n_peaks


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 단계 2: 고해상도 FWHM (전체 78개 영점)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def stage2_highres_fwhm(zeros, g_norm, valid_idx):
    """
    78개 영점 각각에 대해 고해상도 FWHM 측정.
    half_width=0.3, n_points=200 (step=0.003), avoid_delta=0.005
    """
    log("\n" + "=" * 70)
    log("단계 2: 고해상도 FWHM (half_width=0.3, n_points=200, step=0.003)")
    log("=" * 70)
    log("  (모든 영점에 동일 창 적용 — 창 비대칭 문제 회피)")
    log("")

    HALF_WIDTH = 0.3
    N_POINTS = 200
    AVOID_DELTA = 0.02  # 영점 회피 거리: FWHM ≈ 2√2 × 0.02 ≈ 0.057 ("~0.06" 대응)

    # 창 내 영점 사전 색인 (속도 최적화)
    zero_arr = np.array(zeros)

    g_list = []
    fwhm_list = []
    peak_h_list = []
    n_peaks_list = []
    n_total = len(valid_idx)

    t2_start = time.time()

    for count, idx in enumerate(valid_idx):
        t0 = float(zeros[idx])
        g = float(g_norm[idx])
        if np.isnan(g):
            continue

        # 창 내 영점 목록
        t_lo = t0 - HALF_WIDTH
        t_hi = t0 + HALF_WIDTH
        nearby = [float(z) for z in zero_arr if t_lo <= z <= t_hi]

        try:
            ts, kappas = compute_kappa_profile_hr(
                t0, zeros_nearby=nearby,
                half_width=HALF_WIDTH, n_points=N_POINTS,
                avoid_delta=AVOID_DELTA,
            )
            ph, fw, np_ = extract_fwhm(ts, kappas)
            g_list.append(g)
            peak_h_list.append(ph)
            fwhm_list.append(fw)
            n_peaks_list.append(np_)

            if (count + 1) % 10 == 0 or count == 0:
                elapsed = time.time() - t2_start
                eta = elapsed / (count + 1) * (n_total - count - 1)
                log(f"  [{count+1:3d}/{n_total}] t={t0:.3f}, g={g:.4f}, "
                    f"피크높이={ph:.1f}, FWHM={fw:.5f}, 피크수={np_}  "
                    f"({elapsed:.0f}s 경과, ETA {eta:.0f}s)")
        except Exception as e:
            log(f"  [{count+1:3d}/{n_total}] t={t0:.3f} 실패: {e}")

    g_arr = np.array(g_list)
    fwhm_arr = np.array(fwhm_list)
    ph_arr = np.array(peak_h_list)

    log(f"\n  유효 데이터: {len(g_arr)}개 / 전체 {n_total}개")
    log(f"  FWHM: min={np.min(fwhm_arr):.5f}, max={np.max(fwhm_arr):.5f}, "
        f"mean={np.mean(fwhm_arr):.5f}, std={np.std(fwhm_arr):.5f}")
    log(f"  FWHM=0인 영점 수: {np.sum(fwhm_arr == 0.0)}")
    log(f"  [단계 2 소요: {time.time()-t2_start:.0f}초]")

    if len(g_arr) < 5:
        log("  ⚠️ 데이터 부족 — 상관 분석 불가")
        return None

    # 상관 분석
    log("\n--- 고해상도 FWHM 상관 분석 ---")

    # FWHM > 0인 데이터만 (FWHM=0이면 해상도 여전히 불충분)
    valid_fwhm = fwhm_arr > 0
    log(f"  FWHM > 0인 영점: {np.sum(valid_fwhm)}개")

    if np.sum(valid_fwhm) >= 5:
        rho_f, p_f = stats.spearmanr(g_arr[valid_fwhm], fwhm_arr[valid_fwhm])
        log(f"  g_norm vs FWHM (유효 데이터): Spearman ρ = {rho_f:.4f}, p = {p_f:.4e}")
    else:
        rho_f, p_f = np.nan, np.nan
        log("  ⚠️ FWHM > 0 데이터 부족")

    # 전체 데이터 상관 (FWHM=0 포함)
    rho_all, p_all = stats.spearmanr(g_arr, fwhm_arr)
    log(f"  g_norm vs FWHM (전체):    Spearman ρ = {rho_all:.4f}, p = {p_all:.4e}")

    # 피크높이 상관 (비교)
    rho_h, p_h = stats.spearmanr(g_arr, ph_arr)
    log(f"  g_norm vs 피크높이:       Spearman ρ = {rho_h:.4f}, p = {p_h:.4e}")

    # 사분위별
    q25 = np.percentile(g_arr, 25)
    q75 = np.percentile(g_arr, 75)
    close_m = g_arr <= q25
    iso_m = g_arr >= q75
    log(f"\n  하위 25% (g≤{q25:.3f}, n={np.sum(close_m)}):")
    log(f"    FWHM:     {np.mean(fwhm_arr[close_m]):.5f} ± {np.std(fwhm_arr[close_m]):.5f}")
    log(f"    피크높이: {np.mean(ph_arr[close_m]):.1f} ± {np.std(ph_arr[close_m]):.1f}")
    log(f"  상위 25% (g≥{q75:.3f}, n={np.sum(iso_m)}):")
    log(f"    FWHM:     {np.mean(fwhm_arr[iso_m]):.5f} ± {np.std(fwhm_arr[iso_m]):.5f}")
    log(f"    피크높이: {np.mean(ph_arr[iso_m]):.1f} ± {np.std(ph_arr[iso_m]):.1f}")

    # 성공 기준 판정
    log("\n  ★ 성공 기준: |ρ(g_norm, FWHM)| > 0.5")
    if not np.isnan(rho_f):
        if abs(rho_f) > 0.5:
            log(f"  ✅ FWHM 상관 PASS: |ρ| = {abs(rho_f):.4f} > 0.5")
        else:
            log(f"  ❌ FWHM 상관 FAIL: |ρ| = {abs(rho_f):.4f} ≤ 0.5")
    else:
        log("  ⚠️ FWHM 상관 측정 불가")

    return {
        'g_arr': g_arr,
        'fwhm_arr': fwhm_arr,
        'ph_arr': ph_arr,
        'rho_fwhm': rho_f,
        'p_fwhm': p_f,
        'rho_fwhm_all': rho_all,
        'rho_height': rho_h,
        'p_height': p_h,
        'valid_fwhm_count': int(np.sum(valid_fwhm)),
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 단계 3: δ-강건성 (δ∈{0.01, 0.03, 0.05})
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def stage3_delta_robustness(zeros, g_norm, valid_idx):
    """
    δ∈{0.01, 0.03, 0.05}에서 κ(t₀ + δ) 측정.
    g_norm vs κ 피크높이 Spearman ρ 계산.
    이상치 플래그: κ 밖 1/δ² ± 20%.
    """
    log("\n" + "=" * 70)
    log("단계 3: δ-강건성 검증 (δ∈{0.01, 0.03, 0.05})")
    log("=" * 70)
    log("  이상치 기준: κ(t₀+δ) 밖 [1/δ² × 0.8, 1/δ² × 1.2] → 플래그")
    log("")

    DELTAS = [0.01, 0.03, 0.05]
    n_total = len(valid_idx)

    # 각 δ별 결과 저장
    results_by_delta = {d: {'g': [], 'kappa': [], 'outlier_flags': []} for d in DELTAS}

    t3_start = time.time()

    for count, idx in enumerate(valid_idx):
        t0 = float(zeros[idx])
        g = float(g_norm[idx])
        if np.isnan(g):
            continue

        for delta in DELTAS:
            t_eval = t0 + delta
            kappa = compute_kappa_at(t_eval)

            if kappa is None:
                log(f"  WARNING: κ({t0:.3f}+{delta}) = None — 근처 영점 충돌?")
                continue

            # 이상치 플래그: 1/δ² ± 20%
            expected = 1.0 / (delta ** 2)
            outlier = not (0.8 * expected <= kappa <= 1.2 * expected)
            if outlier:
                # 심각한 이상치만 로그 (>10× 차이)
                ratio = kappa / expected
                if ratio > 10 or ratio < 0.1:
                    log(f"  ⚠️ 이상치: t₀={t0:.3f}, δ={delta:.2f}, "
                        f"κ={kappa:.2e}, 예측={expected:.2e} (비율={ratio:.2f})")

            results_by_delta[delta]['g'].append(g)
            results_by_delta[delta]['kappa'].append(kappa)
            results_by_delta[delta]['outlier_flags'].append(outlier)

    log(f"  [측정 완료: {time.time()-t3_start:.1f}초]\n")

    # 각 δ별 상관 분석
    delta_summary = {}

    for delta in DELTAS:
        g_arr = np.array(results_by_delta[delta]['g'])
        k_arr = np.array(results_by_delta[delta]['kappa'])
        flags = np.array(results_by_delta[delta]['outlier_flags'])

        n_outlier = np.sum(flags)
        expected = 1.0 / (delta ** 2)

        log(f"  δ = {delta:.2f}:")
        log(f"    측정 성공: {len(g_arr)}개, 이상치: {n_outlier}개 ({100*n_outlier/max(len(g_arr),1):.1f}%)")
        log(f"    κ 통계: min={np.min(k_arr):.2e}, max={np.max(k_arr):.2e}, "
            f"mean={np.mean(k_arr):.2e} (1/δ²={expected:.2e})")

        if len(g_arr) < 5:
            log(f"    ⚠️ 데이터 부족 — 상관 불가")
            delta_summary[delta] = {'rho': np.nan, 'p': np.nan, 'n': len(g_arr)}
            continue

        rho, p = stats.spearmanr(g_arr, k_arr)
        log(f"    g_norm vs κ: Spearman ρ = {rho:.4f}, p = {p:.4e}")

        # 이상치 제거 후 재검증
        clean_g = g_arr[~flags]
        clean_k = k_arr[~flags]
        if len(clean_g) >= 5:
            rho_c, p_c = stats.spearmanr(clean_g, clean_k)
            log(f"    이상치 제거 후 (n={len(clean_g)}): ρ = {rho_c:.4f}, p = {p_c:.4e}")
        else:
            rho_c, p_c = np.nan, np.nan

        # 성공 기준
        if abs(rho) > 0.7:
            log(f"    ✅ δ={delta:.2f} PASS: |ρ| = {abs(rho):.4f} > 0.7")
        else:
            log(f"    ❌ δ={delta:.2f} FAIL: |ρ| = {abs(rho):.4f} ≤ 0.7")

        delta_summary[delta] = {
            'rho': rho, 'p': p, 'n': len(g_arr),
            'rho_clean': rho_c, 'p_clean': p_c, 'n_clean': len(clean_g),
            'n_outlier': int(n_outlier),
        }
        log("")

    # δ-강건성 종합 판정
    log("--- δ-강건성 종합 ---")
    pass_count = 0
    for delta in DELTAS:
        info = delta_summary.get(delta, {})
        rho = info.get('rho', np.nan)
        if not np.isnan(rho) and abs(rho) > 0.7:
            pass_count += 1
            log(f"  δ={delta:.2f}: PASS (|ρ|={abs(rho):.4f})")
        else:
            log(f"  δ={delta:.2f}: FAIL (|ρ|={abs(rho) if not np.isnan(rho) else 'nan':.4f})")

    if pass_count == 3:
        log(f"  ✅ δ-강건성 PASS: 3/3 통과 → 결과 #16 격상 조건 충족")
    elif pass_count >= 2:
        log(f"  ⚠️ δ-강건성 부분 통과: {pass_count}/3 → 추가 분석 필요")
    else:
        log(f"  ❌ δ-강건성 FAIL: {pass_count}/3 → 결과 #16 조건부 유지")

    return delta_summary


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 최종 요약
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def print_final_summary(stage2_res, stage3_res):
    log("\n" + "=" * 70)
    log("최종 요약 — 결과 #16 격상 기준")
    log("=" * 70)

    # 기준 1: g_norm vs FWHM |ρ| > 0.5
    log("\n[기준 1] g_norm vs FWHM |Spearman ρ| > 0.5")
    if stage2_res is not None:
        rho_f = stage2_res.get('rho_fwhm', np.nan)
        p_f = stage2_res.get('p_fwhm', np.nan)
        n_valid = stage2_res.get('valid_fwhm_count', 0)
        if not np.isnan(rho_f):
            result_str = "✅ PASS" if abs(rho_f) > 0.5 else "❌ FAIL"
            log(f"  {result_str}: ρ = {rho_f:.4f} (p={p_f:.4e}, n_유효={n_valid})")
        else:
            log("  ⚠️ 측정 불가 (FWHM > 0 데이터 부족)")
    else:
        log("  ⚠️ 단계 2 결과 없음")

    # 기준 2: 3개 δ 모두 g_norm vs κ |ρ| > 0.7
    log("\n[기준 2] δ 3개 모두 g_norm vs κ |Spearman ρ| > 0.7")
    if stage3_res:
        pass_all = True
        for d in [0.01, 0.03, 0.05]:
            info = stage3_res.get(d, {})
            rho = info.get('rho', np.nan)
            p = info.get('p', np.nan)
            if not np.isnan(rho) and abs(rho) > 0.7:
                log(f"  ✅ δ={d:.2f}: ρ={rho:.4f} (p={p:.4e})")
            else:
                log(f"  ❌ δ={d:.2f}: ρ={rho if not np.isnan(rho) else 'nan'}")
                pass_all = False
        if pass_all:
            log("  → 기준 2 PASS")
        else:
            log("  → 기준 2 FAIL (일부 δ 미충족)")
    else:
        log("  ⚠️ 단계 3 결과 없음")

    # 최종 판정
    log("\n[최종 판정]")

    crit1 = False
    crit2 = False

    if stage2_res and not np.isnan(stage2_res.get('rho_fwhm', np.nan)):
        crit1 = abs(stage2_res['rho_fwhm']) > 0.5

    if stage3_res:
        crit2 = all(
            not np.isnan(stage3_res.get(d, {}).get('rho', np.nan)) and
            abs(stage3_res.get(d, {}).get('rho', np.nan)) > 0.7
            for d in [0.01, 0.03, 0.05]
        )

    if crit1 and crit2:
        log("  🏆 기준 1+2 모두 충족 → 결과 #16 '양성'으로 격상")
    elif crit1:
        log("  ⚠️ 기준 1만 충족 (FWHM 상관) → 결과 #16 '조건부 양성' 유지 (δ-강건성 미확인)")
    elif crit2:
        log("  ⚠️ 기준 2만 충족 (δ-강건성) → 결과 #16 '조건부 양성' 유지 (FWHM 상관 약함)")
    else:
        log("  ❌ 기준 1+2 모두 미충족 → 결과 #16 '조건부 양성' 유지, 추가 추구 불필요")

    log("\n  사이클 #22 기존 결과 (참고):")
    log("  - FWHM 3.4× (p=0.0097): ✅ 근접쌍 vs 격리 (5 vs 5)")
    log("  - ρ=0.835 (p=2e-21): ✅ g_norm vs κ 피크높이 (δ=0.03, 78개)")
    log("  - |F₂| 수렴: ❌ 예측 반증 (논문 수정 필요)")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    out_dir = os.path.dirname(OUT)
    os.makedirs(out_dir, exist_ok=True)
    f_out = open(OUT, 'w', encoding='utf-8')
    sys.stdout = Tee(sys.__stdout__, f_out)

    log(f"Lehmer κ 고해상도 + δ-강건성 검증 시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"mpmath.mp.dps = {mpmath.mp.dps}")
    log(f"설정: FWHM half_width=0.3, n_points=200 (step=0.003), δ={{0.01,0.03,0.05}}")
    log(f"출력: {OUT}")
    log("")

    # === 단계 1: 영점 수집 ===
    res1 = collect_zeros_and_spacings(200.0)
    if res1 is None:
        log("영점 수집 실패 — 중단")
        f_out.close()
        return
    zeros, spacings, g_norm, valid_idx = res1

    # === 단계 2: 고해상도 FWHM ===
    t2_start = time.time()
    stage2_res = stage2_highres_fwhm(zeros, g_norm, valid_idx)
    log(f"\n  [단계 2 총 소요: {time.time()-t2_start:.0f}초]")

    # === 단계 3: δ-강건성 ===
    t3_start = time.time()
    stage3_res = stage3_delta_robustness(zeros, g_norm, valid_idx)
    log(f"\n  [단계 3 총 소요: {time.time()-t3_start:.0f}초]")

    # === 최종 요약 ===
    print_final_summary(stage2_res, stage3_res)

    total = time.time() - t_start
    log(f"\n총 소요: {total:.1f}초 ({total/60:.1f}분)")
    log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"결과 파일: {OUT}")

    f_out.close()


if __name__ == '__main__':
    main()
