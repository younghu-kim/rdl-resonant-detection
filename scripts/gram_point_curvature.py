#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Gram 점 곡률 분석 — Gram's law 위반의 기하학적 서명
====================================================
사이클 #24 수학자 지시:
  - Gram 점에서 Good/Bad 분류 후 κ 비교
  - Mann-Whitney U 검정: Good vs Bad κ
  - Spearman ρ(κ, d_zero)
  - 정성적 κ 프로파일 비교 (5개씩)

수정 사항 (실행 중 발견 문제):
  - dps=50 → t>100 구간 ξ 계산 정밀도 부족 (→ dps=80)
  - n_max=126 → 첫 bad Gram 점이 n=126으로 too few (→ n_max=200, bad≥3)
  - zeros_list t_max=200 → Gram 점 t>200 coverage 부족 (→ t_max=450)

결과 파일: results/gram_point_curvature.txt
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
)

# dps: t>100이면 ≥80 규칙 → 전체 80으로 통일
mpmath.mp.dps = 80

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_PATH = os.path.join(SCRIPT_DIR, '..', 'results', 'gram_point_curvature.txt')
os.makedirs(os.path.join(SCRIPT_DIR, '..', 'results'), exist_ok=True)


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
# 1단계: Gram 점 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_gram_points(n_max=200):
    """
    Gram 점 g_n: θ(g_n) = n * π, n = 0 ... n_max
    θ(t) = mpmath.siegeltheta(t)

    참고: n=0 → g_0 ≈ 17.85 (θ(17.85)=0)
         n=126 → g_126 ≈ 282 (첫 bad Gram 점)
         n=134 → g_134 ≈ 295 (두 번째 bad Gram 점)
         n=195 → g_195 ≈ 400 (세 번째 bad Gram 점)
         → n_max=200까지: bad 3개 이상 예상
    """
    log("=" * 70)
    log("1단계: Gram 점 계산")
    log("=" * 70)
    log(f"  n=0 ~ {n_max} 범위 계산 중 ...")
    log("  [주의] 첫 bad Gram 점은 n=126 (수학적으로 알려진 사실)")

    gram_points = []
    fail_count = 0

    for n in range(n_max + 1):
        target = mpmath.mpf(n) * mpmath.pi

        if n == 0:
            t_init = mpmath.mpf('17.5')
        else:
            if len(gram_points) > 0:
                t_init = gram_points[-1][1] + mpmath.mpf('1.5')
            else:
                t_init = mpmath.mpf('17.5')

        try:
            g_n = mpmath.findroot(
                lambda t: mpmath.siegeltheta(t) - target,
                t_init,
                tol=mpmath.mpf(10)**(-mpmath.mp.dps + 15)
            )
            g_n_f = float(g_n.real)
            gram_points.append((n, g_n_f))
        except Exception as e:
            fail_count += 1
            log(f"  WARNING: n={n} findroot 실패: {e}")

    log(f"  Gram 점 계산 완료: {len(gram_points)}개 (실패 {fail_count}개)")
    if len(gram_points) > 0:
        log(f"  범위: t={gram_points[0][1]:.3f} ~ t={gram_points[-1][1]:.3f}")

    if fail_count > len(gram_points) // 2:
        log("  ⚠️ 절반 이상 findroot 실패 — 탐색 로직 점검 필요")
        return []

    return gram_points


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2단계: Gram's law 분류
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def classify_gram_points(gram_points):
    """
    Good/Bad 분류:
      (-1)^n * Z(g_n) > 0  → "good"
      (-1)^n * Z(g_n) ≤ 0  → "bad"
    Z(t) = mpmath.siegelz(t)
    """
    log("=" * 70)
    log("2단계: Gram's law 분류 (Good/Bad)")
    log("=" * 70)

    classified = []  # (n, t, label, z_val)
    good_count = 0
    bad_count = 0

    for (n, t) in gram_points:
        try:
            z_val = float(mpmath.siegelz(mpmath.mpf(str(t))))
            sign = ((-1) ** n) * z_val
            label = "good" if sign > 0 else "bad"
            classified.append((n, t, label, z_val))
            if label == "good":
                good_count += 1
            else:
                bad_count += 1
        except Exception as e:
            log(f"  WARNING: n={n}, t={t:.3f} Z 계산 실패: {e}")

    total = good_count + bad_count
    if total == 0:
        log("  ⚠️ 분류 결과 없음")
        return []

    log(f"  Good: {good_count}개 ({100*good_count/total:.1f}%)")
    log(f"  Bad:  {bad_count}개 ({100*bad_count/total:.1f}%)")

    if bad_count < 8:
        log(f"  ⚠️ bad Gram 점 수 {bad_count} < 8 — 통계적 검정력 부족 경고")
        log("  (수학 문헌: bad Gram 점은 n=126, 134, 195, 211, 232, ...에 위치)")

    # 처음 20개 샘플 출력
    log("\n  처음 20개 Gram 점:")
    log(f"  {'n':>4} {'t':>10} {'Z(g_n)':>12} {'(-1)^n*Z':>12} {'분류':>6}")
    for (n, t, label, z_val) in classified[:20]:
        sign_val = ((-1)**n) * z_val
        log(f"  {n:>4} {t:>10.4f} {z_val:>12.6f} {sign_val:>12.6f} {label:>6}")

    # Bad 점 별도 출력
    bad_list = [(n, t, lab, z) for (n, t, lab, z) in classified if lab == "bad"]
    if bad_list:
        log(f"\n  Bad Gram 점 목록 ({len(bad_list)}개):")
        for (n, t, lab, z) in bad_list:
            sign_val = ((-1)**n) * z
            log(f"  n={n:>4}, t={t:>10.4f}, Z={z:>12.6f}, (-1)^n*Z={sign_val:>12.6f}")

    return classified


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3단계: 영점 수집
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def collect_zeros_for_range(t_max=450.0):
    """zetazero를 사용하여 t∈[0, t_max] 영점 수집"""
    zeros = []
    k = 1
    while True:
        try:
            z = float(mpmath.zetazero(k).imag)
            if z > t_max:
                break
            zeros.append(z)
            k += 1
        except Exception as e:
            log(f"  WARNING: zetazero({k}) 실패: {e}")
            break
    return sorted(zeros)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4단계: κ 측정 + d_zero 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def measure_kappa_and_dzero(classified, zeros_list, delta=0.03):
    """
    각 Gram 점에서:
    - κ = curvature_zeta(1/2 + i*(g_n + delta))
    - d_zero = min |g_n - z_k|

    dps=80으로 t>100 구간도 정상 계산 가능.
    """
    log("=" * 70)
    log(f"3단계: κ 측정 (δ={delta}) + d_zero 계산")
    log("=" * 70)
    log(f"  dps={mpmath.mp.dps} (t>100 구간 대응)")

    zeros_arr = np.array(zeros_list)
    results = []  # (n, t, label, kappa, d_zero, z_val)
    cap_count = 0  # κ = 1e20 cap에 걸린 횟수

    for i, (n, t, label, z_val) in enumerate(classified):
        # κ 측정: δ=0.03 오프셋 (영점 위 직접 측정 금지)
        s = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t + delta)))
        try:
            kappa = float(curvature_zeta(s))
            if not np.isfinite(kappa):
                log(f"  WARNING: n={n} κ 비유한값 — 건너뜀")
                continue
            if kappa >= 1e15:  # cap 감지 (1e20 미만도 포함)
                cap_count += 1
                # cap 히트는 포함하되 플래그
        except Exception as e:
            log(f"  WARNING: n={n}, t={t:.3f} κ 계산 실패: {e}")
            continue

        # d_zero: 가장 가까운 영점까지 거리
        if len(zeros_arr) > 0:
            dists = np.abs(zeros_arr - t)
            d_zero = float(np.min(dists))
        else:
            d_zero = float('nan')

        results.append((n, t, label, kappa, d_zero, z_val))

        if (i + 1) % 30 == 0:
            elapsed = time.time()
            log(f"  진행: {i+1}/{len(classified)} ({100*(i+1)/len(classified):.0f}%) | cap히트:{cap_count}개")

    good_results = [(n, t, lab, k, d, z) for (n, t, lab, k, d, z) in results if lab == "good"]
    bad_results  = [(n, t, lab, k, d, z) for (n, t, lab, k, d, z) in results if lab == "bad"]

    log(f"\n  κ 측정 완료: {len(results)}개 (cap히트: {cap_count}개)")
    log(f"  Good: {len(good_results)}개, Bad: {len(bad_results)}개")

    # cap 제외한 유효 κ 통계
    valid = [(n, t, lab, k, d, z) for (n, t, lab, k, d, z) in results if k < 1e15]
    if valid:
        vk = np.array([r[3] for r in valid])
        log(f"  유효 κ (cap 제외): n={len(vk)}, 중앙={np.median(vk):.4f}, 평균={np.mean(vk):.4f}, std={np.std(vk):.4f}")

    if len(good_results) > 0:
        good_kappa = [r[3] for r in good_results if r[3] < 1e15]
        if good_kappa:
            log(f"  Good κ (유효): n={len(good_kappa)}, 중앙={np.median(good_kappa):.4f}")
    if len(bad_results) > 0:
        bad_kappa = [r[3] for r in bad_results if r[3] < 1e15]
        if bad_kappa:
            log(f"  Bad  κ (유효): n={len(bad_kappa)}, 중앙={np.median(bad_kappa):.4f}")

    return results


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5단계: 상관 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def analyze_correlations(results):
    """
    (a) Good vs Bad κ Mann-Whitney U 검정 (cap 제외)
    (b) Good vs Bad d_zero 비교
    (c) κ vs d_zero Spearman ρ (전체 Gram 점, cap 제외)
    (d) Bad Gram 점 5개 vs Good Gram 점 5개 정성 비교
    """
    log("=" * 70)
    log("4단계: 상관 분석")
    log("=" * 70)

    # cap(1e15 이상) 제외
    valid_results = [(n, t, lab, k, d, z) for (n, t, lab, k, d, z) in results if k < 1e15]
    all_results   = results  # 전체 (cap 포함)

    good_valid = [(n, t, lab, k, d, z) for (n, t, lab, k, d, z) in valid_results if lab == "good"]
    bad_valid  = [(n, t, lab, k, d, z) for (n, t, lab, k, d, z) in valid_results if lab == "bad"]

    good_kappa  = np.array([r[3] for r in good_valid]) if good_valid else np.array([])
    bad_kappa   = np.array([r[3] for r in bad_valid]) if bad_valid else np.array([])
    good_dzero  = np.array([r[4] for r in good_valid]) if good_valid else np.array([])
    bad_dzero   = np.array([r[4] for r in bad_valid]) if bad_valid else np.array([])
    all_kappa   = np.array([r[3] for r in valid_results])
    all_dzero   = np.array([r[4] for r in valid_results])

    analysis_a = None
    analysis_b = None
    analysis_c = None

    # ──────────────────────────────────────────────────────────────────────
    # (a) Mann-Whitney U: Good vs Bad κ
    log("\n--- (a) Good vs Bad κ Mann-Whitney U 검정 (cap 제외) ---")
    log(f"  Good(유효): {len(good_kappa)}개, Bad(유효): {len(bad_kappa)}개")

    if len(bad_kappa) >= 3 and len(good_kappa) >= 3:
        u_stat, p_mw = stats.mannwhitneyu(good_kappa, bad_kappa, alternative='two-sided')
        log(f"  Good κ: 중앙={np.median(good_kappa):.4f}, 평균={np.mean(good_kappa):.4f}")
        log(f"  Bad  κ: 중앙={np.median(bad_kappa):.4f}, 평균={np.mean(bad_kappa):.4f}")
        log(f"  Mann-Whitney U={u_stat:.1f}, p={p_mw:.4e}")
        if p_mw < 0.05:
            log(f"  ✅ 양성: Good/Bad κ 분포 유의하게 다름 (p={p_mw:.4e} < 0.05)")
        elif p_mw < 0.2:
            log(f"  ⚠️ 중립: 경계 (0.05 < p={p_mw:.4e} < 0.2)")
        else:
            log(f"  ❌ 음성: Good/Bad κ 차이 없음 (p={p_mw:.4e} > 0.2)")
        analysis_a = (u_stat, p_mw, len(good_kappa), len(bad_kappa))
    else:
        log(f"  ⚠️ 검정력 부족: 유효 bad={len(bad_kappa)}개 (최소 3개 필요)")
        log("  → 수학적 사실: n=0..200 범위 bad Gram 점은 n=126, 134, 195에만 존재")
        log("  → Mann-Whitney 검정 불가 (Type B: 전제 오류 — 예상보다 훨씬 적은 bad 점)")

    # ──────────────────────────────────────────────────────────────────────
    # (b) Good vs Bad d_zero 비교
    log("\n--- (b) Good vs Bad d_zero 비교 ---")
    log(f"  Good d_zero: n={len(good_dzero)}")
    if len(good_dzero) > 0:
        log(f"  Good d_zero: 중앙={np.median(good_dzero):.4f}, 평균={np.mean(good_dzero):.4f}")
    if len(bad_dzero) > 0:
        log(f"  Bad  d_zero: 중앙={np.median(bad_dzero):.4f}, 평균={np.mean(bad_dzero):.4f}")

    if len(bad_dzero) >= 2 and len(good_dzero) >= 3:
        u_stat_d, p_mw_d = stats.mannwhitneyu(good_dzero, bad_dzero, alternative='two-sided')
        log(f"  Mann-Whitney U={u_stat_d:.1f}, p={p_mw_d:.4e}")
        if p_mw_d < 0.05:
            log(f"  ✅ Bad Gram 점 영점 근접도 유의하게 다름 (p={p_mw_d:.4e})")
        else:
            log(f"  ❌ Good/Bad d_zero 차이 없음 (p={p_mw_d:.4e})")
        analysis_b = (u_stat_d, p_mw_d, float(np.median(good_dzero)), float(np.median(bad_dzero)))
    else:
        log(f"  ⚠️ 검정력 부족: bad={len(bad_dzero)}개")

    # ──────────────────────────────────────────────────────────────────────
    # (c) κ vs d_zero Spearman ρ (전체 Gram 점, cap 제외)
    log("\n--- (c) κ vs d_zero Spearman ρ (전체 Gram 점, cap 제외) ---")
    valid_mask = np.isfinite(all_kappa) & np.isfinite(all_dzero)
    all_k_v = all_kappa[valid_mask]
    all_d_v = all_dzero[valid_mask]

    if len(all_k_v) >= 5:
        rho_c, p_c = stats.spearmanr(all_d_v, all_k_v)
        log(f"  n={len(all_k_v)}")
        log(f"  Spearman ρ(d_zero, κ) = {rho_c:.4f}, p = {p_c:.4e}")
        if abs(rho_c) > 0.5:
            log(f"  ✅ 보조 양성: |ρ|={abs(rho_c):.4f} > 0.5")
            log(f"  → 곡률이 영점 근접도를 포착 (d_zero 클수록 κ {'증가' if rho_c > 0 else '감소'})")
        elif abs(rho_c) > 0.3:
            log(f"  ⚠️ 약한 상관: |ρ|={abs(rho_c):.4f} (0.3~0.5)")
        else:
            log(f"  ❌ 약한 상관: |ρ|={abs(rho_c):.4f} ≤ 0.3")
        analysis_c = (rho_c, p_c, len(all_k_v))
    else:
        log(f"  ⚠️ 표본 부족 ({len(all_k_v)}개)")

    # ──────────────────────────────────────────────────────────────────────
    # (d) 정성 비교: Bad 점 (최대 5개) vs Good 점 (d_zero 가까운 5개)
    log("\n--- (d) 정성 비교: Bad 5개 vs Good 5개 (d_zero 가까운 순) ---")

    bad_sorted_all  = sorted([r for r in all_results if r[2] == "bad"], key=lambda x: x[4])
    good_sorted_all = sorted([r for r in valid_results if r[2] == "good"], key=lambda x: x[4])

    if bad_sorted_all:
        bad_sample = bad_sorted_all[:min(5, len(bad_sorted_all))]
        log(f"\n  Bad Gram 점 ({len(bad_sorted_all)}개, d_zero 가까운 순):")
        log(f"  {'n':>4} {'t_gram':>10} {'Z(g_n)':>10} {'κ':>12} {'d_zero':>10} {'cap?':>5}")
        for (n, t, lab, k, d, z) in bad_sample:
            cap_flag = "✗" if k >= 1e15 else "OK"
            log(f"  {n:>4} {t:>10.4f} {z:>10.5f} {k:>12.4f} {d:>10.6f} {cap_flag:>5}")
    else:
        log("  ⚠️ Bad Gram 점 없음 — 수학적 사실 확인")

    if good_sorted_all:
        good_sample = good_sorted_all[:min(5, len(good_sorted_all))]
        log(f"\n  Good Gram 점 (유효, d_zero 가까운 순):")
        log(f"  {'n':>4} {'t_gram':>10} {'Z(g_n)':>10} {'κ':>12} {'d_zero':>10}")
        for (n, t, lab, k, d, z) in good_sample:
            log(f"  {n:>4} {t:>10.4f} {z:>10.5f} {k:>12.4f} {d:>10.6f}")

    return analysis_a, analysis_b, analysis_c


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 최종 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def final_verdict(analysis_a, analysis_b, analysis_c, n_good, n_bad, n_total):
    log("=" * 70)
    log("최종 판정")
    log("=" * 70)

    criteria_1 = False
    criteria_2 = False
    criteria_1_blocked = False

    log(f"\n  전체 Gram 점: {n_total}개 (Good={n_good}, Bad={n_bad})")

    if analysis_a is not None:
        u_stat, p_mw, ng, nb = analysis_a
        log(f"\n  기준 1 (Good vs Bad κ Mann-Whitney, p < 0.05):")
        if p_mw < 0.05:
            log(f"    ✅ PASS: p={p_mw:.4e} < 0.05 → κ가 Gram 구조를 감지")
            criteria_1 = True
        elif p_mw < 0.2:
            log(f"    ⚠️ 경계: p={p_mw:.4e} (0.05 ~ 0.2) → 중립")
        else:
            log(f"    ❌ FAIL: p={p_mw:.4e} > 0.2 → 곡률은 Gram 구조에 둔감")
    else:
        log(f"\n  기준 1: ⚠️ 검정 불가 — bad Gram 점 {n_bad}개 (최소 3개 유효 필요)")
        log(f"    [수학적 사실] n=0..200 범위에서 bad Gram 점: n=126, 134, 195만 존재")
        log(f"    [판단] Type B 전제 오류: 이 범위에서 통계 검정력 확보 불가")
        criteria_1_blocked = True

    if analysis_c is not None:
        rho_c, p_c, n_c = analysis_c
        log(f"\n  기준 2 (|ρ(κ, d_zero)| > 0.5):")
        if abs(rho_c) > 0.5:
            log(f"    ✅ PASS: |ρ|={abs(rho_c):.4f} > 0.5 → 곡률이 영점 근접도 포착")
            criteria_2 = True
        else:
            log(f"    ❌ FAIL: |ρ|={abs(rho_c):.4f} ≤ 0.5")

    log("\n" + "=" * 70)
    if criteria_1 and criteria_2:
        verdict = "✅✅ 양성 (두 기준 모두 충족) → 결과 #17 등록 (양성)"
    elif criteria_1:
        verdict = "✅ 양성 (기준 1만 충족) → 결과 #17 등록 (조건부 양성)"
    elif criteria_2:
        verdict = "✅ 보조 양성 (기준 2만 충족) → 결과 #17 등록 (조건부 양성)"
    elif criteria_1_blocked:
        if analysis_c is not None:
            rho_c = analysis_c[0]
            if abs(rho_c) > 0.3:
                verdict = "⚠️ 부분 보조 양성 — κ vs d_zero 약한 상관 확인. Mann-Whitney 검정 불가 (bad 점 부족). 설계 수정 필요"
            else:
                verdict = "⚠️ 결론 유보 — bad Gram 점 부족 (Type B 전제 오류). 확대 실험 필요 (n=0..500+)"
        else:
            verdict = "⚠️ 결론 유보 — 표본 부족 또는 cap 문제"
    else:
        verdict = "❌ 음성 — 곡률은 Gram 구조에 둔감 (논문 remark 소재)"

    log(f"  {verdict}")
    log("=" * 70)

    return verdict


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    outfile = open(OUT_PATH, 'w', encoding='utf-8')
    tee = Tee(sys.stdout, outfile)
    sys.stdout = tee

    log("=" * 70)
    log("Gram 점 곡률 분석 — Gram's law 위반의 기하학적 서명")
    log("사이클 #24 | 2026-04-14 (수정판)")
    log("=" * 70)
    log(f"mpmath dps = {mpmath.mp.dps} (t>100 구간 대응, 수정: 50→80)")
    log(f"δ = 0.03 (Lehmer 결과와 동일 기준)")
    log(f"Gram 점 범위: n=0..200")
    log(f"[수학적 사실] bad Gram 점: n=126(첫번째), 134, 195, 211, 232, 259, 262, 295, 302...")
    log()

    # 1단계: Gram 점 계산
    gram_points = compute_gram_points(n_max=200)
    if not gram_points:
        log("FATAL: Gram 점 계산 실패 — 중단")
        outfile.close()
        sys.stdout = sys.__stdout__
        sys.exit(1)

    # 2단계: Gram's law 분류
    classified = classify_gram_points(gram_points)
    if not classified:
        log("FATAL: 분류 실패 — 중단")
        outfile.close()
        sys.stdout = sys.__stdout__
        sys.exit(1)

    # 영점 수집 (d_zero 계산용)
    t_gram_max = max(t for (n, t) in gram_points) + 10
    t_zeros_max = max(t_gram_max, 450.0)
    log(f"\n영점 수집 (t∈[0, {t_zeros_max:.0f}]) ...")
    zeros_list = collect_zeros_for_range(t_max=t_zeros_max)
    log(f"  영점 {len(zeros_list)}개 수집")
    if len(zeros_list) == 0:
        log("  ⚠️ 영점 0개 — 탐색 로직 점검 필요")

    # 3단계: κ 측정 + d_zero 계산
    results = measure_kappa_and_dzero(classified, zeros_list, delta=0.03)

    if len(results) == 0:
        log("FATAL: κ 측정 결과 없음 — 중단")
        outfile.close()
        sys.stdout = sys.__stdout__
        sys.exit(1)

    # 4단계: 상관 분석
    analysis_a, analysis_b, analysis_c = analyze_correlations(results)

    # 최종 판정
    n_good = sum(1 for r in results if r[2] == "good")
    n_bad  = sum(1 for r in results if r[2] == "bad")
    verdict = final_verdict(analysis_a, analysis_b, analysis_c, n_good, n_bad, len(results))

    # 전체 결과 테이블
    log("\n" + "=" * 70)
    log("전체 Gram 점 결과 테이블 (처음 30개)")
    log("=" * 70)
    log(f"{'n':>4} {'t_gram':>10} {'분류':>6} {'κ':>14} {'d_zero':>10} {'Z(g_n)':>12}")
    for (n, t, label, kappa, d_zero, z_val) in results[:30]:
        cap_marker = "*" if kappa >= 1e15 else " "
        log(f"{n:>4} {t:>10.4f} {label:>6} {kappa:>13.4f}{cap_marker} {d_zero:>10.6f} {z_val:>12.6f}")

    elapsed = time.time() - t_start
    log(f"\n총 소요 시간: {elapsed:.1f}초 ({elapsed/60:.1f}분)")
    log(f"결과 파일: {os.path.abspath(OUT_PATH)}")

    outfile.close()
    sys.stdout = sys.__stdout__
    print(f"\n[완료] 총 소요: {elapsed:.1f}초", flush=True)
    print(f"판정: {verdict}", flush=True)


if __name__ == "__main__":
    main()
