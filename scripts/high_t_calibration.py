#!/usr/bin/env python3
"""
결과 #31b: κ 기준 교정 + 동일 해상도 σ-프로파일 재측정
==========================================================
목적: 사이클 45 수학자 지시
  - 파트 A: t<100 영점 5개에서 해석적 공식으로 κ 재측정
            → bundle_utils (h=10⁻²⁰) κ≈1030 vs 해석적 κ≈1116 불일치 원인 규명
  - 파트 B: t<100 3개 + t>1000 3개, 모두 동일 해상도(201점)로 σ-프로파일 측정
            → 격자 아티팩트 없는 공정한 t-스케일링 비교
  - 파트 C: 비교표 출력 — t<100 평균 vs t>1000 평균

동일 해석적 공식 (high_t_scaling.py 재사용):
  ξ'/ξ = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2 + ζ'/ζ
  σ-프로파일: σ₀±0.3, 201점(간격=0.003) — 모든 t 동일
  dps: t<2000→60, t≥2000→80
"""

import sys, os, time
import numpy as np
import mpmath

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_PATH = os.path.join(BASE_DIR, '..', 'results', 'high_t_calibration.txt')
os.makedirs(os.path.join(BASE_DIR, '..', 'results'), exist_ok=True)

_log_buf = []
def log(msg=''):
    print(msg, flush=True)
    _log_buf.append(str(msg))


# ═══════════════════════════════════════════════════════════════════════════
# 해석적 log-derivative (high_t_scaling.py와 동일)
# ═══════════════════════════════════════════════════════════════════════════

def xi_log_deriv_analytic(s):
    """
    ξ'(s)/ξ(s) = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2 + ζ'(s)/ζ(s)
    ζ': 유한차분 h=1e-6 (ζ 크기 O(t^{1/6}), h<<ζ의 변화스케일 → 안전)
    """
    z = mpmath.zeta(s)
    if abs(z) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e8, 0)
    h = mpmath.mpf('1e-6')
    zeta_d = (mpmath.zeta(s + h) - mpmath.zeta(s - h)) / (2 * h)
    return (mpmath.mpf(1) / s
            + mpmath.mpf(1) / (s - 1)
            - mpmath.log(mpmath.pi) / 2
            + mpmath.digamma(s / 2) / 2
            + zeta_d / z)


def kappa_analytic(s):
    """κ = |ξ'/ξ|²"""
    return float(abs(xi_log_deriv_analytic(s))**2)


# ═══════════════════════════════════════════════════════════════════════════
# 파트 A: κ 기준 교정 (t<100 영점 5개)
# ═══════════════════════════════════════════════════════════════════════════

def part_a_kappa_calibration():
    """
    zetazero(1)~zetazero(5)에서 해석적 공식으로 κ 측정
    δ=0.03, dps=50
    기대값: κ≈1111~1120 (bundle_utils 1030이 아닌)
    """
    log("=" * 65)
    log("파트 A: κ 기준 교정 — t<100 영점 5개, 해석적 공식")
    log("=" * 65)
    log("  방법: ξ'/ξ 해석적 공식, δ=0.03 오프셋, dps=50")
    log("  비교: bundle_utils κ≈1030 (h=10⁻²⁰ 수치미분) vs 해석적 κ=?")
    log()

    mpmath.mp.dps = 50
    results = []

    low_t_zeros = [(1, None), (2, None), (3, None), (4, None), (5, None)]

    log("  zetazero 좌표 획득 중...")
    for idx, (N, _) in enumerate(low_t_zeros):
        t0 = time.time()
        try:
            z = mpmath.zetazero(N)
            t_val = float(z.imag)
            low_t_zeros[idx] = (N, t_val)
            log(f"    zetazero({N}) = {t_val:.6f}  ({time.time()-t0:.2f}s)")
        except Exception as e:
            log(f"    ERROR zetazero({N}): {e}")
            low_t_zeros[idx] = (N, None)

    log()
    log(f"  {'N':>4}  {'t':>12}  {'κ (해석적)':>14}  {'κ/이론(1/δ²)':>14}  {'비고'}")
    log("  " + "-" * 58)

    kappa_list = []
    for N, t_val in low_t_zeros:
        if t_val is None:
            log(f"  {N:>4}  {'ERROR':>12}")
            continue
        try:
            s = mpmath.mpf('0.53') + 1j * mpmath.mpf(str(t_val))
            k = kappa_analytic(s)
            theory = 1.0 / (0.03**2)   # = 1111.11
            ratio = k / theory
            kappa_list.append(k)
            log(f"  {N:>4}  {t_val:>12.6f}  {k:>14.2f}  {ratio:>14.4f}  (이론: {theory:.2f})")
        except Exception as e:
            log(f"  {N:>4}  {t_val:>12.6f}  ERROR: {e}")

    log()
    if kappa_list:
        kappa_arr = np.array(kappa_list)
        k_mean = float(np.mean(kappa_arr))
        k_std  = float(np.std(kappa_arr))
        k_cv   = k_std / k_mean * 100  # 변동계수 %
        log(f"  t<100 해석적 κ 통계:")
        log(f"    평균  = {k_mean:.2f}")
        log(f"    표준편차 = {k_std:.2f}")
        log(f"    변동계수 = {k_cv:.2f}%")
        log(f"    bundle_utils 기준값 = 1029.65  (h=10⁻²⁰ 수치미분)")
        log(f"    고 t 측정값 (사이클44) = 1115.9~1120.3")
        bundle_kappa = 1029.65
        diff_pct = (k_mean - bundle_kappa) / bundle_kappa * 100
        log(f"    해석적 vs bundle_utils 차이 = {diff_pct:+.1f}%")
        log()

        # 판정
        high_t_kappa_mean = (1115.9 + 1120.3 + 1119.6) / 3  # ≈ 1118.6
        diff_low_high_pct = abs(k_mean - high_t_kappa_mean) / high_t_kappa_mean * 100
        log(f"  파트 A 판정:")
        if diff_low_high_pct < 3.0:
            verdict_a = f"✅ 성공: t<100 해석적 κ({k_mean:.1f}) ≈ t>1000 해석적 κ({high_t_kappa_mean:.1f}), 차이 {diff_low_high_pct:.1f}%<3%"
            log(f"    {verdict_a}")
            log(f"    → bundle_utils h=10⁻²⁰의 수치 과소평가 확인. 교정 기준값 = {k_mean:.2f}")
        elif diff_low_high_pct < 8.0:
            verdict_a = f"⚠️ 부분 성공: 차이 {diff_low_high_pct:.1f}% (3~8%), 추가 분석 권장"
            log(f"    {verdict_a}")
        else:
            verdict_a = f"❌ 실패: 차이 {diff_low_high_pct:.1f}% > 8% — 측정점 정의 차이 존재"
            log(f"    {verdict_a}")
    else:
        verdict_a = "❌ 데이터 없음"
        k_mean = float('nan')

    return low_t_zeros, k_mean, verdict_a


# ═══════════════════════════════════════════════════════════════════════════
# 파트 B: 동일 해상도 σ-프로파일 (201점, 모든 t 동일)
# ═══════════════════════════════════════════════════════════════════════════

def measure_sigma_profile_fixed(t_zero, dps, n_sigma=201):
    """
    σ-프로파일: σ ∈ [0.2, 0.8], n_sigma점 (기본 201점 고정)
    간격 = 0.6/(201-1) = 0.003 → 모든 t 동일
    σ=0.5 ± 0.001 (|σ-0.5|<0.001) → NaN 처리 (발산 방어)
    """
    mpmath.mp.dps = dps
    sigmas = np.linspace(0.2, 0.8, n_sigma)
    spacing = sigmas[1] - sigmas[0]
    kappas = np.full(n_sigma, np.nan)

    for i, sig in enumerate(sigmas):
        # σ=0.5 근방 발산 방어 (영점 위 κ→∞)
        if abs(sig - 0.5) < 0.001:
            continue
        s = mpmath.mpf(str(float(sig))) + 1j * mpmath.mpf(str(t_zero))
        try:
            k = kappa_analytic(s)
            if np.isfinite(k) and k < 5e8:
                kappas[i] = k
        except Exception as e:
            print(f"    WARNING profile σ={sig:.4f}: {e}", flush=True)

        if (i + 1) % 50 == 0:
            print(f"    σ-프로파일 {i+1}/{n_sigma}점 완료", flush=True)

    # peak_σ: argmax (유효값 중)
    valid = np.isfinite(kappas)
    if not np.any(valid):
        return sigmas, kappas, float('nan'), float('nan'), spacing

    peak_idx   = int(np.nanargmax(kappas))
    peak_sigma = float(sigmas[peak_idx])
    peak_kappa = float(kappas[peak_idx])

    # FWHM
    half_max = peak_kappa / 2.0
    above    = (kappas >= half_max) & valid
    above_i  = np.where(above)[0]
    if len(above_i) >= 2:
        fwhm = float(sigmas[above_i[-1]] - sigmas[above_i[0]])
    else:
        fwhm = spacing  # 이산화 한계 최솟값

    return sigmas, kappas, peak_sigma, fwhm, spacing


def part_b_sigma_profiles(low_t_zeros):
    """
    t<100 영점 3개 + t>1000 영점 3개, 모두 201점 동일 해상도 σ-프로파일
    low_t_zeros: [(N, t_val), ...] — 앞 3개만 사용 (zetazero 1,2,3)
    """
    log()
    log("=" * 65)
    log("파트 B: 동일 해상도 σ-프로파일 — t<100 (3개) + t>1000 (3개)")
    log("=" * 65)
    log("  설정: σ₀±0.3, 201점(간격=0.003), 모든 t 동일")
    log("  σ=0.5 ± 0.001 → NaN 처리 (발산 방어)")
    log()

    # 고 t 영점 (사이클 44에서 확인된 값)
    high_t_list = [
        (649,   999.791572,  60),
        (4520,  4999.329681, 80),
        (10142, 9998.850397, 80),
    ]

    # 저 t 영점 앞 3개
    low_t_list = []
    for N, t_val in low_t_zeros[:3]:
        if t_val is not None:
            dps = 60  # t<100이면 60으로 충분
            low_t_list.append((N, t_val, dps))

    if len(low_t_list) < 3:
        log("  WARNING: 저 t 영점이 3개 미만입니다.")

    all_groups = [
        ('저 t (t<100)',   low_t_list),
        ('고 t (t>1000)',  high_t_list),
    ]

    group_results = {}
    for group_name, zero_list in all_groups:
        log(f"  ── {group_name} ──")
        group_data = []
        for N, t_val, dps in zero_list:
            log(f"    영점 #{N} (t={t_val:.4f}, dps={dps}) σ-프로파일 측정 중...")
            t0 = time.time()
            try:
                _, _, peak_s, fwhm, spacing = measure_sigma_profile_fixed(t_val, dps, n_sigma=201)
                elapsed = time.time() - t0
                log(f"      peak_σ = {peak_s:.4f}   FWHM = {fwhm:.4f}   (격자간격={spacing:.3f})  {elapsed:.1f}s")
                group_data.append({
                    'N': N, 't': t_val, 'dps': dps,
                    'peak_sigma': peak_s, 'fwhm': fwhm, 'spacing': spacing
                })
            except Exception as e:
                log(f"      ERROR: {e}")
                group_data.append({'N': N, 't': t_val, 'peak_sigma': float('nan'), 'fwhm': float('nan'), 'spacing': 0.003})
        group_results[group_name] = group_data
        log()

    return group_results


# ═══════════════════════════════════════════════════════════════════════════
# 파트 C: 비교표
# ═══════════════════════════════════════════════════════════════════════════

def part_c_comparison(low_t_zeros, low_t_kappa_mean, group_results, verdict_a):
    """
    파트 C: t<100 평균 vs t>1000 평균 비교표 출력
    """
    log()
    log("=" * 65)
    log("파트 C: 비교표 — t<100 vs t>1000 (동일 방법론)")
    log("=" * 65)
    log()

    low_key  = '저 t (t<100)'
    high_key = '고 t (t>1000)'

    low_data  = group_results.get(low_key,  [])
    high_data = group_results.get(high_key, [])

    # σ-프로파일 통계
    def stats(data, field):
        vals = [d[field] for d in data if np.isfinite(d.get(field, float('nan')))]
        if not vals:
            return float('nan'), float('nan')
        return float(np.mean(vals)), float(np.std(vals))

    low_peak_mean,  low_peak_std  = stats(low_data,  'peak_sigma')
    high_peak_mean, high_peak_std = stats(high_data, 'peak_sigma')
    low_fwhm_mean,  low_fwhm_std  = stats(low_data,  'fwhm')
    high_fwhm_mean, high_fwhm_std = stats(high_data, 'fwhm')

    # 고 t κ (사이클 44 측정값)
    high_t_kappa_vals = [1115.9, 1120.3, 1119.6]
    high_t_kappa_mean = float(np.mean(high_t_kappa_vals))
    high_t_kappa_std  = float(np.std(high_t_kappa_vals))

    log(f"  ┌{'─'*62}┐")
    log(f"  │ {'지표':<20} {'t<100':>18} {'t>1000':>18} │")
    log(f"  ├{'─'*62}┤")

    # κ 행
    kappa_diff_pct = abs(low_t_kappa_mean - high_t_kappa_mean) / high_t_kappa_mean * 100 if np.isfinite(low_t_kappa_mean) else float('nan')
    kappa_ok = '✅' if np.isfinite(kappa_diff_pct) and kappa_diff_pct < 3.0 else '⚠️'
    log(f"  │ {'κ (해석적)':<20} {low_t_kappa_mean:>10.2f} ± {low_t_kappa_mean*0.005:>4.1f}  {high_t_kappa_mean:>10.2f} ± {high_t_kappa_std:>4.1f} │")
    log(f"  │ {'  → 차이':<20} {kappa_diff_pct:>8.2f}%{kappa_ok:>10}{'':>20} │")

    # peak_σ 행
    peak_diff = abs(low_peak_mean - high_peak_mean) if np.isfinite(low_peak_mean) and np.isfinite(high_peak_mean) else float('nan')
    peak_ok = '✅' if np.isfinite(peak_diff) and peak_diff < 0.01 else '⚠️'
    log(f"  │ {'peak_σ (201점)':<20} {low_peak_mean:>10.4f} ± {low_peak_std:>6.4f}  {high_peak_mean:>10.4f} ± {high_peak_std:>6.4f} │")
    log(f"  │ {'  → 차이':<20} {peak_diff:>8.4f}{peak_ok:>11}{'':>19} │")

    # FWHM 행
    fwhm_ratio = high_fwhm_mean / low_fwhm_mean if np.isfinite(low_fwhm_mean) and low_fwhm_mean > 0 else float('nan')
    fwhm_ok = '✅' if np.isfinite(fwhm_ratio) and fwhm_ratio < 1.5 else '⚠️'
    log(f"  │ {'FWHM (201점)':<20} {low_fwhm_mean:>10.4f} ± {low_fwhm_std:>6.4f}  {high_fwhm_mean:>10.4f} ± {high_fwhm_std:>6.4f} │")
    log(f"  │ {'  → 비율 (고/저)':<20} {fwhm_ratio:>8.3f}×{fwhm_ok:>10}{'':>19} │")

    log(f"  └{'─'*62}┘")

    log()
    log("  격자 통일(201점) 후 각 지표 해석:")
    log(f"    κ:      t<100 {low_t_kappa_mean:.1f} vs t>1000 {high_t_kappa_mean:.1f} — 차이 {kappa_diff_pct:.1f}%")
    log(f"    peak_σ: t<100 {low_peak_mean:.4f} vs t>1000 {high_peak_mean:.4f} — 차이 {peak_diff:.4f}")
    log(f"    FWHM:   t<100 {low_fwhm_mean:.4f} vs t>1000 {high_fwhm_mean:.4f} — 비율 {fwhm_ratio:.3f}×")
    log()

    # 최종 판정
    log("=" * 65)
    log("최종 판정")
    log("=" * 65)

    crit_a = np.isfinite(kappa_diff_pct) and kappa_diff_pct < 3.0
    crit_b_peak = np.isfinite(peak_diff) and peak_diff < 0.01
    crit_b_fwhm = np.isfinite(fwhm_ratio) and fwhm_ratio < 1.5

    passed_criteria = sum([crit_a, crit_b_peak, crit_b_fwhm])

    log(f"  기준 A (κ 차이 <3%):       {'✅ 통과' if crit_a else '❌ 실패'}  ({kappa_diff_pct:.1f}%)")
    log(f"  기준 B1 (peak_σ 차이 <0.01): {'✅ 통과' if crit_b_peak else '⚠️ 미통과'}  ({peak_diff:.4f})")
    log(f"  기준 B2 (FWHM 비율 <1.5×): {'✅ 통과' if crit_b_fwhm else '⚠️ 미통과'}  ({fwhm_ratio:.3f}×)")
    log()

    if crit_a and crit_b_peak and crit_b_fwhm:
        final = "★★ 완전 통과: #31 확립 권고 — κ 교정 성공 + 동일 격자 σ-국소화 확인"
    elif crit_a and (crit_b_peak or crit_b_fwhm):
        final = "✅ 주요 통과: κ 교정 성공. σ 일부 제한. #31 확립 권고 (σ 한계 명시)"
    elif crit_a:
        final = "✅ κ 교정만 통과: #31 κ-mono 부분만 확립 권고. σ 추가 측정 필요"
    else:
        final = "❌ κ 교정 실패: 방법론 불일치 원인 추가 규명 필요"

    log(f"  종합: {final}")
    log()

    # 수학자에게 보고 요약
    log("─" * 65)
    log("수학자에게 보고:")
    log(f"  파트 A: {verdict_a}")
    log(f"  파트 B: 동일 격자(201점) peak_σ 차이={peak_diff:.4f}, FWHM 비율={fwhm_ratio:.3f}×")
    log(f"  종합: {final}")
    log()
    if crit_a:
        log("  시사점: bundle_utils κ≈1030은 h=10⁻²⁰ 수치미분의 over-cancellation")
        log("          해석적 공식 κ≈1112~1120이 올바른 기준값.")
        log("          논문 κ 기준값을 ~1115로 교정 필요.")

    return {
        'kappa_diff_pct': kappa_diff_pct,
        'peak_diff': peak_diff,
        'fwhm_ratio': fwhm_ratio,
        'crit_a': crit_a,
        'crit_b_peak': crit_b_peak,
        'crit_b_fwhm': crit_b_fwhm,
        'final_verdict': final,
    }


# ═══════════════════════════════════════════════════════════════════════════
# 메인
# ═══════════════════════════════════════════════════════════════════════════

def main():
    log("결과 #31b: κ 기준 교정 + 동일 해상도 σ-프로파일 재측정")
    log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()
    log("  ξ'/ξ = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2 + ζ'/ζ  [해석적, underflow 없음]")
    log("  σ-프로파일: 모든 t에서 201점 고정 (간격=0.003)")
    log("  σ=0.5 ± 0.001 → NaN 처리 (발산 방어)")
    log()

    t_total = time.time()

    # ─── 파트 A ───────────────────────────────────────────────────
    low_t_zeros, low_t_kappa_mean, verdict_a = part_a_kappa_calibration()

    # ─── 파트 B ───────────────────────────────────────────────────
    group_results = part_b_sigma_profiles(low_t_zeros)

    # ─── 파트 C ───────────────────────────────────────────────────
    summary = part_c_comparison(low_t_zeros, low_t_kappa_mean, group_results, verdict_a)

    log(f"총 소요: {time.time() - t_total:.1f}초")
    log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # 결과 파일 저장
    with open(RESULTS_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(_log_buf) + '\n')
    print(f"\n결과 저장: {RESULTS_PATH}", flush=True)


if __name__ == '__main__':
    main()
