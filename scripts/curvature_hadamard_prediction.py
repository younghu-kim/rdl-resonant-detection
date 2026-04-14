#!/usr/bin/env python3
"""
결과 #33: Hadamard 전개 검증 — Δκ 이론 예측 대 실측 비교
=============================================================
목적: 사이클 47 수학자 지시
  - ξ'/ξ(s) = 1/(s-ρ_n) + G(s) + P_n(s) 분해
  - G(n): 감마+극 기여
  - S_all(n): 제타 영점 기여 (= Re[ζ'/ζ(s) - 1/δ])
    └ S_local(n): ±K 인접 영점 직접 합 (K=5)
    └ S_tail(n): 잔차 = S_all - S_local - S_reflect
    └ S_reflect_n: 반사 영점 ρ̄_n 기여
  - Δκ_decomposed = 2*(G+S_all)/δ + |R|² vs Δκ_measured (#32)
  - 기울기 분해: a_G + a_S_all ≈ 1.59?
  - N=7000 이상치 해부

이론:
  ξ'/ξ(s) = G(s) + ζ'/ζ(s), G(s) = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2
  ζ'/ζ(s) = 1/(s-ρ_n) + P_n(s)  (ρ_n = 가장 가까운 비자명 영점)
  R(s) = ξ'/ξ(s) - 1/δ = G(s) + P_n(s)
  Δκ = 2·Re[R]/δ + |R|²  (δ-독립성: Re[R] = O(δ))
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats as scipy_stats

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_PATH = os.path.join(BASE_DIR, '..', 'results', 'curvature_hadamard_prediction.txt')
os.makedirs(os.path.join(BASE_DIR, '..', 'results'), exist_ok=True)

_log_buf = []
def log(msg=''):
    print(msg, flush=True)
    _log_buf.append(str(msg))

def save():
    with open(RESULTS_PATH, 'w') as f:
        f.write('\n'.join(_log_buf))

# ═══════════════════════════════════════════════════════════════════════════
# 상수
# ═══════════════════════════════════════════════════════════════════════════

ZERO_INDICES = [1, 2, 3, 5, 10, 20, 50, 100, 200, 500,
                649, 1000, 2000, 4520, 5000, 7000, 10142, 12000, 15000, 20000]

DELTA = 0.03   # δ 오프셋

# #32 결과에서 가져온 Δκ_measured (직접 복사)
DELTA_KAPPA_MEASURED = {
    1: 0.4713, 2: 0.8478, 3: 0.6732, 5: 0.8069,
    10: 1.2841, 20: 1.9791, 50: 1.6333, 100: 4.5597,
    200: 5.0646, 500: 4.2606, 649: 4.8267, 1000: 7.8920,
    2000: 6.0463, 4520: 9.1412, 5000: 7.4309, 7000: 14.0846,
    10142: 8.4951, 12000: 8.1503, 15000: 12.0218, 20000: 11.9533,
}


def get_dps(t):
    if t > 30000: return 100
    if t > 5000: return 80
    return max(50, int(30 + t / 200))


# ═══════════════════════════════════════════════════════════════════════════
# 해석적 계산 함수
# ═══════════════════════════════════════════════════════════════════════════

def G_component(s):
    """
    G(s) = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2
    ξ'/ξ(s)의 감마+극+비영점 기여 (복소수 반환)
    """
    return (mpmath.mpf(1)/s
            + mpmath.mpf(1)/(s - 1)
            - mpmath.log(mpmath.pi)/2
            + mpmath.digamma(s/2)/2)


def zeta_log_deriv(s):
    """
    ζ'(s)/ζ(s) 수치 계산 (h=1e-6 유한 차분)
    """
    z = mpmath.zeta(s)
    if abs(z) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e8, 0)
    h = mpmath.mpf('1e-6')
    zd = (mpmath.zeta(s + h) - mpmath.zeta(s - h)) / (2 * h)
    return zd / z


def xi_log_deriv_analytic(s):
    """ξ'/ξ(s) = G(s) + ζ'/ζ(s)"""
    return G_component(s) + zeta_log_deriv(s)


# ═══════════════════════════════════════════════════════════════════════════
# 파트 A: 검증 — Δκ_direct vs Δκ_measured
# ═══════════════════════════════════════════════════════════════════════════

def part_a_verification():
    """
    Δκ_direct = |ξ'/ξ(s)|² - 1/δ² 재계산 → #32와 일치 확인
    """
    log("=" * 70)
    log("파트 A: Δκ_direct 재검증 — #32와 일치 확인")
    log("=" * 70)
    log(f"  목적: 해석적 공식으로 재계산 → 결과 #32 Δκ_measured와 일치 여부")
    log(f"  δ = {DELTA},  1/δ² = {1/DELTA**2:.4f}")
    log()
    log(f"  {'N':>6}  {'t':>12}  {'Δκ_direct':>12}  {'Δκ_measured':>12}  {'오차%':>8}  판정")
    log("  " + "-" * 72)

    results = []
    t_start = time.time()

    for N in ZERO_INDICES:
        dps = 50  # 먼저 t를 얻기 위한 임시값
        mpmath.mp.dps = 60
        try:
            z = mpmath.zetazero(N)
        except Exception as e:
            log(f"  {N:>6}  ERROR: {e}")
            continue
        t_n = float(z.imag)
        dps = get_dps(t_n)
        mpmath.mp.dps = dps

        s = mpmath.mpf('0.5') + mpmath.mpf(str(DELTA)) + 1j * mpmath.mpf(str(t_n))

        try:
            xi_ld = xi_log_deriv_analytic(s)
            kappa_direct = float(abs(xi_ld)**2)
            dk_direct = kappa_direct - 1/DELTA**2
            dk_meas = DELTA_KAPPA_MEASURED[N]
            err_pct = abs(dk_direct - dk_meas) / max(abs(dk_meas), 0.01) * 100

            ok = "✅" if err_pct < 1.0 else ("⚠️" if err_pct < 5.0 else "❌")
            log(f"  {N:>6}  {t_n:>12.4f}  {dk_direct:>12.4f}  {dk_meas:>12.4f}  {err_pct:>8.3f}%  {ok}")
            results.append({
                'N': N, 't': t_n, 'dk_direct': dk_direct,
                'dk_meas': dk_meas, 'err_pct': err_pct,
                'xi_ld': xi_ld, 's': s, 'dps': dps
            })
        except Exception as e:
            log(f"  {N:>6}  {t_n:>12.4f}  ERROR: {e}")

    elapsed = time.time() - t_start
    n_ok = sum(1 for r in results if r['err_pct'] < 1.0)
    log()
    log(f"  검증 통과: {n_ok}/{len(results)}  (오차 < 1.0%)")
    log(f"  파트 A 소요: {elapsed:.1f}초")

    return results


# ═══════════════════════════════════════════════════════════════════════════
# 파트 B: 기여 분해 (G + S_all → S_local + S_tail + S_reflect)
# ═══════════════════════════════════════════════════════════════════════════

def fetch_neighbor_zeros(n, K, t_n):
    """
    zetazero(n-K : n+K) 범위 내 인접 영점 허수부 반환 (n 자신 제외)
    """
    lo = max(1, n - K)
    hi = n + K
    result = []
    for m in range(lo, hi + 1):
        if m == n:
            continue
        try:
            zz = mpmath.zetazero(m)
            result.append(float(zz.imag))
        except Exception as e:
            print(f"    WARNING zetazero({m}): {e}", flush=True)
    return result


def decompose_contributions(partA_results):
    """
    각 영점에서 G_n, S_all_n, S_local_n, S_reflect_n, S_tail_n 계산
    """
    log()
    log("=" * 70)
    log("파트 B: Hadamard 기여 분해")
    log("=" * 70)
    log("  G(n)     = Re[1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2]")
    log("  S_all(n) = Re[ζ'/ζ(s)] - 1/δ  (영점 기여 실부, 극 제거)")
    log("  S_local  = Σ_{|m-n|≤K} δ/(δ²+(t_n-t_m)²),  K=5")
    log("  S_reflect_n = δ/(δ²+(2t_n)²)  [ρ̄_n 기여]")
    log("  S_tail   = S_all - S_local - S_reflect_n")
    log("  Re[R]    = G + S_all  (= 전체 보정 실부)")
    log()
    log("  검증: Δκ_formula = 2·Re[R]/δ + |R|²  ≈  Δκ_measured")
    log()

    # 가장 고-t 인덱스에서는 K를 줄임
    def get_K(n):
        if n >= 10000: return 3
        if n >= 5000: return 5
        if n >= 1000: return 5
        return 10

    header = (f"  {'N':>6}  {'t':>10}  {'G_n':>9}  {'S_all':>9}  "
              f"{'S_loc':>8}  {'S_refl':>8}  {'S_tail':>9}  "
              f"{'Re[R]':>9}  {'Im[R]':>9}  "
              f"{'Δκ_form':>10}  {'Δκ_meas':>10}  {'err%':>7}")
    log(header)
    log("  " + "-" * 115)

    decomp = []
    t_start = time.time()

    for r in partA_results:
        N = r['N']
        t_n = r['t']
        s = r['s']
        dps = r['dps']
        mpmath.mp.dps = dps
        K = get_K(N)

        try:
            # G 기여
            Gval = G_component(s)
            G_n = float(mpmath.re(Gval))
            Im_G_n = float(mpmath.im(Gval))

            # ζ'/ζ 기여 → S_all
            zld = zeta_log_deriv(s)
            Re_zld = float(mpmath.re(zld))
            Im_zld = float(mpmath.im(zld))
            S_all_n = Re_zld - 1.0/DELTA      # Re[ζ'/ζ] - 1/δ (극 제거)
            Im_S_all_n = Im_zld               # Im[ζ'/ζ]

            # S_local: 인접 K개 영점
            t_start_local = time.time()
            neighbors = fetch_neighbor_zeros(N, K, t_n)
            S_local_n = sum(
                DELTA / (DELTA**2 + (t_n - tm)**2)
                for tm in neighbors
            )
            # S_reflect_local: 인접 영점의 반사 기여 (작음)
            S_reflect_local_n = sum(
                DELTA / (DELTA**2 + (t_n + tm)**2)
                for tm in neighbors
            )

            # S_reflect_n: ρ_n 자신의 반사 영점 기여
            S_reflect_n = DELTA / (DELTA**2 + (2*t_n)**2)

            # S_tail: 잔차
            S_tail_n = S_all_n - S_local_n - S_reflect_n - S_reflect_local_n

            # R(s) = G(s) + P_n(s), Re[R] = G_n + S_all_n
            Re_R = G_n + S_all_n
            Im_R = Im_G_n + Im_S_all_n

            # Δκ 공식
            dk_formula = 2*Re_R/DELTA + Re_R**2 + Im_R**2
            dk_meas = r['dk_meas']
            err_pct = abs(dk_formula - dk_meas) / max(abs(dk_meas), 0.01) * 100

            log(f"  {N:>6}  {t_n:>10.2f}"
                f"  {G_n:>9.4f}  {S_all_n:>9.4f}"
                f"  {S_local_n:>8.4f}  {S_reflect_n:>8.6f}  {S_tail_n:>9.4f}"
                f"  {Re_R:>9.6f}  {Im_R:>9.4f}"
                f"  {dk_formula:>10.4f}  {dk_meas:>10.4f}  {err_pct:>7.2f}%")

            decomp.append({
                'N': N, 't': t_n,
                'G_n': G_n, 'Im_G_n': Im_G_n,
                'S_all_n': S_all_n, 'Im_S_all_n': Im_S_all_n,
                'S_local_n': S_local_n, 'S_reflect_n': S_reflect_n,
                'S_reflect_local_n': S_reflect_local_n,
                'S_tail_n': S_tail_n,
                'Re_R': Re_R, 'Im_R': Im_R,
                'dk_formula': dk_formula, 'dk_meas': dk_meas,
                'err_pct': err_pct,
            })

        except Exception as e:
            log(f"  {N:>6}  {t_n:>10.2f}  ERROR: {e}")
            import traceback
            traceback.print_exc()

    elapsed = time.time() - t_start
    log()
    log(f"  파트 B 소요: {elapsed:.1f}초")

    # 정확도 요약
    errs = [d['err_pct'] for d in decomp]
    mean_err = np.mean(errs)
    median_err = np.median(errs)
    n_10 = sum(1 for e in errs if e < 10.0)
    n_20 = sum(1 for e in errs if e < 20.0)
    log()
    log(f"  Δκ_formula vs Δκ_measured 정확도:")
    log(f"    평균 오차: {mean_err:.2f}%")
    log(f"    중앙 오차: {median_err:.2f}%")
    log(f"    < 10%: {n_10}/{len(decomp)}")
    log(f"    < 20%: {n_20}/{len(decomp)}")

    if mean_err < 10.0:
        log(f"  ✅ 강한 양성: 평균 오차 < 10%")
    elif mean_err < 20.0:
        log(f"  ✅ 양성: 평균 오차 < 20%")
    else:
        log(f"  ❌ 음성: 평균 오차 > 20% — 분해 불완전")

    return decomp


# ═══════════════════════════════════════════════════════════════════════════
# 파트 C: 기울기 분해 (vs log(t/2π))
# ═══════════════════════════════════════════════════════════════════════════

def part_c_slope_decomposition(decomp):
    """
    G_n, S_all_n, S_local_n, S_tail_n, Im_R² 각각 vs log(t/2π) 회귀
    이론: Δκ = 1.59·log(t/2π) + const
    어느 항이 이 기울기를 담당하는가?
    """
    log()
    log("=" * 70)
    log("파트 C: 기울기 분해 — 각 기여항의 log(t/2π) 회귀")
    log("=" * 70)
    log("  목표: Δκ = 1.59·log(t/2π) + const 에서 기울기 출처 식별")
    log("  이론 예측: a_G ≈ 0.5, a_S_all ≈ 0.5 → a_total·(2/δ) ≈ 1.59?")
    log()

    log_t = np.array([np.log(d['t'] / (2*np.pi)) for d in decomp])

    def regress(y_arr, label):
        if len(y_arr) < 3:
            return float('nan'), float('nan'), float('nan')
        slope, intercept, r, p, se = scipy_stats.linregress(log_t, y_arr)
        log(f"    {label:30s}: slope={slope:+.4f} ± {se:.4f},  R²={r**2:.4f},  p={p:.2e}")
        return slope, r**2, p

    log("  [Re[R] 기여 분해]")
    G_arr = np.array([d['G_n'] for d in decomp])
    S_all_arr = np.array([d['S_all_n'] for d in decomp])
    S_local_arr = np.array([d['S_local_n'] for d in decomp])
    S_tail_arr = np.array([d['S_tail_n'] for d in decomp])
    S_reflect_arr = np.array([d['S_reflect_n'] for d in decomp])
    Re_R_arr = np.array([d['Re_R'] for d in decomp])
    Im_R_arr = np.array([d['Im_R'] for d in decomp])
    Im_R_sq_arr = Im_R_arr**2

    a_G, r2_G, p_G = regress(G_arr, "G(n) [감마+극]")
    a_S_all, r2_S, p_S = regress(S_all_arr, "S_all(n) [전체 제타 영점]")
    a_S_local, r2_Sl, p_Sl = regress(S_local_arr, "S_local(n) [±K 인접 영점]")
    a_S_tail, r2_St, p_St = regress(S_tail_arr, "S_tail(n) [원거리+상수]")
    a_Re_R, r2_ReR, p_ReR = regress(Re_R_arr, "Re[R](n) = G + S_all")
    log()
    log("  [Im[R] 기여 (Δκ의 |R|² 항)]")
    a_Im_R, r2_ImR, p_ImR = regress(Im_R_arr, "Im[R](n)")
    a_Im_R_sq, r2_ImRsq, p_ImRsq = regress(Im_R_sq_arr, "Im[R]²(n)  [= |R|² 근사]")

    log()
    log("  [Δκ 전체 분해]")
    dk_from_Re = 2 * Re_R_arr / DELTA
    dk_from_ImSq = Im_R_sq_arr
    dk_from_Re_sq = np.array([d['Re_R']**2 for d in decomp])
    a_dk_Re, r2_dkRe, p_dkRe = regress(dk_from_Re, "2·Re[R]/δ  [1차 보정]")
    a_dk_ImSq, r2_dkImSq, p_dkImSq = regress(dk_from_ImSq, "Im[R]²  [2차 보정]")
    a_dk_total, r2_dkT, p_dkT = regress(
        np.array([d['dk_formula'] for d in decomp]), "Δκ_formula (전체)")

    # 이론 해석
    log()
    log("  ─────────────────────────────────────────────────")
    log("  이론적 예측:")
    log("    G(n) ≈ log(t/4π)/2 = log(t/2π)/2 - log(2)/2 → slope ≈ 0.5")
    log("    S_all(n) 밀도 적분: ∫ δ/(δ²+u²) ρ(t) du → ρ(t)·π/2 ... ")
    log("    단, S_all에서 극(1/δ)을 제거했으므로 실질 기여는 작음")
    log("    Im[R]² ∝ [log(t/2π)]^2 또는 a·log(t/2π)?")
    log()
    log(f"  실측 기울기 요약:")
    log(f"    a_G          = {a_G:.4f}  (감마+극)")
    log(f"    a_S_all      = {a_S_all:.4f}  (전체 영점 실부)")
    log(f"    a_Re_R       = {a_Re_R:.4f}  (Re[R] = G + S_all)")
    log(f"    a_Im_R       = {a_Im_R:.4f}  (Im[R])")
    log(f"    a_Im_R_sq    = {a_Im_R_sq:.4f}  (Im[R]²)")
    log(f"    a_dk_Re      = {a_dk_Re:.4f}  (2·Re[R]/δ 기여)")
    log(f"    a_dk_ImSq    = {a_dk_ImSq:.4f}  (Im[R]² 기여)")
    log(f"    a_dk_total   = {a_dk_total:.4f}  (Δκ 전체, 목표: 1.59)")
    log()
    log(f"  기울기 합계 검증: a_dk_Re + a_dk_ImSq = {a_dk_Re + a_dk_ImSq:.4f}")
    log(f"  실측 #32 기울기: 1.5945")

    # 주 기여 식별
    if abs(a_dk_ImSq) > abs(a_dk_Re):
        log()
        log(f"  ★ 핵심 발견: Δκ의 log(t) 기울기는 주로 Im[R]² 항에서 기인.")
        log(f"     a_dk_ImSq = {a_dk_ImSq:.4f} > a_dk_Re = {a_dk_Re:.4f}")
    else:
        log()
        log(f"  ★ 핵심 발견: Δκ의 log(t) 기울기는 주로 2·Re[R]/δ 항에서 기인.")
        log(f"     a_dk_Re = {a_dk_Re:.4f} > a_dk_ImSq = {a_dk_ImSq:.4f}")

    return {
        'a_G': a_G, 'a_S_all': a_S_all, 'a_Re_R': a_Re_R,
        'a_Im_R': a_Im_R, 'a_Im_R_sq': a_Im_R_sq,
        'a_dk_Re': a_dk_Re, 'a_dk_ImSq': a_dk_ImSq, 'a_dk_total': a_dk_total,
    }


# ═══════════════════════════════════════════════════════════════════════════
# 파트 D: N=7000 이상치 해부
# ═══════════════════════════════════════════════════════════════════════════

def part_d_outlier_n7000():
    """
    N=7000 주변 ±10개 영점 구조 분석
    → Δκ=14.08 이상치의 원인 규명 (비정상 근접쌍?)
    """
    log()
    log("=" * 70)
    log("파트 D: N=7000 이상치 해부 (Δκ=14.08, 예측 ~9-10)")
    log("=" * 70)
    log("  목적: 7000번째 영점 근방의 간격 구조로 이상치 원인 규명")
    log()

    n_center = 7000
    window = 10  # 양쪽 10개

    # 영점 수집
    log(f"  zetazero({n_center-window} ~ {n_center+window}) 수집 중...")
    mpmath.mp.dps = 80
    zeros_7000 = []
    t_start = time.time()

    for m in range(n_center - window, n_center + window + 1):
        try:
            zz = mpmath.zetazero(m)
            t_m = float(zz.imag)
            zeros_7000.append((m, t_m))
        except Exception as e:
            log(f"  WARNING zetazero({m}): {e}")

    log(f"  수집 완료: {len(zeros_7000)}개  ({time.time()-t_start:.1f}초)")
    log()

    # 간격 분석
    if len(zeros_7000) < 3:
        log("  ERROR: 영점 부족")
        return

    log(f"  {'m':>6}  {'t_m':>14}  {'gap(m→m+1)':>14}  {'unfolded':>10}  비고")
    log("  " + "-" * 56)

    t_vals = [t for _, t in zeros_7000]
    gaps = np.diff(t_vals)
    # unfolded gap = gap × density = gap × log(t/2π)/(2π)
    mid_t = [(t_vals[i] + t_vals[i+1])/2 for i in range(len(t_vals)-1)]
    density = [np.log(t/(2*np.pi))/(2*np.pi) for t in mid_t]
    unfolded_gaps = [g * d for g, d in zip(gaps, density)]

    min_gap_idx = int(np.argmin(gaps))

    for i, (m, t_m) in enumerate(zeros_7000):
        if i < len(gaps):
            g = gaps[i]
            ug = unfolded_gaps[i]
            note = "★ 최소 간격" if i == min_gap_idx else ""
            note += "  ← N=7000" if m == n_center else ""
            log(f"  {m:>6}  {t_m:>14.6f}  {g:>14.6f}  {ug:>10.4f}  {note}")
        else:
            log(f"  {m:>6}  {t_m:>14.6f}  {'(끝)':>14}")

    log()
    log(f"  간격 통계:")
    log(f"    평균 gap = {np.mean(gaps):.4f}")
    log(f"    최소 gap = {np.min(gaps):.4f}  (인덱스 {zeros_7000[min_gap_idx][0]}-{zeros_7000[min_gap_idx+1][0]})")
    log(f"    최대 gap = {np.max(gaps):.4f}")
    log(f"    std gap  = {np.std(gaps):.4f}")
    log()

    # N=7000에서 ±5 인접 영점이 Δκ에 기여하는 부분
    t_7000_idx = next((i for i, (m, _) in enumerate(zeros_7000) if m == n_center), None)
    if t_7000_idx is not None:
        t_7000 = zeros_7000[t_7000_idx][1]
        s_7000 = mpmath.mpf('0.5') + mpmath.mpf(str(DELTA)) + 1j * mpmath.mpf(str(t_7000))

        log(f"  N=7000 (t={t_7000:.4f}) 에서 인접 영점 기여:")
        log(f"  {'m':>6}  {'t_m':>14}  {'gap |7000-m|':>14}  {'δ/(δ²+gap²)':>14}  기여 순위")
        log("  " + "-" * 60)

        contributions = []
        for m, t_m in zeros_7000:
            if m == n_center:
                continue
            gap = t_7000 - t_m
            contrib_re = DELTA / (DELTA**2 + gap**2)
            contributions.append((m, t_m, gap, contrib_re))

        contributions.sort(key=lambda x: abs(x[3]), reverse=True)

        for rank, (m, t_m, gap, contrib) in enumerate(contributions[:10]):
            log(f"  {m:>6}  {t_m:>14.6f}  {gap:>14.6f}  {contrib:>14.6f}  #{rank+1}")

        # Δκ에서 S_local 기여 추산
        S_local_7000 = sum(DELTA/(DELTA**2 + (t_7000 - tm)**2)
                          for _, tm in zeros_7000 if _ != n_center)
        log()
        log(f"  S_local_7000 (±{window}개 인접 영점 합) = {S_local_7000:.4f}")
        log(f"  이론 Δκ ≈ 2·S_local/δ + Im[R]²:")

        # 필요하다면 Im[R]도 계산
        try:
            mpmath.mp.dps = 80
            xi_ld_7000 = xi_log_deriv_analytic(s_7000)
            R_7000 = xi_ld_7000 - mpmath.mpf(1)/DELTA
            Re_R_7000 = float(mpmath.re(R_7000))
            Im_R_7000 = float(mpmath.im(R_7000))
            dk_7000 = float(abs(xi_ld_7000)**2) - 1/DELTA**2
            log(f"    Re[R] = {Re_R_7000:.6f}")
            log(f"    Im[R] = {Im_R_7000:.4f}")
            log(f"    2·Re[R]/δ = {2*Re_R_7000/DELTA:.4f}")
            log(f"    Im[R]²    = {Im_R_7000**2:.4f}")
            log(f"    Δκ_total  = {dk_7000:.4f}  (#32: 14.0846)")
        except Exception as e:
            log(f"    ERROR: {e}")


# ═══════════════════════════════════════════════════════════════════════════
# 파트 E: Im[R] 성질 상세 분석
# ═══════════════════════════════════════════════════════════════════════════

def part_e_imR_analysis(decomp):
    """
    Im[R(s)]가 Δκ를 지배하는지 검증 (α≈0의 이론적 설명)
    δ-의존성 포함
    """
    log()
    log("=" * 70)
    log("파트 E: Im[R] 분석 — δ-독립성의 이론적 설명")
    log("=" * 70)
    log("  가설: Δκ ≈ Im[R]² (Re[R]/δ 항은 상쇄됨)")
    log("  예측: Im[R]은 δ에 거의 독립 (인접 영점 간격에 의존)")
    log()

    log(f"  {'N':>6}  {'t':>10}  {'Re[R]':>10}  {'Im[R]':>10}"
        f"  {'2ReR/δ':>10}  {'Im²':>10}  {'Δκ':>10}  {'Im²/Δκ%':>9}")
    log("  " + "-" * 90)

    frac_imr = []
    for d in decomp:
        Re_R = d['Re_R']
        Im_R = d['Im_R']
        term_Re = 2*Re_R/DELTA
        term_Im = Im_R**2
        term_ReR2 = Re_R**2
        dk = d['dk_formula']
        frac = term_Im / max(abs(dk), 0.001) * 100
        frac_imr.append(frac)
        log(f"  {d['N']:>6}  {d['t']:>10.2f}  {Re_R:>10.6f}  {Im_R:>10.4f}"
            f"  {term_Re:>10.4f}  {term_Im:>10.4f}  {dk:>10.4f}  {frac:>9.1f}%")

    log()
    log(f"  Im[R]²/Δκ 비율: 평균 {np.mean(frac_imr):.1f}%, 중앙 {np.median(frac_imr):.1f}%")
    log()

    if np.mean(frac_imr) > 70:
        log("  ★★ 핵심 결론: Δκ는 주로 Im[R]²에 의해 결정됨.")
        log("     즉 Δκ ≈ Im[R(s)]² ≈ [Im ζ'/ζ(s)]²")
        log("     Im[ζ'/ζ] = -Σ_m (t_n-t_m)/(δ²+(t_n-t_m)²) → δ-독립 확인")
        log("     기울기 a≈1.59은 Im[R]의 log(t) 성장에서 기인")
    else:
        log("  Im[R]²과 2Re[R]/δ 모두 유의한 기여.")

    # Im[R] vs log(t) 회귀
    log_t = np.array([np.log(d['t']/(2*np.pi)) for d in decomp])
    Im_R_sq = np.array([d['Im_R']**2 for d in decomp])
    slope, intercept, r, p, se = scipy_stats.linregress(log_t, Im_R_sq)
    log()
    log(f"  Im[R]² vs log(t/2π): slope={slope:.4f} ± {se:.4f},  R²={r**2:.4f},  p={p:.2e}")
    log(f"  Δκ ≈ Im[R]² 가정 시: a_Δκ ≈ {slope:.4f} (실측: 1.5945)")


# ═══════════════════════════════════════════════════════════════════════════
# 메인
# ═══════════════════════════════════════════════════════════════════════════

def main():
    t0 = time.time()
    log(f"결과 #33: Hadamard 전개 검증 — Δκ 이론 예측 대 실측 비교")
    log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()
    log(f"  δ = {DELTA},  1/δ² = {1/DELTA**2:.4f}")
    log(f"  영점 인덱스: {ZERO_INDICES}")
    log()

    # 파트 A: 검증
    partA = part_a_verification()
    save()

    if not partA:
        log("ERROR: 파트 A 실패 — 중단")
        save()
        return

    # 파트 B: 분해
    decomp = decompose_contributions(partA)
    save()

    if not decomp:
        log("ERROR: 파트 B 실패 — 중단")
        save()
        return

    # 파트 C: 기울기 분해
    slopes = part_c_slope_decomposition(decomp)
    save()

    # 파트 D: N=7000 이상치
    part_d_outlier_n7000()
    save()

    # 파트 E: Im[R] 분석
    part_e_imR_analysis(decomp)
    save()

    # ─────────────────────────────────────────────
    # 종합 판정
    # ─────────────────────────────────────────────
    errs = [d['err_pct'] for d in decomp]
    mean_err = np.mean(errs)
    a_total = slopes['a_dk_total']

    log()
    log("=" * 70)
    log("★ 종합 판정 — 결과 #33")
    log("=" * 70)
    log()
    log(f"  파트 A: Δκ_direct ≈ Δκ_measured  (동일 공식, 재현 확인)")
    log()
    log(f"  파트 B: Δκ_formula vs Δκ_measured")
    log(f"    평균 오차 {mean_err:.2f}%")
    if mean_err < 10.0:
        log(f"    ✅ 강한 양성: 오차 < 10% → Hadamard 분해 정확")
    elif mean_err < 20.0:
        log(f"    ✅ 양성: 오차 < 20%")
    else:
        log(f"    ❌ 음성: 오차 > 20% — 누락 항 존재")
    log()
    log(f"  파트 C: 기울기 분해")
    log(f"    a_G(감마) = {slopes['a_G']:.4f}")
    log(f"    a_S_all   = {slopes['a_S_all']:.4f}")
    log(f"    a_Im_R²   = {slopes['a_Im_R_sq']:.4f}  (Δκ의 주 기여)")
    log(f"    a_dk_전체 = {a_total:.4f}  (목표: 1.59)")
    diff = abs(a_total - 1.5945)
    if diff < 0.2:
        log(f"    ✅ 기울기 합 1.59 ± 0.2 이내")
    else:
        log(f"    ⚠️ 기울기 합 편차 {diff:.3f}")
    log()
    log(f"  이론-실측 기울기 오차: {diff:.4f} ({diff/1.5945*100:.1f}%)")
    log()
    if mean_err < 10.0 and diff < 0.2:
        log("  ★★ 강한 양성: Hadamard 분해로 Δκ의 수치와 기울기 모두 재현.")
    elif mean_err < 20.0 or diff < 0.3:
        log("  ★ 양성: Hadamard 구조 부분 확인. 추가 항 필요 가능.")
    else:
        log("  ⚠️ 부분 양성: 이론 분해 개선 필요.")

    log()
    log(f"총 소요: {time.time()-t0:.1f}초")
    log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"결과: {RESULTS_PATH}")
    save()


if __name__ == '__main__':
    main()
