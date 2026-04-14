#!/usr/bin/env python3
"""
결과 #34: Δκ 잔차와 GUE 간격 통계의 상관
===========================================
목적: 사이클 48 수학자 지시
  - N=100 등간격 영점에서 Δκ + NNS 수집
  - 잔차 res_n = Δκ_n - (a·log(t/2π) + b) vs NNS_unfolded 상관
  - 2변수 회귀 (3종 f(NNS)) → R² > 0.92 목표
  - Im[R] 분포 vs GUE Wigner surmise KS 검정
  - 이상치 해부: |res| > 2σ → NNS 분포

이론:
  ξ'/ξ(s) = G(s) + ζ'/ζ(s)
  R(s) = ξ'/ξ(s) - 1/δ
  Δκ = |ξ'/ξ|² - 1/δ² = 2·Re[R]/δ + |R|²
  NNS_unfolded = NNS_min × log(t/(2π))/(2π)
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats as scipy_stats
from scipy import integrate as scipy_integrate

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_PATH = os.path.join(BASE_DIR, '..', 'results', 'curvature_gue_residual.txt')
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

DELTA = 0.03   # δ 오프셋
KAPPA_BASE = 1.0 / DELTA**2  # 1/δ² ≈ 1111.11

# 100개 등간격 영점 인덱스
RAW_INDICES = np.linspace(1, 20000, 100).astype(int)
# N<5 (t<33) 제외 — 밀도 근사 불안정
ZERO_INDICES = [int(n) for n in RAW_INDICES if n >= 5]

print(f"총 영점 수: {len(ZERO_INDICES)}개  (5~20000)", flush=True)
print(f"인덱스: {ZERO_INDICES[:5]}...{ZERO_INDICES[-5:]}", flush=True)


def get_dps(t):
    """정밀도 설정"""
    if t > 30000: return 120
    if t > 10000: return 100
    if t > 5000:  return 80
    return max(50, int(30 + t / 200))


# ═══════════════════════════════════════════════════════════════════════════
# 해석적 계산 함수 (bundle_utils.py 함정 회피 — 직접 구현)
# ═══════════════════════════════════════════════════════════════════════════

def G_component(s):
    """
    G(s) = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2
    ξ'/ξ(s)의 감마+극 기여 (복소수 반환)
    """
    return (mpmath.mpf(1)/s
            + mpmath.mpf(1)/(s - 1)
            - mpmath.log(mpmath.pi)/2
            + mpmath.digamma(s/2)/2)


def zeta_log_deriv(s):
    """
    ζ'(s)/ζ(s): h=1e-6 중앙 차분 (mpmath.diff 금지)
    """
    z = mpmath.zeta(s)
    if abs(z) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e8, 0)
    h = mpmath.mpf('1e-6')
    zd = (mpmath.zeta(s + h) - mpmath.zeta(s - h)) / (2 * h)
    return zd / z


def xi_log_deriv(s):
    """ξ'/ξ(s) = G(s) + ζ'/ζ(s) — 해석적 공식"""
    return G_component(s) + zeta_log_deriv(s)


def get_zetazero_t(n, dps=50):
    """n번째 리만 제타 영점의 허수부 반환"""
    with mpmath.workdps(dps):
        z = mpmath.zetazero(n)
        return float(mpmath.im(z))


def compute_kappa_and_R(t, delta=DELTA):
    """
    κ = |ξ'/ξ(s)|², R = ξ'/ξ(s) - 1/δ, Δκ = κ - 1/δ²
    반환: (kappa, delta_kappa, Re_R, Im_R)
    """
    dps = get_dps(t)
    with mpmath.workdps(dps):
        s = mpmath.mpc(0.5 + delta, t)
        xi_ld = xi_log_deriv(s)
        kappa = float(abs(xi_ld)**2)
        # R = ξ'/ξ - 1/δ (1/δ는 실수)
        R = xi_ld - mpmath.mpf(1) / delta
        Re_R = float(mpmath.re(R))
        Im_R = float(mpmath.im(R))
    delta_kappa = kappa - KAPPA_BASE
    return kappa, delta_kappa, Re_R, Im_R


# ═══════════════════════════════════════════════════════════════════════════
# GUE Wigner surmise
# ═══════════════════════════════════════════════════════════════════════════

def gue_pdf(s):
    """GUE Wigner surmise: p(s) = (32/π²)s²exp(-4s²/π)"""
    return (32.0 / np.pi**2) * s**2 * np.exp(-4.0 * s**2 / np.pi)


def gue_cdf_scalar(x):
    """GUE CDF 수치 적분 (스칼라 입력)"""
    if x <= 0:
        return 0.0
    result, _ = scipy_integrate.quad(gue_pdf, 0, x)
    return result

def gue_cdf(x):
    """GUE CDF 수치 적분 — 스칼라 또는 배열 입력 지원"""
    if np.isscalar(x):
        return gue_cdf_scalar(float(x))
    return np.array([gue_cdf_scalar(float(xi)) for xi in x])


# ═══════════════════════════════════════════════════════════════════════════
# 파트 A: 데이터 수집 — 100개 영점
# ═══════════════════════════════════════════════════════════════════════════

def part_a_collect():
    log("=" * 70)
    log("파트 A: 데이터 수집 — 100개 영점 (Δκ + NNS)")
    log("=" * 70)
    log(f"  δ = {DELTA},  1/δ² = {KAPPA_BASE:.4f}")
    log(f"  인덱스: {len(ZERO_INDICES)}개 등간격 (5 ~ 20000)")
    log()

    t_start = time.time()
    data = []
    zero_cache = {}   # 캐시 {n: t_n}
    fail_count = 0

    def cached_zetazero(n, dps=50):
        if n not in zero_cache:
            try:
                zero_cache[n] = get_zetazero_t(n, dps=dps)
            except Exception as e:
                print(f"  WARNING: zetazero({n}) 실패: {e}", flush=True)
                zero_cache[n] = None
        return zero_cache[n]

    log(f"  {'#':>4}  {'N':>6}  {'t_n':>12}  {'Δκ':>10}  {'gap_L':>8}  {'gap_R':>8}  "
        f"{'NNS_unf':>8}  {'Re_R':>8}  {'Im_R':>8}  상태")
    log("  " + "-" * 95)

    for idx, N in enumerate(ZERO_INDICES):
        t0 = time.time()
        try:
            # dps 결정 (t 미리 알 수 없어 임시로 50)
            t_n = cached_zetazero(N, dps=50)
            if t_n is None:
                fail_count += 1
                continue

            # 고 t는 dps 재조정
            dps_needed = get_dps(t_n)
            if dps_needed > 50:
                zero_cache[N] = get_zetazero_t(N, dps=dps_needed)
                t_n = zero_cache[N]

            # 이웃 영점
            t_prev = cached_zetazero(N - 1, dps=min(dps_needed, 80))
            t_next = cached_zetazero(N + 1, dps=min(dps_needed, 80))

            if t_prev is None or t_next is None:
                log(f"  {'':>4}  {N:>6}  {t_n:>12.4f}  {'--':>10}  이웃 영점 실패")
                fail_count += 1
                continue

            gap_left  = t_n - t_prev
            gap_right = t_next - t_n
            NNS_min   = min(gap_left, gap_right)

            # 영점 밀도 D(t) = log(t/(2π))/(2π) → 평균 간격 = 2π/log(t/2π)
            density = np.log(t_n / (2 * np.pi)) / (2 * np.pi)
            NNS_unfolded = NNS_min * density  # = NNS_min × log(t/(2π))/(2π)

            # κ, Δκ, R 계산
            kappa, dk, Re_R, Im_R = compute_kappa_and_R(t_n)

            elapsed = time.time() - t0
            status = "✅" if elapsed < 30 else "⚠️ 느림"

            log(f"  {idx+1:>4}  {N:>6}  {t_n:>12.4f}  {dk:>10.4f}  "
                f"{gap_left:>8.4f}  {gap_right:>8.4f}  {NNS_unfolded:>8.4f}  "
                f"{Re_R:>8.4f}  {Im_R:>8.4f}  {status}")

            data.append({
                'N': N, 't': t_n,
                'kappa': kappa, 'delta_kappa': dk,
                'gap_left': gap_left, 'gap_right': gap_right,
                'NNS_min': NNS_min, 'NNS_unfolded': NNS_unfolded,
                'Re_R': Re_R, 'Im_R': Im_R,
            })

        except Exception as e:
            print(f"  WARNING N={N}: {e}", flush=True)
            fail_count += 1

        # 중간 저장 (10개마다)
        if (idx + 1) % 10 == 0:
            save()
            elapsed_total = time.time() - t_start
            log(f"\n  [{idx+1}/{len(ZERO_INDICES)}] 경과 {elapsed_total:.0f}초\n")

    log()
    log(f"  수집 완료: {len(data)}개 성공, {fail_count}개 실패")
    log(f"  총 소요: {time.time() - t_start:.1f}초")

    if fail_count > len(ZERO_INDICES) // 2:
        log("⚠️ 절반 이상 실패 — 데이터 신뢰도 낮음. 분석 중단.")
        return data

    return data


# ═══════════════════════════════════════════════════════════════════════════
# 파트 B: Δκ 잔차 vs NNS 상관
# ═══════════════════════════════════════════════════════════════════════════

def part_b_residual_nns(data):
    log()
    log("=" * 70)
    log("파트 B: Δκ 잔차 vs NNS 상관")
    log("=" * 70)

    t_vals   = np.array([d['t'] for d in data])
    dk_vals  = np.array([d['delta_kappa'] for d in data])
    nns_unf  = np.array([d['NNS_unfolded'] for d in data])
    nns_min  = np.array([d['NNS_min'] for d in data])

    log_t = np.log(t_vals / (2 * np.pi))

    # 1변수 OLS: Δκ ~ a·log(t/2π) + b
    slope1, intercept1, r1, p1, se1 = scipy_stats.linregress(log_t, dk_vals)
    res = dk_vals - (slope1 * log_t + intercept1)

    log(f"\n  1변수 OLS: Δκ ~ a·log(t/2π) + b")
    log(f"    a = {slope1:.4f} ± {se1:.4f}")
    log(f"    b = {intercept1:.4f}")
    log(f"    R² = {r1**2:.4f}")
    log(f"    p  = {p1:.2e}")
    log()

    # res vs NNS_unfolded
    r_pearson, p_pearson = scipy_stats.pearsonr(res, nns_unf)
    r_spearman, p_spearman = scipy_stats.spearmanr(res, nns_unf)

    log(f"  res_n vs NNS_unfolded:")
    log(f"    Pearson  r = {r_pearson:.4f}  p = {p_pearson:.4e}")
    log(f"    Spearman ρ = {r_spearman:.4f}  p = {p_spearman:.4e}")

    # 이론 예측: res ∝ 1/NNS_min² (Lorentzian²)
    # Lorentzian: δ/(δ²+g²) → Δκ 잔차 ∝ [δ/(δ²+g²)]²
    lorentz2 = (DELTA / (DELTA**2 + nns_min**2))**2
    r_lor, p_lor = scipy_stats.pearsonr(res, lorentz2)
    r_lor_sp, p_lor_sp = scipy_stats.spearmanr(res, lorentz2)

    log()
    log(f"  res_n vs 1/NNS_min² (Lorentzian² — 이론 동기):")
    log(f"    1/NNS² 범위: {1/nns_min.max()**2:.4f} ~ {1/nns_min.min()**2:.1f}")
    log(f"    Pearson  r = {r_lor:.4f}  p = {p_lor:.4e}")
    log(f"    Spearman ρ = {r_lor_sp:.4f}  p = {p_lor_sp:.4e}")

    # 상위/하위 10개 (NNS_unfolded 기준)
    sort_idx = np.argsort(nns_unf)
    log()
    log(f"  NNS_unfolded 하위 10개 (가장 가까운 이웃):")
    log(f"    {'N':>6}  {'t':>12}  {'NNS_unf':>10}  {'NNS_min':>10}  {'res':>10}")
    for i in sort_idx[:10]:
        log(f"    {data[i]['N']:>6}  {data[i]['t']:>12.4f}  {nns_unf[i]:>10.4f}  "
            f"{nns_min[i]:>10.4f}  {res[i]:>10.4f}")
    log()
    log(f"  NNS_unfolded 상위 10개 (가장 먼 이웃):")
    for i in sort_idx[-10:][::-1]:
        log(f"    {data[i]['N']:>6}  {data[i]['t']:>12.4f}  {nns_unf[i]:>10.4f}  "
            f"{nns_min[i]:>10.4f}  {res[i]:>10.4f}")

    # 판정
    log()
    if abs(r_pearson) > 0.5 and p_pearson < 0.01:
        log(f"  ★★ 강한 양성: |r| = {abs(r_pearson):.3f} > 0.5, p = {p_pearson:.2e} < 0.01")
    elif abs(r_pearson) > 0.3 and p_pearson < 0.05:
        log(f"  ★ 양성: |r| = {abs(r_pearson):.3f} > 0.3, p = {p_pearson:.2e} < 0.05")
    else:
        log(f"  음성: |r| = {abs(r_pearson):.3f} ≤ 0.3 또는 p = {p_pearson:.2e} > 0.05")
        log(f"  → 잔차는 NNS와 무관. 다른 원인 탐색 필요.")

    return res, slope1, intercept1, r1**2


# ═══════════════════════════════════════════════════════════════════════════
# 파트 C: 2변수 회귀
# ═══════════════════════════════════════════════════════════════════════════

def part_c_bivariate(data, res, R2_1var):
    log()
    log("=" * 70)
    log("파트 C: 2변수 회귀 — Δκ ~ a·log(t/2π) + b + c·f(NNS)")
    log("=" * 70)

    t_vals   = np.array([d['t'] for d in data])
    dk_vals  = np.array([d['delta_kappa'] for d in data])
    nns_unf  = np.array([d['NNS_unfolded'] for d in data])
    nns_min  = np.array([d['NNS_min'] for d in data])
    log_t    = np.log(t_vals / (2 * np.pi))
    N        = len(data)

    f_candidates = {
        'f₁ = δ²/(δ²+NNS_min²)²  (Lorentzian²)': (DELTA**2 / (DELTA**2 + nns_min**2))**2,
        'f₂ = 1/NNS_unfolded':                    1.0 / nns_unf,
        'f₃ = log(NNS_unfolded)':                  np.log(nns_unf),
    }

    log(f"\n  1변수 R² = {R2_1var:.4f}  (기준)")
    log(f"  {'f':^35}  {'R²_2var':>8}  {'ΔR²':>8}  {'AIC':>10}  {'F-stat':>10}  {'p_F':>10}  판정")
    log("  " + "-" * 100)

    best_r2 = R2_1var
    best_name = "1변수"

    for fname, f_vals in f_candidates.items():
        # NaN/Inf 제거
        mask = np.isfinite(f_vals) & np.isfinite(log_t) & np.isfinite(dk_vals)
        if mask.sum() < 10:
            log(f"  {fname:<35}  유효 데이터 부족 ({mask.sum()}개)")
            continue

        X = np.column_stack([log_t[mask], f_vals[mask], np.ones(mask.sum())])
        y = dk_vals[mask]
        n_obs = mask.sum()
        p_params = 3

        # OLS
        try:
            coeffs, residuals_ss, rank, sv = np.linalg.lstsq(X, y, rcond=None)
        except Exception as e:
            log(f"  {fname:<35}  OLS 실패: {e}")
            continue

        y_pred = X @ coeffs
        SS_res = np.sum((y - y_pred)**2)
        SS_tot = np.sum((y - y.mean())**2)
        R2_2var = 1 - SS_res / SS_tot if SS_tot > 0 else 0.0

        # AIC
        if SS_res > 0 and n_obs > p_params:
            sigma2 = SS_res / n_obs
            AIC = n_obs * np.log(sigma2) + 2 * p_params
        else:
            AIC = np.inf

        # F-test (vs 1변수 모델)
        # H0: c=0 (extra predictor adds nothing)
        X1 = np.column_stack([log_t[mask], np.ones(mask.sum())])
        coeffs1, _, _, _ = np.linalg.lstsq(X1, y, rcond=None)
        y1_pred = X1 @ coeffs1
        SS_res1 = np.sum((y - y1_pred)**2)

        if SS_res1 > SS_res and n_obs > p_params:
            F_stat = ((SS_res1 - SS_res) / 1) / (SS_res / (n_obs - p_params))
            p_F = 1 - scipy_stats.f.cdf(F_stat, 1, n_obs - p_params)
        else:
            F_stat = 0.0
            p_F = 1.0

        delta_R2 = R2_2var - R2_1var

        if R2_2var > best_r2:
            best_r2 = R2_2var
            best_name = fname

        sign = "★★" if R2_2var > 0.92 else ("★" if R2_2var > 0.88 else "")
        log(f"  {fname:<35}  {R2_2var:>8.4f}  {delta_R2:>+8.4f}  {AIC:>10.2f}  "
            f"{F_stat:>10.3f}  {p_F:>10.4f}  {sign}")

    log()
    log(f"  최선 모델: {best_name}  (R² = {best_r2:.4f})")
    if best_r2 > 0.92:
        log(f"  ★★ 강한 양성: R²_2var = {best_r2:.4f} > 0.92")
    elif best_r2 > 0.88:
        log(f"  ★ 양성: R²_2var = {best_r2:.4f} > 0.88")
    else:
        log(f"  음성: R²_2var = {best_r2:.4f} ≤ 0.88 → NNS가 잔차를 충분히 설명하지 못함")

    return best_r2


# ═══════════════════════════════════════════════════════════════════════════
# 파트 D: Im[R] 분포 vs GUE 예측
# ═══════════════════════════════════════════════════════════════════════════

def part_d_imR_gue(data):
    log()
    log("=" * 70)
    log("파트 D: Im[R] 분포 vs GUE Wigner surmise")
    log("=" * 70)

    Im_R_vals = np.array([d['Im_R'] for d in data])
    nns_unf   = np.array([d['NNS_unfolded'] for d in data])
    nns_min   = np.array([d['NNS_min'] for d in data])

    log(f"\n  Im[R] 기초 통계:")
    log(f"    N = {len(Im_R_vals)}")
    log(f"    평균  = {Im_R_vals.mean():.4f}")
    log(f"    표준편차 = {Im_R_vals.std():.4f}")
    log(f"    최솟값 = {Im_R_vals.min():.4f}")
    log(f"    최댓값 = {Im_R_vals.max():.4f}")

    # Im[R]²
    ImR2 = Im_R_vals**2
    log(f"\n  Im[R]² 기초 통계:")
    log(f"    평균  = {ImR2.mean():.4f}")
    log(f"    표준편차 = {ImR2.std():.4f}")

    # 히스토그램 (10 빈)
    log(f"\n  Im[R]² 히스토그램 (10 빈):")
    counts, edges = np.histogram(ImR2, bins=10)
    for i, (c, lo, hi) in enumerate(zip(counts, edges[:-1], edges[1:])):
        bar = '█' * c
        log(f"    [{lo:7.3f}, {hi:7.3f})  {c:>3}  {bar}")

    # KS 검정: NNS_unfolded vs GUE Wigner surmise
    log(f"\n  KS 검정: NNS_unfolded vs GUE Wigner surmise")
    log(f"    GUE: p(s) = (32/π²)s²exp(-4s²/π)")

    # 필터: NNS_unfolded > 0 (음수는 불가)
    mask_pos = nns_unf > 0
    nns_pos = nns_unf[mask_pos]

    ks_stat, ks_p = scipy_stats.kstest(nns_pos, gue_cdf)
    log(f"    KS 통계량 D = {ks_stat:.4f}")
    log(f"    p-value = {ks_p:.4e}")

    if ks_p > 0.05:
        log(f"    → GUE 분포와 일관적 (p={ks_p:.3f} > 0.05)")
    else:
        log(f"    → GUE 분포와 불일치 (p={ks_p:.3f} < 0.05)")

    # Im[R] vs NNS_min 상관 (이론: Im[R] ≈ Σ (γ_n-γ_m)/(δ²+(γ_n-γ_m)²))
    # 주 기여항: 가장 가까운 이웃. 이론적으로 Im[R] ∝ 1/NNS_min (NNS_min >> δ)
    r_imR_nns, p_imR_nns = scipy_stats.pearsonr(np.abs(Im_R_vals), 1.0/nns_min)
    log(f"\n  |Im[R]| vs 1/NNS_min 상관:")
    log(f"    Pearson r = {r_imR_nns:.4f}  p = {p_imR_nns:.4e}")
    if abs(r_imR_nns) > 0.3:
        log(f"    → Im[R]은 NNS에 의해 부분적으로 설명됨 (이론 일관)")
    else:
        log(f"    → Im[R]과 1/NNS_min 상관 약함")

    # NNS_unfolded 히스토그램 vs GUE 이론값
    log(f"\n  NNS_unfolded 히스토그램 vs GUE 이론 (5 빈):")
    bin_edges = np.array([0, 0.4, 0.8, 1.2, 1.8, 3.0])
    emp_counts, _ = np.histogram(nns_pos, bins=bin_edges)
    n_total = len(nns_pos)
    log(f"    {'구간':^20}  {'실측':>6}  {'이론(%)':>8}  {'실측(%)':>8}")
    for i, (lo, hi) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
        theory_prob = gue_cdf(hi) - gue_cdf(lo)
        emp_prob = emp_counts[i] / n_total
        log(f"    [{lo:.1f}, {hi:.1f})  {'':>8}  {emp_counts[i]:>6}  "
            f"{theory_prob*100:>8.2f}%  {emp_prob*100:>8.2f}%")

    return Im_R_vals, nns_unf


# ═══════════════════════════════════════════════════════════════════════════
# 파트 E: 이상치 해부
# ═══════════════════════════════════════════════════════════════════════════

def part_e_outliers(data, res):
    log()
    log("=" * 70)
    log("파트 E: 이상치 해부 — |res| > 2σ")
    log("=" * 70)

    nns_unf = np.array([d['NNS_unfolded'] for d in data])
    sigma = res.std()
    threshold = 2 * sigma

    outlier_mask = np.abs(res) > threshold
    n_outliers = outlier_mask.sum()

    log(f"\n  σ(res) = {sigma:.4f}")
    log(f"  임계값 = 2σ = {threshold:.4f}")
    log(f"  이상치 수: {n_outliers}개 / {len(data)}개 ({n_outliers/len(data)*100:.1f}%)")

    if n_outliers == 0:
        log("  이상치 없음 — 잔차 분포 정규적")
        return

    outlier_nns = nns_unf[outlier_mask]
    non_outlier_nns = nns_unf[~outlier_mask]

    log(f"\n  이상치 NNS_unfolded:")
    log(f"    평균  = {outlier_nns.mean():.4f}")
    log(f"    표준편차 = {outlier_nns.std():.4f}")
    log(f"\n  정상 NNS_unfolded:")
    log(f"    평균  = {non_outlier_nns.mean():.4f}")
    log(f"    표준편차 = {non_outlier_nns.std():.4f}")

    # GUE level repulsion 예측: P(s < 0.5)
    p_s_half = gue_cdf(0.5)
    log(f"\n  GUE level repulsion: P(NNS_unf < 0.5) = {p_s_half:.4f} ({p_s_half*100:.2f}%)")
    n_small_nns = (outlier_nns < 0.5).sum()
    log(f"  이상치 중 NNS_unf < 0.5: {n_small_nns}/{n_outliers}개 ({n_small_nns/max(n_outliers,1)*100:.1f}%)")

    log(f"\n  이상치 목록:")
    log(f"  {'N':>6}  {'t':>12}  {'Δκ':>10}  {'res':>10}  {'NNS_unf':>10}  {'NNS_min':>10}")
    for i, d in enumerate(data):
        if outlier_mask[i]:
            log(f"  {d['N']:>6}  {d['t']:>12.4f}  {d['delta_kappa']:>10.4f}  "
                f"{res[i]:>10.4f}  {nns_unf[i]:>10.4f}  {d['NNS_min']:>10.4f}")

    # Mann-Whitney test: outlier NNS vs non-outlier NNS
    if n_outliers >= 3 and len(non_outlier_nns) >= 3:
        stat_mw, p_mw = scipy_stats.mannwhitneyu(outlier_nns, non_outlier_nns,
                                                   alternative='two-sided')
        log(f"\n  Mann-Whitney test (이상치 vs 정상 NNS):")
        log(f"    U = {stat_mw:.1f},  p = {p_mw:.4e}")
        if p_mw < 0.05:
            log(f"    → 이상치는 정상보다 유의하게 다른 NNS 분포 (p < 0.05)")
        else:
            log(f"    → 이상치 NNS 분포와 정상 NNS 분포 유의 차이 없음")


# ═══════════════════════════════════════════════════════════════════════════
# 종합 판정
# ═══════════════════════════════════════════════════════════════════════════

def final_summary(data, res, R2_1var, R2_2var_best):
    log()
    log("=" * 70)
    log("종합 판정")
    log("=" * 70)

    t_vals  = np.array([d['t'] for d in data])
    nns_unf = np.array([d['NNS_unfolded'] for d in data])

    r_pearson, p_pearson = scipy_stats.pearsonr(res, nns_unf)

    log(f"\n  성공 기준:")
    log(f"    강한 양성: |r| > 0.5 & p < 0.01 & R²_2var > 0.92")
    log(f"    양성:     |r| > 0.3 & p < 0.05 & R²_2var > 0.88")
    log(f"    음성:     |r| < 0.2 또는 p > 0.1")
    log()
    log(f"  실측:")
    log(f"    |r|_Pearson = {abs(r_pearson):.4f}")
    log(f"    p_Pearson   = {p_pearson:.2e}")
    log(f"    R²_1var     = {R2_1var:.4f}")
    log(f"    R²_2var_best = {R2_2var_best:.4f}")
    log(f"    ΔR²         = {R2_2var_best - R2_1var:+.4f}")
    log()

    is_strong = abs(r_pearson) > 0.5 and p_pearson < 0.01 and R2_2var_best > 0.92
    is_pos    = abs(r_pearson) > 0.3 and p_pearson < 0.05 and R2_2var_best > 0.88
    is_neg    = abs(r_pearson) < 0.2 or p_pearson > 0.1

    if is_strong:
        log("  ★★ 강한 양성: 잔차가 NNS로 유의하게 설명됨. GUE-곡률 연결 확인.")
    elif is_pos:
        log("  ★ 양성: 잔차-NNS 상관 존재. 2변수 모델 개선 확인.")
    elif is_neg:
        log("  음성: 잔차는 NNS와 무관. 다른 원인 탐색 필요.")
    else:
        log("  중간: 상관 존재하나 기준 미달.")

    log()
    log(f"  N = {len(data)}개 영점")
    log(f"  t 범위: {t_vals.min():.2f} ~ {t_vals.max():.2f}")
    log(f"  NNS_unfolded 범위: {nns_unf.min():.3f} ~ {nns_unf.max():.3f}")


# ═══════════════════════════════════════════════════════════════════════════
# 메인
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    T_TOTAL = time.time()

    log("결과 #34: Δκ 잔차와 GUE 간격 통계의 상관")
    log(f"실행 시각: 2026-04-15 03:16")
    log(f"δ = {DELTA},  1/δ² = {KAPPA_BASE:.4f}")
    log()

    # 파트 A: 데이터 수집
    data = part_a_collect()
    save()

    if len(data) < 20:
        log("⚠️ 데이터 부족 — 분석 중단")
        save()
        sys.exit(1)

    # 파트 B: 잔차 vs NNS 상관
    res, slope1, intercept1, R2_1var = part_b_residual_nns(data)
    save()

    # 파트 C: 2변수 회귀
    R2_2var_best = part_c_bivariate(data, res, R2_1var)
    save()

    # 파트 D: Im[R] 분포 vs GUE
    Im_R_vals, nns_unf = part_d_imR_gue(data)
    save()

    # 파트 E: 이상치 해부
    part_e_outliers(data, res)
    save()

    # 종합 판정
    final_summary(data, res, R2_1var, R2_2var_best)

    log()
    log(f"총 소요 시간: {time.time() - T_TOTAL:.1f}초")
    save()
    log("완료.")
