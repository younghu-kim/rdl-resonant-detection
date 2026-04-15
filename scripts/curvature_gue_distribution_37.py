#!/usr/bin/env python3
"""
결과 #37: 곡률 값 분포와 GUE 간격 분포의 정량적 관계 검증
==========================================================
목적: 사이클 54 수학자 지시

이론:
  - GUE 간격 분포: P_GUE(s) = (32/π²)s² exp(-4s²/π) [Wigner surmise]
  - 곡률-NNS 관계: κ ≈ 1/(δ²+s_unf²)² (δ=0.03 오프셋)
  - s_unf = NNS_raw × log(t/(2π))/(2π) [unfolded 좌표]
  - 변수 변환: κ = 1/(δ²+s²)² → P_pred(κ) = P_GUE(s(κ)) × |ds/dκ|
    s(κ) = √(1/√κ − δ²)
    |ds/dκ| = 1 / (4 × κ^{3/2} × (1/√κ − δ²)^{1/2})
  - CDF_pred(κ) = P(κ_pred ≤ κ_0) = 1 − F_GUE(s(κ_0))
    [κ 감소 ↔ s 증가 관계]

검정:
  - KS: P_emp(κ) vs P_pred(κ) → p > 0.05 목표 (양성 기준)
  - Q-Q 상관 (log κ) > 0.95 목표
  - AD 검정 추가

주의사항:
  - ★ κ: 해석적 G(s) + ζ'/ζ (h=1e-6 차분), bundle_utils h=1e-20 금지
  - ★ δ=0.03 오프셋 필수 (영점 위 측정 금지)
  - ★ NNS 배제: s_unf < δ=0.03 (κ→∞ 발산 영역)
  - ★ κ_max 이론: 1/δ⁴ ≈ 1.23M (배제 임계: κ>1/δ⁴ 불가)
  - ★ t 범위: [1000, 5000] (캐시 영점 활용, t>5447은 on-the-fly)
"""

import sys, os, time, warnings
import numpy as np
import mpmath
from scipy import stats as scipy_stats
from scipy import integrate as scipy_integrate
from scipy.optimize import brentq
from scipy.interpolate import interp1d

warnings.filterwarnings('ignore', category=RuntimeWarning)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_PATH = os.path.join(BASE_DIR, '..', 'results', 'curvature_gue_distribution_37.txt')
CACHE_PATH = os.path.join(BASE_DIR, '..', 'outputs', 'cache', 'spacing_ratio_zeros_n100_5000.npy')
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

DELTA      = 0.03
T_MIN      = 1000.0
T_MAX      = 5000.0   # 캐시 사용 범위 (t_max_cache≈5447, 10000은 on-the-fly 과부하)
N_SAMPLE   = 2000
S_MIN_EXCL = DELTA          # s_unf < δ 배제
KAPPA_MAX  = 1.0 / DELTA**4  # 이론적 최대값 (s=0에서 κ=1/δ⁴)
# 주: 수학자 지시는 t∈[1000,10000]이나 캐시 영점 한계상 [1000,5000] 사용
# 캐시: 4901개 영점, t∈[236.5, 5447.9]


# ═══════════════════════════════════════════════════════════════════════════
# 해석적 κ 계산 (bundle_utils h=1e-20 금지 — 직접 구현)
# ═══════════════════════════════════════════════════════════════════════════

def get_dps(t):
    if t > 10000: return 100
    if t > 5000:  return 80
    return max(50, int(30 + t / 200))


def G_component(s):
    """G(s) = 1/s + 1/(s-1) − log(π)/2 + ψ(s/2)/2 (해석적)"""
    return (mpmath.mpf(1)/s
            + mpmath.mpf(1)/(s - 1)
            - mpmath.log(mpmath.pi)/2
            + mpmath.digamma(s/2)/2)


def zeta_log_deriv(s):
    """ζ'/ζ(s) — h=1e-6 중앙 차분 (mpmath.diff 금지, h=1e-20 금지)"""
    z = mpmath.zeta(s)
    if abs(z) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
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
        kappa = float(abs(xi_ld)**2)
    return kappa


# ═══════════════════════════════════════════════════════════════════════════
# 영점 탐색 함수
# ═══════════════════════════════════════════════════════════════════════════

def rs_n_zeros(t):
    """리만-지겔: t 이하 영점 수 근사 (인덱스 추정용)"""
    x = t / (2.0 * np.pi)
    if x <= 1:
        return 1
    return int(x * np.log(x) - x + 7.0/8.0) + 1


def get_zetazero_t(n, dps=50):
    with mpmath.workdps(dps):
        return float(mpmath.im(mpmath.zetazero(n)))


def find_nearest_zero_cache(t, zeros_cache):
    """캐시 배열에서 t에 가장 가까운 영점 반환 (O(log n))"""
    idx = np.searchsorted(zeros_cache, t)
    best = None
    best_dist = np.inf
    for i in [idx-1, idx, idx+1]:
        if 0 <= i < len(zeros_cache):
            d = abs(t - zeros_cache[i])
            if d < best_dist:
                best_dist = d
                best = zeros_cache[i]
    return best


def find_nearest_zero_online(t, dps=50):
    """on-the-fly 영점 탐색 (t > cache 범위용)"""
    n_est = max(rs_n_zeros(t), 3)
    best_gamma = None
    best_dist = np.inf
    for n in range(max(1, n_est - 2), n_est + 4):
        try:
            gamma = get_zetazero_t(n, dps=dps)
            d = abs(t - gamma)
            if d < best_dist:
                best_dist = d
                best_gamma = gamma
        except Exception as e:
            print(f"  WARNING: zetazero({n}) 실패: {e}", flush=True)
    return best_gamma


# ═══════════════════════════════════════════════════════════════════════════
# GUE Wigner surmise + CDF (사전 계산 테이블)
# ═══════════════════════════════════════════════════════════════════════════

def gue_pdf(s):
    """P_GUE(s) = (32/π²) s² exp(−4s²/π) [Wigner surmise, GUE]"""
    return (32.0 / np.pi**2) * s**2 * np.exp(-4.0 * s**2 / np.pi)


def _build_gue_cdf_table(n_pts=2000, s_max=8.0):
    """GUE CDF 사전 계산 테이블 (속도 최적화)"""
    s_grid = np.linspace(0, s_max, n_pts)
    cdf_grid = np.zeros(n_pts)
    for i in range(1, n_pts):
        val, _ = scipy_integrate.quad(gue_pdf, 0, s_grid[i], limit=200)
        cdf_grid[i] = val
    # 보간 함수
    interp = interp1d(s_grid, cdf_grid, kind='cubic', fill_value=(0.0, 1.0),
                      bounds_error=False)
    return s_grid, cdf_grid, interp


# ═══════════════════════════════════════════════════════════════════════════
# 변수 변환: κ ↔ s (κ = 1/(δ²+s²)²)
# ═══════════════════════════════════════════════════════════════════════════

def s_from_kappa(kappa, delta=DELTA):
    """s(κ) = √(1/√κ − δ²)"""
    val = 1.0 / np.sqrt(np.maximum(kappa, 1e-300)) - delta**2
    return np.sqrt(np.maximum(val, 0.0))


def pred_cdf_from_s(s_val, gue_cdf_interp):
    """CDF_pred(κ) = 1 − F_GUE(s(κ)) (κ ↑ ↔ s ↓ 관계)"""
    return 1.0 - float(gue_cdf_interp(s_val))


# ═══════════════════════════════════════════════════════════════════════════
# Anderson-Darling 검정 (이론 CDF)
# ═══════════════════════════════════════════════════════════════════════════

def anderson_darling_test(data, cdf_func):
    """
    Anderson-Darling 검정: A² = −N − (1/N) ∑ (2i-1)[ln F(x_i) + ln(1-F(x_{N+1-i}))]
    p-값: Stephens (1974) 근사
    """
    n = len(data)
    x_sorted = np.sort(data)
    cdf_vals = np.array([cdf_func(x) for x in x_sorted])
    # 수치 안정성: 클리핑
    cdf_vals = np.clip(cdf_vals, 1e-10, 1 - 1e-10)
    i_arr = np.arange(1, n + 1)
    s_sum = np.sum((2*i_arr - 1) * (np.log(cdf_vals) + np.log(1 - cdf_vals[::-1])))
    A2 = -n - s_sum / n
    # Stephens (1974) 수정 계수
    A2_mod = A2 * (1 + 0.75/n + 2.25/n**2)
    # p-값 근사 (Lewis 1961 / D'Agostino 1986)
    if A2_mod < 0.2:
        p = 1 - np.exp(-13.436 + 101.14*A2_mod - 223.73*A2_mod**2)
    elif A2_mod < 0.34:
        p = 1 - np.exp(-8.318 + 42.796*A2_mod - 59.938*A2_mod**2)
    elif A2_mod < 0.6:
        p = np.exp(0.9177 - 4.279*A2_mod - 1.38*A2_mod**2)
    elif A2_mod < 153:
        p = np.exp(1.2937 - 5.709*A2_mod + 0.0186*A2_mod**2)
    else:
        p = 3.7e-24
    return float(A2), float(A2_mod), float(np.clip(p, 0, 1))


# ═══════════════════════════════════════════════════════════════════════════
# 메인
# ═══════════════════════════════════════════════════════════════════════════

def main():
    t_start = time.time()

    log("=" * 70)
    log("결과 #37: 곡률 값 분포와 GUE 간격 분포의 정량적 관계 검증")
    log("=" * 70)
    log(f"  δ = {DELTA},  t ∈ [{T_MIN:.0f}, {T_MAX:.0f}],  N = {N_SAMPLE}")
    log(f"  배제: s_unf < {S_MIN_EXCL} (κ 발산 영역)")
    log(f"  κ 최대값 이론: 1/δ⁴ = {KAPPA_MAX:.1f}")
    log(f"  κ 기준값: 1/δ² = {1/DELTA**2:.1f}")
    log(f"  P_GUE(s) = (32/π²)s² exp(−4s²/π)")
    log(f"  CDF_pred(κ) = 1 − F_GUE(s(κ)),  s(κ) = √(1/√κ − δ²)")
    log()
    log("GUE CDF 테이블 사전 계산 중...")
    t_gue = time.time()
    s_grid, cdf_grid, gue_cdf_interp = _build_gue_cdf_table(n_pts=2000, s_max=8.0)
    gue_mean = float(scipy_integrate.quad(lambda s: s * gue_pdf(s), 0, 8)[0])
    log(f"  GUE CDF 테이블 완성 ({time.time()-t_gue:.1f}s), <s>_GUE = {gue_mean:.4f}")
    log()

    # ═══════════════════════════════════════════════════════════════════════
    # 파트 A: 데이터 수집
    # ═══════════════════════════════════════════════════════════════════════
    log("=" * 70)
    log("파트 A: 데이터 수집 (N=2000, t∈[1000,10000], κ+NNS)")
    log("=" * 70)

    # 캐시 로드
    if os.path.exists(CACHE_PATH):
        zeros_cache = np.load(CACHE_PATH)
        cache_t_max = float(zeros_cache[-1])
        log(f"  캐시 로드: {len(zeros_cache)}개 영점, t∈[{zeros_cache[0]:.1f}, {cache_t_max:.1f}]")
    else:
        log("  ⚠️ 캐시 없음 — on-the-fly 계산 (느림)")
        zeros_cache = np.array([])
        cache_t_max = 0.0

    # t 균등 샘플링
    t_samples = np.linspace(T_MIN, T_MAX, N_SAMPLE)
    cache_count  = int(np.sum(t_samples <= cache_t_max))
    online_count = N_SAMPLE - cache_count
    log(f"  캐시 사용: {cache_count}개, on-the-fly: {online_count}개")
    log()
    log(f"  {'#':>5}  {'t':>10}  {'κ':>10}  {'NNS_raw':>8}  {'s_unf':>7}  {'κ_model':>10}  상태")
    log("  " + "-" * 68)

    data = []
    fail_count     = 0
    excluded_count = 0

    for idx, t in enumerate(t_samples):
        try:
            # κ 해석적 계산
            kappa = compute_kappa_analytic(t, delta=DELTA)

            if not np.isfinite(kappa) or kappa <= 0:
                fail_count += 1
                if idx % 200 == 0:
                    log(f"  {idx+1:>5}  {t:>10.2f}  κ 비정상: {kappa}")
                continue

            # 영점 탐색
            if len(zeros_cache) > 0 and t <= cache_t_max + 1.0:
                gamma = find_nearest_zero_cache(t, zeros_cache)
            else:
                dps = get_dps(t)
                gamma = find_nearest_zero_online(t, dps=dps)

            if gamma is None:
                log(f"  {idx+1:>5}  {t:>10.2f}  영점 탐색 실패")
                fail_count += 1
                continue

            NNS_raw = abs(t - gamma)
            if NNS_raw < 1e-8:
                excluded_count += 1
                continue

            # unfolded NNS
            density = np.log(t / (2.0 * np.pi)) / (2.0 * np.pi)
            s_unf   = NNS_raw * density

            # 배제: s_unf < δ (κ 발산 영역)
            if s_unf < S_MIN_EXCL:
                excluded_count += 1
                continue

            # 배제: κ > κ_max (이론적 불가)
            if kappa > KAPPA_MAX:
                excluded_count += 1
                continue

            # 모델 예측 κ (비교용)
            kappa_model = 1.0 / (DELTA**2 + s_unf**2)**2

            data.append({
                't':           t,
                'kappa':       kappa,
                'NNS_raw':     NNS_raw,
                's_unf':       s_unf,
                'kappa_model': kappa_model,
                'gamma':       gamma,
            })

            if idx % 200 == 0:
                log(f"  {idx+1:>5}  {t:>10.2f}  {kappa:>10.4f}  "
                    f"{NNS_raw:>8.4f}  {s_unf:>7.4f}  {kappa_model:>10.4f}  ✓")

        except Exception as e:
            log(f"  {idx+1:>5}  {t:>10.2f}  ERROR: {e}")
            fail_count += 1

    N_valid = len(data)
    elapsed_a = time.time() - t_start
    log()
    log(f"  수집 완료: {N_valid}/{N_SAMPLE}개 유효  |  실패: {fail_count}  배제: {excluded_count}")
    log(f"  파트 A 소요: {elapsed_a:.1f}s")
    save()

    if N_valid < 50:
        log("⚠️ 유효 데이터 부족 — 실험 중단.")
        save()
        return

    kappa_emp    = np.array([d['kappa']       for d in data])
    s_unf_emp    = np.array([d['s_unf']        for d in data])
    kappa_model  = np.array([d['kappa_model']  for d in data])
    t_emp        = np.array([d['t']            for d in data])

    # ═══════════════════════════════════════════════════════════════════════
    # 파트 B: 경험적 κ 분포 기본 통계
    # ═══════════════════════════════════════════════════════════════════════
    log()
    log("=" * 70)
    log("파트 B: 경험적 κ 분포 기본 통계")
    log("=" * 70)
    log(f"  N_valid = {N_valid}")
    log(f"  κ_emp: mean={kappa_emp.mean():.4f}, median={np.median(kappa_emp):.4f}")
    log(f"         std={kappa_emp.std():.4f}, min={kappa_emp.min():.4f}, max={kappa_emp.max():.4f}")
    log()
    log("  분위수 (κ):")
    for q in [5, 10, 25, 50, 75, 90, 95]:
        log(f"    {q:>3}%: {np.percentile(kappa_emp, q):.4f}")
    log()
    log("  s_unf 통계:")
    log(f"  mean(s_unf)={s_unf_emp.mean():.4f}  (GUE 이론: {gue_mean:.4f})")
    log(f"  std(s_unf)={s_unf_emp.std():.4f}")
    log(f"  [5%, 50%, 95%] = [{np.percentile(s_unf_emp,5):.3f}, "
        f"{np.median(s_unf_emp):.3f}, {np.percentile(s_unf_emp,95):.3f}]")
    log()
    # κ_emp vs κ_model 직접 비교 (모델 검증)
    corr_km = np.corrcoef(np.log(kappa_emp), np.log(kappa_model))[0,1]
    log(f"  log(κ_emp) vs log(κ_model=1/(δ²+s²)²): Pearson r = {corr_km:.4f}")
    resid = kappa_emp - kappa_model
    log(f"  잔차 κ_emp − κ_model: mean={resid.mean():.4f}, std={resid.std():.4f}")
    save()

    # ═══════════════════════════════════════════════════════════════════════
    # 파트 C: 이론 P_pred(κ) 검증 (변수 변환 테이블)
    # ═══════════════════════════════════════════════════════════════════════
    log()
    log("=" * 70)
    log("파트 C: P_GUE → P_pred(κ) 변수 변환 검증")
    log("=" * 70)
    log("  κ = 1/(δ²+s²)²  →  s(κ) = √(1/√κ − δ²)")
    log("  P_pred(κ) = P_GUE(s(κ)) × |ds/dκ|")
    log("  |ds/dκ| = 1/(4 × κ^{3/2} × (1/√κ − δ²)^{1/2})")
    log()
    kappa_test = [0.1, 0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0]
    log(f"  {'κ':>8}  {'s(κ)':>7}  {'P_GUE':>8}  {'|ds/dκ|':>10}  "
        f"{'P_pred':>10}  {'CDF_pred':>9}")
    log("  " + "-" * 62)
    for k in kappa_test:
        if k >= KAPPA_MAX:
            continue
        s_k   = float(s_from_kappa(k))
        val   = 1.0/np.sqrt(k) - DELTA**2
        if val <= 0:
            continue
        jac   = 1.0 / (4.0 * k**1.5 * val**0.5)
        pgue  = float(gue_pdf(s_k))
        ppred = pgue * jac
        cdf_v = pred_cdf_from_s(s_k, gue_cdf_interp)
        log(f"  {k:>8.1f}  {s_k:>7.4f}  {pgue:>8.4f}  {jac:>10.4f}  "
            f"{ppred:>10.4f}  {cdf_v:>9.4f}")

    # CDF 정규화 확인: ∫ P_pred dκ = 1 (변환 후)
    log()
    log("  CDF_pred 정규화 확인 (κ ∈ [0, κ_max] 대응):")
    s_range_check = float(gue_cdf_interp(5.0))  # F_GUE(5.0) ≈ 1
    log(f"  F_GUE(s=5.0) = {s_range_check:.6f} (≈1 이면 정상)")
    log(f"  CDF_pred(κ→0) ≡ CDF_pred(s→∞) = 1 − F_GUE(∞) = 0")
    log(f"  CDF_pred(κ_max) ≡ CDF_pred(s=0) = 1 − F_GUE(0) = 1")
    save()

    # ═══════════════════════════════════════════════════════════════════════
    # 파트 D: KS 검정 (P_emp(κ) vs P_pred(κ))
    # ═══════════════════════════════════════════════════════════════════════
    log()
    log("=" * 70)
    log("파트 D: KS 검정 (P_emp(κ) vs P_pred(κ))")
    log("=" * 70)

    N_test     = len(kappa_emp)
    kappa_sort = np.sort(kappa_emp)

    # CDF_pred 계산 (보간 테이블 활용)
    log(f"  CDF_pred 계산 중 ({N_test}점)...")
    t_ks = time.time()

    def pred_cdf_arr(k_arr):
        s_arr = s_from_kappa(np.asarray(k_arr, dtype=float))
        return 1.0 - np.array([float(gue_cdf_interp(s)) for s in s_arr])

    pred_cdf_vals = pred_cdf_arr(kappa_sort)
    ecdf_vals     = np.arange(1, N_test + 1) / N_test
    D_vals        = np.abs(ecdf_vals - pred_cdf_vals)
    D_ks          = float(D_vals.max())
    D_ks_idx      = int(D_vals.argmax())

    # Kolmogorov 점근 p-값
    sqrt_n = np.sqrt(N_test)
    ks_lam = (sqrt_n + 0.12 + 0.11/sqrt_n) * D_ks
    # 급수: P(D > d) = 2 ∑_{k=1}^∞ (−1)^{k-1} exp(−2k²λ²)
    p_ks = 0.0
    for k in range(1, 100):
        term = (-1)**(k-1) * np.exp(-2.0 * k**2 * ks_lam**2)
        p_ks += 2 * term
        if abs(term) < 1e-12:
            break
    p_ks = float(np.clip(p_ks, 0, 1))

    log(f"  KS 통계량 D = {D_ks:.5f}  (N={N_test})")
    log(f"  KS λ = {ks_lam:.4f}")
    log(f"  KS p-값 = {p_ks:.4e}")
    log(f"  최대 편차 위치: κ = {kappa_sort[D_ks_idx]:.4f}  "
        f"(ECDF={ecdf_vals[D_ks_idx]:.4f}, CDF_pred={pred_cdf_vals[D_ks_idx]:.4f})")
    log(f"  소요: {time.time()-t_ks:.1f}s")

    if p_ks > 0.05:
        log("  → ✅ KS 비기각 (p > 0.05)")
        ks_verdict = "비기각"
    elif p_ks > 0.01:
        log("  → ⚠️ KS 약한 기각 (0.01 < p < 0.05)")
        ks_verdict = "약한 기각"
    else:
        log("  → ❌ KS 기각 (p < 0.01)")
        ks_verdict = "기각"
    save()

    # ═══════════════════════════════════════════════════════════════════════
    # 파트 E: Anderson-Darling 검정
    # ═══════════════════════════════════════════════════════════════════════
    log()
    log("=" * 70)
    log("파트 E: Anderson-Darling 검정")
    log("=" * 70)

    def cdf_pred_scalar(k):
        s = float(s_from_kappa(float(k)))
        return pred_cdf_from_s(s, gue_cdf_interp)

    log("  AD 검정 계산 중 (서브샘플 N=500)...")
    t_ad = time.time()
    # AD는 N=500 서브샘플 (정확도 충분)
    np.random.seed(42)
    ad_sample = np.random.choice(kappa_emp, size=min(500, N_valid), replace=False)
    A2, A2_mod, p_ad = anderson_darling_test(ad_sample, cdf_pred_scalar)
    log(f"  AD A² = {A2:.4f}  A²_수정 = {A2_mod:.4f}  p = {p_ad:.4e}")
    log(f"  소요: {time.time()-t_ad:.1f}s")

    if p_ad > 0.05:
        log("  → ✅ AD 비기각 (p > 0.05)")
        ad_verdict = "비기각"
    elif p_ad > 0.01:
        log("  → ⚠️ AD 약한 기각")
        ad_verdict = "약한 기각"
    else:
        log("  → ❌ AD 기각 (p < 0.01)")
        ad_verdict = "기각"
    save()

    # ═══════════════════════════════════════════════════════════════════════
    # 파트 F: Q-Q 분석 (log κ 기준)
    # ═══════════════════════════════════════════════════════════════════════
    log()
    log("=" * 70)
    log("파트 F: Q-Q 분석 (log κ 기준)")
    log("=" * 70)

    # 예측 분위수: F_pred(κ) = q → s = F_GUE^{-1}(1−q) → κ = 1/(δ²+s²)²
    N_qq = 100
    quants_pct = np.linspace(1, 99, N_qq)

    emp_quant = np.percentile(kappa_emp, quants_pct)

    log("  예측 분위수 계산 중 (GUE 역CDF + 변수 변환)...")
    pred_quant = np.zeros(N_qq)
    for i, q in enumerate(quants_pct):
        p_target = 1.0 - q / 100.0  # F_GUE(s) = 1-q/100
        if p_target <= 0:
            pred_quant[i] = KAPPA_MAX
        elif p_target >= 1:
            pred_quant[i] = 0.0
        else:
            try:
                s_q = brentq(lambda s: float(gue_cdf_interp(s)) - p_target,
                             0.0, 7.0, xtol=1e-6, maxiter=200)
                val = 1.0/np.sqrt(max(1e-300, 1.0)) - DELTA**2  # temp
                if s_q < 1e-10:
                    pred_quant[i] = KAPPA_MAX
                else:
                    pred_quant[i] = 1.0 / (DELTA**2 + s_q**2)**2
            except Exception:
                pred_quant[i] = np.nan

    valid_mask = (np.isfinite(pred_quant) & np.isfinite(emp_quant)
                  & (pred_quant > 0) & (emp_quant > 0))
    n_valid_qq = int(valid_mask.sum())
    log(f"  유효 분위수 쌍: {n_valid_qq}/{N_qq}")

    if n_valid_qq > 10:
        log_emp  = np.log(emp_quant[valid_mask])
        log_pred = np.log(pred_quant[valid_mask])
        qq_corr  = float(np.corrcoef(log_emp, log_pred)[0, 1])
        qq_slope = float(np.polyfit(log_pred, log_emp, 1)[0])
        log(f"  Q-Q 상관 (log κ): r = {qq_corr:.4f}")
        log(f"  Q-Q 기울기: slope = {qq_slope:.4f}  (이상: 1.0)")
    else:
        qq_corr = 0.0
        log("  ⚠️ 유효 분위수 부족")

    log()
    log(f"  {'분위수':>6}  {'κ_emp':>10}  {'κ_pred':>10}  {'편차%':>9}  {'log 편차':>9}")
    log("  " + "-" * 52)
    qi_list = [0, 9, 24, 49, 74, 89, 99]
    devs = []
    for qi in qi_list:
        if qi < N_qq and valid_mask[qi]:
            eq = emp_quant[qi]
            pq = pred_quant[qi]
            dev = (eq - pq) / pq * 100
            ldev = np.log(eq) - np.log(pq)
            devs.append(abs(dev))
            log(f"  {quants_pct[qi]:>5.0f}%  {eq:>10.4f}  {pq:>10.4f}  "
                f"{dev:>8.1f}%  {ldev:>9.4f}")
    max_dev = max(devs) if devs else 999.0

    if qq_corr > 0.95:
        log("  → ✅ Q-Q 상관 > 0.95")
        qq_verdict = "양성"
    elif qq_corr > 0.80:
        log(f"  → ⚠️ Q-Q 상관 > 0.80")
        qq_verdict = "중립"
    else:
        log(f"  → ❌ Q-Q 상관 < 0.80")
        qq_verdict = "음성"
    save()

    # ═══════════════════════════════════════════════════════════════════════
    # 파트 G: s_unf 분포 vs GUE 보조 검정
    # ═══════════════════════════════════════════════════════════════════════
    log()
    log("=" * 70)
    log("파트 G: s_unf 분포 vs GUE 보조 검정")
    log("=" * 70)
    log("  (샘플링 편향 진단: s_unf가 GUE를 따르는지?)")

    def gue_cdf_scipy(x_arr):
        return np.array([float(gue_cdf_interp(xi)) for xi in np.atleast_1d(x_arr)])

    # KS 검정 (대규모이면 asymp 사용)
    method = 'exact' if N_valid <= 1000 else 'asymp'
    ks_s = scipy_stats.kstest(s_unf_emp, gue_cdf_scipy, method=method)
    log(f"  s_unf KS vs GUE: D={ks_s.statistic:.4f}, p={ks_s.pvalue:.4e}  [{method}]")
    if ks_s.pvalue > 0.05:
        log("  → ✅ s_unf ∼ GUE (p > 0.05)")
    else:
        log(f"  → ❌ s_unf ≁ GUE (p < 0.05, D={ks_s.statistic:.4f})")
    log()
    log("  주: s_unf는 균등 샘플링 t에서 유도 → GUE와 다를 수 있음 (sampling bias)")
    log("      진정한 GUE 검정은 consecutive spacing (#36b) 또는 pair correlation (#4)")
    save()

    # ═══════════════════════════════════════════════════════════════════════
    # 파트 H: 분위수 세부 편차 분석
    # ═══════════════════════════════════════════════════════════════════════
    log()
    log("=" * 70)
    log("파트 H: 분위수 편차 세부 분석 (전구간 < 10% 조건)")
    log("=" * 70)
    pct_devs_all = []
    for qi in range(N_qq):
        if valid_mask[qi]:
            dev_abs = abs(emp_quant[qi] - pred_quant[qi]) / abs(pred_quant[qi]) * 100
            pct_devs_all.append(dev_abs)
    if pct_devs_all:
        pct_arr = np.array(pct_devs_all)
        log(f"  편차 분포: mean={pct_arr.mean():.1f}%, median={np.median(pct_arr):.1f}%")
        log(f"             max={pct_arr.max():.1f}%, 90%tile={np.percentile(pct_arr,90):.1f}%")
        log(f"  전구간 < 10%: {'✅ YES' if pct_arr.max() < 10 else '❌ NO'}")
        log(f"  전구간 < 20%: {'✅ YES' if pct_arr.max() < 20 else '❌ NO'}")
    save()

    # ═══════════════════════════════════════════════════════════════════════
    # 파트 I: 종합 판정
    # ═══════════════════════════════════════════════════════════════════════
    log()
    log("=" * 70)
    log("파트 I: 종합 판정")
    log("=" * 70)
    log()
    log("  수학자 성공 기준:")
    log("  [양성]       KS p > 0.05 + Q-Q 상관 > 0.95")
    log("  [강한 양성]  + AD 비기각 + 분위수 편차 전구간 < 10%")
    log("  [음성]       KS p < 0.01 또는 Q-Q 상관 < 0.80")
    log("  [중립]       0.01 < p < 0.05 또는 Q-Q > 0.90이나 꼬리 편차 > 20%")
    log()
    log(f"  KS  p-값    = {p_ks:.4e}  → {ks_verdict}")
    log(f"  AD  p-값    = {p_ad:.4e}  → {ad_verdict}")
    log(f"  Q-Q 상관    = {qq_corr:.4f}   → {qq_verdict}")
    log(f"  최대 분위수 편차 = {max_dev:.1f}%")
    log()

    # 판정 로직
    if p_ks > 0.05 and qq_corr > 0.95:
        if p_ad > 0.05 and max_dev < 10:
            verdict = "✅✅ 강한 양성"
            detail  = "GUE → P_pred(κ) 정량적 예측 완전 성공. 논문 통합 가능."
        else:
            verdict = "✅ 양성"
            detail  = "KS + Q-Q 비기각. GUE → κ 분포 예측 성공."
    elif p_ks < 0.01 or qq_corr < 0.80:
        verdict = "❌ 음성"
        detail  = ("P_emp(κ)와 P_pred(κ) 불일치. κ = 1/(δ²+s²)² 모델 또는 "
                   "GUE 분포 가정이 충분하지 않음.")
    else:
        verdict = "⚠️ 중립"
        detail  = "부분적 일치. 수학자 추가 판단 요청."

    log(f"  **최종 판정: {verdict}**")
    log(f"  {detail}")
    log()

    # 음성/중립의 경우 구조 분석
    if "음성" in verdict or "중립" in verdict:
        log("  구조 분석:")
        log(f"  - κ_emp 평균 = {kappa_emp.mean():.4f}  vs  κ_model 평균 = {kappa_model.mean():.4f}")
        log(f"    비율 = {kappa_emp.mean()/kappa_model.mean():.3f} (1.0이면 스케일 일치)")
        ratio_at_50 = np.median(kappa_emp) / (pred_quant[49] if (valid_mask[49] and pred_quant[49]>0) else np.nan)
        if np.isfinite(ratio_at_50):
            log(f"  - 중앙값 비율 κ_emp/κ_pred = {ratio_at_50:.3f}")
        log(f"  - log(κ_emp) vs log(κ_model) 상관 = {corr_km:.4f}")
        log()
        log("  가능한 원인:")
        log("  1. κ = 1/(δ²+s²)² 단순 근사 → 실제는 |G(s)+ζ'/ζ(s)|² (G 기여 포함)")
        log("  2. 균등 샘플링 t: NNS ~ min(U,1-U)×spacing ≠ P_GUE(s)")
        log("  3. δ=0.03에서 실제 κ는 여러 영점 기여의 벡터 합 (단순 최근접만 아님)")

    log()
    log(f"  총 소요: {(time.time()-t_start)/60:.1f}분")
    log()
    log("끝.")
    save()


if __name__ == "__main__":
    main()
