#!/usr/bin/env python3
"""
[C-296] B-57 해소: 고높이 ζ(s) W-의존 A-gap 실험 (v2)

v1 실패: mpmath.zetazero(n>4000)에서 20시간+ 멈춤. N=4000으로 축소.
v2 변경: N=4000, 점진적 캐시 저장, W 최대 1500.

목표: δ(W→∞) = 0 인가? 잔여 상수인가?

방법:
  Phase 1: ζ(s) 영점 캐시 N=4000 확장 (mpmath.zetazero, DPS=30)
  Phase 2: W-의존 ρ_S 측정 + bootstrap 95% CI (1000 리샘플)
    - ζ(s): 캐시에서 로드, 이론적 정규화, trim 20%
    - GUE: N=2000, 30앙상블, 반원 unfolding, trim 20%
    - W = [5, 10, 20, 40, 80, 200, 500, 750, 1000, 1500]
  Phase 3: 외삽 분석
    - (a) δ = a / W^b           → δ→0
    - (b) δ = c + a / W^b       → δ→c (잔여 감쇠?)
    - (c) δ = a · exp(-b·W)     → δ→0
    - AIC/BIC 모델 선택 + c의 95% CI

체크리스트:
  [x] dps=30 for zetazero
  [x] 이론적 d_bar = log(t/(2π))/(2π) for ζ(s)
  [x] GUE: 반원 CDF unfolding → 간격=1
  [x] trim 20%, edge skip
  [x] python -u
  [x] A_bare = S1^2 + 2*H1 (local ±W)
  [x] Spearman only (B-47)
  [x] NaN/Inf 체크
  [x] bootstrap CI 1000 리샘플
  [x] except: pass 금지
  [x] 점진적 캐시 저장 (100개마다)
"""

import sys, os, math, time
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

# ── 설정 ──────────────────────────────────────────────────────────
CACHE_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/outputs/cache/zeta_zeros_t_N4000_c296.npy')
OLD_CACHE_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/outputs/cache/zeta_zeros_t_N5000.npy')
N_ZEROS_TARGET = 4000

TRIM_FRAC = 0.20
W_LIST = [5, 10, 20, 40, 80, 200, 500, 750, 1000, 1500]
GUE_N = 2000
GUE_ENS = 30
GUE_SEED = 42
BOOTSTRAP_N = 1000
BOOTSTRAP_SEED = 12345

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/high_T_w_dependence_c296.txt')

os.makedirs(os.path.dirname(CACHE_PATH), exist_ok=True)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

out_f = open(RESULT_PATH, 'w')
def log(msg=''):
    print(msg, flush=True)
    out_f.write(msg + '\n')
    out_f.flush()


# ══════════════════════════════════════════════════════════════════
# Phase 1: 영점 캐시 확장 (점진적 저장)
# ══════════════════════════════════════════════════════════════════

def build_or_load_cache():
    """N=4000 ζ(s) 영점 캐시. 점진적 저장."""
    if os.path.exists(CACHE_PATH):
        arr = np.load(CACHE_PATH)
        if len(arr) >= N_ZEROS_TARGET:
            log(f"[Phase 1] 캐시 로드: {len(arr)}개 영점, t∈[{arr[0]:.2f}, {arr[-1]:.2f}]")
            return arr
        else:
            log(f"[Phase 1] 기존 캐시 {len(arr)}개 → {N_ZEROS_TARGET}개로 확장")
            existing = list(arr)
    elif os.path.exists(OLD_CACHE_PATH):
        arr = np.load(OLD_CACHE_PATH)
        log(f"[Phase 1] 이전 캐시에서 {len(arr)}개 로드, {N_ZEROS_TARGET}개로 확장")
        existing = list(arr)
    else:
        log(f"[Phase 1] 캐시 없음 → 처음부터 {N_ZEROS_TARGET}개 계산")
        existing = []

    import mpmath
    mpmath.mp.dps = 30

    start_n = len(existing) + 1
    if start_n > N_ZEROS_TARGET:
        return np.array(existing[:N_ZEROS_TARGET])

    t_start = time.time()
    log(f"  mpmath.zetazero({start_n}) ~ zetazero({N_ZEROS_TARGET}) 계산 중...")

    for n in range(start_n, N_ZEROS_TARGET + 1):
        rho = mpmath.zetazero(n)
        existing.append(float(rho.imag))

        # 100개마다 진행 상황 + 점진적 저장
        if n % 100 == 0:
            elapsed = time.time() - t_start
            rate = (n - start_n + 1) / elapsed if elapsed > 0 else 0
            log(f"  {n}/{N_ZEROS_TARGET} ({elapsed:.0f}s, {rate:.1f} zeros/s)")
            # 점진적 저장
            np.save(CACHE_PATH, np.array(existing))

    arr = np.array(existing)
    np.save(CACHE_PATH, arr)
    elapsed = time.time() - t_start
    log(f"  저장 완료: {len(arr)}개 영점, t∈[{arr[0]:.2f}, {arr[-1]:.2f}] ({elapsed:.0f}s)")
    return arr


# ══════════════════════════════════════════════════════════════════
# Phase 2: W-의존 ρ_S 측정 + Bootstrap CI
# ══════════════════════════════════════════════════════════════════

def generate_gue_eigvals(N, seed=None):
    rng = np.random.default_rng(seed)
    real = rng.standard_normal((N, N))
    imag = rng.standard_normal((N, N))
    A = (real + 1j * imag) / np.sqrt(2.0)
    H = (A + A.conj().T) / (2.0 * np.sqrt(N))
    return np.linalg.eigvalsh(H).real


def semicircle_cdf(x):
    x = np.clip(x, -2.0 + 1e-12, 2.0 - 1e-12)
    return (np.arcsin(x / 2.0) + (x / 2.0) * np.sqrt(
        np.maximum(0.0, 1.0 - x**2 / 4.0))) / np.pi + 0.5


def unfold(lambdas):
    return len(lambdas) * semicircle_cdf(lambdas)


def compute_A_bare_vectorized(x, W):
    """벡터화된 A_bare 계산. x: sorted array."""
    n = len(x)
    if n < 2 * W + 1:
        return np.array([]), np.arange(0)

    valid_lo = W
    valid_hi = n - W
    n_valid = valid_hi - valid_lo
    if n_valid < 1:
        return np.array([]), np.arange(0)

    S1 = np.zeros(n_valid)
    H1 = np.zeros(n_valid)

    for j in range(1, W + 1):
        diff_r = x[valid_lo + j : valid_hi + j] - x[valid_lo : valid_hi]
        diff_l = x[valid_lo : valid_hi] - x[valid_lo - j : valid_hi - j]
        diff_r = np.where(np.abs(diff_r) < 1e-15, 1e-15, diff_r)
        diff_l = np.where(np.abs(diff_l) < 1e-15, 1e-15, diff_l)
        S1 += -1.0 / diff_r + 1.0 / diff_l
        H1 += 1.0 / diff_r**2 + 1.0 / diff_l**2

    A = S1**2 + 2.0 * H1
    indices = np.arange(valid_lo, valid_hi)
    return A, indices


def zeta_density(t):
    """ζ(s) 이론적 영점 밀도: d̄(t) = log(t/(2π))/(2π)"""
    if t <= 2 * math.pi:
        return 0.1
    return math.log(t / (2 * math.pi)) / (2 * math.pi)


def get_A_gap_data(x, W, density_func=None):
    """A_bare와 gap_min 데이터 반환 (trim 적용 후)"""
    A, indices = compute_A_bare_vectorized(x, W)
    if len(A) == 0:
        return np.array([]), np.array([])

    gaps_l = x[indices] - x[indices - 1]
    gaps_r = x[indices + 1] - x[indices]
    gap_min = np.minimum(gaps_l, gaps_r)

    if density_func is not None:
        d_bar = np.array([density_func(xi) for xi in x[indices]])
        gap_min = gap_min * d_bar

    # trim 20%
    n_pts = len(A)
    trim_lo = int(n_pts * TRIM_FRAC)
    trim_hi = n_pts - int(n_pts * TRIM_FRAC)
    if trim_hi - trim_lo < 10:
        return np.array([]), np.array([])

    A_t = A[trim_lo:trim_hi]
    gm_t = gap_min[trim_lo:trim_hi]

    mask = np.isfinite(A_t) & np.isfinite(gm_t) & (A_t > 0) & (gm_t > 0)
    return A_t[mask], gm_t[mask]


def bootstrap_spearman(A, gm, n_boot=BOOTSTRAP_N, seed=BOOTSTRAP_SEED):
    """Bootstrap 95% CI for Spearman ρ"""
    rng = np.random.default_rng(seed)
    n = len(A)
    if n < 10:
        return np.nan, np.nan, np.nan, np.nan

    rho_obs, p_obs = stats.spearmanr(A, gm)

    boot_rhos = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        r, _ = stats.spearmanr(A[idx], gm[idx])
        boot_rhos[b] = r

    boot_rhos = boot_rhos[np.isfinite(boot_rhos)]
    if len(boot_rhos) < 100:
        return float(rho_obs), float(p_obs), np.nan, np.nan

    ci_lo = np.percentile(boot_rhos, 2.5)
    ci_hi = np.percentile(boot_rhos, 97.5)
    return float(rho_obs), float(p_obs), float(ci_lo), float(ci_hi)


# ══════════════════════════════════════════════════════════════════
# Phase 3: 외삽 분석
# ══════════════════════════════════════════════════════════════════

def model_a(W, a, b):
    """δ = a / W^b → δ→0"""
    return a / np.power(W, b)

def model_b(W, c, a, b):
    """δ = c + a / W^b → δ→c"""
    return c + a / np.power(W, b)

def model_c(W, a, b):
    """δ = a · exp(-b·W) → δ→0"""
    return a * np.exp(-b * W)

def compute_aic_bic(residuals, k, n):
    ss = np.sum(residuals**2)
    if ss <= 0 or n <= k:
        return np.inf, np.inf
    ll = -n/2 * np.log(2 * np.pi * ss / n) - n/2
    aic = 2 * k - 2 * ll
    bic = k * np.log(n) - 2 * ll
    return aic, bic


def extrapolation_analysis(W_arr, delta_arr, delta_ci_lo, delta_ci_hi):
    """3 모델 적합 + AIC/BIC 비교"""
    log(f"\n{'='*70}")
    log("  [Phase 3] 외삽 분석: δ(W→∞)")
    log(f"{'='*70}\n")

    W = np.array(W_arr, dtype=float)
    d = np.array(delta_arr, dtype=float)
    n = len(W)

    se = (np.array(delta_ci_hi) - np.array(delta_ci_lo)) / (2 * 1.96)
    se = np.where(se > 0, se, 0.005)

    results = {}

    # Model (a): δ = a / W^b
    try:
        popt_a, pcov_a = curve_fit(model_a, W, d, p0=[0.1, 0.5],
                                    sigma=se, absolute_sigma=True,
                                    maxfev=10000)
        pred_a = model_a(W, *popt_a)
        resid_a = d - pred_a
        aic_a, bic_a = compute_aic_bic(resid_a, 2, n)
        se_a = np.sqrt(np.diag(pcov_a)) if pcov_a is not None else [np.nan, np.nan]
        results['a'] = {
            'params': {'a': popt_a[0], 'b': popt_a[1]},
            'se': {'a': se_a[0], 'b': se_a[1]},
            'aic': aic_a, 'bic': bic_a,
            'lim': 0.0,
            'label': 'δ = a/W^b → 0'
        }
        log(f"  모델 (a) δ = a/W^b:")
        log(f"    a = {popt_a[0]:.6f} ± {se_a[0]:.6f}")
        log(f"    b = {popt_a[1]:.6f} ± {se_a[1]:.6f}")
        log(f"    AIC = {aic_a:.2f}, BIC = {bic_a:.2f}")
        log(f"    δ(W→∞) = 0")
    except Exception as e:
        log(f"  모델 (a) 적합 실패: {e}")
        results['a'] = None

    # Model (b): δ = c + a / W^b
    try:
        popt_b, pcov_b = curve_fit(model_b, W, d, p0=[0.01, 0.1, 0.5],
                                    sigma=se, absolute_sigma=True,
                                    maxfev=10000)
        pred_b = model_b(W, *popt_b)
        resid_b = d - pred_b
        aic_b, bic_b = compute_aic_bic(resid_b, 3, n)
        se_b = np.sqrt(np.diag(pcov_b)) if pcov_b is not None else [np.nan]*3

        c_val = popt_b[0]
        c_se = se_b[0]
        c_ci_lo = c_val - 1.96 * c_se
        c_ci_hi = c_val + 1.96 * c_se
        c_includes_zero = (c_ci_lo <= 0 <= c_ci_hi)

        results['b'] = {
            'params': {'c': popt_b[0], 'a': popt_b[1], 'b': popt_b[2]},
            'se': {'c': se_b[0], 'a': se_b[1], 'b': se_b[2]},
            'aic': aic_b, 'bic': bic_b,
            'lim': c_val,
            'c_ci': (c_ci_lo, c_ci_hi),
            'c_includes_zero': c_includes_zero,
            'label': f'δ = c + a/W^b → {c_val:.6f}'
        }
        log(f"\n  모델 (b) δ = c + a/W^b:")
        log(f"    c = {popt_b[0]:.6f} ± {se_b[0]:.6f}  [95%CI: ({c_ci_lo:.6f}, {c_ci_hi:.6f})]")
        log(f"    a = {popt_b[1]:.6f} ± {se_b[1]:.6f}")
        log(f"    b = {popt_b[2]:.6f} ± {se_b[2]:.6f}")
        log(f"    AIC = {aic_b:.2f}, BIC = {bic_b:.2f}")
        log(f"    δ(W→∞) = c = {c_val:.6f}")
        log(f"    c의 95%CI가 0을 포함? {'예 → δ→0 지지' if c_includes_zero else '아니오 → 잔여 감쇠 존재'}")
    except Exception as e:
        log(f"  모델 (b) 적합 실패: {e}")
        results['b'] = None

    # Model (c): δ = a · exp(-b·W)
    try:
        popt_c, pcov_c = curve_fit(model_c, W, d, p0=[0.05, 0.001],
                                    sigma=se, absolute_sigma=True,
                                    maxfev=10000)
        pred_c = model_c(W, *popt_c)
        resid_c = d - pred_c
        aic_c, bic_c = compute_aic_bic(resid_c, 2, n)
        se_c = np.sqrt(np.diag(pcov_c)) if pcov_c is not None else [np.nan, np.nan]
        results['c'] = {
            'params': {'a': popt_c[0], 'b': popt_c[1]},
            'se': {'a': se_c[0], 'b': se_c[1]},
            'aic': aic_c, 'bic': bic_c,
            'lim': 0.0,
            'label': 'δ = a·exp(-bW) → 0'
        }
        log(f"\n  모델 (c) δ = a·exp(-bW):")
        log(f"    a = {popt_c[0]:.6f} ± {se_c[0]:.6f}")
        log(f"    b = {popt_c[1]:.6f} ± {se_c[1]:.6f}")
        log(f"    AIC = {aic_c:.2f}, BIC = {bic_c:.2f}")
        log(f"    δ(W→∞) = 0")
    except Exception as e:
        log(f"  모델 (c) 적합 실패: {e}")
        results['c'] = None

    # AIC/BIC 비교
    log(f"\n  {'─'*60}")
    log(f"  AIC/BIC 모델 비교")
    log(f"  {'─'*60}")
    log(f"  {'모델':<25} {'AIC':>10} {'BIC':>10} {'δ(∞)':>10}")
    log(f"  {'─'*60}")

    for key in ['a', 'b', 'c']:
        if results[key] is not None:
            r = results[key]
            lim_str = f"{r['lim']:.6f}" if r['lim'] != 0.0 else "0"
            log(f"  {r['label']:<25} {r['aic']:>10.2f} {r['bic']:>10.2f} {lim_str:>10}")

    valid = {k: v for k, v in results.items() if v is not None}
    if valid:
        best_aic = min(valid.keys(), key=lambda k: valid[k]['aic'])
        best_bic = min(valid.keys(), key=lambda k: valid[k]['bic'])
        log(f"\n  최적 모델 (AIC): ({best_aic}) {valid[best_aic]['label']}")
        log(f"  최적 모델 (BIC): ({best_bic}) {valid[best_bic]['label']}")

        min_aic = valid[best_aic]['aic']
        for k, v in valid.items():
            delta_aic = v['aic'] - min_aic
            log(f"    ΔAIC({k}) = {delta_aic:.2f}")

    return results


# ══════════════════════════════════════════════════════════════════
# 메인
# ══════════════════════════════════════════════════════════════════

def main():
    t_total = time.time()
    log("=" * 70)
    log("[C-296] B-57 해소: 고높이 ζ(s) W-의존 A-gap 실험 (v2)")
    log(f"  v2: N=4000 (v1 N=5000 mpmath 감속 → 축소)")
    log(f"  Phase 1: 영점 캐시 N={N_ZEROS_TARGET}")
    log(f"  Phase 2: W={W_LIST}, bootstrap {BOOTSTRAP_N}회")
    log(f"  Phase 3: 3-모델 외삽 + AIC/BIC")
    log("=" * 70)

    # ── Phase 1: 영점 캐시 ──────────────────────────────────────
    zeros_all = build_or_load_cache()
    n_total = len(zeros_all)
    log(f"\n  총 {n_total}개 영점, t∈[{zeros_all[0]:.2f}, {zeros_all[-1]:.2f}]")

    # ── Phase 2: W별 ρ 측정 ─────────────────────────────────────
    log(f"\n{'='*70}")
    log("  [Phase 2] W별 ρ_S 측정 + Bootstrap 95% CI")
    log(f"{'='*70}")

    # GUE 앙상블 생성
    log(f"\n[GUE] {GUE_ENS}앙상블 × N={GUE_N} 생성 + unfolding ...")
    t0 = time.time()
    gue_unfolded_list = []
    for ens in range(GUE_ENS):
        seed = GUE_SEED + ens * 7919 + GUE_N
        lambdas = generate_gue_eigvals(GUE_N, seed=seed)
        xi = unfold(lambdas)
        gue_unfolded_list.append(xi)
    log(f"  완료: {time.time()-t0:.1f}s")

    all_results = []

    for W in W_LIST:
        t0 = time.time()
        log(f"\n--- W = {W} ---")

        n_usable = n_total - 2 * W
        if n_usable < 30:
            log(f"  ⚠️ 유효 영점 부족: n_usable={n_usable} < 30. 스킵.")
            continue

        # ζ(s)
        A_z, gm_z = get_A_gap_data(zeros_all, W, density_func=zeta_density)
        if len(A_z) < 10:
            log(f"  ⚠️ ζ(s) 데이터 부족: n={len(A_z)}. 스킵.")
            continue

        rho_z, p_z, ci_lo_z, ci_hi_z = bootstrap_spearman(A_z, gm_z)
        dt_z = time.time() - t0
        log(f"  ζ(s): ρ={rho_z:.6f} [{ci_lo_z:.6f}, {ci_hi_z:.6f}], n={len(A_z)}  ({dt_z:.1f}s)")

        # GUE (앙상블 합산)
        t1 = time.time()
        all_A_g = []
        all_gm_g = []
        for xi in gue_unfolded_list:
            A_g, idx_g = compute_A_bare_vectorized(xi, W)
            if len(A_g) == 0:
                continue
            gaps_l = xi[idx_g] - xi[idx_g - 1]
            gaps_r = xi[idx_g + 1] - xi[idx_g]
            gm_g = np.minimum(gaps_l, gaps_r)
            n_pts = len(A_g)
            tl = int(n_pts * TRIM_FRAC)
            th = n_pts - int(n_pts * TRIM_FRAC)
            if th - tl < 5:
                continue
            A_g_t = A_g[tl:th]
            gm_g_t = gm_g[tl:th]
            mask = np.isfinite(A_g_t) & np.isfinite(gm_g_t) & (A_g_t > 0) & (gm_g_t > 0)
            all_A_g.extend(A_g_t[mask].tolist())
            all_gm_g.extend(gm_g_t[mask].tolist())

        if len(all_A_g) >= 10:
            A_gue = np.array(all_A_g)
            gm_gue = np.array(all_gm_g)
            rho_g, p_g, ci_lo_g, ci_hi_g = bootstrap_spearman(A_gue, gm_gue)
            n_g = len(A_gue)
        else:
            rho_g, p_g, ci_lo_g, ci_hi_g, n_g = np.nan, np.nan, np.nan, np.nan, 0

        dt_g = time.time() - t1
        log(f"  GUE:  ρ={rho_g:.6f} [{ci_lo_g:.6f}, {ci_hi_g:.6f}], n={n_g}  ({dt_g:.1f}s)")

        # δ 계산 + bootstrap CI for δ
        if np.isfinite(rho_z) and np.isfinite(rho_g):
            delta = abs(rho_g) - abs(rho_z)

            rng_d = np.random.default_rng(BOOTSTRAP_SEED + W)
            boot_deltas = np.empty(BOOTSTRAP_N)
            for b in range(BOOTSTRAP_N):
                idx_z = rng_d.integers(0, len(A_z), size=len(A_z))
                idx_g = rng_d.integers(0, len(A_gue), size=len(A_gue))
                r_z, _ = stats.spearmanr(A_z[idx_z], gm_z[idx_z])
                r_g, _ = stats.spearmanr(A_gue[idx_g], gm_gue[idx_g])
                boot_deltas[b] = abs(r_g) - abs(r_z)

            boot_deltas = boot_deltas[np.isfinite(boot_deltas)]
            if len(boot_deltas) >= 100:
                delta_ci_lo = float(np.percentile(boot_deltas, 2.5))
                delta_ci_hi = float(np.percentile(boot_deltas, 97.5))
                delta_se = float(np.std(boot_deltas))
            else:
                delta_ci_lo = delta_ci_hi = delta_se = np.nan
        else:
            delta = np.nan
            delta_ci_lo = delta_ci_hi = delta_se = np.nan

        log(f"  δ(W={W}) = {delta:.6f} [{delta_ci_lo:.6f}, {delta_ci_hi:.6f}] (SE={delta_se:.6f})")

        all_results.append({
            'W': W,
            'rho_zeta': rho_z, 'ci_z': (ci_lo_z, ci_hi_z), 'n_zeta': len(A_z),
            'rho_gue': rho_g, 'ci_g': (ci_lo_g, ci_hi_g), 'n_gue': n_g,
            'delta': delta, 'delta_ci': (delta_ci_lo, delta_ci_hi), 'delta_se': delta_se,
        })

    # ── 최종 비교표 ──────────────────────────────────────────────
    log(f"\n{'='*70}")
    log("  [최종 비교표] W-의존 A-gap 상관 (Bootstrap 95% CI)")
    log(f"{'='*70}\n")

    log(f"  {'W':>5}  {'ρ_ζ':>8}  {'ρ_ζ 95%CI':>20}  {'n_ζ':>5}  "
        f"{'ρ_GUE':>8}  {'n_GUE':>6}  {'δ(W)':>8}  {'δ 95%CI':>22}  {'δ SE':>7}")
    log("  " + "-" * 110)

    for r in all_results:
        ci_z_str = f"[{r['ci_z'][0]:.4f}, {r['ci_z'][1]:.4f}]"
        ci_d_str = f"[{r['delta_ci'][0]:.4f}, {r['delta_ci'][1]:.4f}]"
        log(f"  {r['W']:>5}  {r['rho_zeta']:>8.4f}  {ci_z_str:>20}  {r['n_zeta']:>5}  "
            f"{r['rho_gue']:>8.4f}  {r['n_gue']:>6}  {r['delta']:>8.4f}  {ci_d_str:>22}  {r['delta_se']:>7.4f}")

    # ── δ 안정성 분석 ────────────────────────────────────────────
    log(f"\n{'='*70}")
    log("  [δ 안정성 분석]")
    log(f"{'='*70}\n")

    deltas = [r['delta'] for r in all_results if np.isfinite(r['delta'])]
    Ws = [r['W'] for r in all_results if np.isfinite(r['delta'])]

    if len(deltas) >= 3:
        d_mean = np.mean(deltas)
        d_std = np.std(deltas, ddof=1)
        d_cv = d_std / abs(d_mean) * 100 if abs(d_mean) > 1e-10 else float('inf')
        log(f"  δ 전체: 평균={d_mean:.6f}, σ={d_std:.6f}, CV={d_cv:.1f}%")

        rho_trend, p_trend = stats.spearmanr(Ws, deltas)
        log(f"  W-δ 추세: ρ_S={rho_trend:.4f}, p={p_trend:.4f}")

        sub = [(w, d) for w, d in zip(Ws, deltas) if w >= 40]
        if len(sub) >= 3:
            d_sub = [s[1] for s in sub]
            w_sub = [s[0] for s in sub]
            log(f"  [W≥40] δ 평균={np.mean(d_sub):.6f}, σ={np.std(d_sub, ddof=1):.6f}")
            rs, ps = stats.spearmanr(w_sub, d_sub)
            log(f"  [W≥40] 추세: ρ_S={rs:.4f}, p={ps:.4f}")

    # ── Phase 3: 외삽 ────────────────────────────────────────────
    if len(all_results) >= 5:
        W_fit = [r['W'] for r in all_results]
        d_fit = [r['delta'] for r in all_results]
        d_ci_lo = [r['delta_ci'][0] for r in all_results]
        d_ci_hi = [r['delta_ci'][1] for r in all_results]

        ext_results = extrapolation_analysis(W_fit, d_fit, d_ci_lo, d_ci_hi)
    else:
        log("\n  데이터 부족 — 외삽 분석 스킵")
        ext_results = {}

    # ── 종합 판정 ────────────────────────────────────────────────
    log(f"\n{'='*70}")
    log("  [종합 판정] B-57: δ(W→∞) = ?")
    log(f"{'='*70}\n")

    if len(deltas) >= 5:
        if p_trend < 0.05 and rho_trend < 0:
            log(f"  1. δ(W) 감소 추세 확인: ρ_S={rho_trend:.4f}, p={p_trend:.4f}")
        elif p_trend >= 0.05:
            log(f"  1. δ(W) 추세 비유의: ρ_S={rho_trend:.4f}, p={p_trend:.4f}")

        if ext_results and 'b' in ext_results and ext_results['b'] is not None:
            rb = ext_results['b']
            c_val = rb['params']['c']
            c_ci = rb['c_ci']
            if rb['c_includes_zero']:
                log(f"  2. 모델 (b) c={c_val:.6f}, 95%CI=({c_ci[0]:.6f}, {c_ci[1]:.6f})")
                log(f"     → c의 CI가 0을 포함 → δ→0 지지")
            else:
                log(f"  2. 모델 (b) c={c_val:.6f}, 95%CI=({c_ci[0]:.6f}, {c_ci[1]:.6f})")
                log(f"     → c의 CI가 0을 포함하지 않음 → 잔여 감쇠 c≈{c_val:.4f} 존재")

        last = all_results[-1]
        log(f"\n  3. 최대 W={last['W']}에서:")
        log(f"     δ = {last['delta']:.6f}")
        log(f"     95%CI = [{last['delta_ci'][0]:.6f}, {last['delta_ci'][1]:.6f}]")
        if last['delta_ci'][0] <= 0:
            log(f"     → δ의 CI가 0을 포함 → 최대 W에서도 δ=0 가능")
        else:
            log(f"     → δ의 CI가 0 초과 → 최대 W에서 δ>0 유의")

        log(f"\n  ── 최종 ──")
        if ext_results and 'b' in ext_results and ext_results['b'] is not None:
            rb = ext_results['b']
            if rb['c_includes_zero'] and p_trend < 0.05 and rho_trend < 0:
                log(f"  ★★★★★ δ→0 지지: 감소 추세(p={p_trend:.4f}) + c의 CI가 0 포함")
                log(f"  → B-57 해소: 산술 비국소성 시차는 유한 윈도우 아티팩트")
                log(f"  → GUE universality at d=1 under sufficient window 확립")
            elif not rb['c_includes_zero'] and p_trend < 0.05 and rho_trend < 0:
                log(f"  ★★★★ δ→c={rb['params']['c']:.6f}>0: 감소 추세이나 잔여 감쇠 존재")
                log(f"  → B-57 부분 해소: 진성 산술 감쇠 ≈{rb['params']['c']:.4f} 확정")
            elif p_trend >= 0.05:
                log(f"  ★★★ δ 추세 비유의 + 외삽 불확실")
                log(f"  → B-57 미해소: 더 높은 T 필요")
            else:
                log(f"  ★★★ δ 패턴 복합: 추가 분석 필요")
        else:
            log(f"  ★★★ 외삽 분석 불완전 → 데이터 부족")

    # ── C-295 비교 ───────────────────────────────────────────────
    log(f"\n{'='*70}")
    log("  [C-295 비교] (T=2000, N=1517 → C-296: T≈{:.0f}, N={})".format(
        zeros_all[-1], n_total))
    log(f"{'='*70}\n")

    c295_data = [
        (5, 0.039), (10, 0.039), (20, 0.040), (40, 0.037),
        (80, 0.028), (120, 0.026), (200, 0.018), (300, 0.014)
    ]
    log(f"  {'W':>5}  {'δ(C-295)':>10}  {'δ(C-296)':>10}  {'Δ':>8}")
    log(f"  {'─'*40}")
    for w295, d295 in c295_data:
        d296 = None
        for r in all_results:
            if r['W'] == w295:
                d296 = r['delta']
                break
        if d296 is not None:
            log(f"  {w295:>5}  {d295:>10.4f}  {d296:>10.4f}  {d296-d295:>8.4f}")
        else:
            log(f"  {w295:>5}  {d295:>10.4f}  {'N/A':>10}  {'—':>8}")

    for r in all_results:
        if r['W'] not in [w for w, _ in c295_data]:
            log(f"  {r['W']:>5}  {'(신규)':>10}  {r['delta']:>10.4f}  {'—':>8}")

    log(f"\n  총 소요: {time.time()-t_total:.1f}s")
    log("=" * 70)
    out_f.close()
    print(f"\n✅ 결과 저장: {RESULT_PATH}")


if __name__ == "__main__":
    main()
