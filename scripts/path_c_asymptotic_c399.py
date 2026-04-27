#!/usr/bin/env python3
"""
[C-399] Path C 비율의 T→∞ 점근 안정성 검증

목표: R(T) = Var(2/g²) / (|Cross|+|Resid|) 가 T→∞에서 어떤 행동을 보이는지 확인.
  - R(T) → c > 1 : 해석적 증명 가능
  - R(T) → 1 : 수치적으로만 유효
  - R(T) → ∞ : 더 강한 결과 가능

C-395 기반 파라미터 확장:
  - T_MAX: 2000 → 10000
  - T-대역: 5 → 10 등간격
  - R(T) vs T 선형 회귀

체크리스트:
  [x] python -u
  [x] mpmath dps=50
  [x] trim 20%
  [x] NaN/Inf 체크
  [x] except Exception as e: print(...)
  [x] 이론적 밀도 d(t) = log(t/(2π))/(2π)
  [x] W=100, Gamma 보정 포함 (A_Λ)
  [x] g=0 안전 가드 (1e-15)
"""

import sys, os, time, math
import numpy as np
from scipy import stats

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')

try:
    import mpmath
    mpmath.mp.dps = 50
    print("mpmath OK, dps=50", flush=True)
except Exception as e:
    print(f"FATAL mpmath: {e}", flush=True)
    sys.exit(1)

try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(2048 * 10**6)
    pari.set_real_precision(80)
    print("cypari2 OK", flush=True)
except Exception as e:
    print(f"FATAL cypari2: {e}", flush=True)
    sys.exit(1)

# ── 설정 ──────────────────────────────────────────────────────────
T_MAX = 5000.0
TRIM_FRAC = 0.20
W = 100
SIGMA_C = 0.5
GAMMA_V = [0]
N_BOOTSTRAP = 2000
N_BANDS = 10
RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/path_c_asymptotic_c399.txt'
)

out_f = open(RESULT_PATH, 'w')
def log(msg=''):
    print(msg, flush=True)
    out_f.write(msg + '\n')
    out_f.flush()


def zeta_density(t):
    if t < 2.0:
        return 1.0
    return math.log(t / (2.0 * math.pi)) / (2.0 * math.pi)


def gamma_corrections(t0):
    s = mpmath.mpc(SIGMA_C, t0)
    im_sum = mpmath.mpf(0)
    re_sum = mpmath.mpf(0)
    for mu in GAMMA_V:
        arg = (s + mu) / 2
        psi_val = mpmath.digamma(arg)
        psi1_val = mpmath.psi(1, arg)
        im_sum += mpmath.im(psi_val) / 2
        re_sum += mpmath.re(psi1_val) / 4
    return float(im_sum), float(re_sum)


def get_zeta_zeros():
    """구간별 분할 — T=5000 단일 호출은 30분 이상 걸림."""
    log(f"[ζ(s)] 영점 수집 t∈(0, {T_MAX}] (1000단위 구간별 분할) ...")
    t0 = time.time()
    pari('L_zeta = lfuncreate(1)')

    CHUNK = 1000
    all_zeros = []
    bounds = list(range(0, int(T_MAX), CHUNK)) + [int(T_MAX)]

    for k in range(len(bounds) - 1):
        lo, hi = bounds[k], bounds[k + 1]
        t_chunk = time.time()
        pari(f'Li_ch = lfuninit(L_zeta, [{max(0,lo-10)}, {hi + 10}])')
        pari(f'zv_ch = lfunzeros(Li_ch, [{lo}, {hi}])')
        nc = int(str(pari('#zv_ch')))
        for i in range(1, nc + 1):
            s = str(pari(f'zv_ch[{i}]')).strip().replace(' E', 'e').replace('E ', 'e')
            try:
                t = float(s)
                if t > 0.5:
                    all_zeros.append(t)
            except ValueError:
                pass
        log(f"  구간 [{lo}, {hi}]: {nc}개, 누적 {len(all_zeros)}개 ({time.time()-t_chunk:.0f}s)")

    zeros = np.array(sorted(set(all_zeros)))
    log(f"  총 {len(zeros)}개 영점, t∈[{zeros[0]:.3f}, {zeros[-1]:.3f}], {time.time()-t0:.1f}s")
    return zeros


def precompute_gamma_corrections(zeros):
    log("  Gamma 보정 사전 계산 중 ...")
    t0 = time.time()
    N = len(zeros)
    im_gamma = np.zeros(N)
    re_gamma = np.zeros(N)
    for i in range(N):
        im_g, re_g = gamma_corrections(zeros[i])
        im_gamma[i] = im_g
        re_gamma[i] = re_g
        if (i+1) % 2000 == 0:
            log(f"    {i+1}/{N} ({(i+1)/N*100:.0f}%) — {time.time()-t0:.1f}s")
    log(f"  완료: {N}개, {time.time()-t0:.1f}s")
    return im_gamma, re_gamma


def compute_components(x, W_val, im_gamma_all, re_gamma_all, band_mask=None):
    """영점 배열에서 모든 성분을 계산. trim 적용."""
    N = len(x)
    lo = W_val
    hi = N - W_val - 1
    n_valid = hi - lo
    if n_valid < 30:
        return None

    g = x[lo+1:hi+1] - x[lo:hi]
    g = np.where(np.abs(g) < 1e-15, 1e-15, g)

    delta_L_n = x[lo:hi] - x[lo-1:hi-1]
    delta_L_n = np.where(np.abs(delta_L_n) < 1e-15, 1e-15, delta_L_n)

    delta_R_np1 = x[lo+2:hi+2] - x[lo+1:hi+1]
    delta_R_np1 = np.where(np.abs(delta_R_np1) < 1e-15, 1e-15, delta_R_np1)

    H1_shared = 1.0 / g**2

    H1_tail_n = np.zeros(n_valid)
    H1_tail_np1 = np.zeros(n_valid)
    S1_n = np.zeros(n_valid)
    S1_np1 = np.zeros(n_valid)

    S1_n += -1.0/g + 1.0/delta_L_n
    S1_np1 += -1.0/delta_R_np1 + 1.0/g

    for j in range(2, W_val + 1):
        idx_n = np.arange(lo, hi)
        diff_rj_n = x[idx_n + j] - x[idx_n]
        diff_rj_n = np.where(np.abs(diff_rj_n) < 1e-15, 1e-15, diff_rj_n)
        diff_lj_n = x[idx_n] - x[idx_n - j]
        diff_lj_n = np.where(np.abs(diff_lj_n) < 1e-15, 1e-15, diff_lj_n)

        H1_tail_n += 1.0/diff_rj_n**2 + 1.0/diff_lj_n**2
        S1_n += -1.0/diff_rj_n + 1.0/diff_lj_n

        idx_np1 = np.arange(lo+1, hi+1)
        diff_rj_np1 = x[idx_np1 + j] - x[idx_np1]
        diff_rj_np1 = np.where(np.abs(diff_rj_np1) < 1e-15, 1e-15, diff_rj_np1)
        diff_lj_np1 = x[idx_np1] - x[idx_np1 - j]
        diff_lj_np1 = np.where(np.abs(diff_lj_np1) < 1e-15, 1e-15, diff_lj_np1)

        H1_tail_np1 += 1.0/diff_rj_np1**2 + 1.0/diff_lj_np1**2
        S1_np1 += -1.0/diff_rj_np1 + 1.0/diff_lj_np1

    # Gamma 보정
    im_g_n = im_gamma_all[lo:hi]
    im_g_np1 = im_gamma_all[lo+1:hi+1]
    re_g_n = re_gamma_all[lo:hi]
    re_g_np1 = re_gamma_all[lo+1:hi+1]

    S1L_n = S1_n - im_g_n
    S1L_np1 = S1_np1 - im_g_np1

    H1_total_n = H1_shared + 1.0/delta_L_n**2 + H1_tail_n + re_g_n
    H1_total_np1 = H1_shared + 1.0/delta_R_np1**2 + H1_tail_np1 + re_g_np1

    A_L_n = S1L_n**2 + 2.0 * H1_total_n
    A_L_np1 = S1L_np1**2 + 2.0 * H1_total_np1

    shared_term = 2.0 / g**2
    A_L_prime_n = A_L_n - shared_term
    A_L_prime_np1 = A_L_np1 - shared_term

    t_n = x[lo:hi]

    # trim
    trim_lo = int(n_valid * TRIM_FRAC)
    trim_hi = n_valid - int(n_valid * TRIM_FRAC)
    sl = slice(trim_lo, trim_hi)

    # NaN/Inf 필터
    mask = (np.isfinite(A_L_n[sl]) & np.isfinite(A_L_np1[sl]) &
            np.isfinite(shared_term[sl]) & (A_L_n[sl] > 0) & (A_L_np1[sl] > 0) &
            (np.abs(g[sl]) > 1e-15))

    return {
        'A_L_n': A_L_n[sl][mask],
        'A_L_np1': A_L_np1[sl][mask],
        'A_L_prime_n': A_L_prime_n[sl][mask],
        'A_L_prime_np1': A_L_prime_np1[sl][mask],
        'shared': shared_term[sl][mask],
        'g': g[sl][mask],
        'S1L_sq_n': S1L_n[sl][mask]**2,
        'S1L_sq_np1': S1L_np1[sl][mask]**2,
        't_n': t_n[sl][mask],
        'n': int(np.sum(mask)),
    }


def cov_decompose(data):
    """Cov = Var(shared) + Cross + Cov(residual) 분해."""
    cov_full = np.cov(data['A_L_n'], data['A_L_np1'])[0, 1]
    var_shared = np.var(data['shared'])
    cov_resid = np.cov(data['A_L_prime_n'], data['A_L_prime_np1'])[0, 1]
    cross = cov_full - var_shared - cov_resid
    return cov_full, var_shared, cross, cov_resid


def compute_path_c_ratio(data):
    """Path C 비율 R = Var(2/g²) / (|Cross|+|Resid|)"""
    cov_full, var_shared, cross, cov_resid = cov_decompose(data)
    denom = abs(cross) + abs(cov_resid)
    if denom < 1e-30:
        return float('inf'), cov_full, var_shared, cross, cov_resid
    R = var_shared / denom
    return R, cov_full, var_shared, cross, cov_resid


# ══════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════
log("=" * 75)
log("[C-399] Path C 비율의 T→∞ 점근 안정성 검증")
log("=" * 75)
log(f"  T_MAX={T_MAX}, W={W}, TRIM={TRIM_FRAC}, N_BANDS={N_BANDS}, N_BOOTSTRAP={N_BOOTSTRAP}")
log(f"  밀도: d(t) = log(t/(2π))/(2π) (이론적)")
log()

t_start = time.time()

# ── 영점 수집 ──
zeros = get_zeta_zeros()
n_zeros = len(zeros)
log(f"  영점 수: {n_zeros}")

# ── Gamma 보정 ──
im_gamma_all, re_gamma_all = precompute_gamma_corrections(zeros)

# ── 전체 데이터 계산 ──
log()
log("=" * 75)
log("  [Phase 1] 전체 데이터 3-tier 분해")
log("=" * 75)

data_all = compute_components(zeros, W, im_gamma_all, re_gamma_all)
if data_all is None:
    log("FATAL: 데이터 부족")
    sys.exit(1)

n = data_all['n']
log(f"  n = {n} 쌍 (trim 후)")

R_all, cov_full, var_shared, cross, cov_resid = compute_path_c_ratio(data_all)

log(f"\n  Cov(A_Λ_n, A_Λ_{{n+1}}) = {cov_full:.4e}")
log(f"  Var(2/g²)               = {var_shared:.4e}  ({var_shared/cov_full*100:.1f}%)")
log(f"  Cross                   = {cross:.4e}  ({cross/cov_full*100:.1f}%)")
log(f"  Cov(residual)           = {cov_resid:.4e}  ({cov_resid/cov_full*100:.1f}%)")
log(f"  합계 검증: {(var_shared + cross + cov_resid)/cov_full:.6f} (= 1.0)")
log(f"\n  ★ Path C 비율 R(전체) = {R_all:.4f}")
log(f"  ★ Path C 성립: {'✅ YES' if R_all > 1 else '❌ NO'} (R > 1 필요)")

# ══════════════════════════════════════════════════════════════════
# Phase 2: T-대역별 Path C 비율
# ══════════════════════════════════════════════════════════════════
log()
log("=" * 75)
log(f"  [Phase 2] T-대역별 Path C 비율 ({N_BANDS}개 등간격)")
log("=" * 75)

t_arr = data_all['t_n']
t_min_data = t_arr.min()
t_max_data = t_arr.max()

band_edges = np.linspace(t_min_data, t_max_data, N_BANDS + 1)

log(f"\n  {'대역':>25s} │ {'n':>5s} │ {'Var(sh)%':>9s} │ {'Cross%':>8s} │ {'Resid%':>8s} │ {'R(T)':>7s} │ {'PathC':>5s}")
log("  " + "─" * 90)

band_results = []
for i in range(N_BANDS):
    t_lo, t_hi = band_edges[i], band_edges[i+1]
    if i == N_BANDS - 1:
        mask = (t_arr >= t_lo) & (t_arr <= t_hi)
    else:
        mask = (t_arr >= t_lo) & (t_arr < t_hi)
    n_band = np.sum(mask)

    if n_band < 30:
        log(f"  [{t_lo:8.0f}, {t_hi:8.0f}] │ {n_band:5d} │ (부족 — 최소 30쌍 필요)")
        continue

    # 대역 데이터
    band_data = {
        'A_L_n': data_all['A_L_n'][mask],
        'A_L_np1': data_all['A_L_np1'][mask],
        'A_L_prime_n': data_all['A_L_prime_n'][mask],
        'A_L_prime_np1': data_all['A_L_prime_np1'][mask],
        'shared': data_all['shared'][mask],
        'g': data_all['g'][mask],
        't_n': data_all['t_n'][mask],
        'n': int(n_band),
    }

    R_band, cf, vs, cr, cres = compute_path_c_ratio(band_data)
    t_mid = (t_lo + t_hi) / 2
    path_c_ok = R_band > 1

    band_results.append({
        'i': i+1,
        't_lo': t_lo, 't_hi': t_hi, 't_mid': t_mid,
        'n': int(n_band),
        'R': R_band,
        'cov_full': cf, 'var_shared': vs, 'cross': cr, 'cov_resid': cres,
        'path_c': path_c_ok,
    })

    vs_pct = vs/cf*100 if abs(cf) > 1e-30 else 0
    cr_pct = cr/cf*100 if abs(cf) > 1e-30 else 0
    cres_pct = cres/cf*100 if abs(cf) > 1e-30 else 0
    status = "✅" if path_c_ok else "❌"

    log(f"  [{t_lo:8.0f}, {t_hi:8.0f}] │ {n_band:5d} │ {vs_pct:8.1f}% │ {cr_pct:7.1f}% │ {cres_pct:7.1f}% │ {R_band:7.3f} │ {status}")

n_pass = sum(1 for b in band_results if b['path_c'])
n_total = len(band_results)
log(f"\n  Path C 성립: {n_pass}/{n_total} 대역")

# R 최소/최대/평균
if band_results:
    R_vals = [b['R'] for b in band_results]
    R_min = min(R_vals)
    R_max = max(R_vals)
    R_mean = np.mean(R_vals)
    R_min_band = [b for b in band_results if b['R'] == R_min][0]
    log(f"  R 최소: {R_min:.4f} (대역 {R_min_band['i']}: T∈[{R_min_band['t_lo']:.0f}, {R_min_band['t_hi']:.0f}])")
    log(f"  R 최대: {R_max:.4f}")
    log(f"  R 평균: {R_mean:.4f}")

# ══════════════════════════════════════════════════════════════════
# Phase 3: R(T) vs T 선형 회귀
# ══════════════════════════════════════════════════════════════════
log()
log("=" * 75)
log("  [Phase 3] R(T) vs T 선형 회귀")
log("=" * 75)

if len(band_results) >= 3:
    T_mids = np.array([b['t_mid'] for b in band_results])
    R_array = np.array([b['R'] for b in band_results])

    slope, intercept, r_value, p_value, std_err = stats.linregress(T_mids, R_array)
    log(f"\n  선형 회귀: R(T) = {slope:.6e} × T + {intercept:.4f}")
    log(f"  기울기: {slope:.6e} (±{std_err:.6e})")
    log(f"  기울기 부호: {'≥ 0 (비감소 ✅)' if slope >= 0 else '< 0 (감소 ⚠️)'}")
    log(f"  p-value (기울기=0 검정): {p_value:.4e}")
    log(f"  R²: {r_value**2:.4f}")
    log(f"  절편: {intercept:.4f}")

    # T=10000, T=50000, T=100000에서 예측
    log(f"\n  [외삽 예측]")
    for T_pred in [10000, 50000, 100000]:
        R_pred = slope * T_pred + intercept
        log(f"    R(T={T_pred:>6d}) = {R_pred:.4f} {'✅ > 1' if R_pred > 1 else '⚠️ ≤ 1'}")

    # 2차 회귀도 시도
    if len(band_results) >= 5:
        coeffs = np.polyfit(T_mids, R_array, 2)
        R_quad_pred = np.polyval(coeffs, T_mids)
        ss_res_lin = np.sum((R_array - (slope * T_mids + intercept))**2)
        ss_res_quad = np.sum((R_array - R_quad_pred)**2)
        log(f"\n  [2차 회귀]")
        log(f"  R(T) = {coeffs[0]:.8e}·T² + {coeffs[1]:.6e}·T + {coeffs[2]:.4f}")
        log(f"  2차 계수 부호: {'≥ 0 (볼록 — 발산 경향)' if coeffs[0] >= 0 else '< 0 (오목 — 수렴 경향)'}")
        log(f"  잔차제곱합: 선형={ss_res_lin:.4e}, 2차={ss_res_quad:.4e}")
else:
    log("  대역 부족으로 회귀 분석 불가 (<3)")

# ══════════════════════════════════════════════════════════════════
# Phase 4: Bootstrap 대역별 Path C 비율 안정성
# ══════════════════════════════════════════════════════════════════
log()
log("=" * 75)
log(f"  [Phase 4] Bootstrap Path C 비율 (B={N_BOOTSTRAP})")
log("=" * 75)

rng = np.random.default_rng(42)

# 전체 데이터 bootstrap
AL_n = data_all['A_L_n']
AL_np1 = data_all['A_L_np1']
ALp_n = data_all['A_L_prime_n']
ALp_np1 = data_all['A_L_prime_np1']
sh = data_all['shared']

boot_R = np.zeros(N_BOOTSTRAP)
boot_path_c = np.zeros(N_BOOTSTRAP)

for b in range(N_BOOTSTRAP):
    idx = rng.integers(0, n, n)
    cf = np.cov(AL_n[idx], AL_np1[idx])[0, 1]
    vs = np.var(sh[idx])
    cr_r = np.cov(ALp_n[idx], ALp_np1[idx])[0, 1]
    cr = cf - vs - cr_r
    denom = abs(cr) + abs(cr_r)
    boot_R[b] = vs / denom if denom > 1e-30 else float('inf')
    boot_path_c[b] = 1 if boot_R[b] > 1 else 0

log(f"\n  [전체 Bootstrap]")
log(f"  R 중앙값: {np.median(boot_R):.4f}")
log(f"  R 95% CI: [{np.percentile(boot_R, 2.5):.4f}, {np.percentile(boot_R, 97.5):.4f}]")
log(f"  Path C 성립 빈도: {np.mean(boot_path_c)*100:.1f}% ({int(np.sum(boot_path_c))}/{N_BOOTSTRAP})")

# 대역별 bootstrap
log(f"\n  [대역별 Bootstrap Path C 빈도]")
log(f"  {'대역':>25s} │ {'PathC%':>7s} │ {'R_median':>9s} │ {'R_CI_lo':>8s} │ {'R_CI_hi':>8s}")
log("  " + "─" * 65)

for br in band_results:
    t_lo, t_hi = br['t_lo'], br['t_hi']
    if br['i'] == N_BANDS:
        mask = (t_arr >= t_lo) & (t_arr <= t_hi)
    else:
        mask = (t_arr >= t_lo) & (t_arr < t_hi)
    n_b = np.sum(mask)
    if n_b < 30:
        continue

    b_AL_n = AL_n[mask]
    b_AL_np1 = AL_np1[mask]
    b_ALp_n = ALp_n[mask]
    b_ALp_np1 = ALp_np1[mask]
    b_sh = sh[mask]

    b_R = np.zeros(N_BOOTSTRAP)
    b_pc = np.zeros(N_BOOTSTRAP)
    for b in range(N_BOOTSTRAP):
        idx = rng.integers(0, n_b, n_b)
        cf = np.cov(b_AL_n[idx], b_AL_np1[idx])[0, 1]
        vs = np.var(b_sh[idx])
        cr_r = np.cov(b_ALp_n[idx], b_ALp_np1[idx])[0, 1]
        cr = cf - vs - cr_r
        denom = abs(cr) + abs(cr_r)
        b_R[b] = vs / denom if denom > 1e-30 else float('inf')
        b_pc[b] = 1 if b_R[b] > 1 else 0

    pct = np.mean(b_pc) * 100
    r_med = np.median(b_R)
    r_lo = np.percentile(b_R, 2.5)
    r_hi = np.percentile(b_R, 97.5)
    log(f"  [{t_lo:8.0f}, {t_hi:8.0f}] │ {pct:6.1f}% │ {r_med:9.3f} │ {r_lo:8.3f} │ {r_hi:8.3f}")

# ══════════════════════════════════════════════════════════════════
# Phase 5: 성공 기준 판정
# ══════════════════════════════════════════════════════════════════
log()
log("=" * 75)
log("  [Phase 5] 성공 기준 판정")
log("=" * 75)
log()

# 기준 1: 10개 대역 전부 R > 1
all_pass = all(b['path_c'] for b in band_results)
log(f"  기준 1: 전 대역 R > 1? {'✅ YES' if all_pass else '❌ NO'} ({n_pass}/{n_total})")

# 기준 2: 선형 회귀 기울기 ≥ 0
if len(band_results) >= 3:
    log(f"  기준 2: 기울기 ≥ 0? {'✅ YES' if slope >= 0 else '⚠️ NO'} (기울기={slope:.2e}, p={p_value:.2e})")

# 기준 3: R 최소값 > 1.05
if band_results:
    R_min_check = R_min > 1.05
    log(f"  기준 3: R_min > 1.05? {'✅ YES' if R_min_check else '⚠️ NO'} (R_min={R_min:.4f})")

# 기준 4 (음성): 3개 이상 대역에서 R < 1
n_fail = sum(1 for b in band_results if not b['path_c'])
negative = n_fail >= 3
log(f"  기준 4 (음성): 3개+ 대역 R < 1? {'❌ YES — 전략 재고!' if negative else '✅ NO'} ({n_fail}/{n_total} 실패)")

# ── 종합 판정 ──
log()
if all_pass and (len(band_results) < 3 or slope >= 0) and R_min > 1.05:
    verdict = "강한 양성 ✅ — 해석적 증명 가능"
elif n_pass == n_total and R_min > 1.0:
    verdict = "양성 ✅ — R > 1 전 대역, 마진은 좁음"
elif n_pass >= n_total * 0.7:
    verdict = "조건부 양성 ⚠️ — 대부분 성립, 일부 불안정"
elif negative:
    verdict = "음성 ❌ — Path C 전략 재고 필요"
else:
    verdict = "미결정 — 추가 데이터 필요"

log(f"  ★ 종합 판정: {verdict}")

# ══════════════════════════════════════════════════════════════════
# Phase 6: T-대역별 요약 테이블 (수학자 보고용)
# ══════════════════════════════════════════════════════════════════
log()
log("=" * 75)
log("  [Phase 6] 요약 테이블")
log("=" * 75)
log()

log(f"  | 대역 | T 범위 | n쌍 | Var(sh)% | Cross% | Resid% | R(T) | PathC |")
log(f"  |------|--------|------|----------|--------|--------|------|-------|")
for br in band_results:
    cf = br['cov_full']
    vs_p = br['var_shared']/cf*100 if abs(cf) > 1e-30 else 0
    cr_p = br['cross']/cf*100 if abs(cf) > 1e-30 else 0
    cres_p = br['cov_resid']/cf*100 if abs(cf) > 1e-30 else 0
    log(f"  | {br['i']:4d} | [{br['t_lo']:.0f}, {br['t_hi']:.0f}] | {br['n']:4d} | {vs_p:7.1f}% | {cr_p:6.1f}% | {cres_p:6.1f}% | {br['R']:.3f} | {'✅' if br['path_c'] else '❌'} |")

log()
elapsed = time.time() - t_start
log(f"총 소요: {elapsed:.1f}s")
log("=" * 75)

out_f.close()
print(f"\n결과 저장: {RESULT_PATH}")
