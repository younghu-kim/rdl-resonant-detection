#!/usr/bin/env python3
"""
[C-398] Path C 비율 R(T)의 T→∞ 점근 안정성 검증

목표: ζ(s) T=10000까지 영점을 수집하고, 10개 등간격 T-대역에서
      R(T) = Var(2/g²) / (|Cross|+|Resid|) 를 계산하여 추세 분석.

핵심 질문: R(T)가 T→∞에서:
  (a) 상수 c > 1에 수렴 → 해석적 증명 가능
  (b) 1에 수렴 → Path C는 수치적으로만 유효
  (c) 발산 → 더 강한 결과 가능

최적화: lfunzeros를 구간별로 호출하여 속도 개선 (T=10000 한 번에 → 느림).

체크리스트:
  [x] python -u
  [x] mpmath dps=50
  [x] trim 20%
  [x] NaN/Inf 체크
  [x] 이론적 밀도 d(t) = log(t/(2π))/(2π)
  [x] W=100 (C-395와 동일)
  [x] Gamma 보정 (A_Λ)
  [x] 선형 회귀 + 추세 분석
  [x] lfunzeros 구간 분할 최적화
"""

import sys, os, time, math
import numpy as np
from scipy import stats as sp_stats

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
N_BANDS = 10
N_BOOTSTRAP = 2000
RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/path_c_asymptotic_c398.txt'
)

out_f = open(RESULT_PATH, 'w')
def log(msg=''):
    print(msg, flush=True)
    out_f.write(msg + '\n')
    out_f.flush()


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
    """ζ(s) 영점 수집 — 벡터 일괄 추출 (원소별 PARI 호출 회피)"""
    log(f"[ζ(s)] 영점 수집 t∈(0, {T_MAX}] ...")
    t0 = time.time()

    pari('L_zeta = lfuncreate(1)')
    pari(f'Li_z = lfuninit(L_zeta, [0, {int(T_MAX) + 5}])')
    log(f"  lfuninit 완료 ({time.time()-t0:.0f}s)")

    pari(f'zv = lfunzeros(Li_z, {T_MAX})')
    n = int(str(pari('#zv')))
    log(f"  lfunzeros 완료: {n}개 ({time.time()-t0:.0f}s)")

    # 벡터를 문자열로 일괄 변환하여 파싱 (원소별 호출 대비 수십배 빠름)
    raw = str(pari('zv'))
    log(f"  벡터 문자열 변환 ({time.time()-t0:.0f}s)")
    # 형식: [14.1347..., 21.0220..., ...] 또는 공백 구분
    raw = raw.strip('[]')
    parts = [p.strip() for p in raw.replace(',', ' ').split() if p.strip()]
    zeros = []
    for p in parts:
        p = p.replace(' E', 'e').replace('E ', 'e')
        try:
            t = float(p)
            if t > 0.5:
                zeros.append(t)
        except ValueError:
            pass

    zeros = np.array(sorted(zeros))
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
        if (i + 1) % 2000 == 0:
            log(f"    {i+1}/{N} ({time.time()-t0:.0f}s)")
    log(f"  완료: {time.time()-t0:.1f}s")
    return im_gamma, re_gamma


def compute_components(x, W_val, im_gamma_all, re_gamma_all):
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

    trim_lo = int(n_valid * TRIM_FRAC)
    trim_hi = n_valid - int(n_valid * TRIM_FRAC)
    sl = slice(trim_lo, trim_hi)

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
        't_n': t_n[sl][mask],
        'n': int(np.sum(mask)),
    }


def cov_decompose(A_n, A_np1, Ap_n, Ap_np1, shared):
    cov_full = np.cov(A_n, A_np1)[0, 1]
    var_shared = np.var(shared)
    cov_resid = np.cov(Ap_n, Ap_np1)[0, 1]
    cross = cov_full - var_shared - cov_resid
    return cov_full, var_shared, cross, cov_resid


def path_c_ratio(var_sh, cross, resid):
    denom = abs(cross) + abs(resid)
    if denom < 1e-30:
        return float('inf')
    return var_sh / denom


# ══════════════════════════════════════════════════════════════════
log("=" * 70)
log("[C-398] Path C 비율 R(T)의 T→∞ 점근 안정성 검증")
log("=" * 70)
log(f"  T_MAX={T_MAX}, W={W}, TRIM={TRIM_FRAC}, N_BANDS={N_BANDS}")
log(f"  N_BOOTSTRAP={N_BOOTSTRAP}")
log()

t_start = time.time()

zeros = get_zeta_zeros()
im_gamma_all, re_gamma_all = precompute_gamma_corrections(zeros)

# ── 전체 데이터 기본 분해 ──────────────────────────────────────
log()
log("=" * 70)
log("  [Phase 1] 전체 데이터 기본 분해")
log("=" * 70)

data_all = compute_components(zeros, W, im_gamma_all, re_gamma_all)
if data_all is None:
    log("FATAL: 데이터 부족")
    sys.exit(1)

n_all = data_all['n']
log(f"  n = {n_all} 쌍")

cov_full, var_shared, cross, cov_resid = cov_decompose(
    data_all['A_L_n'], data_all['A_L_np1'],
    data_all['A_L_prime_n'], data_all['A_L_prime_np1'],
    data_all['shared']
)
R_global = path_c_ratio(var_shared, cross, cov_resid)

log(f"\n  Cov(A_Λ_n, A_Λ_{{n+1}}) = {cov_full:.4e}")
log(f"  Var(2/g²)               = {var_shared:.4e}  ({var_shared/cov_full*100:.1f}%)")
log(f"  Cross                   = {cross:.4e}  ({cross/cov_full*100:.1f}%)")
log(f"  Cov(residual)           = {cov_resid:.4e}  ({cov_resid/cov_full*100:.1f}%)")
log(f"  합계 검증: {(var_shared + cross + cov_resid)/cov_full:.6f}")
log(f"  Path C 비율 R = {R_global:.3f}×")

# ── T-대역별 분석 ──────────────────────────────────────────────
log()
log("=" * 70)
log(f"  [Phase 2] T-대역별 Path C 비율 ({N_BANDS}개 대역)")
log("=" * 70)

t_arr = data_all['t_n']
t_min_data = t_arr.min()
t_max_data = t_arr.max()

band_edges = np.linspace(t_min_data, t_max_data, N_BANDS + 1)

header = f"  {'대역':>24s} │ {'n':>5s} │ {'Var%':>7s} │ {'Cross%':>8s} │ {'Resid%':>8s} │ {'R(T)':>7s} │ {'R>1':>4s}"
log(f"\n{header}")
log("  " + "─" * 85)

band_results = []
for i in range(N_BANDS):
    t_lo, t_hi = band_edges[i], band_edges[i+1]
    if i == N_BANDS - 1:
        mask = (t_arr >= t_lo) & (t_arr <= t_hi)
    else:
        mask = (t_arr >= t_lo) & (t_arr < t_hi)
    n_band = np.sum(mask)
    if n_band < 30:
        log(f"  [{t_lo:8.0f}, {t_hi:8.0f}] │ {n_band:5d} │ (부족)")
        continue

    cf, vs, cr, cres = cov_decompose(
        data_all['A_L_n'][mask], data_all['A_L_np1'][mask],
        data_all['A_L_prime_n'][mask], data_all['A_L_prime_np1'][mask],
        data_all['shared'][mask]
    )
    R = path_c_ratio(vs, cr, cres)
    t_mid = (t_lo + t_hi) / 2

    status = "Y" if R > 1 else "N"
    log(f"  [{t_lo:8.0f}, {t_hi:8.0f}] │ {n_band:5d} │ {vs/cf*100:6.1f}% │ {cr/cf*100:+7.1f}% │ {cres/cf*100:+7.1f}% │ {R:6.2f}x │  {status}")

    band_results.append({
        't_lo': t_lo, 't_hi': t_hi, 't_mid': t_mid,
        'n': n_band, 'cov': cf, 'var_sh': vs, 'cross': cr, 'resid': cres,
        'R': R
    })

n_bands_valid = len(band_results)
n_R_above_1 = sum(1 for b in band_results if b['R'] > 1)

log(f"\n  R > 1: {n_R_above_1}/{n_bands_valid} 대역")

# ── R(T) 추세 분석 ─────────────────────────────────────────────
log()
log("=" * 70)
log("  [Phase 3] R(T) 추세 분석")
log("=" * 70)

if n_bands_valid >= 3:
    t_mids = np.array([b['t_mid'] for b in band_results])
    R_vals = np.array([b['R'] for b in band_results])

    # 선형 회귀: R = a·T + b
    slope, intercept, r_val, p_val, std_err = sp_stats.linregress(t_mids, R_vals)
    log(f"\n  선형: R(T) = {slope:+.6f}*T + {intercept:.4f}")
    log(f"    기울기: {slope:+.6f} +/- {std_err:.6f}, r^2={r_val**2:.4f}, p={p_val:.4e}")
    log(f"    추세: {'비감소' if slope >= 0 else '감소'}")

    # log(T) 회귀
    log_t = np.log(t_mids)
    sl_log, int_log, r_log, p_log, se_log = sp_stats.linregress(log_t, R_vals)
    log(f"\n  로그: R(T) = {sl_log:+.4f}*log(T) + {int_log:.4f}")
    log(f"    기울기: {sl_log:+.4f} +/- {se_log:.4f}, r^2={r_log**2:.4f}, p={p_log:.4e}")

    # 1/T 회귀 (수렴 여부)
    inv_t = 1.0 / t_mids
    sl_inv, int_inv, r_inv, p_inv, se_inv = sp_stats.linregress(inv_t, R_vals)
    log(f"\n  1/T: R(T) = {sl_inv:+.1f}*(1/T) + {int_inv:.4f}")
    log(f"    T->inf 극한: {int_inv:.4f}, r^2={r_inv**2:.4f}, p={p_inv:.4e}")

    R_min = R_vals.min()
    R_max = R_vals.max()
    R_mean = R_vals.mean()
    R_std = R_vals.std()
    log(f"\n  통계: min={R_min:.3f}, max={R_max:.3f}, mean={R_mean:.3f}, std={R_std:.3f}")

# ── Bootstrap 안정성 ───────────────────────────────────────────
log()
log("=" * 70)
log(f"  [Phase 4] Bootstrap 대역별 (B={N_BOOTSTRAP})")
log("=" * 70)

rng = np.random.default_rng(42)

for bi, b in enumerate(band_results):
    mask = (t_arr >= b['t_lo'])
    if bi < len(band_results) - 1:
        mask = mask & (t_arr < b['t_hi'])
    else:
        mask = mask & (t_arr <= b['t_hi'])

    An = data_all['A_L_n'][mask]
    Anp1 = data_all['A_L_np1'][mask]
    Apn = data_all['A_L_prime_n'][mask]
    Apnp1 = data_all['A_L_prime_np1'][mask]
    sh = data_all['shared'][mask]
    nb = len(An)

    boot_R = np.zeros(N_BOOTSTRAP)
    for bb in range(N_BOOTSTRAP):
        idx = rng.integers(0, nb, nb)
        cf = np.cov(An[idx], Anp1[idx])[0, 1]
        vs = np.var(sh[idx])
        cres = np.cov(Apn[idx], Apnp1[idx])[0, 1]
        cr = cf - vs - cres
        denom = abs(cr) + abs(cres)
        boot_R[bb] = vs / denom if denom > 1e-30 else 999

    pct = np.mean(boot_R > 1) * 100
    ci_lo = np.percentile(boot_R, 2.5)
    ci_hi = np.percentile(boot_R, 97.5)
    log(f"  [{b['t_lo']:8.0f}, {b['t_hi']:8.0f}] R={b['R']:.2f}x, R>1: {pct:.1f}%, 95%CI: [{ci_lo:.2f}, {ci_hi:.2f}]")

# ── 최종 판정 ─────────────────────────────────────────────────
log()
log("=" * 70)
log("  [Phase 5] 최종 판정")
log("=" * 70)
log()

if n_bands_valid > 0:
    R_vals = np.array([b['R'] for b in band_results])
    R_min = R_vals.min()
    n_fail = np.sum(R_vals < 1)

    log(f"  R(T) 전체: {n_R_above_1}/{n_bands_valid} 대역 R > 1")
    log(f"  R(T) 최소: {R_min:.3f}x")
    log(f"  R(T) 전역: {R_global:.3f}x")

    if n_bands_valid >= 3:
        log(f"  기울기 (선형): {slope:+.6f} (p={p_val:.4e})")
        log(f"  T->inf 극한 (1/T 적합): {int_inv:.4f}")
    log()

    if n_fail == 0 and R_min > 1.05:
        if n_bands_valid >= 3 and slope >= -abs(std_err):
            verdict = "강한 양성"
        else:
            verdict = "양성"
    elif n_fail == 0 and R_min > 1.0:
        verdict = "조건부 양성 (마진 좁음)"
    elif n_fail <= 2:
        verdict = f"약한 양성 ({n_fail}개 대역 R<1)"
    else:
        verdict = f"음성 ({n_fail}개 대역 R<1)"

    log(f"  ** 판정: {verdict}")
    log()

    if n_bands_valid >= 3:
        log("  [해석적 증명 feasibility]")
        if int_inv > 1.0:
            log(f"    T->inf 극한 R = {int_inv:.3f} > 1 — Path C bound 존재 가능")
        else:
            log(f"    T->inf 극한 R = {int_inv:.3f} <= 1 — 유한 T 전용")

        if slope >= 0:
            log(f"    기울기 >= 0 — 안전 마진 비감소")
        elif abs(slope) < std_err:
            log(f"    기울기 ~ 0 (유의하지 않음) — R(T) 거의 상수")
        else:
            log(f"    기울기 < 0 — 안전 마진 감소 (율: {slope:.6f}/T)")

elapsed = time.time() - t_start
log(f"\n총 소요: {elapsed:.1f}s")
log("=" * 70)

out_f.close()
print(f"\n결과 저장: {RESULT_PATH}")
