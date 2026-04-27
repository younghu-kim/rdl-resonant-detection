#!/usr/bin/env python3
"""
[C-397] Path C 교차검증 — Dirichlet L-함수 공분산 3-tier 분해

목표: ζ(s) T=2000에서 확립된 Path C 부등식
  Var(2/g²) > |Cross| + |Resid|
이 Dirichlet L-함수에서도 성립하는지 검증.

대상:
  1. χ mod 5 [1] (원시, order 4, odd, gammaV=[1])
  2. χ mod 7 [1] (원시, order 6, odd, gammaV=[1])

방법: C-395 proof_phase25_cross_sign_c395.py와 동일한 3-tier 분해.
  - Cov(A_Λ_n, A_Λ_{n+1}) = Var(2/g²) + Cross + Residual
  - Path C: Var(2/g²) > |Cross| + |Resid| → Cov > 0 절대 보장
  - Bootstrap 2000회 → Cross > 0 빈도 + 95% CI

밀도: d̄(t) = (1/2π)·log(q·t/(2π))  (conductor q 반영)
Gamma: gammaV에 따라 ψ, ψ₁ 보정

체크리스트:
  [x] python -u
  [x] mpmath dps=50
  [x] trim 20%
  [x] NaN/Inf 체크
  [x] except Exception as e: print(...)
  [x] PARI lfuncreate([znstar(q,1), [chi_idx]])
  [x] 이론적 밀도 사용 (경험적 d_bar 금지)
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
    pari.allocatemem(1024 * 10**6)
    pari.set_real_precision(80)
    print("cypari2 OK", flush=True)
except Exception as e:
    print(f"FATAL cypari2: {e}", flush=True)
    sys.exit(1)

# ── 설정 ──────────────────────────────────────────────────────────
T_MAX = 500.0      # Dirichlet: T=500 → ~400 zeros
TRIM_FRAC = 0.20
W = 50              # W=50 (영점 수 ζ보다 적으므로 축소)
SIGMA_C = 0.5
N_BOOTSTRAP = 2000
RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/path_c_dirichlet_c397.txt'
)

# Dirichlet 지표 설정
CHARACTERS = [
    {
        'label': 'χ mod 5 [1]',
        'q': 5,
        'chi_idx': [1],
        'gammaV': [1],   # odd
        'order': 4,
    },
    {
        'label': 'χ mod 7 [1]',
        'q': 7,
        'chi_idx': [1],
        'gammaV': [1],   # odd
        'order': 6,
    },
]

out_f = open(RESULT_PATH, 'w')
def log(msg=''):
    print(msg, flush=True)
    out_f.write(msg + '\n')
    out_f.flush()


def dirichlet_density(t, q):
    """이론적 영점 밀도: d̄(t) = (1/2π)·log(q·t/(2π))"""
    if t < 2.0:
        return 1.0
    return math.log(q * t / (2.0 * math.pi)) / (2.0 * math.pi)


def gamma_corrections(t0, gammaV):
    """Gamma 보정: ψ(s+μ)/2 합산"""
    s = mpmath.mpc(SIGMA_C, t0)
    im_sum = mpmath.mpf(0)
    re_sum = mpmath.mpf(0)
    for mu in gammaV:
        arg = (s + mu) / 2
        psi_val = mpmath.digamma(arg)
        psi1_val = mpmath.psi(1, arg)
        im_sum += mpmath.im(psi_val) / 2
        re_sum += mpmath.re(psi1_val) / 4
    return float(im_sum), float(re_sum)


def collect_zeros_pari(q, chi_idx, var_prefix):
    """PARI로 Dirichlet L-함수 영점 수집."""
    G_var = f"G_{var_prefix}"
    L_var = f"L_{var_prefix}"
    Li_var = f"Li_{var_prefix}"
    zv_var = f"zv_{var_prefix}"

    chi_str = str(chi_idx)
    pari(f'{G_var} = znstar({q}, 1)')
    pari(f'{L_var} = lfuncreate([{G_var}, {chi_str}])')
    pari(f'{Li_var} = lfuninit({L_var}, [0, {int(T_MAX) + 20}])')
    pari(f'{zv_var} = lfunzeros({Li_var}, {T_MAX})')

    n_z = int(str(pari(f'#{zv_var}')))
    zeros = []
    for i in range(1, n_z + 1):
        try:
            t_str = str(pari(f'{zv_var}[{i}]')).strip().replace(' E', 'e')
            t = float(t_str)
            if t > 2.0:
                zeros.append(t)
        except Exception as e:
            if len(zeros) == 0:
                log(f"  WARNING: 영점 파싱 실패 i={i}: {e}")

    zeros = np.array(sorted(zeros))
    if len(zeros) == 0:
        log("⚠️ 영점 0개 — 탐색 로직 점검 필요")
    return zeros


def precompute_gamma_corrections(zeros, gammaV):
    """Gamma 보정 사전 계산."""
    N = len(zeros)
    im_gamma = np.zeros(N)
    re_gamma = np.zeros(N)
    for i in range(N):
        im_g, re_g = gamma_corrections(zeros[i], gammaV)
        im_gamma[i] = im_g
        re_gamma[i] = re_g
    return im_gamma, re_gamma


def compute_components(x, W_val, im_gamma_all, re_gamma_all):
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
        'H1_tail_n': (H1_tail_n + re_g_n)[sl][mask],
        'H1_tail_np1': (H1_tail_np1 + re_g_np1)[sl][mask],
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


def bootstrap_analysis(data, n_boot=N_BOOTSTRAP):
    """Bootstrap 2000회: Cross > 0 빈도, 95% CI."""
    rng = np.random.default_rng(42)
    n = data['n']
    AL_n = data['A_L_n']
    AL_np1 = data['A_L_np1']
    ALp_n = data['A_L_prime_n']
    ALp_np1 = data['A_L_prime_np1']
    sh = data['shared']

    boot_cross_ratio = np.zeros(n_boot)
    boot_cross_sign = np.zeros(n_boot)
    boot_cov_full = np.zeros(n_boot)
    boot_path_c = np.zeros(n_boot)  # Var(shared) > |Cross|+|Resid|

    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        cf = np.cov(AL_n[idx], AL_np1[idx])[0, 1]
        vs = np.var(sh[idx])
        cr_r = np.cov(ALp_n[idx], ALp_np1[idx])[0, 1]
        cr = cf - vs - cr_r
        boot_cross_ratio[b] = cr / cf if abs(cf) > 1e-30 else 0
        boot_cross_sign[b] = 1 if cr > 0 else 0
        boot_cov_full[b] = cf
        boot_path_c[b] = 1 if vs > abs(cr) + abs(cr_r) else 0

    return {
        'cross_ratio': boot_cross_ratio,
        'cross_sign': boot_cross_sign,
        'cov_full': boot_cov_full,
        'path_c': boot_path_c,
        'pct_cross_pos': np.mean(boot_cross_sign) * 100,
        'ci_lo': np.percentile(boot_cross_ratio, 2.5),
        'ci_hi': np.percentile(boot_cross_ratio, 97.5),
        'pct_path_c': np.mean(boot_path_c) * 100,
    }


def analyze_character(char_info):
    """단일 Dirichlet character에 대한 전체 분석."""
    label = char_info['label']
    q = char_info['q']
    chi_idx = char_info['chi_idx']
    gammaV = char_info['gammaV']

    log()
    log("=" * 70)
    log(f"  [{label}] — conductor q={q}, order={char_info['order']}, gammaV={gammaV}")
    log("=" * 70)

    # 1. 영점 수집
    t0 = time.time()
    log(f"\n  영점 수집 (T={T_MAX}) ...")
    prefix = f"q{q}_c{chi_idx[0]}"
    zeros = collect_zeros_pari(q, chi_idx, prefix)
    log(f"  {len(zeros)}개 영점, t∈[{zeros[0]:.3f}, {zeros[-1]:.3f}], {time.time()-t0:.1f}s")

    if len(zeros) < 100:
        log(f"  ⚠️ 영점 부족 ({len(zeros)} < 100) — 건너뜀")
        return None

    # 2. Gamma 보정
    log(f"  Gamma 보정 사전 계산 (gammaV={gammaV}) ...")
    t1 = time.time()
    im_gamma_all, re_gamma_all = precompute_gamma_corrections(zeros, gammaV)
    log(f"  완료: {time.time()-t1:.1f}s")

    # 3. 성분 계산
    data = compute_components(zeros, W, im_gamma_all, re_gamma_all)
    if data is None:
        log(f"  ⚠️ 데이터 부족 (W={W}에서 30쌍 미만) — 건너뜀")
        return None

    n = data['n']
    log(f"\n  n = {n} 쌍 (W={W}, trim {TRIM_FRAC*100:.0f}%)")

    # 4. 3-tier 분해
    cov_full, var_shared, cross, cov_resid = cov_decompose(data)

    log(f"\n  Cov(A_Λ_n, A_Λ_{{n+1}}) = {cov_full:.4e}")
    log(f"  Var(2/g²)               = {var_shared:.4e}  ({var_shared/cov_full*100:.1f}%)")
    log(f"  Cross                   = {cross:.4e}  ({cross/cov_full*100:.1f}%)")
    log(f"  Cov(residual)           = {cov_resid:.4e}  ({cov_resid/cov_full*100:.1f}%)")
    log(f"  합계 검증: {(var_shared + cross + cov_resid)/cov_full:.6f} (= 1.0)")

    # Path C 부등식
    path_c_holds = var_shared > abs(cross) + abs(cov_resid)
    safety_margin = var_shared / (abs(cross) + abs(cov_resid)) if (abs(cross) + abs(cov_resid)) > 1e-30 else float('inf')
    margin_pct = (var_shared - abs(cross) - abs(cov_resid)) / cov_full * 100 if cov_full > 1e-30 else 0

    log(f"\n  ── Path C 부등식 ──")
    log(f"  Var(2/g²) = {var_shared:.4e}")
    log(f"  |Cross| + |Resid| = {abs(cross) + abs(cov_resid):.4e}")
    if path_c_holds:
        log(f"  → ✅ Var(2/g²) > |Cross|+|Resid| — Path C 성립! (비율 {safety_margin:.2f}×, 안전 마진 {margin_pct:.1f}%)")
    else:
        log(f"  → ❌ Path C 미성립 (비율 {safety_margin:.2f}×)")

    # 5. T-대역별 Cross 부호
    log(f"\n  ── T-대역별 Cross 부호 ──")
    t_arr = data['t_n']
    N_BANDS = 3  # 영점 수가 적으므로 3개 대역
    band_edges = np.linspace(t_arr.min(), t_arr.max(), N_BANDS + 1)

    log(f"  {'대역':>20s} │ {'n':>5s} │ {'Cov':>10s} │ {'Var(sh)':>10s} │ {'Cross':>10s} │ {'Resid':>10s} │ {'PathC':>7s}")
    log("  " + "─" * 85)

    band_results = []
    for i in range(N_BANDS):
        t_lo, t_hi = band_edges[i], band_edges[i+1]
        mask = (t_arr >= t_lo) & (t_arr < t_hi) if i < N_BANDS - 1 else (t_arr >= t_lo)
        n_band = np.sum(mask)
        if n_band < 15:
            log(f"  [{t_lo:.0f}, {t_hi:.0f}] │ {n_band:5d} │ (부족)")
            continue

        cov_f = np.cov(data['A_L_n'][mask], data['A_L_np1'][mask])[0, 1]
        var_s = np.var(data['shared'][mask])
        cov_r = np.cov(data['A_L_prime_n'][mask], data['A_L_prime_np1'][mask])[0, 1]
        cr = cov_f - var_s - cov_r
        pc = var_s > abs(cr) + abs(cov_r)
        band_results.append({'cross_pos': cr > 0, 'path_c': pc})

        pc_str = "✅" if pc else "❌"
        log(f"  [{t_lo:7.0f}, {t_hi:7.0f}] │ {n_band:5d} │ {cov_f:10.3e} │ {var_s:10.3e} │ {cr:10.3e} │ {cov_r:10.3e} │ {pc_str}")

    n_cross_pos = sum(1 for b in band_results if b['cross_pos'])
    n_path_c = sum(1 for b in band_results if b['path_c'])
    log(f"\n  Cross > 0: {n_cross_pos}/{len(band_results)} 대역")
    log(f"  Path C 성립: {n_path_c}/{len(band_results)} 대역")

    # 6. Bootstrap
    log(f"\n  ── Bootstrap (B={N_BOOTSTRAP}) ──")
    boot = bootstrap_analysis(data)

    log(f"  Cross/Cov 비율: 중앙={np.median(boot['cross_ratio']):.4f}")
    log(f"  95% CI: [{boot['ci_lo']:.4f}, {boot['ci_hi']:.4f}]")
    log(f"  Cross > 0 빈도: {boot['pct_cross_pos']:.1f}%")
    log(f"  Path C 성립 빈도: {boot['pct_path_c']:.1f}%")
    log(f"  Cov > 0 빈도: {np.mean(boot['cov_full'] > 0)*100:.1f}%")

    return {
        'label': label,
        'q': q,
        'n_zeros': len(zeros),
        'n_pairs': n,
        'cov_full': cov_full,
        'var_shared': var_shared,
        'cross': cross,
        'cov_resid': cov_resid,
        'path_c_holds': path_c_holds,
        'safety_margin': safety_margin,
        'margin_pct': margin_pct,
        'shared_pct': var_shared / cov_full * 100,
        'cross_pct': cross / cov_full * 100,
        'resid_pct': cov_resid / cov_full * 100,
        'boot_cross_pos': boot['pct_cross_pos'],
        'boot_ci': (boot['ci_lo'], boot['ci_hi']),
        'boot_path_c': boot['pct_path_c'],
        'n_band_cross_pos': n_cross_pos,
        'n_band_total': len(band_results),
        'n_band_path_c': n_path_c,
    }


# ══════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════
log("=" * 70)
log("[C-397] Path C 교차검증 — Dirichlet L-함수 공분산 3-tier 분해")
log("=" * 70)
log(f"  T_MAX={T_MAX}, W={W}, TRIM={TRIM_FRAC}, N_BOOTSTRAP={N_BOOTSTRAP}")
log()

t_global = time.time()

all_results = []
for char_info in CHARACTERS:
    try:
        result = analyze_character(char_info)
        if result is not None:
            all_results.append(result)
    except Exception as e:
        log(f"\n  ❌ {char_info['label']} 분석 실패: {e}")
        import traceback
        traceback.print_exc()

# ══════════════════════════════════════════════════════════════════
# 비교표: ζ(s) vs Dirichlet
# ══════════════════════════════════════════════════════════════════
log()
log("=" * 70)
log("  [비교표] ζ(s) T=2000 vs Dirichlet L-함수")
log("=" * 70)
log()

# C-395 ζ(s) 결과 (하드코딩)
zeta_ref = {
    'label': 'ζ(s) T=2000',
    'n_pairs': 790,
    'shared_pct': 57.7,
    'cross_pct': 36.7,
    'resid_pct': 5.6,
    'path_c_holds': True,
    'safety_margin': 1.37,
    'margin_pct': 15.4,
    'boot_cross_pos': 99.3,
    'boot_ci': (0.1355, 0.4103),
    'boot_path_c': None,  # C-395에서 미계산
}

header = f"  {'L-함수':>20s} │ {'n쌍':>5s} │ {'Var(sh)%':>9s} │ {'Cross%':>8s} │ {'Resid%':>8s} │ {'PathC':>5s} │ {'마진':>6s} │ {'안전마진%':>8s} │ {'Boot>0%':>8s} │ {'95%CI':>16s}"
log(header)
log("  " + "─" * 120)

def fmt_row(r):
    pc_str = "✅" if r['path_c_holds'] else "❌"
    ci_str = f"[{r['boot_ci'][0]:.3f},{r['boot_ci'][1]:.3f}]" if r['boot_ci'] else "N/A"
    margin_str = f"{r['margin_pct']:.1f}%" if r.get('margin_pct') is not None else "N/A"
    return f"  {r['label']:>20s} │ {r.get('n_pairs','?'):>5} │ {r['shared_pct']:>8.1f}% │ {r['cross_pct']:>7.1f}% │ {r['resid_pct']:>7.1f}% │ {pc_str:>5s} │ {r['safety_margin']:>5.2f}× │ {margin_str:>8s} │ {r['boot_cross_pos']:>7.1f}% │ {ci_str:>16s}"

log(fmt_row(zeta_ref))
for r in all_results:
    log(fmt_row(r))

# ══════════════════════════════════════════════════════════════════
# 최종 판정
# ══════════════════════════════════════════════════════════════════
log()
log("=" * 70)
log("  [최종 판정]")
log("=" * 70)
log()

n_tested = len(all_results)
n_path_c_pass = sum(1 for r in all_results if r['path_c_holds'])
n_cross_pos = sum(1 for r in all_results if r['boot_cross_pos'] > 90)

log(f"  테스트된 Dirichlet L-함수: {n_tested}개")
log(f"  Path C 성립: {n_path_c_pass}/{n_tested}")
log(f"  Cross > 0 (bootstrap >90%): {n_cross_pos}/{n_tested}")
log()

if n_path_c_pass == n_tested and n_tested >= 2:
    log("  ★ PATH C 일반화 양성:")
    log("    Var(2/g²) > |Cross| + |Resid| 가 ζ(s)뿐 아니라")
    log("    Dirichlet L-함수에서도 성립.")
    log("    → Proposition을 \"FE+GUE 만족하는 L-함수\"로 일반화 근거 확보.")
    verdict = "POSITIVE"
elif n_path_c_pass > 0:
    log("  ★ 부분 양성:")
    log(f"    {n_path_c_pass}/{n_tested} L-함수에서 Path C 성립.")
    log("    → 일반화에 조건부 근거.")
    verdict = "PARTIAL"
else:
    log("  ★ Path C 미성립:")
    log("    Dirichlet L-함수에서 Path C 부등식 깨짐.")
    log("    → 경계 발견. 어떤 조건에서 깨지는지 분석 필요.")
    verdict = "NEGATIVE"

# 안전 마진 비교
if all_results:
    margins = [r['safety_margin'] for r in all_results]
    log(f"\n  안전 마진 비교:")
    log(f"    ζ(s) T=2000: {zeta_ref['safety_margin']:.2f}×")
    for r in all_results:
        log(f"    {r['label']}: {r['safety_margin']:.2f}×")
    log(f"    최소 마진: {min(margins):.2f}× ({'ζ보다 좁음' if min(margins) < zeta_ref['safety_margin'] else 'ζ와 유사/넓음'})")

log()
log(f"  판정: {verdict}")
elapsed = time.time() - t_global
log(f"\n총 소요: {elapsed:.1f}s")
log("=" * 70)

out_f.close()
print(f"\n결과 저장: {RESULT_PATH}")
