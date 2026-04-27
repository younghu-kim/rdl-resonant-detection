#!/usr/bin/env python3
"""
[C-394] Adjacent Amplitude Correlation — A_Λ Decomposition + lag-structure

C-393 후속: A_bare/d² 대신 A_Λ (Gamma-completed) 사용하여 분해 재현.
목표: 편상관 ρ(A,A|g) +0.118 vs 논문 +0.391 괴리 해소.

A_Λ(t) = (S1_Λ(t))² + 2·H1_Λ(t)
  S1_Λ = S1_bare - Im[(1/2)ψ((s+μ)/2)]
  H1_Λ = H1_bare + Re[(1/4)ψ₁((s+μ)/2)]
  (ζ(s): μ=0, σ_c=1/2)

체크리스트:
  [x] A_Λ = Gamma 보정 포함 (mpmath.digamma, mpmath.psi)
  [x] GUE-정규화 gap: g_GUE = g·d(t)
  [x] d_bar 이론적: log(t/(2π))/(2π)
  [x] trim 20%
  [x] python -u
  [x] Spearman only, p-value 보고
  [x] lag-1~5 자기상관 프로파일 (새 관측)
  [x] 공분산 교차항 세부 분해: Cov(2/g², H1_tail), Cov(2/g², B²)
  [x] mpmath dps=50 (t<2000이므로 충분)
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
T_MAX = 2000.0
TRIM_FRAC = 0.20
W_LIST = [100, 500]
SIGMA_C = 0.5
GAMMA_V = [0]   # ζ(s): Γ_R(s) = π^{-s/2} Γ(s/2)
RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/adj_corr_decomp_ALambda_c394.txt'
)

out_f = open(RESULT_PATH, 'w')
def log(msg=''):
    print(msg, flush=True)
    out_f.write(msg + '\n')
    out_f.flush()


# ── ζ(s) 밀도 ───────────────────────────────────────────────────
def zeta_density(t):
    if t < 2.0:
        return 1.0
    return math.log(t / (2.0 * math.pi)) / (2.0 * math.pi)


# ── Gamma 보정 ──────────────────────────────────────────────────
def gamma_corrections(t0):
    """
    ζ(s)용 Gamma 보정 (gammaV=[0], σ_c=1/2).

    Returns:
      im_gamma: Im[(1/2)ψ((1/2+it₀)/2)] = Im[(1/2)ψ(1/4+it₀/2)]
      re_gamma: Re[(1/4)ψ₁((1/2+it₀)/2)] = Re[(1/4)ψ₁(1/4+it₀/2)]
    """
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


# ── ζ(s) 영점 ────────────────────────────────────────────────────
def get_zeta_zeros():
    log(f"[ζ(s)] 영점 수집 t∈(0, {T_MAX}] ...")
    t0 = time.time()
    pari('L_zeta = lfuncreate(1)')
    pari(f'Li_zeta = lfuninit(L_zeta, [0, {int(T_MAX) + 5}])')
    pari(f'zv_zeta = lfunzeros(Li_zeta, {T_MAX})')
    n = int(str(pari('#zv_zeta')))
    zeros = []
    for i in range(1, n + 1):
        s = str(pari(f'zv_zeta[{i}]')).strip().replace(' E', 'e').replace('E ', 'e')
        try:
            t = float(s)
            if t > 0.5:
                zeros.append(t)
        except ValueError:
            pass
    zeros = np.array(sorted(zeros))
    log(f"  {len(zeros)}개 영점, t∈[{zeros[0]:.3f}, {zeros[-1]:.3f}], {time.time()-t0:.1f}s")
    return zeros


# ── Gamma 보정 사전 계산 ─────────────────────────────────────────
def precompute_gamma_corrections(zeros):
    """모든 영점에 대해 Gamma 보정 사전 계산."""
    log("  Gamma 보정 사전 계산 중 ...")
    t0 = time.time()
    N = len(zeros)
    im_gamma = np.zeros(N)
    re_gamma = np.zeros(N)
    for i in range(N):
        im_g, re_g = gamma_corrections(zeros[i])
        im_gamma[i] = im_g
        re_gamma[i] = re_g
    log(f"  완료: {time.time()-t0:.1f}s")
    return im_gamma, re_gamma


# ══════════════════════════════════════════════════════════════════
# 4-way 분해 (A_Λ 포함)
# ══════════════════════════════════════════════════════════════════
def compute_4way_ALambda(x, W, im_gamma_all, re_gamma_all):
    """
    A_bare 4-way 분해 + Gamma 보정으로 A_Λ 계산.

    A_Λ(n) = (S1(n) - im_γ(n))² + 2·(H1_shared(n) + H1_left(n) + H1_tail(n) + re_γ(n))
    """
    N = len(x)
    if N < 2 * W + 3:
        return None

    lo = W
    hi = N - W - 1
    n_valid = hi - lo
    if n_valid < 30:
        return None

    # 간격
    g = x[lo+1:hi+1] - x[lo:hi]
    g = np.where(np.abs(g) < 1e-15, 1e-15, g)

    delta_L_n = x[lo:hi] - x[lo-1:hi-1]
    delta_L_n = np.where(np.abs(delta_L_n) < 1e-15, 1e-15, delta_L_n)

    delta_R_np1 = x[lo+2:hi+2] - x[lo+1:hi+1]
    delta_R_np1 = np.where(np.abs(delta_R_np1) < 1e-15, 1e-15, delta_R_np1)

    # ── H1 분해 ──
    H1_shared_n = 1.0 / g**2
    H1_shared_np1 = 1.0 / g**2
    H1_left_n = 1.0 / delta_L_n**2
    H1_right_np1 = 1.0 / delta_R_np1**2

    # tail + S1 계산
    H1_tail_n = np.zeros(n_valid)
    H1_tail_np1 = np.zeros(n_valid)
    S1_n = np.zeros(n_valid)
    S1_np1 = np.zeros(n_valid)

    # j=1 기여
    S1_n += -1.0/g + 1.0/delta_L_n
    S1_np1 += -1.0/delta_R_np1 + 1.0/g

    for j in range(2, W + 1):
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

    B_sq_n = S1_n**2
    B_sq_np1 = S1_np1**2

    # ── A_bare ──
    H1_total_n = H1_shared_n + H1_left_n + H1_tail_n
    H1_total_np1 = H1_shared_np1 + H1_right_np1 + H1_tail_np1
    A_bare_n = B_sq_n + 2.0 * H1_total_n
    A_bare_np1 = B_sq_np1 + 2.0 * H1_total_np1

    # ── Gamma 보정 → A_Λ ──
    im_g_n = im_gamma_all[lo:hi]
    im_g_np1 = im_gamma_all[lo+1:hi+1]
    re_g_n = re_gamma_all[lo:hi]
    re_g_np1 = re_gamma_all[lo+1:hi+1]

    S1_Lambda_n = S1_n - im_g_n
    S1_Lambda_np1 = S1_np1 - im_g_np1
    H1_Lambda_n = H1_total_n + re_g_n
    H1_Lambda_np1 = H1_total_np1 + re_g_np1

    A_Lambda_n = S1_Lambda_n**2 + 2.0 * H1_Lambda_n
    A_Lambda_np1 = S1_Lambda_np1**2 + 2.0 * H1_Lambda_np1

    # A_Λ에서 공유항: 여전히 2/g² (Gamma는 공유 아님)
    A_Lambda_prime_n = A_Lambda_n - 2.0 / g**2
    A_Lambda_prime_np1 = A_Lambda_np1 - 2.0 / g**2

    # A_bare에서 공유항 제거 (C-393 대조용)
    A_bare_prime_n = A_bare_n - 2.0 / g**2
    A_bare_prime_np1 = A_bare_np1 - 2.0 / g**2

    t_n = x[lo:hi]
    t_np1 = x[lo+1:hi+1]

    return {
        # A_Λ 시리즈
        'A_Lambda_n': A_Lambda_n, 'A_Lambda_np1': A_Lambda_np1,
        'A_Lambda_prime_n': A_Lambda_prime_n, 'A_Lambda_prime_np1': A_Lambda_prime_np1,
        # A_bare 시리즈 (비교용)
        'A_bare_n': A_bare_n, 'A_bare_np1': A_bare_np1,
        'A_bare_prime_n': A_bare_prime_n, 'A_bare_prime_np1': A_bare_prime_np1,
        # 성분
        'H1_shared_n': H1_shared_n, 'H1_shared_np1': H1_shared_np1,
        'H1_left_n': H1_left_n, 'H1_right_np1': H1_right_np1,
        'H1_tail_n': H1_tail_n, 'H1_tail_np1': H1_tail_np1,
        'B_sq_n': B_sq_n, 'B_sq_np1': B_sq_np1,
        'S1_Lambda_n': S1_Lambda_n, 'S1_Lambda_np1': S1_Lambda_np1,
        're_gamma_n': re_g_n, 're_gamma_np1': re_g_np1,
        'g': g,
        'delta_L_n': delta_L_n, 'delta_R_np1': delta_R_np1,
        't_n': t_n, 't_np1': t_np1,
        'n_valid': n_valid,
    }


def apply_trim(data):
    """trim 적용 (밀도 정규화는 A_Λ에서 불필요)."""
    n_pts = data['n_valid']
    trim_lo = int(n_pts * TRIM_FRAC)
    trim_hi = n_pts - int(n_pts * TRIM_FRAC)
    if trim_hi - trim_lo < 30:
        return None

    sl = slice(trim_lo, trim_hi)
    trimmed = {}
    for key in data:
        if isinstance(data[key], np.ndarray) and len(data[key]) == n_pts:
            trimmed[key] = data[key][sl]
        else:
            trimmed[key] = data[key]

    # NaN/Inf/음수 필터
    mask = np.ones(trimmed['A_Lambda_n'].shape[0], dtype=bool)
    for key in ['A_Lambda_n', 'A_Lambda_np1', 'A_bare_n', 'A_bare_np1', 'g']:
        if key in trimmed:
            mask &= np.isfinite(trimmed[key])
    mask &= (trimmed['A_Lambda_n'] > 0) & (trimmed['A_Lambda_np1'] > 0)
    mask &= (np.abs(trimmed['g']) > 1e-15)

    filtered = {}
    for key in trimmed:
        if isinstance(trimmed[key], np.ndarray) and len(trimmed[key]) == len(mask):
            filtered[key] = trimmed[key][mask]
        else:
            filtered[key] = trimmed[key]
    filtered['n_valid'] = int(np.sum(mask))
    return filtered


# ══════════════════════════════════════════════════════════════════
# 분석
# ══════════════════════════════════════════════════════════════════
def safe_spearmanr(a, b):
    mask = np.isfinite(a) & np.isfinite(b)
    if np.sum(mask) < 10:
        return float('nan'), float('nan')
    return stats.spearmanr(a[mask], b[mask])


def partial_corr_spearman(x, y, z):
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if np.sum(mask) < 10:
        return float('nan'), float('nan')
    x, y, z = x[mask], y[mask], z[mask]
    rho_xy, _ = stats.spearmanr(x, y)
    rho_xz, _ = stats.spearmanr(x, z)
    rho_yz, _ = stats.spearmanr(y, z)
    denom = np.sqrt(max(1e-15, (1 - rho_xz**2) * (1 - rho_yz**2)))
    rho_partial = (rho_xy - rho_xz * rho_yz) / denom
    n = np.sum(mask)
    z_val = 0.5 * np.log((1 + rho_partial) / max(1e-15, 1 - rho_partial))
    se = 1.0 / np.sqrt(max(1, n - 4))
    p_val = 2.0 * (1.0 - stats.norm.cdf(abs(z_val) / se))
    return rho_partial, p_val


def analyze_all(data, W):
    """A_Λ + A_bare 비교 분석."""
    n = data['n_valid']
    if n < 30:
        log(f"  W={W}: 데이터 부족 ({n}개)")
        return None

    # ── A_Λ 분석 ──
    AL_n = data['A_Lambda_n']
    AL_np1 = data['A_Lambda_np1']
    ALp_n = data['A_Lambda_prime_n']
    ALp_np1 = data['A_Lambda_prime_np1']
    g = data['g']

    # A_bare (밀도 정규화 적용)
    d_n = np.array([zeta_density(t) for t in data['t_n']])
    d_np1 = np.array([zeta_density(t) for t in data['t_np1']])
    d_sq_n = np.maximum(d_n**2, 1e-30)
    d_sq_np1 = np.maximum(d_np1**2, 1e-30)
    AB_norm_n = data['A_bare_n'] / d_sq_n
    AB_norm_np1 = data['A_bare_np1'] / d_sq_np1

    # GUE-정규화 gap
    d_mid = (d_n + d_np1) / 2.0
    g_gue = g * d_mid

    log(f"\n  W={W}: n={n} 쌍")
    log("=" * 70)

    # ── Part A: A_Λ 상관 ──
    log("\n  [Part A] A_Λ 기반 분석")
    rho_full_L, p_full_L = safe_spearmanr(AL_n, AL_np1)
    rho_resid_L, p_resid_L = safe_spearmanr(ALp_n, ALp_np1)
    log(f"    ρ_full(A_Λ_n, A_Λ_{{n+1}})       = {rho_full_L:+.4f}  (p={p_full_L:.2e})")
    log(f"    ρ_residual(A'_Λ_n, A'_Λ_{{n+1}}) = {rho_resid_L:+.4f}  (p={p_resid_L:.2e})")
    log(f"    → Δρ = {rho_full_L - rho_resid_L:+.4f}")

    # 편상관 (raw gap)
    rho_part_L, p_part_L = partial_corr_spearman(AL_n, AL_np1, g)
    log(f"    ρ_partial(A_Λ | g)               = {rho_part_L:+.4f}  (p={p_part_L:.2e})")

    # 편상관 (GUE-정규화 gap)
    rho_part_gue, p_part_gue = partial_corr_spearman(AL_n, AL_np1, g_gue)
    log(f"    ρ_partial(A_Λ | g_GUE)           = {rho_part_gue:+.4f}  (p={p_part_gue:.2e})")
    log(f"    (논문 ρ_partial ≈ +0.391 참조)")

    # Pearson Cov 분해 (A_Λ)
    shared_term = 2.0 / g**2
    cov_full_L = np.cov(AL_n, AL_np1)[0, 1]
    cov_shared_L = np.var(shared_term)
    cov_resid_L = np.cov(ALp_n, ALp_np1)[0, 1]
    cov_cross_L = cov_full_L - cov_shared_L - cov_resid_L

    log(f"\n    [Pearson 공분산 분해 — A_Λ]")
    log(f"    Cov(A_Λ)      = {cov_full_L:.4e}")
    log(f"    Var(shared)   = {cov_shared_L:.4e}  ({cov_shared_L/max(abs(cov_full_L),1e-30)*100:.1f}%)")
    log(f"    Cov(residual) = {cov_resid_L:.4e}  ({cov_resid_L/max(abs(cov_full_L),1e-30)*100:.1f}%)")
    log(f"    Cross         = {cov_cross_L:.4e}  ({cov_cross_L/max(abs(cov_full_L),1e-30)*100:.1f}%)")

    # ── Part B: A_bare/d² 상관 (C-393 재현) ──
    log(f"\n  [Part B] A_bare/d² 기반 (C-393 재현)")
    ABp_n = AB_norm_n - 2.0/(g * d_mid)**2  # 밀도 정규화 후의 공유항 제거
    ABp_np1 = AB_norm_np1 - 2.0/(g * d_mid)**2
    rho_full_B, p_full_B = safe_spearmanr(AB_norm_n, AB_norm_np1)
    rho_resid_B, p_resid_B = safe_spearmanr(ABp_n, ABp_np1)
    log(f"    ρ_full(A_bare/d²)  = {rho_full_B:+.4f}  (p={p_full_B:.2e})")
    log(f"    ρ_residual         = {rho_resid_B:+.4f}  (p={p_resid_B:.2e})")

    rho_part_B, p_part_B = partial_corr_spearman(AB_norm_n, AB_norm_np1, g)
    log(f"    ρ_partial(|g)      = {rho_part_B:+.4f}  (p={p_part_B:.2e})")

    # ── Part C: 잔차 원천 분해 (A_Λ 기준) ──
    log(f"\n  [Part C] 잔차 원천 분해 (A_Λ)")
    H1_tail_n = data['H1_tail_n']
    H1_tail_np1 = data['H1_tail_np1']
    B_sq_n = data['B_sq_n']
    B_sq_np1 = data['B_sq_np1']
    # A_Λ에서의 B²는 (S1_Λ)²
    BL_sq_n = data['S1_Lambda_n']**2
    BL_sq_np1 = data['S1_Lambda_np1']**2
    # Gamma 보정 re_γ
    re_g_n = data['re_gamma_n']
    re_g_np1 = data['re_gamma_np1']
    # H1_Λ_tail = H1_tail + re_γ (Gamma는 tail에 합산)
    H1L_tail_n = H1_tail_n + re_g_n
    H1L_tail_np1 = H1_tail_np1 + re_g_np1

    rho_tail_L, p_tail_L = safe_spearmanr(H1L_tail_n, H1L_tail_np1)
    rho_BLsq, p_BLsq = safe_spearmanr(BL_sq_n, BL_sq_np1)
    rho_tail_bare, p_tail_bare = safe_spearmanr(H1_tail_n, H1_tail_np1)
    rho_Bsq_bare, p_Bsq_bare = safe_spearmanr(B_sq_n, B_sq_np1)

    log(f"    ρ(H1_Λ_tail_n, H1_Λ_tail_{{n+1}})  = {rho_tail_L:+.4f}  (p={p_tail_L:.2e})")
    log(f"    ρ((S1_Λ)²_n, (S1_Λ)²_{{n+1}})      = {rho_BLsq:+.4f}  (p={p_BLsq:.2e})")
    log(f"    ρ(H1_tail_n, H1_tail_{{n+1}}) [bare]  = {rho_tail_bare:+.4f}  (p={p_tail_bare:.2e})")
    log(f"    ρ(B²_n, B²_{{n+1}}) [bare]            = {rho_Bsq_bare:+.4f}  (p={p_Bsq_bare:.2e})")

    # ── Part D: 교차항 세부 분해 ──
    log(f"\n  [Part D] 교차항 세부: Cov(2/g², 각 성분)")
    cov_sh_tail = np.cov(shared_term, H1L_tail_n)[0,1] + np.cov(shared_term, H1L_tail_np1)[0,1]
    cov_sh_Bsq = np.cov(shared_term, BL_sq_n)[0,1] + np.cov(shared_term, BL_sq_np1)[0,1]
    log(f"    Cov(2/g², H1_Λ_tail)  = {cov_sh_tail:.4e}")
    log(f"    Cov(2/g², (S1_Λ)²)    = {cov_sh_Bsq:.4e}")
    total_cross = cov_sh_tail + cov_sh_Bsq
    if abs(total_cross) > 1e-30:
        log(f"    H1_Λ_tail 비율: {cov_sh_tail/total_cross*100:.1f}%")
        log(f"    (S1_Λ)² 비율:  {cov_sh_Bsq/total_cross*100:.1f}%")

    # ── Part E: lag-구조 프로파일 (A_Λ) ──
    log(f"\n  [Part E] lag-구조 프로파일 (A_Λ)")
    # AL_n[i] = A_Λ at zero index (lo+i), use single array for autocorrelation
    lag_max = min(5, n // 10)
    for lag in range(1, lag_max + 1):
        if lag < n:
            rho_lag, p_lag = safe_spearmanr(AL_n[:-lag], AL_n[lag:])
            log(f"    lag-{lag}: ρ = {rho_lag:+.4f}  (p={p_lag:.2e})")

    # ── Part F: Gamma 보정 크기 통계 ──
    log(f"\n  [Part F] Gamma 보정 크기")
    log(f"    |re_γ|: mean={np.mean(np.abs(re_g_n)):.4f}, max={np.max(np.abs(re_g_n)):.4f}")
    log(f"    |im_γ|: mean={np.mean(np.abs(data.get('im_gamma_n', np.zeros(1)))):.4f}")
    log(f"    H1_bare mean: {np.mean(data['H1_tail_n'] + data['H1_shared_n'] + data['H1_left_n']):.2f}")
    log(f"    re_γ / H1_bare: {np.mean(np.abs(re_g_n))/np.mean(data['H1_tail_n'] + data['H1_shared_n'] + data['H1_left_n'])*100:.1f}%")

    return {
        'n': n, 'W': W,
        'rho_full_L': rho_full_L, 'rho_resid_L': rho_resid_L,
        'rho_part_L': rho_part_L, 'rho_part_gue': rho_part_gue,
        'rho_full_B': rho_full_B, 'rho_resid_B': rho_resid_B,
        'rho_part_B': rho_part_B,
        'cov_shared_pct': cov_shared_L/max(abs(cov_full_L),1e-30),
        'cov_resid_pct': cov_resid_L/max(abs(cov_full_L),1e-30),
        'cov_cross_pct': cov_cross_L/max(abs(cov_full_L),1e-30),
    }


# ══════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════
log("=" * 70)
log("[C-394] Adjacent Amplitude Correlation — A_Λ Decomposition")
log("=" * 70)
log(f"  T_MAX={T_MAX}, TRIM={TRIM_FRAC}, W_LIST={W_LIST}")
log(f"  Gamma: μ={GAMMA_V}, σ_c={SIGMA_C}")
log()

t_start = time.time()

zeros_zeta = get_zeta_zeros()
im_gamma_all, re_gamma_all = precompute_gamma_corrections(zeros_zeta)

log()
log("=" * 70)
log("  [Phase 1] A_Λ vs A_bare/d² 분해 비교")
log("=" * 70)

all_results = {}
for W in W_LIST:
    data = compute_4way_ALambda(zeros_zeta, W, im_gamma_all, re_gamma_all)
    if data is None:
        log(f"  W={W}: 데이터 부족")
        continue
    # im_gamma도 저장 (Part F용)
    processed = apply_trim(data)
    if processed is None:
        log(f"  W={W}: trim 후 부족")
        continue
    r = analyze_all(processed, W)
    if r:
        all_results[W] = r

# ── 비교 요약표 ──
if len(all_results) >= 1:
    log()
    log("=" * 70)
    log("  [Phase 2] A_Λ vs A_bare/d² 비교 요약")
    log("=" * 70)
    log()
    log(f"{'W':>5} │ {'ρ_full_Λ':>9} {'ρ_res_Λ':>9} {'ρ_part_Λ':>9} {'ρ_part_GUE':>11} │"
        f" {'ρ_full_B':>9} {'ρ_res_B':>9} {'ρ_part_B':>9}")
    log("─" * 95)
    for W in sorted(all_results.keys()):
        r = all_results[W]
        log(f"{W:5d} │ {r['rho_full_L']:+.4f}   {r['rho_resid_L']:+.4f}   "
            f"{r['rho_part_L']:+.4f}   {r['rho_part_gue']:+.4f}     │"
            f" {r['rho_full_B']:+.4f}   {r['rho_resid_B']:+.4f}   {r['rho_part_B']:+.4f}")

# ── 판정 ──
log()
log("=" * 70)
log("  [Phase 3] 판정")
log("=" * 70)

if all_results:
    W_ref = max(all_results.keys())
    r = all_results[W_ref]
    log(f"\n  참조: W={W_ref}")
    log(f"  A_Λ:      ρ_full={r['rho_full_L']:+.4f}, ρ_partial={r['rho_part_L']:+.4f}, ρ_partial_GUE={r['rho_part_gue']:+.4f}")
    log(f"  A_bare/d²: ρ_full={r['rho_full_B']:+.4f}, ρ_partial={r['rho_part_B']:+.4f}")
    log()

    # 논문 +0.391 재현 판정
    paper_val = 0.391
    best_partial = max(r['rho_part_L'], r['rho_part_gue'])
    gap = abs(best_partial - paper_val)

    if gap < 0.05:
        log(f"  ★★★★★ 논문 ρ_partial ≈ +0.391 재현 성공 (차이 {gap:.3f})")
    elif gap < 0.15:
        log(f"  ★★★★ 부분 재현 (차이 {gap:.3f}) — 범위/방법론 차이 가능")
    else:
        log(f"  ★★★ 괴리 지속 (차이 {gap:.3f}) — A_Λ vs A_bare 외 다른 원인 탐색 필요")

    # Cov 분해 비교
    log(f"\n  Cov 분해 (A_Λ): shared={r['cov_shared_pct']*100:.1f}%, "
        f"residual={r['cov_resid_pct']*100:.1f}%, cross={r['cov_cross_pct']*100:.1f}%")

log()
elapsed = time.time() - t_start
log(f"총 소요: {elapsed:.1f}s")
log("=" * 70)

out_f.close()
print(f"\n결과 저장: {RESULT_PATH}")
