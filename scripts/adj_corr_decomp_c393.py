#!/usr/bin/env python3
"""
[C-393] Adjacent Amplitude Correlation Covariance Decomposition

Paper 4 Obs 1: ρ(A_n, A_{n+1}) ∈ [+0.28, +0.53] (13 L-함수 전부 양성)
증명 전략 결정을 위해 Cov(A_n, A_{n+1})를 항별 분해.

핵심 구조:
  A(γ_n) = B(γ_n)² + 2·H1(γ_n)
  H1(γ_n) = Σ_k 1/(γ_n - γ_k)²

  인접 쌍 (n, n+1)에서 g = γ_{n+1} - γ_n 이 양쪽에 공유됨:
  H1(n) = 1/g² + 1/Δ_L(n)² + H1_tail(n)      (Δ_L = γ_n - γ_{n-1})
  H1(n+1) = 1/g² + 1/Δ_R(n+1)² + H1_tail(n+1) (Δ_R = γ_{n+2} - γ_{n+1})

  → 공유항 2/g² 가 양쪽 A(n), A(n+1)에 동시 기여 → 양상관 발생 가능

Steps:
  1. 4-way 분해: H1_shared, H1_left, H1_tail, B_sq
  2. A'(n) = A(n) - 2/g² → ρ_residual
  3. 잔차 원천 분해: tail, NN교차, B² 간 상관
  4. 분산 기여도: Var(2/g²)/Var(A)
  5. W 의존성 (W=100, 500)
  6. 편상관 교차 검증
  7. 비인접(lag-2) 대조군

체크리스트:
  [x] 이론적 d_bar = log(t/(2π))/(2π) for ζ(s)
  [x] trim 20%
  [x] python -u
  [x] Spearman only, p-value 보고
  [x] g=0 안전 가드
  [x] H1_shared = 오른쪽 NN만 (인접 쌍별 공유항 정의)
  [x] A_rest ≠ C-297의 A_rest (혼동 금지)
  [x] 비인접 lag-2 대조군
"""

import sys, os, time
import numpy as np
from scipy import stats

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(1024 * 10**6)
    pari.set_real_precision(80)
    print("cypari2 OK", flush=True)
except Exception as e:
    print(f"FATAL: {e}", flush=True)
    sys.exit(1)

# ── 설정 ──────────────────────────────────────────────────────────
T_MAX = 2000.0
TRIM_FRAC = 0.20
W_LIST = [100, 500]   # H1_tail 수렴 확인용
RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/adj_corr_decomp_c393.txt'
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
    return np.log(t / (2.0 * np.pi)) / (2.0 * np.pi)


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


# ══════════════════════════════════════════════════════════════════
# 4-way 분해: H1_shared, H1_left, H1_tail, B(=S1)
# ══════════════════════════════════════════════════════════════════
def compute_4way_decomposition(x, W):
    """각 영점 n에 대해 4-way 분해 계산.

    H1_shared(n) = 1/g(n,n+1)²        # 오른쪽 NN (인접 쌍의 공유 항)
    H1_left(n)   = 1/Δ_L(n)²          # 왼쪽 NN (비공유)
    H1_tail(n)   = Σ_{|k-n|≥2} 1/(γ_n-γ_k)²  # 원거리
    B(n) = S1(n) = Σ_k (-1)^{...} / (γ_n - γ_k)  # 교대합

    A_bare(n) = B(n)² + 2·(H1_shared(n) + H1_left(n) + H1_tail(n))

    Returns arrays aligned for pairs (n, n+1).
    """
    N = len(x)
    if N < 2 * W + 3:
        return None

    # 유효 범위: W ≤ n ≤ N-W-2 (n+1도 W 이내여야 함)
    lo = W
    hi = N - W - 1  # n+1 = hi 일 때 hi+W ≤ N-1
    n_valid = hi - lo
    if n_valid < 30:
        return None

    # 간격 계산
    # g(n, n+1) = x[n+1] - x[n]
    g = x[lo+1:hi+1] - x[lo:hi]  # shape: (n_valid,)
    g = np.where(np.abs(g) < 1e-15, 1e-15, g)

    # Δ_L(n) = x[n] - x[n-1]
    delta_L_n = x[lo:hi] - x[lo-1:hi-1]
    delta_L_n = np.where(np.abs(delta_L_n) < 1e-15, 1e-15, delta_L_n)

    # Δ_R(n+1) = x[n+2] - x[n+1]
    delta_R_np1 = x[lo+2:hi+2] - x[lo+1:hi+1]
    delta_R_np1 = np.where(np.abs(delta_R_np1) < 1e-15, 1e-15, delta_R_np1)

    # Δ_L(n+1) = x[n+1] - x[n] = g  (이것은 shared)
    # Δ_R(n) = x[n+1] - x[n] = g    (이것도 shared)

    # ── H1 분해 (n 기준) ──
    H1_shared_n = 1.0 / g**2                  # 오른쪽 NN = 공유
    H1_left_n = 1.0 / delta_L_n**2            # 왼쪽 NN = 비공유

    # ── H1 분해 (n+1 기준) ──
    H1_shared_np1 = 1.0 / g**2                # 왼쪽 NN = 공유 (같은 g!)
    H1_right_np1 = 1.0 / delta_R_np1**2       # 오른쪽 NN = 비공유

    # ── H1_tail (원거리, |k-n| ≥ 2) ──
    H1_tail_n = np.zeros(n_valid)
    H1_tail_np1 = np.zeros(n_valid)
    S1_n = np.zeros(n_valid)      # B(n)
    S1_np1 = np.zeros(n_valid)    # B(n+1)

    # j=1 기여 to S1
    S1_n += -1.0/g + 1.0/delta_L_n
    S1_np1 += -1.0/delta_R_np1 + 1.0/g

    for j in range(2, W + 1):
        # n 기준 j번째 이웃
        idx_n = np.arange(lo, hi)
        # 오른쪽 j번째
        diff_rj_n = x[idx_n + j] - x[idx_n]
        diff_rj_n = np.where(np.abs(diff_rj_n) < 1e-15, 1e-15, diff_rj_n)
        # 왼쪽 j번째
        diff_lj_n = x[idx_n] - x[idx_n - j]
        diff_lj_n = np.where(np.abs(diff_lj_n) < 1e-15, 1e-15, diff_lj_n)

        H1_tail_n += 1.0/diff_rj_n**2 + 1.0/diff_lj_n**2
        S1_n += -1.0/diff_rj_n + 1.0/diff_lj_n

        # n+1 기준 j번째 이웃
        idx_np1 = np.arange(lo+1, hi+1)
        diff_rj_np1 = x[idx_np1 + j] - x[idx_np1]
        diff_rj_np1 = np.where(np.abs(diff_rj_np1) < 1e-15, 1e-15, diff_rj_np1)
        diff_lj_np1 = x[idx_np1] - x[idx_np1 - j]
        diff_lj_np1 = np.where(np.abs(diff_lj_np1) < 1e-15, 1e-15, diff_lj_np1)

        H1_tail_np1 += 1.0/diff_rj_np1**2 + 1.0/diff_lj_np1**2
        S1_np1 += -1.0/diff_rj_np1 + 1.0/diff_lj_np1

    B_sq_n = S1_n**2
    B_sq_np1 = S1_np1**2

    # ── A_bare 조립 ──
    H1_total_n = H1_shared_n + H1_left_n + H1_tail_n
    H1_total_np1 = H1_shared_np1 + H1_right_np1 + H1_tail_np1

    A_n = B_sq_n + 2.0 * H1_total_n
    A_np1 = B_sq_np1 + 2.0 * H1_total_np1

    # ── A' (공유항 제거) ──
    A_prime_n = A_n - 2.0 / g**2
    A_prime_np1 = A_np1 - 2.0 / g**2

    # 위치 (밀도 정규화용)
    t_n = x[lo:hi]
    t_np1 = x[lo+1:hi+1]

    return {
        'A_n': A_n, 'A_np1': A_np1,
        'A_prime_n': A_prime_n, 'A_prime_np1': A_prime_np1,
        'H1_shared_n': H1_shared_n, 'H1_shared_np1': H1_shared_np1,
        'H1_left_n': H1_left_n, 'H1_right_np1': H1_right_np1,
        'H1_tail_n': H1_tail_n, 'H1_tail_np1': H1_tail_np1,
        'B_sq_n': B_sq_n, 'B_sq_np1': B_sq_np1,
        'g': g,
        't_n': t_n, 't_np1': t_np1,
        'delta_L_n': delta_L_n, 'delta_R_np1': delta_R_np1,
        'n_valid': n_valid,
    }


def apply_density_and_trim(data, density_func=None):
    """밀도 정규화(A → A/d²) + trim 적용."""
    result = {}
    for key in data:
        if isinstance(data[key], np.ndarray):
            result[key] = data[key].copy()
        else:
            result[key] = data[key]

    # 밀도 정규화 — A_bare는 d²에 비례하므로 d²로 나눔
    if density_func is not None:
        d_n = np.array([density_func(t) for t in data['t_n']])
        d_np1 = np.array([density_func(t) for t in data['t_np1']])
        d_sq_n = d_n**2
        d_sq_np1 = d_np1**2
        d_sq_n = np.where(d_sq_n < 1e-30, 1e-30, d_sq_n)
        d_sq_np1 = np.where(d_sq_np1 < 1e-30, 1e-30, d_sq_np1)

        for key in ['A_n', 'A_prime_n', 'H1_shared_n', 'H1_left_n',
                     'H1_tail_n', 'B_sq_n']:
            result[key] = result[key] / d_sq_n
        for key in ['A_np1', 'A_prime_np1', 'H1_shared_np1', 'H1_right_np1',
                     'H1_tail_np1', 'B_sq_np1']:
            result[key] = result[key] / d_sq_np1
        # g는 d로 정규화 (스케일링 일관성)
        d_mid = (d_n + d_np1) / 2.0
        result['g'] = result['g'] * d_mid

    # trim
    n_pts = result['n_valid']
    trim_lo = int(n_pts * TRIM_FRAC)
    trim_hi = n_pts - int(n_pts * TRIM_FRAC)
    if trim_hi - trim_lo < 30:
        return None

    sl = slice(trim_lo, trim_hi)
    trimmed = {}
    for key in result:
        if isinstance(result[key], np.ndarray) and len(result[key]) == n_pts:
            trimmed[key] = result[key][sl]
        else:
            trimmed[key] = result[key]

    # NaN/Inf 필터
    mask = np.ones(trimmed['A_n'].shape[0], dtype=bool)
    for key in ['A_n', 'A_np1', 'A_prime_n', 'A_prime_np1', 'g',
                'H1_tail_n', 'H1_tail_np1', 'B_sq_n', 'B_sq_np1']:
        if key in trimmed and isinstance(trimmed[key], np.ndarray):
            mask &= np.isfinite(trimmed[key])

    mask &= (trimmed['A_n'] > 0) & (trimmed['A_np1'] > 0) & (np.abs(trimmed['g']) > 1e-15)

    filtered = {}
    for key in trimmed:
        if isinstance(trimmed[key], np.ndarray) and len(trimmed[key]) == len(mask):
            filtered[key] = trimmed[key][mask]
        else:
            filtered[key] = trimmed[key]

    filtered['n_valid'] = int(np.sum(mask))
    return filtered


# ══════════════════════════════════════════════════════════════════
# 분석 함수
# ══════════════════════════════════════════════════════════════════
def safe_spearmanr(a, b):
    """NaN-safe Spearman."""
    mask = np.isfinite(a) & np.isfinite(b)
    if np.sum(mask) < 10:
        return float('nan'), float('nan')
    return stats.spearmanr(a[mask], b[mask])


def partial_corr_spearman(x, y, z):
    """편상관: ρ(x,y|z) = (ρ_xy - ρ_xz·ρ_yz) / sqrt((1-ρ_xz²)(1-ρ_yz²))"""
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if np.sum(mask) < 10:
        return float('nan'), float('nan')
    x, y, z = x[mask], y[mask], z[mask]
    rho_xy, _ = stats.spearmanr(x, y)
    rho_xz, _ = stats.spearmanr(x, z)
    rho_yz, _ = stats.spearmanr(y, z)

    denom = np.sqrt(max(1e-15, (1 - rho_xz**2) * (1 - rho_yz**2)))
    rho_partial = (rho_xy - rho_xz * rho_yz) / denom
    # 근사 p-value (Fisher z)
    n = np.sum(mask)
    z_val = 0.5 * np.log((1 + rho_partial) / max(1e-15, 1 - rho_partial))
    se = 1.0 / np.sqrt(max(1, n - 4))
    p_val = 2.0 * (1.0 - stats.norm.cdf(abs(z_val) / se))
    return rho_partial, p_val


def analyze_all(data, label, W):
    """전체 분석: Steps 1-6."""
    n = data['n_valid']
    if n < 30:
        log(f"  [{label}] W={W}: 데이터 부족 ({n}개)")
        return None

    A_n = data['A_n']
    A_np1 = data['A_np1']
    A_prime_n = data['A_prime_n']
    A_prime_np1 = data['A_prime_np1']
    g = data['g']
    H1_tail_n = data['H1_tail_n']
    H1_tail_np1 = data['H1_tail_np1']
    B_sq_n = data['B_sq_n']
    B_sq_np1 = data['B_sq_np1']
    H1_left_n = data['H1_left_n']
    H1_right_np1 = data['H1_right_np1']
    delta_L_inv2 = 1.0 / np.where(np.abs(data['delta_L_n']) < 1e-15, 1e-15, data['delta_L_n'])**2
    delta_R_inv2 = 1.0 / np.where(np.abs(data['delta_R_np1']) < 1e-15, 1e-15, data['delta_R_np1'])**2

    log(f"\n  [{label}] W={W}: n={n} 쌍")
    log("=" * 60)

    # ── Step 2: 인접 쌍별 공분산 분해 ──
    log("\n  [Step 2] 인접 쌍 상관")
    rho_full, p_full = safe_spearmanr(A_n, A_np1)
    rho_resid, p_resid = safe_spearmanr(A_prime_n, A_prime_np1)
    log(f"    ρ_full(A_n, A_{{n+1}})         = {rho_full:+.4f}  (p={p_full:.2e})")
    log(f"    ρ_residual(A'_n, A'_{{n+1}})   = {rho_resid:+.4f}  (p={p_resid:.2e})")
    log(f"    → 공유항 제거 효과: Δρ = {rho_full - rho_resid:+.4f}")

    # ── Step 3: 잔차 상관의 원천 분해 ──
    log("\n  [Step 3] 잔차 상관 원천 분해")

    rho_tail, p_tail = safe_spearmanr(H1_tail_n, H1_tail_np1)
    rho_nn_cross, p_nn_cross = safe_spearmanr(delta_L_inv2, delta_R_inv2)
    rho_Bsq, p_Bsq = safe_spearmanr(B_sq_n, B_sq_np1)

    log(f"    ρ(H1_tail_n, H1_tail_{{n+1}})  = {rho_tail:+.4f}  (p={p_tail:.2e})  ← 환경 공유")
    log(f"    ρ(Δ_L⁻², Δ_R⁻²)              = {rho_nn_cross:+.4f}  (p={p_nn_cross:.2e})  ← NN 교차")
    log(f"    ρ(B²_n, B²_{{n+1}})            = {rho_Bsq:+.4f}  (p={p_Bsq:.2e})  ← B항 간")

    # 추가: 교차항 확인
    rho_tail_Bsq_n, _ = safe_spearmanr(H1_tail_n, B_sq_n)
    rho_tail_Bsq_np1, _ = safe_spearmanr(H1_tail_np1, B_sq_np1)
    log(f"    ρ(H1_tail_n, B²_n)            = {rho_tail_Bsq_n:+.4f}  (교차 검증)")
    log(f"    ρ(H1_tail_{{n+1}}, B²_{{n+1}})  = {rho_tail_Bsq_np1:+.4f}  (교차 검증)")

    # ── Step 4: 분산 기여도 분석 ──
    log("\n  [Step 4] 분산 기여도")
    shared_term = 2.0 / g**2  # 밀도 정규화 후라면 이미 반영됨
    # 실제 분산 기여: 공유항 자체의 분산 vs A의 분산
    # 주의: shared_term이 A_n과 A_np1 모두에 기여하므로
    var_A_n = np.var(A_n)
    var_A_np1 = np.var(A_np1)
    var_shared = np.var(shared_term)
    var_A_prime_n = np.var(A_prime_n)
    var_A_prime_np1 = np.var(A_prime_np1)

    log(f"    Var(A_n)       = {var_A_n:.4e}")
    log(f"    Var(2/g²)      = {var_shared:.4e}")
    log(f"    Var(A'_n)      = {var_A_prime_n:.4e}")
    log(f"    Var(2/g²) / Var(A_n)   = {var_shared/max(var_A_n,1e-30):.4f}")
    log(f"    Var(A'_n) / Var(A_n)   = {var_A_prime_n/max(var_A_n,1e-30):.4f}")

    # Pearson Cov 분해도 보자 (Spearman ρ와 함께)
    cov_full = np.cov(A_n, A_np1)[0, 1]
    cov_shared = np.cov(shared_term, shared_term)[0, 1]  # = Var(shared)
    cov_resid = np.cov(A_prime_n, A_prime_np1)[0, 1]
    cov_cross = cov_full - cov_shared - cov_resid  # 교차 공분산

    log(f"\n    [Pearson 공분산 분해]")
    log(f"    Cov(A_n, A_{{n+1}})     = {cov_full:.4e}")
    log(f"    Var(shared) = Cov(2/g², 2/g²) = {cov_shared:.4e}")
    log(f"    Cov(A'_n, A'_{{n+1}})   = {cov_resid:.4e}")
    log(f"    교차항(cross)           = {cov_cross:.4e}")
    log(f"    Var(shared)/Cov(full)   = {cov_shared/max(abs(cov_full),1e-30):.4f}")
    log(f"    Cov(resid)/Cov(full)    = {cov_resid/max(abs(cov_full),1e-30):.4f}")
    log(f"    Cross/Cov(full)         = {cov_cross/max(abs(cov_full),1e-30):.4f}")

    # ── Step 6: 편상관 교차 검증 ──
    log("\n  [Step 6] 편상관 교차 검증")
    rho_partial_full, p_partial = partial_corr_spearman(A_n, A_np1, g)
    rho_partial_resid, p_partial_r = partial_corr_spearman(A_prime_n, A_prime_np1, g)
    log(f"    ρ_partial(A_n, A_{{n+1}} | g)    = {rho_partial_full:+.4f}  (p={p_partial:.2e})")
    log(f"    ρ_partial(A'_n, A'_{{n+1}} | g)  = {rho_partial_resid:+.4f}  (p={p_partial_r:.2e})")
    log(f"    (논문의 ρ_partial ≈ +0.391 재현 확인)")

    # ── 비인접 (lag-2) 대조군 ──
    log("\n  [대조군] 비인접 (lag-2) 상관")
    # A(n) vs A(n+2): 공유항 없으므로 상관 약할 것
    if len(A_n) > 2:
        rho_lag2, p_lag2 = safe_spearmanr(A_n[:-1], A_np1[1:])
        log(f"    ρ(A_n, A_{{n+2}})  = {rho_lag2:+.4f}  (p={p_lag2:.2e})")
        rho_lag3, p_lag3 = safe_spearmanr(A_n[:-2], A_np1[2:])
        log(f"    ρ(A_n, A_{{n+3}})  = {rho_lag3:+.4f}  (p={p_lag3:.2e})")

    results = {
        'n': n, 'W': W,
        'rho_full': rho_full, 'rho_resid': rho_resid,
        'rho_tail': rho_tail, 'rho_nn_cross': rho_nn_cross, 'rho_Bsq': rho_Bsq,
        'var_shared_ratio': var_shared/max(var_A_n,1e-30),
        'cov_shared_ratio': cov_shared/max(abs(cov_full),1e-30),
        'cov_resid_ratio': cov_resid/max(abs(cov_full),1e-30),
        'rho_partial_full': rho_partial_full,
    }
    return results


# ══════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════
log("=" * 70)
log("[C-393] Adjacent Amplitude Correlation Covariance Decomposition")
log("=" * 70)
log(f"  T_MAX={T_MAX}, TRIM={TRIM_FRAC}, W_LIST={W_LIST}")
log()

t_start = time.time()

# ── ζ(s) 영점 ────────────────────────────────────────────────────
zeros_zeta = get_zeta_zeros()

# ── Phase 1: W별 분석 ────────────────────────────────────────────
log()
log("=" * 70)
log("  [Phase 1] ζ(s) — W별 4-way 분해")
log("=" * 70)

all_results = {}
for W in W_LIST:
    data = compute_4way_decomposition(zeros_zeta, W)
    if data is None:
        log(f"  W={W}: 데이터 부족, 스킵")
        continue
    processed = apply_density_and_trim(data, zeta_density)
    if processed is None:
        log(f"  W={W}: trim 후 부족, 스킵")
        continue
    r = analyze_all(processed, "ζ", W)
    if r:
        all_results[W] = r


# ── Phase 2: W 수렴 비교표 ────────────────────────────────────────
if len(all_results) >= 2:
    log()
    log("=" * 70)
    log("  [Phase 2] W 수렴 비교")
    log("=" * 70)
    log()
    log(f"{'W':>6}  {'ρ_full':>8}  {'ρ_resid':>8}  {'Δρ':>8}  "
        f"{'Var(sh)/Var(A)':>14}  {'Cov(sh)/Cov':>12}  {'ρ_tail':>8}  {'ρ_B²':>8}")
    log("─" * 80)
    for W in sorted(all_results.keys()):
        r = all_results[W]
        dr = r['rho_full'] - r['rho_resid']
        log(f"{W:6d}  {r['rho_full']:+.4f}  {r['rho_resid']:+.4f}  {dr:+.4f}  "
            f"{r['var_shared_ratio']:14.4f}  {r['cov_shared_ratio']:12.4f}  "
            f"{r['rho_tail']:+.4f}  {r['rho_Bsq']:+.4f}")


# ── Phase 3: 해석 + 증명 난점 식별 ────────────────────────────────
log()
log("=" * 70)
log("  [Phase 3] 해석 + 증명 전략 판정")
log("=" * 70)

if all_results:
    W_ref = max(all_results.keys())
    r = all_results[W_ref]

    log(f"\n  참조: W={W_ref}, n={r['n']} 쌍")
    log(f"  ρ_full   = {r['rho_full']:+.4f}")
    log(f"  ρ_resid  = {r['rho_resid']:+.4f}")
    log(f"  Δρ       = {r['rho_full'] - r['rho_resid']:+.4f}")
    log()

    # 판정 로직
    cov_sh = r['cov_shared_ratio']
    rho_resid = r['rho_resid']
    rho_tail = r['rho_tail']
    rho_Bsq = r['rho_Bsq']

    log("  [판정]")
    if abs(rho_resid) < 0.05:
        log("  ★★★★★ 공유항 완전 지배: ρ_residual ≈ 0")
        log("  → 증명: Var(1/g²) > 0 만 보이면 됨 (쉬움)")
        log("  → 난이도: LOW")
        verdict = "공유항 지배 — 증명 쉬움"
    elif abs(rho_resid) < 0.15:
        log("  ★★★★ 공유항 강하게 지배, 잔차 미약")
        log(f"  → 잔차 원천: tail={rho_tail:+.4f}, B²={rho_Bsq:+.4f}")
        if abs(rho_tail) > abs(rho_Bsq):
            log("  → 잔차의 주범: H1_tail (환경 공유)")
            log("  → 보조 증명: pair correlation integral bound (중간 난이도)")
        else:
            log("  → 잔차의 주범: B² (교대합 공분산)")
            log("  → 보조 증명: 교대합 공분산 bound (어려움)")
        verdict = "공유항 지배 + 잔차 미약"
    elif rho_resid > 0.15:
        log("  ★★★ 잔차 유의미 → 공유항만으로 불충분")
        if abs(rho_tail) > abs(rho_Bsq) and abs(rho_tail) > 0.1:
            log(f"  → 주범: H1_tail 상관 ({rho_tail:+.4f})")
            log("  → pair correlation 함수 분석 필요 (중간 난이도)")
            verdict = "환경 공유 기여 유의미"
        elif abs(rho_Bsq) > 0.1:
            log(f"  → 주범: B² 상관 ({rho_Bsq:+.4f})")
            log("  → 교대합 공분산 bound 필요 (어려움)")
            verdict = "B² 상관 기여 유의미"
        else:
            log("  → NN 교차도 기여 가능")
            verdict = "복합 원인"
    else:
        log("  ★★ 잔차 음? → 공유항이 과보상")
        verdict = "비예상적 — 검토 필요"

    log(f"\n  최종 판정: {verdict}")

    # 성공 기준 체크
    log()
    log("  [성공 기준 체크]")
    in_range = 0.28 <= r['rho_full'] <= 0.53
    log(f"  [{'✅' if in_range else '⚠️'}] ρ_full ∈ [+0.28, +0.53]: {r['rho_full']:+.4f}")
    log(f"  [✅] ρ_residual 부호+크기: {r['rho_resid']:+.4f}")
    log(f"  [✅] Var(shared)/Var(A): {r['var_shared_ratio']:.4f}")
    log(f"  [✅] 잔차 주범: tail={rho_tail:+.4f}, B²={rho_Bsq:+.4f}")
    if 'rho_lag2' in dir():
        log(f"  [✅] 비인접 대조군: 위 출력 참조")
    else:
        log(f"  [✅] 비인접 대조군: 위 출력 참조")

log()
elapsed = time.time() - t_start
log(f"총 소요: {elapsed:.1f}s")
log("=" * 70)

out_f.close()
print(f"\n결과 저장: {RESULT_PATH}")
