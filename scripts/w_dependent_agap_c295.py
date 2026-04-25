#!/usr/bin/env python3
"""
[C-295] W-의존 A-gap 상관 프로파일: ζ(s) vs GUE

질문: δ_arith(W) = |ρ_GUE(W)| - |ρ_ζ(W)|가 W에 독립인가?

방법:
  - ζ(s): T=[10,2000], 이론적 정규화, trim 20%
  - GUE: N=2000, 30앙상블, 반원 unfolding, trim 20%
  - W (=N_MAX) ∈ {5, 10, 20, 40, 80, 120, 200, 300}
  - 각 W에서 ρ_S(A_bare, gap_min_theo) 측정

v2: numpy 벡터화로 O(n×W) 계산 최적화 (Python 루프 제거)

체크리스트:
  [x] dps=80 (t>100, ξ 언더플로 방지)
  [x] 이론적 d_bar = log(t/(2π))/(2π) for ζ(s)
  [x] GUE: 반원 CDF unfolding → 간격=1
  [x] trim 20%, edge skip
  [x] python -u
  [x] A_bare = S1^2 + 2*H1 (local ±W)
  [x] Spearman only (B-47)
"""

import sys, os, math, time
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 80

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
W_LIST = [5, 10, 20, 40, 80, 120, 200, 300]
GUE_N = 2000
GUE_ENS = 30
GUE_SEED = 42

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/w_dependent_agap_c295.txt'
)

out_f = open(RESULT_PATH, 'w')
def log(msg=''):
    print(msg, flush=True)
    out_f.write(msg + '\n')
    out_f.flush()


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


# ── GUE ──────────────────────────────────────────────────────────
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


# ── 벡터화된 A_bare 계산 ─────────────────────────────────────────
def compute_A_bare_vectorized(x, W):
    """x: sorted array of points. Returns A_bare for each interior point.

    벡터화: 오프셋 j=1..W에 대해 numpy 배열 연산으로 S1, H1 누적.
    범위: [W, n-W-1] (양쪽 W개 skip — 윈도우가 완전한 영점만 사용)
    """
    n = len(x)
    if n < 2 * W + 1:
        return np.array([]), np.arange(0)

    # 유효 인덱스: W <= i < n-W (양쪽 W개 여유 확보)
    valid_lo = W
    valid_hi = n - W
    n_valid = valid_hi - valid_lo
    if n_valid < 1:
        return np.array([]), np.arange(0)

    S1 = np.zeros(n_valid)
    H1 = np.zeros(n_valid)

    for j in range(1, W + 1):
        # 오른쪽 이웃: x[i+j] - x[i]
        diff_r = x[valid_lo + j : valid_hi + j] - x[valid_lo : valid_hi]
        # 왼쪽 이웃: x[i] - x[i-j]
        diff_l = x[valid_lo : valid_hi] - x[valid_lo - j : valid_hi - j]

        # 0 방지
        diff_r = np.where(np.abs(diff_r) < 1e-15, 1e-15, diff_r)
        diff_l = np.where(np.abs(diff_l) < 1e-15, 1e-15, diff_l)

        # S1 누적: 1/(x_i - x_{i+j}) + 1/(x_i - x_{i-j})
        S1 += -1.0 / diff_r + 1.0 / diff_l
        # H1 누적
        H1 += 1.0 / diff_r**2 + 1.0 / diff_l**2

    A = S1**2 + 2.0 * H1
    indices = np.arange(valid_lo, valid_hi)
    return A, indices


# ── ρ 측정 ───────────────────────────────────────────────────────
def measure_rho(x, W, density_func=None):
    """x: sorted array. density_func(x_i) → d_bar at x_i. None이면 d_bar=1."""
    n = len(x)
    A, indices = compute_A_bare_vectorized(x, W)
    if len(A) == 0:
        return np.nan, np.nan, 0

    # gap_min 계산
    gaps_l = x[indices] - x[indices - 1]
    gaps_r = x[indices + 1] - x[indices]
    gap_min = np.minimum(gaps_l, gaps_r)

    # 밀도 정규화
    if density_func is not None:
        d_bar = np.array([density_func(xi) for xi in x[indices]])
        gap_min = gap_min * d_bar

    # trim 20%
    n_pts = len(A)
    trim_lo = int(n_pts * TRIM_FRAC)
    trim_hi = n_pts - int(n_pts * TRIM_FRAC)
    if trim_hi - trim_lo < 10:
        return np.nan, np.nan, 0

    A_t = A[trim_lo:trim_hi]
    gm_t = gap_min[trim_lo:trim_hi]

    # NaN/inf 필터
    mask = np.isfinite(A_t) & np.isfinite(gm_t) & (A_t > 0) & (gm_t > 0)
    A_t = A_t[mask]
    gm_t = gm_t[mask]

    if len(A_t) < 10:
        return np.nan, np.nan, 0

    rho, p = stats.spearmanr(A_t, gm_t)
    return float(rho), float(p), len(A_t)


def zeta_density(t):
    """ζ(s) 이론적 영점 밀도"""
    if t <= 2 * math.pi:
        return 0.1
    return math.log(t / (2 * math.pi)) / (2 * math.pi)


# ── 메인 ─────────────────────────────────────────────────────────
def main():
    t_total = time.time()
    log("=" * 70)
    log("[C-295] W-의존 A-gap 상관 프로파일: ζ(s) vs GUE")
    log(f"  ζ(s): T=[10,{T_MAX}], 이론적 정규화, trim {TRIM_FRAC*100:.0f}%")
    log(f"  GUE: N={GUE_N}, {GUE_ENS}앙상블, 반원 unfolding, trim {TRIM_FRAC*100:.0f}%")
    log(f"  W = {W_LIST}")
    log(f"  v2: numpy 벡터화")
    log("=" * 70)

    # 1) ζ(s) 영점 수집
    zeros = get_zeta_zeros()

    # 2) GUE 앙상블 생성 + unfolding (한 번만)
    log(f"\n[GUE] {GUE_ENS}앙상블 생성 + unfolding ...")
    t0 = time.time()
    gue_unfolded_list = []
    for ens in range(GUE_ENS):
        seed = GUE_SEED + ens * 7919 + GUE_N
        lambdas = generate_gue_eigvals(GUE_N, seed=seed)
        xi = unfold(lambdas)
        gue_unfolded_list.append(xi)
    log(f"  완료: {GUE_ENS}개, {time.time()-t0:.1f}s")

    # 3) W별 ρ 측정
    log(f"\n{'='*70}")
    log("  W별 ρ_S 측정")
    log(f"{'='*70}")

    results = []
    for W in W_LIST:
        t0 = time.time()
        log(f"\n--- W = {W} ---")

        # ζ(s)
        rho_z, p_z, n_z = measure_rho(zeros, W, density_func=zeta_density)
        dt_z = time.time() - t0
        log(f"  ζ(s): ρ={rho_z:.6f}, p={p_z:.2e}, n={n_z}  ({dt_z:.1f}s)")

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
            # trim per ensemble
            n_pts = len(A_g)
            tl = int(n_pts * TRIM_FRAC)
            th = n_pts - int(n_pts * TRIM_FRAC)
            A_g_t = A_g[tl:th]
            gm_g_t = gm_g[tl:th]
            mask = np.isfinite(A_g_t) & np.isfinite(gm_g_t) & (A_g_t > 0) & (gm_g_t > 0)
            all_A_g.extend(A_g_t[mask].tolist())
            all_gm_g.extend(gm_g_t[mask].tolist())

        if len(all_A_g) >= 10:
            rho_g, p_g = stats.spearmanr(all_A_g, all_gm_g)
            n_g = len(all_A_g)
        else:
            rho_g, p_g, n_g = np.nan, np.nan, 0
        rho_g, p_g = float(rho_g), float(p_g)

        dt_g = time.time() - t1
        log(f"  GUE:  ρ={rho_g:.6f}, p={p_g:.2e}, n={n_g}  ({dt_g:.1f}s)")

        delta = abs(rho_g) - abs(rho_z) if (np.isfinite(rho_z) and np.isfinite(rho_g)) else np.nan
        log(f"  δ(W) = |ρ_GUE| - |ρ_ζ| = {delta:.6f}")

        results.append({
            'W': W,
            'rho_zeta': rho_z, 'p_zeta': p_z, 'n_zeta': n_z,
            'rho_gue': rho_g, 'p_gue': p_g, 'n_gue': n_g,
            'delta': delta,
        })

    # 4) 최종 비교표
    log(f"\n{'='*70}")
    log("  [최종 비교표] W-의존 A-gap 상관")
    log(f"{'='*70}\n")

    log(f"  {'W':>5}  {'ρ_ζ(s)':>10}  {'p_ζ':>10}  {'n_ζ':>5}  "
        f"{'ρ_GUE':>10}  {'p_GUE':>10}  {'n_GUE':>6}  {'δ(W)':>8}")
    log("  " + "-" * 78)

    for r in results:
        log(f"  {r['W']:>5}  {r['rho_zeta']:>10.6f}  {r['p_zeta']:>10.2e}  {r['n_zeta']:>5}  "
            f"{r['rho_gue']:>10.6f}  {r['p_gue']:>10.2e}  {r['n_gue']:>6}  {r['delta']:>8.4f}")

    # 5) δ 안정성 분석
    log(f"\n{'='*70}")
    log("  [δ 안정성 분석]")
    log(f"{'='*70}\n")

    deltas = [r['delta'] for r in results if np.isfinite(r['delta'])]
    if len(deltas) >= 3:
        d_mean = np.mean(deltas)
        d_std = np.std(deltas, ddof=1)
        d_cv = d_std / abs(d_mean) * 100 if abs(d_mean) > 1e-10 else float('inf')
        d_range = max(deltas) - min(deltas)

        log(f"  δ 평균: {d_mean:.4f}")
        log(f"  δ 표준편차: {d_std:.4f}")
        log(f"  δ CV: {d_cv:.1f}%")
        log(f"  δ 범위: {d_range:.4f}")
        log(f"  δ min/max: {min(deltas):.4f} / {max(deltas):.4f}")
        log("")

        # W≥40 서브셋
        deltas_large = [r['delta'] for r in results
                        if np.isfinite(r['delta']) and r['W'] >= 40]
        if len(deltas_large) >= 3:
            d_mean_L = np.mean(deltas_large)
            d_std_L = np.std(deltas_large, ddof=1)
            d_cv_L = d_std_L / abs(d_mean_L) * 100 if abs(d_mean_L) > 1e-10 else float('inf')
            log(f"  [W≥40 서브셋] δ 평균: {d_mean_L:.4f}, σ: {d_std_L:.4f}, CV: {d_cv_L:.1f}%")

        # 추세 검정
        W_arr = np.array([r['W'] for r in results if np.isfinite(r['delta'])])
        d_arr = np.array(deltas)
        rho_trend, p_trend = stats.spearmanr(W_arr, d_arr)
        log(f"\n  W-δ 추세: ρ_S={rho_trend:.4f}, p={p_trend:.3f}")

        if p_trend < 0.05:
            if rho_trend < 0:
                log("  → δ가 W 증가에 따라 유의미하게 감소 → 유한 윈도우 아티팩트 가능성")
            else:
                log("  → δ가 W 증가에 따라 유의미하게 증가 → 예상 외 패턴")
        else:
            log("  → δ와 W 사이 유의미한 추세 없음 → δ는 W-독립적")

    # 6) ρ 수렴 분석: ρ_ζ(W)와 ρ_GUE(W) 각각의 W-의존성
    log(f"\n{'='*70}")
    log("  [ρ 수렴 분석]")
    log(f"{'='*70}\n")

    rhos_z = [r['rho_zeta'] for r in results if np.isfinite(r['rho_zeta'])]
    rhos_g = [r['rho_gue'] for r in results if np.isfinite(r['rho_gue'])]
    ws = [r['W'] for r in results if np.isfinite(r['rho_zeta'])]

    if len(rhos_z) >= 3:
        rho_z_trend, p_z_trend = stats.spearmanr(ws, rhos_z)
        log(f"  ρ_ζ vs W: 추세 ρ_S={rho_z_trend:.4f}, p={p_z_trend:.3f}")
        log(f"    W=5: ρ_ζ={rhos_z[0]:.6f}, W={ws[-1]}: ρ_ζ={rhos_z[-1]:.6f}, Δ={rhos_z[-1]-rhos_z[0]:.6f}")

    if len(rhos_g) >= 3:
        rho_g_trend, p_g_trend = stats.spearmanr(ws, rhos_g)
        log(f"  ρ_GUE vs W: 추세 ρ_S={rho_g_trend:.4f}, p={p_g_trend:.3f}")
        log(f"    W=5: ρ_GUE={rhos_g[0]:.6f}, W={ws[-1]}: ρ_GUE={rhos_g[-1]:.6f}, Δ={rhos_g[-1]-rhos_g[0]:.6f}")

    # 7) 판정
    log(f"\n{'='*70}")
    log("  [판정]")
    log(f"{'='*70}\n")

    if len(deltas) >= 3:
        if d_cv < 25 and d_range < 0.03:
            log(f"  ★★★★★ δ(W) ≈ {d_mean:.4f} ± {d_std:.4f} (CV={d_cv:.1f}%, range={d_range:.4f})")
            log(f"  → W-독립적 잔여 감쇠. 진정한 산술 효과.")
            log(f"  → 해석: GUE 스펙트럼 강직성 ρ_GUE ≈ {np.mean(rhos_g):.3f}")
            log(f"           ζ(s) 산술 구조 ρ_ζ ≈ {np.mean(rhos_z):.3f}")
            log(f"           δ_arith ≈ {d_mean:.3f} (degree-independent)")
        elif d_cv < 25:
            log(f"  ★★★★ δ(W) ≈ {d_mean:.4f} (CV={d_cv:.1f}%)")
            log(f"  → 실질적으로 안정적. 범위 {d_range:.4f}는 허용 내.")
        elif p_trend < 0.05 and rho_trend < 0:
            log(f"  ★★★ δ(W)가 감소 추세 (ρ={rho_trend:.3f}, p={p_trend:.3f})")
            log(f"  → 부분적 윈도우 아티팩트. 외삽 필요.")
        elif p_trend < 0.05 and rho_trend > 0:
            log(f"  ★★★ δ(W)가 증가 추세 (ρ={rho_trend:.3f}, p={p_trend:.3f})")
            log(f"  → 대형 W에서 GUE-ζ 괴리 확대. 산술 효과 강화.")
        else:
            log(f"  ★★ δ(W) 변동 (CV={d_cv:.1f}%). 추가 분석 필요.")
    else:
        log("  데이터 부족 — 판정 불가")

    log(f"\n  [참조]")
    log(f"  C-282b ζ(s) ρ=-0.929 (N_MAX=200, 이론적)")
    log(f"  C-283 GUE_local ρ=-0.967 (N_MAX=300)")
    log(f"  C-293 GUE_smooth ρ=-0.973 (N_MAX=300, smooth)")
    log(f"  C-291 GL(2) ρ=-0.934 (이론적)")
    log(f"  v1 W=5: δ=0.039, W=10: δ=0.040 (Python loop, 재현 확인용)")

    log(f"\n  총 소요: {time.time()-t_total:.1f}s")
    log("=" * 70)
    out_f.close()
    print(f"\n✅ 결과 저장: {RESULT_PATH}")


if __name__ == "__main__":
    main()
