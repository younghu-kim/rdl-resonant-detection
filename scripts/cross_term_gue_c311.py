#!/usr/bin/env python3
"""
=============================================================================
[C-311] 교차항-GUE pair correlation 정량 분석 (B-61 탐사)
=============================================================================
목표:
  1. E_diag 해석적 검증: E_diag = πN/δ vs 수치 비교
  2. 교차항 쌍별 분해: E_cross = Σ_{i<j} 8πδ/(Δ²+4δ²) (해석적 정확)
  3. GUE pair correlation: R₂_empirical vs R₂_GUE = 1-(sinπx/πx)²
  4. 교차항 스케일링: E_cross ~ B/δ^β → β 측정

수학:
  Hadamard 표현에서 f_n(t) = 1/(δ+i(t-γ_n)), E = ∫|Σ f_n|² dt

  E_diag = Σ_n ∫|f_n|² dt = Σ_n π/δ · [arctan 보정]
  E_cross = Σ_{i≠j} ∫ f_i·f̄_j dt

  핵심 공식 (contour integral):
    ∫_{-∞}^{∞} f_i·f̄_j dt = 2πi/(Δ_{ij}+2iδ)
    where Δ_{ij} = γ_i - γ_j

  따라서 쌍 {i,j} (i<j)의 기여:
    2Re[2πi/(Δ+2iδ) + 2πi/(-Δ+2iδ)] = 8πδ/(Δ²+4δ²)

  E_cross = Σ_{i<j} 8πδ/(Δ_{ij}² + 4δ²)

결과: results/cross_term_gue_c311.txt
=============================================================================
"""

import sys, os, time
import numpy as np
from scipy import stats

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

OUTFILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/cross_term_gue_c311.txt')
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))

def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파라미터
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

N_USE = 400          # 사용할 영점 수 (처음 400개: 가장 밀도 높은 데이터)
N_T_POINTS = 2000    # 수치 적분 격자 (검증용)

# δ (= σ - 1/2) 격자
DELTAS = np.array([0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10,
                   0.15, 0.20, 0.30, 0.50, 1.00])

# GUE pair correlation 비닝
N_BINS_R2 = 80       # R₂ histogram 빈 수
R2_MAX_X = 4.0       # 정규화 간격 최대값

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. 영점 로드
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def load_zeros():
    """캐시에서 ζ 영점 로드 (N4000 우선)"""
    cache_4k = os.path.expanduser(
        '~/Desktop/gdl_unified/outputs/cache/zeta_zeros_t_N4000_c296.npy')
    cache_5k = os.path.expanduser(
        '~/Desktop/gdl_unified/outputs/cache/zeta_zeros_t_N5000.npy')

    if os.path.exists(cache_4k):
        all_zeros = np.load(cache_4k)
        log(f"  캐시 로드: {cache_4k} ({len(all_zeros)}개)")
    elif os.path.exists(cache_5k):
        all_zeros = np.load(cache_5k)
        log(f"  캐시 로드: {cache_5k} ({len(all_zeros)}개)")
    else:
        raise FileNotFoundError("영점 캐시 없음")

    zeros = all_zeros[:N_USE]
    log(f"  사용: {len(zeros)}개, t ∈ [{zeros[0]:.4f}, {zeros[-1]:.4f}]")

    gaps = np.diff(zeros)
    mean_gap = gaps.mean()
    log(f"  평균 간격: {mean_gap:.6f}")

    return zeros, all_zeros


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. E_diag: 해석적 + 이론값
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def E_diag_analytic(zeros, delta, t_min, t_max):
    """
    E_diag = Σ_n [arctan((t_max-γ_n)/δ) - arctan((t_min-γ_n)/δ)] / δ

    각 Lorentzian의 정확한 적분.
    """
    if abs(delta) < 1e-15:
        delta = 1e-15
    total = 0.0
    for g in zeros:
        total += (np.arctan((t_max - g) / delta) - np.arctan((t_min - g) / delta)) / delta
    return total


def E_diag_theory(N, delta):
    """E_diag ≈ πN/δ (무한 적분 구간 근사)"""
    return np.pi * N / delta


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. E_cross: 해석적 (쌍별 정확 공식)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def E_cross_analytic(zeros, delta):
    """
    E_cross = Σ_{i<j} 8πδ / (Δ_{ij}² + 4δ²)

    Δ_{ij} = |γ_i - γ_j|

    해석적 정확: contour 적분에서 도출.
    경계 효과 없음 (무한 적분 구간 사용 시 정확).
    """
    N = len(zeros)
    total = 0.0
    four_d2 = 4.0 * delta * delta
    coeff = 8.0 * np.pi * delta

    for i in range(N):
        for j in range(i + 1, N):
            D2 = (zeros[i] - zeros[j]) ** 2
            total += coeff / (D2 + four_d2)

    return total


def E_cross_analytic_vectorized(zeros, delta):
    """벡터화 버전 (N이 크면 O(N²) 메모리 → N=400이면 160K 쌍, OK)"""
    N = len(zeros)
    coeff = 8.0 * np.pi * delta
    four_d2 = 4.0 * delta**2

    # 상삼각 쌍만
    D = zeros[:, None] - zeros[None, :]  # (N, N)
    D2 = D ** 2
    # 상삼각 마스크
    mask = np.triu(np.ones((N, N), dtype=bool), k=1)
    terms = coeff / (D2[mask] + four_d2)
    return terms.sum(), terms, D[mask]


def E_cross_by_spacing_bin(zeros, delta, n_bins=50, max_spacing=None):
    """
    교차항을 쌍 간격 Δ별로 그룹화.
    C(Δ_bin, δ) = Σ_{pairs in bin} 8πδ/(Δ²+4δ²)
    """
    N = len(zeros)
    gaps = np.diff(zeros)
    mean_gap = gaps.mean()
    if max_spacing is None:
        max_spacing = 10.0 * mean_gap

    coeff = 8.0 * np.pi * delta
    four_d2 = 4.0 * delta**2

    D = np.abs(zeros[:, None] - zeros[None, :])
    mask_upper = np.triu(np.ones((N, N), dtype=bool), k=1)
    spacings = D[mask_upper]
    terms = coeff / (spacings**2 + four_d2)

    bin_edges = np.linspace(0, max_spacing, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_sums = np.zeros(n_bins)
    bin_counts = np.zeros(n_bins, dtype=int)

    indices = np.digitize(spacings, bin_edges) - 1
    for k in range(n_bins):
        in_bin = (indices == k)
        bin_sums[k] = terms[in_bin].sum()
        bin_counts[k] = in_bin.sum()

    return bin_centers, bin_sums, bin_counts


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. E_Had 수치 검증
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def E_hadamard_numeric(ts, zeros, delta):
    """
    수치 적분: E = ∫|Σ_n 1/(δ+i(t-γ_n))|² dt

    반환: E_total (스칼라)
    """
    tau = ts[:, None] - zeros[None, :]  # (N_t, N_z)
    denom = delta**2 + tau**2

    re_sum = np.sum(delta / denom, axis=1)
    im_sum = np.sum(-tau / denom, axis=1)

    integrand = re_sum**2 + im_sum**2
    return np.trapezoid(integrand, ts)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. GUE pair correlation
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def R2_GUE(x):
    """GUE pair correlation: R₂(x) = 1 - (sin πx / πx)²"""
    x = np.asarray(x, dtype=float)
    result = np.ones_like(x)
    nonzero = np.abs(x) > 1e-12
    sinc_val = np.sin(np.pi * x[nonzero]) / (np.pi * x[nonzero])
    result[nonzero] = 1.0 - sinc_val**2
    result[~nonzero] = 0.0  # R₂(0) = 0 (level repulsion)
    return result


def compute_pair_correlation(zeros, n_bins=N_BINS_R2, max_x=R2_MAX_X):
    """
    경험적 pair correlation 함수 계산.

    절차:
    1. 각 영점 γ_n에서 local mean spacing d̄(γ_n) = 2π/log(γ_n/(2π))
    2. 정규화 간격 x_{nm} = |γ_n - γ_m| / d̄(γ_n)
    3. histogram → R₂(x) = (쌍 수) / (기대 쌍 수)
    """
    N = len(zeros)

    # local mean spacing (Riemann-von Mangoldt)
    d_bar = 2.0 * np.pi / np.log(zeros / (2.0 * np.pi))

    # 정규화된 모든 쌍 간격 (중복 방지: i < j, 양쪽 d_bar 평균 사용)
    all_x = []
    for i in range(N):
        for j in range(i + 1, N):
            d_mean = 0.5 * (d_bar[i] + d_bar[j])
            x = abs(zeros[i] - zeros[j]) / d_mean
            if x < max_x:
                all_x.append(x)

    all_x = np.array(all_x)
    log(f"  pair correlation: {len(all_x)} 쌍 (x < {max_x})")

    # histogram
    bin_edges = np.linspace(0, max_x, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    dx = bin_edges[1] - bin_edges[0]

    counts, _ = np.histogram(all_x, bins=bin_edges)

    # 정규화: 기대 쌍 밀도 = N(N-1)/2 · dx / max_x (균일 분포 가정)
    # 실제로는 R₂=1일 때의 기대 수: N · (T/d̄_avg) · dx per unit x
    # 간단히: counts / (N * dx * density)
    # density of pairs at distance x: ~ N/d̄_avg
    d_bar_avg = d_bar.mean()
    expected_per_bin = N * dx / d_bar_avg  # 한 영점당, 한 bin당 기대 쌍 수
    R2_empirical = counts / (expected_per_bin + 1e-30)

    return bin_centers, R2_empirical, counts


def compute_pair_correlation_fast(zeros, n_bins=N_BINS_R2, max_x=R2_MAX_X):
    """
    빠른 버전: 인접 쌍만 사용 (간격 < max_x * max_d_bar)
    먼 쌍은 R₂≈1 기여하므로 가까운 쌍이 핵심.
    """
    N = len(zeros)
    d_bar = 2.0 * np.pi / np.log(zeros / (2.0 * np.pi))
    d_bar_max = d_bar.max()

    # 최대 물리 간격
    max_phys = max_x * d_bar_max

    all_x = []
    for i in range(N):
        for j in range(i + 1, N):
            gap = zeros[j] - zeros[i]
            if gap > max_phys:
                break  # zeros는 정렬되어 있으므로
            d_mean = 0.5 * (d_bar[i] + d_bar[j])
            x = gap / d_mean
            if x < max_x:
                all_x.append(x)

    all_x = np.array(all_x)
    log(f"  pair correlation (fast): {len(all_x)} 쌍 (x < {max_x})")

    bin_edges = np.linspace(0, max_x, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    dx = bin_edges[1] - bin_edges[0]

    counts, _ = np.histogram(all_x, bins=bin_edges)

    # 정규화: 각 bin에서 R₂=1일 때의 기대 수
    # 총 "유효" 영점 쌍: Σ_i (이웃 수) → 근사: N
    # R₂=1일 때: 각 영점은 양쪽으로 x=max_x까지의 이웃 가짐
    # 기대 쌍 수 (R₂=1): N · dx (정규화 간격 단위)
    expected_per_bin = N * dx
    R2_empirical = counts / (expected_per_bin + 1e-30)

    return bin_centers, R2_empirical, counts


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. 멱법칙 피팅
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def fit_power_law(ds, E, label=""):
    """E ~ A / δ^α → log(E) = log(A) - α·log(δ)"""
    mask = (ds > 0) & (E > 0) & np.isfinite(E)
    if mask.sum() < 3:
        return 0, 0, 0, 0
    x = np.log(ds[mask])
    y = np.log(E[mask])
    A_mat = np.column_stack([np.ones(mask.sum()), x])
    coeffs = np.linalg.lstsq(A_mat, y, rcond=None)[0]
    alpha = -coeffs[1]
    A = np.exp(coeffs[0])
    y_pred = A_mat @ coeffs
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0
    return alpha, A, r2, int(mask.sum())


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7. 교차항-R₂ 편차 상관 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def cross_r2_correlation(zeros, delta, r2_centers, r2_empirical, r2_gue):
    """
    교차항 기여도 C(Δ)와 R₂ 편차(R₂_emp - R₂_GUE)의 상관.

    절차:
    1. 각 쌍 (i,j)의 교차항 기여 w_{ij} = 8πδ/(Δ²+4δ²)
    2. 정규화 간격 x_{ij} = Δ_{ij}/d̄ 계산
    3. x_{ij} 빈별로 w 합산 → W(x_bin)
    4. corr(W(x_bin), R₂_emp - R₂_GUE)
    """
    N = len(zeros)
    d_bar = 2.0 * np.pi / np.log(zeros / (2.0 * np.pi))
    coeff = 8.0 * np.pi * delta
    four_d2 = 4.0 * delta**2

    n_bins = len(r2_centers)
    dx = r2_centers[1] - r2_centers[0]
    max_x = r2_centers[-1] + dx / 2

    W_bins = np.zeros(n_bins)
    d_bar_max = d_bar.max()
    max_phys = max_x * d_bar_max

    for i in range(N):
        for j in range(i + 1, N):
            gap = zeros[j] - zeros[i]
            if gap > max_phys:
                break
            d_mean = 0.5 * (d_bar[i] + d_bar[j])
            x = gap / d_mean
            if x >= max_x:
                continue
            w = coeff / (gap**2 + four_d2)
            bk = int(x / dx)
            if 0 <= bk < n_bins:
                W_bins[bk] += w

    # 상관
    deviation = r2_empirical - r2_gue
    # 유효 빈만 (counts > 0)
    valid = (r2_empirical > 0) & np.isfinite(W_bins)
    if valid.sum() < 5:
        return 0.0, 1.0, W_bins

    corr_val, p_val = stats.pearsonr(W_bins[valid], deviation[valid])
    return corr_val, p_val, W_bins


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Main
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t0_global = time.time()

    log("=" * 72)
    log("[C-311] 교차항-GUE pair correlation 정량 분석 (B-61)")
    log("=" * 72)
    log(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"N_USE={N_USE}, N_T_POINTS={N_T_POINTS}, N_DELTAS={len(DELTAS)}")
    log()

    # ─── 1. 영점 로드 ───
    log("[1] 영점 로드")
    zeros, all_zeros = load_zeros()
    N = len(zeros)
    piN = np.pi * N
    gaps = np.diff(zeros)
    mean_gap = gaps.mean()
    t_min, t_max = zeros[0] - 5 * mean_gap, zeros[-1] + 5 * mean_gap  # 경계 여유
    ts = np.linspace(t_min, t_max, N_T_POINTS)

    log(f"  N={N}, πN={piN:.4f}, mean_gap={mean_gap:.6f}")
    log(f"  적분 구간: [{t_min:.2f}, {t_max:.2f}]")
    log()

    # ─── 2. E_diag 해석적 검증 ───
    log("[2] E_diag 해석적 검증")
    log(f"  {'δ':>8s} {'E_diag_exact':>14s} {'πN/δ':>14s} {'비율':>8s} {'E_diag_num':>14s} {'exact/num':>10s}")
    log("  " + "-" * 74)

    E_diag_arr = np.zeros(len(DELTAS))
    E_diag_theory_arr = np.zeros(len(DELTAS))
    E_diag_num_arr = np.zeros(len(DELTAS))

    for k, delta in enumerate(DELTAS):
        ed_exact = E_diag_analytic(zeros, delta, t_min, t_max)
        ed_theory = E_diag_theory(N, delta)
        # 수치 검증 (대각항만)
        tau = ts[:, None] - zeros[None, :]
        denom = delta**2 + tau**2
        diag_integrand = np.sum(delta**2 / denom**2 + tau**2 / denom**2, axis=1)
        # 실제로 |f_n|² = 1/(δ²+τ²), 합은 Σ 1/(δ²+τ²)
        diag_integrand2 = np.sum(1.0 / denom, axis=1)
        ed_num = np.trapezoid(diag_integrand2, ts)

        E_diag_arr[k] = ed_exact
        E_diag_theory_arr[k] = ed_theory
        E_diag_num_arr[k] = ed_num

        ratio1 = ed_exact / ed_theory if ed_theory > 0 else 0
        ratio2 = ed_exact / ed_num if ed_num > 0 else 0
        log(f"  {delta:8.4f} {ed_exact:14.2f} {ed_theory:14.2f} {ratio1:8.5f} {ed_num:14.2f} {ratio2:10.6f}")

    # A_diag / πN
    log()
    log("  [A_diag/πN 분석]")
    for k, delta in enumerate(DELTAS):
        A_diag = E_diag_arr[k] * delta  # E_diag ≈ A/δ → A = E·δ
        ratio = A_diag / piN
        log(f"    δ={delta:.4f}: A_diag/πN = {ratio:.6f}")
    log()

    # ─── 3. 교차항 정확 계산 ───
    log("[3] 교차항 해석적 계산")
    log(f"  {'δ':>8s} {'E_cross':>14s} {'E_diag':>14s} {'E_Had=d+c':>14s} {'E_Had_num':>14s} {'num/(d+c)':>10s} {'cross/diag%':>12s}")
    log("  " + "-" * 88)

    E_cross_arr = np.zeros(len(DELTAS))
    E_had_num_arr = np.zeros(len(DELTAS))

    for k, delta in enumerate(DELTAS):
        t1 = time.time()

        # 해석적 교차항
        ec = E_cross_analytic(zeros, delta)
        E_cross_arr[k] = ec

        # Hadamard 수치 (검증)
        eh_num = E_hadamard_numeric(ts, zeros, delta)
        E_had_num_arr[k] = eh_num

        # 해석적 합 (무한 구간) vs 수치 (유한 구간)
        e_sum = E_diag_theory_arr[k] + ec  # πN/δ + E_cross
        e_sum_exact = E_diag_arr[k] + ec   # arctan 보정 + E_cross

        ratio = eh_num / e_sum_exact if e_sum_exact > 0 else 0
        cross_pct = ec / E_diag_arr[k] * 100 if E_diag_arr[k] > 0 else 0
        dt = time.time() - t1

        log(f"  {delta:8.4f} {ec:14.2f} {E_diag_arr[k]:14.2f} {e_sum_exact:14.2f} {eh_num:14.2f} {ratio:10.6f} {cross_pct:12.2f}% [{dt:.1f}s]")

    log()

    # 교차항 경계 보정 주의
    log("  [참고] 교차항 해석적 공식은 무한 적분 구간 가정.")
    log("  경계 보정: E_cross_finite ≈ E_cross_inf × (유효 비율)")
    log("  num/(d+c) 비율이 ~1이면 보정 불필요 확인.")
    log()

    # ─── 4. 교차항 스케일링 (β 측정) ───
    log("[4] 교차항 스케일링: E_cross ~ B/δ^β")

    alpha_cross, A_cross, r2_cross, n_cross = fit_power_law(DELTAS, E_cross_arr)
    log(f"  전체: β = {alpha_cross:.4f}, R² = {r2_cross:.6f} (n={n_cross})")

    # 소δ (≤0.10)
    mask_sm = DELTAS <= 0.10
    if mask_sm.sum() >= 3:
        a_sm, _, r2_sm, n_sm = fit_power_law(DELTAS[mask_sm], E_cross_arr[mask_sm])
        log(f"  δ≤0.10: β = {a_sm:.4f}, R² = {r2_sm:.6f} (n={n_sm})")

    # 초소δ (≤0.05)
    mask_xs = DELTAS <= 0.05
    if mask_xs.sum() >= 3:
        a_xs, _, r2_xs, n_xs = fit_power_law(DELTAS[mask_xs], E_cross_arr[mask_xs])
        log(f"  δ≤0.05: β = {a_xs:.4f}, R² = {r2_xs:.6f} (n={n_xs})")

    # 대δ (≥0.10)
    mask_lg = DELTAS >= 0.10
    if mask_lg.sum() >= 3:
        a_lg, _, r2_lg, n_lg = fit_power_law(DELTAS[mask_lg], E_cross_arr[mask_lg])
        log(f"  δ≥0.10: β = {a_lg:.4f}, R² = {r2_lg:.6f} (n={n_lg})")

    log()

    # E_Had = E_diag + E_cross 의 α 측정
    E_had_analytic_arr = E_diag_arr + E_cross_arr
    alpha_had, _, r2_had, _ = fit_power_law(DELTAS, E_had_analytic_arr)
    alpha_diag, _, r2_diag, _ = fit_power_law(DELTAS, E_diag_arr)
    log(f"  비교: E_diag α = {alpha_diag:.4f} (R²={r2_diag:.6f})")
    log(f"  비교: E_Had  α = {alpha_had:.4f} (R²={r2_had:.6f})")
    log(f"  판정: β_cross={alpha_cross:.4f} {'< 1 → α=1 보존 ✓' if alpha_cross < 1 else '≥ 1 → 주의!'}")
    log()

    # ─── 5. 교차항 쌍별 분해 (Δ-프로파일) ───
    log("[5] 교차항 Δ-프로파일")
    log("  교차항을 쌍 간격 Δ별로 그룹화 (대표 δ 3개)")

    for delta_repr in [0.01, 0.05, 0.20]:
        centers, sums, counts = E_cross_by_spacing_bin(zeros, delta_repr, n_bins=30, max_spacing=8*mean_gap)
        log(f"\n  δ = {delta_repr}")
        log(f"  {'Δ/d̄':>8s} {'C(Δ)':>12s} {'쌍수':>8s} {'누적%':>8s}")
        total_cross = sums.sum()
        cumul = 0
        for b in range(len(centers)):
            if counts[b] > 0:
                cumul += sums[b]
                pct = cumul / total_cross * 100 if total_cross > 0 else 0
                log(f"  {centers[b]/mean_gap:8.3f} {sums[b]:12.4f} {counts[b]:8d} {pct:8.1f}%")

    log()

    # ─── 6. GUE pair correlation ───
    log("[6] GUE pair correlation 분석")

    # 4000개 영점 사용 (통계 강화)
    zeros_gue = all_zeros[:min(2000, len(all_zeros))]  # 2000개면 충분
    log(f"  R₂ 계산: {len(zeros_gue)}개 영점 사용")

    r2_centers, r2_emp, r2_counts = compute_pair_correlation_fast(zeros_gue)
    r2_gue = R2_GUE(r2_centers)

    log(f"\n  {'x':>6s} {'R₂_emp':>10s} {'R₂_GUE':>10s} {'Δ=emp-GUE':>12s} {'쌍수':>8s}")
    log("  " + "-" * 52)
    for b in range(len(r2_centers)):
        if r2_counts[b] > 0:
            dev = r2_emp[b] - r2_gue[b]
            log(f"  {r2_centers[b]:6.3f} {r2_emp[b]:10.4f} {r2_gue[b]:10.4f} {dev:+12.4f} {r2_counts[b]:8d}")

    # KS test: empirical CDF vs GUE CDF
    log()
    log("  [통계 검정]")

    # 개별 정규화 간격 수집 (인접 쌍만)
    d_bar_gue = 2.0 * np.pi / np.log(zeros_gue / (2.0 * np.pi))
    nn_spacings = np.diff(zeros_gue) / (0.5 * (d_bar_gue[:-1] + d_bar_gue[1:]))

    # GUE 최인접 간격 분포 (Wigner surmise, β=2):
    # p(s) = (32/π²)s²exp(-4s²/π)
    # CDF는 수치 적분으로 계산
    from scipy import integrate as _sci_int

    def gue_pdf(s):
        return (32.0 / np.pi**2) * s**2 * np.exp(-4.0 * s**2 / np.pi)

    def gue_cdf(s):
        if np.isscalar(s):
            val, _ = _sci_int.quad(gue_pdf, 0, s)
            return val
        return np.array([_sci_int.quad(gue_pdf, 0, si)[0] for si in s])

    ks_stat, ks_p = stats.kstest(nn_spacings, gue_cdf)
    log(f"  최인접 간격 KS test (GUE Wigner surmise): D = {ks_stat:.4f}, p = {ks_p:.4f}")
    log(f"  해석: p > 0.05 → GUE와 일치 {'✓' if ks_p > 0.05 else '✗'}")

    # χ² test for R₂ histogram
    valid_bins = r2_counts > 5  # 기대값 > 5인 빈만
    if valid_bins.sum() >= 5:
        # 기대 빈도 = N_pairs · dx · R₂_GUE(x)
        N_gue = len(zeros_gue)
        dx_r2 = r2_centers[1] - r2_centers[0]
        expected_counts = N_gue * dx_r2 * r2_gue
        valid_bins2 = valid_bins & (expected_counts > 5)
        if valid_bins2.sum() >= 5:
            # 기대값 합을 관측값 합에 맞춤 (χ² 요구사항)
            obs = r2_counts[valid_bins2].astype(float)
            exp = expected_counts[valid_bins2]
            exp = exp * (obs.sum() / exp.sum())  # 합 정규화
            chi2, chi2_p = stats.chisquare(obs, exp)
            log(f"  R₂ histogram χ² test: χ² = {chi2:.2f}, p = {chi2_p:.4f} (n_bins={valid_bins2.sum()})")
        else:
            log(f"  R₂ χ² test: 유효 빈 부족 (valid={valid_bins2.sum()})")
    log()

    # ─── 7. 교차항-R₂ 편차 상관 ───
    log("[7] 교차항 기여 vs R₂ 편차 상관")

    # 400개 영점 기준 (교차항 계산용)
    r2_c400, r2_e400, r2_cnt400 = compute_pair_correlation_fast(zeros)
    r2_g400 = R2_GUE(r2_c400)

    for delta_corr in [0.01, 0.05, 0.10]:
        corr, p, W = cross_r2_correlation(zeros, delta_corr, r2_c400, r2_e400, r2_g400)
        log(f"  δ={delta_corr:.2f}: Pearson r = {corr:+.4f}, p = {p:.4f} {'(유의)' if p < 0.05 else '(무의미)'}")
    log()

    # ─── 8. 교차항 부호 전환점 ───
    log("[8] 교차항 부호 전환 분석")
    log("  E_cross는 항상 양수 (모든 항 8πδ/(Δ²+4δ²) > 0)")
    log("  그러나 경계 보정 후 E_Had_num - E_diag_exact 의 부호는 변할 수 있음")
    log()
    log(f"  {'δ':>8s} {'E_Had_num-E_diag':>16s} {'E_cross_inf':>14s} {'잔차':>12s}")
    for k, delta in enumerate(DELTAS):
        diff = E_had_num_arr[k] - E_diag_arr[k]
        resid = diff - E_cross_arr[k]  # 경계 효과
        log(f"  {delta:8.4f} {diff:16.2f} {E_cross_arr[k]:14.2f} {resid:12.2f}")
    log()

    # GL(d) 교차항 부호 전환과의 비교 (수학자 지시)
    log("  [GL(d) 교차항 부호 전환 비교]")
    log("  GL(1) ζ(s): E_cross > 0 항상 (Hadamard 이론)")
    log("  GL(2): Δσ_cross ≈ 0.10에서 cross% 부호 전환 (C-303)")
    log("  GL(4): Δσ_cross ≈ 0.05에서 cross% 부호 전환 (C-308c)")
    log("  → GL(1)은 Gamma 배경 없음. GL(d≥2)의 부호 전환은 Gamma-Had 간섭에서 기인.")
    log()

    # ─── 9. 종합 ───
    log("=" * 72)
    log("[종합] C-311 결과")
    log("=" * 72)

    log()
    log("1. E_diag 검증:")
    # 경계 보정 후 비율
    boundary_ratios = E_diag_arr * DELTAS / piN
    log(f"   A_diag/πN: min={boundary_ratios.min():.5f}, max={boundary_ratios.max():.5f}")
    log(f"   δ≤0.05 평균: {boundary_ratios[DELTAS<=0.05].mean():.6f}")
    log(f"   πN/δ 이론값과 차이: < 2% {'✓' if abs(1 - boundary_ratios[DELTAS<=0.05].mean()) < 0.02 else '✗'}")

    log()
    log("2. 교차항 스케일링:")
    log(f"   β = {alpha_cross:.4f} (전체), R² = {r2_cross:.6f}")
    if mask_sm.sum() >= 3:
        log(f"   β = {a_sm:.4f} (δ≤0.10)")
    if mask_xs.sum() >= 3:
        log(f"   β = {a_xs:.4f} (δ≤0.05)")
    log(f"   판정: β {'< 1 → E_cross는 E_diag보다 약해서 α=1 보존' if alpha_cross < 1 else '≥ 1 → 문제'}")

    log()
    log("3. GUE pair correlation:")
    log(f"   최인접 간격 KS test: D={ks_stat:.4f}, p={ks_p:.4f}")
    log(f"   GUE 일치: {'✓' if ks_p > 0.05 else '✗'}")

    log()
    log("4. 교차항-R₂ 상관:")
    log(f"   (δ=0.05): Pearson r, p 위 참조")

    log()
    elapsed = time.time() - t0_global
    log(f"총 실행 시간: {elapsed:.1f}s ({elapsed/60:.1f}분)")
    log()

    flush_file()
    log(f"결과 저장: {OUTFILE}")


if __name__ == "__main__":
    main()
