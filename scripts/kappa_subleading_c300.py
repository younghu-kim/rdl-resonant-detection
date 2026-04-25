"""
=============================================================================
[C-300] κ 부선도(Sub-Leading) 구조: Hadamard 전개에서의 해석적 예측 vs 수치
=============================================================================
목적:
    κδ² = 1 + c₁δ + c₂δ² + O(δ³) 의 Taylor 계수가 Hadamard 급수의
    S₁ = Σ 1/(γₙ-γₖ), H₁ = Σ 1/(γₙ-γₖ)² 로 예측 가능한지 검증.

    핵심 예측: c₂ ≈ const - 2H₁  →  ρ(c₂, H₁) < -0.9

    이것이 확인되면 Prop 5 (κδ²→1) 와 Prop 12 (A≥2/gap_min²) 를
    단일 전개로 통합하는 정리가 된다.

방법:
    1. ζ(s) 영점 γ₁,...,γ_N (t∈[14,50]) 확보
    2. 각 영점에서 δ ∈ [0.001, 0.15] 다점 κ(δ) 측정
    3. κδ² = 1 + c₁δ + c₂δ² 다항식 적합
    4. Hadamard 급수로 S₁, H₁ 직접 계산
    5. c₂ vs -2H₁ 상관 분석

판정 기준:
    양성: ρ(c₂, H₁) < -0.9 AND c₂ ≈ a - 2H₁ 적합 R²>0.95
    중립: ρ < -0.7 OR R²>0.8
    음성: 상관 미약
=============================================================================
"""

import sys
import time
import numpy as np
import mpmath

# bundle_utils 사용
sys.path.insert(0, '/home/k0who029/Desktop/gdl_unified/scripts')
from bundle_utils import xi_func, connection_zeta, find_zeros_zeta

mpmath.mp.dps = 80  # t<50 이므로 80 충분

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

T_MIN, T_MAX = 14.0, 50.0
DELTA_VALUES = np.array([0.001, 0.002, 0.003, 0.005, 0.007,
                          0.01, 0.015, 0.02, 0.03, 0.05,
                          0.07, 0.1, 0.12, 0.15])
N_NEIGHBORS = 200  # Hadamard 급수 절단 (S₁, H₁ 계산용)
POLY_DEGREE = 2    # κδ² = 1 + c₁δ + c₂δ²

OUT_FILE = '/home/k0who029/Desktop/gdl_unified/results/kappa_subleading_c300.txt'


def compute_kappa(gamma_n, delta):
    """영점 γₙ에서 오프셋 δ만큼 이동한 점의 곡률 κ"""
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(gamma_n + delta))
    L = connection_zeta(s)
    return float(abs(L)**2)


def compute_hadamard_sums(gamma_n, all_zeros, n_idx):
    """Hadamard 급수 S₁, H₁, S₃ 계산 (최근접 N_NEIGHBORS개)"""
    S1 = 0.0  # Σ 1/(γₙ-γₖ)
    H1 = 0.0  # Σ 1/(γₙ-γₖ)²
    S3 = 0.0  # Σ 1/(γₙ-γₖ)³

    for k_idx, gk in enumerate(all_zeros):
        if k_idx == n_idx:
            continue
        diff = gamma_n - gk
        if abs(diff) < 1e-10:
            continue
        S1 += 1.0 / diff
        H1 += 1.0 / diff**2
        S3 += 1.0 / diff**3

    return S1, H1, S3


def main():
    t_start = time.time()

    # 1. 영점 탐색
    print("=" * 70)
    print(f"[C-300] κ 부선도 구조: 해석적 예측 vs 수치 검증")
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"t 구간: [{T_MIN}, {T_MAX}]")
    print("=" * 70)

    print("\n[1] 영점 탐색...")
    zeros = find_zeros_zeta(T_MIN, T_MAX)
    N = len(zeros)
    print(f"  발견 영점: {N}개")
    for i, g in enumerate(zeros):
        print(f"  γ_{i+1} = {g:.6f}")

    # 확장 영점 (S₁, H₁ 정밀 계산용) — 양쪽 확장
    print(f"\n[2] 확장 영점 탐색 (S₁/H₁ 정밀도용)...")
    ext_zeros = find_zeros_zeta(1.0, max(T_MAX + 100, 200))
    print(f"  확장 영점: {len(ext_zeros)}개 (t∈[1, {max(T_MAX+100, 200)}])")

    # 2. 각 영점에서 κδ² 측정
    print(f"\n[3] κδ² 측정 (영점 {N}개 × δ {len(DELTA_VALUES)}개)...")

    kappa_d2_data = np.zeros((N, len(DELTA_VALUES)))

    for i, gn in enumerate(zeros):
        for j, delta in enumerate(DELTA_VALUES):
            kappa = compute_kappa(gn, delta)
            kappa_d2_data[i, j] = kappa * delta**2
        print(f"  γ_{i+1} = {gn:.4f}: κδ²(0.001)={kappa_d2_data[i,0]:.6f}, "
              f"κδ²(0.1)={kappa_d2_data[i,9]:.6f}")

    # 3. 다항식 적합: κδ² = 1 + c₁δ + c₂δ²
    print(f"\n[4] 다항식 적합 (κδ² - 1 = c₁δ + c₂δ²)...")

    c1_fit = np.zeros(N)
    c2_fit = np.zeros(N)
    r2_fit = np.zeros(N)

    for i in range(N):
        y = kappa_d2_data[i, :] - 1.0  # κδ² - 1
        x = DELTA_VALUES

        # 원점을 지나는 2차: y = c₁x + c₂x²
        # 행렬: [x, x²] @ [c₁, c₂]ᵀ = y
        A = np.column_stack([x, x**2])
        coeffs, residuals, _, _ = np.linalg.lstsq(A, y, rcond=None)
        c1_fit[i] = coeffs[0]
        c2_fit[i] = coeffs[1]

        y_pred = A @ coeffs
        ss_res = np.sum((y - y_pred)**2)
        ss_tot = np.sum((y - y.mean())**2)
        r2_fit[i] = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

        print(f"  γ_{i+1}: c₁={c1_fit[i]:+.4f}, c₂={c2_fit[i]:+.4f}, R²={r2_fit[i]:.6f}")

    # 4. Hadamard 급수 S₁, H₁ 직접 계산
    print(f"\n[5] Hadamard 급수 (S₁, H₁) 계산...")

    S1_arr = np.zeros(N)
    H1_arr = np.zeros(N)
    S3_arr = np.zeros(N)
    gap_min_arr = np.zeros(N)

    for i, gn in enumerate(zeros):
        # ext_zeros에서 n_idx 찾기
        n_idx = np.argmin(np.abs(np.array(ext_zeros) - gn))
        S1, H1, S3 = compute_hadamard_sums(gn, ext_zeros, n_idx)
        S1_arr[i] = S1
        H1_arr[i] = H1
        S3_arr[i] = S3

        # gap_min
        sorted_diffs = sorted([abs(gn - g) for g in ext_zeros if abs(gn - g) > 1e-10])
        gap_min_arr[i] = sorted_diffs[0] if sorted_diffs else 1.0

        print(f"  γ_{i+1}: S₁={S1:+.4f}, H₁={H1:.4f}, gap_min={gap_min_arr[i]:.4f}")

    # 5. 핵심 상관 분석
    print(f"\n[6] 핵심 상관 분석")

    # c₂ vs H₁
    from scipy import stats

    r_c2_H1, p_c2_H1 = stats.spearmanr(c2_fit, H1_arr)
    r_c2_H1_p, p_c2_H1_p = stats.pearsonr(c2_fit, H1_arr)

    print(f"  ρ_Spearman(c₂, H₁) = {r_c2_H1:.4f}, p = {p_c2_H1:.2e}")
    print(f"  ρ_Pearson(c₂, H₁)  = {r_c2_H1_p:.4f}, p = {p_c2_H1_p:.2e}")

    # c₂ vs -2H₁ 적합: c₂ = a - 2H₁
    slope_pred = -2.0
    A_lin = np.column_stack([np.ones(N), H1_arr])
    coeffs_lin = np.linalg.lstsq(A_lin, c2_fit, rcond=None)[0]
    intercept_lin, slope_lin = coeffs_lin[0], coeffs_lin[1]
    y_pred_lin = A_lin @ coeffs_lin
    ss_res_lin = np.sum((c2_fit - y_pred_lin)**2)
    ss_tot_lin = np.sum((c2_fit - c2_fit.mean())**2)
    r2_lin = 1.0 - ss_res_lin / ss_tot_lin if ss_tot_lin > 0 else 0.0

    print(f"\n  c₂ = a + b·H₁ 적합:")
    print(f"    a (intercept) = {intercept_lin:.4f}")
    print(f"    b (slope)     = {slope_lin:.4f}  (이론 예측: -2.0)")
    print(f"    R²            = {r2_lin:.6f}")
    print(f"    |b - (-2)| / 2 = {abs(slope_lin - (-2.0)) / 2.0:.4f} (상대 편차)")

    # c₁ vs S₁ (보조 검증)
    r_c1_S1, p_c1_S1 = stats.spearmanr(c1_fit, S1_arr)
    print(f"\n  ρ_Spearman(c₁, S₁) = {r_c1_S1:.4f}, p = {p_c1_S1:.2e}")

    # c₁ vs 2S₁ 적합: c₁ = a' + 2S₁ (이론)
    coeffs_c1 = np.linalg.lstsq(np.column_stack([np.ones(N), S1_arr]),
                                 c1_fit, rcond=None)[0]
    print(f"  c₁ = a' + b'·S₁: b'={coeffs_c1[1]:.4f} (이론 예측: 2.0)")

    # c₂ vs gap_min (간접)
    r_c2_gap, p_c2_gap = stats.spearmanr(c2_fit, gap_min_arr)
    print(f"  ρ_Spearman(c₂, gap_min) = {r_c2_gap:.4f}, p = {p_c2_gap:.2e}")

    # H₁ vs 1/gap_min² (Prop 12 확인)
    r_H1_gap, p_H1_gap = stats.spearmanr(H1_arr, 1.0/gap_min_arr**2)
    print(f"  ρ_Spearman(H₁, 1/gap_min²) = {r_H1_gap:.4f}, p = {p_H1_gap:.2e}")

    # 6. 잔차 분석: c₂ + 2H₁ = |R₀|²-2S₁Im(C₀)+S₁² ≈ const?
    residual = c2_fit + 2.0 * H1_arr  # 이론적으로 이것이 상수
    print(f"\n[7] 잔차 분석: c₂ + 2H₁")
    print(f"  평균: {residual.mean():.4f}")
    print(f"  표준편차: {residual.std():.4f}")
    print(f"  CV: {abs(residual.std()/residual.mean())*100:.1f}%")
    print(f"  범위: [{residual.min():.4f}, {residual.max():.4f}]")

    # c₂ + 2H₁ vs S₁² (이론적으로 상관 있어야)
    r_res_S1sq, p_res_S1sq = stats.spearmanr(residual, S1_arr**2)
    print(f"  ρ(c₂+2H₁, S₁²) = {r_res_S1sq:.4f}, p = {p_res_S1sq:.2e}")

    elapsed = time.time() - t_start

    # 7. 종합 판정
    print(f"\n{'='*70}")
    print("[8] 종합 판정")

    # 기준 1: c₂-H₁ 음상관
    test1 = r_c2_H1 < -0.9
    test1_weak = r_c2_H1 < -0.7
    # 기준 2: 선형 적합 R²
    test2 = r2_lin > 0.95
    test2_weak = r2_lin > 0.80
    # 기준 3: 기울기 ≈ -2
    test3 = abs(slope_lin - (-2.0)) < 0.5
    # 기준 4: 잔차 일관성 (CV < 30%)
    test4 = abs(residual.std() / residual.mean()) < 0.30 if abs(residual.mean()) > 0.01 else True

    if test1 and test2 and test3:
        verdict = "★★★★★ 강양성"
        verdict_detail = "c₂=-2H₁+const 확정. Prop 5+12 통합 정리."
    elif (test1 or test1_weak) and (test2 or test2_weak):
        verdict = "★★★★ 양성"
        verdict_detail = "c₂∝H₁ 확인. 기울기 또는 잔차에 보정 필요."
    elif test1_weak or test2_weak:
        verdict = "★★★ 중립"
        verdict_detail = "상관 존재하나 정량적 일치 부족."
    else:
        verdict = "★★ 음성"
        verdict_detail = "c₂-H₁ 관계 미확인."

    print(f"  판정: {verdict}")
    print(f"  근거: {verdict_detail}")
    print(f"  경과: {elapsed:.1f}s")

    # 결과 파일 저장
    with open(OUT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"[C-300] κ 부선도 구조: Hadamard 전개 검증\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"경과: {elapsed:.1f}s\n")
        f.write(f"t 구간: [{T_MIN}, {T_MAX}], N={N} 영점\n")
        f.write(f"δ 범위: [{DELTA_VALUES[0]}, {DELTA_VALUES[-1]}], {len(DELTA_VALUES)}점\n")
        f.write(f"확장 영점: {len(ext_zeros)}개 (S₁/H₁용)\n")
        f.write("=" * 70 + "\n")

        f.write(f"\n[1] 영점 목록\n")
        for i, g in enumerate(zeros):
            f.write(f"  γ_{i+1} = {g:.6f}\n")

        f.write(f"\n[2] κδ² 측정 (영점×δ)\n")
        f.write(f"  {'δ':>8s}")
        for i in range(N):
            f.write(f"  {'γ_'+str(i+1):>8s}")
        f.write("\n")
        for j, delta in enumerate(DELTA_VALUES):
            f.write(f"  {delta:8.4f}")
            for i in range(N):
                f.write(f"  {kappa_d2_data[i,j]:8.5f}")
            f.write("\n")

        f.write(f"\n[3] 다항식 적합: κδ² - 1 = c₁δ + c₂δ²\n")
        f.write(f"  {'γ':>8s}  {'c₁':>10s}  {'c₂':>10s}  {'R²':>8s}\n")
        for i in range(N):
            f.write(f"  {zeros[i]:8.4f}  {c1_fit[i]:+10.4f}  {c2_fit[i]:+10.4f}  {r2_fit[i]:8.6f}\n")

        f.write(f"\n[4] Hadamard 급수\n")
        f.write(f"  {'γ':>8s}  {'S₁':>10s}  {'H₁':>10s}  {'S₃':>10s}  {'gap_min':>8s}\n")
        for i in range(N):
            f.write(f"  {zeros[i]:8.4f}  {S1_arr[i]:+10.4f}  {H1_arr[i]:10.4f}  {S3_arr[i]:+10.4f}  {gap_min_arr[i]:8.4f}\n")

        f.write(f"\n[5] 핵심 상관 분석\n")
        f.write(f"  ρ_Spearman(c₂, H₁) = {r_c2_H1:.4f}, p = {p_c2_H1:.2e}\n")
        f.write(f"  ρ_Pearson(c₂, H₁)  = {r_c2_H1_p:.4f}, p = {p_c2_H1_p:.2e}\n")
        f.write(f"\n  c₂ = a + b·H₁:\n")
        f.write(f"    a = {intercept_lin:.4f}, b = {slope_lin:.4f} (이론: -2.0), R² = {r2_lin:.6f}\n")
        f.write(f"    |b-(-2)|/2 = {abs(slope_lin+2)/2:.4f}\n")
        f.write(f"\n  ρ_Spearman(c₁, S₁) = {r_c1_S1:.4f}, p = {p_c1_S1:.2e}\n")
        f.write(f"  c₁ = a'+b'·S₁: b' = {coeffs_c1[1]:.4f} (이론: 2.0)\n")
        f.write(f"\n  ρ_Spearman(c₂, gap_min) = {r_c2_gap:.4f}, p = {p_c2_gap:.2e}\n")
        f.write(f"  ρ_Spearman(H₁, 1/gap_min²) = {r_H1_gap:.4f}, p = {p_H1_gap:.2e}\n")

        f.write(f"\n[6] 잔차: c₂ + 2H₁ (이론적으로 ≈ S₁² + |C₀|²)\n")
        f.write(f"  평균 = {residual.mean():.4f}, σ = {residual.std():.4f}, "
                f"CV = {abs(residual.std()/residual.mean())*100 if abs(residual.mean()) > 0.01 else 0:.1f}%\n")
        f.write(f"  범위: [{residual.min():.4f}, {residual.max():.4f}]\n")
        f.write(f"  ρ(c₂+2H₁, S₁²) = {r_res_S1sq:.4f}, p = {p_res_S1sq:.2e}\n")
        for i in range(N):
            f.write(f"  γ_{i+1}: c₂+2H₁ = {residual[i]:.4f}, S₁² = {S1_arr[i]**2:.4f}\n")

        f.write(f"\n[7] 종합 판정: {verdict}\n")
        f.write(f"  근거: {verdict_detail}\n")
        f.write(f"  기준: ρ(c₂,H₁)<-0.9={test1}, R²>0.95={test2}, "
                f"|b+2|<0.5={test3}, 잔차CV<30%={test4}\n")

        f.write(f"\n{'='*70}\n")
        f.write(f"[정리 후보] κ Sub-Leading Expansion Theorem\n\n")
        f.write(f"  Theorem (RH 조건부):\n")
        f.write(f"    ξ(s)의 단순 영점 ρₙ = 1/2+iγₙ에서, δ>0 소:\n")
        f.write(f"    κ(γₙ,δ)·δ² = 1 + c₁δ + c₂δ² + O(δ³)\n")
        f.write(f"    c₁ = 2S₁ + const   (S₁ = Σ_{{k≠n}} 1/(γₙ-γₖ))\n")
        f.write(f"    c₂ = S₁² + |C₀|² - 2H₁   (H₁ = Σ_{{k≠n}} 1/(γₙ-γₖ)²)\n")
        f.write(f"\n  수치 검증 (N={N}, t∈[{T_MIN},{T_MAX}]):\n")
        f.write(f"    ρ(c₂, H₁) = {r_c2_H1:.4f}\n")
        f.write(f"    c₂ = {intercept_lin:.3f} + ({slope_lin:.3f})·H₁, R²={r2_lin:.4f}\n")
        f.write(f"    c₁ = const + ({coeffs_c1[1]:.3f})·S₁\n")

    print(f"\n  결과: {OUT_FILE}")
    print(f"  총 경과: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
