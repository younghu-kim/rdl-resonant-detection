"""
=============================================================================
[Project RDL] 사이클 #256 — A(γ) vs 영점 간격 상관 분석
=============================================================================

목적:
  κ(σ,γ) = 1/(σ-1/2)² + A(γ) + O((σ-1/2)²) 에서
  A(γ) = Im(c₀)² + 2Re(c₁) 의 통계적 성질 탐사.

  구체적으로:
  1. ζ(s) 영점 200개에서 Laurent 계수 c₀, c₁ → A(γ)
  2. 각 영점의 gap_left, gap_right 계산
  3. 상관 분석:
     - Spearman ρ(A, gap_left/right/sum/min)
     - GUE 스케일링 (정규화 간격) 포함
  4. 인접 영점 쌍 A 상관: ρ(A_n, A_{n+1}) — "짝 패턴" 검증
  5. A의 t-의존성: A vs t 회귀
  6. B-42 연결: κ_mid vs (A₁+A₂)/2

방법:
  - mpmath.zetazero(n) 으로 영점 위치 + gap 계산
  - Cauchy 적분으로 c₀, c₁ 추출 (C-255 방법 재사용)
  - scipy.stats.spearmanr 로 상관 분석

결과 파일: results/A_gap_correlation_c256.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

from bundle_utils import xi_func

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

mpmath.mp.dps = 60   # 수학자 지시: C-255와 동일 dps

N_ZEROS = 200
CONTOUR_RADIUS = mpmath.mpf('0.01')
N_CONTOUR_PTS = 128
NU_DIFF_H = mpmath.mpf(1) / mpmath.mpf(10**18)

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/A_gap_correlation_c256.txt'
)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_A_coeffs(t_zero):
    """
    Cauchy 적분으로 c₀, c₁ 계산 (C-255 방법).

    ξ'/ξ(ρ + u) = 1/u + c₀ + c₁u + ...
    g(u) = ξ'/ξ(ρ+u) - 1/u (정칙 부분)
    cₙ = (1/N) Σ g(u_k) · u_k^{-n}

    반환: (A, c0, c1)
      A = Im(c₀)² + 2Re(c₁)  (Cor 4.2)
    """
    rho = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_zero)))
    r = CONTOUR_RADIUS
    N = N_CONTOUR_PTS
    h = NU_DIFF_H

    c0_sum = mpmath.mpc(0, 0)
    c1_sum = mpmath.mpc(0, 0)

    for k in range(N):
        theta = 2 * mpmath.pi * k / N
        u = r * mpmath.exp(mpmath.mpc(0, theta))
        s = rho + u

        # L = ξ'/ξ (중앙차분)
        xi_val = xi_func(s)
        xi_p = xi_func(s + h)
        xi_m = xi_func(s - h)
        L_val = (xi_p - xi_m) / (2 * h * xi_val)

        # g = L - 1/u
        g = L_val - 1 / u

        c0_sum += g            # g · u^0
        c1_sum += g / u        # g · u^{-1}

    c0 = c0_sum / N
    c1 = c1_sum / N

    c0_f = complex(c0)
    c1_f = complex(c1)

    A = c0_f.imag**2 + 2 * c1_f.real
    return A, c0_f, c1_f


def kappa_midpoint(t_mid):
    """
    κ = |ξ'/ξ(1/2 + it_mid)|²
    두 영점 사이 중간점에서의 곡률 (B-42 검증용).
    """
    s = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_mid)))
    h = mpmath.mpf(1) / mpmath.mpf(10**20)
    xi_val = xi_func(s)
    xi_p = xi_func(s + h)
    xi_m = xi_func(s - h)
    L = (xi_p - xi_m) / (2 * h * xi_val)
    return float(abs(L)**2)


def spearman_summary(name, x, y):
    """Spearman 상관 계산 + 결과 문자열 반환"""
    rho, pval = stats.spearmanr(x, y)
    sig = "✅ 유의" if pval < 0.01 else ("⚠️ p<0.05" if pval < 0.05 else "❌ 비유의")
    return rho, pval, f"  ρ({name}) = {rho:+.4f}  (p={pval:.3e})  {sig}"


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    print("=" * 80)
    print("[사이클 #256] A(γ) vs 영점 간격 상관 분석")
    print(f"  dps={mpmath.mp.dps}, N={N_ZEROS}, 반경={CONTOUR_RADIUS}, N_pts={N_CONTOUR_PTS}")
    print("=" * 80)

    # ──────────────────────────────────────────────────────────────────────
    # Step 0: 영점 수집
    # ──────────────────────────────────────────────────────────────────────
    print(f"\n[Step 0] ζ 영점 {N_ZEROS}개 수집 (mpmath.zetazero)...")
    zeros_t = []
    for n in range(1, N_ZEROS + 1):
        t = float(mpmath.zetazero(n).imag)
        zeros_t.append(t)
    zeros_t = np.array(zeros_t)
    print(f"  수집 완료: t ∈ [{zeros_t[0]:.2f}, {zeros_t[-1]:.2f}]")
    print(f"  최소 간격: {np.min(np.diff(zeros_t)):.4f}, 최대 간격: {np.max(np.diff(zeros_t)):.4f}")

    # ──────────────────────────────────────────────────────────────────────
    # Step 1: A(γ) 계산
    # ──────────────────────────────────────────────────────────────────────
    print(f"\n[Step 1] A(γ) 계산 (Cauchy 적분, 200 영점)...")
    A_list = []
    c0_list = []
    c1_list = []
    err_count = 0

    for j in range(N_ZEROS):
        t_n = zeros_t[j]
        if (j + 1) % 20 == 0 or j == 0:
            elapsed = time.time() - t_start
            print(f"  영점 #{j+1}/{N_ZEROS}  t={t_n:.2f}  경과={elapsed:.0f}s", flush=True)

        try:
            A, c0, c1 = compute_A_coeffs(t_n)
            if np.isnan(A) or np.isinf(A):
                print(f"  ⚠️ 영점 #{j+1} A=NaN/Inf — 건너뜀")
                A_list.append(np.nan)
                err_count += 1
            else:
                A_list.append(A)
            c0_list.append(c0)
            c1_list.append(c1)
        except Exception as e:
            print(f"  WARNING: 영점 #{j+1} t={t_n:.2f} 오류: {e}")
            A_list.append(np.nan)
            c0_list.append(complex(0))
            c1_list.append(complex(0))
            err_count += 1

    A_arr = np.array(A_list)
    valid_mask = ~np.isnan(A_arr)
    n_valid = np.sum(valid_mask)
    print(f"\n  계산 완료: {n_valid}/{N_ZEROS} 유효  ({err_count} 오류)")
    if err_count > N_ZEROS // 2:
        print("  ❌ 절반 이상 실패 — 중단")
        return

    if n_valid == 0:
        print("  ⚠️ 유효 데이터 없음 — 중단")
        return

    # Thm 5 검증
    re_c0_arr = np.array([c.real for c in c0_list])
    im_c0_arr = np.array([c.imag for c in c0_list])
    re_c1_arr = np.array([c.real for c in c1_list])
    im_c1_arr = np.array([c.imag for c in c1_list])

    im0_safe = np.where(im_c0_arr != 0, im_c0_arr, np.nan)
    re1_safe = np.where(re_c1_arr != 0, re_c1_arr, np.nan)
    thm5_re_c0 = float(np.nanmax(np.abs(re_c0_arr) / np.abs(im0_safe)))
    thm5_im_c1 = float(np.nanmax(np.abs(im_c1_arr) / np.abs(re1_safe)))
    print(f"  [Thm 5] max|Re(c₀)|/|Im(c₀)| = {thm5_re_c0:.2e}  (목표 < 1e-10)")
    print(f"  [Thm 5] max|Im(c₁)|/|Re(c₁)| = {thm5_im_c1:.2e}  (목표 < 1e-10)")

    # ──────────────────────────────────────────────────────────────────────
    # Step 2: 간격 계산
    # ──────────────────────────────────────────────────────────────────────
    print(f"\n[Step 2] gap_left, gap_right 계산...")

    gaps = np.diff(zeros_t)  # gaps[i] = zeros_t[i+1] - zeros_t[i], len=199

    # 내부 영점 (2..199번째, 인덱스 1..198): 양쪽 간격 모두 존재
    inner_idx = np.arange(1, N_ZEROS - 1)  # 인덱스 1..198 (내부 영점 198개)
    A_inner = A_arr[inner_idx]
    t_inner = zeros_t[inner_idx]

    gap_left  = gaps[inner_idx - 1]   # = t[n] - t[n-1]
    gap_right = gaps[inner_idx]        # = t[n+1] - t[n]
    gap_sum   = gap_left + gap_right
    gap_min   = np.minimum(gap_left, gap_right)

    # GUE 정규화: 평균 간격 ≈ 2π/log(t/2π)
    # normalized_gap = gap × log(t/(2π)) / (2π)
    log_t = np.log(t_inner / (2 * np.pi))
    norm_factor = log_t / (2 * np.pi)
    gap_left_n  = gap_left  * norm_factor
    gap_right_n = gap_right * norm_factor
    gap_sum_n   = gap_sum   * norm_factor
    gap_min_n   = gap_min   * norm_factor

    # 유효 내부 영점 필터
    valid_inner = ~np.isnan(A_inner)
    n_inner_valid = np.sum(valid_inner)
    print(f"  내부 영점: {n_inner_valid}/{len(inner_idx)} 유효")

    A_v = A_inner[valid_inner]
    gl_v = gap_left[valid_inner]
    gr_v = gap_right[valid_inner]
    gs_v = gap_sum[valid_inner]
    gm_v = gap_min[valid_inner]
    gln_v = gap_left_n[valid_inner]
    grn_v = gap_right_n[valid_inner]
    gsn_v = gap_sum_n[valid_inner]
    gmn_v = gap_min_n[valid_inner]

    # ──────────────────────────────────────────────────────────────────────
    # Step 3: Spearman 상관 분석
    # ──────────────────────────────────────────────────────────────────────
    print(f"\n[Step 3] Spearman 상관 분석 (n={n_inner_valid} 내부 영점)...")

    corr_results = []

    # (a) 원시 간격
    for name, arr in [
        ("A, gap_left",   gl_v),
        ("A, gap_right",  gr_v),
        ("A, gap_sum",    gs_v),
        ("A, gap_min",    gm_v),
    ]:
        rho, pval, line = spearman_summary(name, A_v, arr)
        corr_results.append((name, rho, pval))
        print(line)

    # (b) GUE 정규화 간격
    print()
    for name, arr in [
        ("A, gap_left_GUE",  gln_v),
        ("A, gap_right_GUE", grn_v),
        ("A, gap_sum_GUE",   gsn_v),
        ("A, gap_min_GUE",   gmn_v),
    ]:
        rho, pval, line = spearman_summary(name, A_v, arr)
        corr_results.append((name, rho, pval))
        print(line)

    # ──────────────────────────────────────────────────────────────────────
    # Step 4: 인접 영점 쌍 A 상관 — "짝 패턴" 검증
    # ──────────────────────────────────────────────────────────────────────
    print(f"\n[Step 4] 인접 영점 쌍 A 상관 (ρ(A_n, A_{{n+1}}))...")

    valid_full = ~np.isnan(A_arr)
    A_full = A_arr[valid_full]
    t_full = zeros_t[valid_full]
    # 연속 쌍 (A[n], A[n+1]) — valid indices만 사용
    valid_idx = np.where(valid_full)[0]
    # 연속 유효 쌍 찾기
    pair_mask = np.diff(valid_idx) == 1
    An  = A_arr[valid_idx[:-1][pair_mask]]
    An1 = A_arr[valid_idx[1:][pair_mask]]

    rho_pair, pval_pair = stats.spearmanr(An, An1)
    sig_pair = "✅ 유의" if pval_pair < 0.01 else ("⚠️ p<0.05" if pval_pair < 0.05 else "❌ 비유의")
    print(f"  ρ(A_n, A_{{n+1}}) = {rho_pair:+.4f}  (p={pval_pair:.3e})  n={len(An)}쌍  {sig_pair}")

    # 추가: ρ(A_n, A_{n+2}) — 두 칸 건너
    pair2_mask = (valid_idx[2:] - valid_idx[:-2]) == 2
    An_2  = A_arr[valid_idx[:-2][pair2_mask]]
    An2   = A_arr[valid_idx[2:][pair2_mask]]
    rho_pair2, pval_pair2 = stats.spearmanr(An_2, An2) if len(An_2) > 5 else (np.nan, np.nan)
    print(f"  ρ(A_n, A_{{n+2}}) = {rho_pair2:+.4f}  (p={pval_pair2:.3e})  n={len(An_2)}쌍")

    # ──────────────────────────────────────────────────────────────────────
    # Step 5: A vs t 회귀 (t-의존성 추정)
    # ──────────────────────────────────────────────────────────────────────
    print(f"\n[Step 5] A vs t 회귀...")

    t_v = t_full
    A_v_full = A_full

    # 선형: A ~ a + b·t
    slope_lin, intercept_lin, r_lin, p_lin, se_lin = stats.linregress(t_v, A_v_full)
    r2_lin = r_lin**2
    print(f"  선형 A = {intercept_lin:.4f} + {slope_lin:.6f}·t")
    print(f"    R² = {r2_lin:.4f}, p = {p_lin:.3e}")

    # 로그: A ~ a + b·log(t)
    log_t_v = np.log(t_v)
    slope_log, intercept_log, r_log, p_log, se_log = stats.linregress(log_t_v, A_v_full)
    r2_log = r_log**2
    print(f"  로그 A = {intercept_log:.4f} + {slope_log:.4f}·log(t)")
    print(f"    R² = {r2_log:.4f}, p = {p_log:.3e}")

    # 제곱근: A ~ a + b·sqrt(t)
    sqrt_t_v = np.sqrt(t_v)
    slope_sqrt, intercept_sqrt, r_sqrt, p_sqrt, se_sqrt = stats.linregress(sqrt_t_v, A_v_full)
    r2_sqrt = r_sqrt**2
    print(f"  제곱근 A = {intercept_sqrt:.4f} + {slope_sqrt:.4f}·√t")
    print(f"    R² = {r2_sqrt:.4f}, p = {p_sqrt:.3e}")

    # ──────────────────────────────────────────────────────────────────────
    # Step 6: B-42 연결 — κ_mid vs (A₁+A₂)/2
    # ──────────────────────────────────────────────────────────────────────
    print(f"\n[Step 6] B-42 연결: κ_mid vs (A_n + A_{{n+1}})/2...")
    print(f"  (중간점 곡률 계산 중, {N_ZEROS-1}쌍)")

    kappa_mid_list = []
    A_avg_list = []
    gap_list_b42 = []

    skip_b42 = 0
    for i in range(N_ZEROS - 1):
        if np.isnan(A_arr[i]) or np.isnan(A_arr[i + 1]):
            skip_b42 += 1
            continue
        t_mid = (zeros_t[i] + zeros_t[i + 1]) / 2
        try:
            km = kappa_midpoint(t_mid)
            if np.isnan(km) or np.isinf(km):
                skip_b42 += 1
                continue
            kappa_mid_list.append(km)
            A_avg_list.append((A_arr[i] + A_arr[i + 1]) / 2)
            gap_list_b42.append(zeros_t[i + 1] - zeros_t[i])
        except Exception as e:
            print(f"  WARNING: 쌍 ({i+1},{i+2}) κ_mid 오류: {e}")
            skip_b42 += 1

    kappa_mid_arr = np.array(kappa_mid_list)
    A_avg_arr = np.array(A_avg_list)
    gap_b42_arr = np.array(gap_list_b42)
    n_b42 = len(kappa_mid_arr)

    if n_b42 > 5:
        rho_b42_A, pval_b42_A = stats.spearmanr(kappa_mid_arr, A_avg_arr)
        rho_b42_g, pval_b42_g = stats.spearmanr(kappa_mid_arr, gap_b42_arr)
        sig_b42_A = "✅ 유의" if pval_b42_A < 0.01 else "❌ 비유의"
        sig_b42_g = "✅ 유의" if pval_b42_g < 0.01 else "❌ 비유의"
        print(f"  n={n_b42}쌍 (건너뜀: {skip_b42})")
        print(f"  ρ(κ_mid, (A₁+A₂)/2) = {rho_b42_A:+.4f}  p={pval_b42_A:.3e}  {sig_b42_A}")
        print(f"  ρ(κ_mid, gap)        = {rho_b42_g:+.4f}  p={pval_b42_g:.3e}  {sig_b42_g}")
        print(f"  κ_mid 범위: [{kappa_mid_arr.min():.3f}, {kappa_mid_arr.max():.3f}]")
        print(f"  (A₁+A₂)/2 범위: [{A_avg_arr.min():.3f}, {A_avg_arr.max():.3f}]")
    else:
        rho_b42_A = rho_b42_g = pval_b42_A = pval_b42_g = np.nan
        print(f"  ⚠️ B-42 데이터 부족 (n={n_b42})")

    # ──────────────────────────────────────────────────────────────────────
    # Step 7: A 기초 통계
    # ──────────────────────────────────────────────────────────────────────
    print(f"\n[Step 7] A(γ) 기초 통계...")

    A_valid = A_arr[valid_full]
    print(f"  n={len(A_valid)}, mean={np.mean(A_valid):.4f}, std={np.std(A_valid):.4f}")
    print(f"  min={np.min(A_valid):.4f}, max={np.max(A_valid):.4f}")
    print(f"  CV = {np.std(A_valid)/np.mean(A_valid)*100:.1f}%")
    print(f"  분위: Q1={np.percentile(A_valid,25):.4f}, Med={np.median(A_valid):.4f}, Q3={np.percentile(A_valid,75):.4f}")

    # ──────────────────────────────────────────────────────────────────────
    # Step 8: 종합 판정
    # ──────────────────────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("[종합 판정]")
    print("=" * 80)

    elapsed = time.time() - t_start

    # 최강 상관
    best_raw = max(corr_results[:4], key=lambda x: abs(x[1]))
    best_gue = max(corr_results[4:], key=lambda x: abs(x[1]))

    abs_best = max(abs(r[1]) for r in corr_results)
    verdict = "양성" if abs_best > 0.3 and any(r[2] < 0.01 for r in corr_results) else \
              ("음성(독립)" if abs_best < 0.1 else "중간")

    print(f"\n  최강 원시 상관: {best_raw[0]} — ρ={best_raw[1]:+.4f} (p={best_raw[2]:.3e})")
    print(f"  최강 GUE 상관:  {best_gue[0]} — ρ={best_gue[1]:+.4f} (p={best_gue[2]:.3e})")
    print(f"  인접 쌍 상관:   ρ(A_n,A_{{n+1}}) = {rho_pair:+.4f} (p={pval_pair:.3e})")
    print(f"  B-42 연결:      ρ(κ_mid,(A₁+A₂)/2) = {rho_b42_A:+.4f}")
    print(f"  A vs t 성장률:  R²_lin={r2_lin:.3f}, R²_log={r2_log:.3f}, b_log={slope_log:.4f}")
    print(f"\n  판정: {verdict}")
    if verdict == "양성":
        print(f"  → A(γ)가 gap 통계와 유의미한 상관 (B-42 메커니즘 후보)")
    elif verdict == "음성(독립)":
        print(f"  → A(γ)는 gap과 독립적인 새로운 양 (그 자체로 가치 있는 결과)")
    else:
        print(f"  → 보통 수준 상관. 추가 분석 필요.")

    print(f"\n  소요 시간: {elapsed:.1f}초")

    # ──────────────────────────────────────────────────────────────────────
    # 결과 파일 저장
    # ──────────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

    with open(RESULT_PATH, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("[사이클 #256] A(γ) vs 영점 간격 상관 분석\n")
        f.write(f"생성: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"소요: {elapsed:.1f}초  dps={mpmath.mp.dps}  N={N_ZEROS}\n")
        f.write("=" * 80 + "\n\n")

        # --- 설정 ---
        f.write("설정:\n")
        f.write(f"  N_ZEROS={N_ZEROS}, contour_radius={CONTOUR_RADIUS}, N_pts={N_CONTOUR_PTS}\n")
        f.write(f"  t 범위: [{zeros_t[0]:.2f}, {zeros_t[-1]:.2f}]\n")
        f.write(f"  유효: {n_valid}/{N_ZEROS}\n\n")

        # --- Laurent 계수 전수 ---
        f.write("=" * 80 + "\n")
        f.write("Laurent 계수 및 A(γ) 전수\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"  {'#':>4} {'t':>10} {'Re(c₀)':>12} {'Im(c₀)':>12} "
                f"{'Re(c₁)':>12} {'Im(c₁)':>12} {'A(γ)':>10} {'gap_right':>10}\n")
        f.write("  " + "-" * 88 + "\n")
        for j in range(N_ZEROS):
            c0 = c0_list[j]
            c1 = c1_list[j]
            A_j = A_arr[j]
            gr = gaps[j] if j < N_ZEROS - 1 else float('nan')
            f.write(f"  {j+1:4d} {zeros_t[j]:10.4f} {c0.real:12.6f} {c0.imag:12.6f} "
                    f"{c1.real:12.6f} {c1.imag:12.6f} {A_j:10.6f} {gr:10.4f}\n")

        f.write(f"\n  [Thm 5] max|Re(c₀)|/|Im(c₀)| = {thm5_re_c0:.2e}\n")
        f.write(f"  [Thm 5] max|Im(c₁)|/|Re(c₁)| = {thm5_im_c1:.2e}\n")

        # --- A 기초 통계 ---
        f.write("\n" + "=" * 80 + "\n")
        f.write("A(γ) 기초 통계\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"  n={len(A_valid)}, mean={np.mean(A_valid):.6f}, std={np.std(A_valid):.6f}\n")
        f.write(f"  min={np.min(A_valid):.6f}, max={np.max(A_valid):.6f}\n")
        f.write(f"  CV = {np.std(A_valid)/np.mean(A_valid)*100:.2f}%\n")
        f.write(f"  Q1={np.percentile(A_valid,25):.4f}, Med={np.median(A_valid):.4f}, "
                f"Q3={np.percentile(A_valid,75):.4f}\n")

        # --- 상관 분석 ---
        f.write("\n" + "=" * 80 + "\n")
        f.write("Spearman 상관 분석\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"  {'변수쌍':<30} {'ρ':>8} {'p값':>12} {'판정':>12}\n")
        f.write("  " + "-" * 70 + "\n")
        for name, rho, pval in corr_results:
            sig = "유의(p<0.01)" if pval < 0.01 else ("p<0.05" if pval < 0.05 else "비유의")
            f.write(f"  {name:<30} {rho:+8.4f} {pval:12.3e} {sig:>12}\n")

        f.write(f"\n  성공 기준: |ρ| > 0.3 (p<0.01) → 양성 / |ρ| < 0.1 → 음성(독립)\n")

        # --- 인접 쌍 상관 ---
        f.write("\n" + "=" * 80 + "\n")
        f.write("인접 영점 쌍 A 상관\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"  ρ(A_n, A_{{n+1}}) = {rho_pair:+.6f}  (p={pval_pair:.3e}, n={len(An)}쌍)\n")
        f.write(f"  ρ(A_n, A_{{n+2}}) = {rho_pair2:+.6f}  (p={pval_pair2:.3e}, n={len(An_2)}쌍)\n")
        if abs(rho_pair) > 0.3 and pval_pair < 0.01:
            f.write(f"  → '짝 패턴' 확인: 비국소 상관의 미시적 메커니즘\n")
        else:
            f.write(f"  → 인접 쌍 독립\n")

        # --- A vs t 회귀 ---
        f.write("\n" + "=" * 80 + "\n")
        f.write("A vs t 회귀\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"  선형:   A = {intercept_lin:.6f} + {slope_lin:.8f}·t,  R²={r2_lin:.4f}, p={p_lin:.3e}\n")
        f.write(f"  로그:   A = {intercept_log:.6f} + {slope_log:.6f}·log(t),  R²={r2_log:.4f}, p={p_log:.3e}\n")
        f.write(f"  제곱근: A = {intercept_sqrt:.6f} + {slope_sqrt:.6f}·√t,  R²={r2_sqrt:.4f}, p={p_sqrt:.3e}\n")
        best_model = max([("선형", r2_lin), ("로그", r2_log), ("제곱근", r2_sqrt)], key=lambda x: x[1])
        f.write(f"\n  최적 모델: {best_model[0]} (R²={best_model[1]:.4f})\n")

        # --- B-42 연결 ---
        f.write("\n" + "=" * 80 + "\n")
        f.write("B-42 연결: κ_mid vs (A₁+A₂)/2\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"  n={n_b42}쌍\n")
        f.write(f"  ρ(κ_mid, (A₁+A₂)/2) = {rho_b42_A:+.6f}  (p={pval_b42_A:.3e})\n")
        f.write(f"  ρ(κ_mid, gap)        = {rho_b42_g:+.6f}  (p={pval_b42_g:.3e})\n")
        if n_b42 > 5:
            f.write(f"  κ_mid: mean={np.mean(kappa_mid_arr):.4f}, std={np.std(kappa_mid_arr):.4f}\n")
            f.write(f"  (A₁+A₂)/2: mean={np.mean(A_avg_arr):.4f}, std={np.std(A_avg_arr):.4f}\n")
            if abs(rho_b42_A) > 0.3 and pval_b42_A < 0.01:
                f.write(f"  → B-42 메커니즘 후보: κ_mid ↔ A 평균 상관 발견\n")
            else:
                f.write(f"  → B-42: κ_mid ↔ A 상관 미발견\n")

        # --- 종합 판정 ---
        f.write("\n" + "=" * 80 + "\n")
        f.write("종합 판정\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"  최강 원시 상관: {best_raw[0]} — ρ={best_raw[1]:+.4f} (p={best_raw[2]:.3e})\n")
        f.write(f"  최강 GUE 상관:  {best_gue[0]} — ρ={best_gue[1]:+.4f} (p={best_gue[2]:.3e})\n")
        f.write(f"  인접 쌍:        ρ(A_n,A_{{n+1}}) = {rho_pair:+.4f} (p={pval_pair:.3e})\n")
        f.write(f"  B-42:           ρ(κ_mid,(A₁+A₂)/2) = {rho_b42_A:+.4f}\n")
        f.write(f"  A 성장률:       R²_log={r2_log:.3f}, b={slope_log:.4f} (A ~ log(t))\n")
        f.write(f"\n  판정: {verdict}\n")

    print(f"\n결과 저장: {RESULT_PATH}")
    print("완료.")


if __name__ == '__main__':
    main()
