"""
=============================================================================
[Project RDL] 사이클 #254 — σ-국소화 해석적 구조 탐사 (Paper 4 씨앗)
=============================================================================

목적:
  σ=1/2에서 κ(곡률)가 발산하는 메커니즘의 유형을 판별.
  - algebraic singularity: log κ ~ α·log|σ - 1/2| + β
  - exponential/essential: log κ ~ γ / |σ - 1/2|
  - logarithmic: log κ ~ δ·log(log(1/|σ - 1/2|))

방법:
  1. ζ(s) 처음 20개 영점의 t_n에 대해 19개 σ값에서 κ(σ + i·t_n) 측정
  2. log κ vs log|σ - 0.5| 피팅 → α 추정 (algebraic order)
  3. log κ vs 1/|σ - 0.5| 피팅 → γ 추정 (exponential)
  4. FE 대칭 검증: κ(σ, t) vs κ(1-σ, t)
  5. 20개 영점에 걸친 α, γ 분포 → CV 계산

결과 파일: results/sigma_localization_profile_c254.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

from bundle_utils import xi_func, find_zeros_zeta

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

mpmath.mp.dps = 60  # σ≈0.5 근방 소거 오차 방지 (수학자: ≥50)

# 19개 σ값 (수학자 지정)
SIGMAS = [
    0.3, 0.35, 0.40, 0.45, 0.48, 0.49, 0.495, 0.499, 0.4999,
    0.5,
    0.5001, 0.501, 0.505, 0.51, 0.52, 0.55, 0.60, 0.65, 0.70
]

N_ZEROS = 20  # 처음 20개 영점

RESULT_PATH = os.path.expanduser('~/Desktop/gdl_unified/results/sigma_localization_profile_c254.txt')


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수: 임의 s에서 κ 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def curvature_at_s(sigma, t):
    """
    κ(σ + it) = |ξ''/ξ - (ξ'/ξ)²|  (= |L'| where L = ξ'/ξ)

    수학자 주: bundle_utils의 curvature_zeta는 |ξ'/ξ|²를 반환.
    off-critical에서도 동일하게 |L|² = |ξ'/ξ|²를 사용.
    """
    s = mpmath.mpc(mpmath.mpf(str(sigma)), mpmath.mpf(str(t)))
    h = mpmath.mpf(1) / mpmath.mpf(10**20)

    xi_val = xi_func(s)

    # ξ가 영점에 매우 가까우면 κ는 극대값
    if abs(xi_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 15):
        return float('inf')

    # ξ'/ξ (접속)
    xi_plus = xi_func(s + h)
    xi_minus = xi_func(s - h)
    L = (xi_plus - xi_minus) / (2 * h * xi_val)

    return float(abs(L)**2)


def curvature_connection_at_s(sigma, t):
    """
    |L|² 와 함께 L 자체도 반환 (디버깅/분석용)
    """
    s = mpmath.mpc(mpmath.mpf(str(sigma)), mpmath.mpf(str(t)))
    h = mpmath.mpf(1) / mpmath.mpf(10**20)

    xi_val = xi_func(s)

    if abs(xi_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 15):
        return float('inf'), complex(1e10, 0)

    xi_plus = xi_func(s + h)
    xi_minus = xi_func(s - h)
    L = (xi_plus - xi_minus) / (2 * h * xi_val)

    kappa = float(abs(L)**2)
    L_complex = complex(L)
    return kappa, L_complex


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 피팅 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def fit_algebraic(delta_sigma, log_kappa):
    """
    log κ ~ α · log|σ - 1/2| + β 피팅
    σ=0.5 제외, |σ-0.5| > 0인 점만 사용
    반환: α, β, R²
    """
    mask = np.array(delta_sigma) > 0
    x = np.log(np.array(delta_sigma)[mask])
    y = np.array(log_kappa)[mask]

    if len(x) < 3:
        return np.nan, np.nan, np.nan

    # 선형 피팅
    A = np.vstack([x, np.ones_like(x)]).T
    result = np.linalg.lstsq(A, y, rcond=None)
    coeffs = result[0]
    alpha, beta = coeffs[0], coeffs[1]

    y_pred = alpha * x + beta
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return alpha, beta, r_sq


def fit_exponential(delta_sigma, log_kappa):
    """
    log κ ~ γ / |σ - 1/2| + δ 피팅
    반환: γ, δ, R²
    """
    mask = np.array(delta_sigma) > 0
    x = 1.0 / np.array(delta_sigma)[mask]  # 1/|σ - 0.5|
    y = np.array(log_kappa)[mask]

    if len(x) < 3:
        return np.nan, np.nan, np.nan

    A = np.vstack([x, np.ones_like(x)]).T
    result = np.linalg.lstsq(A, y, rcond=None)
    coeffs = result[0]
    gamma, delta = coeffs[0], coeffs[1]

    y_pred = gamma * x + delta
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return gamma, delta, r_sq


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    start_time = time.time()

    print("=" * 80)
    print("[사이클 #254] σ-국소화 해석적 구조 탐사")
    print(f"  dps = {mpmath.mp.dps}")
    print(f"  σ 값: {len(SIGMAS)}개")
    print(f"  영점: 처음 {N_ZEROS}개")
    print(f"  총 측정: {len(SIGMAS)} × {N_ZEROS} = {len(SIGMAS) * N_ZEROS}개")
    print("=" * 80)

    # Step 0: 영점 수집
    print("\n[Step 0] ζ 영점 수집...")
    zeros_t = []
    for n in range(1, N_ZEROS + 1):
        t = float(mpmath.zetazero(n).imag)
        zeros_t.append(t)
        print(f"  영점 #{n}: t = {t:.6f}")

    # Step 1: κ 프로파일 측정
    print(f"\n[Step 1] κ(σ, t_n) 측정 ({len(SIGMAS)} × {N_ZEROS} = {len(SIGMAS)*N_ZEROS}개)...")

    # kappa_matrix[i][j] = κ(σ_i, t_j)
    kappa_matrix = np.zeros((len(SIGMAS), N_ZEROS))

    for j, t_n in enumerate(zeros_t):
        print(f"\n  영점 #{j+1}/{N_ZEROS}: t = {t_n:.4f}")
        for i, sigma in enumerate(SIGMAS):
            try:
                kappa = curvature_at_s(sigma, t_n)
                kappa_matrix[i][j] = kappa

                if kappa == float('inf'):
                    print(f"    σ={sigma:.4f}: κ = ∞ (영점)")
                elif kappa > 1e10:
                    print(f"    σ={sigma:.4f}: κ = {kappa:.4e} (극대)")
                else:
                    print(f"    σ={sigma:.4f}: κ = {kappa:.4e}")
            except Exception as e:
                print(f"    σ={sigma:.4f}: ERROR — {e}")
                kappa_matrix[i][j] = np.nan

    # Step 2: log κ vs σ 분석
    print("\n" + "=" * 80)
    print("[Step 2] 발산 유형 분석: algebraic vs exponential")
    print("=" * 80)

    delta_sigmas = [abs(s - 0.5) for s in SIGMAS]

    alphas = []  # algebraic order for each zero
    gammas = []  # exponential coefficient for each zero
    r2_alg = []
    r2_exp = []

    for j, t_n in enumerate(zeros_t):
        kappa_col = kappa_matrix[:, j]

        # σ=0.5에서의 κ (기준값)
        idx_half = SIGMAS.index(0.5)
        kappa_half = kappa_col[idx_half]

        # σ ≠ 0.5인 점들만 사용
        ds = []
        lk = []
        for i, sigma in enumerate(SIGMAS):
            if abs(sigma - 0.5) > 1e-10 and kappa_col[i] > 0 and not np.isnan(kappa_col[i]) and kappa_col[i] != float('inf'):
                ds.append(abs(sigma - 0.5))
                lk.append(np.log(kappa_col[i]))

        alpha, beta, r2a = fit_algebraic(ds, lk)
        gamma, delta_coeff, r2e = fit_exponential(ds, lk)

        alphas.append(alpha)
        gammas.append(gamma)
        r2_alg.append(r2a)
        r2_exp.append(r2e)

        print(f"\n  영점 #{j+1} (t={t_n:.4f}), κ(1/2) = {kappa_half:.4e}:")
        print(f"    Algebraic: α = {alpha:.4f}, R² = {r2a:.4f}")
        print(f"    Exponential: γ = {gamma:.6f}, R² = {r2e:.4f}")
        winner = "algebraic" if r2a > r2e else "exponential"
        print(f"    → 우세: {winner} (R² 차이: {abs(r2a - r2e):.4f})")

    # Step 2b: α, γ 통계
    alphas_arr = np.array([a for a in alphas if not np.isnan(a)])
    gammas_arr = np.array([g for g in gammas if not np.isnan(g)])

    print(f"\n  α (algebraic order) 통계:")
    print(f"    mean = {np.mean(alphas_arr):.4f}")
    print(f"    std  = {np.std(alphas_arr):.4f}")
    print(f"    CV   = {np.std(alphas_arr)/abs(np.mean(alphas_arr))*100:.1f}%")
    print(f"    min  = {np.min(alphas_arr):.4f}, max = {np.max(alphas_arr):.4f}")

    print(f"\n  γ (exponential coeff) 통계:")
    print(f"    mean = {np.mean(gammas_arr):.6f}")
    print(f"    std  = {np.std(gammas_arr):.6f}")

    print(f"\n  R² 비교 (알제브라 vs 지수):")
    print(f"    mean R²_alg = {np.mean(r2_alg):.4f}")
    print(f"    mean R²_exp = {np.mean(r2_exp):.4f}")
    n_alg_wins = sum(1 for a, e in zip(r2_alg, r2_exp) if a > e)
    print(f"    algebraic 우세: {n_alg_wins}/{N_ZEROS}")

    # Step 3: FE 대칭 검증
    print("\n" + "=" * 80)
    print("[Step 3] FE 대칭 검증: κ(σ) vs κ(1-σ)")
    print("=" * 80)

    # σ < 0.5인 값과 1-σ > 0.5인 값의 쌍
    sigma_pairs = []
    for s in SIGMAS:
        if s < 0.5:
            s_mirror = 1.0 - s
            if s_mirror in SIGMAS or any(abs(s_mirror - ss) < 1e-6 for ss in SIGMAS):
                # 미러 σ가 목록에 있는 경우
                idx_mirror = min(range(len(SIGMAS)), key=lambda k: abs(SIGMAS[k] - s_mirror))
                sigma_pairs.append((SIGMAS.index(s), idx_mirror, s, SIGMAS[idx_mirror]))

    print(f"\n  대칭 쌍 {len(sigma_pairs)}개:")

    symmetry_errors = []
    for j, t_n in enumerate(zeros_t):
        kappa_col = kappa_matrix[:, j]
        kappa_half = kappa_col[SIGMAS.index(0.5)]

        for idx_l, idx_r, s_l, s_r in sigma_pairs:
            k_l = kappa_col[idx_l]
            k_r = kappa_col[idx_r]

            if k_l > 0 and k_r > 0 and kappa_half > 0:
                rel_err = abs(k_l - k_r) / kappa_half
                symmetry_errors.append((j+1, s_l, s_r, k_l, k_r, rel_err))

    if symmetry_errors:
        print(f"\n  {'영점':>4} | {'σ_L':>6} {'σ_R':>6} | {'κ_L':>12} {'κ_R':>12} | {'|Δ|/κ½':>10}")
        print("  " + "-" * 65)
        for zn, sl, sr, kl, kr, re in symmetry_errors:
            print(f"  {zn:4d} | {sl:6.4f} {sr:6.4f} | {kl:12.4e} {kr:12.4e} | {re:10.2e}")

        all_rel = [e[5] for e in symmetry_errors]
        print(f"\n  대칭 오차 통계:")
        print(f"    mean |κ(σ) - κ(1-σ)| / κ(1/2) = {np.mean(all_rel):.2e}")
        print(f"    max                             = {np.max(all_rel):.2e}")
        fe_pass = np.max(all_rel) < 0.01
        print(f"    FE 대칭 판정 (< 1%): {'✅ PASS' if fe_pass else '❌ FAIL'}")
    else:
        print("  ⚠️ 대칭 쌍 없음 — σ 목록에 (σ, 1-σ) 쌍이 부족")
        fe_pass = None

    # Step 4: 종합 판별
    print("\n" + "=" * 80)
    print("[Step 4] 종합 판별")
    print("=" * 80)

    mean_r2_alg = np.mean(r2_alg)
    mean_r2_exp = np.mean(r2_exp)

    if mean_r2_alg > mean_r2_exp + 0.05:
        peak_type = "algebraic singularity"
        confidence = "high" if mean_r2_alg > 0.95 else "moderate"
    elif mean_r2_exp > mean_r2_alg + 0.05:
        peak_type = "essential singularity (exponential)"
        confidence = "high" if mean_r2_exp > 0.95 else "moderate"
    else:
        peak_type = "undetermined (algebraic ≈ exponential)"
        confidence = "low"

    alpha_cv = np.std(alphas_arr) / abs(np.mean(alphas_arr)) * 100 if abs(np.mean(alphas_arr)) > 0 else float('inf')

    print(f"\n  Peak 유형: {peak_type}")
    print(f"  신뢰도: {confidence}")
    print(f"  α (algebraic order): {np.mean(alphas_arr):.4f} ± {np.std(alphas_arr):.4f} (CV={alpha_cv:.1f}%)")
    print(f"  R²_alg = {mean_r2_alg:.4f}, R²_exp = {mean_r2_exp:.4f}")
    if fe_pass is not None:
        print(f"  FE 대칭: {'✅ 확인' if fe_pass else '❌ 위반'}")

    elapsed = time.time() - start_time
    print(f"\n  총 소요 시간: {elapsed:.1f}초")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 결과 저장
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

    with open(RESULT_PATH, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("[사이클 #254] σ-국소화 해석적 구조 탐사 결과\n")
        f.write(f"생성: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"소요: {elapsed:.1f}초, dps={mpmath.mp.dps}\n")
        f.write("=" * 80 + "\n\n")

        # 설정
        f.write("설정:\n")
        f.write(f"  σ 값 ({len(SIGMAS)}개): {SIGMAS}\n")
        f.write(f"  영점 수: {N_ZEROS}\n")
        f.write(f"  총 측정: {len(SIGMAS) * N_ZEROS}\n\n")

        # 원시 데이터: κ 행렬
        f.write("=" * 80 + "\n")
        f.write("원시 데이터: κ(σ, t_n)\n")
        f.write("=" * 80 + "\n\n")

        # 헤더
        f.write(f"{'σ':>8}")
        for j in range(N_ZEROS):
            f.write(f" | {'t_'+str(j+1)+'='+str(round(zeros_t[j],2)):>14}")
        f.write("\n")
        f.write("-" * (8 + 17 * N_ZEROS) + "\n")

        for i, sigma in enumerate(SIGMAS):
            f.write(f"{sigma:8.4f}")
            for j in range(N_ZEROS):
                k = kappa_matrix[i][j]
                if k == float('inf'):
                    f.write(f" | {'∞':>14}")
                elif np.isnan(k):
                    f.write(f" | {'NaN':>14}")
                else:
                    f.write(f" | {k:14.4e}")
            f.write("\n")

        # 영점별 분석
        f.write("\n" + "=" * 80 + "\n")
        f.write("영점별 발산 유형 분석\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"{'영점':>4} {'t_n':>10} {'κ(1/2)':>14} | {'α':>8} {'R²_alg':>8} | {'γ':>10} {'R²_exp':>8} | {'우세':>12}\n")
        f.write("-" * 90 + "\n")

        for j in range(N_ZEROS):
            kh = kappa_matrix[SIGMAS.index(0.5), j]
            a = alphas[j]
            ra = r2_alg[j]
            g = gammas[j]
            re = r2_exp[j]
            winner = "algebraic" if ra > re else "exponential"
            f.write(f"{j+1:4d} {zeros_t[j]:10.4f} {kh:14.4e} | {a:8.4f} {ra:8.4f} | {g:10.6f} {re:8.4f} | {winner:>12}\n")

        # 통계
        f.write(f"\n통계:\n")
        f.write(f"  α (algebraic order): mean={np.mean(alphas_arr):.4f}, std={np.std(alphas_arr):.4f}, CV={alpha_cv:.1f}%\n")
        f.write(f"  γ (exponential coeff): mean={np.mean(gammas_arr):.6f}, std={np.std(gammas_arr):.6f}\n")
        f.write(f"  R²_alg: mean={mean_r2_alg:.4f}\n")
        f.write(f"  R²_exp: mean={mean_r2_exp:.4f}\n")
        f.write(f"  algebraic 우세: {n_alg_wins}/{N_ZEROS}\n\n")

        # FE 대칭
        f.write("=" * 80 + "\n")
        f.write("FE 대칭 검증\n")
        f.write("=" * 80 + "\n\n")

        if symmetry_errors:
            for zn, sl, sr, kl, kr, rel_e in symmetry_errors:
                f.write(f"  영점 #{zn}: σ={sl:.4f} vs σ={sr:.4f}, κ_L={kl:.6e}, κ_R={kr:.6e}, |Δ|/κ½={rel_e:.2e}\n")
            f.write(f"\n  mean |Δ|/κ½ = {np.mean(all_rel):.2e}\n")
            f.write(f"  max  |Δ|/κ½ = {np.max(all_rel):.2e}\n")
            f.write(f"  FE 대칭 판정: {'PASS' if fe_pass else 'FAIL'}\n\n")

        # 종합 판별
        f.write("=" * 80 + "\n")
        f.write("종합 판별\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"  Peak 유형: {peak_type}\n")
        f.write(f"  신뢰도: {confidence}\n")
        f.write(f"  근거:\n")
        f.write(f"    - R²_alg ({mean_r2_alg:.4f}) vs R²_exp ({mean_r2_exp:.4f})\n")
        f.write(f"    - α = {np.mean(alphas_arr):.4f} ± {np.std(alphas_arr):.4f} (CV={alpha_cv:.1f}%)\n")
        f.write(f"    - algebraic 우세: {n_alg_wins}/{N_ZEROS} 영점\n")
        if fe_pass is not None:
            f.write(f"    - FE 대칭: {'확인' if fe_pass else '위반'} (max |Δ|/κ½ = {np.max(all_rel):.2e})\n")

        # 좌우 비대칭 상세 (σ < 0.5 vs σ > 0.5)
        f.write("\n" + "=" * 80 + "\n")
        f.write("좌우 프로파일 상세\n")
        f.write("=" * 80 + "\n\n")

        for j in range(min(5, N_ZEROS)):  # 처음 5개 영점만 상세
            f.write(f"  영점 #{j+1} (t={zeros_t[j]:.4f}):\n")
            f.write(f"    {'σ':>8} {'|σ-½|':>8} {'κ':>14} {'log κ':>10}\n")
            for i, sigma in enumerate(SIGMAS):
                k = kappa_matrix[i][j]
                ds = abs(sigma - 0.5)
                if k > 0 and k != float('inf') and not np.isnan(k):
                    f.write(f"    {sigma:8.4f} {ds:8.4f} {k:14.4e} {np.log(k):10.4f}\n")
                elif k == float('inf'):
                    f.write(f"    {sigma:8.4f} {ds:8.4f} {'∞':>14} {'∞':>10}\n")
            f.write("\n")

    print(f"\n결과 저장: {RESULT_PATH}")
    print("완료.")


if __name__ == '__main__':
    main()
