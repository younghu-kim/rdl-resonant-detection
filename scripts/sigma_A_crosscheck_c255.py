"""
=============================================================================
[Project RDL] 사이클 #255 — σ-프로파일 ↔ A 공식 교차검증
=============================================================================

목적:
  C-254에서 관측된 σ-프로파일 부차항 c(γ)가 Cor 4.2의 진폭 계수
  A(γ) = Im(c₀)² + 2c₁과 정확히 일치하는지 교차검증.

  이론적 근거 (Thm 1 + Thm 5의 따름정리):
    κ(σ,γ) = |ξ'/ξ(σ+iγ)|²
           = 1/(σ-1/2)² + A(γ) + O((σ-1/2)²)
    여기서 A = Im(c₀)² + 2c₁, 1/ε 항은 Re(c₀)=0에 의해 소거.

방법:
  1. ζ(s) 처음 20개 영점에서 c₀, c₁ (Laurent 계수) 직접 계산
  2. A(γ) = Im(c₀)² + 2c₁ 산출
  3. C-254 데이터에서 c(γ) = κ(0.3, γ) - 1/(0.2)² 추출
  4. ρ(c(γ), A(γ)) 상관계수 및 절대 일치도 측정
  5. 추가: O(ε²) 보정항 검증 — c(γ,ε) = A + B·ε² 피팅

결과 파일: results/sigma_A_crosscheck_c255.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

from bundle_utils import xi_func

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

mpmath.mp.dps = 60

N_ZEROS = 20

# C-254와 동일한 σ값 (좌측만, ε 추출용)
SIGMAS_LEFT = [0.3, 0.35, 0.4, 0.45, 0.48, 0.49, 0.495, 0.499, 0.4999]
EPSILONS = [0.5 - s for s in SIGMAS_LEFT]  # [0.2, 0.15, 0.1, 0.05, ...]

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/sigma_A_crosscheck_c255.txt'
)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Laurent 계수 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_laurent_coeffs(t_zero, n_coeffs=2):
    """
    ξ'/ξ(ρ + u) = 1/u + c₀ + c₁u + c₂u² + ... 의 c₀, c₁ 계산.

    방법: 폐곡선 적분 (Cauchy 적분 공식)
      cₙ = (1/(2πi)) ∮ [L(ρ+u) - 1/u] u^{-(n+1)} du

    작은 원 |u| = r 위에서 수치 적분.
    """
    rho = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_zero)))
    r = mpmath.mpf('0.01')  # 작은 반지름
    N_pts = 128  # 적분점 수

    coeffs = []
    for n in range(n_coeffs):
        # cₙ = (1/(2πi)) ∮ g(u) u^{-(n+1)} du  where g(u) = L(ρ+u) - 1/u
        integral = mpmath.mpc(0, 0)
        for k in range(N_pts):
            theta = 2 * mpmath.pi * k / N_pts
            u = r * mpmath.exp(1j * theta)
            s = rho + u

            # L(s) = ξ'/ξ(s) — 수치 미분
            h = mpmath.mpf(1) / mpmath.mpf(10**18)
            xi_val = xi_func(s)
            xi_plus = xi_func(s + h)
            xi_minus = xi_func(s - h)
            L_val = (xi_plus - xi_minus) / (2 * h * xi_val)

            # g(u) = L(ρ+u) - 1/u (정칙 부분)
            g_val = L_val - 1 / u

            # 피적분함수: g(u) · u^{-(n+1)} · du/dθ
            # du = i·r·e^{iθ} dθ = i·u dθ
            integrand = g_val * u**(-(n)) * (1j * u)  # g · u^{-n-1} · i·u = g · u^{-n} · i

            integral += integrand

        # (1/(2πi)) × (2π/N)
        cn = integral / (N_pts)  # (Σ f(θ) Δθ)/(2πi) = (Σ f · 2π/N)/(2πi) = Σ f / (N·i)
        # 보정: integral에 이미 i가 곱해져 있으므로
        cn = integral / (N_pts * 1j)
        # 아니다, 다시 정리:
        # ∮ g(u) u^{-(n+1)} du ≈ Σ g(u_k) u_k^{-(n+1)} Δu_k
        # Δu_k = i u_k Δθ = i u_k (2π/N)
        # 따라서 ∮ ≈ Σ g u_k^{-(n+1)} · i u_k · (2π/N)
        #        = (2πi/N) Σ g u_k^{-n}
        # cₙ = (1/(2πi)) · (2πi/N) Σ g u_k^{-n} = (1/N) Σ g u_k^{-n}

        # 재계산
        cn_sum = mpmath.mpc(0, 0)
        for k in range(N_pts):
            theta = 2 * mpmath.pi * k / N_pts
            u = r * mpmath.exp(1j * theta)
            s = rho + u

            h = mpmath.mpf(1) / mpmath.mpf(10**18)
            xi_val = xi_func(s)
            xi_plus = xi_func(s + h)
            xi_minus = xi_func(s - h)
            L_val = (xi_plus - xi_minus) / (2 * h * xi_val)

            g_val = L_val - 1 / u
            cn_sum += g_val * u**(-n)

        cn = cn_sum / N_pts
        coeffs.append(complex(cn))

    return coeffs


def kappa_at_sigma(sigma, t_zero):
    """σ+i·t_zero에서 κ = |ξ'/ξ|² 계산"""
    s = mpmath.mpc(mpmath.mpf(str(sigma)), mpmath.mpf(str(t_zero)))
    h = mpmath.mpf(1) / mpmath.mpf(10**20)

    xi_val = xi_func(s)
    if abs(xi_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 15):
        return float('inf')

    xi_plus = xi_func(s + h)
    xi_minus = xi_func(s - h)
    L = (xi_plus - xi_minus) / (2 * h * xi_val)

    return float(abs(L)**2)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    start_time = time.time()

    print("=" * 80)
    print("[사이클 #255] σ-프로파일 ↔ A 공식 교차검증")
    print(f"  dps = {mpmath.mp.dps}, 영점 = {N_ZEROS}개")
    print("=" * 80)

    # Step 0: 영점 수집
    print("\n[Step 0] ζ 영점 수집...")
    zeros_t = []
    for n in range(1, N_ZEROS + 1):
        t = float(mpmath.zetazero(n).imag)
        zeros_t.append(t)
    print(f"  {N_ZEROS}개 수집 완료: t ∈ [{zeros_t[0]:.2f}, {zeros_t[-1]:.2f}]")

    # Step 1: Laurent 계수 c₀, c₁ 계산
    print("\n[Step 1] Laurent 계수 c₀, c₁ 계산 (Cauchy 적분)...")

    c0_list = []
    c1_list = []
    A_list = []  # A = Im(c₀)² + 2·Re(c₁)

    for j, t_n in enumerate(zeros_t):
        print(f"  영점 #{j+1}/{N_ZEROS}: t = {t_n:.4f} ... ", end='', flush=True)
        coeffs = compute_laurent_coeffs(t_n, n_coeffs=2)
        c0 = coeffs[0]
        c1 = coeffs[1]
        c0_list.append(c0)
        c1_list.append(c1)

        A = c0.imag**2 + 2 * c1.real
        A_list.append(A)

        print(f"c₀ = {c0.real:+.6f}{c0.imag:+.6f}i, "
              f"c₁ = {c1.real:+.6f}{c1.imag:+.6f}i, "
              f"A = {A:.6f}")

    # Thm 5 검증: Re(c₀) ≈ 0, Im(c₁) ≈ 0
    re_c0 = [abs(c.real) for c in c0_list]
    im_c0 = [abs(c.imag) for c in c0_list]
    re_c1 = [abs(c.real) for c in c1_list]
    im_c1 = [abs(c.imag) for c in c1_list]

    print(f"\n  [Thm 5 검증]")
    print(f"  max |Re(c₀)| / |Im(c₀)| = {max(r/i if i > 0 else 0 for r, i in zip(re_c0, im_c0)):.2e}")
    print(f"  max |Im(c₁)| / |Re(c₁)| = {max(i/r if r > 0 else 0 for r, i in zip(im_c1, re_c1)):.2e}")

    # Step 2: σ-프로파일에서 c(γ) 추출
    print("\n[Step 2] σ-프로파일에서 c(γ) 추출...")

    # 각 영점, 여러 ε에서 c(γ,ε) = κ(σ,γ) - 1/ε² 계산
    c_from_sigma = np.zeros((N_ZEROS, len(EPSILONS)))

    for j, t_n in enumerate(zeros_t):
        print(f"  영점 #{j+1}/{N_ZEROS}: t = {t_n:.4f}", flush=True)
        for i, eps in enumerate(EPSILONS):
            sigma = 0.5 - eps  # 좌측
            kappa = kappa_at_sigma(sigma, t_n)
            c_val = kappa - 1.0 / eps**2
            c_from_sigma[j, i] = c_val

    # 각 영점의 c(γ) 평균 (ε에 무관해야 함)
    c_mean = np.mean(c_from_sigma, axis=1)
    c_std = np.std(c_from_sigma, axis=1)

    print(f"\n  c(γ) 통계 (ε에 걸친 CV):")
    for j, t_n in enumerate(zeros_t):
        cv = c_std[j] / abs(c_mean[j]) * 100 if abs(c_mean[j]) > 1e-10 else float('inf')
        print(f"    #{j+1} t={t_n:.4f}: c = {c_mean[j]:.6f} ± {c_std[j]:.6f} (CV={cv:.2f}%)")

    # Step 3: 교차검증 — c(γ) vs A(γ)
    print("\n" + "=" * 80)
    print("[Step 3] 교차검증: c(γ) vs A(γ)")
    print("=" * 80)

    A_arr = np.array(A_list)

    print(f"\n  {'#':>3} {'t':>10} {'c(γ) σ-scan':>14} {'A(γ) Laurent':>14} {'|Δ|':>10} {'|Δ|/A':>10}")
    print("  " + "-" * 68)

    deltas = []
    rel_deltas = []
    for j in range(N_ZEROS):
        c_val = c_mean[j]
        a_val = A_arr[j]
        delta = abs(c_val - a_val)
        rel_delta = delta / abs(a_val) if abs(a_val) > 1e-10 else float('inf')
        deltas.append(delta)
        rel_deltas.append(rel_delta)
        print(f"  {j+1:3d} {zeros_t[j]:10.4f} {c_val:14.6f} {a_val:14.6f} {delta:10.6f} {rel_delta:10.4e}")

    # 상관계수
    corr = np.corrcoef(c_mean, A_arr)[0, 1]
    mean_rel = np.mean(rel_deltas)
    max_rel = np.max(rel_deltas)

    print(f"\n  ρ(c, A) = {corr:.10f}")
    print(f"  mean |Δ|/|A| = {mean_rel:.4e}")
    print(f"  max  |Δ|/|A| = {max_rel:.4e}")

    # Step 4: O(ε²) 보정항 검증
    print("\n" + "=" * 80)
    print("[Step 4] O(ε²) 보정항 검증")
    print("=" * 80)

    # c(γ,ε) = A + B·ε² 피팅
    for j in range(min(5, N_ZEROS)):
        eps_arr = np.array(EPSILONS)
        c_arr = c_from_sigma[j, :]

        # 선형 피팅: c(ε) = a + b·ε²
        X = np.vstack([np.ones(len(eps_arr)), eps_arr**2]).T
        result = np.linalg.lstsq(X, c_arr, rcond=None)
        a_fit, b_fit = result[0]

        y_pred = a_fit + b_fit * eps_arr**2
        ss_res = np.sum((c_arr - y_pred)**2)
        ss_tot = np.sum((c_arr - np.mean(c_arr))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0

        print(f"\n  영점 #{j+1} (t={zeros_t[j]:.4f}):")
        print(f"    c(ε) = {a_fit:.6f} + {b_fit:.4f}·ε²")
        print(f"    R² = {r2:.8f}")
        print(f"    A(Laurent) = {A_list[j]:.6f}")
        print(f"    |a_fit - A| = {abs(a_fit - A_list[j]):.6e}")

    # Step 5: 종합 판정
    print("\n" + "=" * 80)
    print("[Step 5] 종합 판정")
    print("=" * 80)

    pass_corr = corr > 0.999
    pass_abs = max_rel < 0.01

    status = "양성" if pass_corr and pass_abs else "음성"

    print(f"\n  c(γ) = A(γ) 교차검증: {'✅ PASS' if pass_corr and pass_abs else '❌ FAIL'}")
    print(f"    상관: ρ = {corr:.10f} {'✅' if pass_corr else '❌'} (기준 >0.999)")
    print(f"    절대: max|Δ|/|A| = {max_rel:.4e} {'✅' if pass_abs else '❌'} (기준 <1%)")
    print(f"  판정: {status}")

    if pass_corr and pass_abs:
        print(f"\n  의의: σ-프로파일(C-254)과 A 공식(Cor 4.2)이 정확히 일치.")
        print(f"  이는 Thm 1 (접속 반대칭)과 Thm 5 (Laurent 패리티)의 ")
        print(f"  따름정리로서 이론적으로 예측된 결과를 수치적으로 확인.")
        print(f"  κ(σ,γ) = 1/(σ-1/2)² + A(γ) + O((σ-1/2)²)")
        print(f"  where A = Im(c₀)² + 2c₁, Re(c₀)=0 by Thm 5.")

    elapsed = time.time() - start_time
    print(f"\n  소요 시간: {elapsed:.1f}초")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 결과 저장
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

    with open(RESULT_PATH, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("[사이클 #255] σ-프로파일 ↔ A 공식 교차검증 결과\n")
        f.write(f"생성: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"소요: {elapsed:.1f}초, dps={mpmath.mp.dps}\n")
        f.write("=" * 80 + "\n\n")

        f.write("이론:\n")
        f.write("  Thm 1 (접속 반대칭): L(1-s) = -L(s)\n")
        f.write("  Thm 5 (Laurent 패리티): Re(c₀)=0, Im(c₁)=0\n")
        f.write("  따름정리: κ(σ,γ) = 1/(σ-1/2)² + A(γ) + O((σ-1/2)²)\n")
        f.write("    where A(γ) = Im(c₀)² + 2c₁ (Cor 4.2)\n")
        f.write("    1/(σ-1/2) 항은 Re(c₀)=0에 의해 정확히 소거\n\n")

        f.write("=" * 80 + "\n")
        f.write("Laurent 계수\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"  {'#':>3} {'t':>10} {'Re(c₀)':>12} {'Im(c₀)':>12} {'Re(c₁)':>12} {'Im(c₁)':>12} {'A':>12}\n")
        f.write("  " + "-" * 78 + "\n")
        for j in range(N_ZEROS):
            c0 = c0_list[j]
            c1 = c1_list[j]
            f.write(f"  {j+1:3d} {zeros_t[j]:10.4f} {c0.real:12.6f} {c0.imag:12.6f} "
                    f"{c1.real:12.6f} {c1.imag:12.6f} {A_list[j]:12.6f}\n")

        f.write(f"\n  Thm 5 검증:\n")
        f.write(f"    max |Re(c₀)|/|Im(c₀)| = {max(r/i if i > 0 else 0 for r, i in zip(re_c0, im_c0)):.2e}\n")
        f.write(f"    max |Im(c₁)|/|Re(c₁)| = {max(i/r if r > 0 else 0 for r, i in zip(im_c1, re_c1)):.2e}\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("교차검증: c(γ) vs A(γ)\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"  {'#':>3} {'t':>10} {'c(γ) σ-scan':>14} {'A(γ) Laurent':>14} {'|Δ|':>10} {'|Δ|/A':>10}\n")
        f.write("  " + "-" * 68 + "\n")
        for j in range(N_ZEROS):
            f.write(f"  {j+1:3d} {zeros_t[j]:10.4f} {c_mean[j]:14.6f} {A_arr[j]:14.6f} "
                    f"{deltas[j]:10.6f} {rel_deltas[j]:10.4e}\n")

        f.write(f"\n  ρ(c, A) = {corr:.10f}\n")
        f.write(f"  mean |Δ|/|A| = {mean_rel:.4e}\n")
        f.write(f"  max  |Δ|/|A| = {max_rel:.4e}\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("O(ε²) 보정항 피팅\n")
        f.write("=" * 80 + "\n\n")

        for j in range(min(5, N_ZEROS)):
            eps_arr = np.array(EPSILONS)
            c_arr = c_from_sigma[j, :]
            X = np.vstack([np.ones(len(eps_arr)), eps_arr**2]).T
            result = np.linalg.lstsq(X, c_arr, rcond=None)
            a_fit, b_fit = result[0]
            f.write(f"  영점 #{j+1} (t={zeros_t[j]:.4f}): c(ε) = {a_fit:.6f} + {b_fit:.4f}·ε²\n")
            f.write(f"    A(Laurent) = {A_list[j]:.6f}, |a_fit - A| = {abs(a_fit - A_list[j]):.6e}\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("종합 판정\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"  c(γ) = A(γ) 교차검증: {'PASS' if pass_corr and pass_abs else 'FAIL'}\n")
        f.write(f"    상관: ρ = {corr:.10f} (기준 >0.999)\n")
        f.write(f"    절대: max|Δ|/|A| = {max_rel:.4e} (기준 <1%)\n")
        f.write(f"  판정: {status}\n")

    print(f"\n결과 저장: {RESULT_PATH}")
    print("완료.")


if __name__ == '__main__':
    main()
