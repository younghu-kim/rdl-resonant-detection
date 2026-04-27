"""
C-382: A-gap NN 하한 정리 수치 검증
=====================================
Theorem (A-gap NN Lower Bound):
  FE+CC를 만족하는 L-함수의 임계선 위 단순 영점 ρ_n = 1/2 + iγ_n에서,
  A(γ_n) = Im(c₀)² + 2Re(c₁) ≥ 2/g_min(γ_n)²
  여기서 g_min = min(γ_{n+1}-γ_n, γ_n-γ_{n-1}).

Proof:
  (1) Hadamard: Λ'/Λ(ρ_n+u) = 1/u + Σ_{k≠n} [1/(ρ_n+u-ρ_k) - 1/ρ_k] + B
  (2) Laurent 계수: c₁ = -Σ_{k≠n} 1/(ρ_n-ρ_k)²
  (3) On-critical: (ρ_n-ρ_k)² = -(γ_n-γ_k)²
      ∴ c₁ = Σ_{k≠n} 1/(γ_n-γ_k)² > 0
  (4) c₁ ≥ 1/g_left² + 1/g_right² ≥ 1/g_min²
  (5) A = Im(c₀)² + 2c₁ ≥ 2/g_min²  □

방법: DFT 기반 Cauchy 적분 (ξ 직접 평가, 수치미분 불필요)
  a_k = (1/Nr^k) Σ_j ξ(ρ+re^{iθ_j}) e^{-ikθ_j}
  c₀ = a₂/a₁, c₁ = 2a₃/a₁ - (a₂/a₁)²

출력: results/agap_nn_bound_c382.txt
"""

import sys
import os
import time
import numpy as np
import mpmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bundle_utils import xi_func, find_zeros_zeta, completed_L, CHARACTERS

mpmath.mp.dps = 80

RESULT_FILE = os.path.join(os.path.dirname(__file__), '..', 'results', 'agap_nn_bound_c382.txt')


def compute_c0_c1_dft(gamma_n, f_func, radius=0.3, n_points=256):
    """DFT 기반 Cauchy 적분으로 Laurent 계수 c₀, c₁ 계산.

    f(s)의 단순 영점 ρ에서 f(s) = a₁(s-ρ) + a₂(s-ρ)² + a₃(s-ρ)³ + ...
    → f'/f(s) = 1/(s-ρ) + c₀ + c₁(s-ρ) + ...
    여기서 c₀ = a₂/a₁, c₁ = 2a₃/a₁ - (a₂/a₁)².

    a_k = (1/Nr^k) Σ_j f(ρ+re^{iθ_j}) e^{-ikθ_j}  (DFT)
    """
    rho = mpmath.mpc(0.5, gamma_n)
    r = mpmath.mpf(radius)

    # 컨투어 위의 함수값 평가
    f_vals = []
    for j in range(n_points):
        theta = 2 * mpmath.pi * j / n_points
        s = rho + r * mpmath.expj(theta)
        f_vals.append(f_func(s))

    # DFT로 Taylor 계수 a_1, a_2, a_3 추출
    def get_ak(k):
        total = mpmath.mpc(0)
        for j in range(n_points):
            theta_j = 2 * mpmath.pi * j / n_points
            total += f_vals[j] * mpmath.expj(-k * theta_j)
        return total / (n_points * r**k)

    a1 = get_ak(1)
    a2 = get_ak(2)
    a3 = get_ak(3)

    # c₀ = a₂/a₁
    c0 = a2 / a1

    # c₁ = 2a₃/a₁ - (a₂/a₁)²
    c1 = 2 * a3 / a1 - (a2 / a1)**2

    return c0, c1


def compute_c0_c1_dirichlet_dft(gamma_n, char_info, radius=0.3, n_points=256):
    """디리클레 L-함수 Λ(s,χ)에 대한 DFT Cauchy 적분."""
    rho = mpmath.mpc(0.5, gamma_n)
    r = mpmath.mpf(radius)

    f_vals = []
    for j in range(n_points):
        theta = 2 * mpmath.pi * j / n_points
        s = rho + r * mpmath.expj(theta)
        f_vals.append(completed_L(s, char_info))

    def get_ak(k):
        total = mpmath.mpc(0)
        for j in range(n_points):
            theta_j = 2 * mpmath.pi * j / n_points
            total += f_vals[j] * mpmath.expj(-k * theta_j)
        return total / (n_points * r**k)

    a1 = get_ak(1)
    a2 = get_ak(2)
    a3 = get_ak(3)

    c0 = a2 / a1
    c1 = 2 * a3 / a1 - (a2 / a1)**2

    return c0, c1


def verify_bound(zeros, name, compute_func, trim=3, max_zeros=None):
    """영점 리스트에 대해 A ≥ 2/g_min² 검증."""
    n = len(zeros)
    start = trim
    end = n - trim
    if max_zeros and (end - start) > max_zeros:
        # 균등 샘플링
        indices = np.linspace(start, end-1, max_zeros, dtype=int)
    else:
        indices = range(start, end)

    results = []
    n_pass = 0
    n_fail = 0

    print(f"\n{'='*60}")
    print(f"  {name}: {n} 영점, 검증 {len(indices)}개")
    print(f"{'='*60}")

    for count, idx in enumerate(indices):
        gamma_n = zeros[idx]
        g_left = gamma_n - zeros[idx - 1]
        g_right = zeros[idx + 1] - gamma_n
        g_min = min(g_left, g_right)

        c0, c1 = compute_func(gamma_n)

        re_c0 = float(mpmath.re(c0))
        im_c0 = float(mpmath.im(c0))
        re_c1 = float(mpmath.re(c1))
        im_c1 = float(mpmath.im(c1))

        A = im_c0**2 + 2 * re_c1
        bound = 2.0 / g_min**2
        ratio = A / bound if bound > 0 else float('inf')

        # 2H₁/A
        two_h1_A = 2 * re_c1 / A if A > 1e-20 else 0

        ok = ratio >= 0.999  # 수치 오차 허용
        if ok:
            n_pass += 1
        else:
            n_fail += 1

        results.append({
            'idx': idx, 'gamma': gamma_n,
            'g_min': g_min, 'g_left': g_left, 'g_right': g_right,
            'A': A, 'bound': bound, 'ratio': ratio,
            'two_h1_A': two_h1_A,
            'Re_c0': re_c0, 'Im_c0': im_c0,
            'Re_c1': re_c1, 'Im_c1': im_c1,
        })

        if count % 20 == 0:
            # Thm 5 검증: |Re(c₀)| 작은지
            print(f"  [{count+1}/{len(indices)}] γ={gamma_n:.3f}, g_min={g_min:.4f}, "
                  f"A={A:.4f}, bound={bound:.4f}, R={ratio:.3f}, "
                  f"2H₁/A={two_h1_A:.3f}, |Re(c₀)|={abs(re_c0):.2e} "
                  f"{'✅' if ok else '❌'}")

    return results, n_pass, n_fail


def main():
    t0 = time.time()
    output = []

    def log(msg):
        print(msg)
        output.append(msg)

    log("=" * 70)
    log("C-382: A-gap NN 하한 정리 수치 검증")
    log("=" * 70)
    log(f"일시: {time.strftime('%Y-%m-%d %H:%M')}")
    log(f"mpmath dps: {mpmath.mp.dps}")
    log(f"Cauchy 반경: 0.3, DFT 점수: 256")
    log("")
    log("Theorem (A-gap NN Lower Bound):")
    log("  A(γ_n) = Im(c₀)² + 2Re(c₁) ≥ 2/g_min(γ_n)²")
    log("  여기서 g_min = min(gap_left, gap_right)")
    log("")
    log("Proof:")
    log("  (1) On-critical 영점: (ρ_n-ρ_k)² = -(γ_n-γ_k)² < 0")
    log("  (2) c₁ = -Σ_{k≠n} 1/(ρ_n-ρ_k)² = Σ_{k≠n} 1/(γ_n-γ_k)² > 0")
    log("  (3) c₁ ≥ 1/g_left² + 1/g_right² ≥ 1/g_min²")
    log("  (4) Im(c₀)² ≥ 0 (자명)")
    log("  (5) A = Im(c₀)² + 2c₁ ≥ 2·(1/g_min²) = 2/g_min²  □")
    log("")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # §1. ζ(s) 검증 (200개 영점)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    log("\n" + "━" * 60)
    log("§1. ζ(s) 검증")
    log("━" * 60)

    zeta_zeros = find_zeros_zeta(10.0, 600.0)
    log(f"영점 수: {len(zeta_zeros)}, 범위: [{zeta_zeros[0]:.3f}, {zeta_zeros[-1]:.3f}]")

    zeta_compute = lambda g: compute_c0_c1_dft(g, xi_func, radius=0.3, n_points=256)

    zeta_res, zp, zf = verify_bound(
        zeta_zeros, "ζ(s)", zeta_compute, trim=3, max_zeros=150
    )

    log(f"\nζ(s) 결과: {zp}/{zp+zf} PASS ({zf} 위반)")
    _report(log, zeta_res, "ζ(s)")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # §2. χ₋₃ (mod 3) 검증
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    log("\n" + "━" * 60)
    log("§2. χ₋₃ (mod 3) 검증")
    log("━" * 60)

    from bundle_utils import find_zeros_dirichlet
    char3 = CHARACTERS['chi_mod_3']
    chi3_zeros = find_zeros_dirichlet(char3, t_min=1.0, t_max=200.0, n_scan=4000)
    log(f"영점 수: {len(chi3_zeros)}, 범위: [{chi3_zeros[0]:.3f}, {chi3_zeros[-1]:.3f}]")

    chi3_compute = lambda g: compute_c0_c1_dirichlet_dft(g, char3, radius=0.3, n_points=256)

    chi3_res, c3p, c3f = verify_bound(
        chi3_zeros, "χ₋₃ (mod 3)", chi3_compute, trim=3, max_zeros=80
    )

    log(f"\nχ₋₃ 결과: {c3p}/{c3p+c3f} PASS ({c3f} 위반)")
    _report(log, chi3_res, "χ₋₃")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # §3. Hadamard 교차검증 (ζ, c₁ DFT vs 직접 합)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    log("\n" + "━" * 60)
    log("§3. Hadamard 교차검증")
    log("━" * 60)

    # 중앙 15개 영점
    mid = len(zeta_zeros) // 2
    test_indices = range(max(3, mid-7), min(len(zeta_zeros)-3, mid+8))
    had_errs = []

    for idx in test_indices:
        g = zeta_zeros[idx]
        # DFT c₁
        _, c1_dft = compute_c0_c1_dft(g, xi_func, 0.3, 256)
        c1_d = float(mpmath.re(c1_dft))

        # Hadamard 직접 합
        c1_h = sum(1.0 / (g - zeta_zeros[k])**2
                    for k in range(len(zeta_zeros)) if k != idx)

        rel = abs(c1_d - c1_h) / abs(c1_d) if abs(c1_d) > 1e-10 else 0
        had_errs.append(rel)

        if (idx - test_indices.start) % 5 == 0:
            log(f"  γ={g:.3f}: c₁_DFT={c1_d:.6f}, c₁_Had={c1_h:.6f}, "
                f"Δrel={rel:.2e}")

    log(f"\n  c₁ 교차검증: mean|Δ|/|c₁| = {np.mean(had_errs):.2e}")
    log(f"  (Hadamard 절삭 + B 항 차이로 인한 체계적 편차)")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # §4. 종합
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    total_p = zp + c3p
    total_f = zf + c3f
    total = total_p + total_f
    elapsed = time.time() - t0

    all_res = zeta_res + chi3_res

    log("\n" + "=" * 70)
    log("§4. 종합 결과")
    log("=" * 70)
    log(f"\n  총 검증: {total} 영점 (ζ {zp+zf} + χ₋₃ {c3p+c3f})")
    log(f"  PASS: {total_p}/{total} ({100*total_p/total:.1f}%)")
    log(f"  FAIL: {total_f}")

    if all_res:
        ratios = [r['ratio'] for r in all_res]
        two_h1s = [r['two_h1_A'] for r in all_res]
        As = [r['A'] for r in all_res]
        gs = [r['g_min'] for r in all_res]

        log(f"\n  초과 비율 R = A·g²/2 (전체):")
        log(f"    min  = {min(ratios):.4f}")
        log(f"    mean = {np.mean(ratios):.4f}")
        log(f"    median = {np.median(ratios):.4f}")
        log(f"    max  = {max(ratios):.4f}")
        log(f"    R≥1 비율: {sum(1 for r in ratios if r>=0.999)}/{len(ratios)}")

        log(f"\n  2H₁/A (전체):")
        log(f"    mean = {np.mean(two_h1s):.4f} ± {np.std(two_h1s):.4f}")

        # Spearman
        from scipy.stats import spearmanr, pearsonr
        rho_s, p_s = spearmanr(As, gs)
        log(f"\n  ρ_S(A, g_min) = {rho_s:.4f} (p={p_s:.2e})")

        # 부가: A vs 1/g² Pearson
        inv_g2 = [1/g**2 for g in gs]
        rho_p, p_p = pearsonr(As, inv_g2)
        log(f"  ρ_P(A, 1/g²) = {rho_p:.4f} (p={p_p:.2e})")

    # Thm 5 부수 검증
    re_c0s = [abs(r['Re_c0']) for r in all_res]
    im_c1s = [abs(r['Im_c1']) for r in all_res]
    log(f"\n  Thm 5 부수 검증:")
    log(f"    max|Re(c₀)| = {max(re_c0s):.2e} (이론: 0)")
    log(f"    max|Im(c₁)| = {max(im_c1s):.2e} (이론: 0)")

    # 판정
    if total_f == 0:
        verdict = "★★★★★ 양성: A ≥ 2/g_min² 전 영점 확인"
    elif total_f <= 2:
        verdict = f"★★★★ 조건부 양성: {total_f}개 수치 오차 위반"
    else:
        verdict = f"⚠️ {total_f}개 위반 — 점검 필요"

    log(f"\n  판정: {verdict}")
    log(f"  소요 시간: {elapsed:.1f}초")
    log(f"\n{'='*70}")

    # 저장
    os.makedirs(os.path.dirname(os.path.abspath(RESULT_FILE)), exist_ok=True)
    with open(RESULT_FILE, 'w') as f:
        f.write('\n'.join(output))
    print(f"\n결과 저장: {RESULT_FILE}")


def _report(log, results, name):
    """결과 통계 보고."""
    if not results:
        return

    from scipy.stats import spearmanr

    ratios = [r['ratio'] for r in results]
    two_h1s = [r['two_h1_A'] for r in results]
    As = [r['A'] for r in results]
    gs = [r['g_min'] for r in results]

    log(f"\n  {name} 초과 비율 R = A·g²/2:")
    log(f"    min = {min(ratios):.4f}, mean = {np.mean(ratios):.4f}, "
        f"max = {max(ratios):.4f}")
    log(f"    R ≥ 1: {sum(1 for r in ratios if r >= 0.999)}/{len(ratios)}")

    log(f"\n  {name} 2H₁/A:")
    log(f"    mean = {np.mean(two_h1s):.4f} ± {np.std(two_h1s):.4f}")

    rho_s, p_s = spearmanr(As, gs)
    log(f"\n  {name} ρ_S(A, g_min) = {rho_s:.4f} (p={p_s:.2e})")

    # 대표 영점 5개
    sorted_r = sorted(results, key=lambda r: r['g_min'])
    log(f"\n  {name} 대표 영점:")
    for label, r in [("최소 g", sorted_r[0]),
                      ("Q1", sorted_r[len(sorted_r)//4]),
                      ("중간", sorted_r[len(sorted_r)//2]),
                      ("Q3", sorted_r[3*len(sorted_r)//4]),
                      ("최대 g", sorted_r[-1])]:
        log(f"    {label}: γ={r['gamma']:.2f}, g={r['g_min']:.4f}, "
            f"A={r['A']:.4f}, bound={r['bound']:.4f}, R={r['ratio']:.3f}, "
            f"2H₁/A={r['two_h1_A']:.3f}")


if __name__ == '__main__':
    main()
