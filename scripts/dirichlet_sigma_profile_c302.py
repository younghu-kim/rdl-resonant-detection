"""
=============================================================================
[C-302] Dirichlet E(σ) 프로파일 — Prop 6 GL(1) 보편성 검증
=============================================================================
목적:
    Prop 6 (E(σ) = πN/|σ-1/2| + o) 가 ζ(s) 뿐 아니라 Dirichlet L-함수
    L(s,χ)에서도 성립하는지 검증.

방법:
    χ mod 3, 4, 5 각각에 대해:
    1. E(σ) = ∫|Λ'/Λ(σ+it,χ)|² dt 를 σ∈[0.1,0.9] 20점에서 측정
    2. α 멱법칙 피팅: E(σ) = A/|σ-1/2|^α → α≈1 확인
    3. 계수 확인: A ≈ πN_χ
    4. ζ(s) 결과(C-298)와 비교

판정 기준:
    양성: 3지표 모두 |α-1| < 0.05, A/πN ∈ [0.8, 1.2]
    중립: 일부 지표만 양성
    음성: α≠1 또는 A/πN 크게 벗어남
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))
from bundle_utils import (
    connection_dirichlet, find_zeros_dirichlet, CHARACTERS, completed_L
)

mpmath.mp.dps = 50  # t<40, conductor 작음 → dps=50 충분

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

T_RANGE = (10.0, 40.0)
N_T_POINTS = 300        # t 격자 (Δt=0.1) — 속도 우선, σ 추세 파악에 충분
SIGMA_SAFE_MIN = 0.05   # σ 최소 간격 (|σ-1/2| > 0.05)
N_SIGMA = 12            # σ 격자 수 — 6+6 (좌우 대칭)
FIT_SAFE = 0.03         # 피팅 안전 영역 |σ-1/2| > 0.03

OUT_FILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/dirichlet_sigma_profile_c302.txt'
)


def compute_energy_profile(char_info, t_min, t_max, n_t, sigma_values):
    """각 σ에서 E(σ) = ∫|Λ'/Λ(σ+it)|² dt 계산"""
    ts = np.linspace(t_min, t_max, n_t)
    dt = ts[1] - ts[0]
    energies = np.zeros(len(sigma_values))

    for j, sigma in enumerate(sigma_values):
        integrand = np.zeros(n_t)
        for i, t in enumerate(ts):
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            try:
                L_val = connection_dirichlet(s, char_info)
                integrand[i] = float(abs(L_val)**2)
            except Exception:
                integrand[i] = 0.0

        # 극값 클리핑 (영점 위 또는 매우 가까운 점)
        cap = np.percentile(integrand[integrand < 1e10], 99)
        integrand = np.minimum(integrand, 10 * cap)

        energies[j] = np.trapezoid(integrand, ts)

        if (j + 1) % 5 == 0 or j == 0:
            print(f"    σ={sigma:.3f}: E={energies[j]:.2f}")

    return energies


def fit_power_law(delta_sigma, energies, safe_min=0.02):
    """E = A/Δσ^α 피팅 (안전 영역만)"""
    mask = delta_sigma >= safe_min
    if mask.sum() < 3:
        return 0, 0, 0, 0

    log_x = np.log(delta_sigma[mask])
    log_y = np.log(energies[mask])

    # y = -α·x + log(A)
    A_mat = np.column_stack([np.ones(mask.sum()), log_x])
    coeffs = np.linalg.lstsq(A_mat, log_y, rcond=None)[0]
    log_A, neg_alpha = coeffs[0], coeffs[1]
    alpha = -neg_alpha
    A = np.exp(log_A)

    # R²
    y_pred = A_mat @ coeffs
    ss_res = np.sum((log_y - y_pred)**2)
    ss_tot = np.sum((log_y - log_y.mean())**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0

    return alpha, A, r2, mask.sum()


def main():
    t_start = time.time()

    # σ 격자: 0.1, 0.15, ..., 0.45, 0.55, ..., 0.9 (0.5 제외)
    sigma_left = np.linspace(0.1, 0.47, N_SIGMA // 2)
    sigma_right = np.linspace(0.53, 0.9, N_SIGMA // 2)
    sigmas = np.sort(np.concatenate([sigma_left, sigma_right]))
    delta_sigmas = np.abs(sigmas - 0.5)

    print("=" * 70)
    print("[C-302] Dirichlet E(σ) 프로파일 — Prop 6 보편성 검증")
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"t 구간: [{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점")
    print(f"σ 격자: {len(sigmas)}점, |σ-1/2| ∈ [{delta_sigmas.min():.3f}, {delta_sigmas.max():.3f}]")
    print("=" * 70)

    results = {}

    for char_name, char_info in CHARACTERS.items():
        print(f"\n{'━'*60}")
        print(f"[{char_info['label']}]  (q={char_info['q']})")
        print(f"{'━'*60}")

        # 영점 수
        zeros = find_zeros_dirichlet(char_info, T_RANGE[0], T_RANGE[1])
        N_zeros = len(zeros)
        pi_N = np.pi * N_zeros
        print(f"  영점: {N_zeros}개, πN={pi_N:.4f}")
        for i, z in enumerate(zeros[:5]):
            print(f"    γ_{i+1} = {z:.6f}")
        if N_zeros > 5:
            print(f"    ... (총 {N_zeros}개)")

        # E(σ) 프로파일
        print(f"  E(σ) 계산 ({len(sigmas)}점)...")
        energies = compute_energy_profile(char_info, T_RANGE[0], T_RANGE[1],
                                          N_T_POINTS, sigmas)

        # 대칭성 점검: E(σ) vs E(1-σ)
        n_half = len(sigmas) // 2
        sym_ratios = energies[:n_half] / energies[n_half:][::-1]
        sym_mean = sym_ratios.mean()
        sym_std = sym_ratios.std()
        print(f"  대칭성 E(σ)/E(1-σ): {sym_mean:.4f} ± {sym_std:.4f}")

        # 멱법칙 피팅 (양쪽 합산)
        # 대칭이면 한쪽만 피팅해도 됨. 평균 사용.
        energies_sym = (energies[:n_half] + energies[n_half:][::-1]) / 2
        ds = delta_sigmas[:n_half]

        alpha, A_fit, r2, n_fit = fit_power_law(ds, energies_sym, safe_min=FIT_SAFE)
        print(f"  멱법칙: α={alpha:.4f}, A={A_fit:.4f}, πN={pi_N:.4f}, "
              f"A/πN={A_fit/pi_N:.4f}, R²={r2:.6f} (n_fit={n_fit})")

        # Δσ·E 상수성 (safe 영역)
        mask_safe = ds >= FIT_SAFE
        if mask_safe.sum() > 0:
            ds_E = ds[mask_safe] * energies_sym[mask_safe]
            print(f"  Δσ·E (safe): {ds_E.mean():.4f} ± {ds_E.std():.4f}, "
                  f"CV={ds_E.std()/ds_E.mean()*100:.1f}%")

        results[char_name] = {
            'label': char_info['label'],
            'q': char_info['q'],
            'N_zeros': N_zeros,
            'pi_N': pi_N,
            'alpha': alpha,
            'A_fit': A_fit,
            'r2': r2,
            'n_fit': n_fit,
            'sym_mean': sym_mean,
            'sym_std': sym_std,
            'energies': energies.copy(),
            'energies_sym': energies_sym.copy(),
        }

    elapsed = time.time() - t_start

    # ━━━ 비교 표 ━━━
    print(f"\n{'='*70}")
    print(f"[비교 요약] (총 {elapsed:.1f}s)")
    print(f"{'='*70}")

    # ζ(s) 참조 (C-298)
    print(f"\n  {'지표':<20} {'N':>4} {'α':>8} {'A/πN':>8} {'R²':>10} {'대칭':>8}")
    print(f"  {'-'*60}")
    print(f"  {'ζ(s) (C-298)':<20} {'5':>4} {'1.005':>8} {'0.950':>8} {'0.9999':>10} {'1.000':>8}")

    for name, r in results.items():
        print(f"  {r['label']:<20} {r['N_zeros']:>4} {r['alpha']:>8.4f} "
              f"{r['A_fit']/r['pi_N']:>8.4f} {r['r2']:>10.6f} {r['sym_mean']:>8.4f}")

    # ━━━ 판정 ━━━
    print(f"\n{'='*70}")
    print("[판정]")

    n_pass_alpha = sum(1 for r in results.values() if abs(r['alpha'] - 1) < 0.05)
    n_pass_coeff = sum(1 for r in results.values() if 0.8 <= r['A_fit']/r['pi_N'] <= 1.2)
    n_pass_r2 = sum(1 for r in results.values() if r['r2'] > 0.99)
    n_total = len(results)

    print(f"  |α-1| < 0.05: {n_pass_alpha}/{n_total}")
    print(f"  A/πN ∈ [0.8,1.2]: {n_pass_coeff}/{n_total}")
    print(f"  R² > 0.99: {n_pass_r2}/{n_total}")

    if n_pass_alpha == n_total and n_pass_coeff == n_total:
        verdict = "★★★★★ 강양성 — Prop 6 GL(1) 보편성 확립"
    elif n_pass_alpha >= 2 and n_pass_coeff >= 2:
        verdict = "★★★★ 양성 — 대부분 지표에서 Prop 6 성립"
    elif n_pass_alpha >= 1:
        verdict = "★★★ 중립 — 일부 성립, 일부 불일치"
    else:
        verdict = "★★ 음성 — Prop 6 보편성 미확인"

    print(f"\n  종합: {verdict}")

    # ━━━ 결과 파일 저장 ━━━
    with open(OUT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"[C-302] Dirichlet E(σ) 프로파일 — Prop 6 보편성 검증\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"경과: {elapsed:.1f}s\n")
        f.write(f"t 구간: [{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점\n")
        f.write(f"σ 격자: {len(sigmas)}점\n")
        f.write("=" * 70 + "\n")

        for name, r in results.items():
            f.write(f"\n[{r['label']}] q={r['q']}, N={r['N_zeros']}, πN={r['pi_N']:.4f}\n")
            f.write(f"  α = {r['alpha']:.4f}, A = {r['A_fit']:.4f}, "
                    f"A/πN = {r['A_fit']/r['pi_N']:.4f}, R² = {r['r2']:.6f}\n")
            f.write(f"  대칭: E(σ)/E(1-σ) = {r['sym_mean']:.4f} ± {r['sym_std']:.4f}\n")

            f.write(f"  σ 프로파일:\n")
            f.write(f"    {'σ':>6} {'|σ-½|':>8} {'E(σ)':>12} {'E_sym':>12} {'Δσ·E':>10}\n")
            for j in range(len(sigmas)):
                ds_j = abs(sigmas[j] - 0.5)
                es = r['energies_sym'][j % (len(sigmas)//2)] if j < len(sigmas)//2 else r['energies_sym'][len(sigmas)-1-j]
                f.write(f"    {sigmas[j]:6.3f} {ds_j:8.4f} {r['energies'][j]:12.2f} "
                        f"{es:12.2f} {ds_j*es:10.2f}\n")

        f.write(f"\n{'='*70}\n")
        f.write(f"비교 요약:\n")
        f.write(f"  {'지표':<20} {'N':>4} {'α':>8} {'A/πN':>8} {'R²':>10} {'대칭':>8}\n")
        f.write(f"  {'-'*60}\n")
        f.write(f"  {'ζ(s) (C-298)':<20} {'5':>4} {'1.005':>8} {'0.950':>8} {'0.9999':>10} {'1.000':>8}\n")
        for name, r in results.items():
            f.write(f"  {r['label']:<20} {r['N_zeros']:>4} {r['alpha']:>8.4f} "
                    f"{r['A_fit']/r['pi_N']:>8.4f} {r['r2']:>10.6f} {r['sym_mean']:>8.4f}\n")

        f.write(f"\n종합 판정: {verdict}\n")

    print(f"\n  결과: {OUT_FILE}")
    print(f"  총 경과: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
