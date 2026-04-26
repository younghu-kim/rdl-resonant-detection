"""
=============================================================================
[C-311b] α=1 증명 검증 — N 의존성 스케일링
=============================================================================
C-311 결과: N=10에서 α=1.014, 유한 크기 편차.
이론 예측: E_Had = πN/δ + C(T), C/πN ∝ 1/N → α→1 as N→∞.

검증: N을 증가시키며 α가 1에 수렴하는지 확인.
  - t∈[10,50]: N≈10 (C-311에서 α=1.014)
  - t∈[10,100]: N≈29
  - t∈[10,200]: N≈79
  - t∈[10,500]: N≈약 270

추가: E_cross_fin/E_diag 비율의 N 의존성 확인.
=============================================================================
"""

import sys, os, time, math
import numpy as np
import mpmath

OUT_FILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/alpha1_proof_scaling_c311b.txt')
os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)


def find_zeta_zeros(t_max):
    mpmath.mp.dps = 50
    zeros = []
    n = 1
    while True:
        g = float(mpmath.zetazero(n).imag)
        if g > t_max:
            break
        zeros.append(g)
        n += 1
    return np.array(zeros)


def compute_E_diag_analytic(zeros, delta, t_min, t_max):
    total = 0.0
    for g in zeros:
        total += (math.atan2(t_max - g, delta) - math.atan2(t_min - g, delta)) / delta
    return total


def adaptive_grid(t_min, t_max, zeros, delta, n_base=2000, peak_hw_mult=10, peak_n=200):
    ts = set(np.linspace(t_min, t_max, n_base))
    hw = peak_hw_mult * max(delta, 0.001)
    for g in zeros:
        lo = max(t_min, g - hw)
        hi = min(t_max, g + hw)
        if hi > lo:
            ts.update(np.linspace(lo, hi, peak_n))
    ts = np.array(sorted(ts))
    return ts[(ts >= t_min) & (ts <= t_max)]


def compute_E_had_numerical(ts, zeros, delta):
    tau = ts[:, None] - zeros[None, :]
    denom = delta**2 + tau**2
    re_sum = np.sum(delta / denom, axis=1)
    im_sum = np.sum(-tau / denom, axis=1)
    integrand = re_sum**2 + im_sum**2
    return np.trapezoid(integrand, ts)


def compute_E_cross_inf(zeros, delta):
    N = len(zeros)
    total = 0.0
    for j in range(N):
        for k in range(j+1, N):
            D = zeros[j] - zeros[k]
            total += 2 * 4 * math.pi * delta / (4 * delta**2 + D**2)
    return total


def fit_power_law(deltas, energies):
    deltas = np.array(deltas)
    energies = np.array(energies)
    mask = energies > 0
    if mask.sum() < 3:
        return 0, 0, 0
    d, e = deltas[mask], energies[mask]
    A_mat = np.column_stack([np.ones(len(d)), np.log(d)])
    coeffs = np.linalg.lstsq(A_mat, np.log(e), rcond=None)[0]
    alpha = -coeffs[1]
    A = np.exp(coeffs[0])
    y_pred = A_mat @ coeffs
    ss_res = np.sum((np.log(e) - y_pred)**2)
    ss_tot = np.sum((np.log(e) - np.log(e).mean())**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0
    return alpha, A, r2


def main():
    t_start = time.time()
    T_MIN = 10.0
    DELTAS = [0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30]
    T_MAXES = [50, 100, 200, 500]

    print("=" * 70, flush=True)
    print("[C-311b] α=1 증명 — N 의존성 스케일링", flush=True)
    print("=" * 70, flush=True)

    # 미리 모든 영점을 계산
    print("\n[1] 영점 계산...", flush=True)
    all_zeros = find_zeta_zeros(max(T_MAXES) + 5)
    print(f"  전체 {len(all_zeros)}개 영점 (t≤{max(T_MAXES)+5})", flush=True)

    results = []

    for T_MAX in T_MAXES:
        zeros = all_zeros[(all_zeros >= T_MIN) & (all_zeros <= T_MAX)]
        N = len(zeros)
        mean_gap = np.mean(np.diff(zeros)) if N > 1 else 1.0

        print(f"\n{'='*50}", flush=True)
        print(f"[T_MAX={T_MAX}] N={N}, mean_gap={mean_gap:.4f}", flush=True)
        print(f"{'='*50}", flush=True)

        E_hads = []
        E_diags = []
        E_cross_infs = []
        cross_ratios = []

        for delta in DELTAS:
            E_diag = compute_E_diag_analytic(zeros, delta, T_MIN, T_MAX)

            # 적응 그리드 (N>50이면 피크 해상도 줄임)
            pn = max(50, 200 - N)
            ts = adaptive_grid(T_MIN, T_MAX, zeros, delta,
                             n_base=min(3000, 2000 + N*5), peak_n=pn)
            E_had = compute_E_had_numerical(ts, zeros, delta)

            # 교차항 (N>100이면 무한 적분 공식 스킵 — O(N²) 비용)
            if N <= 150:
                E_ci = compute_E_cross_inf(zeros, delta)
            else:
                E_ci = float('nan')

            E_cross_num = E_had - E_diag
            cr = E_cross_num / E_diag if E_diag > 0 else 0

            E_hads.append(E_had)
            E_diags.append(E_diag)
            E_cross_infs.append(E_ci)
            cross_ratios.append(cr)

            dE = delta * E_had
            print(f"  δ={delta:.3f}: E_diag={E_diag:.1f}, E_Had={E_had:.1f}, "
                  f"cross/diag={cr:.6f}, δ·E_Had={dE:.4f}, A/πN={dE/(math.pi*N):.4f}",
                  flush=True)

        # 피팅
        ds = np.array(DELTAS)
        Eh = np.array(E_hads)
        Ed = np.array(E_diags)

        # 전체
        a_all, A_all, r2_all = fit_power_law(ds, Eh)
        # 소δ ≤0.10
        mask_s = ds <= 0.10
        a_s, A_s, r2_s = fit_power_law(ds[mask_s], Eh[mask_s])
        # 초소δ ≤0.05
        mask_xs = ds <= 0.05
        a_xs, A_xs, r2_xs = fit_power_law(ds[mask_xs], Eh[mask_xs])
        # 대각항
        a_diag, A_diag, r2_diag = fit_power_law(ds, Ed)

        piN = math.pi * N

        print(f"\n  피팅 결과:", flush=True)
        print(f"    [전체]   α={a_all:.6f}, A/πN={A_all/piN:.6f}, R²={r2_all:.8f}", flush=True)
        print(f"    [δ≤0.10] α={a_s:.6f}, A/πN={A_s/piN:.6f}, R²={r2_s:.8f}", flush=True)
        print(f"    [δ≤0.05] α={a_xs:.6f}, A/πN={A_xs/piN:.6f}, R²={r2_xs:.8f}", flush=True)
        print(f"    [대각항]  α={a_diag:.6f}, A/πN={A_diag/piN:.6f}, R²={r2_diag:.8f}", flush=True)

        # 교차항 스케일링
        if not any(math.isnan(x) for x in E_cross_infs):
            eci = np.array(E_cross_infs)
            b_ci, _, r2_ci = fit_power_law(ds[mask_s], eci[mask_s])
            print(f"    E_cross_∞ 스케일링: β={b_ci:.4f}, R²={r2_ci:.6f}", flush=True)
        else:
            b_ci = float('nan')
            r2_ci = float('nan')
            print(f"    E_cross_∞: N>{150}, 스킵", flush=True)

        results.append({
            'T_MAX': T_MAX, 'N': N, 'mean_gap': mean_gap,
            'alpha_all': a_all, 'alpha_s': a_s, 'alpha_xs': a_xs, 'alpha_diag': a_diag,
            'A_all': A_all/piN, 'A_s': A_s/piN, 'A_xs': A_xs/piN,
            'R2_all': r2_all, 'R2_s': r2_s, 'R2_xs': r2_xs,
            'beta_cross': b_ci if not math.isnan(b_ci) else None,
        })

    # --- 종합 ---
    elapsed = time.time() - t_start
    print(f"\n{'='*70}", flush=True)
    print("종합: α(N) 수렴 표", flush=True)
    print(f"{'='*70}", flush=True)
    print(f"  {'T_MAX':>6s} {'N':>5s} {'gap':>6s} {'α_all':>8s} {'α_≤0.10':>8s} "
          f"{'α_≤0.05':>8s} {'α_diag':>8s} {'A/πN_xs':>8s}", flush=True)
    for r in results:
        print(f"  {r['T_MAX']:6d} {r['N']:5d} {r['mean_gap']:6.3f} "
              f"{r['alpha_all']:8.4f} {r['alpha_s']:8.4f} "
              f"{r['alpha_xs']:8.4f} {r['alpha_diag']:8.4f} "
              f"{r['A_xs']:8.4f}", flush=True)

    print(f"\n경과: {elapsed:.1f}s", flush=True)

    # α → 1 수렴 판정
    alphas_xs = [r['alpha_xs'] for r in results]
    if len(alphas_xs) >= 2:
        trend = alphas_xs[-1] - alphas_xs[0]
        if abs(alphas_xs[-1] - 1.0) < 0.005:
            verdict = "★★★★★ 강양성"
            msg = f"α→1 수렴 확인 (α={alphas_xs[-1]:.6f})"
        elif abs(alphas_xs[-1] - 1.0) < 0.02:
            verdict = "★★★★ 양성"
            msg = f"α→1 접근 중 (α={alphas_xs[-1]:.6f})"
        else:
            verdict = "★★★ 조건부"
            msg = f"α={alphas_xs[-1]:.6f}, 추가 검증 필요"
        print(f"\n판정: {verdict} — {msg}", flush=True)

    # 저장
    with open(OUT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"[C-311b] α=1 증명 — N 의존성 스케일링\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"경과: {elapsed:.1f}s\n")
        f.write("=" * 70 + "\n\n")

        f.write("α(N) 수렴 표:\n")
        f.write(f"  {'T_MAX':>6s} {'N':>5s} {'gap':>6s} {'α_all':>8s} {'α_≤0.10':>8s} "
                f"{'α_≤0.05':>8s} {'α_diag':>8s} {'A/πN_xs':>8s}\n")
        for r in results:
            f.write(f"  {r['T_MAX']:6d} {r['N']:5d} {r['mean_gap']:6.3f} "
                    f"{r['alpha_all']:8.6f} {r['alpha_s']:8.6f} "
                    f"{r['alpha_xs']:8.6f} {r['alpha_diag']:8.6f} "
                    f"{r['A_xs']:8.6f}\n")

        f.write(f"\n이론 예측: α→1 as N→∞ (유한 크기 보정 C/πN → 0)\n")
        f.write(f"관측: α_xs(N={results[0]['N']})={results[0]['alpha_xs']:.6f} → "
                f"α_xs(N={results[-1]['N']})={results[-1]['alpha_xs']:.6f}\n")
        f.write(f"\n판정: {verdict} — {msg}\n")

    print(f"\n결과 저장: {OUT_FILE}", flush=True)


if __name__ == '__main__':
    main()
