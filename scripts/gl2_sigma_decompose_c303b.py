"""
=============================================================================
[C-303b] GL(2) E(σ) 분해 — γ-인자 기저선 분리
=============================================================================
목적:
    C-303에서 GL(2) α≈0.66으로 나온 원인 진단.
    E_Λ(σ) = E_L(σ) + E_γ(σ) + E_cross(σ) 로 분해하여
    E_L(σ) = ∫|L'/L(σ+it)|² dt 가 α=1을 보이는지 검증.

    가설: Γ(s)의 digamma ψ(s)가 σ-무관 기저선을 형성하여
    겉보기 α를 낮춤. L'/L만 분리하면 α≈1 회복 예상.
=============================================================================
"""

import sys, os, time, math
import numpy as np
import mpmath

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(500_000_000)
pari.set_real_precision(50)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

SIGMA_CRIT = 1.0
T_RANGE = (5.0, 30.0)
N_T_POINTS = 500          # Δt=0.05 (C-298 동일)
FIT_SAFE = 0.05

OUT_FILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl2_sigma_decompose_c303b.txt'
)


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 곡선 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CURVES = [
    {'name': '11a1', 'coeffs': '[0,-1,1,-10,-20]', 'N_cond': 11, 'eps': 1, 'rank': 0,
     'label': 'GL(2) 11a1 (N=11, ε=+1, rank 0)'},
    {'name': '37a1', 'coeffs': '[0,0,1,-1,0]', 'N_cond': 37, 'eps': -1, 'rank': 1,
     'label': 'GL(2) 37a1 (N=37, ε=-1, rank 1)'},
]


def init_curve(curve):
    name = curve['name']
    coeffs = curve['coeffs']
    N_exp = curve['N_cond']

    pari(f'E_{name} = ellinit({coeffs})')
    actual_N = int(str(pari(f'ellglobalred(E_{name})[1]')))
    if actual_N != N_exp:
        print(f"  ⚠️ [{name}] conductor 불일치: 예상 {N_exp}, 실제 {actual_N}", flush=True)
        return None
    print(f"  [{name}] conductor 검증 ✅ N={actual_N}", flush=True)

    # lfuninit: w=1.0 (넓게)
    T_max = T_RANGE[1] + 5
    pari(f'Li_{name} = lfuninit(E_{name}, [1, 1.0, {int(T_max)}])')

    # 영점
    pari(f'zv_{name} = lfunzeros(Li_{name}, {T_RANGE[1]})')
    n_z = int(str(pari(f'#zv_{name}')))
    zeros = []
    for i in range(1, n_z + 1):
        t = pf(pari(f'zv_{name}[{i}]'))
        if not math.isnan(t) and t >= T_RANGE[0]:
            zeros.append(t)
    zeros = sorted(zeros)
    print(f"  [{name}] 영점: {len(zeros)}개 (t∈[{T_RANGE[0]}, {T_RANGE[1]}])", flush=True)
    return zeros


def compute_decomposed_energy(name, N_cond, sigma, ts):
    """
    E를 3개 성분으로 분해:
      L_conn  = L'/L(s)                  ← 순수 L-함수 접속
      gamma_c = (1/2)log(N)-log(2π)+ψ(s) ← γ-인자 기여
      full    = gamma_c + L_conn         ← 전체 Λ'/Λ

    E_L = ∫|L'/L|² dt, E_γ = ∫|γ-인자|² dt, E_full = ∫|전체|² dt
    """
    const_part = 0.5 * math.log(N_cond) - math.log(2 * math.pi)

    int_L = np.zeros(len(ts))
    int_gamma = np.zeros(len(ts))
    int_full = np.zeros(len(ts))

    mpmath.mp.dps = 30

    for i, t in enumerate(ts):
        try:
            # ψ(s)
            s_mp = mpmath.mpc(sigma, t)
            psi_val = mpmath.digamma(s_mp)
            psi_c = complex(float(mpmath.re(psi_val)), float(mpmath.im(psi_val)))
            gamma_c = const_part + psi_c

            # L(s) and L'(s)
            L_re = pf(pari(f'real(lfun(Li_{name}, {sigma} + {t}*I))'))
            L_im = pf(pari(f'imag(lfun(Li_{name}, {sigma} + {t}*I))'))
            Lp_re = pf(pari(f'real(lfun(Li_{name}, {sigma} + {t}*I, 1))'))
            Lp_im = pf(pari(f'imag(lfun(Li_{name}, {sigma} + {t}*I, 1))'))

            L_val = complex(L_re, L_im)
            Lp_val = complex(Lp_re, Lp_im)

            if abs(L_val) < 1e-30:
                # 영점 근방 → 발산
                int_L[i] = 0
                int_gamma[i] = abs(gamma_c)**2
                int_full[i] = 0
                continue

            LpL = Lp_val / L_val
            full_conn = gamma_c + LpL

            val_L = abs(LpL)**2
            val_g = abs(gamma_c)**2
            val_f = abs(full_conn)**2

            int_L[i] = val_L if np.isfinite(val_L) else 0
            int_gamma[i] = val_g if np.isfinite(val_g) else 0
            int_full[i] = val_f if np.isfinite(val_f) else 0

        except Exception:
            pass

    # 극값 클리핑 (L'/L만)
    finite_L = int_L[int_L < 1e12]
    if len(finite_L) > 10:
        cap = np.percentile(finite_L, 99)
        int_L = np.minimum(int_L, 10 * cap)

    finite_f = int_full[int_full < 1e12]
    if len(finite_f) > 10:
        cap = np.percentile(finite_f, 99)
        int_full = np.minimum(int_full, 10 * cap)

    E_L = np.trapezoid(int_L, ts)
    E_gamma = np.trapezoid(int_gamma, ts)
    E_full = np.trapezoid(int_full, ts)
    E_cross = E_full - E_L - E_gamma  # 교차항 간접 계산

    return E_L, E_gamma, E_full, E_cross


def fit_power_law(delta_sigma, energies, safe_min=0.05):
    mask = delta_sigma >= safe_min
    if mask.sum() < 3:
        return 0, 0, 0, 0
    log_x = np.log(delta_sigma[mask])
    log_y = np.log(np.maximum(energies[mask], 1e-30))
    A_mat = np.column_stack([np.ones(mask.sum()), log_x])
    coeffs = np.linalg.lstsq(A_mat, log_y, rcond=None)[0]
    alpha = -coeffs[1]
    A = np.exp(coeffs[0])
    y_pred = A_mat @ coeffs
    ss_res = np.sum((log_y - y_pred)**2)
    ss_tot = np.sum((log_y - log_y.mean())**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0
    return alpha, A, r2, int(mask.sum())


def main():
    t_start = time.time()

    offsets = np.array([0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5])
    # 양쪽 대칭 측정 후 평균
    ts = np.linspace(T_RANGE[0], T_RANGE[1], N_T_POINTS)

    print("=" * 70, flush=True)
    print("[C-303b] GL(2) E(σ) 분해 — γ-인자 기저선 분리", flush=True)
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print(f"σ_crit = {SIGMA_CRIT}, t ∈ [{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점 (Δt={25/N_T_POINTS:.4f})", flush=True)
    print(f"σ offsets: {offsets}", flush=True)
    print(f"분해: E_full = E_L + E_γ + E_cross", flush=True)
    print("=" * 70, flush=True)

    results = {}

    for curve in CURVES:
        name = curve['name']
        N_cond = curve['N_cond']

        print(f"\n{'━'*60}", flush=True)
        print(f"[{curve['label']}]", flush=True)
        print(f"{'━'*60}", flush=True)

        zeros = init_curve(curve)
        if zeros is None:
            continue
        N_zeros = len(zeros)
        piN = np.pi * N_zeros

        # 양쪽 σ에서 계산 후 평균
        E_L_arr = np.zeros(len(offsets))
        E_gamma_arr = np.zeros(len(offsets))
        E_full_arr = np.zeros(len(offsets))
        E_cross_arr = np.zeros(len(offsets))

        print(f"  E(σ) 분해 계산 ({len(offsets)}×2 σ점 × {N_T_POINTS}t점)...", flush=True)

        for j, off in enumerate(offsets):
            t0 = time.time()
            # σ = σ_c + off
            EL_r, Eg_r, Ef_r, Ec_r = compute_decomposed_energy(name, N_cond, SIGMA_CRIT + off, ts)
            # σ = σ_c - off
            EL_l, Eg_l, Ef_l, Ec_l = compute_decomposed_energy(name, N_cond, SIGMA_CRIT - off, ts)

            E_L_arr[j] = (EL_r + EL_l) / 2
            E_gamma_arr[j] = (Eg_r + Eg_l) / 2
            E_full_arr[j] = (Ef_r + Ef_l) / 2
            E_cross_arr[j] = (Ec_r + Ec_l) / 2

            dt = time.time() - t0
            print(f"    Δσ={off:.2f}: E_L={E_L_arr[j]:.1f}, E_γ={E_gamma_arr[j]:.1f}, "
                  f"E_full={E_full_arr[j]:.1f}, E_cross={E_cross_arr[j]:.1f} ({dt:.1f}s)", flush=True)

        # 멱법칙 피팅
        alpha_full, A_full, r2_full, n_full = fit_power_law(offsets, E_full_arr, FIT_SAFE)
        alpha_L, A_L, r2_L, n_L = fit_power_law(offsets, E_L_arr, FIT_SAFE)
        alpha_g, A_g, r2_g, n_g = fit_power_law(offsets, E_gamma_arr, FIT_SAFE)

        print(f"\n  === 멱법칙 피팅 (Δσ ≥ {FIT_SAFE}) ===", flush=True)
        print(f"  E_full: α={alpha_full:.4f}, A/πN={A_full/piN:.4f}, R²={r2_full:.6f}", flush=True)
        print(f"  E_L:    α={alpha_L:.4f}, A/πN={A_L/piN:.4f}, R²={r2_L:.6f}  ← 핵심", flush=True)
        print(f"  E_γ:    α={alpha_g:.4f}, R²={r2_g:.6f} (σ-무관이면 α≈0)", flush=True)

        # Δσ·E_L 상수성
        mask_safe = offsets >= FIT_SAFE
        if mask_safe.sum() > 0:
            ds_EL = offsets[mask_safe] * E_L_arr[mask_safe]
            print(f"  Δσ·E_L: {ds_EL.mean():.2f} ± {ds_EL.std():.2f}, CV={ds_EL.std()/ds_EL.mean()*100:.1f}%", flush=True)

        results[name] = {
            'label': curve['label'],
            'N_cond': N_cond,
            'rank': curve['rank'],
            'N_zeros': N_zeros,
            'piN': piN,
            'offsets': offsets.copy(),
            'E_L': E_L_arr.copy(),
            'E_gamma': E_gamma_arr.copy(),
            'E_full': E_full_arr.copy(),
            'E_cross': E_cross_arr.copy(),
            'alpha_full': alpha_full,
            'alpha_L': alpha_L,
            'alpha_g': alpha_g,
            'A_L': A_L,
            'r2_full': r2_full,
            'r2_L': r2_L,
        }

    elapsed = time.time() - t_start

    # ━━━ 비교 요약 ━━━
    print(f"\n{'='*70}", flush=True)
    print(f"[비교 요약] 분해 결과 (총 {elapsed:.1f}s)", flush=True)
    print(f"{'='*70}", flush=True)

    print(f"\n  {'곡선':<30} {'α_full':>8} {'α_L':>8} {'α_γ':>8} {'A_L/πN':>8} {'R²_L':>10}", flush=True)
    print(f"  {'-'*74}", flush=True)
    print(f"  {'ζ(s) GL(1) 참조':<30} {'1.005':>8} {'—':>8} {'—':>8} {'0.950':>8} {'0.9999':>10}", flush=True)

    for name, r in results.items():
        a_pn = r['A_L'] / r['piN'] if r['piN'] > 0 else 0
        print(f"  {r['label'][:30]:<30} {r['alpha_full']:>8.4f} {r['alpha_L']:>8.4f} "
              f"{r['alpha_g']:>8.4f} {a_pn:>8.4f} {r['r2_L']:>10.6f}", flush=True)

    # 판정
    print(f"\n{'='*70}", flush=True)
    print("[판정]", flush=True)

    for name, r in results.items():
        a_pn = r['A_L'] / r['piN'] if r['piN'] > 0 else 0
        alpha_dev = abs(r['alpha_L'] - 1.0)
        if alpha_dev < 0.05 and r['r2_L'] > 0.99:
            verdict = "★★★★★ 강양성"
        elif alpha_dev < 0.15 and r['r2_L'] > 0.99:
            verdict = "★★★★ 양성"
        elif alpha_dev < 0.15:
            verdict = "★★★ 중립"
        else:
            verdict = "★★ 음성"
        print(f"  {name}: E_L α={r['alpha_L']:.4f}, A_L/πN={a_pn:.4f}, R²={r['r2_L']:.6f} → {verdict}", flush=True)

    alpha_Ls = [r['alpha_L'] for r in results.values()]
    if len(alpha_Ls) == 2:
        diff = abs(alpha_Ls[0] - alpha_Ls[1])
        print(f"\n  α_L 차이(11a1-37a1): {diff:.4f}", flush=True)

    # γ-인자 오염 정량화
    print(f"\n  [γ-인자 오염 분석]", flush=True)
    for name, r in results.items():
        at_close = r['E_gamma'][0] / r['E_full'][0] * 100  # Δσ 최소에서
        at_far = r['E_gamma'][-1] / r['E_full'][-1] * 100  # Δσ 최대에서
        print(f"  {name}: E_γ/E_full = {at_close:.1f}% (Δσ={r['offsets'][0]}) → {at_far:.1f}% (Δσ={r['offsets'][-1]})", flush=True)

    # ━━━ 결과 파일 ━━━
    with open(OUT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"[C-303b] GL(2) E(σ) 분해 — γ-인자 기저선 분리\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"경과: {elapsed:.1f}s\n")
        f.write(f"σ_crit = {SIGMA_CRIT}, t∈[{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점\n")
        f.write(f"분해: E_full = E_L(|L'/L|²) + E_γ(|ψ+const|²) + E_cross\n")
        f.write("=" * 70 + "\n")

        for name, r in results.items():
            a_pn = r['A_L'] / r['piN'] if r['piN'] > 0 else 0
            f.write(f"\n{'━'*60}\n")
            f.write(f"[{r['label']}]\n")
            f.write(f"  N_zeros={r['N_zeros']}, πN={r['piN']:.4f}\n\n")

            f.write(f"  {'Δσ':>6} {'E_full':>12} {'E_L':>12} {'E_γ':>12} {'E_cross':>12} {'E_γ/E_f%':>10} {'Δσ·E_L':>10}\n")
            for j in range(len(r['offsets'])):
                ds = r['offsets'][j]
                ef = r['E_full'][j]
                el = r['E_L'][j]
                eg = r['E_gamma'][j]
                ec = r['E_cross'][j]
                pct = eg/ef*100 if ef > 0 else 0
                f.write(f"  {ds:6.3f} {ef:12.2f} {el:12.2f} {eg:12.2f} {ec:12.2f} {pct:10.1f} {ds*el:10.2f}\n")

            f.write(f"\n  멱법칙 (E_full): α={r['alpha_full']:.6f}, R²={r['r2_full']:.6f}\n")
            f.write(f"  멱법칙 (E_L):    α={r['alpha_L']:.6f}, A_L={r['A_L']:.4f}, A_L/πN={a_pn:.6f}, R²={r['r2_L']:.6f}\n")
            f.write(f"  멱법칙 (E_γ):    α={r['alpha_g']:.6f}, R²={r['r2_g']:.6f}\n")

        f.write(f"\n{'='*70}\n")
        f.write(f"비교:\n")
        f.write(f"  GL(1) ζ(s):  α_full=1.005, A/πN=0.950, R²=0.9999\n")
        for name, r in results.items():
            a_pn = r['A_L'] / r['piN'] if r['piN'] > 0 else 0
            f.write(f"  {name}: α_full={r['alpha_full']:.4f}, α_L={r['alpha_L']:.4f}, A_L/πN={a_pn:.4f}, R²_L={r['r2_L']:.6f}\n")

        f.write(f"\n가설: E_Λ에서 α<1은 Γ(s) digamma 기저선 오염.\n")
        f.write(f"E_L(|L'/L|²만)에서 α≈1 회복 여부가 Prop 6 보편성 판정 핵심.\n")

    print(f"\n  결과: {OUT_FILE}", flush=True)
    print(f"  총 경과: {elapsed:.1f}s", flush=True)


if __name__ == '__main__':
    main()
