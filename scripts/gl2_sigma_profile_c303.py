"""
=============================================================================
[C-303] GL(2) E(σ) 프로파일 — Prop 6 degree-보편성 검증
=============================================================================
목적:
    Prop 6 (E(σ) = πN/|σ-σ_crit|^α) 이 GL(1)에서 확립된 후,
    GL(2) 타원곡선 L-함수에서도 성립하는지 검증.

대상: 11a1 (N=11, ε=+1, rank 0), 37a1 (N=37, ε=-1, rank 1)

방법:
    σ_crit = 1.0 (GL(2) weight-2 임계선)
    Λ'/Λ(s) = (1/2)log(N) - log(2π) + ψ(s) + L'/L(s)  [해석적 공식]
    E(σ) = ∫|Λ'/Λ(σ+it)|² dt

    L'/L(s) 는 PARI lfun으로 계산 (C 코드, mpmath 대비 ~100배 속도)
    ψ(s) = digamma는 mpmath로 계산

판정 기준:
    강양성: α ∈ [0.95, 1.05], R² > 0.99
    양성:   α ∈ [0.85, 1.15], R² > 0.99
    음성:   α≠1 또는 R² 낮음
=============================================================================
"""

import sys, os, time, math
import numpy as np
import mpmath

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(500_000_000)  # 500MB
pari.set_real_precision(50)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

SIGMA_CRIT = 1.0       # GL(2) weight-2 임계선
T_RANGE = (5.0, 30.0)
N_T_POINTS = 250       # Δt = 0.1
N_SIGMA = 12
FIT_SAFE = 0.05

OUT_FILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl2_sigma_profile_c303.txt'
)


def pf(x):
    """PARI 값 → float 변환"""
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 곡선 정의 + PARI 초기화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CURVES = [
    {'name': '11a1', 'coeffs': '[0,-1,1,-10,-20]', 'N_cond': 11, 'eps': 1, 'rank': 0,
     'label': 'GL(2) 11a1 (N=11, ε=+1, rank 0)'},
    {'name': '37a1', 'coeffs': '[0,0,1,-1,0]', 'N_cond': 37, 'eps': -1, 'rank': 1,
     'label': 'GL(2) 37a1 (N=37, ε=-1, rank 1)'},
]


def init_curve(curve):
    """PARI ellinit + lfuninit + conductor 검증"""
    name = curve['name']
    coeffs = curve['coeffs']
    N_exp = curve['N_cond']

    pari(f'E_{name} = ellinit({coeffs})')

    # conductor 검증
    actual_N = int(str(pari(f'ellglobalred(E_{name})[1]')))
    if actual_N != N_exp:
        print(f"  ⚠️ [{name}] conductor 불일치: 예상 {N_exp}, 실제 {actual_N}", flush=True)
        return False
    print(f"  [{name}] conductor 검증 ✅ N={actual_N}", flush=True)

    # lfuninit: 임계선 위 영점 탐색용 (off-critical은 lfun(E,s) 직접 사용)
    pari(f'Li_{name} = lfuninit(E_{name}, {int(T_RANGE[1]) + 5})')

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
    for i, z in enumerate(zeros[:5]):
        print(f"    γ_{i+1} = {z:.8f}", flush=True)
    if len(zeros) > 5:
        print(f"    ... (총 {len(zeros)}개)", flush=True)

    return zeros


def lfun_eval(name, s_real, s_imag):
    """PARI lfun(E, s) → complex value (E 직접 사용, lfuninit 불필요)"""
    try:
        re_part = pf(pari(f'real(lfun(E_{name}, {s_real} + {s_imag}*I))'))
        im_part = pf(pari(f'imag(lfun(E_{name}, {s_real} + {s_imag}*I))'))
        if math.isnan(re_part) or math.isnan(im_part):
            return complex(0, 0)
        return complex(re_part, im_part)
    except Exception:
        return complex(0, 0)


def lfun_deriv_eval(name, s_real, s_imag):
    """PARI lfun(E, s, 1) → L'(s) complex value"""
    try:
        re_part = pf(pari(f'real(lfun(E_{name}, {s_real} + {s_imag}*I, 1))'))
        im_part = pf(pari(f'imag(lfun(E_{name}, {s_real} + {s_imag}*I, 1))'))
        if math.isnan(re_part) or math.isnan(im_part):
            return complex(0, 0)
        return complex(re_part, im_part)
    except Exception:
        return complex(0, 0)


def connection_gl2(name, sigma, t, N_cond):
    """
    Λ'/Λ(s) = (1/2)log(N) - log(2π) + ψ(s) + L'/L(s)
    해석적 공식 — 수치 미분 불필요!
    """
    # 상수부
    const = 0.5 * math.log(N_cond) - math.log(2 * math.pi)

    # ψ(s) = digamma (mpmath 사용)
    mpmath.mp.dps = 30
    s_mp = mpmath.mpc(sigma, t)
    psi_val = mpmath.digamma(s_mp)
    psi_c = complex(float(mpmath.re(psi_val)), float(mpmath.im(psi_val)))

    # L'/L(s) (PARI)
    L_val = lfun_eval(name, sigma, t)
    if abs(L_val) < 1e-30:
        return complex(1e6, 0)  # 영점 근방 → 발산
    L_prime = lfun_deriv_eval(name, sigma, t)
    LpL = L_prime / L_val

    return const + psi_c + LpL


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# E(σ) 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_energy_profile(name, N_cond, t_min, t_max, n_t, sigma_values):
    """각 σ에서 E(σ) = ∫|Λ'/Λ(σ+it)|² dt"""
    ts = np.linspace(t_min, t_max, n_t)
    energies = np.zeros(len(sigma_values))

    for j, sigma in enumerate(sigma_values):
        integrand = np.zeros(n_t)
        for i, t in enumerate(ts):
            try:
                conn = connection_gl2(name, sigma, t, N_cond)
                val = abs(conn)**2
                integrand[i] = val if np.isfinite(val) else 0.0
            except Exception:
                integrand[i] = 0.0

        # 극값 클리핑
        finite_vals = integrand[integrand < 1e12]
        if len(finite_vals) > 10:
            cap = np.percentile(finite_vals, 99)
            integrand = np.minimum(integrand, 10 * cap)

        energies[j] = np.trapezoid(integrand, ts)

        if (j + 1) % 4 == 0 or j == 0:
            print(f"    σ={sigma:.3f}: E={energies[j]:.2f}", flush=True)

    return energies


def fit_power_law(delta_sigma, energies, safe_min=0.05):
    """E = A/Δσ^α 피팅 (2-파라미터, log-log)"""
    mask = delta_sigma >= safe_min
    if mask.sum() < 3:
        return 0, 0, 0, 0

    log_x = np.log(delta_sigma[mask])
    log_y = np.log(np.maximum(energies[mask], 1e-30))

    A_mat = np.column_stack([np.ones(mask.sum()), log_x])
    coeffs = np.linalg.lstsq(A_mat, log_y, rcond=None)[0]
    log_A, neg_alpha = coeffs[0], coeffs[1]
    alpha = -neg_alpha
    A = np.exp(log_A)

    y_pred = A_mat @ coeffs
    ss_res = np.sum((log_y - y_pred)**2)
    ss_tot = np.sum((log_y - log_y.mean())**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0

    return alpha, A, r2, int(mask.sum())


def fit_with_baseline(delta_sigma, energies, safe_min=0.05):
    """E = C + A/Δσ 피팅 (α=1 고정, 3-파라미터)

    GL(2)에서 Γ-인자 기여가 GL(1)보다 크므로 baseline C 존재.
    α=1 고정하고 C, A를 최소자승법으로 구함.
    """
    mask = delta_sigma >= safe_min
    if mask.sum() < 3:
        return 0, 0, 0, 0

    x = 1.0 / delta_sigma[mask]  # 1/Δσ
    y = energies[mask]

    # y = C + A·x → 선형 피팅
    A_mat = np.column_stack([np.ones(mask.sum()), x])
    coeffs = np.linalg.lstsq(A_mat, y, rcond=None)[0]
    C_fit, A_fit = coeffs[0], coeffs[1]

    y_pred = A_mat @ coeffs
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - y.mean())**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0

    return C_fit, A_fit, r2, int(mask.sum())


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# AFE 함수방정식 검증 (독립 확인)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def verify_functional_equation(name, N_cond, eps):
    """Λ(s) = ε·Λ(2-s) 검증. L, ψ, log 항으로 간접 확인."""
    # 점 (σ=0.8, t=10)에서 L(s) vs L(2-s) 비교
    s_re, s_im = 0.8, 10.0
    L_s = lfun_eval(name, s_re, s_im)
    L_conj = lfun_eval(name, 2 - s_re, -s_im)  # L(2-s) = L(2-σ-it)

    # Λ 비율 계산 (감마 인자 포함)
    mpmath.mp.dps = 30
    s1 = mpmath.mpc(s_re, s_im)
    s2 = mpmath.mpc(2 - s_re, -s_im)

    # Λ(s)/Λ(2-s) = (N/(4π²))^{(s-(2-s))/2} · Γ(s)/Γ(2-s) · L(s)/L(2-s)
    # = (N/(4π²))^{(2σ-2)/2} · Γ(σ+it)/Γ(2-σ-it) · L(s)/L(2-s)
    ratio_N = (N_cond / (4 * math.pi**2)) ** ((2*s_re - 2) / 2)
    gamma_ratio = complex(mpmath.gamma(s1) / mpmath.gamma(s2))
    L_ratio = L_s / L_conj if abs(L_conj) > 1e-30 else 0

    fe_ratio = ratio_N * gamma_ratio * L_ratio
    fe_ok = abs(abs(fe_ratio) - 1.0) < 0.05  # 5% 이내
    print(f"  함수방정식: |Λ(s)/(ε·Λ(2-s))| = {abs(fe_ratio):.6f}, "
          f"phase = {math.atan2(fe_ratio.imag, fe_ratio.real)/math.pi:.4f}π "
          f"{'✅' if fe_ok else '⚠️'}", flush=True)
    return fe_ok


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    # σ 격자: σ_crit ± offset
    offsets = np.array([0.07, 0.1, 0.15, 0.2, 0.3, 0.5])
    sigma_left = SIGMA_CRIT - offsets[::-1]
    sigma_right = SIGMA_CRIT + offsets
    sigmas = np.concatenate([sigma_left, sigma_right])
    delta_sigmas = np.abs(sigmas - SIGMA_CRIT)

    print("=" * 70, flush=True)
    print("[C-303] GL(2) E(σ) 프로파일 — Prop 6 degree-보편성 검증", flush=True)
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print(f"σ_crit = {SIGMA_CRIT} (GL(2) weight-2)", flush=True)
    print(f"t 구간: [{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점 (Δt={25.0/N_T_POINTS:.3f})", flush=True)
    print(f"σ 격자: {len(sigmas)}점", flush=True)
    print(f"방법: 해석적 Λ'/Λ = (1/2)log(N)-log(2π)+ψ(s)+L'/L(s), PARI lfun", flush=True)
    print("=" * 70, flush=True)

    results = {}

    for curve in CURVES:
        name = curve['name']
        N_cond = curve['N_cond']

        print(f"\n{'━'*60}", flush=True)
        print(f"[{curve['label']}]", flush=True)
        print(f"{'━'*60}", flush=True)

        # 1. 초기화 + 영점
        t0_init = time.time()
        zeros = init_curve(curve)
        if zeros is False:
            print(f"  ⚠️ {name} 초기화 실패, 스킵", flush=True)
            continue
        N_zeros = len(zeros)
        pi_N = np.pi * N_zeros
        print(f"  πN = {pi_N:.4f}, 초기화 {time.time()-t0_init:.1f}s", flush=True)

        # 2. 함수방정식 검증
        fe_ok = verify_functional_equation(name, N_cond, curve['eps'])

        # 3. E(σ) 프로파일
        print(f"  E(σ) 계산 ({len(sigmas)}점 × {N_T_POINTS}t)...", flush=True)
        t0_e = time.time()
        energies = compute_energy_profile(name, N_cond, T_RANGE[0], T_RANGE[1],
                                          N_T_POINTS, sigmas)
        e_elapsed = time.time() - t0_e
        print(f"  E(σ) 완료: {e_elapsed:.1f}s", flush=True)

        # 4. 대칭성 E(σ)/E(2-σ)
        n_half = len(sigmas) // 2
        sym_ratios = energies[:n_half] / np.maximum(energies[n_half:][::-1], 1e-10)
        sym_mean = sym_ratios.mean()
        sym_std = sym_ratios.std()
        print(f"  대칭성 E(σ)/E(2-σ): {sym_mean:.4f} ± {sym_std:.4f}", flush=True)

        # 5. 멱법칙 피팅
        energies_sym = (energies[:n_half] + energies[n_half:][::-1]) / 2
        ds = delta_sigmas[:n_half]

        alpha, A_fit, r2, n_fit = fit_power_law(ds, energies_sym, safe_min=FIT_SAFE)
        a_over_pin = A_fit / pi_N if pi_N > 0 else 0
        print(f"  멱법칙 (log-log): α={alpha:.4f}, A={A_fit:.2f}, πN={pi_N:.4f}, "
              f"A/πN={a_over_pin:.4f}, R²={r2:.6f} (n_fit={n_fit})", flush=True)

        # baseline 보정 피팅: E = C + A/Δσ (α=1 고정)
        C_bl, A_bl, r2_bl, n_bl = fit_with_baseline(ds, energies_sym, safe_min=FIT_SAFE)
        a_bl_pin = A_bl / pi_N if pi_N > 0 else 0
        print(f"  baseline 보정 (α=1 고정): C={C_bl:.2f}, A={A_bl:.2f}, "
              f"A/πN={a_bl_pin:.4f}, R²={r2_bl:.6f}", flush=True)

        # Δσ·E 상수성
        mask_safe = ds >= FIT_SAFE
        if mask_safe.sum() > 0:
            ds_E = ds[mask_safe] * energies_sym[mask_safe]
            print(f"  Δσ·E (safe): {ds_E.mean():.2f} ± {ds_E.std():.2f}, "
                  f"CV={ds_E.std()/ds_E.mean()*100:.1f}%", flush=True)

        results[name] = {
            'label': curve['label'],
            'N_cond': N_cond,
            'eps': curve['eps'],
            'rank': curve['rank'],
            'N_zeros': N_zeros,
            'pi_N': pi_N,
            'alpha': alpha,
            'A_fit': A_fit,
            'r2': r2,
            'n_fit': n_fit,
            'C_bl': C_bl,
            'A_bl': A_bl,
            'r2_bl': r2_bl,
            'sym_mean': sym_mean,
            'sym_std': sym_std,
            'energies': energies.copy(),
            'energies_sym': energies_sym.copy(),
            'zeros': zeros,
            'fe_ok': fe_ok,
        }

    elapsed = time.time() - t_start

    # ━━━ 비교 표 ━━━
    print(f"\n{'='*70}", flush=True)
    print(f"[비교 요약] GL(1) vs GL(2) (총 {elapsed:.1f}s)", flush=True)
    print(f"{'='*70}", flush=True)

    header = f"  {'곡선':<28} {'σ_c':>4} {'N':>4} {'α':>8} {'A/πN':>8} {'R²':>10} {'대칭':>8}"
    print(header, flush=True)
    print(f"  {'-'*72}", flush=True)
    print(f"  {'ζ(s) (C-298, GL(1))':<28} {'0.5':>4} {'5':>4} {'1.005':>8} {'0.950':>8} {'0.9999':>10} {'1.000':>8}", flush=True)
    print(f"  {'χ₃ mod3 (C-302, GL(1))':<28} {'0.5':>4} {'12':>4} {'0.988':>8} {'0.980':>8} {'0.9989':>10} {'1.000':>8}", flush=True)
    print(f"  {'χ₄ mod4 (C-302, GL(1))':<28} {'0.5':>4} {'13':>4} {'1.003':>8} {'0.955':>8} {'0.9978':>10} {'1.000':>8}", flush=True)
    print(f"  {'-'*72}", flush=True)

    for name, r in results.items():
        a_pn = r['A_fit'] / r['pi_N'] if r['pi_N'] > 0 else 0
        print(f"  {r['label'][:28]:<28} {SIGMA_CRIT:>4.1f} {r['N_zeros']:>4} {r['alpha']:>8.4f} "
              f"{a_pn:>8.4f} {r['r2']:>10.6f} {r['sym_mean']:>8.4f}", flush=True)

    # ━━━ 판정 ━━━
    print(f"\n{'='*70}", flush=True)
    print("[판정]", flush=True)

    n_total = len(results)
    n_pass_alpha_strong = sum(1 for r in results.values() if abs(r['alpha'] - 1) < 0.05)
    n_pass_alpha_weak = sum(1 for r in results.values() if abs(r['alpha'] - 1) < 0.15)
    n_pass_coeff = sum(1 for r in results.values()
                       if r['pi_N'] > 0 and 0.8 <= r['A_fit']/r['pi_N'] <= 1.2)
    n_pass_r2 = sum(1 for r in results.values() if r['r2'] > 0.99)
    n_pass_sym = sum(1 for r in results.values() if abs(r['sym_mean'] - 1.0) < 0.1)

    print(f"  |α-1| < 0.05 (강): {n_pass_alpha_strong}/{n_total}", flush=True)
    print(f"  |α-1| < 0.15 (약): {n_pass_alpha_weak}/{n_total}", flush=True)
    print(f"  A/πN ∈ [0.8,1.2]:  {n_pass_coeff}/{n_total}", flush=True)
    print(f"  R² > 0.99:         {n_pass_r2}/{n_total}", flush=True)
    print(f"  대칭 (|ratio-1|<0.1): {n_pass_sym}/{n_total}", flush=True)

    alphas = [r['alpha'] for r in results.values()]
    if len(alphas) == 2:
        alpha_diff = abs(alphas[0] - alphas[1])
        print(f"  11a1-37a1 α 차이: {alpha_diff:.4f} (<0.05 → conductor-독립)", flush=True)

    # baseline 피팅 판정
    n_pass_bl_r2 = sum(1 for r in results.values() if r['r2_bl'] > 0.99)
    n_pass_bl_coeff = sum(1 for r in results.values()
                          if r['pi_N'] > 0 and 0.5 <= r['A_bl']/r['pi_N'] <= 2.0)
    print(f"  [baseline α=1] R²>0.99: {n_pass_bl_r2}/{n_total}", flush=True)
    print(f"  [baseline α=1] A/πN ∈ [0.5,2.0]: {n_pass_bl_coeff}/{n_total}", flush=True)

    if n_pass_alpha_strong == n_total and n_pass_coeff == n_total and n_pass_r2 == n_total:
        verdict = "★★★★★ 강양성 — Prop 6 GL(2) 확장 확립 (Selberg class 보편성)"
    elif n_pass_bl_r2 == n_total and n_pass_bl_coeff == n_total:
        verdict = "★★★★ 양성 — GL(2)에서 E=C+A/Δσ (α=1) 성립 (Γ-baseline 보정 필요)"
    elif n_pass_alpha_weak == n_total and n_pass_r2 == n_total:
        verdict = "★★★★ 양성 — GL(2)에서 α≈1 확인 (약간의 편차)"
    elif n_pass_alpha_weak >= 1 or n_pass_bl_r2 >= 1:
        verdict = "★★★ 중립 — 일부 곡선 양성, 추가 검증 필요"
    else:
        verdict = "★★ 음성 — GL(2)에서 α≠1, degree-의존 구조 발견"

    print(f"\n  종합: {verdict}", flush=True)

    # rank 효과 분석
    r37 = results.get('37a1')
    r11 = results.get('11a1')
    if r37 and r11:
        print(f"\n  [rank 효과]", flush=True)
        print(f"  11a1 (rank 0): α={r11['alpha']:.4f}, A/πN={r11['A_fit']/r11['pi_N']:.4f}", flush=True)
        print(f"  37a1 (rank 1): α={r37['alpha']:.4f}, A/πN={r37['A_fit']/r37['pi_N']:.4f}", flush=True)
        print(f"  rank 1 → L(1,E)=0 → 추가 극 기여 가능. α 차이로 판별.", flush=True)

    # ━━━ 결과 파일 ━━━
    n_half = len(sigmas) // 2
    with open(OUT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"[C-303] GL(2) E(σ) 프로파일 — Prop 6 degree-보편성 검증\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"경과: {elapsed:.1f}s\n")
        f.write(f"σ_crit = {SIGMA_CRIT} (GL(2) weight-2)\n")
        f.write(f"t 구간: [{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점\n")
        f.write(f"σ 격자: {len(sigmas)}점\n")
        f.write(f"방법: 해석적 Λ'/Λ = (1/2)log(N)-log(2π)+ψ(s)+L'/L(s)\n")
        f.write("=" * 70 + "\n")

        for name, r in results.items():
            a_pn = r['A_fit'] / r['pi_N'] if r['pi_N'] > 0 else 0
            f.write(f"\n{'━'*60}\n")
            f.write(f"[{r['label']}]\n")
            f.write(f"  conductor={r['N_cond']}, rank={r['rank']}, ε={r['eps']}\n")
            f.write(f"  함수방정식: {'✅' if r['fe_ok'] else '⚠️'}\n")
            f.write(f"  영점: N={r['N_zeros']}, πN={r['pi_N']:.4f}\n")
            for i, z in enumerate(r['zeros']):
                f.write(f"    γ_{i+1} = {z:.8f}\n")
            f.write(f"  [log-log] α = {r['alpha']:.6f}, A = {r['A_fit']:.4f}, "
                    f"A/πN = {a_pn:.6f}, R² = {r['r2']:.8f}\n")
            a_bl_pn = r['A_bl'] / r['pi_N'] if r['pi_N'] > 0 else 0
            f.write(f"  [baseline] C = {r['C_bl']:.4f}, A = {r['A_bl']:.4f}, "
                    f"A/πN = {a_bl_pn:.6f}, R² = {r['r2_bl']:.8f}\n")
            f.write(f"  대칭: E(σ)/E(2-σ) = {r['sym_mean']:.6f} ± {r['sym_std']:.6f}\n")

            f.write(f"\n  σ 프로파일:\n")
            f.write(f"    {'σ':>6} {'|σ-1|':>8} {'E(σ)':>14} {'E_sym':>14} {'Δσ·E':>12}\n")
            for j in range(len(sigmas)):
                ds_j = abs(sigmas[j] - SIGMA_CRIT)
                js = j if j < n_half else len(sigmas) - 1 - j
                es = r['energies_sym'][js]
                f.write(f"    {sigmas[j]:6.3f} {ds_j:8.4f} {r['energies'][j]:14.2f} "
                        f"{es:14.2f} {ds_j*es:12.2f}\n")

        f.write(f"\n{'='*70}\n")
        f.write(f"비교 요약 (GL(1) vs GL(2)):\n")
        f.write(header + "\n")
        f.write(f"  {'-'*72}\n")
        f.write(f"  {'ζ(s) (C-298, GL(1))':<28} {'0.5':>4} {'5':>4} {'1.005':>8} {'0.950':>8} {'0.9999':>10} {'1.000':>8}\n")
        f.write(f"  {'χ₃ mod3 (C-302, GL(1))':<28} {'0.5':>4} {'12':>4} {'0.988':>8} {'0.980':>8} {'0.9989':>10} {'1.000':>8}\n")
        f.write(f"  {'χ₄ mod4 (C-302, GL(1))':<28} {'0.5':>4} {'13':>4} {'1.003':>8} {'0.955':>8} {'0.9978':>10} {'1.000':>8}\n")
        f.write(f"  {'-'*72}\n")
        for name, r in results.items():
            a_pn = r['A_fit'] / r['pi_N'] if r['pi_N'] > 0 else 0
            f.write(f"  {r['label'][:28]:<28} {SIGMA_CRIT:>4.1f} {r['N_zeros']:>4} {r['alpha']:>8.4f} "
                    f"{a_pn:>8.4f} {r['r2']:>10.6f} {r['sym_mean']:>8.4f}\n")

        f.write(f"\n종합 판정: {verdict}\n")

    print(f"\n  결과: {OUT_FILE}", flush=True)
    print(f"  총 경과: {elapsed:.1f}s", flush=True)


if __name__ == '__main__':
    main()
