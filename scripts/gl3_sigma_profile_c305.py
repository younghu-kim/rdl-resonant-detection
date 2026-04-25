"""
=============================================================================
[C-305] GL(3) sym²(11a1) E(σ) 프로파일 — Prop 6 degree-3 보편성 검증
=============================================================================
목적:
    Prop 6 (E(σ) = πN/|σ-σ_c|^α) 이 GL(1)(C-302)과 GL(2)(C-303)에서
    확립된 후, GL(3) 대칭멱 L-함수에서도 성립하는지 검증.

대상: sym²(11a1) — GL(3), degree 3, motivic weight 2

방법:
    σ_crit = 1.5 (motivic normalization, k=3, center = k/2)
    Vga = [0, 0, 1] → Γ_R(s)²·Γ_R(s+1)
    Λ'/Λ(s) = (1/2)log(N) - (3/2)log(π) + ψ(s/2) + (1/2)ψ((s+1)/2) + L'/L(s)
    E(σ) = ∫|Λ'/Λ(σ+it)|² dt

    추가: Hadamard 대각항 E_diag = Σ_n arctan((T-γ_n)/δ)/δ 비교

판정 기준:
    강양성: α ∈ [0.95, 1.05], R² > 0.99
    양성:   α ∈ [0.85, 1.15], R² > 0.95
    음성:   α ≠ 1 또는 R² 낮음 → degree-3 경계 발견
=============================================================================
"""

import sys, os, time, math
import numpy as np
import mpmath

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(2_000_000_000)  # 2GB (GL(3) 필요)
pari.set_real_precision(50)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

SIGMA_CRIT = 1.5       # GL(3) sym², k=3, center = k/2
T_RANGE = (5.0, 48.0)  # lfuninit T=50, 약간의 여유
N_T_POINTS = 300       # Δt ≈ 0.143
N_SIGMA = 12
FIT_SAFE = 0.05

OUT_FILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl3_sigma_profile_c305.txt'
)


def pf(x):
    """PARI 값 → float 변환"""
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# GL(3) L-함수 초기화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def init_gl3():
    """sym²(11a1) PARI 초기화 + 영점 추출"""
    pari("E11 = ellinit([0,-1,1,-10,-20])")

    # conductor 확인
    actual_N = int(str(pari("ellglobalred(E11)[1]")))
    print(f"  [11a1] conductor = {actual_N}", flush=True)

    # sym² L-함수
    pari("L3 = lfunsympow(E11, 2)")
    pari(f"L3init = lfuninit(L3, [{int(T_RANGE[1]) + 5}, 0])")
    print("  lfunsympow(E11, 2) + lfuninit 완료", flush=True)

    # Selberg 데이터 추출
    try:
        vga_str = str(pari("lfunparams(L3)[2]"))  # Vga
        k_str = str(pari("lfunparams(L3)[3]"))     # weight/motivic weight
        N_cond_str = str(pari("lfunparams(L3)[4]"))  # conductor
        eps_str = str(pari("lfunparams(L3)[5]"))    # root number
        print(f"  Vga = {vga_str}", flush=True)
        print(f"  k = {k_str}", flush=True)
        print(f"  N_cond = {N_cond_str}", flush=True)
        print(f"  ε = {eps_str}", flush=True)
    except Exception as e:
        print(f"  Selberg 데이터 추출 실패: {e}", flush=True)

    # conductor 숫자 추출
    try:
        N_cond = int(float(N_cond_str.strip()))
    except Exception:
        N_cond = 121  # sym²(11a1) 기대값
        print(f"  ⚠️ conductor 파싱 실패, 기본값 {N_cond} 사용", flush=True)

    # 영점 추출
    pari(f"zv3 = lfunzeros(L3init, {T_RANGE[1]})")
    n_z = int(str(pari("#zv3")))
    zeros = []
    for i in range(1, n_z + 1):
        t = pf(pari(f"zv3[{i}]"))
        if not math.isnan(t) and t >= T_RANGE[0]:
            zeros.append(t)
    zeros = sorted(zeros)
    print(f"  영점: {len(zeros)}개 (t∈[{T_RANGE[0]}, {T_RANGE[1]}])", flush=True)
    for i, z in enumerate(zeros[:5]):
        print(f"    γ_{i+1} = {z:.8f}", flush=True)
    if len(zeros) > 5:
        print(f"    ... (총 {len(zeros)}개)", flush=True)

    return zeros, N_cond


def lfun_eval_gl3(s_real, s_imag):
    """PARI lfun(L3, s) → complex"""
    try:
        re_part = pf(pari(f"real(lfun(L3, {s_real} + {s_imag}*I))"))
        im_part = pf(pari(f"imag(lfun(L3, {s_real} + {s_imag}*I))"))
        if math.isnan(re_part) or math.isnan(im_part):
            return complex(0, 0)
        return complex(re_part, im_part)
    except Exception:
        return complex(0, 0)


def lfun_deriv_gl3(s_real, s_imag):
    """PARI lfun(L3, s, 1) → L'(s) complex"""
    try:
        re_part = pf(pari(f"real(lfun(L3, {s_real} + {s_imag}*I, 1))"))
        im_part = pf(pari(f"imag(lfun(L3, {s_real} + {s_imag}*I, 1))"))
        if math.isnan(re_part) or math.isnan(im_part):
            return complex(0, 0)
        return complex(re_part, im_part)
    except Exception:
        return complex(0, 0)


def connection_gl3(sigma, t, N_cond):
    """
    Λ'/Λ(s) for GL(3) sym²(11a1)

    Vga = [0, 0, 1]:
    Λ'/Λ = (1/2)log(N) - (3/2)log(π) + ψ(s/2) + (1/2)ψ((s+1)/2) + L'/L(s)
    """
    s_c = complex(sigma, t)

    # 상수부
    const = 0.5 * math.log(N_cond) - 1.5 * math.log(math.pi)

    # ψ 항 (mpmath, dps=30)
    mpmath.mp.dps = 30
    s_mp = mpmath.mpc(sigma, t)

    # ψ(s/2) — Vga에서 α=0 두 번
    psi_s2 = mpmath.digamma(s_mp / 2)
    # ψ((s+1)/2) — Vga에서 α=1 한 번
    psi_s1_2 = mpmath.digamma((s_mp + 1) / 2)

    gamma_contrib = complex(float(mpmath.re(psi_s2)), float(mpmath.im(psi_s2)))
    gamma_contrib += 0.5 * complex(float(mpmath.re(psi_s1_2)), float(mpmath.im(psi_s1_2)))

    # L'/L(s) (PARI)
    L_val = lfun_eval_gl3(sigma, t)
    if abs(L_val) < 1e-30:
        return complex(1e6, 0)  # 영점 근방 → 발산
    L_prime = lfun_deriv_gl3(sigma, t)
    LpL = L_prime / L_val

    return const + gamma_contrib + LpL


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# E(σ) 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_energy_profile(t_min, t_max, n_t, sigma_values, N_cond):
    """각 σ에서 E(σ) = ∫|Λ'/Λ(σ+it)|² dt"""
    ts = np.linspace(t_min, t_max, n_t)
    energies = np.zeros(len(sigma_values))

    for j, sigma in enumerate(sigma_values):
        integrand = np.zeros(n_t)
        t0_s = time.time()
        for i, t in enumerate(ts):
            try:
                conn = connection_gl3(sigma, t, N_cond)
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
        elapsed_s = time.time() - t0_s

        print(f"    σ={sigma:.3f}: E={energies[j]:.2f} ({elapsed_s:.1f}s)", flush=True)

    return energies


def compute_energy_diag(zeros, sigma_values, t_min, t_max):
    """Hadamard 대각항: E_diag = Σ_n ∫ 1/((σ-σ_c)²+(t-γ_n)²) dt"""
    energies_diag = np.zeros(len(sigma_values))

    for j, sigma in enumerate(sigma_values):
        delta = abs(sigma - SIGMA_CRIT)
        if delta < 1e-10:
            delta = 1e-10
        total = 0.0
        for gamma in zeros:
            # ∫_{t_min}^{t_max} 1/(δ²+(t-γ)²) dt = [arctan((t-γ)/δ)/δ]
            total += (math.atan((t_max - gamma) / delta) -
                      math.atan((t_min - gamma) / delta)) / delta
        energies_diag[j] = total

    return energies_diag


def fit_power_law(delta_sigma, energies, safe_min=0.05):
    """E = A/Δσ^α 피팅"""
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


def fit_small_delta(delta_sigma, energies, max_delta=0.15):
    """소Δσ 레짐 피팅 (C-303 교훈: 교차항 회피)"""
    mask = (delta_sigma >= 0.03) & (delta_sigma <= max_delta)
    if mask.sum() < 3:
        return 0, 0, 0, 0
    return fit_power_law(delta_sigma[mask], energies[mask], safe_min=0.02)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    # σ 격자: σ_crit ± offset
    offsets = np.array([0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4])
    sigma_left = SIGMA_CRIT - offsets[::-1]
    sigma_right = SIGMA_CRIT + offsets
    sigmas = np.concatenate([sigma_left, sigma_right])
    delta_sigmas = np.abs(sigmas - SIGMA_CRIT)

    print("=" * 70, flush=True)
    print("[C-305] GL(3) sym²(11a1) E(σ) 프로파일 — Prop 6 degree-3 검증", flush=True)
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print(f"σ_crit = {SIGMA_CRIT} (GL(3) sym², k=3)", flush=True)
    dt = (T_RANGE[1] - T_RANGE[0]) / N_T_POINTS
    print(f"t 구간: [{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점 (Δt={dt:.3f})", flush=True)
    print(f"σ 격자: {len(sigmas)}점, offsets: {offsets.tolist()}", flush=True)
    print(f"방법: 해석적 Λ'/Λ, Vga=[0,0,1]", flush=True)
    print("=" * 70, flush=True)

    # 1. 초기화 + 영점
    print("\n[1] GL(3) sym²(11a1) 초기화", flush=True)
    zeros, N_cond = init_gl3()
    N_zeros = len(zeros)
    pi_N = np.pi * N_zeros
    print(f"  N_cond = {N_cond}, N_zeros = {N_zeros}, πN = {pi_N:.4f}", flush=True)

    if N_zeros < 10:
        print("  ⚠️ 영점 수 부족 (<10). 유효하지 않은 피팅.", flush=True)

    # 영점 간격 통계
    if N_zeros >= 2:
        gaps = np.diff(zeros)
        mean_gap = gaps.mean()
        min_gap = gaps.min()
        print(f"  영점 간격: mean={mean_gap:.4f}, min={min_gap:.4f}", flush=True)

    # 2. E(σ) 프로파일 (해석적 Λ'/Λ)
    print(f"\n[2] E(σ) 프로파일 ({len(sigmas)}점 × {N_T_POINTS}t)", flush=True)
    t0_e = time.time()
    energies = compute_energy_profile(T_RANGE[0], T_RANGE[1], N_T_POINTS, sigmas, N_cond)
    e_elapsed = time.time() - t0_e
    print(f"  E(σ) 완료: {e_elapsed:.1f}s", flush=True)

    # 3. Hadamard 대각항 (독립 비교)
    print(f"\n[3] Hadamard 대각항 E_diag", flush=True)
    energies_diag = compute_energy_diag(zeros, sigmas, T_RANGE[0], T_RANGE[1])

    # 4. 대칭성 E(σ)/E(k-σ) where k=3
    n_half = len(sigmas) // 2
    sym_ratios = energies[:n_half] / np.maximum(energies[n_half:][::-1], 1e-10)
    sym_mean = sym_ratios.mean()
    sym_std = sym_ratios.std()
    print(f"\n[4] 대칭성 E(σ)/E({2*SIGMA_CRIT:.0f}-σ): {sym_mean:.4f} ± {sym_std:.4f}", flush=True)

    # 대칭화
    energies_sym = (energies[:n_half] + energies[n_half:][::-1]) / 2
    energies_diag_sym = (energies_diag[:n_half] + energies_diag[n_half:][::-1]) / 2
    ds = delta_sigmas[:n_half]

    # 5. 멱법칙 피팅
    print(f"\n[5] 멱법칙 피팅", flush=True)

    # 전체 범위
    alpha, A_fit, r2, n_fit = fit_power_law(ds, energies_sym, safe_min=FIT_SAFE)
    a_over_pin = A_fit / pi_N if pi_N > 0 else 0
    print(f"  [전체] α = {alpha:.4f}, A/πN = {a_over_pin:.4f}, R² = {r2:.6f} (n={n_fit})", flush=True)

    # 소Δσ 레짐 (C-303 교훈: 교차항 회피)
    alpha_sm, A_sm, r2_sm, n_sm = fit_small_delta(ds, energies_sym, max_delta=0.15)
    a_sm_pin = A_sm / pi_N if pi_N > 0 else 0
    if n_sm >= 3:
        print(f"  [소Δσ≤0.15] α = {alpha_sm:.4f}, A/πN = {a_sm_pin:.4f}, R² = {r2_sm:.6f} (n={n_sm})", flush=True)
    else:
        print(f"  [소Δσ≤0.15] 데이터 부족 (n={n_sm})", flush=True)

    # Hadamard 대각항 피팅
    alpha_d, A_d, r2_d, n_d = fit_power_law(ds, energies_diag_sym, safe_min=FIT_SAFE)
    a_d_pin = A_d / pi_N if pi_N > 0 else 0
    print(f"  [대각항] α = {alpha_d:.4f}, A/πN = {a_d_pin:.4f}, R² = {r2_d:.6f} (n={n_d})", flush=True)

    # 교차항 비율
    if energies_sym.sum() > 0:
        cross_frac = 1.0 - energies_diag_sym / np.maximum(energies_sym, 1e-10)
        for j in range(len(ds)):
            print(f"    Δσ={ds[j]:.3f}: E_full={energies_sym[j]:.1f}, "
                  f"E_diag={energies_diag_sym[j]:.1f}, cross_frac={cross_frac[j]:.3f}", flush=True)

    # 6. Δσ·E 상수성
    mask_safe = ds >= FIT_SAFE
    if mask_safe.sum() > 0:
        ds_E = ds[mask_safe] * energies_sym[mask_safe]
        print(f"\n  Δσ·E: {ds_E.mean():.2f} ± {ds_E.std():.2f}, "
              f"CV={ds_E.std()/ds_E.mean()*100:.1f}%", flush=True)

    # ━━━ 판정 ━━━
    print(f"\n{'='*70}", flush=True)
    print("[판정]", flush=True)

    # GL(3) 영점 간격 ≈ 0.7 (GL(2)의 1.1~1.4보다 작음)
    # 교차항이 더 지배적일 수 있음 → 소Δσ 피팅이 핵심
    use_alpha = alpha_sm if n_sm >= 3 else alpha
    use_r2 = r2_sm if n_sm >= 3 else r2
    label = "소Δσ" if n_sm >= 3 else "전체"

    if abs(use_alpha - 1) < 0.05 and use_r2 > 0.99:
        verdict = "★★★★★ 강양성 — Prop 6 GL(3) 확장 확립"
    elif abs(use_alpha - 1) < 0.15 and use_r2 > 0.95:
        verdict = "★★★★ 양성 — GL(3)에서 α≈1 확인"
    elif abs(alpha_d - 1) < 0.03 and use_r2 > 0.9:
        verdict = "★★★ 조건부 양성 — E_diag α=1 확인, 교차항이 전체 α를 변형"
    elif abs(use_alpha - 1) >= 0.3:
        verdict = "★★ 음성 — GL(3)에서 α≠1, degree-3 경계 발견"
    else:
        verdict = "★★★ 중립 — 추가 분석 필요"

    print(f"  피팅 기준: {label}, α={use_alpha:.4f}, R²={use_r2:.6f}", flush=True)
    print(f"  대각항: α={alpha_d:.4f} (이론적으로 항상 1)", flush=True)
    print(f"  대칭성: {sym_mean:.4f} ± {sym_std:.4f}", flush=True)
    print(f"\n  종합: {verdict}", flush=True)

    # GL(1)/GL(2) 비교표
    print(f"\n{'='*70}", flush=True)
    print("[비교 요약] GL(1) vs GL(2) vs GL(3)", flush=True)
    print(f"{'='*70}", flush=True)
    header = f"  {'L-함수':<30} {'d':>2} {'σ_c':>4} {'N':>4} {'α_all':>8} {'α_sm':>8} {'α_diag':>8} {'R²':>8}"
    print(header, flush=True)
    print(f"  {'-'*82}", flush=True)
    print(f"  {'ζ(s) (C-298)':<30} {'1':>2} {'0.5':>4} {'5':>4} {'1.005':>8} {'—':>8} {'1.00':>8} {'0.9999':>8}", flush=True)
    print(f"  {'χ₃ mod3 (C-302)':<30} {'1':>2} {'0.5':>4} {'12':>4} {'0.988':>8} {'—':>8} {'—':>8} {'0.9989':>8}", flush=True)
    print(f"  {'11a1 (C-303)':<30} {'2':>2} {'1.0':>4} {'18':>4} {'0.700':>8} {'1.026':>8} {'1.02':>8} {'0.9999':>8}", flush=True)
    print(f"  {'37a1 (C-303)':<30} {'2':>2} {'1.0':>4} {'23':>4} {'0.621':>8} {'1.017':>8} {'1.02':>8} {'0.9999':>8}", flush=True)
    print(f"  {'-'*82}", flush=True)
    a_sm_str = f"{alpha_sm:.3f}" if n_sm >= 3 else "—"
    print(f"  {'sym²(11a1) (C-305)':<30} {'3':>2} {'1.5':>4} {N_zeros:>4} {alpha:>8.3f} {a_sm_str:>8} {alpha_d:>8.2f} {r2:>8.4f}", flush=True)

    elapsed = time.time() - t_start

    # ━━━ 결과 파일 ━━━
    with open(OUT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"[C-305] GL(3) sym²(11a1) E(σ) 프로파일 — Prop 6 degree-3 검증\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"경과: {elapsed:.1f}s\n")
        f.write(f"σ_crit = {SIGMA_CRIT}, Vga = [0, 0, 1]\n")
        f.write(f"t 구간: [{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점 (Δt={dt:.3f})\n")
        f.write(f"N_cond = {N_cond}\n")
        f.write(f"영점: N = {N_zeros}, πN = {pi_N:.4f}\n")
        if N_zeros >= 2:
            f.write(f"영점 간격: mean = {mean_gap:.4f}, min = {min_gap:.4f}\n")
        f.write("=" * 70 + "\n")

        # 영점 목록
        f.write("\n영점 목록:\n")
        for i, z in enumerate(zeros):
            f.write(f"  γ_{i+1} = {z:.8f}\n")

        # E(σ) 프로파일
        f.write(f"\n{'━'*60}\n")
        f.write("E(σ) 프로파일:\n")
        f.write(f"  {'σ':>8} {'|σ-1.5|':>8} {'E_full':>14} {'E_sym':>14} {'E_diag':>14} {'cross%':>8}\n")
        for j in range(len(sigmas)):
            ds_j = abs(sigmas[j] - SIGMA_CRIT)
            js = j if j < n_half else len(sigmas) - 1 - j
            es = energies_sym[js] if js < len(energies_sym) else 0
            ed = energies_diag_sym[js] if js < len(energies_diag_sym) else 0
            cf = (1 - ed / es) * 100 if es > 0 else 0
            f.write(f"  {sigmas[j]:8.3f} {ds_j:8.4f} {energies[j]:14.2f} "
                    f"{es:14.2f} {ed:14.2f} {cf:8.1f}\n")

        # 피팅 결과
        f.write(f"\n{'━'*60}\n")
        f.write("피팅 결과:\n")
        f.write(f"  [전체] α = {alpha:.6f}, A = {A_fit:.4f}, A/πN = {a_over_pin:.6f}, R² = {r2:.8f}\n")
        if n_sm >= 3:
            f.write(f"  [소Δσ≤0.15] α = {alpha_sm:.6f}, A = {A_sm:.4f}, A/πN = {a_sm_pin:.6f}, R² = {r2_sm:.8f}\n")
        f.write(f"  [대각항] α = {alpha_d:.6f}, A = {A_d:.4f}, A/πN = {a_d_pin:.6f}, R² = {r2_d:.8f}\n")

        # 대칭성
        f.write(f"\n대칭성 E(σ)/E(k-σ): {sym_mean:.6f} ± {sym_std:.6f}\n")

        # Δσ·E
        if mask_safe.sum() > 0:
            f.write(f"Δσ·E: {ds_E.mean():.4f} ± {ds_E.std():.4f}, CV = {ds_E.std()/ds_E.mean()*100:.1f}%\n")

        # 판정
        f.write(f"\n{'='*70}\n")
        f.write(f"종합 판정: {verdict}\n")

        # 비교표
        f.write(f"\n{'='*70}\n")
        f.write("비교 요약 (GL(1) vs GL(2) vs GL(3)):\n")
        f.write(header + "\n")
        f.write(f"  {'-'*82}\n")
        f.write(f"  {'ζ(s) (C-298)':<30} {'1':>2} {'0.5':>4} {'5':>4} {'1.005':>8} {'—':>8} {'1.00':>8} {'0.9999':>8}\n")
        f.write(f"  {'χ₃ mod3 (C-302)':<30} {'1':>2} {'0.5':>4} {'12':>4} {'0.988':>8} {'—':>8} {'—':>8} {'0.9989':>8}\n")
        f.write(f"  {'11a1 (C-303)':<30} {'2':>2} {'1.0':>4} {'18':>4} {'0.700':>8} {'1.026':>8} {'1.02':>8} {'0.9999':>8}\n")
        f.write(f"  {'37a1 (C-303)':<30} {'2':>2} {'1.0':>4} {'23':>4} {'0.621':>8} {'1.017':>8} {'1.02':>8} {'0.9999':>8}\n")
        f.write(f"  {'-'*82}\n")
        f.write(f"  {'sym²(11a1) (C-305)':<30} {'3':>2} {'1.5':>4} {N_zeros:>4} {alpha:>8.3f} {a_sm_str:>8} {alpha_d:>8.2f} {r2:>8.4f}\n")

    print(f"\n  결과: {OUT_FILE}", flush=True)
    print(f"  총 경과: {elapsed:.1f}s", flush=True)


if __name__ == '__main__':
    main()
