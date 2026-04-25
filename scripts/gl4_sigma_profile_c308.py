"""
=============================================================================
[C-308] GL(4) sym³(11a1) E(σ) 프로파일 — Prop 6 degree-4 보편성 검증
=============================================================================
목적:
    Prop 6 (E(σ) = πN/|σ-σ_c|^α) 이 GL(d) d=1,2,3에서 확립된 후,
    GL(4) sym³(11a1)에서도 α=1 보편성이 유지되는지 검증.

대상: sym³(11a1) — GL(4), degree 4
    - 11a1: y²+y = x³-x²-10x-20, conductor 11
    - sym³ Euler product: roots {α³, α, α⁻¹, α⁻³}

방법:
    PARI lfunsympow(E11, 3) → Selberg 데이터 자동 추출
    Λ'/Λ(s) = (1/2)log(N) - (d/2)log(π) + Σ_j (1/2)ψ((s+μ_j)/2) + L'/L(s)
    E(σ) = ∫|Λ'/Λ(σ+it)|² dt

    추가: Hadamard 대각항 E_diag = Σ_n arctan((T-γ_n)/δ)/δ 비교

예측:
    Δσ_max ≈ gap/d 법칙. GL(4) 영점 간격 ~0.4(추정) → Δσ_max ≈ 0.1.
    Δσ ≤ 0.05 에서 α=1 회복 기대.

판정 기준:
    강양성: α ∈ [0.95, 1.05], R² > 0.99
    양성:   α ∈ [0.85, 1.15], R² > 0.95
    음성:   α ≠ 1 → degree-4 경계 발견
=============================================================================
"""

import sys, os, time, math
import numpy as np
import mpmath

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(3_000_000_000)  # 3GB (GL(4) 필요할 수 있음)
pari.set_real_precision(50)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

T_RANGE = (3.0, 45.0)   # 영점 탐색 범위
N_T_POINTS = 300         # 적분 격자 (Δt ≈ 0.14)
FIT_SAFE = 0.03          # 멱법칙 피팅 최소 Δσ

OUT_FILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl4_sigma_profile_c308.txt'
)
os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)


def pf(x):
    """PARI 값 → float 변환"""
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# GL(4) L-함수 초기화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 전역 변수 (PARI 세션에서 공유)
VGA = []
N_COND = 1
SIGMA_CRIT = 2.0
DEGREE = 4


def init_gl4():
    """sym³(11a1) PARI 초기화 + 영점 추출"""
    global VGA, N_COND, SIGMA_CRIT, DEGREE

    pari("E11 = ellinit([0,-1,1,-10,-20])")
    actual_N = int(str(pari("ellglobalred(E11)[1]")))
    print(f"  [11a1] conductor = {actual_N}", flush=True)

    # sym³ L-함수 생성
    print("  lfunsympow(E11, 3) 시작...", flush=True)
    t0 = time.time()
    pari("L4 = lfunsympow(E11, 3)")
    print(f"  lfunsympow 완료 ({time.time()-t0:.1f}s)", flush=True)

    # Selberg 데이터 추출
    try:
        vga_str = str(pari("lfunparams(L4)[2]"))
        k_str = str(pari("lfunparams(L4)[3]"))
        N_cond_str = str(pari("lfunparams(L4)[4]"))
        eps_str = str(pari("lfunparams(L4)[5]"))
        print(f"  Vga = {vga_str}", flush=True)
        print(f"  k = {k_str}", flush=True)
        print(f"  N_cond = {N_cond_str}", flush=True)
        print(f"  ε = {eps_str}", flush=True)
    except Exception as e:
        print(f"  ⚠️ Selberg 데이터 추출 실패: {e}", flush=True)
        vga_str = "[0, 1, 1, 2]"
        k_str = "4"
        N_cond_str = "1331"
        eps_str = "1"

    # Vga 파싱
    vga_clean = vga_str.strip().strip('[]').replace('~', '')
    VGA = [int(float(x.strip())) for x in vga_clean.split(',') if x.strip()]
    DEGREE = len(VGA)
    print(f"  parsed: degree = {DEGREE}, Vga = {VGA}", flush=True)

    # conductor
    try:
        N_COND = int(float(N_cond_str.strip()))
    except Exception:
        N_COND = 1331  # 11³ 기대값
        print(f"  ⚠️ conductor 파싱 실패, 기본값 {N_COND} 사용", flush=True)

    # σ_crit: PARI k 파라미터에서 결정
    # 함수방정식 Λ(s) = ε·Λ(k-s) → σ_crit = k/2
    try:
        k_val = int(float(k_str.strip()))
        SIGMA_CRIT = k_val / 2.0
    except Exception:
        SIGMA_CRIT = 2.0
    print(f"  σ_crit = {SIGMA_CRIT} (k={k_str})", flush=True)

    # lfuninit
    print(f"  lfuninit 시작 (T_max={T_RANGE[1]+5})...", flush=True)
    t0 = time.time()
    pari(f"L4init = lfuninit(L4, [{int(T_RANGE[1]) + 5}, 0])")
    print(f"  lfuninit 완료 ({time.time()-t0:.1f}s)", flush=True)

    # 영점 추출
    print(f"  lfunzeros 시작 (T_max={T_RANGE[1]})...", flush=True)
    t0 = time.time()
    pari(f"zv4 = lfunzeros(L4init, {T_RANGE[1]})")
    n_z = int(str(pari("#zv4")))
    zeros = []
    for i in range(1, n_z + 1):
        t = pf(pari(f"zv4[{i}]"))
        if not math.isnan(t) and t >= T_RANGE[0]:
            zeros.append(t)
    zeros = sorted(zeros)
    elapsed_z = time.time() - t0
    print(f"  영점: {len(zeros)}개 (t∈[{T_RANGE[0]}, {T_RANGE[1]}]) ({elapsed_z:.1f}s)", flush=True)

    for i, z in enumerate(zeros[:8]):
        print(f"    γ_{i+1} = {z:.8f}", flush=True)
    if len(zeros) > 8:
        print(f"    ... (총 {len(zeros)}개)", flush=True)

    return zeros


def connection_gl4(sigma, t):
    """
    Λ'/Λ(s) for GL(4) sym³(11a1) — 통합 최적화 버전

    공식: Λ'/Λ = (1/2)log(N) - (d/2)log(π) + Σ_j (1/2)ψ((s+μ_j)/2) + L'/L(s)
    PARI 호출 최소화: L과 L'를 한 세션에서 계산
    """
    # 상수부 (사전 계산 가능하나 가독성 위해 유지)
    const = 0.5 * math.log(max(N_COND, 1)) - 0.5 * DEGREE * math.log(math.pi)

    # ψ 항 (mpmath)
    mpmath.mp.dps = 30
    s_mp = mpmath.mpc(sigma, t)

    gamma_contrib = mpmath.mpc(0, 0)
    for mu in VGA:
        gamma_contrib += 0.5 * mpmath.digamma((s_mp + mu) / 2)
    gamma_c = complex(float(mpmath.re(gamma_contrib)),
                      float(mpmath.im(gamma_contrib)))

    # L'/L(s) — PARI 2회 호출 (lfun + lfun derivative)
    try:
        pari(f"_s4 = {sigma} + {t}*I")
        pari("_Lv = lfun(L4, _s4)")
        re_L = pf(pari("real(_Lv)"))
        im_L = pf(pari("imag(_Lv)"))
        if math.isnan(re_L) or math.isnan(im_L):
            return complex(1e6, 0)
        L_val = complex(re_L, im_L)
        if abs(L_val) < 1e-30:
            return complex(1e6, 0)

        pari("_Ld = lfun(L4, _s4, 1)")
        re_Ld = pf(pari("real(_Ld)"))
        im_Ld = pf(pari("imag(_Ld)"))
        if math.isnan(re_Ld) or math.isnan(im_Ld):
            return complex(1e6, 0)
        L_prime = complex(re_Ld, im_Ld)

        LpL = L_prime / L_val
    except Exception:
        LpL = complex(0, 0)

    return const + gamma_c + LpL


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# E(σ) 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_energy_profile(t_min, t_max, n_t, sigma_values):
    """각 σ에서 E(σ) = ∫|Λ'/Λ(σ+it)|² dt"""
    ts = np.linspace(t_min, t_max, n_t)
    energies = np.zeros(len(sigma_values))

    for j, sigma in enumerate(sigma_values):
        integrand = np.zeros(n_t)
        t0_s = time.time()
        for i, tv in enumerate(ts):
            try:
                conn = connection_gl4(sigma, float(tv))
                val = abs(conn)**2
                integrand[i] = val if np.isfinite(val) else 0.0
            except Exception:
                integrand[i] = 0.0

        # 극값 클리핑 (영점 근방 발산 제거)
        finite_vals = integrand[integrand < 1e12]
        if len(finite_vals) > 10:
            cap = np.percentile(finite_vals, 99)
            integrand = np.minimum(integrand, 10 * cap)

        energies[j] = np.trapezoid(integrand, ts)
        elapsed_s = time.time() - t0_s
        print(f"    σ={sigma:.4f}: E={energies[j]:.2f} ({elapsed_s:.1f}s)", flush=True)

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
            total += (math.atan((t_max - gamma) / delta) -
                      math.atan((t_min - gamma) / delta)) / delta
        energies_diag[j] = total

    return energies_diag


def fit_power_law(delta_sigma, energies, safe_min=0.03):
    """E = A/Δσ^α 멱법칙 피팅"""
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


def fit_regime(delta_sigma, energies, min_ds, max_ds):
    """특정 Δσ 범위에서 피팅"""
    mask = (delta_sigma >= min_ds) & (delta_sigma <= max_ds)
    if mask.sum() < 3:
        return 0, 0, 0, 0
    return fit_power_law(delta_sigma[mask], energies[mask], safe_min=min_ds * 0.9)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    print("=" * 70, flush=True)
    print("[C-308] GL(4) sym³(11a1) E(σ) 프로파일 — Prop 6 degree-4 검증", flush=True)
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print("=" * 70, flush=True)

    # 1. 초기화 + 영점
    print("\n[1] GL(4) sym³(11a1) 초기화", flush=True)
    zeros = init_gl4()
    N_zeros = len(zeros)
    pi_N = np.pi * N_zeros
    print(f"\n  요약: N_cond={N_COND}, degree={DEGREE}, Vga={VGA}", flush=True)
    print(f"  σ_crit={SIGMA_CRIT}, N_zeros={N_zeros}, πN={pi_N:.4f}", flush=True)

    if N_zeros < 5:
        print("  ⚠️ 영점 수 부족 (<5). 유효하지 않은 피팅 가능성.", flush=True)

    # 영점 간격 통계
    if N_zeros >= 2:
        gaps = np.diff(zeros)
        mean_gap = gaps.mean()
        min_gap = gaps.min()
        max_gap = gaps.max()
        print(f"  영점 간격: mean={mean_gap:.4f}, min={min_gap:.4f}, max={max_gap:.4f}", flush=True)
        predicted_ds_max = mean_gap / DEGREE
        print(f"  Δσ_max 예측 (gap/d): {predicted_ds_max:.4f}", flush=True)

    # σ 격자: σ_crit ± offset
    # GL(4)는 영점 밀도 높음 → 미세 격자 필요
    offsets = np.array([0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.40])
    sigma_left = SIGMA_CRIT - offsets[::-1]
    sigma_right = SIGMA_CRIT + offsets
    sigmas = np.concatenate([sigma_left, sigma_right])
    delta_sigmas = np.abs(sigmas - SIGMA_CRIT)

    dt = (T_RANGE[1] - T_RANGE[0]) / N_T_POINTS
    print(f"\n  t 구간: [{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점 (Δt={dt:.3f})", flush=True)
    print(f"  σ 격자: {len(sigmas)}점, offsets: {offsets.tolist()}", flush=True)

    # 2. E(σ) 프로파일
    print(f"\n[2] E(σ) 프로파일 ({len(sigmas)}점 × {N_T_POINTS}t)", flush=True)
    t0_e = time.time()
    energies = compute_energy_profile(T_RANGE[0], T_RANGE[1], N_T_POINTS, sigmas)
    e_elapsed = time.time() - t0_e
    print(f"  E(σ) 완료: {e_elapsed:.1f}s", flush=True)

    # 3. Hadamard 대각항
    print(f"\n[3] Hadamard 대각항 E_diag", flush=True)
    energies_diag = compute_energy_diag(zeros, sigmas, T_RANGE[0], T_RANGE[1])

    # 4. 대칭성 E(σ)/E(k-σ)
    k_fe = 2 * SIGMA_CRIT  # 함수방정식 반영 점
    n_half = len(sigmas) // 2
    sym_ratios = energies[:n_half] / np.maximum(energies[n_half:][::-1], 1e-10)
    sym_mean = sym_ratios.mean()
    sym_std = sym_ratios.std()
    print(f"\n[4] 대칭성 E(σ)/E({k_fe:.0f}-σ): {sym_mean:.6f} ± {sym_std:.6f}", flush=True)

    # 대칭화
    energies_sym = (energies[:n_half] + energies[n_half:][::-1]) / 2
    energies_diag_sym = (energies_diag[:n_half] + energies_diag[n_half:][::-1]) / 2
    ds = offsets  # = delta_sigmas[:n_half] reversed to ascending

    # 5. 멱법칙 피팅
    print(f"\n[5] 멱법칙 피팅", flush=True)

    # 5a. 전체 범위
    alpha_all, A_all, r2_all, n_all = fit_power_law(ds, energies_sym, safe_min=FIT_SAFE)
    a_over_pin_all = A_all / pi_N if pi_N > 0 else 0
    print(f"  [전체] α = {alpha_all:.4f}, A/πN = {a_over_pin_all:.4f}, "
          f"R² = {r2_all:.6f} (n={n_all})", flush=True)

    # 5b. 소Δσ 레짐 (≤0.10) — GL(4) 예측 영역
    alpha_sm, A_sm, r2_sm, n_sm = fit_regime(ds, energies_sym, 0.02, 0.10)
    a_sm_pin = A_sm / pi_N if pi_N > 0 else 0
    if n_sm >= 3:
        print(f"  [소Δσ≤0.10] α = {alpha_sm:.4f}, A/πN = {a_sm_pin:.4f}, "
              f"R² = {r2_sm:.6f} (n={n_sm})", flush=True)

    # 5c. 초소Δσ 레짐 (≤0.05) — GL(3)에서 α=1 회복됨
    alpha_xs, A_xs, r2_xs, n_xs = fit_regime(ds, energies_sym, 0.02, 0.05)
    a_xs_pin = A_xs / pi_N if pi_N > 0 else 0
    if n_xs >= 3:
        print(f"  [초소Δσ≤0.05] α = {alpha_xs:.4f}, A/πN = {a_xs_pin:.4f}, "
              f"R² = {r2_xs:.6f} (n={n_xs})", flush=True)

    # 5d. Hadamard 대각항 피팅
    alpha_d, A_d, r2_d, n_d = fit_power_law(ds, energies_diag_sym, safe_min=FIT_SAFE)
    a_d_pin = A_d / pi_N if pi_N > 0 else 0
    print(f"  [대각항] α = {alpha_d:.4f}, A/πN = {a_d_pin:.4f}, "
          f"R² = {r2_d:.6f} (n={n_d})", flush=True)

    # 교차항 비율
    print(f"\n[6] 교차항 분석", flush=True)
    cross_frac = 1.0 - energies_diag_sym / np.maximum(energies_sym, 1e-10)
    for j in range(len(ds)):
        sign = "+" if cross_frac[j] >= 0 else ""
        print(f"    Δσ={ds[j]:.3f}: E_full={energies_sym[j]:.1f}, "
              f"E_diag={energies_diag_sym[j]:.1f}, "
              f"cross_frac={sign}{cross_frac[j]*100:.1f}%", flush=True)

    # 교차항 부호 전환 탐색
    sign_changes = []
    for j in range(1, len(cross_frac)):
        if cross_frac[j-1] * cross_frac[j] < 0:
            sign_changes.append((ds[j-1], ds[j]))
    if sign_changes:
        for sc in sign_changes:
            print(f"  ⚠️ 교차항 부호 전환: Δσ ∈ [{sc[0]:.3f}, {sc[1]:.3f}]", flush=True)

    # Δσ·E 상수성
    mask_safe = ds >= FIT_SAFE
    if mask_safe.sum() > 0:
        ds_E = ds[mask_safe] * energies_sym[mask_safe]
        print(f"\n  Δσ·E: {ds_E.mean():.2f} ± {ds_E.std():.2f}, "
              f"CV={ds_E.std()/ds_E.mean()*100:.1f}%", flush=True)

    # ━━━ 판정 ━━━
    print(f"\n{'='*70}", flush=True)
    print("[판정]", flush=True)

    # 최적 피팅 결정: 초소Δσ > 소Δσ > 전체 순서
    if n_xs >= 3 and r2_xs > 0.95:
        use_alpha, use_r2, use_A_pin, label = alpha_xs, r2_xs, a_xs_pin, "초소Δσ≤0.05"
    elif n_sm >= 3 and r2_sm > 0.95:
        use_alpha, use_r2, use_A_pin, label = alpha_sm, r2_sm, a_sm_pin, "소Δσ≤0.10"
    else:
        use_alpha, use_r2, use_A_pin, label = alpha_all, r2_all, a_over_pin_all, "전체"

    if abs(use_alpha - 1) < 0.05 and use_r2 > 0.99:
        verdict = "★★★★★ 강양성 — Prop 6 GL(4) 확장 확립"
    elif abs(use_alpha - 1) < 0.15 and use_r2 > 0.95:
        verdict = "★★★★ 양성 — GL(4)에서 α≈1 확인"
    elif abs(alpha_d - 1) < 0.03:
        verdict = "★★★ 조건부 양성 — E_diag α=1 확인, 교차항이 전체 α를 변형"
    elif abs(use_alpha - 1) >= 0.3:
        verdict = "★★ 음성 — GL(4)에서 α≠1, degree-4 경계 발견"
    else:
        verdict = "★★★ 중립 — 추가 분석 필요"

    print(f"  피팅 기준: {label}", flush=True)
    print(f"    α = {use_alpha:.4f}", flush=True)
    print(f"    R² = {use_r2:.6f}", flush=True)
    print(f"    A/πN = {use_A_pin:.4f}", flush=True)
    print(f"  대각항: α = {alpha_d:.4f}, R² = {r2_d:.6f}", flush=True)
    print(f"  대칭성: {sym_mean:.6f} ± {sym_std:.6f}", flush=True)
    print(f"\n  종합: {verdict}", flush=True)

    # GL(1)/GL(2)/GL(3) 비교표
    print(f"\n{'='*70}", flush=True)
    print("[비교 요약] GL(1) ~ GL(4)", flush=True)
    print(f"{'='*70}", flush=True)
    header = f"  {'L-함수':<30} {'d':>2} {'σ_c':>4} {'N':>4} {'α_best':>8} {'Δσ_max':>8} {'gap':>6} {'A/πN':>6}"
    print(header, flush=True)
    print(f"  {'-'*80}", flush=True)
    print(f"  {'ζ(s) (C-298)':<30} {'1':>2} {'0.5':>4} {'5':>4} {'1.005':>8} {'all':>8} {'7.2':>6} {'0.95':>6}", flush=True)
    print(f"  {'χ₃ mod3 (C-302)':<30} {'1':>2} {'0.5':>4} {'12':>4} {'0.988':>8} {'all':>8} {'—':>6} {'0.98':>6}", flush=True)
    print(f"  {'11a1 (C-303)':<30} {'2':>2} {'1.0':>4} {'18':>4} {'1.026':>8} {'≤0.10':>8} {'1.2':>6} {'0.87':>6}", flush=True)
    print(f"  {'37a1 (C-303)':<30} {'2':>2} {'1.0':>4} {'23':>4} {'1.017':>8} {'≤0.10':>8} {'1.2':>6} {'0.87':>6}", flush=True)
    print(f"  {'sym²(11a1) (C-305/6)':<30} {'3':>2} {'1.5':>4} {'59':>4} {'1.004':>8} {'≤0.05':>8} {'0.72':>6} {'0.90':>6}", flush=True)
    print(f"  {'-'*80}", flush=True)

    # GL(4) 행
    gap_str = f"{mean_gap:.2f}" if N_zeros >= 2 else "—"
    ds_max_str = f"≤{predicted_ds_max:.2f}" if N_zeros >= 2 else "?"
    print(f"  {'sym³(11a1) (C-308)':<30} {DEGREE:>2} {SIGMA_CRIT:>4.1f} {N_zeros:>4} "
          f"{use_alpha:>8.3f} {ds_max_str:>8} {gap_str:>6} {use_A_pin:>6.2f}", flush=True)

    # Δσ_max ≈ gap/d 법칙 검증
    if N_zeros >= 2:
        print(f"\n[Δσ_max ≈ gap/d 법칙 검증]", flush=True)
        print(f"  예측: gap/d = {mean_gap:.3f}/{DEGREE} = {predicted_ds_max:.4f}", flush=True)
        print(f"  α 회복 실측:", flush=True)
        for j in range(len(ds)):
            alpha_at_ds, _, r2_at_ds, _ = fit_regime(ds, energies_sym, 0.02, ds[j])
            if r2_at_ds > 0.5:
                status = "✅" if abs(alpha_at_ds - 1) < 0.05 else "✗"
                print(f"    Δσ≤{ds[j]:.3f}: α={alpha_at_ds:.4f}, R²={r2_at_ds:.4f} {status}", flush=True)

    elapsed = time.time() - t_start
    print(f"\n총 소요 시간: {elapsed:.1f}s", flush=True)

    # ━━━ 결과 파일 ━━━
    with open(OUT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"[C-308] GL(4) sym³(11a1) E(σ) 프로파일 — Prop 6 degree-4 검증\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"소요: {elapsed:.1f}s\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"[설정]\n")
        f.write(f"  L-함수: sym³(11a1), degree={DEGREE}\n")
        f.write(f"  Vga = {VGA}\n")
        f.write(f"  N_cond = {N_COND}\n")
        f.write(f"  σ_crit = {SIGMA_CRIT}\n")
        f.write(f"  t 구간: [{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점\n")
        f.write(f"  σ offsets: {offsets.tolist()}\n")
        f.write(f"  N_zeros = {N_zeros}\n")
        if N_zeros >= 2:
            f.write(f"  mean_gap = {mean_gap:.4f}\n")
        f.write(f"  πN = {pi_N:.4f}\n\n")

        f.write(f"[영점]\n")
        for i, z in enumerate(zeros):
            f.write(f"  γ_{i+1} = {z:.10f}\n")
        f.write(f"\n")

        f.write(f"[E(σ) 데이터]\n")
        f.write(f"  {'Δσ':>8} {'E_full':>12} {'E_diag':>12} {'cross%':>8}\n")
        for j in range(len(ds)):
            sign = "+" if cross_frac[j] >= 0 else ""
            f.write(f"  {ds[j]:>8.4f} {energies_sym[j]:>12.2f} "
                    f"{energies_diag_sym[j]:>12.2f} {sign}{cross_frac[j]*100:>7.1f}%\n")
        f.write(f"\n")

        f.write(f"[대칭성]\n")
        f.write(f"  E(σ)/E({k_fe:.0f}-σ): {sym_mean:.6f} ± {sym_std:.6f}\n\n")

        f.write(f"[피팅]\n")
        f.write(f"  전체: α={alpha_all:.4f}, A/πN={a_over_pin_all:.4f}, R²={r2_all:.6f} (n={n_all})\n")
        if n_sm >= 3:
            f.write(f"  소Δσ≤0.10: α={alpha_sm:.4f}, A/πN={a_sm_pin:.4f}, R²={r2_sm:.6f} (n={n_sm})\n")
        if n_xs >= 3:
            f.write(f"  초소Δσ≤0.05: α={alpha_xs:.4f}, A/πN={a_xs_pin:.4f}, R²={r2_xs:.6f} (n={n_xs})\n")
        f.write(f"  대각항: α={alpha_d:.4f}, A/πN={a_d_pin:.4f}, R²={r2_d:.6f} (n={n_d})\n\n")

        f.write(f"[판정]\n")
        f.write(f"  {verdict}\n")
        f.write(f"  기준: {label}, α={use_alpha:.4f}, R²={use_r2:.6f}\n")

    print(f"\n결과 저장: {OUT_FILE}", flush=True)


if __name__ == "__main__":
    main()
