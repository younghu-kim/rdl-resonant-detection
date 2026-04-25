"""
=============================================================================
[C-308b] GL(4) sym³(11a1) E(σ) — Gamma 배경 분리 + 순수 L'/L 에너지
=============================================================================
C-308 결과: E_full ≈ 상수 (α≈0). 원인: 4개 digamma 항이 Gamma 배경 지배.
수정: E_gamma (순수 Gamma 기여) 분리, E_LpL = ∫|L'/L|² dt 별도 계산.

가설: E_LpL(σ) = πN/|σ-σ_c| + lower order (Hadamard 대각항과 일치)
버그 수정: ds 정렬 (C-305와 동일한 delta_sigmas[:n_half] 사용)
=============================================================================
"""

import sys, os, time, math
import numpy as np
import mpmath

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(3_000_000_000)
pari.set_real_precision(50)

T_RANGE = (3.0, 45.0)
N_T_POINTS = 100          # GL(4) PARI 호출 느림 → 100점 (Δt≈0.42, 예상 ~22분)
FIT_SAFE = 0.03

OUT_FILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl4_sigma_profile_c308b.txt')
os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)

VGA = []
N_COND = 1
SIGMA_CRIT = 2.0
DEGREE = 4


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def init_gl4():
    global VGA, N_COND, SIGMA_CRIT, DEGREE
    pari("E11 = ellinit([0,-1,1,-10,-20])")
    print(f"  [11a1] conductor = {int(str(pari('ellglobalred(E11)[1]')))}", flush=True)

    pari("L4 = lfunsympow(E11, 3)")
    print("  lfunsympow(E11, 3) 완료", flush=True)

    # Selberg 데이터 — lfunparams 실패 시 기본값 사용
    try:
        vga_str = str(pari("lfunparams(L4)[2]"))
        k_str = str(pari("lfunparams(L4)[3]"))
        N_cond_str = str(pari("lfunparams(L4)[4]"))
    except Exception:
        vga_str = "[0, 1, 1, 2]"
        k_str = "4"
        N_cond_str = "1331"

    vga_clean = vga_str.strip().strip('[]').replace('~', '')
    VGA = [int(float(x.strip())) for x in vga_clean.split(',') if x.strip()]
    DEGREE = len(VGA)

    try:
        N_COND = int(float(N_cond_str.strip()))
    except Exception:
        N_COND = 1331

    try:
        SIGMA_CRIT = int(float(k_str.strip())) / 2.0
    except Exception:
        SIGMA_CRIT = 2.0

    print(f"  degree={DEGREE}, Vga={VGA}, N_cond={N_COND}, σ_crit={SIGMA_CRIT}", flush=True)

    # lfuninit + zeros
    pari(f"L4init = lfuninit(L4, [{int(T_RANGE[1]) + 5}, 0])")
    pari(f"zv4 = lfunzeros(L4init, {T_RANGE[1]})")
    n_z = int(str(pari("#zv4")))
    zeros = []
    for i in range(1, n_z + 1):
        t = pf(pari(f"zv4[{i}]"))
        if not math.isnan(t) and t >= T_RANGE[0]:
            zeros.append(t)
    zeros = sorted(zeros)
    print(f"  영점: {len(zeros)}개 (t∈[{T_RANGE[0]}, {T_RANGE[1]}])", flush=True)
    return zeros


def gamma_contrib_only(sigma, t):
    """순수 Gamma 기여: (1/2)log(N) - (d/2)log(π) + Σ (1/2)ψ((s+μ)/2)"""
    const = 0.5 * math.log(max(N_COND, 1)) - 0.5 * DEGREE * math.log(math.pi)
    mpmath.mp.dps = 30
    s_mp = mpmath.mpc(sigma, t)
    gc = mpmath.mpc(0, 0)
    for mu in VGA:
        gc += 0.5 * mpmath.digamma((s_mp + mu) / 2)
    return const + complex(float(mpmath.re(gc)), float(mpmath.im(gc)))


def LpL_only(sigma, t):
    """순수 L'/L(s) — PARI 변수 캐싱 (L 1회, L' 1회 평가)"""
    s_r = float(sigma)
    s_i = float(t)
    try:
        pari(f"sv4 = {s_r} + {s_i}*I")
        pari("Lv4c = lfun(L4, sv4)")        # L(s) 1회
        re_L = pf(pari("real(Lv4c)"))
        im_L = pf(pari("imag(Lv4c)"))
        if math.isnan(re_L) or math.isnan(im_L):
            return complex(0, 0)
        L_val = complex(re_L, im_L)
        if abs(L_val) < 1e-30:
            return complex(1e10, 0)

        pari("Ld4c = lfun(L4, sv4, 1)")     # L'(s) 1회
        re_Ld = pf(pari("real(Ld4c)"))
        im_Ld = pf(pari("imag(Ld4c)"))
        if math.isnan(re_Ld) or math.isnan(im_Ld):
            return complex(0, 0)
        return complex(re_Ld, im_Ld) / L_val
    except Exception as e:
        print(f"WARNING: LpL_only σ={sigma:.4f} t={t:.3f} 오류: {e}", flush=True)
        return complex(0, 0)


def connection_gl4(sigma, t):
    """Λ'/Λ(s) = Gamma 기여 + L'/L(s)"""
    return gamma_contrib_only(sigma, t) + LpL_only(sigma, t)


def compute_three_energies(t_min, t_max, n_t, sigma_values):
    """E_full, E_gamma, E_LpL을 한 번에 계산"""
    ts = np.linspace(t_min, t_max, n_t)
    n_sig = len(sigma_values)
    E_full = np.zeros(n_sig)
    E_gamma = np.zeros(n_sig)
    E_LpL = np.zeros(n_sig)

    for j, sigma in enumerate(sigma_values):
        int_full = np.zeros(n_t)
        int_gamma = np.zeros(n_t)
        int_LpL = np.zeros(n_t)
        t0_s = time.time()

        for i, tv in enumerate(ts):
            tv_f = float(tv)
            try:
                gc = gamma_contrib_only(sigma, tv_f)
                lpl = LpL_only(sigma, tv_f)
                full = gc + lpl

                v_full = abs(full)**2
                v_gamma = abs(gc)**2
                v_lpl = abs(lpl)**2

                int_full[i] = v_full if np.isfinite(v_full) else 0.0
                int_gamma[i] = v_gamma if np.isfinite(v_gamma) else 0.0
                int_LpL[i] = v_lpl if np.isfinite(v_lpl) else 0.0
            except Exception as e:
                print(f"WARNING: t={tv_f:.3f} σ={sigma:.4f} 계산 오류: {e}", flush=True)

        # 극값 클리핑 (영점 근방 L'/L 발산)
        for arr in [int_full, int_LpL]:
            fv = arr[arr < 1e12]
            if len(fv) > 10:
                cap = np.percentile(fv, 99)
                arr[:] = np.minimum(arr, 10 * cap)

        E_full[j] = np.trapezoid(int_full, ts)
        E_gamma[j] = np.trapezoid(int_gamma, ts)
        E_LpL[j] = np.trapezoid(int_LpL, ts)

        elapsed_s = time.time() - t0_s
        print(f"    σ={sigma:.4f}: E_full={E_full[j]:.1f}, E_γ={E_gamma[j]:.1f}, "
              f"E_LpL={E_LpL[j]:.1f} ({elapsed_s:.1f}s)", flush=True)

    return E_full, E_gamma, E_LpL


def compute_energy_diag(zeros, sigma_values, t_min, t_max):
    """Hadamard 대각항"""
    ed = np.zeros(len(sigma_values))
    for j, sigma in enumerate(sigma_values):
        delta = abs(sigma - SIGMA_CRIT)
        if delta < 1e-10:
            delta = 1e-10
        for gamma in zeros:
            ed[j] += (math.atan((t_max - gamma) / delta) -
                      math.atan((t_min - gamma) / delta)) / delta
    return ed


def fit_power_law(delta_sigma, energies, safe_min=0.03):
    mask = (delta_sigma >= safe_min) & (energies > 0)
    if mask.sum() < 3:
        return 0, 0, 0, 0
    log_x = np.log(delta_sigma[mask])
    log_y = np.log(energies[mask])
    A_mat = np.column_stack([np.ones(mask.sum()), log_x])
    coeffs = np.linalg.lstsq(A_mat, log_y, rcond=None)[0]
    alpha = -coeffs[1]
    A = np.exp(coeffs[0])
    y_pred = A_mat @ coeffs
    ss_res = np.sum((log_y - y_pred)**2)
    ss_tot = np.sum((log_y - log_y.mean())**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0
    return alpha, A, r2, int(mask.sum())


def fit_regime(ds, E, min_d, max_d):
    mask = (ds >= min_d) & (ds <= max_d) & (E > 0)
    if mask.sum() < 3:
        return 0, 0, 0, 0
    return fit_power_law(ds[mask], E[mask], safe_min=min_d * 0.9)


def main():
    t_start = time.time()
    print("=" * 70, flush=True)
    print("[C-308b] GL(4) sym³(11a1) E(σ) — Gamma 분리 + L'/L 에너지", flush=True)
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print("=" * 70, flush=True)

    # 1. 초기화
    print("\n[1] 초기화", flush=True)
    zeros = init_gl4()
    N_zeros = len(zeros)
    pi_N = np.pi * N_zeros
    print(f"  N_zeros={N_zeros}, πN={pi_N:.4f}", flush=True)

    mean_gap = max_gap = min_gap = 0
    if N_zeros >= 2:
        gaps = np.diff(zeros)
        mean_gap, min_gap, max_gap = gaps.mean(), gaps.min(), gaps.max()
        print(f"  gap: mean={mean_gap:.4f}, min={min_gap:.4f}, max={max_gap:.4f}", flush=True)
        print(f"  Δσ_max 예측 (gap/d): {mean_gap/DEGREE:.4f}", flush=True)

    # σ 격자
    offsets = np.array([0.02, 0.03, 0.04, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.40])
    sigma_left = SIGMA_CRIT - offsets[::-1]   # ascending
    sigma_right = SIGMA_CRIT + offsets        # ascending
    sigmas = np.concatenate([sigma_left, sigma_right])
    delta_sigmas = np.abs(sigmas - SIGMA_CRIT)
    n_half = len(offsets)

    dt = (T_RANGE[1] - T_RANGE[0]) / N_T_POINTS
    print(f"  t: [{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점 (Δt={dt:.3f})", flush=True)
    print(f"  σ: {len(sigmas)}점", flush=True)

    # 2. E(σ) 3분해 (full, gamma, L'/L)
    print(f"\n[2] E(σ) 3분해 ({len(sigmas)}σ × {N_T_POINTS}t)", flush=True)
    t0_e = time.time()
    E_full, E_gamma, E_LpL = compute_three_energies(
        T_RANGE[0], T_RANGE[1], N_T_POINTS, sigmas)
    print(f"  완료: {time.time()-t0_e:.1f}s", flush=True)

    # 3. Hadamard 대각항
    print(f"\n[3] Hadamard 대각항", flush=True)
    E_diag = compute_energy_diag(zeros, sigmas, T_RANGE[0], T_RANGE[1])

    # 4. 대칭화 — ds는 delta_sigmas[:n_half] 사용 (C-305와 동일)
    # sigma_left[0] = σ_c - 0.40, sigma_right[9] = σ_c + 0.40 → 둘 다 δ=0.40
    ds_raw = delta_sigmas[:n_half]  # [0.40, 0.30, ..., 0.02] descending
    sort_idx = np.argsort(ds_raw)  # ascending 정렬 인덱스
    ds = ds_raw[sort_idx]  # ascending [0.02, 0.03, ..., 0.40]
    print(f"  [ds 정렬] ds_raw[0]={ds_raw[0]:.3f} → ds[0]={ds[0]:.3f} (ascending 보장)", flush=True)

    E_full_sym_raw = (E_full[:n_half] + E_full[n_half:][::-1]) / 2
    E_gamma_sym_raw = (E_gamma[:n_half] + E_gamma[n_half:][::-1]) / 2
    E_LpL_sym_raw = (E_LpL[:n_half] + E_LpL[n_half:][::-1]) / 2
    E_diag_sym_raw = (E_diag[:n_half] + E_diag[n_half:][::-1]) / 2

    E_full_sym = E_full_sym_raw[sort_idx]
    E_gamma_sym = E_gamma_sym_raw[sort_idx]
    E_LpL_sym = E_LpL_sym_raw[sort_idx]
    E_diag_sym = E_diag_sym_raw[sort_idx]

    # E_net = E_full - E_gamma (간섭항 포함 순수 영점 기여)
    E_net = np.maximum(E_full_sym - E_gamma_sym, 1.0)

    # 대칭성 검사
    k_fe = 2 * SIGMA_CRIT
    sym_ratios = E_full[:n_half] / np.maximum(E_full[n_half:][::-1], 1e-10)
    sym_mean, sym_std = sym_ratios.mean(), sym_ratios.std()
    print(f"\n[4] 대칭성 E(σ)/E({k_fe:.0f}-σ): {sym_mean:.6f} ± {sym_std:.6f}", flush=True)

    # 5. 멱법칙 피팅 — 5개 에너지 지표
    print(f"\n[5] 멱법칙 피팅", flush=True)

    labels = ["E_full", "E_LpL", "E_net", "E_diag"]
    arrays = [E_full_sym, E_LpL_sym, E_net, E_diag_sym]

    results = {}
    for lbl, arr in zip(labels, arrays):
        # 전체
        a, A, r2, n = fit_power_law(ds, arr, safe_min=FIT_SAFE)
        apin = A / pi_N if pi_N > 0 else 0
        results[f"{lbl}_all"] = (a, apin, r2, n)
        print(f"  [{lbl} 전체] α={a:.4f}, A/πN={apin:.4f}, R²={r2:.6f} (n={n})", flush=True)

        # 소Δσ ≤ 0.10
        a2, A2, r2_2, n2 = fit_regime(ds, arr, 0.02, 0.10)
        apin2 = A2 / pi_N if pi_N > 0 else 0
        if n2 >= 3:
            results[f"{lbl}_sm"] = (a2, apin2, r2_2, n2)
            print(f"  [{lbl} ≤0.10] α={a2:.4f}, A/πN={apin2:.4f}, R²={r2_2:.6f} (n={n2})", flush=True)

        # 초소Δσ ≤ 0.05
        a3, A3, r2_3, n3 = fit_regime(ds, arr, 0.02, 0.05)
        apin3 = A3 / pi_N if pi_N > 0 else 0
        if n3 >= 3:
            results[f"{lbl}_xs"] = (a3, apin3, r2_3, n3)
            print(f"  [{lbl} ≤0.05] α={a3:.4f}, A/πN={apin3:.4f}, R²={r2_3:.6f} (n={n3})", flush=True)

    # 6. 데이터 테이블
    print(f"\n[6] 데이터 테이블", flush=True)
    print(f"  {'Δσ':>6} {'E_full':>10} {'E_γ':>10} {'E_LpL':>10} {'E_net':>10} {'E_diag':>10} {'LpL/diag':>8}", flush=True)
    for j in range(len(ds)):
        ratio = E_LpL_sym[j] / E_diag_sym[j] if E_diag_sym[j] > 0 else 0
        print(f"  {ds[j]:>6.3f} {E_full_sym[j]:>10.1f} {E_gamma_sym[j]:>10.1f} "
              f"{E_LpL_sym[j]:>10.1f} {E_net[j]:>10.1f} {E_diag_sym[j]:>10.1f} {ratio:>8.3f}", flush=True)

    # 교차항 분석 (E_LpL vs E_diag)
    print(f"\n[7] 교차항 (L'/L 내부)", flush=True)
    cross_lpl = 1.0 - E_diag_sym / np.maximum(E_LpL_sym, 1e-10)
    for j in range(len(ds)):
        sign = "+" if cross_lpl[j] >= 0 else ""
        print(f"    Δσ={ds[j]:.3f}: cross_frac={sign}{cross_lpl[j]*100:.1f}%", flush=True)

    # Gamma 지배도
    print(f"\n[8] Gamma 지배도", flush=True)
    gamma_frac = E_gamma_sym / np.maximum(E_full_sym, 1e-10)
    for j in range(len(ds)):
        print(f"    Δσ={ds[j]:.3f}: E_γ/E_full={gamma_frac[j]*100:.1f}%", flush=True)

    # ━━━ 판정 ━━━
    print(f"\n{'='*70}", flush=True)
    print("[판정]", flush=True)

    # E_LpL이 핵심 지표 (순수 L'/L 기여)
    best = results.get("E_LpL_xs", results.get("E_LpL_sm", results.get("E_LpL_all", (0,0,0,0))))
    lpl_alpha, lpl_apin, lpl_r2, lpl_n = best

    # E_diag은 항상 α=1 (수학적 확실)
    diag_all = results.get("E_diag_all", (0,0,0,0))

    best_label = "E_LpL"
    for regime in ["xs", "sm", "all"]:
        key = f"E_LpL_{regime}"
        if key in results and results[key][2] > 0.95:
            lpl_alpha, lpl_apin, lpl_r2, lpl_n = results[key]
            best_label = f"E_LpL ({regime})"
            break

    print(f"  E_LpL 최적: {best_label}", flush=True)
    print(f"    α = {lpl_alpha:.4f}, A/πN = {lpl_apin:.4f}, R² = {lpl_r2:.6f} (n={lpl_n})", flush=True)
    print(f"  E_diag 전체: α = {diag_all[0]:.4f}, R² = {diag_all[2]:.6f}", flush=True)
    print(f"  대칭성: {sym_mean:.6f} ± {sym_std:.6f}", flush=True)
    print(f"  Gamma 지배도: {gamma_frac.mean()*100:.1f}% (평균)", flush=True)

    if abs(lpl_alpha - 1) < 0.05 and lpl_r2 > 0.99:
        verdict = "★★★★★ 강양성 — E_LpL에서 α=1 확인, Prop 6 GL(4) 확장"
    elif abs(lpl_alpha - 1) < 0.15 and lpl_r2 > 0.95:
        verdict = "★★★★ 양성 — E_LpL에서 α≈1, Gamma 분리 필수"
    elif abs(diag_all[0] - 1) < 0.03:
        verdict = "★★★ 조건부 양성 — E_diag α=1 (자명), E_LpL은 교차항으로 왜곡"
    else:
        verdict = "★★ 음성 — GL(4)에서 E(σ) 프레임워크 실패"

    print(f"\n  종합: {verdict}", flush=True)

    # 비교표
    print(f"\n{'='*70}", flush=True)
    print("[비교] GL(1)~GL(4) E(σ) 보편성", flush=True)
    print(f"  {'d':>2} {'L-함수':<25} {'α_Λ':>6} {'α_LpL':>7} {'α_diag':>7} {'gap':>5} {'γ%':>5}", flush=True)
    print(f"  {'-'*60}", flush=True)
    print(f"  {'1':>2} {'ζ, χ':<25} {'1.00':>6} {'—':>7} {'1.00':>7} {'7.2':>5} {'~5':>5}", flush=True)
    print(f"  {'2':>2} {'11a1, 37a1':<25} {'0.70':>6} {'—':>7} {'1.02':>7} {'1.2':>5} {'~20':>5}", flush=True)
    print(f"  {'3':>2} {'sym²(11a1)':<25} {'0.55':>6} {'—':>7} {'1.01':>7} {'0.72':>5} {'~40':>5}", flush=True)
    a_lpl_str = f"{lpl_alpha:.2f}" if lpl_n >= 3 else "?"
    print(f"  {'4':>2} {'sym³(11a1)':<25} {'0.00':>6} {a_lpl_str:>7} {diag_all[0]:>7.2f} "
          f"{mean_gap:>5.2f} {gamma_frac.mean()*100:>5.0f}", flush=True)

    if N_zeros >= 2:
        print(f"\n[Δσ_max ≈ gap/d 검증]", flush=True)
        pred = mean_gap / DEGREE
        print(f"  예측: gap/d = {mean_gap:.3f}/{DEGREE} = {pred:.4f}", flush=True)
        for j in range(len(ds)):
            a_at, _, r2_at, n_at = fit_regime(ds, E_LpL_sym, 0.02, ds[j])
            if n_at >= 3 and r2_at > 0.5:
                st = "✅" if abs(a_at - 1) < 0.10 else "✗"
                print(f"    E_LpL Δσ≤{ds[j]:.3f}: α={a_at:.4f}, R²={r2_at:.4f} {st}", flush=True)

    elapsed = time.time() - t_start
    print(f"\n총 소요: {elapsed:.1f}s", flush=True)

    # ━━━ 결과 파일 ━━━
    with open(OUT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"[C-308b] GL(4) sym³(11a1) E(σ) — Gamma 분리\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"소요: {elapsed:.1f}s\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"[설정]\n")
        f.write(f"  degree={DEGREE}, Vga={VGA}, N_cond={N_COND}, σ_crit={SIGMA_CRIT}\n")
        f.write(f"  t=[{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점\n")
        f.write(f"  N_zeros={N_zeros}, πN={pi_N:.4f}\n")
        if N_zeros >= 2:
            f.write(f"  mean_gap={mean_gap:.4f}\n\n")

        f.write(f"[데이터]\n")
        f.write(f"  {'Δσ':>6} {'E_full':>10} {'E_γ':>10} {'E_LpL':>10} {'E_net':>10} {'E_diag':>10}\n")
        for j in range(len(ds)):
            f.write(f"  {ds[j]:>6.3f} {E_full_sym[j]:>10.1f} {E_gamma_sym[j]:>10.1f} "
                    f"{E_LpL_sym[j]:>10.1f} {E_net[j]:>10.1f} {E_diag_sym[j]:>10.1f}\n")
        f.write(f"\n")

        f.write(f"[피팅]\n")
        for key, val in results.items():
            f.write(f"  {key}: α={val[0]:.4f}, A/πN={val[1]:.4f}, R²={val[2]:.6f} (n={val[3]})\n")
        f.write(f"\n")

        f.write(f"[판정]\n")
        f.write(f"  {verdict}\n")
        f.write(f"  E_LpL: α={lpl_alpha:.4f}, R²={lpl_r2:.6f}\n")
        f.write(f"  E_diag: α={diag_all[0]:.4f}, R²={diag_all[2]:.6f}\n")
        f.write(f"  Gamma 지배도: {gamma_frac.mean()*100:.1f}%\n")

    print(f"결과 저장: {OUT_FILE}", flush=True)


if __name__ == "__main__":
    main()
