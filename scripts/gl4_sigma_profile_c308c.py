"""
=============================================================================
[C-308c] GL(4) sym³(11a1) E(σ) — Hadamard 직접 계산 (PARI lfun 미사용)
=============================================================================
C-308/b 교훈: GL(4) PARI lfun 점당 ~1s → 전체 2시간+. 비실용적.

전략: Hadamard 곱 표현
    Λ'/Λ(s) = b + Σ_ρ 1/(s-ρ) + 1/ρ

    영점 ρ_n = σ_c + iγ_n (임계선 위)에서:
    L'/L 기여 ≈ Σ_n 1/(s-ρ_n) = Σ_n 1/(δ + i(t-γ_n))

    E_Had(σ) = ∫|Σ_n 1/(δ+i(t-γ_n))|² dt  (= E_diag + E_cross)

이 계산은 PARI 불필요, Python-only, O(N_t × N_zeros) per σ → 수초 완료.

추가: Gamma 배경 분리 + E_full = |Gamma + Had|² 직접 비교
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
N_T_POINTS = 500         # 빠르므로 충분한 해상도
FIT_SAFE = 0.03

OUT_FILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl4_sigma_profile_c308c.txt')
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
    """sym³(11a1) PARI: Selberg 데이터 + 영점만 추출. lfun 평가 안함."""
    global VGA, N_COND, SIGMA_CRIT, DEGREE
    pari("E11 = ellinit([0,-1,1,-10,-20])")
    pari("L4 = lfunsympow(E11, 3)")

    # Selberg 데이터
    try:
        vga_str = str(pari("lfunparams(L4)[2]"))
        k_str = str(pari("lfunparams(L4)[3]"))
        N_cond_str = str(pari("lfunparams(L4)[4]"))
    except Exception:
        vga_str, k_str, N_cond_str = "[0, 1, 1, 2]", "4", "1331"

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
    return np.array(zeros)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Hadamard 직접 계산 (Python-only, vectorized)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_E_hadamard_vec(ts, zeros, delta):
    """
    E_Had(σ) = ∫|Σ_n 1/(δ+i(t-γ_n))|² dt

    벡터화: ts shape (N_t,), zeros shape (N_z,)
    """
    # tau[i, n] = t_i - γ_n
    tau = ts[:, None] - zeros[None, :]  # (N_t, N_z)
    denom = delta**2 + tau**2            # (N_t, N_z)

    # Σ_n δ/D_n and Σ_n τ_n/D_n
    re_sum = np.sum(delta / denom, axis=1)   # (N_t,)
    im_sum = np.sum(-tau / denom, axis=1)    # (N_t,)

    integrand = re_sum**2 + im_sum**2  # |f(t)|²
    return np.trapezoid(integrand, ts)


def compute_E_diag(zeros, delta, t_min, t_max):
    """Hadamard 대각항 (해석적): Σ_n [arctan((T-γ)/δ) - arctan((t-γ)/δ)] / δ"""
    if abs(delta) < 1e-12:
        delta = 1e-12
    total = 0.0
    for g in zeros:
        total += (math.atan((t_max - g) / delta) - math.atan((t_min - g) / delta)) / delta
    return total


def compute_E_gamma(ts, sigma, vga, n_cond, degree):
    """Gamma 배경 에너지: ∫|Gamma'/Gamma(σ+it)|² dt"""
    mpmath.mp.dps = 30
    const = 0.5 * math.log(max(n_cond, 1)) - 0.5 * degree * math.log(math.pi)

    integrand = np.zeros(len(ts))
    for i, t in enumerate(ts):
        s_mp = mpmath.mpc(sigma, float(t))
        gc = mpmath.mpc(0, 0)
        for mu in vga:
            gc += 0.5 * mpmath.digamma((s_mp + mu) / 2)
        gamma_c = const + complex(float(mpmath.re(gc)), float(mpmath.im(gc)))
        integrand[i] = abs(gamma_c)**2

    return np.trapezoid(integrand, ts)


def compute_E_full(ts, zeros, delta, sigma, vga, n_cond, degree):
    """
    E_full = ∫|Γ'/Γ + Σ 1/(s-ρ)|² dt  (Gamma + Hadamard 합)
    """
    mpmath.mp.dps = 30
    const = 0.5 * math.log(max(n_cond, 1)) - 0.5 * degree * math.log(math.pi)

    # Hadamard sum (vectorized)
    tau = ts[:, None] - zeros[None, :]
    denom = delta**2 + tau**2
    had_re = np.sum(delta / denom, axis=1)
    had_im = np.sum(-tau / denom, axis=1)

    integrand = np.zeros(len(ts))
    for i, t_val in enumerate(ts):
        s_mp = mpmath.mpc(sigma, float(t_val))
        gc = mpmath.mpc(0, 0)
        for mu in vga:
            gc += 0.5 * mpmath.digamma((s_mp + mu) / 2)
        gamma_c = const + complex(float(mpmath.re(gc)), float(mpmath.im(gc)))

        full_re = gamma_c.real + had_re[i]
        full_im = gamma_c.imag + had_im[i]
        integrand[i] = full_re**2 + full_im**2

    return np.trapezoid(integrand, ts)


def fit_power_law(ds, E, safe_min=0.03):
    mask = (ds >= safe_min) & (E > 0)
    if mask.sum() < 3:
        return 0, 0, 0, 0
    log_x = np.log(ds[mask])
    log_y = np.log(E[mask])
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
    print("[C-308c] GL(4) sym³(11a1) E(σ) — Hadamard 직접 계산", flush=True)
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
    print("=" * 70, flush=True)

    # 1. 초기화 (PARI: 영점 추출만)
    print("\n[1] 초기화", flush=True)
    zeros = init_gl4()
    N_zeros = len(zeros)
    pi_N = np.pi * N_zeros

    mean_gap = min_gap = max_gap = 0
    if N_zeros >= 2:
        gaps = np.diff(zeros)
        mean_gap, min_gap, max_gap = gaps.mean(), gaps.min(), gaps.max()
        print(f"  gap: mean={mean_gap:.4f}, min={min_gap:.4f}, max={max_gap:.4f}", flush=True)
        print(f"  Δσ_max 예측 (gap/d): {mean_gap/DEGREE:.4f}", flush=True)

    ts = np.linspace(T_RANGE[0], T_RANGE[1], N_T_POINTS)

    # σ 격자 — 더 세밀한 격자 (빠르므로)
    offsets = np.array([0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07,
                        0.10, 0.13, 0.15, 0.20, 0.30, 0.40, 0.50])
    n_off = len(offsets)

    print(f"  t: [{T_RANGE[0]}, {T_RANGE[1]}], {N_T_POINTS}점", flush=True)
    print(f"  σ: {2*n_off}점 (σ_crit ± {offsets.tolist()})", flush=True)

    # 2. 4가지 에너지 계산
    print(f"\n[2] 에너지 계산", flush=True)

    E_had = np.zeros(n_off)
    E_diag = np.zeros(n_off)
    E_gamma = np.zeros(n_off)
    E_full = np.zeros(n_off)

    for j, off in enumerate(offsets):
        delta = off
        sigma_r = SIGMA_CRIT + off  # 오른쪽만 계산 (대칭)

        t0 = time.time()
        E_had[j] = compute_E_hadamard_vec(ts, zeros, delta)
        E_diag[j] = compute_E_diag(zeros, delta, T_RANGE[0], T_RANGE[1])
        E_gamma[j] = compute_E_gamma(ts, sigma_r, VGA, N_COND, DEGREE)
        E_full[j] = compute_E_full(ts, zeros, delta, sigma_r, VGA, N_COND, DEGREE)
        elapsed = time.time() - t0

        cross_pct = (E_had[j] - E_diag[j]) / E_had[j] * 100 if E_had[j] > 0 else 0
        gamma_pct = E_gamma[j] / E_full[j] * 100 if E_full[j] > 0 else 0
        print(f"    Δσ={off:.3f}: E_had={E_had[j]:.1f}, E_diag={E_diag[j]:.1f}, "
              f"E_γ={E_gamma[j]:.1f}, E_full={E_full[j]:.1f} "
              f"(cross={cross_pct:+.1f}%, γ%={gamma_pct:.0f}%) [{elapsed:.1f}s]", flush=True)

    # E_net = E_full - E_gamma
    E_net = np.maximum(E_full - E_gamma, 1.0)

    # 3. 멱법칙 피팅
    print(f"\n[3] 멱법칙 피팅", flush=True)
    ds = offsets

    results = {}
    for lbl, arr in [("E_had", E_had), ("E_diag", E_diag), ("E_full", E_full),
                      ("E_gamma", E_gamma), ("E_net", E_net)]:
        a, A, r2, n = fit_power_law(ds, arr)
        apin = A / pi_N if pi_N > 0 else 0
        results[f"{lbl}_all"] = (a, apin, r2, n)
        print(f"  [{lbl} 전체] α={a:.4f}, A/πN={apin:.4f}, R²={r2:.6f} (n={n})", flush=True)

        # 소Δσ ≤ 0.10
        a2, A2, r2_2, n2 = fit_regime(ds, arr, 0.005, 0.10)
        apin2 = A2 / pi_N if pi_N > 0 else 0
        if n2 >= 3:
            results[f"{lbl}_sm"] = (a2, apin2, r2_2, n2)
            print(f"  [{lbl} ≤0.10] α={a2:.4f}, A/πN={apin2:.4f}, R²={r2_2:.6f} (n={n2})", flush=True)

        # 초소Δσ ≤ 0.05
        a3, A3, r2_3, n3 = fit_regime(ds, arr, 0.005, 0.05)
        apin3 = A3 / pi_N if pi_N > 0 else 0
        if n3 >= 3:
            results[f"{lbl}_xs"] = (a3, apin3, r2_3, n3)
            print(f"  [{lbl} ≤0.05] α={a3:.4f}, A/πN={apin3:.4f}, R²={r2_3:.6f} (n={n3})", flush=True)

    # 4. 데이터 테이블
    print(f"\n[4] 데이터", flush=True)
    print(f"  {'Δσ':>6} {'E_had':>10} {'E_diag':>10} {'cross%':>8} "
          f"{'E_γ':>10} {'E_full':>10} {'γ%':>6} {'E_net':>10}", flush=True)
    for j in range(n_off):
        cr = (E_had[j] - E_diag[j]) / E_had[j] * 100 if E_had[j] > 0 else 0
        gp = E_gamma[j] / E_full[j] * 100 if E_full[j] > 0 else 0
        print(f"  {ds[j]:>6.3f} {E_had[j]:>10.1f} {E_diag[j]:>10.1f} {cr:>+8.1f} "
              f"{E_gamma[j]:>10.1f} {E_full[j]:>10.1f} {gp:>6.0f} {E_net[j]:>10.1f}", flush=True)

    # 5. Gamma 지배 임계점 분석
    print(f"\n[5] Gamma 지배 임계점", flush=True)
    gamma_frac = E_gamma / np.maximum(E_full, 1e-10)
    had_frac = E_had / np.maximum(E_full, 1e-10)
    for j in range(n_off):
        dominant = "Gamma" if gamma_frac[j] > 0.5 else "Had"
        print(f"    Δσ={ds[j]:.3f}: γ/full={gamma_frac[j]:.3f}, Had/full={had_frac[j]:.3f} → {dominant}", flush=True)

    # 임계 Δσ: E_had = E_gamma 교차점
    cross_idx = None
    for j in range(1, n_off):
        if E_had[j-1] > E_gamma[j-1] and E_had[j] < E_gamma[j]:
            # 선형 보간
            r = (E_had[j-1] - E_gamma[j-1]) / ((E_had[j-1] - E_gamma[j-1]) - (E_had[j] - E_gamma[j]))
            ds_cross = ds[j-1] + r * (ds[j] - ds[j-1])
            cross_idx = j
            print(f"  ⚠️ E_Had = E_γ 교차점: Δσ ≈ {ds_cross:.4f}", flush=True)
            break

    if cross_idx is None and E_had[0] > E_gamma[0]:
        print(f"  E_Had > E_γ for all tested Δσ", flush=True)
    elif cross_idx is None:
        print(f"  E_Had < E_γ for all tested Δσ (Gamma 완전 지배)", flush=True)

    # ━━━ 판정 ━━━
    print(f"\n{'='*70}", flush=True)
    print("[판정]", flush=True)

    # E_had이 핵심 (순수 영점 기여, Hadamard 곱 표현)
    had_all = results.get("E_had_all", (0,0,0,0))
    had_sm = results.get("E_had_sm", had_all)
    had_xs = results.get("E_had_xs", had_sm)
    diag_all = results.get("E_diag_all", (0,0,0,0))
    full_xs = results.get("E_full_xs", results.get("E_full_sm", results.get("E_full_all", (0,0,0,0))))
    net_xs = results.get("E_net_xs", results.get("E_net_sm", results.get("E_net_all", (0,0,0,0))))

    # 최적 E_had 피팅
    for label, key in [("xs", "E_had_xs"), ("sm", "E_had_sm"), ("all", "E_had_all")]:
        if key in results and results[key][2] > 0.95:
            best_had = results[key]
            best_label = f"E_had ({label})"
            break
    else:
        best_had = had_all
        best_label = "E_had (all)"

    print(f"  {best_label}: α={best_had[0]:.4f}, A/πN={best_had[1]:.4f}, R²={best_had[2]:.6f}", flush=True)
    print(f"  E_diag (all): α={diag_all[0]:.4f}, A/πN={diag_all[1]:.4f}, R²={diag_all[2]:.6f}", flush=True)
    print(f"  E_full (xs): α={full_xs[0]:.4f}, R²={full_xs[2]:.6f}", flush=True)
    print(f"  E_net (xs): α={net_xs[0]:.4f}, R²={net_xs[2]:.6f}", flush=True)
    print(f"  Gamma 지배: {gamma_frac.mean()*100:.1f}% (평균)", flush=True)

    # E_had 교차항 비율
    cross_ratio = (E_had - E_diag) / np.maximum(E_had, 1e-10)
    print(f"  E_had 내 교차항: {cross_ratio.mean()*100:.1f}% (평균)", flush=True)

    if abs(best_had[0] - 1) < 0.05 and best_had[2] > 0.99:
        verdict = "★★★★★ 강양성 — E_Had α=1, Prop 6 GL(4) 확장 확립"
    elif abs(best_had[0] - 1) < 0.15 and best_had[2] > 0.95:
        verdict = "★★★★ 양성 — E_Had α≈1, Gamma 배경 분리 필수"
    elif abs(diag_all[0] - 1) < 0.03:
        if abs(full_xs[0] - 1) < 0.15:
            verdict = "★★★ 조건부 양성 — E_full 소Δσ에서 α≈1, E_diag 항상 α=1"
        else:
            verdict = "★★★ 조건부 양성 — E_diag α=1 (자명), E_full은 Gamma 지배"
    else:
        verdict = "★★ 음성 — GL(4)에서 프레임워크 실패"

    print(f"\n  종합: {verdict}", flush=True)

    # 비교표
    print(f"\n{'='*70}", flush=True)
    print("[비교] GL(1)~GL(4) E(σ) 보편성", flush=True)
    print(f"  {'d':>2} {'L-함수':<25} {'α_Λ':>6} {'α_Had':>7} {'α_diag':>7} {'gap':>5} {'γ%':>5}", flush=True)
    print(f"  {'-'*60}", flush=True)
    print(f"  {'1':>2} {'ζ, χ':<25} {'1.00':>6} {'—':>7} {'1.00':>7} {'7.2':>5} {'~5':>5}", flush=True)
    print(f"  {'2':>2} {'11a1, 37a1':<25} {'0.70':>6} {'—':>7} {'1.02':>7} {'1.2':>5} {'~20':>5}", flush=True)
    print(f"  {'3':>2} {'sym²(11a1)':<25} {'0.55':>6} {'—':>7} {'1.01':>7} {'0.72':>5} {'~40':>5}", flush=True)
    print(f"  {'4':>2} {'sym³(11a1)':<25} {full_xs[0]:>6.2f} {best_had[0]:>7.3f} {diag_all[0]:>7.3f} "
          f"{mean_gap:>5.2f} {gamma_frac.mean()*100:>5.0f}", flush=True)

    # Δσ_max ≈ gap/d 법칙 검증
    if N_zeros >= 2:
        pred = mean_gap / DEGREE
        print(f"\n[Δσ_max ≈ gap/d = {pred:.4f} 검증]", flush=True)
        for j in range(len(ds)):
            a_at, _, r2_at, n_at = fit_regime(ds, E_full, 0.005, ds[j])
            if n_at >= 3:
                st = "✅" if abs(a_at - 1) < 0.10 else "✗"
                print(f"  E_full Δσ≤{ds[j]:.3f}: α={a_at:.4f}, R²={r2_at:.4f} {st}", flush=True)

    elapsed = time.time() - t_start
    print(f"\n총 소요: {elapsed:.1f}s", flush=True)

    # ━━━ 결과 파일 ━━━
    with open(OUT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"[C-308c] GL(4) sym³(11a1) E(σ) — Hadamard 직접\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"소요: {elapsed:.1f}s\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"degree={DEGREE}, Vga={VGA}, N_cond={N_COND}, σ_crit={SIGMA_CRIT}\n")
        f.write(f"N_zeros={N_zeros}, πN={pi_N:.4f}, mean_gap={mean_gap:.4f}\n\n")

        f.write(f"[데이터]\n")
        f.write(f"  {'Δσ':>6} {'E_had':>10} {'E_diag':>10} {'cross%':>8} "
                f"{'E_γ':>10} {'E_full':>10} {'γ%':>6} {'E_net':>10}\n")
        for j in range(n_off):
            cr = (E_had[j] - E_diag[j]) / E_had[j] * 100 if E_had[j] > 0 else 0
            gp = E_gamma[j] / E_full[j] * 100 if E_full[j] > 0 else 0
            f.write(f"  {ds[j]:>6.3f} {E_had[j]:>10.1f} {E_diag[j]:>10.1f} {cr:>+8.1f} "
                    f"{E_gamma[j]:>10.1f} {E_full[j]:>10.1f} {gp:>6.0f} {E_net[j]:>10.1f}\n")
        f.write(f"\n")

        f.write(f"[피팅]\n")
        for key, val in sorted(results.items()):
            f.write(f"  {key}: α={val[0]:.4f}, A/πN={val[1]:.4f}, R²={val[2]:.6f} (n={val[3]})\n")
        f.write(f"\n")

        f.write(f"[판정]\n")
        f.write(f"  {verdict}\n")

    print(f"결과 저장: {OUT_FILE}", flush=True)


if __name__ == "__main__":
    main()
