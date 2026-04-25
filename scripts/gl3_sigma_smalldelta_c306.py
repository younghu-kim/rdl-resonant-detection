"""
=============================================================================
[C-306] GL(3) sym²(11a1) E(σ) 초소Δσ 보충 — α=1 회복 여부 확인
=============================================================================
목적:
    C-305에서 GL(3) 소Δσ(≤0.15) α=0.78. GL(2)에서는 소Δσ α=1.02.
    Δσ = 0.02, 0.03, 0.04를 추가 측정하여:
    (a) α가 1에 수렴하는지, (b) 음의 교차항이 유지되는지 확인.

    또한 C-305의 기존 7점(0.05~0.40)을 재사용하여 통합 피팅.

방법: C-305와 동일한 해석적 Λ'/Λ (PARI lfunsympow)
=============================================================================
"""

import sys, os, time, math
import numpy as np
import mpmath

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(2_000_000_000)
pari.set_real_precision(50)

SIGMA_CRIT = 1.5
T_RANGE = (5.0, 48.0)
N_T_POINTS = 300
N_COND = 121  # sym²(11a1)

OUT_FILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl3_sigma_smalldelta_c306.txt'
)


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


# PARI 초기화
print("=" * 70, flush=True)
print("[C-306] GL(3) sym²(11a1) E(σ) 초소Δσ 보충", flush=True)
print(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)
print("=" * 70, flush=True)

pari("E11 = ellinit([0,-1,1,-10,-20])")
pari("L3 = lfunsympow(E11, 2)")
pari(f"L3init = lfuninit(L3, [{int(T_RANGE[1]) + 5}, 0])")
print("  PARI 초기화 완료", flush=True)

# 영점 추출
pari(f"zv3 = lfunzeros(L3init, {T_RANGE[1]})")
n_z = int(str(pari("#zv3")))
zeros = []
for i in range(1, n_z + 1):
    t = pf(pari(f"zv3[{i}]"))
    if not math.isnan(t) and t >= T_RANGE[0]:
        zeros.append(t)
zeros = sorted(zeros)
N_zeros = len(zeros)
pi_N = np.pi * N_zeros
print(f"  영점: {N_zeros}개, πN = {pi_N:.4f}", flush=True)


def connection_gl3(sigma, t):
    """Λ'/Λ(s) for GL(3) sym²(11a1), Vga=[0,0,1]"""
    const = 0.5 * math.log(N_COND) - 1.5 * math.log(math.pi)
    mpmath.mp.dps = 30
    s_mp = mpmath.mpc(sigma, t)
    psi_s2 = mpmath.digamma(s_mp / 2)
    psi_s1_2 = mpmath.digamma((s_mp + 1) / 2)
    gamma_c = complex(float(mpmath.re(psi_s2)), float(mpmath.im(psi_s2)))
    gamma_c += 0.5 * complex(float(mpmath.re(psi_s1_2)), float(mpmath.im(psi_s1_2)))

    try:
        L_re = pf(pari(f"real(lfun(L3, {sigma} + {t}*I))"))
        L_im = pf(pari(f"imag(lfun(L3, {sigma} + {t}*I))"))
        Lp_re = pf(pari(f"real(lfun(L3, {sigma} + {t}*I, 1))"))
        Lp_im = pf(pari(f"imag(lfun(L3, {sigma} + {t}*I, 1))"))
        L_val = complex(L_re, L_im)
        Lp_val = complex(Lp_re, Lp_im)
        if abs(L_val) < 1e-30:
            return complex(1e6, 0)
        LpL = Lp_val / L_val
    except Exception:
        return complex(1e6, 0)

    return const + gamma_c + LpL


def compute_energy(sigma):
    """E(σ) = ∫|Λ'/Λ(σ+it)|² dt"""
    ts = np.linspace(T_RANGE[0], T_RANGE[1], N_T_POINTS)
    integrand = np.zeros(N_T_POINTS)
    for i, t in enumerate(ts):
        try:
            conn = connection_gl3(sigma, t)
            val = abs(conn)**2
            integrand[i] = val if np.isfinite(val) else 0.0
        except Exception:
            integrand[i] = 0.0

    finite_vals = integrand[integrand < 1e12]
    if len(finite_vals) > 10:
        cap = np.percentile(finite_vals, 99)
        integrand = np.minimum(integrand, 10 * cap)

    return np.trapezoid(integrand, ts)


def compute_energy_diag(delta):
    """E_diag = Σ_n ∫ 1/(δ²+(t-γ_n)²) dt"""
    if delta < 1e-10:
        delta = 1e-10
    total = 0.0
    for gamma in zeros:
        total += (math.atan((T_RANGE[1] - gamma) / delta) -
                  math.atan((T_RANGE[0] - gamma) / delta)) / delta
    return total


# 초소Δσ 측정
new_offsets = [0.02, 0.03, 0.04]
# C-305 기존 결과 (재사용)
c305_data = {
    0.05: 3394.30,
    0.07: 2496.44,
    0.10: 1868.29,
    0.15: 1439.97,
    0.20: 1265.62,
    0.30: 1136.40,
    0.40: 1094.26,
}

print(f"\n[1] 초소Δσ 측정 (Δσ = {new_offsets})", flush=True)
t_total = time.time()

new_results = {}
for offset in new_offsets:
    sigma_left = SIGMA_CRIT - offset
    sigma_right = SIGMA_CRIT + offset

    t0 = time.time()
    E_left = compute_energy(sigma_left)
    E_right = compute_energy(sigma_right)
    elapsed = time.time() - t0

    E_sym = (E_left + E_right) / 2.0
    E_diag = compute_energy_diag(offset)
    cross_frac = 1.0 - E_diag / E_sym if E_sym > 0 else 0
    sym_ratio = E_left / E_right if E_right > 0 else 0

    new_results[offset] = {
        'E_left': E_left, 'E_right': E_right, 'E_sym': E_sym,
        'E_diag': E_diag, 'cross_frac': cross_frac, 'sym_ratio': sym_ratio,
    }

    print(f"  Δσ={offset:.3f}: E_sym={E_sym:.2f}, E_diag={E_diag:.2f}, "
          f"cross={cross_frac:.3f}, sym={sym_ratio:.4f} ({elapsed:.1f}s)", flush=True)

# 통합 데이터
print(f"\n[2] 통합 피팅 (C-305 + C-306)", flush=True)

all_offsets = sorted(list(new_results.keys()) + list(c305_data.keys()))
all_E_sym = []
all_E_diag = []
for d in all_offsets:
    if d in new_results:
        all_E_sym.append(new_results[d]['E_sym'])
        all_E_diag.append(new_results[d]['E_diag'])
    else:
        all_E_sym.append(c305_data[d])
        all_E_diag.append(compute_energy_diag(d))

all_offsets = np.array(all_offsets)
all_E_sym = np.array(all_E_sym)
all_E_diag = np.array(all_E_diag)


def fit_power_law(ds, E, safe_min=0.02):
    mask = ds >= safe_min
    if mask.sum() < 3:
        return 0, 0, 0, 0
    log_x = np.log(ds[mask])
    log_y = np.log(np.maximum(E[mask], 1e-30))
    A_mat = np.column_stack([np.ones(mask.sum()), log_x])
    coeffs = np.linalg.lstsq(A_mat, log_y, rcond=None)[0]
    alpha = -coeffs[1]
    A = np.exp(coeffs[0])
    y_pred = A_mat @ coeffs
    ss_res = np.sum((log_y - y_pred)**2)
    ss_tot = np.sum((log_y - log_y.mean())**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0
    return alpha, A, r2, int(mask.sum())


# 전체 범위
alpha_all, A_all, r2_all, n_all = fit_power_law(all_offsets, all_E_sym)
print(f"  [전체 10점] α = {alpha_all:.4f}, A/πN = {A_all/pi_N:.4f}, "
      f"R² = {r2_all:.6f}", flush=True)

# 소Δσ (≤0.10)
mask_sm = all_offsets <= 0.10
if mask_sm.sum() >= 3:
    alpha_sm, A_sm, r2_sm, n_sm = fit_power_law(
        all_offsets[mask_sm], all_E_sym[mask_sm], safe_min=0.015)
    print(f"  [Δσ≤0.10, {n_sm}점] α = {alpha_sm:.4f}, A/πN = {A_sm/pi_N:.4f}, "
          f"R² = {r2_sm:.6f}", flush=True)

# 초소Δσ (≤0.05)
mask_vsm = all_offsets <= 0.05
if mask_vsm.sum() >= 3:
    alpha_vsm, A_vsm, r2_vsm, n_vsm = fit_power_law(
        all_offsets[mask_vsm], all_E_sym[mask_vsm], safe_min=0.015)
    print(f"  [Δσ≤0.05, {n_vsm}점] α = {alpha_vsm:.4f}, A/πN = {A_vsm/pi_N:.4f}, "
          f"R² = {r2_vsm:.6f}", flush=True)

# 대각항
alpha_d, A_d, r2_d, n_d = fit_power_law(all_offsets, all_E_diag)
print(f"  [대각항 10점] α = {alpha_d:.4f}, A/πN = {A_d/pi_N:.4f}, "
      f"R² = {r2_d:.6f}", flush=True)

# Δσ·E 일람
print(f"\n[3] Δσ·E 상수성 일람", flush=True)
print(f"  {'Δσ':>6} {'E_full':>10} {'E_diag':>10} {'Δσ·E':>10} {'Δσ·E/πN':>10} {'cross%':>8}", flush=True)
for j in range(len(all_offsets)):
    d = all_offsets[j]
    ef = all_E_sym[j]
    ed = all_E_diag[j]
    cf = (1 - ed / ef) * 100 if ef > 0 else 0
    print(f"  {d:6.3f} {ef:10.2f} {ed:10.2f} {d*ef:10.2f} {d*ef/pi_N:10.4f} {cf:8.1f}", flush=True)

# 판정
print(f"\n{'='*70}", flush=True)
print("[판정]", flush=True)

if mask_vsm.sum() >= 3:
    use_alpha = alpha_vsm
    use_r2 = r2_vsm
    label = f"초소Δσ≤0.05 ({n_vsm}점)"
elif mask_sm.sum() >= 3:
    use_alpha = alpha_sm
    use_r2 = r2_sm
    label = f"소Δσ≤0.10 ({n_sm}점)"
else:
    use_alpha = alpha_all
    use_r2 = r2_all
    label = "전체"

if abs(use_alpha - 1) < 0.05 and use_r2 > 0.99:
    verdict = "★★★★★ 강양성 — α=1 회복 확인"
elif abs(use_alpha - 1) < 0.10 and use_r2 > 0.95:
    verdict = "★★★★ 양성 — α≈1 근접"
elif abs(use_alpha - 1) < 0.20:
    verdict = "★★★ 조건부 양성 — α=1에 접근 중"
else:
    verdict = "★★ 음성 — α=1 회복 안 됨. 구조적 차이"

print(f"  기준: {label}", flush=True)
print(f"  α = {use_alpha:.4f}, R² = {use_r2:.6f}", flush=True)
print(f"  대각항: α = {alpha_d:.4f}", flush=True)
print(f"\n  종합: {verdict}", flush=True)

elapsed_total = time.time() - t_total

# 결과 파일
with open(OUT_FILE, 'w') as f:
    f.write("=" * 70 + "\n")
    f.write(f"[C-306] GL(3) sym²(11a1) E(σ) 초소Δσ 보충\n")
    f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"경과: {elapsed_total:.1f}s\n")
    f.write(f"N_zeros = {N_zeros}, πN = {pi_N:.4f}\n")
    f.write("=" * 70 + "\n")

    f.write("\n초소Δσ 측정:\n")
    for d in new_offsets:
        r = new_results[d]
        f.write(f"  Δσ={d:.3f}: E_left={r['E_left']:.2f}, E_right={r['E_right']:.2f}, "
                f"E_sym={r['E_sym']:.2f}, E_diag={r['E_diag']:.2f}, "
                f"cross={r['cross_frac']:.3f}, sym={r['sym_ratio']:.6f}\n")

    f.write(f"\n통합 피팅:\n")
    f.write(f"  [전체 10점] α = {alpha_all:.6f}, A/πN = {A_all/pi_N:.6f}, R² = {r2_all:.8f}\n")
    if mask_sm.sum() >= 3:
        f.write(f"  [Δσ≤0.10] α = {alpha_sm:.6f}, A/πN = {A_sm/pi_N:.6f}, R² = {r2_sm:.8f}\n")
    if mask_vsm.sum() >= 3:
        f.write(f"  [Δσ≤0.05] α = {alpha_vsm:.6f}, A/πN = {A_vsm/pi_N:.6f}, R² = {r2_vsm:.8f}\n")
    f.write(f"  [대각항] α = {alpha_d:.6f}, A/πN = {A_d/pi_N:.6f}, R² = {r2_d:.8f}\n")

    f.write(f"\nΔσ·E 일람:\n")
    f.write(f"  {'Δσ':>6} {'E_full':>10} {'E_diag':>10} {'Δσ·E':>10} {'Δσ·E/πN':>10} {'cross%':>8}\n")
    for j in range(len(all_offsets)):
        d = all_offsets[j]
        ef = all_E_sym[j]
        ed = all_E_diag[j]
        cf = (1 - ed / ef) * 100 if ef > 0 else 0
        f.write(f"  {d:6.3f} {ef:10.2f} {ed:10.2f} {d*ef:10.2f} {d*ef/pi_N:10.4f} {cf:8.1f}\n")

    f.write(f"\n종합 판정: {verdict}\n")

print(f"\n  결과: {OUT_FILE}", flush=True)
print(f"  총 경과: {elapsed_total:.1f}s", flush=True)


if __name__ == '__main__':
    pass  # 모듈 수준에서 이미 실행
