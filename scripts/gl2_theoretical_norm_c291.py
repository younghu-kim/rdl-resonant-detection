#!/usr/bin/env python3
"""
C-291: GL(2) 11a1 이론적 밀도 정규화 — Paper 4 d-감쇠 검증
  목적: C-285의 GL(2) ρ=-0.658이 정규화 아티팩트인지 확인

  C-282b (ζ): 이론적 정규화 d_bar = log(t/(2π))/(2π) → ρ=-0.929
  C-285 (GL2): 경험적 정규화 d_bar = 2/(t_{n+1}-t_{n-1}) → ρ=-0.658

  GL(2) 이론적 밀도: d(t) = (1/π) × log(√N × t / (2π))
  이것으로 재측정하면 δ가 유지되는가?
"""

import sys, os, math, time
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 30  # GL(2) T<500

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(1024 * 10**6)
pari.set_real_precision(100)
print("cypari2 OK", flush=True)

CENTER = 1.0
MU_LIST = [0, 1]
TRIM_FRAC = 0.20
N_MAX = 300
T_MAX = 500.0
T_MIN = 2.0

CURVES = [
    {'name': '11a1', 'coeffs': [0, -1, 1, -10, -20], 'N': 11},
    {'name': '37a1', 'coeffs': [0, 1, 1, -23, -50], 'N': 37},
]

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl2_theoretical_norm_c291.txt'
)

out_f = open(RESULT_PATH, 'w')
def log(msg=''):
    print(msg, flush=True)
    out_f.write(msg + '\n')
    out_f.flush()


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeros_gl2(coeffs, t_max):
    pari(f'E_gl2 = ellinit({coeffs})')
    pari(f'Li_gl2 = lfuninit(E_gl2, [0, {t_max + 2}])')
    pari(f'zv_gl2 = lfunzeros(Li_gl2, {t_max})')
    n = int(str(pari('#zv_gl2')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv_gl2[{i}]'))
        if not math.isnan(t) and t > 0.5:
            zeros.append(t)
    return sorted(zeros)


def compute_A_bare(all_zeros, idx, n_max=N_MAX):
    gamma_0 = all_zeros[idx]
    n_total = len(all_zeros)
    S1 = 0.0
    H1 = 0.0
    for k in range(max(0, idx - n_max), min(n_total, idx + n_max + 1)):
        if k == idx:
            continue
        dg = gamma_0 - all_zeros[k]
        if abs(dg) < 1e-15:
            continue
        S1 += 1.0 / dg
        H1 += 1.0 / (dg * dg)
    return S1 ** 2 + 2.0 * H1


def density_zeta(t):
    """ζ(s) 이론적 밀도: d(t) = log(t/(2π)) / (2π)"""
    if t <= 2 * math.pi:
        return 0.01
    return math.log(t / (2 * math.pi)) / (2 * math.pi)


def density_gl2(t, N_cond):
    """GL(2) 이론적 밀도: d(t) = log(√N × t/(2π)) / π"""
    arg = math.sqrt(N_cond) * t / (2 * math.pi)
    if arg <= 1:
        return 0.01
    return math.log(arg) / math.pi


def analyze_curve(curve_info):
    name = curve_info['name']
    coeffs = curve_info['coeffs']
    N_cond = curve_info['N']

    log(f"\n{'='*60}")
    log(f"  L(s, {name}) — N={N_cond}, T=[{T_MIN},{T_MAX}], N_MAX={N_MAX}")
    log(f"{'='*60}")

    # 영점 수집
    t0 = time.time()
    all_zeros = get_zeros_gl2(coeffs, T_MAX)
    zeros_in = [z for z in all_zeros if z >= T_MIN]
    n = len(zeros_in)
    log(f"  {n}개 영점, t∈[{zeros_in[0]:.3f}, {zeros_in[-1]:.3f}] ({time.time()-t0:.1f}s)")

    # A 계산
    t0 = time.time()
    data = []
    for z in zeros_in:
        idx = all_zeros.index(z)
        a = compute_A_bare(all_zeros, idx)
        if a > 0 and math.isfinite(a):
            data.append({'t': z, 'A_bare': a})
    log(f"  A 계산: {len(data)}/{n} 유효 ({time.time()-t0:.1f}s)")

    n_data = len(data)
    if n_data < 40:
        log(f"  데이터 부족 — 스킵")
        return None

    # 내부 영점 + 간격 (edge 2 제외)
    valid_theo = []
    valid_emp = []
    for pos in range(2, n_data - 2):
        d = data[pos]
        t_n = d['t']
        t_prev = data[pos - 1]['t']
        t_next = data[pos + 1]['t']
        gap_r = t_next - t_n
        gap_l = t_n - t_prev
        if gap_r <= 0 or gap_l <= 0:
            continue

        # 이론적 정규화 (GL(2))
        d_bar_theo = density_gl2(t_n, N_cond)
        # 경험적 정규화 (C-285 방식)
        d_bar_emp = 2.0 / (t_next - t_prev)
        # ζ(s) 이론적 (비교용 — GL(2)에 ζ 밀도 적용하면 안 맞지만 참조용)

        valid_theo.append({
            'A_bare': d['A_bare'],
            'gap_min': min(gap_r, gap_l) * d_bar_theo,
        })
        valid_emp.append({
            'A_bare': d['A_bare'],
            'gap_min': min(gap_r, gap_l) * d_bar_emp,
        })

    # trim 20%
    n_inner = len(valid_theo)
    n_trim = int(n_inner * TRIM_FRAC)
    trim_theo = valid_theo[n_trim: n_inner - n_trim]
    trim_emp = valid_emp[n_trim: n_inner - n_trim]
    n_final = len(trim_theo)

    if n_final < 20:
        log(f"  trim 후 부족 — 스킵")
        return None

    # 이론적 정규화 ρ
    A_arr = np.array([d['A_bare'] for d in trim_theo])
    gm_theo = np.array([d['gap_min'] for d in trim_theo])
    gm_emp = np.array([d['gap_min'] for d in trim_emp])

    mask = np.isfinite(A_arr) & np.isfinite(gm_theo) & (A_arr > 0) & (gm_theo > 0)

    rho_theo, p_theo = stats.spearmanr(A_arr[mask], gm_theo[mask])
    rho_emp, p_emp = stats.spearmanr(A_arr[mask], gm_emp[mask])
    se = 1.0 / np.sqrt(np.sum(mask) - 3)

    log(f"\n  n_trim = {np.sum(mask)}")
    log(f"  ρ_S (이론적 정규화): {rho_theo:+.6f}  p={p_theo:.3e}  SE={se:.4f}")
    log(f"  ρ_S (경험적 정규화): {rho_emp:+.6f}  p={p_emp:.3e}  SE={se:.4f}")
    log(f"  C-285 참조 (경험적): ~-0.658")
    log(f"  Δρ (이론-경험):      {rho_theo - rho_emp:+.4f}")

    return {
        'name': name, 'N': N_cond,
        'n_trim': int(np.sum(mask)),
        'rho_theo': rho_theo, 'rho_emp': rho_emp,
        'p_theo': p_theo, 'se': se,
    }


def main():
    t_start = time.time()
    log("=" * 60)
    log("  C-291 — GL(2) 이론적 vs 경험적 밀도 정규화 비교")
    log("=" * 60)
    log(f"  T=[{T_MIN},{T_MAX}], trim={TRIM_FRAC}, N_MAX={N_MAX}")
    log(f"  GL(2) 이론적 밀도: d(t) = log(√N·t/(2π)) / π")
    log()

    results = []
    for curve in CURVES:
        r = analyze_curve(curve)
        if r:
            results.append(r)

    # 요약표
    log(f"\n{'='*60}")
    log("  요약표")
    log(f"{'='*60}")
    log(f"{'곡선':<8} {'N':>4} {'n':>6} {'ρ_theo':>10} {'ρ_emp':>10} {'Δρ':>8}")
    log("-" * 50)
    for r in results:
        log(f"{r['name']:<8} {r['N']:>4d} {r['n_trim']:>6d} "
            f"{r['rho_theo']:>+10.4f} {r['rho_emp']:>+10.4f} "
            f"{r['rho_theo']-r['rho_emp']:>+8.4f}")

    # ζ(s) 비교
    log(f"\n  비교:")
    log(f"    ζ(s) 이론적: ρ = -0.929  (C-282b)")
    log(f"    ζ(s) 경험적: ρ = -0.636  (C-290)")
    for r in results:
        delta_theo = abs(r['rho_theo']) - 0.929
        delta_emp = abs(r['rho_emp']) - 0.636
        log(f"    {r['name']} δ_theo = {abs(r['rho_theo']):.4f} - 0.929 = {delta_theo:+.4f}")
        log(f"    {r['name']} δ_emp  = {abs(r['rho_emp']):.4f} - 0.636 = {delta_emp:+.4f}")

    # 판정
    log(f"\n{'='*60}")
    log("  최종 판정")
    log(f"{'='*60}")
    if results:
        avg_theo = np.mean([r['rho_theo'] for r in results])
        avg_emp = np.mean([r['rho_emp'] for r in results])
        delta_theo = 0.929 - abs(avg_theo)
        delta_emp = 0.636 - abs(avg_emp)

        if abs(delta_theo) < 0.05:
            log(f"  ★★★★★ d-감쇠 아티팩트 확정. GL(2) 이론적 ρ ≈ ζ(s) 이론적.")
            log(f"  Paper 4의 δ≈0.24는 정규화 불일치의 산물.")
        elif delta_theo > 0.05:
            log(f"  ★★★★ d-감쇠 실재 확인. GL(2) 이론적 ρ < ζ(s) 이론적.")
            log(f"  δ_theo = {delta_theo:.4f} (정규화 통일 후에도 유지)")
        else:
            log(f"  ★★★ 예상 외. GL(2) 이론적 > ζ(s) 이론적?")

        log(f"  avg_theo_GL2 = {avg_theo:+.4f}")
        log(f"  avg_emp_GL2  = {avg_emp:+.4f}")
        log(f"  ζ 이론적     = -0.929")
        log(f"  δ_theo       = {delta_theo:+.4f}")

    log(f"\n  소요: {time.time()-t_start:.1f}s")
    out_f.close()


if __name__ == '__main__':
    main()
