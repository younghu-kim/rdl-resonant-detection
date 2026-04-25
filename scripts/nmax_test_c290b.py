#!/usr/bin/env python3
"""
C-290b: N_MAX 민감도 테스트 — N_MAX=200 vs 300이 ρ에 미치는 영향 확인
빠른 진단용. 전체 ζ(s) 영점에서 A_bare를 N_MAX=200, 300, 500으로 계산 후 ρ 비교.
"""

import sys, os, math, time
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 80

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(1024 * 10**6)
pari.set_real_precision(80)

T_MAX = 2000.0
T_MIN = 10.0
CENTER = 0.5
TRIM_FRAC = 0.20
NMAX_LIST = [100, 200, 300, 500, 1000]


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeta_zeros(t_max):
    print(f"[영점] T_MAX={t_max} ...", flush=True)
    t0 = time.time()
    pari('L_zeta = lfuncreate(1)')
    pari(f'Li_zeta = lfuninit(L_zeta, [0, {int(t_max) + 5}])')
    pari(f'zv_zeta = lfunzeros(Li_zeta, {t_max})')
    n = int(str(pari('#zv_zeta')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv_zeta[{i}]'))
        if not math.isnan(t) and t > 0.5:
            zeros.append(t)
    zeros = sorted(zeros)
    print(f"  {len(zeros)}개 영점 ({time.time()-t0:.1f}s)", flush=True)
    return zeros


def compute_A_bare(all_zeros, idx, n_max):
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


def main():
    t0 = time.time()
    all_zeros = get_zeta_zeros(T_MAX)
    zeros_in = [z for z in all_zeros if z >= T_MIN]
    n = len(zeros_in)
    print(f"T≥{T_MIN}: {n}개", flush=True)

    # 인덱스 매핑
    idx_map = {}
    for z in zeros_in:
        idx_map[z] = all_zeros.index(z)

    # 각 N_MAX에서 A_bare 계산
    results = {nm: [] for nm in NMAX_LIST}

    for nm in NMAX_LIST:
        print(f"\n[N_MAX={nm}] A 계산 중...", flush=True)
        t1 = time.time()
        for z in zeros_in:
            a = compute_A_bare(all_zeros, idx_map[z], nm)
            results[nm].append(a)
        print(f"  완료 ({time.time()-t1:.1f}s)", flush=True)

    # 간격 계산
    gaps_min = []
    valid_mask = []
    for i in range(2, n - 2):
        t_n = zeros_in[i]
        t_prev = zeros_in[i - 1]
        t_next = zeros_in[i + 1]
        gap_r = t_next - t_n
        gap_l = t_n - t_prev
        if gap_r > 0 and gap_l > 0:
            d_bar = 2.0 / (t_next - t_prev)
            gaps_min.append(min(gap_r, gap_l) * d_bar)
            valid_mask.append(i)
        else:
            gaps_min.append(float('nan'))
            valid_mask.append(i)

    # trim
    n_inner = len(valid_mask)
    n_trim = int(n_inner * TRIM_FRAC)
    trimmed_indices = valid_mask[n_trim: n_inner - n_trim]
    trimmed_gaps = gaps_min[n_trim: n_inner - n_trim]

    gm = np.array(trimmed_gaps)
    good = np.isfinite(gm) & (gm > 0)

    print(f"\n{'='*60}")
    print(f"N_MAX 민감도 테스트 (n_trim={np.sum(good)}, T=[{T_MIN},{T_MAX}], trim={TRIM_FRAC})")
    print(f"{'='*60}")
    print(f"{'N_MAX':>6} {'ρ_S':>10} {'p':>14} {'mean_A':>12}")
    print("-" * 50)

    for nm in NMAX_LIST:
        A_arr = np.array([results[nm][i] for i in trimmed_indices])
        mask = good & np.isfinite(A_arr) & (A_arr > 0)
        rho, p = stats.spearmanr(A_arr[mask], gm[mask])
        print(f"{nm:>6d} {rho:>+10.6f} {p:>14.3e} {np.mean(A_arr[mask]):>12.4f}")

    print(f"\n총 소요: {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
