#!/usr/bin/env python3
"""
C-290c: T-대역 A-gap 안정성 — 이론적 밀도 정규화 버전
C-290과 동일한 실험이지만 gap 정규화에 이론적 밀도 사용:
  d_bar = log(t/(2π)) / (2π)  (C-282b와 동일)

이전 C-290은 경험적 밀도 d_bar = 2/(t_{n+1}-t_{n-1}) 사용 → ρ 압축.
이 버전에서 T-대역 안정성이 유지되는지 확인.
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

T_MIN = 10.0
T_MAX = 2000.0
CENTER = 0.5
N_MAX = 300  # C-282b와 동일
TRIM_FRAC = 0.20

BANDS = [(10, 500), (500, 1000), (1000, 1500), (1500, 2000)]

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/zeta_tband_agap_c290c.txt'
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


def get_zeta_zeros(t_max):
    log(f"[영점] T_MAX={t_max} ...")
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
    log(f"  {len(zeros)}개 영점 ({time.time()-t0:.1f}s)")
    return zeros


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


def theoretical_density(t):
    """이론적 영점 밀도: d(t) = log(t/(2π)) / (2π)"""
    if t <= 2 * math.pi:
        return 0.01  # 아주 낮은 t에서 안전값
    return math.log(t / (2 * math.pi)) / (2 * math.pi)


def analyze_band(data_in_band, label, t_lo, t_hi):
    n = len(data_in_band)
    if n < 20:
        log(f"  [{label}] n={n} < 20 — 스킵")
        return None

    # 내부 영점 + 간격 (edge 2 제외)
    valid = []
    for pos in range(2, n - 2):
        d = data_in_band[pos]
        t_n = d['t']
        t_prev = data_in_band[pos - 1]['t']
        t_next = data_in_band[pos + 1]['t']
        gap_r = t_next - t_n
        gap_l = t_n - t_prev
        if gap_r <= 0 or gap_l <= 0:
            continue
        # 이론적 밀도 정규화 (C-282b와 동일)
        d_bar = theoretical_density(t_n)
        if d_bar <= 0:
            continue
        d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
        d['gap_r_gue'] = gap_r * d_bar
        valid.append(d)

    n_inner = len(valid)
    if n_inner < 20:
        return None

    # trim 20%
    n_trim_each = int(n_inner * TRIM_FRAC)
    trimmed = valid[n_trim_each: n_inner - n_trim_each]
    n_trim = len(trimmed)
    if n_trim < 10:
        return None

    A_arr = np.array([d['A_bare'] for d in trimmed])
    gm_arr = np.array([d['gap_min_gue'] for d in trimmed])
    mask = np.isfinite(A_arr) & np.isfinite(gm_arr) & (A_arr > 0)
    A_arr = A_arr[mask]
    gm_arr = gm_arr[mask]
    n_final = len(A_arr)
    if n_final < 10:
        return None

    rho_S, p_val = stats.spearmanr(A_arr, gm_arr)
    se = 1.0 / np.sqrt(n_final - 3) if n_final > 3 else float('nan')
    density = n / (t_hi - t_lo)

    sig = '✅' if p_val < 0.001 else '⚠️'
    log(f"\n  [{label}] T=[{t_lo},{t_hi}]")
    log(f"    n={n}, n_trim={n_final}, density={density:.4f}")
    log(f"    ρ_S(A_bare, gap_min) = {rho_S:+.6f}  p={p_val:.3e}  SE={se:.4f}  {sig}")

    return {
        'band': label, 't_lo': t_lo, 't_hi': t_hi,
        'n_zeros': n, 'n_trim': n_final,
        'density': density,
        'rho_S': rho_S, 'p_val': p_val, 'se': se,
    }


def main():
    t_start = time.time()
    log("=" * 70)
    log("  C-290c — T-대역 A-gap (이론적 밀도 정규화, N_MAX=300)")
    log("=" * 70)
    log(f"  T=[{T_MIN},{T_MAX}], trim={TRIM_FRAC}, N_MAX={N_MAX}")
    log(f"  gap 정규화: d_bar = log(t/(2π))/(2π)  [이론적]")
    log()

    all_zeros = get_zeta_zeros(T_MAX)
    zeros_in = [z for z in all_zeros if z >= T_MIN]
    n = len(zeros_in)
    log(f"  T≥{T_MIN}: {n}개")

    # 인덱스 매핑
    idx_map = {}
    for z in zeros_in:
        idx_map[z] = all_zeros.index(z)

    # A 계산
    log(f"\n[A(γ)] 계산 ({n}개)...")
    t0 = time.time()
    all_data = []
    for i, z in enumerate(zeros_in):
        a = compute_A_bare(all_zeros, idx_map[z])
        if not (math.isnan(a) or math.isinf(a) or a <= 0):
            all_data.append({'t': z, 'A_bare': a})
        if (i + 1) % 200 == 0:
            log(f"  ...{i+1}/{n} ({time.time()-t0:.1f}s)")
    log(f"  완료: {len(all_data)}/{n} ({time.time()-t0:.1f}s)")

    # 대역별 분석
    log(f"\n{'='*70}")
    log("  [대역별] 이론적 밀도 정규화")
    log(f"{'='*70}")

    band_results = []
    for t_lo, t_hi in BANDS:
        bd = [d for d in all_data if t_lo <= d['t'] < t_hi]
        r = analyze_band(bd, f"T[{t_lo},{t_hi}]", t_lo, t_hi)
        if r:
            band_results.append(r)

    # 전체
    log(f"\n{'='*70}")
    log("  [전체] 이론적 밀도 정규화")
    log(f"{'='*70}")
    r_all = analyze_band(all_data, "전체", T_MIN, T_MAX)
    if r_all:
        band_results.append(r_all)

    # 요약표
    log(f"\n{'='*70}")
    log("  요약표 (이론적 밀도 정규화)")
    log(f"{'='*70}")
    log(f"{'대역':<16} {'n_trim':>7} {'density':>9} {'ρ_S':>10} {'SE':>8} {'p':>12}")
    log("-" * 70)
    for r in band_results:
        log(f"{r['band']:<16} {r['n_trim']:>7d} {r['density']:>9.4f} "
            f"{r['rho_S']:>+10.6f} {r['se']:>8.4f} {r['p_val']:>12.3e}")

    # 대역간 상관
    band_only = [r for r in band_results if r['band'] != '전체']
    if len(band_only) >= 3:
        rho_arr = np.array([r['rho_S'] for r in band_only])
        dens_arr = np.array([r['density'] for r in band_only])
        rho_meta, p_meta = stats.spearmanr(dens_arr, rho_arr)
        rho_range = max(rho_arr) - min(rho_arr)
        log(f"\n  ρ(density, ρ_S) = {rho_meta:+.4f}  p={p_meta:.4f}")
        log(f"  ρ 범위: {rho_range:.4f} (max={max(rho_arr):+.4f}, min={min(rho_arr):+.4f})")

        # C-282b 전체 비교
        log(f"\n  C-282b 비교: ρ_전체(이론적)={r_all['rho_S']:+.4f} vs C-282b=-0.929")

    # 판정
    log(f"\n{'='*70}")
    log("  최종 판정")
    log(f"{'='*70}")

    if len(band_only) >= 3:
        rho_arr = np.array([r['rho_S'] for r in band_only])
        rho_range = max(rho_arr) - min(rho_arr)

        if rho_range < 0.05:
            verdict = "★★★★ 양성 — T-대역 무관. 밀도 메커니즘 기각."
        elif rho_range < 0.10:
            verdict = "★★★ 중립 — 약한 대역 의존성."
        else:
            verdict = "★★★ 중립 — 대역 의존성 유의."

        log(f"  {verdict}")
        log(f"  ρ 범위: {rho_range:.4f}")

    log(f"\n  소요: {time.time()-t_start:.1f}s")
    out_f.close()


if __name__ == '__main__':
    main()
