#!/usr/bin/env python3
"""
[C-282c] GL(2) 11a1/37a1 T확장 + trim 효과 검증

목적: C-282b에서 ζ(s)의 ρ가 T=600/no-trim(-0.44) → T=2000/trim(-0.93)로
      급변한 것이 GL(2)에서도 동일한지 확인.
      → YES면 통일표 전면 폐기, NO면 ζ(s) 특수 현상.

설계:
  11a1: T_MAX=500 (C-279는 T=100), trim 0%/20%
  37a1: T_MAX=500, trim 0%/20%
  A_bare (±200, same-side) — C-282b에서 최강으로 확인된 공식

체크리스트:
  [x] python -u, mpmath dps=30, Spearman
  [x] C-279와 동일한 A_bare 공식 (비교 가능)
"""

import sys, os, math, time
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 30

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(1000000000)
    pari.set_real_precision(100)
    print("cypari2 OK")
except Exception as e:
    print(f"FATAL: {e}")
    sys.exit(1)

CENTER = 1.0
T_MAX = 500.0
T_MIN = 2.0
N_MAX = 200
TRIM_FRAC = 0.20

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl2_trim_test_c282c.txt'
)

out_lines = []
def log(msg):
    print(msg, flush=True)
    out_lines.append(msg)


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeros_gl2(coeffs, t_max):
    pari(f'E = ellinit({coeffs})')
    pari(f'Li = lfuninit(E, [0, {int(t_max) + 5}])')
    pari(f'zv = lfunzeros(Li, {t_max})')
    n = int(str(pari('#zv')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv[{i}]'))
        if not math.isnan(t) and t > 0.5:
            zeros.append(t)
    return sorted(zeros)


def gamma_smooth_im(gamma_0, mu_list, N_cond):
    s = mpmath.mpc(CENTER, gamma_0)
    total = mpmath.mpc(0)
    for mu in mu_list:
        total += mpmath.digamma((s + mu) / 2) - mpmath.log(mpmath.pi)
    total /= 2
    total += mpmath.log(N_cond) / 2
    return float(mpmath.im(total))


def compute_A(zeros, idx, mu_list, N_cond, n_max=N_MAX):
    g0 = zeros[idx]
    n_z = len(zeros)
    S1 = 0.0
    H1 = 0.0
    lo = max(0, idx - n_max)
    hi = min(n_z, idx + n_max + 1)
    for k in range(lo, hi):
        if k == idx:
            continue
        dg = g0 - zeros[k]
        if abs(dg) < 1e-15:
            continue
        S1 += 1.0 / dg
        H1 += 1.0 / (dg * dg)

    sm_im = gamma_smooth_im(g0, mu_list, N_cond)
    A_bare = S1 ** 2 + 2.0 * H1
    A_L = (S1 - sm_im) ** 2 + 2.0 * H1

    return {'t': g0, 'A_bare': A_bare, 'A_L': A_L, 'S1': S1, 'H1': H1, 'sm_im': sm_im}


def analyze_curve(name, coeffs, mu_list, N_cond):
    log(f"\n{'='*70}")
    log(f"  {name} (N={N_cond}, T_MAX={T_MAX})")
    log(f"{'='*70}")

    t0 = time.time()
    all_zeros = get_zeros_gl2(coeffs, T_MAX)
    zeros = [z for z in all_zeros if z >= T_MIN]
    log(f"  영점: {len(zeros)}개, t ∈ [{zeros[0]:.2f}, {zeros[-1]:.2f}], {time.time()-t0:.1f}s")

    # A 계산
    t0 = time.time()
    data = []
    for i, z in enumerate(zeros):
        idx = all_zeros.index(z)
        r = compute_A(all_zeros, idx, mu_list, N_cond)
        if r['A_bare'] > 0 and r['A_L'] > 0:
            data.append(r)
    log(f"  A 계산: {len(data)}개 유효, {time.time()-t0:.1f}s")

    # gap_min
    valid = []
    for d in data:
        idx = all_zeros.index(d['t'])
        if idx > 0 and idx < len(all_zeros) - 1:
            gap_r = all_zeros[idx + 1] - all_zeros[idx]
            gap_l = all_zeros[idx] - all_zeros[idx - 1]
            d_bar = 2.0 / (all_zeros[idx + 1] - all_zeros[idx - 1])
            d['gap_min'] = min(gap_r, gap_l) * d_bar
            valid.append(d)
    data = valid
    n = len(data)
    log(f"  gap 계산: {n}개")

    # ── 비교: no-trim vs trim ──
    def report_rho(subset, label):
        A_bare = np.array([d['A_bare'] for d in subset])
        A_L = np.array([d['A_L'] for d in subset])
        gm = np.array([d['gap_min'] for d in subset])
        rho_bare, p_bare = stats.spearmanr(A_bare, gm)
        rho_L, p_L = stats.spearmanr(A_L, gm)
        log(f"  [{label}] n={len(subset)}")
        log(f"    ρ_S(A_bare, gap_min) = {rho_bare:.4f} (p={p_bare:.2e})")
        log(f"    ρ_S(A_L, gap_min)    = {rho_L:.4f} (p={p_L:.2e})")
        return rho_bare, rho_L

    # No trim
    rho_bare_full, rho_L_full = report_rho(data, "전체 (no trim)")

    # C-279 범위 재현: T<=100
    data_t100 = [d for d in data if d['t'] <= 100]
    if len(data_t100) >= 10:
        rho_bare_t100, rho_L_t100 = report_rho(data_t100, "T≤100 (C-279 재현)")

    # 20% trim
    lo = int(n * TRIM_FRAC)
    hi = int(n * (1.0 - TRIM_FRAC))
    data_trim = data[lo:hi]
    rho_bare_trim, rho_L_trim = report_rho(data_trim, f"중앙 {1-2*TRIM_FRAC:.0%} trim")

    # t-bin 분석 (4 bins)
    log(f"\n  --- t-bin 분석 (trim) ---")
    n_trim = len(data_trim)
    n_bins = min(4, n_trim // 15)
    if n_bins >= 2:
        bin_size = n_trim // n_bins
        log(f"  {'bin':>3} {'t_mean':>8} {'ρ_bare':>8} {'ρ_L':>8} {'n':>5}")
        for b in range(n_bins):
            bd = data_trim[b*bin_size:(b+1)*bin_size]
            if len(bd) < 8:
                continue
            t_mean = np.mean([d['t'] for d in bd])
            rb, _ = stats.spearmanr([d['A_bare'] for d in bd], [d['gap_min'] for d in bd])
            rl, _ = stats.spearmanr([d['A_L'] for d in bd], [d['gap_min'] for d in bd])
            log(f"  {b+1:>3} {t_mean:>8.1f} {rb:>8.3f} {rl:>8.3f} {len(bd):>5}")

    # 판정
    log(f"\n  --- {name} 판정 ---")
    delta = abs(rho_bare_trim) - abs(rho_bare_full)
    log(f"  Trim 효과: Δ|ρ_bare| = {delta:+.4f}")
    if delta > 0.10:
        log(f"  → trim 효과 유의 (>0.10). C-279 결과는 과소추정.")
    else:
        log(f"  → trim 효과 미미 (<0.10). C-279 결과 유효.")

    return {
        'name': name,
        'n_full': len(data),
        'rho_bare_full': rho_bare_full,
        'rho_bare_trim': rho_bare_trim,
        'rho_L_full': rho_L_full,
        'rho_L_trim': rho_L_trim,
    }


def main():
    log("=" * 70)
    log("  C-282c — GL(2) T확장 + trim 효과 검증")
    log("=" * 70)
    log(f"  날짜: {time.strftime('%Y-%m-%d %H:%M')}")
    log(f"  T_MAX={T_MAX}, TRIM={TRIM_FRAC}, N_MAX={N_MAX}")

    t_start = time.time()

    curves = [
        {'name': '11a1', 'coeffs': '[0,-1,1,-10,-20]', 'mu': [0, 1], 'N_cond': 11, 'label': 'GL(2) 11a1'},
        {'name': '37a1', 'coeffs': '[0,1,1,-23,-50]', 'mu': [0, 1], 'N_cond': 37, 'label': 'GL(2) 37a1'},
    ]

    results = []
    for c in curves:
        r = analyze_curve(c['label'], c['coeffs'], c['mu'], c['N_cond'])
        if r:
            results.append(r)

    # 최종 요약
    log(f"\n{'='*70}")
    log(f"  최종 요약")
    log(f"{'='*70}")
    for r in results:
        log(f"  {r['name']}:")
        log(f"    전체: ρ_bare={r['rho_bare_full']:.4f}, ρ_L={r['rho_L_full']:.4f}")
        log(f"    trim: ρ_bare={r['rho_bare_trim']:.4f}, ρ_L={r['rho_L_trim']:.4f}")
        delta = abs(r['rho_bare_trim']) - abs(r['rho_bare_full'])
        log(f"    trim 효과: {delta:+.4f}")

    log(f"\n  C-279 참조: 11a1 ρ_bare=-0.423, 37a1 ρ_bare=-0.525")
    log(f"  C-282b ζ(s): no-trim ρ=-0.44, trim ρ=-0.93")

    all_trim_effects = [abs(r['rho_bare_trim']) - abs(r['rho_bare_full']) for r in results]
    avg_effect = np.mean(all_trim_effects)
    log(f"\n  평균 trim 효과: {avg_effect:+.4f}")

    if avg_effect > 0.10:
        log(f"\n  ✅ trim 효과 확인. 통일표 ρ 값은 과소추정. 전면 재측정 필요.")
    else:
        log(f"\n  ❌ GL(2)에서는 trim 효과 미미. ζ(s) 특수 현상 가능.")

    log(f"\n총 소요: {time.time()-t_start:.1f}s")

    os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)
    with open(RESULT_PATH, 'w') as f:
        f.write('\n'.join(out_lines))
    log(f"결과 저장: {RESULT_PATH}")


if __name__ == '__main__':
    main()
