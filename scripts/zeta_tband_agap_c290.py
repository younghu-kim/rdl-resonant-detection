#!/usr/bin/env python3
"""
[사이클 #290] ζ(s) T-대역별 A-gap 안정성 — B-54 밀도 메커니즘 검증
  목표: d-의존 감쇠의 원인이 영점 밀도인지 산술 구조인지 판별

  설계:
    - T=[10, 2000], dps=80
    - 4개 T-대역: [10,500], [500,1000], [1000,1500], [1500,2000]
    - 각 대역 내에서: trim 20%, ρ_S(A_bare, gap_min_GUE) 독립 측정
    - A 계산은 전체 영점 사용 (local ±200)
    - 대역별 밀도 = n_zeros / ΔT 측정
    - ρ vs density 상관 확인

  가설:
    - H0: ρ는 T-대역에 무관 (산술 메커니즘, 밀도와 무관)
    - H1: ρ가 높은 T (높은 밀도) 대역에서 약해짐 (밀도 메커니즘)

  이론적 배경:
    ζ(s)의 영점 밀도 ~ (1/2π)log(T/2π). T=500 → ~0.72/unit, T=2000 → ~0.90/unit
    GL(2)는 같은 T에서 밀도 ~2배. 밀도가 원인이면 ζ(s) 자체 내에서도 감쇠 경사 관찰 가능.

  체크리스트:
    [x] dps=80 (t>100)
    [x] A_bare = S1_bare^2 + 2*H1_bare
    [x] gap_min_GUE = min(gap_r, gap_l) * d_bar
    [x] trim 20% 각 대역 독립 적용
    [x] python -u
    [x] N_MAX=200 (local)
"""

import sys, os, math, time
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 80  # t 최대 2000 → dps=80 필수

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(1024 * 10**6)  # 1 GB
    pari.set_real_precision(80)
    print("cypari2 OK (1 GB, precision=80)")
except Exception as e:
    print(f"FATAL: {e}")
    sys.exit(1)

# ── 설정 ──────────────────────────────────────────────────────────
T_MIN = 10.0
T_MAX = 2000.0
CENTER = 0.5
MU_LIST = [0]
N_COND = 1
N_MAX = 200
TRIM_FRAC = 0.20  # 양쪽 20% (중앙 60%)

BANDS = [
    (10, 500),
    (500, 1000),
    (1000, 1500),
    (1500, 2000),
]

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/zeta_tband_agap_c290.txt'
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
    """PARI lfunzeros로 ζ(s) 영점 수집"""
    log(f"[영점] PARI lfuninit(ζ, [0, {t_max+5}]) ...")
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
    log(f"  완료: {n}개 원시, {len(zeros)}개 유효 (t>0.5)")
    log(f"  t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
    log(f"  소요: {time.time()-t0:.1f}s")
    return zeros


def gamma_smooth_im(gamma_0):
    """ζ(s) Γ smooth part 허수부"""
    s = mpmath.mpc(CENTER, gamma_0)
    total = mpmath.mpc(0)
    for mu in MU_LIST:
        total += mpmath.digamma((s + mu) / 2) - mpmath.log(mpmath.pi)
    total /= 2
    return float(mpmath.im(total))


def compute_A(all_zeros, idx, n_max=N_MAX):
    """영점 idx에서 A_bare, A_L 계산 (local ±n_max)"""
    gamma_0 = all_zeros[idx]
    n_total = len(all_zeros)
    S1_bare = 0.0
    H1_bare = 0.0

    for k in range(max(0, idx - n_max), min(n_total, idx + n_max + 1)):
        if k == idx:
            continue
        dg = gamma_0 - all_zeros[k]
        if abs(dg) < 1e-15:
            continue
        S1_bare += 1.0 / dg
        H1_bare += 1.0 / (dg * dg)

    sm_im = gamma_smooth_im(gamma_0)
    A_bare = S1_bare ** 2 + 2.0 * H1_bare
    S1_L = S1_bare - sm_im
    A_L = S1_L ** 2 + 2.0 * H1_bare

    return {'A_bare': A_bare, 'A_L': A_L, 't': gamma_0}


def analyze_band(data_in_band, band_label, t_lo, t_hi):
    """대역 내 A-gap 상관 분석 (trim 20% 적용)"""
    n = len(data_in_band)
    if n < 20:
        log(f"\n  [{band_label}] n={n} < 20 — 스킵")
        return None

    # 내부 영점 + 간격 계산 (edge 2개씩 제외)
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
        d_bar = 2.0 / (t_next - t_prev)
        d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
        valid.append(d)

    n_inner = len(valid)
    if n_inner < 20:
        log(f"\n  [{band_label}] 내부 영점 {n_inner} < 20 — 스킵")
        return None

    # trim 20%
    n_trim_each = int(n_inner * TRIM_FRAC)
    trimmed = valid[n_trim_each: n_inner - n_trim_each]
    n_trim = len(trimmed)
    if n_trim < 10:
        log(f"\n  [{band_label}] trim 후 {n_trim} < 10 — 스킵")
        return None

    A_arr = np.array([d['A_bare'] for d in trimmed])
    gm_arr = np.array([d['gap_min_gue'] for d in trimmed])

    # NaN 필터
    mask = np.isfinite(A_arr) & np.isfinite(gm_arr) & (A_arr > 0)
    A_arr = A_arr[mask]
    gm_arr = gm_arr[mask]
    n_final = len(A_arr)

    if n_final < 10:
        log(f"\n  [{band_label}] NaN 필터 후 {n_final} < 10 — 스킵")
        return None

    rho_S, p_val = stats.spearmanr(A_arr, gm_arr)
    se = 1.0 / np.sqrt(n_final - 3) if n_final > 3 else float('nan')

    # 밀도 계산
    density = n / (t_hi - t_lo)
    # 이론적 밀도: (1/2π) log(T_mid/2π)
    t_mid = (t_lo + t_hi) / 2
    density_theory = math.log(t_mid / (2 * math.pi)) / (2 * math.pi) if t_mid > 2*math.pi else 0

    result = {
        'band': band_label,
        't_lo': t_lo, 't_hi': t_hi,
        'n_zeros': n,
        'n_trim': n_final,
        'density': density,
        'density_theory': density_theory,
        'rho_S': rho_S,
        'p_val': p_val,
        'se': se,
    }

    sig = '✅' if p_val < 0.001 else ('⚠️ ' if p_val < 0.01 else '❌')
    log(f"\n  [{band_label}] T=[{t_lo},{t_hi}]")
    log(f"    n_zeros={n}, n_trim={n_final}, density={density:.4f} (/unit)")
    log(f"    density_theory={density_theory:.4f}")
    log(f"    ρ_S(A_bare, gap_min) = {rho_S:+.6f}  p={p_val:.3e}  SE={se:.4f}  {sig}")

    return result


def main():
    t_start = time.time()

    log("=" * 70)
    log("  사이클 #290 — ζ(s) T-대역별 A-gap 안정성 (B-54 밀도 메커니즘)")
    log("=" * 70)
    log(f"  T=[{T_MIN},{T_MAX}], trim={TRIM_FRAC*100:.0f}%, N_MAX={N_MAX}")
    log(f"  대역: {BANDS}")
    log()

    # ── 1. 영점 수집 ────────────────────────────────────────────────
    all_zeros = get_zeta_zeros(T_MAX)
    zeros_in_range = [z for z in all_zeros if z >= T_MIN]
    log(f"  T≥{T_MIN}: {len(zeros_in_range)}개")
    dt_vals = [zeros_in_range[i+1] - zeros_in_range[i] for i in range(len(zeros_in_range)-1)]
    log(f"  dt_min={min(dt_vals):.6f}  dt_mean={np.mean(dt_vals):.4f}")

    # ── 2. A 계산 (전체) ────────────────────────────────────────────
    log(f"\n[A(γ)] 계산 ({len(zeros_in_range)}개)...")
    t0 = time.time()
    all_data = []
    for i, z in enumerate(zeros_in_range):
        all_idx = all_zeros.index(z)
        r = compute_A(all_zeros, all_idx)
        if not (math.isnan(r['A_bare']) or math.isinf(r['A_bare']) or r['A_bare'] <= 0):
            all_data.append(r)
        if (i + 1) % 200 == 0:
            log(f"  ...{i+1}/{len(zeros_in_range)} ({time.time()-t0:.1f}s)")

    log(f"  완료: {len(all_data)}/{len(zeros_in_range)} 유효 ({time.time()-t0:.1f}s)")

    # ── 3. 대역별 분석 ──────────────────────────────────────────────
    log(f"\n{'='*70}")
    log("  [대역별 분석] trim 20%")
    log(f"{'='*70}")

    band_results = []
    for t_lo, t_hi in BANDS:
        band_data = [d for d in all_data if t_lo <= d['t'] < t_hi]
        label = f"T[{t_lo},{t_hi}]"
        r = analyze_band(band_data, label, t_lo, t_hi)
        if r is not None:
            band_results.append(r)

    # ── 4. 전체 (비교용) ────────────────────────────────────────────
    log(f"\n{'='*70}")
    log("  [전체 범위] trim 20% (C-282b 비교용)")
    log(f"{'='*70}")
    r_all = analyze_band(all_data, "전체", T_MIN, T_MAX)
    if r_all:
        band_results.append(r_all)

    # ── 5. 요약표 ───────────────────────────────────────────────────
    log(f"\n{'='*70}")
    log("  T-대역 요약표")
    log(f"{'='*70}")
    log(f"{'대역':<16} {'n_trim':>7} {'density':>9} {'ρ_S':>10} {'SE':>8} {'p':>12}")
    log("-" * 70)
    for r in band_results:
        log(f"{r['band']:<16} {r['n_trim']:>7d} {r['density']:>9.4f} "
            f"{r['rho_S']:>+10.6f} {r['se']:>8.4f} {r['p_val']:>12.3e}")

    # ── 6. ρ vs density 상관 (대역 4개) ─────────────────────────────
    band_only = [r for r in band_results if r['band'] != '전체']
    if len(band_only) >= 3:
        rho_arr = np.array([r['rho_S'] for r in band_only])
        dens_arr = np.array([r['density'] for r in band_only])
        rho_meta, p_meta = stats.spearmanr(dens_arr, rho_arr)
        log(f"\n  ρ(density, ρ_S) across bands = {rho_meta:+.4f}  p={p_meta:.4f}")
        log(f"  해석: {'밀도 증가 → ρ 약화 (밀도 메커니즘 지지)' if rho_meta > 0.3 else '밀도-ρ 무상관 (산술 메커니즘 지지)' if abs(rho_meta) < 0.3 else '밀도 증가 → ρ 강화 (예상 외)'}")

        # ρ 범위 (최대-최소)
        rho_range = max(rho_arr) - min(rho_arr)
        log(f"  ρ 범위: {rho_range:.4f} (max={max(rho_arr):+.4f}, min={min(rho_arr):+.4f})")
        log(f"  ρ 범위/전체|ρ|: {rho_range/abs(r_all['rho_S'])*100:.1f}%" if r_all else "")
    else:
        log("\n  대역 수 부족 — ρ-density 상관 불가")

    # ── 7. 판정 ─────────────────────────────────────────────────────
    log(f"\n{'='*70}")
    log("  최종 판정")
    log(f"{'='*70}")

    if len(band_only) >= 3:
        rho_arr = np.array([r['rho_S'] for r in band_only])
        rho_range = max(rho_arr) - min(rho_arr)

        if rho_range < 0.05 and abs(rho_meta) < 0.5:
            verdict = "★★★★ 양성 — ρ T-대역 무관. 밀도 메커니즘 기각. d-감쇠는 산술 구조에 기인."
        elif rho_meta > 0.5 and rho_range > 0.05:
            verdict = "★★★★ 양성 — ρ가 높은 T에서 약화. 밀도 메커니즘 지지. d-감쇠의 부분적 원인."
        elif rho_range > 0.10:
            verdict = "★★★ 중립 — ρ 변동 크나 방향성 불확실. 추가 분석 필요."
        else:
            verdict = "★★★ 중립 — 약한 경향성. 결정적이지 않음."

        log(f"  {verdict}")
        log(f"  근거: ρ_range={rho_range:.4f}, ρ(density,ρ_S)={rho_meta:+.4f}")
    else:
        log("  판정 불가 — 대역 부족")

    elapsed = time.time() - t_start
    log(f"\n  총 소요: {elapsed:.1f}s ({elapsed/60:.1f}분)")

    out_f.close()
    print(f"\n결과 저장: {RESULT_PATH}")


if __name__ == '__main__':
    main()
