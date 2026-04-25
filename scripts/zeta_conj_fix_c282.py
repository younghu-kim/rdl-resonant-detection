#!/usr/bin/env python3
"""
=============================================================================
[사이클 #282] ζ(s) 켤레합 포함 A_L 재측정 — B-50 검증
=============================================================================

목적:
  C-278의 A_L 계산에서 켤레 영점합(Σ 1/(γ₀+γₖ))이 누락되어 있었음.
  C-270에서는 포함 → ρ=-0.80. C-278에서는 미포함 → ρ=-0.44.
  이 차이가 "damping 0.38"이 실재인지 아티팩트인지를 판별.

수정 사항 (C-278 → C-282):
  1. compute_A에 켤레 영점합 추가 (C-270 lines 200-233 로직)
  2. N_MAX=300 제거 → 전체 영점 사용
  3. T=2000 (C-270과 동일)
  4. 60% trim (C-270과 동일)
  5. A_bare_old (C-278 방식, ±300 이웃, 같은편만) 함께 출력하여 직접 비교

성공 기준:
  - ρ(A_L_new, gap_min) > -0.70 → 켤레합 효과 확인 (C-270=-0.80 근접)
  - ρ(A_bare_old, gap_min) ≈ -0.44 → C-278과 일관
  - C-270 결과와 교차 검증

체크리스트:
  [x] CENTER=0.5, mu=[0], N_cond=1
  [x] 켤레합 포함: S1_same + S1_conj (C-270 동일)
  [x] 전체 영점 사용 (N_MAX 제거)
  [x] A_bare_old: C-278 방식 (±300, same-side only)
  [x] A_L_new: 켤레합 포함 primitive
  [x] A_Λ_new: 켤레합 + Gamma 포함 completed
  [x] T=2000, 60% trim
  [x] Spearman (scipy.stats.spearmanr)
  [x] python -u
  [x] NaN/Inf 체크
  [x] 에러 처리 (except: pass 없음)

결과: results/zeta_conj_fix_c282.txt
=============================================================================
"""

import sys
import os
import math
import time

import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 40  # t up to 2000 → 40자리 충분 (C-270에서도 30 사용)

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(512 * 10**6)
    pari.set_real_precision(50)
    print("cypari2 OK (512 MB, precision=50)")
except Exception as e:
    print(f"FATAL: cypari2 로드 실패: {e}")
    sys.exit(1)

# ── 설정 ──────────────────────────────────────────────────────────
T_MAX = 2000.0        # C-270과 동일
CENTER = 0.5
MU_LIST = [0]
N_COND = 1
TRIM_FRAC = 0.20      # 양쪽 20% 제외 → 중앙 60% (C-270과 동일)
N_MAX_OLD = 300       # C-278 방식 비교용

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/zeta_conj_fix_c282.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)


def pf(x):
    """PARI 객체를 float으로 변환"""
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeta_zeros(t_max):
    """PARI lfunzeros로 ζ(s) 영점 수집."""
    print(f"[영점] PARI lfuninit(ζ, [0, {int(t_max)+100}]) ...")
    t0 = time.time()
    pari('L_zeta = lfuncreate(1)')
    pari(f'Li_zeta = lfuninit(L_zeta, [0, {int(t_max) + 100}])')
    pari(f'zv_zeta = lfunzeros(Li_zeta, {t_max})')
    n = int(str(pari('#zv_zeta')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv_zeta[{i}]'))
        if not math.isnan(t) and t > 0.5:
            zeros.append(t)
    zeros = sorted(zeros)
    print(f"  완료: {n}개 원시, {len(zeros)}개 유효 (t>0.5)")
    print(f"  소요: {time.time()-t0:.1f}s")
    return zeros


def compute_A_new(zeros, idx):
    """
    C-270 방식: 같은편 전체 + 켤레 전체 영점합.

    S1_same = Σ_{k≠idx} 1/(γ₀ - γₖ)     (전체 영점)
    S1_conj = Σ_k       1/(γ₀ + γₖ)     (전체 영점, 자기 포함)
    H1_same = Σ_{k≠idx} 1/(γ₀ - γₖ)²
    H1_conj = Σ_k       1/(γ₀ + γₖ)²

    A_L   = (S1_same + S1_conj)² + 2*(H1_same + H1_conj)     [primitive]
    A_Λ   = (S1_total - im_Γ)² + 2*(H1_total + re_Γ')        [completed]
    """
    gamma_0 = zeros[idx]
    n_total = len(zeros)

    # 같은편 영점합 (전체)
    S1_same = 0.0
    H1_same = 0.0
    for k in range(n_total):
        if k == idx:
            continue
        diff = gamma_0 - zeros[k]
        if abs(diff) < 1e-15:
            continue
        S1_same += 1.0 / diff
        H1_same += 1.0 / (diff * diff)

    # 켤레 영점합 (전체, 자기 포함)
    S1_conj = 0.0
    H1_conj = 0.0
    for k in range(n_total):
        s = gamma_0 + zeros[k]
        S1_conj += 1.0 / s
        H1_conj += 1.0 / (s * s)

    S1_total = S1_same + S1_conj
    H1_total = H1_same + H1_conj
    A_L = S1_total ** 2 + 2.0 * H1_total

    # Gamma 보정 (C-270 lines 224-232)
    s_val = mpmath.mpc(CENTER, gamma_0)
    gamma_S = -mpmath.log(mpmath.pi) / 2 + mpmath.digamma(s_val / 2) / 2
    im_gamma_S = float(mpmath.im(gamma_S))
    psi1_val = mpmath.psi(1, s_val / 2)
    re_gamma_H = float(mpmath.re(psi1_val)) / 4.0

    S1_Lambda = S1_total - im_gamma_S
    H1_Lambda = H1_total + re_gamma_H
    A_Lambda = S1_Lambda ** 2 + 2.0 * H1_Lambda

    return {
        'S1_same': S1_same,
        'S1_conj': S1_conj,
        'S1_total': S1_total,
        'H1_same': H1_same,
        'H1_conj': H1_conj,
        'H1_total': H1_total,
        'A_L': A_L,
        'A_Lambda': A_Lambda,
        'im_gamma_S': im_gamma_S,
        're_gamma_H': re_gamma_H,
    }


def compute_A_old(zeros, idx, n_max=N_MAX_OLD):
    """
    C-278 방식: 같은편만, ±n_max 이웃만.
    직접 비교용.
    """
    gamma_0 = zeros[idx]
    n_total = len(zeros)

    S1_bare = 0.0
    H1_bare = 0.0
    for k in range(max(0, idx - n_max), min(n_total, idx + n_max + 1)):
        if k == idx:
            continue
        dg = gamma_0 - zeros[k]
        if abs(dg) < 1e-15:
            continue
        S1_bare += 1.0 / dg
        H1_bare += 1.0 / (dg * dg)

    # B_smooth (C-278과 동일)
    s = mpmath.mpc(CENTER, gamma_0)
    total = mpmath.mpc(0)
    for mu in MU_LIST:
        total += mpmath.digamma((s + mu) / 2) - mpmath.log(mpmath.pi)
    total /= 2
    sm_im = float(mpmath.im(total))

    A_bare = S1_bare ** 2 + 2.0 * H1_bare
    A_L = (S1_bare - sm_im) ** 2 + 2.0 * H1_bare

    return {
        'A_bare_old': A_bare,
        'A_L_old': A_L,
        'S1_bare_old': S1_bare,
        'H1_bare_old': H1_bare,
        'sm_im': sm_im,
    }


def main():
    t_start = time.time()

    print("=" * 80)
    print("  [C-282] ζ(s) 켤레합 포함 A_L 재측정 — B-50 검증")
    print("=" * 80)
    print(f"  T_MAX={T_MAX}, CENTER={CENTER}, TRIM={TRIM_FRAC*2*100:.0f}% 제외")
    print(f"  NEW: 전체 영점 + 켤레합 (C-270 방식)")
    print(f"  OLD: ±{N_MAX_OLD} 이웃, 같은편만 (C-278 방식)")
    print()

    # ── 1. 영점 수집 ────────────────────────────────────────────────
    zeros = get_zeta_zeros(T_MAX)
    n_all = len(zeros)
    print(f"  전체 영점: {n_all}개, t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
    print(f"  소요: {time.time()-t_start:.1f}s")
    print()

    if n_all < 100:
        print("FATAL: 영점 부족 (<100)")
        sys.exit(1)

    # ── 2. A 계산 (NEW + OLD 동시) ──────────────────────────────────
    print(f"[A(γ)] 계산 ({n_all}개 × 전체합)...")
    print(f"  주의: O(n²) 연산 — 시간 소요 예상")
    t0 = time.time()

    data = []
    n_fail = 0

    for i in range(n_all):
        try:
            r_new = compute_A_new(zeros, i)
            r_old = compute_A_old(zeros, i)

            # NaN/Inf 체크
            if (math.isnan(r_new['A_L']) or math.isinf(r_new['A_L']) or
                r_new['A_L'] <= 0):
                n_fail += 1
                continue

            if (math.isnan(r_old['A_bare_old']) or math.isinf(r_old['A_bare_old']) or
                r_old['A_bare_old'] <= 0):
                n_fail += 1
                continue

            entry = {
                't': zeros[i],
                # NEW (C-282: 전체 + 켤레)
                'A_L_new': r_new['A_L'],
                'A_Lambda_new': r_new['A_Lambda'],
                'S1_total': r_new['S1_total'],
                'H1_total': r_new['H1_total'],
                'S1_conj': r_new['S1_conj'],
                'H1_conj': r_new['H1_conj'],
                'im_gamma_S': r_new['im_gamma_S'],
                're_gamma_H': r_new['re_gamma_H'],
                # OLD (C-278: ±300, same-side)
                'A_bare_old': r_old['A_bare_old'],
                'A_L_old': r_old['A_L_old'],
                'S1_bare_old': r_old['S1_bare_old'],
                'H1_bare_old': r_old['H1_bare_old'],
                'sm_im': r_old['sm_im'],
            }
            data.append(entry)

        except Exception as e:
            n_fail += 1
            if n_fail <= 5:
                print(f"  WARNING: i={i}, t={zeros[i]:.3f}: {e}", flush=True)

        if (i + 1) % 200 == 0:
            elapsed = time.time() - t0
            eta = elapsed / (i + 1) * (n_all - i - 1)
            print(f"  {i+1}/{n_all} 완료 ({elapsed:.0f}s, ETA {eta:.0f}s, "
                  f"유효={len(data)}, 실패={n_fail})", flush=True)

    print(f"  A 계산 완료: 유효={len(data)}/{n_all}, 실패={n_fail}, "
          f"소요={time.time()-t0:.1f}s")
    print()

    if len(data) < 100:
        print("FATAL: 유효 데이터 100개 미만")
        sys.exit(1)

    if n_fail > n_all // 2:
        print(f"⚠️ 실패율 {n_fail}/{n_all} — 절반 초과, 신뢰도 낮음")

    # ── 3. 60% trim + gap 계산 ──────────────────────────────────────
    print("[중앙 60% 선택 + gap 계산]")
    N = len(data)
    lo = int(N * TRIM_FRAC)
    hi = int(N * (1.0 - TRIM_FRAC))

    valid = []
    for idx in range(lo, hi):
        if idx <= 0 or idx >= N - 1:
            continue
        d = data[idx]
        t_n = d['t']
        t_prev = data[idx - 1]['t']
        t_next = data[idx + 1]['t']
        gap_r = t_next - t_n
        gap_l = t_n - t_prev
        if gap_r <= 0 or gap_l <= 0:
            continue
        # GUE 정규화: d̄ = log(t/(2π))/(2π)  (C-270과 동일)
        d_bar = np.log(t_n / (2.0 * np.pi)) / (2.0 * np.pi)
        if d_bar <= 0:
            continue
        d['gap_r_gue'] = gap_r * d_bar
        d['gap_l_gue'] = gap_l * d_bar
        d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
        d['d_bar'] = d_bar
        valid.append(d)

    n_valid = len(valid)
    print(f"  내부 영점 (중앙 60%): {n_valid}")
    print(f"  t 범위: [{valid[0]['t']:.3f}, {valid[-1]['t']:.3f}]")
    print()

    if n_valid < 50:
        print("FATAL: 내부 영점 부족 (<50)")
        sys.exit(1)

    # ── 4. 배열 추출 ────────────────────────────────────────────────
    A_L_new_arr     = np.array([d['A_L_new']     for d in valid])
    A_Lambda_new_arr= np.array([d['A_Lambda_new'] for d in valid])
    A_bare_old_arr  = np.array([d['A_bare_old']  for d in valid])
    A_L_old_arr     = np.array([d['A_L_old']     for d in valid])
    gm_arr          = np.array([d['gap_min_gue'] for d in valid])
    gr_arr          = np.array([d['gap_r_gue']   for d in valid])
    t_arr           = np.array([d['t']           for d in valid])

    # NaN 체크
    for name, arr in [('A_L_new', A_L_new_arr), ('A_Λ_new', A_Lambda_new_arr),
                       ('A_bare_old', A_bare_old_arr), ('gap_min', gm_arr)]:
        n_nan = int(np.sum(np.isnan(arr) | np.isinf(arr)))
        if n_nan > 0:
            print(f"  WARNING: {name}에 NaN/Inf {n_nan}개")

    # ── 5. Spearman 상관 ────────────────────────────────────────────
    sig = lambda p: '✅' if p < 0.001 else ('⚠️' if p < 0.01 else '❌')
    GUE_REF = -0.857

    print("=" * 80)
    print(f"  [Spearman 상관] n={n_valid}")
    print("=" * 80)
    print()

    results = {}
    for label, x_arr, y_arr in [
        # OLD 방식 (C-278 재현)
        ("A_bare_old × gap_min  [C-278]", A_bare_old_arr, gm_arr),
        ("A_L_old    × gap_min  [C-278]", A_L_old_arr,    gm_arr),
        # NEW 방식 (C-282 수정)
        ("A_L_new    × gap_min  [C-282]", A_L_new_arr,    gm_arr),
        ("A_Λ_new    × gap_min  [C-282]", A_Lambda_new_arr, gm_arr),
        # gap_right 교차검증
        ("A_L_new    × gap_right[C-282]", A_L_new_arr,    gr_arr),
        ("A_Λ_new    × gap_right[C-282]", A_Lambda_new_arr, gr_arr),
    ]:
        try:
            r, p = stats.spearmanr(x_arr, y_arr)
            results[label] = (r, p)
            damping = abs(GUE_REF) - abs(r)
            print(f"  ρ_S({label}) = {r:+.4f}  p={p:.3e}  {sig(p)}  "
                  f"damping={damping:+.3f}")
        except Exception as e:
            print(f"  WARNING: {label} — {e}")

    print()

    # ── 6. 보조 통계 ────────────────────────────────────────────────
    print("[보조 통계]")

    # 2H₁/A 비율
    H1_total_arr = np.array([d['H1_total'] for d in valid])
    H1_bare_old_arr = np.array([d['H1_bare_old'] for d in valid])

    ratio_new = 2.0 * H1_total_arr / A_L_new_arr
    ratio_old = 2.0 * H1_bare_old_arr / A_bare_old_arr

    print(f"  2H₁/A_L_new  (mean): {np.nanmean(ratio_new):.4f} ± {np.nanstd(ratio_new):.4f}")
    print(f"  2H₁/A_bare_old(mean): {np.nanmean(ratio_old):.4f} ± {np.nanstd(ratio_old):.4f}")
    print(f"  <A_L_new>:     {np.mean(A_L_new_arr):.3f}")
    print(f"  <A_Λ_new>:     {np.mean(A_Lambda_new_arr):.3f}")
    print(f"  <A_bare_old>:  {np.mean(A_bare_old_arr):.3f}")
    print(f"  <A_L_old>:     {np.mean(A_L_old_arr):.3f}")
    print()

    # S1_conj 기여도
    S1_conj_arr = np.array([d['S1_conj'] for d in valid])
    S1_total_arr = np.array([d['S1_total'] for d in valid])
    print(f"  <S1_conj>:     {np.mean(S1_conj_arr):.4f} ± {np.std(S1_conj_arr):.4f}")
    print(f"  <S1_total>:    {np.mean(S1_total_arr):.4f} ± {np.std(S1_total_arr):.4f}")
    print(f"  S1_conj/S1_total: {np.mean(np.abs(S1_conj_arr) / (np.abs(S1_total_arr) + 1e-20)):.4f}")
    print()

    # B_smooth (old) vs im_gamma_S (new Gamma)
    sm_arr = np.array([d['sm_im'] for d in valid])
    im_G_arr = np.array([d['im_gamma_S'] for d in valid])
    rho_bsm_gm, p_bsm_gm = stats.spearmanr(sm_arr, gm_arr)
    print(f"  ρ(B_smooth, gap_min): {rho_bsm_gm:+.4f} (p={p_bsm_gm:.3e})")
    print()

    # ── 7. t-bin 분석 ───────────────────────────────────────────────
    N_BINS = 4
    print(f"[t-bin 분석 ({N_BINS} bins)]")
    print()
    print(f"  {'bin':>4}  {'t̄':>8}  {'n':>5}  "
          f"{'ρ(A_bare_old,gm)':>18}  {'ρ(A_L_old,gm)':>15}  "
          f"{'ρ(A_L_new,gm)':>15}  {'ρ(A_Λ_new,gm)':>15}")
    print("  " + "-" * 100)

    t_quantiles = np.percentile(t_arr, np.linspace(0, 100, N_BINS + 1))

    for b in range(N_BINS):
        if b < N_BINS - 1:
            mask = (t_arr >= t_quantiles[b]) & (t_arr < t_quantiles[b+1])
        else:
            mask = (t_arr >= t_quantiles[b]) & (t_arr <= t_quantiles[b+1])

        n_bin = int(np.sum(mask))
        if n_bin < 10:
            print(f"  {b+1:>4}  n={n_bin} < 10 — 건너뜀")
            continue

        t_mean = float(np.mean(t_arr[mask]))

        r_bo, _ = stats.spearmanr(A_bare_old_arr[mask], gm_arr[mask])
        r_lo, _ = stats.spearmanr(A_L_old_arr[mask],    gm_arr[mask])
        r_ln, _ = stats.spearmanr(A_L_new_arr[mask],    gm_arr[mask])
        r_an, _ = stats.spearmanr(A_Lambda_new_arr[mask], gm_arr[mask])

        print(f"  {b+1:>4}  {t_mean:>8.1f}  {n_bin:>5}  "
              f"{r_bo:>+18.4f}  {r_lo:>+15.4f}  {r_ln:>+15.4f}  {r_an:>+15.4f}")

    print()

    # ── 8. 판정 ─────────────────────────────────────────────────────
    print("=" * 80)
    print("  [판정] — B-50 검증")
    print("=" * 80)
    print()

    r_old_gm = results.get("A_bare_old × gap_min  [C-278]", (np.nan, np.nan))[0]
    r_Lold_gm = results.get("A_L_old    × gap_min  [C-278]", (np.nan, np.nan))[0]
    r_Lnew_gm = results.get("A_L_new    × gap_min  [C-282]", (np.nan, np.nan))[0]
    r_Anew_gm = results.get("A_Λ_new    × gap_min  [C-282]", (np.nan, np.nan))[0]

    print(f"  C-278 A_bare_old gap_min: ρ = {r_old_gm:+.4f}  (기대: ≈ -0.44)")
    print(f"  C-278 A_L_old    gap_min: ρ = {r_Lold_gm:+.4f}  (기대: ≈ -0.44)")
    print(f"  C-282 A_L_new    gap_min: ρ = {r_Lnew_gm:+.4f}  (기대: < -0.70)")
    print(f"  C-282 A_Λ_new    gap_min: ρ = {r_Anew_gm:+.4f}")
    print(f"  C-270 참조:       A_L    gap_min: ρ = -0.800")
    print(f"  C-270 참조:       A_Λ    gap_min: ρ = -0.900")
    print(f"  GUE 참조:                          ρ = {GUE_REF}")
    print()

    # 켤레합 효과
    if not np.isnan(r_Lnew_gm) and not np.isnan(r_Lold_gm):
        delta = abs(r_Lnew_gm) - abs(r_Lold_gm)
        print(f"  켤레합 효과: |ρ_new| - |ρ_old| = {delta:+.4f}")
        if delta > 0.20:
            print(f"  ★★★★ 켤레합 효과 대폭 확인 — damping 0.38은 아티팩트!")
            print(f"         C-278 방법론의 켤레합 누락이 ρ를 50% 약화시킴")
        elif delta > 0.10:
            print(f"  ★★★ 켤레합 효과 중등 — damping 부분적 아티팩트")
        elif delta > 0.05:
            print(f"  ★★ 켤레합 효과 미약 — damping 대부분 실재")
        else:
            print(f"  ★★ 켤레합 효과 없음 — 차이는 T범위/trim 차이에 기인")
        print()

    # C-270 일치 검증
    if not np.isnan(r_Lnew_gm):
        diff_c270 = abs(r_Lnew_gm - (-0.800))
        print(f"  C-270 A_L 대비 차이: |ρ_new - (-0.80)| = {diff_c270:.4f}")
        if diff_c270 < 0.05:
            print(f"  ✅ C-270과 완전 일치 — 방법론 동일 확인")
        elif diff_c270 < 0.10:
            print(f"  ⚠️ C-270과 근사 일치 (trim/T 미세 차이)")
        else:
            print(f"  ❌ C-270과 불일치 — 추가 분석 필요")
        print()

    # ── 9. 결과 파일 저장 ────────────────────────────────────────────
    with open(RESULT_PATH, 'w') as f:
        f.write("# C-282 ζ(s) 켤레합 포함 A_L 재측정 — B-50 검증\n")
        f.write(f"# T_MAX={T_MAX}, CENTER={CENTER}, TRIM={TRIM_FRAC*2*100:.0f}%\n")
        f.write(f"# n_all={n_all}, n_valid={n_valid}\n")
        f.write(f"# NEW: 전체 영점 + 켤레합 (C-270 방식)\n")
        f.write(f"# OLD: ±{N_MAX_OLD} 이웃, 같은편만 (C-278 방식)\n\n")

        f.write("## 핵심 결과\n")
        for label, (r, p) in results.items():
            f.write(f"ρ_S({label}) = {r:+.6f}  p={p:.4e}\n")
        f.write(f"\nGUE_ref = {GUE_REF}\n")
        f.write(f"C-270_ref_A_L = -0.800\n")
        f.write(f"C-270_ref_A_Lambda = -0.900\n\n")

        f.write("## 보조 통계\n")
        f.write(f"2H1/A_L_new_mean  = {np.nanmean(ratio_new):.6f}\n")
        f.write(f"2H1/A_bare_old_mean = {np.nanmean(ratio_old):.6f}\n")
        f.write(f"mean_A_L_new      = {np.mean(A_L_new_arr):.4f}\n")
        f.write(f"mean_A_Lambda_new = {np.mean(A_Lambda_new_arr):.4f}\n")
        f.write(f"mean_A_bare_old   = {np.mean(A_bare_old_arr):.4f}\n")
        f.write(f"mean_A_L_old      = {np.mean(A_L_old_arr):.4f}\n")
        f.write(f"mean_S1_conj      = {np.mean(S1_conj_arr):.6f}\n\n")

        f.write("## 개별 데이터\n")
        f.write("t,A_bare_old,A_L_old,A_L_new,A_Lambda_new,"
                "gap_min_gue,gap_right_gue,S1_conj,H1_conj,sm_im,im_gamma_S\n")
        for d in valid:
            f.write(f"{d['t']:.6f},{d['A_bare_old']:.6f},{d['A_L_old']:.6f},"
                    f"{d['A_L_new']:.6f},{d['A_Lambda_new']:.6f},"
                    f"{d['gap_min_gue']:.6f},{d['gap_r_gue']:.6f},"
                    f"{d['S1_conj']:.6f},{d['H1_conj']:.6f},"
                    f"{d['sm_im']:.6f},{d['im_gamma_S']:.6f}\n")

    print(f"  결과 저장: {RESULT_PATH}")
    print(f"  총 소요: {time.time()-t_start:.1f}s")
    print("=" * 80)
    print("  [C-282 완료]")
    print("=" * 80)


if __name__ == '__main__':
    main()
