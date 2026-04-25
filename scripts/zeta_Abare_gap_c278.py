#!/usr/bin/env python3
"""
[사이클 #278] ζ(s) A_bare gap_min Spearman ρ_S 통합 재측정
  목표: 비교표 방법론 통일 — d=1 ζ(s)에서 GL(4)~GL(6)와 동일한 방법론 적용

  설계:
    - T=[10,600], n≈198 내부 영점
    - A_bare = S₁² + 2H₁  (smooth 보정 없음)
    - A_L = (S₁ - B_smooth)² + 2H₁  (Γ smooth 보정, GL(d)와 동일 공식)
    - B_smooth = Im(Γ_smooth) = Im[(digamma(s/2) - log(π)) / 2] at s=0.5+iγ
                              = Im(digamma(0.25 + iγ/2)) / 2
    - gap_min_GUE, gap_right_GUE: 국소 밀도 정규화 (GL(d)와 동일)
    - 4 Spearman 상관:
        1. ρ_S(A_bare, gap_min_GUE)  ← 핵심: GUE 직접 비교
        2. ρ_S(A_bare, gap_right_GUE) ← 기존 A_Λ gap_right 비교용
        3. ρ_S(A_L, gap_min_GUE)      ← d≥3 비교표 통일
        4. ρ_S(A_L, gap_right_GUE)    ← 교차 검증

  기준 참고:
    - gl6_A_gap_c275.py 구조 기반
    - 동일 방법론: same-side zeros only (GL(d)와 일치)
    - CENTER = 0.5, mu=[0], N_cond=1

  체크리스트:
    [x] CENTER=0.5, mu=[0], N_cond=1
    [x] A_bare = S1_bare^2 + 2*H1_bare (smooth 미보정)
    [x] A_L = (S1_bare - sm_im)^2 + 2*H1 (sm_im = GL(d)와 동일 공식)
    [x] gap_min_GUE = min(gap_r, gap_l) * d_bar (국소밀도)
    [x] Spearman (scipy.stats.spearmanr)
    [x] python -u
    [x] 에러 처리
    [x] NaN 필터
"""

import sys
import os
import math
import time

import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 40  # t up to 600 → 40자리 충분

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(512 * 10**6)  # 512 MB (ζ(s)는 작은 conductor)
    pari.set_real_precision(50)
    print("cypari2 OK (512 MB, precision=50)")
except Exception as e:
    print(f"FATAL: cypari2 로드 실패: {e}")
    sys.exit(1)

# ── 설정 ──────────────────────────────────────────────────────────
T_MIN = 10.0
T_MAX = 600.0
CENTER = 0.5     # ζ(s) critical line: Re(s)=1/2
MU_LIST = [0]    # Γ_R(s) = π^{-s/2} Γ(s/2) → mu=0
N_COND = 1       # ζ(s) conductor=1
N_MAX = 300      # 이웃 영점 최대 탐색 범위 (GL(d)에서 200; ζ는 더 풍부)

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/zeta_Abare_gap_c278.txt'
)


def pf(x):
    """PARI 객체를 float으로 변환"""
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeta_zeros(t_max):
    """PARI lfunzeros로 ζ(s) 영점 수집 (t ∈ (0, t_max])."""
    print(f"[영점] PARI lfuninit(ζ, [0, {t_max+5}]) ...")
    t0 = time.time()
    try:
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
        print(f"  완료: {n}개 원시, {len(zeros)}개 유효 (t>0.5)")
        print(f"  소요: {time.time()-t0:.1f}s")
        return zeros
    except Exception as e:
        print(f"FATAL: 영점 수집 실패: {e}")
        sys.exit(1)


def gamma_smooth_im(gamma_0):
    """
    ζ(s)의 Γ smooth part 허수부.

    Γ_∞'/Γ_∞(s) = (1/2)[ψ(s/2) - log(π)]  at s=CENTER+iγ
    Im[smooth] = Im[(ψ(0.25 + iγ/2) - log(π)) / 2]
               = Im(ψ(0.25 + iγ/2)) / 2   (log(π)는 실수)

    GL(d) 공식과 일치: gamma_smooth(γ, mu=[0], N=1)
      = (ψ((s+0)/2) - log(π))/2 + log(1)/2
      = (ψ(s/2) - log(π))/2
    """
    s = mpmath.mpc(CENTER, gamma_0)
    # GL(d)와 동일 공식: (digamma((s+mu)/2) - log(π))/2 summed over mu
    total = mpmath.mpc(0)
    for mu in MU_LIST:
        total += mpmath.digamma((s + mu) / 2) - mpmath.log(mpmath.pi)
    total /= 2
    # log(N_cond)/2 = 0 (N=1)
    return float(mpmath.im(total))


def compute_A(zeros, all_zeros, idx, n_max=N_MAX):
    """
    영점 idx에서 A_bare, A_L 계산.

    - A_bare = S1_bare^2 + 2*H1_bare  (smooth 보정 없음)
    - A_L    = (S1_bare - sm_im)^2 + 2*H1_bare  (Γ smooth 보정)

    same-side zeros only (GL(d)와 일치).
    """
    gamma_0 = all_zeros[idx]
    n_total = len(all_zeros)

    S1_bare = 0.0
    H1_bare = 0.0
    fail_count = 0

    for k in range(max(0, idx - n_max), min(n_total, idx + n_max + 1)):
        if k == idx:
            continue
        dg = gamma_0 - all_zeros[k]
        if abs(dg) < 1e-15:
            fail_count += 1
            continue
        S1_bare += 1.0 / dg
        H1_bare += 1.0 / (dg * dg)

    if fail_count > 0:
        print(f"  WARNING: 영점 {idx}에서 dg<1e-15 {fail_count}회 — 건너뜀")

    sm_im = gamma_smooth_im(gamma_0)

    A_bare = S1_bare ** 2 + 2.0 * H1_bare
    S1_L = S1_bare - sm_im
    A_L = S1_L ** 2 + 2.0 * H1_bare

    return {
        'S1_bare': S1_bare,
        'H1_bare': H1_bare,
        'sm_im': sm_im,
        'A_bare': A_bare,
        'A_L': A_L,
        'S1_L': S1_L,
    }


def main():
    print("=" * 70)
    print("  사이클 #278 — ζ(s) A_bare gap_min Spearman ρ_S 통합 재측정")
    print("=" * 70)
    print(f"  T범위: [{T_MIN}, {T_MAX}]  CENTER={CENTER}  mu={MU_LIST}  N_cond={N_COND}")
    print()

    # ── 1. 영점 수집 ────────────────────────────────────────────────
    all_zeros = get_zeta_zeros(T_MAX)
    zeros_in_range = [z for z in all_zeros if z >= T_MIN]
    n_all = len(all_zeros)
    n_range = len(zeros_in_range)
    print(f"  T≥{T_MIN}: {n_range}개, t ∈ [{zeros_in_range[0]:.3f}, {zeros_in_range[-1]:.3f}]")
    if n_range < 20:
        print("FATAL: 영점 부족 (<20)")
        sys.exit(1)

    dt_vals = [zeros_in_range[i+1]-zeros_in_range[i] for i in range(len(zeros_in_range)-1)]
    print(f"  dt_min={min(dt_vals):.4f}  dt_max={max(dt_vals):.4f}  dt_mean={np.mean(dt_vals):.4f}")

    # ── 2. A 계산 ───────────────────────────────────────────────────
    print(f"\n[A(γ)] 계산 ({n_range}개 × Σ)...")
    t0 = time.time()
    data = []
    for i, z in enumerate(zeros_in_range):
        all_idx = all_zeros.index(z)
        r = compute_A(zeros_in_range, all_zeros, all_idx)
        if not (math.isnan(r['A_bare']) or math.isinf(r['A_bare'])
                or r['A_bare'] <= 0 or r['A_L'] <= 0):
            r['t'] = z
            data.append(r)
        if (i + 1) % 50 == 0:
            print(f"  ...{i+1}/{n_range} ({time.time()-t0:.1f}s)")

    print(f"  완료: {len(data)}/{n_range} 유효 ({time.time()-t0:.1f}s)")
    if len(data) < 10:
        print("FATAL: 유효 데이터 부족 (<10)")
        sys.exit(1)

    # ── 3. 내부 영점 + 간격 계산 ────────────────────────────────────
    # GL(d)와 동일: data[2:-2] 내부 영점
    inner_data = data[2:-2]
    valid = []
    for pos, d in enumerate(inner_data):
        orig_pos = pos + 2  # data 내 원래 인덱스
        if orig_pos <= 0 or orig_pos >= len(data) - 1:
            continue
        t_n = d['t']
        t_prev = data[orig_pos - 1]['t']
        t_next = data[orig_pos + 1]['t']
        gap_r = t_next - t_n
        gap_l = t_n - t_prev
        if gap_r <= 0 or gap_l <= 0:
            continue
        # GUE 정규화: 국소 밀도 d̄ = 2/(t_{n+1} - t_{n-1})
        d_bar = 2.0 / (t_next - t_prev)
        d['gap_r_gue'] = gap_r * d_bar
        d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
        d['H1_frac_bare'] = (2.0 * d['H1_bare'] / d['A_bare']
                              if d['A_bare'] > 1e-10 else float('nan'))
        d['H1_frac_L'] = (2.0 * d['H1_bare'] / d['A_L']
                           if d['A_L'] > 1e-10 else float('nan'))
        valid.append(d)

    n_valid = len(valid)
    print(f"\n  내부 영점: {n_valid}")
    print(f"  t 범위: [{valid[0]['t']:.3f}, {valid[-1]['t']:.3f}]")
    if n_valid < 10:
        print("FATAL: 내부 영점 부족 (<10)")
        sys.exit(1)

    # ── 4. Spearman 상관 ─────────────────────────────────────────────
    A_bare_arr = np.array([d['A_bare'] for d in valid])
    A_L_arr    = np.array([d['A_L']    for d in valid])
    gr_arr     = np.array([d['gap_r_gue']   for d in valid])
    gm_arr     = np.array([d['gap_min_gue'] for d in valid])

    # NaN 체크
    for name, arr in [('A_bare', A_bare_arr), ('A_L', A_L_arr),
                       ('gap_r', gr_arr), ('gap_min', gm_arr)]:
        n_nan = int(np.sum(np.isnan(arr) | np.isinf(arr)))
        if n_nan > 0:
            print(f"  WARNING: {name}에 NaN/Inf {n_nan}개")

    rho_bare_min,  p_bare_min  = stats.spearmanr(A_bare_arr, gm_arr)
    rho_bare_r,    p_bare_r    = stats.spearmanr(A_bare_arr, gr_arr)
    rho_L_min,     p_L_min     = stats.spearmanr(A_L_arr,    gm_arr)
    rho_L_r,       p_L_r       = stats.spearmanr(A_L_arr,    gr_arr)

    sig = lambda p: '✅' if p < 0.001 else ('⚠️ ' if p < 0.01 else '❌')

    print(f"\n{'='*70}")
    print(f"  [Spearman 상관] n={n_valid}")
    print(f"{'='*70}")
    print(f"  ① ρ_S(A_bare, gap_min_GUE) = {rho_bare_min:+.4f}  p={p_bare_min:.3e}  {sig(p_bare_min)}"
          f"  ← 핵심 (GUE 비교: -0.857)")
    print(f"  ② ρ_S(A_bare, gap_right_GUE)= {rho_bare_r:+.4f}  p={p_bare_r:.3e}  {sig(p_bare_r)}"
          f"  ← 기존 -0.59 비교")
    print(f"  ③ ρ_S(A_L,    gap_min_GUE) = {rho_L_min:+.4f}  p={p_L_min:.3e}  {sig(p_L_min)}"
          f"  ← d≥3 비교표 통일")
    print(f"  ④ ρ_S(A_L,    gap_right_GUE)= {rho_L_r:+.4f}  p={p_L_r:.3e}  {sig(p_L_r)}"
          f"  ← 교차 검증")

    # ── 5. 보조 통계 ────────────────────────────────────────────────
    hf_bare = np.array([d['H1_frac_bare'] for d in valid])
    hf_L    = np.array([d['H1_frac_L']    for d in valid])
    sm_arr  = np.array([d['sm_im'] for d in valid])

    print(f"\n  [보조 통계]")
    print(f"  2H₁/A_bare:  {np.nanmean(hf_bare):.4f} ± {np.nanstd(hf_bare):.4f}")
    print(f"  2H₁/A_L:     {np.nanmean(hf_L):.4f} ± {np.nanstd(hf_L):.4f}")
    print(f"  <A_bare>:    {np.mean(A_bare_arr):.3f}")
    print(f"  <A_L>:       {np.mean(A_L_arr):.3f}")
    print(f"  <B_smooth>:  {np.mean(sm_arr):.4f} ± {np.std(sm_arr):.4f}")
    print(f"  <gap_min_GUE>: {np.mean(gm_arr):.4f} ± {np.std(gm_arr):.4f}")

    # gap_min vs gap_right 차이
    rho_mm_cross, _ = stats.spearmanr(gm_arr, gr_arr)
    print(f"  ρ_S(gap_min, gap_right): {rho_mm_cross:+.4f}")

    # ── 6. 판정 ─────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"  [판정]")
    gue_ref = -0.857  # GUE 기준값 (C-277)

    def star(rho, p):
        if p >= 0.001:
            return "★★ (비유의)"
        if abs(rho) >= 0.70:
            return "★★★★★ (GUE 근접, 감쇠 최소)"
        if abs(rho) >= 0.60:
            return "★★★★ (감쇠 있음)"
        if abs(rho) >= 0.45:
            return "★★★ (감쇠 보편)"
        return "★★★ (약한 상관)"

    print(f"  A_bare gap_min ρ_S = {rho_bare_min:+.4f}  →  {star(rho_bare_min, p_bare_min)}")
    print(f"  GUE 대비 gap_min 감쇠: {abs(gue_ref) - abs(rho_bare_min):.3f} "
          f"({'산술 감쇠 큼' if abs(gue_ref) - abs(rho_bare_min) > 0.2 else '산술 감쇠 작음'})")
    if abs(rho_bare_min) > abs(rho_bare_r):
        print(f"  gap_min이 gap_right보다 강함 ({abs(rho_bare_min):.3f} > {abs(rho_bare_r):.3f})")
    else:
        print(f"  gap_right이 gap_min보다 강함 ({abs(rho_bare_r):.3f} > {abs(rho_bare_min):.3f})")

    # ── 7. 결과 파일 저장 ────────────────────────────────────────────
    os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)
    with open(RESULT_PATH, 'w') as f:
        f.write("# C-278 ζ(s) A_bare gap_min Spearman ρ_S 통합 재측정\n")
        f.write(f"# T=[{T_MIN},{T_MAX}], CENTER={CENTER}, mu={MU_LIST}, N_cond={N_COND}\n")
        f.write(f"# n_all_zeros={n_all}, n_in_range={n_range}, n_valid={n_valid}\n")
        f.write(f"# same-side zeros only (GL(d)와 동일 방법론)\n\n")

        f.write("## 핵심 결과\n")
        f.write(f"rho_S(A_bare, gap_min_GUE)  = {rho_bare_min:+.6f}  p={p_bare_min:.4e}\n")
        f.write(f"rho_S(A_bare, gap_right_GUE)= {rho_bare_r:+.6f}  p={p_bare_r:.4e}\n")
        f.write(f"rho_S(A_L,    gap_min_GUE)  = {rho_L_min:+.6f}  p={p_L_min:.4e}\n")
        f.write(f"rho_S(A_L,    gap_right_GUE)= {rho_L_r:+.6f}  p={p_L_r:.4e}\n\n")

        f.write("## 보조 통계\n")
        f.write(f"n_valid = {n_valid}\n")
        f.write(f"2H1/A_bare_mean = {np.nanmean(hf_bare):.6f}\n")
        f.write(f"2H1/A_L_mean    = {np.nanmean(hf_L):.6f}\n")
        f.write(f"mean_A_bare     = {np.mean(A_bare_arr):.4f}\n")
        f.write(f"mean_A_L        = {np.mean(A_L_arr):.4f}\n")
        f.write(f"mean_Bsmooth    = {np.mean(sm_arr):.6f}\n")
        f.write(f"mean_gap_min_GUE = {np.mean(gm_arr):.6f}\n")
        f.write(f"GUE_ref_rho_S   = {gue_ref}\n")
        f.write(f"arithmetic_damping = {abs(gue_ref) - abs(rho_bare_min):.4f}\n\n")

        f.write("## 개별 데이터 (t, A_bare, A_L, gap_min_GUE, gap_right_GUE, B_smooth)\n")
        f.write("t,A_bare,A_L,gap_min_gue,gap_right_gue,sm_im\n")
        for d in valid:
            f.write(f"{d['t']:.6f},{d['A_bare']:.6f},{d['A_L']:.6f},"
                    f"{d['gap_min_gue']:.6f},{d['gap_r_gue']:.6f},{d['sm_im']:.6f}\n")

    print(f"\n  결과 저장: {RESULT_PATH}")
    print("=" * 70)
    print("  완료.")
    print("=" * 70)


if __name__ == '__main__':
    main()
