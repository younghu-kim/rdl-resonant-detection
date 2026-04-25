#!/usr/bin/env python3
"""
=============================================================================
[사이클 #266] Gamma 잔차 분리: A_Λ vs A_L 정량화
=============================================================================

목적:
  GUE가 A-gap ρ ≈ -0.50 예측, 관측 ρ ≈ -0.578. 잔차 Δ≈0.08의 기원 분리.
  같은 ζ(s) 영점들에서 A_L (zero-sum, Gamma 없음)과 A_Λ (Gamma 포함)을
  동시에 계산하여, Gamma factor의 상관 기여를 정량화.

핵심 관계:
  c₀^Λ = c₀^L + (1/2)ψ(ρ/2) - (1/2)log(π)
  c₁^Λ = c₁^L + (1/4)ψ'(ρ/2)
  A_Λ = Im(c₀^Λ)² + 2Re(c₁^Λ)
  A_L = Im(c₀^L)² + 2Re(c₁^L)   ← 새로 계산
  ΔA = A_Λ - A_L = Gamma contribution

측정 항목 (수학자 보드 지시):
  1. ρ(A_L, gap_min), ρ(A_L, gap_right)
  2. ρ(A_Λ, gap_min), ρ(A_Λ, gap_right)
  3. ρ(ΔA, gap_right) ← 핵심: Gamma 기여분의 gap 상관
  4. ΔA/A_Λ 비율 (Gamma가 A에서 차지하는 비중)
  5. Δρ = ρ(A_Λ) - ρ(A_L): Gamma에 의한 상관 증폭량

결과 파일: results/gamma_residual_c266.txt
=============================================================================
"""

import sys, os, time
import numpy as np
from scipy import stats
import mpmath

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CENTER = 0.5
T_MAX = 2000
TRIM_FRAC = 0.20     # 양쪽 20% 절단 → 중앙 60%

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gamma_residual_c266.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

_log_buf = []
def log(msg=''):
    print(msg, flush=True)
    _log_buf.append(str(msg))

def save():
    with open(RESULT_PATH, 'w') as f:
        f.write('\n'.join(_log_buf))


def main():
    t_start = time.time()
    log(f"[사이클 #266] Gamma 잔차 분리: A_Λ vs A_L")
    log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log("=" * 72)

    # ──────────────────────────────────────────────
    # 1. 영점 수집 (PARI)
    # ──────────────────────────────────────────────
    log()
    log("[1] ζ(s) 영점 수집 (T=2000)...")
    pari = cypari2.Pari()
    pari.allocatemem(128 * 10**6)
    pari.set_real_precision(38)

    pari('Li_z = lfuninit(lfuncreate(1), [0, 2100])')
    pari('zv = lfunzeros(Li_z, 2000)')
    n_z = int(str(pari('#zv')))
    zeros = []
    for i in range(1, n_z + 1):
        t = float(str(pari(f'zv[{i}]')).replace(' E', 'e'))
        if t > 5:
            zeros.append(t)
    zeros = sorted(zeros)
    log(f"  {len(zeros)}개 영점, t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
    log(f"  수집 시간: {time.time()-t_start:.1f}초")

    # ──────────────────────────────────────────────
    # 2. Hadamard 분해: A_L (zero-sum) + A_Λ (Gamma 포함) 동시 계산
    # ──────────────────────────────────────────────
    log()
    log("[2] Hadamard 분해: A_L vs A_Λ 동시 계산...")
    mpmath.mp.dps = 30

    data = []
    for i in range(len(zeros)):
        g0 = zeros[i]

        # --- Zero-sum: 같은 부호 영점 ---
        S1_same = 0.0
        H1_same = 0.0
        for k in range(len(zeros)):
            if k == i:
                continue
            diff = g0 - zeros[k]
            S1_same += 1.0 / diff
            H1_same += 1.0 / diff**2

        # --- 반사 영점 (conjugate γ → -γ) ---
        S1_conj = sum(1.0 / (g0 + zeros[k]) for k in range(len(zeros)))
        H1_conj = sum(1.0 / (g0 + zeros[k])**2 for k in range(len(zeros)))

        # --- A_L: primitive L, Gamma 배제 ---
        S1_L = S1_same + S1_conj
        H1_L = H1_same + H1_conj
        A_L = S1_L**2 + 2 * H1_L

        # --- Gamma 보정 ---
        s = mpmath.mpc(CENTER, g0)
        # c₀ 보정: (1/2)ψ(s/2) - (1/2)log(π)
        gpg = -mpmath.log(mpmath.pi)/2 + mpmath.digamma(s/2)/2
        im_gpg = float(mpmath.im(gpg))  # Im(Gamma contribution to c₀)
        re_gpg = float(mpmath.re(gpg))  # Re(Gamma contribution to c₀)
        # c₁ 보정: (1/4)ψ'(s/2)
        psi1 = mpmath.psi(1, s/2)
        re_correction_c1 = float(mpmath.re(psi1)) / 4.0

        # --- A_Λ: completed Λ, Gamma 포함 ---
        S1_Lambda = S1_L - im_gpg  # Im(c₀^Λ) = S1_L - Im(gamma_shift)
        H1_Lambda = H1_L + re_correction_c1  # Re(c₁^Λ) = H1_L + Re(gamma_c₁)
        A_Lambda = S1_Lambda**2 + 2 * H1_Lambda

        # --- ΔA = A_Λ - A_L ---
        Delta_A = A_Lambda - A_L

        # --- Gaps ---
        if 0 < i < len(zeros) - 1:
            gap_r = zeros[i+1] - g0
            gap_l = g0 - zeros[i-1]
            gap_min = min(gap_r, gap_l)

            # GUE 정규화
            d_bar = np.log(g0 / (2*np.pi)) / (2*np.pi)
            gap_r_gue = gap_r * d_bar
            gap_l_gue = gap_l * d_bar
            gap_min_gue = gap_min * d_bar
        else:
            gap_r = gap_l = gap_min = np.nan
            gap_r_gue = gap_l_gue = gap_min_gue = np.nan

        if A_Lambda > 0 and A_L > 0 and not np.isnan(gap_r_gue):
            data.append({
                't': g0,
                'A_L': A_L,
                'A_Lambda': A_Lambda,
                'Delta_A': Delta_A,
                'S1_L': S1_L,
                'S1_Lambda': S1_Lambda,
                'H1_L': H1_L,
                'H1_Lambda': H1_Lambda,
                'im_gamma_c0': im_gpg,
                're_gamma_c1': re_correction_c1,
                'gap_r_gue': gap_r_gue,
                'gap_min_gue': gap_min_gue,
                'gap_r_raw': gap_r,
                'gap_min_raw': gap_min,
            })

        if (i+1) % 300 == 0:
            log(f"  {i+1}/{len(zeros)} 완료 ({time.time()-t_start:.0f}s)")

    log(f"  전체 유효: {len(data)}")

    # ──────────────────────────────────────────────
    # 3. 중앙 60% 선택
    # ──────────────────────────────────────────────
    N = len(data)
    lo = int(N * TRIM_FRAC)
    hi = int(N * (1.0 - TRIM_FRAC))
    valid = data[lo:hi]
    log(f"  중앙 60% 선택: {len(valid)}개")
    log()

    # ──────────────────────────────────────────────
    # 4. 핵심 상관 분석
    # ──────────────────────────────────────────────
    log("=" * 72)
    log("[4] 핵심 상관 분석: A_L vs A_Λ vs ΔA")
    log("=" * 72)
    log()

    A_L_arr = np.array([d['A_L'] for d in valid])
    A_Lam_arr = np.array([d['A_Lambda'] for d in valid])
    DA_arr = np.array([d['Delta_A'] for d in valid])
    gr_arr = np.array([d['gap_r_gue'] for d in valid])
    gm_arr = np.array([d['gap_min_gue'] for d in valid])
    t_arr = np.array([d['t'] for d in valid])
    im_gc0_arr = np.array([d['im_gamma_c0'] for d in valid])
    re_gc1_arr = np.array([d['re_gamma_c1'] for d in valid])

    sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')
    n = len(valid)
    log(f"  n = {n}")
    log()

    # --- 4.1: A_Λ 상관 (기존 결과 재확인) ---
    log("  [4.1] A_Λ (completed, Gamma 포함)")
    pairs_lam = [
        ("A_Λ, gap_right_GUE", A_Lam_arr, gr_arr),
        ("A_Λ, gap_min_GUE", A_Lam_arr, gm_arr),
    ]
    rho_lam = {}
    for name, x, y in pairs_lam:
        r, p = stats.spearmanr(x, y)
        log(f"    ρ({name:30s}) = {r:+.4f}  (p={p:.2e})  {sig(p)}")
        rho_lam[name] = r
    log()

    # --- 4.2: A_L 상관 (zero-sum, Gamma 배제) ---
    log("  [4.2] A_L (primitive, Gamma 배제)")
    pairs_l = [
        ("A_L, gap_right_GUE", A_L_arr, gr_arr),
        ("A_L, gap_min_GUE", A_L_arr, gm_arr),
    ]
    rho_l = {}
    for name, x, y in pairs_l:
        r, p = stats.spearmanr(x, y)
        log(f"    ρ({name:30s}) = {r:+.4f}  (p={p:.2e})  {sig(p)}")
        rho_l[name] = r
    log()

    # --- 4.3: ΔA = A_Λ - A_L (Gamma 기여분) 상관 ← 핵심 ---
    log("  [4.3] ΔA = A_Λ - A_L (Gamma 기여분) ← 핵심")
    pairs_da = [
        ("ΔA, gap_right_GUE", DA_arr, gr_arr),
        ("ΔA, gap_min_GUE", DA_arr, gm_arr),
    ]
    rho_da = {}
    for name, x, y in pairs_da:
        r, p = stats.spearmanr(x, y)
        log(f"    ρ({name:30s}) = {r:+.4f}  (p={p:.2e})  {sig(p)}")
        rho_da[name] = r
    log()

    # --- 4.4: Gamma 성분 개별 상관 ---
    log("  [4.4] Gamma 성분 개별")
    pairs_gc = [
        ("Im(γ_c₀), gap_right_GUE", im_gc0_arr, gr_arr),
        ("Im(γ_c₀), gap_min_GUE", im_gc0_arr, gm_arr),
        ("Re(γ_c₁), gap_right_GUE", re_gc1_arr, gr_arr),
        ("Re(γ_c₁), gap_min_GUE", re_gc1_arr, gm_arr),
        ("Im(γ_c₀), t", im_gc0_arr, t_arr),
        ("Re(γ_c₁), t", re_gc1_arr, t_arr),
    ]
    for name, x, y in pairs_gc:
        r, p = stats.spearmanr(x, y)
        log(f"    ρ({name:30s}) = {r:+.4f}  (p={p:.2e})  {sig(p)}")
    log()

    # ──────────────────────────────────────────────
    # 5. 정량적 요약
    # ──────────────────────────────────────────────
    log("=" * 72)
    log("[5] 정량적 요약")
    log("=" * 72)
    log()

    # 5.1: Gamma 기여 비율
    ratio_arr = DA_arr / A_Lam_arr
    log("  [5.1] Gamma 기여 비율: ΔA / A_Λ")
    log(f"    평균: {np.mean(ratio_arr):.4f}")
    log(f"    중위: {np.median(ratio_arr):.4f}")
    log(f"    표준편차: {np.std(ratio_arr):.4f}")
    log(f"    범위: [{np.min(ratio_arr):.4f}, {np.max(ratio_arr):.4f}]")
    log()

    # 5.2: 상관 증폭량
    delta_rho_gr = rho_lam["A_Λ, gap_right_GUE"] - rho_l["A_L, gap_right_GUE"]
    delta_rho_gm = rho_lam["A_Λ, gap_min_GUE"] - rho_l["A_L, gap_min_GUE"]
    log("  [5.2] Gamma에 의한 상관 변화 Δρ")
    log(f"    Δρ(gap_right) = ρ(A_Λ) - ρ(A_L) = {rho_lam['A_Λ, gap_right_GUE']:+.4f} - ({rho_l['A_L, gap_right_GUE']:+.4f}) = {delta_rho_gr:+.4f}")
    log(f"    Δρ(gap_min)   = ρ(A_Λ) - ρ(A_L) = {rho_lam['A_Λ, gap_min_GUE']:+.4f} - ({rho_l['A_L, gap_min_GUE']:+.4f}) = {delta_rho_gm:+.4f}")
    log()

    # 5.3: 3요인 분해
    rho_gue = -0.50  # GUE 예측 (C-265)
    rho_obs_gr = rho_lam["A_Λ, gap_right_GUE"]
    rho_gamma_effect = delta_rho_gr
    rho_arithmetic = rho_obs_gr - rho_gue - rho_gamma_effect

    log("  [5.3] 3요인 분해 (gap_right 기준)")
    log(f"    ρ_GUE     = {rho_gue:+.4f}  (GUE 이론, C-265)")
    log(f"    Δρ_Gamma  = {rho_gamma_effect:+.4f}  (Gamma 기여)")
    log(f"    Δρ_arith  = {rho_arithmetic:+.4f}  (산술 잔차)")
    log(f"    합계      = {rho_gue + rho_gamma_effect + rho_arithmetic:+.4f}")
    log(f"    관측      = {rho_obs_gr:+.4f}")
    log()

    # 5.4: GUE 설명 비율 갱신
    if abs(rho_obs_gr) > 0:
        gue_pct = abs(rho_gue / rho_obs_gr) * 100
        gamma_pct = abs(rho_gamma_effect / rho_obs_gr) * 100
        arith_pct = abs(rho_arithmetic / rho_obs_gr) * 100
        log("  [5.4] 기여 비율")
        log(f"    GUE:    {gue_pct:.1f}%")
        log(f"    Gamma:  {gamma_pct:.1f}%")
        log(f"    산술:   {arith_pct:.1f}%")
    log()

    # 5.5: A_L과 A_Λ 기술 통계 비교
    log("  [5.5] A_L vs A_Λ 기술 통계")
    log(f"    {'':>12} {'A_L':>12} {'A_Λ':>12} {'ΔA':>12}")
    log(f"    {'평균':>12} {np.mean(A_L_arr):>12.4f} {np.mean(A_Lam_arr):>12.4f} {np.mean(DA_arr):>12.4f}")
    log(f"    {'중위':>12} {np.median(A_L_arr):>12.4f} {np.median(A_Lam_arr):>12.4f} {np.median(DA_arr):>12.4f}")
    log(f"    {'표준편차':>12} {np.std(A_L_arr):>12.4f} {np.std(A_Lam_arr):>12.4f} {np.std(DA_arr):>12.4f}")
    log()

    # 5.6: ΔA의 높이 의존성
    log("  [5.6] ΔA 높이 의존성")
    r_da_t, p_da_t = stats.spearmanr(DA_arr, t_arr)
    log(f"    ρ(ΔA, t) = {r_da_t:+.4f}  (p={p_da_t:.2e})")
    r_ratio_t, p_ratio_t = stats.spearmanr(ratio_arr, t_arr)
    log(f"    ρ(ΔA/A_Λ, t) = {r_ratio_t:+.4f}  (p={p_ratio_t:.2e})")
    log()

    # 5.7: Pearson 상관 (선형 관계 확인)
    log("  [5.7] Pearson 상관 (선형)")
    r_p_lam, _ = stats.pearsonr(A_Lam_arr, gr_arr)
    r_p_l, _ = stats.pearsonr(A_L_arr, gr_arr)
    r_p_da, _ = stats.pearsonr(DA_arr, gr_arr)
    log(f"    r_Pearson(A_Λ, gap_right) = {r_p_lam:+.4f}")
    log(f"    r_Pearson(A_L, gap_right) = {r_p_l:+.4f}")
    log(f"    r_Pearson(ΔA, gap_right)  = {r_p_da:+.4f}")
    log()

    # ──────────────────────────────────────────────
    # 6. 판정
    # ──────────────────────────────────────────────
    log("=" * 72)
    log("[6] 판정")
    log("=" * 72)
    log()

    # 판정 기준
    da_gap_sig = abs(rho_da.get("ΔA, gap_right_GUE", 0)) > 0.1
    residual_small = abs(rho_arithmetic) < 0.02

    if da_gap_sig and residual_small:
        verdict = "★★★★ Gamma 기여 확립 + 3요인 분해 완결"
    elif da_gap_sig:
        verdict = "★★★ Gamma 기여 확립, 잔차 미해소"
    elif residual_small:
        verdict = "★★★ Gamma 무상관이나 GUE+산술로 설명"
    else:
        verdict = "★★ 부분 양성, 추가 분석 필요"

    log(f"  판정: {verdict}")
    log()
    log(f"  핵심 수치:")
    log(f"    ρ(A_L, gap_right)  = {rho_l.get('A_L, gap_right_GUE', 0):+.4f}  (Gamma 없는 상관)")
    log(f"    ρ(A_Λ, gap_right)  = {rho_lam.get('A_Λ, gap_right_GUE', 0):+.4f}  (Gamma 포함 상관)")
    log(f"    ρ(ΔA, gap_right)   = {rho_da.get('ΔA, gap_right_GUE', 0):+.4f}  (Gamma 기여분)")
    log(f"    Δρ_Gamma           = {delta_rho_gr:+.4f}")
    log(f"    GUE 예측           = {rho_gue:+.4f}")
    log(f"    관측               = {rho_obs_gr:+.4f}")
    log(f"    잔차               = {rho_arithmetic:+.4f}")
    log()
    log(f"  소요: {time.time()-t_start:.1f}초")

    save()
    log()
    log("결과 저장 완료.")


if __name__ == '__main__':
    main()
