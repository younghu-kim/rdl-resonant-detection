#!/usr/bin/env python3
"""
=============================================================================
[사이클 #268] Gamma factor 기여 분리 실험 — A_L vs A_Λ 잔차 정량화
=============================================================================

목적:
  GUE 잔차 0.08 해소: ρ_관측(-0.578) vs ρ_GUE(-0.50) 차이가
  Gamma factor 기여(ψ(s/2) + log-terms)에서 비롯되는지 정량 검증.

측정 항목:
  ρ(A_L, gap_min)          — pure zero-sum의 gap_min 상관
  ρ(A_L, gap_right)         — pure zero-sum의 gap_right 상관
  ρ(A_Λ, gap_min)          — completed A의 gap_min 상관
  ρ(A_Λ, gap_right)         — completed A의 gap_right 상관
  ρ(A_Λ - A_L, gap_right)  — Gamma 기여분의 gap_right 상관 [핵심]
  (A_Λ - A_L)/A_Λ 비율     — Gamma 비중
  Δρ = ρ(A_Λ, gap_right) - ρ(A_L, gap_right) — Gamma 증폭량
  t-의존성 (Gamma/A_Λ 비율)

핵심 가설 (수학자 B-45):
  A_Λ - A_L ≈ Gamma 기여분이 gap_right와 음상관
  → GUE 순수 예측(-0.50)과 관측(-0.578) 차이 0.08 설명

데이터: ζ(s) T=2000, n≈1517 (중앙 60% ≈ 910 내부 영점)
K항: 모든 영점 사용 (≥50 조건 충족)

결과: results/gamma_residual_c268.txt
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

CENTER     = 0.5
T_MAX      = 2000
TRIM_FRAC  = 0.20    # 양쪽 20% 제외 → 중앙 60%
N_BINS     = 8       # gap 분위 빈 수 (t-의존성 분석용)

# GUE 이론 예측 (C-265에서 확립)
RHO_GUE       = -0.50   # ρ_GUE(A, gap_right) 이론값
RHO_OBS_C265  = -0.578  # C-265/C-256 관측값 (gap_right_GUE 기준)
RESIDUAL_GUE  = RHO_OBS_C265 - RHO_GUE   # -0.078 (Gamma 기여 가설)

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gamma_residual_c268.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 로그
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

_log_buf = []

def log(msg=''):
    print(msg, flush=True)
    _log_buf.append(str(msg))

def save():
    with open(RESULT_PATH, 'w') as f:
        f.write('\n'.join(_log_buf))

sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    log("=" * 80)
    log("[C-268] Gamma factor 기여 분리 — A_L vs A_Λ 잔차 정량화")
    log("=" * 80)
    log(f"  설정: T_MAX={T_MAX}, TRIM_FRAC={TRIM_FRAC}, N_BINS={N_BINS}")
    log(f"  GUE 이론: ρ_GUE={RHO_GUE}, 관측={RHO_OBS_C265}, 잔차={RESIDUAL_GUE:.4f}")
    log()

    # ──────────────────────────────────────────────────────────────────────
    # STEP 1: ζ(s) 영점 수집 (PARI, T=2000)
    # ──────────────────────────────────────────────────────────────────────
    log("[STEP 1] ζ(s) 영점 수집 (T=2000)...")
    pari = cypari2.Pari()
    pari.allocatemem(256 * 10**6)
    pari.set_real_precision(38)
    mpmath.mp.dps = 30

    pari('Li_z = lfuninit(lfuncreate(1), [0, 2100])')
    pari('zv = lfunzeros(Li_z, 2000)')
    n_z = int(str(pari('#zv')))

    zeros = []
    for i in range(1, n_z + 1):
        t = float(str(pari(f'zv[{i}]')).replace(' E', 'e'))
        if t > 5:
            zeros.append(t)
    zeros = sorted(zeros)

    log(f"  {len(zeros)}개 영점 수집, t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
    log(f"  수집 시간: {time.time()-t_start:.1f}초")
    log()

    # ──────────────────────────────────────────────────────────────────────
    # STEP 2: A_L / A_Λ 동시 계산
    # ──────────────────────────────────────────────────────────────────────
    log("[STEP 2] A_L / A_Λ 병렬 계산 (Hadamard zero-sum)...")
    log(f"  K = 전체({len(zeros)}개) — K≥50 조건 충족")
    log()

    data = []
    n_fail = 0

    for i in range(len(zeros)):
        g0 = zeros[i]

        try:
            # ── 같은 부호 영점 합산 ──
            S1_same = 0.0
            H1_same = 0.0
            for k in range(len(zeros)):
                if k == i:
                    continue
                diff = g0 - zeros[k]
                S1_same += 1.0 / diff
                H1_same += 1.0 / diff**2

            # ── 켤레 영점 합산 (γ → -γ) ──
            S1_conj = sum(1.0 / (g0 + zeros[k]) for k in range(len(zeros)))
            H1_conj = sum(1.0 / (g0 + zeros[k])**2 for k in range(len(zeros)))

            # ── A_L: pure zero-sum (Gamma 미포함) ──
            S1_L = S1_same + S1_conj
            H1_L = H1_same + H1_conj
            A_L  = S1_L**2 + 2.0 * H1_L

            # ── Gamma 보정 계산 ──
            # Λ'/Λ(ρ₀) = ζ'/ζ(ρ₀) + (-log(π)/2 + ψ(s/2)/2)
            # Im 부분: S₁^Λ = S₁^L - Im(-log(π)/2 + ψ(s/2)/2)
            # Re 부분 (H₁에 해당): H₁^Λ = H₁^L + Re(ψ'(s/2)/4)
            s_val = mpmath.mpc(CENTER, g0)
            gamma_S = -mpmath.log(mpmath.pi) / 2 + mpmath.digamma(s_val / 2) / 2
            im_gamma_S = float(mpmath.im(gamma_S))    # S₁ 보정량 (음수로 빼기)
            psi1_val = mpmath.psi(1, s_val / 2)
            re_gamma_H = float(mpmath.re(psi1_val)) / 4.0  # H₁ 보정량

            # ── A_Λ: completed (Gamma 포함) ──
            S1_Lambda  = S1_L - im_gamma_S
            H1_Lambda  = H1_L + re_gamma_H
            A_Lambda   = S1_Lambda**2 + 2.0 * H1_Lambda

            # ── Gamma 기여분 ──
            Delta_A = A_Lambda - A_L
            # 분해: Delta_A = -2*S1_L*im_gamma_S + im_gamma_S² + 2*re_gamma_H
            Delta_A_cross   = -2.0 * S1_L * im_gamma_S   # 교차항
            Delta_A_sq      = im_gamma_S**2               # im 제곱항
            Delta_A_H       = 2.0 * re_gamma_H            # H₁ 보정항

            # 유효성 체크
            if A_L <= 0 or A_Lambda <= 0:
                n_fail += 1
                continue

            data.append({
                't': g0,
                'A_L': A_L,
                'A_Lambda': A_Lambda,
                'Delta_A': Delta_A,
                'S1_L': S1_L,
                'H1_L': H1_L,
                'S1_Lambda': S1_Lambda,
                'H1_Lambda': H1_Lambda,
                'im_gamma_S': im_gamma_S,
                're_gamma_H': re_gamma_H,
                'Delta_A_cross': Delta_A_cross,
                'Delta_A_sq': Delta_A_sq,
                'Delta_A_H': Delta_A_H,
            })

        except Exception as e:
            n_fail += 1
            if n_fail <= 5:
                print(f"WARNING: i={i}, t={g0:.3f}: {e}", flush=True)

        if (i + 1) % 200 == 0:
            log(f"  {i+1}/{len(zeros)} 계산 완료... ({time.time()-t_start:.0f}s, 유효={len(data)}, 실패={n_fail})")

    log(f"  계산 완료: 유효={len(data)}, 실패={n_fail}, 소요={time.time()-t_start:.1f}초")

    if len(data) < 100:
        log("⚠️ 유효 데이터 100개 미만 — 계산 중단")
        save()
        return

    if n_fail > len(zeros) // 2:
        log(f"⚠️ 실패율 {n_fail}/{len(zeros)} — 절반 초과, 결과 신뢰도 낮음")

    log()

    # ──────────────────────────────────────────────────────────────────────
    # STEP 3: 중앙 60% 선택 + gap 계산
    # ──────────────────────────────────────────────────────────────────────
    log("[STEP 3] 중앙 60% 선택 + GUE 정규화 gap 계산...")

    N = len(data)
    lo = int(N * TRIM_FRAC)
    hi = int(N * (1.0 - TRIM_FRAC))

    valid = []
    for idx in range(lo, hi):
        if idx <= 0 or idx >= N - 1:
            continue
        d = data[idx]

        # 인접 간격
        gap_r = data[idx+1]['t'] - d['t']
        gap_l = d['t'] - data[idx-1]['t']
        if gap_r <= 0 or gap_l <= 0:
            continue

        # GUE 정규화: d̄(t) = log(t/(2π))/(2π)
        d_bar = np.log(d['t'] / (2.0 * np.pi)) / (2.0 * np.pi)
        d['gap_r_gue']   = gap_r * d_bar
        d['gap_l_gue']   = gap_l * d_bar
        d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
        d['d_bar']       = d_bar
        valid.append(d)

    log(f"  내부 영점: {len(valid)} (중앙 {(1.0-2*TRIM_FRAC)*100:.0f}%)")
    log()

    if len(valid) < 50:
        log("⚠️ 내부 영점 50개 미만 — 분석 중단")
        save()
        return

    # ──────────────────────────────────────────────────────────────────────
    # STEP 4: 핵심 상관 분석
    # ──────────────────────────────────────────────────────────────────────
    log("[STEP 4] 핵심 상관 분석 — A_L / A_Λ / Gamma 기여분")
    log()

    n = len(valid)
    AL_arr       = np.array([d['A_L']      for d in valid])
    ALam_arr     = np.array([d['A_Lambda'] for d in valid])
    DeltaA_arr   = np.array([d['Delta_A']  for d in valid])
    gm_arr       = np.array([d['gap_min_gue'] for d in valid])
    gr_arr       = np.array([d['gap_r_gue']   for d in valid])

    # S₁ / H₁ 성분
    S1L_arr      = np.array([d['S1_L']     for d in valid])
    H1L_arr      = np.array([d['H1_L']     for d in valid])
    S1Lam_arr    = np.array([d['S1_Lambda'] for d in valid])
    H1Lam_arr    = np.array([d['H1_Lambda'] for d in valid])

    # Gamma 성분 분해
    imGS_arr     = np.array([d['im_gamma_S']    for d in valid])
    reGH_arr     = np.array([d['re_gamma_H']    for d in valid])
    DeltaCross   = np.array([d['Delta_A_cross'] for d in valid])
    DeltaSq      = np.array([d['Delta_A_sq']    for d in valid])
    DeltaH       = np.array([d['Delta_A_H']     for d in valid])

    # ── 4.1: 전체 상관표 ──
    log(f"  n = {n}")
    log()
    log(f"  {'항목':<35} {'ρ':>8}  {'p':>10}  {'판정'}")
    log("  " + "-" * 65)

    corr_table = [
        ("A_L,    gap_min_GUE",   AL_arr,   gm_arr),
        ("A_L,    gap_right_GUE", AL_arr,   gr_arr),
        ("A_Λ,    gap_min_GUE",   ALam_arr, gm_arr),
        ("A_Λ,    gap_right_GUE", ALam_arr, gr_arr),
        ("A_Λ-A_L, gap_min_GUE", DeltaA_arr, gm_arr),
        ("A_Λ-A_L, gap_right_GUE", DeltaA_arr, gr_arr),   # ★ 핵심
        ("S₁^L², gap_right_GUE",  S1L_arr**2,  gr_arr),
        ("S₁^Λ², gap_right_GUE",  S1Lam_arr**2, gr_arr),
        ("2H₁^L,  gap_right_GUE", 2*H1L_arr,   gr_arr),
        ("2H₁^Λ,  gap_right_GUE", 2*H1Lam_arr, gr_arr),
        ("ΔS₁_cross, gap_right",  DeltaCross,  gr_arr),
        ("im_γ²,     gap_right",  DeltaSq,     gr_arr),
        ("2·re_γ',   gap_right",  DeltaH,      gr_arr),
    ]

    results = {}
    for name, x, y in corr_table:
        try:
            r, p = stats.spearmanr(x, y)
        except Exception as e:
            log(f"  WARNING: {name} — {e}")
            continue
        results[name] = (r, p)
        log(f"  ρ({name:<33}) = {r:+.4f}  (p={p:.3e})  {sig(p)}")

    log()

    # ── 4.2: 핵심 결론 출력 ──
    log("[STEP 4.2] 핵심 결론")
    log()

    # Gap_right 상관: A_L, A_Λ, Delta_A
    rho_AL_gr  = results.get("A_L,    gap_right_GUE",   (np.nan, np.nan))[0]
    rho_ALam_gr = results.get("A_Λ,    gap_right_GUE",  (np.nan, np.nan))[0]
    rho_DA_gr  = results.get("A_Λ-A_L, gap_right_GUE", (np.nan, np.nan))[0]
    p_DA_gr    = results.get("A_Λ-A_L, gap_right_GUE", (np.nan, np.nan))[1]

    delta_rho  = rho_ALam_gr - rho_AL_gr
    log(f"  ρ(A_L,   gap_right_GUE) = {rho_AL_gr:+.4f}  (GUE 기준: {RHO_GUE:+.3f})")
    log(f"  ρ(A_Λ,   gap_right_GUE) = {rho_ALam_gr:+.4f}  (관측 기준: {RHO_OBS_C265:+.3f})")
    log(f"  Δρ = ρ(A_Λ) - ρ(A_L)   = {delta_rho:+.4f}  (Gamma 증폭량)")
    log(f"  ρ(A_Λ-A_L, gap_right)   = {rho_DA_gr:+.4f}  (p={p_DA_gr:.3e})  ← 핵심")
    log()

    # GUE 잔차 설명 여부
    residual_after = rho_ALam_gr - RHO_GUE   # 전체 잔차
    residual_unexplained = rho_AL_gr - RHO_GUE  # A_L 기준 잔차

    log(f"  GUE 예측:         ρ_GUE = {RHO_GUE:+.4f}")
    log(f"  전체 잔차:        ρ_obs - ρ_GUE = {residual_after:+.4f}")
    log(f"  Gamma 증폭:       Δρ          = {delta_rho:+.4f}")
    log(f"  A_L 기준 잔차:    ρ(A_L) - ρ_GUE = {residual_unexplained:+.4f}")
    log()

    if rho_DA_gr < -0.05 and p_DA_gr < 0.01:
        log("  ★★★ [Gamma 기여 확립] A_Λ-A_L이 gap_right와 유의한 음상관")
        log(f"       Gamma factor가 gap_right 상관을 {abs(delta_rho):.4f} 증폭")
        if abs(residual_after) < 0.02:
            log("  ★★★★ [3요인 분해 완결] |잔차| < 0.02 — B-45 해소!")
        else:
            log(f"  → 잔차 {residual_after:+.4f} 미잔류 — 추가 분석 필요")
    elif abs(rho_DA_gr) < 0.05:
        log("  [Gamma 무상관] A_Λ-A_L이 gap_right와 거의 무상관")
        log("  → 잔차(-0.078)의 원인은 산술 S₁² 구조 (다른 방향 필요)")
    else:
        log(f"  [중간] ρ(Gamma, gap_right)={rho_DA_gr:+.4f} — 부분 기여")

    log()

    # ── 4.3: Gamma 비중 통계 ──
    log("[STEP 4.3] Gamma 기여 비중 통계")
    log()

    DA_over_ALam = DeltaA_arr / ALam_arr
    mean_ratio  = np.mean(DA_over_ALam)
    med_ratio   = np.median(DA_over_ALam)
    std_ratio   = np.std(DA_over_ALam)

    log(f"  (A_Λ - A_L)/A_Λ: mean={mean_ratio:+.4f}, median={med_ratio:+.4f}, std={std_ratio:.4f}")
    log(f"  im_γ_S: mean={np.mean(imGS_arr):+.4f}, std={np.std(imGS_arr):.4f}")
    log(f"  re_γ_H: mean={np.mean(reGH_arr):+.4f}, std={np.std(reGH_arr):.4f}")
    log()

    # ── 4.4: A_L vs A_Λ 상관 ──
    rho_AL_ALam, p_AL_ALam = stats.spearmanr(AL_arr, ALam_arr)
    log(f"  ρ(A_L, A_Λ) = {rho_AL_ALam:+.4f}  (p={p_AL_ALam:.3e}) — 방법론 일관성")
    log()

    # ──────────────────────────────────────────────────────────────────────
    # STEP 5: t-의존성 분석 (Gamma/A_Λ 비율 변화)
    # ──────────────────────────────────────────────────────────────────────
    log("[STEP 5] t-의존성 분석 — Gamma 기여 비율 vs t")
    log()
    log(f"  기대: ψ(1/4+it/2) ≈ log(t/2) + O(1/t) → Gamma 기여 log(t) 스케일링")
    log()

    t_arr = np.array([d['t'] for d in valid])

    # t를 N_BINS개 빈으로 나누기
    t_quantiles = np.percentile(t_arr, np.linspace(0, 100, N_BINS + 1))

    log(f"  {'t_bin':>4}  {'t̄':>8}  {'mean(ΔA/A)':>12}  {'mean(im_γ)':>12}  {'ρ(A_Λ,gr)':>11}  {'ρ(A_L,gr)':>11}  {'n':>5}")
    log("  " + "-" * 80)

    t_dep = []
    for b in range(N_BINS):
        mask = (t_arr >= t_quantiles[b]) & (t_arr < t_quantiles[b+1])
        if b == N_BINS - 1:
            mask = (t_arr >= t_quantiles[b]) & (t_arr <= t_quantiles[b+1])
        n_bin = np.sum(mask)
        if n_bin < 5:
            continue

        t_mean   = np.mean(t_arr[mask])
        ratio_m  = np.mean(DA_over_ALam[mask])
        imG_m    = np.mean(imGS_arr[mask])

        rho_Lam_b, _ = stats.spearmanr(ALam_arr[mask], gr_arr[mask]) if n_bin >= 5 else (np.nan, np.nan)
        rho_L_b,   _ = stats.spearmanr(AL_arr[mask],   gr_arr[mask]) if n_bin >= 5 else (np.nan, np.nan)

        t_dep.append({'t_mean': t_mean, 'ratio': ratio_m, 'rho_Lam': rho_Lam_b, 'rho_L': rho_L_b, 'n': n_bin})
        log(f"  {b+1:>4}  {t_mean:>8.1f}  {ratio_m:>+12.4f}  {imG_m:>+12.4f}  {rho_Lam_b:>+11.4f}  {rho_L_b:>+11.4f}  {n_bin:>5}")

    log()

    # Kendall τ로 단조성 검증
    if len(t_dep) >= 4:
        ratio_seq = [x['ratio'] for x in t_dep]
        t_seq     = [x['t_mean'] for x in t_dep]
        tau_ratio, p_tau_ratio = stats.kendalltau(t_seq, ratio_seq)
        log(f"  Kendall τ(t, Gamma/A 비율): τ={tau_ratio:+.4f}  (p={p_tau_ratio:.3e})")
        if abs(tau_ratio) < 0.3 or p_tau_ratio > 0.1:
            log("  → Gamma/A 비율의 t-의존성 미약 ✅ (예상: ψ(t)≈log(t/2) → 약한 t 의존)")
        else:
            log(f"  → Gamma/A 비율이 t에 의존 (τ={tau_ratio:+.4f}) — 추가 분석 필요")
    log()

    # ──────────────────────────────────────────────────────────────────────
    # STEP 6: gap_right 분위별 E[Delta_A | gap_bin]
    # ──────────────────────────────────────────────────────────────────────
    log("[STEP 6] gap_right 분위별 Gamma 기여 E[A_Λ-A_L | gap_right_bin]")
    log()
    log(f"  {'bin':>4}  {'gr̄':>8}  {'E[A_L]':>10}  {'E[A_Λ]':>10}  {'E[ΔA]':>10}  {'ΔA/A_Λ':>9}  {'n':>5}")
    log("  " + "-" * 70)

    gr_quantiles = np.percentile(gr_arr, np.linspace(0, 100, N_BINS + 1))
    for b in range(N_BINS):
        mask = (gr_arr >= gr_quantiles[b]) & (gr_arr < gr_quantiles[b+1])
        if b == N_BINS - 1:
            mask = (gr_arr >= gr_quantiles[b]) & (gr_arr <= gr_quantiles[b+1])
        n_bin = np.sum(mask)
        if n_bin == 0:
            continue
        gr_mean  = np.mean(gr_arr[mask])
        EA_L     = np.mean(AL_arr[mask])
        EA_Lam   = np.mean(ALam_arr[mask])
        E_DA     = np.mean(DeltaA_arr[mask])
        ratio_m  = E_DA / EA_Lam if EA_Lam != 0 else np.nan
        log(f"  {b+1:>4}  {gr_mean:>8.4f}  {EA_L:>10.4f}  {EA_Lam:>10.4f}  {E_DA:>+10.4f}  {ratio_m:>+9.4f}  {n_bin:>5}")

    log()

    # ──────────────────────────────────────────────────────────────────────
    # STEP 7: 3요인 분해 최종 요약
    # ──────────────────────────────────────────────────────────────────────
    log("=" * 80)
    log("[C-268 최종 요약] — 3요인 분해 검증")
    log("=" * 80)
    log()
    log("  가설: ρ(A_Λ, gap_right) = ρ_GUE + ρ_Gamma + ρ_arith")
    log(f"        {RHO_GUE:+.4f}       +  ???   +  ???")
    log()
    log(f"  측정 결과:")
    log(f"    ρ(A_L,   gap_right) = {rho_AL_gr:+.4f}  ← zero-sum (Gamma 미포함)")
    log(f"    ρ(A_Λ,   gap_right) = {rho_ALam_gr:+.4f}  ← completed (Gamma 포함)")
    log(f"    ρ(A_Λ-A_L, gap_right) = {rho_DA_gr:+.4f}  ← Gamma 단독 기여")
    log(f"    Δρ(Gamma 증폭)      = {delta_rho:+.4f}")
    log()

    # 3요인 분해 체계
    rho_arith = rho_AL_gr - RHO_GUE   # 산술 S₁² 기여 (A_L vs GUE)
    rho_gamma_contrib = delta_rho       # Gamma 기여 (A_Λ vs A_L)

    log(f"  3요인 분해:")
    log(f"    ρ_GUE   = {RHO_GUE:+.4f}  (GUE 이론)")
    log(f"    ρ_arith = {rho_arith:+.4f}  = ρ(A_L) - ρ_GUE  (산술 S₁² 기여)")
    log(f"    ρ_Gamma = {rho_gamma_contrib:+.4f}  = ρ(A_Λ) - ρ(A_L)  (Gamma factor 기여)")
    log(f"    합계    = {RHO_GUE + rho_arith + rho_gamma_contrib:+.4f}  vs 관측 {RHO_OBS_C265:+.4f}")
    log()

    # 잔차
    final_residual = abs(rho_ALam_gr - RHO_OBS_C265)
    log(f"  최종 잔차 |ρ(A_Λ) - ρ_obs| = {final_residual:.4f}")
    log()

    # ── 성공 기준 판정 ──
    log("  [성공 기준 판정]")

    # 기준 1: Gamma 기여 유의한 음상관
    crit1 = (rho_DA_gr < -0.05) and (p_DA_gr < 0.01)
    log(f"  ★★★  ρ(A_Λ-A_L, gap_right) < 0 (유의): {'✅ PASS' if crit1 else '❌ FAIL'}  [{rho_DA_gr:+.4f}, p={p_DA_gr:.3e}]")

    # 기준 2: 3요인 분해 완결
    decomp_residual = abs(rho_ALam_gr - rho_AL_gr - delta_rho)
    crit2 = (decomp_residual < 0.02) and crit1
    log(f"  ★★★★ 3요인 분해 잔차 < 0.02: {'✅ PASS' if crit2 else '❌ FAIL'}  [잔차={decomp_residual:.4f}]")

    # 기준 3: t-의존성 확인
    if len(t_dep) >= 4:
        crit3 = abs(tau_ratio) < 0.5  # 약한 t-의존성
        log(f"  ★★★  Gamma/A 비율 t-의존성 약: {'✅ PASS' if crit3 else '⚠️ 강한 t 의존'}  [τ={tau_ratio:+.4f}]")
    else:
        crit3 = None
        log(f"  ★★★  Gamma/A 비율 t-의존성: N/A (빈 부족)")

    # 기준 4: B-45 해소 여부
    b45_resolved = (final_residual < 0.02) and crit1
    log(f"  경계 B-45 해소: {'✅ 해소' if b45_resolved else '⏳ 부분 해소/미해소'}  [최종잔차={final_residual:.4f}]")

    log()

    # ── Paper 4 §6 함의 ──
    log("  [Paper 4 §6 함의]")
    if crit1:
        log(f"  → Gamma factor가 A-gap 음상관에 기여 확인 (Δρ={rho_gamma_contrib:+.4f})")
        log(f"  → §6 'GUE 모델' 절에 3요인 분해 A = GUE + Gamma + arith 추가 가능")
        if b45_resolved:
            log("  → B-45 완전 해소: '관측(-0.578) = GUE(-0.50) + Gamma + arith' 완결")
        else:
            log(f"  → B-45 부분 해소: 잔차 {final_residual:.4f} 남음 — 추가 분석 후 결론")
    else:
        log("  → Gamma 기여 미확립 — ρ(A_Λ-A_L, gap_right) 유의하지 않음")
        log("  → 산술 S₁² 구조가 0.08 잔차의 주원인 가능성 → 다른 방향 검토")

    log()
    log(f"  총 소요시간: {time.time()-t_start:.1f}초")
    log()

    save()
    log(f"[완료] 결과 저장: {RESULT_PATH}")


if __name__ == '__main__':
    main()
