#!/usr/bin/env python3
"""
=============================================================================
[C-375] 교차 가족 A_Λ-gap_min 보편성 검정 — Paper 4 Go/No-Go
=============================================================================

목적:
  5개 L-함수 가족에서 ρ(A_Λ, gap_min_GUE)의 보편성을 검증.
  Paper 4 핵심 가설: "A-gap 음의 상관은 L-함수 보편 법칙".

대상:
  1. ζ(s)             — degree 1, N=1   (기준선, C-270 재현)
  2. χ₋₃ (mod 3)      — degree 1, N=3
  3. 11a1 타원곡선     — degree 2, N=11
  4. Δ (Ramanujan tau) — degree 2, N=1, weight 12
  5. sym²(11a1)        — degree 3, N=121

방법:
  [1] PARI lfunzeros → 영점 수집
  [2] Hadamard zero-sum → S₁, H₁
  [3] Gamma 보정 → S₁^Λ = S₁ - Im(Γ'/Γ), H₁^Λ = H₁ + Re(Γ''/Γ)
  [4] A_Λ = (S₁^Λ)² + 2H₁^Λ
  [5] GUE 정규화: density(t) = (1/2π)[d·log(t/2π) + log(N)]
  [6] gap_min_GUE = gap_min × density(t)
  [7] Spearman ρ(A_Λ, gap_min_GUE) + 95% CI

성공 기준:
  ★★★★★ 5개 ρ < -0.70, 변동 < 0.20  → 보편성 확정
  ★★★★  3-4개 통과                    → 조건부 보편성
  ★★★   2개 이하 또는 파괴            → 보편성 깨짐 (가장 흥미로운 결과)

참조: gap_stability_c270.py (ζ 기준선), dirichlet_Abare_gap_c288.py (d=1 보편성)
주의:
  - GUE 정규화: 이론적 밀도 사용 (경험적 d̄ = 2/Δt 사용 금지)
  - A_Λ 사용 (A_L 아님): C-270에서 A_Λ+gap_min이 t-안정(span=0.16)
  - 중앙 60% 트리밍 (양쪽 20% 제외)

결과: results/cross_family_agap_c375.txt
=============================================================================
"""

import sys
import os
import time
import math

import numpy as np
from scipy import stats
import mpmath

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')

import cypari2

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

TRIM_FRAC = 0.20         # 양쪽 20% → 중앙 60%
MIN_ZEROS = 30           # 최소 영점 수
RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/cross_family_agap_c375.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

mpmath.mp.dps = 30       # ψ 계산용 (t ≤ 2000이면 충분)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 로깅
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

_buf = []

def log(msg=''):
    print(msg, flush=True)
    _buf.append(str(msg))

def save():
    with open(RESULT_PATH, 'w') as f:
        f.write('\n'.join(_buf))

sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 이론적 밀도 (GUE 정규화용)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def theoretical_density(t, degree, conductor):
    """
    L-함수 영점 밀도 (이론적):
      density(t) = (1/2π) [d·log(t/2π) + log(N)]

    유도: von Mangoldt 공식 + Stirling
    """
    arg_t = max(t / (2.0 * math.pi), 1.001)
    return (1.0 / (2.0 * math.pi)) * (
        degree * math.log(arg_t) + math.log(max(conductor, 1))
    )

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Gamma 보정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def gamma_corrections(t0, sigma_c, gammaV):
    """
    Gamma factor corrections for A_Λ.

    Γ_∞'/Γ_∞(s) = Σ_j [-(1/2)log(π) + (1/2)ψ((s+μ_j)/2)]

    Returns:
      im_gamma: Im[Γ_∞'/Γ_∞(σ_c + it₀)]  = Im[Σ_j (1/2)ψ((σ_c+it₀+μ_j)/2)]
      re_gamma: Re[d/ds Γ_∞'/Γ_∞(σ_c+it₀)] = Re[Σ_j (1/4)ψ₁((σ_c+it₀+μ_j)/2)]
    """
    s = mpmath.mpc(sigma_c, t0)
    im_sum = mpmath.mpf(0)
    re_sum = mpmath.mpf(0)
    for mu in gammaV:
        arg = (s + mu) / 2
        psi_val = mpmath.digamma(arg)
        psi1_val = mpmath.psi(1, arg)
        im_sum += mpmath.im(psi_val) / 2
        re_sum += mpmath.re(psi1_val) / 4
    return float(im_sum), float(re_sum)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PARI 유틸리티
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def pari_float(x):
    """PARI 객체를 float으로 안전하게 변환."""
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')

def get_zeros_from_pari(zvec):
    """PARI 영점 벡터를 Python 리스트로 변환."""
    n = int(str(pari('#' + str(zvec))))
    zeros = []
    for i in range(1, n + 1):
        t = pari_float(pari(f'{zvec}[{i}]'))
        if not math.isnan(t) and t > 1.0:
            zeros.append(t)
    return sorted(zeros)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# A_Λ 계산 (Hadamard zero-sum)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_all_A(zeros, sigma_c, gammaV):
    """
    모든 영점에 대해 A_L, A_Λ 계산.
    Returns: list of dicts {'t', 'A_L', 'A_Lambda'}
    """
    n = len(zeros)
    data = []
    n_fail = 0

    for i in range(n):
        t0 = zeros[i]

        try:
            # 같은 부호 영점 합산
            S1_same = 0.0
            H1_same = 0.0
            for k in range(n):
                if k == i:
                    continue
                d = t0 - zeros[k]
                if abs(d) < 1e-15:
                    continue
                S1_same += 1.0 / d
                H1_same += 1.0 / (d * d)

            # 켤레 영점 합산
            S1_conj = 0.0
            H1_conj = 0.0
            for k in range(n):
                denom = t0 + zeros[k]
                S1_conj += 1.0 / denom
                H1_conj += 1.0 / (denom * denom)

            S1_L = S1_same + S1_conj
            H1_L = H1_same + H1_conj
            A_L = S1_L ** 2 + 2.0 * H1_L

            # Gamma 보정
            im_gamma, re_gamma = gamma_corrections(t0, sigma_c, gammaV)

            S1_Lambda = S1_L - im_gamma
            H1_Lambda = H1_L + re_gamma
            A_Lambda = S1_Lambda ** 2 + 2.0 * H1_Lambda

            if A_Lambda <= 0 or A_L <= 0:
                n_fail += 1
                continue
            if math.isnan(A_Lambda) or math.isinf(A_Lambda):
                n_fail += 1
                continue

            data.append({
                't': t0,
                'A_L': A_L,
                'A_Lambda': A_Lambda,
            })

        except Exception as e:
            n_fail += 1
            if n_fail <= 3:
                print(f"  WARNING: i={i}, t={t0:.3f}: {e}", flush=True)

    return data, n_fail

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 단일 L-함수 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def analyze_family(name, zeros, degree, conductor, sigma_c, gammaV):
    """
    하나의 L-함수 가족에 대한 완전한 A_Λ-gap 분석.
    Returns: result dict or None
    """
    n_raw = len(zeros)
    log(f"  영점: {n_raw}개, t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")

    if n_raw < MIN_ZEROS:
        log(f"  ⚠️ 영점 {n_raw}개 < {MIN_ZEROS} — 건너뜀")
        return None

    # A_Λ 계산
    log(f"  A_Λ 계산 (σ_c={sigma_c}, gammaV={gammaV}) ...")
    t0 = time.time()
    data, n_fail = compute_all_A(zeros, sigma_c, gammaV)
    log(f"  유효={len(data)}, 실패={n_fail}, 소요={time.time()-t0:.1f}초")

    if len(data) < MIN_ZEROS:
        log(f"  ⚠️ 유효 데이터 부족 — 건너뜀")
        return None

    # Gap 계산 + 중앙 60% 트리밍
    N = len(data)
    lo = int(N * TRIM_FRAC)
    hi = int(N * (1.0 - TRIM_FRAC))

    valid = []
    for idx in range(lo, hi):
        if idx <= 0 or idx >= N - 1:
            continue
        d = data[idx]
        gap_r = data[idx + 1]['t'] - d['t']
        gap_l = d['t'] - data[idx - 1]['t']
        if gap_r <= 0 or gap_l <= 0:
            continue

        # 이론적 밀도로 GUE 정규화
        dens = theoretical_density(d['t'], degree, conductor)
        if dens <= 0:
            continue

        d['gap_r_gue'] = gap_r * dens
        d['gap_l_gue'] = gap_l * dens
        d['gap_min_gue'] = min(gap_r, gap_l) * dens
        d['density'] = dens
        valid.append(d)

    n_valid = len(valid)
    log(f"  중앙 60%: {N} → {n_valid} (trim {lo}+{N-hi}개)")

    if n_valid < 20:
        log(f"  ⚠️ trim 후 부족 — 건너뜀")
        return None

    # 배열 추출
    t_arr = np.array([d['t'] for d in valid])
    ALa_arr = np.array([d['A_Lambda'] for d in valid])
    AL_arr = np.array([d['A_L'] for d in valid])
    gm_arr = np.array([d['gap_min_gue'] for d in valid])

    # Spearman 상관
    rho_Lambda, p_Lambda = stats.spearmanr(ALa_arr, gm_arr)
    rho_L, p_L = stats.spearmanr(AL_arr, gm_arr)

    # 95% CI (Fisher z-transform)
    def fisher_ci(r, n, alpha=0.05):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            z = np.arctanh(r)
            se = 1.0 / math.sqrt(n - 3)
            z_crit = stats.norm.ppf(1 - alpha / 2)
            lo_z = z - z_crit * se
            hi_z = z + z_crit * se
            return np.tanh(lo_z), np.tanh(hi_z)

    ci_lo, ci_hi = fisher_ci(rho_Lambda, n_valid)

    log()
    log(f"  === {name} 결과 ===")
    log(f"  n = {n_valid}, t ∈ [{t_arr.min():.1f}, {t_arr.max():.1f}]")
    log(f"  ρ(A_Λ, gap_min_GUE) = {rho_Lambda:+.4f}  (p={p_Lambda:.3e})  "
        f"95%CI=[{ci_lo:+.4f}, {ci_hi:+.4f}]  {sig(p_Lambda)}")
    log(f"  ρ(A_L, gap_min_GUE)  = {rho_L:+.4f}  (p={p_L:.3e})  {sig(p_L)}")
    log(f"  <A_Λ> = {np.mean(ALa_arr):.4f}, <gap_min> = {np.mean(gm_arr):.4f}")
    log()

    return {
        'name': name,
        'degree': degree,
        'conductor': conductor,
        'n': n_valid,
        'n_raw': n_raw,
        't_range': (float(t_arr.min()), float(t_arr.max())),
        'rho_Lambda': rho_Lambda,
        'p_Lambda': p_Lambda,
        'ci_lo': ci_lo,
        'ci_hi': ci_hi,
        'rho_L': rho_L,
        'p_L': p_L,
        'mean_A_Lambda': float(np.mean(ALa_arr)),
        'mean_gap_min': float(np.mean(gm_arr)),
    }

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    global pari
    t_start = time.time()

    log("=" * 80)
    log("[C-375] 교차 가족 A_Λ-gap_min 보편성 검정")
    log("=" * 80)
    log(f"  대상: ζ(s), χ₋₃, 11a1, Δ(Ramanujan), sym²(11a1)")
    log(f"  방법: Hadamard zero-sum + Gamma 보정 + 이론적 밀도 GUE 정규화")
    log(f"  Trim: 양쪽 {TRIM_FRAC*100:.0f}% (중앙 {(1-2*TRIM_FRAC)*100:.0f}%)")
    log(f"  시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()

    # ── PARI 초기화 ─────────────────────────────────────────────────────
    pari = cypari2.Pari()
    pari.allocatemem(512 * 10**6)
    pari.set_real_precision(50)
    log("PARI 초기화: 512MB, realprecision=50")
    log()

    results = []

    # ════════════════════════════════════════════════════════════════════
    # (1) ζ(s) — degree 1, N=1
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("[1/5] ζ(s) — degree 1, N=1, gammaV=[0]")
    log("═" * 80)
    t1 = time.time()

    pari('Lz = lfuninit(lfuncreate(1), [0, 2100])')
    pari('zz = lfunzeros(Lz, 2000)')
    zeros_z = get_zeros_from_pari('zz')
    log(f"  PARI 영점 수집 완료 ({time.time()-t1:.1f}s)")

    r = analyze_family('ζ(s)', zeros_z, degree=1, conductor=1,
                       sigma_c=0.5, gammaV=[0])
    if r:
        results.append(r)
    save()

    # ════════════════════════════════════════════════════════════════════
    # (2) χ₋₃ (mod 3) — degree 1, N=3
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("[2/5] χ₋₃ (mod 3) — degree 1, N=3, gammaV=[1]")
    log("═" * 80)
    t1 = time.time()

    pari('Ld3 = lfuninit(lfuncreate(Mod(2,3)), [0, 1100])')
    pari('zd3 = lfunzeros(Ld3, 1000)')
    zeros_d3 = get_zeros_from_pari('zd3')
    log(f"  PARI 영점 수집 완료 ({time.time()-t1:.1f}s)")

    r = analyze_family('χ₋₃ (mod 3)', zeros_d3, degree=1, conductor=3,
                       sigma_c=0.5, gammaV=[1])
    if r:
        results.append(r)
    save()

    # ════════════════════════════════════════════════════════════════════
    # (3) 11a1 — degree 2, N=11
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("[3/5] 11a1 타원곡선 — degree 2, N=11, gammaV=[0,1]")
    log("═" * 80)
    t1 = time.time()

    pari('E11 = ellinit([0,-1,1,-10,-20])')
    pari('Le = lfuninit(lfuncreate(E11), [0, 600])')
    pari('ze = lfunzeros(Le, 500)')
    zeros_e = get_zeros_from_pari('ze')
    log(f"  PARI 영점 수집 완료 ({time.time()-t1:.1f}s)")

    r = analyze_family('11a1', zeros_e, degree=2, conductor=11,
                       sigma_c=1.0, gammaV=[0, 1])
    if r:
        results.append(r)
    save()

    # ════════════════════════════════════════════════════════════════════
    # (4) Ramanujan Δ — degree 2, N=1, weight 12
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("[4/5] Δ (Ramanujan tau) — degree 2, N=1, k=12, gammaV=[0,1]")
    log("═" * 80)
    t1 = time.time()

    log("  τ(n) 사전계산 (n=1..20000) ...")
    pari('Tvec = vector(20000, n, ramanujantau(n))')
    log(f"  τ(n) 완료 ({time.time()-t1:.1f}s)")

    pari('Ldelta = lfuncreate([Tvec, 0, [0,1], 12, 1, 1])')
    pari('Ldinit = lfuninit(Ldelta, [0, 600])')
    pari('zd = lfunzeros(Ldinit, 500)')
    zeros_d = get_zeros_from_pari('zd')
    log(f"  PARI 영점 수집 완료 ({time.time()-t1:.1f}s)")

    # σ_c = k/2 = 6 for Ramanujan Δ
    r = analyze_family('Δ (Ramanujan)', zeros_d, degree=2, conductor=1,
                       sigma_c=6.0, gammaV=[0, 1])
    if r:
        results.append(r)
    save()

    # ════════════════════════════════════════════════════════════════════
    # (5) sym²(11a1) — degree 3, N=121
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("[5/5] sym²(11a1) — degree 3, N=121")
    log("═" * 80)
    t1 = time.time()

    pari('Es2 = ellinit([0,-1,1,-10,-20])')
    pari('Ls2 = lfunsympow(Es2, 2)')

    # gammaV, k 자동 추출
    try:
        k_raw = str(pari('Ls2[4]'))
        k_val = int(round(float(k_raw)))
        sigma_c_s2 = k_val / 2.0
        log(f"  k = {k_val}, σ_c = {sigma_c_s2}")

        # gammaV 추출 시도
        gammaV_raw = str(pari('Ls2[3]'))
        log(f"  gammaV_raw = {gammaV_raw}")
    except Exception as e:
        log(f"  gammaV 자동 추출 실패: {e}")
        log(f"  기본값 사용: gammaV=[0,1,1], k=3, σ_c=1.5")
        k_val = 3
        sigma_c_s2 = 1.5

    # sym²(11a1): gammaV = [0, 1, 1] (standard for sym² of GL(2))
    gammaV_s2 = [0, 1, 1]
    log(f"  사용: gammaV={gammaV_s2}, σ_c={sigma_c_s2}")

    pari('Ls2init = lfuninit(Ls2, [0, 120])')
    pari('zs2 = lfunzeros(Ls2init, 100)')
    zeros_s2 = get_zeros_from_pari('zs2')
    log(f"  PARI 영점 수집 완료 ({time.time()-t1:.1f}s)")

    r = analyze_family('sym²(11a1)', zeros_s2, degree=3, conductor=121,
                       sigma_c=sigma_c_s2, gammaV=gammaV_s2)
    if r:
        results.append(r)
    save()

    # ════════════════════════════════════════════════════════════════════
    # 최종 비교표 + 판정
    # ════════════════════════════════════════════════════════════════════

    log()
    log("=" * 90)
    log("[C-375 비교표] ρ(A_Λ, gap_min_GUE) — 교차 가족 보편성")
    log("=" * 90)
    log()

    hdr = (f"{'L-함수':<20} {'d':>2} {'N':>5} {'n':>6} "
           f"{'ρ(A_Λ,gm)':>11} {'95%CI':>22} {'p':>12} {'ρ(A_L,gm)':>11}")
    log(hdr)
    log("-" * 95)

    for r in results:
        log(f"{r['name']:<20} {r['degree']:>2} {r['conductor']:>5} {r['n']:>6} "
            f"{r['rho_Lambda']:>+11.4f} [{r['ci_lo']:>+.4f},{r['ci_hi']:>+.4f}] "
            f"{r['p_Lambda']:>12.3e} {r['rho_L']:>+11.4f}")

    log()

    # ── 보편성 판정 ──────────────────────────────────────────────────────
    log("=" * 80)
    log("[C-375 판정]")
    log("=" * 80)
    log()

    n_pass = sum(1 for r in results if r['rho_Lambda'] < -0.70)
    n_total = len(results)

    rho_vals = [r['rho_Lambda'] for r in results]
    if len(rho_vals) >= 2:
        rho_range = max(rho_vals) - min(rho_vals)
        rho_mean = np.mean(rho_vals)
        rho_std = np.std(rho_vals)
    else:
        rho_range = float('nan')
        rho_mean = rho_vals[0] if rho_vals else float('nan')
        rho_std = float('nan')

    log(f"  통과 (ρ < -0.70): {n_pass}/{n_total}")
    log(f"  ρ 평균 ± 표준편차: {rho_mean:+.4f} ± {rho_std:.4f}")
    log(f"  ρ 범위 (max-min): {rho_range:.4f}")
    log()

    if n_pass == n_total and rho_range < 0.20:
        verdict = "★★★★★ 강한 보편성"
        log(f"  {verdict}")
        log(f"  모든 가족에서 ρ < -0.70, 변동 {rho_range:.4f} < 0.20")
        log(f"  → Paper 4 본격 착수 권고")
    elif n_pass == n_total and rho_range < 0.40:
        verdict = "★★★★ 보편성 (중등 변동)"
        log(f"  {verdict}")
        log(f"  모든 가족에서 ρ < -0.70이나 변동 {rho_range:.4f} > 0.20")
        log(f"  → Paper 4 진행 가능, degree/conductor 의존성 탐색 필요")
    elif n_pass >= 3:
        verdict = "★★★★ 조건부 보편성"
        log(f"  {verdict}")
        log(f"  {n_pass}/{n_total} 가족 통과")
        failed = [r['name'] for r in results if r['rho_Lambda'] >= -0.70]
        log(f"  실패: {', '.join(failed)}")
        log(f"  → 실패 가족 해부 필요")
    elif n_pass >= 1:
        verdict = "★★★ 부분 보편성"
        log(f"  {verdict}")
        log(f"  {n_pass}/{n_total} 가족만 통과")
        log(f"  → 보편성 가설 약화, degree/conductor 의존성 강함")
    else:
        verdict = "★★ 보편성 기각"
        log(f"  {verdict}")
        log(f"  ρ < -0.70 가족 없음")
        log(f"  → A-gap 보편성 가설 폐기")

    log()

    # ── degree별 분석 ─────────────────────────────────────────────────
    log("[degree별 분석]")
    for d in sorted(set(r['degree'] for r in results)):
        d_results = [r for r in results if r['degree'] == d]
        d_rhos = [r['rho_Lambda'] for r in d_results]
        names = [r['name'] for r in d_results]
        if len(d_rhos) >= 2:
            log(f"  degree {d}: ρ = {np.mean(d_rhos):+.4f} ± {np.std(d_rhos):.4f}  "
                f"({', '.join(names)})")
        else:
            log(f"  degree {d}: ρ = {d_rhos[0]:+.4f}  ({names[0]})")
    log()

    # ── C-270 비교 ────────────────────────────────────────────────────
    log("[C-270 참조 비교 (ζ)]")
    zeta_result = next((r for r in results if 'ζ' in r['name']), None)
    if zeta_result:
        log(f"  C-270 (T=2000, n=910): ρ(A_Λ, gap_min_GUE) = -0.8998")
        log(f"  C-375 (T=2000, n={zeta_result['n']}): ρ(A_Λ, gap_min_GUE) "
            f"= {zeta_result['rho_Lambda']:+.4f}")
        diff = abs(zeta_result['rho_Lambda'] - (-0.8998))
        log(f"  차이: {diff:.4f} {'✅ 재현' if diff < 0.05 else '⚠️ 불일치'}")
    log()

    # ── Paper 4 방향 ─────────────────────────────────────────────────
    log("[Paper 4 방향]")
    if n_pass == n_total:
        log("  보편성 확인 → Paper 4 주제: 'A-gap inequality across L-function families'")
        log("  다음: t-bin 안정성 교차 검증 (각 가족별 8-bin)")
    elif n_pass >= 3:
        log("  부분 보편성 → Paper 4 주제: 'A-gap correlation: universality and its breakdown'")
        log("  다음: 실패 가족의 메커니즘 규명")
    else:
        log("  보편성 기각 → Paper 4 불가. B-60 또는 Selberg 방향 전환.")
    log()

    log(f"  총 소요시간: {time.time()-t_start:.1f}초")
    log()
    save()
    log(f"[완료] {RESULT_PATH}")


if __name__ == '__main__':
    main()
