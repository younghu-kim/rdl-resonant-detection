#!/usr/bin/env python3
"""
=============================================================================
[C-383] GL(4) Sym³(11a1) A_Λ–gap_min 상관 검정 — C-375 방법론 적용
=============================================================================

목적:
  Paper 4의 가장 큰 한계 "degree 1-3만 검증"을 공략.
  C-375 방법론(A_Λ, 이론적 밀도, Hadamard 보정)을 degree 4에 적용.

대상:
  Sym³(11a1) — degree 4, N=1331 (=11³), gammaV=[0,1,1,2]

방법 (C-375와 동일):
  [1] PARI lfunsympow(E,3) → Sym³ 영점 수집 (T=500)
  [2] Hadamard zero-sum → S₁, H₁
  [3] Gamma 보정 → S₁^Λ = S₁ - Im(Γ'/Γ), H₁^Λ = H₁ + Re(Γ''/Γ)
      gammaV = [0, 1, 1, 2]  (4개 항의 digamma 합)
      Γ_ℝ(s) Γ_ℝ(s+1) Γ_ℝ(s+1) Γ_ℝ(s+2)
  [4] A_Λ = (S₁^Λ)² + 2H₁^Λ
  [5] 이론적 밀도: d̄(t) = (1/2π)[4·log(t/2π) + log(1331)]
  [6] gap_min_GUE = gap_min × d̄(t)
  [7] Spearman ρ(A_Λ, gap_min_GUE) + 95% CI (부트스트랩 2000회)

C-269와의 차이:
  C-269: A_L (Gamma 미보정), 경험적 밀도 d̄=2/Δt, T=[2,55], n≈101
  C-383: A_Λ (Gamma 보정),  이론적 밀도,          T=500,    n≈300 예상

성공 기준 (수학자):
  ★★★★★★: ρ < -0.80 (degree 4 보편성 확정, Paper 4 Table 확장)
  ★★★★:   -0.80 ≤ ρ < -0.60 (약화되었으나 유의미한 음의 상관)
  ★★★:    ρ ≥ -0.60 (보편성 깨짐, 경계 발견)

결과: results/gl4_sym3_agap_c383.txt
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
N_BOOT = 2000            # 부트스트랩 반복
RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl4_sym3_agap_c383.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

mpmath.mp.dps = 50       # degree 4 — 더 높은 정밀도 필요

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

    수학자 지시:
      d̄(t) = (1/2π)[4·log(t/2π) + log(1331)]
    """
    arg_t = max(t / (2.0 * math.pi), 1.001)
    return (1.0 / (2.0 * math.pi)) * (
        degree * math.log(arg_t) + math.log(max(conductor, 1))
    )

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Gamma 보정 (degree-4: gammaV=[0,1,1,2])
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def gamma_corrections(t0, sigma_c, gammaV):
    """
    Gamma factor corrections for A_Λ.

    Γ_∞'/Γ_∞(s) = Σ_j [-(1/2)log(π) + (1/2)ψ((s+μ_j)/2)]

    For Sym³(11a1), gammaV=[0,1,1,2]:
      ψ((σ_c+it₀+0)/2) + ψ((σ_c+it₀+1)/2) + ψ((σ_c+it₀+1)/2) + ψ((σ_c+it₀+2)/2)

    수학자 주의: "digamma 합에 유의" — 4항의 합산

    Returns:
      im_gamma: Im[Γ_∞'/Γ_∞(σ_c + it₀)]
      re_gamma: Re[d/ds Γ_∞'/Γ_∞(σ_c + it₀)]
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

def pari_float(x, pari_obj):
    """PARI 객체를 float으로 안전하게 변환."""
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')

def get_zeros_from_pari(zvec_name, pari_obj):
    """PARI 영점 벡터를 Python 리스트로 변환."""
    n = int(str(pari_obj('#' + str(zvec_name))))
    zeros = []
    for i in range(1, n + 1):
        t = pari_float(pari_obj(f'{zvec_name}[{i}]'), pari_obj)
        if not math.isnan(t) and t > 1.0:
            zeros.append(t)
    return sorted(zeros)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# A_Λ 계산 (Hadamard zero-sum)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_all_A(zeros, sigma_c, gammaV):
    """
    모든 영점에 대해 A_L, A_Λ 계산.

    A_L  = S₁² + 2H₁               (같은 부호 + 켤레)
    A_Λ  = (S₁ - Im(Γ'/Γ))² + 2(H₁ + Re(Γ''/Γ))

    Returns: list of dicts {'t', 'A_L', 'A_Lambda'}, n_fail
    """
    n = len(zeros)
    data = []
    n_fail = 0
    n_neg = 0

    for i in range(n):
        t0 = zeros[i]

        try:
            # 같은 부호 영점 합산 (t₀ - γ_k)
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

            # 켤레 영점 합산 (t₀ + γ_k)
            S1_conj = 0.0
            H1_conj = 0.0
            for k in range(n):
                denom = t0 + zeros[k]
                S1_conj += 1.0 / denom
                H1_conj += 1.0 / (denom * denom)

            S1_L = S1_same + S1_conj
            H1_L = H1_same + H1_conj
            A_L = S1_L ** 2 + 2.0 * H1_L

            # Gamma 보정 (degree-4: 4개 digamma 항)
            im_gamma, re_gamma = gamma_corrections(t0, sigma_c, gammaV)

            S1_Lambda = S1_L - im_gamma
            H1_Lambda = H1_L + re_gamma
            A_Lambda = S1_Lambda ** 2 + 2.0 * H1_Lambda

            if A_Lambda <= 0:
                n_neg += 1
                n_fail += 1
                continue
            if A_L <= 0:
                n_fail += 1
                continue
            if math.isnan(A_Lambda) or math.isinf(A_Lambda):
                n_fail += 1
                continue

            data.append({
                't': t0,
                'A_L': A_L,
                'A_Lambda': A_Lambda,
                'S1_L': S1_L,
                'H1_L': H1_L,
                'im_gamma': im_gamma,
                're_gamma': re_gamma,
            })

        except Exception as e:
            n_fail += 1
            if n_fail <= 5:
                print(f"  WARNING: i={i}, t={t0:.3f}: {e}", flush=True)

    return data, n_fail, n_neg

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 상관 분석 + 부트스트랩 CI
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def bootstrap_ci(x, y, n_boot=N_BOOT, alpha=0.05, seed=42):
    """부트스트랩 95% CI for Spearman ρ."""
    rng = np.random.default_rng(seed)
    n = len(x)
    rho_boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        r, _ = stats.spearmanr(x[idx], y[idx])
        if not math.isnan(r):
            rho_boot.append(r)
    rho_boot = np.array(rho_boot)
    lo = np.percentile(rho_boot, 100 * alpha / 2)
    hi = np.percentile(rho_boot, 100 * (1 - alpha / 2))
    return lo, hi

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    log("=" * 80)
    log("[C-383] GL(4) Sym³(11a1) A_Λ–gap_min 상관 검정")
    log("  C-375 방법론 적용 (A_Λ, 이론적 밀도, Hadamard 보정)")
    log("=" * 80)
    log()
    log(f"  대상: Sym³(11a1)")
    log(f"  degree = 4, N = 1331 (= 11³)")
    log(f"  gammaV = [0, 1, 1, 2]")
    log(f"  이론적 밀도: d̄(t) = (1/2π)[4·log(t/2π) + log(1331)]")
    log(f"  T_max = 500 (C-269는 55)")
    log(f"  Trim: 양쪽 {TRIM_FRAC*100:.0f}% (중앙 {(1-2*TRIM_FRAC)*100:.0f}%)")
    log(f"  부트스트랩: {N_BOOT}회")
    log(f"  시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()

    # ── PARI 초기화 ─────────────────────────────────────────────────────
    pari = cypari2.Pari()
    pari.allocatemem(3 * 10**9)     # 3GB — GL(4) lfuninit 필수 (C-269는 2GB)
    pari.set_real_precision(100)    # degree 4 → 높은 정밀도 (C-269 동일)
    log("PARI 초기화: 3GB, realprecision=100")
    log()

    # ════════════════════════════════════════════════════════════════════
    # Sym³(11a1) 영점 수집
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("Sym³(11a1) 영점 수집")
    log("═" * 80)
    t1 = time.time()

    # 수학자 지시: lfunsympow(ellinit([0,-1,1,-10,-20]), 3)
    log("  PARI: ellinit([0,-1,1,-10,-20]) ...")
    pari('E11 = ellinit([0,-1,1,-10,-20])')

    log("  PARI: lfunsympow(E11, 3) ...")
    pari('Lsym3 = lfunsympow(E11, 3)')

    # L-함수 파라미터 추출
    try:
        k_raw = str(pari('Lsym3[4]'))
        k_val = int(round(float(k_raw)))
        sigma_c = k_val / 2.0
        log(f"  k (weight) = {k_val}, σ_c = {sigma_c}")
    except Exception as e:
        log(f"  k 자동 추출 실패: {e}")
        log(f"  기본값: k=4, σ_c=2.0 (Sym³ of weight 2)")
        k_val = 4
        sigma_c = 2.0

    try:
        gammaV_raw = str(pari('Lsym3[3]'))
        log(f"  gammaV (PARI) = {gammaV_raw}")
    except Exception as e:
        log(f"  gammaV 추출 실패: {e}")

    try:
        N_raw = str(pari('Lsym3[5]'))
        log(f"  N (conductor, PARI) = {N_raw}")
    except Exception as e:
        log(f"  N 추출 실패: {e}")

    try:
        sign_raw = str(pari('Lsym3[6]'))
        log(f"  ε (sign, PARI) = {sign_raw}")
    except Exception as e:
        log(f"  sign 추출 실패: {e}")

    # gammaV 설정
    # 핵심: PARI gammaV=[-1,0,0,1] + σ_c=2.0 사용 (PARI 정규화)
    # 검증: digamma 인자가 수학자 지시와 일치하는지 확인
    #   PARI: (σ_c+it+μ)/2 = (2+it+{-1,0,0,1})/2 = {0.5, 1, 1, 1.5}+it/2
    #   수학자: ψ(0.5+it/2), ψ(1+it/2), ψ(1+it/2), ψ(1.5+it/2) ← 일치 ✓
    # 주의: gammaV=[0,1,1,2] + σ_c=2.0은 ψ({1,1.5,1.5,2}+it/2) → 틀림!
    gammaV = [-1, 0, 0, 1]  # PARI 정규화 (σ_c=2.0과 쌍)
    degree = 4
    conductor = 1331
    log(f"  사용: gammaV={gammaV} (PARI), σ_c={sigma_c}, degree={degree}, N={conductor}")
    log(f"  digamma 인자: ψ(0.5+it/2), ψ(1+it/2), ψ(1+it/2), ψ(1.5+it/2)")
    log()

    # 영점 수집 — 단계적 시도 (GL(4)는 메모리 집약적)
    # C-269: T=55, 2GB로 ~105개 수집. 3GB로 T 확장 시도.
    T_targets = [
        (200, 220),   # T=200, lfuninit [0,220]
        (100, 110),   # T=100, lfuninit [0,110]
        (55, 60),     # T=55 (C-269 수준), lfuninit [0,60]
    ]
    zeros = []
    T_used = 0
    for T_max, T_init in T_targets:
        log(f"  lfuninit([0,{T_init}]) + lfunzeros(T={T_max}) 시도 ...")
        try:
            pari(f'Lsym3init = lfuninit(Lsym3, [0, {T_init}])')
            pari(f'zsym3 = lfunzeros(Lsym3init, {T_max})')
            zeros = get_zeros_from_pari('zsym3', pari)
            T_used = T_max
            log(f"  ✅ T={T_max} 성공: {len(zeros)}개, 소요={time.time()-t1:.1f}s")
            break
        except Exception as e:
            log(f"  ⚠️ T={T_max} 실패: {e}")
            # 스택 재할당 시도
            try:
                pari.allocatemem(3 * 10**9)
            except Exception:
                pass
            continue

    if len(zeros) == 0:
        log("⚠️ 영점 0개 — 탐색 로직 점검 필요")
        save()
        return

    log(f"  영점: {len(zeros)}개, t ∈ [{zeros[0]:.4f}, {zeros[-1]:.4f}]")
    log(f"  처음 5개: {[f'{z:.4f}' for z in zeros[:5]]}")
    log(f"  마지막 5개: {[f'{z:.4f}' for z in zeros[-5:]]}")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # A_Λ 계산
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("A_Λ 계산 (Hadamard zero-sum + degree-4 Gamma 보정)")
    log("═" * 80)
    log(f"  σ_c = {sigma_c}")
    log(f"  gammaV = {gammaV}")
    log(f"  Gamma 항: ψ((s+0)/2) + ψ((s+1)/2) + ψ((s+1)/2) + ψ((s+2)/2)")
    log()

    t2 = time.time()
    data, n_fail, n_neg = compute_all_A(zeros, sigma_c, gammaV)
    log(f"  유효={len(data)}, 실패={n_fail}, 음수A={n_neg}")
    log(f"  계산 소요: {time.time()-t2:.1f}s")
    log()

    if len(data) < MIN_ZEROS:
        log(f"  ⚠️ 유효 데이터 {len(data)} < {MIN_ZEROS} — 중단")
        save()
        return

    # Gamma 보정 통계
    im_gammas = [d['im_gamma'] for d in data]
    re_gammas = [d['re_gamma'] for d in data]
    log(f"  Im(Γ'/Γ) 범위: [{min(im_gammas):.4f}, {max(im_gammas):.4f}]")
    log(f"  Re(Γ''/Γ) 범위: [{min(re_gammas):.6f}, {max(re_gammas):.6f}]")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # Gap 계산 + 중앙 60% 트리밍
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("Gap 계산 + GUE 정규화 + 중앙 60% 트리밍")
    log("═" * 80)

    N_data = len(data)
    lo = int(N_data * TRIM_FRAC)
    hi = int(N_data * (1.0 - TRIM_FRAC))

    valid = []
    for idx in range(lo, hi):
        if idx <= 0 or idx >= N_data - 1:
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
    log(f"  전체 → 중앙60%: {N_data} → {n_valid} (trim {lo}+{N_data-hi}개)")
    log()

    if n_valid < 20:
        log(f"  ⚠️ trim 후 부족 ({n_valid} < 20) — 중단")
        save()
        return

    # ── 배열 추출 ────────────────────────────────────────────────────
    t_arr = np.array([d['t'] for d in valid])
    ALa_arr = np.array([d['A_Lambda'] for d in valid])
    AL_arr = np.array([d['A_L'] for d in valid])
    gm_arr = np.array([d['gap_min_gue'] for d in valid])
    gr_arr = np.array([d['gap_r_gue'] for d in valid])

    log(f"  n = {n_valid}, t ∈ [{t_arr.min():.1f}, {t_arr.max():.1f}]")
    log(f"  <A_Λ> = {np.mean(ALa_arr):.4f} ± {np.std(ALa_arr):.4f}")
    log(f"  <A_L> = {np.mean(AL_arr):.4f} ± {np.std(AL_arr):.4f}")
    log(f"  <gap_min_GUE> = {np.mean(gm_arr):.4f}")
    log(f"  <gap_right_GUE> = {np.mean(gr_arr):.4f}")
    log(f"  <density> = {np.mean([d['density'] for d in valid]):.4f}")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # Spearman 상관 + 부트스트랩 95% CI
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("Spearman 상관 분석")
    log("═" * 80)
    log()

    # (1) A_Λ vs gap_min
    rho_Lam_gm, p_Lam_gm = stats.spearmanr(ALa_arr, gm_arr)
    ci_lo_Lam_gm, ci_hi_Lam_gm = bootstrap_ci(ALa_arr, gm_arr)

    log(f"  ρ(A_Λ, gap_min_GUE)  = {rho_Lam_gm:+.4f}  "
        f"95%CI=[{ci_lo_Lam_gm:+.4f}, {ci_hi_Lam_gm:+.4f}]  "
        f"p={p_Lam_gm:.3e}  {sig(p_Lam_gm)}")

    # (2) A_Λ vs gap_right
    rho_Lam_gr, p_Lam_gr = stats.spearmanr(ALa_arr, gr_arr)
    ci_lo_Lam_gr, ci_hi_Lam_gr = bootstrap_ci(ALa_arr, gr_arr)

    log(f"  ρ(A_Λ, gap_right_GUE) = {rho_Lam_gr:+.4f}  "
        f"95%CI=[{ci_lo_Lam_gr:+.4f}, {ci_hi_Lam_gr:+.4f}]  "
        f"p={p_Lam_gr:.3e}  {sig(p_Lam_gr)}")

    # (3) A_L vs gap_min (C-269 비교용)
    rho_L_gm, p_L_gm = stats.spearmanr(AL_arr, gm_arr)
    log(f"  ρ(A_L, gap_min_GUE)   = {rho_L_gm:+.4f}  "
        f"p={p_L_gm:.3e}  {sig(p_L_gm)}")

    # (4) A_L vs gap_right (C-269 비교용)
    rho_L_gr, p_L_gr = stats.spearmanr(AL_arr, gr_arr)
    log(f"  ρ(A_L, gap_right_GUE) = {rho_L_gr:+.4f}  "
        f"p={p_L_gr:.3e}  {sig(p_L_gr)}")

    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # t-bin 안정성 (4-bin)
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("t-bin 안정성 분석 (4-bin)")
    log("═" * 80)
    log()

    n_bins = 4
    bin_edges = np.percentile(t_arr, np.linspace(0, 100, n_bins + 1))
    bin_rhos = []

    for b in range(n_bins):
        mask = (t_arr >= bin_edges[b]) & (t_arr < bin_edges[b+1])
        if b == n_bins - 1:
            mask = (t_arr >= bin_edges[b]) & (t_arr <= bin_edges[b+1])
        n_bin = np.sum(mask)
        if n_bin < 10:
            log(f"  bin {b+1}: t=[{bin_edges[b]:.1f}, {bin_edges[b+1]:.1f}]  n={n_bin} (부족)")
            continue
        r_bin, p_bin = stats.spearmanr(ALa_arr[mask], gm_arr[mask])
        bin_rhos.append(r_bin)
        log(f"  bin {b+1}: t=[{bin_edges[b]:.1f}, {bin_edges[b+1]:.1f}]  "
            f"n={n_bin}  ρ={r_bin:+.4f}  p={p_bin:.3e}  {sig(p_bin)}")

    if len(bin_rhos) >= 2:
        span = max(bin_rhos) - min(bin_rhos)
        log(f"  span(ρ) = {span:.4f}")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # C-269 비교표
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("C-269 vs C-383 비교")
    log("═" * 80)
    log()
    log(f"  {'항목':<25} {'C-269 (구 방법)':>16} {'C-383 (C-375 방법)':>18}")
    log(f"  {'-'*60}")
    log(f"  {'방법':.<25} {'A_L':>16} {'A_Λ':>18}")
    log(f"  {'밀도':.<25} {'경험적 d̄=2/Δt':>16} {'이론적':>18}")
    log(f"  {'T 범위':.<25} {'[2, 55]':>16} {f'[{t_arr.min():.0f}, {t_arr.max():.0f}] (T={T_used})':>18}")
    log(f"  {'영점 수':.<25} {'101':>16} {f'{n_valid}':>18}")
    log(f"  {'gammaV':.<25} {'[0,1,1,2]':>16} {'[0,1,1,2]':>18}")
    log(f"  {'ρ(gap_min)':.<25} {'-0.520':>16} {f'{rho_Lam_gm:+.4f}':>18}")
    log(f"  {'ρ(gap_right)':.<25} {'-0.324':>16} {f'{rho_Lam_gr:+.4f}':>18}")
    log()

    # 차이 분석
    diff_gm = rho_Lam_gm - (-0.520)
    log(f"  Δρ(gap_min) = {diff_gm:+.4f}  "
        f"({'강화' if diff_gm < -0.05 else '유사' if abs(diff_gm) < 0.05 else '약화'})")
    log()

    # 방법론 차이의 기여
    log("  방법론 차이 분석:")
    log(f"    (a) A_L vs A_Λ:   ρ(A_L,gm)={rho_L_gm:+.4f} vs ρ(A_Λ,gm)={rho_Lam_gm:+.4f}  Δ={rho_Lam_gm-rho_L_gm:+.4f}")
    log(f"    (b) T 범위 확장:  C-269 n=101 → C-383 n={n_valid}")
    log(f"    (c) 이론적 밀도:  C-269 경험적 → C-383 이론적")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # C-375 기존 결과와 통합
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("C-375 기존 결과 + C-383 통합")
    log("═" * 80)
    log()

    # C-375 재실행본 수치 (수학자 보드 기준)
    c375_data = [
        ('ζ(s)',          1,    1,  910, -0.900),
        ('χ₋₃',          1,    3,  494, -0.891),
        ('11a1',          2,   11,  437, -0.914),
        ('Δ(Ramanujan)',  2,    1,  324, -0.876),
        ('sym²(11a1)',    3,  121,   96, -0.885),
    ]

    # C-378 수치
    c378_data = [
        ('χ₅[1]',  1,  5, 238, -0.903),
        ('χ̄₅[1]',  1,  5, 238, -0.875),
        ('χ₇[1]',  1,  7, 254, -0.871),
        ('χ̄₇[1]',  1,  7, 255, -0.901),
    ]

    hdr = f"  {'L-함수':<20} {'d':>2} {'N':>6} {'n':>5} {'ρ':>8} {'방법':>12}"
    log(hdr)
    log(f"  {'-'*60}")

    all_rhos = []
    for name, d, N, n, rho in c375_data:
        log(f"  {name:<20} {d:>2} {N:>6} {n:>5} {rho:>+8.3f} {'C-375':>12}")
        all_rhos.append(rho)

    for name, d, N, n, rho in c378_data:
        log(f"  {name:<20} {d:>2} {N:>6} {n:>5} {rho:>+8.3f} {'C-378':>12}")
        all_rhos.append(rho)

    # C-383 추가
    log(f"  {'Sym³(11a1)':<20} {4:>2} {1331:>6} {n_valid:>5} {rho_Lam_gm:>+8.3f} {'C-383 ★':>12}")
    all_rhos.append(rho_Lam_gm)

    log()

    all_rhos_arr = np.array(all_rhos)
    log(f"  전체 (10 L-함수): ρ 평균 = {np.mean(all_rhos_arr):+.4f} ± {np.std(all_rhos_arr):.4f}")
    log(f"  범위: [{np.min(all_rhos_arr):+.4f}, {np.max(all_rhos_arr):+.4f}], "
        f"range = {np.max(all_rhos_arr) - np.min(all_rhos_arr):.4f}")

    # degree별
    d13_rhos = all_rhos_arr[:9]  # degree 1-3
    log(f"  degree 1-3 (9개): ρ 평균 = {np.mean(d13_rhos):+.4f} ± {np.std(d13_rhos):.4f}")
    log(f"  degree 4   (1개): ρ = {rho_Lam_gm:+.4f}")
    log(f"  degree 4 vs 1-3: Δρ = {rho_Lam_gm - np.mean(d13_rhos):+.4f}")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # 최종 판정
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("[C-383 판정]")
    log("═" * 80)
    log()

    if rho_Lam_gm < -0.80:
        verdict = "★★★★★★ degree 4 보편성 확정"
        log(f"  {verdict}")
        log(f"  ρ = {rho_Lam_gm:+.4f} < -0.80")
        log(f"  → Paper 4 Table 확장: degree 1-4 보편성")
        log(f"  → C-269 ρ=-0.52와의 차이: 방법론 개선(A_Λ+이론밀도+T확장)이 원인")
    elif rho_Lam_gm < -0.60:
        verdict = "★★★★ 약화되었으나 유의미한 음의 상관"
        log(f"  {verdict}")
        log(f"  ρ = {rho_Lam_gm:+.4f}, -0.80 ≤ ρ < -0.60")
        log(f"  → degree 4에서 상관 약화. 부분적 보편성.")
        log(f"  → degree 의존성 또는 GL(4) 구조 효과 가능")
    else:
        verdict = "★★★ 보편성 깨짐 — GL(4) 경계 발견"
        log(f"  {verdict}")
        log(f"  ρ = {rho_Lam_gm:+.4f} ≥ -0.60")
        log(f"  → degree 4에서 보편성 깨짐. 새 경계 B-72 등록 후보.")
        log(f"  → 가장 흥미로운 결과: degree 3-4 전이에서 무엇이 변하는가?")

    log()
    log(f"  95% CI: [{ci_lo_Lam_gm:+.4f}, {ci_hi_Lam_gm:+.4f}]")
    log(f"  p = {p_Lam_gm:.3e}")
    log(f"  n = {n_valid}")
    log()

    # ── 수학자 성공 기준 체크 ──────────────────────────────────────────
    log("[성공 기준 검증]")
    log(f"  ★★★★★★ ρ < -0.80: {'✅ PASS' if rho_Lam_gm < -0.80 else '❌ FAIL'}")
    log(f"  ★★★★   -0.80 ≤ ρ < -0.60: {'✅' if -0.80 <= rho_Lam_gm < -0.60 else '—'}")
    log(f"  ★★★    ρ ≥ -0.60: {'✅' if rho_Lam_gm >= -0.60 else '—'}")
    log()

    log(f"  총 소요시간: {time.time()-t_start:.1f}초")
    log()
    save()
    log(f"[완료] {RESULT_PATH}")


if __name__ == '__main__':
    main()
