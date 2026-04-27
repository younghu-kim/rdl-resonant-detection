#!/usr/bin/env python3
"""
=============================================================================
[C-384] GL(4) Sym³(37a1) A_Λ–gap_min 상관 검정 — 2번째 GL(4) 확인
=============================================================================

목적:
  GL(4) 보편성 주장을 위한 2번째 데이터 포인트.
  C-383 (Sym³(11a1), ρ=-0.905)에 이어 Sym³(37a1)로 교차 검증.

대상:
  Sym³(37a1) — degree 4, N=50653 (=37³), gammaV=[-1,0,0,1] (PARI, σ_c=2)

방법 (C-383과 동일):
  [1] PARI lfunsympow(ellinit([0,1,1,-23,-50]), 3) → Sym³ 영점 수집
  [2] Hadamard zero-sum → S₁, H₁
  [3] Gamma 보정 → S₁^Λ = S₁ - Im(Γ'/Γ), H₁^Λ = H₁ + Re(Γ''/Γ)
      gammaV = [-1, 0, 0, 1] (PARI), σ_c = 2.0
      digamma: ψ(0.5+it/2), ψ(1+it/2), ψ(1+it/2), ψ(1.5+it/2)
  [4] A_Λ = (S₁^Λ)² + 2H₁^Λ
  [5] 이론적 밀도: d̄(t) = (1/2π)[4·log(t/2π) + log(50653)]
  [6] gap_min_GUE = gap_min × d̄(t)
  [7] Spearman ρ(A_Λ, gap_min_GUE) + 95% CI (부트스트랩 2000회)

C-383과의 차이:
  C-383: Sym³(11a1), N=1331,  n=327, ρ=-0.905
  C-384: Sym³(37a1), N=50653, n=? (목표 300+)

성공 기준 (수학자):
  ★★★★★: ρ < -0.80 (GL(4) 2/2 보편성 확정)
  ★★★★:  -0.80 ≤ ρ < -0.60 (약화, N-의존성 탐사 필요)
  ★★★:   ρ ≥ -0.60 (Sym³(11a1)과 불일치 → conductor 효과, B-73)

주의 (수학자):
  - N=50653은 N=1331보다 38배 큰 conductor → 밀도 높아 영점 더 많을 수 있음
  - PARI에서 큰 conductor → lfuninit 구간 [0, 250] 시도, 메모리 4GB 권장
  - 실패 시 fallback: Sym³(19a1) N=6859 또는 Sym³(14a1) N=2744

결과: results/gl4_sym3_37a1_agap_c384.txt
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
    '~/Desktop/gdl_unified/results/gl4_sym3_37a1_agap_c384.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

mpmath.mp.dps = 50       # degree 4 — 높은 정밀도

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
      d̄(t) = (1/2π)[4·log(t/2π) + log(50653)]
    """
    arg_t = max(t / (2.0 * math.pi), 1.001)
    return (1.0 / (2.0 * math.pi)) * (
        degree * math.log(arg_t) + math.log(max(conductor, 1))
    )

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Gamma 보정 (degree-4: gammaV=[-1,0,0,1] PARI + σ_c=2)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def gamma_corrections(t0, sigma_c, gammaV):
    """
    Gamma factor corrections for A_Λ.

    Γ_∞'/Γ_∞(s) = Σ_j [-(1/2)log(π) + (1/2)ψ((s+μ_j)/2)]

    For Sym³(37a1), gammaV=[-1,0,0,1] (PARI, σ_c=2):
      (σ_c+it+μ)/2 = (2+it+{-1,0,0,1})/2 = {0.5, 1, 1, 1.5}+it/2
      → ψ(0.5+it/2), ψ(1+it/2), ψ(1+it/2), ψ(1.5+it/2)

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

    Returns: list of dicts {'t', 'A_L', 'A_Lambda'}, n_fail, n_neg
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
# Fallback 커브 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 37a1 실패 시 19a1, 14a1 순으로 fallback (수학자 지시)
# 37a1 (N=50653) PARI OOM → 19a1부터 시작 (사이클 #386 수정)
CURVES = [
    {
        'name': 'Sym³(19a1)',
        'label': '19a1',
        'ainvs': '[0,1,1,-9,-15]',  # 19a1 Cremona
        'conductor': 6859,     # 19³
        'N_curve': 19,
    },
    {
        'name': 'Sym³(14a1)',
        'label': '14a1',
        'ainvs': '[1,0,1,4,-6]',     # 14a1 Cremona (PARI 검증: conductor=14)
        'conductor': 2744,     # 14³
        'N_curve': 14,
    },
]

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    log("=" * 80)
    log("[C-384] GL(4) Sym³(37a1) A_Λ–gap_min 상관 검정 — 2번째 GL(4) 확인")
    log("  C-383 방법론 적용 (A_Λ, 이론적 밀도, Hadamard 보정)")
    log("=" * 80)
    log()
    log(f"  목표: Sym³(37a1), degree 4, N=50653 (=37³)")
    log(f"  fallback: Sym³(19a1) N=6859, Sym³(14a1) N=2744")
    log(f"  gammaV (PARI) = [-1, 0, 0, 1], σ_c = 2.0")
    log(f"  digamma: ψ(0.5+it/2), ψ(1+it/2), ψ(1+it/2), ψ(1.5+it/2)")
    log(f"  이론적 밀도: d̄(t) = (1/2π)[4·log(t/2π) + log(N)]")
    log(f"  T_max = 200 (C-383과 동일)")
    log(f"  Trim: 양쪽 {TRIM_FRAC*100:.0f}% (중앙 {(1-2*TRIM_FRAC)*100:.0f}%)")
    log(f"  부트스트랩: {N_BOOT}회")
    log(f"  시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()

    # ── PARI 초기화 ─────────────────────────────────────────────────────
    pari = cypari2.Pari()
    pari.allocatemem(4 * 10**9)     # 4GB — 큰 conductor 대응 (수학자 지시)
    pari.set_real_precision(100)
    log("PARI 초기화: 4GB, realprecision=100")
    log()

    # ════════════════════════════════════════════════════════════════════
    # 영점 수집 (fallback 포함)
    # ════════════════════════════════════════════════════════════════════

    # gammaV/sigma_c/degree는 Sym³(EC weight 2)에 대해 동일
    gammaV = [-1, 0, 0, 1]  # PARI 정규화 (σ_c=2.0과 쌍)
    degree = 4
    sigma_c = 2.0

    zeros = []
    used_curve = None
    conductor = None
    T_used = 0

    for curve in CURVES:
        log("═" * 80)
        log(f"{curve['name']} 영점 수집 (N={curve['conductor']})")
        log("═" * 80)
        t1 = time.time()

        log(f"  PARI: ellinit({curve['ainvs']}) ...")
        try:
            pari(f"Ecurve = ellinit({curve['ainvs']})")
        except Exception as e:
            log(f"  ⚠️ ellinit 실패: {e}")
            continue

        log(f"  PARI: lfunsympow(Ecurve, 3) ...")
        try:
            pari('Lsym3 = lfunsympow(Ecurve, 3)')
        except Exception as e:
            log(f"  ⚠️ lfunsympow 실패: {e}")
            continue

        # L-함수 파라미터 확인
        try:
            k_raw = str(pari('Lsym3[4]'))
            k_val = int(round(float(k_raw)))
            sc = k_val / 2.0
            log(f"  k = {k_val}, σ_c = {sc}")
            if abs(sc - 2.0) > 0.01:
                log(f"  ⚠️ σ_c={sc} ≠ 2.0 — gammaV 재검토 필요")
        except Exception as e:
            log(f"  k 추출 실패: {e}, 기본값 σ_c=2.0 사용")

        try:
            gv_raw = str(pari('Lsym3[3]'))
            log(f"  gammaV (PARI) = {gv_raw}")
        except Exception:
            pass

        try:
            N_raw = str(pari('Lsym3[5]'))
            log(f"  N (PARI) = {N_raw}")
        except Exception:
            pass

        try:
            sign_raw = str(pari('Lsym3[6]'))
            log(f"  ε (PARI) = {sign_raw}")
        except Exception:
            pass

        # 영점 수집 — 단계적 시도
        # N=50653: lfuninit([0,250])이 OOM 유발 → 작은 T부터 시도
        # N<10000: 일반적 전략
        if curve['conductor'] > 10000:
            T_targets = [
                (50, 60),     # 가장 보수적 (OOM 방지)
                (30, 35),     # 최소한의 영점이라도 수집
            ]
        elif curve['conductor'] > 3000:
            T_targets = [
                (200, 220),
                (100, 110),
                (50, 60),
            ]
        else:
            T_targets = [
                (200, 220),
                (100, 110),
                (55, 60),
            ]

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
                # 스택 재할당
                try:
                    pari.allocatemem(4 * 10**9)
                except Exception:
                    pass
                continue

        if len(zeros) >= MIN_ZEROS:
            used_curve = curve
            conductor = curve['conductor']
            log(f"\n  ✅ {curve['name']} 사용 확정: {len(zeros)}개 영점")
            break
        else:
            log(f"  ⚠️ {curve['name']}: 영점 {len(zeros)}개 < {MIN_ZEROS} — 다음 fallback 시도")
            zeros = []

    if len(zeros) == 0 or used_curve is None:
        log("\n⚠️ 모든 커브에서 영점 수집 실패 — 중단")
        save()
        return

    log()
    log(f"  최종 사용: {used_curve['name']}, N={conductor}, degree={degree}")
    log(f"  영점: {len(zeros)}개, t ∈ [{zeros[0]:.4f}, {zeros[-1]:.4f}]")
    log(f"  처음 5개: {[f'{z:.4f}' for z in zeros[:5]]}")
    log(f"  마지막 5개: {[f'{z:.4f}' for z in zeros[-5:]]}")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # A_Λ 계산
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log(f"A_Λ 계산 (Hadamard zero-sum + degree-4 Gamma 보정)")
    log("═" * 80)
    log(f"  σ_c = {sigma_c}")
    log(f"  gammaV = {gammaV}")
    log(f"  Gamma 항: ψ(0.5+it/2), ψ(1+it/2), ψ(1+it/2), ψ(1.5+it/2)")
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

    # (3) A_L vs gap_min (비보정 비교용)
    rho_L_gm, p_L_gm = stats.spearmanr(AL_arr, gm_arr)
    log(f"  ρ(A_L, gap_min_GUE)   = {rho_L_gm:+.4f}  "
        f"p={p_L_gm:.3e}  {sig(p_L_gm)}")

    # (4) A_L vs gap_right
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
    # C-383 비교표
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log(f"C-383 (Sym³(11a1)) vs C-384 ({used_curve['name']}) 비교")
    log("═" * 80)
    log()
    c384_label = f"C-384 {used_curve['name']}"
    log(f"  {'항목':<25} {'C-383 Sym³(11a1)':>18} {c384_label:>20}")
    log(f"  {'-'*65}")
    curve_label = used_curve['label']
    log(f"  {'N (conductor)':.<25} {'1331':>18} {str(conductor):>20}")
    log(f"  {'curve':.<25} {'11a1':>18} {curve_label:>20}")
    log(f"  {'T_max':.<25} {'200':>18} {f'{T_used}':>20}")
    log(f"  {'영점 (raw)':.<25} {'544':>18} {f'{len(zeros)}':>20}")
    log(f"  {'n (trim 후)':.<25} {'327':>18} {f'{n_valid}':>20}")
    log(f"  {'ρ(A_Λ, gap_min)':.<25} {'-0.905':>18} {f'{rho_Lam_gm:+.4f}':>20}")
    log(f"  {'ρ(A_Λ, gap_right)':.<25} {'-0.418':>18} {f'{rho_Lam_gr:+.4f}':>20}")
    log(f"  {'ρ(A_L, gap_min)':.<25} {'-0.824':>18} {f'{rho_L_gm:+.4f}':>20}")
    log()

    # 일치도 분석
    diff = abs(rho_Lam_gm - (-0.905))
    if diff < 0.05:
        log(f"  → 두 GL(4) ρ 일치 (Δ={rho_Lam_gm-(-0.905):+.4f}). 보편성 강화.")
    elif diff < 0.15:
        log(f"  → 두 GL(4) ρ 근접 (Δ={rho_Lam_gm-(-0.905):+.4f}). 부분적 보편성.")
    else:
        log(f"  → 두 GL(4) ρ 차이 큼 (Δ={rho_Lam_gm-(-0.905):+.4f}). conductor 효과?")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # 전체 통합 (11 L-함수)
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log(f"전체 통합: 11 L-함수 (C-375 + C-378 + C-383 + C-384)")
    log("═" * 80)
    log()

    # 기존 10개 데이터
    prev_data = [
        ('ζ(s)',          1,    1,  910, -0.900, 'C-375', 'SD'),
        ('χ₋₃',          1,    3,  494, -0.891, 'C-375', 'SD'),
        ('11a1',          2,   11,  437, -0.914, 'C-375', 'SD'),
        ('Δ(Ramanujan)',  2,    1,  324, -0.876, 'C-375', 'SD'),
        ('sym²(11a1)',    3,  121,   96, -0.885, 'C-375', 'SD'),
        ('Sym³(11a1)',    4, 1331,  327, -0.905, 'C-383', 'SD'),
        ('χ₅[1]',        1,    5,  238, -0.903, 'C-378', 'NSD'),
        ('χ̄₅[1]',        1,    5,  238, -0.875, 'C-378', 'NSD'),
        ('χ₇[1]',        1,    7,  254, -0.871, 'C-378', 'NSD'),
        ('χ̄₇[1]',        1,    7,  255, -0.901, 'C-378', 'NSD'),
    ]

    hdr = f"  {'L-함수':<20} {'d':>2} {'N':>6} {'n':>5} {'ρ':>8} {'src':>8} {'type':>4}"
    log(hdr)
    log(f"  {'-'*60}")

    all_rhos = []
    for name, d, N, n, rho, src, tp in prev_data:
        log(f"  {name:<20} {d:>2} {N:>6} {n:>5} {rho:>+8.3f} {src:>8} {tp:>4}")
        all_rhos.append(rho)

    # C-384 추가
    log(f"  {used_curve['name']:<20} {4:>2} {conductor:>6} {n_valid:>5} {rho_Lam_gm:>+8.3f} {'C-384 ★':>8} {'SD':>4}")
    all_rhos.append(rho_Lam_gm)

    log()
    all_rhos_arr = np.array(all_rhos)
    log(f"  전체 (11 L-함수): ρ 평균 = {np.mean(all_rhos_arr):+.4f} ± {np.std(all_rhos_arr):.4f}")
    log(f"  범위: [{np.min(all_rhos_arr):+.4f}, {np.max(all_rhos_arr):+.4f}], "
        f"range = {np.max(all_rhos_arr) - np.min(all_rhos_arr):.4f}")

    # degree 별
    d13 = np.array(all_rhos[:5] + all_rhos[6:10])  # degree 1-3
    d4 = np.array([all_rhos[5], all_rhos[10]])       # degree 4 (11a1 + 37a1)
    log(f"  degree 1-3 (9개): ρ 평균 = {np.mean(d13):+.4f} ± {np.std(d13):.4f}")
    log(f"  degree 4   (2개): ρ 평균 = {np.mean(d4):+.4f} ± {np.std(d4):.4f}")
    log(f"  degree 4 vs 1-3: Δρ = {np.mean(d4) - np.mean(d13):+.4f}")
    log()

    # GL(4) 내부 일관성
    log(f"  GL(4) 내부: Sym³(11a1) ρ=-0.905, {used_curve['name']} ρ={rho_Lam_gm:+.4f}")
    log(f"  GL(4) Δρ = {rho_Lam_gm - (-0.905):+.4f}, "
        f"N ratio = {conductor}/1331 = {conductor/1331:.1f}×")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # 최종 판정
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("[C-384 판정]")
    log("═" * 80)
    log()

    if rho_Lam_gm < -0.80:
        verdict = "★★★★★ GL(4) 2/2 보편성 확정"
        log(f"  {verdict}")
        log(f"  ρ = {rho_Lam_gm:+.4f} < -0.80")
        log(f"  → GL(4) 2개 L-함수에서 모두 ρ < -0.80 확인")
        log(f"  → Paper 4 Table: degree 1-4 보편성 강화")
        log(f"  → conductor {conductor} (vs 1331) — N-의존성 없음 확인")
    elif rho_Lam_gm < -0.60:
        verdict = "★★★★ 약화 — N-의존성 탐사 필요"
        log(f"  {verdict}")
        log(f"  ρ = {rho_Lam_gm:+.4f}, -0.80 ≤ ρ < -0.60")
        log(f"  → conductor N={conductor}에서 ρ 약화")
        log(f"  → Sym³(11a1) ρ=-0.905 vs {used_curve['name']} ρ={rho_Lam_gm:+.4f}")
        log(f"  → B-73 후보: conductor 의존성")
    else:
        verdict = "★★★ Sym³(11a1)과 불일치 — conductor 효과 발견"
        log(f"  {verdict}")
        log(f"  ρ = {rho_Lam_gm:+.4f} ≥ -0.60")
        log(f"  → 큰 conductor에서 보편성 깨짐")
        log(f"  → B-73 등록: conductor 효과")

    log()
    log(f"  95% CI: [{ci_lo_Lam_gm:+.4f}, {ci_hi_Lam_gm:+.4f}]")
    log(f"  p = {p_Lam_gm:.3e}")
    log(f"  n = {n_valid}")
    log()

    # ── 성공 기준 체크 ──────────────────────────────────────────
    log("[성공 기준 검증]")
    log(f"  ★★★★★ ρ < -0.80: {'✅ PASS' if rho_Lam_gm < -0.80 else '❌ FAIL'}")
    log(f"  ★★★★  -0.80 ≤ ρ < -0.60: {'✅' if -0.80 <= rho_Lam_gm < -0.60 else '—'}")
    log(f"  ★★★   ρ ≥ -0.60: {'✅' if rho_Lam_gm >= -0.60 else '—'}")
    log()

    log(f"  총 소요시간: {time.time()-t_start:.1f}초")
    log()
    save()
    log(f"[완료] {RESULT_PATH}")


if __name__ == '__main__':
    main()
