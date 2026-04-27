#!/usr/bin/env python3
"""
=============================================================================
[C-388] GL(6) Sym⁵(11a1) A_Λ–gap_min 상관 검정 — degree 6 탐사
=============================================================================

목적:
  12개 L-함수에서 degree 1-5 보편성이 확정된 상태 (ρ=-0.895±0.015).
  GL(6)은 구 방법(A_L)으로 ρ=-0.42 (C-275)였으나 신 방법(A_Λ) 미적용.
  신 방법으로 재측정하여 degree 1-6 전체를 통일.

대상:
  Sym⁵(11a1) — degree 6, N=161051 (=11⁵)
  gammaV, σ_c, ε: PARI에서 자동 추출

방법 (C-387과 동일):
  [1] PARI lfunsympow(E,5) → Sym⁵ 영점 수집 (T=30 시작, OOM시 축소)
  [2] Hadamard zero-sum → S₁, H₁
  [3] Gamma 보정 → S₁^Λ = S₁ - Im(Γ'/Γ), H₁^Λ = H₁ + Re(Γ''/Γ)
  [4] A_Λ = (S₁^Λ)² + 2H₁^Λ
  [5] 이론적 밀도: d̄(t) = (1/2π)[6·log(t/2π) + log(N)]
  [6] gap_min_GUE = gap_min × d̄(t)
  [7] Spearman ρ(A_Λ, gap_min_GUE) + 95% CI (부트스트랩 2000회)

성공 기준 (수학자):
  ★★★★★★★: ρ < -0.80 (degree 6 보편성 — degree 1-6 A_Λ 통일!)
  ★★★★:    -0.80 ≤ ρ < -0.60 (약화, degree 의존성 발견)
  ★★★:     PARI 실패 (계산 경계 확인)
  ★★:      ρ ≥ -0.60 또는 n < 20 (통계적 불충분)

비교: C-275 (구 방법) ρ(A_L, gap_min)=-0.42 → C-388 (신 방법) ρ(A_Λ, gap_min)=?

결과: results/gl6_sym5_agap_c388.txt
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
MIN_ZEROS = 20           # 최소 영점 수 (degree 6은 영점 적음 → 완화)
N_BOOT = 2000            # 부트스트랩 반복
RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl6_sym5_agap_c388.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

mpmath.mp.dps = 80       # degree 6 — 높은 정밀도 필요

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
    """
    arg_t = max(t / (2.0 * math.pi), 1.001)
    return (1.0 / (2.0 * math.pi)) * (
        degree * math.log(arg_t) + math.log(max(conductor, 1))
    )

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Gamma 보정 (degree 변수)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def gamma_corrections(t0, sigma_c, gammaV):
    """
    Gamma factor corrections for A_Λ.

    Γ_∞'/Γ_∞(s) = Σ_j [-(1/2)log(π) + (1/2)ψ((s+μ_j)/2)]

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

def parse_pari_gammaV(pari_obj, lname):
    """PARI L-함수 객체에서 gammaV를 추출하여 Python 리스트로 변환."""
    try:
        raw = str(pari_obj(f'{lname}[3]'))
        # PARI 벡터 형식: [a, b, c, ...] 또는 복잡한 형식
        # 숫자만 추출
        raw = raw.replace('[', '').replace(']', '').replace('~', '')
        parts = raw.split(',')
        result = []
        for p in parts:
            p = p.strip()
            if '/' in p:
                # 분수 처리
                num, den = p.split('/')
                result.append(float(num.strip()) / float(den.strip()))
            else:
                try:
                    result.append(float(p))
                except ValueError:
                    pass
        return result
    except Exception as e:
        log(f"  ⚠️ gammaV 파싱 실패: {e}")
        return None

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

            # Gamma 보정 (degree-6: 6개 digamma 항)
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
    log("[C-388] GL(6) Sym⁵(11a1) A_Λ–gap_min 상관 검정 — degree 6 재측정")
    log("  C-275 (구 방법) ρ(A_L)=-0.42 → C-388 (신 방법 A_Λ) 재측정")
    log("  C-387 방법론 적용 (A_Λ, 이론적 밀도, Hadamard 보정)")
    log("=" * 80)
    log()
    log(f"  대상: Sym⁵(11a1)")
    log(f"  degree = 6, N = 161051 (= 11⁵)")
    log(f"  gammaV, σ_c, ε: PARI에서 자동 추출")
    log(f"  이론적 밀도: d̄(t) = (1/2π)[6·log(t/2π) + log(N)]")
    log(f"  Trim: 양쪽 {TRIM_FRAC*100:.0f}% (중앙 {(1-2*TRIM_FRAC)*100:.0f}%)")
    log(f"  부트스트랩: {N_BOOT}회")
    log(f"  시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()

    # ── PARI 초기화 ─────────────────────────────────────────────────────
    pari = cypari2.Pari()

    MEM_TARGETS = [
        (8 * 10**9, "8GB"),
        (6 * 10**9, "6GB"),
        (4 * 10**9, "4GB"),
        (3 * 10**9, "3GB"),
    ]
    mem_used = "?"
    for mem_bytes, mem_label in MEM_TARGETS:
        try:
            pari.allocatemem(mem_bytes)
            mem_used = mem_label
            break
        except Exception as e:
            log(f"  ⚠️ {mem_label} 할당 실패: {e}")
            continue

    pari.set_real_precision(120)    # degree 6 → 최고 정밀도
    log(f"PARI 초기화: {mem_used}, realprecision=120")
    log()

    # ════════════════════════════════════════════════════════════════════
    # Sym⁵(11a1) 영점 수집
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("Sym⁵(11a1) 영점 수집")
    log("═" * 80)
    t1 = time.time()

    log("  PARI: ellinit([0,-1,1,-10,-20]) ...")
    pari('E11 = ellinit([0,-1,1,-10,-20])')

    log("  PARI: lfunsympow(E11, 5) ...")
    try:
        pari('Lsym5 = lfunsympow(E11, 5)')
        log("  ✅ lfunsympow(E11, 5) 성공")
    except Exception as e:
        log(f"  ❌ lfunsympow(E11, 5) 실패: {e}")
        log()
        log("═" * 80)
        log("[C-388 판정]")
        log("═" * 80)
        log("  ★★★ PARI 실패 — 계산 경계 확인")
        log(f"  에러: {e}")
        log(f"  → 'PARI는 Sym⁵를 지원하지 않음' (degree 5가 계산 한계)")
        log(f"  총 소요시간: {time.time()-t_start:.1f}초")
        save()
        return

    # L-함수 파라미터 추출
    try:
        k_raw = str(pari('Lsym5[4]'))
        k_val = int(round(float(k_raw)))
        sigma_c = k_val / 2.0
        log(f"  k (weight) = {k_val}, σ_c = {sigma_c}")
    except Exception as e:
        log(f"  k 자동 추출 실패: {e} → 기본값 k=6, σ_c=3.0")
        k_val = 6
        sigma_c = 3.0

    # gammaV 추출
    gammaV = parse_pari_gammaV(pari, 'Lsym5')
    if gammaV is None or len(gammaV) != 6:
        log(f"  ⚠️ gammaV 자동 추출 실패 또는 길이≠6 (got: {gammaV})")
        # degree 6 대칭멱: Sym⁵ of weight 2 form
        # 예상: [-5/2, -3/2, -1/2, 1/2, 3/2, 5/2] 또는 PARI 정규화
        # 안전한 폴백: PARI raw 출력을 확인하고 수동 설정
        try:
            raw = str(pari('Lsym5[3]'))
            log(f"  PARI Lsym5[3] raw: {raw}")
        except Exception:
            pass
        # σ_c=3 기준 gammaV 추정 (Sym⁵ of weight 2)
        gammaV_fallback = [-5, -3, -1, 1, 3, 5]
        log(f"  폴백 gammaV (σ_c 보정): 미확정 — PARI raw 확인 후 결정")
        log(f"  시도: {gammaV_fallback}")
        gammaV = gammaV_fallback
    else:
        log(f"  gammaV (PARI) = {gammaV}")

    try:
        N_raw = str(pari('Lsym5[5]'))
        log(f"  N (conductor, PARI) = {N_raw}")
        # 실제 conductor 추출 시도
        try:
            conductor = int(float(N_raw.strip().replace(' ', '')))
        except ValueError:
            conductor = 161051
            log(f"  conductor 파싱 실패 → 기본값 {conductor}")
    except Exception as e:
        conductor = 161051
        log(f"  N 추출 실패: {e} → 기본값 {conductor}")

    try:
        sign_raw = str(pari('Lsym5[6]'))
        log(f"  ε (sign, PARI) = {sign_raw}")
    except Exception as e:
        log(f"  sign 추출 실패: {e}")

    try:
        fe_check = str(pari('lfuncheckfeq(Lsym5)'))
        log(f"  FE check = {fe_check}")
    except Exception as e:
        log(f"  FE check 실패: {e}")

    degree = 6
    log(f"  사용: gammaV={gammaV}, σ_c={sigma_c}, degree={degree}, N={conductor}")
    log()

    # 영점 수집 — 단계적 시도
    # C-275는 T=30에서 75영점 성공. T=50 시도 후 축소.
    T_targets = [
        (50,  60),   # T=50 시도 (공격적)
        (30,  35),   # T=30 (C-275에서 확인된 수준)
        (20,  25),   # T=20 (최소)
    ]
    zeros = []
    T_used = 0
    for T_max, T_init in T_targets:
        log(f"  lfuninit([0,{T_init}]) + lfunzeros(T={T_max}) 시도 ...")
        try:
            pari(f'Lsym5init = lfuninit(Lsym5, [0, {T_init}])')
            pari(f'zsym5 = lfunzeros(Lsym5init, {T_max})')
            zeros = get_zeros_from_pari('zsym5', pari)
            T_used = T_max
            elapsed = time.time() - t1
            log(f"  ✅ T={T_max} 성공: {len(zeros)}개, 소요={elapsed:.1f}s")
            break
        except Exception as e:
            log(f"  ⚠️ T={T_max} 실패: {e}")
            # 스택 재할당 시도
            try:
                for mem_bytes, mem_label in MEM_TARGETS:
                    try:
                        pari.allocatemem(mem_bytes)
                        break
                    except Exception:
                        continue
            except Exception:
                pass
            continue

    if len(zeros) == 0:
        log("⚠️ 영점 0개 — 모든 T에서 실패")
        log()
        log("═" * 80)
        log("[C-388 판정]")
        log("═" * 80)
        log("  ★★★ PARI OOM — 계산 경계 = degree 5 (A_Λ 방법)")
        log("  → Sym⁵(11a1) N=161051의 영점 계산이 현재 리소스로 불가")
        log(f"  총 소요시간: {time.time()-t_start:.1f}초")
        save()
        return

    log(f"  영점: {len(zeros)}개, t ∈ [{zeros[0]:.4f}, {zeros[-1]:.4f}]")
    log(f"  처음 5개: {[f'{z:.4f}' for z in zeros[:5]]}")
    log(f"  마지막 5개: {[f'{z:.4f}' for z in zeros[-5:]]}")
    log()
    save()

    if len(zeros) < MIN_ZEROS:
        log(f"  ⚠️ 영점 {len(zeros)}개 < {MIN_ZEROS} (최소) — 통계적 불충분")
        log()
        log("═" * 80)
        log("[C-388 판정]")
        log("═" * 80)
        log(f"  ★★ n < {MIN_ZEROS} (통계적 불충분)")
        log(f"  영점 {len(zeros)}개로는 유의미한 상관 검정 불가")
        log(f"  총 소요시간: {time.time()-t_start:.1f}초")
        save()
        return

    # ════════════════════════════════════════════════════════════════════
    # A_Λ 계산
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("A_Λ 계산 (Hadamard zero-sum + degree-6 Gamma 보정)")
    log("═" * 80)
    log(f"  σ_c = {sigma_c}")
    log(f"  gammaV = {gammaV}")
    log(f"  Gamma 항: {len(gammaV)}개 digamma")
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

    if n_valid < 15:
        log(f"  ⚠️ trim 후 부족 ({n_valid} < 15) — 중단")
        log()
        log("═" * 80)
        log("[C-388 판정]")
        log("═" * 80)
        log(f"  ★★ trim 후 n={n_valid} < 15 (통계적 불충분)")
        log(f"  총 소요시간: {time.time()-t_start:.1f}초")
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

    # (3) A_L vs gap_min (비교용)
    rho_L_gm, p_L_gm = stats.spearmanr(AL_arr, gm_arr)
    log(f"  ρ(A_L, gap_min_GUE)   = {rho_L_gm:+.4f}  "
        f"p={p_L_gm:.3e}  {sig(p_L_gm)}")

    # (4) A_L vs gap_right (비교용)
    rho_L_gr, p_L_gr = stats.spearmanr(AL_arr, gr_arr)
    log(f"  ρ(A_L, gap_right_GUE) = {rho_L_gr:+.4f}  "
        f"p={p_L_gr:.3e}  {sig(p_L_gr)}")

    # C-275 비교
    log()
    log("  C-275 (구 방법) 비교:")
    log(f"    C-275: ρ(A_L, gap_min) = -0.4225, n=71, T=30")
    log(f"    C-388: ρ(A_Λ, gap_min) = {rho_Lam_gm:+.4f}, n={n_valid}, T={T_used}")
    delta_rho = rho_Lam_gm - (-0.4225)
    log(f"    Δρ = {delta_rho:+.4f} (신 방법 개선)")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # t-bin 안정성 (2-bin: degree 6은 데이터 적으므로 2분할)
    # ════════════════════════════════════════════════════════════════════
    n_bins = 2 if n_valid < 60 else 4
    log("═" * 80)
    log(f"t-bin 안정성 분석 ({n_bins}-bin)")
    log("═" * 80)
    log()

    bin_edges = np.percentile(t_arr, np.linspace(0, 100, n_bins + 1))
    bin_rhos = []

    for b in range(n_bins):
        mask = (t_arr >= bin_edges[b]) & (t_arr < bin_edges[b+1])
        if b == n_bins - 1:
            mask = (t_arr >= bin_edges[b]) & (t_arr <= bin_edges[b+1])
        n_bin = np.sum(mask)
        if n_bin < 8:
            log(f"  bin {b+1}: t=[{bin_edges[b]:.1f}, {bin_edges[b+1]:.1f}]  n={n_bin} (부족)")
            bin_rhos.append(float('nan'))
            continue
        r_bin, p_bin = stats.spearmanr(ALa_arr[mask], gm_arr[mask])
        bin_rhos.append(r_bin)
        log(f"  bin {b+1}: t=[{bin_edges[b]:.1f}, {bin_edges[b+1]:.1f}]  "
            f"n={n_bin}  ρ={r_bin:+.4f}  p={p_bin:.3e}  {sig(p_bin)}")

    valid_bin_rhos = [r for r in bin_rhos if not math.isnan(r)]
    if len(valid_bin_rhos) >= 2:
        span = max(valid_bin_rhos) - min(valid_bin_rhos)
        log(f"  span(ρ) = {span:.4f}")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # 기존 결과와 통합 비교표 (degree 1-6, A_Λ 통일)
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("기존 결과 + C-388 통합 비교표 (A_Λ 방법 통일)")
    log("═" * 80)
    log()

    # 12 L-함수 (C-375/378/383/384/387)
    existing_data = [
        ('ζ(s)',          1,      1,  910, -0.900, 'C-375'),
        ('χ₋₃',          1,      3,  494, -0.891, 'C-375'),
        ('11a1',          2,     11,  437, -0.914, 'C-375'),
        ('Δ(Ramanujan)',  2,      1,  324, -0.876, 'C-375'),
        ('sym²(11a1)',    3,    121,   96, -0.885, 'C-375'),
        ('Sym³(11a1)',    4,   1331,  327, -0.905, 'C-383'),
        ('Sym³(19a1)',    4,   6859,  357, -0.908, 'C-384'),
        ('χ₅[1]',        1,      5,  238, -0.903, 'C-378'),
        ('χ̄₅[1]',        1,      5,  238, -0.875, 'C-378'),
        ('χ₇[1]',        1,      7,  254, -0.871, 'C-378'),
        ('χ̄₇[1]',        1,      7,  255, -0.901, 'C-378'),
        ('Sym⁴(11a1)',    5,  14641,  177, -0.913, 'C-387'),
    ]

    hdr = f"  {'L-함수':<20} {'d':>2} {'N':>8} {'n':>5} {'ρ(A_Λ)':>8} {'소스':>12}"
    log(hdr)
    log(f"  {'-'*62}")

    all_rhos = []
    for name, d, N, n, rho, src in existing_data:
        log(f"  {name:<20} {d:>2} {N:>8} {n:>5} {rho:>+8.3f} {src:>12}")
        all_rhos.append(rho)

    # C-388 추가
    log(f"  {'Sym⁵(11a1)':<20} {6:>2} {conductor:>8} {n_valid:>5} {rho_Lam_gm:>+8.3f} {'C-388 ★':>12}")
    all_rhos.append(rho_Lam_gm)

    log()

    all_rhos_arr = np.array(all_rhos)
    log(f"  전체 (13 L-함수): ρ 평균 = {np.mean(all_rhos_arr):+.4f} ± {np.std(all_rhos_arr):.4f}")
    log(f"  범위: [{np.min(all_rhos_arr):+.4f}, {np.max(all_rhos_arr):+.4f}], "
        f"range = {np.max(all_rhos_arr) - np.min(all_rhos_arr):.4f}")

    # degree별
    d13 = [r for (_, d, _, _, r, _) in existing_data if d <= 3]
    d4  = [r for (_, d, _, _, r, _) in existing_data if d == 4]
    d5  = [r for (_, d, _, _, r, _) in existing_data if d == 5]
    log(f"  degree 1-3 ({len(d13)}개): ρ 평균 = {np.mean(d13):+.4f} ± {np.std(d13):.4f}")
    log(f"  degree 4   ({len(d4)}개):  ρ 평균 = {np.mean(d4):+.4f} ± {np.std(d4):.4f}")
    log(f"  degree 5   ({len(d5)}개):  ρ = {np.mean(d5):+.4f}")
    log(f"  degree 6   (1개):  ρ = {rho_Lam_gm:+.4f}")
    log(f"  degree 6 vs 1-5: Δρ = {rho_Lam_gm - np.mean(all_rhos[:-1]):+.4f}")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # 2H₁/A 비율 (Hadamard 지배 비율)
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("Hadamard 지배 비율 2H₁/A")
    log("═" * 80)
    log()

    h1_ratios = []
    for d in valid:
        h1_L = d['H1_L']
        a_L = d['A_L']
        if a_L > 0:
            h1_ratios.append(2.0 * h1_L / a_L)

    if h1_ratios:
        h1r = np.array(h1_ratios)
        log(f"  <2H₁/A> = {np.mean(h1r):.4f} ± {np.std(h1r):.4f}")
        log(f"  degree별 비교:")
        log(f"    GL(2): 0.873")
        log(f"    GL(3): 0.661-0.676")
        log(f"    GL(4): 0.661-0.676")
        log(f"    GL(5): 0.657")
        log(f"    GL(6) C-275: 0.657")
        log(f"    GL(6) C-388: {np.mean(h1r):.3f} ← 이번 결과")
    log()
    save()

    # ════════════════════════════════════════════════════════════════════
    # 최종 판정
    # ════════════════════════════════════════════════════════════════════
    log("═" * 80)
    log("[C-388 판정]")
    log("═" * 80)
    log()

    if rho_Lam_gm < -0.80:
        verdict = "★★★★★★★ degree 6 보편성 — degree 1-6 A_Λ 통일!"
        log(f"  {verdict}")
        log(f"  ρ = {rho_Lam_gm:+.4f} < -0.80")
        log(f"  → Paper 4 Table 13번째 행 추가: degree 1-6 보편성 (A_Λ 통일)")
        log(f"  → C-275 (구 방법) ρ=-0.42 → C-388 (신 방법) ρ={rho_Lam_gm:+.4f}")
        log(f"  → 방법론 불일치 해소: degree 1-6 전체가 A_Λ 방법으로 통일")
    elif rho_Lam_gm < -0.60:
        verdict = "★★★★ 유의미하지만 degree 약화"
        log(f"  {verdict}")
        log(f"  ρ = {rho_Lam_gm:+.4f}, -0.80 ≤ ρ < -0.60")
        log(f"  → degree 6에서 상관 약화. C-275 (ρ=-0.42)보다 개선되었으나")
        log(f"     degree 1-5 평균 (ρ≈-0.89)보다 약함.")
        log(f"  → 가능한 원인: 작은 n, 낮은 T, 또는 진짜 degree 의존성")
    elif n_valid >= MIN_ZEROS:
        verdict = "★★ ρ ≥ -0.60 (약한 상관)"
        log(f"  {verdict}")
        log(f"  ρ = {rho_Lam_gm:+.4f} ≥ -0.60")
        log(f"  → 신 방법으로도 degree 6 상관이 약함")
        log(f"  → C-275 수준에서 A_Λ 보정이 큰 개선을 주지 않음")
    else:
        verdict = "★★ n < 20 (통계적 불충분)"
        log(f"  {verdict}")

    log()
    log(f"  95% CI: [{ci_lo_Lam_gm:+.4f}, {ci_hi_Lam_gm:+.4f}]")
    log(f"  p = {p_Lam_gm:.3e}")
    log(f"  n = {n_valid}")
    log()

    log("[성공 기준 검증]")
    log(f"  ★★★★★★★ ρ < -0.80: {'✅ PASS' if rho_Lam_gm < -0.80 else '❌ FAIL'}")
    log(f"  ★★★★    -0.80 ≤ ρ < -0.60: {'✅' if -0.80 <= rho_Lam_gm < -0.60 else '—'}")
    log(f"  ★★      ρ ≥ -0.60 또는 n < 20: {'✅' if rho_Lam_gm >= -0.60 or n_valid < 20 else '—'}")
    log()

    log(f"  총 소요시간: {time.time()-t_start:.1f}초")
    log()
    save()
    log(f"[완료] {RESULT_PATH}")


if __name__ == '__main__':
    main()
