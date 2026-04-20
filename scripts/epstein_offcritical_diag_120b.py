#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #120b] Epstein ζ off-critical 진단 — 1D 스캔 + 인수 원리 비교
=============================================================================
목표: Epstein ζ (h>1)에서 off-critical 영점의 존재 여부를 2가지 방법으로 판단:
  (A) |Λ(σ+it)| 1D 스캔: σ=0.3,0.4,0.6,0.7에서 t∈[1,1000] 스캔, 국소 최솟값 추적
  (B) 인수 원리: N(T) vs N_0(T) 비교 (총 영점 수 vs 임계선 영점 수)

대상: x²+5y² (disc=-20, h=2) — 가장 단순한 h>1 형식

결과: results/epstein_offcritical_diag_120b.txt
=============================================================================
"""

import sys, os, time
import numpy as np
from datetime import datetime

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 100)

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'epstein_offcritical_diag_120b.txt')
outf = open(RESULT_FILE, 'w')

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

log("=" * 72)
log("[실험 #120b] Epstein ζ off-critical 진단")
log("=" * 72)
log(f"시작: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log()

M_STR = "[1,0;0,5]"
NAME = "x²+5y² (disc=-20, h=2)"
lfun_obj = pari(f"lfunqf({M_STR})")
log(f"대상: {NAME}")
log()

def Lambda_eval(sigma, t):
    s_str = f"{sigma} + {t}*I" if t >= 0 else f"{sigma} - {abs(t)}*I"
    return complex(pari.lfun(lfun_obj, pari(s_str), 1))


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Part A: |Λ(σ+it)| 1D 스캔
# ━���━━━━━━━━��━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━━━━━━━━━━
log("━━━ Part A: |Λ| 1D 스캔 ━━━")

SIGMAS = [0.3, 0.4, 0.6, 0.7]
T_MAX = 1000
DT = 0.2  # 간격 0.2 → 5000점/라인

for sigma in SIGMAS:
    log(f"\n  σ={sigma}, t∈[1,{T_MAX}], ��t={DT}")
    t0 = time.time()

    ts = np.arange(1.0, T_MAX + DT, DT)
    abs_vals = []
    local_mins = []

    for i, t in enumerate(ts):
        try:
            val = abs(Lambda_eval(sigma, t))
        except Exception:
            val = 1e30
        abs_vals.append(val)

        # 진행률 (200점마다)
        if (i + 1) % 1000 == 0:
            log(f"    진행: {i+1}/{len(ts)} ({time.time()-t0:.1f}초)")

    abs_arr = np.array(abs_vals)

    # 국소 최솟값 찾기 (이웃보다 작은 점)
    for i in range(1, len(abs_arr) - 1):
        if abs_arr[i] < abs_arr[i-1] and abs_arr[i] < abs_arr[i+1]:
            if abs_arr[i] < 0.1:  # 상대적으로 작은 값
                local_mins.append((ts[i], abs_arr[i]))

    # 가장 작은 10개
    local_mins.sort(key=lambda x: x[1])
    log(f"    완료: {time.time()-t0:.1f}초")
    log(f"    전역 최소 |Λ| = {abs_arr.min():.6e} at t={ts[abs_arr.argmin()]:.2f}")
    log(f"    |Λ|<0.1 국소 최솟값: {len(local_mins)}개")
    if local_mins:
        log(f"    최소 10개:")
        for t_val, v_val in local_mins[:10]:
            log(f"      t={t_val:.2f}, |Λ|={v_val:.6e}")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Part B: 인수 원리 비교 — N(T) vs N_0(T)
# ━━━━━━���━━━━━���━━━━━━���━━━━━━━━━━━━━━━━���━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log(f"\n━━━ Part B: N(T) vs N_0(T) ━━━")

# N_0(T) = on-critical 영점 수 (PARI lfunzeros)
T_checks = [100, 200, 500, 1000]

log(f"\n  {'T':>6} {'N_0(T)':>8} {'N_asymp(T)':>12} {'차이':>8}")
log(f"  {'─'*40}")

for T in T_checks:
    try:
        # on-critical 영점 계수
        zeros = pari(f"lfunzeros(lfunqf({M_STR}), {T})")
        n0 = len([float(z) for z in zeros if float(z) > 0.5])
    except Exception as e:
        log(f"  T={T}: lfunzeros 실패 ({e})")
        continue

    # 점근 공식: N(T) ≈ (T/π)·log(√|D|·T/(2πe))
    # Epstein for binary form Q with |disc|=D:
    # Λ(s) = (D/(4π²))^{s/2} Γ(s) Z_Q(s)
    # conductor = D = 20
    # N(T) ≈ T/π · log(T·√D/(2πe)) + O(1)
    D = 20
    N_asymp = T / np.pi * np.log(T * np.sqrt(D) / (2 * np.pi * np.e))

    log(f"  {T:>6} {n0:>8} {N_asymp:>12.1f} {n0 - N_asymp:>8.1f}")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Part C: 정밀 인수 원리 (선택적 — 정확한 N(T))
# ━━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log(f"\n━━━ Part C: 정밀 인수 원리 — 윤곽 적분 ━━━")

def count_zeros_rectangle(T, sigma_left=0.01, sigma_right=0.99, n_points=2000):
    """
    사각형 {sigma_left ≤ σ ≤ sigma_right, 0 ≤ t ≤ T} 내부 영점 수.
    인수 원리: N = (1/2π) × arg(Λ) 변화량 (반시계 윤곽).

    윤곽: 아래(t=0.5) → 오른쪽(σ=sigma_right) → 위(t=T) → 왼쪽(σ=sigma_left)
    """
    # 아래쪽: σ = sigma_left→sigma_right, t = 0.5
    total_phase = 0.0
    prev_phase = None

    contour_points = []
    # Bottom: left→right
    sigmas = np.linspace(sigma_left, sigma_right, n_points // 4)
    for s in sigmas:
        contour_points.append((s, 0.5))
    # Right: bottom→top
    ts = np.linspace(0.5, T, n_points // 4)
    for t in ts:
        contour_points.append((sigma_right, t))
    # Top: right→left
    sigmas_rev = np.linspace(sigma_right, sigma_left, n_points // 4)
    for s in sigmas_rev:
        contour_points.append((s, T))
    # Left: top→bottom
    ts_rev = np.linspace(T, 0.5, n_points // 4)
    for t in ts_rev:
        contour_points.append((sigma_left, t))

    phases = []
    for (s, t) in contour_points:
        try:
            val = Lambda_eval(s, t)
            phases.append(np.angle(val))
        except Exception:
            if phases:
                phases.append(phases[-1])
            else:
                phases.append(0)

    # 위상 변화 누적
    total = 0.0
    for k in range(len(phases) - 1):
        diff = phases[k+1] - phases[k]
        while diff > np.pi:
            diff -= 2 * np.pi
        while diff < -np.pi:
            diff += 2 * np.pi
        total += diff

    return total / (2 * np.pi)


for T in [100, 200]:
    t0 = time.time()
    N_contour = count_zeros_rectangle(T, n_points=4000)
    elapsed = time.time() - t0

    # on-critical count
    zeros = pari(f"lfunzeros(lfunqf({M_STR}), {T})")
    n0 = len([float(z) for z in zeros if float(z) > 0.5])

    log(f"  T={T}: N_contour={N_contour:.2f} (반올림: {round(N_contour)}), N_0={n0}, "
        f"차이={round(N_contour)-n0}, {elapsed:.1f}초")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 종합
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━��━━━━━━━━━━━━━━━━━━
log(f"\n{'='*72}")
log("종합 진단")
log(f"{'='*72}")
log(f"  t∈[1,1000] 4개 σ 라인 1D 스캔 + 인수 원리 비교 완료.")
log(f"  N_contour ≈ N_0이면: 이 T 범위에 off-critical 영점 없음")
log(f"  N_contour > N_0이면: off-critical 영점 존재 (위치는 2D 탐색 필요)")
log(f"  |Λ| 최솟값 패턴으로 off-critical 영점 존재 가능성 판단")

elapsed_total = time.time() - float(os.environ.get('START_TIME', time.time()))
log(f"\n  총 소요: 로그 참조")

outf.close()
