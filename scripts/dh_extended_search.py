#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] DH 함수 Off-critical 영점 확장 탐색
=============================================================================
결과 #30 보충: t∈[60,200], σ∈[0.5,1.2] 범위에서 off-critical 영점 탐색.
DH 첫 off-critical zero가 t>60에 있을 가능성 (Spira 1968 참조).
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

mpmath.mp.dps = 80
START = time.time()

print("=" * 70)
print("[DH 확장 탐색] t∈[1,200], σ∈[0.5,1.2]")
print("=" * 70)
print(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print()

# DH 함수 정의 (본 스크립트에서 재정의)
CHI_MOD5 = [0, 1, 1j, -1j, -1]
CHI_MOD5_BAR = [0, 1, -1j, 1j, -1]

sqrt5 = mpmath.sqrt(5)
kappa_DH = (mpmath.sqrt(10 - 2*sqrt5) - 2) / (sqrt5 - 1)
coeff_chi = (1 - 1j * kappa_DH) / 2
coeff_chi_bar = (1 + 1j * kappa_DH) / 2


def dh_func(s):
    s = mpmath.mpc(s)
    L_chi = mpmath.dirichlet(s, CHI_MOD5)
    L_chi_bar = mpmath.dirichlet(s, CHI_MOD5_BAR)
    return coeff_chi * L_chi + coeff_chi_bar * L_chi_bar


# 전략: σ를 고정하고 t를 스캔하며 |f|의 지역 최솟값 추적
# 더 넓은 σ 범위와 더 높은 t 범위

sigmas_test = [0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10]
t_ranges = [(1, 50, 500), (50, 100, 500), (100, 150, 500), (150, 200, 500)]

all_candidates = []

for sigma in sigmas_test:
    print(f"\n--- σ = {sigma:.2f} ---", flush=True)
    for t_lo, t_hi, n_t in t_ranges:
        ts = np.linspace(t_lo, t_hi, n_t)
        min_absf = float('inf')
        min_t = 0
        local_mins = []

        prev_absf = None
        prev_prev_absf = None
        for j, t in enumerate(ts):
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            absf = float(abs(dh_func(s)))

            # 지역 최솟값 탐지
            if prev_absf is not None and prev_prev_absf is not None:
                if prev_absf < prev_prev_absf and prev_absf < absf and prev_absf < 0.5:
                    local_mins.append((sigma, ts[j-1], prev_absf))

            prev_prev_absf = prev_absf
            prev_absf = absf

            if absf < min_absf:
                min_absf = absf
                min_t = t

        if local_mins:
            for lm in local_mins:
                print(f"  t∈[{t_lo},{t_hi}]: 지역 최솟값 t={lm[1]:.4f}, |f|={lm[2]:.6f}", flush=True)
                all_candidates.append(lm)
        else:
            print(f"  t∈[{t_lo},{t_hi}]: 최소 |f|={min_absf:.6f} at t={min_t:.2f}", flush=True)

print(f"\n\n총 후보: {len(all_candidates)}개", flush=True)
all_candidates.sort(key=lambda x: x[2])

# 상위 후보에서 findroot 시도
off_critical_zeros = []

print("\n--- findroot 정밀화 ---", flush=True)
for sigma_c, t_c, absf_c in all_candidates[:50]:
    if abs(sigma_c - 0.5) < 0.02:
        continue
    if any(abs(sigma_c - z[0]) < 0.01 and abs(t_c - z[1]) < 0.1 for z in off_critical_zeros):
        continue

    try:
        def f_system(sv, tv):
            sc = mpmath.mpf(sv) + 1j * mpmath.mpf(tv)
            fv = dh_func(sc)
            return (mpmath.re(fv), mpmath.im(fv))

        result = mpmath.findroot(f_system, (mpmath.mpf(str(sigma_c)), mpmath.mpf(str(t_c))))
        sr = float(result[0])
        tr = float(result[1])

        sv = mpmath.mpf(str(sr)) + 1j * mpmath.mpf(str(tr))
        res = float(abs(dh_func(sv)))

        if res < 1e-10:
            if not any(abs(sr - z[0]) < 0.01 and abs(tr - z[1]) < 0.1 for z in off_critical_zeros):
                off_critical_zeros.append((sr, tr, res))
                is_off = abs(sr - 0.5) > 0.02
                marker = "★★★ OFF-CRITICAL" if is_off else "on-critical"
                print(f"  {marker}: σ={sr:.6f}, t={tr:.6f}, |f|={res:.2e}", flush=True)
    except Exception as e:
        pass

print(f"\n\n총 영점: {len(off_critical_zeros)}개", flush=True)
genuinely_off = [(s, t, r) for s, t, r in off_critical_zeros if abs(s - 0.5) > 0.02]
print(f"진짜 off-critical (|σ-0.5|>0.02): {len(genuinely_off)}개", flush=True)

for s, t, r in genuinely_off:
    print(f"  σ={s:.6f}, t={t:.6f}, |f|={r:.2e}", flush=True)

elapsed = time.time() - START
print(f"\n소요: {elapsed:.1f}초 ({elapsed/60:.1f}분)")
print(f"종료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
