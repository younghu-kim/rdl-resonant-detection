#!/usr/bin/env python3
"""
ζ(s) 첫 N개 비자명 영점 t_n 캐시 파일 생성.
결과: outputs/cache/zeta_zeros_t_N5000.npy
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))
import mpmath

DPS = 30  # 영점 위치는 30자리면 충분
mpmath.mp.dps = DPS

N_MAX = 2000
CACHE_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/outputs/cache/zeta_zeros_t_N5000.npy')  # same path, reused
os.makedirs(os.path.dirname(CACHE_PATH), exist_ok=True)

if os.path.exists(CACHE_PATH):
    arr = np.load(CACHE_PATH)
    if len(arr) >= N_MAX:
        print(f"캐시 이미 존재: {len(arr)}개. 종료.", flush=True)
        sys.exit(0)
    else:
        print(f"캐시 부족 ({len(arr)}개). 재생성.", flush=True)

print(f"ζ 첫 {N_MAX}개 영점 계산 중 (DPS={DPS})...", flush=True)
t_start = time.time()
t_vals = []
for n in range(1, N_MAX+1):
    rho = mpmath.zetazero(n)
    t_vals.append(float(rho.imag))
    if n % 500 == 0:
        print(f"  {n}/{N_MAX} ({time.time()-t_start:.0f}s)", flush=True)

arr = np.array(t_vals)
np.save(CACHE_PATH, arr)
print(f"저장 완료: {CACHE_PATH}  ({len(arr)}개, {time.time()-t_start:.0f}s)", flush=True)
