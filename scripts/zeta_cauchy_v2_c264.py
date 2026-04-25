#!/usr/bin/env python3
"""
[C-264 교차검증 v5] ζ(s) Cauchy — C-263 방식 그대로 적용
  lfunlambda 직접 사용 (lfuninit 없이), Cauchy 루프 PARI 내부 실행
"""
import sys, math, time
import numpy as np
from scipy import stats

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(256 * 10**6)
pari.set_real_precision(38)

CAUCHY_R = 1e-4
N_PTS = 64
H_DEC = '1e-5'
ABS_THR = '1e-30'

def pari_to_float(x):
    try:
        return float(str(x).replace(' E','e'))
    except:
        return float('nan')

# L-함수 객체
print("[ζ(s)] L-함수 초기화...")
pari('Lz = lfuncreate(1)')
print("  Lz OK")

# 영점 수집
print("[ζ(s)] 영점 수집...")
pari('zvz = lfunzeros(Lz, 550)')
n_z = int(str(pari('#zvz')))
zeros = []
for i in range(1, n_z + 1):
    t = pari_to_float(pari(f'zvz[{i}]'))
    if not math.isnan(t) and t > 5:
        zeros.append(t)
zeros = sorted(zeros)
print(f"  {len(zeros)}개, t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")

# Cauchy integral — C-263 방식 (PARI 내부 루프)
print(f"[Cauchy] R={CAUCHY_R}, N={N_PTS}...")
data = []
t_start = time.time()

for idx, t0 in enumerate(zeros):
    script = (
        f'rho=1/2+{t0:.15f}*I; c0s=0+0.*I; c1s=0+0.*I; cnt=0; '
        f'for(k=0,{N_PTS-1}, '
        f'th=2*Pi*k/{N_PTS}; uu={CAUCHY_R}*exp(th*I); ss=rho+uu; '
        f'Lv=lfunlambda(Lz,ss); '
        f'if(abs(Lv)<{ABS_THR},next); '
        f'Lp=lfunlambda(Lz,ss+{H_DEC}); '
        f'Lm=lfunlambda(Lz,ss-{H_DEC}); '
        f'LpL=(Lp-Lm)/(2*{H_DEC}*Lv); '
        f'gg=LpL-1/uu; c0s=c0s+gg; c1s=c1s+gg/uu; cnt=cnt+1); '
        f'[c0s/cnt, c1s/cnt, cnt]'
    )
    try:
        res = pari(script)
        c0_re = pari_to_float(pari(f'real(({res})[1])'))
        c0_im = pari_to_float(pari(f'imag(({res})[1])'))
        c1_re = pari_to_float(pari(f'real(({res})[2])'))
        n_cnt = pari_to_float(pari(f'({res})[3]'))
        if any(math.isnan(v) for v in [c0_re, c0_im, c1_re]):
            continue
        if n_cnt < N_PTS // 2:
            continue
        A = c0_im**2 + 2*c1_re
        if math.isnan(A) or math.isinf(A) or A <= 0:
            continue
        data.append({'t': t0, 'A': A, 'c0_im': c0_im, 'c1_re': c1_re})
    except:
        continue

    if (idx+1) % 50 == 0:
        elapsed = time.time() - t_start
        print(f"  {idx+1}/{len(zeros)} ({len(data)} ok, {elapsed:.0f}s)")

elapsed = time.time() - t_start
print(f"  {len(data)}/{len(zeros)} 유효 ({elapsed:.0f}s)")

if len(data) < 10:
    print("❌ 유효 데이터 부족")
    sys.exit(1)

# 간격 + GUE 정규화
valid = []
for i in range(2, len(data)-2):
    d = data[i]
    gap_r = data[i+1]['t'] - d['t']
    gap_l = d['t'] - data[i-1]['t']
    if gap_r <= 0 or gap_l <= 0:
        continue
    # GUE: d_bar = (1/π) log(T/(2π)) for ζ(s)
    d_bar = 2.0 / (data[i+1]['t'] - data[i-1]['t'])
    d['gap_r_gue'] = gap_r * d_bar
    d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
    valid.append(d)

print(f"  내부: {len(valid)}")

A_arr = np.array([d['A'] for d in valid])
gr = np.array([d['gap_r_gue'] for d in valid])
gm = np.array([d['gap_min_gue'] for d in valid])

rho_r, p_r = stats.spearmanr(A_arr, gr)
rho_m, p_m = stats.spearmanr(A_arr, gm)
sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')

print(f"\n[Spearman Cauchy/lfunlambda] (n={len(valid)})")
print(f"  ρ(A, gap_right_GUE) = {rho_r:+.4f}  (p={p_r:.3e})  {sig(p_r)}  [C-256: -0.5898]")
print(f"  ρ(A, gap_min_GUE)   = {rho_m:+.4f}  (p={p_m:.3e})  {sig(p_m)}")

if abs(rho_r - (-0.59)) < 0.15:
    print(f"\n✅ 교차검증 통과: ρ(gap_right) = {rho_r:.4f} ≈ -0.59")
else:
    print(f"\n❌ gap_right = {rho_r:.4f} (C-256: -0.5898)")
    if abs(rho_m - (-0.59)) < 0.15:
        print(f"  ⚠️ gap_min은 일치: {rho_m:.4f}")
