#!/usr/bin/env python3
"""
[C-264 교차검증 v2] ζ(s) δ-추출 A(γ) — C-256 Cauchy 결과 재현 확인
  δ-extraction: Λ'/Λ(center+δ+it₀)로 A 추출, ρ(A, gap_right) ≈ -0.59 확인
"""
import sys, math
import numpy as np
from scipy import stats

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(128 * 10**6)
pari.set_real_precision(38)

CENTER = 0.5

# 영점 수집
print("[ζ(s)] 영점 수집...")
pari('Li_z = lfuninit(lfuncreate(1), [0, 600])')
pari('zv = lfunzeros(Li_z, 550)')
n_z = int(str(pari('#zv')))
zeros = []
for i in range(1, n_z + 1):
    t = float(str(pari(f'zv[{i}]')).replace(' E','e'))
    if t > 5:
        zeros.append(t)
zeros = sorted(zeros)
print(f"  {len(zeros)}개")

# δ-extraction
delta = 0.005
print(f"[δ-extraction] δ={delta}...")
data = []
for idx, t0 in enumerate(zeros):
    s_str = f'{CENTER + delta} + {t0}*I'
    try:
        Lv = pari(f'lfunlambda(Li_z, {s_str})')
        Lv_c = complex(Lv)
        if abs(Lv_c) < 1e-30:
            continue
        Ld = pari(f'lfunlambda(Li_z, {s_str}, 1)')
        Ld_c = complex(Ld)
        LpL = Ld_c / Lv_c
        S1 = -LpL.imag
        H1 = (LpL.real - 1.0/delta) / delta
        A = S1**2 + 2*H1
        if A > 0:
            data.append({'t': t0, 'A': A, 'S1': S1, 'H1': H1})
    except:
        continue
    if (idx+1) % 50 == 0:
        print(f"  {idx+1}/{len(zeros)}")

print(f"  {len(data)} 유효")

# 간격 계산
valid = []
for idx in range(2, len(data)-2):
    d = data[idx]
    gap_r = data[idx+1]['t'] - d['t']
    gap_l = d['t'] - data[idx-1]['t']
    if gap_r <= 0 or gap_l <= 0:
        continue
    d_bar = 2.0 / (data[idx+1]['t'] - data[idx-1]['t'])
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

print(f"\n[Spearman δ-extraction] (n={len(valid)})")
print(f"  ρ(A, gap_right_GUE) = {rho_r:+.4f}  (p={p_r:.3e})  {sig(p_r)}  [C-256: -0.5898]")
print(f"  ρ(A, gap_min_GUE)   = {rho_m:+.4f}  (p={p_m:.3e})  {sig(p_m)}")

if abs(rho_r - (-0.59)) < 0.1:
    print(f"\n✅ δ-extraction 교차검증 통과: ρ = {rho_r:.4f} ≈ -0.59")
else:
    print(f"\n❌ δ-extraction ρ = {rho_r:.4f}, C-256 = -0.5898, 차이 = {abs(rho_r-(-0.5898)):.4f}")
