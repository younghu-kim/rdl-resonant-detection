#!/usr/bin/env python3
"""
[C-264 교차검증 v3] ζ(s) Cauchy 방법 — C-256 재현
  C-263 GL(2) 방법을 ζ(s)에 적용. center=0.5
  lfun 직접 사용 (lfuninit 불필요)
"""
import sys, math, cmath
import numpy as np
from scipy import stats

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(256 * 10**6)
pari.set_real_precision(38)

CENTER = 0.5
N_CAUCHY = 64
RADIUS = 1e-4

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

# Cauchy 적분으로 c₀, c₁ 추출
# L'/L(s) = 1/(s-ρ₀) + c₀ + c₁(s-ρ₀) + ...
# c₀ = (1/2πi) ∮ L'/L(s)/(s-ρ₀)⁰ ds (on circle |s-ρ₀|=r)
# c₁ = (1/2πi) ∮ L'/L(s)/(s-ρ₀)¹ ds (minus 1/(s-ρ₀)² from simple pole)
print(f"[Cauchy] N={N_CAUCHY}, r={RADIUS}...")

def eval_lfun(sigma, t):
    """lfunlambda at s = sigma + it, returns complex"""
    s_str = f'{sigma} + {t}*I'
    v = pari(f'lfunlambda(Li_z, {s_str})')
    return complex(v)

data = []
for idx, t0 in enumerate(zeros):
    rho0 = complex(CENTER, t0)
    # Cauchy contour: s = ρ₀ + r·e^{2πiθ}
    vals_LpL = []
    ok = True
    for j in range(N_CAUCHY):
        theta = 2 * math.pi * j / N_CAUCHY
        ds = RADIUS * cmath.exp(1j * theta)
        s = rho0 + ds
        try:
            L_val = eval_lfun(s.real, s.imag)
        except:
            ok = False; break
        if abs(L_val) < 1e-50:
            ok = False; break
        # finite diff for L'/L
        h = 1e-5
        ds2 = complex(h, 0)
        try:
            L_plus = eval_lfun((s+ds2).real, (s+ds2).imag)
        except:
            ok = False; break
        LpL = (L_plus - L_val) / (h * L_val)
        vals_LpL.append((LpL, ds))

    if not ok or len(vals_LpL) != N_CAUCHY:
        continue

    # c₀ = (1/N) Σ [L'/L(s_j) - 1/ds_j]  (remove simple pole)
    c0 = sum((LpL - 1.0/ds) for LpL, ds in vals_LpL) / N_CAUCHY
    # c₁ = (1/N) Σ [L'/L(s_j) - 1/ds_j] / ds_j
    c1 = sum((LpL - 1.0/ds) / ds for LpL, ds in vals_LpL) / N_CAUCHY

    A = c0.imag**2 + 2*c1.real
    if A > 0:
        data.append({'t': t0, 'A': A, 'c0': c0, 'c1': c1})

    if (idx+1) % 50 == 0:
        print(f"  {idx+1}/{len(zeros)} ({len(data)} ok)")

print(f"  {len(data)} 유효")

# 간격 + 상관
valid = []
for i in range(2, len(data)-2):
    d = data[i]
    gap_r = data[i+1]['t'] - d['t']
    gap_l = d['t'] - data[i-1]['t']
    if gap_r <= 0 or gap_l <= 0:
        continue
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

print(f"\n[Spearman Cauchy] (n={len(valid)})")
print(f"  ρ(A, gap_right_GUE) = {rho_r:+.4f}  (p={p_r:.3e})  {sig(p_r)}  [C-256: -0.5898]")
print(f"  ρ(A, gap_min_GUE)   = {rho_m:+.4f}  (p={p_m:.3e})  {sig(p_m)}")

if abs(rho_r - (-0.59)) < 0.15:
    print(f"\n✅ Cauchy 교차검증 통과: ρ(gap_right) = {rho_r:.4f} ≈ -0.59")
elif abs(rho_m - (-0.59)) < 0.15:
    print(f"\n⚠️ gap_min은 일치: ρ(gap_min) = {rho_m:.4f} ≈ -0.59, gap_right = {rho_r:.4f}")
else:
    print(f"\n❌ 불일치: gap_right = {rho_r:.4f}, gap_min = {rho_m:.4f}")
