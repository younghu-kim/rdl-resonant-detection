#!/usr/bin/env python3
"""
[C-264 교차검증] ζ(s) 영점-합 A(γ) vs C-256 Cauchy A(γ) 비교
  영점-합 방법이 ρ(A, gap_right) ≈ -0.59를 재현하는지 확인.
  GL(1) center=0.5, μ=[0]
"""
import sys, os, math
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 30

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(128 * 10**6)
pari.set_real_precision(38)

CENTER = 0.5  # ζ(s) critical line

# 영점 수집 — 넓은 범위
print("[ζ(s)] 영점 수집 (T=2000)...")
pari('Li_z = lfuninit(lfuncreate(1), [0, 2100])')
pari('zv = lfunzeros(Li_z, 2000)')
n_z = int(str(pari('#zv')))
zeros = []
for i in range(1, n_z + 1):
    t = float(str(pari(f'zv[{i}]')).replace(' E','e'))
    if t > 5:
        zeros.append(t)
zeros = sorted(zeros)
print(f"  {len(zeros)}개, t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")

# Gamma smooth for ζ: Im part of (1/2)[ψ(s/2) - log(π)]
def gamma_smooth_im(gamma_0):
    s = mpmath.mpc(CENTER, gamma_0)
    val = (mpmath.digamma(s / 2) - mpmath.log(mpmath.pi)) / 2
    return float(mpmath.im(val))

# A(γ) from zero-sum — conjugate zeros 포함
# Λ'/Λ(1/2+iγ₀) 의 Im = Σ_{γₖ>0} [1/(γ₀-γₖ) + 1/(γ₀+γₖ)] + Gamma smooth
# 1/(γ₀+γₖ) 항 = conjugate zero ρ*=1/2-iγₖ 기여
print("[A(γ)] 영점 합 (conjugate 포함)...")
data = []
for i in range(len(zeros)):
    g0 = zeros[i]
    # 같은 부호 영점: 1/(γ₀-γₖ)
    S1_same = sum(1.0/(g0-zeros[k]) for k in range(len(zeros)) if k!=i)
    H1_same = sum(1.0/(g0-zeros[k])**2 for k in range(len(zeros)) if k!=i)
    # conjugate 영점: 1/(γ₀+γₖ) for all k (including k=i)
    S1_conj = sum(1.0/(g0+zeros[k]) for k in range(len(zeros)))
    H1_conj = sum(1.0/(g0+zeros[k])**2 for k in range(len(zeros)))
    S1_full = S1_same + S1_conj
    H1_full = H1_same + H1_conj
    sm_im = gamma_smooth_im(g0)
    S1_corr = S1_full - sm_im
    A = S1_corr**2 + 2*H1_full
    # 비교용: same-only (이전 방법)
    S1_same_corr = S1_same - sm_im
    A_same = S1_same_corr**2 + 2*H1_same
    if A > 0:
        data.append({'t': g0, 'A': A, 'A_same': A_same, 'S1': S1_corr, 'H1': H1_full})

print(f"  {len(data)} 유효")

# 내부 — 중앙 60%만 사용 (가장자리 truncation 회피)
N = len(data)
lo = int(N * 0.2)
hi = int(N * 0.8)
valid = []
for idx in range(lo, hi):
    if idx <= 0 or idx >= N-1:
        continue
    d = data[idx]
    t_n = d['t']
    gap_r = data[idx+1]['t'] - t_n
    gap_l = t_n - data[idx-1]['t']
    if gap_r <= 0 or gap_l <= 0:
        continue
    d_bar = 2.0 / (data[idx+1]['t'] - data[idx-1]['t'])
    d['gap_r_gue'] = gap_r * d_bar
    d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
    d['S1_plus_2overDR'] = d['S1'] + 2.0/gap_r
    d['H1_frac'] = 2*d['H1']/d['A'] if d['A'] > 0 else float('nan')
    valid.append(d)

print(f"  내부 (중앙 60%): {len(valid)}, t ∈ [{valid[0]['t']:.1f}, {valid[-1]['t']:.1f}]")

A_arr = np.array([d['A'] for d in valid])
A_same = np.array([d['A_same'] for d in valid])
gr = np.array([d['gap_r_gue'] for d in valid])
gm = np.array([d['gap_min_gue'] for d in valid])

rho_r, p_r = stats.spearmanr(A_arr, gr)
rho_same, p_same = stats.spearmanr(A_same, gr)
rho_m, p_m = stats.spearmanr(A_arr, gm)

sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')

print(f"\n[Spearman] (n={len(valid)})")
print(f"  ρ(A_full, gap_right_GUE)  = {rho_r:+.4f}  (p={p_r:.3e})  {sig(p_r)}  [C-256 참조: -0.5898]")
print(f"  ρ(A_same, gap_right_GUE)  = {rho_same:+.4f}  (p={p_same:.3e})  {sig(p_same)}  [conjugate 미포함]")
print(f"  ρ(A_full, gap_min_GUE)    = {rho_m:+.4f}  (p={p_m:.3e})")

# Prop 12
s12 = np.array([d['S1_plus_2overDR'] for d in valid])
hf = np.array([d['H1_frac'] for d in valid])
n_pos = int(np.sum(s12 > 0))
print(f"\n[Prop 12]")
print(f"  S₁+2/Δ_R > 0: {n_pos}/{len(valid)} ({100*n_pos/len(valid):.0f}%)")
print(f"  2H₁/A: {np.nanmean(hf):.3f} ± {np.nanstd(hf):.3f}")

# 판정
if abs(rho_r - (-0.59)) < 0.1:
    print(f"\n✅ 교차검증 통과: 영점-합 ρ = {rho_r:.4f} ≈ Cauchy ρ = -0.5898")
else:
    print(f"\n❌ 교차검증 실패: 영점-합 ρ = {rho_r:.4f} ≠ Cauchy ρ = -0.5898")
    print(f"  차이 = {abs(rho_r - (-0.5898)):.4f} — 영점-합 방법 불일치")
