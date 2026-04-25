#!/usr/bin/env python3
"""
[C-264 교차검증 v4] ζ(s) zero-sum + 완전 Gamma 보정
  A_Λ = (S₁_L + Im(γ'))² + 2(H₁_L + Re(ψ'(ρ₀/2)/4))
  여기서 γ' = γ_R'/γ_R(ρ₀), conjugate zeros 포함
"""
import sys, math
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 30

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
import cypari2
pari = cypari2.Pari()
pari.allocatemem(128 * 10**6)
pari.set_real_precision(38)

CENTER = 0.5

# 영점 수집
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

def gamma_corrections(gamma_0):
    """
    ζ(s)의 Gamma factor: γ_R(s) = π^{-s/2} Γ(s/2)
    γ_R'/γ_R(s) = -log(π)/2 + ψ(s/2)/2
    (γ_R'/γ_R)'(s) = ψ'(s/2)/4

    Returns: Im(γ_R'/γ_R(ρ₀)), Re(ψ'(ρ₀/2)/4)
    """
    s = mpmath.mpc(CENTER, gamma_0)
    # γ_R'/γ_R(s)
    gpg = -mpmath.log(mpmath.pi)/2 + mpmath.digamma(s/2)/2
    im_gpg = float(mpmath.im(gpg))
    # (γ_R'/γ_R)'(s) = ψ'(s/2)/4  (chain rule: d/ds ψ(s/2)/2 = ψ'(s/2)/4)
    psi1 = mpmath.psi(1, s/2)  # ψ'(s/2) = polygamma(1, s/2)
    re_correction = float(mpmath.re(psi1)) / 4.0
    return im_gpg, re_correction

print("[A(γ)] zero-sum + Gamma 보정...")
data = []
for i in range(len(zeros)):
    g0 = zeros[i]
    # L'/L zero-sum: same-sign + conjugate
    S1_same = sum(1.0/(g0-zeros[k]) for k in range(len(zeros)) if k!=i)
    H1_same = sum(1.0/(g0-zeros[k])**2 for k in range(len(zeros)) if k!=i)
    S1_conj = sum(1.0/(g0+zeros[k]) for k in range(len(zeros)))
    H1_conj = sum(1.0/(g0+zeros[k])**2 for k in range(len(zeros)))

    S1_L = S1_same + S1_conj
    H1_L = H1_same + H1_conj

    im_gpg, re_psi1_4 = gamma_corrections(g0)

    # A_Λ = (S1_L + im_gpg)² + 2*(H1_L + re_psi1_4)
    # Note: S₁ = -Im(c₀^L), but Im(γ'/γ) adds to c₀^Λ
    # c₀^Λ = c₀^L + γ'/γ, so Im(c₀^Λ) = Im(c₀^L) + Im(γ'/γ)
    # S₁_Λ = -Im(c₀^Λ) = S1_L - im_gpg (since -Im(c₀^L)=S1_L approximately)
    # Wait: Im(1/(ρ₀-ρₖ)) = -1/(γ₀-γₖ), so -Im(sum)=S1_same+S1_conj=S1_L
    # So c₀^L has Im(c₀^L) = -S1_L (up to constant)
    # c₀^Λ = c₀^L + γ'/γ → Im(c₀^Λ) = -S1_L + im_gpg
    # A_Λ = Im(c₀^Λ)² + 2*Re(c₁^Λ) = (-S1_L + im_gpg)² + 2*(H1_L + re_psi1_4)
    # = (S1_L - im_gpg)² + 2*(H1_L + re_psi1_4)

    S1_Lambda = S1_L - im_gpg
    H1_Lambda = H1_L + re_psi1_4
    A_Lambda = S1_Lambda**2 + 2*H1_Lambda

    # 비교: without Gamma correction
    A_L = S1_L**2 + 2*H1_L
    # 비교: same-only (original bad method)
    S1_same_corr = S1_same - im_gpg
    A_same = S1_same_corr**2 + 2*H1_same

    if A_Lambda > 0:
        data.append({
            't': g0,
            'A_Lambda': A_Lambda,
            'A_L': A_L,
            'A_same': A_same,
            'S1': S1_Lambda, 'H1': H1_Lambda
        })

print(f"  {len(data)} 유효")

# 중앙 60%
N = len(data)
lo = int(N * 0.2)
hi = int(N * 0.8)
valid = []
for idx in range(lo, hi):
    if idx <= 0 or idx >= N-1:
        continue
    d = data[idx]
    gap_r = data[idx+1]['t'] - d['t']
    gap_l = d['t'] - data[idx-1]['t']
    if gap_r <= 0 or gap_l <= 0:
        continue
    d_bar = 2.0 / (data[idx+1]['t'] - data[idx-1]['t'])
    d['gap_r_gue'] = gap_r * d_bar
    d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
    valid.append(d)

print(f"  내부 (중앙 60%): {len(valid)}")

A_lam = np.array([d['A_Lambda'] for d in valid])
A_l = np.array([d['A_L'] for d in valid])
A_s = np.array([d['A_same'] for d in valid])
gr = np.array([d['gap_r_gue'] for d in valid])
gm = np.array([d['gap_min_gue'] for d in valid])

sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')

for name, arr in [('A_Λ(full+Gamma)', A_lam), ('A_L(full,no Gamma)', A_l), ('A(same-only)', A_s)]:
    r1, p1 = stats.spearmanr(arr, gr)
    r2, p2 = stats.spearmanr(arr, gm)
    print(f"\n[{name}] (n={len(valid)})")
    print(f"  ρ(A, gap_right_GUE) = {r1:+.4f}  (p={p1:.3e})  {sig(p1)}")
    print(f"  ρ(A, gap_min_GUE)   = {r2:+.4f}  (p={p2:.3e})  {sig(p2)}")

print(f"\n[C-256 참조] ρ(A, gap_right) = -0.5898")
