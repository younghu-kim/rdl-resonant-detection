#!/usr/bin/env python3
"""
[사이클 #264] GL(3) A-gap 상관 + Prop 12 검증
  핵심 혁신: A(γ) = S₁² + 2H₁ 를 영점 위치만으로 계산
  S₁ = -Σ_{k≠0} 1/(γ₀-γₖ)  + Γ-smooth
  H₁ = Σ_{k≠0} 1/(γ₀-γₖ)²  + Γ'-smooth
  Γ-항은 smooth → Spearman ρ에 거의 무영향 → 영점 합만으로 충분

  PARI는 영점 수집(lfunzeros)에만 사용 — 초고속
  Γ 보정은 mpmath digamma로 추가 (정확도 향상)

체크리스트:
  [x] S₁, H₁ = 영점 합 + mpmath Gamma 보정
  [x] center=1.5
  [x] python -u
"""

import sys, os, time, math
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 30

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(512 * 10**6)
    pari.set_real_precision(38)
    print("cypari2 OK (512 MB)")
except Exception as e:
    print(f"FATAL: {e}")
    sys.exit(1)

T_MAX = 35.0
T_MIN = 3.0
CENTER = 1.5

RESULT_PATH = os.path.expanduser('~/Desktop/gdl_unified/results/gl3_A_gap_c264.txt')

CURVES = [
    {'name': 'sym2_11a1', 'coeffs': '[0,-1,1,-10,-10]', 'N_cond': 121,
     'mu': [1, 1, 2], 'label': 'Sym²(11a1) (N=121, d=3)'},
    {'name': 'sym2_37a1', 'coeffs': '[0,0,1,-1,0]', 'N_cond': 1369,
     'mu': [1, 1, 2], 'label': 'Sym²(37a1) (N=1369, d=3)'},
]


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeros(coeffs, t_max):
    """PARI lfunzeros로 GL(3) 영점 수집."""
    pari(f'E_tmp = ellinit({coeffs})')
    pari('L_tmp = lfunsympow(E_tmp, 2)')
    pari(f'Li_tmp = lfuninit(L_tmp, [0, {t_max + 2}])')
    pari(f'zv = lfunzeros(Li_tmp, {t_max})')
    n = int(str(pari('#zv')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv[{i}]'))
        if not math.isnan(t) and t > 0.5:
            zeros.append(t)
    return sorted(zeros)


def gamma_smooth(gamma_0, mu_list, N_cond):
    """
    Γ 보정: Λ'/Λ의 smooth part at ρ₀ = center + i·γ₀
    smooth = (1/2)log(N) + Γ_∞'/Γ_∞(ρ₀)
    Γ_∞'/Γ_∞(s) = (1/2) Σⱼ [ψ((s+μⱼ)/2) - log(π)]

    반환: (Re, Im) of smooth part
    """
    s = mpmath.mpc(CENTER, gamma_0)
    total = mpmath.mpc(0)
    for mu in mu_list:
        total += mpmath.digamma((s + mu) / 2) - mpmath.log(mpmath.pi)
    total /= 2
    total += mpmath.log(N_cond) / 2
    return float(mpmath.re(total)), float(mpmath.im(total))


def compute_A_from_zeros(zeros, idx, mu_list, N_cond, n_max=200):
    """
    영점 위치에서 직접 A(γ₀) 계산.

    c₀ = Σ_{k≠0} 1/(ρ₀ - ρₖ) + smooth
    ρ₀ - ρₖ = i(γ₀ - γₖ) (on-critical)
    1/(i·Δγ) = -i/Δγ → Re=0, Im=-1/Δγ

    So zero-sum c₀ = -i·S₁ where S₁ = Σ 1/(γ₀-γₖ)
    Re(c₀_zeros) = 0, Im(c₀_zeros) = -S₁

    c₁ = R'(ρ₀) = -Σ 1/(ρ₀-ρₖ)² + smooth'
    1/(i·Δγ)² = -1/Δγ² → Re = -1/Δγ², Im = 0
    R'(ρ₀)_zeros = Σ 1/Δγ² = H₁
    Re(c₁_zeros) = H₁

    A_zeros = Im(c₀)² + 2Re(c₁) = S₁² + 2H₁

    Γ 보정:
    c₀_total = c₀_zeros + c₀_smooth
    Im(c₀_total) = -S₁ + Im(smooth) → S₁_corrected = S₁ - Im(smooth)
    c₁_smooth = d/ds[smooth] at ρ₀ → 작고 smooth → 무시 (첫 근사)
    """
    gamma_0 = zeros[idx]
    n_zeros = len(zeros)

    # 영점 합: S₁ = Σ 1/(γ₀-γₖ), H₁ = Σ 1/(γ₀-γₖ)²
    S1_sum = 0.0
    H1_sum = 0.0
    for k in range(max(0, idx - n_max), min(n_zeros, idx + n_max + 1)):
        if k == idx:
            continue
        dg = gamma_0 - zeros[k]
        if abs(dg) < 1e-15:
            continue
        S1_sum += 1.0 / dg
        H1_sum += 1.0 / (dg * dg)

    # Γ smooth 보정
    sm_re, sm_im = gamma_smooth(gamma_0, mu_list, N_cond)
    # c₀ = -i·S₁ + (sm_re + i·sm_im)
    # Im(c₀) = -S₁ + sm_im
    # Re(c₀) = sm_re (이론: =0 on-critical, Thm 5) → 실제 ≈ (1/2)logN + Re(Γ'/Γ)
    #   → Re(c₀)≠0 here because smooth part includes (1/2)logN
    #   → but for A calculation: A = Im(c₀)² + 2Re(c₁)
    #   → Im(c₀) = -S₁ + sm_im, Re(c₁) = H₁ + Re(sm') ≈ H₁

    S1_corrected = S1_sum - sm_im
    H1_corrected = H1_sum  # c₁ smooth correction small, ignore
    A = S1_corrected ** 2 + 2 * H1_corrected

    return {
        'A': A, 'S1': S1_corrected, 'H1': H1_corrected,
        'S1_bare': S1_sum, 'H1_bare': H1_sum,
        'smooth_im': sm_im, 'smooth_re': sm_re,
    }


def analyze_curve(curve):
    name = curve['name']
    label = curve['label']
    mu = curve['mu']
    N = curve['N_cond']

    print("=" * 70)
    print(f"  {label}")
    print("=" * 70)

    # 영점 수집
    print("[영점]...")
    t0 = time.time()
    all_zeros = get_zeros(curve['coeffs'], T_MAX)
    zeros = [z for z in all_zeros if z >= T_MIN]
    print(f"  {len(zeros)}개 (t≥{T_MIN}), 전체 {len(all_zeros)}개")
    print(f"  t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
    print(f"  dt_min = {min(b-a for a,b in zip(zeros,zeros[1:])):.4f}")
    print(f"  소요: {time.time()-t0:.1f}s")

    if len(zeros) < 8:
        return None

    # A(γ) 계산 (영점 합 + Gamma 보정)
    print(f"\n[A(γ)] 영점 합 계산 ({len(zeros)}개)...")
    t0 = time.time()
    data = []
    for i in range(len(zeros)):
        # all_zeros에서의 실제 인덱스 사용 (가장자리 합 정확도)
        all_idx = all_zeros.index(zeros[i])
        r = compute_A_from_zeros(all_zeros, all_idx, mu, N)
        if r['A'] > 0 and not math.isnan(r['A']):
            r['t'] = zeros[i]
            data.append(r)
    print(f"  완료: {len(data)}/{len(zeros)} ({time.time()-t0:.1f}s)")

    if len(data) < 6:
        return None

    # 내부 영점 + 간격
    inner = data[2:-2]
    valid = []
    for d in inner:
        idx = data.index(d)
        if idx <= 0 or idx >= len(data) - 1:
            continue
        t_n = d['t']
        t_prev = data[idx-1]['t']
        t_next = data[idx+1]['t']
        gap_r = t_next - t_n
        gap_l = t_n - t_prev
        if gap_r <= 0 or gap_l <= 0:
            continue
        d_bar = 2.0 / (t_next - t_prev)
        d['gap_r_gue'] = gap_r * d_bar
        d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
        d['S1_plus_2overDR'] = d['S1'] + 2.0 / gap_r
        d['H1_frac'] = (2.0 * d['H1'] / d['A']) if d['A'] > 1e-10 else float('nan')
        valid.append(d)

    print(f"  내부: {len(valid)}")
    if len(valid) < 5:
        return None

    # Spearman
    A_arr = np.array([d['A'] for d in valid])
    gr = np.array([d['gap_r_gue'] for d in valid])
    gm = np.array([d['gap_min_gue'] for d in valid])

    rho_r, p_r = stats.spearmanr(A_arr, gr)
    rho_m, p_m = stats.spearmanr(A_arr, gm)
    rho_adj, p_adj = stats.spearmanr(A_arr[:-1], A_arr[1:]) if len(A_arr) > 2 else (float('nan'), float('nan'))

    # bare (no Gamma) 비교
    A_bare = np.array([d['S1_bare']**2 + 2*d['H1_bare'] for d in valid])
    rho_bare, p_bare = stats.spearmanr(A_bare, gr)

    sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')
    print(f"\n[Spearman] (n={len(valid)})")
    print(f"  ρ(A, gap_right_GUE)      = {rho_r:+.4f}  (p={p_r:.3e})  {sig(p_r)}")
    print(f"  ρ(A_bare, gap_right_GUE) = {rho_bare:+.4f}  (p={p_bare:.3e})  {sig(p_bare)}  [Γ 무보정]")
    print(f"  ρ(A, gap_min_GUE)        = {rho_m:+.4f}  (p={p_m:.3e})")
    print(f"  ρ(A_n, A_{{n+1}})          = {rho_adj:+.4f}  (p={p_adj:.3e})")

    # Prop 12
    s12 = np.array([d['S1_plus_2overDR'] for d in valid])
    hf = np.array([d['H1_frac'] for d in valid])
    n_pos = int(np.sum(s12 > 0))
    print(f"\n[Prop 12]")
    print(f"  S₁+2/Δ_R > 0: {n_pos}/{len(valid)} ({100*n_pos/len(valid):.0f}%)")
    print(f"  2H₁/A: {np.nanmean(hf):.3f} ± {np.nanstd(hf):.3f}")
    print(f"  <A>: {np.mean(A_arr):.2f}")

    verdict = '✅' if abs(rho_r) > 0.3 and p_r < 0.01 else ('⚠️' if abs(rho_r) > 0.3 and p_r < 0.05 else '❌')
    print(f"\n  [판정] {verdict}  부호={'음' if rho_r < 0 else '양'}")

    return {
        'label': label, 'N': N,
        'n_zeros': len(zeros), 'n_inner': len(valid),
        'rho_right_gue': float(rho_r), 'p_right_gue': float(p_r),
        'rho_bare_gue': float(rho_bare), 'p_bare_gue': float(p_bare),
        'rho_min_gue': float(rho_m), 'p_min_gue': float(p_m),
        'rho_adj': float(rho_adj), 'p_adj': float(p_adj),
        'verdict': verdict,
        'prop12_pos_frac': n_pos / len(valid),
        'prop12_H1_frac': float(np.nanmean(hf)),
        'mean_A': float(np.mean(A_arr)),
    }


# ===== 메인 =====
print("=" * 70)
print("[사이클 #264] GL(3) A-gap + Prop 12")
print(f"  영점 합 방식, T=[{T_MIN},{T_MAX}], center={CENTER}")
print("=" * 70)
print()

all_results = {}
for curve in CURVES:
    result = analyze_curve(curve)
    if result is not None:
        all_results[curve['name']] = result
    print()

# 비교표
print("=" * 70)
print("  전체 비교표: ρ(A, gap_right_GUE)")
print("=" * 70)

REF = [
    ('ζ(s)', 1, 1, 198, -0.5898, 6.1e-20),
    ('χ₃ (mod 3)', 1, 3, 72, -0.5733, 1.4e-7),
    ('χ₄ (mod 4)', 1, 4, 80, -0.5530, 1.0e-7),
    ('χ₅ (mod 5)', 1, 5, 46, -0.6334, 2.3e-6),
    ('E₁₁a₁', 2, 11, 43, -0.6275, 6.7e-6),
    ('E₃₇a₁', 2, 37, 17, -0.8676, 6.4e-6),
]

hdr = f"{'L-함수':<24} {'d':>2} {'N':>5} {'n':>5}  {'ρ':>10}  {'p':>12}  s"
print(hdr)
print("-" * 70)
for lbl, d, N, n, rho, p in REF:
    s = '✅' if abs(rho) > 0.3 and p < 0.01 else '❌'
    print(f"  {lbl:<22} {d:>2} {N:>5} {n:>5}  {rho:>10.4f}  {p:>12.3e}  {s}")
for nm, r in all_results.items():
    print(f"  {r['label']:<22} {3:>2} {r['N']:>5} {r['n_inner']:>5}  "
          f"{r['rho_right_gue']:>10.4f}  {r['p_right_gue']:>12.3e}  {r['verdict']}")

all_rhos = [r[4] for r in REF] + [r['rho_right_gue'] for r in all_results.values()]
rho_mean = np.mean(all_rhos)
rho_std = np.std(all_rhos)
all_neg = all(r < 0 for r in all_rhos)

if len(all_results) >= 1 and all_neg and all(abs(r) < 0.9 for r in all_rhos):
    vt = f"★★★★ degree-독립 보편성 (d=1,2,3): ρ = {rho_mean:.3f} ± {rho_std:.3f}"
elif len(all_results) >= 1 and all_neg:
    vt = f"★★★ GL(3) 확장: ρ = {rho_mean:.3f} ± {rho_std:.3f}"
else:
    vt = "조건부 또는 미완"
print(f"\n  최종: {vt}")

print("\n  Prop 12:")
for nm, r in all_results.items():
    print(f"    {r['label']}: S₁+2/Δ>0={r['prop12_pos_frac']*100:.0f}%, 2H₁/A={r['prop12_H1_frac']:.3f}")

# 저장
with open(RESULT_PATH, 'w') as f:
    f.write("=" * 70 + "\n")
    f.write("[사이클 #264] GL(3) A-gap + Prop 12\n")
    f.write(f"  영점 합 방식, T=[{T_MIN},{T_MAX}], center={CENTER}\n")
    f.write("=" * 70 + "\n\n")
    for nm, r in all_results.items():
        f.write(f"  {r['label']}:\n")
        f.write(f"    영점 {r['n_zeros']}, 내부 {r['n_inner']}\n")
        f.write(f"    ρ(A, gap_right_GUE)      = {r['rho_right_gue']:+.4f} (p={r['p_right_gue']:.3e})\n")
        f.write(f"    ρ(A_bare, gap_right_GUE) = {r['rho_bare_gue']:+.4f} (p={r['p_bare_gue']:.3e})\n")
        f.write(f"    ρ(A, gap_min_GUE)        = {r['rho_min_gue']:+.4f} (p={r['p_min_gue']:.3e})\n")
        f.write(f"    ρ(A_n, A_{{n+1}})          = {r['rho_adj']:+.4f} (p={r['p_adj']:.3e})\n")
        f.write(f"    Prop 12: S₁+2/Δ>0 = {r['prop12_pos_frac']*100:.0f}%\n")
        f.write(f"    Prop 12: 2H₁/A = {r['prop12_H1_frac']:.3f}\n")
        f.write(f"    <A> = {r['mean_A']:.2f}\n\n")
    f.write("=" * 70 + "\n  비교표\n" + "=" * 70 + "\n")
    f.write(hdr + "\n" + "-" * 70 + "\n")
    for lbl, d, N, n, rho, p in REF:
        s = '✅' if abs(rho) > 0.3 and p < 0.01 else '❌'
        f.write(f"  {lbl:<22} {d:>2} {N:>5} {n:>5}  {rho:>10.4f}  {p:>12.3e}  {s}\n")
    for nm, r in all_results.items():
        f.write(f"  {r['label']:<22} {3:>2} {r['N']:>5} {r['n_inner']:>5}  "
                f"{r['rho_right_gue']:>10.4f}  {r['p_right_gue']:>12.3e}  {r['verdict']}\n")
    f.write(f"\n  최종: {vt}\n")
    f.write(f"\n  Prop 12:\n")
    for nm, r in all_results.items():
        f.write(f"    {r['label']}: S₁+2/Δ>0={r['prop12_pos_frac']*100:.0f}%, 2H₁/A={r['prop12_H1_frac']:.3f}\n")
    f.write(f"\n  [방법]\n")
    f.write(f"  영점 합: S₁=Σ1/(γ₀-γₖ), H₁=Σ1/(γ₀-γₖ)², Γ보정=mpmath digamma\n")
    f.write(f"  A = (S₁+Im(Γ'))² + 2H₁ ≈ S₁² + 2H₁ (Γ' smooth)\n")
    f.write(f"  GUE: 국소 밀도 d̄=2/(t_{{n+1}}-t_{{n-1}})\n")
    f.write(f"  ρ(A_bare,gap)도 측정 → Γ 보정 영향 정량화\n")

print(f"\n저장: {RESULT_PATH}")
print("완료.")
