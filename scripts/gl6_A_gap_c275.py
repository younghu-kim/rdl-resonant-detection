#!/usr/bin/env python3
"""
[사이클 #275] GL(6) Sym⁵(11a1) A-gap 상관 + Prop 12 검증
  GL(5) C-273 스크립트 확장. degree 6까지 보편성 검증.

  Sym⁵(11a1): gammaV=[-2,-1,-1,0,0,1], k=6, N=161051, ε=-1
  center = k/2 = 3.0 (PARI 정규화)

  Gamma smooth: (1/2)Σⱼ[ψ((s+μⱼ)/2) - log(π)] + (1/2)log(N)
  μ = [-2, -1, -1, 0, 0, 1]

  B-38 경계: d≥5 계산 시간 경고. #206에서 76영점/49s → 실현 가능 확인.

  체크리스트:
    [x] S₁, H₁ = 영점 합 + mpmath Gamma 보정
    [x] center=3.0 (k=6)
    [x] python -u
    [x] mu = [-2, -1, -1, 0, 0, 1] (gammaV from #206)
    [x] 가장자리 절삭 (내부 영점만)
    [x] n_max=200 (충분한 이웃)
    [x] 메모리 6GB (d=6 확대)
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
    pari.allocatemem(6000000000)  # 6GB (GL(6) 필요)
    pari.set_real_precision(100)
    print("cypari2 OK (6 GB, precision=100)")
except Exception as e:
    print(f"FATAL: {e}")
    sys.exit(1)

T_MAX = 30.0   # #206에서 76영점 확인
T_MIN = 1.0    # 저-t 포함 (conductor 효과 관찰)
CENTER = 3.0   # k=6 → center=k/2=3.0

RESULT_PATH = os.path.expanduser('~/Desktop/gdl_unified/results/gl6_A_gap_c275.txt')


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeros_sym5(coeffs, t_max):
    """PARI lfunzeros로 GL(6) Sym⁵ 영점 수집."""
    pari(f'E_sym5 = ellinit({coeffs})')
    pari('L_sym5 = lfunsympow(E_sym5, 5)')
    pari(f'Li_sym5 = lfuninit(L_sym5, [0, {t_max + 2}])')
    pari(f'zv_sym5 = lfunzeros(Li_sym5, {t_max})')
    n = int(str(pari('#zv_sym5')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv_sym5[{i}]'))
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
    영점 위치에서 직접 A(γ₀) 계산 (zero-sum A_L).

    ρₖ = center + i·γₖ (on-critical)
    ρ₀ - ρₖ = i(γ₀ - γₖ)
    1/(i·Δγ) = -i/Δγ → Re=0, Im=-1/Δγ

    c₀_zeros: Im = -S₁
    c₁_zeros: Re = H₁ = Σ 1/Δγ²

    A = (S₁ - sm_im)² + 2·H₁
    """
    gamma_0 = zeros[idx]
    n_zeros = len(zeros)

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

    sm_re, sm_im = gamma_smooth(gamma_0, mu_list, N_cond)
    S1_corrected = S1_sum - sm_im
    H1_corrected = H1_sum
    A = S1_corrected ** 2 + 2 * H1_corrected

    return {
        'A': A, 'S1': S1_corrected, 'H1': H1_corrected,
        'S1_bare': S1_sum, 'H1_bare': H1_sum,
        'smooth_im': sm_im, 'smooth_re': sm_re,
    }


def analyze_sym5(curve):
    name = curve['name']
    label = curve['label']
    mu = curve['mu']
    N = curve['N_cond']
    degree = curve['degree']

    print("=" * 70)
    print(f"  {label}")
    print("=" * 70)

    # 영점 수집
    print("[영점]...")
    t0 = time.time()
    all_zeros = get_zeros_sym5(curve['coeffs'], T_MAX)
    zeros = [z for z in all_zeros if z >= T_MIN]
    elapsed = time.time() - t0
    print(f"  전체 {len(all_zeros)}개, t≥{T_MIN}: {len(zeros)}개")
    if len(zeros) < 2:
        print("  영점 부족 — 건너뜀")
        return None
    print(f"  t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
    print(f"  dt_min = {min(b-a for a,b in zip(zeros,zeros[1:])):.4f}")
    print(f"  소요: {elapsed:.1f}s")

    if len(zeros) < 8:
        print("  영점 부족 (<8) — 건너뜀")
        return None

    # A(γ) 계산
    print(f"\n[A(γ)] 영점 합 계산 ({len(zeros)}개)...")
    t0 = time.time()
    data = []
    for i in range(len(zeros)):
        all_idx = all_zeros.index(zeros[i])
        r = compute_A_from_zeros(all_zeros, all_idx, mu, N)
        if r['A'] > 0 and not math.isnan(r['A']):
            r['t'] = zeros[i]
            data.append(r)
    print(f"  완료: {len(data)}/{len(zeros)} ({time.time()-t0:.1f}s)")

    if len(data) < 6:
        print("  유효 데이터 부족 (<6) — 건너뜀")
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
        # GUE 정규화: 국소 밀도 d̄ = 2/(t_{n+1} - t_{n-1})
        d_bar = 2.0 / (t_next - t_prev)
        d['gap_r_gue'] = gap_r * d_bar
        d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
        d['S1_plus_2overDR'] = d['S1'] + 2.0 / gap_r
        d['H1_frac'] = (2.0 * d['H1'] / d['A']) if d['A'] > 1e-10 else float('nan')
        valid.append(d)

    print(f"  내부: {len(valid)}")
    if len(valid) < 5:
        print("  내부 영점 부족 (<5) — 건너뜀")
        return None

    # Spearman 상관
    A_arr = np.array([d['A'] for d in valid])
    gr = np.array([d['gap_r_gue'] for d in valid])
    gm = np.array([d['gap_min_gue'] for d in valid])

    rho_r, p_r = stats.spearmanr(A_arr, gr)
    rho_m, p_m = stats.spearmanr(A_arr, gm)
    rho_adj, p_adj = stats.spearmanr(A_arr[:-1], A_arr[1:]) if len(A_arr) > 2 else (float('nan'), float('nan'))

    # bare (Gamma 무보정) 비교
    A_bare = np.array([d['S1_bare']**2 + 2*d['H1_bare'] for d in valid])
    rho_bare_r, p_bare_r = stats.spearmanr(A_bare, gr)
    rho_bare_m, p_bare_m = stats.spearmanr(A_bare, gm)

    sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')
    print(f"\n[Spearman] (n={len(valid)})")
    print(f"  ρ(A_L, gap_right_GUE)      = {rho_r:+.4f}  (p={p_r:.3e})  {sig(p_r)}")
    print(f"  ρ(A_L, gap_min_GUE)        = {rho_m:+.4f}  (p={p_m:.3e})  {sig(p_m)}")
    print(f"  ρ(A_bare, gap_right_GUE)   = {rho_bare_r:+.4f}  (p={p_bare_r:.3e})  {sig(p_bare_r)}")
    print(f"  ρ(A_bare, gap_min_GUE)     = {rho_bare_m:+.4f}  (p={p_bare_m:.3e})  {sig(p_bare_m)}")
    print(f"  ρ(Aₙ, Aₙ₊₁)               = {rho_adj:+.4f}  (p={p_adj:.3e})")

    # Prop 12
    s12 = np.array([d['S1_plus_2overDR'] for d in valid])
    hf = np.array([d['H1_frac'] for d in valid])
    n_pos = int(np.sum(s12 > 0))
    print(f"\n[Prop 12]")
    print(f"  S₁+2/Δ_R > 0: {n_pos}/{len(valid)} ({100*n_pos/len(valid):.0f}%)")
    print(f"  2H₁/A: {np.nanmean(hf):.3f} ± {np.nanstd(hf):.3f}")
    print(f"  <A>: {np.mean(A_arr):.2f}")

    # 판정
    verdict_min = '✅' if abs(rho_m) > 0.3 and p_m < 0.01 else ('⚠️' if abs(rho_m) > 0.2 and p_m < 0.05 else '❌')
    verdict_r = '✅' if abs(rho_r) > 0.3 and p_r < 0.01 else ('⚠️' if abs(rho_r) > 0.2 and p_r < 0.05 else '❌')
    print(f"\n  [판정] gap_min: {verdict_min}  gap_right: {verdict_r}")

    return {
        'label': label, 'N': N, 'degree': degree,
        'n_zeros': len(zeros), 'n_inner': len(valid),
        'rho_right_gue': float(rho_r), 'p_right_gue': float(p_r),
        'rho_min_gue': float(rho_m), 'p_min_gue': float(p_m),
        'rho_bare_right': float(rho_bare_r), 'p_bare_right': float(p_bare_r),
        'rho_bare_min': float(rho_bare_m), 'p_bare_min': float(p_bare_m),
        'rho_adj': float(rho_adj), 'p_adj': float(p_adj),
        'verdict_min': verdict_min, 'verdict_right': verdict_r,
        'prop12_pos_frac': n_pos / len(valid),
        'prop12_H1_frac': float(np.nanmean(hf)),
        'mean_A': float(np.mean(A_arr)),
        'A_values': A_arr.tolist(),
        'gap_min_values': gm.tolist(),
    }


# ===== 메인 =====
print("=" * 70)
print("[사이클 #275] GL(6) Sym⁵ A-gap + Prop 12")
print(f"  영점 합 방식 (A_L), T=[{T_MIN},{T_MAX}], center={CENTER}")
print(f"  B-38 경계 탐사: d=6 계산 가능성 확인")
print("=" * 70)
print()

CURVES = [
    {'name': 'sym5_11a1',
     'coeffs': '[0,-1,1,-10,-20]',
     'N_cond': 161051,  # 11^5
     'mu': [-2, -1, -1, 0, 0, 1],
     'degree': 6,
     'label': 'Sym⁵(11a1) (N=161051, d=6)'},
]

all_results = {}
for curve in CURVES:
    try:
        result = analyze_sym5(curve)
        if result is not None:
            all_results[curve['name']] = result
    except Exception as e:
        print(f"  ⚠️ {curve['name']} 실패: {e}")
        import traceback; traceback.print_exc()
    print()

# 비교표 (GL(1)~GL(6) 전체)
print("=" * 70)
print("  전체 비교표: A-gap 보편성 (GL(1)~GL(6))")
print("=" * 70)

# 기존 결과
REF_A_LAMBDA = [
    ('ζ(s)', 1, 1, 198, -0.5898, 6.1e-20, 'gap_right'),
    ('χ₃ (mod 3)', 1, 3, 72, -0.5733, 1.4e-7, 'gap_right'),
    ('χ₄ (mod 4)', 1, 4, 80, -0.5530, 1.0e-7, 'gap_right'),
    ('χ₅ (mod 5)', 1, 5, 46, -0.6334, 2.3e-6, 'gap_right'),
    ('11a1', 2, 11, 92, -0.5688, 3.3e-9, 'gap_right'),
    ('37a1', 2, 37, 112, -0.5500, 3.4e-10, 'gap_right'),
]

REF_A_L = [
    ('Sym²(11a1)', 3, 121, 96, -0.42, 1.9e-5, 'gap_min'),
    ('Sym²(37a1)', 3, 1369, 48, -0.49, 3.8e-4, 'gap_min'),
    ('Sym³(11a1)', 4, 1331, 101, -0.5199, 2.5e-8, 'gap_min'),
    ('Sym³(37a1)', 4, 50653, 131, -0.5137, 3.5e-10, 'gap_min'),
    ('Sym⁴(11a1)', 5, 14641, 132, -0.4854, 3.7e-9, 'gap_min'),
]

hdr = f"{'L-함수':<24} {'d':>2} {'N':>8} {'n':>5}  {'gap':>9}  {'ρ':>8}  {'p':>12}  s"
print(hdr)
print("-" * 80)

for lbl, d, N, n, rho, p, gap in REF_A_LAMBDA:
    s = '✅' if abs(rho) > 0.3 and p < 0.01 else '❌'
    print(f"  {lbl:<22} {d:>2} {N:>8} {n:>5}  {gap:>9}  {rho:>8.4f}  {p:>12.3e}  {s}")

for lbl, d, N, n, rho, p, gap in REF_A_L:
    s = '✅' if abs(rho) > 0.3 and p < 0.01 else '❌'
    print(f"  {lbl:<22} {d:>2} {N:>8} {n:>5}  {gap:>9}  {rho:>8.4f}  {p:>12.3e}  {s}")

for nm, r in all_results.items():
    s_min = r['verdict_min']
    s_r = r['verdict_right']
    print(f"  {r['label']:<22} {r['degree']:>2} {r['N']:>8} {r['n_inner']:>5}  "
          f"{'gap_min':>9}  {r['rho_min_gue']:>8.4f}  {r['p_min_gue']:>12.3e}  {s_min}")
    print(f"  {'':22} {'':>2} {'':>8} {'':>5}  "
          f"{'gap_right':>9}  {r['rho_right_gue']:>8.4f}  {r['p_right_gue']:>12.3e}  {s_r}")

# 보편성 판정
n_total = len(REF_A_LAMBDA) + len(REF_A_L) + len(all_results)
all_neg = True
for r in REF_A_LAMBDA + REF_A_L:
    if r[4] >= 0:
        all_neg = False
for r in all_results.values():
    if r['rho_min_gue'] >= 0:
        all_neg = False

if all_results and all_neg:
    vt = f"★★★★ GL(1)~GL(6) degree-독립 보편성 ({n_total}/{n_total} 음)"
else:
    vt = "미완 또는 반례 발견"
print(f"\n  최종: {vt}")

# Prop 12 + 2H₁/A degree 경향
print("\n  Prop 12 + 2H₁/A:")
deg_h1a = {
    1: 0.867,   # ζ(s)
    2: 0.82,    # 11a1/37a1
    3: 0.748,   # Sym²
    4: 0.67,    # Sym³
    5: 0.661,   # Sym⁴ (C-273)
}
for d, h in sorted(deg_h1a.items()):
    print(f"    d={d}: 2H₁/A ≈ {h:.3f}")
for nm, r in all_results.items():
    print(f"    d={r['degree']}: 2H₁/A = {r['prop12_H1_frac']:.3f}  "
          f"S₁+2/Δ>0={r['prop12_pos_frac']*100:.0f}%  <A>={r['mean_A']:.2f}")
    print(f"    → degree 경향: {'감소 지속' if r['prop12_H1_frac'] < 0.661 else '포화 확정'}")

# 저장
with open(RESULT_PATH, 'w') as f:
    f.write("=" * 70 + "\n")
    f.write("[사이클 #275] GL(6) Sym⁵ A-gap + Prop 12\n")
    f.write(f"  영점 합 방식 (A_L), T=[{T_MIN},{T_MAX}], center={CENTER}\n")
    f.write(f"  대상: Sym⁵(11a1) d=6 N=161051\n")
    f.write(f"  B-38 경계 탐사: d=6 A-gap 실현 가능성 검증\n")
    f.write("=" * 70 + "\n\n")

    for nm, r in all_results.items():
        f.write(f"  {r['label']}:\n")
        f.write(f"    영점 {r['n_zeros']}, 내부 {r['n_inner']}\n")
        f.write(f"    ρ(A_L, gap_right_GUE) = {r['rho_right_gue']:+.4f} (p={r['p_right_gue']:.3e})\n")
        f.write(f"    ρ(A_L, gap_min_GUE)   = {r['rho_min_gue']:+.4f} (p={r['p_min_gue']:.3e})\n")
        f.write(f"    ρ(A_bare, gap_right)   = {r['rho_bare_right']:+.4f} (p={r['p_bare_right']:.3e})\n")
        f.write(f"    ρ(A_bare, gap_min)     = {r['rho_bare_min']:+.4f} (p={r['p_bare_min']:.3e})\n")
        f.write(f"    ρ(Aₙ, Aₙ₊₁)           = {r['rho_adj']:+.4f} (p={r['p_adj']:.3e})\n")
        f.write(f"    Prop 12: S₁+2/Δ>0={r['prop12_pos_frac']*100:.0f}%, "
                f"2H₁/A={r['prop12_H1_frac']:.3f}, <A>={r['mean_A']:.2f}\n")
        f.write(f"    판정: gap_min {r['verdict_min']}  gap_right {r['verdict_right']}\n")
        f.write(f"\n")

    f.write("\n전체 비교표 (GL(1)~GL(6)):\n")
    f.write(f"{hdr}\n")
    f.write("-" * 80 + "\n")
    for lbl, d, N, n, rho, p, gap in REF_A_LAMBDA:
        s = '✅' if abs(rho) > 0.3 and p < 0.01 else '❌'
        f.write(f"  {lbl:<22} {d:>2} {N:>8} {n:>5}  {gap:>9}  {rho:>8.4f}  {p:>12.3e}  {s}\n")
    for lbl, d, N, n, rho, p, gap in REF_A_L:
        s = '✅' if abs(rho) > 0.3 and p < 0.01 else '❌'
        f.write(f"  {lbl:<22} {d:>2} {N:>8} {n:>5}  {gap:>9}  {rho:>8.4f}  {p:>12.3e}  {s}\n")
    for nm, r in all_results.items():
        s_min = r['verdict_min']
        f.write(f"  {r['label']:<22} {r['degree']:>2} {r['N']:>8} {r['n_inner']:>5}  "
                f"{'gap_min':>9}  {r['rho_min_gue']:>8.4f}  {r['p_min_gue']:>12.3e}  {s_min}\n")

    f.write(f"\n최종: {vt}\n")

    # |ρ| degree 경향 기록
    f.write(f"\n|ρ| degree 경향:\n")
    f.write(f"  d=1 (GL(1)): |ρ| ≈ 0.59 (A_Λ, gap_right)\n")
    f.write(f"  d=2 (GL(2)): |ρ| ≈ 0.56 (A_Λ, gap_right)\n")
    f.write(f"  d=3 (GL(3)): |ρ| ≈ 0.46 (A_L, gap_min)\n")
    f.write(f"  d=4 (GL(4)): |ρ| ≈ 0.52 (A_L, gap_min)\n")
    f.write(f"  d=5 (GL(5)): |ρ| ≈ 0.49 (A_L, gap_min)\n")
    for nm, r in all_results.items():
        f.write(f"  d={r['degree']} (GL({r['degree']})): |ρ| = {abs(r['rho_min_gue']):.4f} (A_L, gap_min)\n")

print(f"\n결과 저장: {RESULT_PATH}")
print("완료.")
