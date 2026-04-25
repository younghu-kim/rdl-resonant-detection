#!/usr/bin/env python3
"""
[사이클 #280] GL(3) Sym²(11a1, 37a1) A_bare gap_min Spearman ρ_S 통합 측정
  목표: Paper 4 비교표 d=3 빈칸 채움 — C-278 ζ(s), C-279 GL(2), GL(4)~GL(6)와 동일 방법론

  설계:
    - Sym²(11a1): N=121, d=3, mu=[1,1,2], center=1.5
    - Sym²(37a1): N=1369, d=3, mu=[1,1,2], center=1.5
    - T=[3.0, 35.0] (C-264와 동일 범위)
    - A_bare = S₁² + 2H₁  (smooth 미보정)
    - A_L = (S₁ - sm_im)² + 2H₁  (Γ smooth 보정)
    - gap_min_GUE, gap_right_GUE: 국소 밀도 정규화
    - 4 Spearman 상관 (A_bare/A_L × gap_min/gap_right)

  GL(3) Γ-인자: Γ_R(s+1)·Γ_R(s+1)·Γ_R(s+2) → mu=[1,1,2]
  smooth = (1/2)Σⱼ[ψ((s+μⱼ)/2) - log(π)] + (1/2)log(N)
  at s = CENTER + i·γ = 1.5 + i·γ

  영점 수집: PARI lfunsympow(E, 2) → Sym² L-함수

  체크리스트:
    [x] CENTER=1.5, mu=[1,1,2], N_cond=121/1369
    [x] lfunsympow(E, 2) — GL(3) Sym² 영점
    [x] A_bare = S1_bare^2 + 2*H1_bare
    [x] A_L = (S1_bare - sm_im)^2 + 2*H1_bare
    [x] gap_min_GUE = min(gap_r, gap_l) * d_bar
    [x] Spearman (scipy.stats.spearmanr)
    [x] python -u
    [x] 2H1/A 비율 출력
    [x] NaN 필터
    [x] arithmetic_damping 계산
"""

import sys, os, math, time
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 30  # GL(3) T<35 → 30자리 충분

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(1000000000)  # 1GB
    pari.set_real_precision(100)
    print("cypari2 OK (1 GB, precision=100)")
except Exception as e:
    print(f"FATAL: {e}")
    sys.exit(1)

T_MAX = 35.0
T_MIN = 3.0
CENTER = 1.5  # GL(3) Sym² → center = (k+1)/2 ≈ 1.5

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl3_Abare_gap_c280.txt'
)


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeros_gl3_sym2(coeffs, t_max):
    """PARI lfunsympow(E, 2)로 Sym²(E) L-함수 영점 수집."""
    pari(f'E_gl3 = ellinit({coeffs})')
    pari('L_gl3 = lfunsympow(E_gl3, 2)')
    pari(f'Li_gl3 = lfuninit(L_gl3, [0, {t_max + 2}])')
    pari(f'zv_gl3 = lfunzeros(Li_gl3, {t_max})')
    n = int(str(pari('#zv_gl3')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv_gl3[{i}]'))
        if not math.isnan(t) and t > 0.5:
            zeros.append(t)
    return sorted(zeros)


def gamma_smooth(gamma_0, mu_list, N_cond):
    """
    Γ 보정: smooth part of Λ'/Λ at ρ₀ = center + i·γ₀
    smooth = (1/2)log(N) + (1/2) Σⱼ [ψ((s+μⱼ)/2) - log(π)]
    반환: (Re, Im)
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
    영점 위치에서 A(γ₀) 계산 — A_bare + A_L 모두 반환.

    S₁ = Σ_{k≠idx} 1/(γ₀ - γ_k)  (bare, same-side)
    H₁ = Σ_{k≠idx} 1/(γ₀ - γ_k)²  (bare)
    A_bare = S₁² + 2H₁
    A_L = (S₁ - sm_im)² + 2H₁  (Γ smooth 보정)
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
    A_bare = S1_sum ** 2 + 2 * H1_sum
    A_L = S1_corrected ** 2 + 2 * H1_sum

    return {
        'A_bare': A_bare, 'A_L': A_L,
        'S1_bare': S1_sum, 'H1': H1_sum,
        'S1_L': S1_corrected,
        'smooth_im': sm_im,
    }


def analyze_gl3(curve):
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
    all_zeros = get_zeros_gl3_sym2(curve['coeffs'], T_MAX)
    zeros = [z for z in all_zeros if z >= T_MIN]
    print(f"  전체 {len(all_zeros)}개, t≥{T_MIN}: {len(zeros)}개")
    if len(zeros) < 8:
        print("⚠️ 영점 0개 또는 부족 — 건너뜀")
        return None
    print(f"  t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
    print(f"  dt_min = {min(b-a for a,b in zip(zeros,zeros[1:])):.4f}")
    print(f"  소요: {time.time()-t0:.1f}s")

    # A(γ) 계산
    print(f"\n[A(γ)] 영점 합 계산 ({len(zeros)}개)...")
    t0 = time.time()
    data = []
    fail_count = 0
    for i in range(len(zeros)):
        try:
            all_idx = all_zeros.index(zeros[i])
            r = compute_A_from_zeros(all_zeros, all_idx, mu, N)
            if r['A_bare'] > 0 and r['A_L'] > 0 and not math.isnan(r['A_bare']) and not math.isnan(r['A_L']):
                r['t'] = zeros[i]
                data.append(r)
        except Exception as e:
            fail_count += 1
            print(f"  WARNING: idx={i} 실패: {e}")
    print(f"  완료: {len(data)}/{len(zeros)}, 실패: {fail_count} ({time.time()-t0:.1f}s)")

    if fail_count > len(zeros) // 2:
        print("  ⚠️ 절반 이상 실패 — 중단")
        return None

    if len(data) < 6:
        print("  유효 데이터 부족")
        return None

    # 내부 영점 + 간격
    inner = data[2:-2]
    valid = []
    for d in inner:
        idx = data.index(d)
        if idx <= 0 or idx >= len(data) - 1:
            continue
        t_n = d['t']
        t_prev = data[idx - 1]['t']
        t_next = data[idx + 1]['t']
        gap_r = t_next - t_n
        gap_l = t_n - t_prev
        if gap_r <= 0 or gap_l <= 0:
            continue
        d_bar = 2.0 / (t_next - t_prev)
        d['gap_r_gue'] = gap_r * d_bar
        d['gap_l_gue'] = gap_l * d_bar
        d['gap_min_gue'] = min(gap_r, gap_l) * d_bar
        d['H1_frac_bare'] = (2.0 * d['H1'] / d['A_bare']) if d['A_bare'] > 1e-10 else float('nan')
        d['H1_frac_L'] = (2.0 * d['H1'] / d['A_L']) if d['A_L'] > 1e-10 else float('nan')
        valid.append(d)

    print(f"  내부: {len(valid)}")
    if len(valid) < 5:
        print("  내부 영점 부족")
        return None

    # Spearman 상관 — 4가지 조합 (C-279와 동일)
    A_bare = np.array([d['A_bare'] for d in valid])
    A_L = np.array([d['A_L'] for d in valid])
    gr = np.array([d['gap_r_gue'] for d in valid])
    gm = np.array([d['gap_min_gue'] for d in valid])

    # NaN/Inf 필터
    mask = (np.isfinite(A_bare) & np.isfinite(A_L) &
            np.isfinite(gr) & np.isfinite(gm))
    A_bare, A_L, gr, gm = A_bare[mask], A_L[mask], gr[mask], gm[mask]

    if len(A_bare) < 5:
        print("  NaN 필터 후 데이터 부족")
        return None

    rho_bare_min, p_bare_min = stats.spearmanr(A_bare, gm)
    rho_bare_right, p_bare_right = stats.spearmanr(A_bare, gr)
    rho_L_min, p_L_min = stats.spearmanr(A_L, gm)
    rho_L_right, p_L_right = stats.spearmanr(A_L, gr)
    rho_adj, p_adj = stats.spearmanr(A_L[:-1], A_L[1:]) if len(A_L) > 2 else (float('nan'), float('nan'))

    sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')
    print(f"\n[Spearman] (n={len(A_bare)})")
    print(f"  ρ(A_bare, gap_min_GUE)   = {rho_bare_min:+.6f}  (p={p_bare_min:.4e})  {sig(p_bare_min)}")
    print(f"  ρ(A_bare, gap_right_GUE) = {rho_bare_right:+.6f}  (p={p_bare_right:.4e})  {sig(p_bare_right)}")
    print(f"  ρ(A_L,    gap_min_GUE)   = {rho_L_min:+.6f}  (p={p_L_min:.4e})  {sig(p_L_min)}")
    print(f"  ρ(A_L,    gap_right_GUE) = {rho_L_right:+.6f}  (p={p_L_right:.4e})  {sig(p_L_right)}")
    print(f"  ρ(Aₙ, Aₙ₊₁)            = {rho_adj:+.6f}  (p={p_adj:.4e})")

    # 2H₁/A
    hf_bare_arr = np.array([d['H1_frac_bare'] for d in valid])[mask]
    hf_L_arr = np.array([d['H1_frac_L'] for d in valid])[mask]
    sm_arr = np.array([d['smooth_im'] for d in valid])[mask]

    print(f"\n[2H₁/A]")
    print(f"  2H₁/A_bare = {np.nanmean(hf_bare_arr):.6f}")
    print(f"  2H₁/A_L    = {np.nanmean(hf_L_arr):.6f}")
    print(f"  <A_bare>    = {np.mean(A_bare):.4f}")
    print(f"  <A_L>       = {np.mean(A_L):.4f}")
    print(f"  <B_smooth>  = {np.mean(sm_arr):.6f}")

    damping = 0.857 - abs(float(rho_bare_min))
    print(f"\n  arithmetic_damping = {damping:.4f}")

    return {
        'label': label, 'name': name, 'N': N,
        'n_zeros': len(zeros), 'n_inner': len(A_bare),
        'rho_bare_min': float(rho_bare_min), 'p_bare_min': float(p_bare_min),
        'rho_bare_right': float(rho_bare_right), 'p_bare_right': float(p_bare_right),
        'rho_L_min': float(rho_L_min), 'p_L_min': float(p_L_min),
        'rho_L_right': float(rho_L_right), 'p_L_right': float(p_L_right),
        'rho_adj': float(rho_adj), 'p_adj': float(p_adj),
        'H1_frac_bare': float(np.nanmean(hf_bare_arr)),
        'H1_frac_L': float(np.nanmean(hf_L_arr)),
        'mean_A_bare': float(np.mean(A_bare)),
        'mean_A_L': float(np.mean(A_L)),
        'mean_Bsmooth': float(np.mean(sm_arr)),
        'mean_gap_min': float(np.mean(gm)),
        'arithmetic_damping': damping,
        'valid_data': [{
            't': d['t'],
            'A_bare': d['A_bare'], 'A_L': d['A_L'],
            'gap_min_gue': d['gap_min_gue'], 'gap_r_gue': d['gap_r_gue'],
            'sm_im': d['smooth_im'],
        } for d in valid],
    }


# ===== 메인 =====
print("=" * 70)
print("[사이클 #280] GL(3) Sym²(11a1, 37a1) A_bare gap_min Spearman ρ_S")
print(f"  T=[{T_MIN},{T_MAX}], center={CENTER}, mu=[1,1,2]")
print("=" * 70)
print()

CURVES = [
    {'name': 'sym2_11a1',
     'coeffs': '[0,-1,1,-10,-10]',
     'N_cond': 121,
     'mu': [1, 1, 2],
     'label': 'Sym²(11a1) (N=121, d=3)'},
    {'name': 'sym2_37a1',
     'coeffs': '[0,0,1,-1,0]',
     'N_cond': 1369,
     'mu': [1, 1, 2],
     'label': 'Sym²(37a1) (N=1369, d=3)'},
]

all_results = {}
t_total = time.time()
for curve in CURVES:
    try:
        result = analyze_gl3(curve)
        if result is not None:
            all_results[curve['name']] = result
    except Exception as e:
        print(f"  ⚠️ {curve['name']} 실패: {e}")
        import traceback
        traceback.print_exc()
    print()

print(f"총 소요: {time.time()-t_total:.1f}s")

# ===== 결과 파일 쓰기 =====
print("\n" + "=" * 70)
print("[결과 저장]")
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

with open(RESULT_PATH, 'w') as f:
    f.write("# C-280 GL(3) Sym²(11a1, 37a1) A_bare gap_min Spearman ρ_S 통합 측정\n")
    f.write(f"# T=[{T_MIN},{T_MAX}], CENTER={CENTER}, mu=[1,1,2]\n")
    f.write(f"# curves: {list(all_results.keys())}\n\n")

    for name, r in all_results.items():
        f.write(f"## {r['label']}\n")
        f.write(f"n_zeros={r['n_zeros']}, n_inner={r['n_inner']}\n\n")
        f.write(f"### 핵심 결과\n")
        f.write(f"rho_S(A_bare, gap_min_GUE)  = {r['rho_bare_min']:.6f}  p={r['p_bare_min']:.4e}\n")
        f.write(f"rho_S(A_bare, gap_right_GUE)= {r['rho_bare_right']:.6f}  p={r['p_bare_right']:.4e}\n")
        f.write(f"rho_S(A_L,    gap_min_GUE)  = {r['rho_L_min']:.6f}  p={r['p_L_min']:.4e}\n")
        f.write(f"rho_S(A_L,    gap_right_GUE)= {r['rho_L_right']:.6f}  p={r['p_L_right']:.4e}\n")
        f.write(f"\n### 보조 통계\n")
        f.write(f"2H1/A_bare = {r['H1_frac_bare']:.6f}\n")
        f.write(f"2H1/A_L    = {r['H1_frac_L']:.6f}\n")
        f.write(f"mean_A_bare = {r['mean_A_bare']:.4f}\n")
        f.write(f"mean_A_L    = {r['mean_A_L']:.4f}\n")
        f.write(f"mean_Bsmooth = {r['mean_Bsmooth']:.6f}\n")
        f.write(f"mean_gap_min = {r['mean_gap_min']:.6f}\n")
        f.write(f"GUE_ref_rho_S = -0.857\n")
        f.write(f"arithmetic_damping = {r['arithmetic_damping']:.4f}\n")

        # 개별 데이터
        f.write(f"\n### 개별 데이터\n")
        f.write("t,A_bare,A_L,gap_min_gue,gap_right_gue,sm_im\n")
        for d in r['valid_data']:
            f.write(f"{d['t']:.6f},{d['A_bare']:.6f},{d['A_L']:.6f},"
                    f"{d['gap_min_gue']:.6f},{d['gap_r_gue']:.6f},{d['sm_im']:.6f}\n")
        f.write("\n")

    # 통합 비교표 (d=1~6 + GUE)
    f.write("\n## 통합 비교표 (d=1~6 + GUE Spearman)\n")
    f.write("| d | L-함수 | ρ_S(A_L, gap_min) | 2H₁/A_L | n | damping |\n")
    f.write("|---|--------|-------------------|---------|---|--------|\n")
    f.write("| 1 | ζ(s) | -0.4371 | 0.666 | 337 | 0.42 |\n")
    f.write("| 2 | 11a1 | -0.4230 | 0.652 | 90 | 0.43 |\n")
    f.write("| 2 | 37a1 | -0.5250 | 0.675 | 110 | 0.33 |\n")
    for name, r in all_results.items():
        f.write(f"| 3 | {r['label'][:20]} | {r['rho_L_min']:.4f} "
                f"| {r['H1_frac_L']:.3f} | {r['n_inner']} "
                f"| {r['arithmetic_damping']:.2f} |\n")
    f.write("| 4 | Sym³(11a1) | -0.5199 | 0.661 | 101 | 0.34 |\n")
    f.write("| 4 | Sym³(37a1) | -0.5137 | 0.676 | 131 | 0.34 |\n")
    f.write("| 5 | Sym⁴(11a1) | -0.4854 | 0.661 | 132 | 0.37 |\n")
    f.write("| 6 | Sym⁵(11a1) | -0.4225 | 0.657 | 71 | 0.43 |\n")
    f.write("| GUE | N=50~2000 | -0.857 | 0.669 | — | 0.00 |\n")

print(f"  저장: {RESULT_PATH}")
print("\n[완료]")
