#!/usr/bin/env python3
"""
[사이클 #279] GL(2) 11a1 + 37a1 A_bare gap_min Spearman ρ_S 통합 측정
  목표: 비교표 d=2 빈칸 채움 — C-278 ζ(s), GL(4)~GL(6)와 동일한 방법론

  설계:
    - 11a1: N=11, ε=+1, mu=[0,1], center=1.0
    - 37a1: N=37, ε=-1, mu=[0,1], center=1.0
    - T_MAX=100 (GL(2)는 영점 풍부, 100 이내에 충분)
    - A_bare = S₁² + 2H₁  (smooth 미보정)
    - A_L = (S₁ - sm_im)² + 2H₁  (Γ smooth 보정)
    - gap_min_GUE, gap_right_GUE: 국소 밀도 정규화
    - 4 Spearman 상관 (A_bare/A_L × gap_min/gap_right)

  GL(2) Γ-인자: Γ_R(s)·Γ_R(s+1) → gammaV=[0,1]
  smooth = (1/2)Σⱼ[ψ((s+μⱼ)/2) - log(π)] + (1/2)log(N)
         = (1/2)[ψ((s)/2) + ψ((s+1)/2) - 2log(π)] + (1/2)log(N)
  at s = CENTER + i·γ = 1 + i·γ

  체크리스트:
    [x] CENTER=1.0, mu=[0,1], gammaV 일치
    [x] A_bare = S1_bare^2 + 2*H1_bare
    [x] A_L = (S1_bare - sm_im)^2 + 2*H1
    [x] gap_min_GUE = min(gap_r, gap_l) * d_bar
    [x] Spearman (scipy.stats.spearmanr)
    [x] python -u
    [x] 2H1/A 비율 출력
    [x] NaN 필터
"""

import sys, os, math, time
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 30  # GL(2) T<100 → 30자리 충분

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

T_MAX = 100.0
T_MIN = 2.0
CENTER = 1.0  # GL(2) weight 2 → center = k/2 = 1

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl2_Abare_gap_c279.txt'
)


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeros_gl2(coeffs, t_max):
    """PARI lfunzeros로 GL(2) 타원곡선 L-함수 영점 수집."""
    pari(f'E_gl2 = ellinit({coeffs})')
    pari(f'Li_gl2 = lfuninit(E_gl2, [0, {t_max + 2}])')
    pari(f'zv_gl2 = lfunzeros(Li_gl2, {t_max})')
    n = int(str(pari('#zv_gl2')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv_gl2[{i}]'))
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
    영점 위치에서 A(γ₀) 계산.
    ρₖ = center + i·γₖ (on-critical)
    S₁ = Σ_{k≠idx} 1/(γ₀ - γ_k)  (bare)
    H₁ = Σ_{k≠idx} 1/(γ₀ - γ_k)²  (bare)
    A_bare = S₁² + 2H₁
    A_L = (S₁ - sm_im)² + 2H₁
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


def analyze_gl2(curve):
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
    all_zeros = get_zeros_gl2(curve['coeffs'], T_MAX)
    zeros = [z for z in all_zeros if z >= T_MIN]
    print(f"  전체 {len(all_zeros)}개, t≥{T_MIN}: {len(zeros)}개")
    if len(zeros) < 8:
        print("  영점 부족 — 건너뜀")
        return None
    print(f"  t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
    print(f"  dt_min = {min(b-a for a,b in zip(zeros,zeros[1:])):.4f}")
    print(f"  소요: {time.time()-t0:.1f}s")

    # A(γ) 계산
    print(f"\n[A(γ)] 영점 합 계산 ({len(zeros)}개)...")
    t0 = time.time()
    data = []
    for i in range(len(zeros)):
        all_idx = all_zeros.index(zeros[i])
        r = compute_A_from_zeros(all_zeros, all_idx, mu, N)
        if r['A_bare'] > 0 and r['A_L'] > 0:
            r['t'] = zeros[i]
            data.append(r)
    print(f"  완료: {len(data)}/{len(zeros)} ({time.time()-t0:.1f}s)")

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

    # Spearman 상관 — 4가지 조합
    A_bare = np.array([d['A_bare'] for d in valid])
    A_L = np.array([d['A_L'] for d in valid])
    gr = np.array([d['gap_r_gue'] for d in valid])
    gm = np.array([d['gap_min_gue'] for d in valid])

    rho_bare_min, p_bare_min = stats.spearmanr(A_bare, gm)
    rho_bare_right, p_bare_right = stats.spearmanr(A_bare, gr)
    rho_L_min, p_L_min = stats.spearmanr(A_L, gm)
    rho_L_right, p_L_right = stats.spearmanr(A_L, gr)
    rho_adj, p_adj = stats.spearmanr(A_L[:-1], A_L[1:]) if len(A_L) > 2 else (float('nan'), float('nan'))

    sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')
    print(f"\n[Spearman] (n={len(valid)})")
    print(f"  ρ(A_bare, gap_min_GUE)   = {rho_bare_min:+.6f}  (p={p_bare_min:.4e})  {sig(p_bare_min)}")
    print(f"  ρ(A_bare, gap_right_GUE) = {rho_bare_right:+.6f}  (p={p_bare_right:.4e})  {sig(p_bare_right)}")
    print(f"  ρ(A_L,    gap_min_GUE)   = {rho_L_min:+.6f}  (p={p_L_min:.4e})  {sig(p_L_min)}")
    print(f"  ρ(A_L,    gap_right_GUE) = {rho_L_right:+.6f}  (p={p_L_right:.4e})  {sig(p_L_right)}")
    print(f"  ρ(Aₙ, Aₙ₊₁)            = {rho_adj:+.6f}  (p={p_adj:.4e})")

    # 2H₁/A
    hf_bare = np.array([d['H1_frac_bare'] for d in valid])
    hf_L = np.array([d['H1_frac_L'] for d in valid])
    print(f"\n[2H₁/A]")
    print(f"  2H₁/A_bare = {np.nanmean(hf_bare):.6f}")
    print(f"  2H₁/A_L    = {np.nanmean(hf_L):.6f}")
    print(f"  <A_bare>    = {np.mean(A_bare):.4f}")
    print(f"  <A_L>       = {np.mean(A_L):.4f}")
    sm_arr = np.array([d['smooth_im'] for d in valid])
    print(f"  <B_smooth>  = {np.mean(sm_arr):.6f}")

    return {
        'label': label, 'name': name, 'N': N,
        'n_zeros': len(zeros), 'n_inner': len(valid),
        'rho_bare_min': float(rho_bare_min), 'p_bare_min': float(p_bare_min),
        'rho_bare_right': float(rho_bare_right), 'p_bare_right': float(p_bare_right),
        'rho_L_min': float(rho_L_min), 'p_L_min': float(p_L_min),
        'rho_L_right': float(rho_L_right), 'p_L_right': float(p_L_right),
        'rho_adj': float(rho_adj), 'p_adj': float(p_adj),
        'H1_frac_bare': float(np.nanmean(hf_bare)),
        'H1_frac_L': float(np.nanmean(hf_L)),
        'mean_A_bare': float(np.mean(A_bare)),
        'mean_A_L': float(np.mean(A_L)),
        'mean_Bsmooth': float(np.mean(sm_arr)),
        'mean_gap_min': float(np.mean(gm)),
        'valid_data': [{
            't': d['t'],
            'A_bare': d['A_bare'], 'A_L': d['A_L'],
            'gap_min_gue': d['gap_min_gue'], 'gap_r_gue': d['gap_r_gue'],
            'sm_im': d['smooth_im'],
        } for d in valid],
    }


# ===== 메인 =====
print("=" * 70)
print("[사이클 #279] GL(2) A_bare gap_min 통합 측정")
print(f"  T=[{T_MIN},{T_MAX}], center={CENTER}, mu=[0,1]")
print("=" * 70)
print()

CURVES = [
    {'name': 'gl2_11a1',
     'coeffs': '[0,-1,1,-10,-20]',
     'N_cond': 11,
     'mu': [0, 1],
     'label': 'L(s, 11a1) (N=11, d=2, ε=+1)'},
    {'name': 'gl2_37a1',
     'coeffs': '[0,0,1,-1,0]',
     'N_cond': 37,
     'mu': [0, 1],
     'label': 'L(s, 37a1) (N=37, d=2, ε=-1)'},
]

all_results = {}
for curve in CURVES:
    try:
        result = analyze_gl2(curve)
        if result is not None:
            all_results[curve['name']] = result
    except Exception as e:
        print(f"  ⚠️ {curve['name']} 실패: {e}")
        import traceback
        traceback.print_exc()
    print()

# ===== 결과 파일 쓰기 =====
print("\n" + "=" * 70)
print("[결과 저장]")
with open(RESULT_PATH, 'w') as f:
    f.write("# C-279 GL(2) A_bare gap_min Spearman ρ_S 통합 측정\n")
    f.write(f"# T=[{T_MIN},{T_MAX}], CENTER={CENTER}, mu=[0,1]\n")
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
        f.write(f"arithmetic_damping = {0.857 - abs(r['rho_bare_min']):.4f}\n")

        # 개별 데이터
        f.write(f"\n### 개별 데이터\n")
        f.write("t,A_bare,A_L,gap_min_gue,gap_right_gue,sm_im\n")
        for d in r['valid_data']:
            f.write(f"{d['t']:.6f},{d['A_bare']:.6f},{d['A_L']:.6f},"
                    f"{d['gap_min_gue']:.6f},{d['gap_r_gue']:.6f},{d['sm_im']:.6f}\n")
        f.write("\n")

    # 요약표
    f.write("\n## 통합 비교표 (d=2 + 기존)\n")
    f.write("| d | L-함수 | A_bare gap_min ρ_S | A_L gap_min ρ_S | 2H₁/A_L | n |\n")
    f.write("|---|--------|-------------------|----------------|---------|---|\n")
    for name, r in all_results.items():
        f.write(f"| 2 | {r['label'][:20]} | {r['rho_bare_min']:.4f} "
                f"| {r['rho_L_min']:.4f} | {r['H1_frac_L']:.3f} | {r['n_inner']} |\n")
    f.write("| 1 | ζ(s) | -0.4484 | -0.4371 | 0.666 | 337 |\n")
    f.write("| 4 | Sym³(11a1) | -0.5660 | -0.5199 | 0.661 | 101 |\n")
    f.write("| 4 | Sym³(37a1) | -0.4900 | -0.5137 | 0.676 | 131 |\n")
    f.write("| 5 | Sym⁴(11a1) | -0.4660 | -0.4854 | 0.661 | 132 |\n")
    f.write("| 6 | Sym⁵(11a1) | -0.5864 | -0.4225 | 0.657 | 71 |\n")
    f.write("| GUE | N=50~2000 | -0.857 | — | 0.669 | — |\n")

print(f"  저장: {RESULT_PATH}")
print("\n[완료]")
