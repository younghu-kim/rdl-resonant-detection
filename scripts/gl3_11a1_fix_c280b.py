#!/usr/bin/env python3
"""
[사이클 #280b] GL(3) Sym²(11a1) 단독 재시도
  이전 시도(c280): set_real_precision(100) → PARI 스택 오버플로우
  수정: precision=38 (C-264와 동일), 동일 방법론 유지

  Sym²(11a1): N=121, d=3, mu=[1,1,2], center=1.5, T=[3.0,35.0]
"""

import sys, os, math, time
import numpy as np
from scipy import stats
import mpmath

mpmath.mp.dps = 30

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(1000000000)  # 1GB
    pari.set_real_precision(38)   # C-264와 동일 — 100은 스택 오버플로우 유발
    print("cypari2 OK (1 GB, precision=38)")
except Exception as e:
    print(f"FATAL: {e}")
    sys.exit(1)

T_MAX = 35.0
T_MIN = 3.0
CENTER = 1.5
RESULT_PATH_11a1 = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl3_11a1_c280b.txt'
)
RESULT_COMBINED = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gl3_Abare_gap_c280.txt'
)


def pf(x):
    s = str(x).strip().replace(' E', 'e').replace('E ', 'e')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeros_sym2_11a1():
    """Sym²(11a1) 영점 수집."""
    pari("E_11a1 = ellinit([0,-1,1,-10,-10])")
    pari("L_sym2_11 = lfunsympow(E_11a1, 2)")
    pari(f"Li_sym2_11 = lfuninit(L_sym2_11, [0, {T_MAX + 2}])")
    pari(f"zv_sym2_11 = lfunzeros(Li_sym2_11, {T_MAX})")
    n = int(str(pari('#zv_sym2_11')))
    zeros = []
    for i in range(1, n + 1):
        t = pf(pari(f'zv_sym2_11[{i}]'))
        if not math.isnan(t) and t > 0.5:
            zeros.append(t)
    return sorted(zeros)


def gamma_smooth(gamma_0):
    """GL(3) Sym²: mu=[1,1,2], N=121"""
    s = mpmath.mpc(CENTER, gamma_0)
    mu_list = [1, 1, 2]
    N_cond = 121
    total = mpmath.mpc(0)
    for mu in mu_list:
        total += mpmath.digamma((s + mu) / 2) - mpmath.log(mpmath.pi)
    total /= 2
    total += mpmath.log(N_cond) / 2
    return float(mpmath.re(total)), float(mpmath.im(total))


def compute_A(zeros, idx, n_max=200):
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
    sm_re, sm_im = gamma_smooth(gamma_0)
    S1_corrected = S1_sum - sm_im
    A_bare = S1_sum ** 2 + 2 * H1_sum
    A_L = S1_corrected ** 2 + 2 * H1_sum
    return {
        'A_bare': A_bare, 'A_L': A_L,
        'S1_bare': S1_sum, 'H1': H1_sum,
        'S1_L': S1_corrected, 'smooth_im': sm_im,
    }


# ===== 메인 =====
print("=" * 70)
print("[사이클 #280b] Sym²(11a1) 단독 재시도 (precision=38)")
print(f"  T=[{T_MIN},{T_MAX}], center={CENTER}, mu=[1,1,2], N=121")
print("=" * 70)

t_start = time.time()

print("[영점]...")
try:
    all_zeros = get_zeros_sym2_11a1()
except Exception as e:
    print(f"FATAL 영점 수집 실패: {e}")
    import traceback; traceback.print_exc()
    sys.exit(1)

zeros = [z for z in all_zeros if z >= T_MIN]
print(f"  전체 {len(all_zeros)}개, t≥{T_MIN}: {len(zeros)}개")
if len(zeros) == 0:
    print("⚠️ 영점 0개 — 탐색 로직 점검 필요")
    sys.exit(1)
print(f"  t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
print(f"  dt_min = {min(b-a for a,b in zip(zeros,zeros[1:])):.4f}")
print(f"  소요: {time.time()-t_start:.1f}s")

# A(γ) 계산
print(f"\n[A(γ)] ({len(zeros)}개)...")
t0 = time.time()
data = []
fail_count = 0
for i in range(len(zeros)):
    try:
        all_idx = all_zeros.index(zeros[i])
        r = compute_A(all_zeros, all_idx)
        if r['A_bare'] > 0 and r['A_L'] > 0 and not math.isnan(r['A_bare']):
            r['t'] = zeros[i]
            data.append(r)
    except Exception as e:
        fail_count += 1
        print(f"  WARNING: idx={i}: {e}")
print(f"  완료: {len(data)}/{len(zeros)}, 실패: {fail_count} ({time.time()-t0:.1f}s)")

if fail_count > len(zeros) // 2:
    print("절반 이상 실패 — 중단")
    sys.exit(1)

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
    print("내부 영점 부족")
    sys.exit(1)

# Spearman
A_bare = np.array([d['A_bare'] for d in valid])
A_L = np.array([d['A_L'] for d in valid])
gr = np.array([d['gap_r_gue'] for d in valid])
gm = np.array([d['gap_min_gue'] for d in valid])

mask = np.isfinite(A_bare) & np.isfinite(A_L) & np.isfinite(gr) & np.isfinite(gm)
A_bare, A_L, gr, gm = A_bare[mask], A_L[mask], gr[mask], gm[mask]

rho_bare_min, p_bare_min = stats.spearmanr(A_bare, gm)
rho_bare_right, p_bare_right = stats.spearmanr(A_bare, gr)
rho_L_min, p_L_min = stats.spearmanr(A_L, gm)
rho_L_right, p_L_right = stats.spearmanr(A_L, gr)

sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')
print(f"\n[Spearman] (n={len(A_bare)})")
print(f"  ρ(A_bare, gap_min_GUE)   = {rho_bare_min:+.6f}  (p={p_bare_min:.4e})  {sig(p_bare_min)}")
print(f"  ρ(A_bare, gap_right_GUE) = {rho_bare_right:+.6f}  (p={p_bare_right:.4e})  {sig(p_bare_right)}")
print(f"  ρ(A_L,    gap_min_GUE)   = {rho_L_min:+.6f}  (p={p_L_min:.4e})  {sig(p_L_min)}")
print(f"  ρ(A_L,    gap_right_GUE) = {rho_L_right:+.6f}  (p={p_L_right:.4e})  {sig(p_L_right)}")

hf_bare = np.array([d['H1_frac_bare'] for d in valid])[mask]
hf_L = np.array([d['H1_frac_L'] for d in valid])[mask]
sm_arr = np.array([d['smooth_im'] for d in valid])[mask]
damping = 0.857 - abs(float(rho_bare_min))

print(f"\n[2H₁/A]")
print(f"  2H₁/A_bare = {np.nanmean(hf_bare):.6f}")
print(f"  2H₁/A_L    = {np.nanmean(hf_L):.6f}")
print(f"  <A_bare>    = {np.mean(A_bare):.4f}")
print(f"  <A_L>       = {np.mean(A_L):.4f}")
print(f"  <B_smooth>  = {np.mean(sm_arr):.6f}")
print(f"  arithmetic_damping = {damping:.4f}")
print(f"\n  총 소요: {time.time()-t_start:.1f}s")

# 결과 저장 (단독 파일)
result = {
    'label': 'Sym²(11a1) (N=121, d=3)',
    'n_zeros': len(zeros), 'n_inner': len(A_bare),
    'rho_bare_min': float(rho_bare_min), 'p_bare_min': float(p_bare_min),
    'rho_bare_right': float(rho_bare_right), 'p_bare_right': float(p_bare_right),
    'rho_L_min': float(rho_L_min), 'p_L_min': float(p_L_min),
    'rho_L_right': float(rho_L_right), 'p_L_right': float(p_L_right),
    'H1_frac_bare': float(np.nanmean(hf_bare)),
    'H1_frac_L': float(np.nanmean(hf_L)),
    'mean_A_bare': float(np.mean(A_bare)),
    'mean_A_L': float(np.mean(A_L)),
    'mean_Bsmooth': float(np.mean(sm_arr)),
    'mean_gap_min': float(np.mean(gm)),
    'arithmetic_damping': damping,
}

with open(RESULT_PATH_11a1, 'w') as f:
    f.write("# C-280b Sym²(11a1) A_bare gap_min Spearman (precision=38 fix)\n")
    f.write(f"n_zeros={result['n_zeros']}, n_inner={result['n_inner']}\n\n")
    f.write(f"rho_S(A_bare, gap_min_GUE)  = {result['rho_bare_min']:.6f}  p={result['p_bare_min']:.4e}\n")
    f.write(f"rho_S(A_bare, gap_right_GUE)= {result['rho_bare_right']:.6f}  p={result['p_bare_right']:.4e}\n")
    f.write(f"rho_S(A_L,    gap_min_GUE)  = {result['rho_L_min']:.6f}  p={result['p_L_min']:.4e}\n")
    f.write(f"rho_S(A_L,    gap_right_GUE)= {result['rho_L_right']:.6f}  p={result['p_L_right']:.4e}\n")
    f.write(f"2H1/A_bare = {result['H1_frac_bare']:.6f}\n")
    f.write(f"2H1/A_L    = {result['H1_frac_L']:.6f}\n")
    f.write(f"mean_A_bare = {result['mean_A_bare']:.4f}\n")
    f.write(f"mean_A_L    = {result['mean_A_L']:.4f}\n")
    f.write(f"mean_Bsmooth = {result['mean_Bsmooth']:.6f}\n")
    f.write(f"GUE_ref_rho_S = -0.857\n")
    f.write(f"arithmetic_damping = {result['arithmetic_damping']:.4f}\n")
    f.write(f"\n### 개별 데이터\n")
    f.write("t,A_bare,A_L,gap_min_gue,gap_right_gue,sm_im\n")
    for d in valid:
        f.write(f"{d['t']:.6f},{d['A_bare']:.6f},{d['A_L']:.6f},"
                f"{d['gap_min_gue']:.6f},{d['gap_r_gue']:.6f},{d['smooth_im']:.6f}\n")

print(f"\n저장: {RESULT_PATH_11a1}")

# 기존 통합 파일에 11a1 결과 추가 (append)
with open(RESULT_COMBINED, 'a') as f:
    f.write(f"\n## Sym²(11a1) (N=121, d=3) [C-280b 보완]\n")
    f.write(f"n_zeros={result['n_zeros']}, n_inner={result['n_inner']}\n\n")
    f.write(f"### 핵심 결과\n")
    f.write(f"rho_S(A_bare, gap_min_GUE)  = {result['rho_bare_min']:.6f}  p={result['p_bare_min']:.4e}\n")
    f.write(f"rho_S(A_bare, gap_right_GUE)= {result['rho_bare_right']:.6f}  p={result['p_bare_right']:.4e}\n")
    f.write(f"rho_S(A_L,    gap_min_GUE)  = {result['rho_L_min']:.6f}  p={result['p_L_min']:.4e}\n")
    f.write(f"rho_S(A_L,    gap_right_GUE)= {result['rho_L_right']:.6f}  p={result['p_L_right']:.4e}\n")
    f.write(f"2H1/A_bare = {result['H1_frac_bare']:.6f}\n")
    f.write(f"2H1/A_L    = {result['H1_frac_L']:.6f}\n")
    f.write(f"arithmetic_damping = {result['arithmetic_damping']:.4f}\n")

print(f"통합 파일에 11a1 결과 추가 완료.")
print("\n[완료]")
