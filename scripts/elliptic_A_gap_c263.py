#!/usr/bin/env python3
"""
[사이클 #263] GL(2) A-gap 상관 교차검증
  타원곡선 L-함수: 11a1 (N=11, rank 0) + 37a1 (N=37, rank 1)
  PARI/GP cypari2 인터페이스 사용

체크리스트 검증:
  [x] Cauchy 적분: lfunlambda 사용 (GL(2) 완성 함수)
  [x] 영점: lfunzeros → s=1+it (number-theoretic 정규화)
  [x] GUE 정규화: gap × (1/π) × log(√N × T/(2π))
  [x] NaN/Inf 체크
  [x] 에러 처리
  [x] 결과 파일에 설정·시드별 결과·통계·판정 포함
"""

import sys
import os
import time
import math
import numpy as np
from scipy import stats

# cypari2: ~/.local/lib에 설치됨
sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
try:
    import cypari2
    pari = cypari2.Pari()
    pari.set_real_precision(50)
    print("cypari2 OK (PARI 50 dps)")
except Exception as e:
    print(f"FATAL: cypari2 로드 실패: {e}")
    sys.exit(1)

# ===== 파라미터 =====
T_MIN = 1.0       # 영점 탐색 하한 (t=0 제외)
T_MAX = 60.0      # 영점 탐색 상한
CAUCHY_R = 0.01   # Cauchy 적분 반지름
N_PTS = 64        # Cauchy 적분 점 수
H_DEC = '0.00001' # 수치 미분 간격 (문자열, PARI 파싱용)
ABS_THR = '1e-20' # Λ 절대값 임계치 (이 이하면 skip)

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/elliptic_A_gap_c263.txt'
)

CURVES = [
    {
        'name': '11a1',
        'coeffs': '[0,-1,1,-10,-10]',
        'N': 11,
        'rank': 0,
        'label': 'E₁₁a₁ (N=11, rank 0, ε=+1)',
        'pari_var': 'L11',
    },
    {
        'name': '37a1',
        'coeffs': '[0,0,1,-1,0]',
        'N': 37,
        'rank': 1,
        'label': 'E₃₇a₁ (N=37, rank 1, ε=-1)',
        'pari_var': 'L37',
    },
]

# ===== 헬퍼 함수 =====

def pari_to_float(x):
    """cypari2 객체 → Python float. 'N E-5' 형식 처리."""
    s = str(x).strip()
    # PARI: '1.234 E-5' 형식 → '1.234e-5'
    s = s.replace(' E', 'e').replace('E ', 'e')
    # 'None' 처리
    if s == 'None' or s == '':
        return float('nan')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def get_zeros_list(coeffs_str, t_min, t_max):
    """타원곡선 계수로부터 영점 t 값 목록 반환."""
    z_pari = pari(f'lfunzeros(lfuncreate(ellinit({coeffs_str})), [{t_min}, {t_max}])')
    n = int(str(pari(f'#lfunzeros(lfuncreate(ellinit({coeffs_str})), [{t_min}, {t_max}])')))
    zeros = []
    for i in range(1, n + 1):
        t = pari_to_float(pari(f'({z_pari})[{i}]'))
        if not math.isnan(t) and t > 0.5:  # t=0 제외
            zeros.append(t)
    return sorted(zeros)


def compute_A_cauchy(t_zero, lf_var, r=CAUCHY_R, n_pts=N_PTS, h_dec=H_DEC):
    """
    Cauchy 적분으로 Laurent 계수 c₀, c₁ 계산.
    Λ'/Λ(ρ+u) = 1/u + c₀ + c₁u + ...
    A(γ) = Im(c₀)² + 2Re(c₁)
    반환: (A, Re(c₀), Im(c₀), Re(c₁)) 또는 None (실패 시)
    """
    script = (
        f'rho=1+{t_zero:.15f}*I; c0s=0+0.*I; c1s=0+0.*I; cnt=0; '
        f'for(k=0,{n_pts-1}, '
        f'th=2*Pi*k/{n_pts}; uu={r}*exp(th*I); ss=rho+uu; '
        f'Lv=lfunlambda({lf_var},ss); '
        f'if(abs(Lv)<{ABS_THR},next); '
        f'Lp=lfunlambda({lf_var},ss+{h_dec}); '
        f'Lm=lfunlambda({lf_var},ss-{h_dec}); '
        f'LpL=(Lp-Lm)/(2*{h_dec}*Lv); '
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
            return None
        if n_cnt < n_pts // 2:  # 절반 이상 실패
            return None
        A = c0_im**2 + 2 * c1_re
        if math.isnan(A) or math.isinf(A):
            return None
        return (A, c0_re, c0_im, c1_re)
    except Exception as e:
        return None


def gue_normalization(gap, t_mid, N_cond):
    """
    GUE 정규화: gap × (1/π) × log(√N × T/(2π))
    """
    log_factor = math.log(math.sqrt(N_cond) * t_mid / (2 * math.pi))
    return gap * log_factor / math.pi


def spearman_with_p(x, y):
    """Spearman 상관 + p값 반환."""
    if len(x) < 5:
        return float('nan'), float('nan')
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p)


# ===== 메인 분석 =====

print("=" * 70)
print("[사이클 #263] GL(2) A-gap 상관 교차검증")
print(f"  T=[{T_MIN},{T_MAX}], CAUCHY_R={CAUCHY_R}, N_PTS={N_PTS}")
print("=" * 70)
print()

# PARI L-함수 객체 초기화
print("[초기화] PARI L-함수 객체 생성...")
for curve in CURVES:
    lv = curve['pari_var']
    coeffs = curve['coeffs']
    pari(f'E_tmp = ellinit({coeffs}); {lv} = lfuncreate(E_tmp)')
    print(f"  {curve['label']}: {lv} OK")
print()

results_all = {}

for curve in CURVES:
    name = curve['name']
    label = curve['label']
    N_cond = curve['N']
    lv = curve['pari_var']
    coeffs = curve['coeffs']

    print("=" * 70)
    print(f"  {label}")
    print("=" * 70)

    # --- Step 0: 영점 탐색 ---
    print(f"[Step 0] 영점 탐색 t∈[{T_MIN},{T_MAX}]...")
    t0 = time.time()
    zeros = get_zeros_list(coeffs, T_MIN, T_MAX)
    print(f"  영점 {len(zeros)}개  t∈[{zeros[0]:.3f}, {zeros[-1]:.3f}]")
    print(f"  최소 간격: {min(b-a for a,b in zip(zeros,zeros[1:])):.4f}")
    print(f"  소요: {time.time()-t0:.1f}s")

    if len(zeros) < 5:
        print(f"  ⚠️ 영점 {len(zeros)}개 — 부족, 스킵")
        continue

    # --- Step 1: A(γ) 계산 ---
    print(f"[Step 1] A(γ) 계산 ({len(zeros)}개 영점)...")
    A_list = []
    t_list = []
    c0re_list = []
    c0im_list = []
    n_fail = 0
    t0 = time.time()

    for i, tz in enumerate(zeros):
        result = compute_A_cauchy(tz, lv)
        if result is None:
            n_fail += 1
            print(f"  ⚠️ 영점 #{i+1} (t={tz:.3f}) — 적분 실패")
        else:
            A, c0re, c0im, c1re = result
            A_list.append(A)
            t_list.append(tz)
            c0re_list.append(c0re)
            c0im_list.append(c0im)

        if (i + 1) % 10 == 0:
            print(f"  [{time.time()-t0:.0f}s] {i+1}/{len(zeros)}  t={tz:.3f}")

    print(f"  계산 완료: {len(A_list)}/{len(zeros)} 유효  ({n_fail} 오류)")
    print(f"  |Im(c₀)| 평균: {np.mean(np.abs(c0im_list)):.4f}")
    print(f"  |Re(c₀)| 평균: {np.mean(np.abs(c0re_list)):.4e}  (비율: {np.mean(np.abs(c0re_list))/max(1e-10,np.mean(np.abs(c0im_list))):.3f})")

    if n_fail > len(zeros) * 0.5:
        print(f"  ⚠️ 실패율 {n_fail}/{len(zeros)} > 50% — 결과 신뢰도 낮음")

    # --- Step 2: 내부 영점 선택 (edge 2개씩 제외) ---
    n_total = len(t_list)
    inner_idx = list(range(2, n_total - 2))
    A_inner = [A_list[i] for i in inner_idx]
    t_inner = [t_list[i] for i in inner_idx]
    print(f"  내부 영점: {len(A_inner)}/{n_total} 유효")

    # --- Step 3: 간격 계산 ---
    # 내부 영점 n에 대해:
    # gap_right[n] = t[n+1] - t[n]  (다음과의 간격)
    # gap_left[n] = t[n] - t[n-1]   (이전과의 간격)
    # gap_min = min(left, right)
    if len(inner_idx) < 5:
        print(f"  ⚠️ 내부 영점 {len(inner_idx)}개 — 부족, 상관 계산 불가")
        continue

    gap_right_raw = []
    gap_left_raw = []
    gap_right_gue = []
    gap_left_gue = []
    gap_min_gue = []

    for pos_in_inner, orig_i in enumerate(inner_idx):
        t_n = t_list[orig_i]
        t_next = t_list[orig_i + 1]
        t_prev = t_list[orig_i - 1]
        gap_r = t_next - t_n
        gap_l = t_n - t_prev
        t_mid_r = (t_n + t_next) / 2
        t_mid_l = (t_prev + t_n) / 2

        gap_right_raw.append(gap_r)
        gap_left_raw.append(gap_l)

        # GUE 정규화
        gue_r = gue_normalization(gap_r, t_mid_r, N_cond)
        gue_l = gue_normalization(gap_l, t_mid_l, N_cond)
        gap_right_gue.append(gue_r)
        gap_left_gue.append(gue_l)
        gap_min_gue.append(min(gue_r, gue_l))

    # --- Step 4: Spearman 상관 ---
    A_arr = np.array(A_inner)
    gr_arr = np.array(gap_right_gue)
    gl_arr = np.array(gap_left_gue)
    gm_arr = np.array(gap_min_gue)
    gr_raw_arr = np.array(gap_right_raw)

    print(f"\n[Step 4] Spearman 상관 (n={len(A_inner)})...")

    rho_right, p_right = spearman_with_p(A_arr, gr_arr)
    rho_left, p_left = spearman_with_p(A_arr, gl_arr)
    rho_min, p_min = spearman_with_p(A_arr, gm_arr)
    rho_raw, p_raw = spearman_with_p(A_arr, gr_raw_arr)
    # 인접 A 상관
    if len(A_arr) > 1:
        rho_adj, p_adj = spearman_with_p(A_arr[:-1], A_arr[1:])
    else:
        rho_adj, p_adj = float('nan'), float('nan')

    def sig(p):
        return '✅ 유의(p<0.01)' if p < 0.01 else ('⚠️ 유의(p<0.05)' if p < 0.05 else '❌ 비유의')

    print(f"  ρ(A, gap_right_raw)  = {rho_raw:+.4f}  (p={p_raw:.3e})  {sig(p_raw)}")
    print(f"  ρ(A, gap_right_GUE) = {rho_right:+.4f}  (p={p_right:.3e})  {sig(p_right)}")
    print(f"  ρ(A, gap_left_GUE)  = {rho_left:+.4f}  (p={p_left:.3e})  {sig(p_left)}")
    print(f"  ρ(A, gap_min_GUE)   = {rho_min:+.4f}  (p={p_min:.3e})  {sig(p_min)}")
    print(f"  ρ(A_n, A_{{n+1}})     = {rho_adj:+.4f}  (p={p_adj:.3e})  {sig(p_adj)}")

    verdict_flag = '✅' if (abs(rho_right) > 0.3 and p_right < 0.01) else '❌'
    sign_match = '음(GL(1)과 일치)' if rho_right < 0 else '양(GL(1)과 반대)'
    print(f"\n  [판정] ρ(A,gap_right_GUE)={rho_right:.4f}  {verdict_flag}  부호={sign_match}")

    results_all[name] = {
        'label': label,
        'N': N_cond,
        'n_zeros': len(zeros),
        'n_inner': len(A_inner),
        'rho_right_gue': rho_right,
        'p_right_gue': p_right,
        'rho_min_gue': rho_min,
        'p_min_gue': p_min,
        'rho_adj': rho_adj,
        'p_adj': p_adj,
        'verdict_flag': verdict_flag,
    }

# ===== 비교표 + 판정 =====
print()
print("=" * 70)
print("  비교표: ρ(A, gap_right_GUE) 요약")
print("=" * 70)
print(f"{'L-함수':<30} {'q':>4} {'n내부':>6}  {'ρ(right_GUE)':>14}  {'p값':>12}  {'유의':>4}")
print("-" * 80)

# GL(1) 참조값
GL1_REF = [
    ('ζ(s) [C-256]', 1, 198, -0.5898, 6.129e-20),
    ('χ₃ (mod 3)', 3, 72, -0.5733, 1.413e-7),
    ('χ₄ (mod 4)', 4, 80, -0.5530, 1.039e-7),
    ('χ₅ (mod 5)', 5, 46, -0.6334, 2.305e-6),
]
for lbl, q, n, rho, p in GL1_REF:
    flag = '✅' if abs(rho) > 0.3 and p < 0.01 else '❌'
    print(f"  {lbl:<28} {q:>4} {n:>6}  {rho:>14.4f}  {p:>12.3e}  {flag}")

print()
for name, r in results_all.items():
    flag = r['verdict_flag']
    print(f"  {r['label']:<28} {r['N']:>4} {r['n_inner']:>6}  {r['rho_right_gue']:>14.4f}  {r['p_right_gue']:>12.3e}  {flag}")

print()
# 판정
n_pos = sum(1 for r in results_all.values() if abs(r['rho_right_gue']) > 0.3 and r['p_right_gue'] < 0.01)
n_curves = len(results_all)
sign_consistent = all(r['rho_right_gue'] < 0 for r in results_all.values())
rho_range = [r['rho_right_gue'] for r in results_all.values() if not math.isnan(r['rho_right_gue'])]
rho_in_range = all(-0.7 <= r <= -0.4 for r in rho_range) if rho_range else False

if n_pos == n_curves and sign_consistent and rho_in_range:
    verdict = "★★★★ degree-독립 보편성 확인 — ρ ∈ [-0.7,-0.4], 음수 일치"
elif n_pos == n_curves and sign_consistent:
    verdict = "★★★ GL(2) 보편성 확인 — |ρ|>0.3 (p<0.01), 음수 일치"
elif n_pos >= 1:
    verdict = f"★★ 조건부 양성 — {n_pos}/{n_curves} 곡선 |ρ|>0.3"
else:
    verdict = "★ 음성 — GL(2) A-gap 상관 없음 (GL(1) 특이 현상)"

print(f"  최종 판정: {verdict}")
print(f"  부호 일치(음수): {sum(1 for r in results_all.values() if r['rho_right_gue']<0)}/{n_curves}")
print()
print(f"  해석:")
if n_pos >= 1:
    rho_vals = [f"{r['rho_right_gue']:.3f}" for r in results_all.values()]
    print(f"  - GL(2) 타원곡선 L-함수에서 ρ ≈ {', '.join(rho_vals)}")
    print(f"  - GL(1) 범위 ρ ≈ -0.55~-0.63과 비교")
    print(f"  - Laurent 진폭 A(γ)의 영점간격 예측: degree-독립 보편성")

# ===== 결과 저장 =====
with open(RESULT_PATH, 'w') as f:
    f.write("=" * 70 + "\n")
    f.write("[사이클 #263] GL(2) A-gap 상관 교차검증 — 최종 결과\n")
    f.write(f"  T=[{T_MIN},{T_MAX}], CAUCHY_R={CAUCHY_R}, N_PTS={N_PTS}\n")
    f.write("=" * 70 + "\n\n")

    for name, r in results_all.items():
        f.write("=" * 70 + "\n")
        f.write(f"  {r['label']}\n")
        f.write("=" * 70 + "\n")
        f.write(f"  영점 수: {r['n_zeros']}, 내부 유효: {r['n_inner']}\n")
        flag_r = '✅' if abs(r['rho_right_gue'])>0.3 and r['p_right_gue']<0.01 else '❌'
        f.write(f"  ρ(A, gap_right_GUE) = {r['rho_right_gue']:+.4f}  (p={r['p_right_gue']:.3e})  {flag_r}\n")
        f.write(f"  ρ(A, gap_min_GUE)   = {r['rho_min_gue']:+.4f}  (p={r['p_min_gue']:.3e})\n")
        f.write(f"  ρ(A_n, A_{{n+1}})    = {r['rho_adj']:+.4f}  (p={r['p_adj']:.3e})\n\n")

    f.write("=" * 70 + "\n")
    f.write("  비교표: ρ(A, gap_right_GUE) 요약\n")
    f.write("=" * 70 + "\n")
    hdr = f"{'L-함수':<30} {'q':>4} {'n내부':>6}  {'ρ(right_GUE)':>14}  {'p값':>12}  유의\n"
    f.write(hdr)
    f.write("-" * 80 + "\n")
    for lbl, q, n, rho, p in GL1_REF:
        flag = '✅' if abs(rho) > 0.3 and p < 0.01 else '❌'
        f.write(f"  {lbl:<28} {q:>4} {n:>6}  {rho:>14.4f}  {p:>12.3e}  {flag}\n")
    f.write("\n")
    for name, r in results_all.items():
        flag = r['verdict_flag']
        f.write(f"  {r['label']:<28} {r['N']:>4} {r['n_inner']:>6}  {r['rho_right_gue']:>14.4f}  {r['p_right_gue']:>12.3e}  {flag}\n")

    f.write("\n")
    f.write("=" * 70 + "\n")
    f.write(f"  최종 판정: {verdict}\n")
    f.write(f"  부호 일치(음수): {sum(1 for r in results_all.values() if r['rho_right_gue']<0)}/{n_curves}\n\n")
    f.write("  [기술적 비고]\n")
    f.write("  - cypari2 (PARI/GP) 인터페이스, lfunlambda + lfunzeros\n")
    f.write("  - Cauchy 적분: Λ'/Λ(1+it), r=0.01, N=64, h=1e-5\n")
    f.write("  - GUE 정규화: gap × (1/π)log(√N·T/(2π))\n")
    f.write("  - 내부 영점: 양끝 2개씩 제외\n")

print(f"\n결과 저장: {RESULT_PATH}")
print("완료.")
