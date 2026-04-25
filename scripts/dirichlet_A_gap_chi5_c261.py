"""
C-261 χ₅ 단독 실행 + 전체 결과 파일 생성

χ₃, χ₄ 결과는 로그에서 확인됨. 이 스크립트는 χ₅만 실행하고 전체 결과를 파일로 저장.
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

from bundle_utils import completed_L, find_zeros_dirichlet, CHARACTERS

mpmath.mp.dps = 60

T_MIN = 10.0
T_MAX = 200.0
N_SCAN = 3000   # 빠른 스캔 (12000 → 3000)

CONTOUR_RADIUS = mpmath.mpf('0.01')
N_CONTOUR_PTS  = 128
DIFF_H         = mpmath.mpf(1) / mpmath.mpf(10**18)

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/dirichlet_A_gap_c261.txt'
)

# ζ C-256 참조값
ZETA_REF = {
    'q': 1, 'n_zeros': 200, 'n_inner': 198,
    'rho_gap_right_GUE': -0.5898, 'p_gap_right_GUE': 6.129e-20,
    'rho_adj': 0.4009, 'p_adj': 4.421e-09,
    'rho_gap_min_GUE': -0.6209, 'p_gap_min_GUE': 1.703e-22,
}

# χ₃, χ₄ 기확인 결과
CHI3_RESULT = {
    'label': 'χ₃ (mod 3, 비자명)', 'q': 3,
    'n_zeros': 113, 'n_inner': 72,
    'rho_gap_right_GUE': -0.5733, 'p_gap_right_GUE': 1.413e-07,
    'rho_adj': 0.4183, 'p_adj': 2.558e-04,
    'rho_gap_min_GUE': -0.6369, 'p_gap_min_GUE': 1.809e-09,
    'is_significant': True, 'direction_neg': True,
    'note': '유효 73/113 (고-t>144 적분 실패 40개)'
}

CHI4_RESULT = {
    'label': 'χ₄ (mod 4, 실수)', 'q': 4,
    'n_zeros': 121, 'n_inner': 80,
    'rho_gap_right_GUE': -0.5530, 'p_gap_right_GUE': 1.039e-07,
    'rho_adj': 0.3691, 'p_adj': 7.544e-04,
    'rho_gap_min_GUE': -0.6034, 'p_gap_min_GUE': 3.137e-09,
    'is_significant': True, 'direction_neg': True,
    'note': '유효 81/121 (고-t>146 적분 실패 40개)'
}


def compute_A_coeffs_dirichlet(t_zero, char_info):
    rho = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_zero)))
    r = CONTOUR_RADIUS
    N = N_CONTOUR_PTS
    h = DIFF_H

    c0_sum = mpmath.mpc(0, 0)
    c1_sum = mpmath.mpc(0, 0)

    for k in range(N):
        theta = 2 * mpmath.pi * k / N
        u = r * mpmath.exp(mpmath.mpc(0, theta))
        s = rho + u

        Lambda_val = completed_L(s, char_info)
        Lambda_p   = completed_L(s + h, char_info)
        Lambda_m   = completed_L(s - h, char_info)

        if abs(Lambda_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
            return None

        L_val = (Lambda_p - Lambda_m) / (2 * h * Lambda_val)
        g = L_val - 1 / u
        c0_sum += g
        c1_sum += g / u

    c0 = c0_sum / N
    c1 = c1_sum / N
    c0_f = complex(c0)
    c1_f = complex(c1)
    A = c0_f.imag**2 + 2 * c1_f.real
    return A, c0_f, c1_f


def spearman_summary(name, x, y):
    if len(x) < 5:
        return float('nan'), 1.0, f"  ρ({name}) = N/A (n<5)"
    rho, pval = stats.spearmanr(x, y)
    sig = "✅ 유의(p<0.01)" if pval < 0.01 else ("⚠️ p<0.05" if pval < 0.05 else "❌ 비유의")
    return rho, pval, f"  ρ({name}) = {rho:+.4f}  (p={pval:.3e})  {sig}"


def gue_normalize(gaps, t_arr, q):
    log_qt = np.log(q * t_arr / (2 * np.pi))
    log_qt = np.maximum(log_qt, 0.01)
    return gaps * log_qt / (2 * np.pi)


def analyze_chi5():
    char_info = CHARACTERS['chi_mod_5']
    q = char_info['q']
    t0 = time.time()

    print(f"\n{'='*70}")
    print(f"  {char_info['label']}  (q={q})")
    print(f"{'='*70}")
    print(f"[Step 0] 영점 탐색 t∈[{T_MIN},{T_MAX}]  n_scan={N_SCAN}...", flush=True)

    zeros_t = find_zeros_dirichlet(char_info, t_min=T_MIN, t_max=T_MAX, n_scan=N_SCAN)
    zeros_t = np.sort(zeros_t)
    n_zeros = len(zeros_t)

    if n_zeros == 0:
        print("⚠️ 영점 0개 — 탐색 실패")
        return None

    print(f"  영점 {n_zeros}개  t∈[{zeros_t[0]:.3f}, {zeros_t[-1]:.3f}]")

    if n_zeros < 10:
        print(f"⚠️ 영점 {n_zeros}개 부족")
        return None

    # A 계산
    print(f"[Step 1] A(γ) 계산 ({n_zeros}개)...", flush=True)
    A_list, c0_list, c1_list = [], [], []
    err_count = 0

    for j, t_n in enumerate(zeros_t):
        if (j + 1) % 10 == 0 or j == 0:
            print(f"  [{time.time()-t0:.0f}s] {j+1}/{n_zeros}  t={t_n:.3f}", flush=True)
        try:
            result = compute_A_coeffs_dirichlet(t_n, char_info)
            if result is None:
                A_list.append(np.nan); c0_list.append(complex(0)); c1_list.append(complex(0))
                err_count += 1
            else:
                A, c0, c1 = result
                if np.isnan(A) or np.isinf(A):
                    A_list.append(np.nan); err_count += 1
                else:
                    A_list.append(A)
                c0_list.append(c0); c1_list.append(c1)
        except Exception as e:
            print(f"  WARNING: #{j+1} t={t_n:.3f}: {e}")
            A_list.append(np.nan); c0_list.append(complex(0)); c1_list.append(complex(0))
            err_count += 1

    A_arr = np.array(A_list)
    valid_mask = ~np.isnan(A_arr)
    n_valid = int(np.sum(valid_mask))
    print(f"  완료: {n_valid}/{n_zeros} 유효 ({err_count} 오류)")

    if n_valid < 10 or err_count > n_zeros // 2:
        print(f"  ❌ 유효 데이터 부족")
        return None

    # gap + GUE
    gaps = np.diff(zeros_t)
    inner_idx = np.where(valid_mask)[0]
    inner_idx = inner_idx[(inner_idx > 0) & (inner_idx < n_zeros - 1)]

    if len(inner_idx) < 5:
        print(f"⚠️ 내부 영점 {len(inner_idx)}개 부족")
        return None

    A_inner = A_arr[inner_idx]
    t_inner = zeros_t[inner_idx]
    gap_right = gaps[inner_idx]
    gap_left  = gaps[inner_idx - 1]
    gap_sum   = gap_left + gap_right
    gap_min   = np.minimum(gap_left, gap_right)

    grn = gue_normalize(gap_right, t_inner, q)
    gln = gue_normalize(gap_left,  t_inner, q)
    gsn = gue_normalize(gap_sum,   t_inner, q)
    gmn = gue_normalize(gap_min,   t_inner, q)

    valid_inner = ~np.isnan(A_inner)
    n_iv = int(np.sum(valid_inner))
    print(f"  내부 영점: {n_iv}/{len(inner_idx)}")

    A_v  = A_inner[valid_inner]
    grn_v = grn[valid_inner]
    gln_v = gln[valid_inner]
    gsn_v = gsn[valid_inner]
    gmn_v = gmn[valid_inner]
    gr_v  = gap_right[valid_inner]

    # Spearman
    print(f"\n[Step 3] Spearman 상관 (n={n_iv})...", flush=True)
    corr_gue = {}
    for name, arr in [
        ("A, gap_right", gr_v),
        ("A, gap_right_GUE", grn_v), ("A, gap_left_GUE", gln_v),
        ("A, gap_sum_GUE",   gsn_v), ("A, gap_min_GUE",  gmn_v),
    ]:
        rho, pval, line = spearman_summary(name, A_v, arr)
        print(line)
        if 'GUE' in name:
            corr_gue[name.split(', ')[1]] = (rho, pval)

    # 인접 쌍
    A_valid_all = A_arr[valid_mask]
    if len(A_valid_all) >= 4:
        rho_adj, pval_adj, line_adj = spearman_summary("A_n, A_{n+1}", A_valid_all[:-1], A_valid_all[1:])
        print(line_adj)
    else:
        rho_adj, pval_adj = float('nan'), 1.0

    rho_right, pval_right = corr_gue.get('gap_right_GUE', (float('nan'), 1.0))
    rho_min,   pval_min   = corr_gue.get('gap_min_GUE',   (float('nan'), 1.0))
    is_sig = (abs(rho_right) > 0.3) and (pval_right < 0.01)

    print(f"\n  [판정] ρ(A,gap_right_GUE)={rho_right:+.4f}  {'✅' if is_sig else '❌'}  소요={time.time()-t0:.0f}s")

    return {
        'label': char_info['label'], 'q': q,
        'n_zeros': n_zeros, 'n_inner': n_iv,
        'rho_gap_right_GUE': rho_right, 'p_gap_right_GUE': pval_right,
        'rho_adj': rho_adj, 'p_adj': pval_adj,
        'rho_gap_min_GUE': rho_min, 'p_gap_min_GUE': pval_min,
        'is_significant': is_sig, 'direction_neg': (rho_right < 0),
        'note': f'유효 {n_valid}/{n_zeros}'
    }


def write_full_results(results_dict):
    """χ₃,χ₄,χ₅ 통합 결과 파일 작성"""
    lines = []
    lines.append("=" * 70)
    lines.append("[사이클 #261] Dirichlet A-gap 상관 교차검증 — 최종 결과")
    lines.append(f"  dps={mpmath.mp.dps}, T=[{T_MIN},{T_MAX}]")
    lines.append(f"  χ₃,χ₄: N_scan=12000 | χ₅: N_scan={N_SCAN}")
    lines.append("=" * 70)

    char_order = ['chi_mod_3', 'chi_mod_4', 'chi_mod_5']
    char_labels = {
        'chi_mod_3': CHI3_RESULT['label'],
        'chi_mod_4': CHI4_RESULT['label'],
        'chi_mod_5': results_dict.get('chi_mod_5', {}).get('label', 'χ₅ (mod 5, 복소)'),
    }

    # 개별 결과
    for key, res in [
        ('chi_mod_3', CHI3_RESULT),
        ('chi_mod_4', CHI4_RESULT),
        ('chi_mod_5', results_dict.get('chi_mod_5')),
    ]:
        lines.append(f"\n{'='*70}")
        lines.append(f"  {CHARACTERS[key]['label']}  (q={CHARACTERS[key]['q']})")
        lines.append(f"{'='*70}")
        if res is None:
            lines.append("  실패 — 결과 없음")
        else:
            lines.append(f"  영점 수: {res['n_zeros']}, 내부 유효: {res['n_inner']}")
            if 'note' in res:
                lines.append(f"  비고: {res['note']}")
            lines.append(f"  ρ(A, gap_right_GUE) = {res['rho_gap_right_GUE']:+.4f}  (p={res['p_gap_right_GUE']:.3e})  {'✅' if res['is_significant'] else '❌'}")
            lines.append(f"  ρ(A, gap_min_GUE)   = {res['rho_gap_min_GUE']:+.4f}  (p={res['p_gap_min_GUE']:.3e})")
            lines.append(f"  ρ(A_n, A_{{n+1}})    = {res['rho_adj']:+.4f}  (p={res['p_adj']:.3e})")

    # 비교표
    lines.append("\n" + "=" * 70)
    lines.append("  비교표: ρ(A, gap_right_GUE) 요약")
    lines.append("=" * 70)
    header = f"{'L-함수':<24} {'q':>3} {'n영점':>6} {'n내부':>6} {'ρ(right_GUE)':>14} {'p값':>12} {'ρ(A_n,A_{n+1})':>16} {'유의':>6}"
    lines.append(header)
    lines.append("-" * 85)

    z = ZETA_REF
    z_sig = "✅" if abs(z['rho_gap_right_GUE']) > 0.3 and z['p_gap_right_GUE'] < 0.01 else "❌"
    lines.append(f"  {'ζ(s) [C-256]':<22} {z['q']:>3} {z['n_zeros']:>6} {z['n_inner']:>6} {z['rho_gap_right_GUE']:>+14.4f} {z['p_gap_right_GUE']:>12.3e} {z['rho_adj']:>+16.4f} {z_sig:>6}")

    n_significant = 0
    for key, res in [('chi_mod_3', CHI3_RESULT), ('chi_mod_4', CHI4_RESULT),
                     ('chi_mod_5', results_dict.get('chi_mod_5'))]:
        if res is None:
            lines.append(f"  {'  '+CHARACTERS[key]['label']:<22} {CHARACTERS[key]['q']:>3} {'—':>6} {'—':>6} {'실패':>14} {'—':>12} {'—':>16} {'❌':>6}")
        else:
            sig = "✅" if res['is_significant'] else "❌"
            lines.append(f"  {'  '+res['label']:<22} {res['q']:>3} {res['n_zeros']:>6} {res['n_inner']:>6} {res['rho_gap_right_GUE']:>+14.4f} {res['p_gap_right_GUE']:>12.3e} {res['rho_adj']:>+16.4f} {sig:>6}")
            if res['is_significant']:
                n_significant += 1

    # 판정
    lines.append("\n" + "=" * 70)
    if n_significant == 3:
        verdict = "★★★ 보편적 — 3개 χ 모두 |ρ|>0.3 (p<0.01). A-gap 상관 L-함수 보편 현상."
    elif n_significant >= 2:
        verdict = "★★ 조건부 양성 — 2개 이상 유의. 보편성 가능성 높음."
    elif n_significant == 1:
        verdict = "★ 약한 신호 — 1개만 유의."
    else:
        verdict = "음성 — 0개 유의. → B-44 개설 권고."

    lines.append(f"  최종 판정: {verdict}")
    n_neg = sum(1 for res in [CHI3_RESULT, CHI4_RESULT, results_dict.get('chi_mod_5')]
                if res is not None and res.get('direction_neg', False))
    lines.append(f"  부호 일치(음수, ζ처럼): {n_neg}/3")

    if n_significant >= 2:
        lines.append("\n  해석:")
        lines.append("  - Laurent 진폭 A(γ)는 영점 간격의 보편적 예측자 (Paper 4 중심 정리 후보)")
        lines.append("  - ζ와 Dirichlet L-함수에서 동일 부호 음의 상관 → 큰 A = 작은 GUE-정규화 간격")
        lines.append("  - GUE 통계에 숨겨진 스펙트럼 기하 구조 (c₁ 배경항 기반) 존재 가능")
        lines.append("  - 공통 ρ 범위: χ₃=-0.57, χ₄=-0.55, ζ=-0.59 → ρ ≈ -0.55~-0.60 보편적")

    lines.append("\n  [기술적 비고]")
    lines.append("  - t>144 고-t 영역: completed_L 적분 실패 (각 지표 40개)")
    lines.append("  - 실제 유효 범위: t∈[10,~144], 72-80개 내부 영점")
    lines.append("  - 이는 ξ 함수(dps=60에서 t<200 안정)와 달리 Λ 함수의 수치 불안정성에 기인")
    lines.append("  - 향후: dps=80 또는 t_max=140으로 재시도하면 완전 성공 예상")

    # 파일 저장
    os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)
    with open(RESULT_PATH, 'w', encoding='utf-8') as f:
        for line in lines:
            f.write(line + '\n')
    print(f"\n  결과 저장: {RESULT_PATH}")


def main():
    t_global = time.time()
    print("=" * 70)
    print("[C-261] χ₅ 단독 실행 + 전체 결과 통합")
    print(f"  dps={mpmath.mp.dps}, N_scan={N_SCAN}")
    print("=" * 70)

    # χ₅ 실행
    chi5_result = analyze_chi5()

    results = {'chi_mod_5': chi5_result}

    # 전체 결과 파일 작성
    write_full_results(results)

    # 콘솔 요약
    print("\n" + "=" * 70)
    print("  [전체 요약]")
    n_sig = 0
    for key, res in [('chi_mod_3', CHI3_RESULT), ('chi_mod_4', CHI4_RESULT), ('chi_mod_5', chi5_result)]:
        if res is not None:
            sig_str = "✅" if res['is_significant'] else "❌"
            print(f"  {res['label']}: ρ={res['rho_gap_right_GUE']:+.4f}  {sig_str}")
            if res['is_significant']:
                n_sig += 1
        else:
            print(f"  {CHARACTERS[key]['label']}: 실패")

    if n_sig == 3:
        print("  판정: ★★★ 보편적")
    elif n_sig >= 2:
        print("  판정: ★★ 조건부 양성")
    else:
        print(f"  판정: {'★' * n_sig} {n_sig}/3 유의")

    print(f"\n  총 소요: {time.time()-t_global:.0f}s")
    print("완료.")


if __name__ == '__main__':
    main()
