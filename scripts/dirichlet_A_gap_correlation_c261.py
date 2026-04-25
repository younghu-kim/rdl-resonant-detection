"""
=============================================================================
[Project RDL] 사이클 #261 — Dirichlet A-gap 상관 교차검증
=============================================================================

목적:
  C-256에서 확립한 ζ(s) 결과 ρ(A, gap_right_GUE) = -0.59가
  Dirichlet L-함수(χ mod 3, 4, 5)에서도 재현되는지 검증.

  핵심 질문: A(γ)-gap 상관이 ζ에만 고유한가, L-함수 보편 현상인가?

방법:
  1. 각 χ에 대해 t∈[10,200]에서 영점 추출
  2. Cauchy 적분으로 Laurent 계수 c₀, c₁ → A(γ) = Im(c₀)² + 2Re(c₁)
  3. GUE 정규화 간격: gap_GUE = gap × log(qT/(2π)) / (2π)  (q=conductor)
  4. Spearman ρ(A, gap_right_GUE) + 인접 쌍 ρ(A_n, A_{n+1})
  5. ζ 대비 비교표

성공 기준:
  - 3개 χ 모두 |ρ| > 0.3 (p<0.01) → ★★★ 보편적
  - 2개 이상 → ★★ 조건부
  - 1개 이하 → 음성 (ζ-고유 현상)

결과 파일: results/dirichlet_A_gap_c261.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

from bundle_utils import (
    completed_L, find_zeros_dirichlet, CHARACTERS,
)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

mpmath.mp.dps = 60   # 수학자 지시: dps=60 충분

T_MIN = 10.0
T_MAX = 200.0
N_SCAN = 12000       # 영점 탐색 격자 수 (t∈[10,200], 간격~0.016)

CONTOUR_RADIUS = mpmath.mpf('0.01')
N_CONTOUR_PTS  = 128
DIFF_H         = mpmath.mpf(1) / mpmath.mpf(10**18)

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/dirichlet_A_gap_c261.txt'
)

# ζ C-256 참조값 (비교표용)
ZETA_REF = {
    'q': 1,
    'n_zeros': 198,   # 내부 영점
    'rho_gap_right_GUE': -0.5898,
    'p_gap_right_GUE': 6.129e-20,
    'rho_adj': 0.4009,
    'p_adj': 4.421e-09,
    'rho_gap_min_GUE': -0.6209,
}

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수: Cauchy 적분으로 Laurent 계수 계산 (Λ 버전)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_A_coeffs_dirichlet(t_zero, char_info):
    """
    Cauchy 적분으로 디리클레 L-함수의 Laurent 계수 추출.

    Λ'/Λ(ρ+u) = 1/u + c₀ + c₁u + ...  (ρ = 1/2 + iγ 영점)
    g(u) = Λ'/Λ(ρ+u) - 1/u  (정칙 부분)
    cₙ = (1/N) Σ g(u_k) · u_k^{-n}

    A(γ) = Im(c₀)² + 2Re(c₁)  (Cor 4.2와 동일 형식)

    반환: (A, c0, c1) 또는 예외 시 None
    """
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

        # Λ'/Λ (중앙차분)
        Lambda_val = completed_L(s, char_info)
        Lambda_p   = completed_L(s + h, char_info)
        Lambda_m   = completed_L(s - h, char_info)

        if abs(Lambda_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
            # 영점 위: 이 점은 건너뜀 (방어 코드)
            return None

        L_val = (Lambda_p - Lambda_m) / (2 * h * Lambda_val)

        # g = L - 1/u
        g = L_val - 1 / u

        c0_sum += g           # g · u^0
        c1_sum += g / u       # g · u^{-1}

    c0 = c0_sum / N
    c1 = c1_sum / N

    c0_f = complex(c0)
    c1_f = complex(c1)

    A = c0_f.imag**2 + 2 * c1_f.real
    return A, c0_f, c1_f


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 유틸리티
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def spearman_summary(name, x, y):
    """Spearman 상관 계산 + 결과 문자열"""
    if len(x) < 5:
        return float('nan'), 1.0, f"  ρ({name}) = N/A (n<5)"
    rho, pval = stats.spearmanr(x, y)
    sig = "✅ 유의(p<0.01)" if pval < 0.01 else ("⚠️ p<0.05" if pval < 0.05 else "❌ 비유의")
    return rho, pval, f"  ρ({name}) = {rho:+.4f}  (p={pval:.3e})  {sig}"


def gue_normalize(gaps, t_arr, q):
    """
    GUE 정규화: gap_GUE = gap × log(qT/(2π)) / (2π)
    N(T, χ) ≈ (T/2π)log(qT/2πe), dN/dT ≈ log(qT/(2π)) / (2π)
    mean_spacing ≈ 2π / log(qT/(2π))
    """
    log_qt = np.log(q * t_arr / (2 * np.pi))
    log_qt = np.maximum(log_qt, 0.01)   # log<0 방어 (극소 t 시)
    norm_factor = log_qt / (2 * np.pi)
    return gaps * norm_factor


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 단일 χ 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def analyze_character(char_key, char_info, lines_out):
    """
    단일 디리클레 지표에 대한 전체 분석 수행.
    결과를 lines_out 리스트에 추가.
    반환: dict (주요 수치) 또는 None (실패)
    """
    q = char_info['q']
    label = char_info['label']
    t0 = time.time()

    print(f"\n{'='*70}")
    print(f"  {label}  (q={q})")
    print(f"{'='*70}")

    lines_out.append(f"\n{'='*70}")
    lines_out.append(f"  {label}  (q={q})")
    lines_out.append(f"{'='*70}")

    # ────────────────────────────────────────────────
    # Step 0: 영점 탐색
    # ────────────────────────────────────────────────
    print(f"[Step 0] 영점 탐색 중 t∈[{T_MIN},{T_MAX}]  (n_scan={N_SCAN})...", flush=True)
    zeros_t = find_zeros_dirichlet(
        char_info, t_min=T_MIN, t_max=T_MAX, n_scan=N_SCAN
    )
    zeros_t = np.sort(zeros_t)
    n_zeros = len(zeros_t)

    if n_zeros == 0:
        msg = f"⚠️ 영점 0개 — 탐색 실패. 분석 건너뜀."
        print(msg)
        lines_out.append(msg)
        return None

    print(f"  영점 {n_zeros}개 발견  t∈[{zeros_t[0]:.3f}, {zeros_t[-1]:.3f}]")
    print(f"  최소 간격: {np.min(np.diff(zeros_t)):.4f}, 평균 간격: {np.mean(np.diff(zeros_t)):.4f}")
    lines_out.append(f"  영점 {n_zeros}개  t∈[{zeros_t[0]:.3f}, {zeros_t[-1]:.3f}]")
    lines_out.append(f"  최소 간격: {np.min(np.diff(zeros_t)):.4f}, 평균 간격: {np.mean(np.diff(zeros_t)):.4f}")

    if n_zeros < 10:
        msg = f"⚠️ 영점 {n_zeros}개 — 분석에 부족. 건너뜀."
        print(msg)
        lines_out.append(msg)
        return None

    # ────────────────────────────────────────────────
    # Step 1: A(γ) 계산 (Cauchy 적분)
    # ────────────────────────────────────────────────
    print(f"[Step 1] A(γ) 계산 중 ({n_zeros}개 영점)...", flush=True)
    A_list  = []
    c0_list = []
    c1_list = []
    err_count = 0

    for j, t_n in enumerate(zeros_t):
        if (j + 1) % 10 == 0 or j == 0:
            elapsed = time.time() - t0
            print(f"  [{elapsed:.0f}s] {j+1}/{n_zeros}  t={t_n:.3f}", flush=True)

        try:
            result = compute_A_coeffs_dirichlet(t_n, char_info)
            if result is None:
                print(f"  ⚠️ 영점 #{j+1} (t={t_n:.3f}) — 적분 실패 (None)")
                A_list.append(np.nan)
                c0_list.append(complex(0))
                c1_list.append(complex(0))
                err_count += 1
            else:
                A, c0, c1 = result
                if np.isnan(A) or np.isinf(A):
                    print(f"  ⚠️ 영점 #{j+1} A=NaN/Inf")
                    A_list.append(np.nan)
                    err_count += 1
                else:
                    A_list.append(A)
                c0_list.append(c0)
                c1_list.append(c1)
        except Exception as e:
            print(f"  WARNING: 영점 #{j+1} t={t_n:.3f} 오류: {e}")
            A_list.append(np.nan)
            c0_list.append(complex(0))
            c1_list.append(complex(0))
            err_count += 1

    A_arr = np.array(A_list)
    valid_mask = ~np.isnan(A_arr)
    n_valid = int(np.sum(valid_mask))

    print(f"  계산 완료: {n_valid}/{n_zeros} 유효  ({err_count} 오류)")
    lines_out.append(f"\n  A 계산 완료: {n_valid}/{n_zeros} 유효  ({err_count} 오류)")

    if err_count > n_zeros // 2:
        msg = f"❌ 절반 이상 실패({err_count}/{n_zeros}) — 분석 중단"
        print(msg)
        lines_out.append(msg)
        return None

    if n_valid < 10:
        msg = f"⚠️ 유효 A 데이터 {n_valid}개 — 분석 불가"
        print(msg)
        lines_out.append(msg)
        return None

    # ────────────────────────────────────────────────
    # Re(c₀) 확인 (복소 지표 검증)
    # ────────────────────────────────────────────────
    re_c0_arr = np.array([c.real for c in c0_list])[valid_mask]
    im_c0_arr = np.array([c.imag for c in c0_list])[valid_mask]
    re_c1_arr = np.array([c.real for c in c1_list])[valid_mask]

    re_c0_mean = np.mean(np.abs(re_c0_arr))
    im_c0_mean = np.mean(np.abs(im_c0_arr))
    c0_complex_ratio = re_c0_mean / (im_c0_mean + 1e-12)

    print(f"  |Re(c₀)| 평균: {re_c0_mean:.4f},  |Im(c₀)| 평균: {im_c0_mean:.4f}  (비율: {c0_complex_ratio:.3f})")
    lines_out.append(f"  |Re(c₀)| 평균: {re_c0_mean:.4f},  |Im(c₀)| 평균: {im_c0_mean:.4f}")
    if re_c0_mean > 0.1 * im_c0_mean:
        lines_out.append(f"  ℹ️  Re(c₀) ≠ 0 확인 (복소 지표 특성 — A 정의는 여전히 유효)")

    # ────────────────────────────────────────────────
    # Step 2: gap 계산 + GUE 정규화
    # ────────────────────────────────────────────────
    gaps = np.diff(zeros_t)  # 길이 n_zeros-1

    # 내부 영점 인덱스: 첫/마지막 제외
    inner_idx = np.where(valid_mask)[0]
    inner_idx = inner_idx[(inner_idx > 0) & (inner_idx < n_zeros - 1)]

    if len(inner_idx) < 5:
        msg = f"⚠️ 내부 영점 {len(inner_idx)}개 — 상관 분석 불가"
        print(msg)
        lines_out.append(msg)
        return None

    A_inner   = A_arr[inner_idx]
    t_inner   = zeros_t[inner_idx]
    gap_left  = gaps[inner_idx - 1]
    gap_right = gaps[inner_idx]
    gap_sum   = gap_left + gap_right
    gap_min   = np.minimum(gap_left, gap_right)

    # GUE 정규화 (conductor q 반영)
    gap_right_n = gue_normalize(gap_right, t_inner, q)
    gap_left_n  = gue_normalize(gap_left,  t_inner, q)
    gap_sum_n   = gue_normalize(gap_sum,   t_inner, q)
    gap_min_n   = gue_normalize(gap_min,   t_inner, q)

    # 유효 필터
    valid_inner = ~np.isnan(A_inner)
    n_inner_valid = int(np.sum(valid_inner))
    print(f"  내부 영점: {n_inner_valid}/{len(inner_idx)} 유효")
    lines_out.append(f"  내부 영점 (상관 분석용): {n_inner_valid}/{len(inner_idx)}")

    A_v    = A_inner[valid_inner]
    grn_v  = gap_right_n[valid_inner]
    gln_v  = gap_left_n[valid_inner]
    gsn_v  = gap_sum_n[valid_inner]
    gmn_v  = gap_min_n[valid_inner]
    gr_v   = gap_right[valid_inner]
    gl_v   = gap_left[valid_inner]

    # ────────────────────────────────────────────────
    # Step 3: Spearman 상관
    # ────────────────────────────────────────────────
    print(f"\n[Step 3] Spearman 상관 분석 (n={n_inner_valid})...", flush=True)
    lines_out.append(f"\n  [Spearman 상관 분석]  n={n_inner_valid}")

    # 원시 간격
    for name, arr in [("A, gap_right", gr_v), ("A, gap_left", gl_v)]:
        rho, pval, line = spearman_summary(name, A_v, arr)
        print(line)
        lines_out.append(line)

    # GUE 정규화 간격
    lines_out.append("")
    corr_gue = {}
    for name, arr in [
        ("A, gap_right_GUE", grn_v),
        ("A, gap_left_GUE",  gln_v),
        ("A, gap_sum_GUE",   gsn_v),
        ("A, gap_min_GUE",   gmn_v),
    ]:
        rho, pval, line = spearman_summary(name, A_v, arr)
        print(line)
        lines_out.append(line)
        corr_gue[name.split(', ')[1]] = (rho, pval)

    # 인접 쌍 상관
    A_valid_all = A_arr[valid_mask]
    if len(A_valid_all) >= 4:
        A_n   = A_valid_all[:-1]
        A_np1 = A_valid_all[1:]
        rho_adj, pval_adj, line_adj = spearman_summary("A_n, A_{n+1}", A_n, A_np1)
        print(line_adj)
        lines_out.append("")
        lines_out.append(line_adj)
    else:
        rho_adj, pval_adj = float('nan'), 1.0
        lines_out.append("  인접 쌍: 데이터 부족")

    # ────────────────────────────────────────────────
    # Step 4: 판정
    # ────────────────────────────────────────────────
    rho_right, pval_right = corr_gue.get('gap_right_GUE', (float('nan'), 1.0))
    rho_min,   pval_min   = corr_gue.get('gap_min_GUE',   (float('nan'), 1.0))

    is_significant = (abs(rho_right) > 0.3) and (pval_right < 0.01)
    direction_match = (rho_right < 0)  # ζ와 같은 부호?

    status = "✅ |ρ|>0.3 p<0.01" if is_significant else "❌ 비유의"
    dir_str = "음(ζ와 일치)" if direction_match else "양(ζ와 반대!)"

    elapsed = time.time() - t0
    summary_line = (
        f"  [판정] ρ(A,gap_right_GUE)={rho_right:+.4f} ({status})  "
        f"부호={dir_str}  소요={elapsed:.0f}s"
    )
    print(summary_line)
    lines_out.append("")
    lines_out.append(summary_line)

    return {
        'label': label,
        'q': q,
        'n_zeros': n_zeros,
        'n_inner': n_inner_valid,
        'rho_gap_right_GUE': rho_right,
        'p_gap_right_GUE': pval_right,
        'rho_gap_min_GUE': rho_min,
        'p_gap_min_GUE': pval_min,
        'rho_adj': rho_adj,
        'p_adj': pval_adj,
        're_c0_mean': re_c0_mean,
        'elapsed': elapsed,
        'is_significant': is_significant,
        'direction_neg': direction_match,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_global = time.time()

    print("=" * 70)
    print("[사이클 #261] Dirichlet A-gap 상관 교차검증")
    print(f"  dps={mpmath.mp.dps}, T=[{T_MIN},{T_MAX}], N_scan={N_SCAN}")
    print(f"  반경={CONTOUR_RADIUS}, N_pts={N_CONTOUR_PTS}")
    print("=" * 70)

    lines_out = []
    lines_out.append("=" * 70)
    lines_out.append("[사이클 #261] Dirichlet A-gap 상관 교차검증")
    lines_out.append(f"  dps={mpmath.mp.dps}, T=[{T_MIN},{T_MAX}], N_scan={N_SCAN}")
    lines_out.append(f"  반경={CONTOUR_RADIUS}, N_pts={N_CONTOUR_PTS}")
    lines_out.append("=" * 70)

    # 3개 지표 분석
    char_order = ['chi_mod_3', 'chi_mod_4', 'chi_mod_5']
    results = {}

    for char_key in char_order:
        char_info = CHARACTERS[char_key]
        res = analyze_character(char_key, char_info, lines_out)
        if res is not None:
            results[char_key] = res
        else:
            results[char_key] = None

    # ────────────────────────────────────────────────
    # 비교표 출력
    # ────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("  비교표: A(γ)-gap GUE 상관 요약")
    print("=" * 70)

    lines_out.append("\n" + "=" * 70)
    lines_out.append("  비교표: A(γ)-gap GUE 상관 요약")
    lines_out.append("=" * 70)

    header = f"{'L-함수':<22} {'q':>3} {'n영점':>6} {'n내부':>6} {'ρ(gap_right_GUE)':>18} {'p값':>12} {'ρ(A_n,A_{n+1})':>16} {'|ρ|>0.3':>8}"
    print(header)
    lines_out.append(header)
    sep = "-" * 90
    print(sep)
    lines_out.append(sep)

    # ζ 참조행
    z = ZETA_REF
    z_sig = "✅" if abs(z['rho_gap_right_GUE']) > 0.3 and z['p_gap_right_GUE'] < 0.01 else "❌"
    z_line = (f"{'ζ(s) [C-256참조]':<22} {z['q']:>3} {z['n_zeros']:>6} {'198':>6} "
              f"{z['rho_gap_right_GUE']:>+18.4f} {z['p_gap_right_GUE']:>12.3e} "
              f"{z['rho_adj']:>+16.4f} {z_sig:>8}")
    print(z_line)
    lines_out.append(z_line)

    n_significant = 0
    n_neg_direction = 0

    for char_key in char_order:
        res = results[char_key]
        if res is None:
            row = f"{'  '+CHARACTERS[char_key]['label']:<22} {CHARACTERS[char_key]['q']:>3} {'—':>6} {'—':>6} {'실패':>18} {'—':>12} {'—':>16} {'❌':>8}"
        else:
            sig = "✅" if res['is_significant'] else "❌"
            row = (f"{'  '+res['label']:<22} {res['q']:>3} {res['n_zeros']:>6} {res['n_inner']:>6} "
                   f"{res['rho_gap_right_GUE']:>+18.4f} {res['p_gap_right_GUE']:>12.3e} "
                   f"{res['rho_adj']:>+16.4f} {sig:>8}")
            if res['is_significant']:
                n_significant += 1
            if res.get('direction_neg', False):
                n_neg_direction += 1
        print(row)
        lines_out.append(row)

    # 보조: gap_min_GUE 상관
    print()
    lines_out.append("")
    print("  [보조] ρ(A, gap_min_GUE):")
    lines_out.append("  [보조] ρ(A, gap_min_GUE):")
    z_min_line = f"    ζ(s): ρ={z['rho_gap_min_GUE']:+.4f}"
    print(z_min_line)
    lines_out.append(z_min_line)
    for char_key in char_order:
        res = results[char_key]
        if res is not None:
            r_min_line = f"    {res['label']}: ρ={res['rho_gap_min_GUE']:+.4f}  (p={res['p_gap_min_GUE']:.3e})"
            print(r_min_line)
            lines_out.append(r_min_line)

    # ────────────────────────────────────────────────
    # 최종 판정
    # ────────────────────────────────────────────────
    print("\n" + "=" * 70)
    lines_out.append("\n" + "=" * 70)

    if n_significant == 3:
        verdict = "★★★ 보편적 — 3개 χ 모두 |ρ|>0.3 (p<0.01). A-gap 상관은 L-함수 보편 현상."
        verdict_code = 3
    elif n_significant >= 2:
        verdict = "★★ 조건부 양성 — 2개 이상 유의. 보편성 가능성 높음."
        verdict_code = 2
    elif n_significant == 1:
        verdict = "★ 약한 신호 — 1개만 유의. ζ-고유 가능성 있음 (→ B-44 검토)."
        verdict_code = 1
    else:
        verdict = "음성 — 0개 유의. ζ-고유 현상. → B-44 개설 권고."
        verdict_code = 0

    dir_note = f"  부호 일치(ζ처럼 음수): {n_neg_direction}/3"

    print(f"  최종 판정: {verdict}")
    print(dir_note)
    lines_out.append(f"  최종 판정: {verdict}")
    lines_out.append(dir_note)

    if n_significant >= 2:
        lines_out.append("\n  해석:")
        lines_out.append("  - Laurent 진폭 A(γ)는 영점 간격의 보편적 예측자 (Paper 4 중심 정리 후보)")
        lines_out.append("  - ζ와 Dirichlet L-함수에서 동일 부호 음의 상관 → 큰 A = 작은 정규화 간격")
        lines_out.append("  - GUE 통계에 숨겨진 스펙트럼 기하 구조 (c₁ 배경항 기반) 존재 가능")
    else:
        lines_out.append("\n  해석:")
        lines_out.append("  - A-gap 상관이 ζ-특이적이면 conductor 의존성 → B-44 ('conductor 의존성')")
        lines_out.append("  - 추가 검증: q=7,11 등 더 큰 conductor에서 재시도 권고")

    elapsed_total = time.time() - t_global
    lines_out.append(f"\n  총 소요 시간: {elapsed_total:.0f}s")
    print(f"\n  총 소요 시간: {elapsed_total:.0f}s")

    # ────────────────────────────────────────────────
    # 결과 파일 저장
    # ────────────────────────────────────────────────
    os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)
    with open(RESULT_PATH, 'w', encoding='utf-8') as f:
        for line in lines_out:
            f.write(line + '\n')

    print(f"\n결과 저장: {RESULT_PATH}")
    print("완료.")


if __name__ == '__main__':
    main()
