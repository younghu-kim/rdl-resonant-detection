"""
=============================================================================
[Project RDL] 사이클 #262 — χ₅ (mod 5, 복소) A-gap 상관 (v3)
=============================================================================

v1: findroot on complex function → 과결정, 영점 0개
v2: scipy minimize + 검증 → 영점 75개 OK, Cauchy 적분 Gamma 언더플로
v3: bare L'/L + 해석적 ψ 보정 → Gamma 인자 회피

핵심 공식:
  Λ'/Λ(s) = L'/L(s) + (1/2)·ψ((s+a)/2) + (1/2)·log(q/π)

L(s,χ)에는 Gamma 인자가 없으므로 |L| ~ O(1), 언더플로 없음.
ψ, log는 해석적 함수이므로 수치적으로 안정.

결과 파일: results/dirichlet_A_gap_chi5_c262.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats, optimize

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

from bundle_utils import completed_L, CHARACTERS

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

mpmath.mp.dps = 50

T_MIN = 10.0
T_MAX = 200.0
N_SCAN = 6000

CONTOUR_RADIUS = mpmath.mpf('0.01')
N_CONTOUR_PTS  = 48    # c₀, c₁에 충분 (32-64)
DIFF_H         = mpmath.mpf(1) / mpmath.mpf(10**15)

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/dirichlet_A_gap_chi5_c262.txt'
)

REF = {
    'zeta':  {'q': 1, 'rho': -0.5898, 'p': 6.1e-20, 'n': 198, 'rho_adj': 0.4009},
    'chi3':  {'q': 3, 'rho': -0.5733, 'p': 1.4e-7,  'n': 72,  'rho_adj': 0.4183},
    'chi4':  {'q': 4, 'rho': -0.5530, 'p': 1.0e-7,  'n': 80,  'rho_adj': float('nan')},
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 영점 탐색: |Λ| 극소 + scipy (v2와 동일)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def find_zeros_complex_char(char_info, t_min, t_max, n_scan):
    ts = np.linspace(t_min, t_max, n_scan)
    abs_vals = np.zeros(n_scan)

    print(f"  [Step A] |Λ| 스캔 ({n_scan}점)...", flush=True)
    t0_scan = time.time()
    for i, t in enumerate(ts):
        s = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t)))
        val = completed_L(s, char_info)
        abs_vals[i] = float(abs(val))
        if (i + 1) % 1000 == 0:
            print(f"    {i+1}/{n_scan} ({time.time()-t0_scan:.0f}s)", flush=True)
    print(f"  [Step A] 완료 ({time.time()-t0_scan:.0f}s)", flush=True)

    # 극소 후보
    median_abs = np.median(abs_vals)
    candidates = []
    for i in range(1, len(abs_vals) - 1):
        if abs_vals[i] < abs_vals[i-1] and abs_vals[i] < abs_vals[i+1]:
            if abs_vals[i] < median_abs * 0.5:
                candidates.append((ts[i], abs_vals[i]))
    print(f"  [Step B] 극소 후보 {len(candidates)}개", flush=True)

    # scipy 정밀화 + 검증
    dt = ts[1] - ts[0]
    def f_abs_float(t_val):
        s = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_val)))
        return float(abs(completed_L(s, char_info)))

    zeros = []
    for t_cand, _ in candidates:
        try:
            res = optimize.minimize_scalar(
                f_abs_float,
                bounds=(t_cand - dt, t_cand + dt),
                method='bounded',
                options={'xatol': 1e-12}
            )
            t_ref = res.x
            s_check = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_ref)))
            val_check = completed_L(s_check, char_info)
            if float(abs(val_check)) < 1e-20:
                if not zeros or abs(t_ref - zeros[-1]) > 0.05:
                    zeros.append(t_ref)
        except Exception:
            pass

    zeros.sort()
    print(f"  [Step C] 확인된 영점: {len(zeros)}개", flush=True)
    return np.array(zeros)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# v3: bare L'/L + 해석적 ψ 보정으로 Λ'/Λ Laurent 계수 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_A_coeffs_v3(t_zero, char_info):
    """
    Gamma 언더플로 회피 Cauchy 적분.

    Λ'/Λ(s) = L'/L(s) + (1/2)·ψ((s+a)/2) + (1/2)·log(q/π)

    L(s,χ) = Σ χ(n)·n^{-s}는 Gamma 인자 없음 → |L| ~ O(1).
    ψ((s+a)/2)는 mpmath.digamma로 안정적 계산.
    """
    rho = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_zero)))
    r = CONTOUR_RADIUS
    N = N_CONTOUR_PTS
    h = DIFF_H

    q = char_info['q']
    a = char_info['a']
    chi_vals = [mpmath.mpc(v) if isinstance(v, complex) else mpmath.mpf(v)
                for v in char_info['chi']]

    log_q_pi = mpmath.log(mpmath.mpf(q) / mpmath.pi) / 2

    c0_sum = mpmath.mpc(0, 0)
    c1_sum = mpmath.mpc(0, 0)

    for k in range(N):
        theta = 2 * mpmath.pi * k / N
        u = r * mpmath.exp(mpmath.mpc(0, theta))
        s = rho + u

        # bare L(s, χ) via mpmath.dirichlet
        L_val = mpmath.dirichlet(s, chi_vals)
        L_p   = mpmath.dirichlet(s + h, chi_vals)
        L_m   = mpmath.dirichlet(s - h, chi_vals)

        if abs(L_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
            return None  # 영점 위 (contour에서는 발생 안 해야 함)

        # L'/L (bare, no Gamma)
        LL = (L_p - L_m) / (2 * h * L_val)

        # 해석적 Gamma 보정: (1/2)ψ((s+a)/2) + (1/2)log(q/π)
        psi_term = mpmath.digamma((s + a) / 2) / 2

        # 전체 Λ'/Λ
        LambdaL = LL + psi_term + log_q_pi

        # g = Λ'/Λ - 1/u  (정칙 부분)
        g = LambdaL - 1 / u

        c0_sum += g
        c1_sum += g / u

    c0 = c0_sum / N
    c1 = c1_sum / N

    c0_f = complex(c0)
    c1_f = complex(c1)

    A = c0_f.imag**2 + 2 * c1_f.real
    return A, c0_f, c1_f


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 유틸리티
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def gue_normalize(gaps, t_arr, q):
    log_qt = np.log(q * t_arr / (2 * np.pi))
    log_qt = np.maximum(log_qt, 0.01)
    return gaps * log_qt / (2 * np.pi)


def spearman_summary(name, x, y):
    if len(x) < 5:
        return float('nan'), 1.0, f"  ρ({name}) = N/A (n<5)"
    rho, pval = stats.spearmanr(x, y)
    sig = "✅ 유의(p<0.01)" if pval < 0.01 else ("⚠️ p<0.05" if pval < 0.05 else "❌ 비유의")
    return rho, pval, f"  ρ({name}) = {rho:+.4f}  (p={pval:.3e})  {sig}"


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t0 = time.time()
    char_info = CHARACTERS['chi_mod_5']
    q = char_info['q']
    lines = []

    print("=" * 70)
    print(f"[사이클 #262] χ₅ (mod 5, 복소) A-gap 상관 (v3)")
    print(f"  dps={mpmath.mp.dps}, T=[{T_MIN},{T_MAX}], N_scan={N_SCAN}")
    print(f"  Cauchy 적분: bare L'/L + 해석적 ψ 보정 (Γ 언더플로 회피)")
    print("=" * 70)

    lines.append("=" * 70)
    lines.append(f"[사이클 #262] χ₅ (mod 5, 복소) A-gap 상관 (v3)")
    lines.append(f"  dps={mpmath.mp.dps}, T=[{T_MIN},{T_MAX}], N_scan={N_SCAN}")
    lines.append(f"  N_contour={N_CONTOUR_PTS}, r={float(CONTOUR_RADIUS)}")
    lines.append(f"  Cauchy: bare L'/L + ψ 보정 (Γ 언더플로 회피)")
    lines.append("=" * 70)

    # ── Step 0: 영점 탐색 ──
    print(f"\n[Step 0] 복소 지표 영점 탐색 중...", flush=True)
    zeros_t = find_zeros_complex_char(char_info, T_MIN, T_MAX, N_SCAN)
    zeros_t = np.sort(zeros_t)
    n_zeros = len(zeros_t)

    if n_zeros < 10:
        msg = f"❌ 영점 {n_zeros}개 — 분석 불가"
        print(msg)
        lines.append(msg)
        with open(RESULT_PATH, 'w') as f:
            f.write('\n'.join(lines))
        return

    print(f"  영점 {n_zeros}개  t∈[{zeros_t[0]:.3f}, {zeros_t[-1]:.3f}]")
    gap_stats = np.diff(zeros_t)
    print(f"  최소 간격: {np.min(gap_stats):.4f}, 평균 간격: {np.mean(gap_stats):.4f}")
    lines.append(f"\n  영점 {n_zeros}개  t∈[{zeros_t[0]:.3f}, {zeros_t[-1]:.3f}]")
    lines.append(f"  최소 간격: {np.min(gap_stats):.4f}, 평균 간격: {np.mean(gap_stats):.4f}")

    # ── Step 1: A(γ) 계산 (v3: bare L'/L) ──
    print(f"\n[Step 1] A(γ) 계산 중 ({n_zeros}개, v3 bare L'/L)...", flush=True)
    A_list, c0_list, c1_list = [], [], []
    err_count = 0

    for j, t_n in enumerate(zeros_t):
        if (j + 1) % 10 == 0 or j == 0:
            elapsed = time.time() - t0
            print(f"  [{elapsed:.0f}s] {j+1}/{n_zeros}  t={t_n:.3f}", flush=True)
        try:
            result = compute_A_coeffs_v3(t_n, char_info)
            if result is None:
                print(f"  ⚠️ #{j+1} (t={t_n:.3f}) — 적분 실패")
                A_list.append(np.nan)
                c0_list.append(complex(0))
                c1_list.append(complex(0))
                err_count += 1
            else:
                A, c0, c1 = result
                if np.isnan(A) or np.isinf(A):
                    A_list.append(np.nan)
                    err_count += 1
                else:
                    A_list.append(A)
                c0_list.append(c0)
                c1_list.append(c1)
        except Exception as e:
            print(f"  WARNING: #{j+1} t={t_n:.3f}: {e}")
            A_list.append(np.nan)
            c0_list.append(complex(0))
            c1_list.append(complex(0))
            err_count += 1

    A_arr = np.array(A_list)
    valid_mask = ~np.isnan(A_arr)
    n_valid = int(np.sum(valid_mask))

    print(f"  완료: {n_valid}/{n_zeros} 유효  ({err_count} 오류)")
    lines.append(f"\n  A 계산: {n_valid}/{n_zeros} 유효  ({err_count} 오류)")

    if n_valid < 10:
        msg = f"❌ 유효 데이터 부족 ({n_valid}개)"
        print(msg)
        lines.append(msg)
        with open(RESULT_PATH, 'w') as f:
            f.write('\n'.join(lines))
        return

    # Re(c₀) 확인
    re_c0_arr = np.array([c.real for c in c0_list])[valid_mask]
    im_c0_arr = np.array([c.imag for c in c0_list])[valid_mask]
    re_c0_mean = np.mean(np.abs(re_c0_arr))
    im_c0_mean = np.mean(np.abs(im_c0_arr))
    print(f"  |Re(c₀)| = {re_c0_mean:.4f}, |Im(c₀)| = {im_c0_mean:.4f}")
    lines.append(f"  |Re(c₀)| = {re_c0_mean:.4f}, |Im(c₀)| = {im_c0_mean:.4f}")

    # Thm 5 (FE+CC): 복소 비자기쌍대에서도 패리티 유지 (C-241 확인)
    # Re(c₀) should be 0 if CC holds for this character
    if re_c0_mean > 0.01:
        lines.append(f"  ⚠️ Re(c₀) ≠ 0: 비자기쌍대 복소 지표. Thm 5 FE+CC 확인 필요.")
    else:
        lines.append(f"  ✅ Re(c₀) ≈ 0: Thm 5 (FE+CC) 일관.")

    # ── Step 2: gap + GUE 정규화 ──
    gaps = np.diff(zeros_t)
    inner_idx = np.where(valid_mask)[0]
    inner_idx = inner_idx[(inner_idx > 0) & (inner_idx < n_zeros - 1)]

    A_inner = A_arr[inner_idx]
    t_inner = zeros_t[inner_idx]
    gap_right = gaps[inner_idx]
    gap_left  = gaps[inner_idx - 1]
    gap_min_v = np.minimum(gap_right, gap_left)

    gap_right_GUE = gue_normalize(gap_right, t_inner, q)
    gap_left_GUE  = gue_normalize(gap_left,  t_inner, q)
    gap_sum_GUE   = gue_normalize(gap_right + gap_left, t_inner, q)
    gap_min_GUE   = gue_normalize(gap_min_v, t_inner, q)

    n_inner = len(inner_idx)
    print(f"  내부 영점: {n_inner}개", flush=True)
    lines.append(f"  내부 영점: {n_inner}개")

    # ── Step 3: Spearman 상관 ──
    print(f"\n[Step 3] Spearman 상관 분석 (n={n_inner})...", flush=True)
    lines.append(f"\n[Step 3] Spearman 상관 분석 (n={n_inner})")

    for name, y_arr in [
        ('A, gap_right', gap_right),
        ('A, gap_left', gap_left),
        ('A, gap_right_GUE', gap_right_GUE),
        ('A, gap_left_GUE', gap_left_GUE),
        ('A, gap_sum_GUE', gap_sum_GUE),
        ('A, gap_min_GUE', gap_min_GUE),
    ]:
        rho_val, pval, msg = spearman_summary(name, A_inner, y_arr)
        print(msg, flush=True)
        lines.append(msg)

    # 인접 쌍 상관
    inner_consec = inner_idx[:-1]
    inner_consec = inner_consec[np.isin(inner_consec + 1, inner_idx)]
    rho_adj = float('nan')
    if len(inner_consec) >= 5:
        A_n   = A_arr[inner_consec]
        A_np1 = A_arr[inner_consec + 1]
        rho_adj, _, msg = spearman_summary('A_n, A_{n+1}', A_n, A_np1)
        print(msg, flush=True)
        lines.append(msg)

    elapsed = time.time() - t0
    rho_main, p_main, _ = spearman_summary('', A_inner, gap_right_GUE)

    # 판정
    if abs(rho_main) > 0.3 and p_main < 0.01:
        sign_match = "음(ζ와 일치)" if rho_main < 0 else "⚠️ 부호 반전"
        verdict = f"✅ |ρ|>0.3 p<0.01  부호={sign_match}"
    else:
        verdict = "❌ |ρ|≤0.3 또는 p≥0.01"

    print(f"\n  [판정] ρ(A,gap_right_GUE)={rho_main:+.4f} ({verdict})  소요={elapsed:.0f}s")
    lines.append(f"\n  [판정] ρ(A,gap_right_GUE)={rho_main:+.4f} ({verdict})  소요={elapsed:.0f}s")

    # ── 비교표 ──
    lines.append(f"\n{'='*70}")
    lines.append("  비교표: ζ vs χ₃ vs χ₄ vs χ₅")
    lines.append(f"{'='*70}")
    lines.append(f"  {'L-함수':<12} {'q':>3}  {'n':>4}  {'ρ(A,gap_GUE)':>14}  {'p-value':>12}")
    lines.append(f"  {'-'*55}")
    for name, ref in REF.items():
        lines.append(f"  {name:<12} {ref['q']:3d}  {ref['n']:4d}  {ref['rho']:+14.4f}  {ref['p']:12.1e}")
    lines.append(f"  {'chi5':<12} {q:3d}  {n_inner:4d}  {rho_main:+14.4f}  {p_main:12.1e}")

    # 최종 판정
    lines.append(f"\n{'='*70}")
    all_rhos = [REF['zeta']['rho'], REF['chi3']['rho'], REF['chi4']['rho'], rho_main]
    n_sig = sum(1 for r in all_rhos if abs(r) > 0.3)
    if n_sig >= 4:
        final = "★★★ 양성 — 4/4 L-함수 보편적 (복소 지표 포함)"
    elif n_sig >= 3:
        final = "★★★ 양성 — 3/4 L-함수 보편적"
    elif n_sig >= 2:
        final = "★★ 조건부 양성"
    else:
        final = "음성"
    lines.append(f"  최종 판정: {final}")
    lines.append(f"  ρ 평균: {np.mean(all_rhos):+.4f}")
    lines.append(f"  총 소요: {elapsed:.0f}s")
    lines.append(f"{'='*70}")

    print(f"\n  최종 판정: {final}")
    print(f"  ρ 평균: {np.mean(all_rhos):+.4f}")

    with open(RESULT_PATH, 'w') as f:
        f.write('\n'.join(lines))
    print(f"\n결과 저장: {RESULT_PATH}")


if __name__ == '__main__':
    main()
