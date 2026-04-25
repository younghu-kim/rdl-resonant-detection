"""
=============================================================================
[Project RDL] 사이클 #263 — GL(2) A(γ) vs 영점 간격 상관 분석
=============================================================================

목적:
  C-261/262의 GL(1) A-gap 보편성 (ρ ≈ -0.57)을 GL(2)로 확장.
  타원곡선 L-함수에서 A(γ) vs gap_right_GUE Spearman 상관 측정.

  대상:
  - 11a1 (N=11, rank 0, ε=+1): ~94 영점 (T<100)
  - 37a1 (N=37, rank 1, ε=-1): rank 1이므로 L(1)=0 (추가 영점)

  방법:
  - PARI lfunzeros로 영점 위치 수집
  - PARI lfunlambda + Cauchy 적분으로 c₀, c₁ → A(γ)
  - GUE 정규화: gap × d̄(T), d̄ = (1/2π) log(N·T²/(4π²)) (GL(2) 밀도)
  - Spearman ρ(A, gap_right_GUE)

  성공 기준:
  - |ρ| > 0.3, p < 0.01 → GL(2) 보편성 확인

  참고: PARI EC L-함수는 weight 2, 임계선 Re(s)=1 (FE: s ↔ 2-s)

결과 파일: results/gl2_A_gap_c263.txt
=============================================================================
"""

import sys, os, time, math, cmath
import numpy as np
from scipy import stats
import mpmath

import cypari2
pari = cypari2.Pari()
pari.allocatemem(10**9)
pari.default("realprecision", 50)

mpmath.mp.dps = 30

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CONTOUR_RADIUS = 0.01
N_CONTOUR_PTS = 48
H_DERIV = 1e-10
T_MAX = 100
CRITICAL_SIGMA = 1.0  # GL(2) weight 2: 임계선 Re(s)=1

CURVES = [
    {"label": "11a1", "coeffs": [0,-1,1,-10,-20], "N": 11},
    {"label": "37a1", "coeffs": [0,0,1,-1,0],     "N": 37},
]

RESULT_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "results", "gl2_A_gap_c263.txt"
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_A_coeffs_pari(lf_obj, t_zero, N_cond, sigma_c=CRITICAL_SIGMA):
    """
    Bare L'/L + ψ 보정 Cauchy 적분으로 c₀, c₁ 추출.

    Λ'/Λ = L'/L + ψ(s) + (1/2)log(N) - log(2π)

    bare L(s)는 Gamma 인자 미포함 → t 큰 영점에서도 수치 안정.
    Gamma 기여는 해석적 ψ(s) 보정.

    반환: (A, c0, c1) 또는 None (실패)
    """
    rho = complex(sigma_c, t_zero)
    R = CONTOUR_RADIUS
    N = N_CONTOUR_PTS
    h = H_DERIV
    log_corr = 0.5 * math.log(N_cond) - math.log(2 * math.pi)

    c0_sum = 0+0j
    c1_sum = 0+0j

    for k in range(N):
        theta = 2 * math.pi * k / N
        u = R * cmath.exp(1j * theta)
        s = rho + u

        s_str = f'{s.real} + {s.imag}*I'
        sp_str = f'{s.real + h} + {s.imag}*I'
        sm_str = f'{s.real - h} + {s.imag}*I'

        try:
            L_s = complex(pari.lfun(lf_obj, s_str))
            L_p = complex(pari.lfun(lf_obj, sp_str))
            L_m = complex(pari.lfun(lf_obj, sm_str))
        except Exception:
            return None

        if abs(L_s) < 1e-30:
            return None

        # L'/L (bare)
        LL_bare = (L_p - L_m) / (2 * h) / L_s

        # ψ(s) digamma correction
        s_mp = mpmath.mpc(s.real, s.imag)
        psi_val = complex(mpmath.digamma(s_mp))

        # Λ'/Λ = L'/L + ψ(s) + log_correction
        LamLam = LL_bare + psi_val + log_corr

        g = LamLam - 1.0 / u
        c0_sum += g
        c1_sum += g / u

    c0 = c0_sum / N
    c1 = c1_sum / N

    A = c0.imag**2 + 2 * c1.real
    return A, c0, c1


def spearman_summary(name, x, y):
    """Spearman 상관 + 결과 문자열"""
    rho, pval = stats.spearmanr(x, y)
    sig = "✅" if pval < 0.01 else ("⚠️" if pval < 0.05 else "❌")
    return rho, pval, f"  ρ({name}) = {rho:+.4f}  (p={pval:.3e})  {sig}"


def analyze_curve(curve_info, out_lines):
    """단일 타원곡선 A-gap 상관 분석"""
    label = curve_info["label"]
    coeffs = curve_info["coeffs"]
    N_cond = curve_info["N"]

    out_lines.append(f"\n{'='*70}")
    out_lines.append(f"  타원곡선: {label}  (N={N_cond})")
    out_lines.append(f"{'='*70}")

    t0 = time.time()

    # 타원곡선 초기화
    E = pari.ellinit(coeffs)
    lf = pari.lfuncreate(E)
    eps = int(pari.ellrootno(E))
    out_lines.append(f"  ε = {eps}")

    # lfuninit for zero-finding (lfunzeros still needs it internally)
    # bare L(s) evaluation uses lfun(lf, s) directly — no linit needed

    # 영점 수집
    zeros_raw = [float(z) for z in pari.lfunzeros(lf, T_MAX)]
    # t > 1 필터 (rank 1 곡선의 s=1 영점 제외)
    zeros_t = np.array([z for z in zeros_raw if z > 1.0])
    n_zeros = len(zeros_t)
    out_lines.append(f"  영점 {n_zeros}개 (t > 1, T < {T_MAX})")
    if n_zeros < 5:
        out_lines.append("  ❌ 영점 부족 — 건너뜀")
        return None

    out_lines.append(f"  t ∈ [{zeros_t[0]:.2f}, {zeros_t[-1]:.2f}]")
    print(f"\n[{label}] 영점 {n_zeros}개 수집, Cauchy 적분 시작...", flush=True)

    # A(γ) 계산
    A_list = []
    c0_list = []
    c1_list = []
    err_count = 0

    for j in range(n_zeros):
        t_n = zeros_t[j]
        if (j + 1) % 20 == 0 or j == 0:
            elapsed = time.time() - t0
            print(f"  [{label}] 영점 #{j+1}/{n_zeros}  t={t_n:.2f}  경과={elapsed:.0f}s", flush=True)

        result = compute_A_coeffs_pari(lf, t_n, N_cond)
        if result is None:
            A_list.append(np.nan)
            c0_list.append(0+0j)
            c1_list.append(0+0j)
            err_count += 1
            continue

        A, c0, c1 = result
        if np.isnan(A) or np.isinf(A) or A < -100 or A > 1000:
            A_list.append(np.nan)
            err_count += 1
        else:
            A_list.append(A)
        c0_list.append(c0)
        c1_list.append(c1)

    A_arr = np.array(A_list)
    valid = ~np.isnan(A_arr)
    n_valid = int(np.sum(valid))
    elapsed = time.time() - t0

    out_lines.append(f"  A(γ) 계산: {n_valid}/{n_zeros} 유효  ({err_count} 오류)  {elapsed:.0f}s")
    print(f"  [{label}] A 계산 완료: {n_valid}/{n_zeros} 유효, {elapsed:.0f}s", flush=True)

    if n_valid < 10:
        out_lines.append("  ❌ 유효 데이터 부족")
        return None

    # Thm 5 검증
    re_c0 = np.array([c.real for c in c0_list])
    im_c0 = np.array([c.imag for c in c0_list])
    re_c1 = np.array([c.real for c in c1_list])
    im_c1 = np.array([c.imag for c in c1_list])

    mask_c0 = valid & (np.abs(im_c0) > 1e-15)
    mask_c1 = valid & (np.abs(re_c1) > 1e-15)

    if np.sum(mask_c0) > 0:
        thm5_c0 = float(np.max(np.abs(re_c0[mask_c0]) / np.abs(im_c0[mask_c0])))
        out_lines.append(f"  [Thm 5] max|Re(c₀)|/|Im(c₀)| = {thm5_c0:.2e}")

    if np.sum(mask_c1) > 0:
        thm5_c1 = float(np.max(np.abs(im_c1[mask_c1]) / np.abs(re_c1[mask_c1])))
        out_lines.append(f"  [Thm 5] max|Im(c₁)|/|Re(c₁)| = {thm5_c1:.2e}")

    # 간격 계산 (내부 영점)
    gaps = np.diff(zeros_t)
    inner_idx = np.arange(1, n_zeros - 1)
    A_inner = A_arr[inner_idx]
    t_inner = zeros_t[inner_idx]

    gap_left  = gaps[inner_idx - 1]
    gap_right = gaps[inner_idx]
    gap_sum   = gap_left + gap_right
    gap_min   = np.minimum(gap_left, gap_right)

    # GL(2) GUE 정규화
    # 밀도: d̄(t) = (1/2π) log(N·t²/(4π²))  (analytic conductor)
    # 또는 간단히: d̄(t) = (1/2π) [log(N) + 2·log(t) - 2·log(2π)]
    log_analytic_cond = np.log(N_cond) + 2*np.log(t_inner) - 2*np.log(2*np.pi)
    d_bar = log_analytic_cond / (2 * np.pi)  # 평균 밀도
    norm_factor = d_bar  # gap_GUE = gap × d̄(t) → 평균 ~1

    gap_left_n  = gap_left  * norm_factor
    gap_right_n = gap_right * norm_factor
    gap_sum_n   = gap_sum   * norm_factor
    gap_min_n   = gap_min   * norm_factor

    # 유효 내부 필터
    valid_inner = ~np.isnan(A_inner) & (d_bar > 0)
    n_inner = int(np.sum(valid_inner))
    out_lines.append(f"  유효 내부 영점: {n_inner}")

    if n_inner < 10:
        out_lines.append("  ❌ 유효 내부 영점 부족")
        return None

    A_v = A_inner[valid_inner]
    grn_v = gap_right_n[valid_inner]
    gln_v = gap_left_n[valid_inner]
    gsn_v = gap_sum_n[valid_inner]
    gmn_v = gap_min_n[valid_inner]
    gr_v = gap_right[valid_inner]

    # GUE 정규화 간격 통계
    mean_grn = np.mean(grn_v)
    std_grn = np.std(grn_v)
    out_lines.append(f"  GUE gap_right: mean={mean_grn:.3f}, std={std_grn:.3f}")

    # Spearman 상관 분석
    out_lines.append(f"\n  [상관 분석] n={n_inner}")

    results = {}
    for name, arr in [
        ("A, gap_right",     gap_right[valid_inner]),
        ("A, gap_right_GUE", grn_v),
        ("A, gap_left_GUE",  gln_v),
        ("A, gap_sum_GUE",   gsn_v),
        ("A, gap_min_GUE",   gmn_v),
    ]:
        rho, pval, line = spearman_summary(name, A_v, arr)
        results[name] = (rho, pval)
        out_lines.append(line)

    # 인접 쌍 상관
    valid_full = ~np.isnan(A_arr)
    valid_idx = np.where(valid_full)[0]
    pair_mask = np.diff(valid_idx) == 1
    An  = A_arr[valid_idx[:-1][pair_mask]]
    An1 = A_arr[valid_idx[1:][pair_mask]]

    if len(An) > 5:
        rho_pair, pval_pair = stats.spearmanr(An, An1)
        sig_p = "✅" if pval_pair < 0.01 else ("⚠️" if pval_pair < 0.05 else "❌")
        out_lines.append(f"\n  ρ(A_n, A_{{n+1}}) = {rho_pair:+.4f}  (p={pval_pair:.3e}, n={len(An)})  {sig_p}")
    else:
        rho_pair = None

    # A 통계
    out_lines.append(f"\n  A 통계: mean={np.mean(A_v):.3f}, std={np.std(A_v):.3f}, "
                     f"min={np.min(A_v):.3f}, max={np.max(A_v):.3f}")

    rho_main = results.get("A, gap_right_GUE", (None, None))[0]
    pval_main = results.get("A, gap_right_GUE", (None, None))[1]

    return {
        "label": label, "N": N_cond, "eps": eps,
        "n_zeros": n_zeros, "n_valid": n_valid, "n_inner": n_inner,
        "rho": rho_main, "pval": pval_main,
        "rho_pair": rho_pair,
        "A_mean": float(np.mean(A_v)), "A_std": float(np.std(A_v)),
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    print("=" * 70)
    print("[사이클 #263] GL(2) A(γ) vs 영점 간격 상관 분석")
    print(f"  PARI realprecision=50, 윤곽 반경={CONTOUR_RADIUS}, N_pts={N_CONTOUR_PTS}")
    print(f"  대상: {', '.join(c['label'] for c in CURVES)}")
    print("=" * 70)

    out = []
    out.append("=" * 70)
    out.append("[사이클 #263] GL(2) A(γ) vs 영점 간격 상관 분석")
    out.append("=" * 70)
    out.append(f"날짜: 2026-04-25")
    out.append(f"PARI realprecision=50, R={CONTOUR_RADIUS}, N_pts={N_CONTOUR_PTS}, T_MAX={T_MAX}")

    all_results = []

    for curve in CURVES:
        try:
            result = analyze_curve(curve, out)
            if result:
                all_results.append(result)
        except Exception as e:
            out.append(f"\n  ❌ {curve['label']} 오류: {e}")
            import traceback
            traceback.print_exc()

    # 종합 요약
    out.append(f"\n{'='*70}")
    out.append("  종합 요약 — GL(2) A-gap 상관")
    out.append(f"{'='*70}")

    if all_results:
        out.append("")
        out.append("| L-함수 | N | ε | n_inner | ρ(A, gap_right_GUE) | p-value | 판정 |")
        out.append("|--------|---|---|---------|---------------------|---------|------|")

        # GL(1) 참조값 추가
        gl1_ref = [
            ("ζ(s)", 1, "+1", 198, -0.5898, 6.1e-20),
            ("χ₃ mod 3", 3, "+1", 72, -0.5733, 1.4e-7),
            ("χ₄ mod 4", 4, "+1", 80, -0.5530, 1.0e-7),
            ("χ₅ mod 5", 5, "복소", 46, -0.6334, 2.3e-6),
        ]
        for name, q, e, n, rho, pv in gl1_ref:
            sig = "✅" if pv < 0.01 else "❌"
            out.append(f"| {name} (GL1) | {q} | {e} | {n} | {rho:+.4f} | {pv:.1e} | {sig} |")

        for r in all_results:
            sig = "✅" if r["pval"] and r["pval"] < 0.01 else ("⚠️" if r["pval"] and r["pval"] < 0.05 else "❌")
            rho_str = f'{r["rho"]:+.4f}' if r["rho"] is not None else "N/A"
            pval_str = f'{r["pval"]:.1e}' if r["pval"] is not None else "N/A"
            out.append(f"| {r['label']} (GL2) | {r['N']} | {r['eps']:+d} | {r['n_inner']} | {rho_str} | {pval_str} | {sig} |")

        # 판정
        gl2_rhos = [r["rho"] for r in all_results if r["rho"] is not None]
        gl2_pvals = [r["pval"] for r in all_results if r["pval"] is not None]

        if gl2_rhos:
            all_sig = all(p < 0.01 for p in gl2_pvals)
            all_neg = all(r < 0 for r in gl2_rhos)
            mean_rho = np.mean(gl2_rhos)

            out.append(f"\n  GL(2) 평균 ρ = {mean_rho:+.4f}")
            out.append(f"  부호 일치: {sum(1 for r in gl2_rhos if r < 0)}/{len(gl2_rhos)} 음수")
            out.append(f"  유의성: {sum(1 for p in gl2_pvals if p < 0.01)}/{len(gl2_pvals)} (p<0.01)")

            if all_sig and all_neg and abs(mean_rho + 0.57) < 0.15:
                verdict = "★★★★ degree-독립 보편성 확립"
            elif all_sig and all_neg:
                verdict = "★★★ GL(2) 보편성 확인"
            elif any(p < 0.01 for p in gl2_pvals):
                verdict = "★★ 부분 양성"
            else:
                verdict = "★ 음성 또는 약한 신호"

            out.append(f"\n  판정: {verdict}")
    else:
        out.append("  ❌ 모든 곡선 실패")

    elapsed_total = time.time() - t_start
    out.append(f"\n  총 소요: {elapsed_total:.0f}s")

    # 결과 저장
    result_text = "\n".join(out)
    with open(RESULT_PATH, 'w') as f:
        f.write(result_text)

    print(f"\n{'='*70}")
    print(result_text)
    print(f"\n결과 저장: {RESULT_PATH}")
    print(f"총 소요: {elapsed_total:.0f}s")


if __name__ == "__main__":
    main()
