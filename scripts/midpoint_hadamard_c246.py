"""
=============================================================================
[C-246] 중간점 곡률의 Hadamard 예측 — 비국소성 메커니즘 규명
=============================================================================
핵심 질문: 중간점 κ_mid의 변동이 Hadamard 곱에서 정량적으로 도출되는가?

배경:
  v2b 실험 (#18)에서 κ_mid vs gap_next: ρ=-0.654 (p=7e-30) 양성.
  NN 상쇄: ρ(κ_mid, gap) = -0.03 확인.
  GUE에서는 비국소 상관 없음 (ρ≈0).
  → "왜?" 에 대한 해석적 답이 없음.

가설 (Hadamard 도출):
  중간점 m = (γ_n + γ_{n+1})/2에서 NN 쌍의 기여가 O(Δ/m)으로 상쇄.
  잔여 곡률은 차차근접 영점 (γ_{n-1}, γ_{n+2})이 지배.
  고t 근사:
    κ_mid ≈ [1/(Δ_{-1}+Δ/2) - 1/(Δ_{+1}+Δ/2)]²
  여기서 Δ = gap, Δ_{-1} = gap_{n-1}, Δ_{+1} = gap_{n+1}.

검증 계획:
  (A) κ_mid_exact  — 해석적 L(s) 공식으로 계산
  (B) κ_NN_removed — Hadamard에서 NN 쌍 제거
  (C) κ_NNN_only   — 차차근접 2영점만의 기여
  (D) κ_asymptotic — 고t 근사 공식 [1/(Δ_{-1}+Δ/2) - 1/(Δ_{+1}+Δ/2)]²
  (E) κ_multi_zero — NN 제외 + 좌우 5영점씩 (10영점 합)

  ρ(A, D) > 0.7이면 ★★★ (공식이 변동의 >50% 설명)
  ρ(A, E) > 0.9이면 ★★★★ (Hadamard 구조가 지배적)

성공 기준:
  1. ρ(κ_exact, κ_asymptotic) > 0.5, p < 0.01 — 공식 양성
  2. ρ(κ_exact, κ_NN_removed) > 0.95 — NN 상쇄 확인
  3. ρ(κ_exact, κ_multi_zero) > ρ(κ_exact, κ_NNN_only) — 다영점 효과
  4. ρ(κ_asymptotic, gap_next) < -0.3 — 공식이 비국소 상관 재현
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats

mpmath.mp.dps = 100

RESULT_FILE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "results", "midpoint_hadamard_c246.txt"
)

N_ZEROS = 300  # γ_1 ~ γ_300 (t ≈ 14 ~ 560)


def log(msg=""):
    print(msg, flush=True)


# ── 해석적 L(s) = ξ'/ξ(s) ──────────────────────────────────────────────────
def connection_analytic(s):
    """
    L(s) = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2 + ζ'/ζ(s)
    ξ 직접 계산 없이 안전.
    """
    zeta_val = mpmath.zeta(s)
    if abs(zeta_val) < mpmath.mpf(10) ** (-mpmath.mp.dps + 15):
        return mpmath.mpc(0, 1e15)  # 영점 근방

    h = mpmath.mpf(10) ** (-20)
    zeta_deriv = (mpmath.zeta(s + h) - mpmath.zeta(s - h)) / (2 * h)
    zeta_ratio = zeta_deriv / zeta_val

    result = 1 / s + 1 / (s - 1) - mpmath.log(mpmath.pi) / 2
    result += mpmath.digamma(s / 2) / 2
    result += zeta_ratio
    return result


def main():
    t_start = time.time()
    log("=" * 70)
    log("[C-246] 중간점 곡률의 Hadamard 예측 검증")
    log(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}")
    log(f"N_zeros: {N_ZEROS}, dps: {mpmath.mp.dps}")
    log("=" * 70)

    # ── 1. 영점 수집 ──
    log("\n1단계: ζ 영점 수집...")
    zeros = []
    for n in range(1, N_ZEROS + 1):
        gamma = mpmath.zetazero(n).imag
        zeros.append(float(gamma))
    zeros = np.array(zeros)
    log(f"  수집: {len(zeros)}개, t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")

    # 간격 계산
    gaps = np.diff(zeros)  # gaps[i] = zeros[i+1] - zeros[i]
    n_pairs = len(gaps)
    log(f"  연속 쌍: {n_pairs}개")

    # ── 2. 중간점 계산 ──
    # 유효 범위: index 1 ~ n_pairs-2 (양쪽에 이웃 필요)
    log("\n2단계: 중간점 곡률 계산 (5개 수준)...")

    results = []
    valid_start = 2   # 좌측에 최소 2개 영점 필요
    valid_end = n_pairs - 2  # 우측에 최소 2개 영점 필요
    n_valid = valid_end - valid_start

    for idx, i in enumerate(range(valid_start, valid_end)):
        m = (zeros[i] + zeros[i + 1]) / 2.0
        gap = gaps[i]                 # Δ = γ_{i+1} - γ_i
        gap_prev = gaps[i - 1]        # Δ_{-1}
        gap_next = gaps[i + 1]        # Δ_{+1}

        # (A) κ_exact: 해석적 공식
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(m))
        L_val = connection_analytic(s)
        kappa_exact = float(abs(L_val) ** 2)

        # (B) κ_NN_removed: NN 쌍 (γ_i, γ_{i+1}) 기여 제거
        rho_n = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(zeros[i]))
        rho_n1 = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(zeros[i + 1]))
        rho_n_bar = mpmath.mpf('0.5') - 1j * mpmath.mpf(str(zeros[i]))
        rho_n1_bar = mpmath.mpf('0.5') - 1j * mpmath.mpf(str(zeros[i + 1]))

        nn_contrib = (1 / (s - rho_n) + 1 / (s - rho_n_bar)
                      + 1 / (s - rho_n1) + 1 / (s - rho_n1_bar))
        L_no_nn = L_val - nn_contrib
        kappa_nn_removed = float(abs(L_no_nn) ** 2)

        # (C) κ_NNN_only: 차차근접 (γ_{i-1}, γ_{i+2})만
        rho_prev = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(zeros[i - 1]))
        rho_next2 = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(zeros[i + 2]))
        rho_prev_bar = mpmath.mpf('0.5') - 1j * mpmath.mpf(str(zeros[i - 1]))
        rho_next2_bar = mpmath.mpf('0.5') - 1j * mpmath.mpf(str(zeros[i + 2]))

        nnn_contrib = (1 / (s - rho_prev) + 1 / (s - rho_prev_bar)
                       + 1 / (s - rho_next2) + 1 / (s - rho_next2_bar))
        kappa_nnn_only = float(abs(nnn_contrib) ** 2)

        # (D) κ_asymptotic: 고t 근사 공식
        term_prev = 1.0 / (gap_prev + gap / 2)
        term_next = 1.0 / (gap_next + gap / 2)
        kappa_asymptotic = (term_prev - term_next) ** 2

        # (E) κ_multi_zero: NN 제외, 좌우 5영점씩 합 (mirror 포함)
        multi_contrib = mpmath.mpc(0, 0)
        for j_offset in [-5, -4, -3, -2, -1, 2, 3, 4, 5, 6]:
            j = i + j_offset
            if 0 <= j < len(zeros):
                rho_j = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(zeros[j]))
                rho_j_bar = mpmath.mpf('0.5') - 1j * mpmath.mpf(str(zeros[j]))
                multi_contrib += 1 / (s - rho_j) + 1 / (s - rho_j_bar)
        kappa_multi = float(abs(multi_contrib) ** 2)

        results.append({
            'm': m, 'gap': gap, 'gap_prev': gap_prev, 'gap_next': gap_next,
            'kappa_exact': kappa_exact,
            'kappa_nn_removed': kappa_nn_removed,
            'kappa_nnn_only': kappa_nnn_only,
            'kappa_asymptotic': kappa_asymptotic,
            'kappa_multi': kappa_multi,
        })

        if (idx + 1) % 50 == 0:
            log(f"  [{idx+1}/{n_valid}] m={m:.2f}, κ_exact={kappa_exact:.4f}, "
                f"κ_asym={kappa_asymptotic:.6f}")

    log(f"\n  유효 데이터: {len(results)}개")

    # ── 3. 상관 분석 ──
    log("\n3단계: Spearman 상관 분석")
    log("=" * 70)

    ke = np.array([r['kappa_exact'] for r in results])
    knr = np.array([r['kappa_nn_removed'] for r in results])
    knnn = np.array([r['kappa_nnn_only'] for r in results])
    kasym = np.array([r['kappa_asymptotic'] for r in results])
    kmulti = np.array([r['kappa_multi'] for r in results])
    gap_arr = np.array([r['gap'] for r in results])
    gap_next_arr = np.array([r['gap_next'] for r in results])
    gap_prev_arr = np.array([r['gap_prev'] for r in results])

    correlations = {}

    pairs = [
        ("κ_exact vs κ_NN_removed", ke, knr, ">0.95", "NN 상쇄 확인"),
        ("κ_exact vs κ_NNN_only", ke, knnn, ">0.5", "NNN 지배 검증"),
        ("κ_exact vs κ_asymptotic", ke, kasym, ">0.5", "고t 공식 검증"),
        ("κ_exact vs κ_multi_zero", ke, kmulti, ">0.9", "Hadamard 지배"),
        ("κ_NN_removed vs κ_multi_zero", knr, kmulti, ">0.9", "다영점 재현"),
        ("κ_asymptotic vs gap_next", kasym, gap_next_arr, "<-0.3", "공식→비국소"),
        ("κ_exact vs gap", ke, gap_arr, "≈0", "NN 상쇄 (새니티)"),
        ("κ_exact vs gap_next", ke, gap_next_arr, "<-0.3", "비국소 (재현)"),
    ]

    for label, x, y, criterion, desc in pairs:
        rho_val, p_val = stats.spearmanr(x, y)
        status = "✅" if _check_criterion(rho_val, criterion) else "❌"
        correlations[label] = (rho_val, p_val, criterion, status, desc)
        log(f"  {status} {label}: ρ={rho_val:+.4f}, p={p_val:.2e} [{criterion}] — {desc}")

    # ── 4. R² 분석 (Pearson) ──
    log("\n4단계: Pearson R² (분산 설명률)")
    pearson_pairs = [
        ("κ_exact ~ κ_NN_removed", ke, knr),
        ("κ_exact ~ κ_NNN_only", ke, knnn),
        ("κ_exact ~ κ_asymptotic", ke, kasym),
        ("κ_exact ~ κ_multi_zero", ke, kmulti),
    ]
    for label, x, y in pearson_pairs:
        r, p = stats.pearsonr(x, y)
        log(f"  {label}: R²={r**2:.4f} (분산 {r**2*100:.1f}% 설명)")

    # ── 5. 기술 통계 ──
    log("\n5단계: 기술 통계")
    for name, arr in [("κ_exact", ke), ("κ_NN_removed", knr),
                       ("κ_NNN_only", knnn), ("κ_asymptotic", kasym),
                       ("κ_multi_zero", kmulti)]:
        log(f"  {name}: mean={np.mean(arr):.4e}, std={np.std(arr):.4e}, "
            f"median={np.median(arr):.4e}")

    # ── 6. t 의존성 점검 ──
    log("\n6단계: t-구간별 상관 (공식 안정성)")
    m_arr = np.array([r['m'] for r in results])
    for t_lo, t_hi in [(14, 100), (100, 200), (200, 350), (350, 560)]:
        mask = (m_arr >= t_lo) & (m_arr < t_hi)
        n_seg = np.sum(mask)
        if n_seg > 10:
            rho_seg, p_seg = stats.spearmanr(ke[mask], kasym[mask])
            rho_nn, p_nn = stats.spearmanr(ke[mask], knr[mask])
            log(f"  t∈[{t_lo},{t_hi}): n={n_seg}, "
                f"ρ(exact,asym)={rho_seg:+.3f} (p={p_seg:.2e}), "
                f"ρ(exact,NN_rm)={rho_nn:+.3f}")

    # ── 7. 최종 판정 ──
    log("\n" + "=" * 70)
    log("최종 판정")
    log("=" * 70)

    pass_count = sum(1 for v in correlations.values() if v[3] == "✅")
    total = len(correlations)
    log(f"  통과: {pass_count}/{total}")

    rho_asym = correlations["κ_exact vs κ_asymptotic"][0]
    rho_multi = correlations["κ_exact vs κ_multi_zero"][0]
    rho_nn_rm = correlations["κ_exact vs κ_NN_removed"][0]

    if rho_nn_rm > 0.95 and rho_multi > 0.9:
        log("  ★★★★ Hadamard 구조가 κ_mid 변동을 지배. NN 상쇄 + 다영점 재현.")
    elif rho_nn_rm > 0.95 and rho_asym > 0.5:
        log("  ★★★ NN 상쇄 확인 + 고t 공식 양성. NNN가 주요 기여.")
    elif rho_nn_rm > 0.9:
        log("  ★★ NN 상쇄 확인. 공식 정량적 정밀도 미달.")
    else:
        log("  ★ 예상과 다른 결과. 추가 분석 필요.")

    if rho_asym > 0.5:
        log(f"\n  → 고t 근사 공식 κ ≈ [1/(Δ₋₁+Δ/2) − 1/(Δ₊₁+Δ/2)]² 양성!")
        log(f"    Proposition 후보: Hadamard 중간점 비국소 정리")
    else:
        log(f"\n  → 고t 근사 공식 정밀도 부족 (ρ={rho_asym:.3f}). 수정 필요.")

    elapsed = time.time() - t_start
    log(f"\n소요 시간: {elapsed:.1f}초")

    # ── 결과 저장 ──
    with open(RESULT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("[C-246] 중간점 곡률의 Hadamard 예측 검증\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"N_zeros: {N_ZEROS}, 유효 쌍: {len(results)}, dps: {mpmath.mp.dps}\n")
        f.write(f"소요 시간: {elapsed:.1f}초\n")
        f.write("=" * 70 + "\n\n")

        f.write("Hadamard 예측 공식:\n")
        f.write("  κ_mid ≈ [1/(Δ_{-1}+Δ/2) - 1/(Δ_{+1}+Δ/2)]²\n")
        f.write("  (NN 쌍 상쇄 후 NNN 기여, 고t 근사)\n\n")

        f.write("Spearman 상관:\n")
        f.write(f"{'상관쌍':<40} {'ρ':>8} {'p':>12} {'기준':>8} {'판정':>4}\n")
        f.write("-" * 75 + "\n")
        for label, (rho_val, p_val, criterion, status, desc) in correlations.items():
            f.write(f"{label:<40} {rho_val:>+8.4f} {p_val:>12.2e} "
                    f"{criterion:>8} {status:>4}\n")

        f.write(f"\nPearson R² (분산 설명률):\n")
        for label, x, y in pearson_pairs:
            r, p = stats.pearsonr(x, y)
            f.write(f"  {label}: R²={r**2:.4f} ({r**2*100:.1f}%)\n")

        f.write(f"\n기술 통계:\n")
        for name, arr in [("κ_exact", ke), ("κ_NN_removed", knr),
                           ("κ_NNN_only", knnn), ("κ_asymptotic", kasym),
                           ("κ_multi_zero", kmulti)]:
            f.write(f"  {name}: mean={np.mean(arr):.4e}, std={np.std(arr):.4e}, "
                    f"median={np.median(arr):.4e}\n")

        f.write(f"\nt-구간별 ρ(exact,asym):\n")
        for t_lo, t_hi in [(14, 100), (100, 200), (200, 350), (350, 560)]:
            mask = (m_arr >= t_lo) & (m_arr < t_hi)
            n_seg = np.sum(mask)
            if n_seg > 10:
                rho_seg, _ = stats.spearmanr(ke[mask], kasym[mask])
                f.write(f"  t∈[{t_lo},{t_hi}): n={n_seg}, ρ={rho_seg:+.3f}\n")

        f.write(f"\n최종 판정: {pass_count}/{total} 통과\n")
        if rho_nn_rm > 0.95 and rho_multi > 0.9:
            f.write("★★★★ Hadamard 구조가 κ_mid 변동을 지배\n")
        elif rho_nn_rm > 0.95 and rho_asym > 0.5:
            f.write("★★★ NN 상쇄 + 고t 공식 양성\n")
        elif rho_nn_rm > 0.9:
            f.write("★★ NN 상쇄 확인. 공식 정밀도 미달\n")
        else:
            f.write("★ 예상 밖 결과\n")

        # 상세 데이터 (처음 30개)
        f.write(f"\n상세 데이터 (처음 30개):\n")
        f.write(f"{'m':>10} {'gap':>8} {'Δ-1':>8} {'Δ+1':>8} "
                f"{'κ_exact':>12} {'κ_NNrm':>12} {'κ_asym':>12} {'κ_multi':>12}\n")
        f.write("-" * 95 + "\n")
        for r in results[:30]:
            f.write(f"{r['m']:>10.3f} {r['gap']:>8.4f} {r['gap_prev']:>8.4f} "
                    f"{r['gap_next']:>8.4f} {r['kappa_exact']:>12.6f} "
                    f"{r['kappa_nn_removed']:>12.6f} {r['kappa_asymptotic']:>12.8f} "
                    f"{r['kappa_multi']:>12.6f}\n")

    log(f"\n결과 저장: {RESULT_FILE}")
    log("완료.")


def _check_criterion(rho, criterion):
    """기준 문자열에 따라 통과/실패 판정"""
    if criterion.startswith(">"):
        return rho > float(criterion[1:])
    elif criterion.startswith("<-"):
        return rho < float(criterion[1:])
    elif criterion.startswith("<"):
        return rho < float(criterion[1:])
    elif criterion == "≈0":
        return abs(rho) < 0.15
    return False


if __name__ == '__main__':
    main()
