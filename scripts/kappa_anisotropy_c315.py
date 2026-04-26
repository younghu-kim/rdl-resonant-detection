#!/usr/bin/env python3
"""
=============================================================================
[C-315] κ-비등방성 검증: σ-방향 vs t-방향 κδ² 비교
=============================================================================
목적: κ-비등방성 명제 수치 검증
  (a) σ-방향: κδ² = 1 + A_σ δ² + O(δ⁴)  [이차, 일차항 소멸]
  (b) t-방향: κδ² = 1 - 2Im(c₀)δ + B_t δ² + O(δ³)  [일차]
  (c) A_σ + B_t = 2Im(c₀)²  [일관성 검증]
  (d) Off-critical (DH): σ-방향에도 일차항 존재

이론 근거:
  Theorem 5 (Laurent 패리티) → Re(c₀) = 0 → σ-방향 일차항 소멸
  이것이 σ-국소화(FWHM < 0.004)의 이론적 메커니즘.
=============================================================================
"""

import sys
import numpy as np
import mpmath
from pathlib import Path

# 정밀도
mpmath.mp.dps = 80

sys.path.insert(0, str(Path(__file__).parent))
from bundle_utils import xi_func, connection_zeta

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DELTAS = [0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02,
          0.03, 0.05, 0.07, 0.1, 0.12, 0.15]
OUTFILE = Path(__file__).parent.parent / "results" / "kappa_anisotropy_c315.txt"


def compute_kappa_sigma(gamma, delta):
    """σ-방향: s = (1/2+δ) + iγ"""
    s = mpmath.mpf('0.5') + mpmath.mpf(str(delta)) + 1j * mpmath.mpf(str(gamma))
    L = connection_zeta(s)
    return float(abs(L)**2)


def compute_kappa_t(gamma, delta):
    """t-방향: s = 1/2 + i(γ+δ)"""
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(gamma + delta))
    L = connection_zeta(s)
    return float(abs(L)**2)


def fit_sigma(deltas, kd2_vals):
    """σ-방향 적합: κδ² - 1 = A_σ δ² + A₄ δ⁴
    일차항이 없어야 함 (Re(c₀)=0 예측)
    """
    y = np.array(kd2_vals) - 1.0
    d = np.array(deltas)
    # 모델 1: y = a₁δ + a₂δ² (일차항 포함)
    X_full = np.column_stack([d, d**2])
    coef_full, res_full, _, _ = np.linalg.lstsq(X_full, y, rcond=None)
    # 모델 2: y = a₂δ² (일차항 강제 0)
    X_quad = d.reshape(-1, 1)**2
    coef_quad, res_quad, _, _ = np.linalg.lstsq(X_quad, y, rcond=None)
    return {
        'c1_full': coef_full[0],  # 일차 계수 (0이어야 함)
        'A_sigma_full': coef_full[1],  # 이차 계수 (모델 1)
        'A_sigma_pure': coef_quad[0],  # 이차 계수 (모델 2)
        'R2_full': 1 - np.sum((y - X_full @ coef_full)**2) / np.sum(y**2),
        'R2_pure': 1 - np.sum((y - X_quad @ coef_quad)**2) / np.sum(y**2),
    }


def fit_t(deltas, kd2_vals):
    """t-방향 적합: κδ² - 1 = -2α δ + B_t δ²"""
    y = np.array(kd2_vals) - 1.0
    d = np.array(deltas)
    X = np.column_stack([d, d**2])
    coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    y_pred = X @ coef
    R2 = 1 - np.sum((y - y_pred)**2) / np.sum(y**2)
    return {
        'minus_2alpha': coef[0],  # -2Im(c₀)
        'B_t': coef[1],
        'R2': R2,
    }


def compute_laurent_c0(gamma, zeros_all):
    """c₀ = ξ''/ξ' / 2 의 직접 계산 대신, Hadamard 급수로 Im(c₀) 추정
    Im(c₀) ≈ -Σ_{j≠0} 1/(γ-γⱼ) + Im(digamma 항)
    """
    # 간단 추정: Hadamard 기여
    s_sum = 0.0
    for g in zeros_all:
        if abs(g - gamma) > 1e-6:
            s_sum += 1.0 / (gamma - g)
    # ψ(s/2) 기여 (순허수 부분)
    s0 = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(gamma))
    psi_val = mpmath.digamma(s0 / 2)
    im_psi = float(psi_val.imag) / 2
    # 1/s + 1/(s-1) 기여
    gs = 1.0 / complex(s0) + 1.0 / complex(s0 - 1)
    im_extra = gs.imag
    return -s_sum + im_psi + im_extra


def main():
    print("=" * 70)
    print("[C-315] κ-비등방성 검증: σ-방향 vs t-방향")
    print("=" * 70)
    print()

    # ζ(s) 영점 (확장 리스트)
    zeros_zeta = [float(mpmath.zetazero(n).imag) for n in range(1, 80)]
    test_zeros = zeros_zeta[:10]  # 처음 10개 검증

    lines = []
    lines.append("=" * 70)
    lines.append("[C-315] κ-비등방성 검증: σ-방향 vs t-방향")
    lines.append(f"날짜: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"영점 수: {len(test_zeros)}, δ 수: {len(DELTAS)}")
    lines.append("=" * 70)
    lines.append("")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # [1] σ-방향 측정
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("[1] σ-방향 κδ² 측정...")
    sigma_data = {}  # gamma -> list of (delta, kd2)
    for gi, g in enumerate(test_zeros):
        sigma_data[g] = []
        for d in DELTAS:
            kappa = compute_kappa_sigma(g, d)
            kd2 = kappa * d * d
            sigma_data[g].append((d, kd2))
        print(f"  γ_{gi+1} = {g:.6f} 완료")

    # σ-방향 적합
    lines.append("[1] σ-방향 적합: κδ² - 1 = c₁δ + A_σδ²")
    lines.append("  이론 예측: c₁ = 0 (Re(c₀)=0에서)")
    lines.append(f"  {'γ':>10s}  {'c₁(일차)':>12s}  {'A_σ(이차)':>12s}  {'R²(full)':>10s}  {'A_σ(pure)':>12s}  {'R²(pure)':>10s}")

    sigma_fits = {}
    for g in test_zeros:
        ds = [x[0] for x in sigma_data[g]]
        kd2s = [x[1] for x in sigma_data[g]]
        fit = fit_sigma(ds, kd2s)
        sigma_fits[g] = fit
        lines.append(f"  {g:10.4f}  {fit['c1_full']:+12.6f}  {fit['A_sigma_full']:+12.4f}  "
                      f"{fit['R2_full']:10.6f}  {fit['A_sigma_pure']:+12.4f}  {fit['R2_pure']:10.6f}")

    # c₁ 통계
    c1_vals = [sigma_fits[g]['c1_full'] for g in test_zeros]
    lines.append(f"\n  c₁ 통계: mean={np.mean(c1_vals):+.6f}, std={np.std(c1_vals):.6f}, "
                 f"|max|={max(abs(v) for v in c1_vals):.6f}")
    lines.append(f"  ★ c₁ ≈ 0 검증: {'PASS' if max(abs(v) for v in c1_vals) < 0.05 else 'FAIL'} "
                 f"(임계값 0.05)")
    lines.append("")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # [2] t-방향 측정
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n[2] t-방향 κδ² 측정...")
    t_data = {}
    for gi, g in enumerate(test_zeros):
        t_data[g] = []
        for d in DELTAS:
            kappa = compute_kappa_t(g, d)
            kd2 = kappa * d * d
            t_data[g].append((d, kd2))
        print(f"  γ_{gi+1} = {g:.6f} 완료")

    # t-방향 적합
    lines.append("[2] t-방향 적합: κδ² - 1 = (-2α)δ + B_tδ²")
    lines.append(f"  {'γ':>10s}  {'-2α(일차)':>12s}  {'B_t(이차)':>12s}  {'R²':>10s}")

    t_fits = {}
    for g in test_zeros:
        ds = [x[0] for x in t_data[g]]
        kd2s = [x[1] for x in t_data[g]]
        fit = fit_t(ds, kd2s)
        t_fits[g] = fit
        lines.append(f"  {g:10.4f}  {fit['minus_2alpha']:+12.6f}  {fit['B_t']:+12.4f}  {fit['R2']:10.6f}")

    # t-방향 일차 계수 통계
    m2a_vals = [t_fits[g]['minus_2alpha'] for g in test_zeros]
    lines.append(f"\n  -2α 통계: mean={np.mean(m2a_vals):+.4f}, std={np.std(m2a_vals):.4f}")
    lines.append(f"  ★ |−2α| >> 0 검증: {'PASS' if min(abs(v) for v in m2a_vals) > 0.5 else 'FAIL'}")
    lines.append("")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # [3] 일관성 검증: A_σ + B_t = 2Im(c₀)²
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n[3] 일관성 검증...")
    lines.append("[3] 일관성: A_σ + B_t = 2α²  (α = Im(c₀) = -(-2α)/2)")
    lines.append(f"  {'γ':>10s}  {'A_σ':>10s}  {'B_t':>10s}  {'A_σ+B_t':>10s}  {'2α²':>10s}  {'비율':>8s}")

    consistency_ratios = []
    for g in test_zeros:
        A_s = sigma_fits[g]['A_sigma_full']
        B_t = t_fits[g]['B_t']
        alpha = -t_fits[g]['minus_2alpha'] / 2  # Im(c₀)
        sum_AB = A_s + B_t
        two_alpha2 = 2 * alpha**2
        ratio = sum_AB / two_alpha2 if abs(two_alpha2) > 1e-10 else float('nan')
        consistency_ratios.append(ratio)
        lines.append(f"  {g:10.4f}  {A_s:+10.4f}  {B_t:+10.4f}  {sum_AB:+10.4f}  {two_alpha2:10.4f}  {ratio:8.4f}")

    valid_ratios = [r for r in consistency_ratios if not np.isnan(r)]
    mean_ratio = np.mean(valid_ratios)
    lines.append(f"\n  비율 통계: mean={mean_ratio:.4f}, std={np.std(valid_ratios):.4f}")
    lines.append(f"  ★ A_σ+B_t = 2α² 검증: {'PASS' if abs(mean_ratio - 1.0) < 0.1 else 'FAIL'} "
                 f"(|mean-1| < 0.1)")
    lines.append("")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # [4] 비등방성 비율 (핵심 지표)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n[4] 비등방성 비율...")
    lines.append("[4] 비등방성: σ-방향 c₁ vs t-방향 c₁")
    lines.append("  σ-방향 c₁ ≈ 0 (이론), t-방향 c₁ = -2α ≠ 0 (이론)")
    lines.append(f"  {'γ':>10s}  {'σ-c₁':>12s}  {'t-c₁':>12s}  {'|σ/t| 비율':>12s}")

    aniso_ratios = []
    for g in test_zeros:
        sc1 = sigma_fits[g]['c1_full']
        tc1 = t_fits[g]['minus_2alpha']
        ratio = abs(sc1 / tc1) if abs(tc1) > 1e-10 else float('nan')
        aniso_ratios.append(ratio)
        lines.append(f"  {g:10.4f}  {sc1:+12.6f}  {tc1:+12.6f}  {ratio:12.6f}")

    valid_aniso = [r for r in aniso_ratios if not np.isnan(r)]
    lines.append(f"\n  |σ-c₁/t-c₁| 비율: mean={np.mean(valid_aniso):.6f}, max={max(valid_aniso):.6f}")
    lines.append(f"  ★ σ-방향 일차항 억제: {'PASS' if np.mean(valid_aniso) < 0.05 else 'FAIL'} "
                 f"(비율 < 0.05)")
    lines.append("")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # [5] FWHM 비교
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n[5] FWHM 추정...")
    lines.append("[5] FWHM 추정 (κδ² = 0.5 지점)")
    lines.append(f"  {'γ':>10s}  {'FWHM_σ':>10s}  {'FWHM_t':>10s}  {'σ/t 비율':>10s}")

    fwhm_ratios = []
    for g in test_zeros:
        # σ-FWHM: 1 + A_σ δ² = 0.5 → δ = √(0.5/|A_σ|)
        A_s = sigma_fits[g]['A_sigma_pure']
        if A_s < 0:
            fwhm_s = 2 * np.sqrt(0.5 / abs(A_s))
        else:
            fwhm_s = float('nan')

        # t-FWHM: 1 - 2|α|δ = 0.5 → δ = 0.25/|α|
        alpha = abs(t_fits[g]['minus_2alpha']) / 2
        if alpha > 1e-10:
            fwhm_t = 2 * 0.5 / (2 * alpha)  # 양쪽이므로 ×2... 아니, t-방향은 비대칭
            # 정확하게: κδ²=0.5는 -2αδ = -0.5 → δ=0.25/α (한 방향)
            fwhm_t = 0.5 / alpha  # 대칭 가정 양쪽 합
        else:
            fwhm_t = float('nan')

        ratio = fwhm_s / fwhm_t if (not np.isnan(fwhm_s) and not np.isnan(fwhm_t)) else float('nan')
        fwhm_ratios.append(ratio)
        lines.append(f"  {g:10.4f}  {fwhm_s:10.4f}  {fwhm_t:10.4f}  {ratio:10.4f}")

    valid_fwhm = [r for r in fwhm_ratios if not np.isnan(r)]
    if valid_fwhm:
        lines.append(f"\n  FWHM_σ/FWHM_t: mean={np.mean(valid_fwhm):.4f}, std={np.std(valid_fwhm):.4f}")
    lines.append("")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # [6] C-300 데이터와 일관성 확인
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    print("\n[6] C-300 t-방향 데이터와 비교...")
    # C-300의 c₁ 적합값 (t-방향)
    c300_c1 = {
        14.1347: -1.1599, 21.0220: -1.5117, 25.0109: -1.1735,
        30.4249: -1.9148, 32.9351: -0.9342, 37.5862: -1.6809,
        40.9187: -1.7201, 43.3271: -0.9026, 48.0052: -2.4057,
        49.7738: -0.8572,
    }
    lines.append("[6] C-300 (t-방향) vs C-315 (t-방향) 일관성")
    lines.append(f"  {'γ':>10s}  {'C-300 c₁':>12s}  {'C-315 -2α':>12s}  {'차이':>10s}")
    for g in test_zeros:
        g_key = min(c300_c1.keys(), key=lambda k: abs(k - g))
        c300 = c300_c1[g_key]
        c315 = t_fits[g]['minus_2alpha']
        diff = c315 - c300
        lines.append(f"  {g:10.4f}  {c300:+12.4f}  {c315:+12.4f}  {diff:+10.4f}")

    lines.append("")

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # [7] 종합 판정
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    lines.append("=" * 70)
    lines.append("[판정]")

    pass_count = 0
    total = 4

    # P1: σ-방향 c₁ ≈ 0
    p1 = max(abs(v) for v in c1_vals) < 0.05
    pass_count += int(p1)
    lines.append(f"  P1. σ-방향 일차항 소멸 (|c₁|<0.05): {'✅ PASS' if p1 else '❌ FAIL'}")

    # P2: t-방향 |c₁| >> 0
    p2 = min(abs(v) for v in m2a_vals) > 0.5
    pass_count += int(p2)
    lines.append(f"  P2. t-방향 일차항 존재 (|c₁|>0.5): {'✅ PASS' if p2 else '❌ FAIL'}")

    # P3: 일관성 A_σ+B_t = 2α²
    p3 = abs(mean_ratio - 1.0) < 0.15
    pass_count += int(p3)
    lines.append(f"  P3. 일관성 (A_σ+B_t)/2α² ≈ 1 (|r-1|<0.15): {'✅ PASS' if p3 else '❌ FAIL'}")

    # P4: 비등방성 비율 σ/t < 0.05
    p4 = np.mean(valid_aniso) < 0.05 if valid_aniso else False
    pass_count += int(p4)
    lines.append(f"  P4. 비등방성 비율 <0.05: {'✅ PASS' if p4 else '❌ FAIL'}")

    lines.append(f"\n  종합: {pass_count}/{total} PASS")
    if pass_count == total:
        lines.append("  ★★★★★ κ-비등방성 명제 완전 검증")
        lines.append("  σ-국소화 메커니즘: FE+CC → Re(c₀)=0 → σ-방향 이차 → 날카로운 피크")
    elif pass_count >= 3:
        lines.append("  ★★★★ 강한 양성 (부분 검증)")
    else:
        lines.append("  검증 실패 — 재분석 필요")

    lines.append("=" * 70)

    # 출력 및 저장
    output = "\n".join(lines)
    print("\n" + output)
    OUTFILE.parent.mkdir(parents=True, exist_ok=True)
    OUTFILE.write_text(output, encoding='utf-8')
    print(f"\n결과 저장: {OUTFILE}")


if __name__ == "__main__":
    main()
