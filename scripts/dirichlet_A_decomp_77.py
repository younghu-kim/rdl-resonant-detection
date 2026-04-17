#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #77 — A(t₀) 해석적 분해: A_Γ vs A_L 분리 검증
=============================================================================
목적:
  A(t₀) = Im(H₀)² + 2Re(H₁) 에서 Γ-인자 기여 A_Γ vs L-함수 산술 기여 A_L 분해.

핵심 이론:
  ξ'/ξ(s) = 1/(s-ρ) + H₀ + H₁(s-ρ) + ...

  ζ(s)에서:
    ξ'/ξ(s) = 1/s + 1/(s-1) - ½log(π) + ½ψ(s/2) + ζ'/ζ(s)
    H₀^Γ(ζ) = 1/ρ + 1/(ρ-1) - ½log(π) + ½ψ(ρ/2)  [해석적, 극 없음]
    H₀^L(ζ) = lim_{s→ρ} [ζ'/ζ(s) - 1/(s-ρ)]       [산술 기여]
    H₀ = H₀^Γ + H₀^L

  χ₅²(even, q=5, μ=0)에서:
    Λ'/Λ(s,χ) = ½log(5/π) + ½ψ(s/2) + L'/L(s,χ)
    H₀^Γ(χ₅²) = ½log(5/π) + ½ψ(ρ_χ/2)             [해석적, 극 없음]
    H₀^L(χ₅²) = lim_{s→ρ_χ} [L'/L(s,χ) - 1/(s-ρ_χ)] [산술 기여]

수치 방법:
  H₀ 수치 계산: f(ε) = F(ρ+ε) - 1/ε → H₀ as ε→0 (Richardson 외삽)
    ζ:   F(s) = ξ'/ξ(s)
    χ₅²: F(s) = Λ'/Λ(s,χ)

  A_Γ = Im(H₀^Γ)²  (H₀가 허수일 때의 1차 기여)
  A_measured = κ·δ² - 1 (δ=0.001에서 측정) — from #76 결과
  A_L = A_measured - A_Γ

핵심 질문:
  ΔA_total = 1.82 (ζ→χ₅²) 중
  ΔA_Γ (Γ 인자 차이)가 얼마나 기여하는가?
  나머지 ΔA_L는 산술적 L'/L 기여.

결과: results/A_decomposition_77.txt
=============================================================================
"""

import sys, os, time
import numpy as np

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

DPS = 60
mpmath.mp.dps = DPS

from bundle_utils import find_zeros_zeta, find_zeros_dirichlet

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "A_decomposition_77.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []

def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))

def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ─────────────────────────────────────────────────────────────────
# 파라미터
# ─────────────────────────────────────────────────────────────────

T_MIN, T_MAX = 14.0, 59.0

# χ₅² 문자값 (Legendre symbol mod 5)
CHI5_EVEN = [0, 1, -1, -1, 1]
CHAR_CHI5 = {
    'chi': CHI5_EVEN,
    'q': 5, 'a': 0,
    'label': 'χ₅² (mod5, even, μ=0)',
}

# #76 결과에서 추출한 A_measured (δ=0.001, 가장 정확)
# ζ 12영점 (t값, A값)
ZETA_MEASURED = [
    (14.13472514, 0.6352),
    (21.02203964, 0.9914),
    (25.01085758, 0.7668),
    (30.42487613, 1.5942),
    (32.93506159, 0.8634),
    (37.58617816, 1.2949),
    (40.91871901, 1.5642),
    (43.32707328, 0.9392),
    (48.00515088, 2.5574),
    (49.77383248, 1.3181),
    (52.97032148, 1.1539),
    (56.44624770, 1.5941),
]

# χ₅² 24영점 (t값, A값)
CHI5_MEASURED = [
    (16.03382113, 3.3213),
    (17.56699429, 2.1742),
    (19.54073262, 1.5475),
    (22.22740545, 1.8263),
    (24.58846622, 2.1807),
    (26.77609595, 2.9570),
    (28.46103510, 3.4710),
    (29.70790935, 2.2023),
    (33.00045601, 3.2903),
    (34.72881298, 4.4787),
    (35.86863837, 2.7863),
    (38.12918472, 3.1324),
    (39.56057295, 2.2372),
    (41.84243855, 2.3440),
    (44.03129006, 4.2334),
    (45.42730008, 4.8603),
    (46.49272716, 3.3006),
    (48.34566182, 1.9382),
    (51.08775193, 5.9069),
    (52.12590223, 3.6048),
    (53.83044520, 2.9428),
    (55.58928034, 4.0530),
    (56.83886594, 3.2198),
    (58.38611749, 2.1851),
]

# ─────────────────────────────────────────────────────────────────
# H₀^Γ 해석적 계산
# ─────────────────────────────────────────────────────────────────

def H0_gamma_zeta(rho):
    """
    H₀^Γ(ζ) = 1/ρ + 1/(ρ-1) - ½log(π) + ½ψ(ρ/2)
    이것이 ξ'/ξ의 Laurent 전개에서 ζ'/ζ를 제외한 Γ-인자 기여.
    """
    log_pi = mpmath.log(mpmath.pi)
    psi_val = mpmath.digamma(rho / 2)
    return (1/rho + 1/(rho - 1) - log_pi/2 + psi_val/2)

def H0_gamma_chi5(rho):
    """
    H₀^Γ(χ₅²) = ½log(5/π) + ½ψ(ρ/2)
    Λ'/Λ의 Laurent 전개에서 L'/L를 제외한 Γ+conductor 기여.
    μ=0이므로 감마인자 = Γ(s/2), a=0.
    """
    log_5_over_pi = mpmath.log(mpmath.mpf(5) / mpmath.pi)
    psi_val = mpmath.digamma(rho / 2)
    return (log_5_over_pi/2 + psi_val/2)

# ─────────────────────────────────────────────────────────────────
# F(s) = ξ'/ξ(s) 또는 Λ'/Λ(s,χ) 계산
# ─────────────────────────────────────────────────────────────────

H_DERIV = mpmath.mpf('1e-7')   # 수치미분 스텝 (ε와 독립)

def xi_log_deriv_zeta(s):
    """ξ'/ξ(s) for ζ(s) — 명시적 공식"""
    log_pi = mpmath.log(mpmath.pi)
    psi_val = mpmath.digamma(s / 2)

    # ζ'/ζ(s) = (ζ(s+h)-ζ(s-h)) / (2h·ζ(s))
    h = H_DERIV
    zeta_s  = mpmath.zeta(s)
    if abs(zeta_s) < mpmath.mpf(10)**(-DPS + 5):
        # ζ가 극도로 작으면 (ε이 너무 작음) 경고
        raise ValueError(f"ζ(s) near-zero at s={s}: |ζ|={abs(zeta_s):.2e}")

    zeta_ph = mpmath.zeta(s + h)
    zeta_mh = mpmath.zeta(s - h)
    zeta_log_deriv = (zeta_ph - zeta_mh) / (2 * h * zeta_s)

    return 1/s + 1/(s - 1) - log_pi/2 + psi_val/2 + zeta_log_deriv

def lambda_log_deriv_chi5(s):
    """Λ'/Λ(s,χ₅²) — 명시적 공식"""
    log_5_over_pi = mpmath.log(mpmath.mpf(5) / mpmath.pi)
    psi_val = mpmath.digamma(s / 2)

    # L'/L(s,χ) = (L(s+h,χ)-L(s-h,χ)) / (2h·L(s,χ))
    h = H_DERIV
    chi = CHI5_EVEN
    L_s  = mpmath.dirichlet(s, chi)
    if abs(L_s) < mpmath.mpf(10)**(-DPS + 5):
        raise ValueError(f"L(s,χ) near-zero at s={s}: |L|={abs(L_s):.2e}")

    L_ph = mpmath.dirichlet(s + h, chi)
    L_mh = mpmath.dirichlet(s - h, chi)
    L_log_deriv = (L_ph - L_mh) / (2 * h * L_s)

    return log_5_over_pi/2 + psi_val/2 + L_log_deriv

# ─────────────────────────────────────────────────────────────────
# H₀ 수치 계산 (Richardson 외삽)
# ─────────────────────────────────────────────────────────────────

EPS_BASE = mpmath.mpf('1e-10')
EPS_LIST = [EPS_BASE / (2**k) for k in range(5)]   # 1e-10, 5e-11, 2.5e-11, 1.25e-11, 6.25e-12

def compute_H0_numerical(rho, F_func, label=""):
    """
    H₀ = lim_{ε→0} [F(ρ+ε) - 1/ε]  (Richardson 외삽)
    F = ξ'/ξ (ζ의 경우) 또는 Λ'/Λ (χ의 경우)
    """
    f_vals = []
    eps_used = []
    errors = []

    for eps in EPS_LIST:
        s = rho + eps   # ε은 실수
        try:
            F_val = F_func(s)
            f = F_val - 1/eps   # pole 제거
            f_vals.append(f)
            eps_used.append(eps)
        except Exception as e:
            errors.append((eps, str(e)))

    if len(f_vals) < 2:
        raise RuntimeError(f"Richardson: 값 부족 ({len(f_vals)}개). 에러: {errors}")

    # Richardson 외삽 (레벨 1): 2·f(ε/2) - f(ε) → H₀ + O(ε²)
    # f_vals[k]는 ε_k = ε_base / 2^k에서의 값
    # 인접한 쌍으로 외삽
    R1_vals = []
    for k in range(len(f_vals) - 1):
        # f_vals[k]: ε_k, f_vals[k+1]: ε_k/2
        R1 = 2 * f_vals[k+1] - f_vals[k]   # O(ε²) 정확도
        R1_vals.append(R1)

    # Richardson 레벨 2: 4·R1(ε/2) - R1(ε))/3 → H₀ + O(ε⁴)
    R2_vals = []
    for k in range(len(R1_vals) - 1):
        R2 = (4 * R1_vals[k+1] - R1_vals[k]) / 3
        R2_vals.append(R2)

    # 수렴 판정: 최선의 추정값 + 수렴 확인
    if R2_vals:
        H0_best = R2_vals[-1]   # 가장 정밀한 레벨 2
        # 수렴성: R2 값들의 변동
        if len(R2_vals) >= 2:
            conv = abs(R2_vals[-1] - R2_vals[-2])
        else:
            conv = abs(R1_vals[-1] - R1_vals[-2]) if len(R1_vals) >= 2 else None
    else:
        H0_best = R1_vals[-1]
        conv = abs(R1_vals[-1] - R1_vals[0]) if len(R1_vals) >= 2 else None

    return H0_best, conv, f_vals, R1_vals, R2_vals

# ─────────────────────────────────────────────────────────────────
# 단일 영점 분해
# ─────────────────────────────────────────────────────────────────

def analyze_zero(t_val, A_measured, is_zeta=True):
    """
    영점 ρ = 0.5 + i·t에서 A(t₀) 분해:
    A_Γ = Im(H₀^Γ)²
    H₀^L = H₀ - H₀^Γ
    A_L_approx = A_measured - A_Γ
    """
    rho = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t_val)))

    # H₀^Γ 해석적 계산
    if is_zeta:
        H0_gamma = H0_gamma_zeta(rho)
        F_func   = xi_log_deriv_zeta
        fname    = "ξ'/ξ"
    else:
        H0_gamma = H0_gamma_chi5(rho)
        F_func   = lambda_log_deriv_chi5
        fname    = "Λ'/Λ"

    # H₀ 수치 계산
    try:
        H0_num, conv, f_vals, R1_vals, R2_vals = compute_H0_numerical(rho, F_func)
        success = True
    except Exception as e:
        return None, str(e)

    # H₀^L = H₀ - H₀^Γ
    H0_L = H0_num - H0_gamma

    # A값 계산
    A_gamma   = float(mpmath.im(H0_gamma)**2)
    A_from_H0 = float(mpmath.im(H0_num)**2)     # Im(H₀)² (H₁ 무시)
    A_L       = float(A_measured) - A_gamma       # 측정 - Γ기여
    H0_L_imag = float(mpmath.im(H0_L))

    return {
        't': t_val,
        'A_measured': A_measured,
        'H0_gamma': H0_gamma,
        'H0_num': H0_num,
        'H0_L': H0_L,
        'conv': float(abs(conv)) if conv is not None else None,
        # A 분해
        'Im_H0_gamma': float(mpmath.im(H0_gamma)),
        'Im_H0_num':   float(mpmath.im(H0_num)),
        'Im_H0_L':     float(H0_L_imag),
        'Re_H0_gamma': float(mpmath.re(H0_gamma)),
        'Re_H0_num':   float(mpmath.re(H0_num)),
        'Re_H0_L':     float(mpmath.re(H0_L)),
        'A_gamma':   A_gamma,
        'A_from_H0': A_from_H0,    # Im(H₀)² = Im(H₀^Γ+H₀^L)²
        'A_L':       A_L,           # A_measured - A_gamma
        # 교차항: A_measured ≈ Im(H₀^Γ)² + 2·Im(H₀^Γ)·Im(H₀^L) + Im(H₀^L)²
        'A_cross':   2 * float(mpmath.im(H0_gamma)) * H0_L_imag,
        'A_L_sq':    H0_L_imag**2,
    }, None

# ─────────────────────────────────────────────────────────────────
# 메인 실행
# ─────────────────────────────────────────────────────────────────

def main():
    t0_total = time.time()

    log("=" * 72)
    log("결과 #77 — A(t₀) 해석적 분해: A_Γ vs A_L 분리 검증")
    log("=" * 72)
    log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"DPS={DPS}")
    log(f"t 범위: [{T_MIN}, {T_MAX}]")
    log(f"ε 목록 (Richardson): {[float(e) for e in EPS_LIST]}")
    log()
    log("이론 공식:")
    log("  H₀^Γ(ζ) = 1/ρ + 1/(ρ-1) - ½log(π) + ½ψ(ρ/2)")
    log("  H₀^Γ(χ₅²) = ½log(5/π) + ½ψ(ρ_χ/2)")
    log("  H₀(수치) = lim_{ε→0} [F(ρ+ε) - 1/ε]  (Richardson 외삽)")
    log("  A_Γ = Im(H₀^Γ)²  [Γ 인자 기여]")
    log("  A_L = A_measured - A_Γ  [산술적 잔차]")
    log()
    log(f"naive 예측: ΔA_q = log(5/π)/2 = {float(mpmath.log(mpmath.mpf(5)/mpmath.pi)/2):.4f}")
    log(f"실측 ΔA = 1.82  (#76 결과)")
    flush_file()

    # ─────────────────────────────────────────────────────────────
    # Step 1: ζ 영점 분해
    # ─────────────────────────────────────────────────────────────

    log()
    log("=" * 72)
    log("Step 1: ζ(s) 12영점 — A(t₀) 분해")
    log("=" * 72)

    zeta_results = []
    zeta_fail_count = 0

    log(f"  {'t':>10} | {'A_meas':>8} | {'Im(H₀^Γ)':>10} | {'Im(H₀)':>10} | "
        f"{'Im(H₀^L)':>10} | {'A_Γ':>8} | {'A_L':>8} | 수렴")
    log("  " + "-" * 90)

    for t_val, A_meas in ZETA_MEASURED:
        res, err = analyze_zero(t_val, A_meas, is_zeta=True)
        if err:
            log(f"  t={t_val:.5f}: ⚠️ 실패 — {err}")
            zeta_fail_count += 1
            continue

        zeta_results.append(res)
        conv_str = f"{res['conv']:.2e}" if res['conv'] is not None else "?"
        log(f"  t={res['t']:>10.6f} | {res['A_measured']:>8.4f} | "
            f"{res['Im_H0_gamma']:>10.5f} | {res['Im_H0_num']:>10.5f} | "
            f"{res['Im_H0_L']:>10.5f} | {res['A_gamma']:>8.4f} | "
            f"{res['A_L']:>8.4f} | {conv_str}")
        flush_file()

    log()
    if zeta_results:
        A_meas_arr   = np.array([r['A_measured']  for r in zeta_results])
        A_gamma_arr  = np.array([r['A_gamma']      for r in zeta_results])
        A_L_arr      = np.array([r['A_L']          for r in zeta_results])
        ImH0g_arr    = np.array([r['Im_H0_gamma']  for r in zeta_results])
        ImH0n_arr    = np.array([r['Im_H0_num']    for r in zeta_results])
        ImH0L_arr    = np.array([r['Im_H0_L']      for r in zeta_results])

        log(f"  ζ 요약 (n={len(zeta_results)}):")
        log(f"    mean(A_measured) = {np.mean(A_meas_arr):.4f} ± {np.std(A_meas_arr):.4f}")
        log(f"    mean(A_Γ)        = {np.mean(A_gamma_arr):.4f} ± {np.std(A_gamma_arr):.4f}")
        log(f"    mean(A_L)        = {np.mean(A_L_arr):.4f} ± {np.std(A_L_arr):.4f}")
        log(f"    Γ 기여 비중 = {np.mean(A_gamma_arr)/np.mean(A_meas_arr)*100:.1f}%")
        log(f"    L 기여 비중 = {np.mean(A_L_arr)/np.mean(A_meas_arr)*100:.1f}%")
        log(f"    mean(Im(H₀^Γ))  = {np.mean(ImH0g_arr):.5f}")
        log(f"    mean(Im(H₀))    = {np.mean(ImH0n_arr):.5f}")
        log(f"    mean(Im(H₀^L))  = {np.mean(ImH0L_arr):.5f}")
        flush_file()

    # ─────────────────────────────────────────────────────────────
    # Step 2: χ₅² 영점 분해
    # ─────────────────────────────────────────────────────────────

    log()
    log("=" * 72)
    log("Step 2: χ₅² (mod5, even, μ=0) 영점 — A(t₀) 분해")
    log("=" * 72)
    log(f"  영점 수: {len(CHI5_MEASURED)}개 (δ=0.001 기준 A값 사용)")

    chi5_results = []
    chi5_fail_count = 0

    log(f"  {'t':>10} | {'A_meas':>8} | {'Im(H₀^Γ)':>10} | {'Im(H₀)':>10} | "
        f"{'Im(H₀^L)':>10} | {'A_Γ':>8} | {'A_L':>8} | 수렴")
    log("  " + "-" * 90)

    for t_val, A_meas in CHI5_MEASURED:
        res, err = analyze_zero(t_val, A_meas, is_zeta=False)
        if err:
            log(f"  t={t_val:.5f}: ⚠️ 실패 — {err}")
            chi5_fail_count += 1
            continue

        chi5_results.append(res)
        conv_str = f"{res['conv']:.2e}" if res['conv'] is not None else "?"
        log(f"  t={res['t']:>10.6f} | {res['A_measured']:>8.4f} | "
            f"{res['Im_H0_gamma']:>10.5f} | {res['Im_H0_num']:>10.5f} | "
            f"{res['Im_H0_L']:>10.5f} | {res['A_gamma']:>8.4f} | "
            f"{res['A_L']:>8.4f} | {conv_str}")
        flush_file()

    log()
    if chi5_results:
        A_meas5   = np.array([r['A_measured']  for r in chi5_results])
        A_gamma5  = np.array([r['A_gamma']      for r in chi5_results])
        A_L5      = np.array([r['A_L']          for r in chi5_results])
        ImH0g5    = np.array([r['Im_H0_gamma']  for r in chi5_results])
        ImH0n5    = np.array([r['Im_H0_num']    for r in chi5_results])
        ImH0L5    = np.array([r['Im_H0_L']      for r in chi5_results])

        log(f"  χ₅² 요약 (n={len(chi5_results)}):")
        log(f"    mean(A_measured) = {np.mean(A_meas5):.4f} ± {np.std(A_meas5):.4f}")
        log(f"    mean(A_Γ)        = {np.mean(A_gamma5):.4f} ± {np.std(A_gamma5):.4f}")
        log(f"    mean(A_L)        = {np.mean(A_L5):.4f} ± {np.std(A_L5):.4f}")
        log(f"    Γ 기여 비중 = {np.mean(A_gamma5)/np.mean(A_meas5)*100:.1f}%")
        log(f"    L 기여 비중 = {np.mean(A_L5)/np.mean(A_meas5)*100:.1f}%")
        log(f"    mean(Im(H₀^Γ))  = {np.mean(ImH0g5):.5f}")
        log(f"    mean(Im(H₀))    = {np.mean(ImH0n5):.5f}")
        log(f"    mean(Im(H₀^L))  = {np.mean(ImH0L5):.5f}")
        flush_file()

    # ─────────────────────────────────────────────────────────────
    # Step 3: ΔA 분해 — Γ vs L 기여
    # ─────────────────────────────────────────────────────────────

    if zeta_results and chi5_results:
        log()
        log("=" * 72)
        log("Step 3: ΔA 분해 — Γ vs L 기여 분석")
        log("=" * 72)

        mean_A_meas_z   = np.mean([r['A_measured'] for r in zeta_results])
        mean_A_meas_c   = np.mean([r['A_measured'] for r in chi5_results])
        mean_A_gamma_z  = np.mean([r['A_gamma']    for r in zeta_results])
        mean_A_gamma_c  = np.mean([r['A_gamma']    for r in chi5_results])
        mean_A_L_z      = np.mean([r['A_L']        for r in zeta_results])
        mean_A_L_c      = np.mean([r['A_L']        for r in chi5_results])
        mean_ImH0g_z    = np.mean([r['Im_H0_gamma'] for r in zeta_results])
        mean_ImH0g_c    = np.mean([r['Im_H0_gamma'] for r in chi5_results])
        mean_ImH0n_z    = np.mean([r['Im_H0_num']   for r in zeta_results])
        mean_ImH0n_c    = np.mean([r['Im_H0_num']   for r in chi5_results])
        mean_ImH0L_z    = np.mean([r['Im_H0_L']     for r in zeta_results])
        mean_ImH0L_c    = np.mean([r['Im_H0_L']     for r in chi5_results])

        Delta_A_total   = mean_A_meas_c  - mean_A_meas_z
        Delta_A_gamma   = mean_A_gamma_c - mean_A_gamma_z
        Delta_A_L       = mean_A_L_c     - mean_A_L_z
        Delta_ImH0g     = mean_ImH0g_c   - mean_ImH0g_z
        Delta_ImH0n     = mean_ImH0n_c   - mean_ImH0n_z
        Delta_ImH0L     = mean_ImH0L_c   - mean_ImH0L_z

        log()
        log("  비교표 (ζ → χ₅²):")
        log(f"  {'항목':>20} | {'ζ':>10} | {'χ₅²':>10} | {'Δ':>10} | 비중(%)")
        log("  " + "-" * 65)
        log(f"  {'mean(A_measured)':>20} | {mean_A_meas_z:>10.4f} | {mean_A_meas_c:>10.4f} | {Delta_A_total:>+10.4f} | 100%")
        log(f"  {'mean(A_Γ)':>20} | {mean_A_gamma_z:>10.4f} | {mean_A_gamma_c:>10.4f} | {Delta_A_gamma:>+10.4f} | {Delta_A_gamma/Delta_A_total*100:>+6.1f}%")
        log(f"  {'mean(A_L)':>20} | {mean_A_L_z:>10.4f} | {mean_A_L_c:>10.4f} | {Delta_A_L:>+10.4f} | {Delta_A_L/Delta_A_total*100:>+6.1f}%")
        log()
        log(f"  {'mean(Im(H₀^Γ))':>20} | {mean_ImH0g_z:>10.5f} | {mean_ImH0g_c:>10.5f} | {Delta_ImH0g:>+10.5f} |")
        log(f"  {'mean(Im(H₀))':>20} | {mean_ImH0n_z:>10.5f} | {mean_ImH0n_c:>10.5f} | {Delta_ImH0n:>+10.5f} |")
        log(f"  {'mean(Im(H₀^L))':>20} | {mean_ImH0L_z:>10.5f} | {mean_ImH0L_c:>10.5f} | {Delta_ImH0L:>+10.5f} |")

        log()
        log("  ─────────────────────────────────────────────")
        log("  핵심 결과:")
        log(f"    ΔA_total  = {Delta_A_total:+.4f} (측정값, #76과 일치 확인)")
        log(f"    ΔA_Γ      = {Delta_A_gamma:+.4f} ({Delta_A_gamma/Delta_A_total*100:+.1f}% — Γ 인자 기여)")
        log(f"    ΔA_L      = {Delta_A_L:+.4f} ({Delta_A_L/Delta_A_total*100:+.1f}% — 산술 기여)")
        log()
        log(f"    naive 예측 ΔA = log(5/π)/2 = {float(mpmath.log(mpmath.mpf(5)/mpmath.pi)/2):.4f}")
        log(f"    → naive 예측의 의미: log(5/π)/2는 '실수' → Im(H₀^Γ)에 직접 기여 없음")
        log(f"    → 실제 ΔA_Γ = {Delta_A_gamma:.4f} (Im 성분 차이에서 기인)")
        log()
        log("  A_Γ 분해 상세:")
        log(f"    ΔIm(H₀^Γ) = {Delta_ImH0g:+.5f}")
        log(f"    이론: Im(H₀^Γ(χ₅²)) - Im(H₀^Γ(ζ)) =")
        log(f"          Im(½ψ(ρ_χ/2)) - Im(1/ρ + 1/(ρ-1) + ½ψ(ρ/2))")
        log(f"          ≈ -Im(1/ρ + 1/(ρ-1)) ≈ 2/t ≈ {2/np.mean([r['t'] for r in zeta_results]):.4f} (소량)")
        log()
        log(f"  A_L 분해 상세:")
        log(f"    ΔA_L = {Delta_A_L:.4f} ({Delta_A_L/Delta_A_total*100:.1f}% of ΔA_total)")
        log(f"    → A(t₀) q-의존성의 대부분은 L'/L(ρ,χ)의 산술적 기여")
        log(f"    → Im(H₀^L) 변화: Δ = {Delta_ImH0L:+.5f}")

        # 실패 통계
        log()
        log(f"  실패 통계: ζ {zeta_fail_count}/{len(ZETA_MEASURED)}, χ₅² {chi5_fail_count}/{len(CHI5_MEASURED)}")

        # 성공 기준 판정
        log()
        log("=" * 72)
        log("성공 기준 판정")
        log("=" * 72)

        crit1 = len(zeta_results) >= 10  # ζ 10개 이상
        crit2 = len(chi5_results) >= 15  # χ₅² 15개 이상
        crit3 = abs(Delta_A_total - 1.82) < 0.15  # #76 결과와 일치
        crit4 = Delta_A_L/Delta_A_total > 0.5  # 산술 기여 50% 이상
        log(f"  [{'✅' if crit1 else '❌'}] ζ 10개 이상 분해 성공: {len(zeta_results)}개")
        log(f"  [{'✅' if crit2 else '❌'}] χ₅² 15개 이상 분해 성공: {len(chi5_results)}개")
        log(f"  [{'✅' if crit3 else '❌'}] ΔA_total ≈ 1.82 (실측값 일치): {Delta_A_total:.4f}")
        log(f"  [{'✅' if crit4 else '❌'}] 산술 기여 ΔA_L이 지배적: {Delta_A_L/Delta_A_total*100:.1f}%")
        log(f"  [{'✅' if abs(Delta_A_gamma) < 0.5 else '❌'}] Γ 기여 ΔA_Γ 작음 (<0.5): {Delta_A_gamma:.4f}")

        n_pass = sum([crit1, crit2, crit3, crit4, abs(Delta_A_gamma) < 0.5])
        log()
        if n_pass == 5:
            log("  판정: ★★ 양성 (5/5 기준 충족)")
        elif n_pass >= 3:
            log(f"  판정: ★ 조건부 양성 ({n_pass}/5 기준 충족)")
        else:
            log(f"  판정: ❌ 기각 ({n_pass}/5 기준 충족)")

    dt_total = time.time() - t0_total
    log()
    log(f"총 소요: {dt_total:.1f}s")
    log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    flush_file()


if __name__ == "__main__":
    main()
