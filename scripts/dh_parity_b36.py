#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 C-242 / B-36] DH Laurent 패리티 — on-critical vs off-critical
=============================================================================
목적:
  패리티 정리 cₙ = (-1)^(n+1) c̄ₙ의 필요조건 검증:
  - σ=1/2 (on-critical) 영점: FE + SR → 1-ρ = ρ̄ → 패리티 PASS 예측
  - σ≠1/2 (off-critical) 영점: 1-ρ ≠ ρ̄ → 패리티 FAIL 예측

  DH는 FE + SR 모두 성립하나, off-critical 영점이 존재.
  → 패리티 깨짐은 SR 실패가 아니라 σ=1/2 이탈이 원인임을 실증.

이론:
  f(s) = α·L(s,χ₅) + ᾱ·L(s,χ̄₅)
  FE: Λ(s) = Λ(1-s)     [ε=1]
  SR: conj(f(s)) = f(s̄) [α·L̄(s,χ₅) = α·L(s̄,χ̄₅) = f(s̄)]

  On-critical: 1-ρ = ρ̄ → FE+SR combine → cₙ = (-1)^(n+1) c̄ₙ
  Off-critical: 1-ρ ≠ ρ̄ → FE and SR reference different points → no self-constraint

방법:
  Λ'/Λ(ρ+δ) = 1/δ + c₀ + c₁δ + c₂δ² + c₃δ³ + ...
  대칭/반대칭 분리 → 다항식 피팅 → c₀-c₃ 추출
  패리티 비율 측정: |Re(c₀)|/|Im(c₀)|, |Im(c₁)|/|Re(c₁)| 등

대상:
  On-critical: 5개 (σ=0.5, t=5.09..27.92) — 실험 #115
  Off-critical: 4개 (σ=0.57..0.81, t=85.7..176.7) — 실험 #115
결과: results/dh_parity_b36.txt
=============================================================================
"""
import sys, os, time
import numpy as np
import mpmath

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'dh_parity_b36.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 78)
log("[실험 C-242 / B-36] DH Laurent 패리티 — on-critical vs off-critical")
log("=" * 78)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()
log("핵심 가설:")
log("  패리티 cₙ = (-1)^(n+1) c̄ₙ 은 σ=1/2 (임계선)을 요구한다.")
log("  DH는 FE+SR 모두 성립하지만 off-critical 영점이 존재.")
log("  → on-critical: PASS, off-critical: FAIL 예측")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# DH 함수 구현 (mpmath, dps=120 for off-critical at t~177)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
mpmath.mp.dps = 120

CHI_MOD5 = [0, 1, 1j, -1j, -1]
CHI_MOD5_BAR = [0, 1, -1j, 1j, -1]

sqrt5 = mpmath.sqrt(5)
kappa_DH = (mpmath.sqrt(10 - 2*sqrt5) - 2) / (sqrt5 - 1)
coeff_chi     = (1 - 1j * kappa_DH) / 2
coeff_chi_bar = (1 + 1j * kappa_DH) / 2

log(f"κ_DH = {float(kappa_DH):.10f}")
log(f"coeff_chi = {complex(coeff_chi)}")
log(f"coeff_chi_bar = {complex(coeff_chi_bar)}")
log(f"mpmath.dps = {mpmath.mp.dps}")
log()

def dh_func(s):
    s = mpmath.mpc(s)
    return coeff_chi * mpmath.dirichlet(s, CHI_MOD5) + \
           coeff_chi_bar * mpmath.dirichlet(s, CHI_MOD5_BAR)

def Lambda_dh(s):
    """완비 DH: Λ(s) = (5/π)^{s/2} · Γ((s+1)/2) · f(s)"""
    s = mpmath.mpc(s)
    return mpmath.power(5/mpmath.pi, s/2) * mpmath.gamma((s+1)/2) * dh_func(s)

def log_deriv_Lambda(s):
    """Λ'/Λ(s) via central difference (h = 10^{-dps/3})"""
    s = mpmath.mpc(s)
    h = mpmath.mpf(10) ** (-(mpmath.mp.dps // 3))
    Lp = Lambda_dh(s + h)
    Lm = Lambda_dh(s - h)
    L0 = Lambda_dh(s)
    deriv = (Lp - Lm) / (2 * h)
    return deriv / L0

def refine_zero(sigma_guess, t_guess):
    """findroot로 영점 정밀화"""
    z0 = mpmath.mpc(sigma_guess, t_guess)
    try:
        root = mpmath.findroot(Lambda_dh, z0, tol=mpmath.mpf(10)**(-mpmath.mp.dps + 10))
        return float(mpmath.re(root)), float(mpmath.im(root))
    except:
        return sigma_guess, t_guess

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 영점 목록 (실험 #115에서)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ON_CRITICAL = [
    (0.5, 5.094160),
    (0.5, 8.939914),
    (0.5, 12.133545),
    (0.5, 14.404003),
    (0.5, 17.130239),
]

OFF_CRITICAL = [
    (0.808517, 85.699348),
    (0.650830, 114.163343),
    (0.574356, 166.479306),
    (0.724258, 176.702461),
]

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# δ 범위 (C-240과 동일)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DELTAS = np.array([0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.007,
                    0.01, 0.012, 0.015, 0.02, 0.025, 0.03, 0.035,
                    0.04, 0.045, 0.05, 0.06])


def extract_laurent_coefficients(sigma0, t0, deltas):
    """
    c₀, c₁, c₂, c₃ 동시 추출 — δ² 다항식 피팅.
    δ는 σ 방향 (실수축).

    Λ'/Λ(σ₀+δ+it₀) = 1/δ + c₀ + c₁δ + c₂δ² + ...

    대칭: sym(δ) = [f(+δ)+f(-δ)]/2 = c₀ + c₂δ² + ...
    반대칭: anti(δ) = [(f(+δ)-f(-δ))/2 - 1/δ]/δ = c₁ + c₃δ² + ...
    """
    sym_vals = []
    anti_vals = []
    good_deltas = []

    for d in deltas:
        try:
            s_plus = mpmath.mpc(sigma0 + d, t0)
            s_minus = mpmath.mpc(sigma0 - d, t0)

            f_plus = log_deriv_Lambda(s_plus)
            f_minus = log_deriv_Lambda(s_minus)

            f_p = complex(f_plus)
            f_m = complex(f_minus)

            sym = (f_p + f_m) / 2.0
            anti_raw = (f_p - f_m) / 2.0
            anti = (anti_raw - 1.0/d) / d

            sym_vals.append(sym)
            anti_vals.append(anti)
            good_deltas.append(d)
        except Exception as e:
            log(f"    ⚠️ δ={d}: {e}")

    if len(good_deltas) < 6:
        return None, None, None, None, False

    d_arr = np.array(good_deltas)
    d2 = d_arr ** 2

    # 대칭: sym = c₀ + c₂δ² + c₄δ⁴
    sym_re = np.array([s.real for s in sym_vals])
    sym_im = np.array([s.imag for s in sym_vals])

    p_sym_re = np.polyfit(d2, sym_re, 2)
    c0_re = p_sym_re[2]
    c2_re = p_sym_re[1]

    p_sym_im = np.polyfit(d2, sym_im, 2)
    c0_im = p_sym_im[2]
    c2_im = p_sym_im[1]

    c0 = complex(c0_re, c0_im)
    c2 = complex(c2_re, c2_im)

    # 반대칭: anti = c₁ + c₃δ²
    anti_re = np.array([a.real for a in anti_vals])
    anti_im = np.array([a.imag for a in anti_vals])

    p_anti_re = np.polyfit(d2, anti_re, 2)
    c1_re = p_anti_re[2]
    c3_re = p_anti_re[1]

    p_anti_im = np.polyfit(d2, anti_im, 2)
    c1_im = p_anti_im[2]
    c3_im = p_anti_im[1]

    c1 = complex(c1_re, c1_im)
    c3 = complex(c3_re, c3_im)

    return c0, c1, c2, c3, True


def check_parity(c, n, label):
    """
    패리티 검사: cₙ = (-1)^(n+1) c̄ₙ
    n 짝수 → c should be pure imaginary (Re=0)
    n 홀수 → c should be pure real (Im=0)

    반환: (ratio, pass_flag)
    ratio = |vanishing part| / |non-vanishing part|
    """
    if n % 2 == 0:
        # 짝수: Re(cₙ) = 0 예측
        vanish = abs(c.real)
        nonvanish = abs(c.imag)
        ratio = vanish / nonvanish if nonvanish > 1e-30 else float('inf')
        return ratio, f"|Re(c{n})|/|Im(c{n})| = {ratio:.3e}"
    else:
        # 홀수: Im(cₙ) = 0 예측
        vanish = abs(c.imag)
        nonvanish = abs(c.real)
        ratio = vanish / nonvanish if nonvanish > 1e-30 else float('inf')
        return ratio, f"|Im(c{n})|/|Re(c{n})| = {ratio:.3e}"


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SR 사전 검증 (DH에서 SR 성립 확인)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 0. SR 사전 검증: conj(Λ(s)) =? Λ(s̄) ━━━")
test_points = [mpmath.mpc(0.3, 10), mpmath.mpc(0.7, 20), mpmath.mpc(0.8, 50)]
for s_test in test_points:
    Ls = Lambda_dh(s_test)
    Lsbar = Lambda_dh(mpmath.conj(s_test))
    conjLs = mpmath.conj(Ls)
    rel = abs(Lsbar - conjLs) / abs(Ls) if abs(Ls) > 0 else 0
    log(f"  s={complex(s_test)}: |Λ(s̄)-conj(Λ(s))|/|Λ(s)| = {float(rel):.3e}  "
        f"{'✅' if float(rel) < 1e-20 else '❌'}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인: 패리티 테스트
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
all_results = []

for group_name, zeros in [("ON-CRITICAL (σ=0.5)", ON_CRITICAL),
                           ("OFF-CRITICAL (σ≠0.5)", OFF_CRITICAL)]:
    log(f"━━━ {group_name} ━━━")
    log()

    for idx, (sigma0, t0) in enumerate(zeros, 1):
        label = f"{'ON' if abs(sigma0 - 0.5) < 0.01 else 'OFF'}#{idx}"
        log(f"  [{label}] σ₀={sigma0:.6f}, t₀={t0:.6f}")

        # 영점 정밀화
        sigma_r, t_r = refine_zero(sigma0, t0)
        resid = abs(complex(Lambda_dh(mpmath.mpc(sigma_r, t_r))))
        log(f"    정밀화: σ={sigma_r:.8f}, t={t_r:.8f}, |Λ|={resid:.3e}")

        # Laurent 계수 추출
        c0, c1, c2, c3 = None, None, None, None
        try:
            c0, c1, c2, c3, ok = extract_laurent_coefficients(sigma_r, t_r, DELTAS)
        except Exception as e:
            log(f"    ❌ 추출 에러: {e}")
            ok = False

        if not ok:
            log(f"    ⚠️ 추출 실패")
            log()
            continue

        # 패리티 비율
        r0, s0 = check_parity(c0, 0, label)
        r1, s1 = check_parity(c1, 1, label)
        r2, s2 = check_parity(c2, 2, label)
        r3, s3 = check_parity(c3, 3, label)

        log(f"    c₀ = ({c0.real:+.6e}) + ({c0.imag:+.6e})i")
        log(f"    c₁ = ({c1.real:+.6e}) + ({c1.imag:+.6e})i")
        log(f"    c₂ = ({c2.real:+.6e}) + ({c2.imag:+.6e})i")
        log(f"    c₃ = ({c3.real:+.6e}) + ({c3.imag:+.6e})i")
        log(f"    패리티: {s0}")
        log(f"            {s1}")
        log(f"            {s2}")
        log(f"            {s3}")

        # on-critical: ratio < 1e-3 → PASS, off-critical: ratio > 0.01 → FAIL(= 패리티 깨짐)
        is_on = abs(sigma_r - 0.5) < 0.01
        c0_ok = r0 < 1e-3
        c1_ok = r1 < 1e-3
        c2_ok = r2 < 1e-2  # c₂ 정밀도 낮음
        c3_ok = r3 < 1e-2
        all_pass = c0_ok and c1_ok and c2_ok and c3_ok

        if is_on:
            verdict = "★ PARITY PASS" if all_pass else "⚠️ UNEXPECTED FAIL"
        else:
            verdict = "★ PARITY FAIL (as predicted)" if not all_pass else "⚠️ UNEXPECTED PASS"

        log(f"    판정: {verdict}")
        log()

        all_results.append({
            'label': label,
            'sigma': sigma_r,
            't': t_r,
            'is_on': is_on,
            'c0': c0, 'c1': c1, 'c2': c2, 'c3': c3,
            'r0': r0, 'r1': r1, 'r2': r2, 'r3': r3,
            'parity_pass': all_pass,
        })

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 통합 요약
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log()
log("=" * 78)
log("통합 요약")
log("=" * 78)
log()

if all_results:
    # 그룹별 테이블
    log(f"{'Label':<8} {'σ':>10} {'t':>10} {'|σ-½|':>8} "
        f"{'r(c₀)':>10} {'r(c₁)':>10} {'r(c₂)':>10} {'r(c₃)':>10} {'Parity':>12}")
    log("-" * 95)

    on_results = [r for r in all_results if r['is_on']]
    off_results = [r for r in all_results if not r['is_on']]

    for r in on_results:
        log(f"{r['label']:<8} {r['sigma']:>10.6f} {r['t']:>10.4f} {abs(r['sigma']-0.5):>8.4f} "
            f"{r['r0']:>10.2e} {r['r1']:>10.2e} {r['r2']:>10.2e} {r['r3']:>10.2e} "
            f"{'PASS' if r['parity_pass'] else 'FAIL':>12}")

    log("-" * 95)

    for r in off_results:
        log(f"{r['label']:<8} {r['sigma']:>10.6f} {r['t']:>10.4f} {abs(r['sigma']-0.5):>8.4f} "
            f"{r['r0']:>10.2e} {r['r1']:>10.2e} {r['r2']:>10.2e} {r['r3']:>10.2e} "
            f"{'PASS' if r['parity_pass'] else 'FAIL':>12}")

    log()

    # 통계
    on_pass = sum(1 for r in on_results if r['parity_pass'])
    off_fail = sum(1 for r in off_results if not r['parity_pass'])

    log(f"On-critical:  {on_pass}/{len(on_results)} PASS  (예측: 전부 PASS)")
    log(f"Off-critical: {off_fail}/{len(off_results)} FAIL  (예측: 전부 FAIL)")
    log()

    # on-critical 평균 비율 (small → 패리티 성립)
    if on_results:
        on_r0 = np.mean([r['r0'] for r in on_results])
        on_r1 = np.mean([r['r1'] for r in on_results])
        log(f"On-critical 평균 패리티 비율: r(c₀)={on_r0:.2e}, r(c₁)={on_r1:.2e}  (≈0 → 패리티)")

    # off-critical 평균 비율 (large → 패리티 깨짐)
    if off_results:
        off_r0 = np.mean([r['r0'] for r in off_results])
        off_r1 = np.mean([r['r1'] for r in off_results])
        log(f"Off-critical 평균 패리티 비율: r(c₀)={off_r0:.2e}, r(c₁)={off_r1:.2e}  (≫0 → 패리티 깨짐)")

    log()

    # 최종 판정
    if on_pass == len(on_results) and off_fail == len(off_results):
        log("★★★★ B-36 해결: 패리티 정리는 σ=1/2를 요구한다.")
        log()
        log("결론:")
        log("  1. DH에서 FE와 SR 모두 성립 — SR 실패는 경계가 아님")
        log("  2. On-critical 영점: 패리티 성립 (FE+SR+σ=½ → cₙ = (-1)^{n+1}c̄ₙ)")
        log("  3. Off-critical 영점: 패리티 깨짐 (1-ρ ≠ ρ̄ → 증명의 핵심 단계 불성립)")
        log("  4. ⟹ 패리티 ⟺ 임계선 위의 영점")
        log("  5. ⟹ GRH 위반 영점은 패리티 깨짐으로 탐지 가능")
    elif on_pass == len(on_results) and off_fail < len(off_results):
        log("⚠️ Off-critical 일부 PASS — 예상 외. 이론 재검토 필요")
    else:
        log("⚠️ On-critical 일부 FAIL — 수치 정밀도 점검 필요")

log()
log(f"총 소요: {time.time()-START:.1f}s")
log("=" * 78)
outf.close()
