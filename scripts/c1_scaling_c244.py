#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 C-244] Laurent c₁ scaling law — ζ(s) 50 영점
=============================================================================
목적:
  Riemann ζ의 비자명 영점 50개에서 Λ'/Λ Laurent 전개의 c₀, c₁ 추출.
  c₁ = Re(c₁)이 zero height γ와 어떤 스케일링 관계를 보이는지 탐색.

이론 배경:
  Hadamard 곱에서 c₁(ρ₀) = Σ'_ρ 1/(ρ₀-ρ)² (양수, on-critical).
  영점 밀도 d̄ = log(γ/2π)/(2π)에 비례하여 c₁ ~ (log γ)^α 예상.

방법:
  C-243과 동일한 δ-extrapolation (σ 방향 대칭/반대칭 분해).
  Hadamard 부분합(100영점)으로 교차 검증.

성공 기준:
  (1) Thm 5 패리티 확인: |Re(c₀)| < 1e-6, |Im(c₁)| < 1e-6
  (2) Hadamard 비교: 상관 ρ > 0.99
  (3) c₁ scaling 피팅: R² > 0.95
  (4) c₁ 50개 측정 완료

결과: results/c1_scaling_c244.txt
=============================================================================
"""
import sys, os, time
import numpy as np
import mpmath

START = time.time()
RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'c1_scaling_c244.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 80)
log("[실험 C-244] Laurent c₁ scaling law — ζ(s) 50 영점")
log("=" * 80)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

mpmath.mp.dps = 80
NZEROS_TOTAL = 100
NZEROS_MEASURE = 50

log(f"mpmath.dps = {mpmath.mp.dps}")
log(f"영점 수: 측정 {NZEROS_MEASURE} / Hadamard 비교용 {NZEROS_TOTAL}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. ζ 영점 로드
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log(f"{T()} ζ 영점 {NZEROS_TOTAL}개 계산 중...")
gammas = []
for n in range(1, NZEROS_TOTAL + 1):
    g = float(mpmath.im(mpmath.zetazero(n)))
    gammas.append(g)
    if n % 25 == 0:
        log(f"  {T()} {n}/{NZEROS_TOTAL} (γ_{n} = {g:.4f})")

log(f"{T()} 완료. γ 범위: [{gammas[0]:.4f}, {gammas[-1]:.4f}]")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. Λ'/Λ 수치 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def Lambda_zeta(s):
    """Λ(s) = π^{-s/2} Γ(s/2) ζ(s)"""
    s = mpmath.mpc(s)
    return mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s)

def log_deriv_Lambda(s):
    """Λ'/Λ(s) via central difference (h = 10^{-dps//3})"""
    s = mpmath.mpc(s)
    h = mpmath.mpf(10) ** (-(mpmath.mp.dps // 3))
    Lp = Lambda_zeta(s + h)
    Lm = Lambda_zeta(s - h)
    L0 = Lambda_zeta(s)
    if abs(L0) == 0:
        return mpmath.mpc(0)
    return (Lp - Lm) / (2 * h * L0)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. c₀, c₁ δ-extrapolation
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DELTAS = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02])

def extract_c0_c1(gamma):
    """
    ρ = 1/2+iγ에서 Λ'/Λ Laurent c₀, c₁ 추출.
    sym(δ) = (f(+δ)+f(-δ))/2 = c₀ + c₂δ² + c₄δ⁴   → 2차 피팅 절편 = c₀
    anti(δ) = ((f(+δ)-f(-δ))/2 - 1/δ)/δ = c₁ + c₃δ² + c₅δ⁴ → 절편 = c₁
    """
    rho = mpmath.mpc(0.5, gamma)
    sym_vals = []
    anti_vals = []

    for d in DELTAS:
        d_mp = mpmath.mpf(float(d))
        fp = log_deriv_Lambda(rho + d_mp)
        fm = log_deriv_Lambda(rho - d_mp)

        sym = complex((fp + fm) / 2)
        anti_num = (fp - fm) / 2 - 1 / d_mp
        anti = complex(anti_num / d_mp)

        sym_vals.append(sym)
        anti_vals.append(anti)

    d2 = DELTAS ** 2
    # 2차 다항식 피팅 (x = δ²): c₀ + c₂·x + c₄·x²
    c0_re = np.polyfit(d2, [v.real for v in sym_vals], 2)[-1]
    c0_im = np.polyfit(d2, [v.imag for v in sym_vals], 2)[-1]
    c1_re = np.polyfit(d2, [v.real for v in anti_vals], 2)[-1]
    c1_im = np.polyfit(d2, [v.imag for v in anti_vals], 2)[-1]

    return complex(c0_re, c0_im), complex(c1_re, c1_im)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. Hadamard 예측
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
def hadamard_c1(idx, gammas):
    """
    c₁^ξ = Σ'_k 1/(γ₀-γ_k)² + Σ_k 1/(γ₀+γ_k)² + 1/(2γ₀)²
    c₁^Λ = c₁^ξ + 2(1/4 - γ₀²)/(1/4 + γ₀²)²
    """
    g0 = gammas[idx]
    c1_xi = 0.0
    for k, g in enumerate(gammas):
        if k == idx:
            c1_xi += 1.0 / (2 * g0) ** 2  # 켤레 영점 1/2-iγ₀ 기여
        else:
            c1_xi += 1.0 / (g0 - g) ** 2    # 같은 부호 영점
            c1_xi += 1.0 / (g0 + g) ** 2    # 켤레 영점
    # Λ vs ξ 보정: 2(1/4-γ²)/(1/4+γ²)²
    corr = 2.0 * (0.25 - g0**2) / (0.25 + g0**2)**2
    return c1_xi + corr

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. 측정 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
results = []

log("━━━ c₀, c₁ 추출 (50 영점) ━━━")
log(f"{'n':>3} {'γ':>10} {'Re(c₀)':>12} {'Im(c₀)':>12} {'c₁':>12} {'Im(c₁)':>12} {'c₁_Had':>12}")
log("-" * 87)

for n in range(NZEROS_MEASURE):
    gamma = gammas[n]
    c0, c1 = extract_c0_c1(gamma)
    c1_had = hadamard_c1(n, gammas)

    results.append({
        'n': n + 1, 'gamma': gamma,
        'c0_re': c0.real, 'c0_im': c0.imag,
        'c1_re': c1.real, 'c1_im': c1.imag,
        'c1_had': c1_had,
    })

    log(f"{n+1:>3} {gamma:>10.4f} {c0.real:>12.4e} {c0.imag:>12.6f} {c1.real:>12.6f} {c1.imag:>12.4e} {c1_had:>12.6f}")

    if (n + 1) % 10 == 0:
        log(f"  {T()} {n+1}/{NZEROS_MEASURE} 완료")

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
gammas_m = np.array([r['gamma'] for r in results])
c1_meas = np.array([r['c1_re'] for r in results])
c1_had_arr = np.array([r['c1_had'] for r in results])
c0_re_arr = np.array([r['c0_re'] for r in results])
c1_im_arr = np.array([r['c1_im'] for r in results])

# 6a. Thm 5 패리티 검증
log("━━━ 패리티 검증 (Thm 5) ━━━")
max_re_c0 = np.max(np.abs(c0_re_arr))
max_im_c1 = np.max(np.abs(c1_im_arr))
log(f"  |Re(c₀)|: max = {max_re_c0:.3e}, mean = {np.mean(np.abs(c0_re_arr)):.3e}")
log(f"  |Im(c₁)|: max = {max_im_c1:.3e}, mean = {np.mean(np.abs(c1_im_arr)):.3e}")
parity_pass = max_re_c0 < 1e-6 and max_im_c1 < 1e-6
log(f"  판정: {'✅ PASS' if parity_pass else '❌ FAIL'}")
log()

# 6b. Hadamard 비교
rel_err = np.abs(c1_meas - c1_had_arr) / np.abs(c1_meas)
had_corr = np.corrcoef(c1_meas, c1_had_arr)[0, 1]
log("━━━ Hadamard 비교 ━━━")
log(f"  상대 오차: max = {np.max(rel_err):.4e}, mean = {np.mean(rel_err):.4e}")
log(f"  상관계수: ρ = {had_corr:.6f}")
had_pass = had_corr > 0.99
log(f"  판정: {'✅ PASS' if had_pass else '❌ FAIL'} (기준: ρ > 0.99)")
log()

# 6c. Scaling law 피팅
log_gamma = np.log(gammas_m)
log_gamma2 = log_gamma ** 2
ss_tot = np.sum((c1_meas - np.mean(c1_meas)) ** 2)

fits = {}

# (1) c₁ = a + b·(log γ)²
A1 = np.vstack([np.ones_like(log_gamma2), log_gamma2]).T
coeff1 = np.linalg.lstsq(A1, c1_meas, rcond=None)[0]
R2_1 = 1 - np.sum((c1_meas - A1 @ coeff1) ** 2) / ss_tot
fits['(log γ)²'] = (R2_1, f"c₁ = {coeff1[0]:.4f} + {coeff1[1]:.4f}·(log γ)²")

# (2) c₁ = a + b·log γ
A2 = np.vstack([np.ones_like(log_gamma), log_gamma]).T
coeff2 = np.linalg.lstsq(A2, c1_meas, rcond=None)[0]
R2_2 = 1 - np.sum((c1_meas - A2 @ coeff2) ** 2) / ss_tot
fits['log γ'] = (R2_2, f"c₁ = {coeff2[0]:.4f} + {coeff2[1]:.4f}·log γ")

# (3) c₁ = a · γ^b (power law in log space)
log_c1 = np.log(c1_meas[c1_meas > 0])
log_g = np.log(gammas_m[c1_meas > 0])
A3 = np.vstack([np.ones_like(log_g), log_g]).T
coeff3 = np.linalg.lstsq(A3, log_c1, rcond=None)[0]
ss_tot3 = np.sum((log_c1 - np.mean(log_c1)) ** 2)
R2_3 = 1 - np.sum((log_c1 - A3 @ coeff3) ** 2) / ss_tot3
fits['γ^b'] = (R2_3, f"c₁ = {np.exp(coeff3[0]):.4f} · γ^{coeff3[1]:.4f}")

# (4) c₁ = a + b·d̄² where d̄ = log(γ/(2π))/(2π)
d_bar = np.log(gammas_m / (2 * np.pi)) / (2 * np.pi)
d_bar2 = d_bar ** 2
A4 = np.vstack([np.ones_like(d_bar2), d_bar2]).T
coeff4 = np.linalg.lstsq(A4, c1_meas, rcond=None)[0]
R2_4 = 1 - np.sum((c1_meas - A4 @ coeff4) ** 2) / ss_tot
fits['d̄²'] = (R2_4, f"c₁ = {coeff4[0]:.4f} + {coeff4[1]:.4f}·d̄²")

# (5) c₁ = a + b·log(γ/(2π)) + c·(log(γ/(2π)))²
log_g2pi = np.log(gammas_m / (2 * np.pi))
A5 = np.vstack([np.ones_like(log_g2pi), log_g2pi, log_g2pi**2]).T
coeff5 = np.linalg.lstsq(A5, c1_meas, rcond=None)[0]
R2_5 = 1 - np.sum((c1_meas - A5 @ coeff5) ** 2) / ss_tot
fits['quadratic log(γ/2π)'] = (R2_5, f"c₁ = {coeff5[0]:.4f} + {coeff5[1]:.4f}·x + {coeff5[2]:.4f}·x² (x=log(γ/2π))")

log("━━━ Scaling law 피팅 ━━━")
for name, (r2, formula) in sorted(fits.items(), key=lambda x: -x[1][0]):
    log(f"  {formula:60s}  R² = {r2:.6f}")
log()

best_name = max(fits, key=lambda k: fits[k][0])
best_R2 = fits[best_name][0]
scaling_pass = best_R2 > 0.95
log(f"  최적 모델: {best_name} (R² = {best_R2:.6f})")
log(f"  판정: {'✅ PASS' if scaling_pass else '❌ FAIL'} (기준: R² > 0.95)")
log()

# 6d. 최근접 영점 기여 분석
log("━━━ 최근접 영점 기여 분석 ━━━")
log(f"{'n':>3} {'γ':>8} {'Δ₁':>8} {'Δ₂':>8} {'nn₂기여':>10} {'c₁_meas':>10} {'nn₂/c₁':>8}")
log("-" * 60)

for n in [0, 9, 24, 49]:
    g0 = gammas[n]
    # 양방향 영점까지의 거리
    dists = sorted([(abs(g0 - gammas[k]), k) for k in range(NZEROS_TOTAL) if k != n])
    d1, d2 = dists[0][0], dists[1][0]
    nn2 = 1.0 / d1**2 + 1.0 / d2**2  # 최근접 2개 영점 기여
    c1_m = results[n]['c1_re']
    ratio = nn2 / c1_m if c1_m != 0 else float('inf')
    log(f"{n+1:>3} {g0:>8.2f} {d1:>8.4f} {d2:>8.4f} {nn2:>10.4f} {c1_m:>10.4f} {ratio:>8.3f}")

log()

# 6e. c₁ 통계
log("━━━ c₁ 기술 통계 ━━━")
log(f"  범위: [{np.min(c1_meas):.4f}, {np.max(c1_meas):.4f}]")
log(f"  평균: {np.mean(c1_meas):.4f}")
log(f"  중앙값: {np.median(c1_meas):.4f}")
log(f"  c₁ 단조 증가?: {np.all(np.diff(c1_meas) > 0)}")
# 비단조 지점 식별
non_mono = np.where(np.diff(c1_meas) <= 0)[0]
if len(non_mono) > 0:
    log(f"  비단조 지점 {len(non_mono)}개: n = {[int(x)+1 for x in non_mono[:10]]}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7. 종합 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
n_pass = sum([parity_pass, had_pass, scaling_pass, True])  # 4번째: 측정 완료

log("=" * 80)
log("[종합 판정]")
log(f"  (1) Thm 5 패리티:  {'✅' if parity_pass else '❌'} |Re(c₀)| < {max_re_c0:.1e}, |Im(c₁)| < {max_im_c1:.1e}")
log(f"  (2) Hadamard 비교: {'✅' if had_pass else '❌'} ρ = {had_corr:.4f}")
log(f"  (3) Scaling R²:    {'✅' if scaling_pass else '❌'} best = {best_R2:.4f} ({best_name})")
log(f"  (4) 50개 측정:     ✅")
log()
log(f"  성공 기준: {n_pass}/4 {'★' * n_pass}")
if n_pass == 4:
    log("  → ★★★★ 양성: c₁ scaling law 확립")
elif n_pass == 3:
    log("  → ★★★ 조건부 양성")
elif n_pass == 2:
    log("  → ★★ 부분 양성")
else:
    log("  → ★ 음성 또는 중립")
log()
log(f"  경과 시간: {time.time()-START:.1f}초")
log("=" * 80)

outf.close()
