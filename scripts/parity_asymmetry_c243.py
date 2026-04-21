#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 C-243 / B-41] Laurent 패리티 짝수/홀수 비대칭 정밀 탐사
=============================================================================
목적:
  C-242에서 발견된 미설명 현상: off-critical 영점에서
  - 짝수차(c₀,c₂) 패리티 위반이 크고 (r ~ 2-3000)
  - 홀수차(c₁,c₃) 패리티 위반이 작다 (r ~ 10⁻⁶-0.03)

  이 비대칭을 c₀-c₇ (8차)까지 확장하여:
  (a) 홀수차 위반이 n→∞에서 지수적으로 감쇄하는지 확인
  (b) r vs |σ-½| power law 관계 추출
  (c) 이것이 FE 단독 잔류 대칭인지 검증

방법:
  Λ'/Λ(ρ+δ) = 1/δ + Σ cₙ δⁿ
  대칭:   sym(δ)  = [f(+δ)+f(-δ)]/2    = c₀ + c₂δ² + c₄δ⁴ + c₆δ⁶
  반대칭: anti(δ) = [(f(+δ)-f(-δ))/2 - 1/δ]/δ = c₁ + c₃δ² + c₅δ⁴ + c₇δ⁶
  → 각각 δ² 기준 3차 다항식 피팅으로 c₀-c₇ 추출

대상:
  On-critical:  5개 (σ=0.5) — 수치 노이즈 바닥 기준
  Off-critical: 4개 (σ≠0.5) — 비대칭 측정

결과: results/parity_asymmetry_c243.txt
=============================================================================
"""
import sys, os, time
import numpy as np
import mpmath

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'parity_asymmetry_c243.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 80)
log("[실험 C-243 / B-41] Laurent 패리티 짝수/홀수 비대칭 정밀 탐사")
log("=" * 80)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# DH 함수 구현 (C-242와 동일, dps=120)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
mpmath.mp.dps = 120

CHI_MOD5 = [0, 1, 1j, -1j, -1]
CHI_MOD5_BAR = [0, 1, -1j, 1j, -1]

sqrt5 = mpmath.sqrt(5)
kappa_DH = (mpmath.sqrt(10 - 2*sqrt5) - 2) / (sqrt5 - 1)
coeff_chi     = (1 - 1j * kappa_DH) / 2
coeff_chi_bar = (1 + 1j * kappa_DH) / 2

log(f"κ_DH = {float(kappa_DH):.10f}")
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
    if abs(L0) == 0:
        return mpmath.mpc(0, 0)
    deriv = (Lp - Lm) / (2 * h)
    return deriv / L0

def refine_zero(sigma_guess, t_guess):
    """findroot로 영점 정밀화"""
    z0 = mpmath.mpc(sigma_guess, t_guess)
    try:
        root = mpmath.findroot(Lambda_dh, z0, tol=mpmath.mpf(10)**(-mpmath.mp.dps + 10))
        return float(mpmath.re(root)), float(mpmath.im(root))
    except Exception as e:
        log(f"    ⚠️ findroot 실패: {e}, 원래 값 사용")
        return sigma_guess, t_guess

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 영점 목록 (C-242와 동일)
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
# δ 범위 — c₆,c₇까지 안정 추출하려면 더 많은 점 필요
# 고차 피팅(δ² 기준 3차)에는 최소 7점 이상, 넉넉히 24점 사용
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DELTAS = np.array([
    0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.004, 0.005,
    0.006, 0.007, 0.008, 0.01, 0.012, 0.015, 0.018, 0.02,
    0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06
])


def extract_laurent_c0_to_c7(sigma0, t0, deltas):
    """
    c₀~c₇ 동시 추출.

    대칭:   sym(δ)  = c₀ + c₂·δ² + c₄·δ⁴ + c₆·δ⁶
    반대칭: anti(δ) = c₁ + c₃·δ² + c₅·δ⁴ + c₇·δ⁶

    각각 x=δ² 기준 3차 다항식으로 피팅:
      sym  = c₆·x³ + c₄·x² + c₂·x + c₀
      anti = c₇·x³ + c₅·x² + c₃·x + c₁
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

            # 대칭 부분: (f(+δ)+f(-δ))/2 = c₀ + c₂δ² + c₄δ⁴ + c₆δ⁶
            sym = (f_p + f_m) / 2.0

            # 반대칭 부분: [(f(+δ)-f(-δ))/2 - 1/δ] / δ = c₁ + c₃δ² + c₅δ⁴ + c₇δ⁶
            anti_raw = (f_p - f_m) / 2.0
            anti = (anti_raw - 1.0/d) / d

            # NaN/Inf 체크
            if np.isnan(sym) or np.isinf(sym) or np.isnan(anti) or np.isinf(anti):
                log(f"    ⚠️ δ={d}: NaN/Inf 검출, 스킵")
                continue

            sym_vals.append(sym)
            anti_vals.append(anti)
            good_deltas.append(d)
        except Exception as e:
            log(f"    ⚠️ δ={d}: {e}")

    if len(good_deltas) < 8:
        log(f"    ⚠️ 유효 δ {len(good_deltas)}개 < 8 — 추출 불가")
        return None, False

    d_arr = np.array(good_deltas)
    x = d_arr ** 2  # x = δ²

    coeffs = {}

    # 대칭 피팅: sym = c₆·x³ + c₄·x² + c₂·x + c₀
    sym_re = np.array([s.real for s in sym_vals])
    sym_im = np.array([s.imag for s in sym_vals])

    p_sym_re = np.polyfit(x, sym_re, 3)  # [c₆_re, c₄_re, c₂_re, c₀_re]
    p_sym_im = np.polyfit(x, sym_im, 3)

    coeffs[0] = complex(p_sym_re[3], p_sym_im[3])  # c₀
    coeffs[2] = complex(p_sym_re[2], p_sym_im[2])  # c₂
    coeffs[4] = complex(p_sym_re[1], p_sym_im[1])  # c₄
    coeffs[6] = complex(p_sym_re[0], p_sym_im[0])  # c₆

    # 반대칭 피팅: anti = c₇·x³ + c₅·x² + c₃·x + c₁
    anti_re = np.array([a.real for a in anti_vals])
    anti_im = np.array([a.imag for a in anti_vals])

    p_anti_re = np.polyfit(x, anti_re, 3)
    p_anti_im = np.polyfit(x, anti_im, 3)

    coeffs[1] = complex(p_anti_re[3], p_anti_im[3])  # c₁
    coeffs[3] = complex(p_anti_re[2], p_anti_im[2])  # c₃
    coeffs[5] = complex(p_anti_re[1], p_anti_im[1])  # c₅
    coeffs[7] = complex(p_anti_re[0], p_anti_im[0])  # c₇

    return coeffs, True


def parity_ratio(c, n):
    """
    패리티 비율 r(cₙ).
    n 짝수 → Re(cₙ)=0 예측 → r = |Re|/|Im|
    n 홀수 → Im(cₙ)=0 예측 → r = |Im|/|Re|
    """
    if n % 2 == 0:
        vanish = abs(c.real)
        nonvanish = abs(c.imag)
    else:
        vanish = abs(c.imag)
        nonvanish = abs(c.real)

    if nonvanish < 1e-50:
        return float('inf')
    return vanish / nonvanish


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실험
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
all_results = []

for group_name, zeros in [("ON-CRITICAL (σ=0.5)", ON_CRITICAL),
                           ("OFF-CRITICAL (σ≠0.5)", OFF_CRITICAL)]:
    log(f"━━━ {group_name} ━━━")
    log()

    for idx, (sigma0, t0) in enumerate(zeros, 1):
        is_on = abs(sigma0 - 0.5) < 0.01
        label = f"{'ON' if is_on else 'OFF'}#{idx}"
        log(f"  [{label}] σ₀={sigma0:.6f}, t₀={t0:.6f}")

        # 영점 정밀화
        sigma_r, t_r = refine_zero(sigma0, t0)
        resid = abs(complex(Lambda_dh(mpmath.mpc(sigma_r, t_r))))
        log(f"    정밀화: σ={sigma_r:.8f}, t={t_r:.8f}, |Λ|={resid:.3e}")

        # Laurent 계수 c₀-c₇ 추출
        coeffs, ok = extract_laurent_c0_to_c7(sigma_r, t_r, DELTAS)

        if not ok or coeffs is None:
            log(f"    ❌ 추출 실패")
            log()
            continue

        # 각 계수 출력
        for n in range(8):
            c = coeffs[n]
            r = parity_ratio(c, n)
            parity_type = "Re=0" if n % 2 == 0 else "Im=0"
            log(f"    c{n} = ({c.real:+.8e}) + ({c.imag:+.8e})i  "
                f"  r(c{n}) = {r:.3e}  [{parity_type}]")

        log()

        all_results.append({
            'label': label,
            'sigma': sigma_r,
            't': t_r,
            'is_on': is_on,
            'sigma_dev': abs(sigma_r - 0.5),
            'coeffs': coeffs,
            'ratios': {n: parity_ratio(coeffs[n], n) for n in range(8)},
        })

    log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 분석 (a): n vs log₁₀(r(cₙ)) 테이블
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("=" * 80)
log("분석 (a): n vs log₁₀(r(cₙ)) — 짝수/홀수 분리")
log("=" * 80)
log()

on_results = [r for r in all_results if r['is_on']]
off_results = [r for r in all_results if not r['is_on']]

# 테이블 헤더
header = f"{'n':>3} {'parity':>8}"
for r in on_results:
    header += f"  {r['label']:>10}"
for r in off_results:
    header += f"  {r['label']:>10}"
log(header)
log("-" * len(header))

for n in range(8):
    row = f"{n:>3} {'even' if n%2==0 else 'odd':>8}"
    for r in on_results:
        val = r['ratios'][n]
        if val == float('inf'):
            row += f"  {'inf':>10}"
        elif val == 0:
            row += f"  {'0':>10}"
        else:
            row += f"  {np.log10(val):>10.2f}"
    for r in off_results:
        val = r['ratios'][n]
        if val == float('inf'):
            row += f"  {'inf':>10}"
        elif val == 0:
            row += f"  {'0':>10}"
        else:
            row += f"  {np.log10(val):>10.2f}"
    log(row)

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 분석 (b): off-critical 평균 r(cₙ) vs n — 감쇄율 α
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("=" * 80)
log("분석 (b): off-critical 평균 r(cₙ) vs n — 감쇄율 추출")
log("=" * 80)
log()

if off_results:
    mean_r_even = {}
    mean_r_odd = {}

    for n in range(8):
        vals = [r['ratios'][n] for r in off_results if r['ratios'][n] != float('inf')]
        if vals:
            mean_val = np.mean(vals)
            if n % 2 == 0:
                mean_r_even[n] = mean_val
            else:
                mean_r_odd[n] = mean_val

    log("짝수차 (even n):")
    for n in sorted(mean_r_even.keys()):
        log(f"  n={n}: mean r = {mean_r_even[n]:.4e}  (log₁₀ = {np.log10(mean_r_even[n]):.2f})")

    log()
    log("홀수차 (odd n):")
    for n in sorted(mean_r_odd.keys()):
        log(f"  n={n}: mean r = {mean_r_odd[n]:.4e}  (log₁₀ = {np.log10(mean_r_odd[n]):.2f})")

    log()

    # 짝수/홀수 비대칭 비율
    if mean_r_even and mean_r_odd:
        overall_even = np.mean(list(mean_r_even.values()))
        overall_odd = np.mean(list(mean_r_odd.values()))
        log(f"전체 평균 — 짝수: {overall_even:.4e}, 홀수: {overall_odd:.4e}")
        if overall_odd > 0:
            log(f"비대칭 비율 (even/odd): {overall_even/overall_odd:.1f}×")
        log()

    # 홀수차 감쇄율 추출: log(r_odd) ~ -α·n + β
    odd_ns = sorted(mean_r_odd.keys())
    if len(odd_ns) >= 3:
        odd_n_arr = np.array(odd_ns, dtype=float)
        odd_logr = np.array([np.log10(mean_r_odd[n]) for n in odd_ns])

        # 선형 피팅: log₁₀(r) = slope·n + intercept
        slope, intercept = np.polyfit(odd_n_arr, odd_logr, 1)

        log(f"홀수차 감쇄 피팅: log₁₀(r_odd) = {slope:.4f}·n + {intercept:.4f}")
        log(f"  → 감쇄율 (per 2 orders): Δlog₁₀(r) = {2*slope:.4f}")
        log(f"  → 지수적 감쇄? slope < 0 → {'YES (α={:.4f})'.format(-slope) if slope < 0 else 'NO'}")
        log()

    # 짝수차 감쇄율도 참고
    even_ns = sorted(mean_r_even.keys())
    if len(even_ns) >= 3:
        even_n_arr = np.array(even_ns, dtype=float)
        even_logr = np.array([np.log10(mean_r_even[n]) for n in even_ns])

        slope_e, intercept_e = np.polyfit(even_n_arr, even_logr, 1)
        log(f"짝수차 감쇄 피팅 (참고): log₁₀(r_even) = {slope_e:.4f}·n + {intercept_e:.4f}")
        log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 분석 (c): off-critical r(c₀) vs |σ-½| — power law β
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("=" * 80)
log("분석 (c): off-critical r(cₙ) vs |σ-½| — power law 추출")
log("=" * 80)
log()

if off_results:
    # 각 영점의 |σ-½|와 r(c₀), r(c₁) 등
    log(f"{'Label':<8} {'|σ-½|':>8} {'r(c₀)':>12} {'r(c₁)':>12} {'r(c₂)':>12} {'r(c₃)':>12}")
    log("-" * 65)
    for r in off_results:
        log(f"{r['label']:<8} {r['sigma_dev']:>8.4f} "
            f"{r['ratios'][0]:>12.4e} {r['ratios'][1]:>12.4e} "
            f"{r['ratios'][2]:>12.4e} {r['ratios'][3]:>12.4e}")
    log()

    # Power law 피팅: log₁₀(r) = β·log₁₀(|σ-½|) + const
    for n_target in [0, 1, 2]:
        devs = np.array([r['sigma_dev'] for r in off_results])
        r_vals = np.array([r['ratios'][n_target] for r in off_results])

        # inf 제거
        mask = np.isfinite(r_vals) & (r_vals > 0) & (devs > 0)
        if np.sum(mask) >= 3:
            log_dev = np.log10(devs[mask])
            log_r = np.log10(r_vals[mask])

            beta, const = np.polyfit(log_dev, log_r, 1)
            log(f"  r(c{n_target}) vs |σ-½|: β = {beta:.3f}  (r ~ |σ-½|^{beta:.3f})")
        else:
            log(f"  r(c{n_target}): 데이터 부족 ({np.sum(mask)}점)")

    log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 분석 (d): on-critical r(cₙ) vs n — 수치 노이즈 바닥
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("=" * 80)
log("분석 (d): on-critical 수치 노이즈 바닥")
log("=" * 80)
log()

if on_results:
    for n in range(8):
        vals = [r['ratios'][n] for r in on_results if np.isfinite(r['ratios'][n])]
        if vals:
            mean_v = np.mean(vals)
            max_v = np.max(vals)
            log(f"  n={n} ({'even' if n%2==0 else 'odd '}): "
                f"mean r = {mean_v:.3e}, max r = {max_v:.3e}  "
                f"{'✅ noise floor' if mean_v < 1e-3 else '⚠️ elevated'}")
        else:
            log(f"  n={n}: 데이터 없음")

    log()

    # 고차(n≥4)에서 노이즈 상승 확인
    on_mean_by_n = {}
    for n in range(8):
        vals = [r['ratios'][n] for r in on_results if np.isfinite(r['ratios'][n])]
        if vals:
            on_mean_by_n[n] = np.mean(vals)

    if on_mean_by_n:
        low_noise = np.mean([on_mean_by_n[n] for n in [0,1,2,3] if n in on_mean_by_n])
        high_noise = np.mean([on_mean_by_n[n] for n in [4,5,6,7] if n in on_mean_by_n])
        log(f"  노이즈 비교: n≤3 평균 = {low_noise:.3e}, n≥4 평균 = {high_noise:.3e}")
        if high_noise > 0 and low_noise > 0:
            log(f"  고차 노이즈 증폭: {high_noise/low_noise:.1f}×")
    log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 종합 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("=" * 80)
log("종합 판정")
log("=" * 80)
log()

success_count = 0
total_criteria = 4

# 기준 1: 72개 값 계산 성공
computed = sum(len(r['ratios']) for r in all_results)
total_expected = 9 * 8  # 9영점 × 8차
log(f"기준 1: {computed}/{total_expected} 값 계산 {'✅' if computed == total_expected else '⚠️'}")
if computed == total_expected:
    success_count += 1

# 기준 2: 짝수/홀수 비대칭 > 100배
if off_results and mean_r_even and mean_r_odd:
    ratio_asym = overall_even / overall_odd if overall_odd > 0 else float('inf')
    log(f"기준 2: 짝수/홀수 비대칭 = {ratio_asym:.1f}× {'✅ > 100' if ratio_asym > 100 else '❌ ≤ 100'}")
    if ratio_asym > 100:
        success_count += 1
else:
    log("기준 2: 데이터 부족 ❌")

# 기준 3: 감쇄율 α > 0
if off_results and len(odd_ns) >= 3:
    alpha = -slope  # slope은 위에서 계산된 값
    log(f"기준 3: 감쇄율 α = {alpha:.4f} {'✅ 양성 (감쇄 확인)' if alpha > 0 else '❌ 비양성'}")
    if alpha > 0:
        success_count += 1
else:
    log("기준 3: 데이터 부족 ❌")

# 기준 4: on-critical 노이즈 바닥이 off-critical보다 6자릿수 아래
if on_results and off_results:
    on_max = max(r['ratios'][0] for r in on_results if np.isfinite(r['ratios'][0]))
    off_min = min(r['ratios'][0] for r in off_results if np.isfinite(r['ratios'][0]))
    if on_max > 0 and off_min > 0:
        sep = np.log10(off_min / on_max)
        log(f"기준 4: on/off 분리도 = {sep:.1f} dex {'✅ ≥ 6' if sep >= 6 else '⚠️ < 6'}")
        if sep >= 6:
            success_count += 1
    else:
        log("기준 4: 값 부족 ❌")

log()
log(f"성공 기준: {success_count}/{total_criteria}")

if success_count == total_criteria:
    log()
    log("★★★★ C-243 양성: 짝수/홀수 패리티 비대칭 확인")
    log("  - 홀수차 위반은 짝수차보다 수 자릿수 약하며 고차에서 감쇄")
    log("  - FE 단독의 잔류 대칭이 홀수차 패리티를 근사적으로 보존")
    log("  → Paper 2에 Observation/Remark 추가 가능")
elif success_count >= 2:
    log()
    log("★★★ C-243 부분 양성: 비대칭은 관측되나 일부 기준 미달")
    log("  → 추가 분석 또는 파라미터 조정 필요")
else:
    log()
    log("★★ C-243 중립/음성: 비대칭이 수치 노이즈 수준이거나 일관성 부족")
    log("  → 수학자에게 보고: 4개 off-critical 영점으로는 통계 불충분 가능성")

log()
log(f"총 소요: {time.time()-START:.1f}s")
log("=" * 80)
outf.close()
