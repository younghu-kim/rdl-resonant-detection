#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 실험 #117 — c₁ 보편 법칙 정밀 검증 + 해석적 유도 (v2)
=============================================================================
목표:
1. 기존 4개 off-critical 영점 (OFF#1-4)에서 c₁ 완전 분석 (Phase 1 — 빠름)
   a. δ-sweep → c₁_fit (절편 회귀)
   b. c₁_analytic = Re(Λ''(ρ)/Λ'(ρ)) 수치 2차 미분
   c. 거울쌍 기여: 1/(2|σ-1/2|) vs 배경
   d. 비교표 + c₁·|σ-1/2| 통계

2. 추가 off-critical 영점 탐색 t∈[200,400] (Phase 2 — 옵션)
   - dps=80, T_STEP=2.0, 시간 예산 600초

핵심 예측: c₁·|σ₀-1/2| → 1 (c₁ ≈ 1/|σ₀-1/2|)
이론: c₁ = Re(Λ''(ρ)/Λ'(ρ)) — FE 대칭이 on-critical에서 이를 0으로 만듦

결과: results/dh_c1_law_117.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

START = time.time()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 출력 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'dh_c1_law_117.txt')

outf = open(RESULT_FILE, 'w')

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"


log("=" * 72)
log("[실험 #117] DH c₁ 보편 법칙 정밀 검증 + 해석적 유도 (v2)")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# DH 함수 구현
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
mpmath.mp.dps = 80  # t<200에는 충분

CHI_MOD5 = [0, 1, 1j, -1j, -1]
CHI_MOD5_BAR = [0, 1, -1j, 1j, -1]

sqrt5 = mpmath.sqrt(5)
kappa_DH = (mpmath.sqrt(10 - 2*sqrt5) - 2) / (sqrt5 - 1)
coeff_chi     = (1 - 1j * kappa_DH) / 2
coeff_chi_bar = (1 + 1j * kappa_DH) / 2


def dh_func(s):
    """Davenport-Heilbronn 함수 f(s)"""
    s = mpmath.mpc(s)
    return coeff_chi * mpmath.dirichlet(s, CHI_MOD5) + \
           coeff_chi_bar * mpmath.dirichlet(s, CHI_MOD5_BAR)


def Lambda_dh(s):
    """완비 DH 함수: Λ(s) = (5/π)^{s/2} · Γ((s+1)/2) · f(s)"""
    s = mpmath.mpc(s)
    return mpmath.power(5/mpmath.pi, s/2) * mpmath.gamma((s+1)/2) * dh_func(s)


def Lambda_deriv1(s, h=None):
    """Λ'(s) — 중앙차분"""
    if h is None:
        h = mpmath.mpf(10) ** (-(mpmath.mp.dps // 3))
    s = mpmath.mpc(s)
    return (Lambda_dh(s+h) - Lambda_dh(s-h)) / (2*h)


def Lambda_deriv2_at_zero(rho, h=None):
    """Λ''(ρ) at zero: (Λ(ρ+h)+Λ(ρ-h)-2Λ(ρ))/h² (Λ(ρ)≈0 이용)"""
    if h is None:
        h = mpmath.mpf(10) ** (-(mpmath.mp.dps // 4))
    rho = mpmath.mpc(rho)
    Lp = Lambda_dh(rho+h)
    Lm = Lambda_dh(rho-h)
    L0 = Lambda_dh(rho)  # ≈ 0
    return (Lp + Lm - 2*L0) / (h**2)


def refine_zero(sigma_guess, t_guess, dps_local=None):
    """findroot로 영점 정밀화"""
    old_dps = mpmath.mp.dps
    if dps_local is not None:
        mpmath.mp.dps = dps_local
    try:
        def f_sys(sv, tv):
            sc = mpmath.mpc(sv, tv)
            fv = dh_func(sc)
            return mpmath.re(fv), mpmath.im(fv)
        result = mpmath.findroot(f_sys, (mpmath.mpf(str(sigma_guess)),
                                          mpmath.mpf(str(t_guess))),
                                  maxsteps=100)
        sr, tr = float(result[0]), float(result[1])
        res = float(abs(dh_func(mpmath.mpc(sr, tr))))
        return sr, tr, res
    except Exception as e:
        return None
    finally:
        mpmath.mp.dps = old_dps


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Phase 1: 기존 4개 영점 정밀화 + c₁ 완전 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("=" * 72)
log("Phase 1: 기존 off-critical 영점 정밀화 + c₁ 완전 분석")
log("=" * 72)
log()

KNOWN_OFF_INIT = [
    (0.808517, 85.699348,  "OFF#1"),
    (0.650830, 114.163343, "OFF#2"),
    (0.574356, 166.479306, "OFF#3"),
    (0.724258, 176.702461, "OFF#4"),
]

log("영점 정밀화 (dps=80→120 for OFF#4):")
log(f"{'이름':<8} {'σ':>12} {'t':>16} {'|f|':>12} {'|σ-0.5|':>10} {'dps':>5}")
log("-" * 65)

all_off_zeros = []  # [(sigma, t, name)]

for sigma0, t0, name in KNOWN_OFF_INIT:
    # OFF#4는 dps=120에서 안정 (이전 확인)
    use_dps = 120 if name == "OFF#4" else 80
    result = refine_zero(sigma0, t0, dps_local=use_dps)
    if result is not None:
        sr, tr, res = result
        d = abs(sr - 0.5)
        log(f"{name:<8} {sr:>12.8f} {tr:>16.8f} {res:>12.3e} {d:>10.6f} {use_dps:>5}")
        if d > 0.01:  # off-critical 확인
            all_off_zeros.append((sr, tr, name))
        else:
            log(f"  ⚠️ {name}: |σ-0.5|={d:.4f} < 0.01 — on-critical 근접, 제외")
    else:
        log(f"  ⚠️ {name}: findroot 실패 — 초기값 사용")
        d = abs(sigma0 - 0.5)
        all_off_zeros.append((sigma0, t0, name))
        log(f"{name:<8} {sigma0:>12.8f} {t0:>16.8f} {'원초기값':>12} {d:>10.6f} {use_dps:>5}")

log()
log(f"영점 확정: {len(all_off_zeros)}개 {T()}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# δ-sweep → c₁_fit
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ δ-sweep → c₁_fit ━━━")
log("κδ² = 1 + c₁·δ + c₂·δ² + ...  →  (κδ²-1)/δ = c₁ + c₂·δ")
log("선형 회귀: (κδ²-1)/δ vs δ, y-절편 = c₁_fit")
log()

DELTA_VALUES = [0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.07]

header_deltas = "  ".join([f"δ={d}" for d in DELTA_VALUES])
log(f"{'이름':<8} {'σ':>10} {'t':>10}  " + "  ".join([f"{'δ='+str(d):>10}" for d in DELTA_VALUES]))
log("-" * (30 + len(DELTA_VALUES)*12))

c1_fit_dict = {}  # name → (c1_fit, r2, pts)

for sigma, t, name in all_off_zeros:
    kd2_row = []

    for dv in DELTA_VALUES:
        s_off = mpmath.mpc(sigma + dv, t)
        L_val = Lambda_dh(s_off)
        thresh = mpmath.mpf(10)**(-mpmath.mp.dps + 15)
        if abs(L_val) < thresh:
            kd2_row.append(float('nan'))
            continue
        h = mpmath.mpf(10)**(-(mpmath.mp.dps//3))
        L_deriv = (Lambda_dh(s_off+h) - Lambda_dh(s_off-h)) / (2*h)
        conn = L_deriv / L_val
        kd2 = float(abs(conn)**2) * dv**2
        kd2_row.append(kd2)

    row_str = "  ".join([f"{k:>10.6f}" if not np.isnan(k) else f"{'N/A':>10}" for k in kd2_row])
    log(f"{name:<8} {sigma:>10.6f} {t:>10.3f}  {row_str}")

    # c₁_fit 계산: (κδ²-1)/δ vs δ 선형 회귀
    valid = [(dv, k) for dv, k in zip(DELTA_VALUES, kd2_row)
             if not np.isnan(k) and k > 1.001]
    if len(valid) >= 4:
        deltas = np.array([d for d, _ in valid])
        kd2s   = np.array([k for _, k in valid])
        y = (kd2s - 1.0) / deltas
        # polyfit: y = c₁ + c₂·δ
        coeffs = np.polyfit(deltas, y, 1)
        c1_fit = coeffs[1]  # y-intercept
        c2_fit = coeffs[0]  # slope (c₂)
        y_pred = np.polyval(coeffs, deltas)
        ss_res = np.sum((y - y_pred)**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r2 = 1.0 - ss_res/ss_tot if ss_tot > 0 else float('nan')
        c1_fit_dict[name] = (c1_fit, r2, valid)
        log(f"  → c₁_fit = {c1_fit:.4f}, c₂_fit = {c2_fit:.4f}, R²={r2:.6f}")
    else:
        c1_fit_dict[name] = (float('nan'), float('nan'), valid)
        log(f"  ⚠️ 유효 포인트 부족 ({len(valid)}개 < 4)")

log()
log(f"{T()}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# c₁_analytic = Re(Λ''(ρ)/Λ'(ρ))
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ c₁_analytic = Re(Λ''(ρ)/Λ'(ρ)) 직접 계산 ━━━")
log("Λ'(ρ) = (Λ(ρ+h)-Λ(ρ-h))/(2h)")
log("Λ''(ρ) = (Λ(ρ+h)+Λ(ρ-h)-2Λ(ρ))/h²  [영점에서 Λ(ρ)≈0]")
log("c₁ = Re(Λ''(ρ)/Λ'(ρ))  → 이론: = Re(2A'(ρ)/A(ρ)) where Λ=(s-ρ)A(s)")
log()

c1_analytic_dict = {}  # name → c1_analytic

log(f"{'이름':<8} {'σ':>10} {'t':>10} {'|Λ_prime|':>14} {'|Λ_dprime|':>14} {'c1_analytic':>13} {'dps':>5}")
log("-" * 75)

for sigma, t, name in all_off_zeros:
    rho = mpmath.mpc(sigma, t)
    try:
        h1 = mpmath.mpf(10)**(-(mpmath.mp.dps//3))
        h2 = mpmath.mpf(10)**(-(mpmath.mp.dps//4))
        L1 = Lambda_deriv1(rho, h=h1)
        L2 = Lambda_deriv2_at_zero(rho, h=h2)
        absL1 = float(abs(L1))
        absL2 = float(abs(L2))
        if absL1 < 1e-100:
            log(f"{name:<8} {sigma:>10.6f} {t:>10.3f} {'Λ_prime≈0':>14} {'N/A':>14} {'N/A':>13}")
            c1_analytic_dict[name] = float('nan')
            continue
        ratio = L2 / L1
        c1_a = float(mpmath.re(ratio))
        c1_analytic_dict[name] = c1_a
        log(f"{name:<8} {sigma:>10.6f} {t:>10.3f} {absL1:>14.4e} {absL2:>14.4e} {c1_a:>13.4f} {mpmath.mp.dps:>5}")
    except Exception as e:
        log(f"{name:<8} {sigma:>10.6f} {t:>10.3f} ERROR: {e}")
        c1_analytic_dict[name] = float('nan')

log()
log(f"{T()}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 비교표 + 거울쌍 기여 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 비교표: c₁ 보편 법칙 검증 ━━━")
log()
log("이론:")
log("  c₁ = Re(Λ''(ρ)/Λ'(ρ)) ≈ 1/|σ₀-1/2|  (c₁·|σ-1/2| ≈ 1)")
log("  거울쌍 기여 (leading): c₁_mirror ≈ 1/(2|σ₀-1/2|)")
log("  배경 기여 = c₁ - c₁_mirror ≈ c₁_mirror (비슷해야)")
log()

col_w = [8, 10, 10, 9, 10, 10, 12, 13, 13, 14, 12]
header = (f"{'이름':<8} {'σ':>10} {'t':>10} {'d=|σ-½|':>9} "
          f"{'c1_fit':>10} {'c1_analy':>10} "
          f"{'1/d':>12} {'1/(2d)':>13} {'c1fit*d':>13} {'c1analy*d':>14} "
          f"{'fit/analy':>12}")
log(header)
log("-" * len(header))

results_table = []
for sigma, t, name in all_off_zeros:
    d = abs(sigma - 0.5)
    c1f, r2f, _ = c1_fit_dict.get(name, (float('nan'), float('nan'), []))
    c1a = c1_analytic_dict.get(name, float('nan'))
    inv_d  = 1.0/d if d > 0 else float('nan')
    inv_2d = 1.0/(2*d) if d > 0 else float('nan')
    prod_f = c1f * d if not np.isnan(c1f) else float('nan')
    prod_a = c1a * d if not np.isnan(c1a) else float('nan')

    if not np.isnan(c1f) and not np.isnan(c1a) and abs(c1a) > 0:
        ratio_fa = c1f / c1a
    else:
        ratio_fa = float('nan')

    def f(x, w=10):
        return f"{x:{w}.4f}" if not np.isnan(x) else f"{'N/A':>{w}}"

    log(f"{name:<8} {sigma:>10.6f} {t:>10.3f} {d:>9.4f} "
        f"{f(c1f):>10} {f(c1a):>10} "
        f"{f(inv_d, 12):>12} {f(inv_2d, 13):>13} {f(prod_f, 13):>13} {f(prod_a, 14):>14} "
        f"{f(ratio_fa, 12):>12}")

    results_table.append({
        'name': name, 'sigma': sigma, 't': t, 'd': d,
        'c1_fit': c1f, 'c1_analytic': c1a,
        'inv_d': inv_d, 'inv_2d': inv_2d,
        'prod_fit': prod_f, 'prod_analytic': prod_a,
        'ratio_fa': ratio_fa
    })

log()

# 통계
valid_prod_fit = [r['prod_fit'] for r in results_table if not np.isnan(r['prod_fit'])]
valid_prod_ana = [r['prod_analytic'] for r in results_table if not np.isnan(r['prod_analytic'])]
valid_ratio    = [r['ratio_fa'] for r in results_table if not np.isnan(r['ratio_fa'])]
valid_c1f      = [r['c1_fit'] for r in results_table if not np.isnan(r['c1_fit'])]
valid_c1a      = [r['c1_analytic'] for r in results_table if not np.isnan(r['c1_analytic'])]

log("【통계 요약】")
if valid_prod_fit:
    log(f"  c₁_fit · |σ-1/2|:      mean={np.mean(valid_prod_fit):.4f}, "
        f"std={np.std(valid_prod_fit):.4f}, N={len(valid_prod_fit)}")
if valid_prod_ana:
    log(f"  c₁_analytic · |σ-1/2|: mean={np.mean(valid_prod_ana):.4f}, "
        f"std={np.std(valid_prod_ana):.4f}, N={len(valid_prod_ana)}")
if valid_ratio:
    log(f"  c₁_fit / c₁_analytic:  mean={np.mean(valid_ratio):.4f}, "
        f"std={np.std(valid_ratio):.4f}")
log()

# 거울쌍 기여 분해
log("【거울쌍 기여 분해】")
log(f"  거울쌍 기여 이론: c₁_mirror = 1/(2|σ-1/2|) (leading term)")
log(f"  배경 기여: c₁_background = c₁ - c₁_mirror")
log()
log(f"{'이름':<8} {'c1_fit':>10} {'c1_mirror':>12} {'background':>12} {'비율':>10}")
log("-" * 55)

bg_list = []
for r in results_table:
    c1f = r['c1_fit']
    inv_2d = r['inv_2d']
    if not np.isnan(c1f) and not np.isnan(inv_2d):
        bg = c1f - inv_2d
        bg_ratio = bg / inv_2d if inv_2d > 0 else float('nan')
        bg_list.append(bg)
        log(f"{r['name']:<8} {c1f:>10.4f} {inv_2d:>12.4f} {bg:>12.4f} {bg_ratio:>10.4f}")
    else:
        log(f"{r['name']:<8} {'N/A':>10} {'N/A':>12} {'N/A':>12} {'N/A':>10}")

if bg_list:
    log(f"\n  배경 기여 평균: {np.mean(bg_list):.4f} ± {np.std(bg_list):.4f}")
    log(f"  (거울쌍과 배경이 비슷하면 c₁ ≈ 2×c₁_mirror = 1/|σ-1/2|)")

log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# log-log 스케일링 분석: c₁ ∝ |σ-1/2|^β
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ log-log 스케일링: c₁ ∝ |σ-1/2|^β ━━━")
log("예측: β = -1 (c₁ = A/|σ-1/2|, A≈1)")
log()

log_d_vals = []
log_c1f_vals = []
log_c1a_vals = []

for r in results_table:
    d = r['d']
    c1f = r['c1_fit']
    c1a = r['c1_analytic']
    if d > 0:
        if not np.isnan(c1f) and c1f > 0:
            log_d_vals.append(np.log(d))
            log_c1f_vals.append(np.log(c1f))
        if not np.isnan(c1a) and c1a > 0:
            log_c1a_vals.append(np.log(c1a))

if len(log_d_vals) >= 3:
    ld = np.array(log_d_vals)
    lc1f = np.array(log_c1f_vals)
    # log(c1) = β·log(d) + log(A)
    cf_fit = np.polyfit(ld, lc1f, 1)
    beta_f  = cf_fit[0]
    A_f     = np.exp(cf_fit[1])
    # R²
    y_pred = np.polyval(cf_fit, ld)
    r2_f = 1.0 - np.sum((lc1f - y_pred)**2) / np.sum((lc1f - np.mean(lc1f))**2)
    log(f"  c₁_fit 스케일링: β = {beta_f:.4f}, A = {A_f:.4f}, R² = {r2_f:.4f}")
    log(f"  이론: β=-1, A=1")

    if len(log_c1a_vals) >= 3:
        lc1a = np.array(log_c1a_vals)
        ca_fit = np.polyfit(ld[:len(lc1a)], lc1a, 1)
        beta_a = ca_fit[0]
        A_a    = np.exp(ca_fit[1])
        y_pred_a = np.polyval(ca_fit, ld[:len(lc1a)])
        r2_a = 1.0 - np.sum((lc1a - y_pred_a)**2) / np.sum((lc1a - np.mean(lc1a))**2)
        log(f"  c₁_analytic 스케일링: β = {beta_a:.4f}, A = {A_a:.4f}, R² = {r2_a:.4f}")
else:
    log(f"  ⚠️ log-log 분석 데이터 부족 ({len(log_d_vals)}개)")

log()
log(f"{T()}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Phase 2: 추가 영점 탐색 t∈[200,400] (시간 예산: 600초)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("=" * 72)
log("Phase 2: 추가 off-critical 영점 탐색 t∈[200,400] (예산 600s)")
log("=" * 72)
log()

phase2_start = time.time()
PHASE2_BUDGET = 600  # 10분
SCAN_DPS = 80        # 탐색 시 정밀도
SCAN_T_STEP = 2.0    # 코스한 t 스텝 (빠름)
SCAN_SIGMAS = [0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85]

log(f"설정: dps={SCAN_DPS}, T_STEP={SCAN_T_STEP}, σ={SCAN_SIGMAS}")
log(f"예상 격자점: {len(SCAN_SIGMAS)} × {int((400-200)/SCAN_T_STEP)} = "
    f"{len(SCAN_SIGMAS) * int((400-200)/SCAN_T_STEP)}")
log()

new_candidates = []
mpmath.mp.dps = SCAN_DPS

for sigma in SCAN_SIGMAS:
    if time.time() - phase2_start > PHASE2_BUDGET:
        log(f"⏱ 예산 초과 ({PHASE2_BUDGET}s) — 탐색 중단")
        break

    t_vals = np.arange(200.0, 400.0 + SCAN_T_STEP, SCAN_T_STEP)
    prev2 = prev1 = prev1_t = None
    n_found = 0

    log(f"  σ={sigma:.2f} 스캔... {T()}", )

    for t in t_vals:
        if time.time() - phase2_start > PHASE2_BUDGET:
            break
        try:
            sv = mpmath.mpc(sigma, t)
            absf = float(abs(dh_func(sv)))
        except Exception as e:
            print(f"WARNING: σ={sigma}, t={t:.1f}: {e}", flush=True)
            absf = float('inf')

        # 지역 최솟값 탐지
        if (prev2 is not None and prev1 is not None and
                prev1 < prev2 and prev1 < absf and prev1 < 0.3):
            new_candidates.append((sigma, prev1_t, prev1))
            n_found += 1

        prev2, prev1, prev1_t = prev1, absf, t

    log(f" → {n_found}개 후보")

log()
log(f"총 후보: {len(new_candidates)}개")
new_candidates.sort(key=lambda x: x[2])

# findroot 정밀화
log()
log("findroot 정밀화 (dps=80, 상위 30개 후보):")
new_off_zeros = []
findroot_tries = 0
findroot_fails = 0

for sigma_c, t_c, absf_c in new_candidates[:30]:
    if time.time() - phase2_start > PHASE2_BUDGET:
        log("⏱ 예산 초과 — findroot 중단")
        break

    # 기존 영점과 중복 체크
    all_known = [(s, t) for s, t, _ in all_off_zeros] + new_off_zeros
    if any(abs(t_c - t2) < 1.5 and abs(sigma_c - s2) < 0.05
           for s2, t2 in all_known):
        continue
    if abs(sigma_c - 0.5) < 0.01:
        continue

    findroot_tries += 1
    result = refine_zero(sigma_c, t_c, dps_local=80)
    if result is None:
        findroot_fails += 1
        continue

    sr, tr, res = result
    if (abs(sr - 0.5) > 0.02 and res < 1e-8 and 195 < tr < 405 and
            not any(abs(tr - t2) < 1.0 for _, t2 in new_off_zeros)):
        new_off_zeros.append((sr, tr))
        log(f"  ★ 신규: σ={sr:.6f}, t={tr:.6f}, |f|={res:.3e}, |σ-0.5|={abs(sr-0.5):.4f}")

log()
log(f"findroot: {findroot_tries}시도, {findroot_fails}실패, "
    f"{len(new_off_zeros)}개 신규 영점 발견")
log(f"Phase 2 소요: {time.time()-phase2_start:.1f}초")
log()

# 신규 영점 c₁ 분석 (빠른 버전)
if new_off_zeros:
    log("신규 영점 c₁ 분석 (dps=80, δ=0.01~0.05):")
    log(f"{'영점':<20} {'d=|σ-0.5|':>10} {'c1_fit':>10} {'c1·d':>10} {'c1_analytic':>13}")
    log("-" * 70)

    mpmath.mp.dps = 80

    for sr, tr in new_off_zeros:
        d = abs(sr - 0.5)
        # 빠른 δ-sweep (3개 δ만)
        pts = []
        for dv in [0.01, 0.02, 0.03]:
            s_off = mpmath.mpc(sr + dv, tr)
            L_val = Lambda_dh(s_off)
            if abs(L_val) < mpmath.mpf(10)**(-60):
                continue
            h = mpmath.mpf(10)**(-(mpmath.mp.dps//3))
            L_d = (Lambda_dh(s_off+h) - Lambda_dh(s_off-h)) / (2*h)
            kd2 = float(abs(L_d/L_val)**2) * dv**2
            pts.append((dv, kd2))

        if len(pts) >= 2:
            deltas = np.array([d for d, _ in pts])
            y = np.array([(k-1)/d for d, k in pts])
            c1f_new = np.mean(y)  # 단순 평균 (포인트 적음)
        else:
            c1f_new = float('nan')

        # c₁_analytic
        rho = mpmath.mpc(sr, tr)
        try:
            h1 = mpmath.mpf(10)**(-(mpmath.mp.dps//3))
            h2 = mpmath.mpf(10)**(-(mpmath.mp.dps//4))
            L1 = (Lambda_dh(rho+h1) - Lambda_dh(rho-h1)) / (2*h1)
            L2 = (Lambda_dh(rho+h2) + Lambda_dh(rho-h2) - 2*Lambda_dh(rho)) / (h2**2)
            c1a_new = float(mpmath.re(L2/L1)) if float(abs(L1)) > 1e-80 else float('nan')
        except Exception as e:
            print(f"WARNING: c1_analytic 계산 오류: {e}", flush=True)
            c1a_new = float('nan')

        prod = c1f_new * d if not np.isnan(c1f_new) else float('nan')
        c1a_str = f"{c1a_new:.4f}" if not np.isnan(c1a_new) else "N/A"
        c1f_str = f"{c1f_new:.4f}" if not np.isnan(c1f_new) else "N/A"
        prod_str = f"{prod:.4f}" if not np.isnan(prod) else "N/A"
        log(f"  σ={sr:.6f}, t={tr:.3f}  {d:>10.4f} {c1f_str:>10} {prod_str:>10} {c1a_str:>13}")

        # 전체 목록에 추가
        name_new = f"NEW#{len(all_off_zeros)+1}"
        all_off_zeros.append((sr, tr, name_new))
        c1_fit_dict[name_new] = (c1f_new, float('nan'), pts)
        c1_analytic_dict[name_new] = c1a_new

    log()
else:
    log("새 영점 발견 없음. 기존 4개 영점으로 분석 완료.")
    log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 최종 비교표 (전체)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("=" * 72)
log("최종 비교표 — 전체 off-critical 영점")
log("=" * 72)
log()

products_all = []
for sigma, t, name in all_off_zeros:
    d = abs(sigma - 0.5)
    c1f_val, _, _ = c1_fit_dict.get(name, (float('nan'), float('nan'), []))
    c1a_val = c1_analytic_dict.get(name, float('nan'))
    inv_d = 1.0/d if d > 0 else float('nan')

    prod_f = c1f_val * d if not np.isnan(c1f_val) else float('nan')
    prod_a = c1a_val * d if not np.isnan(c1a_val) else float('nan')

    rel_err_str = "N/A"
    if not np.isnan(c1f_val) and not np.isnan(c1a_val) and abs(c1a_val) > 0:
        rel_err = abs(c1f_val - c1a_val) / abs(c1a_val)
        rel_err_str = f"{rel_err:.4f}"

    c1f_str = f"{c1f_val:.4f}" if not np.isnan(c1f_val) else "N/A"
    c1a_str = f"{c1a_val:.4f}" if not np.isnan(c1a_val) else "N/A"
    inv_d_str = f"{inv_d:.4f}" if not np.isnan(inv_d) else "N/A"
    pf_str = f"{prod_f:.4f}" if not np.isnan(prod_f) else "N/A"
    pa_str = f"{prod_a:.4f}" if not np.isnan(prod_a) else "N/A"

    log(f"{name:<8} σ={sigma:.6f} t={t:>8.3f} d={d:.4f} "
        f"c1f={c1f_str:>8} c1a={c1a_str:>8} "
        f"1/d={inv_d_str:>7} c1f*d={pf_str:>7} c1a*d={pa_str:>7} "
        f"rel_err={rel_err_str}")

    if not np.isnan(prod_f):
        products_all.append(prod_f)

log()
if products_all:
    log(f"c₁_fit · |σ-1/2| — 전체 통계:")
    log(f"  N = {len(products_all)}")
    log(f"  mean = {np.mean(products_all):.4f}  (이론: 1.0)")
    log(f"  std  = {np.std(products_all):.4f}")
    log(f"  min  = {np.min(products_all):.4f}")
    log(f"  max  = {np.max(products_all):.4f}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 성공 기준 + 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 성공 기준 체크 ━━━")
log()

valid_c1f_main = [c1_fit_dict[n][0] for _, _, n in all_off_zeros[:4]
                  if not np.isnan(c1_fit_dict.get(n, (float('nan'),))[0])]
valid_c1a_main = [c1_analytic_dict[n] for _, _, n in all_off_zeros[:4]
                  if not np.isnan(c1_analytic_dict.get(n, float('nan')))]

sc1 = len(all_off_zeros) >= 7
sc2_count = sum(1 for _, _, n in all_off_zeros[:4]
                if not np.isnan(c1_fit_dict.get(n, (float('nan'),))[0]) and
                   not np.isnan(c1_analytic_dict.get(n, float('nan'))) and
                   abs(c1_fit_dict[n][0]) > 0 and abs(c1_analytic_dict[n]) > 0 and
                   abs(c1_fit_dict[n][0] - c1_analytic_dict[n]) / abs(c1_analytic_dict[n]) < 0.05)
sc2 = sc2_count >= 3

sc3_mean = np.mean(products_all) if products_all else float('nan')
sc3 = not np.isnan(sc3_mean) and abs(sc3_mean - 1.0) < 0.25

sc4 = len(log_d_vals) >= 3
sc5 = len(all_off_zeros) >= 4

log(f"SC1. 총 off-critical 영점 ≥7개: {'✅' if sc1 else f'❌ ({len(all_off_zeros)}개)'}")
log(f"SC2. c₁_fit vs c₁_analytic rel_err<5% (≥3 영점): {'✅' if sc2 else f'❌ ({sc2_count}개)'}")
log(f"SC3. c₁·|σ-1/2| mean∈[0.75,1.25]: {'✅' if sc3 else f'❌ (mean={sc3_mean:.4f})'}")
log(f"SC4. log-log β≈-1 분석 완료 (≥3 포인트): {'✅' if sc4 else f'❌ ({len(log_d_vals)}포인트)'}")
log(f"SC5. 비교표 출력 (≥4 영점): {'✅' if sc5 else f'❌ ({len(all_off_zeros)}개)'}")
log()

sc_pass = sum([sc1, sc2, sc3, sc4, sc5])
if sc_pass >= 4:
    verdict = "★★★ 강양성 — c₁ 보편 법칙 수치 확립"
elif sc_pass >= 3:
    verdict = "★★ 양성 — c₁ 법칙 부분 확인 (추가 영점 필요)"
elif sc_pass >= 2:
    verdict = "★ 조건부 — 주요 측정 완료, 통계 강화 필요"
else:
    verdict = "⚠️ 불충분"

log(f"성공 기준 {sc_pass}/5 통과")
log(f"실험 #117 판정: {verdict}")
log()

log("【수학적 의미】")
if sc3:
    log(f"  c₁·|σ₀-1/2| ≈ {sc3_mean:.3f} (≈1)")
    log("  → 완비 DH 함수 곡률 c₁이 임계선까지 거리의 역수")
    log("  → c₁ = Re(Λ''(ρ)/Λ'(ρ)) ≈ 1/|σ₀-1/2| 관계 확인")
    log("  → 해석: 거울쌍 1-ρ̄이 기여 + 배경 기여 ≈ 각 1/(2|σ-1/2|)")
else:
    log(f"  c₁·|σ₀-1/2| ≈ {sc3_mean:.3f} (이론 1.0에서 벗어남)")
    log("  → c₁ 보편 법칙 수치 확인 부분적")

log()
log("【on-critical 비교 (참조)】")
log("  On-critical (σ=0.5): c₁ = Re(Λ''(ρ)/Λ'(ρ)) = 0 (FE 대칭에 의해)")
log("  → κδ² = 1 + O(δ²) (c₁=0이므로 선형항 없음)")
log("  → 이것이 Theorem 4의 핵심: FE ↔ c₁=0 ↔ κδ² 스케일링 급변")

log()
log("=" * 72)
total = time.time() - START
log(f"총 소요: {total:.1f}초 ({total/60:.1f}분)")
log(f"결과 파일: {RESULT_FILE}")
log(f"종료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

outf.close()
print(f"\n결과 저장: {RESULT_FILE}")
