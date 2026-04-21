#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 B-35] A = Im(c₀)² + 2c₁ 공식 degree≥2 수치 검증
=============================================================================
목적:
  Slope Universality Theorem (thm:slopeuniv)의 A 공식을 degree 1 이상에서 검증.
  #214/B-34는 degree 1 (ζ, 디리클레)에서만 검증. GL(2)/GL(3)/GL(4)에서 미확인.

실험 설계:
  - GL(2) degree 2: 11a1, 37a1, Ramanujan Δ
  - GL(3) degree 3: sym²(11a1)
  - GL(4) degree 4: sym³(Δ)
  각 L-함수에서 3-5영점, Laurent 계수 c₀, c₁ 추출 → A_pred vs A_meas 비교.

결과: results/A_formula_highdeg_b35.txt
=============================================================================
"""
import sys, os, time
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'A_formula_highdeg_b35.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 B-35] A = Im(c₀)² + 2c₁ 공식 degree≥2 검증")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PARI 초기화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import cypari2
gp = cypari2.Pari()
gp.allocatemem(1000 * 1024 * 1024)
gp("default(realprecision, 150)")
log(f"{T()} PARI OK, realprecision=150")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# L-함수 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

LFUNCS = [
    # GL(1) 대조군
    {
        "name": "ζ(s)",
        "setup": "Lcur = lfuncreate(1)",
        "degree": 1,
        "desc": "Riemann zeta",
        "n_zeros": 5,
        "t_max": 50.0,
    },
    # GL(2) degree 2 — a-invariants 직접 사용 (elldata 불필요)
    {
        "name": "L(s,11a1)",
        "setup": 'E11 = ellinit([0,-1,1,-10,-20]); Lcur = lfuncreate(E11)',
        "degree": 2,
        "desc": "Elliptic curve 11a1 (N=11)",
        "n_zeros": 5,
        "t_max": 30.0,
    },
    {
        "name": "L(s,37a1)",
        "setup": 'E37 = ellinit([0,0,1,-1,0]); Lcur = lfuncreate(E37)',
        "degree": 2,
        "desc": "Elliptic curve 37a1 (N=37, rank 1)",
        "n_zeros": 5,
        "t_max": 30.0,
    },
    # GL(3) degree 3
    {
        "name": "L(s,sym²(11a1))",
        "setup": 'E11 = ellinit([0,-1,1,-10,-20]); Lcur = lfunsympow(E11, 2)',
        "degree": 3,
        "desc": "Symmetric square of 11a1",
        "n_zeros": 5,
        "t_max": 30.0,
    },
    # GL(4) degree 4 — sym³(11a1)
    {
        "name": "L(s,sym³(11a1))",
        "setup": 'E11 = ellinit([0,-1,1,-10,-20]); Lcur = lfunsympow(E11, 3)',
        "degree": 4,
        "desc": "Symmetric cube of 11a1",
        "n_zeros": 5,
        "t_max": 15.0,
    },
    # GL(5) degree 5 — sym⁴(11a1)
    {
        "name": "L(s,sym⁴(11a1))",
        "setup": 'E11 = ellinit([0,-1,1,-10,-20]); Lcur = lfunsympow(E11, 4)',
        "degree": 5,
        "desc": "Symmetric fourth power of 11a1",
        "n_zeros": 5,
        "t_max": 15.0,
    },
]

C0_DELTAS = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01])
SLOPE_DELTAS = np.array([0.01, 0.015, 0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2])

# 현재 L-함수의 임계선 중심 (L[4]/2로 결정)
SIGMA_C = 0.5  # 기본값, L-함수별 갱신

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def get_log_deriv(sigma, t0):
    """Λ'/Λ(σ+it₀) via lfunlambda (Lcur 직접 사용, lfuninit 불필요)"""
    gp(f"s_eval = {sigma:.20f} + I*{t0:.20f}")
    gp("Lv = lfunlambda(Lcur, s_eval)")
    gp("dLv = lfunlambda(Lcur, s_eval, 1)")
    gp("r_val = dLv / Lv")
    return complex(float(gp("real(r_val)")), float(gp("imag(r_val)")))

def extract_c0_symmetric(t0):
    """c₀ 추출 — (f(σ_c+δ) + f(σ_c-δ))/2 ≈ c₀ + O(δ²)"""
    c0_estimates = []
    for d in C0_DELTAS:
        f_plus = get_log_deriv(SIGMA_C + d, t0)
        f_minus = get_log_deriv(SIGMA_C - d, t0)
        c0_est = (f_plus + f_minus) / 2.0
        c0_estimates.append(c0_est)
    c0_val = np.mean(c0_estimates[:3])
    return c0_val

def extract_c1_symmetric(t0):
    """c₁ 추출 — [(f(σ_c+δ) - f(σ_c-δ))/2 - 1/δ] / δ ≈ c₁ + O(δ²)"""
    c1_estimates = []
    for d in C0_DELTAS:
        f_plus = get_log_deriv(SIGMA_C + d, t0)
        f_minus = get_log_deriv(SIGMA_C - d, t0)
        antisym = (f_plus - f_minus) / 2.0
        c1_est = (antisym - 1.0/d) / d
        c1_estimates.append(c1_est)
    c1_val = np.mean(c1_estimates[:3])
    return c1_val

def measure_slope(t0):
    """κδ² log-log slope 측정 + A_meas"""
    kd2_vals = []
    for d in SLOPE_DELTAS:
        ratio = get_log_deriv(SIGMA_C + d, t0)
        kappa = abs(ratio)**2
        kd2_vals.append(kappa * d**2)

    valid_d = []
    valid_lkd2 = []
    for d, kd2 in zip(SLOPE_DELTAS, kd2_vals):
        if kd2 > 1.001:
            valid_d.append(np.log(d))
            valid_lkd2.append(np.log(kd2 - 1.0))

    if len(valid_d) < 4:
        return None, None, None

    coeffs_full = np.polyfit(valid_d, valid_lkd2, 1)
    slope_full = coeffs_full[0]
    pred_full = np.polyval(coeffs_full, valid_d)
    ss_res = np.sum((np.array(valid_lkd2) - pred_full)**2)
    ss_tot = np.sum((np.array(valid_lkd2) - np.mean(valid_lkd2))**2)
    r2_full = 1.0 - ss_res/ss_tot if ss_tot > 1e-30 else 0.0

    A_meas = np.exp(coeffs_full[1])
    return slope_full, r2_full, A_meas

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

all_results = []

for lf in LFUNCS:
    log(f"{T()} ━━━ {lf['name']} (d={lf['degree']}) ━━━")
    log(f"  {lf['desc']}")

    try:
        # L-함수 초기화
        setup_cmd = lf['setup'].strip().replace('\n', '; ')
        gp(setup_cmd)

        # 임계선 중심 추출: k = Lcur[4], σ_c = k/2
        k_val = float(str(gp("Lcur[4]")))
        SIGMA_C = k_val / 2.0
        log(f"  k = {k_val}, σ_c = {SIGMA_C}")

        # 영점 탐색
        gp(f"zeros = lfunzeros(Lcur, {lf['t_max']:.1f})")
        n_zeros_found = int(gp("length(zeros)"))

        if n_zeros_found == 0:
            log(f"  ⚠️ 영점 없음 (t < {lf['t_max']})")
            continue

        zeros_to_use = min(lf['n_zeros'], n_zeros_found)
        zero_list = []
        for j in range(1, zeros_to_use + 1):
            zt = float(gp(f"zeros[{j}]"))
            zero_list.append(zt)

        log(f"  영점 {n_zeros_found}개, 사용: {[f'{z:.4f}' for z in zero_list]}")
        log()

        for idx, t0 in enumerate(zero_list):
            # c₀ 추출
            c0 = extract_c0_symmetric(t0)

            # c₁ 추출
            c1 = extract_c1_symmetric(t0)

            # A_pred
            A_pred_im = c0.imag**2
            A_pred_re = 2.0 * c1.real
            A_pred = A_pred_im + A_pred_re

            # slope + A_meas
            slope, r2, A_meas = measure_slope(t0)

            if slope is None:
                log(f"  [영점 #{idx+1}] t₀ = {t0:.6f} — slope 측정 실패")
                continue

            A_err = abs(A_pred - A_meas) / abs(A_meas) * 100 if A_meas else float('nan')

            log(f"  [영점 #{idx+1}] t₀ = {t0:.6f}")
            log(f"    c₀ = ({c0.real:.6e}) + ({c0.imag:.6f})i")
            log(f"    |Re(c₀)| = {abs(c0.real):.6e}")
            log(f"    c₁ = ({c1.real:.6f}) + ({c1.imag:.6e})i")
            log(f"    |Im(c₁)| = {abs(c1.imag):.6e}")
            log(f"    A_pred = Im(c₀)²+2Re(c₁) = {A_pred_im:.4f}+{A_pred_re:.4f} = {A_pred:.4f}")
            log(f"    slope = {slope:.4f}  R² = {r2:.8f}")
            log(f"    A_meas = {A_meas:.4f}")
            log(f"    A 오차 = {A_err:.1f}%")

            # 판정
            re_c0_ok = abs(c0.real) < 1e-6
            im_c1_ok = abs(c1.imag) < 1e-6
            slope_ok = abs(slope - 2.0) < 0.05
            a_ok = A_err < 5.0

            status = "★★★ PASS" if (re_c0_ok and im_c1_ok and slope_ok and a_ok) else "⚠️ CHECK"
            log(f"    판정: {status}")
            log()

            all_results.append({
                "name": lf['name'],
                "degree": lf['degree'],
                "t0": t0,
                "Re_c0": c0.real,
                "Im_c0": c0.imag,
                "Re_c1": c1.real,
                "Im_c1": c1.imag,
                "A_pred": A_pred,
                "A_meas": A_meas,
                "A_err": A_err,
                "slope": slope,
                "r2": r2,
                "pass": re_c0_ok and im_c1_ok and slope_ok and a_ok,
            })

    except Exception as e:
        log(f"  ❌ 오류: {e}")
        log()
        continue

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 종합 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("종합 분석")
log("=" * 72)
log()

# Degree별 통계
for deg in sorted(set(r['degree'] for r in all_results)):
    deg_results = [r for r in all_results if r['degree'] == deg]
    n = len(deg_results)
    n_pass = sum(1 for r in deg_results if r['pass'])
    slopes = [r['slope'] for r in deg_results]
    a_errs = [r['A_err'] for r in deg_results]
    re_c0s = [abs(r['Re_c0']) for r in deg_results]

    log(f"━━━ Degree {deg} (GL({deg})) ━━━")
    log(f"  N = {n} 영점, PASS = {n_pass}/{n}")
    log(f"  max|Re(c₀)| = {max(re_c0s):.3e}")
    log(f"  mean slope = {np.mean(slopes):.4f} ± {np.std(slopes):.4f}")
    log(f"  mean A_err = {np.mean(a_errs):.1f}% (max {max(a_errs):.1f}%)")
    log()

# 전체 결과표
log("━━━ 상세 결과표 ━━━")
header = f"{'L-함수':30s} {'d':>2s} {'t₀':>10s} {'|Re(c₀)|':>12s} {'Im(c₀)':>10s} {'Re(c₁)':>10s} {'|Im(c₁)|':>10s} {'slope':>7s} {'A_pred':>8s} {'A_meas':>8s} {'err%':>6s} {'P/F':>4s}"
log(header)
log("-" * len(header))
for r in all_results:
    pf = "PASS" if r['pass'] else "FAIL"
    log(f"{r['name']:30s} {r['degree']:>2d} {r['t0']:>10.4f} {abs(r['Re_c0']):>12.3e} {r['Im_c0']:>10.4f} {r['Re_c1']:>10.4f} {abs(r['Im_c1']):>10.3e} {r['slope']:>7.4f} {r['A_pred']:>8.4f} {r['A_meas']:>8.4f} {r['A_err']:>5.1f}% {pf:>4s}")

log()
total_pass = sum(1 for r in all_results if r['pass'])
total = len(all_results)
log(f"전체: {total_pass}/{total} PASS")

# 최종 판정
all_pass = total_pass == total
max_a_err = max(r['A_err'] for r in all_results) if all_results else 999
log()
log("=" * 72)
log("최종 판정")
log("=" * 72)

if all_pass and max_a_err < 5.0:
    log(f"★★★ 강양성 — A=Im(c₀)²+2c₁ 공식 degree 1-{max(r['degree'] for r in all_results)} 검증 완료")
    log(f"  {total}영점 전원 PASS. max A_err = {max_a_err:.1f}%.")
elif total_pass >= total * 0.8:
    log(f"★★ 양성 — 대부분 통과 ({total_pass}/{total}), max A_err = {max_a_err:.1f}%")
else:
    log(f"★ 조건부/음성 — {total_pass}/{total} PASS, max A_err = {max_a_err:.1f}%")

log()
elapsed = time.time() - START
log(f"총 소요시간: {elapsed:.1f}초")
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
outf.close()
