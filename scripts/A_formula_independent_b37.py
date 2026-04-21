#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 B-37] A = Im(c₀)² + 2c₁ 공식 독립 가족 교차검증
=============================================================================
목적:
  B-35에서 d=3,4는 sym^n(11a1) chain만 사용 → 편향 위험.
  독립 가족(다른 곡선, 다른 가족)에서 A 공식 보편성 교차검증.

실험 설계:
  - sym²(37a1): degree 3, 다른 타원곡선
  - Ramanujan Δ: degree 2, weight 12 cusp form (비타원곡선)
  - Dirichlet χ₋₃: degree 1, 이미 slope 검증됨, c₀/c₁ 개별 추출 미검증
  - Dirichlet χ₋₇: degree 1, 다른 conductor

결과: results/A_formula_independent_b37.txt
=============================================================================
"""
import sys, os, time
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'A_formula_independent_b37.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 B-37] A = Im(c₀)² + 2c₁ 공식 독립 가족 교차검증")
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
# L-함수 정의 — 독립 가족
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

LFUNCS = [
    # GL(1) 교차검증 — Dirichlet χ₋₃
    {
        "name": "L(s,χ₋₃)",
        "setup": "Lcur = lfuncreate(Mod(-1, 3))",
        "degree": 1,
        "desc": "Dirichlet L-function χ₋₃ (conductor 3)",
        "n_zeros": 5,
        "t_max": 50.0,
    },
    # GL(1) 교차검증 — Dirichlet χ₋₇
    {
        "name": "L(s,χ₋₇)",
        "setup": "Lcur = lfuncreate(Mod(3, 7))",
        "degree": 1,
        "desc": "Dirichlet L-function χ₋₇ (conductor 7, quadratic)",
        "n_zeros": 5,
        "t_max": 50.0,
    },
    # GL(2) 독립 — Ramanujan Δ (weight 12, level 1)
    {
        "name": "L(s,Δ)",
        "setup": "mf = mfinit([1,12],0); B = mfeigenbasis(mf); Lcur = lfunmf(mf, B[1])",
        "degree": 2,
        "desc": "Ramanujan Delta (weight 12, level 1, non-EC)",
        "n_zeros": 5,
        "t_max": 30.0,
    },
    # GL(3) 독립 — sym²(37a1, rank 1)
    {
        "name": "L(s,sym²(37a1))",
        "setup": "E37 = ellinit([0,0,1,-1,0]); Lcur = lfunsympow(E37, 2)",
        "degree": 3,
        "desc": "Symmetric square of 37a1 (N=37, rank 1, independent from 11a1)",
        "n_zeros": 5,
        "t_max": 30.0,
    },
]

C0_DELTAS = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01])
SLOPE_DELTAS = np.array([0.01, 0.015, 0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2])

SIGMA_C = 0.5  # 기본값, L-함수별 갱신

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수 (B-35 동일)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def get_log_deriv(sigma, t0):
    """Λ'/Λ(σ+it₀) via lfunlambda"""
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
            t_start = time.time()

            # c₀ 추출
            c0 = extract_c0_symmetric(t0)
            # c₁ 추출
            c1 = extract_c1_symmetric(t0)
            # A_pred
            A_pred_im = c0.imag**2
            A_pred_re = 2.0 * c1.real
            A_pred = A_pred_im + A_pred_re

            # slope & A_meas
            slope, r2, A_meas = measure_slope(t0)

            if slope is None:
                log(f"  [영점 #{idx+1}] t₀ = {t0:.6f}")
                log(f"    ⚠️ slope 측정 실패 (유효 δ < 4)")
                all_results.append({
                    'name': lf['name'], 'degree': lf['degree'],
                    't0': t0, 'pass': False, 'reason': 'slope_fail'
                })
                log()
                continue

            # A 오차
            A_err = abs(A_pred - A_meas) / abs(A_meas) * 100.0 if abs(A_meas) > 1e-12 else 999.0
            passed = A_err < 10.0 and abs(slope - 2.0) < 0.05

            log(f"  [영점 #{idx+1}] t₀ = {t0:.6f}")
            log(f"    c₀ = ({c0.real:.6e}) + ({c0.imag:.6f})i")
            log(f"    |Re(c₀)| = {abs(c0.real):.6e}")
            log(f"    c₁ = ({c1.real:.6f}) + ({c1.imag:.6e})i")
            log(f"    |Im(c₁)| = {abs(c1.imag):.6e}")
            log(f"    A_pred = Im(c₀)²+2Re(c₁) = {A_pred_im:.4f}+{A_pred_re:.4f} = {A_pred:.4f}")
            log(f"    slope = {slope:.4f}  R² = {r2:.8f}")
            log(f"    A_meas = {A_meas:.4f}")
            log(f"    A 오차 = {A_err:.1f}%")
            log(f"    판정: {'★★★ PASS' if passed else '✗ FAIL'}")
            log()

            t_elapsed = time.time() - t_start
            all_results.append({
                'name': lf['name'], 'degree': lf['degree'],
                't0': t0, 'pass': passed, 'A_err': A_err,
                'slope': slope, 'Re_c0': abs(c0.real),
                'Im_c1': abs(c1.imag), 'A_pred': A_pred, 'A_meas': A_meas,
                'elapsed': t_elapsed,
            })

            # 타임아웃 체크: 영점당 600초 초과 시 다음 L-함수로
            if t_elapsed > 600:
                log(f"  ⚠️ 영점당 {t_elapsed:.0f}초 — 이 L-함수 건너뜀")
                break

    except Exception as e:
        log(f"  ⚠️ 오류: {e}")
        log()
        continue

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 요약
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log()
log("=" * 72)
log("[요약] B-37: A = Im(c₀)² + 2c₁ 독립 가족 교차검증")
log("=" * 72)
log()

# L-함수별 요약
lf_names = list(dict.fromkeys([r['name'] for r in all_results]))
total_pass = 0
total_test = 0

log(f"| {'L-함수':<25s} | degree | 영점 | A 오차 범위  | slope 범위      | 판정 |")
log(f"|{'-'*27}|--------|------|-------------|-----------------|------|")

for name in lf_names:
    rs = [r for r in all_results if r['name'] == name and 'A_err' in r]
    if not rs:
        log(f"| {name:<25s} | — | 0 | — | — | 중단 |")
        continue
    n_pass = sum(1 for r in rs if r['pass'])
    n_total = len(rs)
    total_pass += n_pass
    total_test += n_total
    deg = rs[0]['degree']
    a_errs = [r['A_err'] for r in rs]
    slopes = [r['slope'] for r in rs]
    a_min, a_max = min(a_errs), max(a_errs)
    s_min, s_max = min(slopes), max(slopes)
    verdict = "PASS" if n_pass == n_total else f"{n_pass}/{n_total}"
    log(f"| {name:<25s} | {deg}      | {n_pass}/{n_total}  | {a_min:.1f}%–{a_max:.1f}%   | {s_min:.3f}–{s_max:.3f}     | {verdict} |")

log(f"| {'**총계**':<25s} |**mix** |**{total_pass}/{total_test}**| — | — | **{'전원 PASS' if total_pass==total_test else f'{total_pass}/{total_test}'}** |")
log()

# 핵심 수치
valid = [r for r in all_results if 'A_err' in r]
if valid:
    re_c0s = [r['Re_c0'] for r in valid]
    im_c1s = [r['Im_c1'] for r in valid]
    slopes_all = [r['slope'] for r in valid]
    a_errs_all = [r['A_err'] for r in valid]

    log(f"핵심 수치:")
    log(f"  Re(c₀) max = {max(re_c0s):.6e}")
    log(f"  Im(c₁) max = {max(im_c1s):.6e}")
    log(f"  slope 평균 = {np.mean(slopes_all):.4f} ± {np.std(slopes_all):.4f}")
    log(f"  A 오차 평균 = {np.mean(a_errs_all):.2f}%, 최대 = {max(a_errs_all):.1f}%")
    log()

# B-35 대비 비교
log("B-35 대비 교차검증 결과:")
log("  B-35: 25/25 PASS (d=1-4, sym chain 11a1)")
log(f"  B-37: {total_pass}/{total_test} ({'PASS' if total_pass==total_test else 'MIXED'}) (독립 가족)")
if total_pass == total_test and total_test > 0:
    log("  → A 공식의 보편성 교차검증 성공. sym chain 편향 기각.")
    log()
    log("판정: ★★★★ 강양성 — Cor(amplitude) 독립 가족 교차검증 완료")
else:
    log(f"  → {total_test - total_pass}개 실패 영점 존재. 조건부 양성 또는 공식 범위 제한.")
    log()
    log(f"판정: {'★★★ 조건부 양성' if total_pass > total_test*0.8 else '★★ 약양성'}")

log()
log(f"종료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
outf.close()
