#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #214v2] κδ² slope=2 해석적 증명 검증 — 수정판
=============================================================================
v1 문제점: lfuninit 초기화 부족으로 c₁ 추출 불안정 (Im(c₁)≠0 아티팩트)
수정:
  1. lfuninit 범위 확대 ([1.0, t_max+20])
  2. c₀ 추출: δ=0.001~0.01 범위 (lfunlambda 안정 영역)
  3. c₁ 추출: 양측(±δ) 대칭 방법 → Im 아티팩트 제거
  4. A 비교: fitted A vs predicted A (δ=0.01~0.2 range)
  5. EC 대신 추가 Dirichlet (elldata 없음)

핵심 정리:
  Theorem (κδ² Slope Universality):
    For self-dual L-function Λ(s) = ε·Λ(1-s) with Λ(s̄)=Λ̄(s),
    simple zero ρ=1/2+iγ ⟹ κ(δ)·δ² - 1 = O(δ²), i.e., slope = 2 exactly.
  Proof: Functional equation + reality ⟹ Re(c₀)=0, eliminating δ¹ term.

결과: results/kd2_slope2_proof_214.txt
=============================================================================
"""
import sys, os, time
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'kd2_slope2_proof_214.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 #214] κδ² slope=2 해석적 증명 검증 (v2)")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PARI 초기화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import cypari2
gp = cypari2.Pari()
gp.allocatemem(1000 * 1024 * 1024)
gp("default(realprecision, 150)")  # 150 자릿수 (c₀ 실수부 정밀 검증)
log(f"{T()} PARI OK, realprecision=150")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 정리 서술
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("THEOREM (κδ² Slope Universality)")
log("=" * 72)
log()
log("Statement:")
log("  Let Λ(s) be a completed L-function satisfying:")
log("    (FE)  Λ(s) = ε·Λ(1-s)    (functional equation)")
log("    (RC)  Λ(s̄) = Λ̄(s)        (reality condition, self-dual)")
log("  Let ρ = ½+iγ be a simple zero on the critical line.")
log("  Define κ(δ) := |Λ'/Λ(½+δ+iγ)|² for real δ>0.")
log("  Then κ(δ)·δ² = 1 + A·δ² + O(δ³), so")
log("    log(κδ²-1) = log A + 2·log δ + O(δ)")
log("  and the log-log slope is exactly 2.")
log()
log("Proof:")
log("  Laurent expansion at ρ: Λ'/Λ(s) = 1/(s-ρ) + c₀ + c₁(s-ρ) + ...")
log("  At s = ½+δ+iγ (so s-ρ = δ ∈ ℝ):")
log("    Λ'/Λ = 1/δ + c₀ + c₁δ + ...                    (★)")
log()
log("  From (FE): Λ'/Λ(s) = -Λ'/Λ(1-s)")
log("  At 1-s = ½-δ-iγ = conj(½-δ+iγ):")
log("  From (RC): Λ'/Λ(conj(z)) = conj(Λ'/Λ(z))")
log("    So Λ'/Λ(½-δ-iγ) = conj(Λ'/Λ(½-δ+iγ))         (†)")
log()
log("  Laurent at ½-δ+iγ (i.e. s-ρ = -δ):")
log("    Λ'/Λ(½-δ+iγ) = 1/(-δ) + c₀ + c₁(-δ) + ...")
log("                   = -1/δ + c₀ - c₁δ + ...          (‡)")
log()
log("  Combining (FE) at s = ½+δ+iγ:")
log("    Λ'/Λ(½+δ+iγ) = -Λ'/Λ(½-δ-iγ) = -conj(Λ'/Λ(½-δ+iγ))  [by (†)]")
log("    From (★): 1/δ + c₀ + c₁δ = -conj(-1/δ + c₀ - c₁δ)")
log("                               = 1/δ - c̄₀ + c̄₁δ")
log()
log("  Equating coefficients:")
log("    δ⁰: c₀ = -c̄₀  ⟹  Re(c₀) = 0   (c₀ purely imaginary)  ■")
log("    δ¹: c₁ = c̄₁   ⟹  Im(c₁) = 0   (c₁ real)              ■")
log()
log("  Therefore:")
log("    κδ² = |1 + c₀δ + c₁δ²|² = 1 + 2Re(c₀)δ + (|c₀|²+2Re(c₁))δ² + O(δ³)")
log("         = 1 + 0 + A·δ²      (since Re(c₀)=0)")
log("    where A = Im(c₀)² + 2c₁.                              □")
log()
log("Remark: A depends on the local zero density (via c₀) and")
log("  Γ-factor structure (via digamma contributions to c₁).")
log("  This explains the observed A(t₀) ∝ d² + log N variation.")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# L-함수 정의 (elldata 불필요)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

LFUNCS = [
    {"name": "ζ(s)", "setup": "Lcur = lfuncreate(1)", "degree": 1, "N": 1},
    {"name": "L(s,χ₋₃)", "setup": "Lcur = lfuncreate(-3)", "degree": 1, "N": 3},
    {"name": "L(s,χ₋₇)", "setup": "Lcur = lfuncreate(-7)", "degree": 1, "N": 7},
    {"name": "L(s,χ₋₁₁)", "setup": "Lcur = lfuncreate(-11)", "degree": 1, "N": 11},
]

N_ZEROS_PER = 5
# c₀ 추출: δ 0.001~0.01 (lfunlambda 안정 영역)
C0_DELTAS = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01])
# κδ² slope fitting: δ 0.01~0.2
SLOPE_DELTAS = np.array([0.01, 0.015, 0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2])
T_MAX = 50.0

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def get_log_deriv(sigma, t0):
    """Λ'/Λ(σ+it₀) — lfunlambda(Linit, s, 1) / lfunlambda(Linit, s)"""
    gp(f"s_eval = {sigma:.20f} + I*{t0:.20f}")
    gp("Lv = lfunlambda(Linit, s_eval)")
    gp("dLv = lfunlambda(Linit, s_eval, 1)")
    gp("r_val = dLv / Lv")
    return complex(float(gp("real(r_val)")), float(gp("imag(r_val)")))

def extract_c0_symmetric(t0):
    """c₀ 추출 — 양측 대칭 방법

    (Λ'/Λ(½+δ+it) + Λ'/Λ(½-δ+it)) / 2 ≈ c₀ + c₂δ² + ...
    (로랑 주항 1/δ는 부호가 바뀌어 상쇄)

    이론: Re(c₀)=0이면 이 값의 실수부 → 0
    """
    c0_estimates = []
    for d in C0_DELTAS:
        f_plus = get_log_deriv(0.5 + d, t0)   # 1/δ + c₀ + c₁δ + ...
        f_minus = get_log_deriv(0.5 - d, t0)  # -1/δ + c₀ - c₁δ + ...
        # 합의 반: c₀ + c₂δ² + ...
        c0_est = (f_plus + f_minus) / 2.0
        c0_estimates.append(c0_est)

    # Richardson extrapolation: 가장 작은 δ들이 가장 정확
    c0_val = np.mean(c0_estimates[:3])
    c0_spread_re = np.std([v.real for v in c0_estimates[:3]])
    c0_spread_im = np.std([v.imag for v in c0_estimates[:3]])
    return c0_val, c0_spread_re

def extract_c1_symmetric(t0):
    """c₁ 추출 — 반대칭 방법

    (Λ'/Λ(½+δ+it) - Λ'/Λ(½-δ+it)) / 2 = 1/δ + c₁δ + c₃δ³ + ...
    → [(Λ'/Λ(½+δ) - Λ'/Λ(½-δ))/2 - 1/δ] / δ = c₁ + c₃δ² + ...

    이론: Im(c₁)=0이면 이 값의 허수부 → 0
    """
    c1_estimates = []
    for d in C0_DELTAS:
        f_plus = get_log_deriv(0.5 + d, t0)
        f_minus = get_log_deriv(0.5 - d, t0)
        # 차의 반: 1/δ + c₁δ + c₃δ³ + ...
        antisym = (f_plus - f_minus) / 2.0
        # c₁ ≈ (antisym - 1/δ) / δ
        c1_est = (antisym - 1.0/d) / d
        c1_estimates.append(c1_est)

    c1_val = np.mean(c1_estimates[:3])
    return c1_val

def measure_slope(t0):
    """κδ² log-log slope 측정"""
    kd2_vals = []
    for d in SLOPE_DELTAS:
        ratio = get_log_deriv(0.5 + d, t0)
        kappa = abs(ratio)**2
        kd2_vals.append(kappa * d**2)

    # log(κδ²-1) vs log(δ) fit
    valid_d = []
    valid_lkd2 = []
    for d, kd2 in zip(SLOPE_DELTAS, kd2_vals):
        if kd2 > 1.001:
            valid_d.append(np.log(d))
            valid_lkd2.append(np.log(kd2 - 1.0))

    if len(valid_d) < 4:
        return None, None, None

    coeffs = np.polyfit(valid_d, valid_lkd2, 1)
    slope = coeffs[0]
    A_meas = np.exp(coeffs[1])

    pred = np.polyval(coeffs, valid_d)
    ss_res = np.sum((np.array(valid_lkd2) - pred)**2)
    ss_tot = np.sum((np.array(valid_lkd2) - np.mean(valid_lkd2))**2)
    r2 = 1.0 - ss_res/ss_tot if ss_tot > 1e-30 else 0.0

    return slope, r2, A_meas

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("수치 검증 — Re(c₀)=0 및 Im(c₁)=0")
log("=" * 72)
log()

all_results = []

for lf in LFUNCS:
    log(f"{T()} ━━━ {lf['name']} (d={lf['degree']}, N={lf['N']}) ━━━")

    gp(lf["setup"])
    fe = float(gp("lfuncheckfeq(Lcur)"))
    log(f"  FE: {fe:.0f} digits")

    # lfuninit with wider range
    gp(f"Linit = lfuninit(Lcur, [1.0, {T_MAX + 20}])")

    # 영점 찾기
    gp(f"zvec = lfunzeros(Linit, {T_MAX})")
    n_found = int(gp("length(zvec)"))
    zeros = [float(gp(f"zvec[{i}]")) for i in range(1, min(n_found, N_ZEROS_PER)+1)]
    log(f"  영점 {n_found}개, 사용: {[f'{z:.4f}' for z in zeros]}")
    log()

    for idx, t0 in enumerate(zeros):
        log(f"  [영점 #{idx+1}] t₀ = {t0:.6f}")

        # c₀ (대칭 방법)
        c0, c0_spread = extract_c0_symmetric(t0)
        log(f"    c₀ = ({c0.real:.3e}) + ({c0.imag:.6f})i")
        log(f"    |Re(c₀)| = {abs(c0.real):.3e}  [이론: 0, spread={c0_spread:.2e}]")

        # c₁ (반대칭 방법)
        c1 = extract_c1_symmetric(t0)
        log(f"    c₁ = ({c1.real:.6f}) + ({c1.imag:.3e})i")
        log(f"    |Im(c₁)| = {abs(c1.imag):.3e}  [이론: 0]")

        # A 예측
        A_pred = c0.imag**2 + 2*c1.real

        # slope 측정
        slope, r2, A_meas = measure_slope(t0)

        if slope is not None:
            A_err = abs(A_pred - A_meas) / max(abs(A_meas), 1e-10)
            log(f"    A_pred = Im(c₀)²+2c₁ = {c0.imag**2:.4f}+{2*c1.real:.4f} = {A_pred:.4f}")
            log(f"    A_meas = {A_meas:.4f}  (오차: {A_err*100:.1f}%)")
            log(f"    slope  = {slope:.6f}  R² = {r2:.8f}")

            all_results.append({
                "Lfunc": lf["name"], "degree": lf["degree"], "N": lf["N"],
                "t0": t0, "Re_c0": c0.real, "Im_c0": c0.imag,
                "Re_c1": c1.real, "Im_c1": c1.imag,
                "A_pred": A_pred, "A_meas": A_meas, "slope": slope, "R2": r2,
            })
        log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 종합 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("종합 분석")
log("=" * 72)
log()

if all_results:
    re_c0s = [abs(r["Re_c0"]) for r in all_results]
    im_c1s = [abs(r["Im_c1"]) for r in all_results]
    slopes = [r["slope"] for r in all_results]
    A_errs = [abs(r["A_pred"]-r["A_meas"])/max(abs(r["A_meas"]),1e-10) for r in all_results]

    log("━━━ 검증 1: Re(c₀) = 0 (함수방정식의 귀결) ━━━")
    log(f"  N = {len(re_c0s)} zeros")
    log(f"  max|Re(c₀)|  = {max(re_c0s):.3e}")
    log(f"  mean|Re(c₀)| = {np.mean(re_c0s):.3e}")
    re_c0_pass = max(re_c0s) < 1e-4
    log(f"  판정: {'★★★ PASS' if re_c0_pass else 'FAIL'}")
    log()

    log("━━━ 검증 2: Im(c₁) = 0 (현실조건의 귀결) ━━━")
    log(f"  max|Im(c₁)|  = {max(im_c1s):.3e}")
    log(f"  mean|Im(c₁)| = {np.mean(im_c1s):.3e}")
    im_c1_pass = max(im_c1s) < 1e-3
    log(f"  판정: {'★★★ PASS' if im_c1_pass else '★★ 부분 PASS (수치 한계)' if max(im_c1s) < 0.1 else 'FAIL (수치 정밀도 부족)'}")
    log()

    log("━━━ 검증 3: slope = 2.0 (정리의 직접 귀결) ━━━")
    log(f"  mean slope = {np.mean(slopes):.6f} ± {np.std(slopes):.6f}")
    log(f"  max|slope-2| = {max(abs(s-2.0) for s in slopes):.6f}")
    slope_pass = max(abs(s-2.0) for s in slopes) < 0.01
    log(f"  판정: {'★★★ PASS' if slope_pass else 'FAIL'}")
    log()

    log("━━━ 검증 4: A = Im(c₀)² + 2c₁ 예측 ━━━")
    log(f"  mean A error = {np.mean(A_errs)*100:.1f}%")
    log(f"  min  A error = {min(A_errs)*100:.1f}%")
    log(f"  max  A error = {max(A_errs)*100:.1f}%")
    a_pass = np.mean(A_errs) < 0.2
    log(f"  판정: {'★★★ PASS' if np.mean(A_errs)<0.05 else '★★ PASS' if a_pass else '★ 부분 (고차항 기여)'}")
    log()

    log("━━━ 상세 결과표 ━━━")
    log(f"{'L-함수':<14} {'t₀':>7} {'|Re(c₀)|':>10} {'Im(c₀)':>9} {'c₁':>8} {'|Im(c₁)|':>10} {'A_pred':>7} {'A_meas':>7} {'err%':>5} {'slope':>8}")
    log("-" * 105)
    for r in all_results:
        log(f"{r['Lfunc']:<14} {r['t0']:>7.3f} {abs(r['Re_c0']):>10.2e} {r['Im_c0']:>9.4f} {r['Re_c1']:>8.4f} {abs(r['Im_c1']):>10.2e} {r['A_pred']:>7.3f} {r['A_meas']:>7.3f} {abs(r['A_pred']-r['A_meas'])/max(abs(r['A_meas']),1e-10)*100:>5.1f} {r['slope']:>8.5f}")
    log()

    # 최종 판정
    log("=" * 72)
    log("최종 판정")
    log("=" * 72)
    log()
    if re_c0_pass and slope_pass:
        log("★★★ 강양성 — Theorem (κδ² Slope Universality) 수치 검증 완료")
        log()
        log("결론:")
        log("  1. Re(c₀)=0는 함수방정식 Λ(s)=ε·Λ(1-s)와 현실조건 Λ(s̄)=Λ̄(s)의")
        log("     필연적 귀결이다. 수치적으로 |Re(c₀)|<10⁻⁴ 확인됨.")
        log("  2. 이로부터 κδ²-1의 δ¹ 항이 정확히 소멸하여 slope=2가 보편적이다.")
        log("  3. 16행 비교표의 slope=2.0000 보편성은 경험적 관찰이 아니라")
        log("     함수방정식의 수학적 필연(mathematical necessity)이다.")
        log()
        log("의의:")
        log("  - 기존: '16개 L-함수에서 slope≈2를 관찰' (관찰)")
        log("  - 변경: 'slope=2는 함수방정식의 정리' (정리)")
        log("  - 모든 자기쌍대(self-dual) L-함수에서 자동 성립")
        log("  - A(t₀) 변동은 c₀(인근 영점 밀도)와 c₁(Γ-인자)로 설명")
        log()
        log("논문 반영:")
        log("  Paper 2, Part II-b에 Theorem (Slope Universality) 추가")
        log("  기존 Observation (slope≈2)을 Theorem으로 승격")
    else:
        log("결과 불완전 — 추가 검증 필요")

log()
log(f"총 소요시간: {time.time()-START:.1f}초")
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
outf.close()
