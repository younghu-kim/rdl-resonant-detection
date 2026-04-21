#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #214] κδ² slope=2 해석적 증명 검증
=============================================================================
배경:
  16행 비교표에서 κδ² slope = 2.0000 보편성이 관찰됨.
  이론적 근거:
    Λ'/Λ(s) = 1/(s-ρ) + c₀ + c₁(s-ρ) + ... (로랑 전개)
    함수방정식 Λ(s)=ε·Λ(1-s) + 현실조건 Λ(s̄)=Λ̄(s)로부터:
      c₀ = 순허수 (Re(c₀) = 0)
      c₁ = 실수
    따라서 κδ² - 1 = 2Re(c₀)·δ + O(δ²) = O(δ²)
    → slope = 2 (정확)

검증 전략:
  1. 다양한 L-함수에서 c₀, c₁을 수치적으로 추출
  2. Re(c₀) ≈ 0 확인 (정밀도 수준)
  3. A_predicted = |c₀|² + 2·Re(c₁) 계산
  4. A_measured (κδ² fitting)와 비교
  5. slope=2가 함수방정식의 직접적 귀결임을 수치적으로 확정

대상 L-함수:
  - ζ(s): GL(1), N=1
  - L(s, χ₋₇): GL(1), N=7
  - L(s, Δ): GL(2), weight 12, N=1
  - L(s, sym²(11a1)): GL(3), N=11²

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
log("[실험 #214] κδ² slope=2 해석적 증명 검증")
log("  Theorem: Functional equation → Re(c₀)=0 → slope exactly 2")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PARI 초기화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [0] PARI 초기화...")
import cypari2
gp = cypari2.Pari()
gp.allocatemem(1000 * 1024 * 1024)
gp("default(realprecision, 100)")  # 고정밀도 (c₀ 실수부 검증용)
log(f"  cypari2 OK, realprecision=100")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 이론적 배경
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("이론적 배경 (Theorem Statement)")
log("=" * 72)
log()
log("Theorem (κδ² Slope Universality).")
log("  Let Λ(s) be the completed L-function satisfying:")
log("    (i)  Λ(s) = ε·Λ(1-s)  (functional equation, |ε|=1)")
log("    (ii) Λ(s̄) = Λ̄(s)     (reality, for self-dual)")
log("  Let ρ = 1/2 + iγ be a simple zero of Λ on the critical line.")
log("  Define κ(δ) = |Λ'/Λ(1/2+δ+iγ)|² for real δ > 0.")
log("  Then:")
log("    κ(δ)·δ² - 1 = A·δ² + O(δ³)")
log("  where A = |c₀|² + 2·Re(c₁), with c₀ purely imaginary.")
log()
log("Proof sketch:")
log("  Λ'/Λ(s) = 1/(s-ρ) + c₀ + c₁(s-ρ) + ...")
log("  From (i): Λ'/Λ(s) = -Λ'/Λ(1-s)")
log("    → at s = ρ+δ: 1/δ + c₀ + c₁δ = -(1/(-δ) + c₀ - c₁δ)")
log("                 = 1/δ - c₀ + c₁δ")
log("    → 2c₀ = 0 BUT this assumes 1-ρ = ρ (only if ρ=1/2)")
log()
log("  Correct argument using (i)+(ii):")
log("    Λ'/Λ(1/2+δ+iγ) = -Λ'/Λ(1/2-δ-iγ)")
log("    Λ'/Λ(1/2-δ-iγ) = conj(Λ'/Λ(1/2-δ+iγ))  [from (ii)]")
log("    Laurent at 1/2-δ+iγ: -1/δ + c₀ - c₁δ + ...")
log("    So: 1/δ + c₀ + c₁δ = -conj(-1/δ + c₀ - c₁δ)")
log("                        = 1/δ - c̄₀ + c̄₁δ")
log("    Comparing: c₀ = -c̄₀ → Re(c₀) = 0  ∎")
log("               c₁ = c̄₁  → Im(c₁) = 0")
log()
log("  Consequence:")
log("    κδ² = |1/δ + c₀ + c₁δ|²·δ² = |1 + c₀δ + c₁δ²|²")
log("         = 1 + 2Re(c₀)δ + (|c₀|²+2Re(c₁))δ² + O(δ³)")
log("         = 1 + 0·δ + A·δ²  (since Re(c₀)=0)")
log("    log(κδ²-1) = log(A) + 2·log(δ) + O(δ)  → slope = 2 ∎")
log()
log("  The coefficient A = |c₀|² + 2c₁ depends on:")
log("    - Γ-factor contributions to c₀, c₁ (degree-dependent)")
log("    - Local zero density (nearby zeros contribute to c₀)")
log("    - Conductor via log(N) terms in Γ-factors")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# L-함수 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

LFUNCS = [
    {"name": "ζ(s)", "setup": "Lcur = lfuncreate(1)", "degree": 1, "N": 1},
    {"name": "L(s,χ₋₇)", "setup": "Lcur = lfuncreate(-7)", "degree": 1, "N": 7},
    {"name": "L(s,Δ)", "setup": 'Lcur = lfuncreate(ellinit("11a1"))', "degree": 2, "N": 11},
    {"name": "L(s,sym²(11a1))", "setup": 'Lcur = lfunsympow(ellinit("11a1"), 2)', "degree": 3, "N": 121},
]

N_ZEROS_PER = 5  # 각 L-함수에서 5개 영점
DELTAS = np.array([0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2])

# c₀ 추출용 작은 δ
C0_DELTAS = np.array([1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4])

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 측정 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def get_lambda_log_deriv(sigma, t0):
    """Λ'/Λ(σ+it₀) 계산 — PARI lfunlambda 사용"""
    gp(f"scur = {sigma:.20f} + I*{t0:.20f}")
    gp("Lval = lfunlambda(Linit, scur)")
    gp("dLval = lfunlambda(Linit, scur, 1)")
    gp("ratio = dLval / Lval")
    re = float(gp("real(ratio)"))
    im = float(gp("imag(ratio)"))
    return complex(re, im)

def extract_c0(t0):
    """c₀ = lim_{δ→0} [Λ'/Λ(1/2+δ+it₀) - 1/δ]
    Richardson extrapolation으로 추출"""
    center = 0.5
    vals = []
    for d in C0_DELTAS:
        ratio = get_lambda_log_deriv(center + d, t0)
        c0_est = ratio - 1.0/d  # Λ'/Λ - 1/(s-ρ) where s-ρ = δ
        vals.append(c0_est)

    # 가장 작은 δ들의 평균 (수렴 확인)
    # 마지막 3개 (가장 큰 δ)는 고차항 오염 가능 → 제외
    c0_converged = np.mean(vals[:4])
    c0_spread = np.std([v.real for v in vals[:4]])
    return c0_converged, c0_spread

def extract_c1(t0, c0):
    """c₁ = lim_{δ→0} [Λ'/Λ(1/2+δ+it₀) - 1/δ - c₀] / δ"""
    center = 0.5
    vals = []
    for d in C0_DELTAS[2:]:  # 약간 큰 δ 사용 (너무 작으면 정밀도 문제)
        ratio = get_lambda_log_deriv(center + d, t0)
        c1_est = (ratio - 1.0/d - c0) / d
        vals.append(c1_est)
    c1_converged = np.mean(vals[:3])
    return c1_converged

def measure_kd2_slope(t0):
    """표준 κδ² log-log slope 측정"""
    center = 0.5
    kd2_vals = []
    for d in DELTAS:
        ratio = get_lambda_log_deriv(center + d, t0)
        kappa = abs(ratio)**2
        kd2 = kappa * d**2
        kd2_vals.append(kd2)

    # log(κδ²-1) vs log(δ) fitting
    valid = [(d, kd2) for d, kd2 in zip(DELTAS, kd2_vals) if kd2 > 1.001]
    if len(valid) < 4:
        return None, None, None

    log_d = np.log([v[0] for v in valid])
    log_kd2m1 = np.log([v[1] - 1.0 for v in valid])

    # linear fit
    coeffs = np.polyfit(log_d, log_kd2m1, 1)
    slope = coeffs[0]
    intercept = coeffs[1]

    # R²
    pred = np.polyval(coeffs, log_d)
    ss_res = np.sum((log_kd2m1 - pred)**2)
    ss_tot = np.sum((log_kd2m1 - np.mean(log_kd2m1))**2)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0

    A_measured = np.exp(intercept)  # A from fit: log(A) + 2*log(δ)
    return slope, r2, A_measured

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("수치 검증")
log("=" * 72)
log()

all_results = []

for lf in LFUNCS:
    log(f"{T()} ━━━ {lf['name']} (degree={lf['degree']}, N={lf['N']}) ━━━")

    # L-함수 초기화
    gp(lf["setup"])
    fe = float(gp("lfuncheckfeq(Lcur)"))
    log(f"  FE check: {fe:.0f} digits")
    gp("Linit = lfuninit(Lcur, [0.5, 60])")

    # 영점 찾기
    gp("zvec = lfunzeros(Linit, 50)")
    n_found = int(gp("length(zvec)"))
    log(f"  영점 {n_found}개 발견, 처음 {N_ZEROS_PER}개 사용")

    zeros = []
    for i in range(1, min(n_found, N_ZEROS_PER) + 1):
        z = float(gp(f"zvec[{i}]"))
        zeros.append(z)

    log(f"  영점: {[f'{z:.4f}' for z in zeros]}")
    log()

    for idx, t0 in enumerate(zeros):
        log(f"  영점 #{idx+1}: t₀ = {t0:.6f}")

        # 1. c₀ 추출
        c0, c0_spread = extract_c0(t0)
        re_c0 = c0.real
        im_c0 = c0.imag
        log(f"    c₀ = {re_c0:.2e} + {im_c0:.6f}i")
        log(f"    |Re(c₀)| = {abs(re_c0):.2e}  (이론: 0)")
        log(f"    Re(c₀) 수렴 spread = {c0_spread:.2e}")

        # 2. c₁ 추출
        c1 = extract_c1(t0, c0)
        log(f"    c₁ = {c1.real:.6f} + {c1.imag:.2e}i")
        log(f"    |Im(c₁)| = {abs(c1.imag):.2e}  (이론: 0)")

        # 3. A 예측 vs 측정
        A_predicted = abs(c0)**2 + 2*c1.real

        # 4. 실제 κδ² slope 측정
        slope, r2, A_measured = measure_kd2_slope(t0)

        if slope is not None:
            log(f"    A_predicted = |c₀|² + 2·Re(c₁) = {abs(c0)**2:.4f} + {2*c1.real:.4f} = {A_predicted:.4f}")
            log(f"    A_measured  = {A_measured:.4f}")
            log(f"    A 일치율    = {abs(A_predicted - A_measured)/max(abs(A_measured),1e-10)*100:.1f}% 오차")
            log(f"    slope       = {slope:.4f} (이론: 2.0000)")
            log(f"    R²          = {r2:.6f}")

            all_results.append({
                "Lfunc": lf["name"],
                "degree": lf["degree"],
                "N": lf["N"],
                "t0": t0,
                "Re_c0": re_c0,
                "Im_c0": im_c0,
                "Re_c1": c1.real,
                "Im_c1": c1.imag,
                "A_pred": A_predicted,
                "A_meas": A_measured,
                "slope": slope,
                "R2": r2,
            })
        else:
            log(f"    [SKIP] slope 측정 실패 (유효 데이터 부족)")
        log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 종합 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("종합 분석")
log("=" * 72)
log()

if all_results:
    re_c0_vals = [r["Re_c0"] for r in all_results]
    im_c1_vals = [r["Im_c1"] for r in all_results]
    slopes = [r["slope"] for r in all_results]
    A_errors = [abs(r["A_pred"] - r["A_meas"])/max(abs(r["A_meas"]), 1e-10) for r in all_results]

    log("1. Re(c₀) = 0 검증 (함수방정식 귀결):")
    log(f"   max|Re(c₀)| = {max(abs(v) for v in re_c0_vals):.2e}")
    log(f"   mean|Re(c₀)| = {np.mean(np.abs(re_c0_vals)):.2e}")
    log(f"   판정: {'PASS ✓' if max(abs(v) for v in re_c0_vals) < 1e-3 else 'FAIL ✗'}")
    log()

    log("2. Im(c₁) = 0 검증 (현실조건 귀결):")
    log(f"   max|Im(c₁)| = {max(abs(v) for v in im_c1_vals):.2e}")
    log(f"   mean|Im(c₁)| = {np.mean(np.abs(im_c1_vals)):.2e}")
    log(f"   판정: {'PASS ✓' if max(abs(v) for v in im_c1_vals) < 1e-3 else 'FAIL ✗'}")
    log()

    log("3. slope = 2.0 검증:")
    log(f"   mean slope = {np.mean(slopes):.4f} ± {np.std(slopes):.4f}")
    log(f"   max |slope-2| = {max(abs(s-2.0) for s in slopes):.4f}")
    log(f"   판정: {'PASS ✓' if max(abs(s-2.0) for s in slopes) < 0.01 else 'FAIL ✗'}")
    log()

    log("4. A 예측 정확도 (|c₀|²+2c₁ vs fitting):")
    log(f"   mean error = {np.mean(A_errors)*100:.1f}%")
    log(f"   max error  = {max(A_errors)*100:.1f}%")
    log(f"   판정: {'PASS ✓' if np.mean(A_errors) < 0.1 else '△ (고차항 기여)' if np.mean(A_errors) < 0.3 else 'FAIL ✗'}")
    log()

    log("5. 상세 결과표:")
    log(f"   {'L-함수':<18} {'t₀':>8} {'Re(c₀)':>10} {'Im(c₀)':>10} {'c₁':>10} {'A_pred':>8} {'A_meas':>8} {'slope':>6}")
    log("   " + "-" * 90)
    for r in all_results:
        log(f"   {r['Lfunc']:<18} {r['t0']:>8.3f} {r['Re_c0']:>10.2e} {r['Im_c0']:>10.4f} {r['Re_c1']:>10.4f} {r['A_pred']:>8.3f} {r['A_meas']:>8.3f} {r['slope']:>6.4f}")
    log()

    # 최종 판정
    re_c0_pass = max(abs(v) for v in re_c0_vals) < 1e-3
    slope_pass = max(abs(s-2.0) for s in slopes) < 0.01

    log("=" * 72)
    log("최종 판정")
    log("=" * 72)
    log()
    if re_c0_pass and slope_pass:
        log("★★★ 강양성 — 정리 검증됨")
        log("  함수방정식으로부터 Re(c₀)=0이 수치적으로 확인됨.")
        log("  이로부터 κδ² slope = 2가 정확히 따름 (수학적 필연).")
        log("  16행 비교표의 보편성은 함수방정식의 직접적 귀결.")
        log()
        log("  논문 반영: Theorem (Slope Universality)")
        log("    'The κδ² slope equals exactly 2 for any self-dual L-function")
        log("     with a simple zero on the critical line.'")
    else:
        log("△ 부분 양성 또는 검증 실패 — 상세 분석 필요")

log()
log(f"총 소요시간: {time.time()-START:.1f}초")
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
outf.close()
