#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 실험 #116 — DH κδ² δ-sweep: Theorem 4 스케일링 수치 검증
=============================================================================
목표: Theorem 4 (κδ² critical-line discrimination)의 핵심 예측 검증
  - On-critical  (σ=1/2): κδ² = 1 + O(δ²) → log|κδ²-1| vs log(δ) 기울기 ≈ 2
  - Off-critical (σ≠1/2): κδ² = 1 + cδ + O(δ²) → 기울기 ≈ 1

DH 함수: f(s) = [(1-iκ)/2]·L(s,χ) + [(1+iκ)/2]·L(s,χ̄)  (χ mod 5)
영점: #115에서 확인된 10 on-critical + 4 off-critical

δ = {0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1} (8개)
각 영점, 각 δ에서 κδ² 측정 → log-log 기울기 추출

결과: results/dh_delta_sweep_116.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

mpmath.mp.dps = 80
START = time.time()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 출력
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'dh_delta_sweep_116.txt')

outf = open(RESULT_FILE, 'w')

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()


log("=" * 72)
log("[실험 #116] DH κδ² δ-sweep — Theorem 4 스케일링 수치 검증")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"mpmath.dps = {mpmath.mp.dps}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. DH 함수 구현 (#115와 동일)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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
    pref = mpmath.power(5 / mpmath.pi, s / 2)
    gam  = mpmath.gamma((s + 1) / 2)
    return pref * gam * dh_func(s)


def Lambda_deriv(s, h=None):
    """Λ'(s) 수치 미분 (중앙차분)"""
    if h is None:
        h = mpmath.mpf(10) ** (-min(30, mpmath.mp.dps - 10))
    s = mpmath.mpc(s)
    return (Lambda_dh(s + h) - Lambda_dh(s - h)) / (2 * h)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. 영점 목록 (#115에서 확인)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 1. 영점 목록 (from #115) ━━━")

ON_ZEROS = [
    (0.5, 5.094160),
    (0.5, 8.939914),
    (0.5, 12.133545),
    (0.5, 14.404003),
    (0.5, 17.130239),
    (0.5, 19.308800),
    (0.5, 22.159708),
    (0.5, 23.345370),
    (0.5, 26.094967),
    (0.5, 27.923799),
]

OFF_ZEROS = [
    (0.80851718, 85.69934849),
    (0.65083008, 114.16334273),
    (0.57435605, 166.47930591),
    (0.72425769, 176.70246124),
]

log(f"On-critical: {len(ON_ZEROS)}개")
log(f"Off-critical: {len(OFF_ZEROS)}개")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. δ-sweep 측정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 2. δ-sweep: κδ² 측정 ━━━")
log()

DELTAS = [0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1]

def measure_kappa_delta2(sigma0, t0, delta):
    """영점 (σ₀, t₀)에서 δ 오프셋으로 κδ² 측정"""
    s = mpmath.mpc(sigma0 + delta, t0)
    Lam = Lambda_dh(s)
    Lam_prime = Lambda_deriv(s)
    if abs(Lam) < mpmath.mpf(10)**(-60):
        return None  # 영점 너무 가까움
    ratio = Lam_prime / Lam
    kappa = float(abs(ratio)**2)
    return kappa * delta**2


# --- On-critical ---
log("[ On-critical zeros ]")
header = f"{'영점':<8} {'t':>10}"
for d in DELTAS:
    header += f"  δ={d:<6}"
log(header)
log("-" * (20 + 10 * len(DELTAS)))

on_data = []  # (zero_idx, delta, kd2)

for i, (sig, t0) in enumerate(ON_ZEROS):
    row = f"ON#{i+1:<3} {t0:>10.3f}"
    for d in DELTAS:
        kd2 = measure_kappa_delta2(sig, t0, d)
        if kd2 is not None:
            on_data.append((i, d, kd2))
            row += f"  {kd2:>8.6f}"
        else:
            row += f"  {'N/A':>8}"
    log(row)

log()

# --- Off-critical ---
log("[ Off-critical zeros ]")
header = f"{'영점':<8} {'σ':>8} {'t':>10}"
for d in DELTAS:
    header += f"  δ={d:<6}"
log(header)
log("-" * (28 + 10 * len(DELTAS)))

off_data = []  # (zero_idx, delta, kd2)

for i, (sig, t0) in enumerate(OFF_ZEROS):
    row = f"OFF#{i+1:<2} {sig:>8.4f} {t0:>10.3f}"
    for d in DELTAS:
        kd2 = measure_kappa_delta2(sig, t0, d)
        if kd2 is not None:
            off_data.append((i, d, kd2))
            row += f"  {kd2:>8.6f}"
        else:
            row += f"  {'N/A':>8}"
    log(row)

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. 스케일링 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 3. 스케일링 분석: log|κδ²-1| vs log(δ) ━━━")
log()
log("Theorem 4 예측:")
log("  On-critical:  κδ² = 1 + c₂δ² + ...  →  log|κδ²-1| ~ 2·log(δ)")
log("  Off-critical: κδ² = 1 + c₁δ + ...   →  log|κδ²-1| ~ 1·log(δ)")
log()

from scipy import stats as scipy_stats

# --- 영점별 기울기 추출 ---
log("[ 영점별 log-log 기울기 ]")
log(f"{'영점':<8} {'유형':<12} {'기울기':>8} {'R²':>8} {'절편':>10} {'판정'}")
log("-" * 60)

on_slopes = []
for zi in range(len(ON_ZEROS)):
    pts = [(d, kd2) for (z, d, kd2) in on_data if z == zi and abs(kd2 - 1.0) > 1e-15]
    if len(pts) < 3:
        continue
    log_d = np.array([np.log(d) for d, kd2 in pts])
    log_dev = np.array([np.log(abs(kd2 - 1.0)) for d, kd2 in pts])
    slope, intercept, r_value, p_value, std_err = scipy_stats.linregress(log_d, log_dev)
    on_slopes.append(slope)
    sig, t0 = ON_ZEROS[zi]
    ok = "✅" if abs(slope - 2.0) < 0.5 else "❌"
    log(f"ON#{zi+1:<3} {'on-critical':<12} {slope:>8.3f} {r_value**2:>8.4f} {intercept:>10.3f} {ok}")

off_slopes = []
for zi in range(len(OFF_ZEROS)):
    pts = [(d, kd2) for (z, d, kd2) in off_data if z == zi and abs(kd2 - 1.0) > 1e-15]
    if len(pts) < 3:
        continue
    log_d = np.array([np.log(d) for d, kd2 in pts])
    log_dev = np.array([np.log(abs(kd2 - 1.0)) for d, kd2 in pts])
    slope, intercept, r_value, p_value, std_err = scipy_stats.linregress(log_d, log_dev)
    off_slopes.append(slope)
    sig, t0 = OFF_ZEROS[zi]
    ok = "✅" if abs(slope - 1.0) < 0.5 else "❌"
    log(f"OFF#{zi+1:<2} {'off-critical':<12} {slope:>8.3f} {r_value**2:>8.4f} {intercept:>10.3f} {ok}")

log()

# --- 앙상블 통계 ---
log("━━━ 4. 앙상블 통계 ━━━")
log()

if on_slopes:
    on_arr = np.array(on_slopes)
    log(f"On-critical  기울기: {on_arr.mean():.3f} ± {on_arr.std():.3f}  "
        f"(이론 예측: 2.0)  N={len(on_arr)}")
else:
    log("On-critical: 데이터 부족")

if off_slopes:
    off_arr = np.array(off_slopes)
    log(f"Off-critical 기울기: {off_arr.mean():.3f} ± {off_arr.std():.3f}  "
        f"(이론 예측: 1.0)  N={len(off_arr)}")
else:
    log("Off-critical: 데이터 부족")

log()

# --- 분리도 (discrimination power) ---
log("━━━ 5. δ별 on/off 분리도 ━━━")
log()
log(f"{'δ':>8} {'⟨κδ²⟩_on':>12} {'std_on':>10} {'⟨κδ²⟩_off':>12} {'std_off':>10} {'분리':>8}")
log("-" * 70)

for d in DELTAS:
    on_vals = [kd2 for (z, dd, kd2) in on_data if abs(dd - d) < 1e-10]
    off_vals = [kd2 for (z, dd, kd2) in off_data if abs(dd - d) < 1e-10]
    if on_vals and off_vals:
        on_mean = np.mean(on_vals)
        on_std  = np.std(on_vals)
        off_mean = np.mean(off_vals)
        off_std  = np.std(off_vals)
        # gap = min(off) - max(on)
        gap = min(off_vals) - max(on_vals)
        sep = "✅" if gap > 0 else "❌"
        log(f"{d:>8.3f} {on_mean:>12.6f} {on_std:>10.6f} {off_mean:>12.6f} {off_std:>10.6f} {sep:>8} gap={gap:.6f}")

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. Theorem 4 세부 검증: c₁ 추출 (off-critical)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 6. Off-critical c₁ 추출: κδ² ≈ 1 + c₁·δ ━━━")
log()
log("이론: off-critical에서 κδ² - 1 ≈ c₁·δ (선형). c₁은 |σ-1/2|에 의존.")
log()
log(f"{'영점':<8} {'|σ-1/2|':>10} {'c₁ (fit)':>10} {'R² (linear)':>12}")
log("-" * 50)

for zi in range(len(OFF_ZEROS)):
    pts = [(d, kd2) for (z, d, kd2) in off_data if z == zi]
    if len(pts) < 3:
        continue
    deltas_arr = np.array([d for d, _ in pts])
    dev_arr = np.array([kd2 - 1.0 for _, kd2 in pts])
    # 선형 피팅: (κδ²-1) = c₁·δ + c₂·δ²
    # 단순 선형: (κδ²-1)/δ = c₁ + c₂·δ
    ratio_arr = dev_arr / deltas_arr
    slope_c, intercept_c, r_val, _, _ = scipy_stats.linregress(deltas_arr, ratio_arr)
    c1 = intercept_c
    sig, t0 = OFF_ZEROS[zi]
    dist = abs(sig - 0.5)
    log(f"OFF#{zi+1:<2} {dist:>10.4f} {c1:>10.4f} {r_val**2:>12.4f}")

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 7. 최종 판정 ━━━")
log()

on_ok = len(on_slopes) > 0 and abs(np.mean(on_slopes) - 2.0) < 0.5
off_ok = len(off_slopes) > 0 and abs(np.mean(off_slopes) - 1.0) < 0.5

if on_ok and off_ok:
    verdict = "★★★ 강양성"
    log(f"판정: {verdict}")
    log("On-critical 기울기 ≈ 2 (O(δ²) 확인) + Off-critical 기울기 ≈ 1 (O(δ) 확인)")
    log("Theorem 4의 스케일링 예측 수치적으로 검증됨.")
elif on_ok or off_ok:
    verdict = "★★ 부분 양성"
    log(f"판정: {verdict}")
    if on_ok:
        log("On-critical 기울기 ≈ 2 확인, off-critical 예측과 불일치.")
    else:
        log("Off-critical 기울기 ≈ 1 확인, on-critical 예측과 불일치.")
else:
    verdict = "음성"
    log(f"판정: {verdict}")
    log("두 스케일링 모두 예측과 불일치.")

log()
elapsed = time.time() - START
log(f"소요 시간: {elapsed:.1f}초")
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
outf.close()
print(f"\n결과 저장: {RESULT_FILE}")
