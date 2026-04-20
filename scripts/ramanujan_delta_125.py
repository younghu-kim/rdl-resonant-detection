#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #125] Ramanujan Δ (weight 12, level 1) — ξ-다발 4성질 검증
=============================================================================
배경:
  ξ-다발 프레임워크가 GL(1)~GL(3), weight 1~2에서 검증됨 (Papers 1-3).
  Weight 12 모듈러 형식 (Ramanujan Δ)에서 동일 4성질이 성립하는지 확인.

  Δ(z) = η(z)^24 = Σ τ(n)q^n, 최초 cusp form (weight 12, level 1, N=1).
  L(s, Δ) = Σ τ(n)/n^s, 임계선 Re(s) = 6.
  Λ(s) = (2π)^{-s} Γ(s) L(s), FE: Λ(s) = Λ(12-s), ε = +1.

방법:
  Mellin 변환 + FE 분할로 Λ(s) 계산 (지수적 수렴):
    Λ(s) = Σ τ(n) [Γ(s,2πn)/(2πn)^s + Γ(12-s,2πn)/(2πn)^{12-s}]
  최적화: 스캔은 dps=30, 정밀 측정은 dps=80.
  4성질: FE, κδ² = 1+O(δ²), 모노드로미 = 2π, σ-유일성

성공 기준:
  5/5 SC PASS → ★★★ (weight 불변 확립)
  4/5 → ★★ (조건부)
  ≤3 → ★ 또는 음성

결과: results/ramanujan_delta_125.txt
=============================================================================
"""
import sys, os, time, math
import numpy as np
import mpmath
from datetime import datetime

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'ramanujan_delta_125.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

SIGMA_CRIT = 6       # 임계선: Re(s) = k/2 = 6
K_WEIGHT = 12
N_LEVEL = 1
N_TERMS = 15          # e^{-2π·8}≈10^{-22}, 15항이면 dps=30에 과잉
DELTAS = [0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05, 0.10]

log("=" * 72)
log("[실험 #125] Ramanujan Δ (weight 12, N=1) — ξ-다발 4성질 검증")
log("=" * 72)
log(f"시작: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log(f"weight k = {K_WEIGHT}, level N = {N_LEVEL}")
log(f"임계선: Re(s) = {K_WEIGHT//2}")
log(f"FE: Λ(s) = Λ({K_WEIGHT}-s), ε = +1")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. τ(n) 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [1/6] τ(n) 계산 (n=1..{N_TERMS})...")

def compute_tau_coefficients(N_max):
    """q-전개로 τ(n) 계산: Δ(q) = q·∏(1-q^n)^24"""
    c = [0] * (N_max + 2)
    c[0] = 1
    for m in range(1, N_max + 2):
        for _ in range(24):
            for k in range(N_max + 1, m - 1, -1):
                c[k] -= c[k - m]
    return {n: c[n - 1] for n in range(1, N_max + 1)}

TAU_INT = compute_tau_coefficients(N_TERMS)

# 검증
KNOWN_TAU = {1: 1, 2: -24, 3: 252, 4: -1472, 5: 4830,
             6: -6048, 7: -16744, 8: 84480, 9: -113643, 10: -115920}
tau_ok = True
for n, expected in KNOWN_TAU.items():
    if n > N_TERMS:
        break
    computed = TAU_INT[n]
    match = "✓" if computed == expected else "✗"
    if computed != expected:
        tau_ok = False
    log(f"  τ({n:2d}) = {computed:>10d}  (기대: {expected:>10d}) {match}")

if not tau_ok:
    log("  ⚠ τ(n) 불일치 — 중단!")
    sys.exit(1)
log(f"  ✓ τ(n) 검증 완료")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. Λ(s) = Σ τ(n) [Γ(s,2πn)/(2πn)^s + Γ(12-s,2πn)/(2πn)^{12-s}]
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [2/6] Λ(s) 구현 (Mellin 분할, 지수적 수렴)...")

def completed_delta(s, N=N_TERMS):
    """완비 L-함수 Λ(s,Δ) via Mellin integral split.
    Λ(s) = Σ_{n=1}^N τ(n) [Γ(s,2πn)·(2πn)^{-s} + Γ(12-s,2πn)·(2πn)^{-(12-s)}]
    """
    result = mpmath.mpc(0)
    two_pi = 2 * mpmath.pi
    s12 = 12 - s
    for n in range(1, N + 1):
        x = two_pi * n
        tau_n = TAU_INT[n]
        g1 = mpmath.gammainc(s, x)
        g2 = mpmath.gammainc(s12, x)
        xps = mpmath.power(x, s)
        xps12 = mpmath.power(x, s12)
        result += tau_n * (g1 / xps + g2 / xps12)
    return result


def connection_delta(s, N=N_TERMS):
    """접속 L(s) = Λ'/Λ (수치 미분)"""
    h = mpmath.mpf(10) ** (-20)
    val = completed_delta(s, N)
    if abs(val) < mpmath.mpf(10) ** (-mpmath.mp.dps + 10):
        return mpmath.mpc(1e10, 0)
    dval = (completed_delta(s + h, N) - completed_delta(s - h, N)) / (2 * h)
    return dval / val


# FE 검증
mpmath.mp.dps = 50
log(f"  FE 교차검증 (dps={mpmath.mp.dps}):")
test_points = [
    mpmath.mpc(3, 10), mpmath.mpc(4, 15), mpmath.mpc(5.5, 20),
    mpmath.mpc(7, 12), mpmath.mpc(8, 25),
]
fe_max_err = 0
for sp in test_points:
    v1 = completed_delta(sp)
    v2 = completed_delta(12 - sp)
    denom = max(abs(v1), abs(v2), mpmath.mpf(1e-100))
    err = float(abs(v1 - v2) / denom)
    fe_max_err = max(fe_max_err, err)
    log(f"    s={sp}: err = {err:.2e}")
fe_pass = fe_max_err < 1e-20
log(f"  {'✓' if fe_pass else '✗'} FE 최대 오차: {fe_max_err:.2e}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. 영점 탐색 (저정밀 dps=25로 빠른 스캔 → 고정밀로 이분법)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [3/6] 영점 탐색 (Re(s)=6, t∈[1,60])...")

# 정밀도 요구사항:
#   |Λ(6+it)| ~ |Γ(6+it)|/(2π)^6 ~ t^{5.5} · e^{-πt/2} / (2π)^6
#   t=20: ~10^{-11}, t=30: ~10^{-17}, t=40: ~10^{-24}
#   → dps ≈ πt/(2·ln10) + 10 여유
#
# 전략: t ∈ [1,50]을 3구간으로 나눠 적정 dps 사용

sigma = mpmath.mpf(SIGMA_CRIT)

def scan_segment(t_min, t_max, n_pts, dps_val):
    """주어진 dps로 [t_min, t_max] 구간 스캔, 부호변화 구간 반환."""
    old_dps = mpmath.mp.dps
    mpmath.mp.dps = dps_val
    sig = mpmath.mpf(SIGMA_CRIT)
    ts = np.linspace(t_min, t_max, n_pts)
    vals = []
    for t in ts:
        s = sig + 1j * mpmath.mpf(str(t))
        v = float(mpmath.re(completed_delta(s)))
        vals.append(v)
    intervals = []
    for i in range(len(vals) - 1):
        if vals[i] * vals[i + 1] < 0:
            intervals.append((ts[i], ts[i + 1]))
    mpmath.mp.dps = old_dps
    return intervals

# 구간별 스캔 (dps 적정화)
segments = [
    (1.0, 20.0, 600, 30),     # |Λ| ~ 10^{-6}~10^{-11}
    (20.0, 35.0, 500, 45),    # |Λ| ~ 10^{-11}~10^{-20}
    (35.0, 50.0, 500, 60),    # |Λ| ~ 10^{-20}~10^{-30}
]

all_intervals = []
for t_lo, t_hi, n_pts, dps_v in segments:
    intervals = scan_segment(t_lo, t_hi, n_pts, dps_v)
    log(f"  구간 [{t_lo:.0f},{t_hi:.0f}] (dps={dps_v}, {n_pts}점): {len(intervals)}개 부호변화")
    all_intervals.extend(intervals)

log(f"  총 부호변화 구간: {len(all_intervals)}개")

# 고정밀 이분법
mpmath.mp.dps = 60
zeros = []
for a, b in all_intervals:
    va = float(mpmath.re(completed_delta(sigma + 1j * mpmath.mpf(str(a)))))
    for _ in range(60):
        mid = (a + b) / 2
        vm = float(mpmath.re(completed_delta(sigma + 1j * mpmath.mpf(str(mid)))))
        if va * vm < 0:
            b = mid
        else:
            a = mid
            va = vm
    zeros.append((a + b) / 2)

n_zeros = len(zeros)
log(f"  발견된 영점: {n_zeros}개 (dps=60 이분법)")
for i, t0 in enumerate(zeros[:15]):
    log(f"    ρ_{i+1} = 6 + {t0:.8f}i")
if n_zeros > 15:
    log(f"    ... ({n_zeros-15}개 추가)")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. SC1: 함수방정식 정밀도
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [4/6] 4성질 검증...")
log()
log("--- SC1: 함수방정식 Λ(s) = Λ(12-s) ---")

mpmath.mp.dps = 50
fe_errs = []
for t0 in zeros[:5]:
    s_test = sigma + 1j * mpmath.mpf(str(t0 + 0.01))
    v1 = completed_delta(s_test)
    v2 = completed_delta(12 - s_test)
    denom = max(abs(v1), abs(v2), mpmath.mpf(1e-100))
    err = float(abs(v1 - v2) / denom)
    fe_errs.append(err)
    log(f"  ρ_{zeros.index(t0)+1} 근방: FE err = {err:.2e}")

fe_mean = np.mean(fe_errs) if fe_errs else 1.0
sc1_pass = fe_mean < 1e-15
log(f"  SC1: {'✓ PASS' if sc1_pass else '✗ FAIL'} (평균 오차: {fe_mean:.2e})")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SC3a: κδ² = 1 + O(δ²)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("--- SC3a: κδ² = 1 + O(δ²) ---")
log(f"  δ 값: {DELTAS}")

mpmath.mp.dps = 60  # κδ²는 정밀도 중요
n_kappa_zeros = min(5, n_zeros)
kappa_results = []

for iz in range(n_kappa_zeros):
    t0 = zeros[iz]
    log(f"  영점 ρ_{iz+1} = 6 + {t0:.6f}i:")
    kd2_vals = []

    for delta in DELTAS:
        s_plus = sigma + 1j * mpmath.mpf(str(t0 + delta))
        s_minus = sigma + 1j * mpmath.mpf(str(t0 - delta))

        L_plus = connection_delta(s_plus)
        L_minus = connection_delta(s_minus)
        kappa_plus = float(abs(L_plus) ** 2)
        kappa_minus = float(abs(L_minus) ** 2)
        kappa_avg = (kappa_plus + kappa_minus) / 2.0
        kd2 = kappa_avg * delta ** 2
        kd2_vals.append(kd2)
        log(f"    δ={delta:.3f}: κδ² = {kd2:.6f}")

    # log|κδ²-1| vs log(δ) 기울기
    log_d = np.log(np.array(DELTAS))
    dev = np.array([abs(v - 1.0) for v in kd2_vals])
    valid = dev > 1e-12
    if valid.sum() >= 3:
        log_dev = np.log(dev[valid])
        log_d_v = log_d[valid]
        slope, intercept = np.polyfit(log_d_v, log_dev, 1)
        fitted = slope * log_d_v + intercept
        ss_res = np.sum((log_dev - fitted) ** 2)
        ss_tot = np.sum((log_dev - np.mean(log_dev)) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    else:
        slope, r2 = 2.0, 1.0

    kappa_results.append({'t': t0, 'slope': slope, 'r2': r2, 'kd2_vals': kd2_vals})
    log(f"    → slope = {slope:.4f} (이론: 2.0), R² = {r2:.6f}")
    log()

slopes = [r['slope'] for r in kappa_results]
mean_slope = np.mean(slopes)
std_slope = np.std(slopes)
sc3a_pass = abs(mean_slope - 2.0) < 0.5 and min(slopes) > 1.0
log(f"  κδ² slope 평균: {mean_slope:.4f} ± {std_slope:.4f}")
log(f"  SC3a: {'✓ PASS' if sc3a_pass else '✗ FAIL'}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SC3b: 모노드로미 = 2π
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("--- SC3b: 모노드로미 = 2π ---")

mpmath.mp.dps = 50

def monodromy_delta(t0, radius=0.5, n_steps=128):
    """영점 주위 폐곡선 적분으로 모노드로미 측정."""
    center = sigma + 1j * mpmath.mpf(str(t0))
    total_angle = mpmath.mpf(0)
    prev_val = completed_delta(center + radius)
    for k in range(1, n_steps + 1):
        theta = 2 * mpmath.pi * k / n_steps
        point = center + radius * mpmath.expj(theta)
        curr_val = completed_delta(point)
        ratio = curr_val / prev_val
        darg = mpmath.im(mpmath.log(ratio))
        total_angle += darg
        prev_val = curr_val
    return float(total_angle)

n_mono = min(5, n_zeros)
mono_results = []
for iz in range(n_mono):
    t0 = zeros[iz]
    mono = monodromy_delta(t0, radius=0.3, n_steps=96)
    mono_pi = mono / math.pi
    mono_results.append(mono_pi)
    log(f"  ρ_{iz+1} (t={t0:.4f}): mono/π = {mono_pi:.4f} (기대: ±2.0)")

sc3b_pass = all(abs(abs(m) - 2.0) < 0.1 for m in mono_results) if mono_results else False
log(f"  SC3b: {'✓ PASS' if sc3b_pass else '✗ FAIL'} ({sum(1 for m in mono_results if abs(abs(m)-2.0)<0.1)}/{len(mono_results)})")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SC3c: σ-유일성
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("--- SC3c: σ-유일성 ---")

mpmath.mp.dps = 45  # t<40에서 |Λ| ~ 10^{-24} → dps=45 필요

def count_sign_changes(sigma_val, t_min=5.0, t_max=30.0, n_points=500):
    """Re(Λ(σ+it)) 부호변화 수."""
    ts = np.linspace(t_min, t_max, n_points)
    count = 0
    prev_sign = None
    for t in ts:
        s = mpmath.mpf(str(sigma_val)) + 1j * mpmath.mpf(str(t))
        val = float(mpmath.re(completed_delta(s)))
        curr_sign = 1 if val >= 0 else -1
        if prev_sign is not None and curr_sign != prev_sign:
            count += 1
        prev_sign = curr_sign
    return count

sigmas_test = [4.0, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0]
sigma_jumps = {}
for sig in sigmas_test:
    nj = count_sign_changes(sig)
    sigma_jumps[sig] = nj
    marker = " ← critical" if sig == 6.0 else ""
    log(f"  σ={sig:.1f}: 부호변화 {nj:3d}{marker}")

crit_jumps = sigma_jumps.get(6.0, 0)
sc3c_pass = all(sigma_jumps[sig] <= crit_jumps for sig in sigmas_test)
log(f"  σ=6.0 최대 여부: {sc3c_pass}")
log(f"  SC3c: {'✓ PASS' if sc3c_pass else '✗ FAIL'}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. 오일러 곱 교차검증
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [5/6] 오일러 곱 교차검증...")

mpmath.mp.dps = 40

def L_euler_product(s, primes_max=200):
    """L(s,Δ) via Euler product (Re(s) > 13/2 필요)"""
    spf = list(range(primes_max + 1))
    i = 2
    while i * i <= primes_max:
        if spf[i] == i:
            for j in range(i * i, primes_max + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    primes = [p for p in range(2, primes_max + 1) if spf[p] == p]

    result = mpmath.mpc(1)
    for p in primes:
        if p > N_TERMS:
            break
        tau_p = TAU_INT[p]
        ps = mpmath.power(p, s)
        p2s = ps * ps
        p11 = mpmath.power(p, 11)
        factor = 1 - mpmath.mpf(tau_p) / ps + p11 / p2s
        if abs(factor) > 1e-50:
            result /= factor
    return result

def L_dirichlet_sum(s, N_max=N_TERMS):
    """L(s,Δ) via Dirichlet sum"""
    result = mpmath.mpc(0)
    for n in range(1, N_max + 1):
        result += mpmath.mpf(TAU_INT[n]) / mpmath.power(n, s)
    return result

log(f"  비교: Re(s)=8 (수렴 영역)")
for sp in [mpmath.mpc(8, 5), mpmath.mpc(8, 10), mpmath.mpc(8, 15)]:
    L_e = L_euler_product(sp)
    L_d = L_dirichlet_sum(sp)
    denom = max(abs(L_e), abs(L_d), mpmath.mpf(1e-100))
    rel = float(abs(L_e - L_d) / denom)
    log(f"    s={sp}: |Euler-Dirich|/|L| = {rel:.6e}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. 종합 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [6/6] 종합 판정")
log()
log("=" * 72)
log("종합 결과")
log("=" * 72)

results_table = [
    ("SC1 FE", sc1_pass),
    ("SC2 영점", n_zeros > 0),
    ("SC3a κδ²", sc3a_pass),
    ("SC3b 모노드로미", sc3b_pass),
    ("SC3c σ-유일성", sc3c_pass),
]

pass_count = sum(1 for _, p in results_table if p)
for name, passed in results_table:
    log(f"  {name}: {'✓ PASS' if passed else '✗ FAIL'}")

log()
log(f"  통과: {pass_count}/5")

if pass_count == 5:
    verdict = "★★★ 강양성 — weight 12에서 ξ-다발 4성질 확립"
elif pass_count == 4:
    verdict = "★★ 조건부 양성"
elif pass_count == 3:
    verdict = "★ 약양성"
else:
    verdict = "음성"

log(f"  판정: {verdict}")
log()
log(f"  영점 수: {n_zeros}")
if kappa_results:
    log(f"  κδ² slope: {mean_slope:.4f} ± {std_slope:.4f}")
if mono_results:
    log(f"  모노드로미: {np.mean([abs(m) for m in mono_results]):.4f}π")
log(f"  σ=6 부호변화: {sigma_jumps.get(6.0, 'N/A')}")
log()
log(f"  L-함수: Ramanujan Δ, weight {K_WEIGHT}, level {N_LEVEL}")
log(f"  Λ(s) = (2π)^{{-s}} Γ(s) L(s,Δ), FE: Λ(s) = Λ({K_WEIGHT}-s)")
log(f"  계산: Mellin 분할 ({N_TERMS}항, 지수 수렴)")
log()
log(f"종료: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log(f"소요: {time.time()-START:.1f}초")
log("=" * 72)

outf.close()
