#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] #79 — GL(4) sym³(Δ) ★ PARI lfun 직접 호출 4성질 검증
=============================================================================
#78 실패 원인: PARI lfun을 호출하지 않고 기존 Mellin-Barnes AFE 코드 재사용.
이 스크립트는 cypari2를 통해 PARI의 lfun 모듈을 직접 호출한다.

★ 핵심 수정 (v79b):
  - lfuninit 제거: lfuninit이 L(s)를 손상시킴 (FE=-1 vs direct FE=-11)
  - Ldata를 PARI 글로벌 변수로 저장: 매 호출마다 재파싱 방지
  - N_COEFF=4000: FE=-11 확인된 설정

PARI 파라미터 (저자 확인):
  gammaV = [-1, 0, 0, 1]
  conductor = 144 = 12²
  root_number = 1 (self-dual)

결과: results/gl4_sym3_pari_79.txt
=============================================================================
"""

import sys, os, time
import numpy as np
from scipy.special import loggamma, digamma as scipy_digamma

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

# ─── 출력 ─────────────────────────────────────────────────────────────────
OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "results", "gl4_sym3_pari_79.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))
def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ─── 파라미터 ─────────────────────────────────────────────────────────────
N_COEFF     = 4000      # Dirichlet 계수 개수 (FE=-11 확인된 설정)
DELTA_KAPPA = 0.001     # κ_near 측정 δ
MONO_R      = 0.3       # 모노드로미 반지름
MONO_N      = 64        # 모노드로미 적분 단계
T_RANGE     = [0, 50]   # 영점 탐색 범위

# ─── PARI 초기화 ──────────────────────────────────────────────────────────
import cypari2
gp = cypari2.Pari()
gp.allocatemem(512 * 1024 * 1024)  # 512MB

log("=" * 72)
log("결과 #79 — GL(4) sym³(Δ) PARI lfun 직접 호출 4성질 검증 (v79b)")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"★ lfuninit 제거 — 직접 lfun(Ldata, s) 호출만 사용")
log(f"gammaV=[-1,0,0,1], conductor=144, root_number=1")
log(f"N_COEFF={N_COEFF}, δ_κ={DELTA_KAPPA}, mono_r={MONO_R}")
log()
flush_file()

# ─── [1] Ramanujan τ(n) 계산 ──────────────────────────────────────────────
log("[1] Ramanujan τ(n) 계산 (PARI ramanujantau)")
t0 = time.time()

# PARI로 τ(n) 직접 계산
tau = [0] * (N_COEFF + 1)  # tau[0] unused
for n in range(1, N_COEFF + 1):
    tau[n] = int(gp.ramanujantau(n))
    if n % 500 == 0:
        log(f"  τ: n={n}/{N_COEFF}")
        flush_file()

log(f"  τ(1)={tau[1]}, τ(2)={tau[2]}, τ(3)={tau[3]}, τ(5)={tau[5]}")
log(f"  ✅ 완료 ({time.time()-t0:.1f}s)")
log()
flush_file()

# ─── [2] sym³ Dirichlet 계수 계산 ────────────────────────────────────────
log("[2] sym³(Δ) Dirichlet 계수 c(n) 계산 (해석적 정규화)")
t0 = time.time()

def sieve_primes(limit):
    is_p = [True] * (limit + 1)
    is_p[0] = is_p[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_p[i]:
            for j in range(i*i, limit+1, i):
                is_p[j] = False
    return [i for i in range(2, limit+1) if is_p[i]]

primes = sieve_primes(N_COEFF)
log(f"  소수 {len(primes)}개 (≤{N_COEFF})")

# 해석적 정규화: t_p = τ(p) / p^{11/2}
# sym³ Euler 곱: e₁ = t³-2t, 4항 재귀
cpk = {}
for p in primes:
    tp = tau[p] / (p ** 5.5)
    e1 = tp**3 - 2.0 * tp
    e2 = (tp**2 - 2.0) * (tp**2 - 1.0)
    cpk[(p, 0)] = 1.0
    cpk[(p, 1)] = e1
    pk = p; k = 1
    while pk * p <= N_COEFF:
        pk *= p; k += 1
        c1 = cpk.get((p, k-1), 0.0)
        c2 = cpk.get((p, k-2), 0.0)
        c3 = cpk.get((p, k-3), 0.0)
        c4 = cpk.get((p, k-4), 0.0)
        cpk[(p, k)] = e1*c1 - e2*c2 + e1*c3 - c4

cn = np.zeros(N_COEFF + 1, dtype=np.float64)
cn[1] = 1.0
for n in range(2, N_COEFF + 1):
    temp = n; result = 1.0
    for p in primes:
        if p * p > temp:
            break
        if temp % p == 0:
            k = 0
            while temp % p == 0:
                k += 1; temp //= p
            result *= cpk.get((p, k), 0.0)
    if temp > 1:
        result *= cpk.get((temp, 1), 0.0)
    cn[n] = result

# 검증
tp2 = tau[2] / (2 ** 5.5)
expected_c2 = tp2**3 - 2*tp2
log(f"  c(1)={cn[1]:.6f}, c(2)={cn[2]:.6f} (기대 {expected_c2:.6f})")
log(f"  max|c_p| (p≤100): {max(abs(cn[p]) for p in primes if p <= 100):.4f} (≤4 기대)")
log(f"  ✅ 완료 ({time.time()-t0:.1f}s)")
log()
flush_file()

# ─── [3] PARI L-함수 구성 ────────────────────────────────────────────────
log("[3] PARI lfuncreate — L(s, sym³Δ) 구성")
t0 = time.time()

# 계수를 PARI 벡터로 전달 → 글로벌 변수에 저장 (매번 재파싱 방지)
cn_str = "[" + ",".join(f"{cn[i]:.15e}" for i in range(1, N_COEFF + 1)) + "]"
gp(f"global(Ldata); Ldata = lfuncreate([{cn_str}, 0, [-1,0,0,1], 144, 1, 0])")
log(f"  Ldata 생성 완료 (PARI 글로벌 변수)")
flush_file()

# FE 검증
try:
    fe_score = float(gp("lfuncheckfeq(Ldata)"))
    log(f"  lfuncheckfeq = {fe_score:.1f} (목표: ≤ -8)")
except Exception as e:
    fe_score = 0.0
    log(f"  ❌ lfuncheckfeq 실패: {e}")

log(f"  ✅ 완료 ({time.time()-t0:.1f}s)")
log()
flush_file()

# ─── [4] 함수방정식 검증 ──────────────────────────────────────────────────
log("=" * 72)
log("[4] ★ 함수방정식 검증 (PARI lfuncheckfeq)")
log("=" * 72)

if fe_score <= -8:
    log(f"  FE = {fe_score:.1f}자리 정밀도")
    log(f"  ✅ PASS (≥8자리)")
elif fe_score <= -3:
    log(f"  FE = {fe_score:.1f}자리 정밀도")
    log(f"  ⚠️ 부분 통과 (3~8자리)")
else:
    log(f"  FE = {fe_score:.1f}자리 정밀도")
    log(f"  ❌ FAIL")

log()
flush_file()

# ─── [5] 영점 탐색 (PARI lfunzeros) ──────────────────────────────────────
log("=" * 72)
log("[5] ★ 영점 탐색 (PARI lfunzeros)")
log("=" * 72)
t0 = time.time()

try:
    zeros_pari = gp(f"lfunzeros(Ldata, {T_RANGE[1]})")
    zeros = [float(z) for z in zeros_pari]
    log(f"  영점 {len(zeros)}개 발견 (t∈[0,{T_RANGE[1]}])")
    log(f"  t₁ = {zeros[0]:.6f}" if zeros else "  영점 없음")
    for i, z in enumerate(zeros[:10]):
        log(f"  영점 #{i+1:2d}: t = {z:.6f}")
    if len(zeros) > 10:
        log(f"  ... (총 {len(zeros)}개)")
    log(f"  ✅ 완료 ({time.time()-t0:.1f}s)")
except Exception as e:
    zeros = []
    log(f"  ❌ lfunzeros 실패: {e}")

log()
flush_file()

# ─── [6] L(s) 계산 함수 (★ lfuninit 없이 직접 호출) ─────────────────────
log("[6] L(s) 직접 호출 테스트 (lfuninit 사용 안 함)")
t0 = time.time()

def pari_lfun_value(s_real, s_imag):
    """PARI로 L(s) 값 계산. ★ Ldata 직접 호출 (lfuninit 미사용)."""
    try:
        if abs(s_imag) < 1e-15:
            s_str = f"{s_real:.15f}"
        else:
            s_str = f"{s_real:.15f} + {s_imag:.15f}*I"
        result = gp(f"lfun(Ldata, {s_str})")
        # PARI 결과 → Python complex
        r = float(gp(f"real({result})"))
        i = float(gp(f"imag({result})"))
        return complex(r, i)
    except Exception:
        return complex('nan')

# L(s) 작동 검증 — 여러 점
test_points = [(0.6, 5.0), (0.51, 5.0), (2.0, 5.0)]
for sr, si in test_points:
    v = pari_lfun_value(sr, si)
    log(f"  L({sr}+{si}i) = {v.real:.6e} + {v.imag:.6e}i  |L|={abs(v):.6e}")
log(f"  ✅ 완료 ({time.time()-t0:.1f}s)")
log()
flush_file()

# ─── [7] ξ'/ξ 계산 함수 ─────────────────────────────────────────────────
# Λ(s) = N^{s/2} · ∏ Γ_R(s + μ_j) · L(s)
# Λ'/Λ = log(N)/2 + γ'/γ + L'/L
# γ'/γ = ∑ [-log(π)/2 + ψ((s+μ)/2)/2]  (ψ = digamma)
# L'/L: PARI L(s) 수치 미분

MU_VEC = [-1.0, 0.0, 0.0, 1.0]
LOG_N_HALF = 0.5 * np.log(144.0)

def gamma_prime_over_gamma(s):
    """γ'/γ(s) 해석적 계산 (digamma 사용)"""
    s = complex(s)
    result = 0.0j
    for mu in MU_VEC:
        sm = s + mu
        result += -0.5 * np.log(np.pi) + 0.5 * scipy_digamma(sm / 2)
    result += LOG_N_HALF
    return result

def log_gamma_factor(s):
    """log(N^{s/2} · ∏ Γ_R(s + μ_j)) 계산"""
    s = complex(s)
    result = 0.0j
    for mu in MU_VEC:
        sm = s + mu
        result += -sm / 2 * np.log(np.pi) + loggamma(sm / 2)
    result += s / 2 * np.log(144.0)
    return result

def Lambda_value(s_real, s_imag):
    """Λ(s) = N^{s/2} · γ(s) · L(s)"""
    s = complex(s_real, s_imag)
    L_val = pari_lfun_value(s_real, s_imag)
    if np.isnan(L_val.real):
        return complex('nan')
    log_gf = log_gamma_factor(s)
    return np.exp(log_gf) * L_val

def xi_prime_over_xi(s_real, s_imag, h=1e-6):
    """ξ'/ξ(s) = γ'/γ + L'/L + log(N)/2, L'/L은 수치 미분.
    ★ h=1e-6 (이전 1e-7보다 안정적, L~10^260에서도 11자릿수 정확)"""
    L_plus = pari_lfun_value(s_real + h, s_imag)
    L_minus = pari_lfun_value(s_real - h, s_imag)
    L_center = pari_lfun_value(s_real, s_imag)
    if abs(L_center) < 1e-300 or np.isnan(L_center.real):
        return complex('nan')
    # L'/L 수치 미분
    dLds = (L_plus - L_minus) / (2 * h)
    LpL = dLds / L_center
    # γ'/γ 해석적 (이미 log(N)/2 포함)
    dgdg = gamma_prime_over_gamma(complex(s_real, s_imag))
    return dgdg + LpL

# ─── [8] κ_near 측정 ─────────────────────────────────────────────────────
log("=" * 72)
log("[8] κ_near 측정 (PARI L(s) 기반, δ=0.001)")
log("=" * 72)
t0_sec = time.time()

# t<25 영점만 사용 (정밀도 높은 영역)
zeros_low = [z for z in zeros if z > 0.5 and z < 25]
log(f"  대상 영점: {len(zeros_low)}개 (t∈[0.5, 25])")

kappa_results = []
A_results = []

for idx, t_zero in enumerate(zeros_low):
    sigma = 0.5 + DELTA_KAPPA
    conn = xi_prime_over_xi(sigma, t_zero)
    if np.isnan(conn.real):
        log(f"  [{idx+1:02d}] t₀={t_zero:.6f}: ❌ 계산 실패")
        continue
    kappa = abs(conn) ** 2
    A = kappa - 1.0 / (DELTA_KAPPA ** 2)
    kappa_results.append(kappa)
    A_results.append(A)
    kd2 = kappa * DELTA_KAPPA**2
    log(f"  [{idx+1:02d}] t₀={t_zero:.6f}: κ={kappa:.2f}, A={A:.4f}, κδ²={kd2:.6f}")
    flush_file()

if A_results:
    A_arr = np.array(A_results)
    mean_A = np.mean(A_arr)
    std_A = np.std(A_arr)
    cv_A = abs(std_A / mean_A * 100) if abs(mean_A) > 1e-10 else float('inf')
    mean_kd2 = np.mean(np.array(kappa_results) * DELTA_KAPPA**2)
    log()
    log(f"  mean(A) = {mean_A:.4f}")
    log(f"  std(A)  = {std_A:.4f}")
    log(f"  CV(A)   = {cv_A:.1f}%")
    log(f"  mean(κδ²) = {mean_kd2:.6f}")
    if cv_A < 5:
        log(f"  ✅ κ_near CV < 5%")
    elif cv_A < 30:
        log(f"  ⚠️ κ_near CV = {cv_A:.1f}% (5~30% 범위)")
    else:
        log(f"  ❌ κ_near CV = {cv_A:.1f}% (> 30%)")
else:
    mean_A = float('nan')
    cv_A = float('inf')
    log("  ❌ κ_near 계산 실패")

log(f"  소요: {time.time()-t0_sec:.1f}s")
log()
flush_file()

# ─── [9] δ 독립성 체크 ───────────────────────────────────────────────────
log("[9] δ 독립성 체크 (3개 영점 × 4개 δ)")
deltas = [0.0001, 0.001, 0.01, 0.1]
check_zeros = zeros_low[:3] if len(zeros_low) >= 3 else zeros_low

for t_zero in check_zeros:
    vals = []
    for d in deltas:
        sigma = 0.5 + d
        conn = xi_prime_over_xi(sigma, t_zero)
        if not np.isnan(conn.real):
            kap = abs(conn) ** 2
            A_val = kap - 1.0 / (d ** 2)
            vals.append(f"δ={d}: A={A_val:.3f}")
        else:
            vals.append(f"δ={d}: FAIL")
    log(f"  t₀={t_zero:.5f}: {', '.join(vals)}")

log()
flush_file()

# ─── [10] 모노드로미 ─────────────────────────────────────────────────────
log("=" * 72)
log("[10] 모노드로미 (PARI L(s) 기반 폐곡선 적분)")
log("=" * 72)
t0_sec = time.time()

# 근접 영점 쌍 감지: 간격 < 2*MONO_R인 영점은 반지름 축소
def get_mono_radius(t_zero, all_zeros, default_r=MONO_R):
    """근접 영점이 있으면 반지름을 영점 간격의 40%로 축소"""
    min_gap = float('inf')
    for z in all_zeros:
        if z != t_zero:
            gap = abs(z - t_zero)
            if gap < min_gap:
                min_gap = gap
    if min_gap < 2 * default_r:
        return min_gap * 0.4
    return default_r

mono_results = []
mono_zeros = zeros_low[:20] if len(zeros_low) > 20 else zeros_low

for idx, t_zero in enumerate(mono_zeros):
    radius = get_mono_radius(t_zero, zeros, MONO_R)
    center_sigma = 0.5
    phase_acc = 0.0
    prev_val = None
    ok = True

    for k in range(MONO_N + 1):
        theta = 2 * np.pi * k / MONO_N
        pt_sigma = center_sigma + radius * np.cos(theta)
        pt_t = t_zero + radius * np.sin(theta)
        val = Lambda_value(pt_sigma, pt_t)
        if np.isnan(val.real) or abs(val) < 1e-300:
            ok = False
            break
        if prev_val is not None:
            phase_acc += np.angle(val / prev_val)
        prev_val = val

    if ok:
        mono_pi = abs(phase_acc) / np.pi
        mono_results.append(mono_pi)
        status = "✅" if abs(mono_pi - 2.0) < 0.15 else "⚠️"
        log(f"  [{idx+1:02d}] t₀={t_zero:.6f}: mono/π={mono_pi:.4f} (r={radius:.3f}) {status}")
    else:
        log(f"  [{idx+1:02d}] t₀={t_zero:.6f}: ❌ 계산 실패")
    flush_file()

if mono_results:
    mono_arr = np.array(mono_results)
    n_pass = sum(1 for m in mono_results if abs(m - 2.0) < 0.15)
    log()
    log(f"  mean(mono/π) = {np.mean(mono_arr):.4f}")
    log(f"  mono/π ≈ 2.0: {n_pass}/{len(mono_results)} ({n_pass/len(mono_results)*100:.0f}%)")
    if n_pass / len(mono_results) >= 0.8:
        log(f"  ✅ 모노드로미 PASS")
    else:
        log(f"  ⚠️ 모노드로미 부분 통과")
else:
    n_pass = 0
    log("  ❌ 모노드로미 전체 실패")

log(f"  소요: {time.time()-t0_sec:.1f}s")
log()
flush_file()

# ─── [11] σ-유일성 ────────────────────────────────────────────────────────
log("=" * 72)
log("[11] σ-유일성 검사 (κ(σ=0.5)/κ(σ=0.45) ratio)")
log("=" * 72)
t0_sec = time.time()

sigma_zeros = zeros_low[:15] if len(zeros_low) > 15 else zeros_low
ratios = []

for idx, t_zero in enumerate(sigma_zeros):
    k_half = None
    k_off = None
    for sigma_test in [0.5, 0.45]:
        delta = sigma_test - 0.5 if sigma_test != 0.5 else DELTA_KAPPA
        if sigma_test == 0.5:
            delta = DELTA_KAPPA
            sigma = 0.5 + delta
        else:
            sigma = sigma_test
            delta = abs(sigma - 0.5)
        conn = xi_prime_over_xi(sigma, t_zero)
        if not np.isnan(conn.real):
            kap = abs(conn) ** 2
            if sigma_test == 0.5:
                k_half = kap
            else:
                k_off = kap

    if k_half is not None and k_off is not None and k_off > 0:
        ratio = k_half / k_off
        ratios.append(ratio)
        status = "✅" if ratio > 10 else "❌"
        log(f"  [{idx+1:02d}] t₀={t_zero:.6f}: ratio={ratio:.1f} {status}")
    else:
        log(f"  [{idx+1:02d}] t₀={t_zero:.6f}: ❌ 계산 실패")
    flush_file()

if ratios:
    mean_ratio = np.mean(ratios)
    n_pass_sigma = sum(1 for r in ratios if r > 10)
    log()
    log(f"  mean(ratio) = {mean_ratio:.1f}")
    log(f"  ratio > 10: {n_pass_sigma}/{len(ratios)} ({n_pass_sigma/len(ratios)*100:.0f}%)")
    if n_pass_sigma / len(ratios) >= 0.8:
        log(f"  ✅ σ-유일성 PASS")
    else:
        log(f"  ❌ σ-유일성 FAIL")
else:
    mean_ratio = 0
    n_pass_sigma = 0
    log("  ❌ σ-유일성 전체 실패")

log(f"  소요: {time.time()-t0_sec:.1f}s")
log()
flush_file()

# ─── [12] κ_near(d) 단조증가 비교 ────────────────────────────────────────
log("=" * 72)
log("[12] κ_near(d) 단조증가 비교 — B-12")
log("=" * 72)

prev_data = {
    1: ("ζ(s)", 1.2727),
    2: ("GL(2) avg", 3.93),
    3: ("GL(3) sym²", 12.79),
}

log(f"  d=1 (ζ)       : A = {prev_data[1][1]:.4f}")
log(f"  d=2 (GL2)     : A = {prev_data[2][1]:.4f}")
log(f"  d=3 (GL3)     : A = {prev_data[3][1]:.4f}")
log(f"  d=4 (sym³Δ)   : A = {mean_A:.4f} ★")

if not np.isnan(mean_A):
    if mean_A > prev_data[3][1]:
        log(f"  ✅ 단조증가 확인: {prev_data[3][1]:.2f} < {mean_A:.2f}")
        b12_pass = True
    else:
        log(f"  ❌ 단조증가 위반: {prev_data[3][1]:.2f} ≥ {mean_A:.2f}")
        b12_pass = False
else:
    log(f"  ❌ mean(A) 계산 불가")
    b12_pass = False

gaps = []
prev_A = 0
for d in [1, 2, 3]:
    gaps.append(prev_data[d][1] - prev_A)
    prev_A = prev_data[d][1]
if not np.isnan(mean_A):
    gaps.append(mean_A - prev_data[3][1])
    log(f"  gap: 1→2={gaps[1]:.2f}, 2→3={gaps[2]:.2f}, 3→4={gaps[3]:.2f}")

log()
flush_file()

# ─── [13] 종합 ────────────────────────────────────────────────────────────
log("=" * 72)
log("종합 판정 — GL(4) sym³(Δ) #79 (PARI lfun 직접 호출)")
log("=" * 72)
log()

results = {}

# P0: FE
fe_pass = fe_score <= -8
results['FE'] = fe_pass
log(f"  [P0] 함수방정식 (PARI lfuncheckfeq)")
log(f"       FE = {fe_score:.1f}자리")
log(f"       {'✅ PASS' if fe_pass else '⚠️' if fe_score <= -3 else '❌ FAIL'}")
log()

# P1: 영점
zeros_pass = len(zeros) >= 30
results['zeros'] = zeros_pass
log(f"  [P1] 영점 (PARI lfunzeros)")
log(f"       {len(zeros)}개 (t∈[0,{T_RANGE[1]}])")
log(f"       t₁ = {zeros[0]:.6f}" if zeros else "       없음")
log(f"       {'✅ PASS' if zeros_pass else '❌ FAIL'}")
log()

# P2: κ_near
kappa_pass = cv_A < 5 if not np.isnan(mean_A) else False
results['kappa'] = kappa_pass
log(f"  [P2] κ_near (δ={DELTA_KAPPA})")
log(f"       mean(A) = {mean_A:.4f}, CV = {cv_A:.1f}%")
log(f"       {'✅ PASS' if kappa_pass else '⚠️ CV 범위 확인' if cv_A < 30 else '❌ FAIL'}")
log()

# P3: 모노드로미
mono_pass = n_pass / len(mono_results) >= 0.8 if mono_results else False
results['mono'] = mono_pass
if mono_results:
    log(f"  [P3] 모노드로미")
    log(f"       mono/π ≈ 2.0: {n_pass}/{len(mono_results)}")
    log(f"       {'✅ PASS' if mono_pass else '⚠️ 부분 통과'}")
else:
    log(f"  [P3] 모노드로미: ❌ FAIL")
log()

# P4: σ-유일성
sigma_pass = n_pass_sigma / len(ratios) >= 0.8 if ratios else False
results['sigma'] = sigma_pass
log(f"  [P4] σ-유일성")
log(f"       mean(ratio) = {mean_ratio:.1f}")
log(f"       {'✅ PASS' if sigma_pass else '❌ FAIL'}")
log()

# P5: B-12
results['b12'] = b12_pass
log(f"  [P5] κ_near(d) 단조증가 (B-12)")
log(f"       {'✅ PASS' if b12_pass else '❌ FAIL'}")
log()

total_pass = sum(results.values())
total = len(results)
log(f"  통과: {total_pass}/{total}")
log()

# 수치 요약
log("─" * 72)
log("수치 요약")
log("─" * 72)
log(f"  FE: {fe_score:.1f}자리 (PARI lfuncheckfeq)")
log(f"  영점: {len(zeros)}개, t₁={zeros[0]:.6f}" if zeros else "  영점: 없음")
log(f"  mean(A): {mean_A:.4f} (d=4)")
log(f"  CV(A): {cv_A:.1f}%")
log(f"  mean(κδ²): {mean_kd2:.6f}" if A_results else "  mean(κδ²): N/A")
log(f"  mean(mono/π): {np.mean(mono_arr):.4f}" if mono_results else "  mean(mono/π): N/A")
log(f"  mean(ratio): {mean_ratio:.1f}")
log(f"  완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 72)
log(f"결과: {OUTFILE}")

flush_file()
print(f"\n결과 저장: {OUTFILE}")
