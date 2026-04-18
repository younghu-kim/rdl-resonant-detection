#!/usr/bin/env python3
"""
#78v2 — GL(4) sym³(Δ) PARI lfun 기반 4성질 검증
핵심 수정: t_VEC float (not closure) + N=25000 → FE=-11

알고리즘:
1. PARI mfcoefs로 tau(n) 빠르게 계산
2. sym3 analytic 계수: a_an(n) = a_arith(n) / n^{33/2}  (n^16.5)
3. t_VEC로 PARI lfuncreate → gammaV=[-1,0,0,1], N_cond=144, eps=1
4. lfuncheckfeq → FE
5. lfunzeros → PARI 공인 영점
6. κ_near: ξ'/ξ = (1/2)log(N) + Σ ψ((s+μⱼ)/2)/2 + L'(s)/L(s)
7. 모노드로미: 폐곡선 적분
8. σ-유일성: κ(0.5)/κ(0.45)
"""
import sys, os, time, math, cmath
sys.path.insert(0, os.path.dirname(__file__))

import cypari2
import mpmath
import numpy as np

# ─── 설정 ────────────────────────────────────────────────────────────────────
N_COEFF    = 25000     # Dirichlet 계수 수
GAMMAVRV   = [-1.0, 0.0, 0.0, 1.0]   # gammaV (analytic, Σ=0)
N_COND     = 144       # conductor = 12²
EPS_ROOT   = 1         # root number
DELTA_K    = 0.001     # κ 측정 시 σ-이동 (실수)
MONO_R     = 0.3       # 모노드로미 반지름
MONO_STEPS = 128       # 폐곡선 단계 수
T_ZERO_MAX = 50.0      # 영점 탐색 범위 최대

RESULT_FILE = os.path.join(os.path.dirname(__file__), '..', 'results', 'gl4_sym3_pari_78v2.txt')

t_START = time.time()

lines = []
def log(s=""):
    print(s, flush=True)
    lines.append(s)

def elapsed():
    return f"{time.time()-t_START:.1f}초"

log("=" * 72)
log("결과 #78v2 — GL(4) sym³(Δ) PARI lfun (t_VEC float, N=25000)")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"핵심: t_VEC float 계수 + gammaV=[-1,0,0,1] → FE=-11 달성")
log(f"N_COEFF={N_COEFF}, gammaV={GAMMAVRV}, N_cond={N_COND}")

# ─── PARI 초기화 ────────────────────────────────────────────────────────────
pari = cypari2.Pari()
pari.allocatemem(512 * 1024 * 1024)
log(f"PARI {pari.version()}")

# ─── [1] tau(n) via mfcoefs ──────────────────────────────────────────────────
log(f"\n[1] τ(n) 계산 (mfcoefs, n=1..{N_COEFF})")
pari('mf12 = mfinit([1,12],1)')
pari('df12 = mfeigenbasis(mf12)[1]')
pari(f'tau_c = mfcoefs(df12, {N_COEFF})')
tau = [0] * (N_COEFF + 1)
for n in range(1, N_COEFF + 1):
    tau[n] = int(pari(f'tau_c[{n+1}]'))
log(f"  완료 ({elapsed()})")
log(f"  τ(2)={tau[2]}, τ(3)={tau[3]}, τ(5)={tau[5]}")
assert tau[2] == -24, f"τ(2) 검증 실패: {tau[2]}"

# ─── [2] sym3 analytic 계수 계산 ─────────────────────────────────────────────
log(f"\n[2] sym³ analytic 계수 계산 (n=1..{N_COEFF})")

# 소수 체
sieve = bytearray([1]) * (N_COEFF + 1)
sieve[0] = sieve[1] = 0
for i in range(2, int(N_COEFF**0.5) + 1):
    if sieve[i]:
        sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
primes = [i for i in range(2, N_COEFF + 1) if sieve[i]]
log(f"  소수 {len(primes)}개 (≤{N_COEFF})")

# 최소 소인수 체
spf = [0] * (N_COEFF + 1)
for p in primes:
    for n in range(p, N_COEFF + 1, p):
        if spf[n] == 0:
            spf[n] = p

# sym3 L-인수 계수: p^k별 계산
# α,β = roots of x²-τ(p)x+p^{11}=0 (arithmetic)
# sym3 roots: α³, α²β, αβ², β³
# a_arith(p^k) = [T^k in 1/prod(1-r_j T)]
sym3_ppow = {}  # (p,k) -> arith coeff (integer)

def sym3_local_factor(p, maxk):
    ap = tau[p]
    disc = ap * ap - 4 * p**11
    sq = cmath.sqrt(disc)
    al = (ap + sq) / 2
    be = (ap - sq) / 2
    r = [al**3, al**2 * be, al * be**2, be**3]
    c = [0 + 0j] * (maxk + 1)
    c[0] = 1.0
    for ri in r:
        nc = [0 + 0j] * (maxk + 1)
        nc[0] = c[0]
        for j in range(1, maxk + 1):
            nc[j] = c[j] + ri * nc[j - 1]
        c = nc
    return [round(c[k].real) for k in range(maxk + 1)]

for p in primes:
    maxk = int(math.log(N_COEFF) / math.log(p))
    if maxk == 0:
        continue
    lf = sym3_local_factor(p, maxk)
    for k in range(1, maxk + 1):
        if p**k <= N_COEFF:
            sym3_ppow[(p, k)] = lf[k]

# 모든 n에 대해 sym3 analytic 계수 계산 (multiplicative)
sym3_an = np.zeros(N_COEFF + 1)
sym3_an[1] = 1.0

for n in range(2, N_COEFF + 1):
    m = n
    fac = {}
    while m > 1:
        p = spf[m]
        while m % p == 0:
            fac[p] = fac.get(p, 0) + 1
            m //= p
    val = 1
    for p, k in fac.items():
        val *= sym3_ppow.get((p, k), 0)
    # analytic normalization: divide by n^{33/2} = n^{16.5}
    sym3_an[n] = val / (n ** 16.5)

log(f"  완료 ({elapsed()})")
t2 = tau[2] / (2 ** 5.5)
log(f"  a_an(2)={sym3_an[2]:.8f}, 기대={t2**3-2*t2:.8f}")
log(f"  a_an(3)={sym3_an[3]:.8f}")
assert abs(sym3_an[2] - (t2**3 - 2*t2)) < 1e-6, "a_an(2) 검증 실패"

# ─── [3] PARI L-함수 구성 ─────────────────────────────────────────────────────
log(f"\n[3] PARI lfuncreate (t_VEC float, gammaV=[-1,0,0,1])")
fvec = '[' + ','.join(f'{x:.15e}' for x in sym3_an[1:N_COEFF+1]) + ']'
pari(f'sym3fv = {fvec}')
pari(f'sym3Lraw = [sym3fv, 0, [-1,0,0,1], {N_COND}, {EPS_ROOT}, [], []]')
pari('sym3Lobj = lfuncreate(sym3Lraw)')
log(f"  lfuncreate OK ({elapsed()})")

# ─── [4] 함수방정식 검증 ──────────────────────────────────────────────────────
log(f"\n[4] 함수방정식 검증 (lfuncheckfeq)")
fe_raw = pari('lfuncheckfeq(sym3Lobj)')
fe_val = int(str(fe_raw))
log(f"  FE = {fe_val} 자릿수 (기대 ≤ -8, 즉 |FE| ≥ 8)")
log(f"  ({elapsed()})")
fe_pass = fe_val <= -8

if fe_pass:
    log(f"  FE 판정: ✅ 통과 ({abs(fe_val)}자릿수)")
else:
    log(f"  FE 판정: ❌ 실패 (score={fe_val})")

# ─── [5] 영점 발견 (numpy vectorized Λ(1/2+it) 스캔) ───────────────────────
log(f"\n[5] 영점 발견 (numpy Λ(1/2+it) 부호 변환, t∈[0,{T_ZERO_MAX}])")
# PARI lfunzeros = poles=[] 버그. numpy로 Dirichlet 합 vectorize.
# Λ(1/2+it) = N^{(1/2+it)/2} * Π Γ_R(1/2+it+μⱼ) * Σ a_an(n)/n^{1/2+it}

# 사전 계산
N_ZERO = min(N_COEFF, 5000)  # Dirichlet 항 수
ns = np.arange(1, N_ZERO + 1, dtype=np.float64)
an_v = sym3_an[1:N_ZERO + 1].copy()   # a(1)..a(N_ZERO)
log_ns = np.log(ns)
an_over_sqrtn = an_v / np.sqrt(ns)    # a(n)/n^{1/2}
log(f"  N_ZERO={N_ZERO} ({elapsed()})")

def lambda_fast(t):
    """Λ(1/2+it) via numpy (float64). 자기대칭 → 실수부 사용."""
    # L(1/2+it) = Σ a(n)/n^{1/2+it} = Σ (a(n)/√n) * exp(-it*log(n))
    L_val = np.dot(an_over_sqrtn, np.exp(-1j * t * log_ns))
    # Gamma factors: Π Γ_R(s+μⱼ) with s=1/2+it
    s = 0.5 + 1j * t
    gam = complex(N_COND) ** (s / 2)
    for mu in GAMMAVRV:
        sm = s + mu
        # Γ_R(sm) = π^{-sm/2} * Γ(sm/2)
        try:
            gam *= (math.pi ** (-sm / 2)) * complex(mpmath.gamma(sm / 2))
        except Exception:
            return None
    return gam * L_val

# Z-함수 스캔 (부호 변환 탐지)
DT = 0.015
T_START = 0.1
t_grid = np.arange(T_START, T_ZERO_MAX + DT, DT)
n_grid = len(t_grid)
log(f"  Z 스캔 ({n_grid}점, dt={DT})...")

zeros_list = []
Z_prev = None
t_prev = None

for i, t in enumerate(t_grid):
    Lv = lambda_fast(t)
    if Lv is None:
        Z_prev = None; continue
    # 자기대칭 L-함수: Im(Λ) ≈ 0, 실수부가 Z-함수
    Zv = Lv.real
    if Z_prev is not None and Z_prev * Zv < 0:
        # 이분법 세밀화
        ta, Za = t_prev, Z_prev
        tb, Zb = t, Zv
        for _ in range(30):
            tm = (ta + tb) / 2
            Lm = lambda_fast(tm)
            if Lm is None: break
            Zm = Lm.real
            if Za * Zm < 0:
                tb, Zb = tm, Zm
            else:
                ta, Za = tm, Zm
        zeros_list.append((ta + tb) / 2)
    Z_prev = Zv
    t_prev = t
    if (i % 500) == 0:
        log(f"  스캔 {i+1}/{n_grid} t={t:.2f} ({elapsed()})")

zeros_list = sorted(zeros_list)
n_zeros = len(zeros_list)
log(f"  영점 {n_zeros}개 발견 ({elapsed()})")

if zeros_list:
    t1 = zeros_list[0]
    log(f"  t₁ = {t1:.5f} (기대 0.3239)")
    log(f"  처음 10개: {[f'{z:.4f}' for z in zeros_list[:10]]}")
else:
    log("  ⚠️ 영점 없음!")
    t1 = 999.0

zeros_pass = (n_zeros >= 30) and (t1 < 1.0)
log(f"  영점 판정: {'✅ 통과' if zeros_pass else '❌ 실패'} (≥30개, t₁<1.0)")

# ─── L-함수 평가 함수 (numpy vectorized L'/L) ────────────────────────────────
# numpy 사전 계산 (κ 측정용, 더 많은 항 사용)
N_KAPPA = min(N_COEFF, 5000)
ns_k = np.arange(1, N_KAPPA + 1, dtype=np.float64)
log_ns_k = np.log(ns_k)
an_k = sym3_an[1:N_KAPPA + 1].copy()

def L_and_deriv(sigma, t):
    """
    Returns (L(s), L'(s)) where s = sigma + it
    L(s) = Σ a(n)/n^s
    L'(s) = -Σ a(n)*log(n)/n^s
    """
    n_pow = ns_k ** (-sigma)
    phase = np.exp(-1j * t * log_ns_k)
    weights = an_k * n_pow * phase
    L_val = np.sum(weights)
    Ld_val = -np.sum(weights * log_ns_k)
    return L_val, Ld_val

def log_xi_deriv(s_real, s_imag):
    """
    ξ'/ξ(s) = (1/2)log(N) + Σⱼ ψ((s+μⱼ)/2)/2 + L'(s)/L(s)
    numpy L'/L + mpmath digamma
    """
    s = complex(s_real, s_imag)
    # (1/2)log(N)
    part_N = 0.5 * math.log(N_COND)

    # Σⱼ ψ((s+μⱼ)/2)/2
    part_gamma = 0.0 + 0j
    for mu in GAMMAVRV:
        arg = (s + mu) / 2
        try:
            psi_val = complex(mpmath.digamma(arg))
            part_gamma += psi_val / 2
        except Exception as e:
            print(f"WARNING: digamma({arg}): {e}")

    # L'(s)/L(s) via numpy
    try:
        Ls, Ld = L_and_deriv(s_real, s_imag)
        if abs(Ls) < 1e-50:
            return None
        part_L = Ld / Ls
    except Exception as e:
        print(f"WARNING: L_and_deriv: {e}")
        return None

    return part_N + part_gamma + part_L

def kappa_at(s_real, s_imag):
    """κ(s) = |ξ'/ξ(s)|²"""
    xi_d = log_xi_deriv(s_real, s_imag)
    if xi_d is None:
        return None
    return abs(xi_d) ** 2

# ─── [6] κ_near 측정 ──────────────────────────────────────────────────────────
log(f"\n[6] κ_near 측정 (δ={DELTA_K}, 실수 σ-이동)")
log(f"  대상: {min(n_zeros, 30)}개 영점")

kappa_results = []
A_results = []
kappa_delta2 = []
sigma_ratios = []

# κ(σ=0.45+it₀) for σ-uniqueness
SIGMA_SHIFT = 0.05  # σ = 0.5 - 0.05 = 0.45 vs σ = 0.5

for idx, t0 in enumerate(zeros_list[:30]):
    s_on = 0.5 + DELTA_K  # σ = 0.5+δ on critical line region
    kap = kappa_at(s_on, t0)
    if kap is None:
        log(f"  [{idx+1:02d}] t₀={t0:.5f}: κ 계산 실패")
        continue

    A_val = kap * (DELTA_K ** 2) - 1.0
    kd2 = kap * (DELTA_K ** 2)

    # σ-uniqueness: κ at σ=0.5 vs σ=0.45
    kap_off = kappa_at(0.5 - SIGMA_SHIFT, t0)
    ratio = kap / kap_off if (kap_off and kap_off > 0) else None

    kappa_results.append(kap)
    A_results.append(A_val)
    kappa_delta2.append(kd2)
    if ratio:
        sigma_ratios.append(ratio)

    ratio_str = f"{ratio:.1f}" if ratio else "N/A"
    log(f"  [{idx+1:02d}] t₀={t0:.5f}: κ={kap:.4f}, A={A_val:.4f}, κδ²={kd2:.6f}, ratio={ratio_str}")

if A_results:
    mean_A = np.mean(A_results)
    std_A = np.std(A_results)
    cv_A = std_A / abs(mean_A) * 100 if abs(mean_A) > 0 else float('inf')
    mean_kd2 = np.mean(kappa_delta2)
    log(f"\n  성공: {len(A_results)}개")
    log(f"  mean(A) = {mean_A:.4f}")
    log(f"  std(A)  = {std_A:.4f}")
    log(f"  CV(A)   = {cv_A:.2f}% (< 5% 기대)")
    log(f"  mean(κδ²) = {mean_kd2:.6f} ([0.99, 1.15] 기대)")
    kappa_pass = (cv_A < 5.0) and (0.99 <= mean_kd2 <= 1.15)
    log(f"  κ_near 판정: {'✅ 통과' if kappa_pass else '⚠️ 검토 필요'}")
else:
    log("  κ 계산 전부 실패")
    kappa_pass = False
    mean_A = 0.0
    cv_A = 999.0

# ─── [7] 모노드로미 (폐곡선 적분) ─────────────────────────────────────────────
log(f"\n[7] 모노드로미 (반지름={MONO_R}, {MONO_STEPS}단계)")

def lambda_at(sigma, t):
    """Λ(σ+it) = N^{s/2} * Π Γ_R(s+μⱼ) * L(s) via numpy"""
    s = complex(sigma, t)
    # L(s)
    n_pow = ns_k ** (-sigma)
    phase = np.exp(-1j * t * log_ns_k)
    L_val = np.dot(an_k * n_pow, phase)
    # Gamma factors
    gam = complex(N_COND) ** (s / 2)
    for mu in GAMMAVRV:
        sm = s + mu
        try:
            gam *= (math.pi ** (-sm / 2)) * complex(mpmath.gamma(sm / 2))
        except Exception:
            return None
    return gam * L_val

def monodromy_sym3(t0, radius=MONO_R, steps=MONO_STEPS):
    """폐곡선 arg(Λ) 누적 — numpy Λ 사용"""
    thetas = [2 * math.pi * k / steps for k in range(steps + 1)]
    arg_prev = None
    total_arg = 0.0
    ok_count = 0
    for theta in thetas:
        s_r = 0.5 + radius * math.cos(theta)
        s_i = t0 + radius * math.sin(theta)
        lval = lambda_at(s_r, s_i)
        if lval is None or abs(lval) < 1e-200:
            continue
        arg_cur = math.atan2(lval.imag, lval.real)
        if arg_prev is not None:
            diff = arg_cur - arg_prev
            while diff > math.pi: diff -= 2 * math.pi
            while diff < -math.pi: diff += 2 * math.pi
            total_arg += diff
        arg_prev = arg_cur
        ok_count += 1
    return total_arg / math.pi  # π 단위

mono_results = []
mono_ok = 0
for idx, t0 in enumerate(zeros_list[:20]):
    mono = monodromy_sym3(t0)
    is_ok = abs(abs(mono) - 2.0) < 0.15
    if is_ok:
        mono_ok += 1
    mono_results.append(mono)
    mark = "✅" if is_ok else "⚠️"
    log(f"  [{idx+1:02d}] t₀={t0:.5f}: mono/π = {mono:.4f} {mark}")

if mono_results:
    mean_mono = np.mean(np.abs(mono_results))
    frac_ok = sum(1 for m in mono_results if abs(abs(m)-2.0)<0.15) / len(mono_results)
    log(f"\n  성공: {len(mono_results)}개")
    log(f"  mean(|mono/π|) = {mean_mono:.4f}")
    log(f"  |mono/π-2|<0.15: {frac_ok*100:.0f}%")
    mono_pass = (frac_ok >= 0.8)
    log(f"  모노드로미 판정: {'✅ 통과' if mono_pass else '⚠️'}")
else:
    mono_pass = False
    mean_mono = 0.0
    frac_ok = 0.0

# ─── [8] σ-유일성 요약 ────────────────────────────────────────────────────────
log(f"\n[8] σ-유일성 (κ(σ=0.5)/κ(σ=0.45))")
if sigma_ratios:
    mean_ratio = np.mean(sigma_ratios)
    frac_ratio = sum(1 for r in sigma_ratios if r > 10) / len(sigma_ratios)
    log(f"  mean(ratio) = {mean_ratio:.2f} (>10 기대)")
    log(f"  ratio>10: {frac_ratio*100:.0f}%")
    sigma_pass = (frac_ratio >= 0.9) and (mean_ratio > 10)
    log(f"  σ-유일성 판정: {'✅ 통과' if sigma_pass else '❌'}")
else:
    sigma_pass = False
    mean_ratio = 0.0

# ─── [9] κ_near(d) 비교 — B-12 ───────────────────────────────────────────────
log(f"\n[9] κ_near(d) 단조증가 비교 — B-12")
d_table = [
    ("d=1 (ζ)",    1.2727),
    ("d=2 (GL2)",  3.9300),
    ("d=3 (GL3)", 12.7900),
    ("d=4 (sym³Δ)", mean_A),
]
for label, a_val in d_table:
    log(f"  {label:<20}: A = {a_val:.4f}")
b12_pass = mean_A > 12.79
log(f"  B-12 판정: {'✅ d=4 > d=3' if b12_pass else '❌ 단조증가 위반'}")

# ─── 최종 요약 ────────────────────────────────────────────────────────────────
log("\n" + "=" * 72)
log("성공 기준 총정리")
log("=" * 72)
checks = [
    ("FE ≤ -8자릿수", fe_pass, f"{abs(fe_val)}자릿수" if fe_pass else f"score={fe_val}"),
    ("영점 ≥30개, t₁<1.0", zeros_pass, f"{n_zeros}개, t₁={t1:.4f}"),
    ("κ_near CV<5%", kappa_pass, f"CV={cv_A:.1f}%"),
    ("mono/π≈2.0", mono_pass, f"mean={mean_mono:.3f}, {frac_ok*100:.0f}%OK"),
    ("σ-유일성 ratio>10", sigma_pass, f"mean={mean_ratio:.1f}"),
    ("κ(d=4)>κ(d=3)", b12_pass, f"A={mean_A:.4f} vs 12.79"),
]
n_pass = sum(1 for _, p, _ in checks if p)
for name, passed, detail in checks:
    mark = "✅" if passed else "❌"
    log(f"  {mark} {name:<30} → {detail}")
log(f"\n  통과: {n_pass}/{len(checks)}")

log("\n" + "-" * 72)
log("수치 요약")
log("-" * 72)
log(f"  FE: {fe_val}자릿수")
log(f"  영점: {n_zeros}개, t₁={t1:.6f}")
log(f"  mean(A): {mean_A:.4f} (d=4)")
log(f"  CV(A): {cv_A:.2f}%")
log(f"  mean(mono/π): {mean_mono:.4f}")
log(f"  mean(ratio): {mean_ratio:.2f}")

log("\n" + "-" * 72)
log("경계 갱신")
log("-" * 72)
if b12_pass:
    log(f"  B-12: ★★ 확립 — d=1,2,3,4 모두 확인")
if sigma_pass:
    log(f"  B-05: ✅ d=4 확인")
log(f"  총 소요: {time.time()-t_START:.1f}초")
log(f"  완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 72)
log(f"결과: {RESULT_FILE}")

# 파일 저장
os.makedirs(os.path.dirname(RESULT_FILE), exist_ok=True)
with open(RESULT_FILE, 'w') as f:
    f.write('\n'.join(lines) + '\n')

print(f"\n[저장 완료] {RESULT_FILE}", flush=True)
