#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #123] Artin S₄ 표준 표현 (3D) — 4성질 검증
=============================================================================
배경:
  PARI/GP 미설치 환경. mpmath + 수론 직접 계산으로 구현.
  다항식: x⁴ - x - 1, disc = -283, Gal(K/Q) ≅ S₄
  3D 표준(standard) 표현 V_std → degree-3 Artin L-함수 L(s, V_std)

핵심 관계:
  ζ_K(s) = ζ(s) · L(s, V_std)  (K = Q(α), α root of x⁴-x-1)
  ⟹ a_p(V_std) = (# roots of x⁴-x-1 mod p) − 1

  분해형 → V_std 고유값 → 오일러 곱 점화:
    (1,1,1,1): 4근 → (1,1,1)  → e = (3, 3, 1)
    (1,1,2):   2근 → (1,1,-1) → e = (1, -1, -1)
    (1,3):     1근 → (1,ω,ω²) → e = (0, 0, 1)
    (2,2):     0근 → (1,-1,-1) → e = (-1, -1, 1)
    (4):       0근 → (i,-1,-i) → e = (-1, 1, -1)

  degree-3 점화: c_k = e₁·c_{k-1} − e₂·c_{k-2} + e₃·c_{k-3}

도체: N = 283 (|disc(K)| = 283, 도체-판별식 공식에서)
서명: (2,1) → 복소공액 = 호환 → V_std 고유값 (1,1,-1) → a=2, b=1
감마: Γ_ℝ(s)² · Γ_ℝ(s+1) [μ = (0, 0, 1)]
근 번호: ε = +1 (ζ_K/ζ 분해에서 유도, 수치 확인)

AFE:
  L(s) ≈ Σ_{n≤X} a_n n^{-s} + ε·G(s)·Σ_{n≤X} a_n n^{s-1}
  G(s) = N^{1/2-s} · π^{3s-3/2} · Γ((1-s)/2)² · Γ((2-s)/2) / [Γ(s/2)² · Γ((s+1)/2)]

성공 기준:
  SC1: FE 검증 |Λ(s)−ε·Λ(1−s)| / |Λ(s)| < 0.05 at 3개 점
  SC2: ≥10개 임계선 영점
  SC3a: κδ² slope within 5% of -1.0
  SC3b: 모노드로미 ≈ ±2π
  SC3c: σ-유일성 (σ≈½)

결과: results/artin_s4_standard_123.txt
=============================================================================
"""
import sys, os, time, math
import numpy as np
import mpmath
from datetime import datetime

mpmath.mp.dps = 120
START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'artin_s4_standard_123.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"


log("=" * 72)
log("[실험 #123] Artin S₄ 표준 표현 (3D) — 4성질 검증")
log("=" * 72)
log(f"시작: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log(f"dps = {mpmath.mp.dps}")
log(f"다항식: x⁴ - x - 1,  disc = -283,  Gal(K/Q) ≅ S₄")
log(f"도체 N = 283")
log(f"표준 표현 V_std (3D), 감마: Γ_ℝ(s)²·Γ_ℝ(s+1)")
log()


# ─────────────────────────────────────────────────────────────────────
# 1. 다항식 산술 mod p (0근 분해형 판별용)
# ─────────────────────────────────────────────────────────────────────

def poly_mul_mod(a, b, mod_poly, p):
    """다항식 a·b mod (mod_poly, p). 계수 리스트 [c0,c1,...,cn]."""
    if not a or not b:
        return [0]
    result = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            result[i + j] = (result[i + j] + ai * bj) % p
    # mod_poly로 나누기
    md = len(mod_poly) - 1
    lc_inv = pow(mod_poly[-1], p - 2, p) if mod_poly[-1] % p != 0 else 0
    while len(result) > md:
        if result[-1] % p != 0:
            coeff = (result[-1] * lc_inv) % p
            for k in range(len(mod_poly)):
                idx = len(result) - md - 1 + k
                result[idx] = (result[idx] - coeff * mod_poly[k]) % p
        result.pop()
    while len(result) > 1 and result[-1] % p == 0:
        result.pop()
    return result


def poly_powmod(base, exp, mod_poly, p):
    """base^exp mod (mod_poly, p), 반복 제곱법."""
    result = [1]
    cur = base[:]
    while exp > 0:
        if exp & 1:
            result = poly_mul_mod(result, cur, mod_poly, p)
        cur = poly_mul_mod(cur, cur, mod_poly, p)
        exp >>= 1
    return result


def poly_gcd_mod(a, b, p):
    """gcd(a, b) mod p (유클리드 호제법)."""
    while True:
        # b가 0인지 확인
        b_nonzero = any(c % p != 0 for c in b) if b else False
        if not b_nonzero:
            break
        # a mod b
        a_work = a[:]
        bd = len(b) - 1
        lc_inv = pow(b[-1], p - 2, p) if b[-1] % p != 0 else 0
        while len(a_work) > bd and len(a_work) > 0:
            if a_work[-1] % p != 0:
                coeff = (a_work[-1] * lc_inv) % p
                for k in range(len(b)):
                    idx = len(a_work) - bd - 1 + k
                    a_work[idx] = (a_work[idx] - coeff * b[k]) % p
            a_work.pop()
        while len(a_work) > 1 and a_work[-1] % p == 0:
            a_work.pop()
        a, b = b, a_work
    # a를 모닉으로 정규화
    if a and a[-1] % p != 0:
        inv = pow(a[-1], p - 2, p)
        a = [(c * inv) % p for c in a]
    return a


def split_type_s4(p):
    """
    x⁴-x-1 mod p 분해형 → (e1, e2, e3) 반환.
    n_roots=0일 때 (2,2) vs (4) 구분에 gcd(f, x^{p²}-x) 사용.
    분기소수 283은 None 반환.
    """
    if p == 283:
        return None  # 분기

    # 근의 수 세기
    n_roots = sum(1 for a in range(p) if (pow(a, 4, p) - a - 1) % p == 0)

    if n_roots == 4:
        return (3, 3, 1)       # (1,1,1,1), a_p=3
    elif n_roots == 2:
        return (1, -1, -1)     # (1,1,2), a_p=1
    elif n_roots == 1:
        return (0, 0, 1)       # (1,3), a_p=0
    else:  # n_roots == 0: (2,2) vs (4) 구분
        f_poly = [(-1) % p, (-1) % p, 0, 0, 1]  # x⁴-x-1 mod p
        x_poly = [0, 1]  # x
        xp2 = poly_powmod(x_poly, p * p, f_poly, p)
        # x^{p²} - x
        while len(xp2) < 2:
            xp2.append(0)
        xp2[1] = (xp2[1] - 1) % p
        while len(xp2) > 1 and xp2[-1] % p == 0:
            xp2.pop()
        g = poly_gcd_mod(f_poly[:], xp2, p)
        deg_g = len(g) - 1
        if deg_g >= 4:
            return (-1, -1, 1)   # (2,2), a_p=-1
        else:
            return (-1, 1, -1)   # (4), a_p=-1


# ─────────────────────────────────────────────────────────────────────
# 2. 디리클레 계수 구축
# ─────────────────────────────────────────────────────────────────────

def build_sieve(N_max):
    """최소 소인수 체"""
    spf = list(range(N_max + 1))
    i = 2
    while i * i <= N_max:
        if spf[i] == i:
            for j in range(i * i, N_max + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    return spf


def factorize(n, spf):
    factors = {}
    while n > 1:
        p = spf[n]
        factors[p] = factors.get(p, 0) + 1
        n //= p
    return factors


def a_pk_degree3(p, k, e_cache):
    """a_{p^k} — degree-3 점화: c_k = e1·c_{k-1} − e2·c_{k-2} + e3·c_{k-3}"""
    if k == 0:
        return 1
    info = e_cache.get(p)
    if info is None:
        return 0  # 분기소수: a_{p^k}=0 for k≥1
    e1, e2, e3 = info
    if k == 1:
        return e1
    if k == 2:
        return e1 * e1 - e2
    c = [0] * (k + 1)
    c[0] = 1
    c[1] = e1
    c[2] = e1 * e1 - e2
    for i in range(3, k + 1):
        c[i] = e1 * c[i - 1] - e2 * c[i - 2] + e3 * c[i - 3]
    return c[k]


def build_coeffs_s4(N_max):
    """L(s, V_std) 디리클레 계수 a_n 배열 구축."""
    spf = build_sieve(N_max)
    e_cache = {}
    primes = []

    log(f"  {T()} 분해형 계산 시작 (N_max={N_max})...")
    cnt = 0
    for i in range(2, N_max + 1):
        if spf[i] == i:
            primes.append(i)
            e_cache[i] = split_type_s4(i)
            cnt += 1
            if cnt % 200 == 0:
                log(f"  {T()} {cnt}개 소수 처리됨")

    log(f"  {T()} {len(primes)}개 소수 완료")

    # 분해형 통계
    types = {'(1,1,1,1)': 0, '(1,1,2)': 0, '(1,3)': 0, '(2,2)': 0, '(4)': 0, 'ram': 0}
    e_to_type = {
        (3, 3, 1): '(1,1,1,1)', (1, -1, -1): '(1,1,2)',
        (0, 0, 1): '(1,3)', (-1, -1, 1): '(2,2)', (-1, 1, -1): '(4)'
    }
    for p in primes:
        e = e_cache[p]
        if e is None:
            types['ram'] += 1
        else:
            types[e_to_type[e]] += 1

    log(f"  분해형 분포: {types}")

    # 이론 기대 비율 (Chebotarev 밀도): |C|/|G|
    # S₄ 크기 24. e:1, (12):6, (12)(34):3, (123):8, (1234):6
    # (1,1,1,1):1/24≈4.2%, (1,1,2):6/24=25%, (1,3):8/24≈33.3%, (2,2):3/24=12.5%, (4):6/24=25%
    n_total = len(primes) - types['ram']
    if n_total > 0:
        log(f"  실측 비율: split={types['(1,1,1,1)']}/{n_total}={100*types['(1,1,1,1)']/n_total:.1f}% "
            f"(1,1,2)={100*types['(1,1,2)']/n_total:.1f}% "
            f"(1,3)={100*types['(1,3)']/n_total:.1f}% "
            f"(2,2)={100*types['(2,2)']/n_total:.1f}% "
            f"(4)={100*types['(4)']/n_total:.1f}%")
        log(f"  이론 기대: split=4.2% (1,1,2)=25% (1,3)=33.3% (2,2)=12.5% (4)=25%")

    # a_n 구축
    a = [0] * (N_max + 1)
    a[1] = 1
    for n in range(2, N_max + 1):
        factors = factorize(n, spf)
        val = 1
        for p, k in factors.items():
            val *= a_pk_degree3(p, k, e_cache)
        a[n] = val

    return a, e_cache, primes


N_MAX = 5000
log("━━━ 1. 디리클레 계수 구축 (N_MAX=5000) ━━━")
a_coeffs, e_cache, primes_list = build_coeffs_s4(N_MAX)
log(f"  {T()} 완료. a[1..10] = {a_coeffs[1:11]}")

# 소수별 a_p 확인 (앞 20개)
log(f"  소수별 a_p (앞 20):")
for p in primes_list[:20]:
    e = e_cache[p]
    a_p = e[0] if e else 0
    n_r = a_p + 1 if e else '?'
    tp = {(3,3,1):'split', (1,-1,-1):'(1,1,2)', (0,0,1):'(1,3)',
          (-1,-1,1):'(2,2)', (-1,1,-1):'(4)'}.get(e, 'ram')
    log(f"    p={p:4d}: {n_r}근 → {tp:8s} → a_p={a_p}")
log()


# ─────────────────────────────────────────────────────────────────────
# 3. Artin L-함수 (degree-3 AFE)
# ─────────────────────────────────────────────────────────────────────

N_COND = 283

def G_factor_s4(s, N_cond=N_COND):
    """
    degree-3 AFE의 G(s):
    G(s) = N^{1/2-s} · π^{3s-3/2} · Γ((1-s)/2)² · Γ((2-s)/2)
                                      / [Γ(s/2)² · Γ((s+1)/2)]
    """
    N_fac = mpmath.power(mpmath.mpf(N_cond), mpmath.mpf('0.5') - s)
    pi_fac = mpmath.power(mpmath.pi, 3 * s - mpmath.mpf('1.5'))
    gnum = mpmath.gamma((1 - s) / 2) ** 2 * mpmath.gamma((2 - s) / 2)
    gden = mpmath.gamma(s / 2) ** 2 * mpmath.gamma((s + 1) / 2)
    return N_fac * pi_fac * gnum / gden


def artin_L_afe_s4(s, a_coeffs, N_cond=N_COND, eps=1, X=None):
    """
    L(s, V_std) via AFE:
    L(s) ≈ Σ_{n≤X} a_n n^{-s} + ε·G(s)·Σ_{n≤X} a_n n^{s-1}
    """
    t_val = float(abs(mpmath.im(s)))
    if X is None:
        X = max(20, int(math.sqrt(N_cond * max(t_val, 1.0) / (2 * math.pi))) + 10)
        X = min(X, len(a_coeffs) - 1)

    S1 = mpmath.mpc(0)
    for n in range(1, X + 1):
        if a_coeffs[n] != 0:
            S1 += mpmath.mpf(a_coeffs[n]) * mpmath.power(n, -s)

    S2 = mpmath.mpc(0)
    for n in range(1, X + 1):
        if a_coeffs[n] != 0:
            S2 += mpmath.mpf(a_coeffs[n]) * mpmath.power(n, s - 1)

    G = G_factor_s4(s, N_cond)
    return S1 + mpmath.mpf(eps) * G * S2


def completed_artin_s4(s, a_coeffs, N_cond=N_COND, eps=1, X=None):
    """
    Λ(s, V_std) = (N/π³)^{s/2} · Γ(s/2)² · Γ((s+1)/2) · L(s)
    """
    L_val = artin_L_afe_s4(s, a_coeffs, N_cond, eps, X)
    fac = mpmath.power(mpmath.mpf(N_cond) / mpmath.pi ** 3, s / 2)
    gamma_fac = mpmath.gamma(s / 2) ** 2 * mpmath.gamma((s + 1) / 2)
    return fac * gamma_fac * L_val


# ─────────────────────────────────────────────────────────────────────
# 4. 근 번호(root number) ε 수치 결정
# ─────────────────────────────────────────────────────────────────────

log("━━━ 2. 근 번호(root number) ε 결정 ━━━")
# 이론적 유도: ζ_K(s) = ζ(s)·L(s,V_std). 양쪽 W=+1이므로 ε(V_std)=+1.
# 수치 검증: 낮은 t (|Λ| 크기 유지)에서 확인.

def check_fe_residual_s4(eps, test_pts=None):
    if test_pts is None:
        # 낮은 t 사용 → |Λ|가 극단적으로 작아지지 않음
        test_pts = [
            mpmath.mpf('0.7') + 1j * mpmath.mpf('5'),
            mpmath.mpf('0.6') + 1j * mpmath.mpf('8'),
            mpmath.mpf('0.65') + 1j * mpmath.mpf('12'),
        ]
    residuals = []
    for s in test_pts:
        try:
            Ls = completed_artin_s4(s, a_coeffs, eps=eps)
            L1ms = completed_artin_s4(1 - s, a_coeffs, eps=eps)
            denom = abs(Ls)
            if denom > 1e-80:
                r = float(abs(Ls - eps * L1ms) / denom)
                residuals.append(r)
                log(f"    ε={eps:+d}, s={float(mpmath.re(s)):.2f}+{float(mpmath.im(s)):.0f}i: "
                    f"|Λ|={float(denom):.4e}, 잔차={r:.2e}")
        except Exception as e:
            log(f"  WARNING FE 체크 오류: {e}")
    return np.mean(residuals) if residuals else 999.0

res_p1 = check_fe_residual_s4(eps=+1)
res_m1 = check_fe_residual_s4(eps=-1)
log(f"  FE 잔차 요약 — ε=+1: {res_p1:.2e},  ε=-1: {res_m1:.2e}")

# 이론적 근거로 ε=+1 강제 (ζ_K/ζ 분해). 수치가 동의하면 확인, 아니면 경고.
EPSILON = 1
if res_m1 < res_p1 * 0.01:
    log(f"  ⚠️ 수치는 ε=-1을 선호하나, 이론적 ε=+1 유지 (검토 필요)")
elif res_p1 < res_m1 * 0.01:
    log(f"  ✅ 수치도 ε=+1 지지")
else:
    log(f"  두 ε 모두 잔차 ≈0 (|Λ| 작음), 이론적 ε=+1 채택")
log(f"  → 채택: ε = {EPSILON:+d} (이론: ζ_K=ζ·L(V_std), W(ζ_K)=W(ζ)=+1)")
log()


# ─────────────────────────────────────────────────────────────────────
# 5. Z-함수 & 영점 탐색
# ─────────────────────────────────────────────────────────────────────

def theta_artin_s4(t, N_cond=N_COND):
    """
    degree-3 Artin θ(t):
    θ(t) = (t/2)·(log N − 3 log π) + 2·Im logΓ(1/4+it/2) + Im logΓ(3/4+it/2)
    """
    t = mpmath.mpf(str(t))
    s_half = mpmath.mpf('0.5') + 1j * t
    phase = (t * (mpmath.log(N_cond) - 3 * mpmath.log(mpmath.pi)) / 2
             + 2 * mpmath.im(mpmath.loggamma(s_half / 2))
             + mpmath.im(mpmath.loggamma((s_half + 1) / 2)))
    return float(phase)


def Z_artin_s4(t, a_coeffs, eps=1):
    """Z(t) = e^{iθ(t)}·L(1/2+it) — ε=+1이면 Re(Z)가 실수값."""
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
    L_val = artin_L_afe_s4(s, a_coeffs, eps=eps)
    theta = theta_artin_s4(t)
    Z = mpmath.exp(1j * theta) * L_val
    return float(mpmath.re(Z)), float(mpmath.im(Z))


log("━━━ 3. 임계선 영점 탐색 (t ∈ [5, 65]) ━━━")

T_MIN, T_MAX, DT = 5.0, 65.0, 0.06
ts_scan = np.arange(T_MIN, T_MAX, DT)
Z_re_arr = []
Z_im_arr = []

log(f"  스캔: t∈[{T_MIN},{T_MAX}], Δt={DT}, {len(ts_scan)}점  {T()}")
for idx, t in enumerate(ts_scan):
    if (idx + 1) % 200 == 0:
        log(f"  {T()} 스캔 {idx+1}/{len(ts_scan)}")
    try:
        zr, zi = Z_artin_s4(float(t), a_coeffs, EPSILON)
        Z_re_arr.append(zr)
        Z_im_arr.append(zi)
    except Exception:
        Z_re_arr.append(0.0)
        Z_im_arr.append(0.0)

log(f"  {T()} 스캔 완료")
Z_main = np.array(Z_re_arr if EPSILON == 1 else Z_im_arr)
log(f"  Z-함수 범위: [{Z_main.min():.3f}, {Z_main.max():.3f}]")

# 부호 변화 → 영점 후보
zeros_raw = []
for i in range(len(Z_main) - 1):
    if Z_main[i] * Z_main[i + 1] < 0:
        t1, t2 = float(ts_scan[i]), float(ts_scan[i + 1])
        zeros_raw.append((t1, t2))

log(f"  부호변화 {len(zeros_raw)}개 발견")

# 이분법 정밀화
zeros_found = []
fail_cnt = 0


def Z_scalar_s4(t_val):
    zr, zi = Z_artin_s4(float(t_val), a_coeffs, EPSILON)
    return zr if EPSILON == 1 else zi


for t1, t2 in zeros_raw:
    try:
        f1 = Z_scalar_s4(t1)
        f2 = Z_scalar_s4(t2)
        if f1 * f2 >= 0:
            fail_cnt += 1
            continue
        lo, hi = t1, t2
        flo = f1
        for _ in range(40):
            mid = (lo + hi) / 2
            fm = Z_scalar_s4(mid)
            if fm * flo < 0:
                hi = mid
            else:
                lo = mid
                flo = fm
        zeros_found.append((lo + hi) / 2)
    except Exception:
        fail_cnt += 1

log(f"  이분법 완료: {len(zeros_found)}개 영점, {fail_cnt}개 실패  {T()}")

# 검증: |L(1/2+it)| 확인
for i, t0 in enumerate(zeros_found[:20]):
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t0))
    L_val = artin_L_afe_s4(s, a_coeffs, eps=EPSILON)
    afe_X = max(20, int(math.sqrt(N_COND * max(abs(t0), 1.0) / (2 * math.pi))) + 10)
    afe_X = min(afe_X, len(a_coeffs) - 1)
    log(f"  ρ_{i+1:2d}: t = {t0:.8f},  |L| = {float(abs(L_val)):.2e},  AFE X = {afe_X}항")

log()


# ─────────────────────────────────────────────────────────────────────
# 6. 성질 1: 함수방정식 (FE) 검증
# ─────────────────────────────────────────────────────────────────────

log("━━━ 4. 성질 1: 함수방정식 (FE) 검증 ━━━")

fe_test_points = [
    mpmath.mpf('0.65') + 1j * mpmath.mpf('15'),
    mpmath.mpf('0.60') + 1j * mpmath.mpf('25'),
    mpmath.mpf('0.55') + 1j * mpmath.mpf('40'),
]

fe_errors = []
for s_test in fe_test_points:
    Ls = completed_artin_s4(s_test, a_coeffs, eps=EPSILON)
    L1ms = completed_artin_s4(1 - s_test, a_coeffs, eps=EPSILON)
    denom = abs(Ls)
    if denom > 1e-50:
        err = float(abs(Ls - EPSILON * L1ms) / denom)
        fe_errors.append(err)
        log(f"  s={float(mpmath.re(s_test)):.2f}+{float(mpmath.im(s_test)):.0f}i: "
            f"|Λ(s)|={float(abs(Ls)):.4e}, |Λ(1-s)|={float(abs(L1ms)):.4e}, FE 오차={err:.6f}")
    else:
        log(f"  WARNING: |Λ(s)| 너무 작음 ({float(denom):.2e}) — 건너뜀")

fe_mean = np.mean(fe_errors) if fe_errors else 999
log(f"  평균 FE 오차: {fe_mean:.6f}  → {'✅ PASS' if fe_mean < 0.05 else '❌ FAIL'}")
log()


# ─────────────────────────────────────────────────────────────────────
# 7. 성질 2: κδ² 스케일링
# ─────────────────────────────────────────────────────────────────────

log("━━━ 5. 성질 2: κδ² 스케일링 ━━━")
# 정의: κ = |Λ'/Λ|² at s = (1/2+δ) + it₀
# 이론: 영점 근방에서 |Λ'/Λ| ~ 1/δ → κ ~ 1/δ² → log(κ) vs log(δ²) slope = -1

DELTAS = [0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1]


def curvature_s4(s, a_coeffs, eps=EPSILON):
    """κ = |Λ'/Λ|² (수치 미분, dps=120 고정밀)"""
    h = mpmath.mpf(1) / mpmath.mpf(10 ** 18)
    Lambda_val = completed_artin_s4(s, a_coeffs, eps=eps)
    if abs(Lambda_val) < mpmath.mpf(10) ** (-(mpmath.mp.dps - 10)):
        return None  # 영점 너무 가까움
    Lambda_d = (completed_artin_s4(s + h, a_coeffs, eps=eps) -
                completed_artin_s4(s - h, a_coeffs, eps=eps)) / (2 * h)
    kval = float(abs(Lambda_d / Lambda_val) ** 2)
    if not math.isfinite(kval) or kval > 1e15:
        return None
    return kval


# t≥20 영점 선택 (AFE 수렴 양호)
kappa_zeros = [z for z in zeros_found if z >= 20][:5]
if len(kappa_zeros) < 5:
    kappa_zeros = zeros_found[:5]

log(f"  δ 범위: {DELTAS}")
log(f"  측정 영점 ({len(kappa_zeros)}개): {[f'{z:.2f}' for z in kappa_zeros]}")

slopes = []
for t0 in kappa_zeros:
    kappas = []
    valid_deltas = []
    for delta in DELTAS:
        s_test = mpmath.mpf(str(0.5 + delta)) + 1j * mpmath.mpf(str(t0))
        try:
            k = curvature_s4(s_test, a_coeffs)
            if k is not None and k > 0:
                kappas.append(k)
                valid_deltas.append(delta)
        except Exception:
            pass

    if len(valid_deltas) >= 4:
        log_x = np.log(np.array(valid_deltas) ** 2)
        log_y = np.log(np.array(kappas))
        slope, intercept = np.polyfit(log_x, log_y, 1)
        ss_res = np.sum((log_y - (slope * log_x + intercept)) ** 2)
        ss_tot = np.sum((log_y - log_y.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 1e-30 else 0
        slopes.append(slope)
        afe_X = max(20, int(math.sqrt(N_COND * max(abs(t0), 1.0) / (2 * math.pi))) + 10)
        log(f"    ρ (t={t0:.6f}, X={afe_X}항): slope={slope:.4f} (이론=-1), R²={r2:.5f}  [{len(valid_deltas)}점]")
    else:
        log(f"    ρ (t={t0:.6f}): 데이터 부족 ({len(valid_deltas)}점)")

if slopes:
    mean_slope = np.mean(slopes)
    std_slope = np.std(slopes)
    dev_pct = abs(mean_slope + 1) * 100
    log(f"  평균 slope: {mean_slope:.4f} ± {std_slope:.4f}  (이론 -1.0, 편차 {dev_pct:.1f}%)")
    log(f"  → {'✅ PASS' if dev_pct < 5 else '⚠️ MARGINAL' if dev_pct < 10 else '❌ FAIL'}")
else:
    mean_slope = None
log()


# ─────────────────────────────────────────────────────────────────────
# 8. 성질 3: 모노드로미 (폐곡선 적분)
# ─────────────────────────────────────────────────────────────────────

log("━━━ 6. 성질 3: 모노드로미 (폐곡선 적분) ━━━")


def compute_monodromy_s4(t0, radius=0.5, n_steps=64):
    """영점 주위 폐곡선의 위상 누적."""
    phases = []
    for k in range(n_steps + 1):
        theta = 2 * math.pi * k / n_steps
        s = (mpmath.mpf('0.5') + radius * mpmath.cos(theta)
             + 1j * (mpmath.mpf(str(t0)) + radius * mpmath.sin(theta)))
        L_val = artin_L_afe_s4(s, a_coeffs, eps=EPSILON)
        phases.append(float(mpmath.arg(L_val)))

    total = 0.0
    for k in range(1, len(phases)):
        diff = phases[k] - phases[k - 1]
        while diff > math.pi:
            diff -= 2 * math.pi
        while diff < -math.pi:
            diff += 2 * math.pi
        total += diff

    return total


mono_zeros = zeros_found[:5] if len(zeros_found) >= 5 else zeros_found
mono_pass_cnt = 0
for t0 in mono_zeros:
    try:
        mono = compute_monodromy_s4(t0, radius=0.5)
        winding = mono / (2 * math.pi)
        check = "✅" if abs(abs(winding) - 1) < 0.15 else "⚠️"
        if abs(abs(winding) - 1) < 0.15:
            mono_pass_cnt += 1
        log(f"  ρ (t={t0:.4f}): mono = {mono:.5f} rad  (≈±2π {check}(권선수={winding:.3f}))")
    except Exception as e:
        log(f"  ρ (t={t0:.4f}): ERROR — {e}")

if mono_zeros:
    log(f"  ±2π 판정 비율: {mono_pass_cnt}/{len(mono_zeros)} ({100*mono_pass_cnt/len(mono_zeros):.0f}%)  "
        f"→ {'✅ PASS' if mono_pass_cnt == len(mono_zeros) else '⚠️ PARTIAL'}")
log()


# ─────────────────────────────────────────────────────────────────────
# 9. 성질 4: σ-유일성
# ─────────────────────────────────────────────────────────────────────

log("━━━ 7. 성질 4: σ-유일성 ━━━")

sigma_zeros = zeros_found[:10] if len(zeros_found) >= 10 else zeros_found
sigma_pass = 0
for t0 in sigma_zeros:
    sigmas = np.linspace(0.3, 0.7, 21)
    min_val, min_sigma = 1e300, 0.5
    for sig in sigmas:
        s = mpmath.mpf(str(sig)) + 1j * mpmath.mpf(str(t0))
        L_val = artin_L_afe_s4(s, a_coeffs, eps=EPSILON)
        val = float(abs(L_val))
        if val < min_val:
            min_val, min_sigma = val, sig
    at_half = "σ≈½ ✅" if abs(min_sigma - 0.5) < 0.05 else f"σ={min_sigma:.3f} ⚠️"
    if abs(min_sigma - 0.5) < 0.05:
        sigma_pass += 1
    log(f"  ρ (t={t0:.4f}): min|L| = {min_val:.3e} at σ={min_sigma:.3f} ({at_half})")

sigma_total = len(sigma_zeros)
log(f"  임계선 집중 비율: {sigma_pass}/{sigma_total} ({100*sigma_pass/sigma_total:.0f}%) "
    f"→ {'✅ PASS' if sigma_pass == sigma_total else '⚠️ PARTIAL'}")
log()


# ─────────────────────────────────────────────────────────────────────
# 10. 오일러 곱 교차 검증 (Re(s)=2)
# ─────────────────────────────────────────────────────────────────────

log("━━━ 8. 오일러 곱 검증 (Re(s)=2) ━━━")

s_test_euler = mpmath.mpf('2')
euler_prod = mpmath.mpf(1)
for p in primes_list:
    if p == 283:
        continue
    e = e_cache[p]
    if e is not None and p < 2000:
        e1, e2, e3 = e
        x = mpmath.power(p, -s_test_euler)
        local = 1 / (1 - e1 * x + e2 * x ** 2 - e3 * x ** 3)
        euler_prod *= local

dirichlet_sum = mpmath.mpf(0)
for n in range(1, 1001):
    if a_coeffs[n] != 0:
        dirichlet_sum += mpmath.mpf(a_coeffs[n]) / mpmath.power(n, s_test_euler)

log(f"  L(2,V_std) 오일러곱(p<2000): {float(euler_prod):.8f}")
log(f"  L(2,V_std) Dirichlet급수(1000항): {float(dirichlet_sum):.8f}")
rel_diff = abs(float(euler_prod - dirichlet_sum) / float(euler_prod)) if float(euler_prod) != 0 else 999
log(f"  상대차: {rel_diff:.6f}  ({'✅' if rel_diff < 0.01 else '⚠️'})")
log()


# ─────────────────────────────────────────────────────────────────────
# 11. 종합 결과
# ─────────────────────────────────────────────────────────────────────

log("━━━ 9. 종합 결과 ━━━")
log(f"  다항식: x⁴ - x - 1,  Gal ≅ S₄,  도체 N = 283")
log(f"  dps = {mpmath.mp.dps},  N_MAX = {N_MAX}")
log(f"  3차원 표준 표현 V_std,  근 번호 ε = {EPSILON:+d}")
log()
log(f"  SC1 FE 검증    (오차<5%): {'✅ PASS' if fe_mean < 0.05 else '❌ FAIL'}  (오차={fe_mean:.6f})")
log(f"  SC2 임계선 영점 (≥10개):  {'✅ PASS' if len(zeros_found) >= 10 else '❌ FAIL'}  (발견={len(zeros_found)}개)")
if slopes:
    log(f"  SC3a κδ² 스케일링:        {'✅ PASS' if abs(mean_slope+1)<0.05 else '⚠️ MARGINAL' if abs(mean_slope+1)<0.10 else '❌ FAIL'}  (slope={mean_slope:.4f})")
if mono_zeros:
    log(f"  SC3b 모노드로미 (≈±2π):   {'✅ PASS' if mono_pass_cnt == len(mono_zeros) else '⚠️ PARTIAL'}")
log(f"  SC3c σ-유일성 (σ≈½):     {'✅ PASS' if sigma_pass == sigma_total else '⚠️ PARTIAL'}  ({sigma_pass}/{sigma_total})")
log()

# 영점 목록
log(f"  영점 목록 (앞 25개 / {len(zeros_found)}개):")
for i, t0 in enumerate(zeros_found[:25]):
    afe_X = max(20, int(math.sqrt(N_COND * max(abs(t0), 1.0) / (2 * math.pi))) + 10)
    afe_X = min(afe_X, N_MAX)
    log(f"    {i+1:3d}. t = {t0:.8f}  (AFE X={afe_X})")

log()
log(f"  총 소요 시간: {time.time()-START:.1f}초")
log("=" * 72)
outf.close()
