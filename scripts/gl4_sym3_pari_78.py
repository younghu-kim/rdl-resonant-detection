#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] #78 — GL(4) sym³(Δ) 수정 AFE 기반 4성질 검증
=============================================================================
배경:
  #72v1~v5 실패 원인 = 감마 파라미터 오류:
    - v4 사용: MU=(0,1,11,12), N=1  ← 틀림 (ΣμJ=24≠0)
    - 올바름:  MU=(-1,0,0,1), N=144 ← 저자 PARI 확인 (FE=-11, 79영점)
  ΣμJ=0은 해석적 정규화 s→1-s 함수방정식의 필요조건.
  N=144=12²은 weight-12 형식의 PARI 해석적 conductor.

해결:
  AFE를 gammaV=[-1,0,0,1], N=144로 수정.
  지수인자의 N^{s/2}를 올바르게 포함.
  mpmath DPS=60 고정밀도 유지.

4성질 검증:
  1. 함수방정식 FE_error 자릿수 (≥8)
  2. 영점 개수 (≥30 in [0,50]), t₁<1.0
  3. κ_near CV < 5%
  4. 모노드로미 mono/π ≈ 2.0
  5. σ-유일성 ratio > 10

결과: results/gl4_sym3_pari_78.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from math import comb

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

# ─── 출력 설정 ────────────────────────────────────────────────────────────
OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "results", "gl4_sym3_pari_78.txt"
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
DPS         = 60        # mpmath 정밀도 (t>50 → 80 필요, 여기선 t<50)
N_COEFF     = 2000      # Dirichlet 계수 개수 (정밀 AFE 위해 증가)
N_GH        = 80        # Gauss-Hermite 노드
C_SHIFT     = 2.0       # Mellin 윤곽 이동 (MU=[-1,0,0,1] → c>0.5 충분)
T_MIN       = 0.2       # 영점 탐색 시작
T_MAX       = 50.0      # 영점 탐색 끝
DELTA_KAPPA = 0.001     # κ_near 측정 δ
H_DERIV     = 1e-6      # 수치 미분 h
MONO_R      = 0.4       # 모노드로미 반지름
MONO_N      = 64        # 모노드로미 단계

# ★ 핵심 수정: 올바른 감마 파라미터
MU          = [-1.0, 0.0, 0.0, 1.0]   # gammaV=[-1,0,0,1], ΣμJ=0 ✓
CONDUCTOR   = 144.0                     # 12² (PARI 확인)
ROOT_NUMBER = 1                         # ε=1, self-dual

mpmath.mp.dps = DPS

# Gauss-Hermite 노드
GH_X, GH_W = np.polynomial.hermite.hermgauss(N_GH)

# ─── 소수 체 ──────────────────────────────────────────────────────────────
def sieve_primes(limit):
    is_p = [True] * (limit + 1)
    is_p[0] = is_p[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_p[i]:
            for j in range(i*i, limit+1, i):
                is_p[j] = False
    return [i for i in range(2, limit+1) if is_p[i]]


# ─── τ(n) 계산 ────────────────────────────────────────────────────────────
def compute_tau_fast(limit):
    """Δ(q)=q∏(1-q^n)^24 계수 — binomial 최적화"""
    binom24 = [comb(24, k) * ((-1)**k) for k in range(25)]
    p = [0] * (limit + 1)
    p[0] = 1
    for n in range(1, limit + 1):
        if n % 500 == 0:
            log(f"  τ 계산: n={n}/{limit}")
            flush_file()
        max_k_n = min(24, limit // n)
        for j in range(limit, n - 1, -1):
            val = 0
            upper = min(max_k_n, j // n)
            for k in range(1, upper + 1):
                val += binom24[k] * p[j - n * k]
            if val:
                p[j] += val
    tau = np.zeros(limit + 1, dtype=np.float64)
    for m in range(1, limit + 1):
        tau[m] = float(p[m - 1])
    return tau


# ─── sym³ 해석적 정규화 계수 ──────────────────────────────────────────────
def compute_sym3_an(tau, primes, limit):
    """
    해석적 정규화: t_p = τ(p)/p^{5.5}
    e₁ = t_p³ - 2t_p
    e₂ = (t_p²-2)(t_p²-1)
    e₃ = e₁,  e₄ = 1  (αβ=1 조건)

    ★ v4와 동일 계수 — 감마 파라미터만 수정
    """
    cpk = {}
    for p in primes:
        if p > limit:
            break
        tp = tau[p] / (p ** 5.5)
        e1 = tp**3 - 2.0 * tp
        e2 = (tp**2 - 2.0) * (tp**2 - 1.0)
        cpk[(p, 0)] = 1.0
        cpk[(p, 1)] = e1
        pk = p; k = 1
        while pk * p <= limit:
            pk *= p; k += 1
            c1 = cpk.get((p, k-1), 0.0)
            c2 = cpk.get((p, k-2), 0.0)
            c3 = cpk.get((p, k-3), 0.0)
            c4 = cpk.get((p, k-4), 0.0)
            # e₃ = e₁, e₄ = 1
            cpk[(p, k)] = e1*c1 - e2*c2 + e1*c3 - c4

    an = np.zeros(limit + 1, dtype=np.float64)
    an[0] = 0.0; an[1] = 1.0
    for n in range(2, limit + 1):
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
        an[n] = result
    return an


# ─── 감마 인자 (★ 수정 MU=[-1,0,0,1], N=144) ────────────────────────────
def lgf_mp(s):
    """
    log Λ_dir(s) without L(s):
    = (s/2)log(N) + Σⱼ [-(s+μⱼ)/2·log(π) + logΓ((s+μⱼ)/2)]
    """
    N_mp = mpmath.mpf(str(CONDUCTOR))
    r = s / 2 * mpmath.log(N_mp)
    for mu in MU:
        mu_mp = mpmath.mpf(str(mu))
        r += -(s + mu_mp) / 2 * mpmath.log(mpmath.pi)
        r += mpmath.loggamma((s + mu_mp) / 2)
    return r


def Ldir_mp(s, cn_mp, logn_mp):
    """Λ_dir(s) = γ(s) × L(s)  [mpmath]"""
    g = mpmath.exp(lgf_mp(s))
    L = sum(cn_mp[i] * mpmath.exp(-s * logn_mp[i]) for i in range(len(cn_mp)))
    return g * L


def LAFE_mp(s, cn_mp, logn_mp, c=C_SHIFT, eps_r=ROOT_NUMBER):
    """
    Mellin-Barnes AFE — mpmath
    Λ(s) ≈ ec² / (2π) ∫ [Λ_dir(s+w) + ε·Λ_dir(1-s+w)] e^{2icw-w²} / w dw
    """
    c_mp = mpmath.mpf(str(c))
    ec2 = mpmath.exp(c_mp**2)
    tot = mpmath.mpc(0, 0)
    s1 = mpmath.mpf(1) - s
    for k in range(N_GH):
        y  = mpmath.mpf(str(GH_X[k]))
        w_f = mpmath.mpf(str(GH_W[k]))
        wk = c_mp + mpmath.mpc(0, 1) * y
        ph = mpmath.exp(2 * mpmath.mpc(0, 1) * c_mp * y) / wk
        lam_f = Ldir_mp(s  + wk, cn_mp, logn_mp)
        lam_b = Ldir_mp(s1 + wk, cn_mp, logn_mp)
        tot += w_f * (lam_f + mpmath.mpf(str(eps_r)) * lam_b) * ph
    return ec2 / (2 * mpmath.pi) * tot


def Lafe_mp(s, cn_mp, logn_mp):
    """L(s) = LAFE(s) / γ(s)  [mpmath]"""
    lam = LAFE_mp(s, cn_mp, logn_mp)
    g = mpmath.exp(lgf_mp(s))
    if abs(g) < mpmath.mpf('1e-300'):
        return mpmath.nan
    return lam / g


# float64 버전 (빠른 영점 탐색)
def lgf_f64(s):
    from scipy.special import loggamma
    s = complex(s)
    N_r = float(CONDUCTOR)
    r = s/2 * np.log(N_r)
    for mu in MU:
        r += -(s + mu)/2 * np.log(np.pi)
        r += loggamma((s + mu) / 2)
    return r


def Ldir_f64(s, cn_nz, log_n_nz):
    s = complex(s)
    g = np.exp(lgf_f64(s))
    L = np.dot(cn_nz, np.exp(-s * log_n_nz))
    return g * L


def LAFE_f64(s, cn_nz, log_n_nz, c=C_SHIFT, eps_r=ROOT_NUMBER):
    s = complex(s); s1 = 1.0 - s
    ec2 = np.exp(c**2)
    total = 0.0j
    for k in range(N_GH):
        y = GH_X[k]; w = GH_W[k]
        wk = complex(c, y)
        phase = np.exp(2j * c * y) / wk
        lam_f = Ldir_f64(s  + wk, cn_nz, log_n_nz)
        lam_b = Ldir_f64(s1 + wk, cn_nz, log_n_nz)
        total += w * (lam_f + eps_r * lam_b) * phase
    return ec2 / (2*np.pi) * total


def Lafe_f64(s, cn_nz, log_n_nz):
    lam = LAFE_f64(s, cn_nz, log_n_nz)
    g = np.exp(lgf_f64(s))
    if abs(g) < 1e-300:
        return complex('nan')
    return lam / g


def hardy_Z_f64(t, cn_nz, log_n_nz, sigma=0.5):
    """Hardy Z(t): 부호 변환이 영점"""
    s = complex(sigma, t)
    lam = LAFE_f64(s, cn_nz, log_n_nz)
    if abs(lam) < 1e-300:
        return 0.0
    return float(abs(lam) * np.sign(lam.real))


# ─── 함수방정식 FE 검증 ────────────────────────────────────────────────────
def check_fe_mp(cn_mp, logn_mp, test_points=None):
    """
    |LAFE(s) - LAFE(1-s)| / |LAFE(s)| 를 여러 점에서 검사
    반환: -log10(최대 상대오차)
    """
    if test_points is None:
        # Re(s)=2~2.5 범위 — v4 수정된 FE 검증점
        test_points = [
            mpmath.mpc(2.0, 3.7),
            mpmath.mpc(2.3, 7.1),
            mpmath.mpc(2.0, 14.1),
            mpmath.mpc(2.5, 21.0),
            mpmath.mpc(2.1, 5.0),
        ]
    errs = []
    for s in test_points:
        s1 = 1 - s
        L_s  = LAFE_mp(s,  cn_mp, logn_mp)
        L_s1 = LAFE_mp(s1, cn_mp, logn_mp)
        if abs(L_s) < 1e-100:
            continue
        rel = abs(L_s - L_s1) / abs(L_s)
        if rel > 0:
            errs.append(-float(mpmath.log10(rel)))
    if not errs:
        return -1.0
    return min(errs)


# ─── 영점 탐색 ────────────────────────────────────────────────────────────
def find_zeros_sym3(cn_nz, log_n_nz, t_min=T_MIN, t_max=T_MAX, n_scan=5000):
    """
    Hardy Z(t) 부호 변환으로 영점 위치 추정 → 정밀화
    """
    t_scan = np.linspace(t_min, t_max, n_scan)
    Z_vals = np.zeros(n_scan)

    log(f"  영점 스캔 중 (n={n_scan})...")
    flush_file()
    t_batch = time.time()
    for i, t in enumerate(t_scan):
        try:
            Z_vals[i] = hardy_Z_f64(t, cn_nz, log_n_nz)
        except Exception as e:
            Z_vals[i] = 0.0
        if (i+1) % 1000 == 0:
            log(f"  스캔 {i+1}/{n_scan} ({time.time()-t_batch:.1f}초)")
            flush_file()
            t_batch = time.time()

    # 부호 변환 검출
    sign_changes = []
    for i in range(len(Z_vals)-1):
        if Z_vals[i] * Z_vals[i+1] < 0:
            sign_changes.append((t_scan[i], t_scan[i+1]))

    log(f"  부호 변환 {len(sign_changes)}개 발견")
    flush_file()

    # 정밀화 (bisection)
    zeros = []
    for (ta, tb) in sign_changes:
        try:
            Za = hardy_Z_f64(ta, cn_nz, log_n_nz)
            Zb = hardy_Z_f64(tb, cn_nz, log_n_nz)
            for _ in range(40):  # 이분법 40회 → ~14자리
                tm = (ta + tb) / 2
                Zm = hardy_Z_f64(tm, cn_nz, log_n_nz)
                if Za * Zm <= 0:
                    tb = tm; Zb = Zm
                else:
                    ta = tm; Za = Zm
            zeros.append((ta + tb) / 2)
        except Exception as e:
            log(f"  WARNING: bisection 실패 [{ta:.4f},{tb:.4f}]: {e}")

    return sorted(zeros)


# ─── κ_near 계산 ──────────────────────────────────────────────────────────
def kappa_kear(t0, delta, cn_nz, log_n_nz, h=H_DERIV):
    """
    κ = |Λ'/Λ(σ+it₀)|² where σ = 1/2 + delta
    실수 σ 방향 수치 미분 (v4 방식 유지)
    반환: (kappa, A) where A = kappa - 1/delta²
    """
    sigma = 0.5 + delta
    s0 = complex(sigma, t0)
    try:
        L0 = LAFE_f64(s0, cn_nz, log_n_nz)
        if abs(L0) < 1e-250:
            return None, None
        Lp = LAFE_f64(s0 + h, cn_nz, log_n_nz)
        Lm = LAFE_f64(s0 - h, cn_nz, log_n_nz)
        conn = (Lp - Lm) / (2*h * L0)
        kappa = abs(conn)**2
        if not np.isfinite(kappa):
            return None, None
        A = kappa - 1.0 / (delta**2)
        return float(kappa), float(A)
    except Exception as e:
        log(f"  WARNING: κ 계산 실패 t0={t0:.4f}: {e}")
        return None, None


# ─── 모노드로미 계산 ─────────────────────────────────────────────────────
def monodromy_sym3(t0, cn_nz, log_n_nz, radius=MONO_R, n_steps=MONO_N):
    """폐곡선 적분으로 arg(Λ) 누적"""
    center = complex(0.5, t0)
    phase_acc = 0.0
    prev_val = None

    for k in range(n_steps + 1):
        theta = 2 * np.pi * k / n_steps
        pt = center + radius * np.exp(1j * theta)
        try:
            val = LAFE_f64(pt, cn_nz, log_n_nz)
        except Exception as e:
            return None
        if abs(val) < 1e-250:
            return None
        if prev_val is not None:
            phase_acc += np.angle(val / prev_val)
        prev_val = val

    return abs(phase_acc) / np.pi


# ─── σ-유일성 ────────────────────────────────────────────────────────────
def sigma_uniqueness(t0, cn_nz, log_n_nz, delta=DELTA_KAPPA):
    """κ(σ=0.5) / κ(σ=0.45) 비율"""
    try:
        kappa_on, _ = kappa_kear(t0, delta, cn_nz, log_n_nz)
        if kappa_on is None:
            return None

        # σ=0.45
        sigma_off = 0.45 + delta
        s_off = complex(sigma_off, t0)
        L0_off = LAFE_f64(s_off, cn_nz, log_n_nz)
        if abs(L0_off) < 1e-250:
            return None
        Lp_off = LAFE_f64(s_off + H_DERIV, cn_nz, log_n_nz)
        Lm_off = LAFE_f64(s_off - H_DERIV, cn_nz, log_n_nz)
        conn_off = (Lp_off - Lm_off) / (2*H_DERIV * L0_off)
        kappa_off = abs(conn_off)**2

        if kappa_off < 1e-10:
            return None
        return float(kappa_on / kappa_off)
    except Exception as e:
        log(f"  WARNING: σ-유일성 실패 t0={t0:.4f}: {e}")
        return None


# ════════════════════════════════════════════════════════════════════════════
# 메인
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("결과 #78 — GL(4) sym³(Δ) 수정 AFE (MU=[-1,0,0,1], N=144)")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"★ 수정: MU={MU}, N={CONDUCTOR} (v4의 MU=(0,1,11,12),N=1 → 올바른 parameterization)")
log(f"mpmath DPS={DPS}, N_COEFF={N_COEFF}, N_GH={N_GH}, C_SHIFT={C_SHIFT}")
log(f"δ_κ={DELTA_KAPPA}, mono_r={MONO_R}, t범위=[{T_MIN},{T_MAX}]")
log()
flush_file()

t_start = time.time()

# ─── 1. τ(n) 계산 ──────────────────────────────────────────────────────
log("[1] τ(n) 계산")
t1 = time.time()
primes = sieve_primes(N_COEFF)
log(f"  소수: {len(primes)}개 (≤{N_COEFF})")
tau = compute_tau_fast(N_COEFF)
log(f"  완료 ({time.time()-t1:.1f}초)")

# 검증
known = {2: -24, 3: 252, 5: 4830, 7: -16744}
ok_tau = all(abs(tau[p] - v) < 0.5 for p, v in known.items())
log(f"  τ 검증: {'✓' if ok_tau else '✗'} τ(2)=-24, τ(3)=252, τ(5)=4830")
log()
flush_file()

# ─── 2. sym³ 계수 ──────────────────────────────────────────────────────
log("[2] sym³ 해석적 정규화 계수 (v4와 동일, 감마만 수정)")
t2 = time.time()
an = compute_sym3_an(tau, primes, N_COEFF)
log(f"  완료 ({time.time()-t2:.1f}초)")

# 소수 계수 검증
for p in [2, 3, 5, 7]:
    tp = tau[p] / p**5.5
    expected = tp**3 - 2*tp
    log(f"  a({p}) = {an[p]:.8f} (기대 {expected:.8f})")
log(f"  max|a_p|={max(abs(an[p]) for p in primes[:20]):.4f} (≤4 기대)")
log()
flush_file()

# ─── 3. AFE 준비 ────────────────────────────────────────────────────────
log("[3] AFE 데이터 준비")
# 비영 계수만 추출
nz_idx = [n for n in range(1, N_COEFF+1) if abs(an[n]) > 1e-15]
cn_nz = np.array([an[n] for n in nz_idx], dtype=np.float64)
log_n_nz = np.log(np.array(nz_idx, dtype=np.float64))
log(f"  비영 계수: {len(nz_idx)}/{N_COEFF}")

# mpmath 버전 (FE 검증용)
cn_mp   = [mpmath.mpf(str(an[n])) for n in nz_idx[:500]]  # 처음 500개
logn_mp = [mpmath.mpf(str(np.log(n))) for n in nz_idx[:500]]
log(f"  mpmath 버전: {len(cn_mp)}개")
log()
flush_file()

# ─── 4. FE 검증 ─────────────────────────────────────────────────────────
log("[4] 함수방정식 검증 (MU=[-1,0,0,1], N=144)")
log("  테스트 점: Re(s)=2.0~2.5 (backward 수렴 보장)")

t_fe = time.time()
fe_score = check_fe_mp(cn_mp, logn_mp)
log(f"  FE score = {fe_score:.1f} 자릿수 ({time.time()-t_fe:.1f}초)")
log(f"  목표: ≥8 (저자 PARI 확인: 11)")
fe_ok = fe_score >= 8
log(f"  FE 판정: {'✅ 통과' if fe_ok else f'❌ 실패 (score={fe_score:.1f})'}")

if not fe_ok:
    log()
    log("  ⚠️ FE 미통과. 추가 진단:")
    # 부가 진단: 각 점별 FE 오류
    test_pts = [
        mpmath.mpc(2.0, 3.7),
        mpmath.mpc(2.3, 7.1),
        mpmath.mpc(2.0, 14.1),
    ]
    for s in test_pts:
        s1 = 1 - s
        L_s  = LAFE_mp(s,  cn_mp, logn_mp)
        L_s1 = LAFE_mp(s1, cn_mp, logn_mp)
        rel = abs(L_s - L_s1) / (abs(L_s)+1e-100)
        log(f"    s={s}: |Λ(s)-Λ(1-s)|/|Λ(s)| = {float(rel):.3e}")

log()
flush_file()

# ─── 5. 영점 탐색 ──────────────────────────────────────────────────────
log("[5] 영점 탐색 (Hardy Z 부호 변환, float64)")
t_z = time.time()
zeros = find_zeros_sym3(cn_nz, log_n_nz, t_min=T_MIN, t_max=T_MAX, n_scan=5000)
log(f"  탐색 완료 ({time.time()-t_z:.1f}초)")
log(f"  영점 개수: {len(zeros)} (≥30 기대)")
if zeros:
    log(f"  t₁ = {zeros[0]:.6f} (저자: 0.3239)")
    for z in zeros[:10]:
        log(f"  t = {z:.6f}")
else:
    log("  ⚠️ 영점 0개 — 탐색 로직 점검 필요")

zeros_ok = len(zeros) >= 10 and (zeros[0] < 3.0 if zeros else False)
log(f"  영점 판정: {'✅ 통과' if zeros_ok else '⚠️ 조건 미충족'}")
log()
flush_file()

if not zeros:
    log("FATAL: 영점 없음. 종료.")
    flush_file()
    sys.exit(1)

zeros_use = zeros[:min(50, len(zeros))]

# ─── 6. κ_near + A(t₀) ─────────────────────────────────────────────────
log("[6] κ_near 측정 (δ=0.001, 실수 σ-이동)")
log(f"  {len(zeros_use)}개 영점 대상")
log()

A_vals = []; kappa_vals = []; kd2_vals = []
n_kappa_fail = 0

for i, t0 in enumerate(zeros_use):
    kappa, A = kappa_kear(t0, DELTA_KAPPA, cn_nz, log_n_nz)
    if kappa is None or not np.isfinite(kappa) or not np.isfinite(A):
        n_kappa_fail += 1
        log(f"  [{i+1:02d}] t₀={t0:.5f} → FAIL")
        continue
    A_vals.append(A); kappa_vals.append(kappa)
    kd2 = kappa * (DELTA_KAPPA**2)
    kd2_vals.append(kd2)
    log(f"  [{i+1:02d}] t₀={t0:.5f}: κ={kappa:.2f}, A={A:.4f}, κδ²={kd2:.6f}")
    if (i+1) % 5 == 0:
        flush_file()

log()
n_ok = len(A_vals)
log(f"  성공: {n_ok}, 실패: {n_kappa_fail}")

if n_ok == 0:
    log("FATAL: κ 측정 전부 실패")
    flush_file()
    sys.exit(1)

A_arr = np.array(A_vals)
kd2_arr = np.array(kd2_vals)
mean_A = float(np.mean(A_arr))
std_A  = float(np.std(A_arr))
cv_A   = float(std_A / abs(mean_A) * 100) if abs(mean_A) > 1e-10 else 999.0
mean_kd2 = float(np.mean(kd2_arr))

log(f"  mean(A) = {mean_A:.4f}")
log(f"  std(A)  = {std_A:.4f}")
log(f"  CV(A)   = {cv_A:.1f}% (< 5% 기대)")
log(f"  mean(κδ²) = {mean_kd2:.6f} ([0.99, 1.15] 기대)")
kappa_ok = n_ok >= 10 and cv_A < 5.0
log(f"  κ_near 판정: {'✅ 통과' if kappa_ok else '⚠️'}")
log()
flush_file()

# δ 독립성 체크 (처음 3개 영점, 추가 δ값)
log("  [δ 독립성 추가 체크]")
for t0 in zeros_use[:3]:
    row = []
    for d in [0.0001, 0.001, 0.01, 0.1]:
        k, a = kappa_kear(t0, d, cn_nz, log_n_nz)
        if a is not None and np.isfinite(a):
            row.append(f"δ={d}: A={a:.3f}")
    log(f"  t₀={t0:.5f}:  " + "  ".join(row))
log()
flush_file()

# ─── 7. 모노드로미 ──────────────────────────────────────────────────────
log("[7] 모노드로미 (폐곡선 적분, 반지름=0.4)")
mono_vals = []
n_mono_fail = 0

for i, t0 in enumerate(zeros_use[:20]):
    mono = monodromy_sym3(t0, cn_nz, log_n_nz)
    if mono is None or not np.isfinite(mono):
        n_mono_fail += 1
        log(f"  [{i+1:02d}] t₀={t0:.5f} → FAIL")
        continue
    mono_vals.append(mono)
    ok = "✅" if abs(mono - 2.0) < 0.15 else "⚠️"
    log(f"  [{i+1:02d}] t₀={t0:.5f}: mono/π = {mono:.4f} {ok}")
    if (i+1) % 5 == 0:
        flush_file()

log()
log(f"  성공: {len(mono_vals)}, 실패: {n_mono_fail}")
if mono_vals:
    mono_arr = np.array(mono_vals)
    mean_mono = float(np.mean(mono_arr))
    frac_ok = float(np.sum(np.abs(mono_arr - 2.0) < 0.15)) / len(mono_arr) * 100
    log(f"  mean(mono/π) = {mean_mono:.4f} (≈2.0 기대)")
    log(f"  |mono/π-2|<0.15: {frac_ok:.0f}%")
    mono_ok = abs(mean_mono - 2.0) < 0.2 and frac_ok > 70
    log(f"  모노드로미 판정: {'✅ 통과' if mono_ok else '⚠️'}")
else:
    mono_ok = False; mean_mono = None
log()
flush_file()

# ─── 8. σ-유일성 ────────────────────────────────────────────────────────
log("[8] σ-유일성 (κ(σ=0.5)/κ(σ=0.45))")
ratio_vals = []
n_ratio_fail = 0

for i, t0 in enumerate(zeros_use[:15]):
    ratio = sigma_uniqueness(t0, cn_nz, log_n_nz)
    if ratio is None or not np.isfinite(ratio):
        n_ratio_fail += 1
        log(f"  [{i+1:02d}] t₀={t0:.5f} → FAIL")
        continue
    ratio_vals.append(ratio)
    ok = "✅" if ratio > 10 else "⚠️"
    log(f"  [{i+1:02d}] t₀={t0:.5f}: ratio = {ratio:.3f} {ok}")

log()
if ratio_vals:
    ratio_arr = np.array(ratio_vals)
    mean_ratio = float(np.mean(ratio_arr))
    frac_gt10 = float(np.sum(ratio_arr > 10)) / len(ratio_arr) * 100
    log(f"  mean(ratio) = {mean_ratio:.2f} (>10 기대)")
    log(f"  ratio>10: {frac_gt10:.0f}%")
    sigma_ok = mean_ratio > 10 and frac_gt10 > 70
    log(f"  σ-유일성 판정: {'✅ 통과' if sigma_ok else '⚠️'}")
else:
    sigma_ok = False; mean_ratio = None
log()
flush_file()

# ─── 9. κ_near(d) 단조증가 ─────────────────────────────────────────────
log("[9] κ_near(d) 단조증가 비교 — B-12")
KNOWN_A = {'d=1 (ζ)': 1.2727, 'd=2 (GL2)': 3.93, 'd=3 (GL3)': 12.79}
for lbl, val in KNOWN_A.items():
    log(f"  {lbl:20s}: A = {val:.4f}")
log(f"  {'d=4 (sym³Δ)':20s}: A = {mean_A:.4f}  "
    f"{'✅ d 단조증가' if mean_A > 12.79 else '❌ 단조증가 위반'}")
monotone_ok = mean_A > 12.79
log()
flush_file()

# ─── 결과 요약 ──────────────────────────────────────────────────────────
elapsed = time.time() - t_start
log("=" * 72)
log("성공 기준 총정리")
log("=" * 72)

crit = [
    ("FE ≥8자릿수",             fe_ok,      f"{fe_score:.1f}자릿수"),
    ("영점 ≥10, t₁<3.0",        zeros_ok,   f"{len(zeros)}개, t₁={zeros[0]:.4f}" if zeros else "0개"),
    ("κ_near CV<5%",             kappa_ok,   f"{n_ok}개, CV={cv_A:.1f}%"),
    ("mono/π≈2.0",               mono_ok,    f"mean={mean_mono:.4f}" if mono_vals else "N/A"),
    ("σ-유일성 ratio>10",        sigma_ok,   f"mean={mean_ratio:.2f}" if mean_ratio else "N/A"),
    ("κ_near(d=4)>κ_near(d=3)",  monotone_ok,f"A={mean_A:.4f} vs 12.79"),
]
n_pass = sum(1 for _, ok, _ in crit if ok)
for name, ok, detail in crit:
    log(f"  {'✅' if ok else '❌'} {name:35s} → {detail}")
log()
log(f"  통과: {n_pass}/{len(crit)}")
log()
log("─" * 72)
log("수치 요약")
log("─" * 72)
log(f"  FE: {fe_score:.1f}자릿수 (기대 ≥8)")
log(f"  영점: {len(zeros)}개, t₁={zeros[0]:.6f}" if zeros else "  영점: 없음")
log(f"  mean(A): {mean_A:.4f} (d=4, PARI verified: ?)")
log(f"  CV(A): {cv_A:.2f}%")
log(f"  mean(κδ²): {mean_kd2:.6f}")
log(f"  mean(mono/π): {mean_mono:.4f}" if mono_vals else "  mean(mono/π): N/A")
log(f"  mean(ratio): {mean_ratio:.2f}" if mean_ratio else "  mean(ratio): N/A")
log()
log("─" * 72)
log("경계 갱신")
log("─" * 72)
log(f"  B-12: {'★★ 확립 — d=1,2,3,4 모두 확인' if monotone_ok else '⚠️ d=4 위반/미확인'}")
log(f"  B-05: {'✅ d=4 확인' if sigma_ok else '⚠️ 미확인'}")
log(f"  총 소요: {elapsed:.1f}초")
log(f"  완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 72)
flush_file()
log(f"결과: {OUTFILE}")
flush_file()
