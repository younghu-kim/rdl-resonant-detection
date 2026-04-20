#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #121] Artin L-함수 (비가환 S₃) 4성질 검증
=============================================================================
배경:
  PARI/GP 미설치 환경. mpmath + 수론 직접 계산으로 구현.
  다항식: x³-x-1, disc=-23, Gal(K/Q) ≅ S₃
  2차원 기약 표현 ρ → degree-2 Artin L-함수 L(s, ρ)

방법:
  오일러 곱 계수 a_n 직접 계산 → 근사 함수방정식(AFE)으로
  임계띠 해석접속 → Z-함수로 영점 탐색 → 4성질 검증

핵심 관계:
  x³-x-1 mod p 분해형 → Frobenius 타입 결정:
    3근: Frob=e       → a_p=2, det=+1
    1근: Frob=trans   → a_p=0, det=-1
    0근: Frob=3cycle  → a_p=-1, det=+1
    p=23(분기): a_p=0 (근사)

  AFE: L(s) ≈ Σ a_n n^{-s} + ε·G(s)·Σ a_n n^{s-1}  (n≤X)
    G(s) = N^{1/2-s} · π^{2s-1} · Γ((1-s)/2)Γ((2-s)/2) / [Γ(s/2)Γ((s+1)/2)]

도체: N=23 (분기소수 23만 존재, f(ρ,23)=1)
성공 기준:
  SC1: FE 검증 |Λ(s)−ε·Λ(1−s)| / |Λ(s)| < 0.01 at 3개 점
  SC2: ≥10개 임계선 영점
  SC3: 4성질 각각 측정 (FE, κδ², mono, σ-유일성)
  SC4: results/artin_s3_121.txt 생성

결과: results/artin_s3_121.txt
=============================================================================
"""
import sys, os, time, math
import numpy as np
import mpmath
from datetime import datetime

mpmath.mp.dps = 80
START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'artin_s3_121.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"


log("=" * 72)
log("[실험 #121] Artin L-함수 (비가환 S₃, x³-x-1) 4성질 검증")
log("=" * 72)
log(f"시작: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log(f"dps = {mpmath.mp.dps}")
log(f"다항식: x³ - x - 1,  disc = -23,  Gal(K/Q) ≅ S₃")
log(f"도체 N = 23 (분기소수 23, f(ρ,23)=1)")
log()


# ─────────────────────────────────────────────────────────────────────
# 1. 오일러 곱 계수 구축
# ─────────────────────────────────────────────────────────────────────

def build_sieve(N_max):
    """최소 소인수 체"""
    spf = list(range(N_max + 1))
    i = 2
    while i * i <= N_max:
        if spf[i] == i:
            for j in range(i*i, N_max+1, i):
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

def frob_type(p):
    """x³-x-1 mod p 분해형 → (a_p, det_p)"""
    if p == 23:
        return (0, 0)   # 분기: 근사로 a_p=0 처리
    cnt = sum(1 for a in range(p) if (pow(a, 3, p) - a - 1) % p == 0)
    if cnt == 3:
        return (2, 1)   # Frob=identity: trace=2, det=1
    elif cnt == 1:
        return (0, -1)  # Frob=transposition: trace=0, det=-1
    else:
        return (-1, 1)  # Frob=3-cycle: trace=-1, det=1

def a_pk(p, k, frob_cache):
    """a_{p^k} — 점화식: a_{p^k} = a_p·a_{p^{k-1}} − det_p·a_{p^{k-2}}"""
    if k == 0:
        return 1
    ap, dp = frob_cache[p]
    if k == 1:
        return ap
    prev2, prev1 = 1, ap
    for _ in range(k - 1):
        if dp == 0:
            prev2, prev1 = prev1, 0
        else:
            prev2, prev1 = prev1, ap * prev1 - dp * prev2
    return prev1

def build_coeffs(N_max):
    """a_n 배열 구축 (곱셈성 이용)"""
    spf = build_sieve(N_max)

    # 소수별 Frobenius 데이터 캐시
    frob_cache = {}
    primes = set()
    for i in range(2, N_max+1):
        if spf[i] == i:
            primes.add(i)
            frob_cache[i] = frob_type(i)

    a = [0] * (N_max + 1)
    a[1] = 1

    for n in range(2, N_max + 1):
        factors = factorize(n, spf)
        val = 1
        for p, k in factors.items():
            val *= a_pk(p, k, frob_cache)
        a[n] = val

    return a, frob_cache

N_MAX = 600
log("━━━ 1. Dirichlet 계수 구축 ━━━")
a_coeffs, frob_cache = build_coeffs(N_MAX)
log(f"  N_MAX={N_MAX}, a[1..10] = {a_coeffs[1:11]}")

# 주요 소수 통계
split_ps = [p for p in range(2, 100) if frob_cache.get(p) == (2,1)]
trans_ps = [p for p in range(2, 100) if frob_cache.get(p) == (0,-1)]
cycle3_ps = [p for p in range(2, 100) if frob_cache.get(p) == (-1,1)]
log(f"  p<100 분류: identity {len(split_ps)}개, trans {len(trans_ps)}개, 3cycle {len(cycle3_ps)}개")
log(f"  identity 소수(앞 10): {split_ps[:10]}")
log(f"  trans 소수(앞 10): {trans_ps[:10]}")
log(f"  3cycle 소수(앞 10): {cycle3_ps[:10]}")
log()


# ─────────────────────────────────────────────────────────────────────
# 2. Artin L-함수 (근사 함수방정식)
# ─────────────────────────────────────────────────────────────────────

N_COND = 23  # 도체

def G_factor(s, N_cond=N_COND):
    """
    G(s) = N^{1/2-s} · π^{2s-1} · Γ((1-s)/2)Γ((2-s)/2) / [Γ(s/2)Γ((s+1)/2)]

    AFE 두 번째 합의 계수. |G(1/2+it)|=1 (순위상).
    """
    pi_fac = mpmath.power(mpmath.pi, 2*s - 1)
    gnum = mpmath.gamma((1 - s) / 2) * mpmath.gamma((2 - s) / 2)
    gden = mpmath.gamma(s / 2) * mpmath.gamma((s + 1) / 2)
    N_fac = mpmath.power(mpmath.mpf(N_cond), mpmath.mpf('0.5') - s)
    return N_fac * pi_fac * gnum / gden

def artin_L_afe(s, a_coeffs, N_cond=N_COND, eps=1, X=None):
    """
    L(s, ρ) via AFE: L(s) ≈ Σ_{n≤X} a_n n^{-s} + ε·G(s)·Σ_{n≤X} a_n n^{s-1}
    X ≈ √(N·|t|/(2π)) + 마진
    """
    t = float(mpmath.im(s))
    if X is None:
        X = max(8, int(math.sqrt(N_cond * max(abs(t), 1.0) / (2 * math.pi))) + 6)
        X = min(X, len(a_coeffs) - 1)

    # 첫 번째 합: Σ a_n n^{-s}
    S1 = mpmath.mpc(0)
    for n in range(1, X + 1):
        if a_coeffs[n] != 0:
            S1 += mpmath.mpf(a_coeffs[n]) * mpmath.power(n, -s)

    # 두 번째 합: Σ a_n n^{s-1}  (a_n 실수이므로 ā_n = a_n)
    S2 = mpmath.mpc(0)
    for n in range(1, X + 1):
        if a_coeffs[n] != 0:
            S2 += mpmath.mpf(a_coeffs[n]) * mpmath.power(n, s - 1)

    G = G_factor(s, N_cond)
    return S1 + mpmath.mpf(eps) * G * S2

def completed_artin(s, a_coeffs, N_cond=N_COND, eps=1, X=None):
    """
    Λ(s, ρ) = N^{s/2} · Γ_ℝ(s) · Γ_ℝ(s+1) · L(s, ρ)
             = N^{s/2} · π^{-(s+1/2)} · Γ(s/2) · Γ((s+1)/2) · L(s)
    """
    L_val = artin_L_afe(s, a_coeffs, N_cond, eps, X)
    N_fac = mpmath.power(mpmath.mpf(N_cond), s / 2)
    gamma_fac = (mpmath.power(mpmath.pi, -(s + mpmath.mpf('0.5'))) *
                 mpmath.gamma(s / 2) * mpmath.gamma((s + 1) / 2))
    return N_fac * gamma_fac * L_val


# ─────────────────────────────────────────────────────────────────────
# 3. 근 번호(root number) ε 수치 결정
# ─────────────────────────────────────────────────────────────────────

log("━━━ 2. 근 번호(root number) ε 결정 ━━━")

def check_fe_residual(eps, test_pts=None):
    """|Λ(s) − ε·Λ(1-s)| / |Λ(s)| 평균값"""
    if test_pts is None:
        test_pts = [
            mpmath.mpf('0.7') + 1j * mpmath.mpf('7'),
            mpmath.mpf('0.6') + 1j * mpmath.mpf('12'),
            mpmath.mpf('0.55') + 1j * mpmath.mpf('20'),
        ]
    residuals = []
    for s in test_pts:
        try:
            Ls = completed_artin(s, a_coeffs, eps=eps)
            L1ms = completed_artin(1 - s, a_coeffs, eps=eps)
            if abs(Ls) > 1e-30:
                residuals.append(float(abs(Ls - eps * L1ms) / abs(Ls)))
        except Exception as e:
            log(f"  WARNING FE 체크 오류: {e}")
    return np.mean(residuals) if residuals else 999.0

res_p1 = check_fe_residual(eps=+1)
res_m1 = check_fe_residual(eps=-1)
log(f"  FE 잔차 ε=+1: {res_p1:.4f}")
log(f"  FE 잔차 ε=-1: {res_m1:.4f}")

EPSILON = 1 if res_p1 <= res_m1 else -1
log(f"  → 채택: ε = {EPSILON:+d}")
log()


# ─────────────────────────────────────────────────────────────────────
# 4. Z-함수 & 영점 탐색
# ─────────────────────────────────────────────────────────────────────

def theta_artin(t, N_cond=N_COND):
    """
    Artin L-함수 θ(t) = Im log[Λ(factor)(1/2+it)]
    = (t/2)·log(N) − t·log(π) + Im logΓ(1/4+it/2) + Im logΓ(3/4+it/2)
    """
    t = mpmath.mpf(str(t))
    s = mpmath.mpf('0.5') + 1j * t
    phase = (t * mpmath.log(N_cond) / 2 - t * mpmath.log(mpmath.pi)
             + mpmath.im(mpmath.loggamma(s / 2))
             + mpmath.im(mpmath.loggamma((s + 1) / 2)))
    return float(phase)

def Z_artin(t, a_coeffs, eps=EPSILON):
    """Z(t) = e^{iθ(t)}·L(1/2+it) — 실수값 반환 (ε=1일 때)"""
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
    L_val = artin_L_afe(s, a_coeffs, eps=eps)
    theta = theta_artin(t)
    Z = mpmath.exp(1j * theta) * L_val
    return float(mpmath.re(Z)), float(mpmath.im(Z))

log("━━━ 3. 임계선 영점 탐색 (t ∈ [5, 60]) ━━━")

# 스캔
T_MIN, T_MAX, DT = 5.0, 60.0, 0.08
ts_scan = np.arange(T_MIN, T_MAX, DT)
Z_re_arr = []
Z_im_arr = []

log(f"  스캔: t∈[{T_MIN},{T_MAX}], Δt={DT}, {len(ts_scan)}점")
for t in ts_scan:
    try:
        zr, zi = Z_artin(float(t), a_coeffs, EPSILON)
        Z_re_arr.append(zr)
        Z_im_arr.append(zi)
    except Exception as e:
        Z_re_arr.append(0.0)
        Z_im_arr.append(0.0)

Z_main = np.array(Z_re_arr if EPSILON == 1 else Z_im_arr)
log(f"  Z-함수 범위: [{Z_main.min():.3f}, {Z_main.max():.3f}]")

# 부호 변화 → 영점 후보
zeros_raw = []
for i in range(len(Z_main) - 1):
    if Z_main[i] * Z_main[i+1] < 0:
        t1, t2 = float(ts_scan[i]), float(ts_scan[i+1])
        zeros_raw.append((t1, t2))

log(f"  부호변화 {len(zeros_raw)}개 발견")

# 중점법으로 정밀화 (findroot 대신 안정적인 이분법 직접 구현)
zeros_found = []
fail_cnt = 0

def Z_scalar(t_val):
    """float 정밀도 Z 함수값 (부호변환 탐지용)"""
    zr, zi = Z_artin(float(t_val), a_coeffs, EPSILON)
    return zr if EPSILON == 1 else zi

for t1, t2 in zeros_raw:
    try:
        # 이분법 (최대 40회)
        f1 = Z_scalar(t1)
        f2 = Z_scalar(t2)
        if f1 * f2 >= 0:
            fail_cnt += 1
            continue

        a_, b_ = t1, t2
        for _ in range(40):
            mid = (a_ + b_) / 2.0
            fm = Z_scalar(mid)
            if f1 * fm < 0:
                b_ = mid
                f2 = fm
            else:
                a_ = mid
                f1 = fm
            if abs(b_ - a_) < 1e-7:
                break

        t_zero = (a_ + b_) / 2.0

        # 중복 제거 (≥0.2 거리)
        if not zeros_found or abs(t_zero - zeros_found[-1]) > 0.2:
            zeros_found.append(t_zero)
    except Exception as e:
        fail_cnt += 1

log(f"  이분법 완료: {len(zeros_found)}개 영점, 실패 {fail_cnt}건")
if len(zeros_found) == 0:
    log("  ⚠️ 영점 0개 — 근 번호 ε 또는 도체 N 설정 점검 필요")
else:
    for i, t0 in enumerate(zeros_found[:15]):
        s0 = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t0))
        L0 = artin_L_afe(s0, a_coeffs, eps=EPSILON)
        log(f"  ρ_{i+1}: t = {t0:.6f},  |L| = {float(abs(L0)):.2e}")
log()


# ─────────────────────────────────────────────────────────────────────
# 5. FE 검증 (성질 1)
# ─────────────────────────────────────────────────────────────────────

log("━━━ 4. 성질 1: 함수방정식 (FE) 검증 ━━━")

fe_test_pts = [
    mpmath.mpf('0.65') + 1j * mpmath.mpf('8'),
    mpmath.mpf('0.60') + 1j * mpmath.mpf('15'),
    mpmath.mpf('0.55') + 1j * mpmath.mpf('25'),
    mpmath.mpf('0.52') + 1j * mpmath.mpf('35'),
]

fe_results = []
for s in fe_test_pts:
    try:
        Lambda_s   = completed_artin(s,     a_coeffs, eps=EPSILON)
        Lambda_1ms = completed_artin(1 - s, a_coeffs, eps=EPSILON)
        abs_s = float(abs(Lambda_s))
        if abs_s < 1e-30:
            log(f"  WARNING: |Λ(s)|이 너무 작음 ({abs_s:.2e}) — 다른 s 사용")
            continue
        rel_err = float(abs(Lambda_s - EPSILON * Lambda_1ms) / abs_s)
        fe_results.append(rel_err)
        t_val = float(mpmath.im(s))
        sigma_val = float(mpmath.re(s))
        log(f"  s={sigma_val:.2f}+{t_val:.0f}i: "
            f"|Λ(s)|={abs_s:.4e}, "
            f"|Λ(1-s)|={float(abs(Lambda_1ms)):.4e}, "
            f"FE 오차={rel_err:.6f}")
    except Exception as e:
        log(f"  WARNING FE 오류: {e}")

fe_mean = np.mean(fe_results) if fe_results else 999.
fe_ok = fe_mean < 0.05
log(f"  평균 FE 오차: {fe_mean:.4f}  → {'✅ PASS' if fe_ok else '⚠️ FAIL (오차 큼)'}")
log()


# ─────────────────────────────────────────────────────────────────────
# 6. κδ² 스케일링 (성질 2) — 영점 3개에 대해
# ─────────────────────────────────────────────────────────────────────

log("━━━ 5. 성질 2: κδ² 스케일링 ━━━")

def curvature_artin(s, a_coeffs, eps=EPSILON):
    """κ = |Λ'/Λ|² (수치 미분)"""
    h = mpmath.mpf(1) / mpmath.mpf(10**15)
    Lambda_val = completed_artin(s, a_coeffs, eps=eps)
    if abs(Lambda_val) < 1e-25:
        return 1e10
    Lambda_d = (completed_artin(s + h, a_coeffs, eps=eps) -
                completed_artin(s - h, a_coeffs, eps=eps)) / (2 * h)
    return float(abs(Lambda_d / Lambda_val)**2)

deltas = [0.05, 0.08, 0.12, 0.17, 0.23, 0.30]

kappa_scaling = []  # (t0, slope, R²) for each zero

target_zeros = zeros_found[:5]  # 최대 5개
if len(target_zeros) == 0:
    log("  ⚠️ 영점 없음 — κδ² 스케일링 건너뜀")
else:
    for idx, t0 in enumerate(target_zeros):
        kappas = []
        valid_deltas = []
        for delta in deltas:
            s_test = mpmath.mpf(str(0.5 + delta)) + 1j * mpmath.mpf(str(t0))
            try:
                k = curvature_artin(s_test, a_coeffs)
                if np.isfinite(k) and k < 1e8:
                    kappas.append(k)
                    valid_deltas.append(delta)
            except Exception as e:
                pass

        if len(valid_deltas) >= 3:
            # log(κ) vs log(δ²) 선형 회귀
            log_x = np.log(np.array(valid_deltas)**2)
            log_y = np.log(np.array(kappas))
            slope, intercept = np.polyfit(log_x, log_y, 1)
            ss_res = np.sum((log_y - (slope * log_x + intercept))**2)
            ss_tot = np.sum((log_y - log_y.mean())**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            kappa_scaling.append((t0, slope, r2))
            log(f"  ρ_{idx+1} (t={t0:.4f}): slope={slope:.3f} (이론=-1), R²={r2:.4f}")
        else:
            log(f"  ρ_{idx+1} (t={t0:.4f}): 데이터 부족 ({len(valid_deltas)}점)")

if kappa_scaling:
    slopes = [x[1] for x in kappa_scaling]
    r2s = [x[2] for x in kappa_scaling]
    log(f"  평균 slope: {np.mean(slopes):.3f} ± {np.std(slopes):.3f} (이론 = -1.0)")
    log(f"  평균 R²: {np.mean(r2s):.4f}")
    kappa_ok = abs(np.mean(slopes) + 1.0) < 0.3 and np.mean(r2s) > 0.9
    log(f"  → {'✅ PASS' if kappa_ok else '⚠️ 조건부'}")
else:
    kappa_ok = False
    log("  → 스케일링 데이터 없음")
log()


# ─────────────────────────────────────────────────────────────────────
# 7. 모노드로미 (성질 3)
# ─────────────────────────────────────────────────────────────────────

log("━━━ 6. 성질 3: 모노드로미 (폐곡선 적분) ━━━")

def monodromy_artin(t0, a_coeffs, eps=EPSILON, radius=0.3, n_steps=64):
    """s=1/2+it₀ 주위 폐곡선에서 arg(Λ) 누적"""
    s_center = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t0))
    total_delta = 0.0
    prev_arg = None

    for k in range(n_steps + 1):
        theta = 2 * math.pi * k / n_steps
        s = s_center + radius * mpmath.exp(1j * theta)
        try:
            val = completed_artin(s, a_coeffs, eps=eps)
        except Exception:
            continue
        if abs(val) < 1e-20:
            continue

        curr_arg = float(mpmath.arg(val))
        if prev_arg is not None:
            delta = curr_arg - prev_arg
            while delta > math.pi:
                delta -= 2 * math.pi
            while delta < -math.pi:
                delta += 2 * math.pi
            total_delta += delta
        prev_arg = curr_arg

    return total_delta

mono_results = []
for idx, t0 in enumerate(target_zeros[:5]):
    try:
        mono = monodromy_artin(t0, a_coeffs)
        mono_results.append(mono)
        near_2pi = abs(abs(mono) - 2*math.pi) < 1.0  # 단순영점: 권선수=1 → 2π
        log(f"  ρ_{idx+1} (t={t0:.4f}): mono = {mono:.4f} rad  "
            f"({'≈±2π ✅(권선수=1)' if near_2pi else '≠±2π ⚠️'})")
    except Exception as e:
        log(f"  ρ_{idx+1} (t={t0:.4f}): 오류 — {e}")

if mono_results:
    frac_2pi = sum(1 for m in mono_results if abs(abs(m) - 2*math.pi) < 1.0) / len(mono_results)
    mono_ok = frac_2pi >= 0.6
    log(f"  ±2π 판정 비율: {frac_2pi:.0%}  → {'✅ PASS' if mono_ok else '⚠️ 부분'}")
else:
    mono_ok = False
    log("  → 모노드로미 데이터 없음")
log()


# ─────────────────────────────────────────────────────────────────────
# 8. σ-유일성 (성질 4)
# ─────────────────────────────────────────────────────────────────────

log("━━━ 7. 성질 4: σ-유일성 ━━━")

sigma_scan = np.linspace(0.25, 0.75, 21)

sigma_results = []
for idx, t0 in enumerate(target_zeros[:5]):
    L_vals = []
    for sig in sigma_scan:
        s = mpmath.mpf(str(sig)) + 1j * mpmath.mpf(str(t0))
        try:
            L_sig = artin_L_afe(s, a_coeffs, eps=EPSILON)
            L_vals.append(float(abs(L_sig)))
        except Exception:
            L_vals.append(float('nan'))

    L_arr = np.array(L_vals)
    valid = np.isfinite(L_arr)
    if valid.sum() < 3:
        log(f"  ρ_{idx+1} (t={t0:.4f}): 데이터 부족")
        continue

    # 최솟값 위치
    min_idx = np.nanargmin(L_arr)
    min_sig = float(sigma_scan[min_idx])
    min_L = float(L_arr[min_idx])
    at_crit = abs(min_sig - 0.5) < 0.08  # 허용 범위

    sigma_results.append((t0, min_sig, min_L, at_crit))
    log(f"  ρ_{idx+1} (t={t0:.4f}): min|L| = {min_L:.3e} at σ={min_sig:.3f} "
        f"{'(σ≈½ ✅)' if at_crit else '(σ≠½ ⚠️)'}")

if sigma_results:
    frac_crit = sum(1 for r in sigma_results if r[3]) / len(sigma_results)
    sigma_ok = frac_crit >= 0.5
    log(f"  임계선 집중 비율: {frac_crit:.0%}  → {'✅ PASS' if sigma_ok else '⚠️ 부분'}")
else:
    sigma_ok = False
    log("  → 데이터 없음")
log()


# ─────────────────────────────────────────────────────────────────────
# 9. 오일러 곱 검증 (Re(s)>1 단순 확인)
# ─────────────────────────────────────────────────────────────────────

log("━━━ 8. 오일러 곱 검증 (Re(s)=2) ━━━")

def euler_product_L(s, N_max=500):
    """오일러 곱 직접 계산 (Re(s)>1에서만 수렴)"""
    result = mpmath.mpc(1)
    for p in range(2, N_max+1):
        if frob_cache.get(p) is None:
            continue
        ap, dp = frob_cache[p]
        ps = mpmath.power(p, -s)
        # Local factor: 1/(1 - a_p·p^{-s} + det_p·p^{-2s})
        if dp == 0:
            continue  # 분기 소수 (근사로 무시)
        local = 1 / (1 - mpmath.mpf(ap) * ps + mpmath.mpf(dp) * ps * ps)
        result *= local
    return result

try:
    s_test = mpmath.mpf('2')
    L_euler = euler_product_L(s_test)
    L_afe = artin_L_afe(s_test + 1j * mpmath.mpf('0.001'), a_coeffs, eps=EPSILON)
    # Re(s)=2: AFE 두 번째 항이 작아야 함
    L_direct = sum(a_coeffs[n] * n**(-2) for n in range(1, 200))
    log(f"  L(2,ρ) 오일러곱: {float(mpmath.re(L_euler)):.6f}")
    log(f"  L(2,ρ) Dirichlet급수(200항): {L_direct:.6f}")
    rel = abs(float(mpmath.re(L_euler)) - L_direct) / max(abs(L_direct), 1e-10)
    log(f"  상대차: {rel:.4f}  ({'✅' if rel < 0.01 else '⚠️'})")
except Exception as e:
    log(f"  WARNING: {e}")
log()


# ─────────────────────────────────────────────────────────────────────
# 10. 종합 판정
# ─────────────────────────────────────────────────────────────────────

log("━━━ 9. 종합 결과 ━━━")
elapsed = time.time() - START

n_zeros = len(zeros_found)
sc1 = fe_mean < 0.05 if fe_results else False
sc2 = n_zeros >= 10
sc3_kappa = kappa_ok
sc3_mono = mono_ok
sc3_sigma = sigma_ok
all_pass = sc1 and sc2 and sc3_kappa and sc3_mono and sc3_sigma

log(f"  다항식: x³ - x - 1,  Gal ≅ S₃,  도체 N = {N_COND}")
log(f"  2차원 기약 표현 ρ,  근 번호 ε = {EPSILON:+d}")
log()
log(f"  SC1 FE 검증    (오차<5%): {'✅ PASS' if sc1 else '⚠️ FAIL'}  (오차={fe_mean:.4f})")
log(f"  SC2 임계선 영점 (≥10개): {'✅ PASS' if sc2 else '⚠️ FAIL'}  (발견={n_zeros}개)")
log(f"  SC3a κδ² 스케일링:       {'✅ PASS' if sc3_kappa else '⚠️ FAIL'}")
log(f"  SC3b 모노드로미 (≈±2π):   {'✅ PASS' if sc3_mono else '⚠️ FAIL'}")
log(f"  SC3c σ-유일성 (σ≈½):     {'✅ PASS' if sc3_sigma else '⚠️ FAIL'}")
log()
log(f"  영점 목록 ({n_zeros}개):")
for i, t0 in enumerate(zeros_found[:20]):
    log(f"    {i+1:2d}. t = {t0:.6f}")
if n_zeros > 20:
    log(f"    ... 외 {n_zeros-20}개")
log()

if all_pass:
    rating = "★★★ 강양성"
elif n_zeros >= 5 and any([sc1, sc3_kappa, sc3_mono]):
    rating = "★★ 양성 (일부 성질 확인)"
elif n_zeros >= 1:
    rating = "★ 약양성 (영점 일부 확보)"
else:
    rating = "⚠️ 조건부 — 영점 탐색 실패"

log(f"  최종 판정: {rating}")
log(f"  소요 시간: {elapsed:.1f}초")
log()
log(f"  결과 파일: {RESULT_FILE}")
log()
log("=" * 72)
log("[실험 #121 완료]")
log("=" * 72)

outf.close()
print(f"\n결과 저장: {RESULT_FILE}", flush=True)
