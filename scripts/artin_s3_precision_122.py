#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #122] Artin S₃ κδ² 고정밀 재측정 — AFE artifact vs structural 규명
=============================================================================
배경:
  #121에서 κδ² slope=-0.901 (이론 -1.0에서 10% 이탈).
  AFE X≈9~12항 (t=5~11 낮은-t 영점)으로 인한 유한항 artifact 의심.
  높은-t 영점(t=40~60, X≈18~21항) + N_MAX 확장 + dps 증가로 artifact 여부 판별.

방법:
  #121과 동일 구조, 정밀도 파라미터만 강화:
  - N_MAX: 600 → 5000
  - dps: 80 → 120
  - δ 범위: [0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10] (더 세밀)
  - κδ² 측정: 낮은-t(t=5~15) 5개 + 높은-t(t=40~60) 5개 분리 측정
  - 교차 검증: 낮은-t vs 높은-t slope 비교

성공 기준:
  1. 높은-t(40~60) slope ≤ -0.98 → AFE artifact 확정, ★★★ 승격
  2. 낮은-t vs 높은-t slope 차이 유의 → artifact 메커니즘 확인
  3. 실패: 높은-t에서도 slope > -0.95 → structural 이탈 → B-29 신설

결과: results/artin_s3_precision_122.txt
=============================================================================
"""
import sys, os, time, math
import numpy as np
import mpmath
from datetime import datetime

mpmath.mp.dps = 120          # 80 → 120
START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'artin_s3_precision_122.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"


log("=" * 72)
log("[실험 #122] Artin S₃ κδ² 고정밀 재측정 — artifact vs structural 판별")
log("=" * 72)
log(f"시작: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log(f"dps = {mpmath.mp.dps}  (이전: 80)")
log(f"N_MAX = 5000  (이전: 600)")
log(f"다항식: x³ - x - 1,  disc = -23,  Gal(K/Q) ≅ S₃")
log(f"도체 N = 23")
log(f"목표: 낮은-t(5~15) vs 높은-t(40~60) κδ² slope 비교로 AFE artifact 여부 확정")
log()


# ─────────────────────────────────────────────────────────────────────
# 1. 오일러 곱 계수 구축 (N_MAX=5000)
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
    log(f"  체 구성 중 (N_max={N_max})...")
    spf = build_sieve(N_max)
    log(f"  {T()} 체 완료. 소수 계수 계산 중...")

    # 소수별 Frobenius 데이터 캐시
    frob_cache = {}
    for i in range(2, N_max+1):
        if spf[i] == i:  # 소수
            frob_cache[i] = frob_type(i)

    log(f"  {T()} Frobenius 캐시: {len(frob_cache)}개 소수")

    a = [0] * (N_max + 1)
    a[1] = 1

    for n in range(2, N_max + 1):
        factors = factorize(n, spf)
        val = 1
        for p, k in factors.items():
            val *= a_pk(p, k, frob_cache)
        a[n] = val
        if n % 500000 == 0:
            log(f"  {T()} 계수 구축 진행: n={n}/{N_max}")

    return a, frob_cache

N_MAX = 5000
log("━━━ 1. Dirichlet 계수 구축 (N_MAX=5000) ━━━")
log(f"  {T()} 시작...")
a_coeffs, frob_cache = build_coeffs(N_MAX)
log(f"  {T()} 완료. a[1..10] = {a_coeffs[1:11]}")

# 주요 소수 통계
split_ps = [p for p in range(2, 100) if frob_cache.get(p) == (2,1)]
trans_ps = [p for p in range(2, 100) if frob_cache.get(p) == (0,-1)]
cycle3_ps = [p for p in range(2, 100) if frob_cache.get(p) == (-1,1)]
log(f"  p<100 분류: identity {len(split_ps)}개, trans {len(trans_ps)}개, 3cycle {len(cycle3_ps)}개")
log(f"  identity 소수(앞 5): {split_ps[:5]}")
log()


# ─────────────────────────────────────────────────────────────────────
# 2. Artin L-함수 (근사 함수방정식)
# ─────────────────────────────────────────────────────────────────────

N_COND = 23  # 도체

def G_factor(s, N_cond=N_COND):
    """
    G(s) = N^{1/2-s} · π^{2s-1} · Γ((1-s)/2)Γ((2-s)/2) / [Γ(s/2)Γ((s+1)/2)]
    """
    pi_fac = mpmath.power(mpmath.pi, 2*s - 1)
    gnum = mpmath.gamma((1 - s) / 2) * mpmath.gamma((2 - s) / 2)
    gden = mpmath.gamma(s / 2) * mpmath.gamma((s + 1) / 2)
    N_fac = mpmath.power(mpmath.mpf(N_cond), mpmath.mpf('0.5') - s)
    return N_fac * pi_fac * gnum / gden

def artin_L_afe(s, a_coeffs, N_cond=N_COND, eps=1, X=None):
    """
    L(s, ρ) via AFE: Σ_{n≤X} a_n n^{-s} + ε·G(s)·Σ_{n≤X} a_n n^{s-1}
    X ≈ √(N·|t|/(2π)) + 6, capped at len(a_coeffs)-1
    """
    t = float(mpmath.im(s))
    if X is None:
        X = max(8, int(math.sqrt(N_cond * max(abs(t), 1.0) / (2 * math.pi))) + 6)
        X = min(X, len(a_coeffs) - 1)

    S1 = mpmath.mpc(0)
    for n in range(1, X + 1):
        if a_coeffs[n] != 0:
            S1 += mpmath.mpf(a_coeffs[n]) * mpmath.power(n, -s)

    S2 = mpmath.mpc(0)
    for n in range(1, X + 1):
        if a_coeffs[n] != 0:
            S2 += mpmath.mpf(a_coeffs[n]) * mpmath.power(n, s - 1)

    G = G_factor(s, N_cond)
    return S1 + mpmath.mpf(eps) * G * S2

def completed_artin(s, a_coeffs, N_cond=N_COND, eps=1, X=None):
    """
    Λ(s, ρ) = N^{s/2} · π^{-(s+1/2)} · Γ(s/2) · Γ((s+1)/2) · L(s)
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
# 4. Z-함수 & 영점 탐색 (t∈[5,60])
# ─────────────────────────────────────────────────────────────────────

def theta_artin(t, N_cond=N_COND):
    """θ(t) = (t/2)·log(N) − t·log(π) + Im logΓ(1/4+it/2) + Im logΓ(3/4+it/2)"""
    t = mpmath.mpf(str(t))
    s = mpmath.mpf('0.5') + 1j * t
    phase = (t * mpmath.log(N_cond) / 2 - t * mpmath.log(mpmath.pi)
             + mpmath.im(mpmath.loggamma(s / 2))
             + mpmath.im(mpmath.loggamma((s + 1) / 2)))
    return float(phase)

def Z_artin(t, a_coeffs, eps=None):
    if eps is None:
        eps = EPSILON
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
    L_val = artin_L_afe(s, a_coeffs, eps=eps)
    theta = theta_artin(t)
    Z = mpmath.exp(1j * theta) * L_val
    return float(mpmath.re(Z)), float(mpmath.im(Z))

log("━━━ 3. 임계선 영점 탐색 (t ∈ [5, 65]) ━━━")

T_MIN, T_MAX, DT = 5.0, 65.0, 0.06  # 더 촘촘한 스캔, t=65까지
ts_scan = np.arange(T_MIN, T_MAX, DT)
Z_re_arr = []
Z_im_arr = []

log(f"  스캔: t∈[{T_MIN},{T_MAX}], Δt={DT}, {len(ts_scan)}점  {T()}")
for idx_t, t in enumerate(ts_scan):
    try:
        zr, zi = Z_artin(float(t), a_coeffs, EPSILON)
        Z_re_arr.append(zr)
        Z_im_arr.append(zi)
    except Exception as e:
        Z_re_arr.append(0.0)
        Z_im_arr.append(0.0)
    if (idx_t+1) % 200 == 0:
        log(f"  {T()} 스캔 {idx_t+1}/{len(ts_scan)}")

Z_main = np.array(Z_re_arr if EPSILON == 1 else Z_im_arr)
log(f"  {T()} 스캔 완료. Z-함수 범위: [{Z_main.min():.3f}, {Z_main.max():.3f}]")

# 부호 변화 → 영점 후보
zeros_raw = []
for i in range(len(Z_main) - 1):
    if Z_main[i] * Z_main[i+1] < 0:
        t1, t2 = float(ts_scan[i]), float(ts_scan[i+1])
        zeros_raw.append((t1, t2))

log(f"  부호변화 {len(zeros_raw)}개 발견")

# 이분법으로 정밀화
zeros_found = []
fail_cnt = 0

def Z_scalar(t_val):
    zr, zi = Z_artin(float(t_val), a_coeffs, EPSILON)
    return zr if EPSILON == 1 else zi

for t1, t2 in zeros_raw:
    try:
        f1 = Z_scalar(t1)
        f2 = Z_scalar(t2)
        if f1 * f2 >= 0:
            fail_cnt += 1
            continue

        a_, b_ = t1, t2
        for _ in range(50):  # 더 많은 반복 (50회)
            mid = (a_ + b_) / 2.0
            fm = Z_scalar(mid)
            if f1 * fm < 0:
                b_ = mid
                f2 = fm
            else:
                a_ = mid
                f1 = fm
            if abs(b_ - a_) < 1e-9:  # 더 높은 정밀도
                break

        t_zero = (a_ + b_) / 2.0

        if not zeros_found or abs(t_zero - zeros_found[-1]) > 0.15:
            zeros_found.append(t_zero)
    except Exception as e:
        fail_cnt += 1
        log(f"  WARNING 이분법 오류: {e}")

log(f"  이분법 완료: {len(zeros_found)}개 영점, 실패 {fail_cnt}건  {T()}")
if len(zeros_found) == 0:
    log("  ⚠️ 영점 0개 — 근 번호 ε 또는 도체 N 설정 점검 필요")
else:
    for i, t0 in enumerate(zeros_found[:20]):
        s0 = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t0))
        L0 = artin_L_afe(s0, a_coeffs, eps=EPSILON)
        # AFE에서 사용된 X 값 출력
        X_afe = max(8, int(math.sqrt(N_COND * max(abs(t0), 1.0) / (2 * math.pi))) + 6)
        X_afe = min(X_afe, len(a_coeffs) - 1)
        log(f"  ρ_{i+1:2d}: t = {t0:.8f},  |L| = {float(abs(L0)):.2e},  AFE X = {X_afe}항")
log()


# ─────────────────────────────────────────────────────────────────────
# 5. FE 검증 (성질 1)
# ─────────────────────────────────────────────────────────────────────

log("━━━ 4. 성질 1: 함수방정식 (FE) 검증 ━━━")

fe_test_pts = [
    mpmath.mpf('0.65') + 1j * mpmath.mpf('8'),
    mpmath.mpf('0.60') + 1j * mpmath.mpf('15'),
    mpmath.mpf('0.55') + 1j * mpmath.mpf('25'),
    mpmath.mpf('0.52') + 1j * mpmath.mpf('45'),
]

fe_results = []
for s in fe_test_pts:
    try:
        Lambda_s   = completed_artin(s,     a_coeffs, eps=EPSILON)
        Lambda_1ms = completed_artin(1 - s, a_coeffs, eps=EPSILON)
        abs_s = float(abs(Lambda_s))
        if abs_s < 1e-30:
            log(f"  WARNING: |Λ(s)|이 너무 작음 ({abs_s:.2e}) — 건너뜀")
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
log(f"  평균 FE 오차: {fe_mean:.4f}  → {'✅ PASS' if fe_ok else '⚠️ FAIL'}")
log()


# ─────────────────────────────────────────────────────────────────────
# 6. κδ² 스케일링 — 핵심 실험 (낮은-t vs 높은-t 비교)
# ─────────────────────────────────────────────────────────────────────

log("━━━ 5. 성질 2: κδ² 스케일링 — 낮은-t vs 높은-t 비교 ━━━")
log("  목적: AFE 유한항 artifact 규명")
log()

def curvature_artin(s, a_coeffs, eps=None):
    """κ = |Λ'/Λ|² (수치 미분, dps=120 고정밀)"""
    if eps is None:
        eps = EPSILON
    h = mpmath.mpf(1) / mpmath.mpf(10**18)  # 더 정밀한 미분 간격
    Lambda_val = completed_artin(s, a_coeffs, eps=eps)
    if abs(Lambda_val) < mpmath.mpf(10)**(-(mpmath.mp.dps - 10)):
        return None  # 영점 근처 — 건너뜀
    Lambda_d = (completed_artin(s + h, a_coeffs, eps=eps) -
                completed_artin(s - h, a_coeffs, eps=eps)) / (2 * h)
    kval = float(abs(Lambda_d / Lambda_val)**2)
    if not math.isfinite(kval) or kval > 1e12:
        return None
    return kval

# δ 범위: 더 세밀화 (수학자 지시)
deltas = [0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10]
log(f"  δ 범위: {deltas}")
log()

def measure_kappa_scaling(zeros_list, label):
    """주어진 영점 목록에 대해 κδ² slope 측정"""
    results = []
    if len(zeros_list) == 0:
        log(f"  [{label}] 영점 없음 — 건너뜀")
        return results

    log(f"  [{label}] 영점 {len(zeros_list)}개 측정 중...")
    for idx, t0 in enumerate(zeros_list):
        # AFE X 값 출력 (artifact 원인 추적)
        X_afe = max(8, int(math.sqrt(N_COND * max(abs(t0), 1.0) / (2 * math.pi))) + 6)
        X_afe = min(X_afe, len(a_coeffs) - 1)

        kappas = []
        valid_deltas = []
        for delta in deltas:
            s_test = mpmath.mpf(str(0.5 + delta)) + 1j * mpmath.mpf(str(t0))
            try:
                k = curvature_artin(s_test, a_coeffs)
                if k is not None and k > 0:
                    kappas.append(k)
                    valid_deltas.append(delta)
            except Exception as e:
                log(f"  WARNING κ 계산 오류 (t={t0:.4f}, δ={delta}): {e}")

        if len(valid_deltas) >= 4:
            log_x = np.log(np.array(valid_deltas)**2)
            log_y = np.log(np.array(kappas))
            slope, intercept = np.polyfit(log_x, log_y, 1)
            ss_res = np.sum((log_y - (slope * log_x + intercept))**2)
            ss_tot = np.sum((log_y - log_y.mean())**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            results.append((t0, slope, r2, X_afe))
            log(f"    ρ (t={t0:.6f}, X={X_afe}항): slope={slope:.4f} (이론=-1), R²={r2:.5f}  [{len(valid_deltas)}점]")
        else:
            log(f"    ρ (t={t0:.6f}, X={X_afe}항): 유효 데이터 부족 ({len(valid_deltas)}점)")

    if results:
        slopes = [r[1] for r in results]
        r2s = [r[2] for r in results]
        log(f"  [{label}] 평균 slope: {np.mean(slopes):.4f} ± {np.std(slopes):.4f}  "
            f"(이론 -1.0, 편차 {abs(np.mean(slopes)+1)*100:.1f}%)")
        log(f"  [{label}] 평균 R²: {np.mean(r2s):.5f}")
    return results

# 낮은-t 영점: t < 20
low_t_zeros = [t for t in zeros_found if t < 20.0][:5]
# 높은-t 영점: t >= 40
high_t_zeros = [t for t in zeros_found if t >= 40.0][:5]

log(f"  낮은-t 영점 (t<20): {len(low_t_zeros)}개 → {[f'{t:.2f}' for t in low_t_zeros]}")
log(f"  높은-t 영점 (t≥40): {len(high_t_zeros)}개 → {[f'{t:.2f}' for t in high_t_zeros]}")
log()

low_t_results  = measure_kappa_scaling(low_t_zeros,  "낮은-t (t<20)")
log()
high_t_results = measure_kappa_scaling(high_t_zeros, "높은-t (t≥40)")
log()

# 교차 검증 판정
log("  ─── 교차 검증 요약 ───")
if low_t_results and high_t_results:
    low_slopes  = [r[1] for r in low_t_results]
    high_slopes = [r[1] for r in high_t_results]
    low_mean    = np.mean(low_slopes)
    high_mean   = np.mean(high_slopes)
    diff        = high_mean - low_mean

    log(f"  낮은-t(t<20) slope:  {low_mean:.4f} ± {np.std(low_slopes):.4f}")
    log(f"  높은-t(t≥40) slope:  {high_mean:.4f} ± {np.std(high_slopes):.4f}")
    log(f"  차이 (높음-낮음):    {diff:+.4f}")

    if high_mean <= -0.98:
        verdict_kappa = "★★★ AFE artifact 확정 — 높은-t에서 slope→이론값"
        kappa_ok = True
    elif high_mean <= -0.95:
        verdict_kappa = "★★ 부분 개선 — 높은-t slope 개선되었으나 이론 미달"
        kappa_ok = True
    else:
        verdict_kappa = "⚠️ structural 이탈 — 높은-t에서도 slope > -0.95 → B-29 신설 필요"
        kappa_ok = False

    if diff < -0.05:
        log(f"  → 낮은-t와 높은-t 간 유의한 slope 차이: AFE 유한항 artifact 메커니즘 확인")
    else:
        log(f"  → 낮은-t와 높은-t 간 slope 차이 미미: structural 이탈 가능성")

    log(f"  판정: {verdict_kappa}")
elif high_t_results:
    high_slopes = [r[1] for r in high_t_results]
    high_mean   = np.mean(high_slopes)
    kappa_ok    = high_mean <= -0.95
    verdict_kappa = f"높은-t slope = {high_mean:.4f}"
    log(f"  높은-t slope: {high_mean:.4f}")
elif low_t_results:
    low_slopes = [r[1] for r in low_t_results]
    low_mean   = np.mean(low_slopes)
    kappa_ok   = abs(low_mean + 1.0) < 0.3
    verdict_kappa = f"낮은-t만 측정: slope = {low_mean:.4f}"
    log(f"  낮은-t slope만 측정됨: {low_mean:.4f}")
else:
    kappa_ok = False
    verdict_kappa = "데이터 없음"
    log("  κδ² 측정 실패 — 영점 부족")

log()


# ─────────────────────────────────────────────────────────────────────
# 7. 모노드로미 (성질 3)
# ─────────────────────────────────────────────────────────────────────

log("━━━ 6. 성질 3: 모노드로미 (폐곡선 적분) ━━━")

def monodromy_artin(t0, a_coeffs, eps=None, radius=0.3, n_steps=64):
    if eps is None:
        eps = EPSILON
    s_center = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t0))
    total_delta = 0.0
    prev_arg = None

    for k in range(n_steps + 1):
        theta_k = 2 * math.pi * k / n_steps
        s = s_center + radius * mpmath.exp(1j * theta_k)
        try:
            val = completed_artin(s, a_coeffs, eps=eps)
        except Exception as e:
            log(f"  WARNING mono step {k}: {e}")
            continue
        if abs(val) < mpmath.mpf(10)**(-(mpmath.mp.dps - 10)):
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

# 낮은-t 5개 + 높은-t 1개 모노드로미
mono_targets = zeros_found[:5] if zeros_found else []
if high_t_zeros:
    mono_targets = list(zeros_found[:4]) + [high_t_zeros[0]] if len(zeros_found) >= 4 else zeros_found[:5]

mono_results = []
for idx, t0 in enumerate(mono_targets):
    try:
        mono = monodromy_artin(t0, a_coeffs)
        mono_results.append(mono)
        near_2pi = abs(abs(mono) - 2*math.pi) < 1.0
        log(f"  ρ (t={t0:.4f}): mono = {mono:.5f} rad  "
            f"({'≈±2π ✅(권선수=1)' if near_2pi else '≠±2π ⚠️'})")
    except Exception as e:
        log(f"  ρ (t={t0:.4f}): 오류 — {e}")

if mono_results:
    frac_2pi = sum(1 for m in mono_results if abs(abs(m) - 2*math.pi) < 1.0) / len(mono_results)
    mono_ok = frac_2pi >= 0.6
    log(f"  ±2π 판정 비율: {frac_2pi:.0%}  → {'✅ PASS' if mono_ok else '⚠️ 부분'}")
else:
    mono_ok = False
    log("  → 모노드로미 데이터 없음")
log()


# ─────────────────────────────────────────────────────────────────────
# 8. σ-유일성 (성질 4) — 전체 영점
# ─────────────────────────────────────────────────────────────────────

log("━━━ 7. 성질 4: σ-유일성 (전체 영점) ━━━")

sigma_scan = np.linspace(0.25, 0.75, 21)
sigma_check_zeros = zeros_found[:10] if len(zeros_found) >= 10 else zeros_found

n_crit = 0
n_total = 0
for idx, t0 in enumerate(sigma_check_zeros):
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
        continue

    min_idx = np.nanargmin(L_arr)
    min_sig = float(sigma_scan[min_idx])
    min_L = float(L_arr[min_idx])
    at_crit = abs(min_sig - 0.5) < 0.08
    n_total += 1
    if at_crit:
        n_crit += 1
    log(f"  ρ (t={t0:.4f}): min|L| = {min_L:.3e} at σ={min_sig:.3f} "
        f"{'(σ≈½ ✅)' if at_crit else '(σ≠½ ⚠️)'}")

if n_total > 0:
    frac_crit = n_crit / n_total
    sigma_ok = frac_crit >= 0.5
    log(f"  임계선 집중 비율: {frac_crit:.0%} ({n_crit}/{n_total}) → {'✅ PASS' if sigma_ok else '⚠️ FAIL'}")
else:
    sigma_ok = False
    log("  데이터 없음")
log()


# ─────────────────────────────────────────────────────────────────────
# 9. 오일러 곱 검증 (Re(s)>1)
# ─────────────────────────────────────────────────────────────────────

log("━━━ 8. 오일러 곱 검증 (Re(s)=2) ━━━")

def euler_product_L(s, N_max=2000):
    result = mpmath.mpc(1)
    for p in sorted(frob_cache.keys()):
        if p > N_max:
            break
        ap, dp = frob_cache[p]
        if dp == 0:
            continue
        ps = mpmath.power(p, -s)
        local = 1 / (1 - mpmath.mpf(ap) * ps + mpmath.mpf(dp) * ps * ps)
        result *= local
    return result

try:
    s_test = mpmath.mpf('2')
    L_euler = euler_product_L(s_test, N_max=2000)
    L_direct = sum(a_coeffs[n] * n**(-2) for n in range(1, min(1000, N_MAX)))
    log(f"  L(2,ρ) 오일러곱(p<2000): {float(mpmath.re(L_euler)):.8f}")
    log(f"  L(2,ρ) Dirichlet급수(1000항): {L_direct:.8f}")
    rel = abs(float(mpmath.re(L_euler)) - L_direct) / max(abs(L_direct), 1e-10)
    log(f"  상대차: {rel:.6f}  ({'✅' if rel < 0.01 else '⚠️'})")
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
log(f"  dps = {mpmath.mp.dps},  N_MAX = {N_MAX}")
log(f"  2차원 기약 표현 ρ,  근 번호 ε = {EPSILON:+d}")
log()
log(f"  SC1 FE 검증    (오차<5%): {'✅ PASS' if sc1 else '⚠️ FAIL'}  (오차={fe_mean:.4f})")
log(f"  SC2 임계선 영점 (≥10개): {'✅ PASS' if sc2 else '⚠️ FAIL'}  (발견={n_zeros}개)")
log(f"  SC3a κδ² 스케일링:       {'✅ PASS' if sc3_kappa else '⚠️ FAIL'}  → {verdict_kappa}")
log(f"  SC3b 모노드로미 (≈±2π):  {'✅ PASS' if sc3_mono else '⚠️ FAIL'}")
log(f"  SC3c σ-유일성 (σ≈½):    {'✅ PASS' if sc3_sigma else '⚠️ FAIL'}")
log()

# κδ² 상세 요약
log("  ─── κδ² 교차 검증 상세 ───")
if low_t_results:
    low_mean_s = np.mean([r[1] for r in low_t_results])
    log(f"  낮은-t(t<20) 평균 slope: {low_mean_s:.4f}  (AFE X ≈ {int(np.mean([r[3] for r in low_t_results]))}항)")
if high_t_results:
    high_mean_s = np.mean([r[1] for r in high_t_results])
    log(f"  높은-t(t≥40) 평균 slope: {high_mean_s:.4f}  (AFE X ≈ {int(np.mean([r[3] for r in high_t_results]))}항)")
if low_t_results and high_t_results:
    log(f"  slope 개선: {low_mean_s:.4f} → {high_mean_s:.4f}  (Δ = {high_mean_s-low_mean_s:+.4f})")
    if high_mean_s <= -0.98:
        log(f"  ✅ AFE artifact 확정 — 높은-t에서 이론값 회복")
    elif high_mean_s <= -0.95:
        log(f"  ⚠️ 부분 개선 — 항수 증가로 slope 개선되었으나 아직 이론 미달")
    else:
        log(f"  🔴 structural 이탈 — AFE 항수와 무관하게 slope > -0.95")
log()

log(f"  영점 목록 (앞 25개 / {n_zeros}개):")
for i, t0 in enumerate(zeros_found[:25]):
    X_afe = max(8, int(math.sqrt(N_COND * max(abs(t0), 1.0) / (2 * math.pi))) + 6)
    X_afe = min(X_afe, len(a_coeffs) - 1)
    tag = " ← 낮은-t" if t0 < 20 else (" ← 높은-t" if t0 >= 40 else "")
    log(f"    {i+1:2d}. t = {t0:.8f}  (AFE X={X_afe}){tag}")
if n_zeros > 25:
    log(f"    ... 외 {n_zeros-25}개")
log()

# 최종 등급
if all_pass and high_t_results and np.mean([r[1] for r in high_t_results]) <= -0.98:
    rating = "★★★ 강양성 — κδ² 이론값 회복 확인, Artin S₃ 결과 ★★★ 승격"
elif sc2 and sc3_mono and sc3_sigma and sc3_kappa:
    rating = "★★ 양성 (높은-t slope 개선, 추가 확인 권고)"
elif sc2 and (sc3_mono or sc3_sigma):
    rating = "★★ 조건부 양성"
elif n_zeros >= 5:
    rating = "★ 약양성"
else:
    rating = "⚠️ 영점 탐색 실패"

log(f"  최종 판정: {rating}")
log(f"  소요 시간: {elapsed:.1f}초")
log()
log(f"  결과 파일: {RESULT_FILE}")
log()
log("=" * 72)
log("[실험 #122 완료]")
log("=" * 72)

outf.close()
print(f"\n결과 저장: {RESULT_FILE}", flush=True)
