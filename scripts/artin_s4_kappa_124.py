#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #124] Artin S₄ κδ² 교차검증 — 낮은-t / 높은-t / 인접 제외
=============================================================================
배경:
  #123에서 slope=-0.9586 (편차 4.1%). 측정 대역 t=20~23, 인접 영점(0.28 간격) 포함.
  → 대역 분리 + 인접 제외로 (A) 간섭 효과 vs (B) degree-3 구조적 한계 판별.

핵심 변경 사항 (#123 대비):
  - dps: 120 → 150
  - δ 범위: [0.005,...,0.1] → [0.003, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05]
  - κδ² 그룹 분리:
      그룹 A: 낮은-t (t<15), 인접 간격>0.5, 최대 고립도 5개
      그룹 B: 높은-t (t≥40), 인접 간격>0.5, 최대 고립도 5개
      그룹 C: 전체 (인접 간격≤0.3 제외)
  - SC1/SC3b/SC3c 제외 (이미 #123에서 확인), κδ² 집중

성공 기준:
  - 낮은-t slope ≤ -0.98 또는 높은-t slope ≤ -0.98 → ★★★ 승격
  - 두 대역 모두 slope > -0.95 → B-30 "degree-3 AFE 구조적 한계" 신설

결과: results/artin_s4_kappa_124.txt
=============================================================================
"""
import sys, os, time, math
import numpy as np
import mpmath
from datetime import datetime

mpmath.mp.dps = 150
START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'artin_s4_kappa_124.txt')
outf = open(RESULT_FILE, 'w', buffering=1)


def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()


def T():
    return f"[{time.time()-START:.1f}s]"


log("=" * 72)
log("[실험 #124] Artin S₄ κδ² 교차검증 — 낮은-t / 높은-t / 인접 제외")
log("=" * 72)
log(f"시작: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log(f"dps = {mpmath.mp.dps}")
log(f"다항식: x⁴ - x - 1,  disc = -283,  Gal(K/Q) ≅ S₄")
log(f"도체 N = 283")
log(f"표준 표현 V_std (3D), 감마: Γ_ℝ(s)²·Γ_ℝ(s+1)")
log(f"δ 범위: [0.003, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05]")
log()


# ─────────────────────────────────────────────────────────────────────
# 1. 다항식 산술 mod p
# ─────────────────────────────────────────────────────────────────────

def poly_mul_mod(a, b, mod_poly, p):
    if not a or not b:
        return [0]
    result = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            result[i + j] = (result[i + j] + ai * bj) % p
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
    result = [1]
    cur = base[:]
    while exp > 0:
        if exp & 1:
            result = poly_mul_mod(result, cur, mod_poly, p)
        cur = poly_mul_mod(cur, cur, mod_poly, p)
        exp >>= 1
    return result


def poly_gcd_mod(a, b, p):
    while True:
        b_nonzero = any(c % p != 0 for c in b) if b else False
        if not b_nonzero:
            break
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
    if a and a[-1] % p != 0:
        inv = pow(a[-1], p - 2, p)
        a = [(c * inv) % p for c in a]
    return a


def split_type_s4(p):
    if p == 283:
        return None
    n_roots = sum(1 for a in range(p) if (pow(a, 4, p) - a - 1) % p == 0)
    if n_roots == 4:
        return (3, 3, 1)
    elif n_roots == 2:
        return (1, -1, -1)
    elif n_roots == 1:
        return (0, 0, 1)
    else:
        f_poly = [(-1) % p, (-1) % p, 0, 0, 1]
        x_poly = [0, 1]
        xp2 = poly_powmod(x_poly, p * p, f_poly, p)
        while len(xp2) < 2:
            xp2.append(0)
        xp2[1] = (xp2[1] - 1) % p
        while len(xp2) > 1 and xp2[-1] % p == 0:
            xp2.pop()
        g = poly_gcd_mod(f_poly[:], xp2, p)
        deg_g = len(g) - 1
        if deg_g >= 4:
            return (-1, -1, 1)
        else:
            return (-1, 1, -1)


# ─────────────────────────────────────────────────────────────────────
# 2. 디리클레 계수 구축
# ─────────────────────────────────────────────────────────────────────

def build_sieve(N_max):
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
    if k == 0:
        return 1
    info = e_cache.get(p)
    if info is None:
        return 0
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
log()


# ─────────────────────────────────────────────────────────────────────
# 3. Artin L-함수 (degree-3 AFE)
# ─────────────────────────────────────────────────────────────────────

N_COND = 283
EPSILON = 1  # 이론: ζ_K=ζ·L(V_std), W=+1 강제


def G_factor_s4(s, N_cond=N_COND):
    N_fac = mpmath.power(mpmath.mpf(N_cond), mpmath.mpf('0.5') - s)
    pi_fac = mpmath.power(mpmath.pi, 3 * s - mpmath.mpf('1.5'))
    gnum = mpmath.gamma((1 - s) / 2) ** 2 * mpmath.gamma((2 - s) / 2)
    gden = mpmath.gamma(s / 2) ** 2 * mpmath.gamma((s + 1) / 2)
    return N_fac * pi_fac * gnum / gden


def artin_L_afe_s4(s, a_coeffs, N_cond=N_COND, eps=1, X=None):
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
    L_val = artin_L_afe_s4(s, a_coeffs, N_cond, eps, X)
    fac = mpmath.power(mpmath.mpf(N_cond) / mpmath.pi ** 3, s / 2)
    gamma_fac = mpmath.gamma(s / 2) ** 2 * mpmath.gamma((s + 1) / 2)
    return fac * gamma_fac * L_val


# ─────────────────────────────────────────────────────────────────────
# 4. Z-함수 & 영점 탐색
# ─────────────────────────────────────────────────────────────────────

def theta_artin_s4(t, N_cond=N_COND):
    t = mpmath.mpf(str(t))
    s_half = mpmath.mpf('0.5') + 1j * t
    phase = (t * (mpmath.log(N_cond) - 3 * mpmath.log(mpmath.pi)) / 2
             + 2 * mpmath.im(mpmath.loggamma(s_half / 2))
             + mpmath.im(mpmath.loggamma((s_half + 1) / 2)))
    return float(phase)


def Z_artin_s4(t, a_coeffs, eps=1):
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
    L_val = artin_L_afe_s4(s, a_coeffs, eps=eps)
    theta = theta_artin_s4(t)
    Z = mpmath.exp(1j * theta) * L_val
    return float(mpmath.re(Z)), float(mpmath.im(Z))


log("━━━ 2. 임계선 영점 탐색 (t ∈ [5, 65]) ━━━")

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
    except Exception as e:
        log(f"  WARNING 스캔 t={t:.2f}: {e}")
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
if len(zeros_raw) == 0:
    log("⚠️ 영점 0개 — 탐색 로직 점검 필요")


def Z_scalar_s4(t_val):
    zr, zi = Z_artin_s4(float(t_val), a_coeffs, EPSILON)
    return zr if EPSILON == 1 else zi


# 이분법 정밀화
zeros_found = []
fail_cnt = 0

for t1, t2 in zeros_raw:
    try:
        f1 = Z_scalar_s4(t1)
        f2 = Z_scalar_s4(t2)
        if f1 * f2 >= 0:
            fail_cnt += 1
            continue
        lo, hi = t1, t2
        flo = f1
        for _ in range(50):
            mid = (lo + hi) / 2
            fm = Z_scalar_s4(mid)
            if fm * flo < 0:
                hi = mid
            else:
                lo = mid
                flo = fm
        zeros_found.append((lo + hi) / 2)
    except Exception as e:
        log(f"  WARNING 이분법 실패 ({t1:.2f}~{t2:.2f}): {e}")
        fail_cnt += 1

log(f"  이분법 완료: {len(zeros_found)}개 영점, {fail_cnt}개 실패  {T()}")

if fail_cnt > len(zeros_raw) / 2:
    log(f"  ⚠️ 이분법 실패율 {fail_cnt}/{len(zeros_raw)} — 중단")
    outf.close()
    sys.exit(1)

zeros_found.sort()

# 인접 간격 계산
n_z = len(zeros_found)
min_spacings = []
for i, t0 in enumerate(zeros_found):
    prev_gap = t0 - zeros_found[i - 1] if i > 0 else 999.0
    next_gap = zeros_found[i + 1] - t0 if i < n_z - 1 else 999.0
    min_spacings.append(min(prev_gap, next_gap))

log(f"  영점 목록 (전체 {n_z}개):")
for i, (t0, sp) in enumerate(zip(zeros_found, min_spacings)):
    afe_X = max(20, int(math.sqrt(N_COND * max(abs(t0), 1.0) / (2 * math.pi))) + 10)
    afe_X = min(afe_X, N_MAX)
    flag = " ⚠️인접" if sp < 0.5 else ""
    log(f"  {i+1:3d}. t = {t0:.8f}  간격={sp:.3f}  AFE X={afe_X}{flag}")
log()


# ─────────────────────────────────────────────────────────────────────
# 5. 그룹 선택
# ─────────────────────────────────────────────────────────────────────

log("━━━ 3. κδ² 측정 그룹 선택 ━━━")

# 그룹 A: 낮은-t (t<15), 고립도 순 정렬, 상위 5개 (최소 간격>0.5 우선)
low_t_candidates = [(t0, sp) for t0, sp in zip(zeros_found, min_spacings) if t0 < 15.0]
low_t_candidates.sort(key=lambda x: -x[1])  # 고립도 내림차순
group_A = [t0 for t0, sp in low_t_candidates[:5]]

# 그룹 B: 높은-t (t≥40), 고립도 순 정렬, 상위 5개 (최소 간격>0.5 우선)
high_t_candidates = [(t0, sp) for t0, sp in zip(zeros_found, min_spacings) if t0 >= 40.0]
high_t_candidates.sort(key=lambda x: -x[1])
group_B = [t0 for t0, sp in high_t_candidates[:5]]

# 그룹 C: 전체 (인접 간격≤0.3 제외)
group_C = [t0 for t0, sp in zip(zeros_found, min_spacings) if sp > 0.3]

log(f"  그룹 A (낮은-t, t<15): {len(group_A)}개")
for t0, sp in low_t_candidates[:5]:
    log(f"    t={t0:.6f}  간격={sp:.3f}")

log(f"  그룹 B (높은-t, t≥40): {len(group_B)}개")
for t0, sp in high_t_candidates[:5]:
    log(f"    t={t0:.6f}  간격={sp:.3f}")

log(f"  그룹 C (전체 비-인접, 간격>0.3): {len(group_C)}개")

if len(group_A) == 0:
    log("⚠️ 낮은-t 영점 없음 — 범위 조정 필요")
if len(group_B) == 0:
    log("⚠️ 높은-t 영점 없음 — 범위 조정 필요")
log()


# ─────────────────────────────────────────────────────────────────────
# 6. κδ² 측정 함수
# ─────────────────────────────────────────────────────────────────────

# #124 핵심: 더 세밀한 δ + dps=150
DELTAS = [0.003, 0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05]

log(f"━━━ 4. κδ² 스케일링 측정 (δ={DELTAS}) ━━━")
log(f"  dps={mpmath.mp.dps} (150으로 강화)")
log()


def curvature_s4(s, a_coeffs, eps=EPSILON):
    """κ = |Λ'/Λ|² at 점 s. dps=150 사용."""
    h = mpmath.mpf(1) / mpmath.mpf(10 ** 25)  # dps=150에 맞춰 h 더 작게
    Lambda_val = completed_artin_s4(s, a_coeffs, eps=eps)
    if abs(Lambda_val) < mpmath.mpf(10) ** (-(mpmath.mp.dps - 10)):
        return None  # 영점 너무 가까움
    Lambda_d = (completed_artin_s4(s + h, a_coeffs, eps=eps) -
                completed_artin_s4(s - h, a_coeffs, eps=eps)) / (2 * h)
    ratio = abs(Lambda_d / Lambda_val)
    kval = float(ratio ** 2)
    if not math.isfinite(kval) or kval > 1e20:
        return None
    return kval


def measure_kappa_group(group_zeros, group_name):
    """지정 영점 그룹에 대해 κδ² 측정 후 slope 반환."""
    log(f"  [{group_name}] {len(group_zeros)}개 영점 측정 시작 {T()}")
    slopes = []
    r2_vals = []

    for t0 in group_zeros:
        kappas = []
        valid_deltas = []
        for delta in DELTAS:
            s_test = mpmath.mpf(str(0.5 + delta)) + 1j * mpmath.mpf(str(t0))
            try:
                k = curvature_s4(s_test, a_coeffs)
                if k is not None and k > 0:
                    kappas.append(k)
                    valid_deltas.append(delta)
                else:
                    log(f"    ρ(t={t0:.4f}), δ={delta}: κ=None (영점 근접?)")
            except Exception as e:
                log(f"    WARNING ρ(t={t0:.4f}), δ={delta}: {e}")

        if len(valid_deltas) >= 4:
            log_x = np.log(np.array(valid_deltas) ** 2)
            log_y = np.log(np.array(kappas))
            slope, intercept = np.polyfit(log_x, log_y, 1)
            ss_res = np.sum((log_y - (slope * log_x + intercept)) ** 2)
            ss_tot = np.sum((log_y - log_y.mean()) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 1e-30 else 0
            slopes.append(slope)
            r2_vals.append(r2)
            afe_X = max(20, int(math.sqrt(N_COND * max(abs(t0), 1.0) / (2 * math.pi))) + 10)
            afe_X = min(afe_X, N_MAX)
            dev = abs(slope + 1) * 100
            log(f"    ρ(t={t0:.6f}, X={afe_X}항): slope={slope:.4f}  R²={r2:.5f}  편차={dev:.1f}%  [{len(valid_deltas)}점]")
        else:
            log(f"    ρ(t={t0:.6f}): 유효 δ 부족 ({len(valid_deltas)}점) — 건너뜀")

    if slopes:
        mean_s = np.mean(slopes)
        std_s = np.std(slopes)
        mean_r2 = np.mean(r2_vals)
        dev_pct = abs(mean_s + 1) * 100
        if dev_pct <= 2.0:
            verdict = "✅ ★★★ 기준 달성 (≤2%)"
        elif dev_pct <= 5.0:
            verdict = "✅ PASS (≤5%)"
        elif dev_pct <= 10.0:
            verdict = "⚠️ MARGINAL (≤10%)"
        else:
            verdict = "❌ FAIL (>10%)"
        log(f"  [{group_name}] 평균 slope: {mean_s:.4f} ± {std_s:.4f}  편차={dev_pct:.1f}%  R²={mean_r2:.5f}")
        log(f"  [{group_name}] 판정: {verdict}")
        log()
        return mean_s, std_s, dev_pct, slopes
    else:
        log(f"  [{group_name}] ⚠️ 측정 실패 — 영점 없음")
        log()
        return None, None, None, []


# ─────────────────────────────────────────────────────────────────────
# 7. 그룹별 κδ² 측정
# ─────────────────────────────────────────────────────────────────────

log("━━━ 4-A. 그룹 A: 낮은-t (t<15) ━━━")
res_A = measure_kappa_group(group_A, "낮은-t A")

log("━━━ 4-B. 그룹 B: 높은-t (t≥40) ━━━")
res_B = measure_kappa_group(group_B, "높은-t B")

log("━━━ 4-C. 그룹 C: 전체 비-인접 (간격>0.3) ━━━")
# 그룹 C에서 최대 20개만 측정 (시간 절약)
group_C_sample = group_C[:20]
res_C = measure_kappa_group(group_C_sample, "전체 C")


# ─────────────────────────────────────────────────────────────────────
# 8. 종합 결과 및 판정
# ─────────────────────────────────────────────────────────────────────

log("=" * 72)
log("━━━ 5. 종합 결과 & 판정 ━━━")
log()
log(f"  실험: #124 — S₄ κδ² 교차검증 (낮은-t / 높은-t / 인접 제외)")
log(f"  dps={mpmath.mp.dps}, δ={DELTAS}")
log(f"  영점 총 {n_z}개 발견 (t∈[5,65])")
log()
log("  ┌─────────────────────────────────────────────────────────────┐")
log("  │ 그룹           │ mean slope  │ ±std    │ 편차   │ 판정     │")
log("  ├─────────────────────────────────────────────────────────────┤")

def fmt_row(name, res):
    ms, ss, dp, sl = res
    if ms is None:
        return f"  │ {name:<14s} │ {'N/A':>11s} │ {'N/A':>7s} │ {'N/A':>6s} │ {'⚠️':>8s} │"
    judge = "★★★" if dp <= 2.0 else ("PASS" if dp <= 5.0 else ("MARGIN" if dp <= 10.0 else "FAIL"))
    return (f"  │ {name:<14s} │ {ms:>+11.4f} │ ±{ss:.4f} │ {dp:>5.1f}% │ {judge:>8s} │")

log(fmt_row("낮은-t (A)", res_A))
log(fmt_row("높은-t (B)", res_B))
log(fmt_row("전체비인접 (C)", res_C))
log("  └─────────────────────────────────────────────────────────────┘")
log()

# 최종 판정
star_upgrade = False
b30_structural = False

ms_A = res_A[0]
ms_B = res_B[0]
ms_C = res_C[0]

if ms_A is not None and ms_A <= -0.98:
    star_upgrade = True
    log(f"  ✅ 낮은-t slope {ms_A:.4f} ≤ -0.98 → ★★★ 승격 조건 충족")
if ms_B is not None and ms_B <= -0.98:
    star_upgrade = True
    log(f"  ✅ 높은-t slope {ms_B:.4f} ≤ -0.98 → ★★★ 승격 조건 충족")
if ms_C is not None and ms_C <= -0.98:
    star_upgrade = True
    log(f"  ✅ 전체(비인접) slope {ms_C:.4f} ≤ -0.98 → ★★★ 승격 조건 충족")

# B-30: 두 대역 모두 > -0.95
A_weak = (ms_A is None or ms_A > -0.95)
B_weak = (ms_B is None or ms_B > -0.95)
if A_weak and B_weak:
    b30_structural = True
    log(f"  ⚠️ 두 대역 모두 slope > -0.95 → B-30 'degree-3 AFE 구조적 한계' 신설 필요")

log()
if star_upgrade:
    log("  ★★★ 최종 판정: S₄ ★★★ 승격 → Paper 3 강양성 소재 확보")
elif b30_structural:
    log("  ★★ 최종 판정: S₄ ★★ 확정 + B-30 신설 (degree-3 구조적 한계)")
else:
    ms_best = max(filter(lambda x: x is not None, [ms_A, ms_B, ms_C]), default=None)
    if ms_best is not None:
        dp = abs(ms_best + 1) * 100
        log(f"  ★★☆ 최종 판정: 최선 slope={ms_best:.4f} (편차 {dp:.1f}%) — 추가 분석 필요")
    else:
        log("  ❌ 판정 불가 — 모든 그룹 측정 실패")

log()
log(f"  총 소요 시간: {time.time()-START:.1f}초")
log("=" * 72)
outf.close()
