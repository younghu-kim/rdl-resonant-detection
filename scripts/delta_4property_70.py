#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #70 — Ramanujan Δ L-function (N=1, weight=12) 4성질 검증
=============================================================================
대상: L(Δ, s) — Ramanujan Δ modular form의 L-함수
  - weight k=12, level N=1, self-dual (ε=+1)
  - 감마 인자: (2π)^{-s} Γ(s + 11/2)
  - 함수방정식: Λ(s) = Λ(1-s)

Fourier 계수: Ramanujan τ(n) 함수
  - Δ(q) = q ∏_{n≥1} (1-q^n)^{24} = Σ τ(n) q^n
  - Hecke 정규화: a(n) = τ(n) / n^{11/2}
  - Ramanujan 추측 (Deligne 정리): |a(p)| ≤ 2

σ-유일성 판별 실험 (경계 B-10):
  - N=1이면서 weight=12
  - PASS → 가설 A (N=1 → PASS) 확정
  - FAIL → 가설 B (weight=0 → PASS) 확정
  - 어느 쪽이든 정보적 결과

결과 파일: results/delta_4property_70.txt
=============================================================================
"""

import sys, os, time
import numpy as np
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

DPS = 50
mpmath.mp.dps = DPS

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "delta_4property_70.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

# ━━━━━━━━━━━ 로깅 ━━━━━━━━━━━
lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))
def flush():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━━━━━━ 파라미터 ━━━━━━━━━━━
WEIGHT = 12
HALF_WEIGHT_MINUS_1 = mpmath.mpf(11) / 2  # (k-1)/2 = 11/2
EPSILON = 1    # self-dual, (-1)^{k/2} = (-1)^6 = 1
N_COND = 1     # conductor

# Gauss-Hermite 직교
N_HERMITE = 50
_herm_nodes, _herm_weights = np.polynomial.hermite.hermgauss(N_HERMITE)
HERM_NODES = [mpmath.mpf(float(x)) for x in _herm_nodes]
HERM_WEIGHTS = [mpmath.mpf(float(w)) for w in _herm_weights]

CONTOUR_C = mpmath.mpf(2)  # AFE contour shift

# κ 스윕 파라미터
T_MIN, T_MAX = 2.0, 55.0   # 넓은 범위 (첫 영점 γ₁≈9.22)
DT_SCAN = 0.35              # 영점 스캔 간격
DELTA = 0.03                # σ 오프셋
MONO_R = 0.4                # 모노드로미 반지름
MONO_N = 48                 # 모노드로미 단계

N_MAX_COEFF = 500           # Dirichlet series 항 수

# ━━━━━━━━━━━ Ramanujan τ(n) 계산 ━━━━━━━━━━━

def compute_tau(N_max):
    """
    Δ(q) = q ∏_{n≥1} (1-q^n)^{24} = Σ_{n≥1} τ(n) q^n

    q-product expansion으로 τ(n) 계산 (n=1..N_max).
    정수 연산으로 정확히 계산.
    """
    # f(q) = ∏_{n=1}^{N_max} (1-q^n)^{24}
    # Δ(q) = q · f(q), 따라서 τ(n) = f의 q^{n-1} 계수

    coeffs = [0] * (N_max + 1)
    coeffs[0] = 1

    for n in range(1, N_max + 1):
        # (1-q^n)^{24}를 곱함
        # 효율: (1-x)^{24} = Σ_{k=0}^{24} C(24,k)(-1)^k x^k
        from math import comb
        binom_coeffs = [((-1)**k) * comb(24, k) for k in range(25)]

        # 역순으로 갱신 (in-place)
        for m in range(N_max, 0, -1):
            for k in range(1, min(24, m // n) + 1):
                coeffs[m] += binom_coeffs[k] * coeffs[m - k * n]

    tau = {}
    for n in range(1, N_max + 1):
        tau[n] = coeffs[n - 1]

    return tau


def verify_tau(tau):
    """τ(n) 검증: 알려진 값과 비교"""
    known = {1: 1, 2: -24, 3: 252, 4: -1472, 5: 4830,
             6: -6048, 7: -16744, 8: 84480, 9: -113643, 10: -115920}
    ok = True
    for n, expected in known.items():
        if n in tau:
            if tau[n] != expected:
                log(f"  ❌ τ({n}) = {tau[n]}, 기대값 = {expected}")
                ok = False
            else:
                log(f"  ✅ τ({n}) = {tau[n]}")
        else:
            log(f"  ⚠️ τ({n}) 없음")
            ok = False
    return ok


def verify_hecke(tau):
    """Hecke 곱셈성/재귀 검증"""
    checks = []

    # τ(mn) = τ(m)τ(n) for gcd(m,n)=1
    from math import gcd
    pairs = [(2, 3), (2, 5), (3, 5), (2, 7), (3, 7), (2, 11)]
    for m, n in pairs:
        if m*n <= max(tau.keys()):
            lhs = tau[m * n]
            rhs = tau[m] * tau[n]
            diff = abs(lhs - rhs)
            checks.append((f"τ({m}·{n})=τ({m})τ({n})", diff))
            log(f"  곱셈성: τ({m*n})={lhs}, τ({m})·τ({n})={rhs}, diff={diff}")

    # τ(p^{k+1}) = τ(p)·τ(p^k) - p^{11}·τ(p^{k-1})
    for p in [2, 3, 5]:
        for k in range(1, 4):
            pk1 = p**(k+1)
            pk = p**k
            pk_1 = p**(k-1)
            if pk1 <= max(tau.keys()):
                lhs = tau[pk1]
                rhs = tau[p] * tau[pk] - (p**11) * tau[pk_1]
                diff = abs(lhs - rhs)
                checks.append((f"τ({p}^{k+1})재귀", diff))
                log(f"  Hecke 재귀: τ({pk1})={lhs}, τ({p})τ({pk})-{p}^11·τ({pk_1})={rhs}, diff={diff}")

    all_ok = all(d == 0 for _, d in checks)
    return all_ok


# ━━━━━━━━━━━ 감마 인자 (Ramanujan Δ) ━━━━━━━━━━━

def gamma_factor_delta(s):
    """
    Δ L-함수 감마 인자 (analytic normalization):
    γ(s) = (2π)^{-s} · Γ(s + 11/2)

    함수방정식: Λ(s) = γ(s)·L(s) = Λ(1-s)
    여기서 L(s) = Σ a(n)/n^s, a(n) = τ(n)/n^{11/2}
    """
    return mpmath.power(2 * mpmath.pi, -s) * mpmath.gamma(s + HALF_WEIGHT_MINUS_1)


def dirichlet_series(w, a_coeffs, N_max):
    """D(w) = Σ_{n=1}^{N_max} a(n)/n^w"""
    total = mpmath.mpc(0)
    for n in range(1, N_max + 1):
        if n not in a_coeffs:
            continue
        an = a_coeffs[n]
        if an == 0:
            continue
        total += an * mpmath.power(n, -w)
    return total


def Lambda_AFE(s, a_coeffs, N_max):
    """
    Λ(s) via Mellin-Barnes AFE + Gauss-Hermite.
    Λ(s) = (2π)^{-s} Γ(s+11/2) L(s)
    함수방정식: Λ(s) = Λ(1-s) (ε=+1)
    """
    c = CONTOUR_C
    s_mp = mpmath.mpc(s)
    s1_mp = 1 - s_mp

    prefactor = mpmath.exp(c ** 2) / (2 * mpmath.pi)

    total = mpmath.mpc(0)
    for k in range(N_HERMITE):
        v_k = HERM_NODES[k]
        wt_k = HERM_WEIGHTS[k]

        iv_k = mpmath.mpc(0, v_k)
        w_shift = c + iv_k

        exp_phase = mpmath.exp(2 * mpmath.mpc(0, 1) * c * v_k)

        # A_k(s): gamma at s+w_shift
        sw = s_mp + w_shift
        gamma_sw = gamma_factor_delta(sw)
        A_s = gamma_sw * exp_phase / w_shift
        D_s = dirichlet_series(sw, a_coeffs, N_max)

        # A_k(1-s)
        s1w = s1_mp + w_shift
        gamma_s1w = gamma_factor_delta(s1w)
        A_1s = gamma_s1w * exp_phase / w_shift
        D_1s = dirichlet_series(s1w, a_coeffs, N_max)

        total += wt_k * (A_s * D_s + EPSILON * A_1s * D_1s)

    return prefactor * total


def curvature(t, a_coeffs, N_max, delta=DELTA, h=1e-6):
    """κ(s) = |Λ'/Λ|² at s = 0.5+δ+it"""
    s = mpmath.mpc(0.5 + delta, t)
    h_mp = mpmath.mpf(h)
    try:
        L0 = Lambda_AFE(s, a_coeffs, N_max)
        if abs(L0) < mpmath.mpf(10) ** (-DPS + 15):
            return float('inf')
        Lp = Lambda_AFE(s + h_mp, a_coeffs, N_max)
        Lm = Lambda_AFE(s - h_mp, a_coeffs, N_max)
        conn = (Lp - Lm) / (2 * h_mp * L0)
        k = float(abs(conn) ** 2)
        return k if np.isfinite(k) else float('inf')
    except Exception as e:
        print(f"  WARNING curvature at t={t}: {e}", flush=True)
        return 0.0


def monodromy(t_center, a_coeffs, N_max, radius=MONO_R, n_steps=MONO_N):
    """폐곡선 적분으로 모노드로미 측정"""
    center = mpmath.mpc(0.5, t_center)
    phase_accum = mpmath.mpf(0)
    prev_val = None

    for j in range(n_steps + 1):
        theta = 2 * mpmath.pi * j / n_steps
        pt = center + radius * mpmath.exp(1j * theta)
        try:
            val = Lambda_AFE(pt, a_coeffs, N_max)
            if abs(val) < mpmath.mpf(10) ** (-DPS + 15):
                return None
        except Exception as e:
            return None

        if prev_val is not None:
            ratio = val / prev_val
            if abs(ratio) < mpmath.mpf(10) ** (-DPS + 15):
                return None
            phase_accum += mpmath.im(mpmath.log(ratio))
        prev_val = val

    return float(abs(phase_accum) / mpmath.pi)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#                   메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 70)
log("결과 #70 — Ramanujan Δ L-함수 (N=1, weight=12) 4성질 검증")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"DPS={DPS}, N_hermite={N_HERMITE}, c={float(CONTOUR_C)}")
log(f"weight k={WEIGHT}, N={N_COND}, ε={EPSILON}")
log(f"감마: (2π)^{{-s}} Γ(s + 11/2)")
log(f"함수방정식: Λ(s) = Λ(1-s)")
log(f"Dirichlet series 항 수: N_max={N_MAX_COEFF}")
log()

t_total = time.time()

# ━━━━━━ Step 0: τ(n) 계산 및 검증 ━━━━━━

log("[Step 0] Ramanujan τ(n) 계산 (q-product expansion)")
t0 = time.time()
tau = compute_tau(N_MAX_COEFF)
dt0 = time.time() - t0
log(f"  τ(n) 계산 완료: {len(tau)}개 ({dt0:.1f}s)")
log()

log("  τ(n) 알려진 값 검증:")
tau_ok = verify_tau(tau)
if not tau_ok:
    log("  ❌ τ(n) 검증 실패 — 중단")
    flush()
    sys.exit(1)
log(f"  ✅ τ(n) 10개 검증 통과")
log()

log("  Hecke 곱셈성/재귀 검증:")
hecke_ok = verify_hecke(tau)
if not hecke_ok:
    log("  ❌ Hecke 검증 실패 — 중단")
    flush()
    sys.exit(1)
log(f"  ✅ Hecke 곱셈성/재귀 전부 통과")
log()

# Hecke 정규화 a(n) = τ(n) / n^{11/2}
log("  Hecke 정규화: a(n) = τ(n) / n^{11/2}")
a_coeffs = {}
for n in range(1, N_MAX_COEFF + 1):
    if tau[n] != 0:
        a_coeffs[n] = mpmath.mpf(tau[n]) / mpmath.power(n, HALF_WEIGHT_MINUS_1)

# Ramanujan 추측 검증: |a(p)| ≤ 2 for primes
primes_small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
log("  Ramanujan-Petersson 추측 검증 (|a(p)| ≤ 2):")
for p in primes_small:
    ap = float(abs(a_coeffs[p]))
    ok = ap <= 2.0
    log(f"    a({p}) = {float(a_coeffs[p]):+.8f}, |a({p})| = {ap:.8f} {'✅' if ok else '❌'}")
log()
flush()

# ━━━━━━ Step 1: 함수방정식 검증 ━━━━━━

log("[Step 1] 함수방정식 Λ(s) = Λ(1-s) on critical line")
fe_pts = [
    mpmath.mpc(0.5, 5),
    mpmath.mpc(0.5, 10),
    mpmath.mpc(0.5, 15),
    mpmath.mpc(0.5, 25),
    mpmath.mpc(0.5, 40),
]
fe_ok = True
fe_max_rel = 0.0

for sp in fe_pts:
    t_s = time.time()
    L_s = Lambda_AFE(sp, a_coeffs, N_MAX_COEFF)
    L_1s = Lambda_AFE(1 - sp, a_coeffs, N_MAX_COEFF)
    dt_s = time.time() - t_s

    if abs(L_s) > mpmath.mpf(10) ** (-DPS + 15):
        rel = float(abs(L_s - EPSILON * L_1s) / abs(L_s))
    else:
        rel = 0.0  # 영점 근방

    fe_max_rel = max(fe_max_rel, rel)
    ok = rel < 1e-6
    if not ok:
        fe_ok = False
    log(f"  s={mpmath.nstr(sp, 6)}: |Λ|={float(abs(L_s)):.4e}, rel={rel:.2e} {'✅' if ok else '❌'} ({dt_s:.1f}s)")
    flush()

if fe_ok:
    log(f"  ✅ 함수방정식 통과 (max_rel={fe_max_rel:.2e})")
else:
    log(f"  ❌ 함수방정식 실패 (max_rel={fe_max_rel:.2e})")
    log(f"  경고: 감마 인자 (2π)^{{-s}} Γ(s+11/2) 재확인 필요.")
log()
flush()

# ━━━━━━ Step 2: Λ(1/2+it) 스캔 — 영점 탐색 ━━━━━━

log(f"[Step 2] Λ(1/2+it) 스캔 (t ∈ [{T_MIN}, {T_MAX}])")
log(f"  기대: 첫 영점 γ₁ ≈ 9.22, 12개 이상 확보 목표")
scan_ts = np.arange(T_MIN, T_MAX + DT_SCAN, DT_SCAN)
scan_vals = []
t_scan = time.time()

for i, t in enumerate(scan_ts):
    val = Lambda_AFE(mpmath.mpc(0.5, t), a_coeffs, N_MAX_COEFF)
    re_val = float(mpmath.re(val))
    scan_vals.append(re_val)

    if (i + 1) % 20 == 0:
        elapsed = time.time() - t_scan
        eta = elapsed / (i + 1) * (len(scan_ts) - i - 1)
        log(f"  [{i+1}/{len(scan_ts)}] t={t:.1f}: Re(Λ)={re_val:+.4e} ({elapsed:.0f}s, ETA {eta:.0f}s)")
        flush()

scan_arr = np.array(scan_vals)
dt_scan_elapsed = time.time() - t_scan
log(f"  스캔 완료: {dt_scan_elapsed:.1f}s ({dt_scan_elapsed/60:.1f}분)")

# ━━━━━━ Step 3: 영점 이분법 정밀화 ━━━━━━

zero_approx = []
for i in range(len(scan_arr) - 1):
    if scan_arr[i] * scan_arr[i + 1] < 0:
        t_lo, t_hi = float(scan_ts[i]), float(scan_ts[i + 1])
        val_lo = scan_arr[i]
        try:
            for _ in range(25):  # 25회 이분법
                t_mid = (t_lo + t_hi) / 2
                val_mid = float(mpmath.re(Lambda_AFE(mpmath.mpc(0.5, t_mid), a_coeffs, N_MAX_COEFF)))
                if not np.isfinite(val_mid):
                    break
                if val_mid * val_lo < 0:
                    t_hi = t_mid
                else:
                    t_lo = t_mid
                    val_lo = val_mid
            z_t = (t_lo + t_hi) / 2
            zero_approx.append(z_t)
        except Exception as e:
            log(f"  WARNING bisection at t≈{(t_lo+t_hi)/2:.3f}: {e}")
            zero_approx.append((t_lo + t_hi) / 2)

n_zeros = len(zero_approx)
log(f"\n  발견 영점: {n_zeros}개")
if n_zeros == 0:
    log(f"  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
    flush()
else:
    for j, z in enumerate(zero_approx):
        log(f"    #{j+1}: γ = {z:.10f}")

# LMFDB 알려진 영점 비교 (처음 몇 개)
known_zeros = [9.2224, 13.9076, 17.4428, 19.6565, 22.3350, 25.2748]
log(f"\n  LMFDB 알려진 영점 (근사)과 비교:")
for kz in known_zeros:
    match = min(zero_approx, key=lambda z: abs(z - kz)) if zero_approx else None
    if match and abs(match - kz) < 1.0:
        log(f"    γ≈{kz}: 매칭 → {match:.10f} (차이 {abs(match-kz):.6f})")
    else:
        log(f"    γ≈{kz}: 매칭 실패 ⚠️")

log()
flush()

# ━━━━━━ Step 4: κ near/far 측정 + 모노드로미 ━━━━━━

near_med = 0
far_med = 0
ratio = 0
near_cv = 999
mono_pass = 0
total_mono = 0
near_fin = []
kappa_near = []
mono_results = []

if n_zeros >= 2:
    log(f"[Step 3] κ near/far 비율 + 모노드로미")
    n_tp = min(12, n_zeros)  # 최대 12 TP (≥6 목표)

    for j, z_t in enumerate(zero_approx[:n_tp]):
        t_k = time.time()
        k_val = curvature(z_t, a_coeffs, N_MAX_COEFF, delta=DELTA)
        m_val = monodromy(z_t, a_coeffs, N_MAX_COEFF)
        dt_k = time.time() - t_k
        kappa_near.append(k_val)
        mono_results.append(m_val)
        m_str = f"{m_val:.4f}" if m_val is not None else "FAIL"
        log(f"    TP #{j+1} γ={z_t:.8f}: κ={k_val:.2f}, mono/π={m_str} ({dt_k:.1f}s)")
        flush()

    # FP (영점 사이 중간점)
    kappa_far = []
    fp_ts = []
    for j in range(len(zero_approx) - 1):
        mid = (zero_approx[j] + zero_approx[j + 1]) / 2
        fp_ts.append(mid)
    if zero_approx[0] > T_MIN + 1:
        fp_ts.insert(0, T_MIN + 0.5)

    n_fp = min(8, len(fp_ts))
    for j, t in enumerate(fp_ts[:n_fp]):
        t_k = time.time()
        k_val = curvature(t, a_coeffs, N_MAX_COEFF, delta=DELTA)
        dt_k = time.time() - t_k
        kappa_far.append(k_val)
        log(f"    FP #{j+1} t={t:.6f}: κ={k_val:.4f} ({dt_k:.1f}s)")
        flush()

    # κ 통계
    near_fin = [k for k in kappa_near if np.isfinite(k) and k < 1e15 and k > 0]
    far_fin = [k for k in kappa_far if np.isfinite(k) and k < 1e15 and k > 0]

    if near_fin and far_fin:
        near_med = float(np.median(near_fin))
        far_med = float(np.median(far_fin))
        ratio = near_med / far_med if far_med > 0 else float('inf')
        near_cv = float(np.std(near_fin) / np.mean(near_fin) * 100) if np.mean(near_fin) > 0 else 999
    else:
        near_med = 0
        far_med = 0
        ratio = 0
        near_cv = 999

    log(f"\n  κ near median = {near_med:.2f} (n={len(near_fin)}, CV={near_cv:.1f}%)")
    log(f"  κ far median  = {far_med:.4f} (n={len(far_fin)})")
    log(f"  κ ratio = {ratio:.1f}×")

    # 모노드로미 요약
    mono_pass = sum(1 for m in mono_results if m is not None and 1.5 < m < 2.5)
    total_mono = len([m for m in mono_results if m is not None])
    log(f"  모노드로미: {mono_pass}/{total_mono} (mono/π≈2.0)")

    # 상세 모노드로미
    for j, m in enumerate(mono_results):
        m_str = f"{m:.4f}" if m is not None else "FAIL"
        log(f"    TP #{j+1}: mono/π = {m_str}")

    log()
    flush()

# ━━━━━━ Step 5: σ-유일성 테스트 (≥12 영점) ━━━━━━

log("[Step 4] σ-유일성 테스트 (7σ × 영점, 목표 ≥12개)")
log("  이것이 B-10 판별의 핵심:")
log("  PASS → 가설 A (N=1→PASS) 확정")
log("  FAIL → 가설 B (weight=0→PASS) 확정")
log()
sigma_pass_count = 0
sigma_fail_count = 0
sigma_test_count = 0
sigma_details = []

if n_zeros >= 2:
    sigma_vals = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
    n_sigma_test = min(max(12, n_zeros), n_zeros)  # 최소 12개
    test_zeros = zero_approx[:n_sigma_test]

    for z_idx, z_t in enumerate(test_zeros):
        log(f"\n  영점 #{z_idx+1} γ={z_t:.6f}:")
        kappas_by_sigma = []
        for sigma in sigma_vals:
            k_val = curvature(z_t, a_coeffs, N_MAX_COEFF, delta=sigma - 0.5)
            kappas_by_sigma.append(k_val)
            log(f"    σ={sigma:.2f}: κ={k_val:.4e}")
            flush()

        # σ=0.5에서의 κ가 최대인지 확인
        idx_05 = sigma_vals.index(0.5)
        k_05 = kappas_by_sigma[idx_05]

        # κ(0.5) vs κ(0.45) 비교 (인접 σ)
        idx_045 = sigma_vals.index(0.45)
        k_045 = kappas_by_sigma[idx_045]
        if k_045 > 0 and np.isfinite(k_045) and np.isfinite(k_05):
            k_ratio_local = k_05 / k_045
        else:
            k_ratio_local = float('inf')

        is_max = all(k_05 >= k for k in kappas_by_sigma if np.isfinite(k))
        sigma_test_count += 1
        if is_max:
            sigma_pass_count += 1
        else:
            sigma_fail_count += 1

        sigma_details.append({
            'idx': z_idx + 1,
            'gamma': z_t,
            'k_05': k_05,
            'k_045': k_045,
            'k_ratio': k_ratio_local,
            'is_max': is_max,
            'kappas': list(zip(sigma_vals, kappas_by_sigma)),
        })

        log(f"    σ=0.5 최대? {'✅ YES' if is_max else '❌ NO'} (κ(0.5)={k_05:.4e}, κ(0.5)/κ(0.45)={k_ratio_local:.2e})")

    log(f"\n  σ-유일성 결과: {sigma_pass_count}/{sigma_test_count}")
    log(f"  PASS: {sigma_pass_count}, FAIL: {sigma_fail_count}")

    # B-10 판별 결과
    if sigma_test_count >= 6:
        pass_rate = sigma_pass_count / sigma_test_count
        if pass_rate >= 0.8:
            b10_verdict = "PASS → 가설 A (N=1→PASS) 확정"
        elif pass_rate <= 0.2:
            b10_verdict = "FAIL → 가설 B (weight=0→PASS) 확정"
        else:
            b10_verdict = f"혼합 (통과율 {pass_rate:.0%}) → 추가 분석 필요"
        log(f"\n  ★ B-10 판별: {b10_verdict}")
    log()
    flush()
else:
    log("  ⚠️ 영점 부족 — σ-유일성 스킵")
    flush()

# ━━━━━━ 최종 요약 ━━━━━━

total_time = time.time() - t_total
log("=" * 70)
log("최종 요약")
log("=" * 70)
log(f"대상: Ramanujan Δ L-함수")
log(f"weight k={WEIGHT}, N={N_COND}, ε={EPSILON}")
log(f"감마: (2π)^{{-s}} Γ(s + 11/2)")
log(f"Fourier 계수: τ(n) (q-product expansion, {len(tau)}개)")
log(f"Hecke 정규화: a(n) = τ(n)/n^{{11/2}}")
log(f"")
log(f"1. 함수방정식: {'✅ PASS' if fe_ok else '❌ FAIL'} (max_rel={fe_max_rel:.2e})")
log(f"2. 영점 발견: {n_zeros}개 (t ∈ [{T_MIN}, {T_MAX}])")

if n_zeros >= 2 and near_fin:
    log(f"3. κ_near median = {near_med:.2f} (n={len(near_fin)}, CV={near_cv:.1f}%)")
    log(f"   κ_far median  = {far_med:.4f}")
    log(f"   κ ratio = {ratio:.1f}×")

    # 비교
    diff_maass = abs(near_med - 1114.9) / 1114.9 * 100
    diff_gl3 = abs(near_med - 1125.16) / 1125.16 * 100
    log(f"   GL(2) Maass 평균 κ_near = 1114.9 → 차이 {diff_maass:.2f}%")
    log(f"   GL(3) sym² κ_near = 1125.16 → 차이 {diff_gl3:.2f}%")

    if diff_maass < 1.0:
        log(f"   → κ_near ≈ 1115: GL(2) degree-의존 상수와 일치 (holomorphic w=12 ≈ Maass w=0)")
    elif diff_gl3 < 1.0:
        log(f"   → κ_near ≈ 1125: GL(3)과 일치 → weight-의존 가능")
    elif near_med > 1125.16:
        log(f"   → κ_near > 1125: weight가 클수록 κ_near 증가? → κ_near(d, k) 2변수 의존")
    else:
        log(f"   → 제3의 값: degree, weight, 감마 인자 구조 등 복합 의존")

    log(f"4. 모노드로미: {mono_pass}/{total_mono} (mono/π≈2.0)")
else:
    log(f"3. κ_near: 영점 부족으로 측정 불가")
    log(f"4. 모노드로미: 영점 부족")

log(f"5. σ-유일성: {sigma_pass_count}/{sigma_test_count}")
if sigma_test_count >= 6:
    pass_rate = sigma_pass_count / sigma_test_count
    if pass_rate >= 0.8:
        log(f"   ★★★ B-10 판별: PASS → 가설 A (N=1 → σ-유일성 PASS) 확정")
    elif pass_rate <= 0.2:
        log(f"   ★★★ B-10 판별: FAIL → 가설 B (weight=0 → σ-유일성 PASS) 확정")
    else:
        log(f"   ⚠️ B-10 판별: 혼합 결과 (통과율 {pass_rate:.0%})")

log(f"")
log(f"총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log(f"")

# TP별 상세 테이블
if n_zeros >= 2 and kappa_near:
    log("TP별 상세:")
    log(f"{'TP':>4} | {'γ':>12} | {'κ_near':>10} | {'mono/π':>8}")
    log("-" * 50)
    for j in range(len(kappa_near)):
        m_str = f"{mono_results[j]:.4f}" if mono_results[j] is not None else "FAIL"
        log(f"  #{j+1:>2} | {zero_approx[j]:>12.6f} | {kappa_near[j]:>10.2f} | {m_str:>8}")

log()

# σ-유일성 상세 테이블
if sigma_details:
    log("σ-유일성 상세:")
    log(f"{'#':>4} | {'γ':>10} | {'κ(0.5)':>12} | {'κ(0.45)':>12} | {'비율':>10} | 판정")
    log("-" * 70)
    for d in sigma_details:
        k05_str = f"{d['k_05']:.4e}" if np.isfinite(d['k_05']) else "inf"
        k045_str = f"{d['k_045']:.4e}" if np.isfinite(d['k_045']) else "inf"
        r_str = f"{d['k_ratio']:.2e}" if np.isfinite(d['k_ratio']) else "inf"
        log(f"  #{d['idx']:>2} | {d['gamma']:>10.4f} | {k05_str:>12} | {k045_str:>12} | {r_str:>10} | {'✅ PASS' if d['is_max'] else '❌ FAIL'}")

log()

# 3점 κ_near 비교
if n_zeros >= 2 and near_fin:
    log("κ_near 3점 비교 (degree × weight):")
    log(f"  | L-함수 | degree | N | weight | κ_near | CV |")
    log(f"  |--------|--------|---|--------|--------|----|")
    log(f"  | GL(1) ζ(s) | 1 | 1 | — | (미측정) | — |")
    log(f"  | GL(2) Maass (평균) | 2 | 1 | 0 | 1114.9 | ~0.15% |")
    log(f"  | GL(2) Δ (w=12) | 2 | 1 | 12 | {near_med:.2f} | {near_cv:.1f}% |")
    log(f"  | GL(3) sym² (6곡선) | 3 | >1 | — | 1125.16 | 0.15% |")
    log()

# σ-유일성 종합 갱신
if sigma_test_count >= 6:
    log("σ-유일성 종합 (갱신):")
    log(f"  | L-함수 | degree | N | weight | σ-유일성 | n(영점) |")
    log(f"  |--------|--------|---|--------|----------|---------|")
    log(f"  | ζ(s) | 1 | 1 | — | PASS | 385 |")
    log(f"  | Maass R=9.53 | 2 | 1 | 0 | PASS | 3 |")
    log(f"  | Maass R=13.78 | 2 | 1 | 0 | PASS | 24 |")
    delta_verdict = "PASS" if sigma_pass_count >= sigma_test_count * 0.8 else "FAIL" if sigma_pass_count <= sigma_test_count * 0.2 else "MIXED"
    log(f"  | Δ (w=12) | 2 | 1 | 12 | **{delta_verdict}** | {sigma_test_count} |")
    log(f"  | 11a1 | 2 | 11 | 2 | FAIL | — |")
    log(f"  | 37a1 | 2 | 37 | 2 | FAIL | — |")
    log(f"  | sym²(11a1) | 3 | 121 | — | FAIL | — |")
    log()

    if delta_verdict == "PASS":
        log("★★★ 결론: 가설 A (N=1 → σ-유일성 PASS) 확정")
        log("  N=1인 모든 테스트 케이스(ζ, Maass×2, Δ)가 PASS")
        log("  N>1인 모든 테스트 케이스(11a1, 37a1, sym²)가 FAIL")
        log("  weight(0 vs 12)는 σ-유일성에 영향 없음")
    elif delta_verdict == "FAIL":
        log("★★★ 결론: 가설 B (weight=0 → σ-유일성 PASS) 확정")
        log("  Maass(w=0, N=1)만 PASS, Δ(w=12, N=1)는 FAIL")
        log("  N=1은 필요조건이나 충분조건은 아님")
        log("  weight=0 (또는 Maass 고유 구조)이 핵심 요인")
    else:
        log("⚠️ 결론: 혼합 결과 — 단순한 가설 A/B로 설명 불가")
        log("  추가 분석 필요 (weight 연속 변화? 감마 인자 구조?)")

flush()
log("\n[완료]")
flush()
