#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #203] GL(4) sym³(Δ) — κδ² log-log slope 측정
=============================================================================
배경:
  #102: GL(4) sym³(Δ) ξ-bundle κ 검증 (gammaV=[5.5,6.5,16.5,17.5], N=1)
        lfuninit/lfunzeros 오작동 확인 → lfun(Ldata,s) 직접 사용
  #202: GL(3) sym²(11a1) κδ² slope=2.0000±0.0001 (★★★ 최고 정밀도)

  이번 실험: GL(4) degree-4 L-함수에서 κδ² log-log slope 표준 측정.
  양성이면 10행 비교표 → "GL(1)–GL(4) degree-보편 slope=2.0" 확립.

대상: sym³(Δ), conductor N=1, degree 4, motivic weight 33
  - center = 0.5 (PARI s↔1-s 정규화)
  - 감마: gammaV = [5.5, 6.5, 16.5, 17.5]
  - Euler 계수: a(p) = τ(p)³ - 2p¹¹·τ(p)
  - **lfuninit 사용 금지** (대형 감마이동 오작동)
  - **lfunzeros 사용 금지** (수동 영점 탐색 사용)

방법론: σ-방향 κδ² log-log slope
  Λ(s) = N^(s/2) · Π Γ_R(s+μᵢ) · L(s)
  κ(δ) = |Λ'(center+δ+it₀)/Λ(center+δ+it₀)|²
        = |L'/L(s) + conn_Γ(s)|²
  conn_Γ(s) = (1/2)Σψ((s+μ_j)/2) - 2·log(π)   [d=4, so -d/2·logπ = -2logπ]
  log|κδ²-1| vs log(δ) → slope=2.0 이론값

δ 범위: [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]
영점: 5개 (간격>1.0). #102에서 t₁≈4.1559 확인. lfun 수동 탐색.

+ σ-유일성 (center±0.2, 부호변환 수)
+ 모노드로미 (5영점, 폐곡선 적분, 반지름 0.01, 64분할)

결과: results/gl4_sym3_kd2_203.txt
=============================================================================
"""
import sys, os, time, math, cmath
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'gl4_sym3_kd2_203.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 #203] GL(4) sym³(Δ) — κδ² log-log slope 측정")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"L-함수: sym³(Δ), N=1, degree 4, motivic weight 33")
log(f"center = 0.5 (PARI s↔1-s 정규화)")
log(f"gammaV = [5.5, 6.5, 16.5, 17.5]")
log(f"δ 범위: [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]")
log(f"방법: lfun(Ldata, s) 직접 사용 (NO lfuninit, NO lfunzeros)")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파라미터
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
CENTER = 0.5
GAMMA_V = [5.5, 6.5, 16.5, 17.5]
N_COEFF = 30000
DELTAS = [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]
T_SCAN_MAX = 25.0
H_DERIV = 1e-8  # 수치 미분 간격

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. PARI 초기화 + sym³(Δ) 계수 생성
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [1/7] PARI 초기화 + sym³(Δ) 계수...")

import cypari2
gp = cypari2.Pari()
gp.allocatemem(6000000000)  # 6GB
gp("default(realprecision, 38)")  # #102에서 검증된 값 — 100은 대형 감마이동에서 오작동
log("  PARI: 6GB, realprecision=38 ✅ (#102 호환)")

# τ(n) 생성
log(f"  τ(n) 생성 (n=1..{N_COEFF})...")
t0_s = time.time()

tau = [0] * (N_COEFF + 1)
for n in range(1, N_COEFF + 1):
    tau[n] = int(gp.ramanujantau(n))
    if n % 10000 == 0:
        log(f"    τ: n={n}/{N_COEFF}")

log(f"  τ(n) 완료 ({time.time()-t0_s:.1f}s)")
log(f"  τ(2)={tau[2]}, τ(3)={tau[3]}, τ(5)={tau[5]}")

# 소수 체
sieve = [True] * (N_COEFF + 1)
sieve[0] = sieve[1] = False
for i in range(2, int(N_COEFF**0.5) + 1):
    if sieve[i]:
        for j in range(i*i, N_COEFF+1, i):
            sieve[j] = False
primes = [i for i in range(2, N_COEFF+1) if sieve[i]]
log(f"  소수: {len(primes)}개")

# sym³ Euler 계수
# 정규화: t_p = τ(p) / p^(11/2)
# #102 방식: e1 = tp³ - 2tp, e2 = (tp²-2)(tp²-1)
# Euler 다항식 재귀 (degree 4)
log(f"  sym³ 계수 생성...")
t0_s = time.time()

cpk = {}
for p in primes:
    tp = tau[p] / (p ** 5.5)
    e1 = tp**3 - 2.0*tp
    e2 = (tp**2 - 2.0) * (tp**2 - 1.0)
    cpk[(p, 0)] = 1.0
    cpk[(p, 1)] = e1
    pk = p
    k = 1
    while pk * p <= N_COEFF:
        pk *= p
        k += 1
        cpk[(p, k)] = (e1 * cpk.get((p, k-1), 0.0)
                       - e2 * cpk.get((p, k-2), 0.0)
                       + e1 * cpk.get((p, k-3), 0.0)
                       - cpk.get((p, k-4), 0.0))

cn = np.zeros(N_COEFF + 1)
cn[1] = 1.0
for n in range(2, N_COEFF + 1):
    temp = n
    result = 1.0
    for p in primes:
        if p * p > temp:
            break
        if temp % p == 0:
            k = 0
            while temp % p == 0:
                k += 1
                temp //= p
            result *= cpk.get((p, k), 0.0)
    if temp > 1:
        result *= cpk.get((temp, 1), 0.0)
    cn[n] = result

log(f"  c(2)={cn[2]:.6f}, c(3)={cn[3]:.6f}, c(5)={cn[5]:.6f}")
log(f"  계수 생성 완료 ({time.time()-t0_s:.1f}s)")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. PARI L-함수 생성 (NO lfuninit)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [2/7] PARI Ldata 생성 (NO lfuninit)...")
t0_s = time.time()

cn_str = "[" + ",".join(f"{cn[i]:.15e}" for i in range(1, N_COEFF+1)) + "]"
gp(f"global(Ldata); Ldata = lfuncreate([{cn_str}, 0, [5.5,6.5,16.5,17.5], 1, 1, 0])")

fe_score = float(gp("lfuncheckfeq(Ldata)"))
log(f"  FE check = {fe_score:.1f}  {'✅' if fe_score <= -8 else '⚠️'}")

# L(1/2) — root number 확인
try:
    L_half_re = float(gp("real(lfun(Ldata, 0.5))"))
    L_half_im = float(gp("imag(lfun(Ldata, 0.5))"))
    L_half = complex(L_half_re, L_half_im)
    log(f"  L(1/2) = {L_half.real:.6e} + {L_half.imag:.6e}i")
    has_central_zero = abs(L_half) < 1e-10
    if has_central_zero:
        log(f"  → L(1/2) ≈ 0: root number = -1")
    else:
        log(f"  → L(1/2) ≠ 0: root number = +1")
except Exception as e:
    log(f"  L(1/2) 계산 실패: {e}")
    has_central_zero = False

log(f"  소요: {time.time()-t0_s:.1f}s")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. 수동 영점 탐색 (lfunzeros 대신)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [3/7] 수동 영점 탐색 (|L(0.5+it)| 최솟값)...")
t0_s = time.time()

def eval_absL(t):
    """L(0.5+it)의 절댓값 (#102 방식: complex(gp(...)) 직접)"""
    try:
        v = complex(gp(f"lfun(Ldata, 0.5 + I*{t:.12f})"))
        return abs(v)
    except Exception:
        return 1e30

def refine_zero_golden(t_lo, t_hi, tol=1e-9, max_iter=80):
    """골든 섹션 서치로 |L| 최솟값 위치 정밀 결정"""
    phi = (1 + 5**0.5) / 2
    for _ in range(max_iter):
        if t_hi - t_lo < tol:
            break
        t1 = t_hi - (t_hi - t_lo) / phi
        t2 = t_lo + (t_hi - t_lo) / phi
        if eval_absL(t1) < eval_absL(t2):
            t_hi = t2
        else:
            t_lo = t1
    t_mid = (t_lo + t_hi) / 2
    return t_mid, eval_absL(t_mid)

# 거친 스캔 (dt=0.05)
log("  [1] 거친 스캔 (dt=0.05, t∈[0.2, 25])...")
dt_scan = 0.05
t_scan = np.arange(0.2, T_SCAN_MAX + dt_scan, dt_scan)
absL_scan = []
n_scan = len(t_scan)
for idx, t in enumerate(t_scan):
    val = eval_absL(t)
    absL_scan.append(val)
    if (idx + 1) % 100 == 0:
        min_so_far = min(absL_scan)
        log(f"    스캔 진행: {idx+1}/{n_scan}, min|L|={min_so_far:.4e}")
absL_scan = np.array(absL_scan)
log(f"    스캔 완료: min|L|={absL_scan.min():.4e}, max|L|={absL_scan.max():.4e}")

# 극솟값 후보 (|L| < 1.0 — GL(4)는 변동폭이 클 수 있으므로 여유있게)
candidate_intervals = []
for i in range(1, len(absL_scan) - 1):
    if (absL_scan[i] < absL_scan[i-1] and
        absL_scan[i] < absL_scan[i+1] and
        absL_scan[i] < 1.0):
        candidate_intervals.append((t_scan[max(0, i-1)], t_scan[min(len(t_scan)-1, i+1)]))

log(f"    극솟값 후보: {len(candidate_intervals)}개")

# 정밀 보정
log("  [2] 정밀 보정 (골든 섹션)...")
confirmed_zeros = []
for lo, hi in candidate_intervals:
    t0_val, absL0 = refine_zero_golden(lo, hi)
    if absL0 < 0.1:  # #102 임계값
        confirmed_zeros.append(t0_val)
        log(f"    ✓ t={t0_val:.8f}, |L|={absL0:.4e}")
    elif absL0 < 0.5:
        # 약 영점: GL(4)에서는 수렴이 느릴 수 있음
        confirmed_zeros.append(t0_val)
        log(f"    ~ t={t0_val:.8f}, |L|={absL0:.4e} (약 영점 — 주의)")
    else:
        log(f"    ✗ t={t0_val:.8f}, |L|={absL0:.4e} (영점 아님)")

if has_central_zero:
    confirmed_zeros = [0.0] + confirmed_zeros

all_zeros = np.array(sorted(confirmed_zeros))
N_ALL = len(all_zeros)

log()
log(f"  확인 영점: N={N_ALL}개")
if N_ALL > 0:
    for i, z in enumerate(all_zeros):
        log(f"    γ_{i} = {z:.8f}")
log(f"  소요: {time.time()-t0_s:.1f}s")
log()

if N_ALL < 3:
    log("❌ 영점 부족 (3개 미만) — 중단")
    log("   원인: lfun(Ldata,s) 직접 평가 범위 내 영점이 부족하거나 계수 오류")
    outf.close()
    sys.exit(1)

# 간격>1.0인 5영점 선택
selected_zeros = []
for z in all_zeros:
    if z < 0.1:  # 중앙 영점 제외 (κ 측정에 부적합)
        continue
    if not selected_zeros or (z - selected_zeros[-1]) > 1.0:
        selected_zeros.append(z)
    if len(selected_zeros) >= 5:
        break

log(f"  선택 영점 (간격>1.0): {len(selected_zeros)}개")
for i, z in enumerate(selected_zeros):
    log(f"    ρ_{i+1} = {CENTER} + {z:.8f}i")
log()

if len(selected_zeros) < 3:
    log("❌ 간격 조건 만족 영점 3개 미만 — 중단")
    outf.close()
    sys.exit(1)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. κδ² log-log slope 측정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [4/7] SC3a: κδ² log-log slope 측정...")
log(f"  δ 값: {DELTAS}")
log(f"  κ = |L'/L(s) + conn_Γ(s)|²")
log(f"  conn_Γ(s) = (1/2)Σψ((s+μ_j)/2) - 2·log(π)")
log(f"  L'/L(s) = 수치 미분: (log L(s+h) - log L(s-h))/(2h), h={H_DERIV}")
log()

def eval_L_complex(s_re, s_im):
    """L(s) 복소수 값 반환 (#102 방식)"""
    try:
        v = complex(gp(f"lfun(Ldata, {s_re:.15f} + I*{s_im:.15f})"))
        return v
    except Exception as e:
        log(f"    WARNING: L({s_re:.6f}+{s_im:.6f}i) 계산 실패: {e}")
        return None

def conn_gamma(s_re, s_im):
    """conn_Γ(s) = (1/2)Σψ((s+μ_j)/2) - (d/2)·log(π)  [d=4]"""
    total_re = 0.0
    total_im = 0.0
    for mu in GAMMA_V:
        arg_re = (s_re + mu) / 2.0
        arg_im = s_im / 2.0
        try:
            psi_re = float(gp(f"real(psi({arg_re:.15f} + I*{arg_im:.15f}))"))
            psi_im = float(gp(f"imag(psi({arg_re:.15f} + I*{arg_im:.15f}))"))
            total_re += psi_re
            total_im += psi_im
        except Exception as e:
            log(f"    WARNING: ψ({arg_re:.4f}+{arg_im:.4f}i) 실패: {e}")
            return None
    # (1/2)Σψ - 2·log(π)
    conn_re = 0.5 * total_re - 2.0 * math.log(math.pi)
    conn_im = 0.5 * total_im
    return complex(conn_re, conn_im)

def compute_kappa(s_re, s_im):
    """κ(s) = |L'/L(s) + conn_Γ(s)|²"""
    h = H_DERIV
    L_plus = eval_L_complex(s_re + h, s_im)
    L_minus = eval_L_complex(s_re - h, s_im)
    L_center = eval_L_complex(s_re, s_im)

    if L_plus is None or L_minus is None or L_center is None:
        return None
    if abs(L_center) < 1e-100:
        return None

    try:
        log_Lp = cmath.log(L_plus)
        log_Lm = cmath.log(L_minus)
    except ValueError:
        return None

    LpL = (log_Lp - log_Lm) / (2.0 * h)

    cg = conn_gamma(s_re, s_im)
    if cg is None:
        return None

    total = LpL + cg
    kappa = abs(total) ** 2
    return kappa

kappa_results = []

for iz, t0 in enumerate(selected_zeros):
    t_start = time.time()
    log(f"  영점 ρ_{iz+1} = {CENTER} + {t0:.6f}i:")
    kd2_vals = []

    for delta in DELTAS:
        sigma = CENTER + delta
        kappa = compute_kappa(sigma, t0)

        if kappa is None:
            log(f"    δ={delta:.4f}: 계산 실패")
            continue
        if kappa < 0 or not np.isfinite(kappa):
            log(f"    δ={delta:.4f}: κ 이상값 ({kappa})")
            continue

        kd2 = kappa * delta ** 2
        kd2_vals.append((delta, kd2))
        log(f"    δ={delta:.4f}: κ={kappa:.4f}, κδ²={kd2:.6f}")

    # log|κδ²-1| vs log(δ) 기울기
    slope = float('nan')
    r2 = float('nan')
    if len(kd2_vals) >= 4:
        deltas_arr = np.array([x[0] for x in kd2_vals])
        kd2_arr = np.array([x[1] for x in kd2_vals])
        dev = np.abs(kd2_arr - 1.0)
        valid = dev > 1e-12
        if valid.sum() >= 4:
            log_d = np.log(deltas_arr[valid])
            log_dev = np.log(dev[valid])
            slope, intercept = np.polyfit(log_d, log_dev, 1)
            fitted = slope * log_d + intercept
            ss_res = np.sum((log_dev - fitted) ** 2)
            ss_tot = np.sum((log_dev - np.mean(log_dev)) ** 2)
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    dt_zero = time.time() - t_start
    kappa_results.append({'t': t0, 'slope': slope, 'r2': r2})
    log(f"    → slope = {slope:.4f} (이론: 2.0), R² = {r2:.6f} ({dt_zero:.1f}s)")
    log()

slopes = [r['slope'] for r in kappa_results if np.isfinite(r['slope'])]
r2s = [r['r2'] for r in kappa_results if np.isfinite(r['r2'])]
if slopes:
    mean_slope = np.mean(slopes)
    std_slope = np.std(slopes)
    min_r2 = min(r2s)
    sc3a_pass = abs(mean_slope - 2.0) < 0.05 and min_r2 > 0.999
else:
    mean_slope, std_slope = float('nan'), float('nan')
    min_r2 = float('nan')
    sc3a_pass = False

log(f"  κδ² slope 평균: {mean_slope:.4f} ± {std_slope:.4f}")
log(f"  R² 최솟값: {min_r2:.6f}")
log(f"  SC3a: {'✓ PASS' if sc3a_pass else '✗ FAIL'}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. 모노드로미 (폐곡선 적분)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [5/7] SC3b: 모노드로미 (반지름=0.01, 64분할)...")

def monodromy_gl4(t_center, radius=0.01, n_steps=64):
    """영점 주위 폐곡선 적분으로 winding number 측정"""
    s0_re = CENTER + radius
    s0_im = t_center
    prev_val = eval_L_complex(s0_re, s0_im)
    if prev_val is None or abs(prev_val) < 1e-200:
        return None

    total_angle = 0.0
    for k_step in range(1, n_steps + 1):
        theta = 2 * math.pi * k_step / n_steps
        s_re = CENTER + radius * math.cos(theta)
        s_im = t_center + radius * math.sin(theta)
        cur_val = eval_L_complex(s_re, s_im)

        if cur_val is None or abs(cur_val) < 1e-200:
            return None

        ratio = cur_val / prev_val
        darg = cmath.phase(ratio)
        total_angle += darg
        prev_val = cur_val

    return total_angle

n_mono = min(5, len(selected_zeros))
mono_results = []
for iz in range(n_mono):
    t0 = selected_zeros[iz]
    t_start = time.time()
    mono = monodromy_gl4(t0, radius=0.01, n_steps=64)
    dt_m = time.time() - t_start

    if mono is not None:
        mono_pi = mono / math.pi
        mono_results.append(mono_pi)
        log(f"  ρ_{iz+1} (t={t0:.4f}): mono/π = {mono_pi:.4f} (기대: ±2.0) ({dt_m:.1f}s)")
    else:
        log(f"  ρ_{iz+1} (t={t0:.4f}): 측정 실패 ({dt_m:.1f}s)")

sc3b_pass = (len(mono_results) >= 3 and
             all(abs(abs(m) - 2.0) < 0.1 for m in mono_results))
n_pass_mono = sum(1 for m in mono_results if abs(abs(m) - 2.0) < 0.1)
log(f"  SC3b: {'✓ PASS' if sc3b_pass else '✗ FAIL'} ({n_pass_mono}/{len(mono_results)} = 2π)")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. σ-유일성
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [6/7] SC3c: σ-유일성...")
log(f"  σ = center ± {{0, 0.05, 0.1, 0.15, 0.2}}")

SIGMA_OFFSETS = [-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2]
sigmas_test = [CENTER + off for off in SIGMA_OFFSETS]

# 부호 변환 수 비교
log(f"  부호 변환 수 카운트 (t∈[1, {T_SCAN_MAX}], dt=0.2)...")

def count_sign_changes(sigma_val, t_vals):
    """L(sigma+it)의 실수부 부호 변환 횟수"""
    signs = []
    for t in t_vals:
        v = eval_L_complex(sigma_val, t)
        if v is not None:
            signs.append(1 if v.real >= 0 else -1)
    count = 0
    for i in range(1, len(signs)):
        if signs[i] != signs[i-1]:
            count += 1
    return count

t_grid = np.arange(1.0, T_SCAN_MAX, 0.2)
sigma_counts = {}

for sig in sigmas_test:
    n_change = count_sign_changes(sig, t_grid)
    sigma_counts[sig] = n_change
    marker = " ← center" if abs(sig - CENTER) < 0.01 else ""
    log(f"    σ={sig:.2f}: 부호변환={n_change}{marker}")

center_count = sigma_counts.get(CENTER, 0)
max_sigma = max(sigma_counts, key=sigma_counts.get)
max_count = sigma_counts[max_sigma]

sc3c_pass = (sigma_counts.get(CENTER, 0) == max_count)
if sc3c_pass:
    log(f"  → ✓ PASS (center에서 최대 부호변환: {center_count})")
else:
    log(f"  → ✗ FAIL (최대 σ={max_sigma:.2f}, count={max_count} vs center={center_count})")
log(f"  SC3c: {'✓ PASS' if sc3c_pass else '✗ FAIL'}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7. FE 검증 + 종합 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [7/7] 종합 판정...")

# FE는 lfuncheckfeq 값으로 판정
sc1_pass = fe_score <= -8
log(f"  SC1 (FE): {'✓ PASS' if sc1_pass else '✗ FAIL'} (lfuncheckfeq={fe_score:.1f})")
log()

# 종합
total_time = time.time() - START
log("=" * 72)
log("종합 판정")
log("=" * 72)
log(f"  L-함수: sym³(Δ), degree 4, N=1")
log(f"  center: {CENTER}")
log(f"  gammaV: {GAMMA_V}")
log(f"  영점: {N_ALL}개 발견, {len(selected_zeros)}개 선택")
log()
log(f"  SC1 (FE):       {'✓ PASS' if sc1_pass else '✗ FAIL'} (feq={fe_score:.1f})")
log(f"  SC2 (영점):     {len(selected_zeros)}개 선택 / {N_ALL}개 발견")
log(f"  SC3a (κδ²):     {'✓ PASS' if sc3a_pass else '✗ FAIL'} (slope={mean_slope:.4f}±{std_slope:.4f}, R²≥{min_r2:.6f})")
log(f"  SC3b (mono):    {'✓ PASS' if sc3b_pass else '✗ FAIL'} ({n_pass_mono}/{len(mono_results)})")
log(f"  SC3c (σ-uniq):  {'✓ PASS' if sc3c_pass else '✗ FAIL'} (center={center_count}, max={max_count})")
log()

n_pass_total = sum([sc1_pass, len(selected_zeros) >= 5, sc3a_pass, sc3b_pass, sc3c_pass])
if n_pass_total >= 4:
    rating = "★★★ 강양성"
elif n_pass_total >= 3:
    rating = "★★☆ 양성"
elif n_pass_total >= 2:
    rating = "★☆☆ 조건부"
else:
    rating = "☆☆☆ 음성"

log(f"  PASS: {n_pass_total}/5")
log(f"  판정: {rating}")
log()

# 영점별 상세
log("─" * 72)
log("  영점별 slope 상세")
log("─" * 72)
for r in kappa_results:
    log(f"    t={r['t']:.6f}: slope={r['slope']:.4f}, R²={r['r2']:.6f}")
log()

# 비교표
log("─" * 72)
log("  10행 비교표 (현재 9행 + GL(4) 추가)")
log("─" * 72)
log(f"  {'L-함수':<25} {'degree':>6} {'slope':>10} {'±σ':>8} {'mono':>6}")
log(f"  {'─'*25} {'─'*6} {'─'*10} {'─'*8} {'─'*6}")
log(f"  {'ζ(s)':<25} {'1':>6} {'2.0*':>10} {'—':>8} {'2π':>6}")
log(f"  {'Artin S₃':<25} {'2':>6} {'-0.993**':>10} {'0.001':>8} {'2π':>6}")
log(f"  {'EC 11a1':<25} {'2':>6} {'2.0*':>10} {'—':>8} {'2π':>6}")
log(f"  {'Maass R=9.53 (Odd)':<25} {'2':>6} {'2.0003':>10} {'0.0003':>8} {'2π':>6}")
log(f"  {'Maass R=13.78 (Even)':<25} {'2':>6} {'1.9999':>10} {'0.0008':>8} {'2π':>6}")
log(f"  {'Δ (w=12)':<25} {'2':>6} {'2.0008':>10} {'0.0006':>8} {'2π':>6}")
log(f"  {'Δ·E₄ (w=16)':<25} {'2':>6} {'1.9989':>10} {'0.0035':>8} {'2π':>6}")
log(f"  {'Δ·E₈ (w=20)':<25} {'2':>6} {'1.9984':>10} {'0.0055':>8} {'2π':>6}")
log(f"  {'sym²(11a1)':<25} {'3':>6} {'2.0000':>10} {'0.0001':>8} {'2π':>6}")
log(f"  {'sym³(Δ) [이번]':<25} {'4':>6} {f'{mean_slope:.4f}':>10} {f'{std_slope:.4f}':>8} {'2π' if sc3b_pass else 'FAIL':>6}")
log()

if np.isfinite(mean_slope):
    if abs(mean_slope - 2.0) < 0.05:
        log(f"  ★ GL(4)에서 slope=2.0 확인! → 10행 비교표 확립")
        log(f"    GL(1)×1 + GL(2)×7 + GL(3)×1 + GL(4)×1 = 10행")
        log(f"    degree {{1,2,3,4}} 모두 slope=2.0 → degree-보편성 강력 증거")
    elif abs(mean_slope - 2.0) < 0.5:
        log(f"  ~ 근사적 일치 (추가 검증 필요)")
    else:
        log(f"  ✗ GL(4) 이탈 발견")
log()

log(f"총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log("=" * 72)

outf.close()
print(f"\n결과 저장: {RESULT_FILE}")
