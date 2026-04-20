#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #206] GL(6) sym⁵(11a1) — degree-6 κδ² log-log slope 검증
=============================================================================
배경:
  #205: GL(5) sym⁴(11a1) 4/5 PASS (★★★). slope=1.9999±0.0004.
  Paper 2 완결 선언 (degree 1-5, 12행 비교표).
  이번 실험: degree 6 돌파 시도. B-31 해결 여부 결정.

대상: sym⁵(11a1), conductor N=161051 (=11⁵), degree 6
  - 11a1: y² - y = x³ - x² (Cremona label)
  - PARI lfunsympow(E, 5) 직접 사용
  - gammaV = [-2, -1, -1, 0, 0, 1] (PARI L[2] 슬롯에서 확인)
  - center = 3.0 (arithmetic normalisation, critical line Re(s)=k/2=3.0, k=6)
  - root number ε = -1 → L(3.0)=0 (중앙 영점)

방법론: SC3a κδ² + SC3b 모노드로미 + SC3c σ-유일성
  동일 프로토콜: δ=[0.001..0.03], log|κδ²-1| vs log(δ), slope=2.0 이론값

사전 검증 (사이클 #217):
  - PARI lfunsympow(11a1,5) 성공 (0.0s)
  - FE check = -55 (55자리)
  - 150개 영점 in [0,50] (134s, 2GB 메모리)
  - first nontrivial zero ≈ 1.964

결과: results/gl6_sym5_kd2_206.txt
=============================================================================
"""
import sys, os, time, math, cmath
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'gl6_sym5_kd2_206.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 #206] GL(6) sym⁵(11a1) — degree-6 κδ² log-log slope 검증")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"L-함수: sym⁵(11a1), N=161051, degree 6")
log(f"center = 3.0 (arithmetic normalisation, k=6)")
log(f"gammaV = [-2, -1, -1, 0, 0, 1] (PARI lfunsympow)")
log(f"δ 범위: [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]")
log(f"방법: PARI lfunsympow(E,5) + lfun 직접 평가")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파라미터
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
CENTER = 3.0
DELTAS = [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]
T_SCAN_MAX = 30.0
H_DERIV = 1e-8
GAMMA_V = [-2.0, -1.0, -1.0, 0.0, 0.0, 1.0]
DEGREE = 6

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. PARI 초기화 + sym⁵(11a1) L-함수 생성
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [1/7] PARI 초기화 + sym⁵(11a1)...")

import cypari2
gp = cypari2.Pari()
gp.allocatemem(4000000000)
gp("default(realprecision, 57)")
log("  PARI: 4GB, realprecision=57")

gp('E = ellinit([0, -1, 1, 0, 0])')
cond = int(gp('ellglobalred(E)[1]'))
log(f"  E = 11a1: conductor={cond}")
assert cond == 11

gp('global(L5); L5 = lfunsympow(E, 5)')
log(f"  L5 = lfunsympow(E, 5)")

n_cond = int(gp('L5[5]'))
root_num = int(gp('L5[6]'))
log(f"  conductor = {n_cond}")
log(f"  root number = {root_num}")

fe_score = float(gp("lfuncheckfeq(L5)"))
log(f"  lfuncheckfeq = {fe_score:.1f}  {'✅' if fe_score <= -8 else '⚠️'}")

try:
    L_center_re = float(gp(f"real(lfun(L5, {CENTER}))"))
    L_center_im = float(gp(f"imag(lfun(L5, {CENTER}))"))
    log(f"  L(center={CENTER}) = {L_center_re:.6e} + {L_center_im:.6e}i")
except Exception as e:
    log(f"  L(center) 계산 실패: {e}")

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. 영점 탐색
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [2/7] 영점 탐색 (t∈[0, {T_SCAN_MAX}])...")
t0_s = time.time()

zeros_raw = gp(f"lfunzeros(L5, {T_SCAN_MAX})")
all_zeros = [float(z) for z in zeros_raw]
N_ALL = len(all_zeros)

log(f"  영점: {N_ALL}개 ({time.time()-t0_s:.1f}s)")
for i, z in enumerate(all_zeros[:10]):
    log(f"    γ_{i} = {z:.10f}")
if N_ALL > 10:
    log(f"    ... (나머지 {N_ALL-10}개)")
log()

if N_ALL < 3:
    log("❌ 영점 부족 — 중단")
    outf.close()
    sys.exit(1)

selected_zeros = []
for z in all_zeros:
    if z < 0.5:
        continue
    if not selected_zeros or (z - selected_zeros[-1]) > 1.0:
        selected_zeros.append(z)
    if len(selected_zeros) >= 5:
        break

log(f"  선택 (간격>1.0): {len(selected_zeros)}개")
for i, z in enumerate(selected_zeros):
    log(f"    ρ_{i+1} = {CENTER} + {z:.10f}i")
log()

if len(selected_zeros) < 3:
    log("❌ 선택 영점 3개 미만 — 중단")
    outf.close()
    sys.exit(1)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. 공통 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def eval_L_complex(s_re, s_im):
    try:
        v = complex(gp(f"lfun(L5, {s_re:.15f} + I*{s_im:.15f})"))
        return v
    except Exception:
        return None

def conn_gamma(s_re, s_im):
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
        except Exception:
            return None
    conn_re = 0.5 * math.log(n_cond) + 0.5 * total_re - (DEGREE / 2.0) * math.log(math.pi)
    conn_im = 0.5 * total_im
    return complex(conn_re, conn_im)

def compute_kappa(s_re, s_im):
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
    return abs(total) ** 2

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. κδ² log-log slope
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [4/7] SC3a: κδ² log-log slope...")
log(f"  δ: {DELTAS}")
log(f"  conn(s) = (1/2)log({n_cond}) + (1/2)Σψ((s+μ)/2) - 3·log(π)")
log()

kappa_results = []

for iz, t0 in enumerate(selected_zeros):
    t_start = time.time()
    log(f"  ρ_{iz+1} = {CENTER} + {t0:.6f}i:")
    kd2_vals = []
    for delta in DELTAS:
        sigma = CENTER + delta
        kappa = compute_kappa(sigma, t0)
        if kappa is None or kappa < 0 or not np.isfinite(kappa):
            log(f"    δ={delta:.4f}: 실패")
            continue
        kd2 = kappa * delta ** 2
        kd2_vals.append((delta, kd2))
        log(f"    δ={delta:.4f}: κ={kappa:.4f}, κδ²={kd2:.6f}")

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
    log(f"    → slope = {slope:.4f}, R² = {r2:.6f} ({dt_zero:.1f}s)")
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

log(f"  slope 평균: {mean_slope:.4f} ± {std_slope:.4f}")
log(f"  R² 최솟값: {min_r2:.6f}")
log(f"  SC3a: {'✓ PASS' if sc3a_pass else '✗ FAIL'}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. 모노드로미
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [5/7] SC3b: 모노드로미 (반지름=0.01, 64분할)...")

def monodromy(t_center, radius=0.01, n_steps=64):
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
    mono = monodromy(t0, radius=0.01, n_steps=64)
    dt_m = time.time() - t_start
    if mono is not None:
        mono_pi = mono / math.pi
        mono_results.append(mono_pi)
        log(f"  ρ_{iz+1} (t={t0:.4f}): mono/π = {mono_pi:.4f} ({dt_m:.1f}s)")
    else:
        log(f"  ρ_{iz+1} (t={t0:.4f}): 실패 ({dt_m:.1f}s)")

sc3b_pass = (len(mono_results) >= 3 and
             all(abs(abs(m) - 2.0) < 0.1 for m in mono_results))
n_pass_mono = sum(1 for m in mono_results if abs(abs(m) - 2.0) < 0.1)
log(f"  SC3b: {'✓ PASS' if sc3b_pass else '✗ FAIL'} ({n_pass_mono}/{len(mono_results)} = 2π)")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. σ-유일성
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [6/7] SC3c: σ-유일성...")

SIGMA_OFFSETS = [-0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2]
sigmas_test = [CENTER + off for off in SIGMA_OFFSETS]

log(f"  부호 변환 카운트 (t∈[1, {T_SCAN_MAX}], dt=0.2)...")

def count_sign_changes(sigma_val, t_vals):
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
    log(f"    σ={sig:.2f}: {n_change}{marker}")

center_count = sigma_counts.get(CENTER, 0)
max_sigma = max(sigma_counts, key=sigma_counts.get)
max_count = sigma_counts[max_sigma]

sc3c_pass = (sigma_counts.get(CENTER, 0) == max_count)
if sc3c_pass:
    log(f"  → ✓ PASS (center 최대: {center_count})")
else:
    log(f"  → ✗ FAIL (max σ={max_sigma:.2f}: {max_count} vs center: {center_count})")
log(f"  SC3c: {'✓ PASS' if sc3c_pass else '✗ FAIL'}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7. 종합
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [7/7] 종합 판정...")

sc1_pass = fe_score <= -8

total_time = time.time() - START
log("=" * 72)
log("종합 판정")
log("=" * 72)
log(f"  L-함수: sym⁵(11a1), degree 6, N={n_cond}")
log(f"  center: {CENTER}, ε={root_num}")
log(f"  gammaV: {GAMMA_V}")
log(f"  영점: {N_ALL}개 / 선택 {len(selected_zeros)}개")
log()
log(f"  SC1 (FE):       {'✓ PASS' if sc1_pass else '✗ FAIL'} ({fe_score:.1f})")
log(f"  SC2 (영점):     {len(selected_zeros)}개 / {N_ALL}개")
log(f"  SC3a (κδ²):     {'✓ PASS' if sc3a_pass else '✗ FAIL'} ({mean_slope:.4f}±{std_slope:.4f}, R²≥{min_r2:.6f})")
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

log("─" * 72)
log("  영점별 slope")
log("─" * 72)
for r in kappa_results:
    log(f"    t={r['t']:.6f}: slope={r['slope']:.4f}, R²={r['r2']:.6f}")
log()

log("─" * 72)
log("  모노드로미")
log("─" * 72)
for iz, m in enumerate(mono_results):
    log(f"    ρ_{iz+1}: mono/π = {m:.6f}")
log()

log(f"  소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log("=" * 72)
outf.close()
