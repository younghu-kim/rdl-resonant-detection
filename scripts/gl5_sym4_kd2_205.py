#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #205] GL(5) sym⁴(11a1) — degree-5 κδ² log-log slope 검증
=============================================================================
배경:
  #85: GL(5) sym⁴(11a1) 6/6 PASS (구 4성질 방식). κδ²≈1 단일 δ에서 확인.
  #200-#204: Paper 2 표준 — κδ² log-log slope 방식 (multi-δ).
  이번 실험: GL(5)를 Paper 2 표준으로 재측정 → 12행 degree 1-5 비교표.

대상: sym⁴(11a1), conductor N=14641, degree 5
  - 11a1: y² - y = x³ - x² (Cremona label)
  - PARI lfunsympow(E, 4) 직접 사용 (#85에서 검증됨)
  - gammaV = [-2, -1, 0, 0, 1], k=5
  - center = 2.5 (arithmetic normalisation, critical line Re(s)=k/2=2.5)

방법론: SC3a κδ² + SC3b 모노드로미 + SC3c σ-유일성
  동일 프로토콜: δ=[0.001..0.03], log|κδ²-1| vs log(δ), slope=2.0 이론값

목적: degree 1-5 보편성 완성. 12행 비교표. Paper 2 완결 근거.

결과: results/gl5_sym4_kd2_205.txt
=============================================================================
"""
import sys, os, time, math, cmath
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'gl5_sym4_kd2_205.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 #205] GL(5) sym⁴(11a1) — degree-5 κδ² log-log slope 검증")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"L-함수: sym⁴(11a1), N=14641, degree 5")
log(f"center = 2.5 (arithmetic normalisation, k=5)")
log(f"gammaV = [-2, -1, 0, 0, 1] (PARI lfunsympow)")
log(f"δ 범위: [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]")
log(f"방법: PARI lfunsympow(E,4) + lfun 직접 평가")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파라미터
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
CENTER = 2.5  # arithmetic normalisation: critical line at Re(s) = k/2 = 2.5
DELTAS = [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]
T_SCAN_MAX = 20.0
H_DERIV = 1e-8  # 수치 미분 간격

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. PARI 초기화 + sym⁴(11a1) L-함수 생성
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [1/7] PARI 초기화 + sym⁴(11a1)...")

import cypari2
gp = cypari2.Pari()
gp.allocatemem(4000000000)  # 4GB
gp("default(realprecision, 38)")
log("  PARI: 4GB, realprecision=38")

# 11a1 타원곡선 생성
gp('E = ellinit([0, -1, 1, 0, 0])')  # 11a1 Cremona model
cond = int(gp('ellglobalred(E)[1]'))
log(f"  E = 11a1: y² - y = x³ - x², conductor={cond}")
assert cond == 11, f"conductor mismatch: {cond}"

# sym⁴ L-함수
gp('global(L4); L4 = lfunsympow(E, 4)')
log(f"  L4 = lfunsympow(E, 4)")

# 기본 검증
gamma_v = str(gp('L4[3]'))
n_sym4 = int(gp('L4[5]'))
log(f"  gammaV = {gamma_v}")
log(f"  conductor = {n_sym4}")

# FE 검증
fe_score = float(gp("lfuncheckfeq(L4)"))
log(f"  lfuncheckfeq = {fe_score:.1f}  {'✅' if fe_score <= -8 else '⚠️'}")

# L(center) 확인
try:
    L_center_re = float(gp("real(lfun(L4, 2.5))"))
    L_center_im = float(gp("imag(lfun(L4, 2.5))"))
    L_center_val = complex(L_center_re, L_center_im)
    log(f"  L(center=2.5) = {L_center_val.real:.6e} + {L_center_val.imag:.6e}i")
    has_central_zero = abs(L_center_val) < 1e-10
    if has_central_zero:
        log(f"  → L(2.5) ≈ 0: 중앙 영점 존재")
    else:
        log(f"  → L(2.5) ≠ 0")
except Exception as e:
    log(f"  L(center) 계산 실패: {e}")
    has_central_zero = False

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. 영점 탐색
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [2/7] 영점 탐색 (lfunzeros, t∈[0, {T_SCAN_MAX}])...")
t0_s = time.time()

zeros_raw = gp(f"lfunzeros(L4, {T_SCAN_MAX})")
all_zeros = [float(z) for z in zeros_raw]
N_ALL = len(all_zeros)

log(f"  영점: {N_ALL}개 발견 ({time.time()-t0_s:.1f}s)")
for i, z in enumerate(all_zeros[:10]):
    log(f"    γ_{i} = {z:.10f}")
if N_ALL > 10:
    log(f"    ... (나머지 {N_ALL-10}개 생략)")
log()

if N_ALL < 3:
    log("❌ 영점 부족 (3개 미만) — 중단")
    outf.close()
    sys.exit(1)

# 간격>1.0인 5영점 선택
selected_zeros = []
for z in all_zeros:
    if z < 0.1:
        continue
    if not selected_zeros or (z - selected_zeros[-1]) > 1.0:
        selected_zeros.append(z)
    if len(selected_zeros) >= 5:
        break

log(f"  선택 영점 (간격>1.0): {len(selected_zeros)}개")
for i, z in enumerate(selected_zeros):
    log(f"    ρ_{i+1} = {CENTER} + {z:.10f}i")
log()

if len(selected_zeros) < 3:
    log("❌ 간격 조건 만족 영점 3개 미만 — 중단")
    outf.close()
    sys.exit(1)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. 공통 함수 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def eval_L_complex(s_re, s_im):
    """L(s) 복소수 값 반환"""
    try:
        v = complex(gp(f"lfun(L4, {s_re:.15f} + I*{s_im:.15f})"))
        return v
    except Exception:
        return None

def conn_gamma(s_re, s_im):
    """Λ'/Λ의 감마 + conductor 기여: (1/2)log(N) + (1/2)Σψ((s+μ_j)/2) - (d/2)·log(π)"""
    GAMMA_V = [-2.0, -1.0, 0.0, 0.0, 1.0]  # degree 5
    N_COND = n_sym4
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
    # (1/2)log(N) + (1/2)Σψ - (d/2)·log(π)  [d=5]
    conn_re = 0.5 * math.log(N_COND) + 0.5 * total_re - 2.5 * math.log(math.pi)
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

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. κδ² log-log slope 측정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [4/7] SC3a: κδ² log-log slope 측정...")
log(f"  δ 값: {DELTAS}")
log(f"  κ = |Λ'/Λ(s)|² = |L'/L(s) + conn(s)|²")
log(f"  conn(s) = (1/2)log(N) + (1/2)Σψ((s+μ_j)/2) - (5/2)·log(π), N={n_sym4}, μ∈[-2,-1,0,0,1]")
log(f"  L'/L(s) = 수치 미분: (log L(s+h) - log L(s-h))/(2h), h={H_DERIV}")
log()

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

def monodromy(t_center, radius=0.01, n_steps=64):
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
    mono = monodromy(t0, radius=0.01, n_steps=64)
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
# 7. 종합 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [7/7] 종합 판정...")

sc1_pass = fe_score <= -8
log(f"  SC1 (FE): {'✓ PASS' if sc1_pass else '✗ FAIL'} (lfuncheckfeq={fe_score:.1f})")
log()

# 종합
total_time = time.time() - START
log("=" * 72)
log("종합 판정")
log("=" * 72)
log(f"  L-함수: sym⁴(11a1), degree 5, N={n_sym4}")
log(f"  center: {CENTER}")
log(f"  gammaV: {gamma_v}")
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

# κδ² 상세 (비교표용)
log("─" * 72)
log("  κδ² t-values (비교표 삽입용)")
log("─" * 72)
for iz, t0 in enumerate(selected_zeros):
    log(f"    ρ_{iz+1}: t = {t0:.4f}")
log()

log(f"  소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log("=" * 72)
outf.close()
