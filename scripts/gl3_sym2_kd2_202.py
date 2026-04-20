#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #202] GL(3) sym²(11a1) — κδ² log-log slope 측정
=============================================================================
배경:
  #88: sym²(11a1) ξ-bundle κ 측정 (3영점×7δ). κδ²≈1.000 at 작은 δ 확인.
        δ≥0.05에서 이탈 (A(t₀)≈6~12로 큼).
  #200-#201: GL(2) 8행 비교표 완성. κδ² slope=2.0 보편.

  이번 실험: GL(3) degree-3 L-함수에서 κδ² log-log slope 표준 측정.
  양성이면 9행 비교표 → "Selberg class degree-보편 slope=2.0" 확립.

대상: sym²(E, 11a1), conductor N=121, degree 3
  - center = 1.5 (FE: Λ(s) = Λ(3-s))
  - 감마: Γ_R(s+1)·Γ_R(s+1)·Γ_R(s+2) = Γ_R(s+1)²·Γ_R(s+2)

방법론: σ-방향 κδ² log-log slope (GL(2) 표준과 동일)
  κ(δ) = |Λ'(center+δ+it₀)/Λ(center+δ+it₀)|²
  log|κδ²-1| vs log(δ) → slope=2.0 이론값

δ 범위: [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]
  (#88에서 δ≥0.05 이탈 확인 → 0.03 상한)

영점 5개 (간격>1.0, #63에서 확인):
  γ≈3.899, 6.189, 8.654, 10.132, 12.114

+ 모노드로미 (5영점, 폐곡선 적분)
+ σ-유일성 (center±0.2, 5σ값)

결과: results/gl3_sym2_kd2_202.txt
=============================================================================
"""
import sys, os, time, math
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'gl3_sym2_kd2_202.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 #202] GL(3) sym²(11a1) — κδ² log-log slope 측정")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"L-함수: sym²(E, 11a1), N=121, degree 3")
log(f"center = 1.5 (FE: Λ(s) = Λ(3-s))")
log(f"δ 범위: [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]")
log(f"영점: 5개 (간격>1.0)")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. PARI 초기화 + L-함수 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [1/6] PARI 초기화...")

import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000 * 1024 * 1024)  # 2GB
gp("default(realprecision, 100)")

# 타원곡선 11a1 + sym² L-함수
gp("E11 = ellinit([0,-1,1,-10,-20])")
gp("L11 = lfunsympow(E11, 2)")

# k, center 자동 추출
k_raw = str(gp("L11[4]"))
k = int(round(float(k_raw)))
CENTER = k / 2.0
log(f"  k = {k}, center = {CENTER}")
assert abs(CENTER - 1.5) < 1e-9, f"center={CENTER} ≠ 1.5"
log(f"  ✓ center=1.5 확인")

# lfuninit
log(f"  lfuninit([0, 35]) 시작...")
t_init = time.time()
gp("L11i = lfuninit(L11, [0, 35])")
log(f"  lfuninit 완료 ({time.time()-t_init:.1f}s)")

# 영점 추출 (충분히 많이)
gp("zvec = lfunzeros(L11i, 30)")
n_total = int(gp("length(zvec)"))
all_zeros = []
for i in range(1, n_total + 1):
    z = float(gp(f"zvec[{i}]"))
    all_zeros.append(z)

log(f"  영점 수 (t∈[0,30]): {n_total}개")
log(f"  처음 10개: {[f'{z:.3f}' for z in all_zeros[:10]]}")

# 간격>1.0인 5영점 선택
selected_zeros = []
for z in all_zeros:
    if not selected_zeros or (z - selected_zeros[-1]) > 1.0:
        selected_zeros.append(z)
    if len(selected_zeros) >= 5:
        break

log(f"  선택 5영점 (간격>1.0): {[f'{z:.6f}' for z in selected_zeros]}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. κδ² log-log slope 측정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DELTAS = [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]

log(f"{T()} [2/6] SC3a: κδ² log-log slope 측정...")
log(f"  δ 값: {DELTAS}")
log(f"  σ-방향: s = center + δ + i·t₀ = {CENTER} + δ + i·t₀")
log()

kappa_results = []

for iz, t0 in enumerate(selected_zeros):
    t_start = time.time()
    log(f"  영점 ρ_{iz+1} = {CENTER} + {t0:.6f}i:")
    kd2_vals = []

    for delta in DELTAS:
        sigma = CENTER + delta
        # κ = |Λ'(s)/Λ(s)|²
        gp(f"scur = {sigma:.12f} + I*{t0:.10f}")
        gp(f"L0val = lfunlambda(L11i, scur)")
        gp(f"Lpval = lfunlambda(L11i, scur, 1)")
        abs_L0 = float(gp("abs(L0val)"))
        abs_Lp = float(gp("abs(Lpval)"))

        if abs_L0 < 1e-200:
            log(f"    δ={delta:.3f}: |Λ(s)| 너무 작음 ({abs_L0:.2e})")
            continue

        kappa = (abs_Lp / abs_L0) ** 2
        kd2 = kappa * delta ** 2
        kd2_vals.append((delta, kd2))
        log(f"    δ={delta:.3f}: κ={kappa:.4f}, κδ²={kd2:.6f}")

    # log|κδ²-1| vs log(δ) 기울기
    if len(kd2_vals) >= 3:
        deltas_arr = np.array([x[0] for x in kd2_vals])
        kd2_arr = np.array([x[1] for x in kd2_vals])
        dev = np.abs(kd2_arr - 1.0)
        valid = dev > 1e-12
        if valid.sum() >= 3:
            log_d = np.log(deltas_arr[valid])
            log_dev = np.log(dev[valid])
            slope, intercept = np.polyfit(log_d, log_dev, 1)
            fitted = slope * log_d + intercept
            ss_res = np.sum((log_dev - fitted) ** 2)
            ss_tot = np.sum((log_dev - np.mean(log_dev)) ** 2)
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
        else:
            slope, r2 = float('nan'), float('nan')
    else:
        slope, r2 = float('nan'), float('nan')

    dt_zero = time.time() - t_start
    kappa_results.append({'t': t0, 'slope': slope, 'r2': r2})
    log(f"    → slope = {slope:.4f} (이론: 2.0), R² = {r2:.6f} ({dt_zero:.1f}s)")
    log()

slopes = [r['slope'] for r in kappa_results if np.isfinite(r['slope'])]
if slopes:
    mean_slope = np.mean(slopes)
    std_slope = np.std(slopes)
    sc3a_pass = abs(mean_slope - 2.0) < 0.5 and min(slopes) > 1.0
else:
    mean_slope, std_slope = float('nan'), float('nan')
    sc3a_pass = False

log(f"  κδ² slope 평균: {mean_slope:.4f} ± {std_slope:.4f}")
log(f"  SC3a: {'✓ PASS' if sc3a_pass else '✗ FAIL'}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. 모노드로미 (폐곡선 적분)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [3/6] SC3b: 모노드로미...")

def monodromy_gl3(t_center, radius=0.3, n_steps=96):
    """영점 주위 폐곡선 적분으로 모노드로미 측정 (PARI lfunlambda)."""
    import cmath
    # 시작점
    s0_re = CENTER + radius
    s0_im = t_center
    gp(f"scur = {s0_re:.12f} + I*{s0_im:.10f}")
    gp(f"prev_val = lfunlambda(L11i, scur)")
    prev_re = float(gp("real(prev_val)"))
    prev_im = float(gp("imag(prev_val)"))
    prev_val = complex(prev_re, prev_im)

    total_angle = 0.0
    for k_step in range(1, n_steps + 1):
        theta = 2 * math.pi * k_step / n_steps
        s_re = CENTER + radius * math.cos(theta)
        s_im = t_center + radius * math.sin(theta)
        gp(f"scur = {s_re:.12f} + I*{s_im:.10f}")
        gp(f"cur_val = lfunlambda(L11i, scur)")
        cur_re = float(gp("real(cur_val)"))
        cur_im = float(gp("imag(cur_val)"))
        cur_val = complex(cur_re, cur_im)

        if abs(cur_val) < 1e-200:
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
    mono = monodromy_gl3(t0, radius=0.3, n_steps=96)
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
log(f"  SC3b: {'✓ PASS' if sc3b_pass else '✗ FAIL'} ({n_pass_mono}/{len(mono_results)})")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. σ-유일성
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [4/6] SC3c: σ-유일성...")

SIGMA_DELTA = 0.001  # 영점 근방 오프셋
sigmas_test = [CENTER - 0.2, CENTER - 0.1, CENTER, CENTER + 0.1, CENTER + 0.2]
log(f"  σ 값: {sigmas_test}, δ={SIGMA_DELTA}")

sigma_pass_count = 0
n_sigma_zeros = min(5, len(selected_zeros))

for iz in range(n_sigma_zeros):
    t0 = selected_zeros[iz]
    log(f"  영점 ρ_{iz+1} (t={t0:.4f}):")
    kappa_vals = {}

    for sig in sigmas_test:
        s_re = sig + SIGMA_DELTA
        gp(f"scur = {s_re:.12f} + I*{t0:.10f}")
        gp(f"L0val = lfunlambda(L11i, scur)")
        gp(f"Lpval = lfunlambda(L11i, scur, 1)")
        abs_L0 = float(gp("abs(L0val)"))
        abs_Lp = float(gp("abs(Lpval)"))
        if abs_L0 < 1e-200:
            kappa_vals[sig] = float('inf')
            log(f"    σ={sig:.1f}: κ=∞ (영점 근방)")
        else:
            kappa = (abs_Lp / abs_L0) ** 2
            kappa_vals[sig] = kappa
            marker = " ← center" if abs(sig - CENTER) < 0.01 else ""
            log(f"    σ={sig:.1f}: κ={kappa:.4f}{marker}")

    # center에서 최대인지
    center_kappa = kappa_vals.get(CENTER, float('nan'))
    max_sigma = max(kappa_vals, key=lambda s: kappa_vals[s] if not math.isinf(kappa_vals[s]) else 1e300)
    is_pass = abs(max_sigma - CENTER) < 0.05

    if is_pass:
        log(f"    → ✓ PASS (center에서 최대)")
        sigma_pass_count += 1
    else:
        log(f"    → ✗ FAIL (최대 σ={max_sigma:.1f})")
    log()

sc3c_pass = sigma_pass_count >= 3
log(f"  σ-유일성: {sigma_pass_count}/{n_sigma_zeros} PASS")
log(f"  SC3c: {'✓ PASS' if sc3c_pass else '✗ FAIL'}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. FE 검증
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [5/6] SC1: 함수방정식 검증 (Λ(s) = Λ(3-s))...")

fe_test_points = [(1.5, 5.0), (1.5, 10.0), (1.5, 20.0), (1.3, 15.0), (1.7, 8.0)]
fe_max_err = 0
fe_ok = True

for s_re, s_im in fe_test_points:
    s1_re = 3.0 - s_re  # 3-s의 실수부
    gp(f"scur = {s_re:.12f} + I*{s_im:.10f}")
    gp(f"Ls = lfunlambda(L11i, scur)")
    gp(f"scur2 = {s1_re:.12f} + I*{s_im:.10f}")
    gp(f"L1s = lfunlambda(L11i, scur2)")

    abs_Ls = float(gp("abs(Ls)"))
    abs_L1s = float(gp("abs(L1s)"))
    gp("diff = abs(Ls - L1s)")
    abs_diff = float(gp("diff"))

    denom = max(abs_Ls, abs_L1s, 1e-100)
    rel_err = abs_diff / denom
    fe_max_err = max(fe_max_err, rel_err)

    ok = rel_err < 1e-6
    if not ok:
        fe_ok = False
    log(f"  s=({s_re},{s_im}): |Λ|={abs_Ls:.4e}, rel_err={rel_err:.2e} {'✓' if ok else '✗'}")

sc1_pass = fe_ok
log(f"  SC1: {'✓ PASS' if sc1_pass else '✗ FAIL'} (max_err={fe_max_err:.2e})")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. 종합 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

total_time = time.time() - START
log("=" * 72)
log("종합 판정")
log("=" * 72)
log(f"  SC1 (FE):       {'✓ PASS' if sc1_pass else '✗ FAIL'} (max_err={fe_max_err:.2e})")
log(f"  SC2 (영점):     {len(selected_zeros)}개 선택 / {n_total}개 발견")
log(f"  SC3a (κδ²):     {'✓ PASS' if sc3a_pass else '✗ FAIL'} (slope={mean_slope:.4f}±{std_slope:.4f})")
log(f"  SC3b (mono):    {'✓ PASS' if sc3b_pass else '✗ FAIL'} ({n_pass_mono}/{len(mono_results)})")
log(f"  SC3c (σ-uniq):  {'✓ PASS' if sc3c_pass else '✗ FAIL'} ({sigma_pass_count}/{n_sigma_zeros})")
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
log(f"  κδ² slope = {mean_slope:.4f} ± {std_slope:.4f}")
if np.isfinite(mean_slope):
    if abs(mean_slope - 2.0) < 0.1:
        log(f"  → GL(3)에서도 slope=2.0 확인! Degree-보편성 확립")
        log(f"  → 9행 비교표 (GL(1)×1 + GL(2)×7 + GL(3)×1)")
    elif abs(mean_slope - 2.0) < 0.5:
        log(f"  → 근사적 일치 (추가 검증 필요)")
    else:
        log(f"  → GL(3) 이탈 발견 → B-31 (degree 의존 경계)")
log()

# GL(2) 비교
log("─" * 72)
log("  GL(2) vs GL(3) 비교")
log("─" * 72)
log(f"  GL(2) Maass Odd R=9.53:   slope = 2.0003 ± 0.0003 (5영점)")
log(f"  GL(2) Maass Even R=13.78: slope = 1.9999 ± 0.0008 (5영점)")
log(f"  GL(2) Δ (w=12):           slope = 2.0008 ± 0.0006 (5영점)")
log(f"  GL(3) sym²(11a1):         slope = {mean_slope:.4f} ± {std_slope:.4f} ({len(slopes)}영점)")
log()

log(f"총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log("=" * 72)

outf.close()
print(f"\n결과 저장: {RESULT_FILE}")
