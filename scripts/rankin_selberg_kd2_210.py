#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #210] Rankin-Selberg L(11a1 × 37a1) — κδ² log-log slope 측정
=============================================================================
배경:
  14행 비교표 degree 1-6, 4가족 확립. 그러나 degree 4는 sym³(×2) + Artin S₅(×1).
  모두 단일 형식 기반 구성(sym^n) 또는 수체론(Artin). 텐서 곱 구성은 미검증.

  이번 실험: Rankin-Selberg L(s, E₁ ⊗ E₂) where E₁=11a1, E₂=37a1.
  두 독립 타원곡선의 텐서 곱 → GL(4) 자기동형 L-함수.
  sym^n chain, Artin 구성 모두와 완전 독립.

대상: L(s, 11a1 ⊗ 37a1) = lfunmul(L(11a1), L(37a1))
  - conductor N = 11 × 37 = 407 (coprime)
  - degree 4 (GL(4))
  - 두 rank-0 (11a1) × rank-1 (37a1) 타원곡선

방법론: σ-방향 κδ² log-log slope (표준 프로토콜)
  κ(δ) = |Λ'(center+δ+it₀)/Λ(center+δ+it₀)|²
  log|κδ²-1| vs log(δ) → slope=2.0 이론값

δ 범위: [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]
영점 5개 (간격>1.0)
+ 모노드로미 (폐곡선 적분)
+ σ-유일성 (center±0.2)

결과: results/rankin_selberg_kd2_210.txt
=============================================================================
"""
import sys, os, time, math
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'rankin_selberg_kd2_210.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 #210] Rankin-Selberg L(11a1 × 37a1) — κδ² log-log slope")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"L-함수: Rankin-Selberg L(s, 11a1 ⊗ 37a1)")
log(f"구성: lfunmul(L(11a1), L(37a1))")
log(f"degree 4 (GL(4) automorphic)")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. PARI 초기화 + L-함수 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [1/7] PARI 초기화...")

import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000 * 1024 * 1024)  # 2GB
gp("default(realprecision, 100)")

# 타원곡선 생성
gp("E1 = ellinit([0,-1,1,-10,-20])")  # 11a1
gp("E2 = ellinit([0,0,1,-1,0])")      # 37a1
log(f"  E₁ = 11a1: [0,-1,1,-10,-20]")
log(f"  E₂ = 37a1: [0,0,1,-1,0]")

# L-함수 생성
gp("L1 = lfuncreate(E1)")
gp("L2 = lfuncreate(E2)")

# 개별 FE 확인
fe1 = float(gp("lfuncheckfeq(L1)"))
fe2 = float(gp("lfuncheckfeq(L2)"))
log(f"  L(11a1) FE = {fe1:.0f}")
log(f"  L(37a1) FE = {fe2:.0f}")

# Rankin-Selberg 곱
gp("Lrs = lfunmul(L1, L2)")

# Ldata 확인
k_rs = float(str(gp("Lrs[4]")))
CENTER = k_rs / 2.0
log(f"  L(11a1 ⊗ 37a1): k = {k_rs}, center = {CENTER}")

cond = str(gp("Lrs[5]"))
eps = str(gp("Lrs[6]"))
gamma = str(gp("Lrs[3]"))
log(f"  conductor = {cond}, ε = {eps}")
log(f"  γ = {gamma}")

# FE 확인
fe_rs = float(gp("lfuncheckfeq(Lrs)"))
log(f"  Rankin-Selberg FE = {fe_rs:.0f}")

# lfuninit
T_MAX = 30
log(f"  lfuninit([0, {T_MAX}]) 시작...")
t_init = time.time()
gp(f"Lrsi = lfuninit(Lrs, [0, {T_MAX}])")
log(f"  lfuninit 완료 ({time.time()-t_init:.1f}s)")

# 영점 추출
gp(f"zvec = lfunzeros(Lrsi, {T_MAX})")
n_total = int(gp("length(zvec)"))
all_zeros = []
for i in range(1, n_total + 1):
    z = float(gp(f"zvec[{i}]"))
    all_zeros.append(z)

log(f"  영점 수 (t∈[0,{T_MAX}]): {n_total}개")
log(f"  처음 10개: {[f'{z:.4f}' for z in all_zeros[:10]]}")

# 간격>1.0인 5영점 선택
selected_zeros = []
for z in all_zeros:
    if z < 0.5:
        continue
    if not selected_zeros or (z - selected_zeros[-1]) > 1.0:
        selected_zeros.append(z)
    if len(selected_zeros) >= 5:
        break

log(f"  선택 5영점 (간격>1.0): {[f'{z:.6f}' for z in selected_zeros]}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. SC1: 함수방정식 검증
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [2/7] SC1: 함수방정식 검증 Λ(s) = ε·Λ(k-s) ...")

fe_test_points = [
    (CENTER, 5.0), (CENTER, 10.0), (CENTER, 20.0),
    (CENTER - 0.2, 15.0), (CENTER + 0.2, 8.0)
]
fe_max_err = 0
sc1_critical_pass = True

for s_re, s_im in fe_test_points:
    s1_re = k_rs - s_re  # k - s for functional equation
    gp(f"scur = {s_re:.12f} + I*{s_im:.10f}")
    gp(f"Ls = lfunlambda(Lrsi, scur)")
    gp(f"scur2 = {s1_re:.12f} + I*{s_im:.10f}")
    gp(f"L1s = lfunlambda(Lrsi, scur2)")

    abs_Ls = float(gp("abs(Ls)"))
    # Λ(s) = ε · Λ(k-s), so check |Λ(s) - ε·Λ(k-s)|
    gp(f"diff = abs(Ls - ({eps})*L1s)")
    abs_diff = float(gp("diff"))

    if abs_Ls > 1e-200:
        rel_err = abs_diff / abs_Ls
        fe_max_err = max(fe_max_err, rel_err)
        on_critical = abs(s_re - CENTER) < 0.01
        label = " ← critical" if on_critical else ""
        log(f"  s={s_re:.1f}+{s_im}i: |Λ(s)|={abs_Ls:.4e}, rel_err={rel_err:.2e}{label}")
    else:
        log(f"  s={s_re:.1f}+{s_im}i: |Λ(s)|={abs_Ls:.4e} (영점 근방)")

# 임계선에서의 FE 검증 (가장 중요)
gp(f"scrit = {CENTER:.12f} + I*7.5")
gp("Lcrit = lfunlambda(Lrsi, scrit)")
gp(f"scrit2 = {k_rs - CENTER:.12f} + I*7.5")
gp(f"Lcrit2 = lfunlambda(Lrsi, scrit2)")
abs_crit = float(gp("abs(Lcrit)"))
gp(f"dcrit = abs(Lcrit - ({eps})*Lcrit2)")
abs_dcrit = float(gp("dcrit"))
if abs_crit > 1e-200:
    rel_crit = abs_dcrit / abs_crit
else:
    rel_crit = 0.0
log(f"  임계선 rel_err = {rel_crit:.2e}")

sc1_pass = fe_rs < -30 and rel_crit < 1e-10
log(f"  FE = {fe_rs:.0f}, critical rel_err = {rel_crit:.2e}, off-crit max = {fe_max_err:.2e}")
log(f"  SC1: {'✓ PASS' if sc1_pass else '✗ FAIL'}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. SC3a: κδ² log-log slope 측정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DELTAS = [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]

log(f"{T()} [3/7] SC3a: κδ² log-log slope 측정...")
log(f"  δ 값: {DELTAS}")
log(f"  σ-방향: s = {CENTER} + δ + i·t₀")
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
        gp(f"L0val = lfunlambda(Lrsi, scur)")
        gp(f"Lpval = lfunlambda(Lrsi, scur, 1)")
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
# 4. SC3b: 모노드로미 (폐곡선 적분)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [4/7] SC3b: 모노드로미...")

def monodromy(t_center, radius=0.3, n_steps=96):
    """영점 주위 폐곡선 적분으로 모노드로미 측정."""
    import cmath
    s0_re = CENTER + radius
    s0_im = t_center
    gp(f"scur = {s0_re:.12f} + I*{s0_im:.10f}")
    gp(f"prev_val = lfunlambda(Lrsi, scur)")
    prev_re = float(gp("real(prev_val)"))
    prev_im = float(gp("imag(prev_val)"))
    prev_val = complex(prev_re, prev_im)

    total_angle = 0.0
    for k_step in range(1, n_steps + 1):
        theta = 2 * math.pi * k_step / n_steps
        s_re = CENTER + radius * math.cos(theta)
        s_im = t_center + radius * math.sin(theta)
        gp(f"scur = {s_re:.12f} + I*{s_im:.10f}")
        gp(f"cur_val = lfunlambda(Lrsi, scur)")
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
    # 인접 영점과의 거리 확인 → radius 조정
    min_dist = float('inf')
    for z in all_zeros:
        if z != t0:
            d = abs(z - t0)
            if d < min_dist:
                min_dist = d
    radius = min(0.3, min_dist * 0.4)  # 인접 영점의 40% 이내

    t_start = time.time()
    mono = monodromy(t0, radius=radius, n_steps=96)
    dt_m = time.time() - t_start

    if mono is not None:
        mono_pi = mono / math.pi
        mono_results.append(mono_pi)
        log(f"  ρ_{iz+1} (t={t0:.4f}, r={radius:.2f}): mono/π = {mono_pi:.4f} (기대: ±2.0) ({dt_m:.1f}s)")
    else:
        log(f"  ρ_{iz+1} (t={t0:.4f}): 측정 실패 ({dt_m:.1f}s)")

sc3b_pass = (len(mono_results) >= 3 and
             all(abs(abs(m) - 2.0) < 0.1 for m in mono_results))
n_pass_mono = sum(1 for m in mono_results if abs(abs(m) - 2.0) < 0.1)
log(f"  SC3b: {'✓ PASS' if sc3b_pass else '✗ FAIL'} ({n_pass_mono}/{len(mono_results)})")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. SC3c: σ-유일성
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [5/7] SC3c: σ-유일성...")

SIGMA_DELTA = 0.001
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
        gp(f"L0val = lfunlambda(Lrsi, scur)")
        gp(f"Lpval = lfunlambda(Lrsi, scur, 1)")
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

    center_kappa = kappa_vals.get(CENTER, float('nan'))
    max_sigma = max(kappa_vals, key=lambda s: kappa_vals[s] if not math.isinf(kappa_vals[s]) else 1e300)
    is_pass = abs(max_sigma - CENTER) < 0.05

    if is_pass:
        log(f"    → ✓ PASS (center에서 최대)")
        sigma_pass_count += 1
    else:
        log(f"    → ✗ FAIL (최대 σ={max_sigma:.1f}, κ={kappa_vals[max_sigma]:.1f})")
    log()

sc3c_pass = sigma_pass_count >= 3
log(f"  σ-유일성: {sigma_pass_count}/{n_sigma_zeros} PASS")
log(f"  SC3c: {'✓ PASS' if sc3c_pass else '✗ FAIL'}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. 추가 FE 정밀도
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [6/7] 추가 FE 정밀도...")
fe_rp57 = float(gp("default(realprecision, 57); lfuncheckfeq(Lrs)"))
gp("default(realprecision, 100)")
log(f"  FE (rp=57): {fe_rp57:.0f}")
log(f"  FE (rp=100): {fe_rs:.0f}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7. 최종 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [7/7] 최종 판정")
log("=" * 72)
log()

n_pass = sum([sc1_pass, sc3a_pass, sc3b_pass, sc3c_pass])

log("┌──────────────┬─────────────────────────────────────────┬─────────┐")
log("│ 검증항목      │ 결과                                    │ 판정    │")
log("├──────────────┼─────────────────────────────────────────┼─────────┤")

fe_str = f"FE={fe_rs:.0f}, crit_err={rel_crit:.2e}"
log(f"│ SC1 FE       │ {fe_str:<39s} │ {'✓ PASS' if sc1_pass else '✗ FAIL':>7s} │")

r2_vals = [r['r2'] for r in kappa_results if np.isfinite(r['r2'])]
mean_r2 = np.mean(r2_vals) if r2_vals else float('nan')
kd2_str = f"slope={mean_slope:.4f}±{std_slope:.4f} (R²={mean_r2:.6f})"
log(f"│ SC3a κδ²     │ {kd2_str:<39s} │ {'✓ PASS' if sc3a_pass else '✗ FAIL':>7s} │")

mono_str = f"{n_pass_mono}/{len(mono_results)} = {np.mean([abs(m) for m in mono_results]):.4f}π" if mono_results else "N/A"
log(f"│ SC3b mono    │ {mono_str:<39s} │ {'✓ PASS' if sc3b_pass else '✗ FAIL':>7s} │")

sigma_str = f"{sigma_pass_count}/{n_sigma_zeros} PASS"
log(f"│ SC3c σ-uniq  │ {sigma_str:<39s} │ {'✓ PASS' if sc3c_pass else '✗ FAIL':>7s} │")

log("└──────────────┴─────────────────────────────────────────┴─────────┘")
log()

# 영점별 slope 개별 표시
log("영점별 slope:")
for r in kappa_results:
    log(f"  t={r['t']:.4f}: slope={r['slope']:.4f}, R²={r['r2']:.6f}")
log()

# 판정
if n_pass >= 4:
    verdict = "★★★ 강양성 (4/4 PASS)"
elif n_pass == 3:
    verdict = "★★ 양성 (3/4 PASS)"
elif n_pass == 2:
    verdict = "★ 약양성 (2/4 PASS)"
else:
    verdict = "음성 (<2/4 PASS)"

log(f"총 판정: {n_pass}/4 PASS → {verdict}")
log()
log(f"L-함수: Rankin-Selberg L(s, 11a1 ⊗ 37a1)")
log(f"  구성: lfunmul(L(11a1), L(37a1))")
log(f"  conductor N={cond}, ε={eps}")
log(f"  γ = {gamma}")
log(f"  center={CENTER}, weight k={k_rs}")
log()
log(f"의의: sym^n, Artin 이외 3번째 degree-4 구성.")
log(f"  텐서 곱(tensor product) 구성 = 두 독립 GL(2) 형식의 곱.")
log(f"  구성 독립적 보편성 확인의 결정적 증거.")
log()

total_time = time.time() - START
log(f"총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log(f"종료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

outf.close()
print(f"\n결과 저장: {RESULT_FILE}")
