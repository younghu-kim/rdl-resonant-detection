#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #209] Artin S₅ degree-4 — κδ² log-log slope 측정
=============================================================================
배경:
  13행 비교표 degree 1-6 완성. 그러나 degree≥3 항목이 모두 sym^n chain.
  Red Team 지적: "chain bias" — sym^n 이외 구성에서도 slope=2.0 성립하는가?

  이번 실험: x^5-x-1의 Galois group S₅에서 유래하는 irreducible degree-4
  Artin L-함수. ζ_K(s)/ζ(s) = L(s, ρ₄)로 구성.
  sym^n chain과 완전히 독립된 degree-4 L-함수.

대상: Artin ρ₄ (standard rep of S₅), K = Q(α), α^5-α-1=0
  - conductor N = 2869 (= disc(K))
  - degree 4, gammaV = [0, 0, 1, 1]
  - center = 0.5 (FE: Λ(s) = Λ(1-s))
  - epsilon = +1
  - FE check = -393

방법론: σ-방향 κδ² log-log slope (표준 프로토콜)
  κ(δ) = |Λ'(center+δ+it₀)/Λ(center+δ+it₀)|²
  log|κδ²-1| vs log(δ) → slope=2.0 이론값

δ 범위: [0.001, 0.002, 0.003, 0.005, 0.008, 0.01, 0.02, 0.03]
영점 5개 (간격>1.0)
+ 모노드로미 (폐곡선 적분)
+ σ-유일성 (center±0.2)

결과: results/artin_s5_kd2_209.txt
=============================================================================
"""
import sys, os, time, math
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'artin_s5_kd2_209.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 #209] Artin S₅ degree-4 (irreducible) — κδ² log-log slope")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"L-함수: Artin ρ₄ (standard rep of S₅)")
log(f"구성: ζ_K/ζ, K = Q(α), α^5-α-1=0")
log(f"N=2869, degree 4, gammaV=[0,0,1,1], center=0.5, ε=+1")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. PARI 초기화 + L-함수 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━��━━━━━━━━━��━━━━━���━━━

log(f"{T()} [1/7] PARI 초기화...")

import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000 * 1024 * 1024)  # 2GB
gp("default(realprecision, 100)")

# Number field K = Q(α), α^5-α-1=0
gp("K = nfinit(x^5 - x - 1)")
disc = str(gp("K.disc"))
log(f"  disc(K) = {disc}")

# Galois group
gal = str(gp("polgalois(x^5 - x - 1)"))
log(f"  Gal(K̃/Q) = {gal}")

# Dedekind zeta
gp("LK = lfuncreate(K)")
fe_K = float(gp("lfuncheckfeq(LK)"))
log(f"  ζ_K FE = {fe_K:.0f}")

# Riemann zeta
gp("Lzeta = lfuncreate(1)")

# Artin L-function = ζ_K / ζ
gp("Lartin = lfundiv(LK, Lzeta)")

# Verify
CENTER = 0.5
k_art = int(round(float(str(gp("Lartin[4]")))))
log(f"  k = {k_art}, center = {CENTER}")
assert k_art == 1, f"k={k_art} ≠ 1"

fe_art = float(gp("lfuncheckfeq(Lartin)"))
log(f"  Artin L FE = {fe_art:.0f}")

cond = str(gp("Lartin[5]"))
eps = str(gp("Lartin[6]"))
gamma = str(gp("Lartin[3]"))
log(f"  conductor = {cond}, �� = {eps}, γ = {gamma}")

# lfuninit
log(f"  lfuninit([0, 35]) 시작...")
t_init = time.time()
gp("Larti = lfuninit(Lartin, [0, 35])")
log(f"  lfuninit 완료 ({time.time()-t_init:.1f}s)")

# 영점 추출
gp("zvec = lfunzeros(Larti, 30)")
n_total = int(gp("length(zvec)"))
all_zeros = []
for i in range(1, n_total + 1):
    z = float(gp(f"zvec[{i}]"))
    all_zeros.append(z)

log(f"  영점 수 (t∈[0,30]): {n_total}개")
log(f"  처음 10개: {[f'{z:.4f}' for z in all_zeros[:10]]}")

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
# 2. SC1: 함수방정식 검증
# ━━━━━━━━━��━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━

log(f"{T()} [2/7] SC1: 함수방정식 검증 (Λ(s) = Λ(1-s))...")

fe_test_points = [(0.5, 5.0), (0.5, 10.0), (0.5, 20.0), (0.3, 15.0), (0.7, 8.0)]
fe_max_err = 0
fe_ok = True

for s_re, s_im in fe_test_points:
    s1_re = 1.0 - s_re
    gp(f"scur = {s_re:.12f} + I*{s_im:.10f}")
    gp(f"Ls = lfunlambda(Larti, scur)")
    gp(f"scur2 = {s1_re:.12f} + I*{s_im:.10f}")
    gp(f"L1s = lfunlambda(Larti, scur2)")

    abs_Ls = float(gp("abs(Ls)"))
    gp("diff = abs(Ls - L1s)")
    abs_diff = float(gp("diff"))

    if abs_Ls > 1e-200:
        rel_err = abs_diff / abs_Ls
        fe_max_err = max(fe_max_err, rel_err)
        log(f"  s={s_re}+{s_im}i: |Λ(s)|={abs_Ls:.4e}, rel_err={rel_err:.2e}")
    else:
        log(f"  s={s_re}+{s_im}i: |Λ(s)|={abs_Ls:.4e} (영점 근방)")

sc1_pass = fe_art < -30 and fe_max_err < 1e-20
log(f"  FE = {fe_art:.0f}, max_rel_err = {fe_max_err:.2e}")
log(f"  SC1: {'✓ PASS' if sc1_pass else '✗ FAIL'}")
log()


# ━��━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━━━━━━━━━━━━━�����━━━━━━━━━━━━━━━━━━━━━
# 3. SC3a: κδ² log-log slope 측정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━━━━━━━━━━

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
        gp(f"L0val = lfunlambda(Larti, scur)")
        gp(f"Lpval = lfunlambda(Larti, scur, 1)")
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


# ━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━���━━━
# 4. SC3b: 모노드로미 (폐곡선 ���분)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━���━━━━━━━���━━━━

log(f"{T()} [4/7] SC3b: 모노드로미...")

def monodromy(t_center, radius=0.3, n_steps=96):
    """영점 주위 폐곡선 적분으로 모노드로미 측정."""
    import cmath
    s0_re = CENTER + radius
    s0_im = t_center
    gp(f"scur = {s0_re:.12f} + I*{s0_im:.10f}")
    gp(f"prev_val = lfunlambda(Larti, scur)")
    prev_re = float(gp("real(prev_val)"))
    prev_im = float(gp("imag(prev_val)"))
    prev_val = complex(prev_re, prev_im)

    total_angle = 0.0
    for k_step in range(1, n_steps + 1):
        theta = 2 * math.pi * k_step / n_steps
        s_re = CENTER + radius * math.cos(theta)
        s_im = t_center + radius * math.sin(theta)
        gp(f"scur = {s_re:.12f} + I*{s_im:.10f}")
        gp(f"cur_val = lfunlambda(Larti, scur)")
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


# ━━━━━━━━━━━━��━━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━��━━━━━━━━���━━━━━━━━━━━━━
# 5. SC3c: σ-유일성
# ━━━━━━━━━━━━━━━━━━━���━━━━━━━���━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [5/7] SC3c: σ-유일성...")

SIGMA_DELTA = 0.001
sigmas_test = [CENTER - 0.2, CENTER - 0.1, CENTER, CENTER + 0.1, CENTER + 0.2]
log(f"  σ 값: {sigmas_test}, δ={SIGMA_DELTA}")

sigma_pass_count = 0
n_sigma_zeros = min(5, len(selected_zeros))

for iz in range(n_sigma_zeros):
    t0 = selected_zeros[iz]
    log(f"  ��점 ρ_{iz+1} (t={t0:.4f}):")
    kappa_vals = {}

    for sig in sigmas_test:
        s_re = sig + SIGMA_DELTA
        gp(f"scur = {s_re:.12f} + I*{t0:.10f}")
        gp(f"L0val = lfunlambda(Larti, scur)")
        gp(f"Lpval = lfunlambda(Larti, scur, 1)")
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


# ━━━━━━━━━━━��━━━━━━━━━━━━━━━━���━━━━━━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━━
# 6. 추가 FE 정밀도 (realprecision 영향)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━���━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [6/7] 추가 FE 정밀도 검증...")
fe_rp57 = float(gp("default(realprecision, 57); lfuncheckfeq(Lartin)"))
gp("default(realprecision, 100)")
log(f"  FE (rp=57): {fe_rp57:.0f}")
log(f"  FE (rp=100): {fe_art:.0f}")
log()


# ━━━━━━━━━��━━━━━━━━��━━━━━━━━━━��━━━━��━━━━━━━━━━━━━━━━━━━━��━━━━━━��━━━━━
# 7. 최종 판정
# ━━━━━━━━━━━━━━━━��━━━━━━━━━━���━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [7/7] 최종 판정")
log("=" * 72)
log()

n_pass = sum([sc1_pass, sc3a_pass, sc3b_pass, sc3c_pass])

log("┌──────────────┬─────────────────────────────────────────┬─────────┐")
log("│ 검증항목      │ 결과                                    │ 판정    │")
log("├──────────────┼─────────────────────────────────────────┼─────────┤")

fe_str = f"FE={fe_art:.0f} (rp=100) / {fe_rp57:.0f} (rp=57)"
log(f"│ SC1 FE       │ {fe_str:<39s} │ {'✓ PASS' if sc1_pass else '✗ FAIL':>7s} │")

kd2_str = f"slope={mean_slope:.4f}±{std_slope:.4f} (R²={np.mean([r['r2'] for r in kappa_results if np.isfinite(r['r2'])]):.6f})"
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
log(f"L-함수: Artin ρ₄ (S₅ standard, x^5-x-1)")
log(f"  구성: ζ_K/ζ, irreducible degree 4")
log(f"  conductor N={cond}, ε={eps}, γ={gamma}")
log(f"  center={CENTER}, weight k={k_art}")
log()
log(f"의의: sym^n chain 이외 최초 degree-4 κδ² 측정.")
log(f"  Red Team 'chain bias' 지적에 대한 직접 대응.")
log(f"  Artin L-함수 = 수체론에서 유래, 타원곡선 무관.")
log()

total_time = time.time() - START
log(f"총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log(f"종료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

outf.close()
print(f"\n결과 저장: {RESULT_FILE}")
