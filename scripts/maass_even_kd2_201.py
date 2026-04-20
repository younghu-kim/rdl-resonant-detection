#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #201] Maass Cusp Form (R=13.7798, Even, SL₂(Z)) — κδ² 4성질 검증
=============================================================================
배경:
  #200: R=9.5337 (odd) → ★★★ (slope=2.0003±0.0003, 역대 최고)
  비평가 지적: "단일 R 검증의 한계"

  이번 실험: 제2 Maass eigenvalue + EVEN 대칭.
  R-독립성 + 짝/홀 대칭 독립성 동시 검증.

대상: 첫 번째 EVEN Maass cusp form on SL₂(Z)
  - 스펙트럼 파라미터 R = 13.77975135189073894...
  - 대칭: EVEN (symmetryclass=0)
  - Conductor N = 1, ε = +1
  - 임계선: Re(s) = 1/2

감마 인자 (EVEN Maass form):
  γ(s) = π^{-s} · Γ((s+iR)/2) · Γ((s-iR)/2)
  Λ(s) = γ(s) · L(s)
  함수방정식: Λ(s) = Λ(1-s)

성공 기준:
  κδ² slope = 2.0 ± 0.1 → R-독립성 + 대칭-독립성 ★★★
  slope 불일치 → B-34 (대칭 의존성) 또는 B-35 (R 의존성)

결과: results/maass_even_kd2_201.txt
=============================================================================
"""
import sys, os, time, math, re
import numpy as np
import mpmath
from datetime import datetime

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'maass_even_kd2_201.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파라미터
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━

SIGMA_CRIT = 0.5  # 임계선: Re(s) = 1/2
DPS = 80          # mpmath 정밀도
R_SPECTRAL_STR = "13.7797513518907389442436732815177125971551325687934870692523882216144503335399700941578316095574275767193896975161287443999283110906237284182748522174244756946463880387122207565011726416547421745062599576096354495336996316021746715513298423864094750660773386134821234354065878318107420092697820008456367397046441444204100974295282406307605104847401124353638446013573382146593753596314736943835033155144106627743098436557742134452522088050388124100292975260253857291002701097628060896803728109964163361464036671036092421008204067682445588106173547212003702918902288023376081045377518378746424327842479500483091850356527573612744036562484423094963126545854160767240415502000035902749088024314316144959870379241402279351119682904039793294356303348804835992575434979100056720560323274201995268866661405302361985953046712894731291318630314322461657254703080705963277904288425653715207499454589776535490253467341697549840245849797609566718798252558999366074290484446985173388003448336962339609979135710858933263619745081881025774487437724121972112794320610"
EPSILON = 1       # root number
N_HERMITE = 80    # Gauss-Hermite 직교점 수
CONTOUR_C = 3.0   # AFE contour shift
DELTAS = [0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05, 0.10]
T_MIN, T_MAX = 3.0, 50.0  # 영점 탐색 범위

mpmath.mp.dps = DPS
R_SPECTRAL = mpmath.mpf(R_SPECTRAL_STR)

log("=" * 72)
log("[실험 #201] Maass Cusp Form (R=13.7798, Even, SL₂(Z)) — κδ² 4성질 ��증")
log("=" * 72)
log(f"시작: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log(f"R = {mpmath.nstr(R_SPECTRAL, 20)}")
log(f"대칭: EVEN (symmetryclass=0)")
log(f"임계선: Re(s) = {SIGMA_CRIT}")
log(f"DPS = {DPS}, N_Hermite = {N_HERMITE}, c = {CONTOUR_C}")
log(f"δ 값: {DELTAS}")
log(f"FE: Λ(s) = Λ(1-s), ε = +1")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. Fourier 계수 로드 (Strömbergsson 데이터)
# ��━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━━

log(f"{T()} [1/7] Fourier 계수 로드...")

COEFF_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "maass_coeff13_raw.txt")
if not os.path.exists(COEFF_FILE):
    log(f"❌ 계수 ���일 없음: {COEFF_FILE}")
    outf.close()
    sys.exit(1)

with open(COEFF_FILE, 'r') as f:
    text = f.read()

# Parse coefficients C[n]
coeff_matches = re.findall(r'C\[(\d+)\]=\s*([-+]?\d+\.\d+)', text)
CN = {}
for n_str, v_str in coeff_matches:
    n = int(n_str)
    CN[n] = mpmath.mpf(v_str[:88])

N_MAX = max(CN.keys())
log(f"  계수 수: {len(CN)} (n=1..{N_MAX})")
log(f"  a(1) = {float(CN[1]):.10f}")
log(f"  a(2) = {float(CN[2]):.15f}")
log(f"  a(3) = {float(CN[3]):.15f}")

# Hecke 곱셈성 검증
a6_check = CN[2] * CN[3]
a6_diff = abs(float(CN[6]) - float(a6_check))
log(f"  Hecke: a(6)={float(CN[6]):.10f}, a(2)*a(3)={float(a6_check):.10f}, diff={a6_diff:.2e}")

a4_check = CN[2]**2 - 1
a4_diff = abs(float(CN[4]) - float(a4_check))
log(f"  Hecke: a(4)={float(CN[4]):.10f}, a(2)²-1={float(a4_check):.10f}, diff={a4_diff:.2e}")

a9_check = CN[3]**2 - 1
a9_diff = abs(float(CN[9]) - float(a9_check))
log(f"  Hecke: a(9)={float(CN[9]):.10f}, a(3)²-1={float(a9_check):.10f}, diff={a9_diff:.2e}")

hecke_ok = a6_diff < 1e-30 and a4_diff < 1e-30 and a9_diff < 1e-30
log(f"  {'✓' if hecke_ok else '✗'} Hecke 관계 검증 {'통과' if hecke_ok else '실패'}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━��━━━━━━���━━━━━━━━━━━���━━━━━━━━━━━━━━
# 2. Λ(s) 구현 (AFE + Gauss-Hermite)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━��━━━

log(f"{T()} [2/7] Λ(s) 구현 (AFE, Gauss-Hermite N={N_HERMITE})...")

_herm_nodes, _herm_weights = np.polynomial.hermite.hermgauss(N_HERMITE)
HERM_NODES = [mpmath.mpf(float(x)) for x in _herm_nodes]
HERM_WEIGHTS = [mpmath.mpf(float(w)) for w in _herm_weights]

iR = mpmath.mpc(0, R_SPECTRAL)
c_afe = mpmath.mpf(CONTOUR_C)


def gamma_factor_maass(s):
    """EVEN Maass form 감마 인자:
    γ(s) = π^{-s} · Γ((s+iR)/2) · Γ((s-iR)/2)
    (ODD: π^{-(s+1)} · Γ((s+1+iR)/2) · Γ((s+1-iR)/2))
    """
    return (mpmath.power(mpmath.pi, -s)
            * mpmath.gamma((s + iR) / 2)
            * mpmath.gamma((s - iR) / 2))


def dirichlet_series(w, N_use=None):
    """D(w) = Σ_{n=1}^{N_use} a(n)/n^w"""
    if N_use is None:
        N_use = N_MAX
    total = mpmath.mpc(0)
    for n in range(1, N_use + 1):
        if n not in CN:
            continue
        total += CN[n] * mpmath.power(n, -w)
    return total


def Lambda_AFE(s, N_use=None):
    """완비 L-함수 Λ(s) via AFE + Gauss-Hermite."""
    if N_use is None:
        N_use = N_MAX
    s_mp = mpmath.mpc(s)
    s1_mp = 1 - s_mp

    prefactor = mpmath.exp(c_afe ** 2) / (2 * mpmath.pi)

    total = mpmath.mpc(0)
    for k in range(N_HERMITE):
        v_k = HERM_NODES[k]
        wt_k = HERM_WEIGHTS[k]
        iv_k = mpmath.mpc(0, v_k)
        w_shift = c_afe + iv_k

        sw = s_mp + w_shift
        gamma_sw = gamma_factor_maass(sw)
        D_s = dirichlet_series(sw, N_use)
        term_s = gamma_sw * D_s / w_shift

        s1w = s1_mp + w_shift
        gamma_s1w = gamma_factor_maass(s1w)
        D_1s = dirichlet_series(s1w, N_use)
        term_1s = gamma_s1w * D_1s / w_shift

        total += wt_k * (term_s + EPSILON * term_1s)

    return prefactor * total


def connection_maass(s, N_use=None):
    """접속 L(s) = Λ'/Λ (수치 미분, h=10^{-20})"""
    h = mpmath.mpf(10) ** (-20)
    val = Lambda_AFE(s, N_use)
    if abs(val) < mpmath.mpf(10) ** (-DPS + 15):
        return mpmath.mpc(1e10, 0)
    dval = (Lambda_AFE(s + h, N_use) - Lambda_AFE(s - h, N_use)) / (2 * h)
    return dval / val


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. 함수방정식 검증
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [3/7] 함수방정식 Λ(s) = Λ(1-s) 검증...")

fe_test_points = [
    mpmath.mpc(0.5, 5),
    mpmath.mpc(0.5, 10),
    mpmath.mpc(0.5, 20),
    mpmath.mpc(0.5, 30),
    mpmath.mpc(0.3, 15),
]
fe_max_err = 0
fe_ok = True

for sp in fe_test_points:
    t_s = time.time()
    L_s = Lambda_AFE(sp)
    L_1s = Lambda_AFE(1 - sp)
    dt_s = time.time() - t_s

    denom = max(abs(L_s), abs(L_1s), mpmath.mpf(1e-100))
    rel_err = float(abs(L_s - EPSILON * L_1s) / denom)
    fe_max_err = max(fe_max_err, rel_err)

    ok = rel_err < 1e-6
    if not ok:
        fe_ok = False
    log(f"  s={mpmath.nstr(sp, 6)}: |Λ|={float(abs(L_s)):.4e}, rel_err={rel_err:.2e} {'✓' if ok else '✗'} ({dt_s:.1f}s)")

sc1_pass = fe_ok
log(f"  SC1: {'✓ PASS' if sc1_pass else '✗ FAIL'} (max_err={fe_max_err:.2e})")
log()


# ━━━��━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━━━���━━━━
# 4. 영점 탐색 (σ=1/2, t∈[3,50])
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [4/7] 영점 탐색 (σ={SIGMA_CRIT}, t∈[{T_MIN},{T_MAX}])...")

sigma = mpmath.mpf(SIGMA_CRIT)
DT_SCAN = 0.3

scan_ts = np.arange(T_MIN, T_MAX + DT_SCAN, DT_SCAN)
scan_vals = []
t_scan_start = time.time()

for i, t_val in enumerate(scan_ts):
    val = Lambda_AFE(mpmath.mpc(0.5, t_val))
    re_val = float(mpmath.re(val))
    scan_vals.append(re_val)

    if (i + 1) % 30 == 0:
        elapsed = time.time() - t_scan_start
        eta = elapsed / (i + 1) * (len(scan_ts) - i - 1)
        log(f"  [{i+1}/{len(scan_ts)}] t={t_val:.1f}: Re(Λ)={re_val:+.4e} ({elapsed:.0f}s, ETA {eta:.0f}s)")

dt_scan_total = time.time() - t_scan_start
log(f"  스캔 완료: {dt_scan_total:.1f}s ({dt_scan_total/60:.1f}분)")

scan_arr = np.array(scan_vals)
intervals = []
for i in range(len(scan_arr) - 1):
    if scan_arr[i] * scan_arr[i + 1] < 0:
        intervals.append((float(scan_ts[i]), float(scan_ts[i + 1])))

log(f"  부호변화 구간: {len(intervals)}개")

# 고정밀 이분법 (40회)
zeros = []
for a, b in intervals:
    va = float(mpmath.re(Lambda_AFE(mpmath.mpc(0.5, a))))
    for _ in range(40):
        mid = (a + b) / 2
        vm = float(mpmath.re(Lambda_AFE(mpmath.mpc(0.5, mid))))
        if not np.isfinite(vm):
            break
        if va * vm < 0:
            b = mid
        else:
            a = mid
            va = vm
    zeros.append((a + b) / 2)

n_zeros = len(zeros)
log(f"  발견된 ��점: {n_zeros}개")
for i, t0 in enumerate(zeros[:15]):
    log(f"    ρ_{i+1} = 0.5 + {t0:.10f}i")
if n_zeros > 15:
    log(f"    ... ({n_zeros - 15}개 추가)")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━━━━━
# 5. SC3a: κδ² log-log slope
# ━━━━━━━━━━━━��━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━���━

log(f"{T()} [5/7] SC3a: κδ² log-log slope 측정...")
log(f"  δ 값: {DELTAS}")

n_kappa_zeros = min(5, n_zeros)
kappa_results = []

for iz in range(n_kappa_zeros):
    t0 = zeros[iz]
    t_start = time.time()
    log(f"  영점 ρ_{iz+1} = 0.5 + {t0:.6f}i:")
    kd2_vals = []

    for delta in DELTAS:
        s_plus = sigma + 1j * mpmath.mpf(str(t0 + delta))
        s_minus = sigma + 1j * mpmath.mpf(str(t0 - delta))

        L_plus = connection_maass(s_plus)
        L_minus = connection_maass(s_minus)
        kappa_plus = float(abs(L_plus) ** 2)
        kappa_minus = float(abs(L_minus) ** 2)
        kappa_avg = (kappa_plus + kappa_minus) / 2.0
        kd2 = kappa_avg * delta ** 2
        kd2_vals.append(kd2)
        log(f"    δ={delta:.3f}: κ_avg={kappa_avg:.4f}, κδ²={kd2:.6f}")

    log_d = np.log(np.array(DELTAS))
    dev = np.array([abs(v - 1.0) for v in kd2_vals])
    valid = dev > 1e-12
    if valid.sum() >= 3:
        log_dev = np.log(dev[valid])
        log_d_v = log_d[valid]
        slope, intercept = np.polyfit(log_d_v, log_dev, 1)
        fitted = slope * log_d_v + intercept
        ss_res = np.sum((log_dev - fitted) ** 2)
        ss_tot = np.sum((log_dev - np.mean(log_dev)) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    else:
        slope, r2 = float('nan'), float('nan')

    dt_zero = time.time() - t_start
    kappa_results.append({'t': t0, 'slope': slope, 'r2': r2, 'kd2_vals': kd2_vals})
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


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━��━
# 6. SC3b: 모노드로미 = ±2π
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━

log(f"{T()} [6/7] SC3b: 모���드로미...")

def monodromy_maass(t_center, radius=0.4, n_steps=96):
    """영점 주위 폐곡선 적분으로 모노드로미 측정."""
    center = mpmath.mpc(0.5, t_center)
    total_angle = mpmath.mpf(0)
    prev_val = Lambda_AFE(center + radius)

    for k in range(1, n_steps + 1):
        theta = 2 * mpmath.pi * k / n_steps
        point = center + radius * mpmath.expj(theta)
        curr_val = Lambda_AFE(point)
        if abs(curr_val) < mpmath.mpf(10) ** (-DPS + 15):
            return None
        ratio = curr_val / prev_val
        darg = mpmath.im(mpmath.log(ratio))
        total_angle += darg
        prev_val = curr_val

    return float(total_angle)

n_mono = min(5, n_zeros)
mono_results = []
for iz in range(n_mono):
    t0 = zeros[iz]
    t_start = time.time()
    mono = monodromy_maass(t0, radius=0.3, n_steps=96)
    dt_m = time.time() - t_start

    if mono is not None:
        mono_pi = mono / math.pi
        mono_results.append(mono_pi)
        log(f"  ρ_{iz+1} (t={t0:.4f}): mono/π = {mono_pi:.4f} (기대: ±2.0) ({dt_m:.1f}s)")
    else:
        log(f"  ρ_{iz+1} (t={t0:.4f}): 측정 실패 ({dt_m:.1f}s)")

sc3b_pass = (len(mono_results) >= 3 and
             all(abs(abs(m) - 2.0) < 0.1 for m in mono_results))
n_pass = sum(1 for m in mono_results if abs(abs(m) - 2.0) < 0.1)
log(f"  SC3b: {'✓ PASS' if sc3b_pass else '✗ FAIL'} ({n_pass}/{len(mono_results)})")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━━━
# 7. SC3c: σ-유일성
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━��━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [7/7] SC3c: σ-유일성...")

def count_sign_changes(sigma_val, t_min=5.0, t_max=40.0, n_points=400):
    """Re(Λ(σ+it)) 부호변화 수."""
    ts = np.linspace(t_min, t_max, n_points)
    count = 0
    prev_sign = None
    for t_val in ts:
        s = mpmath.mpf(str(sigma_val)) + 1j * mpmath.mpf(str(t_val))
        val = float(mpmath.re(Lambda_AFE(s)))
        curr_sign = 1 if val >= 0 else -1
        if prev_sign is not None and curr_sign != prev_sign:
            count += 1
        prev_sign = curr_sign
    return count

sigmas_test = [0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9]
sigma_jumps = {}
for sig in sigmas_test:
    t_s = time.time()
    nj = count_sign_changes(sig)
    dt_s = time.time() - t_s
    sigma_jumps[sig] = nj
    marker = " ← critical" if sig == 0.5 else ""
    log(f"  σ={sig:.1f}: 부호변화 {nj:3d}{marker} ({dt_s:.1f}s)")

crit_jumps = sigma_jumps.get(0.5, 0)
sc3c_pass = all(sigma_jumps[sig] <= crit_jumps for sig in sigmas_test)
log(f"  σ=0.5 최대 여부: {sc3c_pass}")
log(f"  SC3c: {'✓ PASS' if sc3c_pass else '✗ FAIL (GL(2) B-01 패턴 예상)'}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 종합 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━━━━━━━━���━━

total_time = time.time() - START
log("=" * 72)
log("종합 판정")
log("=" * 72)
log(f"  SC1 (FE):       {'✓ PASS' if sc1_pass else '✗ FAIL'} (max_err={fe_max_err:.2e})")
log(f"  SC2 (영점):     {n_zeros}개 발견 (t∈[{T_MIN},{T_MAX}])")
log(f"  SC3a (κδ²):     {'✓ PASS' if sc3a_pass else '✗ FAIL'} (slope={mean_slope:.4f}±{std_slope:.4f})")
log(f"  SC3b (mono):    {'✓ PASS' if sc3b_pass else '✗ FAIL'} ({n_pass}/{len(mono_results)})")
log(f"  SC3c (σ-uniq):  {'✓ PASS' if sc3c_pass else '✗ FAIL'}")
log()

n_pass_total = sum([sc1_pass, n_zeros >= 5, sc3a_pass, sc3b_pass, sc3c_pass])
if n_pass_total >= 4:
    rating = "★★★ 강양성"
elif n_pass_total >= 3:
    rating = "★★☆ 양성"
elif n_pass_total >= 2:
    rating = "★☆☆ 조건부"
else:
    rating = "☆��☆ 음성"

log(f"  PASS: {n_pass_total}/5")
log(f"  판정: {rating}")
log()
log(f"  κδ² slope = {mean_slope:.4f} ± {std_slope:.4f}")
if np.isfinite(mean_slope):
    if abs(mean_slope - 2.0) < 0.1:
        log(f"  → R-독립성 확인! R=13.78(even) ≈ R=9.53(odd)")
        log(f"  → 짝/홀 대칭 독립성 확인")
        log(f"  → Selberg 클래스 degree-2 보편성 강화")
    elif abs(mean_slope - 2.0) < 0.5:
        log(f"  → 근사적 일치 (추가 검증 필요)")
    else:
        log(f"  → R-의존성 또는 대칭-의존성 발견 → B-34/B-35")
log()

# #200 vs #201 비교
log("─" * 72)
log("  #200 (R=9.53, Odd)  vs  #201 (R=13.78, Even) 비교")
log("─" * 72)
log(f"  #200 slope: 2.0003 ± 0.0003 (5영점)")
log(f"  #201 slope: {mean_slope:.4f} ± {std_slope:.4f} ({n_kappa_zeros}영점)")
if np.isfinite(mean_slope):
    diff = abs(mean_slope - 2.0003)
    log(f"  차이: {diff:.4f}")
    log(f"  R-독립성: {'✓ 확인' if diff < 0.05 else '✗ 차이 발견'}")
log()

log(f"총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log("=" * 72)

outf.close()
print(f"\n결과 저장: {RESULT_FILE}")
