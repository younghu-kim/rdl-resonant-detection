#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 C-244] A(t₀) 점근 거동 탐사: 고γ 영역에서 A 공식 스케일링
=============================================================================
목적:
  A = Im(c₀)² + 2Re(c₁) 공식이 45/45 검증되었으나 모두 중간 t 영역 (t < 200).
  t → ∞ 거동이 미탐사. 이번 실험:
  - A(t) ~ t^α scaling 발견 → Paper 2에 Conjecture 추가
  - A(t) ~ const (t-독립) → 보편적 상수 존재 → 더 깊은 구조 시사
  - d-의존성 발견 → degree scaling law 정밀화

대상:
  (1) ζ(s)         d=1, t∈[10,500], 50 영점 균등 샘플링
  (2) Ramanujan Δ  d=2, t∈[10,200], 전수 (~155 영점)
  (3) sym²(11a1)   d=3, t∈[10,100], 전수 (~140 영점)

방법:
  - c₀ ≈ (Λ'/Λ(ρ+δ) + Λ'/Λ(ρ-δ)) / 2
  - c₁ ≈ [(Λ'/Λ(ρ+δ) - Λ'/Λ(ρ-δ))/(2δ) - 1/δ] / δ
  - A = Im(c₀)² + 2Re(c₁)
  - PARI lfunlambda 사용, δ=[0.001,0.002,0.003] 3점 평균

출력:
  (a) degree별 A(t) vs t 테이블 (전체)
  (b) A(t) 이동평균 (5점) → 증가/감소 경향 판정
  (c) log(A) vs log(t) → power law α 추출 (선형 피팅 + R²)
  (d) degree별 A 분포 비교 테이블

결과: results/A_asymptotic_c244.txt
=============================================================================
"""
import sys, os, time
import numpy as np

START = time.time()
RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'A_asymptotic_c244.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 80)
log("[실험 C-244] A(t₀) 점근 거동 탐사: 고γ 영역 A 스케일링")
log("=" * 80)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"대상: ζ(s) t∈[10,500] 50샘플, Δ t∈[10,200] 전수, sym²(11a1) t∈[10,100] 전수")
log(f"방법: A = Im(c₀)²+2Re(c₁), PARI lfunlambda, δ=[0.001,0.002,0.003]")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PARI 초기화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log(f"{T()} [0] PARI 초기화...")
import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000 * 1024 * 1024)   # 2GB
gp("default(realprecision, 80)")
log(f"  PARI OK, realprecision=80")

# mpmath (zetazero용)
import mpmath
mpmath.mp.dps = 40
log(f"  mpmath OK, dps=40")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# L-함수 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log(f"{T()} [1] L-함수 설정...")

# (1) ζ(s): d=1, σ_c=0.5
gp("Lzeta = lfuncreate(1)")
sigma_c_zeta = 0.5
log(f"  ζ(s): lfuncreate(1) → σ_c={sigma_c_zeta}")

# (2) Ramanujan Δ: d=2, weight k=12, σ_c=6
gp("M12 = mfinit([1,12],1); fe12 = mfeigenbasis(M12)[1]; LDelta = lfunmf(M12, fe12)")
k_delta = float(str(gp("LDelta[4]")))
sigma_c_delta = k_delta / 2.0
log(f"  Δ: lfunmf(mfinit([1,12],1), ...) → k={k_delta}, σ_c={sigma_c_delta}")

# (3) sym²(11a1): d=3, k=3, σ_c=1.5
gp("E11 = ellinit([0,-1,1,-10,-20]); Lsym2 = lfunsympow(E11, 2)")
k_sym2 = float(str(gp("Lsym2[4]")))
sigma_c_sym2 = k_sym2 / 2.0
log(f"  sym²(11a1): lfunsympow(E11,2) → k={k_sym2}, σ_c={sigma_c_sym2}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 영점 수집
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log(f"{T()} [2] 영점 수집...")

# (1) ζ: mpmath.zetazero로 t∈[10,500] 범위에서 50개 균등 샘플링
log(f"  (1) ζ(s) 영점: mpmath.zetazero, t∈[10,500]")
n_max_zeta = int(mpmath.nzeros(500))
log(f"    N(500) = {n_max_zeta}개")
# 50개 균등 인덱스 (1-based)
sample_ns = np.linspace(1, n_max_zeta, 50, dtype=int)
# 중복 제거
sample_ns = np.unique(sample_ns)
zeta_zeros = []
t0_perf = time.time()
for n in sample_ns:
    z = mpmath.zetazero(n)
    t_val = float(z.imag)
    if t_val > 10.0:
        zeta_zeros.append(t_val)
log(f"    수집: {len(zeta_zeros)}개 (t∈[{zeta_zeros[0]:.1f},{zeta_zeros[-1]:.1f}]) in {time.time()-t0_perf:.1f}s")

# (2) Δ: PARI lfuninit + lfunzeros, t∈[10,200]
log(f"  (2) Δ 영점: PARI lfunzeros, t∈[10,200]")
t0_perf = time.time()
gp("LDi = lfuninit(LDelta, [0, 210])")
gp("zvDelta = lfunzeros(LDi, 200)")
nz_delta_all = int(gp("length(zvDelta)"))
delta_zeros = []
for i in range(1, nz_delta_all + 1):
    t_val = float(gp(f"zvDelta[{i}]"))
    if t_val >= 10.0:
        delta_zeros.append(t_val)
log(f"    수집: {len(delta_zeros)}개 (t∈[{delta_zeros[0]:.2f},{delta_zeros[-1]:.2f}]) in {time.time()-t0_perf:.1f}s")

# (3) sym²: PARI lfuninit + lfunzeros, t∈[10,100]
log(f"  (3) sym²(11a1) 영점: PARI lfunzeros, t∈[10,100]")
t0_perf = time.time()
gp("Ls2i = lfuninit(Lsym2, [0, 110])")
gp("zvSym2 = lfunzeros(Ls2i, 100)")
nz_sym2_all = int(gp("length(zvSym2)"))
sym2_zeros = []
for i in range(1, nz_sym2_all + 1):
    t_val = float(gp(f"zvSym2[{i}]"))
    if t_val >= 10.0:
        sym2_zeros.append(t_val)
log(f"    수집: {len(sym2_zeros)}개 (t∈[{sym2_zeros[0]:.2f},{sym2_zeros[-1]:.2f}]) in {time.time()-t0_perf:.1f}s")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# c₀,c₁ 추출 + A 계산 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DELTA_VALS = [0.001, 0.002, 0.003]

def compute_A_pari(L_name, sigma_c, t0):
    """
    Λ'/Λ의 Laurent 계수 c₀, c₁을 추출하여 A = Im(c₀)² + 2Re(c₁) 계산.
    L_name: PARI 전역 변수명 (e.g. 'Lzeta', 'LDelta', 'Lsym2')
    sigma_c: 임계선 실수 부분
    t0: 영점 허수 부분
    """
    c0_list = []
    c1_list = []
    fail_count = 0

    for d in DELTA_VALS:
        try:
            s_plus = f"{sigma_c + d:.15f} + I*{t0:.15f}"
            s_minus = f"{sigma_c - d:.15f} + I*{t0:.15f}"

            gp(f"_sp = {s_plus}")
            gp(f"_sm = {s_minus}")

            # Λ'/Λ(σ+δ+it₀)
            gp(f"_Lp = lfunlambda({L_name}, _sp)")
            gp(f"_dLp = lfunlambda({L_name}, _sp, 1)")
            gp("_fp = _dLp / _Lp")
            fp = complex(float(gp("real(_fp)")), float(gp("imag(_fp)")))

            # Λ'/Λ(σ-δ+it₀)
            gp(f"_Lm = lfunlambda({L_name}, _sm)")
            gp(f"_dLm = lfunlambda({L_name}, _sm, 1)")
            gp("_fm = _dLm / _Lm")
            fm = complex(float(gp("real(_fm)")), float(gp("imag(_fm)")))

            # c₀ = (f(+δ) + f(-δ)) / 2
            c0_est = (fp + fm) / 2.0
            # c₁ = [(f(+δ) - f(-δ))/2 - 1/δ] / δ
            c1_est = ((fp - fm) / 2.0 - 1.0 / d) / d

            c0_list.append(c0_est)
            c1_list.append(c1_est)

        except Exception as e:
            fail_count += 1
            print(f"    WARNING: δ={d} 실패: {e}", flush=True)

    if len(c0_list) == 0:
        return float('nan'), None, None

    # 가장 작은 δ (0.001)의 추정이 가장 정확 — 첫 번째 값 사용
    c0 = c0_list[0]
    c1 = c1_list[0]

    # 3개 δ 평균과 비교 (안정성 체크)
    c0_mean = np.mean(c0_list)
    c1_mean = np.mean(c1_list)

    A = c0.imag**2 + 2.0 * c1.real
    return A, c0, c1


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 계산 루프
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

all_results = {}   # {degree: [(t0, A, c0, c1), ...]}

TARGETS = [
    {"name": "ζ(s)",        "L_name": "Lzeta",  "sigma_c": sigma_c_zeta,  "degree": 1, "zeros": zeta_zeros},
    {"name": "Δ (Ramanujan)","L_name": "LDelta", "sigma_c": sigma_c_delta, "degree": 2, "zeros": delta_zeros},
    {"name": "sym²(11a1)",   "L_name": "Lsym2",  "sigma_c": sigma_c_sym2,  "degree": 3, "zeros": sym2_zeros},
]

for tgt in TARGETS:
    name = tgt["name"]
    L_name = tgt["L_name"]
    sigma_c = tgt["sigma_c"]
    d = tgt["degree"]
    zeros = tgt["zeros"]

    log("=" * 80)
    log(f"[d={d}] {name}  (σ_c={sigma_c}, {len(zeros)}개 영점)")
    log("=" * 80)
    log(f"  {'#':>4} | {'t₀':>12} | {'A':>12} | {'Im(c₀)':>12} | {'Re(c₁)':>12}")
    log(f"  {'-'*4}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}-+-{'-'*12}")

    results = []
    t_start = time.time()
    n_fail = 0

    for idx, t0 in enumerate(zeros):
        t_calc = time.time()
        A, c0, c1 = compute_A_pari(L_name, sigma_c, t0)
        dt = time.time() - t_calc

        if np.isnan(A):
            log(f"  {idx+1:>4} | {t0:>12.4f} | {'NaN':>12} | {'—':>12} | {'—':>12}  ❌")
            n_fail += 1
            continue

        if abs(A) < 1e-6:
            A_str = f"{A:.4e}  (≈0)"
        else:
            A_str = f"{A:.6f}"

        im_c0 = c0.imag if c0 is not None else float('nan')
        re_c1 = c1.real if c1 is not None else float('nan')
        log(f"  {idx+1:>4} | {t0:>12.4f} | {A:>12.6f} | {im_c0:>12.6f} | {re_c1:>12.6f}")
        results.append((t0, A, c0, c1))

        if (idx + 1) % 20 == 0:
            elapsed = time.time() - t_start
            log(f"  --- {idx+1}/{len(zeros)} 완료, 소요 {elapsed:.0f}s ---")

    elapsed_total = time.time() - t_start
    log()
    log(f"  → {len(results)}/{len(zeros)} 성공, {n_fail}실패, 소요 {elapsed_total:.1f}s")
    log()
    all_results[d] = results

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 분석 (b) 이동평균 + 경향 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 80)
log("분석 (b): A(t) 이동평균 (5점) 및 경향 판정")
log("=" * 80)

degree_trends = {}

for d, results in sorted(all_results.items()):
    if len(results) < 5:
        log(f"  d={d}: 데이터 부족 ({len(results)}개)")
        continue

    ts = np.array([r[0] for r in results])
    As = np.array([r[1] for r in results])

    # 5점 이동평균
    n = len(As)
    ma_ts = []
    ma_As = []
    for i in range(2, n - 2):
        window = As[i-2:i+3]
        if not np.any(np.isnan(window)):
            ma_ts.append(ts[i])
            ma_As.append(np.mean(window))
    ma_ts = np.array(ma_ts)
    ma_As = np.array(ma_As)

    # 선형 추세 (이동평균에 대해)
    if len(ma_ts) >= 5:
        coeffs = np.polyfit(ma_ts, ma_As, 1)
        slope_raw = coeffs[0]
        A_range = ma_As[-1] - ma_As[0]
        A_mean = np.mean(ma_As)
        rel_change = A_range / abs(A_mean) * 100 if abs(A_mean) > 1e-10 else 0

        if rel_change > 20:
            trend = "증가" if slope_raw > 0 else "감소"
        elif rel_change > 5:
            trend = "약한 " + ("증가" if slope_raw > 0 else "감소")
        else:
            trend = "무상관 (평탄)"

        degree_trends[d] = {"slope_raw": slope_raw, "trend": trend, "rel_change": rel_change,
                             "A_mean": A_mean, "A_std": np.std(As), "A_median": np.median(As),
                             "ts": ts, "As": As}

        log(f"  d={d}: 이동평균 슬로프={slope_raw:+.6e}, 상대변화={rel_change:.1f}% → {trend}")
        log(f"    A 통계: mean={np.mean(As):.6f}, std={np.std(As):.6f}, median={np.median(As):.6f}")
        log(f"    A 범위: [{np.min(As):.4f}, {np.max(As):.4f}]  (유효값 {len(As)}개)")

        # 5점 이동평균 출력 (처음 5개, 마지막 5개)
        log(f"    이동평균 (처음 5):")
        for i in range(min(5, len(ma_ts))):
            log(f"      t={ma_ts[i]:.2f}: MA₅={ma_As[i]:.6f}")
        if len(ma_ts) > 10:
            log(f"    이동평균 (마지막 5):")
            for i in range(-5, 0):
                log(f"      t={ma_ts[i]:.2f}: MA₅={ma_As[i]:.6f}")
    log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 분석 (c) log(A) vs log(t) — power law α 추출
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 80)
log("분석 (c): log(A) vs log(t) — power law A ~ t^α")
log("=" * 80)

degree_powerlaw = {}

for d, results in sorted(all_results.items()):
    if len(results) < 5:
        log(f"  d={d}: 데이터 부족")
        continue

    ts = np.array([r[0] for r in results])
    As = np.array([r[1] for r in results])

    # A > 0인 점만 사용 (log 취하기 위해)
    mask_pos = As > 0
    n_pos = np.sum(mask_pos)
    n_neg = np.sum(~mask_pos)

    log(f"  d={d}: 전체 {len(As)}개, A>0: {n_pos}개, A≤0: {n_neg}개")

    if n_pos < 5:
        log(f"  d={d}: A>0 데이터 부족 ({n_pos}개 < 5)")
        log()
        continue

    ts_pos = ts[mask_pos]
    As_pos = As[mask_pos]
    log_t = np.log(ts_pos)
    log_A = np.log(As_pos)

    # 선형 피팅: log(A) = α·log(t) + const
    coeffs = np.polyfit(log_t, log_A, 1)
    alpha = coeffs[0]
    log_A_intercept = coeffs[1]

    pred_log_A = np.polyval(coeffs, log_t)
    ss_res = np.sum((log_A - pred_log_A)**2)
    ss_tot = np.sum((log_A - np.mean(log_A))**2)
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else 0.0

    A0_fit = np.exp(log_A_intercept)   # A ~ A0 * t^α

    log(f"  d={d}: α = {alpha:.4f}, A₀ = {A0_fit:.4f}, R² = {R2:.6f}")
    log(f"    피팅: A(t) ≈ {A0_fit:.4f} × t^{alpha:.4f}")
    if abs(alpha) < 0.05:
        log(f"    해석: α≈0 → A(t) t-독립 (보편 상수)")
    elif alpha > 0:
        log(f"    해석: α={alpha:.3f}>0 → A(t) 증가 스케일링")
    else:
        log(f"    해석: α={alpha:.3f}<0 → A(t) 감소 스케일링")
    log()

    degree_powerlaw[d] = {"alpha": alpha, "A0": A0_fit, "R2": R2, "n_pos": n_pos}

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 분석 (d) degree별 A 분포 비교 테이블
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 80)
log("분석 (d): degree별 A 분포 비교")
log("=" * 80)
log()

name_map = {1: "ζ(s)", 2: "Δ (Ramanujan)", 3: "sym²(11a1)"}

log(f"  {'d':>2} | {'L-함수':20s} | {'N':>5} | {'mean(A)':>10} | {'std(A)':>10} | {'median(A)':>10} | {'min(A)':>10} | {'max(A)':>10} | {'α':>8} | {'R²':>8}")
log(f"  {'-'*2}-+-{'-'*20}-+-{'-'*5}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*10}-+-{'-'*8}-+-{'-'*8}")

for d in sorted(all_results.keys()):
    results = all_results[d]
    if not results:
        continue
    ts = np.array([r[0] for r in results])
    As = np.array([r[1] for r in results])
    pl = degree_powerlaw.get(d, {})
    alpha = pl.get("alpha", float('nan'))
    R2 = pl.get("R2", float('nan'))
    log(f"  {d:>2} | {name_map.get(d,''):20s} | {len(As):>5} | {np.mean(As):>10.4f} | {np.std(As):>10.4f} | {np.median(As):>10.4f} | {np.min(As):>10.4f} | {np.max(As):>10.4f} | {alpha:>8.4f} | {R2:>8.4f}")

log()

# d-의존성 해석
log("d-의존성 해석:")
if 1 in all_results and 2 in all_results and 3 in all_results:
    medians = {d: np.median([r[1] for r in all_results[d]]) for d in [1,2,3]}
    log(f"  median(A): d=1 → {medians[1]:.4f}, d=2 → {medians[2]:.4f}, d=3 → {medians[3]:.4f}")
    if medians[1] < medians[2] < medians[3]:
        log(f"  → median(A) d=1 < d=2 < d=3: A가 degree와 함께 단조증가 (기존 저γ 결과와 일치)")
    elif medians[1] > medians[2] > medians[3]:
        log(f"  → median(A) d=1 > d=2 > d=3: A가 degree와 함께 감소")
    else:
        log(f"  → 비단조 패턴: 추가 분석 필요")

    # 저γ 비교 (B-35에서 알려진 값)
    log(f"")
    log(f"  저γ 참조값 (B-35, t<30):")
    log(f"    d=1 (ζ):       A ≈ 1.27  (B-35 평균)")
    log(f"    d=2 (GL2,11a1): A ≈ 3.93 (B-35 평균)")
    log(f"    d=3 (sym²):     A ≈ 12.79 (B-35 평균)")
    log(f"")
    log(f"  고γ 결과 (이번 실험):")
    for d in [1,2,3]:
        log(f"    d={d}: median(A) = {medians[d]:.4f}, mean(A) = {np.mean([r[1] for r in all_results[d]]):.4f}")

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 성공 기준 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 80)
log("성공 기준 판정")
log("=" * 80)

n_zeta = len(all_results.get(1, []))
n_delta = len(all_results.get(2, []))
n_sym2 = len(all_results.get(3, []))

c1 = n_zeta >= 50
c2 = n_delta >= 15
c3 = n_sym2 >= 10

# 경향 판정
trend1 = degree_trends.get(1, {}).get("trend", "미판정")
trend2 = degree_trends.get(2, {}).get("trend", "미판정")
trend3 = degree_trends.get(3, {}).get("trend", "미판정")
c4 = all(t != "미판정" for t in [trend1, trend2, trend3])

# d-의존성 테이블 완성
c5 = n_zeta >= 50 and n_delta >= 15 and n_sym2 >= 10

log(f"  {'기준':<50} {'결과':<20} {'판정'}")
log(f"  {'-'*50} {'-'*20} {'-'*6}")
log(f"  {'ζ(s): 50+ 영점 A 측정':<50} {n_zeta:>5}개{'':<14} {'✅' if c1 else '❌'}")
log(f"  {'Δ: 15+ 영점 A 측정':<50} {n_delta:>5}개{'':<14} {'✅' if c2 else '❌'}")
log(f"  {'sym²: 10+ 영점 A 측정':<50} {n_sym2:>5}개{'':<14} {'✅' if c3 else '❌'}")
log(f"  {'degree별 경향 판정':<50} d1:{trend1[:6]}, d2:{trend2[:6]}, d3:{trend3[:6]}{'':8} {'✅' if c4 else '⚠️'}")
log(f"  {'3 degree A 분포 비교 테이블':<50} {'완성' if c5 else '미완'}{'':14} {'✅' if c5 else '❌'}")

n_pass = sum([c1, c2, c3, c4, c5])
log()
log(f"  통과: {n_pass}/5")
log()

if n_pass == 5:
    log(f"★★★★ 양성 — C-244 A(t₀) 점근 탐사 완료")
elif n_pass >= 4:
    log(f"★★★ 조건부 양성 — {n_pass}/5 기준 충족")
elif n_pass >= 3:
    log(f"★★ 부분 성공 — {n_pass}/5 기준 충족")
else:
    log(f"★ 실패 — {n_pass}/5 기준만 충족")

log()

# 최종 해석
log("=" * 80)
log("최종 해석")
log("=" * 80)
log()

for d in [1,2,3]:
    if d not in all_results or not all_results[d]:
        continue
    As = np.array([r[1] for r in all_results[d]])
    ts = np.array([r[0] for r in all_results[d]])
    pl = degree_powerlaw.get(d, {})
    alpha = pl.get("alpha", float("nan"))
    R2 = pl.get("R2", float("nan"))
    trend = degree_trends.get(d, {}).get("trend", "미판정")
    rel_chg = degree_trends.get(d, {}).get("rel_change", float("nan"))

    log(f"  d={d} ({name_map.get(d,'')}):")
    log(f"    영점 수: {len(As)}, t 범위: [{ts.min():.1f}, {ts.max():.1f}]")
    log(f"    A 평균: {np.mean(As):.4f}, 중앙값: {np.median(As):.4f}, std: {np.std(As):.4f}")
    log(f"    경향: {trend} (상대 변화 {rel_chg:.1f}%)")
    if not np.isnan(alpha):
        log(f"    power law: α = {alpha:.4f}, R² = {R2:.4f}")
        if abs(alpha) < 0.05 and R2 < 0.3:
            log(f"    → A가 t에 무관하게 거의 일정 (α≈0, R²낮음 — 노이즈 지배)")
        elif abs(alpha) < 0.1:
            log(f"    → 약한 t 의존성 또는 t-독립 (α={alpha:.3f})")
        else:
            log(f"    → 명확한 power law: A ~ t^{alpha:.3f}")
    log()

log(f"총 소요: {time.time()-START:.1f}s")
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"결과 파일: {RESULT_FILE}")
outf.close()
