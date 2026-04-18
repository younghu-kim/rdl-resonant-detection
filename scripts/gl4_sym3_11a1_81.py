#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] #81 — GL(4) sym³(11a1) PARI 4성질 검증 + A(t₀) 측정
=============================================================================
목적: A(d=4)=286.26의 급증이 degree 효과인지 weight 효과인지 분리.
  - sym³(Δ)     : d=4, w=11 → A=286.26 (★★★ 강양성, #80)
  - sym³(11a1)  : d=4, w=1  → A=?  (이번 실험)
    - A ≈ 286 → degree 효과 (weight 무관)
    - A ≈ 13  → weight 효과 (d=3과 동일)
    - 중간값   → 혼합 의존

L-함수: lfunsympow(ellinit([0,-1,1,-10,-20]), 3)
  gammaV=[-1,0,0,1], k=4, N=1331 (=11³), ε=+1

4성질: FE, 영점, κ_near (δ=0.001), 모노드로미
핵심 측정: A(t₀) mean, std, CV — sym³(Δ) 직접 비교

결과: results/gl4_sym3_11a1_81.txt
=============================================================================
"""

import sys, os, time
import statistics

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

# ─── 출력 ─────────────────────────────────────────────────────────────────
OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "results", "gl4_sym3_11a1_81.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))
def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ─── 파라미터 ─────────────────────────────────────────────────────────────
DELTA_KAPPA = 0.001     # κ_near 측정 δ
T_RANGE     = [0, 55]   # 영점 탐색 범위 (pari_80.txt와 동일)
T_KAPPA_MIN = 0.5       # κ 측정 최소 t
T_KAPPA_MAX = 40        # κ 측정 최대 t

# 비교 기준 (기존 결과)
D1_A = 1.2727    # ζ (d=1, w=-)
D2_A = 3.9300    # GL(2) Δ avg (d=2)
D3_A = 12.7900   # sym²(11a1) (d=3, w=1)
D4_DELTA_A = 286.26  # sym³(Δ) (d=4, w=11, #80)

# ─── PARI 초기화 ──────────────────────────────────────────────────────────
import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000000000)           # 2GB (수학자 지시)
gp("default(realprecision, 100)")    # 100자리 정밀도

log("=" * 72)
log("결과 #81 — GL(4) sym³(11a1) PARI 4성질 검증 + A(t₀) 측정")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"★ 핵심: A(d=4)의 degree/weight 기여 분리")
log(f"  sym³(Δ)   d=4, w=11: A={D4_DELTA_A:.2f} (기준, #80)")
log(f"  sym³(11a1) d=4, w=1: A=? (이번 실험)")
log(f"  만약 A≈286 → degree 효과 / A≈13 → weight 효과")
log(f"PARI 파라미터: gammaV=[-1,0,0,1], k=4, N=1331, ε=+1")
log(f"lfuninit 범위: [0, {T_RANGE[1]}], δ={DELTA_KAPPA}")
log(f"PARI realprecision=100, 메모리=2GB")
log()
flush_file()

# =====================================================================
# [1] L-함수 구성: lfunsympow(11a1, 3)
# =====================================================================
log("[1] L-함수 구성 — lfunsympow(11a1, 3)")
t0 = time.time()

gp("E = ellinit([0,-1,1,-10,-20])")
gp("L81 = lfunsympow(E, 3)")

# 파라미터 확인
params = str(gp("lfunparams(L81)"))
log(f"  11a1 = ellinit([0,-1,1,-10,-20])")
log(f"  L81 = lfunsympow(E, 3)")
log(f"  params = {params}")
log(f"  기대: gammaV=[-1,0,0,1], k=4, N=1331, ε=+1")

# lfuninit (영점 계산을 위해)
gp(f"Linit = lfuninit(L81, [{T_RANGE[0]}, {T_RANGE[1]}])")
log(f"  lfuninit([0, {T_RANGE[1]}]) 완료 ({time.time()-t0:.1f}s)")
log()
flush_file()

# =====================================================================
# [2] ★ 함수방정식 검증 (FE)
# =====================================================================
log("=" * 72)
log("[2] ★ 함수방정식 검증 (lfuncheckfeq)")
log("=" * 72)
t0 = time.time()

fe = float(gp("lfuncheckfeq(Linit)"))
log(f"  FE = {fe:.1f}자리 정밀도")
fe_pass = fe <= -8
log(f"  판정: {'✅ PASS' if fe_pass else '❌ FAIL'} (≤-8자리)")
log(f"  ({time.time()-t0:.1f}s)")
log()
flush_file()

# =====================================================================
# [3] ★ 영점 탐색 (lfunzeros)
# =====================================================================
log("=" * 72)
log("[4] ★ 영점 탐색 (lfunzeros)")
log("=" * 72)
t0 = time.time()

z = gp(f"lfunzeros(Linit, {T_RANGE[1]})")
nz = len(z)
zeros = [float(z[i]) for i in range(nz)]

# 비자명 영점 (ε=+1이므로 t=0 trivial zero 없음)
zeros_nontrivial = [t for t in zeros if t > T_KAPPA_MIN]

log(f"  총 영점: {nz}개 (t in [0, {T_RANGE[1]}])")
log(f"  비자명 영점: {len(zeros_nontrivial)}개 (t > {T_KAPPA_MIN})")
if zeros_nontrivial:
    log(f"  t₁ = {zeros_nontrivial[0]:.6f}")

for i, t_z in enumerate(zeros[:20]):
    log(f"  영점 #{i+1:2d}: t = {t_z:.6f}")
if nz > 20:
    log(f"  ... (총 {nz}개)")

zeros_pass = len(zeros_nontrivial) >= 20
log(f"  판정: {'✅ PASS' if zeros_pass else '❌ FAIL'} (≥20개)")
log(f"  ({time.time()-t0:.1f}s)")
log()
flush_file()

# =====================================================================
# [4] Hardy Z 영점 검증 (진짜 영점인지 확인)
# =====================================================================
log("[5] Hardy Z 영점 검증")
t0 = time.time()

z_check_success = 0
for i in range(min(5, len(zeros_nontrivial))):
    t_zero = zeros_nontrivial[i]
    try:
        Z_val = float(gp(f"lfunhardy(Linit, {t_zero})"))
        log(f"  Z({t_zero:.6f}) = {Z_val:.4e}")
        z_check_success += 1
    except Exception as e:
        log(f"  Z({t_zero:.6f}) = ❌ 오류: {str(e)[:50]}")

log(f"  ({time.time()-t0:.1f}s)")
log()
flush_file()

# =====================================================================
# [5] ★ 모노드로미 — 단순 영점 확인 (부호 변환)
# =====================================================================
log("=" * 72)
log("[6] ★ 모노드로미 — 단순 영점 확인 (Hardy Z 부호 변환)")
log("=" * 72)
log("  단순 영점 ↔ Z(t₀-ε)·Z(t₀+ε) < 0 (부호 변환)")
log("  단순 영점 ↔ 모노드로미 ≈ ±π (위상 점프 1회)")
t0 = time.time()

mono_eps = 0.05
mono_results = []
mono_fail_cnt = 0

check_zeros_mono = zeros_nontrivial[:20]
for i, t_zero in enumerate(check_zeros_mono):
    try:
        Za = float(gp(f"lfunhardy(Linit, {t_zero - mono_eps})"))
        Zb = float(gp(f"lfunhardy(Linit, {t_zero + mono_eps})"))
        sign_change = (Za < 0) != (Zb < 0)
        mono_results.append(sign_change)
        icon = "✅" if sign_change else "❌"
        log(f"  [{i+1:02d}] t0={t_zero:.6f}: Z(-ε)={Za:.3e}, Z(+ε)={Zb:.3e} → {icon} {'단순' if sign_change else '비단순?'}")
    except Exception as e:
        log(f"  [{i+1:02d}] t0={t_zero:.6f}: ❌ 오류: {str(e)[:50]}")
        mono_results.append(False)
        mono_fail_cnt += 1

n_simple = sum(mono_results)
n_checked = len(mono_results)
mono_pass = n_checked > 0 and n_simple == n_checked
log(f"  단순 영점: {n_simple}/{n_checked}")
log(f"  판정: {'✅ PASS' if mono_pass else '⚠️ 일부 비단순'}")
log(f"  ({time.time()-t0:.1f}s)")
log()
flush_file()

# =====================================================================
# [6] ★★ κ_near 측정 (Hardy Z 기반, δ=0.001) + A(t₀) 핵심 측정
# =====================================================================
log("=" * 72)
log(f"[7] ★★ κ_near 측정 (Hardy Z 기반, δ={DELTA_KAPPA}) + A(t₀)")
log("=" * 72)
log(f"  핵심: A(t₀) = κ - 1/δ², 비교: sym³(Δ) A={D4_DELTA_A:.2f}, sym²(11a1) A={D3_A}")
t0_time = time.time()

target_zeros = [t for t in zeros_nontrivial if t <= T_KAPPA_MAX]
log(f"  대상 영점: {len(target_zeros)}개 (t ∈ [{T_KAPPA_MIN}, {T_KAPPA_MAX}])")
log()

kappa_results = []
kappa_success = 0
kappa_fail = 0

for i, t_zero in enumerate(target_zeros):
    try:
        # Z'(t₀) — 5점 중심 차분 (4차 정확도)
        h = 1e-6
        Zp2 = float(gp(f"lfunhardy(Linit, {t_zero + 2*h})"))
        Zp1 = float(gp(f"lfunhardy(Linit, {t_zero + h})"))
        Zm1 = float(gp(f"lfunhardy(Linit, {t_zero - h})"))
        Zm2 = float(gp(f"lfunhardy(Linit, {t_zero - 2*h})"))
        Zprime = (-Zp2 + 8*Zp1 - 8*Zm1 + Zm2) / (12 * h)

        # Z(t₀ + δ)
        Z_delta = float(gp(f"lfunhardy(Linit, {t_zero + DELTA_KAPPA})"))

        if abs(Z_delta) < 1e-200:
            log(f"  [{i+1:02d}] t0={t_zero:.6f}: ❌ |Z(t₀+δ)| 너무 작음")
            kappa_fail += 1
            continue

        if abs(Zprime) < 1e-200:
            log(f"  [{i+1:02d}] t0={t_zero:.6f}: ❌ |Z'(t₀)| 너무 작음 (고차 영점?)")
            kappa_fail += 1
            continue

        kappa = (Zprime / Z_delta) ** 2
        kd2 = kappa * DELTA_KAPPA**2
        A = kappa - 1.0 / (DELTA_KAPPA**2)

        kappa_results.append({
            't0': t_zero,
            'kappa': kappa,
            'kd2': kd2,
            'A': A,
        })
        kappa_success += 1
        log(f"  [{i+1:02d}] t0={t_zero:.6f}: κ={kappa:.2f}, A={A:.4f}, κδ²={kd2:.6f}")

    except Exception as e:
        log(f"  [{i+1:02d}] t0={t_zero:.6f}: ❌ 오류: {str(e)[:60]}")
        kappa_fail += 1

log()
log(f"  성공: {kappa_success}, 실패: {kappa_fail}")

mean_A, std_A, cv_A, mean_kd2, median_A = 0, 0, 0, 0, 0
if kappa_results:
    As = [r['A'] for r in kappa_results]
    kd2s = [r['kd2'] for r in kappa_results]

    mean_A = statistics.mean(As)
    std_A = statistics.stdev(As) if len(As) > 1 else 0.0
    cv_A = abs(std_A / mean_A * 100) if abs(mean_A) > 1e-12 else float('inf')
    mean_kd2 = statistics.mean(kd2s)
    median_A = statistics.median(As)

    log(f"\n  ★★ A(t₀) 결과 ★★")
    log(f"  mean(A)   = {mean_A:.4f}")
    log(f"  median(A) = {median_A:.4f}")
    log(f"  std(A)    = {std_A:.4f}")
    log(f"  CV(A)     = {cv_A:.1f}%")
    log(f"  mean(κδ²) = {mean_kd2:.6f} (기대: [0.99, 1.15])")
    log()
    log(f"  ★ 비교 (d=4):")
    log(f"    sym³(Δ)    d=4, w=11: A={D4_DELTA_A:.2f}")
    log(f"    sym³(11a1) d=4, w=1 : A={mean_A:.4f}")
    if abs(mean_A - D4_DELTA_A) < abs(mean_A - D3_A):
        log(f"    → A(sym³(11a1)) ≈ A(sym³(Δ)) → 'degree 효과' 지지")
    elif abs(mean_A - D3_A) < abs(mean_A - D4_DELTA_A):
        log(f"    → A(sym³(11a1)) ≈ A(sym²,d=3) → 'weight 효과' 지지")
    else:
        log(f"    → 중간값 → 혼합 의존")

kd2_pass = bool(kappa_results) and 0.99 <= mean_kd2 <= 1.15
log(f"\n  κδ² 판정: {'✅ PASS' if kd2_pass else '⚠️ FAIL'}")
log(f"  ({time.time()-t0_time:.1f}s)")
log()
flush_file()

# =====================================================================
# [7] ★ δ 독립성 체크 (3 영점 × 4δ)
# =====================================================================
log("=" * 72)
log("[8] ★ δ 독립성 체크 (3영점 × 4δ)")
log("=" * 72)
t0 = time.time()

check_zeros_delta = zeros_nontrivial[:3]
delta_list = [0.0001, 0.001, 0.01, 0.1]

for t_zero in check_zeros_delta:
    try:
        h = 1e-6
        Zp2 = float(gp(f"lfunhardy(Linit, {t_zero + 2*h})"))
        Zp1 = float(gp(f"lfunhardy(Linit, {t_zero + h})"))
        Zm1 = float(gp(f"lfunhardy(Linit, {t_zero - h})"))
        Zm2 = float(gp(f"lfunhardy(Linit, {t_zero - 2*h})"))
        Zprime = (-Zp2 + 8*Zp1 - 8*Zm1 + Zm2) / (12 * h)

        parts = []
        As_delta = []
        for delta in delta_list:
            try:
                Z_d = float(gp(f"lfunhardy(Linit, {t_zero + delta})"))
                if abs(Z_d) > 1e-200 and abs(Zprime) > 1e-200:
                    kd2 = (Zprime / Z_d)**2 * delta**2
                    A_d = (Zprime / Z_d)**2 - 1.0/delta**2
                    parts.append(f"δ={delta}: κδ²={kd2:.6f}, A={A_d:.3f}")
                    As_delta.append(A_d)
            except Exception as e:
                parts.append(f"δ={delta}: ❌")
        log(f"  t0={t_zero:.5f}:")
        for p in parts:
            log(f"    {p}")
        if len(As_delta) >= 2:
            spread = (max(As_delta) - min(As_delta)) / (abs(statistics.mean(As_delta)) + 1e-12) * 100
            log(f"    → A spread: {spread:.2f}% (< 5%이면 δ-독립 ✅)")
    except Exception as e:
        log(f"  t0={t_zero:.5f}: ❌ 오류: {str(e)[:60]}")

log(f"  ({time.time()-t0:.1f}s)")
log()
flush_file()

# =====================================================================
# [8] B-12 단조증가 비교표
# =====================================================================
log("=" * 72)
log("[9] B-12 단조증가 비교표 — degree/weight 분리")
log("=" * 72)

log(f"  {'d':2s}  {'L-함수':20s}  {'weight':7s}  {'A(t₀) mean':12s}  {'판정':6s}")
log(f"  {'-'*2}  {'-'*20}  {'-'*7}  {'-'*12}  {'-'*6}")
log(f"  {'1':2s}  {'ζ(s)':20s}  {'—':7s}  {D1_A:12.4f}  {'기준'}")
log(f"  {'2':2s}  {'GL(2) avg':20s}  {'0/12':7s}  {D2_A:12.4f}  {'기준'}")
log(f"  {'3':2s}  {'sym²(11a1)':20s}  {'1':7s}  {D3_A:12.4f}  {'기준'}")
log(f"  {'4':2s}  {'sym³(Δ)':20s}  {'11':7s}  {D4_DELTA_A:12.2f}  {'#80 ✅'}")
if kappa_results:
    d4_11a1_icon = '★'
    log(f"  {'4':2s}  {'sym³(11a1)':20s}  {'1':7s}  {mean_A:12.4f}  {d4_11a1_icon} (#81)")
    log()
    log(f"  ─── weight 효과 vs degree 효과 판별 ───")
    log(f"  d=4, w=1  sym³(11a1): A = {mean_A:.4f}")
    log(f"  d=4, w=11 sym³(Δ)  : A = {D4_DELTA_A:.2f}")
    log(f"  d=3, w=1  sym²(11a1): A = {D3_A:.4f}")
    diff_same_w = abs(mean_A - D3_A)
    diff_same_d = abs(mean_A - D4_DELTA_A)
    log(f"  |A(sym³,w=1) - A(d=3,w=1)| = {diff_same_w:.4f}")
    log(f"  |A(sym³,w=1) - A(d=4,w=11)| = {diff_same_d:.4f}")
    if diff_same_w < diff_same_d:
        log(f"  → weight 효과 지지: A는 weight(=1)에 따라 결정")
        log(f"    A(d=4,w=11)=286.26의 급증은 weight-11 효과")
    elif diff_same_d < diff_same_w:
        log(f"  → degree 효과 지지: A는 degree(=4)에 따라 결정")
        log(f"    A는 weight 무관, d=4 구조가 원인")
    else:
        log(f"  → 혼합 의존: degree + weight 모두 기여")

log()
flush_file()

# =====================================================================
# 총정리
# =====================================================================
log("=" * 72)
log("성공 기준 총정리")
log("=" * 72)

# 성공 기준 평가
crit_fe = fe_pass
crit_zeros = zeros_pass
crit_kd2 = kd2_pass
crit_key = bool(kappa_results and abs(mean_A - D4_DELTA_A) > 1.0)  # 구별 가능

log(f"  {'✅' if crit_fe else '❌'} FE ≥ 8자리: {fe:.1f}자리")
log(f"  {'✅' if crit_zeros else '❌'} 영점 ≥ 20개 (t∈[0,55]): {len(zeros_nontrivial)}개")
log(f"  {'✅' if crit_kd2 else '❌'} κδ² ∈ [0.99, 1.15]: mean={mean_kd2:.6f}")
log(f"  {'✅' if mono_pass else '❌'} 모노드로미 단순영점: {n_simple}/{n_checked}")
log(f"  {'✅' if crit_key else '⚠️'} A(t₀) degree/weight 판별 가능")

pass_count = sum([crit_fe, crit_zeros, crit_kd2, mono_pass])
log()
log(f"  4성질 통과: {pass_count}/4")
log()

# ─── 수치 요약 ───────────────────────────────────────────────────────────
log("─" * 72)
log("수치 요약")
log("─" * 72)
log(f"  L-함수: lfunsympow(11a1, 3)")
log(f"  파라미터: gammaV=[-1,0,0,1], k=4, N=1331, ε=+1")
log(f"  FE: {fe:.1f}자리")
log(f"  영점: {nz}개 (총), {len(zeros_nontrivial)}개 (비자명)")
if zeros_nontrivial:
    log(f"  t₁ = {zeros_nontrivial[0]:.6f}")
log(f"  모노드로미: {n_simple}/{n_checked} 단순")
if kappa_results:
    log(f"  A(t₀) mean   = {mean_A:.4f}")
    log(f"  A(t₀) median = {median_A:.4f}")
    log(f"  A(t₀) std    = {std_A:.4f}")
    log(f"  A(t₀) CV     = {cv_A:.1f}%")
    log(f"  κδ² mean     = {mean_kd2:.6f}")
    log()
    log(f"  ★★ 핵심 결론 ★★")
    log(f"  sym³(11a1) d=4, w=1: A = {mean_A:.4f}")
    log(f"  sym³(Δ)    d=4, w=11: A = {D4_DELTA_A:.2f}")
    log(f"  sym²(11a1) d=3, w=1: A = {D3_A:.4f}")
    if abs(mean_A - D3_A) < abs(mean_A - D4_DELTA_A):
        log(f"  → ★ weight 효과: A(d=4,w=11) 급증은 weight-11 때문")
    elif abs(mean_A - D4_DELTA_A) < abs(mean_A - D3_A):
        log(f"  → ★ degree 효과: A는 weight 무관, degree-4 구조 원인")
    else:
        log(f"  → ★ 혼합 의존: degree + weight 복합")
log()
log(f"  완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 72)
log(f"결과: {os.path.abspath(OUTFILE)}")
log()

flush_file()
print(f"\n결과 저장: {OUTFILE}")
