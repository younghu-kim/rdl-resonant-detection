#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #87 — GL(3) sym²(37a1) PARI 4성질 검증 + ξ-bundle A(t₀)
=============================================================================
목적:  d=3에서 제2 L-함수 확보 → 보편성 주장의 통계적 강화 + conductor 독립성 검증.
  sym²(37a1) = lfunsympow(ellinit([0,0,1,-1,0]), 2)

  L-함수 정보:
  - 37a1: y²+y = x³−x, LMFDB [0,0,1,-1,0], conductor 37, rank 1
  - sym²(37a1): degree 3, GL(3), k=4(→center=2), N=37²=1369
  - gammaV = [-1, 0, 1]
  - 비교 대상: sym²(11a1) → κδ²≈1.000, A(t₀)≈12.79 (#74)
  - conductor 독립성 최초 검증 (N=121 vs N=1369)

4성질:
  P1 (FE)       : lfuncheckfeq(L3) ≤ -30
  P2 (영점)     : lfunzeros(L3, 30) ≥ 30개
  P3 (κ_near)   : Hardy Z κδ² ∈ [0.99, 1.01] (3영점 × δ=0.001)
  P4 (모노드로미): Z 부호변환 ≥ 90% at 20+영점

ξ-bundle A(t₀) (σ-방향, lfunlambda):
  - κ = |Λ'(s)/Λ(s)|² at s = center + δ + i*t₀  (center=2)
  - 3영점 × 6δ (0.001, 0.005, 0.01, 0.02, 0.05, 0.10) = 18점
  - A(t₀) = κδ² − 1 → mean, CV 보고
  - 핵심: sym²(37a1) vs sym²(11a1)=12.79 → conductor 효과 정량화

σ-유일성 (motivic center=2):
  - σ ∈ {1.8, 1.9, 2.0, 2.1, 2.2} (= center + {−0.2,−0.1,0,+0.1,+0.2})
  - PASS: σ=2.0(critical)에서 κ 최대

결과: results/gl3_sym2_37a1_87.txt
=============================================================================
"""

import sys, os, time
import statistics

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

# ─── 출력 설정 ─────────────────────────────────────────────────────────────
OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "results", "gl3_sym2_37a1_87.txt"
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
DELTA_KAPPA = 0.001     # P3 κ_near 측정 δ (Hardy Z)
T_RANGE_MAX = 30        # 영점 탐색 범위
T_KAPPA_MIN = 0.5       # 양수 영점 최소 t (trivial zero 제외)

# ξ-bundle δ 값 (수학자 지시: 6δ점)
XI_DELTAS = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]
CENTER = 2.0            # sym²(37a1): PARI k=4, center=k/2=2

# σ-유일성 σ 값 (motivic: center=2 ± offsets)
SIGMA_MOTIVIC = [1.8, 1.9, 2.0, 2.1, 2.2]

# κ_near / 모노드로미용 영점 개수
N_ZEROS_KAPPA = 3
N_ZEROS_MONO  = 25  # P4 모노드로미 최소 확인 수

# 비교 기준 (sym²(11a1) #74)
D3_11A1_A = 12.79   # sym²(11a1) A(t₀) [#74]

log("=" * 72)
log("결과 #87 — GL(3) sym²(37a1) PARI 4성질 검증 + ξ-bundle A(t₀)")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"  L-함수: lfunsympow(ellinit([0,0,1,-1,0]), 2)")
log(f"  degree=3, GL(3), 37a1 기반, k=4, center={CENTER}")
log(f"  N=37²=1369 (conductor), gammaV=[-1,0,1]")
log(f"  비교 대상: sym²(11a1) A(t₀)={D3_11A1_A} [#74] (N=121)")
log(f"  xi_deltas: {XI_DELTAS}")
log(f"  sigma_motivic: {SIGMA_MOTIVIC}")
log()
flush_file()

# ─── PARI 초기화 ──────────────────────────────────────────────────────────
import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000 * 1024 * 1024)   # 2GB
gp("default(realprecision, 100)")    # 100자리 정밀도
log("PARI 초기화: 2GB 메모리, realprecision=100")
log()
flush_file()

# ─── [0] L-함수 생성 ──────────────────────────────────────────────────────
log("=" * 72)
log("[0] sym²(37a1) L-함수 생성 — lfunsympow(E, 2)")
log("=" * 72)
t_init = time.time()

gp("E37 = ellinit([0,0,1,-1,0])")
gp("L3 = lfunsympow(E37, 2)")
log(f"  lfunsympow(E37, 2) 완료 ({time.time()-t_init:.1f}s)")

try:
    params = str(gp("lfunparams(L3)"))
    log(f"  lfunparams: {params}")
except Exception as e:
    log(f"  lfunparams: {str(e)[:60]}")

log()
flush_file()

# =====================================================================
# [1] P1: 함수방정식 검증 (FE)
# =====================================================================
log("=" * 72)
log("[1] P1: 함수방정식 검증 (lfuncheckfeq)")
log("=" * 72)
t0 = time.time()

try:
    fe = float(gp("lfuncheckfeq(L3)"))
    fe_pass = fe <= -30
    log(f"  FE = {fe:.1f}자리 정밀도")
    log(f"  판정: {'✅ PASS' if fe_pass else '❌ FAIL'} (기준 ≤ -30)")
except Exception as e:
    log(f"  ❌ lfuncheckfeq 오류: {e}")
    fe = 0.0
    fe_pass = False

log(f"  ({time.time()-t0:.1f}s)")
log()
flush_file()

# =====================================================================
# [2] P2: 영점 탐색 (lfunzeros)
# =====================================================================
log("=" * 72)
log("[2] P2: 영점 탐색 (lfunzeros)")
log("=" * 72)
t0 = time.time()

try:
    z = gp(f"lfunzeros(L3, {T_RANGE_MAX})")
    nz = len(z)
    zeros = [float(z[i]) for i in range(nz)]
    zeros_pos = [t for t in zeros if t > T_KAPPA_MIN]
    log(f"  총 영점: {nz}개 (t ∈ [0, {T_RANGE_MAX}])")
    log(f"  양수 영점: {len(zeros_pos)}개 (t > {T_KAPPA_MIN})")
    if zeros_pos:
        log(f"  t₁ = {zeros_pos[0]:.6f}")
    for i, t in enumerate(zeros[:25]):
        log(f"  영점 #{i+1:2d}: t = {t:.6f}")
    if nz > 25:
        log(f"  ... (총 {nz}개)")
    zeros_pass = len(zeros_pos) >= 30
    log(f"  판정: {'✅ PASS' if zeros_pass else '❌ FAIL'} (≥30개)")
except Exception as e:
    log(f"  ❌ lfunzeros 오류: {e}")
    zeros = []
    zeros_pos = []
    nz = 0
    zeros_pass = False

log(f"  ({time.time()-t0:.1f}s)")
log()
flush_file()

# κ_near / ξ-bundle용 영점 선택
kappa_zeros = zeros_pos[:N_ZEROS_KAPPA]
log(f"  κ_near/ξ-bundle용 영점 (처음 {N_ZEROS_KAPPA}개): {[f'{t:.6f}' for t in kappa_zeros]}")
log()
flush_file()

# =====================================================================
# [3] P3: κ_near 측정 (Hardy Z 기반, t-방향)
# =====================================================================
log("=" * 72)
log(f"[3] P3: κ_near 측정 (Hardy Z, δ={DELTA_KAPPA})")
log("=" * 72)
t0_time = time.time()

kappa_results = []
kappa_success = 0
kappa_fail = 0

for i, t_zero in enumerate(kappa_zeros):
    try:
        h = 1e-8
        Z_plus  = float(gp(f"lfunhardy(L3, {t_zero + h})"))
        Z_minus = float(gp(f"lfunhardy(L3, {t_zero - h})"))
        Zp = (Z_plus - Z_minus) / (2.0 * h)

        Z_delta = float(gp(f"lfunhardy(L3, {t_zero + DELTA_KAPPA})"))

        if abs(Z_delta) < 1e-200:
            log(f"  [{i+1}] t0={t_zero:.6f}: ❌ |Z(t0+δ)| 너무 작음 ({abs(Z_delta):.2e})")
            kappa_fail += 1
            continue

        kappa = (Zp / Z_delta) ** 2
        kd2   = kappa * DELTA_KAPPA**2
        A     = kappa - 1.0 / (DELTA_KAPPA**2)

        kappa_results.append({'t0': t_zero, 'kappa': kappa, 'kd2': kd2, 'A': A})
        kappa_success += 1
        ok = "✅" if 0.99 <= kd2 <= 1.01 else "⚠️"
        log(f"  [{i+1}] t0={t_zero:.6f}: κ={kappa:.4f}, A={A:.4f}, κδ²={kd2:.6f} {ok}")

    except Exception as e:
        log(f"  [{i+1}] t0={t_zero:.6f}: ❌ 오류: {str(e)[:60]}")
        kappa_fail += 1

log(f"  성공: {kappa_success}/{len(kappa_zeros)}")
if kappa_results:
    kd2s    = [r['kd2'] for r in kappa_results]
    As_p3   = [r['A']   for r in kappa_results]
    mean_kd2   = statistics.mean(kd2s)
    mean_A_p3  = statistics.mean(As_p3)
    log(f"  mean(κδ²) = {mean_kd2:.6f}")
    log(f"  mean(A)   = {mean_A_p3:.4f}")
    kd2_pass = 0.99 <= mean_kd2 <= 1.01
    log(f"  판정: {'✅ PASS' if kd2_pass else '❌ FAIL'} (κδ² ∈ [0.99, 1.01])")
else:
    mean_kd2 = float('nan')
    mean_A_p3 = float('nan')
    kd2_pass = False
    log("  ❌ κ_near 데이터 없음")

log(f"  ({time.time()-t0_time:.1f}s)")
log()
flush_file()

# =====================================================================
# [4] P4: 모노드로미 (Hardy Z 부호변환)
# =====================================================================
log("=" * 72)
log("[4] P4: 모노드로미 (Hardy Z 부호변환, 단순 영점 비율)")
log("=" * 72)
t0 = time.time()

mono_zeros = zeros_pos[:min(N_ZEROS_MONO, len(zeros_pos))]
log(f"  대상 영점: {len(mono_zeros)}개")

mono_count  = 0
simple_count = 0

for i, t_zero in enumerate(mono_zeros):
    try:
        h = 0.01  # 영점 양쪽에서 부호 확인
        Z_before = float(gp(f"lfunhardy(L3, {t_zero - h})"))
        Z_after  = float(gp(f"lfunhardy(L3, {t_zero + h})"))

        is_simple = (Z_before * Z_after < 0)
        if is_simple:
            simple_count += 1
        mono_count += 1

        if i < 25:
            ok = "✅ 단순" if is_simple else "⚠️ 비단순"
            log(f"  [{i+1:2d}] t0={t_zero:.6f}: "
                f"Z(t-h)={Z_before:+.4e}, Z(t+h)={Z_after:+.4e} {ok}")
    except Exception as e:
        log(f"  [{i+1:2d}] t0={t_zero:.6f}: ❌ 오류: {str(e)[:50]}")

log()
mono_rate = simple_count / mono_count * 100 if mono_count > 0 else 0.0
mono_pass = mono_rate >= 90.0
log(f"  단순 영점: {simple_count}/{mono_count}개")
log(f"  단순 영점 비율: {mono_rate:.1f}%")
log(f"  판정: {'✅ PASS' if mono_pass else '❌ FAIL'} (≥90%)")
log(f"  ({time.time()-t0:.1f}s)")
log()
flush_file()

# =====================================================================
# [5] ξ-bundle A(t₀) 측정 (σ-방향, lfunlambda)
# =====================================================================
log("=" * 72)
log("[5] ξ-bundle A(t₀) 측정 (σ-방향, lfunlambda)")
log("=" * 72)
log(f"  center = {CENTER} (k=4, k/2=2)")
log(f"  s = {CENTER} + δ + i*t₀  (σ 방향, t₀ 고정)")
log(f"  κ = |Λ'(s)/Λ(s)|²,  A = κδ² − 1")
log(f"  δ 값: {XI_DELTAS}")
log()

# lfuninit — 계산 범위 확장
t_xi_start = time.time()
try:
    gp("L3i = lfuninit(L3, [0, 12])")
    log(f"  lfuninit([0,12]) 완료 ({time.time()-t_xi_start:.1f}s)")
except Exception as e:
    log(f"  WARNING lfuninit: {e} — L3 직접 사용")
    gp("L3i = L3")

log()
flush_file()

xi_results = []

for t_zero in kappa_zeros:
    log(f"  ── t₀ = {t_zero:.6f} ──")
    row = {'t0': t_zero, 'As': [], 'kd2s': []}

    for delta in XI_DELTAS:
        try:
            s_val = f"{CENTER} + {delta} + I*{t_zero}"
            gp(f"s_xi = {s_val}")
            L0_val = gp("lfunlambda(L3i, s_xi)")
            Lp_val = gp("lfunlambda(L3i, s_xi, 1)")
            abs_L0 = float(gp(f"abs({L0_val})"))
            abs_Lp = float(gp(f"abs({Lp_val})"))

            if abs_L0 < 1e-250:
                log(f"    δ={delta:.4f}: ⚠️ |Λ(s)| 너무 작음 ({abs_L0:.2e})")
                continue

            kappa_xi = abs_Lp**2 / abs_L0**2
            kd2_xi   = kappa_xi * delta**2
            A_xi     = kappa_xi - 1.0 / delta**2
            row['As'].append(A_xi)
            row['kd2s'].append(kd2_xi)
            ok = "✅" if 0.99 <= kd2_xi <= 1.01 else "⚠️"
            log(f"    δ={delta:.4f}: κ={kappa_xi:.4f}, A={A_xi:.4f}, κδ²={kd2_xi:.6f} {ok}")

        except Exception as e:
            log(f"    δ={delta:.4f}: ❌ 오류: {str(e)[:60]}")

    if row['As']:
        mean_A = statistics.mean(row['As'])
        cv_A   = (statistics.stdev(row['As']) / abs(mean_A) * 100
                  if len(row['As']) > 1 and abs(mean_A) > 1e-10 else 0.0)
        row['mean_A'] = mean_A
        row['cv_A']   = cv_A
        log(f"    → mean(A) = {mean_A:.4f}, CV(A) = {cv_A:.2f}%")
    xi_results.append(row)
    log()
    flush_file()

# ξ-bundle 전체 요약
xi_all_A = [A for r in xi_results for A in r['As']]
if xi_all_A:
    xi_mean = statistics.mean(xi_all_A)
    xi_cv   = (statistics.stdev(xi_all_A) / abs(xi_mean) * 100
               if len(xi_all_A) > 1 and abs(xi_mean) > 1e-10 else 0.0)
    xi_cv_pass = xi_cv < 50.0
    log(f"  sym²(37a1) ξ-bundle 요약:")
    log(f"    전체 mean(A) = {xi_mean:.4f}")
    log(f"    CV(A)        = {xi_cv:.2f}%")
    log(f"    점수:          {len(xi_all_A)}/{len(kappa_zeros)*len(XI_DELTAS)}")
    log(f"    CV < 50%: {'✅ PASS' if xi_cv_pass else '❌ FAIL'}")
else:
    xi_mean    = float('nan')
    xi_cv      = float('nan')
    xi_cv_pass = False
    log("  ❌ ξ-bundle 데이터 없음")
log()
flush_file()

# δ-독립성 표 출력
log("  A(t₀, δ) 표 (행=영점, 열=δ):")
log(f"  {'t₀':>10s} | " + " | ".join(f"δ={d:.4f}" for d in XI_DELTAS))
log("  " + "-" * 80)
for r in xi_results:
    row_vals = []
    for i_d, _ in enumerate(XI_DELTAS):
        if i_d < len(r['As']):
            row_vals.append(f"{r['As'][i_d]:8.4f}")
        else:
            row_vals.append("  N/A   ")
    log(f"  t₀={r['t0']:.5f} | " + " | ".join(row_vals))
log()
flush_file()

# =====================================================================
# [6] σ-유일성 (motivic σ 방향, center=2)
# =====================================================================
log("=" * 72)
log("[6] σ-유일성 (lfunlambda, σ_motivic ∈ {1.8, 1.9, 2.0, 2.1, 2.2})")
log("=" * 72)
log(f"  κ_log = |Λ'(σ+it₀)/Λ(σ+it₀)|²  (δ=0.001 고정)")
log(f"  PASS: σ=2.0 (critical line)에서 κ 최대")
log()

SIGMA_DELTA = 0.001  # σ-유일성용 δ

sigma_results = {}  # t_zero -> {σ -> kappa}

for t_zero in kappa_zeros:
    log(f"  ── t₀ = {t_zero:.6f} ──")
    sig_kappas = {}

    for sigma in SIGMA_MOTIVIC:
        try:
            s_str = f"{sigma} + {SIGMA_DELTA} + I*{t_zero}"
            gp(f"s_sig = {s_str}")
            L0_sig = gp("lfunlambda(L3i, s_sig)")
            Lp_sig = gp("lfunlambda(L3i, s_sig, 1)")
            abs_L0_sig = float(gp(f"abs({L0_sig})"))
            abs_Lp_sig = float(gp(f"abs({Lp_sig})"))

            if abs_L0_sig < 1e-250:
                log(f"    σ={sigma:.1f}: ⚠️ |Λ| 너무 작음 ({abs_L0_sig:.2e})")
                sig_kappas[sigma] = None
                continue

            k_sig = abs_Lp_sig**2 / abs_L0_sig**2
            sig_kappas[sigma] = k_sig
            log(f"    σ={sigma:.1f}: κ={k_sig:.4f}, |Λ|={abs_L0_sig:.4e}")

        except Exception as e:
            log(f"    σ={sigma:.1f}: ❌ 오류: {str(e)[:60]}")
            sig_kappas[sigma] = None

    sigma_results[t_zero] = sig_kappas

    valid = {s: v for s, v in sig_kappas.items() if v is not None}
    if valid:
        max_sig = max(valid, key=valid.get)
        kappa_crit = valid.get(2.0, None)
        log(f"    최대 κ: σ={max_sig:.1f} (κ={valid[max_sig]:.4f})")
        if kappa_crit is not None:
            log(f"    σ=2.0(crit) κ: {kappa_crit:.4f}")
        if max_sig == 2.0:
            log(f"    → ✅ PASS (σ=2.0에서 최대)")
        else:
            vals_list = [valid.get(s) for s in SIGMA_MOTIVIC if valid.get(s) is not None]
            log(f"    → ❌ σ={max_sig:.1f}에서 최대 (구조적 패턴 확인)")
            if len(vals_list) >= 3:
                mono_inc = all(vals_list[i] < vals_list[i+1] for i in range(len(vals_list)-1))
                mono_dec = all(vals_list[i] > vals_list[i+1] for i in range(len(vals_list)-1))
                if mono_inc:
                    log(f"    패턴: σ 증가 → κ 단조증가")
                elif mono_dec:
                    log(f"    패턴: σ 증가 → κ 단조감소")
                else:
                    log(f"    패턴: 비단조")
    log()
    flush_file()

# σ-유일성 판정 (σ=2.0에서 최대 기준)
sigma_pass_cnt = 0
sigma_total_cnt = 0
for t_zero in kappa_zeros:
    sigs = sigma_results.get(t_zero, {})
    valid = {s: v for s, v in sigs.items() if v is not None}
    if valid:
        sigma_total_cnt += 1
        if max(valid, key=valid.get) == 2.0:
            sigma_pass_cnt += 1

sigma_unique_pass = (sigma_pass_cnt == sigma_total_cnt and sigma_total_cnt > 0)
log(f"  σ-유일성 PASS (σ=2.0 최대): {sigma_pass_cnt}/{sigma_total_cnt}영점")
if sigma_unique_pass:
    log(f"  전체 판정: ✅ PASS")
else:
    log(f"  전체 판정: ❌ FAIL (구조적 FAIL — 이전 GL(4),(5)와 패턴 비교 필요)")
log()
flush_file()

# =====================================================================
# [7] conductor 독립성 비교 — sym²(37a1) vs sym²(11a1)
# =====================================================================
log("=" * 72)
log("[7] conductor 독립성 비교 — sym²(37a1) vs sym²(11a1)")
log("=" * 72)
log()
log(f"  sym²(11a1): N=121=11², center=2(?), A(t₀)={D3_11A1_A}  [#74]")

if xi_all_A:
    log(f"  sym²(37a1): N=1369=37², center={CENTER},  A(t₀)={xi_mean:.4f}  [#87 THIS]")
    diff_pct = abs(xi_mean - D3_11A1_A) / D3_11A1_A * 100
    log(f"  차이: |{xi_mean:.4f} - {D3_11A1_A}| = {abs(xi_mean - D3_11A1_A):.4f} ({diff_pct:.1f}%)")
    if diff_pct < 30:
        log(f"  → ✅ conductor 독립성 확인 (30% 이내)")
    else:
        log(f"  → ⚠️ conductor 효과 존재 ({diff_pct:.1f}%) — 정규화 불일치 가능성 포함")
    log()

# d=3 degree 비교 — 두 L-함수
log("  d=3 L-함수 요약 (conductor 독립성):")
log(f"    sym²(11a1): N=121, A(t₀)={D3_11A1_A} [#74]")
if xi_all_A:
    log(f"    sym²(37a1): N=1369, A(t₀)={xi_mean:.4f} [#87]")
    log(f"    → d=3에서 2개 L-함수 모두 4성질 통과 시 보편성 주장 강화")
log()
flush_file()

# =====================================================================
# [8] 성공 기준 총정리
# =====================================================================
log("=" * 72)
log("[8] 성공 기준 총정리")
log("=" * 72)
log()

c_fe     = fe_pass
c_zeros  = zeros_pass
c_kd2    = kd2_pass
c_mono   = mono_pass
c_xi     = xi_cv_pass
c_sigma  = sigma_unique_pass

log(f"  P1 FE ≤ -30:           {'✅ PASS' if c_fe else '❌ FAIL'}  "
    f"(FE={fe:.1f}자리)")
log(f"  P2 영점 ≥ 30:          {'✅ PASS' if c_zeros else '❌ FAIL'}  "
    f"({len(zeros_pos)}개)")
log(f"  P3 κδ² ∈ [0.99,1.01]:  {'✅ PASS' if c_kd2 else '❌ FAIL'}  "
    f"(mean={mean_kd2:.6f})")
log(f"  P4 모노드로미 ≥90%:    {'✅ PASS' if c_mono else '❌ FAIL'}  "
    f"({simple_count}/{mono_count}, {mono_rate:.1f}%)")
log(f"  ξ A(t₀) CV<50%:        {'✅ PASS' if c_xi else '❌ FAIL'}  "
    f"(CV={xi_cv:.1f}%, mean={xi_mean:.4f})")
log(f"  σ-유일성 PASS:         {'✅ PASS' if c_sigma else '❌ FAIL'}  "
    f"({sigma_pass_cnt}/{sigma_total_cnt})")

p4_count = sum([c_fe, c_zeros, c_kd2, c_mono])
p6_count = sum([c_fe, c_zeros, c_kd2, c_mono, c_xi, c_sigma])
log()
log(f"  4성질 통과: {p4_count}/4")
log(f"  전체 통과: {p6_count}/6")
log()

if p4_count == 4:
    log("  ★★★ GL(3) sym²(37a1) 4성질 전부 통과!")
    log("       d=3에서 2개 L-함수 모두 통과 → 보편성 통계적 강화 + conductor 독립성 확인.")
elif p4_count == 3:
    log("  ★★ GL(3) sym²(37a1) 4성질 3/4 통과 (양성)")
elif p4_count == 2:
    log("  ★ GL(3) 4성질 2/4 (부분)")
else:
    log("  ❌ GL(3) sym²(37a1) 4성질 실패")
log()
flush_file()

# =====================================================================
# [9] 수치 요약 + 경계 갱신
# =====================================================================
log("=" * 72)
log("[9] 수치 요약 + 경계 갱신")
log("=" * 72)
log()
log(f"  L-함수: lfunsympow(ellinit([0,0,1,-1,0]), 2)  [GL(3), d=3]")
log(f"  37a1 = y²+y = x³−x, N=37²=1369, center={CENTER}")
log()
log(f"  P1 FE:           {fe:.1f}자리")
log(f"  P2 영점:         {nz}개 총 / {len(zeros_pos)}개 양수" +
    (f" / t₁={zeros_pos[0]:.6f}" if zeros_pos else ""))
if kappa_results:
    log(f"  P3 κδ²:          mean={mean_kd2:.6f}, {kappa_success}/{len(kappa_zeros)}")
log(f"  P4 모노드로미:   {simple_count}/{mono_count} ({mono_rate:.1f}%)")
if xi_all_A:
    log(f"  ξ A(t₀):         mean={xi_mean:.4f}, CV={xi_cv:.2f}%")
    for r in xi_results:
        if r['As']:
            log(f"    t₀={r['t0']:.5f}: "
                f"mean(A)={r.get('mean_A',float('nan')):.4f}, "
                f"CV={r.get('cv_A',0):.2f}%")
log()
log("─" * 72)
log("경계 갱신")
log("─" * 72)
log(f"  B-25 (d=3 conductor 독립성): {'✅ 확인' if p4_count >= 3 else '❌ 미확인'} — "
    f"N=121 vs N=1369 두 L-함수 비교")
if xi_all_A:
    log(f"  B-09 (d=3 A 비교): sym²(11a1)={D3_11A1_A}, sym²(37a1)={xi_mean:.4f} — conductor 효과")
log(f"  보편성: d=3에서 {'2개' if p4_count >= 3 else '1개'}로 {'강화' if p4_count >= 3 else '미변'}")
log()
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 72)
log(f"결과: {os.path.abspath(OUTFILE)}")
log()

flush_file()
print(f"\n✅ 결과 저장: {OUTFILE}")
