#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] #82 — δ-독립성 진단 + Richardson 외삽
=============================================================================
목표: sym³(Δ) + sym³(11a1)에서 A(δ) = A_true + c/δ 패턴 검증
     Hardy Z 기반 κ에서 1/δ 보정항 c가 L-함수 의존적인지 보편적인지 판별

핵심 질문:
  #81에서 sym³(11a1): A ≈ 706/δ + c/δ + A_true 패턴 발견 (c≈0.0705)
  sym³(Δ)에서도 같은 패턴이면 → Hardy Z 방식의 보편 보정항
  다르면 → weight/conductor에 의존하는 새 수학

데이터 소스:
  - sym³(11a1): #81 결과 하드코딩 (3영점 × 4δ)
  - sym³(Δ): 이 스크립트에서 측정 (3영점 × 4δ)
  - ζ(s): 보너스 (3영점 × 4δ)

결과: results/delta_independence_82.txt
=============================================================================
"""

import sys, os, time
import statistics

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

# ─── 출력 설정 ────────────────────────────────────────────────────────────
OUTFILE = os.path.expanduser(
    "~/Desktop/gdl_unified/results/delta_independence_82.txt"
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
DELTAS = [0.0001, 0.001, 0.01, 0.1]
N_COEFF = 4000

# sym³(Δ) PARI 파라미터 (gammaV=[−11,−10,0,1], k=34, N=1, ε=-1)
SYM3D_ZEROS = [4.155866, 5.549122, 8.111776]  # pari_80.txt에서
SYM3D_GAMMA_V = [-11, -10, 0, 1]
SYM3D_K = 34
SYM3D_N = 1
SYM3D_EPS = -1

# sym³(11a1) — #81 결과 하드코딩 (3영점 × 4δ)
# [8]절 데이터
SYM3_11A1_DATA = {
    2.320021: {0.0001: 706.714, 0.001: 71.800, 0.01: 8.321,  0.1: 2.037},
    3.591881: {0.0001: -2674.139, 0.001: -265.133, 0.01: -24.234, 0.1: -0.125},
    4.622627: {0.0001: -5723.701, 0.001: -569.187, 0.01: -53.749, 0.1: -2.270},
}

log("=" * 72)
log("결과 #82 — δ-독립성 진단 + Richardson 외삽")
log("  sym³(Δ) vs sym³(11a1) 1/δ 보정항 비교")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"δ값: {DELTAS}")
log(f"sym³(Δ) 영점: {SYM3D_ZEROS}")
log()
flush_file()

# ─── Richardson 외삽 함수 ─────────────────────────────────────────────────
def richardson_fit(delta_A_dict):
    """
    A(δ) = A_true + c/δ 모델 피팅 (최소제곱)
    delta_A_dict: {δ: A_measured}
    반환: (c, A_true, residuals, fit_quality)
    """
    deltas = sorted(delta_A_dict.keys())
    A_vals = [delta_A_dict[d] for d in deltas]

    if len(deltas) < 2:
        return None, None, None, "데이터 부족"

    # 선형 회귀: A = A_true + c * (1/δ)
    # y = A, x = 1/δ → 선형 회귀
    xs = [1.0/d for d in deltas]
    ys = A_vals
    n = len(xs)
    sx = sum(xs); sy = sum(ys)
    sxx = sum(x*x for x in xs); sxy = sum(x*y for x, y in zip(xs, ys))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-15:
        return None, None, None, "denom≈0"
    c = (n * sxy - sx * sy) / denom
    A_true = (sy - c * sx) / n

    # 잔차
    preds = [A_true + c / d for d in deltas]
    residuals = [abs(p - a) / (abs(a) + 1e-10) * 100 for p, a in zip(preds, A_vals)]

    return c, A_true, residuals, "OK"

def two_point_richardson(d1, A1, d2, A2):
    """2점 Richardson: c = (A1-A2)/(1/d1-1/d2), A_true = A2 - c/d2"""
    inv_d1 = 1.0/d1; inv_d2 = 1.0/d2
    if abs(inv_d1 - inv_d2) < 1e-15:
        return None, None
    c = (A1 - A2) / (inv_d1 - inv_d2)
    A_true = A2 - c * inv_d2
    return c, A_true

# =====================================================================
# [1] sym³(11a1) — #81 데이터 Richardson 외삽
# =====================================================================
log("=" * 72)
log("[1] sym³(11a1) — #81 데이터 Richardson 외삽")
log("=" * 72)
log("  (원본 데이터: #81 [8]절, 재측정 없음)")
log()

sym3_11a1_results = {}
for t0, data in SYM3_11A1_DATA.items():
    log(f"  t₀ = {t0:.6f}")
    for d in DELTAS:
        log(f"    δ={d}: A={data[d]:.3f}")

    c, A_true, residuals, status = richardson_fit(data)
    if status == "OK":
        log(f"  → Richardson 피팅 (4점):")
        log(f"    c = {c:.6f}")
        log(f"    A_true = {A_true:.4f}")
        log(f"    잔차: {[f'{r:.1f}%' for r in residuals]}")
        max_res = max(residuals)
        log(f"    최대 잔차: {max_res:.1f}% ({'✅ <5%' if max_res < 5 else '⚠️ ≥5%'})")
    else:
        log(f"  → 피팅 실패: {status}")
        c, A_true = None, None

    # 2점 검증 (δ=0.0001, 0.001)
    c2, At2 = two_point_richardson(0.0001, data[0.0001], 0.001, data[0.001])
    log(f"  → 2점 Richardson (δ=0.0001,0.001): c={c2:.6f}, A_true={At2:.4f}")
    log(f"     검증 δ=0.01: 예측={At2+c2/0.01:.3f}, 실제={data[0.01]:.3f}")
    log(f"     검증 δ=0.1:  예측={At2+c2/0.1:.3f}, 실제={data[0.1]:.3f}")

    sym3_11a1_results[t0] = {'c': c2, 'A_true': At2, 'A_obs': data}
    log()

flush_file()

# =====================================================================
# [2] sym³(Δ) — PARI lfuninit + 3영점 × 4δ 측정
# =====================================================================
log("=" * 72)
log("[2] sym³(Δ) — PARI lfun 초기화 + δ-스윕")
log("=" * 72)
log()

# PARI 초기화
import cypari2
gp = cypari2.Pari()
gp.allocatemem(1024 * 1024 * 1024)  # 1GB
gp("default(realprecision, 100)")

log("  PARI 초기화 완료 (realprecision=100, 1GB)")

# sym³(Δ) 계수 생성
log("  direuler로 sym³(Δ) 계수 생성 중...")
t_build = time.time()

gp(f"""
an_d = direuler(p=2, {N_COEFF},
  my(t=ramanujantau(p), q=p^11,
     e1=t*(t^2-2*q),
     e2=q*(t^2-2*q)*(t^2-q),
     e3=q^3*e1,
     e4=q^6);
  1/(1 - e1*X + e2*X^2 - e3*X^3 + e4*X^4)
)
""")

gv_str = str(SYM3D_GAMMA_V).replace(' ', '')
gp(f"Ld_d = lfuncreate([an_d, 0, {gv_str}, {SYM3D_K}, {SYM3D_N}, {SYM3D_EPS}, []])")

c2_check = int(gp("an_d[2]"))
log(f"  계수 생성 완료 ({time.time()-t_build:.1f}s)")
log(f"  c(2) = {c2_check} (기대 84480)")
log()
flush_file()

# 3영점 × 4δ 측정
sym3d_results = {}
h_diff = 1e-8  # 수치 미분 스텝

log("  sym³(Δ) 영점별 δ-스윕 시작:")
log()

for t_zero in SYM3D_ZEROS:
    log(f"  t₀ = {t_zero:.6f}")
    t_loop = time.time()

    # Z'(t₀) 수치 미분 (공통)
    try:
        Z_plus_h  = float(gp(f"lfunhardy(Ld_d, {t_zero + h_diff})"))
        Z_minus_h = float(gp(f"lfunhardy(Ld_d, {t_zero - h_diff})"))
        Zp = (Z_plus_h - Z_minus_h) / (2 * h_diff)
        log(f"    Z'(t₀) = {Zp:.6e}")
    except Exception as e:
        log(f"    ❌ Z' 계산 실패: {e}")
        sym3d_results[t_zero] = None
        continue

    delta_A = {}
    for delta in DELTAS:
        try:
            Z_delta = float(gp(f"lfunhardy(Ld_d, {t_zero + delta})"))
            if abs(Z_delta) < 1e-200:
                log(f"    δ={delta}: ❌ |Z(t₀+δ)| 너무 작음")
                continue
            kappa = (Zp / Z_delta) ** 2
            kd2 = kappa * delta**2
            A = kappa - 1.0 / (delta**2)
            delta_A[delta] = A
            log(f"    δ={delta}: κδ²={kd2:.6f}, A={A:.3f}")
        except Exception as e:
            log(f"    δ={delta}: ❌ 오류: {e}")

    log(f"    ({time.time()-t_loop:.1f}s)")

    # Richardson 외삽
    if len(delta_A) >= 3:
        c, A_true, residuals, status = richardson_fit(delta_A)
        if status == "OK":
            log(f"  → Richardson 피팅 ({len(delta_A)}점):")
            log(f"    c = {c:.6f}")
            log(f"    A_true = {A_true:.4f}")
            log(f"    잔차: {[f'{r:.1f}%' for r in residuals]}")
            max_res = max(residuals)
            log(f"    최대 잔차: {max_res:.1f}% ({'✅ <5%' if max_res < 5 else '⚠️ ≥5%'})")
        else:
            c, A_true = None, None
            log(f"  → 피팅 실패: {status}")

        # 2점 Richardson (δ=0.0001, 0.001)
        if 0.0001 in delta_A and 0.001 in delta_A:
            c2, At2 = two_point_richardson(0.0001, delta_A[0.0001], 0.001, delta_A[0.001])
            log(f"  → 2점 Richardson (δ=0.0001,0.001): c={c2:.6f}, A_true={At2:.4f}")
            for d_verify in [0.01, 0.1]:
                if d_verify in delta_A:
                    pred = At2 + c2/d_verify
                    log(f"     검증 δ={d_verify}: 예측={pred:.3f}, 실제={delta_A[d_verify]:.3f}, 오차={abs(pred-delta_A[d_verify])/max(abs(delta_A[d_verify]),1e-10)*100:.1f}%")
            sym3d_results[t_zero] = {'c': c2, 'A_true': At2, 'A_obs': delta_A}
        else:
            sym3d_results[t_zero] = {'c': c, 'A_true': A_true, 'A_obs': delta_A}
    else:
        log(f"  ⚠️ δ 데이터 부족 ({len(delta_A)}점)")
        sym3d_results[t_zero] = {'c': None, 'A_true': None, 'A_obs': delta_A}

    log()
    flush_file()

# =====================================================================
# [3] 보너스 — ζ(s) 3영점 × 4δ (d=1 기준선)
# =====================================================================
log("=" * 72)
log("[3] 보너스 — ζ(s) δ-스윕 (d=1 기준)")
log("=" * 72)

# ζ(s) 첫 3영점
ZETA_ZEROS = [14.134725, 21.022040, 25.010858]
log(f"  영점: {ZETA_ZEROS}")
log()

# ζ(s) PARI L-function (GL(1), 자명)
gp("Ld_zeta = lfuncreate(1)")  # Riemann zeta

zeta_results = {}
for t_zero in ZETA_ZEROS:
    log(f"  t₀ = {t_zero:.6f}")
    try:
        Z_plus_h  = float(gp(f"lfunhardy(Ld_zeta, {t_zero + h_diff})"))
        Z_minus_h = float(gp(f"lfunhardy(Ld_zeta, {t_zero - h_diff})"))
        Zp = (Z_plus_h - Z_minus_h) / (2 * h_diff)
        log(f"    Z'(t₀) = {Zp:.6e}")
    except Exception as e:
        log(f"    ❌ Z' 계산 실패: {e}")
        zeta_results[t_zero] = None
        continue

    delta_A = {}
    for delta in DELTAS:
        try:
            Z_delta = float(gp(f"lfunhardy(Ld_zeta, {t_zero + delta})"))
            if abs(Z_delta) < 1e-200:
                log(f"    δ={delta}: ❌ |Z(t₀+δ)| 너무 작음")
                continue
            kappa = (Zp / Z_delta) ** 2
            kd2 = kappa * delta**2
            A = kappa - 1.0 / (delta**2)
            delta_A[delta] = A
            log(f"    δ={delta}: κδ²={kd2:.6f}, A={A:.3f}")
        except Exception as e:
            log(f"    δ={delta}: ❌ 오류: {e}")

    # 2점 Richardson
    if 0.0001 in delta_A and 0.001 in delta_A:
        c2, At2 = two_point_richardson(0.0001, delta_A[0.0001], 0.001, delta_A[0.001])
        log(f"  → 2점 Richardson: c={c2:.6f}, A_true={At2:.4f}")
        for d_verify in [0.01, 0.1]:
            if d_verify in delta_A:
                pred = At2 + c2/d_verify
                log(f"     검증 δ={d_verify}: 예측={pred:.3f}, 실제={delta_A[d_verify]:.3f}")
        zeta_results[t_zero] = {'c': c2, 'A_true': At2, 'A_obs': delta_A}
    else:
        zeta_results[t_zero] = {'c': None, 'A_true': None, 'A_obs': delta_A}
    log()

flush_file()

# =====================================================================
# [4] 비교표 + 핵심 판정
# =====================================================================
log("=" * 72)
log("[4] ★★ 핵심 비교표 + Richardson 판정")
log("=" * 72)
log()

log("  ─── c값 비교 (1/δ 보정항) ───")
log(f"  {'L-함수':<22} {'영점 t₀':<12} {'c (Richardson)':<18} {'A_true':<12} {'판정'}")
log(f"  {'─'*22} {'─'*12} {'─'*18} {'─'*12} {'─'*6}")

all_c_zeta = []; all_c_11a1 = []; all_c_delta = []
all_At_zeta = []; all_At_11a1 = []; all_At_delta = []

for t0, res in zeta_results.items():
    if res and res['c'] is not None:
        c_str = f"{res['c']:.6f}"; At_str = f"{res['A_true']:.4f}"
        log(f"  {'ζ(s)':<22} {t0:<12.6f} {c_str:<18} {At_str:<12}")
        all_c_zeta.append(res['c']); all_At_zeta.append(res['A_true'])

for t0, res in sym3_11a1_results.items():
    if res and res['c'] is not None:
        c_str = f"{res['c']:.6f}"; At_str = f"{res['A_true']:.4f}"
        log(f"  {'sym³(11a1)':<22} {t0:<12.6f} {c_str:<18} {At_str:<12} (#81)")
        all_c_11a1.append(res['c']); all_At_11a1.append(res['A_true'])

for t0, res in sym3d_results.items():
    if res and res['c'] is not None:
        c_str = f"{res['c']:.6f}"; At_str = f"{res['A_true']:.4f}"
        log(f"  {'sym³(Δ)':<22} {t0:<12.6f} {c_str:<18} {At_str:<12}")
        all_c_delta.append(res['c']); all_At_delta.append(res['A_true'])

log()
log("  ─── L-함수별 c 평균 ───")

def safe_mean(lst):
    if not lst: return None
    return sum(lst)/len(lst)
def safe_std(lst):
    if len(lst) < 2: return None
    m = safe_mean(lst)
    return (sum((x-m)**2 for x in lst)/len(lst))**0.5

c_zeta_mean   = safe_mean(all_c_zeta)
c_11a1_mean   = safe_mean(all_c_11a1)
c_delta_mean  = safe_mean(all_c_delta)
At_zeta_mean  = safe_mean(all_At_zeta)
At_11a1_mean  = safe_mean(all_At_11a1)
At_delta_mean = safe_mean(all_At_delta)

if c_zeta_mean is not None:
    log(f"  ζ(s)      [d=1]: c_mean={c_zeta_mean:.6f}, A_true_mean={At_zeta_mean:.4f}")
if c_11a1_mean is not None:
    log(f"  sym³(11a1)[d=4,w=1]: c_mean={c_11a1_mean:.6f}, A_true_mean={At_11a1_mean:.4f}")
if c_delta_mean is not None:
    log(f"  sym³(Δ)  [d=4,w=11]: c_mean={c_delta_mean:.6f}, A_true_mean={At_delta_mean:.4f}")

log()
log("  ─── 판정 ───")

# c 비교: 보편적(비슷) vs L-함수 의존적(다름)
judgment_lines = []

if c_11a1_mean is not None and c_delta_mean is not None:
    c_ratio = c_delta_mean / c_11a1_mean if abs(c_11a1_mean) > 1e-10 else float('inf')
    log(f"  c(sym³Δ)/c(sym³(11a1)) = {c_ratio:.4f}")
    if abs(c_ratio - 1.0) < 0.2:
        verdict = "✅ 보편적 (비율≈1): Hardy Z의 구조적 보정항"
        judgment_lines.append("Hardy Z 1/δ 보정항은 L-함수 독립 (보편 상수)")
    elif abs(c_ratio) > 5:
        verdict = "⚠️ L-함수 의존: c(Δ) >> c(11a1) → weight/conductor 의존"
        judgment_lines.append("1/δ 보정항은 weight 또는 conductor에 의존 (새 수학)")
    else:
        verdict = f"⚠️ 부분 의존 (비율={c_ratio:.2f}): 추가 확인 필요"
        judgment_lines.append(f"c 비율={c_ratio:.2f} — 불명확")
    log(f"  {verdict}")

log()
log("  ─── A_true 비교 (δ→0 외삽값) ───")
if At_11a1_mean is not None:
    log(f"  sym³(11a1) A_true ≈ {At_11a1_mean:.4f}")
    log(f"    (δ=0.001 관측값 71.8 → 외삽 {At_11a1_mean:.4f}으로 대폭 감소)")
if At_delta_mean is not None:
    log(f"  sym³(Δ)  A_true ≈ {At_delta_mean:.4f}")
    log(f"    (δ=0.001 관측값 286.26 → 외삽 {At_delta_mean:.4f})")
if At_11a1_mean is not None and At_delta_mean is not None:
    ratio_At = At_delta_mean / At_11a1_mean if abs(At_11a1_mean) > 1e-4 else float('inf')
    log(f"  A_true(Δ)/A_true(11a1) = {ratio_At:.4f}")
    if ratio_At > 5:
        log(f"  → weight 효과 유지: A_true(Δ) >> A_true(11a1)")
    elif abs(ratio_At - 1.0) < 0.5:
        log(f"  → weight 효과 사라짐: A_true 유사 (δ=0.001 차이는 보정항 효과)")
    else:
        log(f"  → 중간: A_true 비율={ratio_At:.2f}")

log()
log("  ─── #73 'c₁=0 정리'와의 정합성 ───")
log("  #73: ξ-bundle 기반 κ에서 1/δ 항 = 0 (함수방정식+Schwarz 반사)")
log("  #82: Hardy Z 기반 κ에서 1/δ 항 ≠ 0 (A ∝ 1/δ 관측)")
log("  → 두 κ의 정의 차이가 핵심:")
log("    ξ-bundle κ: 번들 곡률 (σ-방향 섭동)")
log("    Hardy Z κ: (Z'(t₀)/Z(t₀+δ))² (t-방향 섭동)")
log("  → Hardy Z의 1/δ 항은 Z의 테일러 전개에서 자연스럽게 나타남:")
log("    Z(t₀+δ) ≈ Z'(t₀)δ + Z''(t₀)δ²/2 + ...")
log("    κ = (Z'/Z(t₀+δ))² ≈ 1/δ² - Z''/Z'·(1/δ) + (Z''/Z')²/4 + ...")
log("    → c = -Z''(t₀)/Z'(t₀) · 1 (1차 계수)")
log("  → c값은 영점에서의 Z''/Z' 비율을 반영 (L-함수 의존적)")

log()
log("=" * 72)
log("성공 기준 총정리")
log("=" * 72)

total_points = sum(len(r['A_obs']) for r in sym3d_results.values() if r is not None)
log(f"  sym³(Δ) 측정점: {total_points}개 (목표 ≥12)")
log(f"  판정: {'✅ PASS' if total_points >= 12 else '⚠️ 부족'}")

richardson_ok = sum(1 for r in sym3d_results.values() if r and r['c'] is not None)
log(f"  Richardson 피팅 성공: {richardson_ok}/3 영점")

log()
for jl in judgment_lines:
    log(f"  ★ {jl}")

log()
log(f"완료 시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
flush_file()
print(f"\n결과 저장: {OUTFILE}")
