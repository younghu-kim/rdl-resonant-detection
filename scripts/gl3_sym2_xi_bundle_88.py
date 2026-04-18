#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #88 — GL(3) sym²(11a1) + sym²(37a1) ξ-bundle A(t₀) 재측정
center 자동 감지 (k = L[4], center = k/2)
=============================================================================
목적:
  #87에서 center=2.0 오류로 실패한 sym²(37a1) ξ-bundle A(t₀) 재측정.
  sym²(11a1)도 동일 방법으로 재측정하여 conductor 독립성 비교.
  PARI Ldata[4]에서 k 자동 추출 → center = k/2.

핵심:
  - k = L[4] (PARI Ldata 4번째 성분)
  - center = k/2 (sym²(E): k=3 → center=1.5 예상)
  - κ = |Λ'(s)/Λ(s)|² at s = center + δ + i*t₀ (σ-방향)
  - A(t₀) = κ - 1/δ² (δ-독립이면 ξ-bundle 정리)

σ-유일성:
  s = σ + δ + i*t₀, δ=0.001 고정
  σ ∈ [center-0.2, center-0.1, center, center+0.1, center+0.2]
  → center(=1.5)에서 κ 최대 기대 (GRH: 영점은 Re(s)=center에만 존재)

이전 관련 실험:
  - #87: sym²(37a1) 4성질 4/4 통과. ξ-bundle center=2.0 오류.
  - #74: sym²(11a1) A=12.79 (mpmath, Hardy Z 방법, #83에서 기준값으로 사용)
  - #83: sym³ ξ-bundle (GL4 template) — 이번 스크립트 기반

결과: results/gl3_sym2_xi_bundle_88.txt
=============================================================================
"""

import sys
import os
import time
import math
import statistics

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "results", "gl3_sym2_xi_bundle_88.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []

def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))

def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

log("=" * 72)
log("결과 #88 — GL(3) sym²(11a1) + sym²(37a1) ξ-bundle A(t₀) 재측정")
log("center 자동 감지 (k = L[4], center = k/2)")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()
log("★ 핵심 수정: #87에서 k=4 하드코딩 오류 → L[4]로 k 자동 추출")
log("  sym²(E): k=3, center=1.5 (sym³에서 k=4, center=2.0 가져온 오류)")
log()
flush_file()

# ─── PARI 초기화 ──────────────────────────────────────────────────────────
import cypari2
gp = cypari2.Pari()
gp.allocatemem(2000 * 1024 * 1024)  # 2GB
gp("default(realprecision, 100)")
log("PARI 초기화: 2GB 메모리, realprecision=100")
log()
flush_file()

# ─── 파라미터 ─────────────────────────────────────────────────────────────
# ξ-bundle δ 값 (수학자 지시: 7개)
DELTAS = [0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1]

# σ-유일성 δ 고정값
SIGMA_DELTA = 0.001

# 이전 데이터 (비교용)
D1_A = 1.27    # ζ (d=1) #73
D2_A = 3.93    # GL(2) (d=2) #74
D3_REF = 12.79  # GL(3) sym²(11a1) #74 (mpmath, Hardy Z 방법)


# ─── 유틸리티 함수 ─────────────────────────────────────────────────────────
def extract_pari_float(gp, expr):
    """PARI 표현식을 Python float로 변환 (강건한 파싱)"""
    try:
        return float(gp(expr))
    except Exception as e:
        log(f"  ⚠️ PARI float 추출 실패: {expr} → {e}")
        return float('nan')


def get_zeros(gp, linit_name, n_max=3):
    """lfuninit된 L-함수에서 첫 n_max개 영점 추출
    주의: PARI 변수명에 _ 금지 — zvec 사용
    """
    try:
        gp(f"zvec = lfunzeros({linit_name}, 30)")
        n_total = int(gp("length(zvec)"))
        n_use = min(n_total, n_max)
        zeros = []
        for i in range(1, n_use + 1):
            z = float(gp(f"zvec[{i}]"))
            zeros.append(z)
        return zeros, n_total
    except Exception as e:
        log(f"  ❌ 영점 추출 실패: {e}")
        return [], 0


def compute_kappa(gp, linit_name, sigma, t0):
    """κ = |Λ'(s)/Λ(s)|² 계산, s = sigma + i*t0
    주의: PARI 변수명에 _ 금지 — scur 사용
    """
    try:
        gp(f"scur = {sigma:.12f} + I*{t0:.10f}")
        gp(f"L0val = lfunlambda({linit_name}, scur)")
        gp(f"Lpval = lfunlambda({linit_name}, scur, 1)")
        abs_L0 = float(gp("abs(L0val)"))
        abs_Lp = float(gp("abs(Lpval)"))
        if abs_L0 < 1e-200:
            return None, f"|Λ(s)| 너무 작음 ({abs_L0:.2e})"
        kappa = (abs_Lp / abs_L0) ** 2
        return kappa, None
    except Exception as e:
        return None, str(e)[:80]


def measure_xi_bundle(gp, linit_name, center, zeros, deltas, label):
    """
    ξ-bundle κ 측정: 각 t₀에서 δ별 κ, κδ², A 계산
    s = center + δ + i*t₀

    Returns: list of dicts {'t0', 'data': [(delta, kappa, kd2, A)], 'mean_A', 'cv_A'}
    """
    log(f"  {label}: center={center:.4f}")
    log(f"  δ 값: {deltas}")
    log()

    all_results = []

    for t_zero in zeros:
        log(f"  ── t₀ = {t_zero:.6f} ──")
        row_data = []

        for delta in deltas:
            sigma = center + delta
            kappa, err = compute_kappa(gp, linit_name, sigma, t_zero)
            if kappa is None:
                log(f"    δ={delta:.4f}: ❌ {err}")
                continue
            kd2 = kappa * delta**2
            A = kappa - 1.0 / delta**2
            row_data.append((delta, kappa, kd2, A))
            ok = "✅" if 0.99 <= kd2 <= 1.01 else "⚠️"
            log(f"    δ={delta:.4f}: κ={kappa:.4f}, A={A:.6f}, κδ²={kd2:.6f} {ok}")

        # 통계
        As = [x[3] for x in row_data]
        if As:
            mean_A = statistics.mean(As)
            cv_A = (statistics.stdev(As) / abs(mean_A) * 100) if len(As) > 1 and abs(mean_A) > 1e-10 else 0.0
        else:
            mean_A = float('nan')
            cv_A = float('nan')

        log(f"    → 유효 점수: {len(row_data)}/{len(deltas)}, mean(A) = {mean_A:.4f}, CV = {cv_A:.2f}%")
        all_results.append({'t0': t_zero, 'data': row_data, 'mean_A': mean_A, 'cv_A': cv_A})
        log()
        flush_file()

    return all_results


def measure_sigma_uniqueness(gp, linit_name, center, zeros, sigma_delta, label):
    """
    σ-유일성: s = σ + δ + i*t₀, σ ∈ [center-0.2, ..., center+0.2], δ=0.001 고정
    → center에서 κ 최대이면 PASS
    """
    sigmas = [round(center + offset, 4) for offset in [-0.2, -0.1, 0.0, 0.1, 0.2]]
    log(f"  {label} σ-유일성: σ = {sigmas}, δ={sigma_delta}")
    log()

    pass_count = 0
    for t_zero in zeros:
        log(f"  ── t₀ = {t_zero:.6f} ──")
        kappa_vals = {}

        for sigma in sigmas:
            s = sigma + sigma_delta  # 영점 바로 옆 (δ offset)
            kappa, err = compute_kappa(gp, linit_name, s, t_zero)
            if kappa is None:
                log(f"    σ={sigma:.2f}+δ: ❌ {err}")
                if err and "작음" in err:
                    # Λ(s)=0에 너무 가까움 → κ → ∞ 처리
                    kappa_vals[sigma] = float('inf')
                continue
            kappa_vals[sigma] = kappa
            log(f"    σ={sigma:.2f}+δ: κ={kappa:.6f}")

        if not kappa_vals:
            log(f"    → ⚠️ 데이터 없음")
            log()
            flush_file()
            continue

        # center에서 최대인지 확인
        center_kappa = kappa_vals.get(center, float('nan'))
        if math.isnan(center_kappa) and center not in kappa_vals:
            # center 데이터 없음
            log(f"    → ⚠️ center σ={center} 데이터 없음")
            log()
            flush_file()
            continue

        max_sigma = max(kappa_vals, key=lambda s: kappa_vals[s] if not math.isinf(kappa_vals[s]) else 1e300)
        max_kappa = kappa_vals[max_sigma]

        # center가 최대인지 (±0.05 허용)
        is_pass = abs(max_sigma - center) < 0.05

        if is_pass:
            log(f"    → ✅ PASS: center(σ={center:.1f})에서 최대, κ={center_kappa:.4f}")
            pass_count += 1
        else:
            log(f"    → ❌ FAIL: 최대가 σ={max_sigma:.2f}(κ={max_kappa:.4f}), center κ={center_kappa:.4f}")

        # κ 비율 출력
        if center_kappa not in (float('inf'), float('nan')) and max_kappa not in (float('inf'),):
            try:
                ratio = max_kappa / center_kappa
                log(f"       max/center 비율: {ratio:.4f}")
            except Exception:
                pass

        log()
        flush_file()

    return pass_count, len(zeros)


# ════════════════════════════════════════════════════════════════════════════
# [1] sym²(11a1) 초기화
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[1] sym²(11a1) 초기화 — lfunsympow(E, 2), N=121")
log("=" * 72)
t_start = time.time()

gp("E11 = ellinit([0,-1,1,-10,-20])")  # 11a1
gp("L11 = lfunsympow(E11, 2)")

# k 자동 추출 (Ldata[4] = weight k)
k11_raw = str(gp("L11[4]"))
k11 = int(round(float(k11_raw)))
center11 = k11 / 2.0
log(f"  L11[4] (k) = {k11_raw} → k={k11}, center={center11}")

# 검증
if k11 == 3 and abs(center11 - 1.5) < 1e-9:
    log(f"  ✅ center=1.5 확인 (예상 k=3)")
else:
    log(f"  ⚠️ k={k11} 예상과 다름! 계속 진행 (center={center11})")

# lfuninit — 충분히 넓은 범위 [0, 35] 사용
log(f"  lfuninit([0, 35]) 시작...")
gp("L11i = lfuninit(L11, [0, 35])")
log(f"  lfuninit 완료 ({time.time()-t_start:.1f}s)")

# 영점 추출
zeros11, n_zeros11 = get_zeros(gp, "L11i", n_max=3)
log(f"  영점 수 (t∈[0,30]): {n_zeros11}개")
log(f"  처음 3개: {[f'{z:.6f}' for z in zeros11]}")

if len(zeros11) < 3:
    log("  ⚠️ 영점 3개 미만 — lfuninit 범위 문제 가능성")

log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [2] sym²(11a1) ξ-bundle κ 측정
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[2] sym²(11a1) ξ-bundle κ — 3영점 × 7δ = 21점")
log("=" * 72)

results11 = measure_xi_bundle(gp, "L11i", center11, zeros11, DELTAS, "sym²(11a1)")

all_A_11 = [x[3] for r in results11 for x in r['data']]
all_kd2_11 = [x[2] for r in results11 for x in r['data']]
all_delta_kd2_11 = [(x[0], x[2]) for r in results11 for x in r['data']]

if all_A_11:
    mean_A_11 = statistics.mean(all_A_11)
    cv_A_11 = (statistics.stdev(all_A_11) / abs(mean_A_11) * 100) if len(all_A_11) > 1 and abs(mean_A_11) > 1e-10 else 0.0
    log(f"  sym²(11a1) 전체 요약: mean(A) = {mean_A_11:.4f}, CV = {cv_A_11:.2f}%")
    log(f"  점수: {len(all_A_11)}/21")
else:
    mean_A_11 = float('nan')
    cv_A_11 = float('nan')
    log("  ⚠️ sym²(11a1) 데이터 없음")

log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [3] sym²(11a1) σ-유일성
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[3] sym²(11a1) σ-유일성 검사")
log("=" * 72)

sigma_pass11, sigma_total11 = measure_sigma_uniqueness(
    gp, "L11i", center11, zeros11, SIGMA_DELTA, "sym²(11a1)"
)
log(f"  sym²(11a1) σ-유일성 결과: {sigma_pass11}/{sigma_total11}")
log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [4] sym²(37a1) 초기화
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[4] sym²(37a1) 초기화 — lfunsympow(E, 2), N=1369")
log("=" * 72)
t_start2 = time.time()

gp("E37 = ellinit([0,0,1,-1,0])")  # 37a1
gp("L37 = lfunsympow(E37, 2)")

# k 자동 추출
k37_raw = str(gp("L37[4]"))
k37 = int(round(float(k37_raw)))
center37 = k37 / 2.0
log(f"  L37[4] (k) = {k37_raw} → k={k37}, center={center37}")

# 검증
if k37 == 3 and abs(center37 - 1.5) < 1e-9:
    log(f"  ✅ center=1.5 확인 (예상 k=3)")
else:
    log(f"  ⚠️ k={k37} 예상과 다름! 계속 진행 (center={center37})")

# lfuninit
log(f"  lfuninit([0, 35]) 시작...")
gp("L37i = lfuninit(L37, [0, 35])")
log(f"  lfuninit 완료 ({time.time()-t_start2:.1f}s)")

# 영점 추출
zeros37, n_zeros37 = get_zeros(gp, "L37i", n_max=3)
log(f"  영점 수 (t∈[0,30]): {n_zeros37}개")
log(f"  처음 3개: {[f'{z:.6f}' for z in zeros37]}")

# #87 결과와 비교
log(f"  #87 영점 참조: [2.158694, 3.217147, 3.978335]")
if zeros37:
    diff = abs(zeros37[0] - 2.158694)
    log(f"  t₁ 비교: {zeros37[0]:.6f} vs 2.158694 (차이={diff:.6f})")
    if diff < 0.01:
        log(f"  ✅ 일치 (lfuninit 정상)")
    else:
        log(f"  ⚠️ 차이 큼 — 영점 계산 이슈 가능")

if len(zeros37) < 3:
    log("  ⚠️ 영점 3개 미만 — lfuninit 범위 문제. 대안: [0, 50] 시도")
    log("  lfuninit([0, 50]) 재시도...")
    gp("L37i = lfuninit(L37, [0, 50])")
    zeros37, n_zeros37 = get_zeros(gp, "L37i", n_max=3)
    log(f"  재시도 후 영점: {[f'{z:.6f}' for z in zeros37]}")

log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [5] sym²(37a1) ξ-bundle κ 측정
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[5] sym²(37a1) ξ-bundle κ — 3영점 × 7δ = 21점")
log("=" * 72)

results37 = measure_xi_bundle(gp, "L37i", center37, zeros37, DELTAS, "sym²(37a1)")

all_A_37 = [x[3] for r in results37 for x in r['data']]
all_kd2_37 = [x[2] for r in results37 for x in r['data']]
all_delta_kd2_37 = [(x[0], x[2]) for r in results37 for x in r['data']]

if all_A_37:
    mean_A_37 = statistics.mean(all_A_37)
    cv_A_37 = (statistics.stdev(all_A_37) / abs(mean_A_37) * 100) if len(all_A_37) > 1 and abs(mean_A_37) > 1e-10 else 0.0
    log(f"  sym²(37a1) 전체 요약: mean(A) = {mean_A_37:.4f}, CV = {cv_A_37:.2f}%")
    log(f"  점수: {len(all_A_37)}/21")
else:
    mean_A_37 = float('nan')
    cv_A_37 = float('nan')
    log("  ⚠️ sym²(37a1) 데이터 없음")

log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [6] sym²(37a1) σ-유일성
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[6] sym²(37a1) σ-유일성 검사")
log("=" * 72)

sigma_pass37, sigma_total37 = measure_sigma_uniqueness(
    gp, "L37i", center37, zeros37, SIGMA_DELTA, "sym²(37a1)"
)
log(f"  sym²(37a1) σ-유일성 결과: {sigma_pass37}/{sigma_total37}")
log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [7] 핵심 비교표 + conductor 독립성
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[7] 핵심 비교표 + conductor 독립성")
log("=" * 72)
log()

log(f"  d=1 (ζ)                       A_mean = {D1_A:.4f}  [#73]")
log(f"  d=2 (GL2)                     A_mean = {D2_A:.4f}  [#74]")
log(f"  d=3 sym²(11a1) N=121    (#88) A_mean = {mean_A_11:.4f}")
log(f"  d=3 sym²(37a1) N=1369   (#88) A_mean = {mean_A_37:.4f}")
log(f"  d=3 ref sym²(11a1)      (#74) A       = {D3_REF:.4f}  (mpmath/Hardy Z)")
log()

# conductor 독립성
if not (math.isnan(mean_A_11) or math.isnan(mean_A_37)):
    pair_mean = (abs(mean_A_11) + abs(mean_A_37)) / 2.0
    if pair_mean > 1e-10:
        diff_pct = abs(mean_A_37 - mean_A_11) / pair_mean * 100.0
        log(f"  conductor 비교:")
        log(f"    |A(37a1) - A(11a1)| / mean = {diff_pct:.1f}%")
        if diff_pct < 200:
            log(f"    → ✅ conductor-약의존 시사 (<200%)")
        else:
            log(f"    → ⚠️ conductor-강의존 (≥200%)")
    else:
        log(f"  ⚠️ A_mean이 0에 너무 가까움")
else:
    log(f"  ⚠️ conductor 비교 불가 (nan)")

log()

# degree 단조성 (기존 d=1,2 기준)
d1d2 = D2_A > D1_A
d2d3_11 = (not math.isnan(mean_A_11)) and mean_A_11 > D2_A
d2d3_37 = (not math.isnan(mean_A_37)) and mean_A_37 > D2_A
log(f"  degree 단조성 (d=1→2): {'✅' if d1d2 else '❌'} ({D1_A:.2f} → {D2_A:.2f})")
log(f"  d=2→3 sym²(11a1): {'✅' if d2d3_11 else '❌'} ({D2_A:.2f} → {mean_A_11:.4f})")
log(f"  d=2→3 sym²(37a1): {'✅' if d2d3_37 else '❌'} ({D2_A:.2f} → {mean_A_37:.4f})")
log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [8] δ-독립성 표 출력
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[8] δ-독립성 표 — A(t₀, δ)")
log("=" * 72)
log()

def print_delta_table(results, label, deltas):
    """δ별 A값 표 출력"""
    log(f"  {label} A(t₀, δ):")
    header = f"  {'t₀':>10s} | " + " | ".join(f"δ={d:.4f}" for d in deltas)
    log(header)
    log("  " + "-" * (len(header) - 2))
    for r in results:
        delta_to_A = {x[0]: x[3] for x in r['data']}
        row_vals = []
        for d in deltas:
            if d in delta_to_A:
                row_vals.append(f"{delta_to_A[d]:9.4f}")
            else:
                row_vals.append("      N/A")
        log(f"  t₀={r['t0']:.5f} | " + " | ".join(row_vals))
    log()

print_delta_table(results11, "sym²(11a1)", DELTAS)
print_delta_table(results37, "sym²(37a1)", DELTAS)
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [9] 성공 기준 평가
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[9] 성공 기준 평가 (수학자 지시 기준)")
log("=" * 72)
log()

# 기준 1: sym²(37a1) κδ² ∈ [0.99, 1.01] at δ≤0.01
small_delta_kd2_37 = [(d, kd2) for (d, kd2) in all_delta_kd2_37 if d <= 0.01]
c1_ok = len(small_delta_kd2_37) > 0 and all(0.99 <= kd2 <= 1.01 for _, kd2 in small_delta_kd2_37)
c1_n = len(small_delta_kd2_37)
log(f"  ① sym²(37a1) κδ² ∈ [0.99,1.01] (δ≤0.01): {'✅' if c1_ok else '❌'} (n={c1_n})")
if small_delta_kd2_37:
    for d, kd2 in small_delta_kd2_37[:6]:
        log(f"     δ={d:.4f}: κδ²={kd2:.6f} {'✅' if 0.99<=kd2<=1.01 else '❌'}")

# 기준 2: sym²(11a1) 동일
small_delta_kd2_11 = [(d, kd2) for (d, kd2) in all_delta_kd2_11 if d <= 0.01]
c2_ok = len(small_delta_kd2_11) > 0 and all(0.99 <= kd2 <= 1.01 for _, kd2 in small_delta_kd2_11)
c2_n = len(small_delta_kd2_11)
log(f"  ② sym²(11a1) κδ² ∈ [0.99,1.01] (δ≤0.01): {'✅' if c2_ok else '❌'} (n={c2_n})")
if small_delta_kd2_11:
    for d, kd2 in small_delta_kd2_11[:6]:
        log(f"     δ={d:.4f}: κδ²={kd2:.6f} {'✅' if 0.99<=kd2<=1.01 else '❌'}")

# 기준 3: CV(A) < 50%
c3_37 = (not math.isnan(cv_A_37)) and cv_A_37 < 50.0
c3_11 = (not math.isnan(cv_A_11)) and cv_A_11 < 50.0
log(f"  ③ CV(A) < 50%:")
log(f"     sym²(37a1): {'✅' if c3_37 else '❌'} CV={cv_A_37:.2f}%")
log(f"     sym²(11a1): {'✅' if c3_11 else '❌'} CV={cv_A_11:.2f}%")

# 기준 4: σ-유일성 3/3 PASS
c4_37 = (sigma_total37 > 0) and (sigma_pass37 == sigma_total37)
c4_11 = (sigma_total11 > 0) and (sigma_pass11 == sigma_total11)
log(f"  ④ σ-유일성 3/3:")
log(f"     sym²(37a1): {'✅' if c4_37 else '❌'} ({sigma_pass37}/{sigma_total37})")
log(f"     sym²(11a1): {'✅' if c4_11 else '❌'} ({sigma_pass11}/{sigma_total11})")

# 기준 5: conductor 독립성
if not (math.isnan(mean_A_11) or math.isnan(mean_A_37)):
    pair_mean = (abs(mean_A_11) + abs(mean_A_37)) / 2.0
    if pair_mean > 1e-10:
        diff_pct = abs(mean_A_37 - mean_A_11) / pair_mean * 100.0
        c5_ok = diff_pct < 200.0
    else:
        c5_ok = False
        diff_pct = float('nan')
    log(f"  ⑤ conductor 독립성 (<200%): {'✅' if c5_ok else '❌'} ({diff_pct:.1f}%)")
else:
    c5_ok = False
    log(f"  ⑤ conductor 독립성: ❌ 데이터 부족")

# 전체 집계
criteria = {
    "sym²(37a1) κδ²": c1_ok,
    "sym²(11a1) κδ²": c2_ok,
    "CV(A)<50% 양쪽": c3_37 and c3_11,
    "σ-유일성 양쪽": c4_37 and c4_11,
    "conductor 독립성": c5_ok,
}
n_pass = sum(criteria.values())
log()
log(f"  통과: {n_pass}/5")
log()

# 종합 판정
if n_pass >= 4:
    verdict = "★★★ 강양성"
elif n_pass >= 3:
    verdict = "★★ 양성"
elif n_pass >= 2:
    verdict = "★ 조건부 양성"
else:
    verdict = "✗ 음성"

log(f"  종합 판정: {verdict}")
log()
flush_file()

# ════════════════════════════════════════════════════════════════════════════
# [10] 최종 요약
# ════════════════════════════════════════════════════════════════════════════
log("=" * 72)
log("[10] 최종 요약")
log("=" * 72)
log()

log(f"  ★ k 자동 감지 결과:")
log(f"    sym²(11a1): L11[4]={k11}, center={center11}")
log(f"    sym²(37a1): L37[4]={k37}, center={center37}")
log()
log(f"  ★ ξ-bundle A(t₀):")
log(f"    sym²(11a1): mean={mean_A_11:.4f}, CV={cv_A_11:.2f}%, n={len(all_A_11)}/21")
log(f"    sym²(37a1): mean={mean_A_37:.4f}, CV={cv_A_37:.2f}%, n={len(all_A_37)}/21")
log()
log(f"  ★ σ-유일성:")
log(f"    sym²(11a1): {sigma_pass11}/{sigma_total11}")
log(f"    sym²(37a1): {sigma_pass37}/{sigma_total37}")
log()
log(f"  ★ 성공 기준: {n_pass}/5 → {verdict}")
log()
log(f"  ★ 이전 오류 수정: #87 center=2.0 → #88 center={center37:.1f}")
log()

# per-zero A 값 출력
log(f"  sym²(11a1) per-zero A(t₀):")
for r in results11:
    log(f"    t₀={r['t0']:.5f}: mean(A)={r['mean_A']:.4f}, CV={r['cv_A']:.2f}%")
log()
log(f"  sym²(37a1) per-zero A(t₀):")
for r in results37:
    log(f"    t₀={r['t0']:.5f}: mean(A)={r['mean_A']:.4f}, CV={r['cv_A']:.2f}%")
log()

log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()
flush_file()

print(f"\n✅ 결과 저장: {OUTFILE}")
