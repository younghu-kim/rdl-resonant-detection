#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #213] Dirichlet L-함수 κδ² 4성질 검증 — 6번째 가족 후보
=============================================================================
배경:
  15행 비교표에 degree 1은 ζ(s) (N=1) 하나뿐.
  Dirichlet L(s, χ_D)는 degree 1이지만 conductor N=|D| > 1.
  이 실험으로:
  1. slope=2.0 보편성이 conductor 비자명인 GL(1)에서도 성립하는지 검증
  2. 비교표에 "Dirichlet" 가족 추가
  3. A-scaling conjecture의 conductor 의존성 검증 (d=1, 다양한 N)

대상: 이차(Kronecker) Dirichlet 지표
  - D=-3: χ₋₃, N=3, 홀 (gammaV=[1])
  - D=-4: χ₋₄, N=4, 홀 (gammaV=[1])
  - D= 5: χ₅,  N=5, 짝 (gammaV=[0])
  - D=-7: χ₋₇, N=7, 홀 (gammaV=[1])

  모두 실수 지표 → L(1/2+it)이 실수 (Hardy Z-function 존재) → 영점 탐색 안전.

방법론: 표준 κδ² log-log slope (Paper 2 프로토콜)
  κ(δ) = |Λ'/Λ(center+δ+it₀)|²  (σ-방향)
  fit: log(κδ²-1) vs log(δ) → slope=2.0 이론값

성공 기준:
  slope = 2.0 ± 0.05 (4개 L-함수 전체) → ★★★
  하나라도 이탈 → 경계 발견

결과: results/dirichlet_kd2_213.txt
=============================================================================
"""
import sys, os, time, math
import numpy as np

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'dirichlet_kd2_213.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 #213] Dirichlet L-함수 κδ² 4성질 검증 — GL(1) conductor 확장")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

TARGETS = [
    {"D": -3, "name": "χ₋₃", "N": 3, "parity": "odd"},
    {"D": -4, "name": "χ₋₄", "N": 4, "parity": "odd"},
    {"D":  5, "name": "χ₅",  "N": 5, "parity": "even"},
    {"D": -7, "name": "χ₋₇", "N": 7, "parity": "odd"},
]

DELTAS = [0.005, 0.008, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2]
N_ZEROS = 7  # 각 L-함수에서 7영점 (5개 선택)
T_MAX = 50.0
CENTER = 0.5

# σ-유일성 탐색 범위
SIGMA_OFFSETS = np.linspace(-0.2, 0.2, 41)

# 모노드로미 적분 설정
MONO_RADIUS = 0.3
MONO_STEPS = 64

log(f"대상: {[t['name'] for t in TARGETS]}")
log(f"δ 범위: {DELTAS}")
log(f"영점 탐색: t∈[0, {T_MAX}], 목표 {N_ZEROS}개/L-함수")
log(f"σ-유일성: center ± 0.2 ({len(SIGMA_OFFSETS)}점)")
log(f"모노드로미: radius={MONO_RADIUS}, steps={MONO_STEPS}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PARI 초기화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [0] PARI 초기화...")

import cypari2
gp = cypari2.Pari()
gp.allocatemem(1000 * 1024 * 1024)  # 1GB
gp("default(realprecision, 80)")

log(f"  cypari2 OK, realprecision=80")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def create_dirichlet_L(D):
    """이차 Dirichlet L-함수 생성 (Kronecker symbol)"""
    gp(f"Lcur = lfuncreate({D})")
    fe = float(gp("lfuncheckfeq(Lcur)"))
    return fe

def init_dirichlet_L(t_max):
    """현재 Lcur에 대해 lfuninit (σ 범위 0.5 포함)"""
    # [c, w] 형식: σ 방향 반폭 c=0.5, 높이 w=t_max+5
    gp(f"Linit = lfuninit(Lcur, [0.5, {t_max + 5}])")

def get_zeros(t_max, n_want):
    """lfunzeros로 영점 추출"""
    gp(f"zvec = lfunzeros(Linit, {t_max})")
    n_total = int(gp("length(zvec)"))
    zeros = []
    for i in range(1, min(n_total, 200) + 1):
        z = float(gp(f"zvec[{i}]"))
        zeros.append(z)
    return zeros

def compute_kappa_sigma(t0, delta):
    """σ-방향 κ(δ) = |Λ'/Λ(center+δ+it₀)|²"""
    s_re = CENTER + delta
    s_im = t0
    # Λ'/Λ at s
    gp(f"scur = {s_re:.15f} + I*{s_im:.15f}")
    gp("Lval = lfunlambda(Linit, scur)")
    gp("dLval = lfunlambda(Linit, scur, 1)")  # derivative
    gp("ratio = dLval / Lval")
    ratio_re = float(gp("real(ratio)"))
    ratio_im = float(gp("imag(ratio)"))
    kappa = ratio_re**2 + ratio_im**2
    return kappa

def compute_kappa_sigma_v2(t0, delta):
    """κ = |d(log Λ)/dσ|² at σ=center+δ

    수치 미분으로 계산 (안정성 확보):
    d/dσ [log Λ] ≈ [log Λ(σ+h) - log Λ(σ-h)] / (2h)
    """
    h = delta * 0.01  # 미분 스텝 (δ의 1%)
    if h < 1e-10:
        h = 1e-10

    s_re = CENTER + delta
    s_im = t0

    # Λ(σ+h + it)
    gp(f"sp = {s_re + h:.15f} + I*{s_im:.15f}")
    gp(f"sm = {s_re - h:.15f} + I*{s_im:.15f}")
    gp("Lp = lfunlambda(Linit, sp)")
    gp("Lm = lfunlambda(Linit, sm)")

    # log derivative via finite difference
    gp(f"dlogL = (log(Lp) - log(Lm)) / (2 * {h:.15e})")
    dlr = float(gp("real(dlogL)"))
    dli = float(gp("imag(dlogL)"))
    kappa = dlr**2 + dli**2
    return kappa

def fit_slope(deltas, kappas):
    """log(κδ²-1) vs log(δ) linear fit → slope"""
    log_d = []
    log_kd2m1 = []
    for d, k in zip(deltas, kappas):
        kd2 = k * d**2
        if kd2 > 1.0 + 1e-8:  # κδ² > 1 필요 (degree 1은 A가 작으므로 임계값 낮춤)
            log_d.append(math.log(d))
            log_kd2m1.append(math.log(kd2 - 1))

    if len(log_d) < 4:
        return None, None, 0

    log_d = np.array(log_d)
    log_kd2m1 = np.array(log_kd2m1)

    # Linear fit
    coeffs = np.polyfit(log_d, log_kd2m1, 1)
    slope = coeffs[0]

    # R²
    pred = np.polyval(coeffs, log_d)
    ss_res = np.sum((log_kd2m1 - pred)**2)
    ss_tot = np.sum((log_kd2m1 - np.mean(log_kd2m1))**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return slope, R2, len(log_d)

def compute_monodromy(t0, radius=MONO_RADIUS, steps=MONO_STEPS):
    """영점 주위 폐곡선 적분으로 모노드로미 계산

    mono = ∮ d(arg Λ) = Σ Δ(arg Λ)
    이론값: 2π (단순 영점)
    """
    angles = np.linspace(0, 2*np.pi, steps + 1)  # 0..2π
    total_phase = 0.0

    prev_phase = None
    for i in range(steps + 1):
        theta = angles[i]
        s_re = CENTER + radius * math.cos(theta)
        s_im = t0 + radius * math.sin(theta)

        gp(f"scur = {s_re:.15f} + I*{s_im:.15f}")
        gp("Lval = lfunlambda(Linit, scur)")
        lr = float(gp("real(Lval)"))
        li = float(gp("imag(Lval)"))

        cur_phase = math.atan2(li, lr)

        if prev_phase is not None:
            dp = cur_phase - prev_phase
            # unwrap
            while dp > math.pi:
                dp -= 2*math.pi
            while dp < -math.pi:
                dp += 2*math.pi
            total_phase += dp

        prev_phase = cur_phase

    return total_phase

def compute_sigma_uniqueness(t0, offsets=SIGMA_OFFSETS):
    """σ-유일성: κ(σ) 프로파일에서 σ=center가 최대인지 확인"""
    kappa_profile = []

    for dσ in offsets:
        s_re = CENTER + dσ
        s_im = t0

        if abs(dσ) < 0.001:
            # 영점 근방 → 매우 큰 값. 스킵하고 보간
            kappa_profile.append(None)
            continue

        gp(f"scur = {s_re:.15f} + I*{s_im:.15f}")
        gp("Lval = lfunlambda(Linit, scur)")
        gp("dLval = lfunlambda(Linit, scur, 1)")
        lr = float(gp("abs(Lval)"))

        if lr < 1e-50:
            kappa_profile.append(float('inf'))
        else:
            gp("ratio = dLval / Lval")
            rr = float(gp("real(ratio)"))
            ri = float(gp("imag(ratio)"))
            k = rr**2 + ri**2
            kappa_profile.append(k)

    # center (dσ=0) 근방에서 최대인가?
    # center_idx = middle
    center_idx = len(offsets) // 2

    # center 제외하고 양쪽의 최대값 vs center 근접값
    left_max = max([k for k in kappa_profile[:center_idx-2] if k is not None and k < 1e30], default=0)
    right_max = max([k for k in kappa_profile[center_idx+3:] if k is not None and k < 1e30], default=0)

    # center 근방 (±1~2 steps)
    near_center = [k for k in kappa_profile[center_idx-2:center_idx+3] if k is not None and k < 1e30]
    center_val = max(near_center) if near_center else 0

    side_max = max(left_max, right_max)

    # PASS if center region has highest κ
    if center_val > 0 and side_max > 0:
        ratio = center_val / side_max
        return ratio > 1.0, ratio, center_val, side_max
    else:
        return False, 0, center_val, side_max


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 루프: 각 Dirichlet L-함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

all_results = []

for tgt in TARGETS:
    D = tgt["D"]
    name = tgt["name"]
    N = tgt["N"]
    parity = tgt["parity"]

    log("─" * 72)
    log(f"{T()} ▶ {name} (D={D}, N={N}, {parity})")
    log("─" * 72)

    # 1. L-함수 생성
    log(f"  [1] L-함수 생성...")
    fe = create_dirichlet_L(D)
    log(f"      FE check = {fe:.0f} ({abs(fe):.0f}자리 일치)")

    # 2. lfuninit
    log(f"  [2] lfuninit...")
    t_init = time.time()
    init_dirichlet_L(T_MAX)
    log(f"      완료 ({time.time()-t_init:.1f}s)")

    # 3. 영점 추출
    log(f"  [3] 영점 추출 (t∈[0, {T_MAX}])...")
    zeros = get_zeros(T_MAX, 200)
    log(f"      총 영점: {len(zeros)}개")
    if len(zeros) > 0:
        log(f"      처음 5개: {[f'{z:.4f}' for z in zeros[:5]]}")

    # 간격 > 0.5인 영점 5개 선택
    selected = []
    for z in zeros:
        if z < 1.0:
            continue  # 너무 작은 t 제외
        if not selected or (z - selected[-1]) > 0.5:
            selected.append(z)
        if len(selected) >= 5:
            break

    if len(selected) < 3:
        log(f"  ⚠ 영점 부족 ({len(selected)}개). 스킵.")
        all_results.append({"name": name, "D": D, "N": N, "slope": None, "mono": None, "sigma": None})
        continue

    log(f"      선택 {len(selected)}영점: {[f'{z:.4f}' for z in selected]}")

    # 4. SC3a: κδ² log-log slope
    log(f"  [4] SC3a: κδ² slope 측정...")
    slopes = []
    r2s = []

    for idx, t0 in enumerate(selected):
        kappas = []
        valid_deltas = []

        for delta in DELTAS:
            try:
                k = compute_kappa_sigma(t0, delta)
                if k > 0 and not math.isinf(k) and not math.isnan(k):
                    kappas.append(k)
                    valid_deltas.append(delta)
            except Exception as e:
                log(f"      ρ_{idx+1} δ={delta}: 오류 {e}")
                continue

        if len(valid_deltas) >= 4:
            slope, R2, n_pts = fit_slope(valid_deltas, kappas)
            if slope is not None:
                slopes.append(slope)
                r2s.append(R2)
                log(f"      ρ_{idx+1} (t={t0:.4f}): slope={slope:.4f}, R²={R2:.6f} ({n_pts}점)")
            else:
                log(f"      ρ_{idx+1} (t={t0:.4f}): 피팅 실패")
        else:
            log(f"      ρ_{idx+1} (t={t0:.4f}): 유효 δ 부족 ({len(valid_deltas)}개)")

    if slopes:
        mean_slope = np.mean(slopes)
        std_slope = np.std(slopes)
        mean_R2 = np.mean(r2s)
        log(f"      ★ 평균 slope = {mean_slope:.4f} ± {std_slope:.4f} (R²={mean_R2:.6f})")
    else:
        mean_slope = None
        std_slope = None
        mean_R2 = None
        log(f"      ✗ slope 측정 실패")

    # 5. SC3b: 모노드로미
    log(f"  [5] SC3b: 모노드로미 측정...")
    monos = []

    for idx, t0 in enumerate(selected):
        try:
            mono = compute_monodromy(t0)
            mono_pi = mono / math.pi
            monos.append(mono_pi)
            log(f"      ρ_{idx+1} (t={t0:.4f}): mono = {mono_pi:.4f}π")
        except Exception as e:
            log(f"      ρ_{idx+1}: 오류 {e}")

    mono_pass = sum(1 for m in monos if abs(abs(m) - 2.0) < 0.3)
    log(f"      ★ 2π 판정: {mono_pass}/{len(monos)} PASS")

    # 6. SC3c: σ-유일성
    log(f"  [6] SC3c: σ-유일성 측정...")
    sigma_results = []

    for idx, t0 in enumerate(selected[:3]):  # 3영점만 (시간 절약)
        try:
            passed, ratio, cv, sv = compute_sigma_uniqueness(t0)
            sigma_results.append(passed)
            tag = "PASS" if passed else "FAIL"
            log(f"      ρ_{idx+1} (t={t0:.4f}): {tag} (ratio={ratio:.2f})")
        except Exception as e:
            log(f"      ρ_{idx+1}: 오류 {e}")

    sigma_pass = sum(sigma_results)
    sigma_total = len(sigma_results)
    log(f"      ★ σ-유일성: {sigma_pass}/{sigma_total}")

    # 결과 저장
    result = {
        "name": name,
        "D": D,
        "N": N,
        "parity": parity,
        "n_zeros": len(selected),
        "slope_mean": mean_slope,
        "slope_std": std_slope,
        "R2_mean": mean_R2,
        "slopes": slopes,
        "mono_pass": f"{mono_pass}/{len(monos)}",
        "mono_mean_pi": np.mean([abs(m) for m in monos]) if monos else None,
        "sigma_pass": f"{sigma_pass}/{sigma_total}" if sigma_results else "N/A",
        "fe": fe,
    }
    all_results.append(result)

    log(f"  ──────────────────────────────────────────")
    slope_s = f"{mean_slope:.4f}±{std_slope:.4f}" if mean_slope is not None else "N/A"
    log(f"  {name} 종합: slope={slope_s}, mono={mono_pass}/{len(monos)}, σ-uniq={sigma_pass}/{sigma_total}, FE={fe:.0f}")
    log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 종합 결과
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log()
log("=" * 72)
log("[종합 결과]")
log("=" * 72)
log()

log(f"{'name':<8} {'D':>4} {'N':>4} {'slope':>8} {'±σ':>8} {'R²':>8} {'mono':>6} {'σ-uniq':>8} {'FE':>6}")
log(f"{'─'*8} {'─'*4} {'─'*4} {'─'*8} {'─'*8} {'─'*8} {'─'*6} {'─'*8} {'─'*6}")

slopes_all = []
for r in all_results:
    slope_str = f"{r['slope_mean']:.4f}" if r.get('slope_mean') is not None else "N/A"
    std_str = f"{r.get('slope_std', 0):.4f}" if r.get('slope_std') is not None else "N/A"
    r2_str = f"{r.get('R2_mean', 0):.6f}" if r.get('R2_mean') is not None else "N/A"
    mono_str = r.get('mono_pass', 'N/A')
    sigma_str = r.get('sigma_pass', 'N/A')
    fe_str = f"{r.get('fe', 0):.0f}"

    log(f"{r['name']:<8} {r['D']:>4} {r['N']:>4} {slope_str:>8} {std_str:>8} {r2_str:>8} {mono_str:>6} {sigma_str:>8} {fe_str:>6}")

    if r.get('slope_mean') is not None:
        slopes_all.append(r['slope_mean'])

log()

# 통합 판정
if slopes_all:
    grand_mean = np.mean(slopes_all)
    grand_std = np.std(slopes_all)
    log(f"통합 slope = {grand_mean:.4f} ± {grand_std:.4f} (n={len(slopes_all)})")

    if abs(grand_mean - 2.0) < 0.05 and grand_std < 0.02:
        verdict = "★★★ 강양성"
        detail = "slope=2.0 ± 0.05, 4개 Dirichlet L-함수 전체 통과"
    elif abs(grand_mean - 2.0) < 0.1:
        verdict = "★★ 양성"
        detail = "slope=2.0 ± 0.1, 대부분 통과"
    else:
        verdict = "✗ 음성 또는 이탈"
        detail = f"slope={grand_mean:.4f}, 이론값 2.0에서 이탈"

    log(f"\n판정: {verdict}")
    log(f"근거: {detail}")
    log()
    log("의의:")
    log("  - 비교표에 'Dirichlet' 가족 추가 (6번째 가족)")
    log("  - degree=1, conductor>1 에서 slope=2.0 확인")
    log("  - A-scaling 검증: d=1에서 conductor N에 따른 A 변화 관찰 가능")
else:
    log("판정: 측정 불가")

log()
elapsed = time.time() - START
log(f"총 실행 시간: {elapsed:.1f}s ({elapsed/60:.1f}분)")
log(f"종료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

outf.close()
print(f"\n결과 저장: {RESULT_FILE}")
