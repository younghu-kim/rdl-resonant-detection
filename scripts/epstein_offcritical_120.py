#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #120] Epstein ζ off-critical 영점 탐색 + κδ²/c₁ 측정
=============================================================================
목표: B-28 경계 탐사 — h>1 Epstein ζ에서 off-critical 영점을 찾고,
      κδ² + c₁ 측정으로 Theorem 4 검증 표본 확대 (N=4→N=7+).

배경:
  - Epstein ζ (h>1)은 FE 있으나 오일러 곱 없음
  - RH 실패가 알려져 있음 (Davenport-Heilbronn 논법 확장)
  - #114에서 on-critical 4성질 ★★★ 확인 완료
  - DH off-critical 4개 (#115-#117)로 Theorem 4 검증 (N=4, ★★)

방법:
  1단계: PARI lfunqf로 Epstein ζ의 Λ(s) 평가
  2단계: 2D 격자 탐색 (σ∈[0.05,0.95] × t∈[1,200]) → |Λ| 최솟값 후보
  3단계: scipy minimize로 정밀화 → |Λ| < 1e-10 확인
  4단계: 확인된 off-critical 영점에서 κδ², c₁, 모노드로미 측정

사용 인프라: system python3 + cypari2 + scipy (qrop_env 아님)

형식:
  [1,0;0,5]  (disc=-20, h=2, non-Euler) — x²+5y²
  [2,1;1,3]  (disc=-23, h=3, non-Euler) — 2x²+xy+3y²
  [1,0;0,6]  (disc=-24, h=2, non-Euler) — x²+6y²

성공 기준:
  SC1: ≥3개 off-critical 영점 발견 (|σ-0.5|>0.01)
  SC2: 발견 시 κδ² ≠ 1 (off-critical), c₁ ≠ 0 측정
  SC3: DH off-critical과 정량 비교

결과: results/epstein_offcritical_120.txt
=============================================================================
"""

import sys, os, time
import numpy as np
from datetime import datetime
from scipy.optimize import minimize

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 100)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 출력 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'epstein_offcritical_120.txt')

outf = open(RESULT_FILE, 'w')

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

log("=" * 72)
log("[실험 #120] Epstein ζ off-critical 영점 탐색 + κδ²/c₁ 측정")
log("=" * 72)
log(f"시작: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log(f"PARI realprecision = 100")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. Epstein 형식 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 1. Epstein 형식 정의 ━━━")

FORMS = [
    {"M": "[1,0;0,5]",  "name": "x²+5y²",     "disc": -20, "h": 2},
    {"M": "[2,1;1,3]",  "name": "2x²+xy+3y²", "disc": -23, "h": 3},
    {"M": "[1,0;0,6]",  "name": "x²+6y²",     "disc": -24, "h": 2},
]

for f in FORMS:
    log(f"  {f['name']}: disc={f['disc']}, h={f['h']}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. 핵심 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def make_lfun(M_str):
    """PARI lfunqf L-함수 객체 생성."""
    return pari(f"lfunqf({M_str})")


def Lambda_eval(lfun_obj, sigma, t):
    """Λ(σ+it) 평가. 복소수 반환."""
    s_str = f"{sigma} + {t}*I" if t >= 0 else f"{sigma} - {abs(t)}*I"
    val = pari.lfun(lfun_obj, pari(s_str), 1)  # flag=1 → Λ 함수
    return complex(val)


def abs_Lambda(lfun_obj, sigma, t):
    """|Λ(σ+it)|"""
    return abs(Lambda_eval(lfun_obj, sigma, t))


def Lambda_deriv(lfun_obj, sigma, t, h=1e-8):
    """Λ'(s) 수치 미분 (중앙차분, σ 방향)."""
    lp = Lambda_eval(lfun_obj, sigma + h, t)
    lm = Lambda_eval(lfun_obj, sigma - h, t)
    return (lp - lm) / (2 * h)


def Lambda_deriv2(lfun_obj, sigma, t, h=1e-6):
    """Λ''(s) 수치 미분 (중앙차분, σ 방향)."""
    lp = Lambda_eval(lfun_obj, sigma + h, t)
    lc = Lambda_eval(lfun_obj, sigma, t)
    lm = Lambda_eval(lfun_obj, sigma - h, t)
    return (lp - 2*lc + lm) / (h * h)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. 2D 격자 탐색
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 탐색 범위
SIGMA_RANGE = np.linspace(0.05, 0.95, 37)   # Δσ=0.025
T_MAX = 200
T_RANGE = np.arange(1.0, T_MAX + 0.5, 0.5)  # Δt=0.5
OFF_CRITICAL_THRESH = 0.01  # |σ-0.5| > 이 값이면 off-critical

# 정밀화 기준
REFINE_THRESH = 0.01    # |Λ| < 이 값인 후보만 정밀화
ZERO_THRESH = 1e-8      # 정밀화 후 |Λ| < 이 값이면 영점 확인


def grid_search(lfun_obj, name):
    """2D 격자 탐색으로 off-critical 영점 후보 찾기."""
    log(f"\n{'='*60}")
    log(f"[격자 탐색] {name}: σ∈[0.05,0.95]×t∈[1,{T_MAX}]")
    log(f"  격자: {len(SIGMA_RANGE)}σ × {len(T_RANGE)}t = {len(SIGMA_RANGE)*len(T_RANGE)} 점")
    log(f"{'='*60}")

    t0 = time.time()

    # 격자 평가
    grid_vals = np.zeros((len(SIGMA_RANGE), len(T_RANGE)))
    for i, sigma in enumerate(SIGMA_RANGE):
        for j, t in enumerate(T_RANGE):
            try:
                grid_vals[i, j] = abs_Lambda(lfun_obj, sigma, t)
            except Exception:
                grid_vals[i, j] = 1e30
        if (i+1) % 10 == 0:
            log(f"  σ 진행: {i+1}/{len(SIGMA_RANGE)} ({time.time()-t0:.1f}초)")

    elapsed_grid = time.time() - t0
    log(f"  격자 평가 완료: {elapsed_grid:.1f}초")

    # Step 2: 후보 추출 — |Λ| < REFINE_THRESH 이고 |σ-0.5|>0.01
    candidates = []
    for i, sigma in enumerate(SIGMA_RANGE):
        if abs(sigma - 0.5) < OFF_CRITICAL_THRESH:
            continue  # on-critical은 건너뜀
        for j, t in enumerate(T_RANGE):
            if grid_vals[i, j] < REFINE_THRESH:
                candidates.append((sigma, t, grid_vals[i, j]))

    log(f"  후보 (|Λ|<{REFINE_THRESH}, |σ-0.5|>{OFF_CRITICAL_THRESH}): {len(candidates)}개")

    # Step 3: 후보 클러스터링 (인접한 후보 병합)
    if not candidates:
        log("  ❌ off-critical 영점 후보 없음")
        return []

    # 가장 작은 |Λ| 기준 정렬
    candidates.sort(key=lambda x: x[2])

    # 가까운 후보 병합 (Δσ<0.1, Δt<2)
    clusters = []
    used = set()
    for idx, (s, t, v) in enumerate(candidates):
        if idx in used:
            continue
        cluster = [(s, t, v)]
        used.add(idx)
        for jdx, (s2, t2, v2) in enumerate(candidates):
            if jdx in used:
                continue
            if abs(s - s2) < 0.1 and abs(t - t2) < 2.0:
                cluster.append((s2, t2, v2))
                used.add(jdx)
        # 클러스터에서 |Λ| 최소점
        best = min(cluster, key=lambda x: x[2])
        clusters.append(best)

    log(f"  클러스터링 후: {len(clusters)}개 독립 후보")
    for c in clusters[:20]:  # 최대 20개 표시
        log(f"    σ={c[0]:.4f}, t={c[1]:.2f}, |Λ|={c[2]:.3e}")

    return clusters


def refine_zero(lfun_obj, sigma0, t0, name=""):
    """후보 영점을 scipy로 정밀화."""
    def objective(x):
        return abs_Lambda(lfun_obj, x[0], x[1])

    try:
        res = minimize(objective, [sigma0, t0], method='Nelder-Mead',
                       options={'xatol': 1e-12, 'fatol': 1e-15, 'maxiter': 2000})
        sigma_opt, t_opt = res.x
        val_opt = res.fun

        # 임계대 내부 확인
        if sigma_opt < 0.01 or sigma_opt > 0.99:
            return None
        if t_opt < 0.5:
            return None
        # 실제 영점인지 확인
        if val_opt > ZERO_THRESH:
            return None
        # off-critical인지 확인
        if abs(sigma_opt - 0.5) < OFF_CRITICAL_THRESH:
            return None

        return {'sigma': sigma_opt, 't': t_opt, 'abs_Lambda': val_opt}
    except Exception as e:
        log(f"    정밀화 실패 ({name}): {e}")
        return None


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. κδ² / c₁ / 모노드로미 측정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

DELTAS = [0.01, 0.02, 0.03, 0.05, 0.07, 0.1]
CONTOUR_RADIUS = [0.05, 0.15, 0.5]
CONTOUR_STEPS = 64


def measure_kappa_delta2(lfun_obj, sigma0, t0, delta):
    """κδ² = |Λ'(σ₀+δ+it₀)/Λ(σ₀+δ+it₀)|² · δ²"""
    s_sigma = sigma0 + delta
    lam_val = Lambda_eval(lfun_obj, s_sigma, t0)
    if abs(lam_val) < 1e-50:
        return None
    lam_deriv = Lambda_deriv(lfun_obj, s_sigma, t0)
    kappa = abs(lam_deriv / lam_val) ** 2
    return kappa * delta**2


def measure_c1_fit(lfun_obj, sigma0, t0, deltas=DELTAS):
    """κδ² = 1 + c₁δ + c₂δ² + ... 에서 c₁ 추출 (선형 회귀)."""
    x_vals = []
    y_vals = []
    for d in deltas:
        kd2 = measure_kappa_delta2(lfun_obj, sigma0, t0, d)
        if kd2 is not None:
            x_vals.append(d)
            y_vals.append(kd2 - 1.0)  # κδ² - 1 = c₁δ + c₂δ² + ...

    if len(x_vals) < 3:
        return None, None, []

    x = np.array(x_vals)
    y = np.array(y_vals)

    # 2차 회귀: y = c₁·x + c₂·x²
    A = np.column_stack([x, x**2])
    coeffs, residuals, rank, sv = np.linalg.lstsq(A, y, rcond=None)
    c1_fit = coeffs[0]
    c2_fit = coeffs[1]

    # R² 계산
    y_pred = A @ coeffs
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2) + 1e-30
    r2 = 1 - ss_res / ss_tot

    kd2_data = list(zip(x_vals, [y + 1.0 for y in y_vals]))
    return c1_fit, r2, kd2_data


def measure_c1_analytic(lfun_obj, sigma0, t0):
    """c₁ = Re(Λ''/Λ')(ρ) 해석적 공식 (수치 미분)."""
    ld = Lambda_deriv(lfun_obj, sigma0, t0, h=1e-6)
    ld2 = Lambda_deriv2(lfun_obj, sigma0, t0, h=1e-5)
    if abs(ld) < 1e-50:
        return None
    ratio = ld2 / ld
    return ratio.real


def measure_monodromy(lfun_obj, sigma0, t0, radius):
    """모노드로미: σ₀+it₀ 주위 반지름 r 원을 따른 arg(Λ) 누적."""
    phases = []
    for k in range(CONTOUR_STEPS + 1):
        theta = 2 * np.pi * k / CONTOUR_STEPS
        s_sigma = sigma0 + radius * np.cos(theta)
        s_t = t0 + radius * np.sin(theta)
        try:
            val = Lambda_eval(lfun_obj, s_sigma, s_t)
            phases.append(np.angle(val))
        except Exception:
            return None

    # 연속 위상 차분 합산
    total = 0.0
    for k in range(len(phases) - 1):
        diff = phases[k+1] - phases[k]
        # branch cut 보정
        while diff > np.pi:
            diff -= 2 * np.pi
        while diff < -np.pi:
            diff += 2 * np.pi
        total += diff

    return total / np.pi  # π 단위


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. 메인 루프
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

all_off_zeros = []  # 전체 결과 수집
t_start = time.time()

for form in FORMS:
    M_str = form["M"]
    name = form["name"]
    disc = form["disc"]
    h = form["h"]

    log(f"\n{'#'*72}")
    log(f"# {name} (disc={disc}, h={h})")
    log(f"{'#'*72}")

    # L-함수 생성
    lfun_obj = make_lfun(M_str)

    # FE 확인
    fe_check = float(pari(f"lfuncheckfeq(lfunqf({M_str}))"))
    log(f"  FE 검증: {fe_check:.0f} (< -30이면 PASS)")

    # on-critical 영점 참조
    on_zeros_pari = pari(f"lfunzeros(lfunqf({M_str}), {T_MAX})")
    on_zeros = [float(z) for z in on_zeros_pari if float(z) > 0.5]
    log(f"  on-critical 영점 (t∈[0,{T_MAX}]): {len(on_zeros)}개")

    # 2D 격자 탐색
    candidates = grid_search(lfun_obj, name)

    if not candidates:
        log(f"\n  [{name}] off-critical 영점: 0개 발견")
        continue

    # 정밀화
    log(f"\n  ━━━ 정밀화 (Nelder-Mead) ━━━")
    confirmed = []
    for ci, (s0, t0, v0) in enumerate(candidates[:30]):  # 최대 30개 정밀화
        result = refine_zero(lfun_obj, s0, t0, f"{name} 후보 {ci+1}")
        if result is not None:
            # 거울쌍 중복 제거 (σ > 0.5만 기록)
            sigma = result['sigma']
            t_val = result['t']
            if sigma < 0.5:
                sigma = 1.0 - sigma  # 거울쌍으로 전환

            # 기존 확인된 영점과 중복 확인
            is_dup = False
            for prev in confirmed:
                if abs(prev['sigma'] - sigma) < 0.001 and abs(prev['t'] - t_val) < 0.1:
                    is_dup = True
                    break

            if not is_dup:
                result['sigma'] = sigma
                result['t'] = t_val
                confirmed.append(result)
                log(f"    ✅ OFF#{len(confirmed)}: σ={sigma:.8f}, t={t_val:.6f}, "
                    f"|Λ|={result['abs_Lambda']:.3e}, |σ-0.5|={abs(sigma-0.5):.6f}")

    log(f"\n  [{name}] off-critical 영점: {len(confirmed)}개 확인")

    if not confirmed:
        continue

    # κδ², c₁, 모노드로미 측정
    log(f"\n  ━━━ κδ² / c₁ / 모노드로미 측정 ━━━")

    for zi, zero in enumerate(confirmed):
        sigma0 = zero['sigma']
        t0 = zero['t']
        tag = f"EPST-{name[:5]}-OFF#{zi+1}"

        log(f"\n  [{tag}] σ={sigma0:.8f}, t={t0:.6f}")
        log(f"    |σ-0.5| = {abs(sigma0-0.5):.6f}")

        # κδ² 다중 δ 측정
        log(f"    δ    κδ²")
        log(f"    {'─'*25}")
        for delta in DELTAS:
            kd2 = measure_kappa_delta2(lfun_obj, sigma0, t0, delta)
            if kd2 is not None:
                status = "✅" if abs(kd2 - 1.0) < 0.01 else "❌"
                log(f"    {delta:.3f}  {kd2:.6f}  {status}")

        # c₁ 피팅
        c1_fit, r2, kd2_data = measure_c1_fit(lfun_obj, sigma0, t0)
        c1_analytic = measure_c1_analytic(lfun_obj, sigma0, t0)

        if c1_fit is not None:
            log(f"    c₁ (fit)     = {c1_fit:.6f}")
            log(f"    c₁ (analytic)= {c1_analytic:.6f}" if c1_analytic is not None else "    c₁ (analytic)= N/A")
            if c1_analytic is not None and abs(c1_analytic) > 1e-10:
                rel_err = abs(c1_fit - c1_analytic) / abs(c1_analytic)
                log(f"    rel_err      = {rel_err:.6f}")
            log(f"    R²           = {r2:.6f}")

            # c₁·|σ-0.5| 계산
            dist = abs(sigma0 - 0.5)
            c1_product = c1_fit * dist if c1_fit is not None else None
            c1a_product = c1_analytic * dist if c1_analytic is not None else None
            log(f"    c₁·|σ-0.5| (fit)     = {c1_product:.6f}" if c1_product is not None else "")
            log(f"    c₁·|σ-0.5| (analytic)= {c1a_product:.6f}" if c1a_product is not None else "")

        # 모노드로미 (3반경)
        log(f"    모노드로미 (π 단위):")
        for r in CONTOUR_RADIUS:
            mono = measure_monodromy(lfun_obj, sigma0, t0, r)
            if mono is not None:
                expected = 2.0  # 단순 영점
                status = "✅" if abs(abs(mono) - expected) < 0.1 else "❌"
                log(f"      r={r:.2f}: {mono:.4f}π  {status}")

        # 결과 저장
        zero['c1_fit'] = c1_fit
        zero['c1_analytic'] = c1_analytic
        zero['form'] = name
        zero['disc'] = disc
        all_off_zeros.append(zero)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. DH off-critical과 비교 + 종합 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log(f"\n{'='*72}")
log("종합 판정")
log(f"{'='*72}")

# DH off-critical 참조 (기존 #115-#117 데이터)
DH_OFFCRIT = [
    {"sigma": 0.808517, "t": 85.699, "c1_fit": 3.76, "dist": 0.308517},
    {"sigma": 0.650830, "t": 114.163, "c1_fit": 6.91, "dist": 0.150830},
    {"sigma": 0.574356, "t": 166.479, "c1_fit": 13.49, "dist": 0.074356},
    {"sigma": 0.724258, "t": 176.702, "c1_fit": 5.60, "dist": 0.224258},
]

n_found = len(all_off_zeros)
log(f"\n  Epstein off-critical 영점: {n_found}개 발견")
log(f"  DH off-critical (참조): {len(DH_OFFCRIT)}개")
log(f"  합계: {n_found + len(DH_OFFCRIT)}개")

if n_found > 0:
    log(f"\n  c₁·|σ-0.5| 비교:")
    log(f"  {'L-함수':<20} {'σ':>8} {'t':>10} {'|σ-½|':>8} {'c₁':>10} {'c₁·|σ-½|':>10}")
    log(f"  {'─'*72}")

    # DH 데이터
    for d in DH_OFFCRIT:
        product = d['c1_fit'] * d['dist']
        log(f"  {'DH':<20} {d['sigma']:>8.6f} {d['t']:>10.3f} {d['dist']:>8.6f} "
            f"{d['c1_fit']:>10.4f} {product:>10.4f}")

    # Epstein 데이터
    for z in all_off_zeros:
        dist = abs(z['sigma'] - 0.5)
        c1 = z.get('c1_fit', None)
        product = c1 * dist if c1 is not None else None
        log(f"  {z['form']:<20} {z['sigma']:>8.6f} {z['t']:>10.3f} {dist:>8.6f} "
            f"{c1:>10.4f} {product:>10.4f}" if c1 is not None and product is not None
            else f"  {z['form']:<20} {z['sigma']:>8.6f} {z['t']:>10.3f} {dist:>8.6f} {'N/A':>10} {'N/A':>10}")

# SC 판정
log(f"\n  성공 기준 판정:")
sc1 = n_found >= 3
log(f"  SC1 (≥3 off-critical): {'✅ PASS' if sc1 else '❌ FAIL'} ({n_found}개)")

if n_found > 0:
    kd2_all_off = [abs(z.get('c1_fit', 0)) > 0.1 for z in all_off_zeros if z.get('c1_fit') is not None]
    sc2 = len(kd2_all_off) > 0 and all(kd2_all_off)
    log(f"  SC2 (κδ²≠1, c₁≠0): {'✅ PASS' if sc2 else '❌ FAIL'}")
    sc3 = True  # 비교표 제공됨
    log(f"  SC3 (DH 비교):      ✅ PASS (비교표 제공)")
else:
    log(f"  SC2: N/A (영점 없음)")
    log(f"  SC3: N/A (영점 없음)")

# 판정
if n_found >= 3:
    log(f"\n  ★★★ 양성 — Epstein off-critical {n_found}개 발견. Theorem 4 표본 확대.")
elif n_found >= 1:
    log(f"\n  ★★ 조건부 양성 — {n_found}개 발견. 추가 탐색 필요.")
else:
    log(f"\n  음성 — t∈[1,{T_MAX}]에서 off-critical 영점 미발견. 탐색 범위 확대 필요.")

elapsed_total = time.time() - t_start
log(f"\n  총 소요: {elapsed_total:.1f}초 ({elapsed_total/60:.1f}분)")

outf.close()
log(f"\n결과 저장: {RESULT_FILE}")
