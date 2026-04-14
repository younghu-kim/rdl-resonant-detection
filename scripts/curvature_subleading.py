#!/usr/bin/env python3
"""
결과 #32: κ 차선도 구조 분석 — 곡률의 산술적 내용
========================================================
목적: 사이클 46 수학자 지시
  - Δκ ≡ κ - 1/δ² 의 t-의존 구조 측정
  - Δκ ∝ log(t/2π) 검증 (R² > 0.9이면 "곡률이 영점 밀도 인코딩")
  - δ-의존성: Δκ ∝ 1/δ^α (α ≈ 1 예상)
  - 이웃 영점 간격과의 상관 (Pearson r + p-value)

파트 A: 20개 영점에서 해석적 κ 측정, Δκ = κ - 1/δ² 추출
파트 B: Δκ vs log(t/2π) 선형 회귀
파트 C: δ-의존성 (δ = 0.01, 0.02, 0.03, 0.05, 0.1)
파트 D: 이웃 영점 간격과의 상관

이론:
  ξ'/ξ = 1/(s-ρ) + R(s)  (ρ 최근접 영점)
  κ = |1/δ + R|² = 1/δ² + (2/δ)Re[R] + |R|²
  차선도 ≈ (2/δ)·Re[R(s)], R ~ log(t) 스케일링 예상
"""

import sys, os, time
import numpy as np
import mpmath
from scipy import stats as scipy_stats

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_PATH = os.path.join(BASE_DIR, '..', 'results', 'curvature_subleading.txt')
os.makedirs(os.path.join(BASE_DIR, '..', 'results'), exist_ok=True)

_log_buf = []
def log(msg=''):
    print(msg, flush=True)
    _log_buf.append(str(msg))


# ═══════════════════════════════════════════════════════════════════════════
# 해석적 log-derivative (high_t_calibration.py와 동일)
# ═══════════════════════════════════════════════════════════════════════════

def xi_log_deriv_analytic(s):
    """
    ξ'(s)/ξ(s) = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2 + ζ'(s)/ζ(s)
    ζ'(s)/ζ(s): 유한차분 h=1e-6 (안전 범위)
    """
    z = mpmath.zeta(s)
    if abs(z) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e8, 0)
    h = mpmath.mpf('1e-6')
    zeta_d = (mpmath.zeta(s + h) - mpmath.zeta(s - h)) / (2 * h)
    return (mpmath.mpf(1) / s
            + mpmath.mpf(1) / (s - 1)
            - mpmath.log(mpmath.pi) / 2
            + mpmath.digamma(s / 2) / 2
            + zeta_d / z)


def kappa_analytic(s):
    """κ = |ξ'/ξ|²"""
    return float(abs(xi_log_deriv_analytic(s))**2)


def get_dps(t):
    """t 값에 따른 dps 자동 조절"""
    if t > 30000:
        return 100
    elif t > 5000:
        return 80
    else:
        return max(50, int(30 + t / 200))


# ═══════════════════════════════════════════════════════════════════════════
# 파트 A: Δκ vs t 측정 (20개 영점)
# ═══════════════════════════════════════════════════════════════════════════

# 수학자 지정 영점 인덱스
ZERO_INDICES = [1, 2, 3, 5, 10, 20, 50, 100, 200, 500,
                649, 1000, 2000, 4520, 5000, 7000, 10142, 12000, 15000, 20000]

DELTA = 0.03  # 오프셋: σ = 0.5 + δ


def part_a_delta_kappa():
    """
    20개 영점에서 해석적 κ 측정, Δκ = κ - 1/δ² 추출
    """
    log("=" * 70)
    log("파트 A: Δκ vs t 측정 — 20개 영점, δ=0.03")
    log("=" * 70)
    log(f"  이론값: 1/δ² = 1/{DELTA}² = {1/DELTA**2:.4f}")
    log(f"  영점 인덱스: {ZERO_INDICES}")
    log(f"  dps 자동 조절: max(50, 30+t/200), t>5000→80, t>30000→100")
    log()
    log(f"  {'N':>6}  {'t':>12}  {'κ':>12}  {'1/δ²':>10}  {'Δκ':>10}  {'Δκ/κ%':>8}  dps  시간(s)")
    log("  " + "-" * 78)

    results = []
    t_start = time.time()

    for N in ZERO_INDICES:
        t0 = time.time()
        # dps는 zetazero 계산에도 충분히 높게 설정
        # 먼저 낮은 dps로 t 추정 후 조절
        try:
            mpmath.mp.dps = 50  # zetazero는 50으로 충분
            z = mpmath.zetazero(N)
            t_val = float(z.imag)
        except Exception as e:
            log(f"  {N:>6}  ERROR zetazero: {e}")
            continue

        dps = get_dps(t_val)
        mpmath.mp.dps = dps

        try:
            s = mpmath.mpf('0.5') + mpmath.mpf(str(DELTA)) + mpmath.mpc(0, mpmath.mpf(str(t_val)))
            k = kappa_analytic(s)
            theory = 1.0 / (DELTA ** 2)
            delta_k = k - theory
            ratio_pct = delta_k / k * 100

            elapsed = time.time() - t0
            log(f"  {N:>6}  {t_val:>12.4f}  {k:>12.2f}  {theory:>10.4f}  {delta_k:>10.4f}  {ratio_pct:>8.3f}  {dps:>3}  {elapsed:.1f}")

            results.append({
                'N': N, 't': t_val, 'kappa': k,
                'theory': theory, 'delta_kappa': delta_k,
                'dps': dps, 'elapsed': elapsed
            })
        except Exception as e:
            log(f"  {N:>6}  {t_val:>12.4f}  ERROR: {e}")

    log()
    log(f"  파트 A 완료: {len(results)}/20 영점, 총 {time.time()-t_start:.1f}초")

    if len(results) == 0:
        log("  ⚠️ 결과 없음 — 탐색 로직 점검 필요")

    neg_count = sum(1 for r in results if r['delta_kappa'] < 0)
    if neg_count > 0:
        log(f"  ⚠️ Δκ < 0 이 {neg_count}개 발생 — 이론적으로 드묾, 확인 필요")

    return results


# ═══════════════════════════════════════════════════════════════════════════
# 파트 B: Δκ vs log(t/2π) 선형 회귀
# ═══════════════════════════════════════════════════════════════════════════

def part_b_linear_fit(part_a_results):
    """
    Δκ = a · log(t/2π) + b 선형 회귀
    R² > 0.9 → "곡률이 영점 밀도를 인코딩" (양성)
    """
    log()
    log("=" * 70)
    log("파트 B: Δκ vs log(t/2π) 선형 회귀")
    log("=" * 70)
    log("  이론: Δκ ≈ (2/δ)·Re[R(s)], R ~ log(t) → Δκ = a·log(t/2π) + b")
    log()

    if len(part_a_results) < 4:
        log("  ⚠️ 데이터 부족 (< 4점) — 회귀 스킵")
        return None

    t_arr   = np.array([r['t'] for r in part_a_results])
    dk_arr  = np.array([r['delta_kappa'] for r in part_a_results])
    x_arr   = np.log(t_arr / (2 * np.pi))  # log(t/2π)

    # scipy 선형 회귀
    slope, intercept, r_value, p_value, std_err = scipy_stats.linregress(x_arr, dk_arr)
    r2 = r_value ** 2

    log(f"  {'log(t/2π)':>14}  {'Δκ':>10}  {'예측 Δκ':>10}  {'잔차':>10}")
    log("  " + "-" * 52)
    for i, r in enumerate(part_a_results):
        x  = x_arr[i]
        dk = dk_arr[i]
        pred = slope * x + intercept
        resid = dk - pred
        log(f"  {x:>14.4f}  {dk:>10.4f}  {pred:>10.4f}  {resid:>10.4f}  (N={r['N']}, t={r['t']:.2f})")

    log()
    log(f"  선형 회귀 결과:")
    log(f"    기울기 a  = {slope:.4f}  (표준오차 {std_err:.4f})")
    log(f"    절편 b    = {intercept:.4f}")
    log(f"    R²        = {r2:.4f}")
    log(f"    p-value   = {p_value:.2e}")
    log()

    # 판정
    if r2 > 0.9 and slope > 0:
        verdict = "★ 양성: Δκ ∝ log(t/2π) 확인 (R²>{:.2f}) → 곡률이 영점 밀도를 인코딩".format(r2)
        sign = "✅"
    elif r2 > 0.7 and slope > 0:
        verdict = "⚠️ 약한 양성: R²={:.3f} (0.7~0.9) — 경향 존재, 강도 미달".format(r2)
        sign = "⚠️"
    elif r2 < 0.5:
        verdict = "❌ 음성: R²={:.3f} < 0.5 — Δκ의 log(t) 의존성 불명확".format(r2)
        sign = "❌"
    else:
        verdict = "⚠️ 중간: R²={:.3f} (0.5~0.7), p={:.2e}".format(r2, p_value)
        sign = "⚠️"

    log(f"  판정: {sign} {verdict}")

    return {
        'slope': slope,
        'intercept': intercept,
        'r2': r2,
        'p_value': p_value,
        'std_err': std_err,
        'verdict': verdict,
        'positive': r2 > 0.9 and slope > 0,
    }


# ═══════════════════════════════════════════════════════════════════════════
# 파트 C: δ-의존성
# ═══════════════════════════════════════════════════════════════════════════

DELTA_LIST = [0.01, 0.02, 0.03, 0.05, 0.1]
PART_C_ZEROS = [1, 100, 649, 4520, 10142]  # 수학자 지정


def part_c_delta_dependence():
    """
    δ = 0.01, 0.02, 0.03, 0.05, 0.1 에서 5개 영점의 Δκ 측정
    Δκ(δ) ∝ 1/δ^α → α 추정 (이론: α ≈ 1 for leading correction)
    """
    log()
    log("=" * 70)
    log("파트 C: δ-의존성 — Δκ(δ) ∝ 1/δ^α")
    log("=" * 70)
    log(f"  δ 목록: {DELTA_LIST}")
    log(f"  영점: N = {PART_C_ZEROS}")
    log()

    # 영점 t값 미리 획득
    mpmath.mp.dps = 80  # 고 t 영점에 충분
    zero_t = {}
    for N in PART_C_ZEROS:
        try:
            z = mpmath.zetazero(N)
            zero_t[N] = float(z.imag)
            log(f"  zetazero({N}) = {zero_t[N]:.6f}")
        except Exception as e:
            log(f"  ERROR zetazero({N}): {e}")
            zero_t[N] = None

    log()

    # δ별 Δκ 표 출력
    header = f"  {'N':>6}  {'t':>12}"
    for d in DELTA_LIST:
        header += f"  {'δ='+str(d):>12}"
    log(header)
    log("  " + "-" * (18 + 14 * len(DELTA_LIST)))

    results_c = {}  # delta -> list of (N, t, delta_kappa)
    for d in DELTA_LIST:
        results_c[d] = []

    for N in PART_C_ZEROS:
        t_val = zero_t.get(N)
        if t_val is None:
            continue

        dps = get_dps(t_val)
        mpmath.mp.dps = dps

        row = f"  {N:>6}  {t_val:>12.4f}"
        for d in DELTA_LIST:
            try:
                s = mpmath.mpf('0.5') + mpmath.mpf(str(d)) + mpmath.mpc(0, mpmath.mpf(str(t_val)))
                k = kappa_analytic(s)
                theory = 1.0 / (d ** 2)
                dk = k - theory
                row += f"  {dk:>12.4f}"
                results_c[d].append({'N': N, 't': t_val, 'delta_kappa': dk, 'theory': theory})
            except Exception as e:
                row += f"  {'ERROR':>12}"
                log(f"    WARNING δ={d}, N={N}: {e}")
        log(row)

    log()

    # α 추정: log(Δκ) = -α·log(δ) + const → 선형 회귀
    log("  α 추정 (각 영점별 log-log 회귀: log|Δκ| vs log(1/δ)):")
    log(f"  이론 예측: α ≈ 1 (선도 보정 (2/δ)·Re[R])")
    log()

    alpha_list = []
    for N in PART_C_ZEROS:
        t_val = zero_t.get(N)
        if t_val is None:
            continue

        log_inv_delta = []
        log_dk = []
        for d in DELTA_LIST:
            dk_list = [r['delta_kappa'] for r in results_c[d] if r['N'] == N]
            if dk_list:
                dk = dk_list[0]
                if dk > 0:
                    log_inv_delta.append(np.log(1.0 / d))
                    log_dk.append(np.log(dk))

        if len(log_inv_delta) >= 3:
            try:
                sl, ic, rv, pv, se = scipy_stats.linregress(log_inv_delta, log_dk)
                alpha_list.append(sl)
                log(f"    N={N:>6} (t={t_val:.2f}): α={sl:.4f} ± {se:.4f}, R²={rv**2:.4f}")
            except Exception as e:
                log(f"    N={N:>6} WARNING 회귀 실패: {e}")
        else:
            log(f"    N={N:>6}: 유효 데이터 부족 ({len(log_inv_delta)}점)")

    if alpha_list:
        alpha_mean = float(np.mean(alpha_list))
        alpha_std  = float(np.std(alpha_list))
        log()
        log(f"  평균 α = {alpha_mean:.4f} ± {alpha_std:.4f}")
        if abs(alpha_mean - 1.0) < 0.2:
            log(f"  → ✅ α ≈ 1 확인 — 선도 보정 (2/δ)·Re[R] 구조와 일치")
        else:
            log(f"  → ⚠️ α = {alpha_mean:.3f} ≠ 1 — 이론 예측과 편차 존재")
    else:
        alpha_mean = float('nan')
        alpha_std  = float('nan')

    return results_c, alpha_mean, alpha_std


# ═══════════════════════════════════════════════════════════════════════════
# 파트 D: 이웃 영점 간격과의 상관
# ═══════════════════════════════════════════════════════════════════════════

def part_d_gap_correlation(part_a_results):
    """
    각 20개 영점에서 gap = γ_{n+1} - γ_n 계산
    Δκ vs gap: Pearson 상관계수 + p-value
    """
    log()
    log("=" * 70)
    log("파트 D: 이웃 영점 간격과의 상관 (Δκ vs gap)")
    log("=" * 70)
    log("  이론: 좁은 gap → 큰 Δκ (이웃 영점 기여 증가)")
    log("  방법: gap = zetazero(N+1) - zetazero(N) (순방향 간격)")
    log()

    if len(part_a_results) < 4:
        log("  ⚠️ 데이터 부족 — 스킵")
        return None

    log(f"  {'N':>6}  {'t':>12}  {'gap':>10}  {'Δκ':>10}")
    log("  " + "-" * 44)

    gap_list = []
    dk_list  = []

    for r in part_a_results:
        N     = r['N']
        t_val = r['t']
        dk    = r['delta_kappa']

        # gap = γ_{N+1} - γ_N
        try:
            mpmath.mp.dps = max(40, get_dps(t_val))
            z_next = mpmath.zetazero(N + 1)
            gap = float(z_next.imag) - t_val
            gap_list.append(gap)
            dk_list.append(dk)
            log(f"  {N:>6}  {t_val:>12.4f}  {gap:>10.4f}  {dk:>10.4f}")
        except Exception as e:
            log(f"  {N:>6}  {t_val:>12.4f}  ERROR gap: {e}")

    log()

    if len(gap_list) < 4:
        log("  ⚠️ 유효 gap 데이터 부족 — Pearson 상관 스킵")
        return None

    gap_arr = np.array(gap_list)
    dk_arr  = np.array(dk_list)

    # Pearson 상관
    r_val, p_val = scipy_stats.pearsonr(gap_arr, dk_arr)
    log(f"  Pearson r(gap, Δκ) = {r_val:.4f},  p = {p_val:.3e}")
    log()

    # 이론 예측: 반상관 (gap 작을수록 Δκ 클수록 r < 0)
    if r_val < -0.3 and p_val < 0.05:
        verdict_d = f"✅ 반상관 확인: r={r_val:.3f}, p={p_val:.3e} < 0.05 → 좁은 간격에서 Δκ 증가"
        sign = "✅"
    elif r_val < 0 and p_val < 0.05:
        verdict_d = f"⚠️ 약한 반상관: r={r_val:.3f}, p={p_val:.3e} < 0.05"
        sign = "⚠️"
    elif p_val >= 0.05:
        verdict_d = f"❌ 유의하지 않음: r={r_val:.3f}, p={p_val:.3e} ≥ 0.05"
        sign = "❌"
    else:
        verdict_d = f"⚠️ 정상관 (예상과 반대): r={r_val:.3f}, p={p_val:.3e}"
        sign = "⚠️"

    log(f"  판정: {sign} {verdict_d}")

    # unfolded 간격 상관도 확인 (선택)
    # 평균 간격으로 나눠 unfolded gap
    # N→∞에서 gap_mean ≈ 2π/log(t/2π)
    t_arr = np.array([r['t'] for r in part_a_results[:len(gap_list)]])
    mean_spacing = 2 * np.pi / np.log(t_arr / (2 * np.pi))
    unfolded_gap = gap_arr / mean_spacing

    r_uf, p_uf = scipy_stats.pearsonr(unfolded_gap, dk_arr)
    log()
    log(f"  Unfolded gap 기준: r(unfolded_gap, Δκ) = {r_uf:.4f}, p = {p_uf:.3e}")
    if r_uf < 0 and p_uf < 0.05:
        log("    → ✅ Unfolded 기준에서도 반상관 유의")
    else:
        log("    → (unfolded에서 유의성 변화 — raw gap이 더 적절할 수 있음)")

    return {
        'r_pearson': r_val,
        'p_pearson': p_val,
        'r_unfolded': r_uf,
        'p_unfolded': p_uf,
        'verdict': verdict_d,
        'positive': r_val < -0.3 and p_val < 0.05,
    }


# ═══════════════════════════════════════════════════════════════════════════
# 종합 판정
# ═══════════════════════════════════════════════════════════════════════════

def final_summary(part_a_results, res_b, alpha_mean, res_d):
    log()
    log("=" * 70)
    log("★ 종합 판정 — 결과 #32")
    log("=" * 70)
    log()

    # 파트 A 요약
    if part_a_results:
        t_arr  = np.array([r['t'] for r in part_a_results])
        dk_arr = np.array([r['delta_kappa'] for r in part_a_results])
        mono   = np.all(np.diff(dk_arr[np.argsort(t_arr)]) >= -10)  # 단조 여부 (허용 오차 10)
        log(f"  파트 A: {len(part_a_results)}개 영점 측정 완료")
        log(f"    Δκ 범위: [{dk_arr.min():.3f}, {dk_arr.max():.3f}]")
        log(f"    t와 단조 증가: {'예 ✅' if mono else '아니오 ⚠️'}")
        log()

    # 파트 B 요약
    if res_b:
        log(f"  파트 B: Δκ = {res_b['slope']:.4f}·log(t/2π) + {res_b['intercept']:.4f}")
        log(f"    R² = {res_b['r2']:.4f},  p = {res_b['p_value']:.2e}")
        log(f"    판정: {res_b['verdict']}")
    else:
        log("  파트 B: 데이터 부족으로 회귀 불가")
    log()

    # 파트 C 요약
    if not np.isnan(alpha_mean):
        log(f"  파트 C: δ-의존성 지수 α = {alpha_mean:.4f}")
        if abs(alpha_mean - 1.0) < 0.2:
            log("    → ✅ α ≈ 1 — 선도 보정 (2/δ)·Re[R] 구조 지지")
        else:
            log(f"    → ⚠️ α = {alpha_mean:.3f} — 이론값 1.0에서 편차")
    else:
        log("  파트 C: α 추정 불가")
    log()

    # 파트 D 요약
    if res_d:
        log(f"  파트 D: Pearson r(gap, Δκ) = {res_d['r_pearson']:.4f}, p = {res_d['p_pearson']:.3e}")
        log(f"    판정: {res_d['verdict']}")
    else:
        log("  파트 D: gap 상관 분석 불가")
    log()

    # 최종
    log("─" * 70)
    log("최종 결론:")

    b_pos = res_b is not None and res_b.get('positive', False)
    c_pos = not np.isnan(alpha_mean) and abs(alpha_mean - 1.0) < 0.2
    d_pos = res_d is not None and res_d.get('positive', False)

    if b_pos and c_pos and d_pos:
        conclusion = ("★★ 강한 양성: Δκ ∝ log(t/2π) (R²>0.9) + α≈1 + 이웃 간격 반상관 — "
                      "곡률이 단순 기하(1/δ²)를 넘어 산술 정보(영점 밀도)를 인코딩")
    elif b_pos and (c_pos or d_pos):
        conclusion = ("★ 양성: Δκ ∝ log(t/2π) 확인 + 부가 검증 부분 통과 — "
                      "곡률의 차선도항이 산술 구조 반영 (일부 조건)")
    elif b_pos:
        conclusion = ("✅ 부분 양성: Δκ ∝ log(t/2π) (R²>0.9) — "
                      "로그 스케일링 확인; δ/간격 구조 추가 검증 권장")
    elif res_b is not None and res_b['r2'] > 0.7:
        conclusion = ("⚠️ 약한 양성: 로그 경향 존재 (R²>0.7), 강도 미달 — "
                      "차선도 산술 구조 단정 불가")
    else:
        conclusion = ("❌ 음성: Δκ의 log(t) 의존성 불명확 — "
                      "Leading 곡률 κ≈1/δ²의 자명성 인정. "
                      "논문: 'Leading curvature is geometric (1/δ²); sub-leading corrections are O(log t)'")

    log(f"  {conclusion}")
    log()

    if res_b is not None:
        r2 = res_b['r2']
        if r2 > 0.9:
            log("  논문 서술 (양성): '곡률의 차선도항 Δκ = κ - 1/δ²는 log(t/2π)에 선형 비례")
            log("    (R²={:.3f}), 이는 ξ-다발이 리만 제타의 산술적 성질(영점 밀도)을".format(r2))
            log("    기하학적으로 반영함을 시사한다.'")
        else:
            log("  논문 서술 (음성): '선도 곡률 κ ≈ 1/δ²는 영점으로부터의 기하학적 거리 효과이며,")
            log(f"    차선도 보정 Δκ는 O(log t) 이하의 작은 산술적 기여를 포함한다.'")

    return conclusion


# ═══════════════════════════════════════════════════════════════════════════
# 메인
# ═══════════════════════════════════════════════════════════════════════════

def main():
    log("결과 #32: κ 차선도 구조 분석 — 곡률의 산술적 내용")
    log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()
    log("  Δκ ≡ κ - 1/δ² (δ=0.03 기준)")
    log("  이론: Δκ ≈ (2/δ)·Re[R(s)], R(s) = ξ'/ξ - 1/(s-ρ) ~ O(log t)")
    log()

    t_total = time.time()

    # ─── 파트 A ────────────────────────────────────────────────────────────
    part_a_results = part_a_delta_kappa()

    # ─── 파트 B ────────────────────────────────────────────────────────────
    res_b = part_b_linear_fit(part_a_results)

    # ─── 파트 C ────────────────────────────────────────────────────────────
    results_c, alpha_mean, alpha_std = part_c_delta_dependence()

    # ─── 파트 D ────────────────────────────────────────────────────────────
    res_d = part_d_gap_correlation(part_a_results)

    # ─── 종합 판정 ─────────────────────────────────────────────────────────
    conclusion = final_summary(part_a_results, res_b, alpha_mean, res_d)

    log()
    log(f"총 소요: {time.time() - t_total:.1f}초")
    log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # 결과 파일 저장
    with open(RESULTS_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(_log_buf) + '\n')
    print(f"\n결과 저장: {RESULTS_PATH}", flush=True)


if __name__ == '__main__':
    main()
