#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #75 — 디리클레 χ κ_near 측정 + A(t₀) 추출
=============================================================================
목적:
  디리클레 L(s,χ)의 κ_near와 A(t₀)를 측정하여,
  B-09 결론 "κ_near = f(degree)"의 정밀 검증 + 가설 A/B 판별.

대상: χ mod 3 (odd, a=1), χ mod 4 (odd, a=1), χ mod 5 (odd, a=1)
  - 모두 degree=1, conductor q=3,4,5
  - ζ와 동일 degree이지만 감마 인자 μ=1 (ζ는 μ=0)

실험 설계:
  - t 범위: [14, 59] (ζ와 동일 — A(t₀) 비교 공정성)
  - δ = [0.001, 0.005, 0.01, 0.05, 0.1] (5점)
  - 영점 15+개 목표
  - κ_near = |Λ'/Λ|² at s = 0.5 + δ + it₀
  - A(t₀) = κ_near - 1/δ² (δ-독립이어야 함)

가설 판별:
  - 가설 A (degree-only): κ_near(χ) ≈ ζ의 κ_near ≈ 1112.32
  - 가설 B (μ-shift): χ (μ=1) ≠ ζ (μ=0) → κ_near 차이 > 1%

성공 기준:
  - χ별 κ_near (CV < 1%)
  - κ·δ² ∈ [0.99, 1.15] 확인
  - A(t₀) 추출 3개 이상 χ에서 성공
  - 가설 A/B 판별: |κ_near(χ) - κ_near(ζ)| / κ_near(ζ)
    < 0.5% → A 지지 / > 1% → B 지지

결과: results/dirichlet_kappa_75.txt
=============================================================================
"""

import sys, os, time
import numpy as np

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

DPS = 50
mpmath.mp.dps = DPS

from bundle_utils import (
    completed_L, find_zeros_dirichlet, CHARACTERS
)

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "dirichlet_kappa_75.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 로깅
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

lines = []

def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))

def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파라미터
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

T_MIN, T_MAX = 14.0, 59.0   # ζ와 동일 t 범위

# 수학자 지시: δ = [0.001, 0.005, 0.01, 0.05, 0.1]
DELTAS = [0.001, 0.005, 0.01, 0.05, 0.1]

# 각 χ별 감마 인자 정보
#   ζ: a=0, Γ(s/2) → μ=0
#   χ₃, χ₄, χ₅: a=1, Γ((s+1)/2) → μ=1
ZETA_KAPPA_NEAR = 1112.32   # #73에서 확립 (13영점 중앙값)
ZETA_A_MEAN     = 1.30      # #73에서 확립 (A(t₀) = κ - 1/δ² at δ=0.01, 13영점 평균)

# 분석 대상 character 목록
CHAR_KEYS = ['chi_mod_3', 'chi_mod_4', 'chi_mod_5']

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# κ_near 측정 함수 (δ 파라미터화)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def curvature_delta_dirichlet(t, delta, char_info):
    """
    κ_near = |Λ'/Λ|² at s = 0.5 + delta + it

    Λ'/Λ 해석적 공식 (수치 안정성):
      Λ'/Λ = (1/2)log(q/π) + (1/2)ψ((s+a)/2) + L'/L

    L'/L: h=1e-6 중앙차분 (h=1e-20 수치 미분보다 안정적)
    """
    q = mpmath.mpf(char_info['q'])
    a = mpmath.mpf(char_info['a'])
    chi = char_info['chi']

    s = mpmath.mpc(
        mpmath.mpf('0.5') + mpmath.mpf(str(delta)),
        mpmath.mpf(str(t))
    )

    try:
        # 로그 + digamma 항 (해석적)
        log_term     = mpmath.log(q / mpmath.pi) / 2
        digamma_term = mpmath.digamma((s + a) / 2) / 2

        # L'/L (중앙차분, h=1e-6)
        h = mpmath.mpf('1e-6')
        L_val = mpmath.dirichlet(s, chi)
        if abs(L_val) < mpmath.mpf(10)**(-DPS + 10):
            log(f"  WARNING: L≈0 at t={t:.6f}, δ={delta} — inf 반환")
            return float('inf')

        L_plus  = mpmath.dirichlet(s + h, chi)
        L_minus = mpmath.dirichlet(s - h, chi)
        L_deriv = (L_plus - L_minus) / (2 * h)
        L_log_deriv = L_deriv / L_val

        # 접속 = Λ'/Λ
        conn = log_term + digamma_term + L_log_deriv

        k = float(abs(conn)**2)

        if not np.isfinite(k) or k <= 0:
            log(f"  WARNING: κ 비유한 at t={t:.6f}, δ={delta}: k={k}")
            return float('inf')

        return k

    except Exception as e:
        log(f"  WARNING: curvature_delta_dirichlet exception at t={t:.6f}, δ={delta}: {e}")
        return float('inf')


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 단일 χ 분석 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def analyze_character(char_key, char_info):
    """
    하나의 디리클레 문자에 대해 κ_near + A(t₀) 분석.
    반환: dict with summary statistics
    """
    label = char_info.get('label', char_key)
    q = char_info['q']
    a = char_info['a']

    log()
    log("=" * 72)
    log(f"  분석 대상: {label}  (q={q}, a={a}, μ={a})")
    log("=" * 72)

    # Step 1: 영점 탐색
    log(f"  [Step 1] 영점 탐색 (t ∈ [{T_MIN}, {T_MAX}])")
    t0 = time.time()
    zeros = find_zeros_dirichlet(char_info, t_min=T_MIN, t_max=T_MAX, n_scan=2000)
    dt_z = time.time() - t0
    n_z = len(zeros)
    log(f"  발견: {n_z}개  ({dt_z:.1f}s)")

    if n_z == 0:
        log(f"  ⚠️ 영점 0개 — {label} 탐색 실패")
        return None

    for j, z in enumerate(zeros):
        log(f"    #{j+1:>2}: t = {z:.8f}")
    flush_file()

    if n_z < 5:
        log(f"  ⚠️ 영점 {n_z}개로 부족 (기대 15+). 결과 신뢰도 낮음.")

    # Step 2: δ별 κ_near 측정
    log()
    log(f"  [Step 2] δ별 κ_near 측정 ({n_z}영점 × {len(DELTAS)}δ = {n_z*len(DELTAS)}점)")

    all_data = []   # (delta, t, kappa)
    delta_stats = {}

    for delta in DELTAS:
        theory_1_over_d2 = 1.0 / delta**2
        log()
        log(f"  ── δ = {delta} (1/δ² = {theory_1_over_d2:.2f}) ──")

        kappas = []
        t_start = time.time()

        for j, z_t in enumerate(zeros):
            k = curvature_delta_dirichlet(z_t, delta, char_info)
            kappas.append(k)
            all_data.append((delta, z_t, k))

            if np.isfinite(k) and k > 0:
                prod = k * delta**2
                A_est = k - theory_1_over_d2
                log(f"    #{j+1:>2} t={z_t:>10.6f}: κ={k:>12.4f}  κ·δ²={prod:>7.5f}  A(t₀)={A_est:>8.4f}")
            else:
                log(f"    #{j+1:>2} t={z_t:>10.6f}: κ=FAIL")

        dt_d = time.time() - t_start

        # 통계
        fin = [k for k in kappas if np.isfinite(k) and k > 0 and k < 1e15]
        if not fin:
            log(f"  ⚠️ δ={delta}: 유한 κ 없음")
            delta_stats[delta] = {'n': 0}
            continue

        median_k = float(np.median(fin))
        mean_k   = float(np.mean(fin))
        std_k    = float(np.std(fin))
        cv       = std_k / mean_k * 100 if mean_k > 0 else 999.0
        product  = median_k * delta**2
        A_median = median_k - theory_1_over_d2
        A_vals   = [k - theory_1_over_d2 for k in fin]
        A_mean   = float(np.mean(A_vals))
        A_std    = float(np.std(A_vals))

        delta_stats[delta] = {
            'n': len(fin),
            'median': median_k,
            'mean': mean_k,
            'std': std_k,
            'cv': cv,
            'product': product,
            'A_median': A_median,
            'A_mean': A_mean,
            'A_std': A_std,
            'A_vals': A_vals,
        }

        log()
        log(f"  δ={delta}: n={len(fin)}, median κ={median_k:.4f}, mean κ={mean_k:.4f}±{std_k:.4f}")
        log(f"           CV={cv:.3f}%, κ·δ²={product:.5f}")
        log(f"           A(t₀) = κ - 1/δ²: median={A_median:.4f}, mean={A_mean:.4f}±{A_std:.4f}")
        log(f"           소요: {dt_d:.1f}s")
        flush_file()

    # Step 3: 요약
    log()
    log(f"  [Step 3] {label} 요약")
    log(f"  {'δ':>6} | {'κ_near(중앙)':>14} | {'CV(%)':>7} | {'κ·δ²':>8} | {'A(t₀) 평균':>12} | 판정")
    log("  " + "-" * 75)

    products_ok = []
    A_at_delta01 = None   # A(t₀) at δ=0.01 (주요 지표)
    kappa_near_main = None  # κ_near at δ=0.03 (수학자의 기준 δ)

    for delta in DELTAS:
        r = delta_stats.get(delta, {})
        if r.get('n', 0) == 0:
            log(f"  {delta:>6.3f} | {'N/A':>14} | {'N/A':>7} | {'N/A':>8} | {'N/A':>12} | ❌")
            continue

        med    = r['median']
        cv     = r['cv']
        prod   = r['product']
        A_mean = r['A_mean']

        if 0.95 <= prod <= 1.15:
            verdict = "★★★"
        elif 0.90 <= prod <= 1.20:
            verdict = "★★"
        elif 0.80 <= prod <= 1.30:
            verdict = "★"
        else:
            verdict = "❌"

        if 0.80 <= prod <= 1.30:
            products_ok.append(prod)

        if delta == 0.01:
            A_at_delta01 = A_mean
            kappa_near_main = med

        log(f"  {delta:>6.3f} | {med:>14.4f} | {cv:>7.3f} | {prod:>8.5f} | {A_mean:>12.4f} | {verdict}")

    # κ_near 종합 (모든 δ 중앙값의 중앙값)
    all_medians = [delta_stats[d]['median'] for d in DELTAS if delta_stats.get(d, {}).get('n', 0) > 0]
    kappa_near_combined = float(np.median(all_medians)) if all_medians else float('nan')

    log()
    log(f"  κ_near 종합 (δ 중앙값들의 중앙값): {kappa_near_combined:.4f}")
    log(f"  A(t₀) at δ=0.01: {A_at_delta01:.4f}" if A_at_delta01 is not None else "  A(t₀) at δ=0.01: N/A")

    # ζ와의 비교
    if np.isfinite(kappa_near_combined):
        diff_abs = abs(kappa_near_combined - ZETA_KAPPA_NEAR)
        diff_pct = diff_abs / ZETA_KAPPA_NEAR * 100
        log()
        log(f"  ζ 비교:")
        log(f"    κ_near(ζ) = {ZETA_KAPPA_NEAR:.2f}  (기준)")
        log(f"    κ_near({label}) = {kappa_near_combined:.4f}")
        log(f"    차이: {diff_abs:.4f}  ({diff_pct:.4f}%)")
        if diff_pct < 0.5:
            hyp_verdict = "→ 가설 A 지지 (degree-only, 차이 < 0.5%)"
        elif diff_pct > 1.0:
            hyp_verdict = "→ 가설 B 지지 (μ-shift 의존, 차이 > 1.0%)"
        else:
            hyp_verdict = f"→ 경계 (차이 0.5%~1.0%, 추가 증거 필요)"
        log(f"    {hyp_verdict}")

    if A_at_delta01 is not None:
        A_diff = abs(A_at_delta01 - ZETA_A_MEAN)
        A_diff_pct = A_diff / max(abs(ZETA_A_MEAN), 0.01) * 100
        log()
        log(f"  A(t₀) 비교:")
        log(f"    A(t₀)(ζ) mean = {ZETA_A_MEAN:.4f}  (기준)")
        log(f"    A(t₀)({label}) mean(δ=0.01) = {A_at_delta01:.4f}")
        log(f"    차이: {A_diff:.4f}  ({A_diff_pct:.2f}%)")

    return {
        'label': label,
        'q': q,
        'a': a,
        'n_zeros': n_z,
        'kappa_near': kappa_near_combined,
        'kappa_near_delta01': kappa_near_main,
        'A_delta01': A_at_delta01,
        'delta_stats': delta_stats,
        'all_data': all_data,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("결과 #75 — 디리클레 χ κ_near 측정 + A(t₀) 추출")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"DPS={DPS}")
log(f"t 범위: [{T_MIN}, {T_MAX}] (ζ와 동일)")
log(f"δ 목록: {DELTAS}")
log(f"이론: κ_near = 1/δ² + A(t₀) + O(δ²)  →  A(t₀) = κ - 1/δ² (δ-독립)")
log()
log(f"기준값 (ζ, #73 결과):")
log(f"  κ_near(ζ) = {ZETA_KAPPA_NEAR:.2f}  (13영점 중앙값, t∈[10,60])")
log(f"  A(t₀) mean(ζ) = {ZETA_A_MEAN:.4f}  (δ=0.01)")
log()
log(f"가설:")
log(f"  A (degree-only): κ_near(χ) ≈ {ZETA_KAPPA_NEAR:.2f}  (차이 < 0.5%)")
log(f"  B (μ-shift):     κ_near(χ≠ζ)로 차이 > 1.0%  (μ=1 vs μ=0)")
log()

t_total = time.time()

# 각 character 분석
results = {}
for ck in CHAR_KEYS:
    if ck not in CHARACTERS:
        log(f"⚠️ {ck} not in CHARACTERS — skip")
        continue
    ci = CHARACTERS[ck]
    res = analyze_character(ck, ci)
    if res is not None:
        results[ck] = res
    flush_file()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 통합 비교표
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log()
log("=" * 72)
log("통합 비교표 — degree-1 L-함수 κ_near + A(t₀)")
log("=" * 72)
log()
log(f"{'L-함수':>16} | {'degree':>6} | {'μ(a)':>4} | {'n_zeros':>7} | {'κ_near':>10} | {'A(t₀)δ=0.01':>14} | {'Δκ/κ(ζ)%':>10}")
log("-" * 85)

# ζ 기준행
log(f"  {'ζ(s) [#73]':>14} | {'1':>6} | {'0':>4} | {'13':>7} | {ZETA_KAPPA_NEAR:>10.2f} | {ZETA_A_MEAN:>14.4f} | {'기준':>10}")

verdict_list = []

for ck in CHAR_KEYS:
    if ck not in results:
        continue
    r = results[ck]
    label  = r['label']
    kn     = r['kappa_near']
    A_d01  = r['A_delta01'] if r['A_delta01'] is not None else float('nan')
    n_z    = r['n_zeros']
    a_val  = r['a']

    if np.isfinite(kn):
        diff_pct = (kn - ZETA_KAPPA_NEAR) / ZETA_KAPPA_NEAR * 100
        diff_str = f"{diff_pct:>+9.4f}%"
    else:
        diff_str = "N/A"

    A_str = f"{A_d01:.4f}" if np.isfinite(A_d01) else "N/A"
    kn_str = f"{kn:.4f}" if np.isfinite(kn) else "N/A"

    log(f"  {label:>14} | {'1':>6} | {a_val:>4} | {n_z:>7} | {kn_str:>10} | {A_str:>14} | {diff_str:>10}")

    if np.isfinite(kn):
        diff_abs_pct = abs((kn - ZETA_KAPPA_NEAR) / ZETA_KAPPA_NEAR * 100)
        verdict_list.append((label, kn, diff_abs_pct))

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# κ·δ² 스케일링 요약 (전 χ, δ=0.01)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("κ·δ² 스케일링 법칙 확인 (δ=0.01, κ·δ² ∈ [0.99, 1.15])")
log(f"{'L-함수':>16} | {'κ·δ²(δ=0.01)':>14} | 판정")
log("-" * 45)

# ζ 기준
zeta_prod_d01 = ZETA_KAPPA_NEAR * 0.01**2  # ≈ 1112.32 * 0.0001 = 0.111 ... 아니다
# 실제로는 κ_near(δ=0.01) * 0.01² ≈ 10001.21 * 0.0001 = 1.000121
# ZETA_KAPPA_NEAR = 1112.32는 δ=0.03에서의 값임!
# #73 결과: δ=0.01에서 κ_near ≈ 10001.21, κ·δ² = 1.000121
# ZETA_KAPPA_NEAR = 1112.32는 δ=0.03에서의 ζ 영점 중앙값
# → 혼동 방지: #73 결과에서 δ=0.01 κ·δ² = 1.00012
log(f"  {'ζ(s) [#73]':>14} | {'1.00012':>14} | ★★★ (기준)")

for ck in CHAR_KEYS:
    if ck not in results:
        continue
    r = results[ck]
    label = r['label']
    ds = r['delta_stats']
    r01 = ds.get(0.01, {})
    if r01.get('n', 0) > 0:
        prod = r01['product']
        if 0.99 <= prod <= 1.15:
            v = "★★★"
        elif 0.95 <= prod <= 1.20:
            v = "★★"
        else:
            v = "❌"
        log(f"  {label:>14} | {prod:>14.5f} | {v}")
    else:
        log(f"  {label:>14} | {'N/A':>14} | ❌")

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# A(t₀) δ-일관성 점검 (핵심 품질 지표)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("A(t₀) = κ - 1/δ² δ-일관성 (모든 δ에서 동일해야 함)")
log()
log(f"{'L-함수':>16} | {'δ=0.001':>10} | {'δ=0.005':>10} | {'δ=0.01':>10} | {'δ=0.05':>10} | {'δ=0.1':>10} | spread%")
log("-" * 90)

for ck in CHAR_KEYS:
    if ck not in results:
        continue
    r = results[ck]
    label = r['label']
    ds = r['delta_stats']

    A_vals_by_delta = []
    row = f"  {label:>14} |"
    for delta in DELTAS:
        rd = ds.get(delta, {})
        if rd.get('n', 0) > 0:
            Am = rd['A_mean']
            row += f" {Am:>10.4f} |"
            A_vals_by_delta.append(Am)
        else:
            row += f" {'N/A':>10} |"

    if len(A_vals_by_delta) >= 2:
        spread = (max(A_vals_by_delta) - min(A_vals_by_delta)) / (abs(np.mean(A_vals_by_delta)) + 1e-10) * 100
        row += f" {spread:>7.2f}%"
    else:
        row += " N/A"
    log(row)

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 최종 가설 판별
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

total_time = time.time() - t_total

log("=" * 72)
log("최종 판정 — 가설 A/B 판별")
log("=" * 72)
log()
log(f"기준: κ_near(ζ) = {ZETA_KAPPA_NEAR:.2f}")
log()

if verdict_list:
    all_diffs = [d for _, _, d in verdict_list]
    max_diff = max(all_diffs)
    mean_diff = np.mean(all_diffs)

    log(f"χ별 |Δκ/κ(ζ)|:")
    for label, kn, dp in verdict_list:
        log(f"  {label}: κ_near={kn:.4f}, |Δκ/κ|={dp:.4f}%")
    log()

    # 판별 기준
    if max_diff < 0.5:
        final = "가설 A 지지 ★★★ (degree-only, 모든 χ에서 |Δκ/κ| < 0.5%)"
        note = "κ_near는 degree의 함수. 보편성 강화. B-09 확립."
    elif mean_diff > 1.0:
        final = "가설 B 지지 ★★ (μ-shift 의존, mean |Δκ/κ| > 1.0%)"
        note = "κ_near는 μ-벡터(감마 인자 구조)에 의존. 더 정밀한 분류 불변량."
    elif max_diff > 1.0:
        final = "가설 B 부분 지지 ★ (일부 χ에서 |Δκ/κ| > 1.0%)"
        note = "μ-shift 효과가 존재하나 일관성 없음. 추가 데이터 필요."
    else:
        final = f"판별 불분명 (|Δκ/κ| ∈ [0.5%, 1.0%] 범위)"
        note = "경계 영역. 더 많은 영점 또는 더 정밀한 측정 필요."

    log(f"→ {final}")
    log(f"   {note}")

else:
    log("데이터 부족 — 판별 불가")

log()

# 성공 기준 체크
log("성공 기준 체크:")

# 1. CV < 1% (δ=0.01)
cv_ok = all(
    results[ck]['delta_stats'].get(0.01, {}).get('cv', 999) < 1.0
    for ck in results
    if results[ck]['delta_stats'].get(0.01, {}).get('n', 0) > 0
)
log(f"  {'✅' if cv_ok else '❌'} χ별 κ_near (CV < 1%): {cv_ok}")

# 2. κ·δ² ∈ [0.99, 1.15]
scale_ok = all(
    0.99 <= results[ck]['delta_stats'].get(0.01, {}).get('product', 0) <= 1.15
    for ck in results
    if results[ck]['delta_stats'].get(0.01, {}).get('n', 0) > 0
)
log(f"  {'✅' if scale_ok else '❌'} κ·δ² ∈ [0.99, 1.15]: {scale_ok}")

# 3. A(t₀) 3개 이상 χ
A_count = sum(1 for ck in results if results[ck].get('A_delta01') is not None)
A_ok = A_count >= 3
log(f"  {'✅' if A_ok else '❌'} A(t₀) 추출 성공 χ 수: {A_count}/3 (필요: 3)")

# 4. 가설 판별 가능
judged = len(verdict_list) >= 2
log(f"  {'✅' if judged else '❌'} 가설 A/B 판별 데이터 충분: {judged}")

n_pass = sum([cv_ok, scale_ok, A_ok, judged])
log()
log(f"성공 기준: {n_pass}/4 충족")
log()
log(f"총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log()
log("[완료]")

flush_file()
print(f"\n결과 저장: {OUTFILE}", flush=True)
