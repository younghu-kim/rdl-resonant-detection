#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #76 — 짝수(even) 디리클레 문자 κ_near 측정: μ vs q 분리
=============================================================================
목적:
  χ₅² (mod 5, even, μ=0, Legendre symbol) + ζ 재측정으로
  A(t₀)의 μ 기여 vs conductor q 기여를 완전 분리.

핵심 논리:
  - ζ(s):    μ=0, q=1  → A_Λ(ζ)    (기준)
  - χ₅²(s): μ=0, q=5  → A_Λ(χ₅²) (순수 q 효과: A_Λ(χ₅²) - A_Λ(ζ))
  - χ₅¹(s): μ=1, q=5  → A_Λ(χ₅¹)=3.14 [#75 기결과]
  - χ₃(s):  μ=1, q=3  → A_Λ(χ₃)=2.25   [#75 기결과]

  q 기여 = A(χ₅², μ=0) - A(ζ, μ=0)  [같은 μ, 다른 q]
  μ 기여 = A(χ₅¹, μ=1) - A(χ₅², μ=0) [같은 q, 다른 μ]

χ₅² 문자값 (Legendre symbol mod 5):
  QR mod 5 = {1, 4}, NQR = {2, 3}
  chi=[0, 1, -1, -1, 1]  (indices 0~4)
  χ₅²(-1) = χ₅²(4) = 1 → even ✓, a=0, μ=0

ζ 재측정:
  동일 Λ'/Λ 공식 사용 (chi=[1], q=1, a=0)
  → A(t₀) 비교가 동일 방법 기반으로 이뤄짐 (방법 통일)
  #73 A(t₀)=1.30은 ξ'/ξ 기반이므로 직접 비교 불가 → 재측정 필수

실험 설계:
  - t 범위: [14, 59] (ζ·χ₅¹ 동일 범위)
  - δ = [0.001, 0.005, 0.01, 0.05, 0.1] (5점)
  - 총 ~130점 (ζ 13점×5 + χ₅² ~15점×5)

성공 기준:
  - χ₅²(even) A(t₀) 추출 성공 (CV < 1%, spread < 0.5%)
  - μ vs q 분리 정량화
  - κ·δ² ∈ [0.99, 1.15]
  - A(t₀) δ-독립 (spread < 0.5%)

결과: results/dirichlet_kappa_76.txt
=============================================================================
"""

import sys, os, time
import numpy as np

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

DPS = 50
mpmath.mp.dps = DPS

from bundle_utils import (
    find_zeros_zeta, find_zeros_dirichlet, CHARACTERS
)

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "dirichlet_kappa_76.txt"
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

T_MIN, T_MAX = 14.0, 59.0   # ζ·χ₅¹ 동일 t 범위

DELTAS = [0.001, 0.005, 0.01, 0.05, 0.1]

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 문자 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# ζ(s): q=1, a=0, chi=[1] (자명 지표 mod 1)
CHAR_ZETA = {
    'chi': [1],
    'q': 1, 'a': 0,
    'label': 'ζ(s) (q=1, μ=0)',
    'mu_label': 0,
    'is_zeta': True,
}

# χ₅² (mod 5, Legendre symbol, even, μ=0)
# QR mod 5 = {1,4}, NQR = {2,3}
# chi(0)=0, chi(1)=1, chi(2)=-1, chi(3)=-1, chi(4)=1
# χ₅²(-1) = χ₅²(4) = 1 → even ✓
CHAR_CHI5_EVEN = {
    'chi': [0, 1, -1, -1, 1],
    'q': 5, 'a': 0,
    'label': 'χ₅² (mod5, even, μ=0)',
    'mu_label': 0,
    'is_zeta': False,
}

# #75 기결과 (동일 Λ'/Λ 방법)
PREV_RESULTS = {
    'chi3':  {'label': 'χ₃ (mod3, odd, μ=1)',  'q': 3, 'a': 1, 'A_d01': 2.2451, 'n_zeros': 20},
    'chi4':  {'label': 'χ₄ (mod4, odd, μ=1)',  'q': 4, 'a': 1, 'A_d01': 2.7688, 'n_zeros': 22},
    'chi5':  {'label': 'χ₅¹ (mod5, odd, μ=1)', 'q': 5, 'a': 1, 'A_d01': 3.1407, 'n_zeros': 24},
}

# 수학자 이론 예측 (μ=0, q=5):
# A_predict ≈ A(ζ) + log(5/π)/2 ≈ 1.30 + 0.23 ≈ 1.53
# (log(5/π)/2 ≈ log(1.592)/2 ≈ 0.232)
# 단, 1.30은 ξ'/ξ 기반이므로, 실측 A_Λ(ζ)와 비교해야 함
LOG_5_OVER_PI_HALF = float(mpmath.log(mpmath.mpf(5) / mpmath.pi) / 2)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# κ_near 측정 함수 (Λ'/Λ 기반 — #75와 동일 방법)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def curvature_delta_dirichlet(t, delta, char_info):
    """
    κ_near = |Λ'/Λ|² at s = 0.5 + delta + it  (Λ'/Λ 공식)

    Λ'/Λ = log(q/π)/2 + ψ((s+a)/2)/2 + L'/L(s)

    NOTE: chi=[1] (q=1) 시 L(s,χ₁) = ζ(s) — 자명 지표
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
        L_val   = mpmath.dirichlet(s, chi)
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
# 단일 L-함수 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def analyze_Lfunc(char_key, char_info):
    """
    하나의 L-함수에 대해 κ_near + A(t₀) 분석.
    ζ: find_zeros_zeta 사용 (정확), χ: find_zeros_dirichlet 사용.
    """
    label = char_info.get('label', char_key)
    q = char_info['q']
    a = char_info['a']
    is_zeta = char_info.get('is_zeta', False)

    log()
    log("=" * 72)
    log(f"  분석 대상: {label}  (q={q}, a={a}, μ={a})")
    log("=" * 72)

    # Step 1: 영점 탐색
    log(f"  [Step 1] 영점 탐색 (t ∈ [{T_MIN}, {T_MAX}])")
    t0 = time.time()

    if is_zeta:
        # ζ: mpmath.zetazero 사용 (정확)
        zeros = find_zeros_zeta(T_MIN, T_MAX)
    else:
        # 디리클레 χ: 부호 변화 스캔
        zeros = find_zeros_dirichlet(char_info, t_min=T_MIN, t_max=T_MAX, n_scan=2000)

    dt_z = time.time() - t0
    n_z = len(zeros)
    log(f"  발견: {n_z}개  ({dt_z:.1f}s)")

    if n_z == 0:
        log(f"  ⚠️ 영점 0개 — {label} 탐색 실패")
        return None

    for j, z in enumerate(zeros):
        log(f"    #{j+1:>2}: t = {float(z):.8f}")
    flush_file()

    if n_z < 5:
        log(f"  ⚠️ 영점 {n_z}개로 부족 (기대 10+). 결과 신뢰도 낮음.")

    # Step 2: δ별 κ_near 측정
    log()
    log(f"  [Step 2] δ별 κ_near 측정 ({n_z}영점 × {len(DELTAS)}δ = {n_z*len(DELTAS)}점)")

    all_data = []
    delta_stats = {}

    for delta in DELTAS:
        theory_1_over_d2 = 1.0 / delta**2
        log()
        log(f"  ── δ = {delta} (1/δ² = {theory_1_over_d2:.2f}) ──")

        kappas = []
        t_start = time.time()

        for j, z_t in enumerate(zeros):
            k = curvature_delta_dirichlet(float(z_t), delta, char_info)
            kappas.append(k)
            all_data.append((delta, float(z_t), k))

            if np.isfinite(k) and k > 0:
                prod = k * delta**2
                A_est = k - theory_1_over_d2
                log(f"    #{j+1:>2} t={float(z_t):>10.6f}: κ={k:>12.4f}  κ·δ²={prod:>7.5f}  A(t₀)={A_est:>8.4f}")
            else:
                log(f"    #{j+1:>2} t={float(z_t):>10.6f}: κ=FAIL")

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
        A_vals   = [k - theory_1_over_d2 for k in fin]
        A_median = float(np.median(A_vals))
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

    A_at_delta01 = None
    products_ok  = []

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

        log(f"  {delta:>6.3f} | {med:>14.4f} | {cv:>7.3f} | {prod:>8.5f} | {A_mean:>12.4f} | {verdict}")

    log()
    log(f"  A(t₀) at δ=0.01: {A_at_delta01:.4f}" if A_at_delta01 is not None else "  A(t₀) at δ=0.01: N/A")

    return {
        'label': label,
        'q': q,
        'a': a,
        'n_zeros': n_z,
        'A_delta01': A_at_delta01,
        'delta_stats': delta_stats,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("결과 #76 — 짝수(even) 디리클레 문자 κ_near: μ vs q 분리")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"DPS={DPS}")
log(f"t 범위: [{T_MIN}, {T_MAX}]")
log(f"δ 목록: {DELTAS}")
log()
log("핵심 목표:")
log("  ζ(μ=0, q=1)  vs  χ₅²(μ=0, q=5) → 순수 q 기여 분리")
log("  χ₅²(μ=0, q=5) vs χ₅¹(μ=1, q=5) → 순수 μ 기여 분리")
log()
log(f"χ₅² 문자값 검증:")
log(f"  chi=[0, 1, -1, -1, 1]  (Legendre symbol mod 5)")
log(f"  χ₅²(1)=1, χ₅²(2)=-1, χ₅²(3)=-1, χ₅²(4)=1, χ₅²(5)=0")
log(f"  χ₅²(-1) = χ₅²(4) = 1 → even ✓, a=0, μ=0")
log()
log(f"#75 기결과 (χ₅¹, χ₃, χ₄ — 동일 Λ'/Λ 방법):")
for k, v in PREV_RESULTS.items():
    log(f"  {v['label']}: A(t₀)|δ=0.01 = {v['A_d01']:.4f}  ({v['n_zeros']}영점)")
log()
log(f"이론 예측:")
log(f"  log(5/π)/2 = {LOG_5_OVER_PI_HALF:.4f}")
log(f"  A_predict(χ₅², μ=0, q=5) ≈ A_Λ(ζ) + log(5/π)/2")
log(f"  (단, A_Λ(ζ) ≠ A_ξ(ζ)=1.30 — Λ'/Λ vs ξ'/ξ 차이로 재측정 필수)")
log()

t_total = time.time()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 측정 대상 목록
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

TARGETS = [
    ('zeta',       CHAR_ZETA),
    ('chi5_even',  CHAR_CHI5_EVEN),
]

results = {}
for key, ci in TARGETS:
    res = analyze_Lfunc(key, ci)
    if res is not None:
        results[key] = res
    flush_file()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# A(t₀) δ-일관성 점검
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log()
log("=" * 72)
log("A(t₀) = κ - 1/δ² δ-일관성 (전 δ에서 동일해야 함)")
log("=" * 72)
log()
log(f"{'L-함수':>26} | {'δ=0.001':>10} | {'δ=0.005':>10} | {'δ=0.01':>10} | {'δ=0.05':>10} | {'δ=0.1':>10} | spread%")
log("-" * 105)

for key, _ in TARGETS:
    if key not in results:
        continue
    r = results[key]
    label = r['label']
    ds = r['delta_stats']

    A_vals_by_delta = []
    row = f"  {label:>24} |"
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
# κ·δ² 스케일링 요약
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("κ·δ² 스케일링 법칙 확인 (κ·δ² ∈ [0.99, 1.15])")
log("=" * 72)
log()
log(f"{'L-함수':>26} | {'δ=0.001':>10} | {'δ=0.005':>10} | {'δ=0.01':>10} | {'δ=0.05':>10} | {'δ=0.1':>10}")
log("-" * 95)

for key, _ in TARGETS:
    if key not in results:
        continue
    r = results[key]
    label = r['label']
    ds = r['delta_stats']

    row = f"  {label:>24} |"
    for delta in DELTAS:
        rd = ds.get(delta, {})
        if rd.get('n', 0) > 0:
            prod = rd['product']
            row += f" {prod:>10.5f} |"
        else:
            row += f" {'N/A':>10} |"
    log(row)

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심: μ vs q 분리 비교표 (동일 δ=0.01)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("μ vs q 분리 비교표 (동일 δ=0.01)")
log("=" * 72)
log()
log("수학자 지시: δ=0.01 비교: ζ(μ=0,q=1) vs χ₅²(μ=0,q=5) vs χ₅¹(μ=1,q=5) vs χ₃(μ=1,q=3)")
log()
log(f"{'L-함수':>28} | {'μ':>3} | {'q':>3} | {'A(t₀)|δ=0.01':>14} | {'n_zeros':>7}")
log("-" * 70)

# 4종 L-함수 A(t₀) 수집
A_zeta = results.get('zeta', {}).get('A_delta01')
A_chi5_even = results.get('chi5_even', {}).get('A_delta01')
A_chi5_odd  = PREV_RESULTS['chi5']['A_d01']   # #75
A_chi3      = PREV_RESULTS['chi3']['A_d01']   # #75

n_zeta    = results.get('zeta', {}).get('n_zeros', '?')
n_chi5ev  = results.get('chi5_even', {}).get('n_zeros', '?')

def fmt_A(v):
    return f"{v:.4f}" if (v is not None and np.isfinite(v)) else "N/A"

log(f"  {'ζ(s)':>26} | {0:>3} | {1:>3} | {fmt_A(A_zeta):>14} | {n_zeta:>7}  [#76 신규]")
log(f"  {'χ₅² (mod5, even, μ=0)':>26} | {0:>3} | {5:>3} | {fmt_A(A_chi5_even):>14} | {n_chi5ev:>7}  [#76 신규]")
log(f"  {'χ₅¹ (mod5, odd, μ=1)':>26} | {1:>3} | {5:>3} | {fmt_A(A_chi5_odd):>14} | {PREV_RESULTS['chi5']['n_zeros']:>7}  [#75]")
log(f"  {'χ₃ (mod3, odd, μ=1)':>26} | {1:>3} | {3:>3} | {fmt_A(A_chi3):>14} | {PREV_RESULTS['chi3']['n_zeros']:>7}  [#75]")

log()

# μ vs q 분리 계산
log("분리 계산:")
if A_zeta is not None and np.isfinite(A_zeta) and A_chi5_even is not None and np.isfinite(A_chi5_even):
    q_contrib = A_chi5_even - A_zeta
    log(f"  q 기여 (순수, μ=0 고정): A(χ₅², q=5) - A(ζ, q=1) = {A_chi5_even:.4f} - {A_zeta:.4f} = {q_contrib:+.4f}")
    log(f"    이론 예측: log(5/π)/2 = {LOG_5_OVER_PI_HALF:+.4f}")
    q_theory_match = abs(q_contrib - LOG_5_OVER_PI_HALF) < 0.3
    log(f"    이론 일치: {'✅' if q_theory_match else '⚠️ 편차 있음'} (|측정 - 이론| = {abs(q_contrib - LOG_5_OVER_PI_HALF):.4f})")

    mu_contrib = A_chi5_odd - A_chi5_even
    log(f"  μ 기여 (순수, q=5 고정): A(χ₅¹, μ=1) - A(χ₅², μ=0) = {A_chi5_odd:.4f} - {A_chi5_even:.4f} = {mu_contrib:+.4f}")

    log()
    log(f"  총 기여 (q + μ 합산): {q_contrib:+.4f} + {mu_contrib:+.4f} = {q_contrib + mu_contrib:+.4f}")
    log(f"  직접 측정 (μ·q 동시 변화): A(χ₅¹) - A(ζ) = {A_chi5_odd - A_zeta:+.4f}")
    log(f"  합산 일관성 확인: {q_contrib + mu_contrib:.4f} vs {A_chi5_odd - A_zeta:.4f}")

    log()
    # μ vs q 상대 기여
    total_diff = abs(A_chi5_odd - A_zeta)
    if total_diff > 0.01:
        q_frac = abs(q_contrib) / total_diff * 100
        mu_frac = abs(mu_contrib) / total_diff * 100
        log(f"  μ vs q 상대 기여:")
        log(f"    q 효과 비중: {q_frac:.1f}%")
        log(f"    μ 효과 비중: {mu_frac:.1f}%")
else:
    log("  ⚠️ A_zeta 또는 A_chi5_even 측정 실패 — 분리 불가")

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 최종 성공 기준 체크
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

total_time = time.time() - t_total

log("=" * 72)
log("성공 기준 체크")
log("=" * 72)
log()

# 1. χ₅²(even) A(t₀) 추출 성공
chi5ev_ok = (
    'chi5_even' in results and
    results['chi5_even'].get('A_delta01') is not None and
    results['chi5_even']['delta_stats'].get(0.01, {}).get('cv', 999) < 1.0
)
log(f"  {'✅' if chi5ev_ok else '❌'} χ₅²(even) A(t₀) 추출 성공 (CV < 1%)")

# 2. κ·δ² ∈ [0.99, 1.15]
scale_ok = all(
    0.99 <= results[k]['delta_stats'].get(0.01, {}).get('product', 0) <= 1.15
    for k in results
    if results[k]['delta_stats'].get(0.01, {}).get('n', 0) > 0
)
log(f"  {'✅' if scale_ok else '❌'} κ·δ² ∈ [0.99, 1.15]")

# 3. ζ A(t₀) 재측정 성공
zeta_ok = (
    'zeta' in results and
    results['zeta'].get('A_delta01') is not None
)
log(f"  {'✅' if zeta_ok else '❌'} ζ(s) A(t₀) 재측정 성공")

# 4. μ vs q 분리 정량화
sep_ok = (
    A_zeta is not None and np.isfinite(A_zeta) and
    A_chi5_even is not None and np.isfinite(A_chi5_even)
)
log(f"  {'✅' if sep_ok else '❌'} μ vs q 분리 정량화")

# 5. A(t₀) δ-독립 (spread < 0.5%)
spread_ok = True
for key, _ in TARGETS:
    if key not in results:
        spread_ok = False
        continue
    ds = results[key]['delta_stats']
    A_list = [ds[d]['A_mean'] for d in DELTAS if ds.get(d, {}).get('n', 0) > 0]
    if len(A_list) >= 2:
        sp = (max(A_list) - min(A_list)) / (abs(np.mean(A_list)) + 1e-10) * 100
        if sp > 0.5:
            spread_ok = False
            log(f"    ⚠️ {results[key]['label']} spread = {sp:.2f}% > 0.5%")
log(f"  {'✅' if spread_ok else '❌'} A(t₀) δ-독립 (spread < 0.5%)")

n_pass = sum([chi5ev_ok, scale_ok, zeta_ok, sep_ok, spread_ok])
log()
log(f"성공 기준: {n_pass}/5 충족")
log()
log(f"총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log()
log("[완료]")

flush_file()
print(f"\n결과 저장: {OUTFILE}", flush=True)
