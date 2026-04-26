#!/usr/bin/env python3
"""
[사이클 #359] C-353 — κ_bg Analytical Full Formula (ψ-보정) 검증

목표:
  C-352에서 κ_bg 이론이 중립(CV=48%)이었던 원인 분석.
  단순 공식 κ_bg = (log(N/4π²) + γ)² + π²/4 는 ψ(s) 항 누락.

  Full formula:
    κ_gamma(σ, t) = |g(σ+it)|²
    g(s) = (1/2)·log(N/(4π²)) + ψ(s)

  여기서 ψ는 digamma 함수. EC의 경우 Λ(s) = (√N/2π)^s · Γ(s) · L(s)이므로
  Λ'/Λ(s) = log(√N/2π) + ψ(s) + L'/L(s)
           = (1/2)log(N/4π²) + ψ(s) + L'/L(s)

  배경점(영점과 먼 곳)에서 L'/L ≈ bounded → κ_gamma가 κ_bg의 주도항.

프로토콜:
  1. C-352와 동일 8곡선, σ=1.03, t∈[2,50], dt=0.1
  2. 각 t에서:
     (a) κ_full = |Λ'/Λ|² (PARI, C-352와 동일)
     (b) κ_gamma = |(1/2)log(N/4π²) + ψ(1.03+it)|² (mpmath)
  3. median(κ_gamma) vs median(κ_full) → ratio 비교
  4. C-352의 단순 공식과 비교: CV 개선 측정

성공 기준:
  - ψ-보정 후 ratio CV < 20% → ★★★★ 양성 (이론 정량적 양성)
  - ψ-보정 후 ratio 1.0±0.3 전곡선 → ★★★★★ 양성
  - CV 여전히 >30% → 중립 (L'/L 기여가 무시 불가)

결과: results/ec_kappa_bg_analytical_c353.txt
"""

import sys, os, time, math
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

import mpmath
mpmath.mp.dps = 50

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 곡선 정보 (C-352와 동일) + C-352 경험적 κ_median 값
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CURVES = [
    {"label": "11a1",    "N": 11,    "rank": 0, "kappa_emp_median": 7.4360},
    {"label": "43a1",    "N": 43,    "rank": 0, "kappa_emp_median": 10.7309},
    {"label": "37a1",    "N": 37,    "rank": 1, "kappa_emp_median": 10.5140},
    {"label": "79a1",    "N": 79,    "rank": 1, "kappa_emp_median": 12.1244},
    {"label": "389a1",   "N": 389,   "rank": 2, "kappa_emp_median": 20.9702},
    {"label": "571a1",   "N": 571,   "rank": 2, "kappa_emp_median": 20.9483},
    {"label": "5077a1",  "N": 5077,  "rank": 3, "kappa_emp_median": 34.3023},
    {"label": "11197a1", "N": 11197, "rank": 3, "kappa_emp_median": 39.0306},
]

SIGMA = 1.03
T_MIN = 2.0
T_MAX = 50.0
DT = 0.1
EULER_GAMMA = float(mpmath.euler)

# 결과 파일
RESULT_FILE = os.path.join(os.path.dirname(__file__), "..", "results", "ec_kappa_bg_analytical_c353.txt")

lines = []
def log(msg):
    print(msg, flush=True)
    lines.append(msg)

def save():
    os.makedirs(os.path.dirname(RESULT_FILE), exist_ok=True)
    with open(RESULT_FILE, "w") as f:
        f.write("\n".join(lines) + "\n")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def gamma_contribution(sigma, t, N):
    """
    감마 인자 기여 g(s) = (1/2)log(N/(4π²)) + ψ(s)
    EC: Λ(s) = (√N/(2π))^s · Γ(s) · L(s)
    d/ds log Λ = log(√N/(2π)) + ψ(s) + L'/L(s)
              = (1/2)log(N/(4π²)) + ψ(s) + L'/L(s)

    Returns |g(s)|² = κ_gamma(s)
    """
    s = mpmath.mpc(sigma, t)
    log_term = mpmath.log(N / (4 * mpmath.pi**2)) / 2
    psi_val = mpmath.digamma(s)
    g = log_term + psi_val
    return float(abs(g)**2)


def kappa_gamma_median(N, sigma=SIGMA, t_min=T_MIN, t_max=T_MAX, dt=DT):
    """t∈[t_min, t_max]에서 κ_gamma의 median 계산."""
    t_arr = np.arange(t_min, t_max + dt/2, dt)
    kg_arr = np.array([gamma_contribution(sigma, float(t), N) for t in t_arr])
    return float(np.median(kg_arr)), kg_arr, t_arr


def kappa_bg_simple(N):
    """C-352의 단순 공식 (t-독립)."""
    log_term = math.log(N / (4 * math.pi**2)) + EULER_GAMMA
    return log_term**2 + (math.pi**2) / 4


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 추가 분석: 배경 vs 영점 근방 분리
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_kappa_gamma_profile(N, sigma=SIGMA, t_points=[5, 10, 15, 20, 25, 30, 35, 40, 45]):
    """대표 t값에서 κ_gamma 프로파일 계산."""
    results = []
    for t in t_points:
        kg = gamma_contribution(sigma, t, N)
        # 점근 전개 비교: ψ(σ+it) ≈ log(t) for large t
        approx = (math.log(N/(4*math.pi**2))/2 + math.log(t))**2 + (math.pi/2)**2
        results.append((t, kg, approx))
    return results


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    log("=" * 70)
    log("[Project RDL] C-353 — κ_bg Analytical Full Formula (ψ-보정) 검증")
    log("=" * 70)
    log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"σ={SIGMA}, t∈[{T_MIN},{T_MAX}], dt={DT}")
    log(f"공식: κ_gamma = |(1/2)log(N/4π²) + ψ(σ+it)|²")
    log("")

    # ── 1. 각 곡선별 κ_gamma median 계산 ──
    log("=" * 70)
    log("[Part 1] κ_gamma median 계산 (8곡선)")
    log("=" * 70)

    comparison = []

    for c in CURVES:
        label = c["label"]
        N = c["N"]
        rank = c["rank"]
        k_emp = c["kappa_emp_median"]

        # ψ-보정 (full)
        kg_med, kg_arr, t_arr = kappa_gamma_median(N)

        # C-352 단순 공식
        k_simple = kappa_bg_simple(N)

        ratio_full = k_emp / kg_med
        ratio_simple = k_emp / k_simple

        log(f"\n  [{label}] N={N}, rank={rank}")
        log(f"    κ_empirical (C-352 median) = {k_emp:.4f}")
        log(f"    κ_gamma (ψ-full, median)   = {kg_med:.4f}  → ratio = {ratio_full:.4f}")
        log(f"    κ_simple (C-352 공식)       = {k_simple:.4f}  → ratio = {ratio_simple:.4f}")
        log(f"    개선: |1-ratio| = {abs(1-ratio_full):.3f} (full) vs {abs(1-ratio_simple):.3f} (simple)")

        comparison.append({
            "label": label, "N": N, "rank": rank,
            "k_emp": k_emp, "kg_med": kg_med, "k_simple": k_simple,
            "ratio_full": ratio_full, "ratio_simple": ratio_simple,
        })

    # ── 2. 비교표 ──
    log(f"\n\n{'=' * 70}")
    log("[Part 2] 비교표")
    log("=" * 70)

    log(f"\n  {'곡선':>10s}  {'N':>6s}  {'rank':>4s}  {'κ_emp':>10s}  {'κ_gamma':>10s}  {'κ_simple':>10s}  {'r_full':>8s}  {'r_simple':>8s}")
    log(f"  {'-'*80}")

    ratios_full = []
    ratios_simple = []

    for c in comparison:
        log(f"  {c['label']:>10s}  {c['N']:6d}  {c['rank']:4d}  {c['k_emp']:10.4f}  {c['kg_med']:10.4f}  {c['k_simple']:10.4f}  {c['ratio_full']:8.4f}  {c['ratio_simple']:8.4f}")
        ratios_full.append(c['ratio_full'])
        ratios_simple.append(c['ratio_simple'])

    rf = np.array(ratios_full)
    rs = np.array(ratios_simple)

    log(f"\n  ψ-full:   mean={rf.mean():.4f}, std={rf.std():.4f}, CV={rf.std()/rf.mean()*100:.1f}%")
    log(f"  simple:   mean={rs.mean():.4f}, std={rs.std():.4f}, CV={rs.std()/rs.mean()*100:.1f}%")
    log(f"  CV 개선: {rs.std()/rs.mean()*100:.1f}% → {rf.std()/rf.mean()*100:.1f}%")

    # ── 3. t-의존 프로파일 (대표 곡선) ──
    log(f"\n\n{'=' * 70}")
    log("[Part 3] t-의존 프로파일: κ_gamma(t) vs 점근 근사")
    log("=" * 70)
    log("  점근 근사: κ_asymp = [(1/2)log(Nt/(4π²))]² + (π/2)²")
    log("  (ψ(σ+it) ≈ log(t) for large t)")

    for c in [CURVES[0], CURVES[-1]]:  # 11a1 (smallest N) and 11197a1 (largest N)
        label = c["label"]
        N = c["N"]
        log(f"\n  [{label}] N={N}")
        profile = compute_kappa_gamma_profile(N)
        log(f"    {'t':>5s}  {'κ_gamma':>10s}  {'κ_asymp':>10s}  {'ratio':>8s}")
        for t, kg, ka in profile:
            log(f"    {t:5.0f}  {kg:10.4f}  {ka:10.4f}  {kg/ka:8.4f}")

    # ── 4. L'/L 잔차 추정 ──
    log(f"\n\n{'=' * 70}")
    log("[Part 4] L'/L 잔차 추정")
    log("=" * 70)
    log("  κ_full = |g + L'/L|² = |g|² + 2Re(g·conj(L'/L)) + |L'/L|²")
    log("  → κ_full/κ_gamma - 1 = 2Re(g·conj(L'/L))/|g|² + |L'/L|²/|g|²")
    log("  ratio ≡ κ_emp/κ_gamma:")
    log("    ratio > 1 → L'/L 양의 기여 (영점 근방에서 증폭)")
    log("    ratio < 1 → L'/L 음의 기여 (영점 간 소거)")
    log("    ratio ≈ 1 → L'/L 기여 무시 가능 (감마 인자 지배)")

    log(f"\n  곡선별 잔차:")
    for c in comparison:
        residual_pct = (c['ratio_full'] - 1) * 100
        log(f"    {c['label']:>10s} (N={c['N']:>5d}): ratio={c['ratio_full']:.4f}, 잔차={residual_pct:+.1f}%")

    # ── 5. Conductor 스케일링 법칙 ──
    log(f"\n\n{'=' * 70}")
    log("[Part 5] Conductor 스케일링 법칙")
    log("=" * 70)

    log_N = np.log(np.array([c["N"] for c in comparison]))
    ratios = np.array([c["ratio_full"] for c in comparison])

    # Linear fit: ratio = a + b/log(N)
    inv_logN = 1.0 / log_N
    A = np.vstack([np.ones_like(inv_logN), inv_logN]).T
    coeff, res, _, _ = np.linalg.lstsq(A, ratios, rcond=None)
    a_fit, b_fit = coeff

    log(f"  모델: ratio = a + b/log(N)")
    log(f"  a = {a_fit:.4f} (→ 1이면 N→∞에서 감마 지배)")
    log(f"  b = {b_fit:.4f} (유한 N 보정)")

    fitted = a_fit + b_fit / log_N
    ss_res = np.sum((ratios - fitted)**2)
    ss_tot = np.sum((ratios - ratios.mean())**2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    log(f"  R² = {R2:.4f}")

    log(f"\n  곡선별 피팅:")
    log(f"  {'곡선':>10s}  {'N':>6s}  {'ratio':>8s}  {'fitted':>8s}  {'residual':>8s}")
    for i, c in enumerate(comparison):
        log(f"  {c['label']:>10s}  {c['N']:6d}  {c['ratio_full']:8.4f}  {fitted[i]:8.4f}  {c['ratio_full']-fitted[i]:+8.4f}")

    # ── 6. 종합 판정 ──
    cv_full = rf.std() / rf.mean() * 100
    cv_simple = rs.std() / rs.mean() * 100

    log(f"\n\n{'=' * 70}")
    log("[종합 판정]")
    log("=" * 70)

    log(f"\n  [1] ψ-보정 효과")
    log(f"      CV: {cv_simple:.1f}% (simple) → {cv_full:.1f}% (ψ-full)")

    if cv_full < 10:
        log(f"      ★★★★★ — 양성: ψ-보정으로 CV<10%, 감마 인자가 κ_bg 지배")
    elif cv_full < 20:
        log(f"      ★★★★ — 양성: ψ-보정으로 CV<20%, 감마 인자 준지배")
    elif cv_full < 30:
        log(f"      ★★★ — 조건부 양성: CV<30%, L'/L 잔차 보정 필요")
    else:
        log(f"      ★★ — 중립: CV≥30%, 감마 인자 단독으로 불충분")

    log(f"\n  [2] Conductor 스케일링")
    log(f"      ratio ≈ {a_fit:.3f} + {b_fit:.3f}/log(N), R²={R2:.3f}")
    if abs(a_fit - 1.0) < 0.15 and R2 > 0.8:
        log(f"      → N→∞에서 ratio→{a_fit:.3f}≈1: 감마 인자 점근 지배 확인")

    log(f"\n  [3] 물리적 해석")
    log(f"      Λ'/Λ = g(s) + L'/L(s)")
    log(f"      g(s) = Gamma contribution = (1/2)log(N/4π²) + ψ(s)")
    log(f"      L'/L(s) = zero sum = Σ 1/(s-ρ)")
    log(f"      Large N: |g| >> |L'/L| → κ_bg ≈ |g|² (감마 지배)")
    log(f"      Small N: |g| ~ |L'/L| → 교차항 2Re(g·conj(L'/L)) 무시 불가")

    elapsed = time.time() - t_start
    log(f"\n총 소요: {elapsed:.1f}s ({elapsed/60:.1f}분)")
    log(f"종료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    save()
    log(f"\n결과 저장: {RESULT_FILE}")


if __name__ == "__main__":
    main()
