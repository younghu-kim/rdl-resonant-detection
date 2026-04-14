"""
결과 #36: 간격 비율(spacing ratio) 분포의 GUE 검증
=======================================================
r_n = min(s_n, s_{n+1}) / max(s_n, s_{n+1}),  s_n = gamma_{n+1} - gamma_n

Atas et al., PRL 2013: unfolding 불필요.
GUE 예측: P(r) = (81√3/4π) · (r+r²)² / (1+r+r²)⁴
GOE 예측: P(r) = (27/4) · (r+r²) / (1+r+r²)^{5/2}
Poisson: P(r) = 2/(1+r)²

KS 검정: GUE vs GOE vs Poisson
성공 기준:
  양성: KS p_GUE > 0.05 + p_Poisson < 0.01
  강한 양성: + p_GOE < 0.01

영점: mpmath.zetazero, N=100 ~ 5000 (총 4901개)
"""

import mpmath
import numpy as np
from scipy import stats
import time
import os

mpmath.mp.dps = 25  # 충분한 정밀도

OUTPUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "results", "spacing_ratio_gue.txt")

# ──────────────────────────────────────────────
# GUE / GOE / Poisson P(r) 이론 밀도함수
# ──────────────────────────────────────────────
def p_gue(r):
    """GUE spacing ratio density (Atas et al. 2013)"""
    C_GUE = 81 * np.sqrt(3) / (4 * np.pi)
    return C_GUE * (r + r**2)**2 / (1 + r + r**2)**4

def p_goe(r):
    """GOE spacing ratio density"""
    C_GOE = 27.0 / 4.0
    return C_GOE * (r + r**2) / (1 + r + r**2)**2.5

def p_poisson(r):
    """Poisson spacing ratio density"""
    return 2.0 / (1 + r)**2

def build_cdf(pdf_func, n_points=10000):
    """r ∈ [0,1] 에서 PDF를 수치 적분하여 CDF 반환"""
    r_grid = np.linspace(0, 1, n_points + 1)
    p_vals = pdf_func(r_grid)
    cdf_vals = np.cumsum(p_vals) * (1.0 / n_points)
    cdf_vals = cdf_vals / cdf_vals[-1]  # 정규화
    return r_grid, cdf_vals

def cdf_from_grid(r_grid, cdf_vals, x):
    """CDF 보간"""
    return np.interp(x, r_grid, cdf_vals)

def ks_test_against_theory(data, pdf_func, n_points=10000):
    """KS 통계량 + p-값 계산 (scipy 없이 직접)"""
    r_grid, cdf_vals = build_cdf(pdf_func, n_points)
    # scipy KS: ECDF vs 이론 CDF
    n = len(data)
    data_sorted = np.sort(data)
    ecdf = np.arange(1, n + 1) / n
    theory_cdf = cdf_from_grid(r_grid, cdf_vals, data_sorted)
    D = np.max(np.abs(ecdf - theory_cdf))
    # Kolmogorov 분포 p-값
    p_val = stats.kstest(data, lambda x: cdf_from_grid(r_grid, cdf_vals, x)).pvalue
    return D, p_val

# ──────────────────────────────────────────────
# 메인
# ──────────────────────────────────────────────
def main():
    t_start = time.time()

    lines = []
    def pr(s=""):
        print(s, flush=True)
        lines.append(s)

    pr("결과 #36: 간격 비율(spacing ratio) 분포의 GUE 검증")
    pr(f"실행 시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    pr(f"N 범위: 100 ~ 5000 (총 4901개 영점)")
    pr()

    # ──────────────────────────────────────────
    # 파트 A: 영점 수집
    # ──────────────────────────────────────────
    pr("=" * 70)
    pr("파트 A: Riemann 영점 수집 (N=100 ~ 5000)")
    pr("=" * 70)

    N_min = 100   # N≥100부터 사용 (수학자 지시)
    N_max = 5000
    N_total = N_max - N_min + 1  # 4901개

    t_zeros = []
    t_collect_start = time.time()

    for n in range(N_min, N_max + 1):
        if (n - N_min) % 1000 == 0:
            elapsed = time.time() - t_collect_start
            pr(f"  [{n - N_min}/{N_total}] n={n}, 경과={elapsed:.0f}초")
        try:
            z = mpmath.zetazero(n)
            t_zeros.append(float(z.imag))
        except Exception as e:
            pr(f"  WARNING: n={n} 실패 — {e}")

    t_zeros = np.array(t_zeros)
    pr(f"  수집 완료: {len(t_zeros)}개")
    pr(f"  t 범위: {t_zeros[0]:.2f} ~ {t_zeros[-1]:.2f}")
    pr(f"  소요 시간: {time.time() - t_collect_start:.1f}초")
    pr()

    # ──────────────────────────────────────────
    # 파트 B: 간격 s_n = gamma_{n+1} - gamma_n
    # ──────────────────────────────────────────
    pr("=" * 70)
    pr("파트 B: 간격 s_n 계산 + spacing ratio r_n")
    pr("=" * 70)

    spacings = np.diff(t_zeros)  # s_n = gamma_{n+1} - gamma_n
    N_s = len(spacings)

    # r_n = min(s_n, s_{n+1}) / max(s_n, s_{n+1}) for n=0..N_s-2
    r_vals = []
    for i in range(N_s - 1):
        s1 = spacings[i]
        s2 = spacings[i + 1]
        if s1 <= 0 or s2 <= 0:
            pr(f"  WARNING: 비정상 간격 i={i}, s1={s1}, s2={s2}")
            continue
        r = min(s1, s2) / max(s1, s2)
        r_vals.append(r)

    r_vals = np.array(r_vals)
    pr(f"  s_n 개수: {N_s}")
    pr(f"  r_n 개수: {len(r_vals)}")
    pr(f"  r 범위: [{r_vals.min():.4f}, {r_vals.max():.4f}]")
    pr(f"  r 평균: {r_vals.mean():.4f}  (GUE 예측: ~0.5307, Poisson: ~0.3863)")
    pr(f"  r 중앙값: {np.median(r_vals):.4f}")

    # 이론 평균값 검증
    r_mean_gue = 4 - 2 * np.sqrt(3)  # ≈ 0.5359 (Atas et al.)
    r_mean_poisson = 2 * np.log(2) - 1  # ≈ 0.3863
    pr(f"  GUE 이론 평균: {r_mean_gue:.4f}")
    pr(f"  Poisson 이론 평균: {r_mean_poisson:.4f}")
    pr()

    # ──────────────────────────────────────────
    # 파트 C: N에 따른 r̄ 수렴 (N=2000, 5000, 10000, 20000 부분집합)
    # ──────────────────────────────────────────
    pr("=" * 70)
    pr("파트 C: N 부분집합별 r̄ 수렴 확인")
    pr("=" * 70)
    pr(f"  {'N_used':>8}  {'r̄':>8}  {'GUE 편차':>10}  {'Poisson 편차':>12}")
    pr(f"  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*12}")

    for N_use in [1000, 2000, 3000, len(r_vals)]:
        if N_use > len(r_vals):
            N_use = len(r_vals)
        r_sub = r_vals[:N_use]
        r_m = r_sub.mean()
        dev_gue = abs(r_m - r_mean_gue)
        dev_poi = abs(r_m - r_mean_poisson)
        pr(f"  {N_use:>8d}  {r_m:>8.4f}  {dev_gue:>10.4f}  {dev_poi:>12.4f}")
    pr()

    # ──────────────────────────────────────────
    # 파트 D: KS 검정 — GUE vs GOE vs Poisson
    # ──────────────────────────────────────────
    pr("=" * 70)
    pr("파트 D: KS 검정 — GUE / GOE / Poisson")
    pr("=" * 70)

    # 전체 데이터로 KS 검정
    # r_vals는 [0,1] 범위 — 이론 CDF 필요
    r_grid_gue, cdf_gue = build_cdf(p_gue)
    r_grid_goe, cdf_goe = build_cdf(p_goe)
    r_grid_poi, cdf_poi = build_cdf(p_poisson)

    # scipy kstest 사용 (이론 CDF 보간)
    D_gue, p_gue_ks = stats.kstest(
        r_vals,
        lambda x: np.interp(x, r_grid_gue, cdf_gue)
    )
    D_goe, p_goe_ks = stats.kstest(
        r_vals,
        lambda x: np.interp(x, r_grid_goe, cdf_goe)
    )
    D_poi, p_poi_ks = stats.kstest(
        r_vals,
        lambda x: np.interp(x, r_grid_poi, cdf_poi)
    )

    pr(f"  분포       KS통계량 D     p-값         판정")
    pr(f"  {'-'*55}")
    gue_verdict = "✅ 비기각 (GUE 적합)" if p_gue_ks > 0.05 else "⚠️ 기각"
    goe_verdict = "❌ 기각 (GOE 부적합)" if p_goe_ks < 0.01 else "⚠️ 미기각"
    poi_verdict = "❌ 기각 (Poisson 부적합)" if p_poi_ks < 0.01 else "⚠️ 미기각"
    pr(f"  GUE        {D_gue:.6f}      {p_gue_ks:.4e}    {gue_verdict}")
    pr(f"  GOE        {D_goe:.6f}      {p_goe_ks:.4e}    {goe_verdict}")
    pr(f"  Poisson    {D_poi:.6f}      {p_poi_ks:.4e}    {poi_verdict}")
    pr()

    # 부분집합별 KS 검정 (N=2000, 5000, 10000, 20000)
    pr("  부분집합별 GUE KS 검정:")
    pr(f"  {'N_used':>8}  {'D_GUE':>10}  {'p_GUE':>12}  {'판정':}")
    pr(f"  {'-'*8}  {'-'*10}  {'-'*12}  {'-'*15}")
    for N_use in [1000, 2000, 3000, len(r_vals)]:
        if N_use > len(r_vals):
            N_use = len(r_vals)
        r_sub = r_vals[:N_use]
        D_s, p_s = stats.kstest(
            r_sub,
            lambda x: np.interp(x, r_grid_gue, cdf_gue)
        )
        v = "비기각" if p_s > 0.05 else "기각"
        pr(f"  {N_use:>8d}  {D_s:>10.6f}  {p_s:>12.4e}  {v}")
    pr()

    # ──────────────────────────────────────────
    # 파트 E: 기술 통계 + 분위수 비교
    # ──────────────────────────────────────────
    pr("=" * 70)
    pr("파트 E: 기술 통계 + 이론 분위수 비교")
    pr("=" * 70)

    percentiles = [10, 25, 50, 75, 90]
    pr(f"  {'분위수':>6}  {'실측':>8}  {'GUE':>8}  {'GOE':>8}  {'Poisson':>9}")
    pr(f"  {'-'*6}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*9}")

    # 각 이론 분위수 계산 (역 CDF 보간)
    for pct in percentiles:
        q_obs = np.percentile(r_vals, pct)
        q_gue = np.interp(pct / 100, cdf_gue, r_grid_gue)
        q_goe = np.interp(pct / 100, cdf_goe, r_grid_goe)
        q_poi = np.interp(pct / 100, cdf_poi, r_grid_poi)
        pr(f"  {pct:>6}%  {q_obs:>8.4f}  {q_gue:>8.4f}  {q_goe:>8.4f}  {q_poi:>9.4f}")
    pr()

    # ──────────────────────────────────────────
    # 파트 F: 종합 판정
    # ──────────────────────────────────────────
    pr("=" * 70)
    pr("★ 종합 판정")
    pr("=" * 70)
    pr()

    positive_gue = p_gue_ks > 0.05
    reject_poisson = p_poi_ks < 0.01
    reject_goe = p_goe_ks < 0.01

    pr(f"  성공 기준 체크:")
    pr(f"    GUE 비기각 (p_GUE > 0.05):   {'✅' if positive_gue else '❌'}  p={p_gue_ks:.4e}")
    pr(f"    Poisson 기각 (p_poi < 0.01): {'✅' if reject_poisson else '❌'}  p={p_poi_ks:.4e}")
    pr(f"    GOE 기각 (p_goe < 0.01):     {'✅' if reject_goe else '❌'}  p={p_goe_ks:.4e}")
    pr()

    if positive_gue and reject_poisson and reject_goe:
        verdict = "★★ 강한 양성: GUE 비기각 + Poisson·GOE 모두 기각"
    elif positive_gue and reject_poisson:
        verdict = "★ 양성: GUE 비기각 + Poisson 기각 (GOE 미기각)"
    elif not positive_gue and reject_poisson:
        verdict = "⚠️ 부분 양성: GUE 기각되었으나 Poisson도 기각"
    else:
        verdict = "❌ 음성: GUE 기각"

    pr(f"  {verdict}")
    pr()
    pr(f"  r̄_실측 = {r_vals.mean():.4f}  (GUE: {r_mean_gue:.4f}, Poisson: {r_mean_poisson:.4f})")
    pr(f"  r̄ - r̄_GUE = {r_vals.mean() - r_mean_gue:+.4f}")
    pr()
    pr(f"총 소요 시간: {time.time() - t_start:.1f}초")
    pr(f"결과 파일: results/spacing_ratio_gue.txt")

    # 결과 저장
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    with open(OUTPUT, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n결과 저장 완료: {OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
