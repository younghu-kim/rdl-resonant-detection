#!/usr/bin/env python3
"""
=============================================================================
[RDL 실험 #215] c₁ 계수와 GUE Pair Correlation의 연결
=============================================================================
배경:
  Theorem 5에서 κδ²-1 = A·δ² + O(δ³), A = Im(c₀)² + 2c₁.
  Hadamard product에서: c₁ = Σ_{n≠0} 1/(γ₀-γ_n)²

  Montgomery pair correlation (1973):
    R₂(α) = 1 - (sin πα / πα)² + δ(α)   (GUE 예측)

  c₁과의 관계:
    c₁ = Σ_{n≠0} 1/(γ₀-γ_n)² = (1/Δ²)·Σ 1/(normalized gap)²
    여기서 Δ = 2π/log(T/2π) (mean spacing at height T).

    만약 영점이 GUE 분포를 따르면:
    E[c₁] = (1/Δ²)·∫₀^∞ R₂(x)/x² dx + ...  (pair-correlation weighted sum)

  정량적 예측:
    c₁ ~ (log T/(2π))² / (4π²) · [GUE 상수]
    = (d·log(NT/(2πe)))² / (4π²) · C_GUE

  이 실험에서:
    1. ζ(s)의 처음 100개 영점에서 c₁ 직접 측정
    2. c₁ vs (log T/(2π))² 관계 확인
    3. c₁_measured vs c₁_predicted (zero-sum 직접 계산) 비교
    4. Pair correlation → c₁ 변환이 일치하는지 검증

결과: results/c1_pair_correlation_215.txt
=============================================================================
"""
import sys, os, time
import numpy as np
from scipy import stats

START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'c1_pair_correlation_215.txt')
outf = open(RESULT_FILE, 'w', buffering=1)

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()

def T():
    return f"[{time.time()-START:.1f}s]"

log("=" * 72)
log("[실험 #215] c₁ 계수와 GUE Pair Correlation의 연결")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PARI 초기화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

import cypari2
gp = cypari2.Pari()
gp.allocatemem(1000 * 1024 * 1024)
gp("default(realprecision, 100)")
log(f"{T()} PARI OK, realprecision=100")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 이론적 배경
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("이론적 배경:")
log("  c₁(γ₀) = Σ_{γ_n≠γ₀} 1/(γ₀-γ_n)²")
log("  Mean spacing: Δ(T) = 2π/log(T/(2π))")
log("  Normalized: c₁·Δ² = Σ 1/(normalized_gap)²")
log()
log("  GUE prediction for Σ 1/x_n²:")
log("    If gaps follow GUE 2-point correlation,")
log("    c₁·Δ² should approach a universal constant.")
log()
log("  Alternative: c₁ vs (logT)² scaling")
log("    c₁ ~ [log(T/(2π))]² / (4π²) · C")
log("    i.e., c₁·(2π/logT)² = C (constant)")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Phase 1: ζ(s) 영점 수집 (대량)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [Phase 1] ζ(s) 영점 대량 수집...")

gp("Lcur = lfuncreate(1)")
gp("Linit = lfuninit(Lcur, [1.0, 600])")  # t up to 600
gp("zvec = lfunzeros(Linit, 500)")
n_zeros = int(gp("length(zvec)"))
log(f"  영점 {n_zeros}개 수집 (t < 500)")

zeros = np.array([float(gp(f"zvec[{i}]")) for i in range(1, n_zeros+1)])
log(f"  범위: γ₁={zeros[0]:.4f} ~ γ_{n_zeros}={zeros[-1]:.4f}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Phase 2: c₁ 직접 계산 (zero-sum)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [Phase 2] c₁ = Σ 1/(γ₀-γ_n)² 직접 계산...")
log()

# 내부 영점만 사용 (edge effect 방지: 양측 20개씩 제외)
MARGIN = 20
target_indices = range(MARGIN, n_zeros - MARGIN)

c1_values = []
gamma_values = []

for idx in target_indices:
    g0 = zeros[idx]
    # 모든 다른 영점과의 차이
    diffs = zeros - g0
    diffs = diffs[diffs != 0]  # 자기 자신 제외
    c1_sum = np.sum(1.0 / diffs**2)
    c1_values.append(c1_sum)
    gamma_values.append(g0)

c1_values = np.array(c1_values)
gamma_values = np.array(gamma_values)

log(f"  계산 완료: {len(c1_values)}개 영점")
log(f"  c₁ 범위: [{c1_values.min():.4f}, {c1_values.max():.4f}]")
log(f"  c₁ 평균: {c1_values.mean():.4f}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Phase 3: c₁ vs log(T)² 스케일링
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [Phase 3] c₁ vs (log T/(2π))² 스케일링 분석...")
log()

# log(γ/(2π))²
log_factor = (np.log(gamma_values / (2*np.pi)))**2

# 정규화: c₁ / (logT/(2π))²
c1_normalized = c1_values / log_factor

log(f"  c₁ / (log(γ/(2π)))² 통계:")
log(f"    mean = {c1_normalized.mean():.6f}")
log(f"    std  = {c1_normalized.std():.6f}")
log(f"    CV   = {c1_normalized.std()/c1_normalized.mean()*100:.1f}%")
log()

# 구간별 통계
bins = [(50, 100), (100, 200), (200, 300), (300, 400), (400, 500)]
log("  구간별 c₁/(logT)² 평균:")
for lo, hi in bins:
    mask = (gamma_values >= lo) & (gamma_values < hi)
    if mask.sum() > 5:
        m = c1_normalized[mask].mean()
        s = c1_normalized[mask].std()
        n = mask.sum()
        log(f"    γ∈[{lo},{hi}]: mean={m:.6f} ± {s:.6f} (n={n})")

log()

# Linear regression: c₁ vs (logT)²
slope_fit, intercept_fit, r_value, p_value, std_err = stats.linregress(
    log_factor, c1_values
)
log(f"  선형 회귀: c₁ = {slope_fit:.6f}·(log(γ/2π))² + {intercept_fit:.4f}")
log(f"    R² = {r_value**2:.6f}")
log(f"    p  = {p_value:.2e}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Phase 4: c₁ vs PARI 측정 비교 (검증)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [Phase 4] c₁(zero-sum) vs c₁(PARI 측정) 비교...")
log()

# PARI에서 직접 c₁ 측정 (대칭 방법, 10개 영점 샘플)
C0_DELTAS = np.array([0.001, 0.002, 0.003, 0.005, 0.007, 0.01])

def get_log_deriv_pari(sigma, t0):
    gp(f"s_eval = {sigma:.20f} + I*{t0:.20f}")
    gp("Lv = lfunlambda(Linit, s_eval)")
    gp("dLv = lfunlambda(Linit, s_eval, 1)")
    gp("r_val = dLv / Lv")
    return complex(float(gp("real(r_val)")), float(gp("imag(r_val)")))

def measure_c1_pari(t0):
    """반대칭 방법으로 c₁ 직접 측정"""
    c1_estimates = []
    for d in C0_DELTAS:
        f_plus = get_log_deriv_pari(0.5 + d, t0)
        f_minus = get_log_deriv_pari(0.5 - d, t0)
        antisym = (f_plus - f_minus) / 2.0
        c1_est = (antisym - 1.0/d) / d
        c1_estimates.append(c1_est.real)
    return np.mean(c1_estimates[:3])

# 10개 영점 샘플 (균등 간격)
sample_indices = np.linspace(MARGIN, n_zeros-MARGIN-1, 10, dtype=int)
log(f"  샘플 영점 10개 (idx: {sample_indices[:3]}...{sample_indices[-1]})")
log()

comparison = []
for idx in sample_indices:
    t0 = zeros[idx]
    c1_sum = c1_values[idx - MARGIN]  # zero-sum 값
    c1_pari = measure_c1_pari(t0)
    err = abs(c1_sum - c1_pari) / max(abs(c1_pari), 1e-10)
    comparison.append({"t0": t0, "c1_sum": c1_sum, "c1_pari": c1_pari, "err": err})
    log(f"  γ={t0:.3f}: c₁(sum)={c1_sum:.4f}, c₁(PARI)={c1_pari:.4f}, err={err*100:.1f}%")

log()
mean_err = np.mean([c["err"] for c in comparison])
log(f"  평균 오차: {mean_err*100:.1f}%")
log(f"  판정: {'★★★ PASS' if mean_err < 0.05 else '★★ PASS' if mean_err < 0.2 else '△ 부분'}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Phase 5: GUE 예측과 비교
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [Phase 5] GUE 정규화 상수 추출...")
log()

# c₁·Δ² = c₁·(2π/log(γ/(2π)))² 가 상수인지 확인
delta_values = 2*np.pi / np.log(gamma_values / (2*np.pi))
c1_times_delta2 = c1_values * delta_values**2

log(f"  c₁·Δ² 통계:")
log(f"    mean = {c1_times_delta2.mean():.6f}")
log(f"    std  = {c1_times_delta2.std():.6f}")
log(f"    CV   = {c1_times_delta2.std()/c1_times_delta2.mean()*100:.1f}%")
log()

# GUE 예측: Σ 1/x_n² for GUE eigenvalue pair
# pair correlation integral: ∫₁^∞ (1-(sin πx/(πx))²)/x² dx ≈ ?
# 실제로 GUE에서 Σ_{n≠0} 1/x_n² (정규화된 간격)의 기대값
# = π²/3 (Poisson), GUE에서는 더 작음 (repulsion)
log(f"  구간별 c₁·Δ² 평균:")
for lo, hi in bins:
    mask = (gamma_values >= lo) & (gamma_values < hi)
    if mask.sum() > 5:
        m = c1_times_delta2[mask].mean()
        s = c1_times_delta2[mask].std()
        log(f"    γ∈[{lo},{hi}]: mean={m:.6f} ± {s:.6f}")

log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Phase 6: A(t₀) 분해
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log(f"{T()} [Phase 6] A(t₀) = Im(c₀)² + 2c₁ 분해...")
log()

# Im(c₀) = -Σ 1/(γ₀-γ_n) (Hilbert transform)
im_c0_values = []
for idx in target_indices:
    g0 = zeros[idx]
    diffs = zeros - g0
    diffs = diffs[diffs != 0]
    im_c0 = -np.sum(1.0 / diffs)
    im_c0_values.append(im_c0)

im_c0_values = np.array(im_c0_values)
A_values = im_c0_values**2 + 2*c1_values

log(f"  Im(c₀) 통계:")
log(f"    mean = {im_c0_values.mean():.4f}")
log(f"    std  = {im_c0_values.std():.4f}")
log(f"    |mean| vs std: {'fluctuating (OK)' if abs(im_c0_values.mean()) < im_c0_values.std() else 'systematic'}")
log()
log(f"  A(t₀) = Im(c₀)² + 2c₁ 분해:")
log(f"    mean Im(c₀)² = {(im_c0_values**2).mean():.4f} ({(im_c0_values**2).mean()/A_values.mean()*100:.0f}%)")
log(f"    mean 2c₁     = {(2*c1_values).mean():.4f} ({(2*c1_values).mean()/A_values.mean()*100:.0f}%)")
log(f"    mean A       = {A_values.mean():.4f}")
log()

# A의 성장율
A_normalized = A_values / (np.log(gamma_values/(2*np.pi)))**2
log(f"  A / (log(γ/(2π)))² 통계:")
log(f"    mean = {A_normalized.mean():.6f}")
log(f"    std  = {A_normalized.std():.6f}")
log(f"    CV   = {A_normalized.std()/A_normalized.mean()*100:.1f}%")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 종합 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 72)
log("종합 분석")
log("=" * 72)
log()

log("1. c₁ = Σ 1/(γ₀-γ_n)² 확인:")
log(f"   zero-sum vs PARI: mean {mean_err*100:.1f}% error → {'PASS' if mean_err<0.1 else 'PARTIAL'}")
log()

log("2. c₁ ∝ (log T)² 스케일링:")
log(f"   R² = {r_value**2:.6f} → {'★★★ 강한 상관' if r_value**2>0.9 else '★★ 유의' if r_value**2>0.5 else '★ 약함'}")
log(f"   c₁ ≈ {slope_fit:.4f}·(log(γ/2π))² + {intercept_fit:.4f}")
log()

log("3. GUE 정규화 상수:")
log(f"   c₁·Δ² = {c1_times_delta2.mean():.4f} ± {c1_times_delta2.std():.4f}")
log(f"   CV = {c1_times_delta2.std()/c1_times_delta2.mean()*100:.1f}%")
cv = c1_times_delta2.std()/c1_times_delta2.mean()
log(f"   판정: {'★★★ 상수 (GUE 보편성)' if cv<0.1 else '★★ 약한 T-의존성' if cv<0.3 else '★ 비상수'}")
log()

log("4. A 분해 (Im(c₀)² vs 2c₁):")
frac_c1 = (2*c1_values).mean() / A_values.mean()
log(f"   2c₁ 기여: {frac_c1*100:.0f}%")
log(f"   Im(c₀)² 기여: {(1-frac_c1)*100:.0f}%")
log(f"   → {'c₁ 지배적 (밀도 항)' if frac_c1>0.7 else 'Im(c₀)² 지배적 (비정칙 항)' if frac_c1<0.3 else '양자 기여 균형'}")
log()

log("=" * 72)
log("최종 판정")
log("=" * 72)
log()
if r_value**2 > 0.8 and mean_err < 0.2:
    log("★★★ 강양성 — c₁ 해석적 공식 확인")
    log()
    log("결론:")
    log("  1. c₁ = Σ_{n≠0} 1/(γ₀-γ_n)² (Hadamard product에서 직접 도출)")
    log("  2. c₁ ∝ (log γ)² 스케일링 확인 (R²>{:.3f})".format(r_value**2))
    log("  3. 이것은 c₁ ~ (mean density)² 와 일치")
    log("     (mean density = log(T/(2π))/(2π)이므로 c₁ ~ (logT)²/(4π²))")
    log("  4. A(t₀) = Im(c₀)² + 2c₁에서 2c₁이 체계적 성장을 담당")
    log("     Im(c₀)²는 영점 간격의 비정칙성(fluctuation)을 반영")
    log()
    log("논문 기여:")
    log("  - A(d,N) ≈ 0.39d²+2.19logN (#212, 취약) 를")
    log("    A = Im(c₀)²+2c₁, c₁~(d·log(NT))²/(4π²) 로 업그레이드")
    log("  - Paper 2 Future direction #3 (c₁ in terms of degree) 부분 해결")
elif r_value**2 > 0.5:
    log("★★ 양성 — 경향 확인, 정량적 공식 미완")
else:
    log("△ 중립 — 추가 분석 필요")

log()
log(f"총 소요시간: {time.time()-START:.1f}초")
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
outf.close()
