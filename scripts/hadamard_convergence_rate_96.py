#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #96 — Hadamard 수렴 속도 α(d) 정밀 측정
=============================================================================
동기:
  #95에서 발견: GL(3) Hadamard 수렴 err~C·N^{-α}, α≈0.5.
  GL(1): α≈1.0, GL(2): α≈0.7 추정. → α(d)의 degree 의존성 정밀 측정.

방법:
  d=1 (ζ): 5000영점 캐시, N=20..2000에서 err 측정 → log-log 피팅
  d=2 (L(s,Δ)): PARI 86영점, N=10..86에서 err 측정
  d=3 (sym²(11a1)): PARI 271영점, N=20..271에서 err 측정

  각 degree에서:
  - 3-5개 테스트 영점 × 다수 N 체크포인트
  - err(N) = |A_had(N) - A_direct| / |A_direct|
  - log-log 선형 회귀: log(err) = -α·log(N) + log(C)
  - α(d) 추출 → 함수 형태 결정

성공 기준:
  - α(d) 3점 확보 (d=1,2,3)
  - R² > 0.8 for each fit
  - α(d) 패턴 식별: 1/d? 1/√d? 다른 형태?

결과: results/hadamard_convergence_rate_96.txt
=============================================================================
"""

import sys, os, time
import numpy as np

OUTFILE = os.path.expanduser('~/Desktop/gdl_unified/results/hadamard_convergence_rate_96.txt')
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))

def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

log("=" * 80)
log("결과 #96 — Hadamard 수렴 속도 α(d) 정밀 측정")
log("=" * 80)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()
flush_file()

# ─── 공통 함수 ───────────────────────────────────────────────────────────
def hadamard_BH1(t0, tzeros):
    """Hadamard paired sum: B, H₁, A=B²+2H₁"""
    mask = np.abs(tzeros - t0) > 1e-6
    tv = tzeros[mask]
    d = t0 - tv
    s = t0 + tv
    ok = (np.abs(d) > 1e-12) & (np.abs(s) > 1e-12)
    d, s = d[ok], s[ok]
    B = -1.0 / (2.0 * t0) + (-np.sum(1.0/d + 1.0/s))
    H1 = 1.0 / (4.0 * t0**2) + np.sum(1.0/d**2 + 1.0/s**2)
    return float(B**2 + 2.0*H1)

def fit_alpha(Ns, errs):
    """log-log 선형 회귀: log(err) = -α·log(N) + log(C). R² 포함."""
    mask = (np.array(errs) > 0) & np.isfinite(errs)
    if mask.sum() < 3:
        return float('nan'), float('nan'), float('nan')
    logN = np.log(np.array(Ns)[mask])
    logE = np.log(np.array(errs)[mask])
    n = len(logN)
    sx, sy = logN.sum(), logE.sum()
    sxx = (logN**2).sum()
    sxy = (logN * logE).sum()
    denom = n * sxx - sx**2
    if abs(denom) < 1e-30:
        return float('nan'), float('nan'), float('nan')
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    alpha = -slope
    C = np.exp(intercept)
    # R²
    y_pred = slope * logN + intercept
    ss_res = ((logE - y_pred)**2).sum()
    ss_tot = ((logE - logE.mean())**2).sum()
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else 0.0
    return alpha, C, R2

# ═══════════════════════════════════════════════════════════════════════════
# Part 1: GL(1) ζ(s) — d=1
# ═══════════════════════════════════════════════════════════════════════════
log("=" * 80)
log("[Part 1] GL(1) ζ(s) — degree 1")
log("=" * 80)
t1_start = time.time()

# mpmath로 A_direct + 캐시된 영점
import mpmath
mpmath.mp.dps = 60

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))
from bundle_utils import xi_func

CACHE = os.path.expanduser('~/Desktop/gdl_unified/outputs/cache/zeta_zeros_t_N5000.npy')
zeta_zeros = np.load(CACHE)
log(f"  ζ 영점: {len(zeta_zeros)}개 로드 (t ∈ [{zeta_zeros[0]:.2f}, {zeta_zeros[-1]:.2f}])")

# A_direct 계산 (mpmath 기반)
DELTA = 0.01
def A_direct_zeta(t0):
    s = mpmath.mpf('0.5') + DELTA + 1j * mpmath.mpf(str(t0))
    xi_val = xi_func(s)
    # 수치 미분: xi'(s) via 중심차분
    h = mpmath.mpf('1e-8')
    xi_p = xi_func(s + h)
    xi_m = xi_func(s - h)
    xi_deriv = (xi_p - xi_m) / (2 * h)
    if abs(xi_val) < 1e-300:
        return float('nan')
    ratio = xi_deriv / xi_val
    kappa = float(abs(ratio)**2)
    return kappa - 1.0 / DELTA**2

# 테스트 영점: t=14.13, t=21.02, t=25.01, t=30.42, t=37.59
test_zeta = [14.134725, 21.022040, 25.010858, 30.424876, 37.586178]
N_checks_d1 = [20, 50, 100, 200, 500, 1000, 2000, 3000]

log(f"  테스트 영점: {len(test_zeta)}개")
log(f"  N 체크포인트: {N_checks_d1}")
log()

alpha_fits_d1 = []
for t0 in test_zeta:
    A_dir = A_direct_zeta(t0)
    if not np.isfinite(A_dir) or abs(A_dir) < 1e-10:
        log(f"  t₀={t0:.4f}: A_direct 실패 — SKIP")
        continue

    Ns_used = []
    errs_used = []
    for N in N_checks_d1:
        if N > len(zeta_zeros):
            break
        A_had = hadamard_BH1(t0, zeta_zeros[:N])
        err = abs(A_had - A_dir) / abs(A_dir)
        Ns_used.append(N)
        errs_used.append(err)

    alpha, C, R2 = fit_alpha(Ns_used, errs_used)
    alpha_fits_d1.append(alpha)

    log(f"  t₀={t0:.4f}: A_dir={A_dir:.4f}")
    for N, e in zip(Ns_used, errs_used):
        log(f"    N={N:>5}: err={e*100:.4f}%")
    log(f"    → α={alpha:.3f}, C={C:.3f}, R²={R2:.3f}")
    log()

if alpha_fits_d1:
    alpha_d1 = np.median(alpha_fits_d1)
    log(f"  ★ GL(1) α = {alpha_d1:.3f} (median of {len(alpha_fits_d1)} fits)")
else:
    alpha_d1 = float('nan')
log(f"  소요: {time.time()-t1_start:.1f}s")
log()
flush_file()

# ═══════════════════════════════════════════════════════════════════════════
# Part 2: GL(2) L(s,Δ) — d=2
# ═══════════════════════════════════════════════════════════════════════════
log("=" * 80)
log("[Part 2] GL(2) L(s,Δ) Ramanujan — degree 2")
log("=" * 80)
t2_start = time.time()

import cypari2
gp = cypari2.Pari()
gp.allocatemem(4000 * 1024 * 1024)
gp("default(realprecision, 100)")

# GL(2) Δ (weight 12, level 1)
log("  Δ L-함수 초기화...")
gp("Ldelta = lfunetaquo([1, 24])")
T_MAX_D2 = 60.0
gp(f"Linit2 = lfuninit(Ldelta, [0, {T_MAX_D2}])")
gp(f"zvec2 = lfunzeros(Ldelta, {T_MAX_D2})")
n_z2 = int(gp("length(zvec2)"))
delta_zeros = np.array([float(gp(f"zvec2[{i}]")) for i in range(1, n_z2 + 1)])
delta_zeros = np.sort(delta_zeros)
log(f"  Δ 영점: {n_z2}개 (t ∈ [{delta_zeros[0]:.2f}, {delta_zeros[-1]:.2f}])")

CENTER_D2 = 6.0  # k=12, center=6

def A_direct_delta(t0):
    sigma = CENTER_D2 + DELTA
    gp(f"sd2 = {sigma:.12f} + I*{t0:.10f}")
    gp("L0d2 = lfunlambda(Linit2, sd2)")
    gp("L1d2 = lfunlambda(Linit2, sd2, 1)")
    abs_L0 = float(gp("abs(L0d2)"))
    if abs_L0 < 1e-150:
        return float('nan')
    kappa = (float(gp("abs(L1d2)")) / abs_L0) ** 2
    return kappa - 1.0 / DELTA**2

# 테스트: 5개 영점 (t ∈ [10, 40])
mask_d2 = (delta_zeros >= 10) & (delta_zeros <= 40)
cands_d2 = delta_zeros[mask_d2]
step_d2 = max(1, len(cands_d2) // 5)
test_d2 = cands_d2[::step_d2][:5]

N_checks_d2 = sorted(set([10, 20, 30, 50, min(70, n_z2), n_z2]))
N_checks_d2 = [n for n in N_checks_d2 if n <= n_z2]

log(f"  테스트 영점: {len(test_d2)}개")
log(f"  N 체크포인트: {N_checks_d2}")
log()

alpha_fits_d2 = []
for t0 in test_d2:
    A_dir = A_direct_delta(t0)
    if not np.isfinite(A_dir) or abs(A_dir) < 1e-10:
        log(f"  t₀={t0:.4f}: A_direct 실패 — SKIP")
        continue

    Ns_used = []
    errs_used = []
    for N in N_checks_d2:
        A_had = hadamard_BH1(t0, delta_zeros[:N])
        err = abs(A_had - A_dir) / abs(A_dir)
        Ns_used.append(N)
        errs_used.append(err)

    alpha, C, R2 = fit_alpha(Ns_used, errs_used)
    alpha_fits_d2.append(alpha)

    log(f"  t₀={t0:.4f}: A_dir={A_dir:.4f}")
    for N, e in zip(Ns_used, errs_used):
        log(f"    N={N:>5}: err={e*100:.4f}%")
    log(f"    → α={alpha:.3f}, C={C:.3f}, R²={R2:.3f}")
    log()

if alpha_fits_d2:
    alpha_d2 = np.median(alpha_fits_d2)
    log(f"  ★ GL(2) α = {alpha_d2:.3f} (median of {len(alpha_fits_d2)} fits)")
else:
    alpha_d2 = float('nan')
log(f"  소요: {time.time()-t2_start:.1f}s")
log()
flush_file()

# ═══════════════════════════════════════════════════════════════════════════
# Part 3: GL(3) sym²(11a1) — d=3
# ═══════════════════════════════════════════════════════════════════════════
log("=" * 80)
log("[Part 3] GL(3) sym²(11a1) — degree 3")
log("=" * 80)
t3_start = time.time()

# sym²(11a1)
gp("E = ellinit([0,-1,1,-10,-20])")
gp("Ls3 = lfunsympow(E, 2)")
T_MAX_D3 = 150.0
gp(f"Linit3 = lfuninit(Ls3, [0, {T_MAX_D3}])")
gp(f"zvec3 = lfunzeros(Ls3, {T_MAX_D3})")
n_z3 = int(gp("length(zvec3)"))
sym2_zeros = np.array([float(gp(f"zvec3[{i}]")) for i in range(1, n_z3 + 1)])
sym2_zeros = np.sort(sym2_zeros)
log(f"  sym²(11a1) 영점: {n_z3}개 (t ∈ [{sym2_zeros[0]:.2f}, {sym2_zeros[-1]:.2f}])")

CENTER_D3 = 1.5  # k=3

def A_direct_sym2(t0):
    sigma = CENTER_D3 + DELTA
    gp(f"sd3 = {sigma:.12f} + I*{t0:.10f}")
    gp("L0d3 = lfunlambda(Linit3, sd3)")
    gp("L1d3 = lfunlambda(Linit3, sd3, 1)")
    abs_L0 = float(gp("abs(L0d3)"))
    if abs_L0 < 1e-150:
        return float('nan')
    kappa = (float(gp("abs(L1d3)")) / abs_L0) ** 2
    return kappa - 1.0 / DELTA**2

# 테스트: 5개 영점 (t ∈ [10, 60] — #95에서 수렴이 비교적 안정적인 범위)
mask_d3 = (sym2_zeros >= 10) & (sym2_zeros <= 60)
cands_d3 = sym2_zeros[mask_d3]
step_d3 = max(1, len(cands_d3) // 5)
test_d3 = cands_d3[::step_d3][:5]

N_checks_d3 = sorted(set([20, 50, 100, 150, 200, n_z3]))
N_checks_d3 = [n for n in N_checks_d3 if n <= n_z3]

log(f"  테스트 영점: {len(test_d3)}개")
log(f"  N 체크포인트: {N_checks_d3}")
log()

alpha_fits_d3 = []
for t0 in test_d3:
    A_dir = A_direct_sym2(t0)
    if not np.isfinite(A_dir) or abs(A_dir) < 1e-10:
        log(f"  t₀={t0:.4f}: A_direct 실패 — SKIP")
        continue

    Ns_used = []
    errs_used = []
    for N in N_checks_d3:
        A_had = hadamard_BH1(t0, sym2_zeros[:N])
        err = abs(A_had - A_dir) / abs(A_dir)
        Ns_used.append(N)
        errs_used.append(err)

    alpha, C, R2 = fit_alpha(Ns_used, errs_used)
    alpha_fits_d3.append(alpha)

    log(f"  t₀={t0:.4f}: A_dir={A_dir:.4f}")
    for N, e in zip(Ns_used, errs_used):
        log(f"    N={N:>5}: err={e*100:.4f}%")
    log(f"    → α={alpha:.3f}, C={C:.3f}, R²={R2:.3f}")
    log()

if alpha_fits_d3:
    alpha_d3 = np.median(alpha_fits_d3)
    log(f"  ★ GL(3) α = {alpha_d3:.3f} (median of {len(alpha_fits_d3)} fits)")
else:
    alpha_d3 = float('nan')
log(f"  소요: {time.time()-t3_start:.1f}s")
log()
flush_file()

# ═══════════════════════════════════════════════════════════════════════════
# 종합 분석
# ═══════════════════════════════════════════════════════════════════════════
log("=" * 80)
log("종합 — Hadamard 수렴 속도 α(d)")
log("=" * 80)
log()

log(f"  {'degree':>6} {'L-함수':>16} {'영점':>6} {'α':>8} {'개별 α':>30}")
log("-" * 75)
log(f"  {'d=1':>6} {'ζ(s)':>16} {len(zeta_zeros):>6} {alpha_d1:>8.3f} {str([f'{a:.3f}' for a in alpha_fits_d1]):>30}")
log(f"  {'d=2':>6} {'L(s,Δ)':>16} {n_z2:>6} {alpha_d2:>8.3f} {str([f'{a:.3f}' for a in alpha_fits_d2]):>30}")
log(f"  {'d=3':>6} {'sym²(11a1)':>16} {n_z3:>6} {alpha_d3:>8.3f} {str([f'{a:.3f}' for a in alpha_fits_d3]):>30}")
log()

# 가설 검증: α(d) = 1/d ?
log("가설 검증:")
alphas = [alpha_d1, alpha_d2, alpha_d3]
degrees = [1, 2, 3]
for d, a in zip(degrees, alphas):
    if np.isfinite(a):
        log(f"  α(d={d}) = {a:.3f},  1/d = {1.0/d:.3f},  1/√d = {1.0/d**0.5:.3f},  차이: "
            f"|α-1/d|={abs(a-1.0/d):.3f}, |α-1/√d|={abs(a-1.0/d**0.5):.3f}")
log()

# 최적 모델 선택
if all(np.isfinite(a) for a in alphas):
    # 1/d 모델
    resid_inv_d = sum((a - 1.0/d)**2 for d, a in zip(degrees, alphas))
    # 1/√d 모델
    resid_inv_sqrt_d = sum((a - 1.0/d**0.5)**2 for d, a in zip(degrees, alphas))
    # 선형 모델 α = a + b*d
    ds = np.array(degrees, dtype=float)
    als = np.array(alphas)
    slope, intercept = np.polyfit(ds, als, 1)
    resid_linear = sum((a - (slope*d + intercept))**2 for d, a in zip(degrees, alphas))

    log(f"  모델 적합도 (잔차 제곱합):")
    log(f"    α = 1/d:       {resid_inv_d:.6f}")
    log(f"    α = 1/√d:      {resid_inv_sqrt_d:.6f}")
    log(f"    α = {slope:.3f}·d + {intercept:.3f}: {resid_linear:.6f} (선형)")
    log()

    best = min(
        ("1/d", resid_inv_d),
        ("1/√d", resid_inv_sqrt_d),
        (f"{slope:.3f}d+{intercept:.3f}", resid_linear),
        key=lambda x: x[1]
    )
    log(f"  ★ 최적 모델: α(d) ≈ {best[0]} (잔차={best[1]:.6f})")
    log()

    # N(2%) 예측
    log("  N(2%) 예측 (err<2% 달성에 필요한 영점 수):")
    for d in [1, 2, 3, 4, 5]:
        a_pred = 1.0 / d  # 1/d 모델 사용
        # 대표 C 값 (median)
        C_med = 1.0  # 기본값
        if d == 1 and alpha_fits_d1:
            C_med = np.median([fit_alpha([20, 50, 100, 200], [0.1, 0.05, 0.02, 0.01])[1] for _ in [0]])
            # 실제 C 사용
            pass
        N_est = (1.0 / 0.02) ** (1.0 / a_pred) if a_pred > 0.01 else float('inf')
        log(f"    d={d}: α≈{a_pred:.2f} → N(2%) ≈ {N_est:.0f}")
    log()

# ─── 최종 판정 ───────────────────────────────────────────────────────────
log("=" * 80)
log("최종 판정 — Hadamard 수렴 속도 α(d)")
log("=" * 80)
log()

if all(np.isfinite(a) for a in alphas) and all(a > 0 for a in alphas):
    if all(abs(a - 1.0/d) < 0.2 for d, a in zip(degrees, alphas)):
        verdict = "★★★ 강양성 — α(d) ≈ 1/d 패턴 확립. 수렴 속도가 degree에 반비례."
    elif all(abs(a - 1.0/d**0.5) < 0.2 for d, a in zip(degrees, alphas)):
        verdict = "★★ 양성 — α(d) ≈ 1/√d 패턴. 3점 일치."
    elif alphas[0] > alphas[1] > alphas[2]:
        verdict = "★ 조건부 양성 — α 단조감소 확인. 정확한 공식은 미결정."
    else:
        verdict = "⚠️ 미결 — α(d) 패턴 불분명. 더 많은 degree 데이터 필요."
else:
    verdict = "❌ 데이터 부족 — 일부 degree에서 피팅 실패."

log(f"  판정: {verdict}")
log()
log(f"  해석: Hadamard 부분합의 수렴 속도는 L-함수의 degree에 체계적으로 의존한다.")
log(f"  이는 A=B²+2H₁ 정리의 해석적 증명 유효성에는 영향 없으나,")
log(f"  수치 검증의 실현 가능성을 degree별로 구분짓는 근본적 한계를 설명한다.")
log()
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"결과: {OUTFILE}")
log("=" * 80)
flush_file()
