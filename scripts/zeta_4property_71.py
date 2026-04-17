#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #71 — GL(1) ζ(s) 4성질 검증 + κ_near(d=1) 측정
=============================================================================
대상: Riemann zeta function ζ(s)
  - degree d=1, conductor N=1, trivial character
  - 감마 인자: π^{-s/2} Γ(s/2)
  - 함수방정식: Λ(s) = π^{-s/2} Γ(s/2) ζ(s) → Λ(s) = Λ(1-s)
  - ξ(s) = (1/2)s(s-1)Λ(s) — 정함수

Dirichlet 계수: c(n) = 1 for all n (완전 곱셈적)

목적:
  1. κ_near(d=1) 측정 — κ_near(d) 함수 형태 결정 핵심 데이터
     현재: GL(2)=1114.9, GL(3)=1125.2 → GL(1) 추가로 3점 확보
  2. FE 형식적 확인 (ζ라면 당연 ~0)
  3. mono/π = 2.0 확인 (≥80% TP)
  4. σ-유일성 PASS 확인 (N=1이므로 5번째 PASS 예상)
  5. κ_near 비교표: GL(1) vs GL(2) 1114.9 vs GL(3) 1125.2

검증 범위: t ∈ [10, 60] (γ₁≈14.13 ~ γ₅₀≈52.97)
DPS=50, c=2.0

결과 파일: results/zeta_4property_71.txt
=============================================================================
"""

import sys, os, time
import numpy as np
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

DPS = 50
mpmath.mp.dps = DPS

from bundle_utils import (
    xi_func, connection_zeta, curvature_zeta,
    monodromy_contour, find_zeros_zeta,
)

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "zeta_4property_71.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

# ━━━━━━━━━━━ 로깅 ━━━━━━━━━━━
lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))
def flush():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━━━━━━ 파라미터 ━━━━━━━━━━━
T_MIN, T_MAX = 10.0, 60.0
DELTA = 0.03          # κ near 오프셋
MONO_R = 0.4          # 모노드로미 반지름
MONO_N = 64           # 모노드로미 단계

# Odlyzko 고정밀 영점 (처음 50개, 검증용)
KNOWN_ZEROS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
    52.970321,
]

# ━━━━━━━━━━━ 곡률 함수 (δ 오프셋 지원) ━━━━━━━━━━━
def curvature_at(t, delta=DELTA):
    """ξ'/ξ at s = 0.5+δ+it → κ = |conn|²"""
    s = mpmath.mpc(mpmath.mpf('0.5') + mpmath.mpf(delta), t)
    try:
        xi_val = xi_func(s)
        if abs(xi_val) < mpmath.mpf(10)**(-DPS + 10):
            return float('inf')
        h = mpmath.mpf(1) / mpmath.mpf(10**20)
        xi_d = (xi_func(s + h) - xi_func(s - h)) / (2 * h)
        conn = xi_d / xi_val
        k = float(abs(conn)**2)
        return k if np.isfinite(k) else float('inf')
    except Exception as e:
        print(f"  WARNING curvature at t={t}, δ={delta}: {e}", flush=True)
        return 0.0


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#                   메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 70)
log("결과 #71 — GL(1) ζ(s) 4성질 검증 + κ_near(d=1) 측정")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"DPS={DPS}")
log(f"대상: Riemann zeta function ζ(s)")
log(f"degree=1, conductor N=1, trivial character")
log(f"감마: π^{{-s/2}} Γ(s/2)")
log(f"ξ(s) = (1/2)s(s-1)π^{{-s/2}}Γ(s/2)ζ(s)")
log(f"t 범위: [{T_MIN}, {T_MAX}]")
log(f"δ = {DELTA}, mono_r = {MONO_R}, mono_n = {MONO_N}")
log()

t_total = time.time()

# ━━━━━━ Step 1: 함수방정식 검증 ━━━━━━

log("[Step 1] 함수방정식 ξ(s) = ξ(1-s)")
fe_pts = [
    mpmath.mpc(0.5, 10),
    mpmath.mpc(0.5, 20),
    mpmath.mpc(0.5, 30),
    mpmath.mpc(0.5, 40),
    mpmath.mpc(0.5, 50),
]
fe_ok = True
fe_max_rel = 0.0

for sp in fe_pts:
    xi_s = xi_func(sp)
    xi_1s = xi_func(1 - sp)

    if abs(xi_s) > mpmath.mpf(10)**(-DPS + 15):
        rel = float(abs(xi_s - xi_1s) / abs(xi_s))
    else:
        rel = 0.0

    fe_max_rel = max(fe_max_rel, rel)
    ok = rel < 1e-8
    if not ok:
        fe_ok = False
    log(f"  s={mpmath.nstr(sp, 6)}: |ξ|={float(abs(xi_s)):.4e}, rel={rel:.2e} {'✅' if ok else '❌'}")

if fe_ok:
    log(f"  ✅ 함수방정식 통과 (max_rel={fe_max_rel:.2e})")
else:
    log(f"  ❌ 함수방정식 실패 (max_rel={fe_max_rel:.2e})")
log()
flush()

# ━━━━━━ Step 2: 영점 탐색 (mpmath.zetazero 사용) ━━━━━━

log(f"[Step 2] 영점 탐색 (t ∈ [{T_MIN}, {T_MAX}])")
t_zs = time.time()
zeros = find_zeros_zeta(T_MIN, T_MAX)
dt_zs = time.time() - t_zs
n_zeros = len(zeros)

log(f"  발견 영점: {n_zeros}개 ({dt_zs:.1f}s)")
if n_zeros == 0:
    log("  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
    flush()
    sys.exit(1)

for j, z in enumerate(zeros[:30]):
    log(f"    #{j+1}: γ = {z:.10f}")
if n_zeros > 30:
    log(f"    ... ({n_zeros - 30}개 추가)")

# Odlyzko 영점과 매칭
log(f"\n  Odlyzko 영점 매칭:")
for kz in KNOWN_ZEROS:
    dists = np.abs(zeros - kz)
    best_idx = np.argmin(dists)
    best_err = dists[best_idx]
    log(f"    known γ={kz:.6f} → found γ={zeros[best_idx]:.10f} (오차 {best_err:.6f})")

log()
flush()

# ━━━━━━ Step 3: κ near/far 측정 + 모노드로미 ━━━━━━

log(f"[Step 3] κ near/far + 모노드로미 (≥20 TP 목표)")

# TP: 영점에서 κ 측정
n_tp = min(30, n_zeros)  # 수학자 요구: ≥20 TP
kappa_near = []
mono_results = []

for j in range(n_tp):
    z_t = zeros[j]
    t_k = time.time()
    k_val = curvature_at(z_t, delta=DELTA)
    m_val = monodromy_contour(z_t, radius=MONO_R, n_steps=MONO_N, func='zeta')
    m_pi = abs(m_val) / np.pi if m_val is not None else None
    dt_k = time.time() - t_k
    kappa_near.append(k_val)
    mono_results.append(m_pi)
    m_str = f"{m_pi:.4f}" if m_pi is not None else "FAIL"
    log(f"    TP #{j+1:>2} γ={z_t:.8f}: κ={k_val:.2f}, mono/π={m_str} ({dt_k:.1f}s)")
    if (j + 1) % 5 == 0:
        flush()

# FP: 영점 사이 중간점
kappa_far = []
fp_ts = []
for j in range(min(n_tp - 1, len(zeros) - 1)):
    mid = (zeros[j] + zeros[j + 1]) / 2
    fp_ts.append(mid)

n_fp = min(15, len(fp_ts))
for j in range(n_fp):
    t = fp_ts[j]
    t_k = time.time()
    k_val = curvature_at(t, delta=DELTA)
    dt_k = time.time() - t_k
    kappa_far.append(k_val)
    log(f"    FP #{j+1:>2} t={t:.6f}: κ={k_val:.4f} ({dt_k:.1f}s)")
    if (j + 1) % 5 == 0:
        flush()

# κ 통계
near_fin = [k for k in kappa_near if np.isfinite(k) and 0 < k < 1e15]
far_fin = [k for k in kappa_far if np.isfinite(k) and 0 < k < 1e15]

if near_fin and far_fin:
    near_med = float(np.median(near_fin))
    near_mean = float(np.mean(near_fin))
    near_std = float(np.std(near_fin))
    near_cv = near_std / near_mean * 100 if near_mean > 0 else 999
    far_med = float(np.median(far_fin))
    ratio = near_med / far_med if far_med > 0 else float('inf')
else:
    near_med = near_mean = near_std = 0
    near_cv = 999
    far_med = 0
    ratio = 0

log(f"\n  κ_near median = {near_med:.2f} (n={len(near_fin)}, CV={near_cv:.1f}%)")
log(f"  κ_near mean   = {near_mean:.2f} ± {near_std:.2f}")
log(f"  κ_far median  = {far_med:.4f} (n={len(far_fin)})")
log(f"  κ ratio = {ratio:.1f}×")

# 모노드로미 요약
mono_pass = sum(1 for m in mono_results if m is not None and 1.5 < m < 2.5)
total_mono = len([m for m in mono_results if m is not None])
log(f"  모노드로미: {mono_pass}/{total_mono} (mono/π≈2.0)")

log()
flush()

# ━━━━━━ Step 4: σ-유일성 테스트 ━━━━━━

log(f"[Step 4] σ-유일성 테스트 ({n_zeros}개 영점)")
sigma_vals = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60]
n_sigma_test = min(n_zeros, 30)  # 최대 30개 영점 테스트
sigma_pass_count = 0
sigma_test_count = 0

for z_idx in range(n_sigma_test):
    z_t = zeros[z_idx]
    kappas_by_sigma = []
    for sigma in sigma_vals:
        k_val = curvature_at(z_t, delta=sigma - 0.5)
        kappas_by_sigma.append(k_val)

    # σ=0.5에서의 κ가 최대인지 확인
    idx_05 = sigma_vals.index(0.50)
    k_05 = kappas_by_sigma[idx_05]

    # κ(0.5)/κ(0.45) 비율 (σ-유일성 날카로움 정량)
    idx_045 = sigma_vals.index(0.45)
    k_045 = kappas_by_sigma[idx_045]
    sharpness = k_05 / k_045 if k_045 > 0 and np.isfinite(k_045) else float('inf')

    is_max = all(k_05 >= k for k in kappas_by_sigma if np.isfinite(k))
    sigma_test_count += 1
    if is_max:
        sigma_pass_count += 1

    kstr = ", ".join(f"{k:.2e}" for k in kappas_by_sigma)
    log(f"  #{z_idx+1:>2} γ={z_t:.6f}: {'✅' if is_max else '❌'} κ(0.5)={k_05:.4e}, κ(0.5)/κ(0.45)={sharpness:.2e}")
    if (z_idx + 1) % 10 == 0:
        flush()

log(f"\n  σ-유일성 결과: {sigma_pass_count}/{sigma_test_count}")
log()
flush()

# ━━━━━━ Step 5: 저γ 편향 분석 ━━━━━━

log("[Step 5] 저γ κ_near 편향 분석")
if len(near_fin) >= 10:
    low_gamma = [kappa_near[j] for j in range(len(kappa_near)) if zeros[j] < 20 and np.isfinite(kappa_near[j]) and kappa_near[j] < 1e15]
    high_gamma = [kappa_near[j] for j in range(len(kappa_near)) if zeros[j] >= 20 and np.isfinite(kappa_near[j]) and kappa_near[j] < 1e15]
    if low_gamma and high_gamma:
        low_med = float(np.median(low_gamma))
        high_med = float(np.median(high_gamma))
        log(f"  γ<20: median κ={low_med:.2f} (n={len(low_gamma)})")
        log(f"  γ≥20: median κ={high_med:.2f} (n={len(high_gamma)})")
        log(f"  차이: {abs(high_med - low_med):.2f} ({abs(high_med - low_med)/high_med*100:.2f}%)")
    else:
        log(f"  데이터 부족 (low={len(low_gamma)}, high={len(high_gamma)})")
else:
    log(f"  데이터 부족 (n={len(near_fin)})")
log()
flush()

# ━━━━━━ 최종 요약 ━━━━━━

total_time = time.time() - t_total

log("=" * 70)
log("최종 요약 — GL(1) ζ(s) 4성질 검증")
log("=" * 70)
log(f"대상: Riemann zeta function ζ(s)")
log(f"degree=1, conductor N=1, trivial character")
log(f"감마: π^{{-s/2}} Γ(s/2)")
log(f"")
log(f"1. 함수방정식: {'✅ PASS' if fe_ok else '❌ FAIL'} (max_rel={fe_max_rel:.2e})")
log(f"2. 영점 발견: {n_zeros}개 (t ∈ [{T_MIN}, {T_MAX}])")
log(f"3. κ_near median = {near_med:.2f} (n={len(near_fin)}, CV={near_cv:.1f}%)")
log(f"   κ_near mean   = {near_mean:.2f} ± {near_std:.2f}")
log(f"   κ_far median  = {far_med:.4f}")
log(f"   κ ratio = {ratio:.1f}×")
log(f"4. 모노드로미: {mono_pass}/{total_mono} (mono/π≈2.0)")
log(f"5. σ-유일성: {sigma_pass_count}/{sigma_test_count}")
log(f"")

# κ_near 비교표 (핵심!)
log("=" * 70)
log("κ_near degree별 비교표")
log("=" * 70)
log(f"  GL(1) ζ(s)     degree=1: κ_near = {near_med:.2f} (n={len(near_fin)}, CV={near_cv:.1f}%) ← 이번 실험")
log(f"  GL(2) Maass avg degree=2: κ_near = 1114.9 ± 0.5 (n=20, 2 forms)")
log(f"  GL(2) Δ (w=12)  degree=2: κ_near = 1114.13 ± 1.4 (n=12)")
log(f"  GL(2) 합산       degree=2: κ_near ≈ 1114.6 (n=34, 3 forms)")
log(f"  GL(3) sym²       degree=3: κ_near = 1125.16 ± 1.7 (n=~60, 6곡선)")
log(f"")

# 추세 분석
if near_med > 0:
    diff_gl2 = abs(near_med - 1114.6) / 1114.6 * 100
    diff_gl3 = abs(near_med - 1125.16) / 1125.16 * 100
    log(f"  GL(1) vs GL(2) 차이: {abs(near_med - 1114.6):.2f} ({diff_gl2:.2f}%)")
    log(f"  GL(1) vs GL(3) 차이: {abs(near_med - 1125.16):.2f} ({diff_gl3:.2f}%)")
    log(f"")

    if near_med < 1114.6:
        log(f"  → κ_near(1) < κ_near(2) < κ_near(3): 단조증가 f(d)")
        if abs(near_med - (1114.6 - 10.3)) < 5:
            log(f"  → 등차 가설: Δκ ≈ 10.3/degree step")
        else:
            gap_12 = 1114.6 - near_med
            gap_23 = 1125.16 - 1114.6
            log(f"  → gap(1→2) = {gap_12:.2f}, gap(2→3) = {gap_23:.2f}")
    elif near_med > 1125.16:
        log(f"  → κ_near(1) > κ_near(3) > κ_near(2): 비단조. 감마 인자 구조 의존?")
    elif abs(near_med - 1114.6) < 3:
        log(f"  → κ_near(1) ≈ κ_near(2): degree 독립? 또는 d≥1에서 포화?")
    else:
        log(f"  → κ_near(1) ∈ (κ_near(2), κ_near(3)): 비단조, 복합 구조")

log(f"")
log(f"총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log(f"")

# TP별 상세 테이블
log("TP별 상세:")
log(f"{'TP':>4} | {'γ':>12} | {'κ_near':>10} | {'mono/π':>8}")
log("-" * 50)
for j in range(len(kappa_near)):
    m_str = f"{mono_results[j]:.4f}" if mono_results[j] is not None else "FAIL"
    log(f"  #{j+1:>2} | {zeros[j]:>12.6f} | {kappa_near[j]:>10.2f} | {m_str:>8}")

log()

# σ-유일성 현황표
log("σ-유일성 현황표 (N=1 ⇔ PASS):")
log(f"  | L-함수 | degree | N | σ-유일성 | n |")
log(f"  |--------|--------|---|----------|---|")
log(f"  | ζ(s)           | 1 | 1 | {'PASS' if sigma_pass_count == sigma_test_count else 'FAIL'} ({sigma_pass_count}/{sigma_test_count}) | {n_zeros} | ← 이번")
log(f"  | Maass R=9.53   | 2 | 1 | PASS (3/3) | 18 |")
log(f"  | Maass R=13.78  | 2 | 1 | PASS (24/24) | 24 |")
log(f"  | Δ (w=12)       | 2 | 1 | PASS (23/23) | 23 |")
log(f"  | 11a1           | 2 | 11 | FAIL | — |")
log(f"  | 37a1           | 2 | 37 | FAIL | — |")
log(f"  | sym²(11a1)     | 3 | 121 | FAIL | — |")

log()
log("[완료]")
flush()
