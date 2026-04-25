#!/usr/bin/env python3
"""
[C-293] GUE 이론적 정규화 확인 — 기준점 재설정
  목적: C-292에서 L-함수 이론적 d_bar → ρ↑ 확인.
  GUE도 이론적(스무딩) d_bar 적용 시 ρ가 달라지는가?
  GUE는 unfolded → 이론적 d_bar≈1≈경험적 d_bar → ρ 변화 없을 것으로 예상.

  방법:
    N=1000 GUE 행렬, 30 앙상블
    trim 20%, N_MAX=300
    경험적 d_bar vs 스무딩(W=100) d_bar
    C-283 로컬 결과 (-0.967) 재현 확인

결과: results/gue_theoretical_norm_c293.txt
"""

import sys, os, math, time
import numpy as np
from scipy import stats

np.random.seed(42)

N_MAT = 1000    # 행렬 크기
N_ENS = 30      # 앙상블 수
N_MAX = 300     # A_bare 이웃 수
TRIM_FRAC = 0.20
W_SMOOTH = 100  # 스무딩 윈도우

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/gue_theoretical_norm_c293.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)
out_f = open(RESULT_PATH, 'w')


def log(msg=''):
    print(msg, flush=True)
    out_f.write(str(msg) + '\n')
    out_f.flush()


def generate_gue(N, rng):
    """N×N GUE 행렬 (unfolded via Wigner semicircle)"""
    A = rng.standard_normal((N, N))
    H = (A + A.T) / (2 * math.sqrt(2 * N))
    evals = np.linalg.eigvalsh(H)
    # Wigner semicircle 보정 (unfolding)
    evals_sorted = np.sort(evals)
    # CDF: 적분 된 반원 분포 → 균일 분포로 변환
    # CDF of semicircle: F(x) = (1/2) + (x*sqrt(1-x^2) + arcsin(x)) / pi, x ∈ [-1,1]
    # map eigenvalues in [-1,1] to [0,N] uniform
    x = np.clip(evals_sorted, -1.0 + 1e-10, 1.0 - 1e-10)
    F = 0.5 + (x * np.sqrt(1 - x**2) + np.arcsin(x)) / np.pi
    unfolded = F * N  # uniform in [0, N]
    return unfolded


def compute_A_bare(evs, idx, n_max=N_MAX):
    n = len(evs)
    gamma_0 = evs[idx]
    S1 = 0.0
    H1 = 0.0
    for k in range(max(0, idx - n_max), min(n, idx + n_max + 1)):
        if k == idx:
            continue
        dg = gamma_0 - evs[k]
        if abs(dg) < 1e-15:
            continue
        S1 += 1.0 / dg
        H1 += 1.0 / (dg * dg)
    return S1 * S1 + 2.0 * H1


def smooth_density(evs, idx, W=W_SMOOTH):
    n = len(evs)
    lo = max(0, idx - W // 2)
    hi = min(n - 1, idx + W // 2)
    if hi - lo < 4 or abs(evs[hi] - evs[lo]) < 1e-10:
        return float('nan')
    return (hi - lo) / (evs[hi] - evs[lo])


log("=" * 60)
log("  [C-293] GUE 이론적 정규화 확인")
log("=" * 60)
log(f"  N_MAT={N_MAT}, N_ENS={N_ENS}, N_MAX={N_MAX}, trim={TRIM_FRAC*100:.0f}%")
log(f"  W_SMOOTH={W_SMOOTH}")
log(f"  참조: C-283 GUE_local ρ=-0.967")
log()

rng = np.random.default_rng(42)
t_start = time.time()

rho_smooth_list = []
rho_emp_list = []

for ens_idx in range(N_ENS):
    evs = generate_gue(N_MAT, rng)
    n = len(evs)

    # A_bare 계산
    data = []
    for i in range(n):
        A = compute_A_bare(evs, i)
        if math.isnan(A) or math.isinf(A) or A <= 0:
            continue
        data.append({'idx': i, 'A_bare': A})

    n_data = len(data)
    lo = int(n_data * TRIM_FRAC)
    hi = int(n_data * (1.0 - TRIM_FRAC))
    trimmed = data[lo:hi]

    valid_s = []
    valid_e = []
    for d in trimmed:
        idx = d['idx']
        if idx <= 0 or idx >= n - 1:
            continue
        t_n    = evs[idx]
        t_prev = evs[idx - 1]
        t_next = evs[idx + 1]
        gap_r = t_next - t_n
        gap_l = t_n - t_prev
        if gap_r <= 0 or gap_l <= 0:
            continue
        gap_min = min(gap_r, gap_l)
        d_bar_s = smooth_density(evs, idx, W=W_SMOOTH)
        d_bar_e = 2.0 / (t_next - t_prev)
        if math.isnan(d_bar_s) or d_bar_s <= 0:
            continue
        valid_s.append({'A_bare': d['A_bare'], 'gap': gap_min * d_bar_s})
        valid_e.append({'A_bare': d['A_bare'], 'gap': gap_min * d_bar_e})

    if len(valid_s) < 50:
        continue

    A_s = np.array([d['A_bare'] for d in valid_s])
    g_s = np.array([d['gap'] for d in valid_s])
    A_e = np.array([d['A_bare'] for d in valid_e])
    g_e = np.array([d['gap'] for d in valid_e])

    mask_s = np.isfinite(A_s) & np.isfinite(g_s) & (A_s > 0)
    mask_e = np.isfinite(A_e) & np.isfinite(g_e) & (A_e > 0)

    if np.sum(mask_s) < 20:
        continue

    rho_s, _ = stats.spearmanr(A_s[mask_s], g_s[mask_s])
    rho_e, _ = stats.spearmanr(A_e[mask_e], g_e[mask_e])
    rho_smooth_list.append(rho_s)
    rho_emp_list.append(rho_e)

    if (ens_idx + 1) % 10 == 0:
        log(f"  앙상블 {ens_idx+1}/{N_ENS}: ρ_smooth={rho_s:+.4f}, ρ_emp={rho_e:+.4f}  "
            f"({time.time()-t_start:.0f}s)")

log()
log("=" * 60)
log("  [결과]")
log("=" * 60)

if rho_smooth_list:
    arr_s = np.array(rho_smooth_list)
    arr_e = np.array(rho_emp_list)

    log(f"  앙상블 수: {len(arr_s)}")
    log(f"  ρ_smooth: {np.mean(arr_s):+.6f} ± {np.std(arr_s):.6f}  (95%CI: [{np.percentile(arr_s,2.5):.4f}, {np.percentile(arr_s,97.5):.4f}])")
    log(f"  ρ_emp:    {np.mean(arr_e):+.6f} ± {np.std(arr_e):.6f}  (95%CI: [{np.percentile(arr_e,2.5):.4f}, {np.percentile(arr_e,97.5):.4f}])")
    log(f"  Δρ (smooth-emp): {np.mean(arr_s)-np.mean(arr_e):+.4f}")
    log()
    log(f"  C-283 참조 GUE_local: ρ=-0.967")
    log(f"  재현 여부 (ρ_emp vs -0.967): Δ={abs(np.mean(arr_e))-0.967:+.4f}")
    log()

    delta_s = np.mean(arr_s) - np.mean(arr_e)
    if abs(delta_s) < 0.01:
        log("  ★★★★★ GUE 이론적=경험적 (차이<0.01): 예상대로 d_bar≈1이라 동일")
        log("  → GUE 기준점 유지: ρ_GUE ≈ -0.967 (d_bar 방식 무관)")
    elif abs(delta_s) < 0.05:
        log("  ★★★ GUE 미세 차이 (<0.05): 스무딩 효과 있으나 작음")
    else:
        log(f"  ⚠️ GUE 유의미한 차이 ({delta_s:.3f}): 예상 외")

    # 새로운 arithmetic damping (이론적 정규화 기준)
    rho_GUE_theo = float(np.mean(arr_s))
    rho_L_theo   = -0.931  # GL(1)~GL(2) 이론적 평균 (C-291)
    delta_arith  = abs(rho_GUE_theo) - abs(rho_L_theo)
    log()
    log(f"  [신 산술 감쇠 (이론적 정규화, degree-무관)]")
    log(f"  ρ_GUE_theo = {rho_GUE_theo:+.4f}")
    log(f"  ρ_L_theo   = {rho_L_theo:+.4f}  (GL(1)~GL(2) 이론적 평균)")
    log(f"  δ_arith    = {delta_arith:+.4f}  (= Paper 4의 d=1 δ≈0.038과 일치 여부)")
    if 0.02 < delta_arith < 0.10:
        log(f"  → 새로운 해석: d-무관 산술 감쇠 δ≈{delta_arith:.3f}")
        log(f"     'degree-independent arithmetic damping' → Paper 4 핵심 결론 유지 가능")
    else:
        log(f"  → 예상 범위 밖: 추가 분석 필요")
else:
    log("  결과 없음")

elapsed = time.time() - t_start
log()
log(f"  총 소요: {elapsed:.1f}s")
out_f.close()
log(f"\n결과 저장: {RESULT_PATH}")
log("[C-293 완료]")
