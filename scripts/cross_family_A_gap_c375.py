#!/usr/bin/env python3
"""
=============================================================================
[사이클 #375] 교차 가족 A_Λ–gap_min 보편성 검정 (Paper 4 씨앗)
=============================================================================

목적:
  5개 L-함수 가족에서 ρ(A_Λ, gap_min_GUE)의 보편성을 검정한다.

대상:
  1. ζ(s) — degree 1, N=1 (기준선)
  2. χ₋₃ (mod 3) — degree 1, N=3
  3. 11a1 타원곡선 — degree 2, N=11
  4. Ramanujan Δ — degree 2, N=1, weight 12
  5. sym²(11a1) — degree 3, N=121

방법:
  - PARI lfunzeros()로 영점 수집
  - Hadamard zero-sum으로 A_Λ(t₀) 계산
    A_Λ = Im(c₀)² + 2·Re(c₁)
    S₁^Λ = S₁_L - Im[smooth], H₁^Λ = H₁_L + Re[smooth']
  - gap_min = min(gap_left, gap_right)
  - GUE 정규화: gap_min_GUE = gap_min × d̄(t₀)
    d̄(t) = (1/2π) log(N · t^d / (2π)^d) [이론적 밀도]
  - Spearman ρ(A_Λ, gap_min_GUE) 계산
  - 중앙 60% 트리밍

성공 기준:
  1. 5개 L-함수 전부에서 ρ 측정 완료 (n≥50)
  2. 보편/비보편 판정 + ρ값 ± 95% CI
  3. 비보편이면: 어떤 가족에서, 왜 깨지는지 분석

결과: results/cross_family_A_gap_c375.txt
=============================================================================
"""

import sys, os, time, math
import numpy as np
from scipy import stats
import mpmath

sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import cypari2

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

TRIM_FRAC = 0.20  # 양쪽 20% 제외 → 중앙 60%
N_BOOTSTRAP = 2000  # 95% CI 부트스트랩 횟수

RESULT_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/results/cross_family_A_gap_c375.txt'
)
os.makedirs(os.path.dirname(RESULT_PATH), exist_ok=True)

# L-함수 가족 정의
# center: 임계선 실수부 (PARI 영점의 기준)
# mu: 스펙트럼 파라미터 (Γ_R(s+μ_j) 인자들)
# degree: L-함수 차수
# N_cond: 산술적 도체 (conductor)
FAMILIES = [
    {
        'name': 'zeta',
        'label': 'ζ(s) — degree 1, N=1',
        'degree': 1,
        'N_cond': 1,
        'center': 0.5,
        'mu': [0],
        'T_max': 2000,
        'pari_init': [
            "L_z = lfuncreate(1)",
            "Li_z = lfuninit(L_z, [0, 2100])",
            "zv = lfunzeros(Li_z, 2000)",
        ],
        'zeros_var': 'zv',
    },
    {
        'name': 'chi_m3',
        'label': 'L(s,χ₋₃) — degree 1, N=3',
        'degree': 1,
        'N_cond': 3,
        'center': 0.5,
        'mu': [1],  # χ₋₃(-1) = -1, 홀수 캐릭터 → a=1
        'T_max': 500,
        'pari_init': [
            "L_chi = lfuncreate(-3)",
            "Li_chi = lfuninit(L_chi, [0, 520])",
            "zv = lfunzeros(Li_chi, 500)",
        ],
        'zeros_var': 'zv',
    },
    {
        'name': 'ec_11a1',
        'label': 'L(s,E_11a1) — degree 2, N=11',
        'degree': 2,
        'N_cond': 11,
        'center': 1.0,
        'mu': [0, 1],
        'T_max': 200,
        'pari_init': [
            "E_ec = ellinit([0,-1,1,-10,-20])",
            "L_ec = lfuncreate(E_ec)",
            "Li_ec = lfuninit(L_ec, [0, 220])",
            "zv = lfunzeros(Li_ec, 200)",
        ],
        'zeros_var': 'zv',
    },
    {
        'name': 'ramanujan_delta',
        'label': 'L(s,Δ) — degree 2, N=1, weight 12',
        'degree': 2,
        'N_cond': 1,
        'center': 6.0,  # weight 12 → center = k/2 = 6
        'mu': [0, 1],   # Γ(s) = Γ_R(s)Γ_R(s+1)
        'T_max': 200,
        'pari_init': [
            "L_delta = lfunetaquo([1,24])",
            "Li_delta = lfuninit(L_delta, [0, 220])",
            "zv = lfunzeros(Li_delta, 200)",
        ],
        'zeros_var': 'zv',
    },
    {
        'name': 'sym2_11a1',
        'label': 'L(s,sym²(11a1)) — degree 3, N=121',
        'degree': 3,
        'N_cond': 121,
        'center': 1.5,
        'mu': [1, 1, 2],
        'T_max': 100,
        'pari_init': [
            "E_s2 = ellinit([0,-1,1,-10,-20])",
            "L_s2 = lfunsympow(E_s2, 2)",
            "Li_s2 = lfuninit(L_s2, [0, 120])",
            "zv = lfunzeros(Li_s2, 100)",
        ],
        'zeros_var': 'zv',
    },
]

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 로그
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

_log_buf = []

def log(msg=''):
    print(msg, flush=True)
    _log_buf.append(str(msg))

def save():
    with open(RESULT_PATH, 'w') as f:
        f.write('\n'.join(_log_buf))

sig = lambda p: '✅' if p < 0.01 else ('⚠️' if p < 0.05 else '❌')


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 이론적 영점 밀도 (d̄)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def theoretical_density(t, degree, N_cond):
    """
    이론적 영점 밀도 d̄(t) = (1/2π) log(N · t^d / (2π)^d).

    d = degree, N = conductor.
    경험적 2/(Δt) 방식 사용 금지 — MEMORY.md 경고.
    """
    if t <= 0:
        return 0.0
    val = N_cond * (t / (2.0 * math.pi)) ** degree
    if val <= 0:
        return 0.0
    return math.log(val) / (2.0 * math.pi)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Gamma 보정 (일반)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def gamma_correction(gamma_0, center, mu_list):
    """
    일반 L-함수에 대한 Gamma 보정.

    smooth(ρ₀) = log(N)/2 + Σ_j [-log(π)/2 + ψ((ρ₀+μ_j)/2)/2]
    - Im(smooth) → S₁ 보정
    - Re(smooth') → H₁ 보정

    반환: (im_smooth, re_smooth_deriv)
    """
    s = mpmath.mpc(center, gamma_0)

    im_smooth = 0.0
    re_smooth_deriv = 0.0

    for mu in mu_list:
        arg = (s + mu) / 2
        # S₁ 보정: Im[ψ(arg)/2]
        psi_val = mpmath.digamma(arg)
        im_smooth += float(mpmath.im(psi_val)) / 2.0
        # H₁ 보정: Re[ψ₁(arg)/4]
        psi1_val = mpmath.psi(1, arg)
        re_smooth_deriv += float(mpmath.re(psi1_val)) / 4.0

    return im_smooth, re_smooth_deriv


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# A_Λ 계산 (Hadamard zero-sum)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_A_from_zeros(zeros, idx, center, mu_list):
    """
    영점 위치에서 직접 A_Λ(γ₀) 계산.

    S₁_L = Σ_{k≠i} 1/(γ₀-γ_k) + Σ_k 1/(γ₀+γ_k)  [같은 부호 + 켤레]
    H₁_L = Σ_{k≠i} 1/(γ₀-γ_k)² + Σ_k 1/(γ₀+γ_k)²

    c₀_zeros → Im = -S₁_L, c₁_zeros → Re = +H₁_L
    c₀ = c₀_zeros + smooth → Im(c₀) = -S₁_L + Im(smooth)
    → S₁_Λ = S₁_L - Im(smooth)  [부호 반전 후 제곱하므로 동치]

    c₁ = -Σ 1/(ρ₀-ρ)² + smooth'
    c₁_zeros = +H₁_L  [부호 유도: 1/(ρ₀-ρ)² = -1/Δγ², c₁ = -Σ(-1/Δγ²) = +H₁_L]
    Re(c₁) = H₁_L + Re(smooth')
    → H₁_Λ = H₁_L + Re(smooth')

    A_Λ = S₁_Λ² + 2·H₁_Λ
    """
    g0 = zeros[idx]
    n = len(zeros)

    # 같은 부호 영점 합산
    S1_same = 0.0
    H1_same = 0.0
    for k in range(n):
        if k == idx:
            continue
        diff = g0 - zeros[k]
        if abs(diff) < 1e-15:
            continue
        S1_same += 1.0 / diff
        H1_same += 1.0 / (diff * diff)

    # 켤레 영점 합산: ρ̄_k = center - iγ_k → ρ₀ - ρ̄_k = i(γ₀+γ_k)
    S1_conj = 0.0
    H1_conj = 0.0
    for k in range(n):
        denom = g0 + zeros[k]
        if abs(denom) < 1e-15:
            continue
        S1_conj += 1.0 / denom
        H1_conj += 1.0 / (denom * denom)

    S1_L = S1_same + S1_conj
    H1_L = H1_same + H1_conj

    # Gamma 보정
    im_smooth, re_smooth_deriv = gamma_correction(g0, center, mu_list)

    # A_Λ
    S1_Lambda = S1_L - im_smooth
    H1_Lambda = H1_L + re_smooth_deriv
    A_Lambda = S1_Lambda ** 2 + 2.0 * H1_Lambda

    return A_Lambda


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Bootstrap 95% CI for Spearman ρ
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def bootstrap_ci(x, y, n_bootstrap=N_BOOTSTRAP, ci=0.95):
    """Spearman ρ의 부트스트랩 95% CI."""
    n = len(x)
    rng = np.random.default_rng(42)
    rhos = []
    for _ in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        r, _ = stats.spearmanr(x[idx], y[idx])
        if not np.isnan(r):
            rhos.append(r)
    if len(rhos) < 100:
        return np.nan, np.nan
    rhos = np.sort(rhos)
    alpha = (1 - ci) / 2
    lo = rhos[int(alpha * len(rhos))]
    hi = rhos[int((1 - alpha) * len(rhos))]
    return lo, hi


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 단일 가족 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def analyze_family(pari, fam, t_global_start):
    """단일 L-함수 가족에 대한 A_Λ–gap_min 상관 분석."""
    name = fam['name']
    label = fam['label']
    degree = fam['degree']
    N_cond = fam['N_cond']
    center = fam['center']
    mu_list = fam['mu']
    T_max = fam['T_max']

    log(f"\n{'='*75}")
    log(f"  {label}")
    log(f"{'='*75}")
    t0 = time.time()

    # ──────────────────────────────────────────────────
    # STEP 1: 영점 수집 (PARI)
    # ──────────────────────────────────────────────────
    log(f"[{name}] STEP 1: 영점 수집 (T={T_max})...")

    try:
        for cmd in fam['pari_init']:
            pari(cmd)
    except Exception as e:
        log(f"  ❌ PARI 초기화 실패: {e}")
        return None

    zv = fam['zeros_var']
    try:
        n_z = int(str(pari(f'#{zv}')))
    except Exception as e:
        log(f"  ❌ 영점 개수 조회 실패: {e}")
        return None

    zeros = []
    for i in range(1, n_z + 1):
        try:
            t_str = str(pari(f'{zv}[{i}]')).strip().replace(' E', 'e')
            t = float(t_str)
            if t > 2.0:  # 매우 작은 영점 제외
                zeros.append(t)
        except Exception as e:
            if len(zeros) == 0:
                log(f"  WARNING: 영점 파싱 실패 i={i}: {e}")

    zeros = sorted(zeros)
    if len(zeros) == 0:
        log(f"  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
        return None

    log(f"  {len(zeros)}개 영점, t ∈ [{zeros[0]:.3f}, {zeros[-1]:.3f}]")
    log(f"  소요: {time.time()-t0:.1f}초")

    if len(zeros) < 20:
        log(f"  ⚠️ 영점 {len(zeros)}개 < 20 — 분석 중단")
        return None

    # ──────────────────────────────────────────────────
    # STEP 2: A_Λ 계산 (zero-sum + Gamma)
    # ──────────────────────────────────────────────────
    log(f"[{name}] STEP 2: A_Λ 계산 ({len(zeros)}개 영점)...")

    mpmath.mp.dps = 30

    data = []
    n_fail = 0
    n_neg = 0

    for i in range(len(zeros)):
        try:
            A = compute_A_from_zeros(zeros, i, center, mu_list)
            if np.isnan(A) or np.isinf(A):
                n_fail += 1
                continue
            if A <= 0:
                n_neg += 1
                continue
            data.append({'t': zeros[i], 'A': A, 'idx': i})
        except Exception as e:
            n_fail += 1
            if n_fail <= 3:
                log(f"  WARNING: i={i}, t={zeros[i]:.3f}: {e}")

        if (i + 1) % 200 == 0:
            log(f"  {i+1}/{len(zeros)} 완료 (유효={len(data)}, "
                f"음수={n_neg}, 실패={n_fail}, "
                f"{time.time()-t0:.0f}s)")

    log(f"  완료: 유효={len(data)}, 음수={n_neg}, 실패={n_fail}")

    if len(data) < 20:
        log(f"  ⚠️ 유효 데이터 {len(data)}개 < 20 — 분석 중단")
        return None

    if n_fail > len(zeros) // 2:
        log(f"  ⚠️ 실패율 {n_fail}/{len(zeros)} > 50% — 신뢰도 낮음")

    # ──────────────────────────────────────────────────
    # STEP 3: 중앙 60% + gap 계산 (이론적 밀도)
    # ──────────────────────────────────────────────────
    log(f"[{name}] STEP 3: 중앙 60% 선택 + GUE 정규화...")

    N = len(data)
    lo = int(N * TRIM_FRAC)
    hi = int(N * (1.0 - TRIM_FRAC))

    valid = []
    for idx in range(lo, hi):
        if idx <= 0 or idx >= N - 1:
            continue
        d = data[idx]
        t_n = d['t']
        t_prev = data[idx - 1]['t']
        t_next = data[idx + 1]['t']
        gap_r = t_next - t_n
        gap_l = t_n - t_prev
        if gap_r <= 0 or gap_l <= 0:
            continue

        # 이론적 밀도 (경험적 2/(Δt) 사용 금지!)
        d_bar = theoretical_density(t_n, degree, N_cond)
        if d_bar <= 0:
            continue

        gap_min = min(gap_r, gap_l)
        d['gap_min_gue'] = gap_min * d_bar
        d['gap_r_gue'] = gap_r * d_bar
        d['d_bar'] = d_bar
        valid.append(d)

    log(f"  내부 영점: {len(valid)} (중앙 60%)")

    if len(valid) < 15:
        log(f"  ⚠️ 내부 영점 {len(valid)}개 < 15 — 분석 중단")
        return None

    # 배열 추출
    t_arr = np.array([d['t'] for d in valid])
    A_arr = np.array([d['A'] for d in valid])
    gm_arr = np.array([d['gap_min_gue'] for d in valid])
    gr_arr = np.array([d['gap_r_gue'] for d in valid])

    # ──────────────────────────────────────────────────
    # STEP 4: Spearman 상관 분석
    # ──────────────────────────────────────────────────
    log(f"[{name}] STEP 4: Spearman 상관 (n={len(valid)})...")
    log()

    rho_gm, p_gm = stats.spearmanr(A_arr, gm_arr)
    rho_gr, p_gr = stats.spearmanr(A_arr, gr_arr)

    # Bootstrap 95% CI
    ci_lo_gm, ci_hi_gm = bootstrap_ci(A_arr, gm_arr)
    ci_lo_gr, ci_hi_gr = bootstrap_ci(A_arr, gr_arr)

    log(f"  ρ(A_Λ, gap_min_GUE)   = {rho_gm:+.4f}  "
        f"(p={p_gm:.3e})  95%CI=[{ci_lo_gm:+.4f}, {ci_hi_gm:+.4f}]  {sig(p_gm)}")
    log(f"  ρ(A_Λ, gap_right_GUE) = {rho_gr:+.4f}  "
        f"(p={p_gr:.3e})  95%CI=[{ci_lo_gr:+.4f}, {ci_hi_gr:+.4f}]  {sig(p_gr)}")
    log()

    # 기본 통계
    log(f"  A_Λ: mean={np.mean(A_arr):.2f}, median={np.median(A_arr):.2f}, "
        f"std={np.std(A_arr):.2f}")
    log(f"  gap_min_GUE: mean={np.mean(gm_arr):.4f}, median={np.median(gm_arr):.4f}")
    log(f"  gap_right_GUE: mean={np.mean(gr_arr):.4f}, median={np.median(gr_arr):.4f}")
    log(f"  d̄ range: [{min(d['d_bar'] for d in valid):.4f}, "
        f"{max(d['d_bar'] for d in valid):.4f}]")
    log(f"  t range: [{t_arr.min():.2f}, {t_arr.max():.2f}]")
    log()

    elapsed = time.time() - t0
    log(f"  소요: {elapsed:.1f}초")

    return {
        'name': name,
        'label': label,
        'degree': degree,
        'N_cond': N_cond,
        'n_zeros': len(zeros),
        'n_valid': len(valid),
        'rho_gm': float(rho_gm),
        'p_gm': float(p_gm),
        'ci_lo_gm': float(ci_lo_gm),
        'ci_hi_gm': float(ci_hi_gm),
        'rho_gr': float(rho_gr),
        'p_gr': float(p_gr),
        'ci_lo_gr': float(ci_lo_gr),
        'ci_hi_gr': float(ci_hi_gr),
        'mean_A': float(np.mean(A_arr)),
        'median_A': float(np.median(A_arr)),
        'mean_gm': float(np.mean(gm_arr)),
        't_range': (float(t_arr.min()), float(t_arr.max())),
        'elapsed': elapsed,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_global = time.time()

    log("=" * 80)
    log("[C-375] 교차 가족 A_Λ–gap_min 보편성 검정")
    log("=" * 80)
    log(f"  가족 수: {len(FAMILIES)}")
    log(f"  트리밍: 양쪽 {TRIM_FRAC*100:.0f}% 제외 (중앙 {(1-2*TRIM_FRAC)*100:.0f}%)")
    log(f"  밀도: 이론적 d̄(t) = (1/2π) log(N·t^d/(2π)^d)")
    log(f"  부트스트랩: {N_BOOTSTRAP}회")
    log()

    # PARI 초기화
    pari = cypari2.Pari()
    pari.allocatemem(512 * 10**6)
    pari.set_real_precision(38)
    log("PARI 초기화 완료 (512 MB)")
    log()

    results = []
    for fam in FAMILIES:
        try:
            r = analyze_family(pari, fam, t_global)
            if r is not None:
                results.append(r)
            else:
                log(f"\n  ⚠️ {fam['name']}: 분석 실패")
        except Exception as e:
            log(f"\n  ❌ {fam['name']}: 예외 — {e}")
            import traceback
            traceback.print_exc()
        save()

    # ══════════════════════════════════════════════════════════════════════
    # 최종 비교 + 보편성 판정
    # ══════════════════════════════════════════════════════════════════════

    log()
    log("=" * 80)
    log("[C-375 최종 비교] 교차 가족 보편성 판정")
    log("=" * 80)
    log()

    if len(results) == 0:
        log("❌ 유효한 결과 없음 — 전체 실패")
        save()
        return

    # 비교표
    log(f"  {'가족':<28} {'d':>3} {'N':>6} {'n':>5} "
        f"{'ρ(A_Λ,gm)':>11} {'95% CI':>20} {'p':>10}")
    log("  " + "-" * 95)

    rho_list = []
    for r in results:
        rho_list.append(r['rho_gm'])
        log(f"  {r['label']:<28} {r['degree']:>3} {r['N_cond']:>6} {r['n_valid']:>5} "
            f"{r['rho_gm']:>+11.4f} [{r['ci_lo_gm']:+.4f}, {r['ci_hi_gm']:+.4f}] "
            f"{r['p_gm']:>10.3e}")

    log()
    log("  --- gap_right 참조 ---")
    for r in results:
        log(f"  {r['label']:<28} "
            f"ρ(A_Λ,gr)={r['rho_gr']:+.4f} "
            f"[{r['ci_lo_gr']:+.4f}, {r['ci_hi_gr']:+.4f}]  "
            f"p={r['p_gr']:.3e}")

    log()

    # 보편성 판정
    rho_arr = np.array(rho_list)
    n_strong = sum(1 for r in rho_list if r < -0.70)
    n_moderate = sum(1 for r in rho_list if r < -0.50)
    n_total = len(rho_list)
    rho_range = float(np.max(rho_arr) - np.min(rho_arr))
    rho_mean = float(np.mean(rho_arr))
    rho_std = float(np.std(rho_arr))

    log("[보편성 판정]")
    log()
    log(f"  ρ 전체: mean={rho_mean:+.4f}, std={rho_std:.4f}, "
        f"range={rho_range:.4f}")
    log(f"  ρ < -0.70: {n_strong}/{n_total}")
    log(f"  ρ < -0.50: {n_moderate}/{n_total}")
    log(f"  max-min: {rho_range:.4f} (강한 보편성 기준: < 0.20)")
    log()

    # 판정 분기
    if n_strong == n_total and rho_range < 0.20:
        verdict = "★★★★★ 강한 보편성"
        log(f"  {verdict}: 5개 가족 전부 ρ < -0.70 + 변동 < 0.20")
        log("  → Paper 4 본격 착수. A-gap 상관은 L-함수 보편 법칙.")
    elif n_strong >= 3 and rho_range < 0.30:
        verdict = "★★★★ 보편성 양성"
        log(f"  {verdict}: {n_strong}/5 강한 + 변동 < 0.30")
        log("  → Paper 4 진행. 일부 가족에서 약화 원인 탐색 필요.")
    elif n_moderate >= 3:
        verdict = "★★★ 조건부 양성"
        log(f"  {verdict}: {n_moderate}/5 중등 (ρ < -0.50)")
        # 깨지는 가족 식별
        weak = [r for r in results if r['rho_gm'] > -0.50]
        if weak:
            log(f"  깨지는 가족:")
            for w in weak:
                log(f"    {w['label']}: ρ={w['rho_gm']:+.4f} (n={w['n_valid']})")
            log("  → degree/conductor 의존성 탐색 필요.")
    elif n_moderate >= 1:
        verdict = "★★ 약한 양성"
        log(f"  {verdict}: {n_moderate}/5 중등. 보편성 부분적.")
    else:
        verdict = "★ 음성"
        log(f"  {verdict}: 과반 ρ > -0.50. 보편성 가설 기각.")

    log()

    # degree 별 패턴
    log("[degree 별 패턴]")
    log()
    for d_val in sorted(set(r['degree'] for r in results)):
        family_d = [r for r in results if r['degree'] == d_val]
        rhos = [r['rho_gm'] for r in family_d]
        names = [r['name'] for r in family_d]
        log(f"  degree {d_val}: ρ = {', '.join(f'{r:+.4f}' for r in rhos)} "
            f"({', '.join(names)})")
    log()

    # conductor 의존성
    log("[conductor 의존성]")
    log()
    Ns = np.array([r['N_cond'] for r in results], dtype=float)
    rhos_for_N = np.array([r['rho_gm'] for r in results])
    if len(Ns) >= 4:
        tau_N, p_N = stats.kendalltau(np.log(Ns + 1), rhos_for_N)
        log(f"  Kendall τ(log(N), ρ): τ={tau_N:+.4f}  (p={p_N:.3e})")
    for r in sorted(results, key=lambda x: x['N_cond']):
        log(f"  N={r['N_cond']:>5}: ρ={r['rho_gm']:+.4f}  ({r['name']})")
    log()

    # 종합 요약
    log("=" * 80)
    log(f"[C-375 최종 판정]: {verdict}")
    log("=" * 80)
    log()
    log(f"  가족 수: {n_total}/5 분석 완료")
    log(f"  ρ(A_Λ, gap_min_GUE) 평균: {rho_mean:+.4f} ± {rho_std:.4f}")
    log(f"  ρ 범위: [{np.min(rho_arr):+.4f}, {np.max(rho_arr):+.4f}]")
    log(f"  보편성: {verdict}")
    log()

    # 결과별 분기 명시
    if n_strong == n_total:
        log("  [다음 단계]: Paper 4 본격 착수")
        log("    - t-안정성 교차 검증 (각 가족별)")
        log("    - GUE 이론과의 정밀 비교")
    elif n_moderate >= 3:
        log("  [다음 단계]: 깨지는 가족 해부")
        log("    - degree 의존성 vs conductor 의존성 분리")
        log("    - 영점 수 증대 (통계력 확보)")
    else:
        log("  [다음 단계]: A-gap 보편성 가설 재검토")
        log("    - gap_right vs gap_min 비교")
        log("    - t-범위 의존성 확인")

    log()
    log(f"  총 소요시간: {time.time()-t_global:.1f}초 ({(time.time()-t_global)/60:.1f}분)")
    log()

    save()
    log(f"[완료] 결과: {RESULT_PATH}")


if __name__ == '__main__':
    main()
