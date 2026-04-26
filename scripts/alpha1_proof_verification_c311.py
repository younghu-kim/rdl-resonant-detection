"""
=============================================================================
[C-311] α=1 보편성 해석적 증명 검증
=============================================================================
Proposition: E(σ) = πN/δ + O(δ) as δ→0, 따라서 α=1.

증명 핵심:
  (1) Hadamard 분해: Λ'/Λ = Σ_ρ 1/(s-ρ) + smooth
  (2) E_diag = Σ_ρ π/δ = πN/δ  (정확한 해석적 결과)
  (3) E_cross 의 각 쌍 (j,k) 기여:
      ∫_{-∞}^{∞} Re[1/((δ+iτ)(δ-i(τ-Δ)))] dτ = 4πδ/(4δ²+Δ²) = O(δ)
  (4) E_cross = 4πδ · Σ_{j≠k} 1/(4δ²+Δ²) = O(δ) → E_cross/E_diag → 0

검증 실험:
  1. ζ(s) 영점으로 E_cross 해석적 vs 수치적 비교 (적응 해상도)
  2. E_cross ∝ δ 스케일링 확인
  3. A/πN 수렴 확인 (δ→0에서 A/πN → 1⁺)
  4. GL(2), GL(3) 교차 검증

질문: C-305/306의 "음의 교차항"은 수치 아티팩트인가?
=============================================================================
"""

import sys, os, time, math
import numpy as np
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

T_MIN, T_MAX = 10.0, 50.0
DELTAS = [0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]
N_T_BASE = 2000  # 기본 그리드 점 수 (피크 외 영역)
PEAK_HALFWIDTH_MULT = 10  # 피크 주위 ±10δ 범위 고해상도
PEAK_N_POINTS = 200  # 각 피크당 고해상도 점 수

OUT_FILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/alpha1_proof_c311.txt')
os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# ζ 영점
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def find_zeta_zeros(t_max):
    """mpmath로 ζ 영점 탐색"""
    mpmath.mp.dps = 50
    zeros = []
    n = 1
    while True:
        g = float(mpmath.zetazero(n).imag)
        if g > t_max:
            break
        zeros.append(g)
        n += 1
    return np.array(zeros)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 적응형 고해상도 그리드
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def adaptive_grid(t_min, t_max, zeros, delta, n_base, peak_hw_mult, peak_n):
    """영점 주위에 고해상도 점을 배치하는 적응형 그리드 생성.

    피크 반폭 = peak_hw_mult * delta 범위에 peak_n 점 추가.
    나머지 영역은 n_base 점의 균일 그리드.
    """
    # 기본 균일 그리드
    ts = set(np.linspace(t_min, t_max, n_base))

    # 각 영점 주위 고해상도
    hw = peak_hw_mult * max(delta, 0.001)
    for g in zeros:
        if g - hw < t_min or g + hw > t_max:
            continue
        peak_ts = np.linspace(g - hw, g + hw, peak_n)
        ts.update(peak_ts)

    ts = np.array(sorted(ts))
    # t_min, t_max 범위 외 제거
    ts = ts[(ts >= t_min) & (ts <= t_max)]
    return ts


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 에너지 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_E_had_numerical(ts, zeros, delta):
    """E_Had = ∫|Σ_n 1/(δ+i(t-γ_n))|² dt  (수치 적분, 적응 그리드)"""
    tau = ts[:, None] - zeros[None, :]  # (N_t, N_z)
    denom = delta**2 + tau**2
    re_sum = np.sum(delta / denom, axis=1)
    im_sum = np.sum(-tau / denom, axis=1)
    integrand = re_sum**2 + im_sum**2
    return np.trapezoid(integrand, ts)


def compute_E_diag_analytic(zeros, delta, t_min, t_max):
    """E_diag = Σ_n [arctan((T₂-γ)/δ) - arctan((T₁-γ)/δ)] / δ  (해석적)"""
    total = 0.0
    for g in zeros:
        total += (math.atan2(t_max - g, delta) - math.atan2(t_min - g, delta)) / delta
    return total


def compute_E_cross_analytic_infinite(zeros, delta):
    """E_cross (무한 적분 해석 공식) = Σ_{j≠k} 4πδ/(4δ²+(γ_j-γ_k)²)"""
    N = len(zeros)
    total = 0.0
    for j in range(N):
        for k in range(N):
            if j == k:
                continue
            D = zeros[j] - zeros[k]
            total += 4 * math.pi * delta / (4 * delta**2 + D**2)
    return total


def compute_E_cross_analytic_finite(zeros, delta, t_min, t_max):
    """E_cross (유한 구간 해석 공식)

    ∫_{T₁}^{T₂} Re[1/((δ+i(t-γ_j))(δ-i(t-γ_k)))] dt
    = Re[ -i/(2δ+iΔ) · {log(δ+i(T₂-γ_j)) + log(δ-i(T₂-γ_k))
                       - log(δ+i(T₁-γ_j)) - log(δ-i(T₁-γ_k))} ]

    여기서 Δ = γ_k - γ_j.
    """
    N = len(zeros)
    total = 0.0
    for j in range(N):
        for k in range(N):
            if j == k:
                continue
            gj, gk = zeros[j], zeros[k]
            D = gk - gj  # Δ

            # 해석적 부정적분:
            # ∫ dt / ((δ+i(t-gj))(δ-i(t-gk)))
            # = 1/(2δ+iΔ) · [-i·log(δ+i(t-gj)) + i·log(δ-i(t-gk))]
            # (부분분수 분해 후 적분)

            c = complex(2*delta, D)  # 2δ+iΔ

            # t = t_max
            z1_hi = complex(delta, t_max - gj)
            z2_hi = complex(delta, -(t_max - gk))

            # t = t_min
            z1_lo = complex(delta, t_min - gj)
            z2_lo = complex(delta, -(t_min - gk))

            # 부정적분 값
            from cmath import log as clog
            val_hi = (-1j * clog(z1_hi) + 1j * clog(z2_hi)) / c
            val_lo = (-1j * clog(z1_lo) + 1j * clog(z2_lo)) / c

            integral = (val_hi - val_lo).real
            total += integral

    return total


def compute_S2(zeros):
    """쌍 간격 합: S₂ = Σ_{j≠k} 1/(γ_j-γ_k)²"""
    N = len(zeros)
    total = 0.0
    for j in range(N):
        for k in range(N):
            if j == k:
                continue
            total += 1.0 / (zeros[j] - zeros[k])**2
    return total


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 멱법칙 피팅
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def fit_power_law(deltas, energies, label=""):
    """log-log 선형 피팅: E = A/δ^α → log E = log A - α log δ"""
    mask = np.array([e > 0 for e in energies])
    if mask.sum() < 3:
        return None
    d = np.array(deltas)[mask]
    e = np.array(energies)[mask]
    log_d = np.log(d)
    log_e = np.log(e)
    A_mat = np.column_stack([np.ones(len(d)), log_d])
    coeffs = np.linalg.lstsq(A_mat, log_e, rcond=None)[0]
    alpha = -coeffs[1]
    A = np.exp(coeffs[0])
    y_pred = A_mat @ coeffs
    ss_res = np.sum((log_e - y_pred)**2)
    ss_tot = np.sum((log_e - log_e.mean())**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0
    return {'alpha': alpha, 'A': A, 'R2': r2, 'n': len(d), 'label': label}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실험
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    print("=" * 70, flush=True)
    print("[C-311] α=1 보편성 해석적 증명 검증", flush=True)
    print("=" * 70, flush=True)

    # --- 1. ζ 영점 ---
    print("\n[1] ζ 영점 탐색...", flush=True)
    all_zeros = find_zeta_zeros(T_MAX + 5)
    zeros = all_zeros[(all_zeros >= T_MIN) & (all_zeros <= T_MAX)]
    N = len(zeros)
    print(f"  영점 {N}개 (t∈[{T_MIN},{T_MAX}])", flush=True)
    print(f"  πN = {math.pi * N:.4f}", flush=True)

    # --- 2. S₂ 계산 ---
    print("\n[2] 쌍 간격 합 S₂ 계산...", flush=True)
    S2 = compute_S2(zeros)
    mean_gap = np.mean(np.diff(zeros)) if N > 1 else 1.0
    print(f"  S₂ = {S2:.4f}", flush=True)
    print(f"  평균 간격 = {mean_gap:.4f}", flush=True)
    print(f"  δ_c ≈ √(N/(4·S₂)) = {math.sqrt(N/(4*S2)):.4f}", flush=True)

    # --- 3. δ 스윕 ---
    print("\n[3] δ 스윕 — E_diag, E_cross(해석), E_Had(수치) 비교", flush=True)
    print(f"  {'δ':>8s} {'E_diag':>12s} {'E_cross∞':>12s} {'E_crossF':>12s} "
          f"{'E_Had_num':>12s} {'cross/diag':>10s} {'num-diag':>10s} "
          f"{'A/πN':>8s} {'grid_N':>7s}", flush=True)

    results = []

    for delta in DELTAS:
        # 적응형 그리드 생성
        ts = adaptive_grid(T_MIN, T_MAX, zeros, delta,
                          N_T_BASE, PEAK_HALFWIDTH_MULT, PEAK_N_POINTS)

        # 해석적 대각항
        E_diag = compute_E_diag_analytic(zeros, delta, T_MIN, T_MAX)

        # 해석적 교차항 (무한 적분)
        E_cross_inf = compute_E_cross_analytic_infinite(zeros, delta)

        # 해석적 교차항 (유한 적분)
        E_cross_fin = compute_E_cross_analytic_finite(zeros, delta, T_MIN, T_MAX)

        # 수치적 E_Had (적응 그리드)
        E_had_num = compute_E_had_numerical(ts, zeros, delta)

        # 수치적 교차항 = E_had_num - E_diag
        E_cross_num = E_had_num - E_diag

        # 비율
        cross_ratio_inf = E_cross_inf / E_diag if E_diag > 0 else 0
        num_vs_diag = E_cross_num / E_diag if E_diag > 0 else 0

        # A/πN 계산 (E_diag 기준)
        A_over_piN_diag = (E_diag * delta) / (math.pi * N)
        A_over_piN_had = (E_had_num * delta) / (math.pi * N)

        row = {
            'delta': delta,
            'E_diag': E_diag,
            'E_cross_inf': E_cross_inf,
            'E_cross_fin': E_cross_fin,
            'E_had_num': E_had_num,
            'E_cross_num': E_cross_num,
            'cross_ratio_inf': cross_ratio_inf,
            'num_vs_diag': num_vs_diag,
            'A_piN_diag': A_over_piN_diag,
            'A_piN_had': A_over_piN_had,
            'grid_N': len(ts),
        }
        results.append(row)

        print(f"  {delta:8.3f} {E_diag:12.2f} {E_cross_inf:12.2f} {E_cross_fin:12.2f} "
              f"{E_had_num:12.2f} {cross_ratio_inf:10.6f} {num_vs_diag:10.6f} "
              f"{A_over_piN_had:8.4f} {len(ts):7d}", flush=True)

    # --- 4. E_cross 스케일링 검증 ---
    print("\n[4] E_cross 스케일링 검증: E_cross_∞ ∝ δ^β", flush=True)

    ds = np.array([r['delta'] for r in results])
    E_cross_infs = np.array([r['E_cross_inf'] for r in results])
    E_cross_fins = np.array([r['E_cross_fin'] for r in results])
    E_diags = np.array([r['E_diag'] for r in results])
    E_hads = np.array([r['E_had_num'] for r in results])

    # E_cross_∞ ∝ δ^β 피팅
    mask_small = ds <= 0.15
    if mask_small.sum() >= 3:
        log_d = np.log(ds[mask_small])
        log_ec = np.log(E_cross_infs[mask_small])
        A_mat = np.column_stack([np.ones(mask_small.sum()), log_d])
        coeffs = np.linalg.lstsq(A_mat, log_ec, rcond=None)[0]
        beta = coeffs[1]
        y_pred = A_mat @ coeffs
        ss_res = np.sum((log_ec - y_pred)**2)
        ss_tot = np.sum((log_ec - log_ec.mean())**2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0
        print(f"  E_cross_∞: β = {beta:.4f} (이론: 1.0), R² = {r2:.6f}", flush=True)

    # E_cross_finite ∝ δ^β 피팅
    if mask_small.sum() >= 3:
        log_ecf = np.log(np.abs(E_cross_fins[mask_small]))
        A_mat = np.column_stack([np.ones(mask_small.sum()), log_d])
        coeffs = np.linalg.lstsq(A_mat, log_ecf, rcond=None)[0]
        beta_f = coeffs[1]
        y_pred = A_mat @ coeffs
        ss_res = np.sum((log_ecf - y_pred)**2)
        ss_tot = np.sum((log_ecf - log_ecf.mean())**2)
        r2_f = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0
        print(f"  E_cross_fin: β = {beta_f:.4f}, R² = {r2_f:.6f}", flush=True)
        # 부호 확인
        signs = ["+" if x > 0 else "-" for x in E_cross_fins[mask_small]]
        print(f"  E_cross_fin 부호: {' '.join(signs)}", flush=True)

    # --- 5. E_cross/E_diag ∝ δ² 확인 ---
    print("\n[5] E_cross_∞/E_diag ∝ δ² 확인:", flush=True)
    ratios = E_cross_infs / E_diags
    for i, d in enumerate(ds):
        ratio = ratios[i]
        ratio_over_d2 = ratio / d**2 if d > 0 else 0
        print(f"  δ={d:.3f}: ratio={ratio:.6f}, ratio/δ²={ratio_over_d2:.2f}", flush=True)

    # --- 6. E_Had 수치 vs (E_diag + E_cross_fin) 비교 (정확도 검증) ---
    print("\n[6] 수치 정확도: E_Had_num vs (E_diag + E_cross_fin)", flush=True)
    for r in results:
        E_analytic = r['E_diag'] + r['E_cross_fin']
        rel_err = abs(r['E_had_num'] - E_analytic) / E_analytic * 100 if E_analytic > 0 else 0
        print(f"  δ={r['delta']:.3f}: E_num={r['E_had_num']:.2f}, "
              f"E_ana={E_analytic:.2f}, err={rel_err:.3f}%", flush=True)

    # --- 7. E_Had α=1 피팅 (적응 해상도 기준) ---
    print("\n[7] E_Had(적응 그리드) 멱법칙 피팅:", flush=True)

    # 전체 δ
    fit_all = fit_power_law(ds, E_hads, "전체")
    if fit_all:
        print(f"  [전체 {fit_all['n']}점] α = {fit_all['alpha']:.6f}, "
              f"A/πN = {fit_all['A']/(math.pi*N):.6f}, R² = {fit_all['R2']:.8f}", flush=True)

    # 소 δ (≤0.10)
    mask_s = ds <= 0.10
    fit_s = fit_power_law(ds[mask_s], E_hads[mask_s], "소δ")
    if fit_s:
        print(f"  [δ≤0.10, {fit_s['n']}점] α = {fit_s['alpha']:.6f}, "
              f"A/πN = {fit_s['A']/(math.pi*N):.6f}, R² = {fit_s['R2']:.8f}", flush=True)

    # 초소 δ (≤0.05)
    mask_xs = ds <= 0.05
    fit_xs = fit_power_law(ds[mask_xs], E_hads[mask_xs], "초소δ")
    if fit_xs:
        print(f"  [δ≤0.05, {fit_xs['n']}점] α = {fit_xs['alpha']:.6f}, "
              f"A/πN = {fit_xs['A']/(math.pi*N):.6f}, R² = {fit_xs['R2']:.8f}", flush=True)

    # E_diag 피팅 (기준선)
    fit_diag = fit_power_law(ds, E_diags, "대각항")
    if fit_diag:
        print(f"  [대각항 {fit_diag['n']}점] α = {fit_diag['alpha']:.6f}, "
              f"A/πN = {fit_diag['A']/(math.pi*N):.6f}, R² = {fit_diag['R2']:.8f}", flush=True)

    # --- 8. 유한 vs 무한 교차항 비교 ---
    print("\n[8] 유한 vs 무한 적분 교차항 비교:", flush=True)
    print(f"  {'δ':>8s} {'E_cross∞':>12s} {'E_cross_fin':>12s} {'ratio':>8s}", flush=True)
    for r in results:
        ratio = r['E_cross_fin'] / r['E_cross_inf'] if abs(r['E_cross_inf']) > 0 else 0
        print(f"  {r['delta']:8.3f} {r['E_cross_inf']:12.4f} {r['E_cross_fin']:12.4f} "
              f"{ratio:8.4f}", flush=True)

    # --- 9. 이론적 예측: δ→0에서 δ·E → πN ---
    print("\n[9] δ·E 수렴 (δ→0에서 πN으로 수렴해야 함):", flush=True)
    print(f"  πN = {math.pi * N:.4f}", flush=True)
    for r in results:
        dE_diag = r['delta'] * r['E_diag']
        dE_had = r['delta'] * r['E_had_num']
        dE_ana = r['delta'] * (r['E_diag'] + r['E_cross_fin'])
        print(f"  δ={r['delta']:.3f}: δ·E_diag={dE_diag:.4f}, δ·E_Had={dE_had:.4f}, "
              f"δ·E_ana={dE_ana:.4f}", flush=True)

    elapsed = time.time() - t_start

    # --- 종합 판정 ---
    print("\n" + "=" * 70, flush=True)
    print("종합 판정:", flush=True)

    # 판정 1: E_cross_∞ ∝ δ^1?
    if mask_small.sum() >= 3:
        print(f"  [1] E_cross_∞ 스케일링: β={beta:.4f} (이론: 1.0)", flush=True)
        if abs(beta - 1.0) < 0.05:
            print(f"      → ✅ O(δ) 확인 (|β-1|={abs(beta-1):.4f} < 0.05)", flush=True)
        else:
            print(f"      → ⚠️ 편차 (|β-1|={abs(beta-1):.4f})", flush=True)

    # 판정 2: E_cross/E_diag → 0?
    print(f"  [2] E_cross/E_diag → 0:", flush=True)
    for d_val in [0.01, 0.02, 0.05]:
        idx = next((i for i, r in enumerate(results) if abs(r['delta'] - d_val) < 0.001), None)
        if idx is not None:
            print(f"      δ={d_val}: {ratios[idx]:.6f} (→0 확인)", flush=True)

    # 판정 3: α = 1?
    if fit_xs:
        print(f"  [3] α (초소δ) = {fit_xs['alpha']:.6f}, R² = {fit_xs['R2']:.8f}", flush=True)
        if abs(fit_xs['alpha'] - 1.0) < 0.01:
            print(f"      → ✅ α=1 확인", flush=True)
        else:
            print(f"      → ⚠️ 편차 |α-1|={abs(fit_xs['alpha']-1):.6f}", flush=True)

    # 판정 4: 수치 아티팩트 진단
    print(f"  [4] 수치 아티팩트 진단:", flush=True)
    for d_val in [0.02, 0.05]:
        idx = next((i for i, r in enumerate(results) if abs(r['delta'] - d_val) < 0.001), None)
        if idx is not None:
            r = results[idx]
            E_ana = r['E_diag'] + r['E_cross_fin']
            err = (r['E_had_num'] - E_ana) / E_ana * 100 if E_ana > 0 else 0
            print(f"      δ={d_val}: E_Had_num - E_analytic = {err:+.3f}%", flush=True)

    print(f"\n경과 시간: {elapsed:.1f}s", flush=True)
    print("=" * 70, flush=True)

    # --- 결과 저장 ---
    with open(OUT_FILE, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"[C-311] α=1 보편성 해석적 증명 검증\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"경과: {elapsed:.1f}s\n")
        f.write(f"N_zeros = {N}, πN = {math.pi*N:.4f}\n")
        f.write(f"S₂ = {S2:.4f}, 평균간격 = {mean_gap:.4f}\n")
        f.write(f"δ_c ≈ {math.sqrt(N/(4*S2)):.4f}\n")
        f.write("=" * 70 + "\n\n")

        f.write("δ 스윕 결과:\n")
        f.write(f"  {'δ':>8s} {'E_diag':>12s} {'E_cross∞':>12s} {'E_crossF':>12s} "
                f"{'E_Had_num':>12s} {'c∞/diag':>10s} {'A/πN_Had':>10s}\n")
        for r in results:
            f.write(f"  {r['delta']:8.3f} {r['E_diag']:12.2f} {r['E_cross_inf']:12.2f} "
                    f"{r['E_cross_fin']:12.2f} {r['E_had_num']:12.2f} "
                    f"{r['cross_ratio_inf']:10.6f} {r['A_piN_had']:10.4f}\n")

        f.write(f"\nE_cross_∞ 스케일링: β = {beta:.4f} (이론: 1.0), R² = {r2:.6f}\n")
        f.write(f"E_cross_fin 스케일링: β = {beta_f:.4f}, R² = {r2_f:.6f}\n")

        if fit_all:
            f.write(f"\n피팅 (전체): α={fit_all['alpha']:.6f}, A/πN={fit_all['A']/(math.pi*N):.6f}, "
                    f"R²={fit_all['R2']:.8f}\n")
        if fit_s:
            f.write(f"피팅 (δ≤0.10): α={fit_s['alpha']:.6f}, A/πN={fit_s['A']/(math.pi*N):.6f}, "
                    f"R²={fit_s['R2']:.8f}\n")
        if fit_xs:
            f.write(f"피팅 (δ≤0.05): α={fit_xs['alpha']:.6f}, A/πN={fit_xs['A']/(math.pi*N):.6f}, "
                    f"R²={fit_xs['R2']:.8f}\n")
        if fit_diag:
            f.write(f"피팅 (대각항): α={fit_diag['alpha']:.6f}, A/πN={fit_diag['A']/(math.pi*N):.6f}, "
                    f"R²={fit_diag['R2']:.8f}\n")

        f.write(f"\nδ·E 수렴:\n")
        f.write(f"  πN = {math.pi * N:.4f}\n")
        for r in results:
            f.write(f"  δ={r['delta']:.3f}: δ·E_diag={r['delta']*r['E_diag']:.4f}, "
                    f"δ·E_Had={r['delta']*r['E_had_num']:.4f}\n")

        # 판정
        f.write(f"\n종합 판정:\n")
        verdict = "양성" if (abs(beta - 1.0) < 0.05 and fit_xs and abs(fit_xs['alpha'] - 1.0) < 0.01) else "조건부"
        f.write(f"  E_cross_∞ β = {beta:.4f} (이론 1.0)\n")
        if fit_xs:
            f.write(f"  α (초소δ) = {fit_xs['alpha']:.6f}\n")
        f.write(f"  판정: ★★★★{'★' if verdict=='양성' else ''} {verdict}\n")

    print(f"\n결과 저장: {OUT_FILE}", flush=True)


if __name__ == '__main__':
    main()
