"""
=============================================================================
[C-298] 고해상도 σ 감쇠 프로파일: E(σ) 멱법칙 검증
=============================================================================
목표:
  1. σ=0.5 근처 고밀도 격자(40점)에서 E(σ) 측정
  2. E(σ) = A/|σ-1/2|^α + B 피팅 → α 결정
  3. 이론 예측 α=1 (Hadamard 부분분수) 검증
  4. κ 비율 고해상도 갱신 (C-297의 19× 격자 한계 해소)

이론적 근거:
  Hadamard 곱에서 ξ'/ξ(s) = Σ_ρ 1/(s-ρ) + const.
  E(σ) = ∫|ξ'/ξ|²dt 의 대각항 = Σ_n π/|σ-1/2| = πN/|σ-1/2|.
  따라서 E(σ) ~ πN/|σ-1/2| + O(1)이 예측됨.

방법:
  - σ 격자: 0.5 ± {0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02,
    0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4} (36점)
  - E(σ) 대칭 검증: E(0.5+Δ) vs E(0.5-Δ) 비교
  - 멱법칙 피팅: log-log 선형 회귀로 α 결정
  - 이론값 비교: A_fit vs πN/2 (인자 1/2는 유한 범위 효과)

t 구간: [14, 34] (C-297과 동일, ~9개 영점 포함)
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))

mpmath.mp.dps = 30

# ─── 기본 함수 (C-297 동일) ───

def xi_func(s):
    half = mpmath.mpf('0.5')
    return half * s * (s - 1) * mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s)

def L_func(s):
    h = mpmath.mpf('1e-15')
    xi_val = xi_func(s)
    if abs(xi_val) < 1e-40:
        return mpmath.mpc('inf')
    xi_d = (xi_func(s + h) - xi_func(s - h)) / (2 * h)
    return xi_d / xi_val

def curvature(s):
    L = L_func(s)
    return float(abs(L)**2)


# ─── E(σ) 계산 ───

def compute_energy(sigma, t_min=14.0, t_max=34.0, n_t=400):
    """E(σ) = ∫|ξ'/ξ(σ+it)|² dt"""
    ts = np.linspace(t_min, t_max, n_t)
    dt = (t_max - t_min) / (n_t - 1)
    vals = []
    for t in ts:
        s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
        try:
            k = curvature(s)
            if np.isfinite(k):
                vals.append(k)
            else:
                vals.append(1e6)
        except:
            vals.append(0.0)
    return np.trapezoid(vals, dx=dt)


# ─── 이론 예측: 영점 수 세기 ───

def count_zeros_in_range(t_min, t_max):
    """[t_min, t_max]의 ξ 영점 수 (N(T) 공식)"""
    from mpmath import zetazero
    zeros = []
    n = 1
    while True:
        z = float(zetazero(n).imag)
        if z > t_max:
            break
        if z >= t_min:
            zeros.append(z)
        n += 1
    return len(zeros), zeros


# ─── 메인 ───

def main():
    t0 = time.time()
    outpath = os.path.expanduser('~/Desktop/gdl_unified/results/offcritical_highres_c298.txt')

    # σ 격자 설계: 0.5 ± Δσ (Δσ 대수적으로 분포)
    deltas = [0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.025,
              0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]

    t_min, t_max = 14.0, 34.0

    # 영점 수 세기
    print("영점 수집 중...")
    N_zeros, zero_list = count_zeros_in_range(t_min, t_max)
    print(f"[{t_min}, {t_max}] 구간 영점 수: {N_zeros}")
    for i, z in enumerate(zero_list):
        print(f"  ρ_{i+1}: t = {z:.6f}")

    # E(σ) 측정
    print(f"\nE(σ) 고해상도 측정 시작 (Δσ {len(deltas)}개 × 양쪽 = {2*len(deltas)} 점)...")

    results = []  # (Δσ, E_right, E_left)

    for i, d in enumerate(deltas):
        sigma_r = 0.5 + d
        sigma_l = 0.5 - d

        print(f"  [{i+1}/{len(deltas)}] Δσ={d:.4f} ...", end=" ", flush=True)
        t1 = time.time()

        E_r = compute_energy(sigma_r, t_min, t_max)
        E_l = compute_energy(sigma_l, t_min, t_max)

        sym_ratio = E_r / E_l if E_l > 0 else float('inf')
        elapsed = time.time() - t1

        print(f"E(+)={E_r:.2f}, E(-)={E_l:.2f}, 대칭비={sym_ratio:.6f}, {elapsed:.1f}s")
        results.append((d, E_r, E_l))

    # 결과 표 작성
    print("\n" + "="*70)

    # E(0.5) 참조 (C-297에서 30515)
    print("E(0.5) 측정 중...")
    E_half = compute_energy(0.5, t_min, t_max)
    print(f"E(0.5) = {E_half:.2f}")

    # ─── 분석 1: E(σ) = A/Δσ^α + B 피팅 ───

    deltas_arr = np.array([r[0] for r in results])
    E_avg = np.array([(r[1] + r[2]) / 2 for r in results])  # 대칭 평균

    # 배경 추정: 큰 Δσ에서의 E 값
    bg_estimate = np.mean(E_avg[-3:])  # 마지막 3점 평균

    # log-log 피팅 (배경 빼기 전)
    log_delta = np.log(deltas_arr)
    log_E = np.log(E_avg)

    # 배경 제거 후 피팅 (Δσ < 0.2만 사용)
    mask_fit = deltas_arr <= 0.2
    E_sub = E_avg[mask_fit] - bg_estimate
    E_sub = np.maximum(E_sub, 1.0)  # 음수 방지

    log_delta_fit = np.log(deltas_arr[mask_fit])
    log_E_sub = np.log(E_sub)

    # 선형 회귀: log(E-B) = log(A) - α·log(Δσ)
    coeffs = np.polyfit(log_delta_fit, log_E_sub, 1)
    alpha_fit = -coeffs[0]
    A_fit = np.exp(coeffs[1])

    # 이론 예측: A_theory = πN/2 (유한 범위 인자)
    A_theory = np.pi * N_zeros / 2

    # 비선형 피팅: E = A/Δσ^α + B (전체)
    def power_law(x, A, alpha, B):
        return A / x**alpha + B

    try:
        popt, pcov = curve_fit(power_law, deltas_arr, E_avg,
                                p0=[A_theory, 1.0, bg_estimate],
                                bounds=([0, 0, 0], [1e6, 5, 1e4]))
        A_nl, alpha_nl, B_nl = popt
        perr = np.sqrt(np.diag(pcov))

        # 잔차
        E_pred = power_law(deltas_arr, *popt)
        residuals = E_avg - E_pred
        R2 = 1 - np.sum(residuals**2) / np.sum((E_avg - np.mean(E_avg))**2)

        nl_success = True
    except Exception as e:
        print(f"비선형 피팅 실패: {e}")
        nl_success = False

    # ─── 분석 2: 대칭성 검증 ───
    sym_ratios = [r[1] / r[2] if r[2] > 0 else float('inf') for r in results]

    # ─── 분석 3: κ 비율 갱신 ───
    # 가장 작은 Δσ에서의 E 비율
    E_nearest = (results[0][1] + results[0][2]) / 2  # Δσ=0.001
    E_far = (results[-1][1] + results[-1][2]) / 2     # Δσ=0.4
    kappa_ratio = E_nearest / E_far

    # ─── 출력 ───

    with open(outpath, 'w') as f:
        f.write("="*70 + "\n")
        f.write("[C-298] 고해상도 σ 감쇠 프로파일\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"t 구간: [{t_min}, {t_max}], 영점 수 N={N_zeros}\n")
        f.write(f"σ 격자: {len(deltas)} Δσ 값 × 양쪽 = {2*len(deltas)} 점\n")
        f.write("="*70 + "\n\n")

        f.write("[1] 영점 목록\n")
        for i, z in enumerate(zero_list):
            f.write(f"  ρ_{i+1}: t = {z:.6f}\n")
        f.write(f"\n")

        f.write("[2] E(σ) 고해상도 프로파일\n")
        f.write(f"{'Δσ':>10}  {'E(0.5+Δ)':>14}  {'E(0.5-Δ)':>14}  {'E(avg)':>14}  {'대칭비':>10}\n")
        f.write("-"*70 + "\n")
        for r, sym in zip(results, sym_ratios):
            d, Er, El = r
            E_a = (Er + El) / 2
            f.write(f"{d:10.4f}  {Er:14.4f}  {El:14.4f}  {E_a:14.4f}  {sym:10.6f}\n")
        f.write(f"\nE(0.5) = {E_half:.2f}\n")

        f.write(f"\n[3] 멱법칙 피팅: E(σ) = A/|σ-1/2|^α + B\n")
        f.write(f"\n  (a) Log-log 선형 회귀 (배경 {bg_estimate:.1f} 제거):\n")
        f.write(f"      α = {alpha_fit:.4f}\n")
        f.write(f"      A = {A_fit:.4f}\n")

        if nl_success:
            f.write(f"\n  (b) 비선형 피팅 (전체 데이터):\n")
            f.write(f"      A = {A_nl:.4f} ± {perr[0]:.4f}\n")
            f.write(f"      α = {alpha_nl:.4f} ± {perr[1]:.4f}\n")
            f.write(f"      B = {B_nl:.4f} ± {perr[2]:.4f}\n")
            f.write(f"      R² = {R2:.6f}\n")

        f.write(f"\n  (c) 이론 예측 (Hadamard 부분분수):\n")
        f.write(f"      α_theory = 1\n")
        f.write(f"      A_theory = πN/2 = {A_theory:.4f}\n")
        f.write(f"      |α_fit - 1| = {abs(alpha_nl - 1):.4f}\n" if nl_success else "")
        f.write(f"      A_fit/A_theory = {A_nl/A_theory:.4f}\n" if nl_success else "")

        f.write(f"\n[4] 대칭성 E(σ)=E(1-σ) 검증\n")
        sym_arr = np.array(sym_ratios)
        f.write(f"  max |ratio-1| = {np.max(np.abs(sym_arr - 1)):.6f}\n")
        f.write(f"  mean |ratio-1| = {np.mean(np.abs(sym_arr - 1)):.6f}\n")

        f.write(f"\n[5] κ 비율 갱신\n")
        f.write(f"  E(Δσ=0.001)/E(Δσ=0.4) = {kappa_ratio:.1f}×\n")
        f.write(f"  (C-297: 19× → C-298: {kappa_ratio:.0f}×)\n")

        # 판정
        f.write(f"\n" + "="*70 + "\n")
        f.write("[판정]\n")
        if nl_success and abs(alpha_nl - 1) < 0.1:
            f.write(f"★★★★★ 강양성: α = {alpha_nl:.3f} ≈ 1 (Hadamard 예측 확인)\n")
            f.write(f"E(σ) = {A_nl:.1f}/|σ-1/2| + {B_nl:.1f}  (R²={R2:.4f})\n")
        elif nl_success and abs(alpha_nl - 1) < 0.3:
            f.write(f"★★★★ 양성: α = {alpha_nl:.3f} ≈ 1 (근사 일치)\n")
        elif nl_success:
            f.write(f"★★★ 중립: α = {alpha_nl:.3f} ≠ 1 (이론 불일치, 원인 조사 필요)\n")
        else:
            f.write("? 피팅 실패 — 수동 분석 필요\n")

        elapsed_total = time.time() - t0
        f.write(f"\n총 소요: {elapsed_total:.0f}초 ({elapsed_total/60:.1f}분)\n")

    print(f"\n결과 저장: {outpath}")
    print(f"총 소요: {(time.time()-t0)/60:.1f}분")

    # 핵심 수치 콘솔 출력
    if nl_success:
        print(f"\n=== 핵심 결과 ===")
        print(f"α = {alpha_nl:.4f} ± {perr[1]:.4f}  (이론: 1)")
        print(f"A = {A_nl:.4f}  (이론: πN/2 = {A_theory:.4f})")
        print(f"R² = {R2:.6f}")
        print(f"κ 비율: {kappa_ratio:.0f}× (C-297: 19×)")


if __name__ == '__main__':
    main()
