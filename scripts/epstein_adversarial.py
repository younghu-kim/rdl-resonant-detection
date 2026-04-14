#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #29 — Epstein 제타 함수 적대적 검증
=============================================================================
Epstein 제타 Z_Q(s) = Σ' Q(m,n)^{-s} (Q(m,n) = m² + λn²)

RH를 위반하는 함수에서 다발 서명이 σ=1/2 영점과 어떻게 다른지 검증.

λ=1: Z_Q = 4ζ(s)L(s,χ_{-4}) — RH 만족 (대조군)
λ=5: class number h(-20)=2 → 개별 Epstein zeta에 σ≠1/2 영점 가능

분석적 접속: θ 함수 Mellin 적분, t=1 분할
  Λ(s) = π^{-s}Γ(s)Z_Q(s) = I₊ + I₋/√λ + 1/(√λ(s-1)) - 1/s
  I₊ = ∫₁^∞ (θ_Q(t)-1) t^{s-1} dt
  I₋ = ∫₁^∞ (θ_{Q⁻¹}(t)-1) t^{-s} dt

진단: 곡률 κ, 모노드로미 (폐곡선 적분), σ-국소화 FWHM.
출력: results/epstein_adversarial.txt
참고: Davenport & Heilbronn (1936), Bombieri & Hejhal
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

mpmath.mp.dps = 80

print("=" * 70)
print("[결과 #29] Epstein 제타 함수 적대적 검증")
print("=" * 70)
print(f"시작 시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"mpmath dps = {mpmath.mp.dps}")
print()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. Epstein 제타 함수 구현 (θ 함수 Mellin 적분)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def theta3(q, N_terms=40):
    """야코비 세타 함수 θ₃(0,q) = 1 + 2Σ_{n=1}^N q^{n²}"""
    result = mpmath.mpf(1)
    for n in range(1, N_terms + 1):
        term = q ** (n * n)
        if abs(term) < mpmath.mpf(10) ** (-mpmath.mp.dps + 5):
            break
        result += 2 * term
    return result


def theta_Q(t, lam, N_terms=40):
    """θ_Q(t) = Σ exp(-πt(m² + λn²)) = θ₃(e^{-πt}) · θ₃(e^{-πλt})"""
    q1 = mpmath.exp(-mpmath.pi * t)
    q2 = mpmath.exp(-mpmath.pi * lam * t)
    return theta3(q1, N_terms) * theta3(q2, N_terms)


def completed_epstein(s, lam):
    """
    완비 Epstein: Λ(s) = π^{-s}Γ(s)Z_Q(s)

    θ 함수 Mellin 적분, t=1 분할:
      Λ(s) = I₊ + I₋/√λ + 1/(√λ(s-1)) - 1/s
      I₊ = ∫₁^∞ (θ_Q(t)-1) t^{s-1} dt       [Q(m,n) = m² + λn²]
      I₋ = ∫₁^∞ (θ_{Q⁻¹}(t)-1) t^{-s} dt    [Q⁻¹(m,n) = m² + n²/λ]
    """
    s = mpmath.mpc(s)
    lam_mp = mpmath.mpf(str(lam))
    sqrt_lam = mpmath.sqrt(lam_mp)
    lam_inv = 1 / lam_mp

    def integrand_plus(t):
        return (theta_Q(t, lam_mp) - 1) * t ** (s - 1)

    def integrand_minus(t):
        return (theta_Q(t, lam_inv) - 1) * t ** (-s)

    I_plus = mpmath.quad(integrand_plus, [1, mpmath.inf])
    I_minus = mpmath.quad(integrand_minus, [1, mpmath.inf])

    Lambda_s = I_plus + I_minus / sqrt_lam + 1 / (sqrt_lam * (s - 1)) - 1 / s
    return Lambda_s


def epstein_zeta(s, lam):
    """Z_Q(s) = π^s / Γ(s) · Λ(s)"""
    s = mpmath.mpc(s)
    Lambda_s = completed_epstein(s, lam)
    return mpmath.power(mpmath.pi, s) / mpmath.gamma(s) * Lambda_s


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. 검증: λ=1에서 Z_Q vs 4ζ(s)L(s,χ_{-4})
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

print("━━━ 검증: λ=1에서 Z_{m²+n²}(s) vs 4·ζ(s)·L(s,χ_{-4}) ━━━")

chi_m4 = [0, 1, 0, -1]  # χ_{-4}: Kronecker symbol (-4/n)

verification_ok = True
for s_test in [2, 3, 4, 1.5, 2.5]:
    s_mp = mpmath.mpf(str(s_test))
    Z_computed = epstein_zeta(s_mp, 1)
    Z_exact = 4 * mpmath.zeta(s_mp) * mpmath.dirichlet(s_mp, chi_m4)
    rel_err = float(abs(Z_computed - Z_exact) / abs(Z_exact))
    status = "✅" if rel_err < 1e-10 else "❌"
    if rel_err >= 1e-10:
        verification_ok = False
    print(f"  s={s_test}: Z_comp={float(Z_computed.real):.12f}, "
          f"Z_exact={float(Z_exact.real):.12f}, err={rel_err:.2e} {status}")

if not verification_ok:
    print("❌ 검증 실패 — 구현 오류. 중단.")
    sys.exit(1)

# 임계선 위 검증 (s = 0.5 + it)
print("\n  임계선 위 검증 (복소수):")
for t_test in [5.0, 10.0, 14.134]:
    s_mp = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_test))
    Z_comp = epstein_zeta(s_mp, 1)
    Z_exact = 4 * mpmath.zeta(s_mp) * mpmath.dirichlet(s_mp, chi_m4)
    rel_err = float(abs(Z_comp - Z_exact) / abs(Z_exact)) if abs(Z_exact) > 1e-50 else float(abs(Z_comp))
    status = "✅" if rel_err < 1e-6 else "❌"
    print(f"  s=0.5+{t_test}i: |Z_comp|={float(abs(Z_comp)):.8f}, "
          f"|Z_exact|={float(abs(Z_exact)):.8f}, rel_err={rel_err:.2e} {status}")

print()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. 다발 진단 함수 (Epstein zeta 전용)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _near_zero_epstein(val):
    """dps 기반 영점 판정"""
    return abs(val) < mpmath.mpf(10) ** (-mpmath.mp.dps + 10)


def connection_epstein(s, lam):
    """접속 L = Λ'/Λ (수치 미분, h=10^{-20})"""
    h = mpmath.mpf(1) / mpmath.mpf(10**20)
    Lambda_val = completed_epstein(s, lam)
    if _near_zero_epstein(Lambda_val):
        return mpmath.mpc(1e10, 0)
    Lambda_d = (completed_epstein(s + h, lam) - completed_epstein(s - h, lam)) / (2 * h)
    return Lambda_d / Lambda_val


def curvature_epstein(s, lam):
    """곡률 κ = |Λ'/Λ|²"""
    L = connection_epstein(s, lam)
    return float(abs(L) ** 2)


def monodromy_epstein(sigma, t, lam, radius=0.3, n_steps=64):
    """
    점 s=σ+it 주위 반지름 radius 원에서 모노드로미 (폐곡선 적분).
    영점이 원 안에 있으면 ≈±2π, 없으면 ≈0.
    eps 차분 방식 금지 — 반드시 arg 누적.
    """
    s_center = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
    total_delta = 0.0
    prev_arg = None

    for k in range(n_steps + 1):
        theta = 2 * np.pi * k / n_steps
        s = s_center + radius * mpmath.exp(1j * theta)
        val = completed_epstein(s, lam)

        if _near_zero_epstein(val):
            continue

        curr_arg = float(mpmath.arg(val))

        if prev_arg is not None:
            delta = curr_arg - prev_arg
            while delta > np.pi:
                delta -= 2 * np.pi
            while delta < -np.pi:
                delta += 2 * np.pi
            total_delta += delta

        prev_arg = curr_arg

    return total_delta


def sigma_localization(t, lam, sigma_range=(0.2, 0.8), n_sigma=40):
    """σ-국소화: 고정 t에서 κ(σ) 프로파일의 FWHM 측정"""
    sigmas = np.linspace(sigma_range[0], sigma_range[1], n_sigma)
    kappas = []

    for sigma in sigmas:
        s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
        try:
            k = curvature_epstein(s, lam)
            kappas.append(min(k, 1e8) if np.isfinite(k) else 1e8)
        except Exception as e:
            print(f"    WARNING: σ={sigma:.3f} 곡률 계산 실패: {e}")
            kappas.append(0.0)

    kappas = np.array(kappas)
    peak_idx = np.argmax(kappas)
    peak_val = kappas[peak_idx]
    peak_sigma = sigmas[peak_idx]

    if peak_val <= 0:
        return peak_sigma, float('nan'), sigmas, kappas

    half_max = peak_val / 2

    # 좌측 half-max (선형 보간)
    left = sigma_range[0]
    for i in range(peak_idx, -1, -1):
        if kappas[i] < half_max:
            if i < peak_idx and kappas[i + 1] != kappas[i]:
                frac = (half_max - kappas[i]) / (kappas[i + 1] - kappas[i])
                left = sigmas[i] + frac * (sigmas[i + 1] - sigmas[i])
            break

    # 우측 half-max
    right = sigma_range[1]
    for i in range(peak_idx, len(sigmas)):
        if kappas[i] < half_max:
            if i > peak_idx and kappas[i - 1] != kappas[i]:
                frac = (half_max - kappas[i]) / (kappas[i - 1] - kappas[i])
                right = sigmas[i] - frac * (sigmas[i] - sigmas[i - 1])
            break

    fwhm = right - left
    return peak_sigma, fwhm, sigmas, kappas


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. 영점 탐색 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def find_zeros_on_line(lam, t_min=2.0, t_max=35.0, n_scan=400):
    """임계선 σ=1/2 위에서 |Z_Q(1/2+it)| 최소화로 영점 탐색"""
    print(f"  임계선 영점 탐색: λ={lam}, t∈[{t_min},{t_max}], n_scan={n_scan}")
    t0 = time.time()

    ts = np.linspace(t_min, t_max, n_scan)
    abs_vals = np.zeros(n_scan)

    for i, t in enumerate(ts):
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        try:
            val = epstein_zeta(s, lam)
            abs_vals[i] = float(abs(val))
        except Exception as e:
            abs_vals[i] = 1e10
            print(f"    WARNING: t={t:.4f} 계산 실패: {e}")
        if (i + 1) % 100 == 0:
            elapsed = time.time() - t0
            eta = elapsed / (i + 1) * (n_scan - i - 1)
            print(f"    스캔 진행: {i+1}/{n_scan} ({elapsed:.0f}s, ETA {eta:.0f}s)")

    # 극소점 찾기 (|Z| < threshold에서 findroot 정밀화)
    zeros = []
    findroot_fail = 0
    for i in range(1, len(abs_vals) - 1):
        if abs_vals[i] < abs_vals[i - 1] and abs_vals[i] < abs_vals[i + 1]:
            if abs_vals[i] < 2.0:
                t_approx = ts[i]
                # Re(Z)=0 기반 findroot
                try:
                    def f_re(t_var):
                        sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                        return mpmath.re(epstein_zeta(sv, lam))

                    t_refined = float(mpmath.findroot(f_re, mpmath.mpf(str(t_approx))))

                    # |Z| 확인
                    s_check = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_refined))
                    Z_check = epstein_zeta(s_check, lam)
                    residual = float(abs(Z_check))

                    if residual < 0.5:
                        if not zeros or abs(t_refined - zeros[-1][1]) > 0.2:
                            zeros.append((0.5, t_refined, residual))
                            print(f"    영점: t={t_refined:.6f}, |Z|={residual:.2e}")
                    else:
                        # Im(Z)=0도 시도
                        def f_im(t_var):
                            sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                            return mpmath.im(epstein_zeta(sv, lam))
                        t_refined2 = float(mpmath.findroot(f_im, mpmath.mpf(str(t_approx))))
                        s_check2 = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_refined2))
                        Z_check2 = epstein_zeta(s_check2, lam)
                        residual2 = float(abs(Z_check2))
                        if residual2 < 0.5 and (not zeros or abs(t_refined2 - zeros[-1][1]) > 0.2):
                            zeros.append((0.5, t_refined2, residual2))
                            print(f"    영점: t={t_refined2:.6f}, |Z|={residual2:.2e}")

                except Exception as e:
                    findroot_fail += 1
                    # 근사값이라도 기록
                    if abs_vals[i] < 0.3:
                        zeros.append((0.5, t_approx, abs_vals[i]))
                        print(f"    영점 (근사): t={t_approx:.4f}, |Z|={abs_vals[i]:.2e}")

    if findroot_fail > n_scan // 4:
        print(f"  ⚠️ findroot 실패 {findroot_fail}회 — 과반 초과, 탐색 로직 점검 필요")
    if len(zeros) == 0:
        print(f"  ⚠️ 영점 0개 — 탐색 범위/해상도 점검 필요")

    elapsed = time.time() - t0
    print(f"  → 임계선 영점 {len(zeros)}개 발견 ({elapsed:.0f}s)")
    return zeros


def find_zeros_off_line(lam, sigma_range=(0.25, 0.75), t_range=(2.0, 35.0),
                        n_sigma=12, n_t=120):
    """
    σ≠1/2에서 영점 탐색: |Z_Q(σ+it)| 2D 격자 스캔 후 Newton 정밀화.
    함수방정식으로 σ와 1-σ 대칭이 아닌(Q≠Q* when λ≠1) 경우 비대칭 분포 가능.
    """
    print(f"  격자 스캔: λ={lam}, σ∈{sigma_range}, t∈{t_range}, grid={n_sigma}×{n_t}")
    t0 = time.time()

    sigmas = np.linspace(sigma_range[0], sigma_range[1], n_sigma)
    ts = np.linspace(t_range[0], t_range[1], n_t)

    grid = np.full((n_sigma, n_t), 1e10)

    for i, sigma in enumerate(sigmas):
        for j, t in enumerate(ts):
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
            try:
                val = epstein_zeta(s, lam)
                grid[i, j] = float(abs(val))
            except Exception as e:
                grid[i, j] = 1e10
        elapsed = time.time() - t0
        eta = elapsed / (i + 1) * (n_sigma - i - 1)
        print(f"    σ={sigma:.3f} 완료 ({i+1}/{n_sigma}, {elapsed:.0f}s, ETA {eta:.0f}s)")

    # 2D 극소점 찾기 (σ=0.5 ± 0.03 제외)
    zeros = []
    findroot_fail = 0
    for i in range(1, n_sigma - 1):
        if abs(sigmas[i] - 0.5) < 0.04:
            continue
        for j in range(1, n_t - 1):
            if (grid[i, j] < grid[i - 1, j] and grid[i, j] < grid[i + 1, j] and
                    grid[i, j] < grid[i, j - 1] and grid[i, j] < grid[i, j + 1]):
                if grid[i, j] < 3.0:
                    sigma_approx, t_approx = sigmas[i], ts[j]

                    # 2D Newton: findroot(Re(Z), Im(Z)) = (0, 0)
                    try:
                        def f_sys(sigma_v, t_v):
                            sv = mpmath.mpf(sigma_v) + 1j * mpmath.mpf(t_v)
                            Z = epstein_zeta(sv, lam)
                            return (mpmath.re(Z), mpmath.im(Z))

                        result = mpmath.findroot(
                            f_sys,
                            (mpmath.mpf(str(sigma_approx)), mpmath.mpf(str(t_approx)))
                        )
                        sigma_ref = float(result[0])
                        t_ref = float(result[1])

                        s_check = mpmath.mpf(str(sigma_ref)) + 1j * mpmath.mpf(str(t_ref))
                        Z_check = epstein_zeta(s_check, lam)
                        residual = float(abs(Z_check))

                        if residual < 1e-5 and abs(sigma_ref - 0.5) > 0.01:
                            if not any(abs(sigma_ref - z[0]) < 0.01 and abs(t_ref - z[1]) < 0.2 for z in zeros):
                                zeros.append((sigma_ref, t_ref, residual))
                                print(f"    ⚠️ 임계선 밖 영점: σ={sigma_ref:.6f}, t={t_ref:.6f}, |Z|={residual:.2e}")
                    except Exception as e:
                        findroot_fail += 1
                        if grid[i, j] < 0.5:
                            zeros.append((sigma_approx, t_approx, grid[i, j]))
                            print(f"    후보 (미정밀): σ={sigma_approx:.3f}, t={t_approx:.3f}, |Z|={grid[i,j]:.2e}")

    if findroot_fail > 0:
        print(f"  findroot 실패: {findroot_fail}회")
    if len(zeros) == 0:
        print(f"  임계선 밖 영점 미발견 (탐색 범위: σ∈{sigma_range}, t∈{t_range})")

    elapsed = time.time() - t0
    print(f"  → 임계선 밖 영점 {len(zeros)}개 발견 ({elapsed:.0f}s)")
    return zeros


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. 메인 실험
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

results_all = {}
experiment_start = time.time()

for lam_label, lam_val in [("lambda_1_control", 1), ("lambda_5_experiment", 5)]:
    human_label = "λ=1 (대조군)" if lam_val == 1 else "λ=5 (실험군)"
    print(f"\n{'='*70}")
    print(f"[{human_label}] Q(m,n) = m² + {lam_val}n²")
    print(f"{'='*70}")

    results = {'lambda': lam_val, 'label': human_label, 'on_line': [], 'off_line': []}

    # 5a. 임계선 영점 탐색
    print(f"\n--- 임계선 (σ=1/2) 영점 탐색 ---")
    on_line_zeros = find_zeros_on_line(lam_val, t_min=2.0, t_max=35.0, n_scan=350)

    # 5b. 임계선 밖 영점 탐색
    print(f"\n--- 임계선 밖 영점 탐색 ---")
    off_line_zeros = find_zeros_off_line(lam_val,
                                          sigma_range=(0.25, 0.75),
                                          t_range=(2.0, 35.0),
                                          n_sigma=12, n_t=100)

    # 5c. 각 영점에서 다발 진단
    all_zeros = [(s, t, r, 'on') for s, t, r in on_line_zeros] + \
                [(s, t, r, 'off') for s, t, r in off_line_zeros]

    if len(all_zeros) == 0:
        print(f"  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
        results_all[lam_label] = results
        continue

    # 최대 8개 영점만 진단 (시간 절약)
    all_zeros_sorted = sorted(all_zeros, key=lambda x: x[1])[:8]

    print(f"\n--- 다발 진단 ({len(all_zeros_sorted)}개 영점) ---")
    for sigma_z, t_z, residual, location in all_zeros_sorted:
        print(f"\n  영점: σ={sigma_z:.6f}, t={t_z:.6f} ({location}-line), |Z|={residual:.2e}")
        diag = {
            'sigma': sigma_z, 't': t_z, 'location': location,
            'residual': residual,
            'kappa': float('nan'), 'mono': float('nan'),
            'peak_sigma': float('nan'), 'fwhm': float('nan'),
        }

        # 곡률 (영점 위 직접 측정 금지 — δ=0.03 오프셋)
        s_offset = mpmath.mpf(str(sigma_z)) + 1j * mpmath.mpf(str(t_z + 0.03))
        try:
            kappa = curvature_epstein(s_offset, lam_val)
            diag['kappa'] = kappa
            print(f"    곡률 κ (δ=0.03): {kappa:.4f}")
        except Exception as e:
            print(f"    곡률 계산 실패: {e}")

        # 모노드로미 (폐곡선 적분 — eps 차분 금지)
        try:
            mono = monodromy_epstein(sigma_z, t_z, lam_val, radius=0.3, n_steps=64)
            diag['mono'] = mono
            print(f"    모노드로미: {mono:.4f} (≈{mono/np.pi:.2f}π)")
        except Exception as e:
            print(f"    모노드로미 실패: {e}")

        # σ-국소화 (FWHM)
        try:
            peak_sigma, fwhm, _, _ = sigma_localization(t_z, lam_val,
                                                         sigma_range=(0.2, 0.8), n_sigma=30)
            diag['peak_sigma'] = peak_sigma
            diag['fwhm'] = fwhm
            print(f"    σ-국소화: peak_σ={peak_sigma:.4f}, FWHM={fwhm:.4f}")
        except Exception as e:
            print(f"    σ-국소화 실패: {e}")

        if location == 'on':
            results['on_line'].append(diag)
        else:
            results['off_line'].append(diag)

    results_all[lam_label] = results


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. 결과 정리 및 저장
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

total_elapsed = time.time() - experiment_start

print("\n" + "=" * 70)
print("최종 결과 요약")
print("=" * 70)

out = []
out.append("=" * 70)
out.append("[결과 #29] Epstein 제타 함수 적대적 검증")
out.append(f"실행 시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
out.append(f"총 실행 시간: {total_elapsed:.0f}s ({total_elapsed/60:.1f}분)")
out.append(f"mpmath dps = {mpmath.mp.dps}")
out.append("=" * 70)
out.append("")
out.append("설정:")
out.append("  Q(m,n) = m² + λn²")
out.append("  λ=1 (대조군): Z = 4ζ(s)L(s,χ_{-4}), RH 만족")
out.append("  λ=5 (실험군): h(-20)=2, 개별 Epstein zeta에 RH 위반 가능")
out.append("  탐색 범위: σ∈[0.25,0.75], t∈[2,35]")
out.append("  다발 진단: κ (곡률), mono (모노드로미), FWHM (σ-국소화)")
out.append("")

for lam_label, results in results_all.items():
    human = results.get('label', lam_label)
    out.append(f"--- {human} ---")
    out.append(f"임계선 영점 진단: {len(results['on_line'])}개")
    out.append(f"임계선 밖 영점 진단: {len(results['off_line'])}개")

    for tag, zeros_list in [("임계선", results['on_line']), ("임계선 밖", results['off_line'])]:
        if zeros_list:
            out.append(f"\n  [{tag} 영점]")
            out.append(f"  {'σ':>8s} {'t':>10s} {'κ':>12s} {'mono/π':>8s} {'peak_σ':>8s} {'FWHM':>8s}")
            for d in zeros_list:
                mono_pi = d['mono'] / np.pi if np.isfinite(d['mono']) else float('nan')
                out.append(
                    f"  {d['sigma']:8.4f} {d['t']:10.4f} {d['kappa']:12.2f} "
                    f"{mono_pi:8.2f} {d['peak_sigma']:8.4f} {d['fwhm']:8.4f}"
                )
            kappas = [d['kappa'] for d in zeros_list if np.isfinite(d['kappa'])]
            monos = [abs(d['mono'] / np.pi) for d in zeros_list if np.isfinite(d['mono'])]
            fwhms = [d['fwhm'] for d in zeros_list if np.isfinite(d['fwhm'])]
            if kappas:
                out.append(f"  κ 평균: {np.mean(kappas):.2f} ± {np.std(kappas):.2f}")
            if monos:
                out.append(f"  |mono/π| 평균: {np.mean(monos):.2f} ± {np.std(monos):.2f}")
            if fwhms:
                out.append(f"  FWHM 평균: {np.mean(fwhms):.4f} ± {np.std(fwhms):.4f}")
    out.append("")

# 판정
out.append("=" * 70)
out.append("비교 분석 및 판정")
out.append("=" * 70)

lam5 = results_all.get("lambda_5_experiment", {'on_line': [], 'off_line': []})
lam1 = results_all.get("lambda_1_control", {'on_line': [], 'off_line': []})

has_off = len(lam5.get('off_line', [])) > 0
has_on_5 = len(lam5.get('on_line', [])) > 0
has_on_1 = len(lam1.get('on_line', [])) > 0

if has_off and has_on_5:
    out.append(f"\nσ≠1/2 영점 {len(lam5['off_line'])}개 발견 → 적대적 검증 가능")

    on_k = [d['kappa'] for d in lam5['on_line'] if np.isfinite(d['kappa'])]
    off_k = [d['kappa'] for d in lam5['off_line'] if np.isfinite(d['kappa'])]
    on_m = [abs(d['mono'] / np.pi) for d in lam5['on_line'] if np.isfinite(d['mono'])]
    off_m = [abs(d['mono'] / np.pi) for d in lam5['off_line'] if np.isfinite(d['mono'])]
    on_f = [d['fwhm'] for d in lam5['on_line'] if np.isfinite(d['fwhm'])]
    off_f = [d['fwhm'] for d in lam5['off_line'] if np.isfinite(d['fwhm'])]

    out.append("\n비교표:")
    out.append(f"  {'지표':>12s} | {'σ=1/2':>15s} | {'σ≠1/2':>15s} | {'비율':>8s}")
    out.append(f"  {'-'*12}-+-{'-'*15}-+-{'-'*15}-+-{'-'*8}")
    if on_k and off_k:
        r = np.mean(off_k) / np.mean(on_k) if np.mean(on_k) > 0 else float('inf')
        out.append(f"  {'κ 평균':>12s} | {np.mean(on_k):15.2f} | {np.mean(off_k):15.2f} | {r:8.2f}")
    if on_m and off_m:
        out.append(f"  {'|mono/π|':>12s} | {np.mean(on_m):15.2f} | {np.mean(off_m):15.2f} | {'—':>8s}")
    if on_f and off_f:
        r = np.mean(off_f) / np.mean(on_f) if np.mean(on_f) > 0 else float('inf')
        out.append(f"  {'FWHM':>12s} | {np.mean(on_f):15.4f} | {np.mean(off_f):15.4f} | {r:8.2f}")

    # 판별 기준: κ 차이, mono 편차, FWHM 확대
    discriminated = False
    reasons = []
    if on_k and off_k and abs(np.mean(on_k) - np.mean(off_k)) > 0.5 * np.mean(on_k):
        discriminated = True
        reasons.append("κ 차이 유의")
    if on_m and off_m:
        on_near_pi = sum(1 for m in on_m if abs(m - 1.0) < 0.3) / len(on_m)
        off_near_pi = sum(1 for m in off_m if abs(m - 1.0) < 0.3) / len(off_m)
        if on_near_pi > off_near_pi + 0.3:
            discriminated = True
            reasons.append("mono ±π 준수율 차이")
    if on_f and off_f and np.mean(off_f) > 1.3 * np.mean(on_f):
        discriminated = True
        reasons.append("FWHM 확대")

    if discriminated:
        out.append(f"\n→ 양성: 다발 서명이 σ=1/2 영점과 σ≠1/2 영점에서 정량적으로 구분됨")
        out.append(f"  근거: {', '.join(reasons)}")
        out.append(f"  해석: 프레임워크의 falsifiability 입증 — RH 위반 영점은 다른 기하학적 서명")
    else:
        out.append(f"\n→ 음성: 다발 서명이 유사 — 프레임워크가 σ=1/2 vs σ≠1/2를 구분하지 못함")
        out.append(f"  해석: 프레임워크 한계로 정직하게 기재")

elif has_on_5 and not has_off:
    out.append(f"\nσ≠1/2 영점 미발견 (탐색: σ∈[0.25,0.75], t∈[2,35], 12×100 격자)")
    out.append("→ 부분 음성: 탐색 범위 내 임계선 밖 영점 없음")
    out.append("  가능 원인: (a) 영점이 더 높은 t에 존재, (b) σ 범위 부족, (c) λ=5 형식의 특성")

    # 대조군 vs 실험군 on-line 비교 가능
    if has_on_1:
        out.append("\n대조군(λ=1) vs 실험군(λ=5) 임계선 영점 비교:")
        for tag, data in [("λ=1", lam1['on_line']), ("λ=5", lam5['on_line'])]:
            k = [d['kappa'] for d in data if np.isfinite(d['kappa'])]
            m = [abs(d['mono'] / np.pi) for d in data if np.isfinite(d['mono'])]
            f = [d['fwhm'] for d in data if np.isfinite(d['fwhm'])]
            if k:
                out.append(f"  {tag}: κ={np.mean(k):.2f}±{np.std(k):.2f}, "
                           f"|mono/π|={np.mean(m):.2f}±{np.std(m):.2f}, "
                           f"FWHM={np.mean(f):.4f}±{np.std(f):.4f}")
else:
    out.append("\n영점 탐색 실패 — 판정 불가")

out.append("")
out.append(f"총 실행 시간: {total_elapsed:.0f}s ({total_elapsed/60:.1f}분)")
out.append(f"완료 시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")

output_text = "\n".join(out)
print(output_text)

# 결과 파일 저장
results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(results_dir, exist_ok=True)
output_path = os.path.join(results_dir, 'epstein_adversarial.txt')
with open(output_path, 'w') as f:
    f.write(output_text)
print(f"\n결과 저장: {output_path}")
print(f"완료 시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
