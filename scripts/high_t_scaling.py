#!/usr/bin/env python3
"""
결과 #31: 고 t-범위 핵심 결과 스케일링 검증
==============================================
t ≈ 1000, 5000, 10000에서 Riemann ξ 다발 성질(κ, 모노드로미, σ-프로파일) 측정.
비교 기준: t<100에서의 확립값 (κ~1030, mono=2π, peak_σ~0.494, FWHM~0.016)

핵심 설계 원칙:
  - mpmath.diff(xi, s) 금지 (t>300에서 실패)
  - bundle_utils.connection_zeta 금지 (h=10^{-20} underflow for t>1000)
  - ξ 직접 계산 금지 (t≥1000에서 Γ-인수 underflow ~10^{-3000})
  → 해석적 log-derivative 공식: ξ'/ξ = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2 + ζ'(s)/ζ(s)
  → 모노드로미: arg(ξ)를 log-sum으로 계산 (ξ 직접 계산 불필요)
  → ζ(s), ζ'(s): 직접 계산 가능 (Riemann-Siegel, underflow 없음)

σ-프로파일 표준 (사이클 44 교훈):
  - σ₀±0.3, 201점(t<2000) / 101점(t<6000) / 51점(t≥6000)
  - 간격 ≈0.003
"""

import sys, os, time
import numpy as np
import mpmath

BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_PATH = os.path.join(BASE_DIR, '..', 'results', 'high_t_scaling.txt')
os.makedirs(os.path.join(BASE_DIR, '..', 'results'), exist_ok=True)

_log_buf = []
def log(msg=''):
    print(msg, flush=True)
    _log_buf.append(str(msg))

# ═══════════════════════════════════════════════════════════════════════════
# 핵심 공식: 해석적 log-derivative
# ═══════════════════════════════════════════════════════════════════════════

def xi_log_deriv_analytic(s):
    """
    ξ'(s)/ξ(s) = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2 + ζ'(s)/ζ(s)

    이 공식은 ξ를 직접 계산하지 않아 고 t에서 underflow 없음.
    ζ(s)는 Riemann-Siegel로 직접 계산 → 크기 O(t^{1/6}), underflow 없음.
    ζ'(s)는 h=1e-6 유한차분 → dps≥50에서 정확.
    """
    z = mpmath.zeta(s)
    # 영점 위 발산 방어
    if abs(z) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e8, 0)
    # ζ'(s): 유한차분 (h=1e-6, ζ 크기 O(t^{1/6}) → 차분은 O(h·ζ')로 안전)
    h = mpmath.mpf('1e-6')
    zeta_d = (mpmath.zeta(s + h) - mpmath.zeta(s - h)) / (2 * h)
    return (mpmath.mpf(1) / s
            + mpmath.mpf(1) / (s - 1)
            - mpmath.log(mpmath.pi) / 2
            + mpmath.digamma(s / 2) / 2
            + zeta_d / z)


def kappa_analytic(s):
    """곡률 κ = |ξ'/ξ|² (해석적 공식)"""
    return float(abs(xi_log_deriv_analytic(s))**2)


def xi_log_arg_analytic(s):
    """
    arg(ξ(s)) — log-sum 방식:
      arg(ξ) = arg(s) + arg(s-1) + Re(-s/2·logπ이미부) + Im(logΓ(s/2)) + arg(ζ(s))
    ξ 직접 계산 없이 각 항 분리 → Γ 인수 underflow 완전 회피.
    영점 위(ζ=0)에서는 nan 반환.
    """
    t_val = float(mpmath.im(s))
    a_s   = float(mpmath.arg(s))
    a_sm1 = float(mpmath.arg(s - 1))
    # π^{-s/2}: arg = -Im(s)/2 · log(π)
    a_pi  = -t_val / 2.0 * float(mpmath.log(mpmath.pi))
    # Γ(s/2): loggamma의 허수부 = arg(Γ)
    a_gam = float(mpmath.im(mpmath.loggamma(s / 2)))
    z = mpmath.zeta(s)
    if abs(z) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return float('nan')
    a_zeta = float(mpmath.arg(z))
    return a_s + a_sm1 + a_pi + a_gam + a_zeta


# ═══════════════════════════════════════════════════════════════════════════
# 측정 함수
# ═══════════════════════════════════════════════════════════════════════════

def measure_kappa(t_zero, delta=0.03):
    """κ at σ=0.5+δ, t=t_zero (δ=0.03 표준 오프셋)"""
    s = mpmath.mpf(str(0.5 + delta)) + 1j * mpmath.mpf(str(t_zero))
    return kappa_analytic(s)


def measure_monodromy(t_zero, radius=0.2, n_steps=64):
    """
    모노드로미: 반지름 radius의 폐곡선 적분, log-arg 누적
    [주의] bundle_utils.monodromy_contour = ξ 직접 계산 → 고 t 실패
    이 함수는 xi_log_arg_analytic 사용 → underflow 없음
    """
    s_center = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero))
    total_delta = 0.0
    prev_arg = None
    skipped = 0

    for k in range(n_steps + 1):
        theta = 2.0 * np.pi * k / n_steps
        s = s_center + radius * mpmath.exp(1j * mpmath.mpf(str(theta)))
        try:
            curr_arg = xi_log_arg_analytic(s)
        except Exception as e:
            print(f"    WARNING mono step {k}: {e}", flush=True)
            skipped += 1
            continue

        if curr_arg is None or (isinstance(curr_arg, float) and np.isnan(curr_arg)):
            skipped += 1
            continue

        if prev_arg is not None:
            delta = curr_arg - prev_arg
            while delta >  np.pi: delta -= 2 * np.pi
            while delta < -np.pi: delta += 2 * np.pi
            total_delta += delta
        prev_arg = curr_arg

    return total_delta, skipped


def measure_sigma_profile(t_zero):
    """
    σ-프로파일: σ ∈ [0.2, 0.8], n_sigma 점
    t<2000: 201점(간격0.003), t<6000: 101점(간격0.006), t≥6000: 51점(간격0.012)
    영점 위(σ=0.5 근방) 발산점은 NaN 처리 후 제외
    """
    if t_zero < 2000:
        n_sigma = 201
    elif t_zero < 6000:
        n_sigma = 101
    else:
        n_sigma = 51

    sigmas = np.linspace(0.2, 0.8, n_sigma)
    kappas = np.full(n_sigma, np.nan)

    for i, sig in enumerate(sigmas):
        s = mpmath.mpf(str(float(sig))) + 1j * mpmath.mpf(str(t_zero))
        try:
            k = kappa_analytic(s)
            if np.isfinite(k) and k < 5e8:   # 발산(영점 위) 제외
                kappas[i] = k
        except Exception as e:
            print(f"    WARNING profile σ={sig:.4f}: {e}", flush=True)

        if (i + 1) % 25 == 0:
            print(f"    σ-프로파일 {i+1}/{n_sigma} 완료", flush=True)

    # peak 및 FWHM
    valid = np.isfinite(kappas)
    if not np.any(valid):
        return sigmas, kappas, 0.5, float('nan'), n_sigma

    peak_idx   = int(np.nanargmax(kappas))
    peak_sigma = float(sigmas[peak_idx])
    peak_kappa = float(kappas[peak_idx])

    half_max = peak_kappa / 2.0
    above    = (kappas >= half_max) & valid
    above_i  = np.where(above)[0]

    if len(above_i) >= 2:
        fwhm = float(sigmas[above_i[-1]] - sigmas[above_i[0]])
    else:
        fwhm = float(sigmas[1] - sigmas[0])   # 이산화 한계

    return sigmas, kappas, peak_sigma, fwhm, n_sigma


# ═══════════════════════════════════════════════════════════════════════════
# 단일 영점 진단
# ═══════════════════════════════════════════════════════════════════════════

def diagnose_zero(N, t_zero, dps):
    """영점 #N (t=t_zero)에서 κ, 모노드로미, σ-프로파일 측정"""
    log(f"\n{'='*60}")
    log(f"영점 #{N}:  t = {t_zero:.6f}   dps = {dps}")
    log(f"{'='*60}")
    mpmath.mp.dps = dps

    res = {'N': N, 't': t_zero, 'dps': dps}
    t_start = time.time()

    # ── 1. κ (δ=0.03 오프셋) ──────────────────────────────────────
    log("  [1] κ (σ=0.53 = 0.5+δ, δ=0.03)...")
    try:
        k = measure_kappa(t_zero, delta=0.03)
        res['kappa'] = k
        log(f"      κ = {k:.2f}   (기준 ~1030)")
    except Exception as e:
        log(f"      ERROR κ: {e}")
        res['kappa'] = float('nan')

    # ── 2. 모노드로미 ─────────────────────────────────────────────
    log("  [2] 모노드로미 (r=0.2, n=64)...")
    try:
        t_m = time.time()
        mono, skipped = measure_monodromy(t_zero, radius=0.2, n_steps=64)
        mono_pi = mono / np.pi
        res['mono']    = mono
        res['mono_pi'] = mono_pi
        log(f"      mono = {mono:.4f} ({mono_pi:.4f}π)  skipped={skipped}  {time.time()-t_m:.1f}s")
    except Exception as e:
        log(f"      ERROR 모노드로미: {e}")
        res['mono']    = float('nan')
        res['mono_pi'] = float('nan')

    # ── 3. σ-프로파일 ─────────────────────────────────────────────
    log("  [3] σ-프로파일 (σ₀±0.3)...")
    try:
        t_p = time.time()
        sigmas, kappas, peak_s, fwhm, n_sig = measure_sigma_profile(t_zero)
        spacing = 0.6 / (n_sig - 1)
        res['peak_sigma'] = peak_s
        res['fwhm']       = fwhm
        log(f"      peak_σ = {peak_s:.4f}   (기준 ~0.494)")
        log(f"      FWHM   = {fwhm:.4f}   (격자간격 {spacing:.4f})  {time.time()-t_p:.1f}s")
    except Exception as e:
        log(f"      ERROR σ-프로파일: {e}")
        res['peak_sigma'] = float('nan')
        res['fwhm']       = float('nan')

    res['elapsed'] = time.time() - t_start
    log(f"  영점 #{N} 총 소요: {res['elapsed']:.1f}s")
    return res


# ═══════════════════════════════════════════════════════════════════════════
# 메인
# ═══════════════════════════════════════════════════════════════════════════

def main():
    log("결과 #31: 고 t-범위 핵심 결과 스케일링 검증")
    log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log()
    log("방법:")
    log("  ξ'/ξ = 1/s + 1/(s-1) - log(π)/2 + ψ(s/2)/2 + ζ'(s)/ζ(s)  [해석적]")
    log("  arg(ξ): log-sum [Γ underflow 회피]")
    log("  ζ, ζ': Riemann-Siegel (직접, underflow 없음)")
    log()

    # 기준값 (t<100, 확립)
    BASELINE = {
        'kappa':      1029.65,
        'mono_pi':    2.0000,
        'peak_sigma': 0.4937,
        'fwhm':       0.0163,
    }

    # ─── 영점 좌표 획득 ────────────────────────────────────────────
    # von Mangoldt: N(T) ≈ T/(2π)*(log(T/(2π))-1) + 7/8
    # N≈649→t≈1000, N≈4520→t≈5000, N≈10142→t≈10000
    TARGET_Ns = [
        (649,   60),   # t ≈ 1000
        (4520,  80),   # t ≈ 5000
        (10142, 80),   # t ≈ 10000
    ]

    log("1단계: zetazero 좌표 획득 (dps=50)")
    log("-" * 50)
    zero_info = []
    for N, dps in TARGET_Ns:
        log(f"  zetazero({N}) 계산 중...")
        mpmath.mp.dps = 50
        t0 = time.time()
        try:
            z = mpmath.zetazero(N)
            t_val = float(z.imag)
            elapsed = time.time() - t0
            log(f"  → zetazero({N}) = {t_val:.6f}  ({elapsed:.1f}s)")
            zero_info.append((N, t_val, dps))
        except Exception as e:
            log(f"  ERROR zetazero({N}): {e}")

    if not zero_info:
        log("ERROR: 영점을 하나도 획득하지 못했습니다.")
        return

    log()
    log("2단계: 측정")
    log("-" * 50)

    all_results = []
    for N, t_val, dps in zero_info:
        res = diagnose_zero(N, t_val, dps)
        all_results.append(res)

    # ─── 결과 표 ──────────────────────────────────────────────────
    log()
    log("=" * 75)
    log("최종 결과 표")
    log("=" * 75)
    log(f"{'N':>6}  {'t':>12}  {'κ':>10}  {'κ비':>6}  {'mono/π':>8}  {'peak_σ':>8}  {'FWHM':>7}  {'판정':>8}")
    log("-" * 75)

    judgments = []
    for res in all_results:
        kappa   = res.get('kappa',      float('nan'))
        mpi     = res.get('mono_pi',    float('nan'))
        peak_s  = res.get('peak_sigma', float('nan'))
        fwhm    = res.get('fwhm',       float('nan'))
        k_ratio = kappa / BASELINE['kappa'] if np.isfinite(kappa) else float('nan')

        ok_mono  = np.isfinite(mpi)    and (1.8 < mpi    < 2.2)
        ok_peak  = np.isfinite(peak_s) and (0.47 < peak_s < 0.53)
        ok_kappa = np.isfinite(kappa)  and (400  < kappa  < 4000)
        passed   = sum([ok_mono, ok_peak, ok_kappa])

        if   passed == 3: judge = '✅강양성'
        elif passed == 2: judge = '✅양성'
        elif passed == 1: judge = '⚠️약성'
        else:             judge = '❌음성'
        judgments.append((passed, judge))

        log(f"{res['N']:>6}  {res['t']:>12.4f}  {kappa:>10.1f}  "
            f"{k_ratio:>5.2f}×  {mpi:>8.4f}  {peak_s:>8.4f}  "
            f"{fwhm:>7.4f}  {judge}")

    log()
    log("기준값 (t<100):")
    log(f"  κ={BASELINE['kappa']:.2f}, mono={BASELINE['mono_pi']:.4f}π, "
        f"peak_σ={BASELINE['peak_sigma']:.4f}, FWHM={BASELINE['fwhm']:.4f}")
    log()

    n_full = sum(1 for p, _ in judgments if p == 3)
    n_good = sum(1 for p, _ in judgments if p >= 2)
    n_zero = len(judgments)

    if n_full == n_zero:
        final = f"★ 강한 양성: {n_full}/{n_zero} 완전 통과 → 고 t 스케일링 확립"
    elif n_good >= 2:
        final = f"✅ 양성: {n_good}/{n_zero} 통과 → t-의존성 추가 조사 가능"
    elif n_good == 1:
        final = f"⚠️ 중립: 1/{n_zero} 통과 → 고 t 한계 관찰"
    else:
        final = f"❌ 음성: 0/{n_zero} → 프레임워크 고 t 한계"

    log(f"종합 판정: {final}")
    log()
    log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # ─── 파일 저장 ─────────────────────────────────────────────────
    with open(RESULTS_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(_log_buf) + '\n')

    print(f"\n결과 저장: {RESULTS_PATH}", flush=True)


if __name__ == '__main__':
    main()
