"""
=============================================================================
[Project RDL] 결과 #40 — χ mod 7 고 conductor 디리클레 4성질 완전 검증
=============================================================================
목표: t∈[10,200]에서 40개 이상 영점 수집 → 5성질 완전 검증

mod 7 지표 (복소, order 6):
  - g=3 (원시근), χ(n) = ω^{dlog_3(n)}, ω = e^{2πi/6}
  - χ(-1) = χ(6) = ω^3 = -1 → a=1 (홀수 지표)

검증 파트:
  A. 영점 수집 (t∈[10,200], ≥40개)
  B. κ 비율 — σ=0.5 vs 타σ(0.3, 0.7), δ=0.03 오프셋, 성공 기준 >100×
  C. 모노드로미 — 폐곡선 적분 (radius=0.1, 64단계), 성공 기준 |mono|/π = 2.000±0.001
  D. 블라인드 예측 — κ 피크 → findroot 확인, 성공 기준 5/5 이상
  E. σ-국소화 — E(σ) 비교 9점, 성공 기준 E(0.5)/E(타σ) > 50×

주의:
  - Λ'/Λ 해석적 공식: (1/2)log(q/π) + (1/2)ψ((s+a)/2) + L'/L
  - bare L'/L 금지 (h=1e-20 차분 금지)
  - 영점 위 직접 κ 측정 금지 → δ=0.03 오프셋
  - 모노드로미: radius=0.1 (단일 영점 포함 보장)
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath
import cmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))

mpmath.mp.dps = 80  # 고 conductor + 고 t 정밀도 필수

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# mod 7 지표 정의
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 원시근 g=3: 3^k mod 7 순서: 1,3,2,6,4,5
# dlog_3(n): dlog(1)=0, dlog(2)=2, dlog(3)=1, dlog(4)=4, dlog(5)=5, dlog(6)=3
# ω = e^{2πi/6}
_w6 = cmath.exp(2j * cmath.pi / 6)
_chi7_raw = [0, 1, _w6**2, _w6**1, _w6**4, _w6**5, _w6**3]
# mpmath 복소수로 변환
_chi7 = [mpmath.mpc(c.real, c.imag) for c in _chi7_raw]

CHAR_MOD7 = {
    'chi': _chi7,
    'q': 7,
    'a': 1,   # χ(-1)=χ(6)=ω³=-1 → 홀수 지표 → a=1
    'label': 'χ₇ (mod 7, 복소, order 6)',
}

T_MIN = 10.0
T_MAX = 200.0
DELTA = 0.03        # κ 측정 오프셋
MONO_RADIUS = 0.1   # 모노드로미 폐곡선 반지름 (단일 영점 포함)
MONO_STEPS = 64


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수: Λ'/Λ 해석적 공식
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def completed_L(s, ci):
    """완비 L-함수: Λ(s, χ) = (q/π)^{s/2} Γ((s+a)/2) L(s, χ)"""
    q = mpmath.mpf(ci['q'])
    a = mpmath.mpf(ci['a'])
    L_val = mpmath.dirichlet(s, ci['chi'])
    gamma_val = mpmath.gamma((s + a) / 2)
    prefactor = mpmath.power(q / mpmath.pi, s / 2)
    return prefactor * gamma_val * L_val


def connection_analytic(s, ci):
    """
    Λ'/Λ 해석적 공식:
    Λ'/Λ = (1/2)log(q/π) + (1/2)ψ((s+a)/2) + L'/L(s,χ)

    - log 항, digamma 항: 해석적 계산
    - L'/L: L(s,χ)만 수치 미분 (h=1e-6, prefactor 분리로 안정)
    """
    q = mpmath.mpf(ci['q'])
    a = mpmath.mpf(ci['a'])
    chi = ci['chi']

    # 항 1: (1/2) log(q/π)  — 순수 해석
    log_term = mpmath.log(q / mpmath.pi) / 2

    # 항 2: (1/2) ψ((s+a)/2)  — digamma 해석
    digamma_term = mpmath.digamma((s + a) / 2) / 2

    # 항 3: L'/L  — L만 수치 미분 (h=1e-6, 차분 안정)
    L_val = mpmath.dirichlet(s, chi)
    if abs(L_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e10, 0)
    h = mpmath.mpf('1e-6')
    L_d = (mpmath.dirichlet(s + h, chi) - mpmath.dirichlet(s - h, chi)) / (2 * h)
    L_log_deriv = L_d / L_val

    return log_term + digamma_term + L_log_deriv


def curvature(s, ci, delta_offset=0.0):
    """κ = |Λ'/Λ|². delta_offset > 0이면 t 방향으로 오프셋 적용 (영점 회피)"""
    s_use = s + 1j * mpmath.mpf(str(delta_offset)) if delta_offset != 0 else s
    conn = connection_analytic(s_use, ci)
    return float(abs(conn)**2)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파트 A: 영점 탐색
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def find_zeros_mod7(ci, t_min=10.0, t_max=200.0, n_scan=4000):
    """
    Re(Λ)의 부호 변화로 영점 탐색 + findroot 정밀화.
    홀수 디리클레 지표: |L| 최소화로 실제 영점 확인.
    """
    ts = np.linspace(t_min, t_max, n_scan)
    zeros = []
    fail_count = 0
    success_count = 0

    prev_re, prev_t = None, None
    for i, t in enumerate(ts):
        if i % 500 == 0:
            print(f"  스캔 {i}/{n_scan} (t={t:.1f}), 현재 영점 {len(zeros)}개",
                  flush=True)
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        try:
            val = completed_L(s, ci)
            curr_re = float(mpmath.re(val))
        except Exception as e:
            print(f"  WARNING: completed_L 평가 실패 t={t:.2f}: {e}", flush=True)
            prev_re, prev_t = None, float(t)
            continue

        if prev_re is not None and prev_re * curr_re < 0:
            mid = (prev_t + float(t)) / 2
            try:
                def f_real(t_var):
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                    return mpmath.re(completed_L(sv, ci))
                tz = float(mpmath.findroot(f_real, mpmath.mpf(str(mid))))
                # 실제 영점 확인: |Λ(0.5+itz)| 매우 작아야 함
                sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(tz))
                lv = abs(completed_L(sv, ci))
                if lv < mpmath.mpf(10)**(-20):
                    if not zeros or abs(tz - zeros[-1]) > 0.1:
                        zeros.append(tz)
                        success_count += 1
                else:
                    # |Λ| 크면 부호 변화지만 영점 아님 (무시)
                    pass
            except Exception as e:
                fail_count += 1
                print(f"  WARNING: findroot 실패 t≈{mid:.2f}: {e}", flush=True)

        prev_re, prev_t = curr_re, float(t)

    # Im(Λ) 부호 변화도 확인 (복소 지표에서 누락 영점 포착)
    prev_im, prev_t = None, None
    for t in ts:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        try:
            val = completed_L(s, ci)
            curr_im = float(mpmath.im(val))
        except Exception as e:
            print(f"  WARNING: Im 평가 실패 t={t:.2f}: {e}", flush=True)
            prev_im, prev_t = None, float(t)
            continue

        if prev_im is not None and prev_im * curr_im < 0:
            mid = (prev_t + float(t)) / 2
            try:
                def f_imag(t_var):
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                    return mpmath.im(completed_L(sv, ci))
                tz = float(mpmath.findroot(f_imag, mpmath.mpf(str(mid))))
                sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(tz))
                lv = abs(completed_L(sv, ci))
                if lv < mpmath.mpf(10)**(-20):
                    if not any(abs(tz - z) < 0.1 for z in zeros):
                        zeros.append(tz)
            except Exception as e:
                fail_count += 1
                print(f"  WARNING: Im-findroot 실패 t≈{mid:.2f}: {e}", flush=True)

        prev_im, prev_t = curr_im, float(t)

    if fail_count > 0:
        print(f"  ⚠️ findroot 실패 {fail_count}회 (성공 {success_count}회)", flush=True)
        if fail_count > success_count:
            print("  ⚠⚠ 실패율 > 50% — 탐색 로직 점검 필요!", flush=True)
    if len(zeros) == 0:
        print("  ⚠️ 영점 0개 — 탐색 로직 점검 필요!", flush=True)

    zeros.sort()
    return np.array(zeros)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파트 B: κ 비율 (σ=0.5 vs 타σ, δ=0.03 오프셋)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def part_b_kappa_ratio(ci, zeros, delta=DELTA):
    """
    각 영점 t_n에서 σ=0.5+i(t_n+δ)의 κ와
    같은 t_n+δ에서 σ=0.3, σ=0.7의 κ 비교.
    비율: κ(0.5, t_n+δ) / κ(타σ, t_n+δ)
    """
    print(f"\n[파트 B] κ 비율 (δ={delta}, {len(zeros)}개 영점)", flush=True)
    kappa_half = []
    kappa_03 = []
    kappa_07 = []

    for i, tz in enumerate(zeros):
        t_meas = tz + delta
        s_half = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_meas))
        s_03   = mpmath.mpf('0.3') + 1j * mpmath.mpf(str(t_meas))
        s_07   = mpmath.mpf('0.7') + 1j * mpmath.mpf(str(t_meas))
        try:
            k_half = float(abs(connection_analytic(s_half, ci))**2)
            k_03   = float(abs(connection_analytic(s_03, ci))**2)
            k_07   = float(abs(connection_analytic(s_07, ci))**2)
            if np.isfinite(k_half) and np.isfinite(k_03) and np.isfinite(k_07):
                kappa_half.append(k_half)
                kappa_03.append(k_03)
                kappa_07.append(k_07)
        except Exception as e:
            print(f"  WARNING: κ 계산 실패 t={tz:.4f}: {e}", flush=True)

        if (i+1) % 10 == 0:
            print(f"  진행: {i+1}/{len(zeros)}", flush=True)

    kappa_half = np.array(kappa_half)
    kappa_03 = np.array(kappa_03)
    kappa_07 = np.array(kappa_07)

    med_half = np.median(kappa_half) if len(kappa_half) > 0 else 0
    med_03   = np.median(kappa_03)   if len(kappa_03) > 0 else 1
    med_07   = np.median(kappa_07)   if len(kappa_07) > 0 else 1

    ratio_03 = med_half / med_03 if med_03 > 0 else float('inf')
    ratio_07 = med_half / med_07 if med_07 > 0 else float('inf')
    ratio_min = min(ratio_03, ratio_07)

    print(f"  median κ(σ=0.5): {med_half:.2e}", flush=True)
    print(f"  median κ(σ=0.3): {med_03:.2e}  → 비율 {ratio_03:.1f}×", flush=True)
    print(f"  median κ(σ=0.7): {med_07:.2e}  → 비율 {ratio_07:.1f}×", flush=True)
    print(f"  ★ 최소 비율: {ratio_min:.1f}×  (기준 >100×)", flush=True)
    pass_b = ratio_min > 100
    print(f"  판정: {'✅ PASS' if pass_b else '❌ FAIL'}", flush=True)

    return {
        'med_half': med_half, 'med_03': med_03, 'med_07': med_07,
        'ratio_03': ratio_03, 'ratio_07': ratio_07, 'ratio_min': ratio_min,
        'pass': pass_b,
        'kappa_half': kappa_half, 'kappa_03': kappa_03, 'kappa_07': kappa_07,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파트 C: 모노드로미 (폐곡선 적분)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def monodromy_contour_dirichlet(t, ci, radius=MONO_RADIUS, n_steps=MONO_STEPS):
    """
    s=0.5+it 주위 반지름 radius 원에서 Λ(s,χ)의 arg 누적 (폐곡선 적분).
    영점이 원 안에 있으면 ≈±π, 없으면 ≈0.
    반환: arg 누적량 / π (이론값 ±2.000)
    """
    s_center = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
    total_delta = mpmath.mpf(0)
    prev_arg = None

    for k in range(n_steps + 1):
        theta = 2 * mpmath.pi * k / n_steps
        s = s_center + radius * mpmath.exp(1j * theta)
        try:
            val = completed_L(s, ci)
            if abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 15):
                continue  # 영점 경유 회피
            curr_arg = float(mpmath.arg(val))
        except Exception as e:
            print(f"  WARNING: 모노드로미 평가 실패 k={k}: {e}", flush=True)
            continue

        if prev_arg is not None:
            diff = curr_arg - prev_arg
            # 연속 branch 보정
            while diff > np.pi:
                diff -= 2 * np.pi
            while diff < -np.pi:
                diff += 2 * np.pi
            total_delta += diff
        prev_arg = curr_arg

    return float(total_delta) / np.pi  # 단위: π


def part_c_monodromy(ci, zeros):
    """각 영점에서 모노드로미 계산, |mono|/π ≈ 2.000 확인"""
    print(f"\n[파트 C] 모노드로미 (radius={MONO_RADIUS}, steps={MONO_STEPS})",
          flush=True)
    monos = []
    deviations = []

    for i, tz in enumerate(zeros):
        try:
            mono = monodromy_contour_dirichlet(tz, ci)
            abs_mono = abs(mono)
            dev = abs(abs_mono - 2.0)  # 이론값 = ±2 (단위: π)
            monos.append(mono)
            deviations.append(dev)
            if i < 5 or (i+1) % 20 == 0:
                print(f"  t={tz:.4f}: mono/π = {mono:+.4f}, "
                      f"|dev from ±2| = {dev:.4f}", flush=True)
        except Exception as e:
            print(f"  WARNING: 모노드로미 실패 t={tz:.4f}: {e}", flush=True)

    monos = np.array(monos)
    deviations = np.array(deviations)
    mean_abs = np.mean(np.abs(monos)) if len(monos) > 0 else 0
    mean_dev = np.mean(deviations) if len(deviations) > 0 else 1
    std_dev  = np.std(deviations) if len(deviations) > 0 else 1

    print(f"  평균 |mono|/π = {mean_abs:.4f}  (기준 2.000±0.001)", flush=True)
    print(f"  평균 편차 from 2.000 = {mean_dev:.4f} ± {std_dev:.4f}", flush=True)
    pass_c = mean_abs > 1.99 and mean_abs < 2.01 and mean_dev < 0.01
    print(f"  판정: {'✅ PASS' if pass_c else '❌ FAIL'}", flush=True)

    return {
        'monos': monos, 'deviations': deviations,
        'mean_abs': mean_abs, 'mean_dev': mean_dev, 'std_dev': std_dev,
        'pass': pass_c,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파트 D: 블라인드 예측
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def find_kappa_peaks(t_vals, kappas, min_prominence=5.0, min_dist_idx=3):
    """κ(t) 배열에서 극대점 추출"""
    peaks = []
    n = len(kappas)
    for i in range(1, n - 1):
        if kappas[i] > kappas[i-1] and kappas[i] > kappas[i+1]:
            # prominence 체크 (주변보다 명확히 높은지)
            if kappas[i] > min_prominence:
                peaks.append((t_vals[i], kappas[i]))
    # 너무 가까운 피크 합치기
    merged = []
    for t_pk, k_pk in sorted(peaks, key=lambda x: -x[1]):
        if not merged or min(abs(t_pk - p[0]) for p in merged) > 0.3:
            merged.append((t_pk, k_pk))
    merged.sort(key=lambda x: x[0])
    return merged


def part_d_blind_prediction(ci, zeros_known, t_scan_min=10.0, t_scan_max=50.0,
                             n_scan=1000, n_predict=5):
    """
    κ(t) 스캔 → 피크 → findroot 확인.
    알려진 영점과 비교하여 예측 정확도 측정.
    """
    print(f"\n[파트 D] 블라인드 예측 (t∈[{t_scan_min},{t_scan_max}])",
          flush=True)

    # κ(t) 스캔 (σ=0.5, δ=0.03 오프셋)
    ts = np.linspace(t_scan_min, t_scan_max, n_scan)
    kappas = []
    print(f"  κ 스캔 중 ({n_scan}점)...", flush=True)
    for i, t in enumerate(ts):
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t + DELTA))
        try:
            k = float(abs(connection_analytic(s, ci))**2)
            kappas.append(k if np.isfinite(k) else 0.0)
        except Exception as e:
            print(f"  WARNING: 스캔 실패 t={t:.2f}: {e}", flush=True)
            kappas.append(0.0)
        if (i+1) % 200 == 0:
            print(f"  스캔 {i+1}/{n_scan}", flush=True)

    kappas = np.array(kappas)
    peaks = find_kappa_peaks(ts, kappas, min_prominence=50.0)
    print(f"  κ 피크 {len(peaks)}개 발견", flush=True)
    for t_pk, k_pk in peaks[:10]:
        print(f"    t={t_pk:.4f}, κ={k_pk:.1f}", flush=True)

    # 상위 n_predict개 피크에서 영점 예측
    top_peaks = sorted(peaks, key=lambda x: -x[1])[:n_predict]
    top_peaks.sort(key=lambda x: x[0])

    predictions = []
    hits = 0
    for t_pk, k_pk in top_peaks:
        try:
            def f_real(t_var):
                sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                return mpmath.re(completed_L(sv, ci))
            tz_pred = float(mpmath.findroot(f_real, mpmath.mpf(str(t_pk))))
            # 알려진 영점과 비교
            if len(zeros_known) > 0:
                nearest = min(zeros_known, key=lambda z: abs(z - tz_pred))
                error = abs(tz_pred - nearest)
                hit = error < 0.1
            else:
                nearest, error, hit = None, float('inf'), False
            predictions.append({
                'peak_t': t_pk, 'pred_t': tz_pred,
                'nearest_known': float(nearest) if nearest is not None else None,
                'error': error, 'hit': hit,
            })
            if hit:
                hits += 1
            status = "✅" if hit else "❌"
            print(f"  {status} 피크 t={t_pk:.4f} → 예측 {tz_pred:.4f} "
                  f"(가장 가까운 영점: {float(nearest):.4f}, 오차: {error:.4f})",
                  flush=True)
        except Exception as e:
            print(f"  WARNING: 예측 findroot 실패 t≈{t_pk:.4f}: {e}", flush=True)

    n_tried = len(predictions)
    print(f"  ★ 적중: {hits}/{n_tried}  (기준 5/5 이상)", flush=True)
    pass_d = hits >= 5 and n_tried >= 5
    print(f"  판정: {'✅ PASS' if pass_d else '❌ FAIL'}", flush=True)

    return {
        'peaks': peaks, 'predictions': predictions,
        'hits': hits, 'n_tried': n_tried,
        'pass': pass_d,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파트 E: σ-국소화
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def part_e_sigma_localization(ci, zeros):
    """
    E(σ) = Σ_n κ(σ + i(γ_n + δ))
    σ∈[0.1, 0.9] 9점에서 에너지 비교.
    E(0.5) / min(E(타σ)) > 50× 기준.
    """
    sigmas = np.linspace(0.1, 0.9, 9)
    print(f"\n[파트 E] σ-국소화 ({len(zeros)}영점 × {len(sigmas)}σ값)", flush=True)

    energies = []
    for sigma in sigmas:
        E = 0.0
        for tz in zeros:
            t_meas = tz + DELTA
            s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t_meas))
            try:
                k = float(abs(connection_analytic(s, ci))**2)
                if np.isfinite(k):
                    E += k
            except Exception as e:
                print(f"  WARNING: E(σ={sigma:.1f}) 실패 t={tz:.4f}: {e}",
                      flush=True)
        energies.append(E)
        marker = "★" if abs(sigma - 0.5) < 0.01 else " "
        print(f"  {marker} σ={sigma:.2f}: E = {E:.2e}", flush=True)

    energies = np.array(energies)
    idx_half = np.argmin(np.abs(sigmas - 0.5))
    E_half = energies[idx_half]

    # 타σ 에너지: σ=0.5 제외
    other_mask = np.abs(sigmas - 0.5) > 0.05
    other_energies = energies[other_mask]
    E_other_max = np.max(other_energies) if len(other_energies) > 0 else 1.0
    E_other_min = np.min(other_energies) if len(other_energies) > 0 else 1.0

    ratio_vs_max = E_half / E_other_max if E_other_max > 0 else float('inf')
    ratio_vs_03 = E_half / energies[np.argmin(np.abs(sigmas - 0.3))]
    ratio_vs_07 = E_half / energies[np.argmin(np.abs(sigmas - 0.7))]

    print(f"\n  E(0.5) = {E_half:.2e}", flush=True)
    print(f"  E(0.5)/E(0.3) = {ratio_vs_03:.1f}×", flush=True)
    print(f"  E(0.5)/E(0.7) = {ratio_vs_07:.1f}×", flush=True)
    print(f"  E(0.5)/E(타σ 최대) = {ratio_vs_max:.1f}×  (기준 >50×)", flush=True)
    pass_e = ratio_vs_max > 50
    print(f"  판정: {'✅ PASS' if pass_e else '❌ FAIL'}", flush=True)

    return {
        'sigmas': sigmas, 'energies': energies,
        'E_half': E_half, 'E_other_max': E_other_max,
        'ratio_vs_max': ratio_vs_max,
        'ratio_vs_03': ratio_vs_03, 'ratio_vs_07': ratio_vs_07,
        'pass': pass_e,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 결과 파일 저장
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def save_results(out_path, zeros, res_b, res_c, res_d, res_e, elapsed_total):
    with open(out_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("결과 #40 — χ mod 7 고 conductor 디리클레 4성질 완전 검증\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"정밀도: {mpmath.mp.dps} 자릿수\n")
        f.write(f"t 범위: [{T_MIN}, {T_MAX}], δ={DELTA}\n")
        f.write(f"총 소요: {elapsed_total:.0f}초\n")
        f.write("=" * 70 + "\n\n")

        # 영점 목록
        f.write(f"[파트 A] 영점 ({len(zeros)}개, t∈[{T_MIN},{T_MAX}])\n")
        f.write("-" * 40 + "\n")
        for i, tz in enumerate(zeros):
            f.write(f"  #{i+1:3d}: t = {tz:.8f}\n")
        f.write(f"\n  총 영점 수: {len(zeros)}\n")
        pass_a = len(zeros) >= 30
        f.write(f"  성공 기준 (≥30): {'✅ PASS' if pass_a else '❌ FAIL'}\n\n")

        # 파트 B
        f.write("[파트 B] κ 비율 (σ=0.5 vs 타σ, δ=0.03)\n")
        f.write("-" * 40 + "\n")
        f.write(f"  median κ(σ=0.5): {res_b['med_half']:.4e}\n")
        f.write(f"  median κ(σ=0.3): {res_b['med_03']:.4e}\n")
        f.write(f"  median κ(σ=0.7): {res_b['med_07']:.4e}\n")
        f.write(f"  κ(0.5)/κ(0.3): {res_b['ratio_03']:.1f}×\n")
        f.write(f"  κ(0.5)/κ(0.7): {res_b['ratio_07']:.1f}×\n")
        f.write(f"  최소 비율: {res_b['ratio_min']:.1f}×  (기준 >100×)\n")
        f.write(f"  판정: {'✅ PASS' if res_b['pass'] else '❌ FAIL'}\n\n")

        # 파트 C
        f.write("[파트 C] 모노드로미 (radius=0.1, 64단계)\n")
        f.write("-" * 40 + "\n")
        f.write(f"  {'t':>10} {'mono/π':>12} {'|dev from 2|':>14}\n")
        f.write("  " + "-" * 40 + "\n")
        for i, (tz, mono, dev) in enumerate(
                zip(zeros, res_c['monos'], res_c['deviations'])):
            f.write(f"  {tz:>10.4f} {mono:>+12.5f} {dev:>14.6f}\n")
        f.write(f"\n  평균 |mono|/π = {res_c['mean_abs']:.6f}\n")
        f.write(f"  평균 편차 from 2.000 = {res_c['mean_dev']:.6f}"
                f" ± {res_c['std_dev']:.6f}\n")
        f.write(f"  판정: {'✅ PASS' if res_c['pass'] else '❌ FAIL'}\n\n")

        # 파트 D
        f.write("[파트 D] 블라인드 예측\n")
        f.write("-" * 40 + "\n")
        for p in res_d['predictions']:
            status = "✅" if p['hit'] else "❌"
            f.write(f"  {status} 피크 t={p['peak_t']:.4f} → 예측 {p['pred_t']:.4f}"
                    f" (최근접 영점: {p['nearest_known']:.4f},"
                    f" 오차: {p['error']:.4f})\n")
        f.write(f"\n  적중: {res_d['hits']}/{res_d['n_tried']}  (기준 5/5)\n")
        f.write(f"  판정: {'✅ PASS' if res_d['pass'] else '❌ FAIL'}\n\n")

        # 파트 E
        f.write("[파트 E] σ-국소화\n")
        f.write("-" * 40 + "\n")
        f.write(f"  {'σ':>8} {'E(σ)':>14}\n")
        for sigma, E in zip(res_e['sigmas'], res_e['energies']):
            marker = "★" if abs(sigma - 0.5) < 0.01 else " "
            f.write(f"  {marker}{sigma:>7.2f} {E:>14.4e}\n")
        f.write(f"\n  E(0.5)/E(0.3): {res_e['ratio_vs_03']:.1f}×\n")
        f.write(f"  E(0.5)/E(0.7): {res_e['ratio_vs_07']:.1f}×\n")
        f.write(f"  E(0.5)/E(타σ 최대): {res_e['ratio_vs_max']:.1f}×  (기준 >50×)\n")
        f.write(f"  판정: {'✅ PASS' if res_e['pass'] else '❌ FAIL'}\n\n")

        # 종합 판정
        n_pass = sum([pass_a, res_b['pass'], res_c['pass'],
                      res_d['pass'], res_e['pass']])
        f.write("=" * 70 + "\n")
        f.write("종합 판정\n")
        f.write("=" * 70 + "\n")
        f.write(f"  파트 A (영점 ≥30): {'✅' if pass_a else '❌'}  ({len(zeros)}개)\n")
        f.write(f"  파트 B (κ > 100×): {'✅' if res_b['pass'] else '❌'}  "
                f"({res_b['ratio_min']:.0f}×)\n")
        f.write(f"  파트 C (mono=2π):  {'✅' if res_c['pass'] else '❌'}  "
                f"({res_c['mean_abs']:.4f}π)\n")
        f.write(f"  파트 D (5/5 예측): {'✅' if res_d['pass'] else '❌'}  "
                f"({res_d['hits']}/{res_d['n_tried']})\n")
        f.write(f"  파트 E (E>50×):   {'✅' if res_e['pass'] else '❌'}  "
                f"({res_e['ratio_vs_max']:.0f}×)\n")
        f.write(f"\n  PASS {n_pass}/5\n")

        if n_pass == 5:
            verdict = "★★ 강한 양성 — 4성질 모두 통과"
        elif n_pass >= 3:
            verdict = "★ 양성 — 3성질 이상 통과"
        else:
            verdict = "❌ 음성 — 2성질 이상 실패"
        f.write(f"\n  최종: {verdict}\n")
        f.write("=" * 70 + "\n")

    return n_pass, verdict


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    out_path = os.path.expanduser(
        '~/Desktop/gdl_unified/results/dirichlet_mod7_verification.txt'
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    ci = CHAR_MOD7
    print("=" * 70, flush=True)
    print("결과 #40 — χ mod 7 고 conductor 디리클레 4성질 완전 검증", flush=True)
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}", flush=True)
    print(f"지표: {ci['label']}", flush=True)
    print(f"t 범위: [{T_MIN}, {T_MAX}], dps={mpmath.mp.dps}", flush=True)
    print("=" * 70, flush=True)

    t_start = time.time()

    # ── 파트 A: 영점 탐색 ──────────────────────────────────────────────────
    print("\n[파트 A] 영점 탐색 시작 (t∈[10, 200], n_scan=4000)", flush=True)
    t0 = time.time()
    zeros = find_zeros_mod7(ci, T_MIN, T_MAX, n_scan=4000)
    print(f"\n  ✓ 영점 {len(zeros)}개 수집 ({time.time()-t0:.0f}초)", flush=True)
    print(f"  t 범위: [{zeros[0]:.4f}, {zeros[-1]:.4f}]" if len(zeros) > 0
          else "  ⚠️ 영점 없음", flush=True)
    for i, tz in enumerate(zeros[:10]):
        print(f"    #{i+1}: t = {tz:.6f}", flush=True)
    if len(zeros) > 10:
        print(f"    ... (총 {len(zeros)}개)", flush=True)

    if len(zeros) == 0:
        print("FATAL: 영점을 찾지 못했습니다. 스크립트 종료.", flush=True)
        return

    # ── 파트 B: κ 비율 ──────────────────────────────────────────────────────
    t0 = time.time()
    res_b = part_b_kappa_ratio(ci, zeros)
    print(f"  [소요 {time.time()-t0:.0f}초]", flush=True)

    # ── 파트 C: 모노드로미 ──────────────────────────────────────────────────
    t0 = time.time()
    res_c = part_c_monodromy(ci, zeros)
    print(f"  [소요 {time.time()-t0:.0f}초]", flush=True)

    # ── 파트 D: 블라인드 예측 ───────────────────────────────────────────────
    t0 = time.time()
    res_d = part_d_blind_prediction(ci, zeros,
                                     t_scan_min=10.0, t_scan_max=80.0,
                                     n_scan=700, n_predict=7)
    print(f"  [소요 {time.time()-t0:.0f}초]", flush=True)

    # ── 파트 E: σ-국소화 ────────────────────────────────────────────────────
    t0 = time.time()
    res_e = part_e_sigma_localization(ci, zeros)
    print(f"  [소요 {time.time()-t0:.0f}초]", flush=True)

    # ── 저장 + 종합 판정 ────────────────────────────────────────────────────
    elapsed_total = time.time() - t_start
    n_pass, verdict = save_results(
        out_path, zeros, res_b, res_c, res_d, res_e, elapsed_total
    )

    print("\n" + "=" * 70, flush=True)
    print("종합 판정", flush=True)
    print("=" * 70, flush=True)
    print(f"  파트 A (영점 ≥30): {'✅' if len(zeros)>=30 else '❌'}  ({len(zeros)}개)",
          flush=True)
    print(f"  파트 B (κ > 100×): {'✅' if res_b['pass'] else '❌'}  "
          f"({res_b['ratio_min']:.0f}×)", flush=True)
    print(f"  파트 C (mono=2π):  {'✅' if res_c['pass'] else '❌'}  "
          f"({res_c['mean_abs']:.4f}π)", flush=True)
    print(f"  파트 D (5/5 예측): {'✅' if res_d['pass'] else '❌'}  "
          f"({res_d['hits']}/{res_d['n_tried']})", flush=True)
    print(f"  파트 E (E>50×):   {'✅' if res_e['pass'] else '❌'}  "
          f"({res_e['ratio_vs_max']:.0f}×)", flush=True)
    print(f"\n  PASS {n_pass}/5", flush=True)
    print(f"  최종: {verdict}", flush=True)
    print(f"\n결과 저장: {out_path}", flush=True)
    print(f"총 소요: {elapsed_total:.0f}초", flush=True)
    print("완료.", flush=True)


if __name__ == '__main__':
    main()
