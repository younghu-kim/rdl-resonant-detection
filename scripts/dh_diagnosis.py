#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] DH 함수 다발 진단 — on-critical vs off-critical 비교
=============================================================================
확장 탐색에서 발견한 off-critical 영점 4개 + on-critical 영점 10개에 대해
다발 진단 수행: κ, 모노드로미, peak_σ, FWHM.

결과를 results/davenport_heilbronn_adversarial.txt에 추가(append).
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

mpmath.mp.dps = 80
START = time.time()

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
RESULT_FILE = os.path.join(RESULT_DIR, 'davenport_heilbronn_adversarial.txt')

outf = open(RESULT_FILE, 'a')  # 추가 모드

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()


log()
log("=" * 70)
log("[결과 #30 보충] Off-critical 영점 다발 진단")
log("=" * 70)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()


# DH 함수 정의
CHI_MOD5 = [0, 1, 1j, -1j, -1]
CHI_MOD5_BAR = [0, 1, -1j, 1j, -1]
sqrt5 = mpmath.sqrt(5)
kappa_DH = (mpmath.sqrt(10 - 2*sqrt5) - 2) / (sqrt5 - 1)
coeff_chi = (1 - 1j * kappa_DH) / 2
coeff_chi_bar = (1 + 1j * kappa_DH) / 2


def dh_func(s):
    s = mpmath.mpc(s)
    return coeff_chi * mpmath.dirichlet(s, CHI_MOD5) + coeff_chi_bar * mpmath.dirichlet(s, CHI_MOD5_BAR)


def dh_derivative(s, h=None):
    if h is None:
        h = mpmath.mpf(1) / mpmath.mpf(10**20)
    s = mpmath.mpc(s)
    return (dh_func(s + h) - dh_func(s - h)) / (2 * h)


def dh_monodromy(s_center, radius=0.5, n_steps=64):
    """폐곡선 적분으로 모노드로미 계산"""
    total_delta = 0.0
    prev_arg = None
    for k in range(n_steps + 1):
        theta = 2 * np.pi * k / n_steps
        s = s_center + radius * mpmath.exp(1j * theta)
        val = dh_func(s)
        if abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
            continue
        curr_arg = float(mpmath.arg(val))
        if prev_arg is not None:
            delta = curr_arg - prev_arg
            while delta > np.pi: delta -= 2 * np.pi
            while delta < -np.pi: delta += 2 * np.pi
            total_delta += delta
        prev_arg = curr_arg
    return total_delta


# 영점 목록
on_critical = [
    (0.5, 5.094160),
    (0.5, 8.939914),
    (0.5, 12.133545),
    (0.5, 14.404003),
    (0.5, 17.130239),
    (0.5, 19.308800),
    (0.5, 22.159708),
    (0.5, 23.345370),
    (0.5, 26.094967),
    (0.5, 27.923799),
]

off_critical = [
    (0.808517, 85.699348),
    (0.650830, 114.163343),
    (0.574356, 166.479306),
    (0.724258, 176.702461),
]

DELTA = 0.03

def diagnose(sigma_0, t_0, label):
    """영점 다발 진단: κ, mono, peak_σ, FWHM"""
    log(f"\n--- {label}: σ={sigma_0:.6f}, t={t_0:.6f} ---")

    s_center = mpmath.mpf(str(sigma_0)) + 1j * mpmath.mpf(str(t_0))

    # 1. 곡률 (δ 오프셋)
    s_off = mpmath.mpf(str(sigma_0 + DELTA)) + 1j * mpmath.mpf(str(t_0))
    f_val = dh_func(s_off)
    if abs(f_val) > mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        conn = dh_derivative(s_off) / f_val
        kappa = float(abs(conn)**2)
    else:
        kappa = float('inf')
    log(f"  κ = {kappa:.2f}")

    # 2. 모노드로미
    mono = dh_monodromy(s_center, radius=0.5, n_steps=64)
    mono_pi = mono / np.pi
    log(f"  mono = {mono:.6f} ({mono_pi:.4f}π)")

    # 3. σ-국소화 프로파일
    n_prof = 81
    sig_lo = max(0.05, sigma_0 - 0.6)
    sig_hi = min(1.5, sigma_0 + 0.6)
    sigmas = np.linspace(sig_lo, sig_hi, n_prof)
    kappa_prof = np.zeros(n_prof)

    for i, sig in enumerate(sigmas):
        s_test = mpmath.mpf(str(sig)) + 1j * mpmath.mpf(str(t_0))
        f_v = dh_func(s_test)
        if abs(f_v) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
            kappa_prof[i] = 1e10
            continue
        f_d = dh_derivative(s_test)
        kappa_prof[i] = float(abs(f_d / f_v)**2)

    # 발산 cap
    kappa_cap = np.where(kappa_prof > 1e9, np.nan, kappa_prof)
    valid = ~np.isnan(kappa_cap)

    if np.sum(valid) > 5:
        peak_idx = np.nanargmax(kappa_cap)
        peak_sigma = sigmas[peak_idx]
        peak_kappa = kappa_cap[peak_idx]

        half_max = peak_kappa / 2
        above = (kappa_cap >= half_max) & valid
        above_idx = np.where(above)[0]
        fwhm = sigmas[above_idx[-1]] - sigmas[above_idx[0]] if len(above_idx) >= 2 else 0.015
    else:
        peak_sigma = sigma_0
        fwhm = 0.0

    log(f"  peak_σ = {peak_sigma:.4f} (영점 σ₀={sigma_0:.4f}, 차이={abs(peak_sigma-sigma_0):.4f})")
    log(f"  FWHM = {fwhm:.4f}")

    return {
        'sigma_0': sigma_0, 't_0': t_0,
        'kappa': kappa, 'mono_pi': mono_pi,
        'peak_sigma': peak_sigma, 'fwhm': fwhm,
        'label': label,
    }


# 진단 실행
log("===== On-critical 영점 (σ=1/2) 진단 =====")
on_diag = []
for i, (sig, t) in enumerate(on_critical):
    d = diagnose(sig, t, f"ON #{i+1}")
    on_diag.append(d)

log()
log("===== Off-critical 영점 (σ≠1/2) 진단 =====")
off_diag = []
for i, (sig, t) in enumerate(off_critical):
    d = diagnose(sig, t, f"OFF #{i+1}")
    off_diag.append(d)


# 통계
log()
log("=" * 70)
log("최종 비교 테이블")
log("=" * 70)
log()

def stats(vals):
    a = np.array(vals)
    return np.mean(a), np.std(a)

on_k = [d['kappa'] for d in on_diag if d['kappa'] < 1e9]
on_m = [d['mono_pi'] for d in on_diag]
on_p = [d['peak_sigma'] for d in on_diag]
on_f = [d['fwhm'] for d in on_diag]

off_k = [d['kappa'] for d in off_diag if d['kappa'] < 1e9]
off_m = [d['mono_pi'] for d in off_diag]
off_p = [d['peak_sigma'] for d in off_diag]
off_f = [d['fwhm'] for d in off_diag]

log(f"{'메트릭':<12} {'On-critical (N=10)':<28} {'Off-critical (N=4)':<28} {'판별'}")
log("-" * 90)

mk_on, sk_on = stats(on_k) if on_k else (0,0)
mk_off, sk_off = stats(off_k) if off_k else (0,0)
log(f"{'κ':<12} {mk_on:>10.2f} ± {sk_on:>8.2f}     {mk_off:>10.2f} ± {sk_off:>8.2f}     {'—'}")

mm_on, sm_on = stats(on_m)
mm_off, sm_off = stats(off_m)
log(f"{'mono/π':<12} {mm_on:>10.4f} ± {sm_on:>8.4f}     {mm_off:>10.4f} ± {sm_off:>8.4f}     {'동일' if abs(mm_on-mm_off)<0.1 else '★ 다름'}")

mp_on, sp_on = stats(on_p)
mp_off, sp_off = stats(off_p)
sep = abs(mp_off - mp_on)
log(f"{'peak_σ':<12} {mp_on:>10.4f} ± {sp_on:>8.4f}     {mp_off:>10.4f} ± {sp_off:>8.4f}     {'★★ 분리! (Δ='+f'{sep:.4f})' if sep > 0.03 else '겹침'}")

mf_on, sf_on = stats(on_f)
mf_off, sf_off = stats(off_f)
ratio = mf_off / mf_on if mf_on > 0 else 0
log(f"{'FWHM':<12} {mf_on:>10.4f} ± {sf_on:>8.4f}     {mf_off:>10.4f} ± {sf_off:>8.4f}     {'★★ 넓어짐! (×'+f'{ratio:.1f})' if ratio > 1.3 else '유사'}")


# 개별 테이블
log()
log("개별 영점 진단 상세:")
log(f"{'#':<5} {'유형':<8} {'σ₀':<10} {'t₀':<12} {'κ':<12} {'mono/π':<10} {'peak_σ':<10} {'|Δσ|':<8} {'FWHM':<10}")
log("-" * 90)

all_d = [(d, 'ON') for d in on_diag] + [(d, 'OFF') for d in off_diag]
for i, (d, typ) in enumerate(all_d):
    k_str = f"{d['kappa']:.2f}" if d['kappa'] < 1e9 else "∞"
    delta_sig = abs(d['peak_sigma'] - d['sigma_0'])
    log(f"{i+1:<5} {typ:<8} {d['sigma_0']:<10.6f} {d['t_0']:<12.6f} {k_str:<12} {d['mono_pi']:<10.4f} {d['peak_sigma']:<10.4f} {delta_sig:<8.4f} {d['fwhm']:<10.4f}")


# 판정
log()
log("=" * 70)
log("최종 판정")
log("=" * 70)
log()

if sep > 0.03 and ratio > 1.3:
    verdict = "강한 양성"
elif sep > 0.03 or ratio > 1.3:
    verdict = "양성"
else:
    verdict = "음성"

log(f"1. Off-critical 영점 발견: 4개 (σ∈[0.57, 0.81], t∈[85, 177])")
log(f"   - σ=0.808517, t=85.699 (Spira 1968과 일치)")
log(f"   - σ=0.650830, t=114.163")
log(f"   - σ=0.574356, t=166.479")
log(f"   - σ=0.724258, t=176.702")
log()
log(f"2. peak_σ 분리: {sep:.4f} {'(>0.03 ★★)' if sep > 0.03 else '(<0.03)'}")
log(f"3. FWHM 비율: {ratio:.2f} {'(>1.3 ★★)' if ratio > 1.3 else '(<1.3)'}")
log()
log(f"판정: **{verdict}**")
log()

if verdict in ["강한 양성", "양성"]:
    log("해석:")
    log("  다발 프레임워크의 곡률 프로파일(peak_σ, FWHM)은")
    log("  on-critical 영점(σ=1/2)과 off-critical 영점(σ≠1/2)을")
    log("  정량적으로 구별한다.")
    log()
    log("  → peak_σ가 실제 영점 위치로 이동: 프레임워크의 σ-국소화 예측 확인")
    log("  → FWHM 증가: off-critical 영점에서 곡률 집중이 약화")
    log("  → 이것은 프레임워크의 falsifiability를 확립:")
    log("    'RH를 위반하는 함수에서는 다발 서명이 달라진다'")
else:
    log("해석: 프레임워크가 on-critical과 off-critical을 구별하지 못함.")
    log("  이는 프레임워크의 한계로 논문에 기재해야 함.")

log()
log(f"DH 함수 주의: 오일러 곱 없음 (L-함수 아님). 함수방정식만 만족.")
log(f"총 소요: {(time.time()-START)/60:.1f}분")
log(f"종료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

outf.close()
print(f"\n결과 추가: {RESULT_FILE}")
