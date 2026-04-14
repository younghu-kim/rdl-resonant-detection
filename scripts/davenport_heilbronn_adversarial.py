#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #30 — Davenport-Heilbronn 함수 적대적 검증
=============================================================================
Davenport-Heilbronn (DH) 함수:
  f(s) = [(1-iκ)/2] · L(s,χ) + [(1+iκ)/2] · L(s,χ̄)

  χ: mod 5 원시 지표 (차수 4)
    χ(0)=0, χ(1)=1, χ(2)=i, χ(3)=-i, χ(4)=-1
  κ = (√(10-2√5) - 2) / (√5 - 1)

  f(s)는 함수방정식을 만족하지만 오일러 곱이 없음 (L-함수 아님).
  **σ>1/2에 무한히 많은 영점 존재** (Davenport & Heilbronn, 1936).

목표: off-critical 영점의 다발 서명(κ, mono, peak_σ, FWHM)이
       on-critical 영점과 정량적으로 다른지 검증 → falsifiability 확립.

진단:
  - 접속: L_f = f'/f (대수 미분)
  - 곡률: κ = |f'/f|² at δ=0.03 오프셋
  - 모노드로미: ∮ d(arg f) 폐곡선 적분 (radius=0.5, 64단계)
  - σ-국소화: κ(σ) 프로파일 → FWHM, peak_σ

출력: results/davenport_heilbronn_adversarial.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

mpmath.mp.dps = 80

START_TIME = time.time()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 출력 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'davenport_heilbronn_adversarial.txt')

outf = open(RESULT_FILE, 'w')

def log(msg=''):
    """콘솔 + 파일 동시 출력"""
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()


log("=" * 70)
log("[결과 #30] Davenport-Heilbronn 함수 적대적 검증")
log("=" * 70)
log(f"시작 시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"mpmath dps = {mpmath.mp.dps}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. Davenport-Heilbronn 함수 구현
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("━━━ 1. DH 함수 구현 ━━━")

# χ mod 5 (차수 4, 홀수 지표: χ(-1)=χ(4)=-1 → a=1)
CHI_MOD5 = [0, 1, 1j, -1j, -1]  # χ(0), χ(1), χ(2), χ(3), χ(4)
CHI_MOD5_BAR = [0, 1, -1j, 1j, -1]  # χ̄

# κ 상수
sqrt5 = mpmath.sqrt(5)
kappa_DH = (mpmath.sqrt(10 - 2*sqrt5) - 2) / (sqrt5 - 1)
log(f"κ (DH 상수) = {float(kappa_DH):.10f}")

# DH 계수
coeff_chi = (1 - 1j * kappa_DH) / 2      # χ 항 계수
coeff_chi_bar = (1 + 1j * kappa_DH) / 2   # χ̄ 항 계수

log(f"coeff_chi     = {complex(coeff_chi)}")
log(f"coeff_chi_bar = {complex(coeff_chi_bar)}")
log()


def dh_func(s):
    """
    Davenport-Heilbronn 함수:
    f(s) = [(1-iκ)/2] · L(s,χ) + [(1+iκ)/2] · L(s,χ̄)
    """
    s = mpmath.mpc(s)
    L_chi = mpmath.dirichlet(s, CHI_MOD5)
    L_chi_bar = mpmath.dirichlet(s, CHI_MOD5_BAR)
    return coeff_chi * L_chi + coeff_chi_bar * L_chi_bar


def dh_func_derivative(s, h=None):
    """f'(s) 수치 미분 (중심차분)"""
    if h is None:
        h = mpmath.mpf(1) / mpmath.mpf(10**20)
    s = mpmath.mpc(s)
    return (dh_func(s + h) - dh_func(s - h)) / (2 * h)


def dh_connection(s, delta=0.03):
    """
    접속 L_f = f'/f
    δ 오프셋: 영점 위 직접 계산 방지
    """
    s = mpmath.mpc(s)
    s_off = s + mpmath.mpf(str(delta))  # σ 방향 오프셋
    f_val = dh_func(s_off)
    if abs(f_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e10, 0)
    f_deriv = dh_func_derivative(s_off)
    return f_deriv / f_val


def dh_curvature(s, delta=0.03):
    """곡률 κ = |f'/f|² at δ 오프셋"""
    L = dh_connection(s, delta)
    return float(abs(L)**2)


def dh_monodromy(s_center, radius=0.5, n_steps=64):
    """
    s_center 주위 반지름 radius 원에서 모노드로미.
    ∮ d(arg f) 폐곡선 적분.
    영점이 원 안에 있으면 ≈±2π.
    """
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
            while delta > np.pi:
                delta -= 2 * np.pi
            while delta < -np.pi:
                delta += 2 * np.pi
            total_delta += delta

        prev_arg = curr_arg

    return total_delta


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. 함수방정식 수치 검증
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("━━━ 2. 함수방정식 검증 ━━━")
log()

# 완비 DH: Φ(s) = (5/π)^{s/2} Γ((s+1)/2) f(s)  [a=1, 홀수 지표]
# 함수방정식: Φ(s) = ε · Φ̄(1-s̄) 형태
# 여기서 ε = i·(-1)^a · (5)^{1/2} / τ(χ)·...
# 실제로는 f(s)와 f(1-s) 관계를 수치적으로 확인

def completed_dh(s):
    """완비 DH: Φ(s) = (5/π)^{s/2} Γ((s+1)/2) f(s)"""
    s = mpmath.mpc(s)
    prefactor = mpmath.power(5 / mpmath.pi, s / 2)
    gamma_val = mpmath.gamma((s + 1) / 2)
    return prefactor * gamma_val * dh_func(s)


# 함수방정식 검증: Φ(s) / Φ(1-conj(s)) 의 비율이 상수
test_points = [
    mpmath.mpc(0.3, 5.0),
    mpmath.mpc(0.7, 10.0),
    mpmath.mpc(0.5, 15.0),
    mpmath.mpc(0.2, 20.0),
]

log("함수방정식 검증: Φ(s) / conj(Φ(1-conj(s))) 비율")
ratios = []
for s in test_points:
    Phi_s = completed_dh(s)
    s_dual = 1 - mpmath.conj(s)
    Phi_dual = completed_dh(s_dual)
    Phi_dual_conj = mpmath.conj(Phi_dual)

    if abs(Phi_dual_conj) > mpmath.mpf('1e-50'):
        ratio = Phi_s / Phi_dual_conj
        ratios.append(ratio)
        log(f"  s={complex(s):.1f}: Φ(s)/conj(Φ(1-s̄)) = {complex(ratio):.10f}")
    else:
        log(f"  s={complex(s):.1f}: Φ(1-s̄) ≈ 0, 스킵")

if len(ratios) >= 2:
    # 비율이 모두 같은 상수인지 확인
    r0 = ratios[0]
    all_same = all(abs(r - r0) < mpmath.mpf('1e-10') for r in ratios)
    log(f"비율 일관성: {'✅ 일치' if all_same else '⚠️ 불일치'} (ε ≈ {complex(r0):.6f})")
else:
    log("⚠️ 검증 포인트 부족")

log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. 영점 탐색
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("━━━ 3. 영점 탐색 ━━━")
log()

# 3a. 임계선 영점 (σ=1/2)
log("--- 3a. 임계선 영점 (σ=1/2, t∈[1,60]) ---")
t_scan = np.linspace(1.0, 60.0, 3000)  # 밀도 높은 스캔
on_critical_zeros = []

prev_val = None
prev_t = None
findroot_fail = 0
findroot_total = 0

for t in t_scan:
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
    val = dh_func(s)
    abs_val = float(abs(val))

    if prev_val is not None:
        # |f| 최소화 방식: |f|가 지역 최소인 곳 탐색
        # Re(f)와 Im(f) 부호 변화 모두 확인
        re_change = float(mpmath.re(prev_val)) * float(mpmath.re(val)) < 0
        im_change = float(mpmath.im(prev_val)) * float(mpmath.im(val)) < 0

        if re_change or im_change:
            findroot_total += 1
            try:
                def f_abs_min(t_var):
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                    fv = dh_func(sv)
                    return mpmath.re(fv)

                t_root_re = float(mpmath.findroot(f_abs_min, (prev_t, t)))
                sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_root_re))
                f_at_root = dh_func(sv)

                if abs(f_at_root) < mpmath.mpf('1e-8'):
                    if not any(abs(t_root_re - z[1]) < 0.05 for z in on_critical_zeros):
                        on_critical_zeros.append((0.5, t_root_re, float(abs(f_at_root))))
                else:
                    # Re=0인 곳에서 |f|가 작은지 확인, Im도 시도
                    try:
                        def f_im(t_var):
                            sv2 = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                            return mpmath.im(dh_func(sv2))

                        t_root_im = float(mpmath.findroot(f_im, (prev_t, t)))
                        sv2 = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_root_im))
                        f_at_root2 = dh_func(sv2)

                        if abs(f_at_root2) < mpmath.mpf('1e-8'):
                            if not any(abs(t_root_im - z[1]) < 0.05 for z in on_critical_zeros):
                                on_critical_zeros.append((0.5, t_root_im, float(abs(f_at_root2))))
                    except Exception:
                        pass
            except Exception as e:
                findroot_fail += 1

    prev_val = val
    prev_t = t

# 보충: |f| 직접 최소화 방식 (부호 변화 놓친 영점 포착)
log(f"부호 변화 기반 탐색: {len(on_critical_zeros)}개, findroot 실패: {findroot_fail}/{findroot_total}")

# |f| 스캔으로 추가 영점 탐색
abs_vals = []
for t in t_scan:
    s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
    abs_vals.append(float(abs(dh_func(s))))
abs_vals = np.array(abs_vals)

# 지역 최솟값에서 |f| < 1e-4인 곳 추가 탐색
for i in range(1, len(abs_vals) - 1):
    if abs_vals[i] < abs_vals[i-1] and abs_vals[i] < abs_vals[i+1] and abs_vals[i] < 1e-4:
        t_cand = t_scan[i]
        if not any(abs(t_cand - z[1]) < 0.05 for z in on_critical_zeros):
            try:
                def f_complex_re(t_var):
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                    return mpmath.re(dh_func(sv))
                def f_complex_im(t_var):
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                    return mpmath.im(dh_func(sv))

                # 2D findroot: Re(f)=0, Im(f)=0 동시
                t_root = float(mpmath.findroot(
                    lambda tv: abs(dh_func(mpmath.mpf('0.5') + 1j * mpmath.mpf(tv))),
                    t_cand
                ))
                sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_root))
                f_check = dh_func(sv)
                if abs(f_check) < mpmath.mpf('1e-8'):
                    if not any(abs(t_root - z[1]) < 0.05 for z in on_critical_zeros):
                        on_critical_zeros.append((0.5, t_root, float(abs(f_check))))
            except Exception as e:
                pass

on_critical_zeros.sort(key=lambda x: x[1])
log(f"임계선 영점 총 {len(on_critical_zeros)}개 발견")
for i, (sig, t, res) in enumerate(on_critical_zeros[:15]):
    log(f"  [{i+1}] σ={sig:.4f}, t={t:.6f}, |f|={res:.2e}")
if len(on_critical_zeros) > 15:
    log(f"  ... ({len(on_critical_zeros)-15}개 더)")
log()


# 3b. Off-critical 영점 (σ∈[0.5,1.0], t∈[1,60])
log("--- 3b. Off-critical 영점 탐색 (σ∈[0.5,1.0], t∈[1,60]) ---")
log("격자: n_sigma=40, n_t=300 (총 12,000점)")

n_sigma = 40
n_t = 300
sigma_range = np.linspace(0.52, 1.0, n_sigma)  # σ>0.5만 (임계선은 이미 탐색)
t_range = np.linspace(1.0, 60.0, n_t)

off_critical_candidates = []
grid_scan_count = 0

log("격자 스캔 시작...")
t_grid_start = time.time()

for i, sigma in enumerate(sigma_range):
    for j, t in enumerate(t_range):
        s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
        val = dh_func(s)
        abs_f = float(abs(val))
        grid_scan_count += 1

        if abs_f < 0.1:  # 후보점 수집 (임계값 넉넉하게)
            off_critical_candidates.append((sigma, t, abs_f))

    if (i + 1) % 10 == 0:
        elapsed = time.time() - t_grid_start
        log(f"  σ={sigma:.3f} 완료 ({i+1}/{n_sigma}), 후보 {len(off_critical_candidates)}개, 경과 {elapsed:.0f}초")

log(f"격자 스캔 완료: {grid_scan_count}점, 후보 {len(off_critical_candidates)}개")
log(f"격자 스캔 소요: {time.time()-t_grid_start:.1f}초")

# 후보를 |f| 기준 정렬 후 상위 점에서 findroot 정밀화
off_critical_candidates.sort(key=lambda x: x[2])
off_critical_zeros = []

log(f"상위 후보 100개에서 findroot 정밀화...")

for sigma_c, t_c, abs_c in off_critical_candidates[:100]:
    # 임계선 근처(|σ-0.5|<0.02)는 제외
    if abs(sigma_c - 0.5) < 0.02:
        continue

    # 이미 찾은 영점 근처 제외
    if any(abs(sigma_c - z[0]) < 0.01 and abs(t_c - z[1]) < 0.05 for z in off_critical_zeros):
        continue

    try:
        # 2D findroot: Re(f(σ+it))=0, Im(f(σ+it))=0 동시
        def f_system(sigma_var, t_var):
            sv = mpmath.mpf(sigma_var) + 1j * mpmath.mpf(t_var)
            fv = dh_func(sv)
            return (mpmath.re(fv), mpmath.im(fv))

        result = mpmath.findroot(f_system, (mpmath.mpf(str(sigma_c)), mpmath.mpf(str(t_c))))
        sigma_root = float(result[0])
        t_root = float(result[1])

        # 유효성 검사
        if sigma_root < 0.1 or sigma_root > 1.5 or t_root < 0.5 or t_root > 65:
            continue

        sv = mpmath.mpf(str(sigma_root)) + 1j * mpmath.mpf(str(t_root))
        f_check = dh_func(sv)
        residual = float(abs(f_check))

        if residual < 1e-15:
            # 중복 제거
            if not any(abs(sigma_root - z[0]) < 0.01 and abs(t_root - z[1]) < 0.05
                      for z in off_critical_zeros):
                off_critical_zeros.append((sigma_root, t_root, residual))
                log(f"  ★ Off-critical 영점: σ={sigma_root:.6f}, t={t_root:.6f}, |f|={residual:.2e}")
    except Exception as e:
        pass

# σ<0.5 쪽도 탐색 (함수방정식에 의한 대칭 확인)
log()
log("--- 3b-2. σ<0.5 쪽 추가 탐색 (σ∈[0.0,0.48]) ---")
sigma_range_low = np.linspace(0.0, 0.48, 20)
t_range_low = np.linspace(1.0, 60.0, 200)

for sigma in sigma_range_low:
    for t in t_range_low:
        s = mpmath.mpf(str(sigma)) + 1j * mpmath.mpf(str(t))
        val = dh_func(s)
        abs_f = float(abs(val))

        if abs_f < 0.1:
            try:
                def f_sys2(sv, tv):
                    sc = mpmath.mpf(sv) + 1j * mpmath.mpf(tv)
                    fv = dh_func(sc)
                    return (mpmath.re(fv), mpmath.im(fv))

                result = mpmath.findroot(f_sys2, (mpmath.mpf(str(sigma)), mpmath.mpf(str(t))))
                sr = float(result[0])
                tr = float(result[1])

                if sr < -0.5 or sr > 1.5 or tr < 0.5 or tr > 65:
                    continue

                sv = mpmath.mpf(str(sr)) + 1j * mpmath.mpf(str(tr))
                res = float(abs(dh_func(sv)))

                if res < 1e-15 and abs(sr - 0.5) > 0.02:
                    if not any(abs(sr - z[0]) < 0.01 and abs(tr - z[1]) < 0.05
                             for z in off_critical_zeros):
                        off_critical_zeros.append((sr, tr, res))
                        log(f"  ★ Off-critical 영점 (σ<0.5): σ={sr:.6f}, t={tr:.6f}, |f|={res:.2e}")
            except Exception:
                pass

off_critical_zeros.sort(key=lambda x: x[1])
log()
log(f"Off-critical 영점 총 {len(off_critical_zeros)}개 발견")
for i, (sig, t, res) in enumerate(off_critical_zeros):
    log(f"  [{i+1}] σ={sig:.6f}, t={t:.6f}, |f|={res:.2e}")

if len(off_critical_zeros) == 0:
    log("⚠️ Off-critical 영점 0개 — 탐색 로직 점검 필요")
    log("  DH 함수의 off-critical zeros는 수학적으로 보장되나,")
    log("  t∈[1,60] 범위에서 출현하지 않았을 수 있음.")
    log("  Titchmarsh (1986), Balanzario (2007) 참조: 첫 off-critical zero는 t≈6-30 범위 보고")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. 다발 진단
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("━━━ 4. 다발 진단 ━━━")
log()

DELTA = 0.03  # 오프셋

def diagnose_zero(sigma_0, t_0, label, delta=DELTA):
    """
    영점 (σ₀, t₀) 근처에서 다발 진단 수행.
    Returns dict with κ, mono, peak_σ, FWHM
    """
    log(f"--- 진단: {label} (σ={sigma_0:.6f}, t={t_0:.6f}) ---")

    s_center = mpmath.mpf(str(sigma_0)) + 1j * mpmath.mpf(str(t_0))

    # 4a. 곡률 κ (δ 오프셋)
    s_off = mpmath.mpf(str(sigma_0 + delta)) + 1j * mpmath.mpf(str(t_0))
    f_val = dh_func(s_off)
    f_deriv = dh_func_derivative(s_off)

    if abs(f_val) > mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        conn = f_deriv / f_val
        kappa = float(abs(conn)**2)
    else:
        kappa = float('inf')

    log(f"  κ (δ={delta}) = {kappa:.2f}")

    # 4b. 모노드로미 (폐곡선 적분)
    mono = dh_monodromy(s_center, radius=0.5, n_steps=64)
    mono_pi = mono / np.pi
    log(f"  mono = {mono:.6f} ({mono_pi:.4f}π)")

    # 4c. σ-국소화 프로파일 (peak_σ, FWHM)
    n_sigma_prof = 61
    sigma_min = max(0.05, sigma_0 - 0.5)
    sigma_max = min(1.5, sigma_0 + 0.5)
    sigmas = np.linspace(sigma_min, sigma_max, n_sigma_prof)
    kappa_profile = np.zeros(n_sigma_prof)

    for i, sig in enumerate(sigmas):
        s_test = mpmath.mpf(str(sig)) + 1j * mpmath.mpf(str(t_0))
        f_val_test = dh_func(s_test)

        if abs(f_val_test) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
            kappa_profile[i] = 1e10  # 영점 직접 → 매우 큰 값
            continue

        f_d = dh_func_derivative(s_test)
        conn_test = f_d / f_val_test
        kappa_profile[i] = float(abs(conn_test)**2)

    # 발산 제거 (1e9 이상 cap)
    kappa_capped = np.where(kappa_profile > 1e9, np.nan, kappa_profile)
    valid = ~np.isnan(kappa_capped)

    if np.sum(valid) > 5:
        peak_idx = np.nanargmax(kappa_capped)
        peak_sigma = sigmas[peak_idx]
        peak_kappa = kappa_capped[peak_idx]

        # FWHM 계산
        half_max = peak_kappa / 2
        above_half = kappa_capped >= half_max
        above_indices = np.where(above_half & valid)[0]

        if len(above_indices) >= 2:
            fwhm = sigmas[above_indices[-1]] - sigmas[above_indices[0]]
        else:
            fwhm = sigmas[1] - sigmas[0]  # 최소 해상도
    else:
        peak_sigma = sigma_0
        peak_kappa = kappa
        fwhm = 0.0

    log(f"  peak_σ = {peak_sigma:.4f} (영점 σ₀={sigma_0:.4f})")
    log(f"  FWHM = {fwhm:.4f}")
    log(f"  |peak_σ - σ₀| = {abs(peak_sigma - sigma_0):.4f}")
    log()

    return {
        'sigma_0': sigma_0,
        't_0': t_0,
        'kappa': kappa,
        'mono': mono,
        'mono_pi': mono_pi,
        'peak_sigma': peak_sigma,
        'fwhm': fwhm,
        'label': label,
    }


# 4-1. On-critical 영점 진단
log("===== On-critical 영점 다발 진단 =====")
log()

on_crit_diag = []
max_diag = min(10, len(on_critical_zeros))

if max_diag == 0:
    log("⚠️ 임계선 영점 없음 — 진단 불가")
else:
    for i in range(max_diag):
        sig, t, res = on_critical_zeros[i]
        d = diagnose_zero(sig, t, f"on-crit #{i+1}")
        on_crit_diag.append(d)

# 4-2. Off-critical 영점 진단
log("===== Off-critical 영점 다발 진단 =====")
log()

off_crit_diag = []

if len(off_critical_zeros) == 0:
    log("⚠️ Off-critical 영점 없음 — 진단 불가")
else:
    for i, (sig, t, res) in enumerate(off_critical_zeros):
        d = diagnose_zero(sig, t, f"off-crit #{i+1}")
        off_crit_diag.append(d)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. 비교 테이블
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log()
log("━━━ 5. 비교 테이블 ━━━")
log()

def stats(values):
    """평균 ± 표준편차"""
    if len(values) == 0:
        return "N/A", "N/A"
    arr = np.array(values)
    return f"{np.mean(arr):.4f}", f"{np.std(arr):.4f}"

if on_crit_diag:
    on_kappas = [d['kappa'] for d in on_crit_diag if d['kappa'] < 1e9]
    on_monos = [d['mono_pi'] for d in on_crit_diag]
    on_peaks = [d['peak_sigma'] for d in on_crit_diag]
    on_fwhms = [d['fwhm'] for d in on_crit_diag]

    log(f"On-critical (N={len(on_crit_diag)}):")
    log(f"  κ 평균: {stats(on_kappas)[0]} ± {stats(on_kappas)[1]}")
    log(f"  mono/π: {stats(on_monos)[0]} ± {stats(on_monos)[1]}")
    log(f"  peak_σ: {stats(on_peaks)[0]} ± {stats(on_peaks)[1]}")
    log(f"  FWHM:   {stats(on_fwhms)[0]} ± {stats(on_fwhms)[1]}")
else:
    log("On-critical: 데이터 없음")

log()

if off_crit_diag:
    off_kappas = [d['kappa'] for d in off_crit_diag if d['kappa'] < 1e9]
    off_monos = [d['mono_pi'] for d in off_crit_diag]
    off_peaks = [d['peak_sigma'] for d in off_crit_diag]
    off_fwhms = [d['fwhm'] for d in off_crit_diag]

    log(f"Off-critical (N={len(off_crit_diag)}):")
    log(f"  κ 평균: {stats(off_kappas)[0]} ± {stats(off_kappas)[1]}")
    log(f"  mono/π: {stats(off_monos)[0]} ± {stats(off_monos)[1]}")
    log(f"  peak_σ: {stats(off_peaks)[0]} ± {stats(off_peaks)[1]}")
    log(f"  FWHM:   {stats(off_fwhms)[0]} ± {stats(off_fwhms)[1]}")

    log()
    log("--- 핵심 비교 ---")
    log(f"{'메트릭':<12} {'On-critical':<25} {'Off-critical':<25} {'차이'}")
    log("-" * 80)

    if on_kappas and off_kappas:
        log(f"{'κ':<12} {np.mean(on_kappas):>10.2f} ± {np.std(on_kappas):>8.2f}   {np.mean(off_kappas):>10.2f} ± {np.std(off_kappas):>8.2f}   {'—'}")

    if on_monos and off_monos:
        on_m = np.mean(on_monos)
        off_m = np.mean(off_monos)
        log(f"{'mono/π':<12} {on_m:>10.4f} ± {np.std(on_monos):>8.4f}   {off_m:>10.4f} ± {np.std(off_monos):>8.4f}   {'동일' if abs(on_m - off_m) < 0.1 else '다름'}")

    if on_peaks and off_peaks:
        on_p = np.mean(on_peaks)
        off_p = np.mean(off_peaks)
        log(f"{'peak_σ':<12} {on_p:>10.4f} ± {np.std(on_peaks):>8.4f}   {off_p:>10.4f} ± {np.std(off_peaks):>8.4f}   {'★ 분리!' if abs(on_p - off_p) > 0.03 else '겹침'}")

    if on_fwhms and off_fwhms:
        on_f = np.mean(on_fwhms)
        off_f = np.mean(off_fwhms)
        log(f"{'FWHM':<12} {on_f:>10.4f} ± {np.std(on_fwhms):>8.4f}   {off_f:>10.4f} ± {np.std(off_fwhms):>8.4f}   {'★ 넓어짐!' if off_f > on_f * 1.3 else '유사'}")
else:
    log("Off-critical: 데이터 없음 — 비교 불가")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. 개별 영점 상세 테이블
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log()
log("━━━ 6. 개별 영점 진단 상세 ━━━")
log()
log(f"{'#':<5} {'유형':<12} {'σ₀':<10} {'t₀':<12} {'κ':<12} {'mono/π':<10} {'peak_σ':<10} {'FWHM':<10}")
log("-" * 85)

all_diag = [(d, 'ON') for d in on_crit_diag] + [(d, 'OFF') for d in off_crit_diag]
for i, (d, typ) in enumerate(all_diag):
    kappa_str = f"{d['kappa']:.2f}" if d['kappa'] < 1e9 else "∞"
    log(f"{i+1:<5} {typ:<12} {d['sigma_0']:<10.6f} {d['t_0']:<12.6f} {kappa_str:<12} {d['mono_pi']:<10.4f} {d['peak_sigma']:<10.4f} {d['fwhm']:<10.4f}")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7. 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log()
log("━━━ 7. 판정 ━━━")
log()

n_on = len(on_critical_zeros)
n_off = len(off_critical_zeros)

if n_off > 0 and off_crit_diag:
    # 핵심 판별: peak_σ와 FWHM 분리
    peak_sep = abs(np.mean(off_peaks) - np.mean(on_peaks)) if on_peaks else 0
    fwhm_ratio = np.mean(off_fwhms) / np.mean(on_fwhms) if on_fwhms and np.mean(on_fwhms) > 0 else 0

    if peak_sep > 0.03 and fwhm_ratio > 1.3:
        verdict = "강한 양성"
        reason = f"peak_σ 분리 = {peak_sep:.4f} (>0.03), FWHM 비율 = {fwhm_ratio:.2f} (>1.3)"
    elif peak_sep > 0.03 or fwhm_ratio > 1.3:
        verdict = "양성"
        reason = f"peak_σ 분리 = {peak_sep:.4f}, FWHM 비율 = {fwhm_ratio:.2f} — 한 지표 이상 분리"
    else:
        verdict = "음성"
        reason = f"peak_σ 분리 = {peak_sep:.4f} (<0.03), FWHM 비율 = {fwhm_ratio:.2f} (<1.3) — 프레임워크가 구분 못 함"
elif n_off == 0:
    verdict = "미결"
    reason = f"Off-critical 영점 미발견 (on-critical {n_on}개). 탐색 범위 확대 또는 구현 재검증 필요."
else:
    verdict = "미결"
    reason = "진단 데이터 불충분"

log(f"On-critical 영점: {n_on}개")
log(f"Off-critical 영점: {n_off}개")
log(f"판정: **{verdict}**")
log(f"근거: {reason}")

log()
log(f"DH 함수 주의: 오일러 곱 없음 (L-함수 아님). 함수방정식만 만족.")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 8. 설정 기록
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

elapsed = time.time() - START_TIME
log("━━━ 설정 ━━━")
log(f"mpmath.dps = {mpmath.mp.dps}")
log(f"χ mod 5: {CHI_MOD5}")
log(f"κ (DH 상수) = {float(kappa_DH):.10f}")
log(f"δ 오프셋 = {DELTA}")
log(f"임계선 스캔: t∈[1,60], n=3000")
log(f"격자 스캔: σ∈[0.52,1.0]×t∈[1,60], {n_sigma}×{n_t}={n_sigma*n_t}점")
log(f"σ<0.5 추가: σ∈[0.0,0.48]×t∈[1,60], 20×200=4000점")
log(f"모노드로미: 폐곡선 적분 (radius=0.5, 64단계)")
log(f"σ-국소화 프로파일: 61점")
log(f"총 소요 시간: {elapsed:.1f}초 ({elapsed/60:.1f}분)")
log(f"종료 시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")

outf.close()
print(f"\n결과 저장: {RESULT_FILE}")
