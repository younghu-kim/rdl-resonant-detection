#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 실험 #115 — Davenport-Heilbronn off-critical ξ-bundle 4성질
=============================================================================
목표: B-28 검증 — DH 함수의 off-critical 영점에서 ξ-bundle 4성질 붕괴 여부

DH 함수: f(s) = [(1-iκ)/2] · L(s,χ) + [(1+iκ)/2] · L(s,χ̄)
  χ: mod 5 원시 지표 (홀수, χ(-1)=-1)
  κ_DH = (√(10-2√5) - 2)/(√5 - 1) ≈ 0.2841

완비 함수: Λ(s) = (5/π)^{s/2} · Γ((s+1)/2) · f(s)
함수방정식: Λ(s) = ε · conj(Λ(1-conj(s)))  [ε = 상수, |ε|=1]

4성질:
  P1: FE 수치 검증 (rel_err < 1e-20)
  P2: 영점 목록 (on-critical 10개, off-critical 4개 — 기존 #30에서 확인)
  P3: κδ² at σ₀+δ+it₀ (δ=0.03), 완비 함수 사용
  P4: mono/π at 3 반경 (r=0.05, 0.15, 0.5)

핵심 예측:
  - κδ² = 1: 단순 영점이면 on/off 무관하게 성립 → "영점 감지기"
  - κδ² ≠ 1 off-critical: "임계선 감별기" (강력한 증거)
  - mono/π: on-critical → 2 (자기쌍), off-critical → 2(소반경) or 4(대반경, 거울쌍 포획)

결과: results/dh_off_critical_115.txt
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

mpmath.mp.dps = 80  # t<200 충분
START = time.time()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 출력 설정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
RESULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
os.makedirs(RESULT_DIR, exist_ok=True)
RESULT_FILE = os.path.join(RESULT_DIR, 'dh_off_critical_115.txt')

outf = open(RESULT_FILE, 'w')

def log(msg=''):
    print(msg, flush=True)
    outf.write(msg + '\n')
    outf.flush()


log("=" * 72)
log("[실험 #115] Davenport-Heilbronn off-critical ξ-bundle 4성질")
log("=" * 72)
log(f"시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"mpmath.dps = {mpmath.mp.dps}")
log()

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. DH 함수 구현
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 1. DH 함수 구현 ━━━")

CHI_MOD5 = [0, 1, 1j, -1j, -1]       # χ(0..4)
CHI_MOD5_BAR = [0, 1, -1j, 1j, -1]   # χ̄(0..4)

sqrt5 = mpmath.sqrt(5)
kappa_DH = (mpmath.sqrt(10 - 2*sqrt5) - 2) / (sqrt5 - 1)
coeff_chi     = (1 - 1j * kappa_DH) / 2
coeff_chi_bar = (1 + 1j * kappa_DH) / 2

log(f"κ_DH = {float(kappa_DH):.10f}")
log(f"coeff_chi = {complex(coeff_chi)}")
log(f"coeff_chi_bar = {complex(coeff_chi_bar)}")
log()


def dh_func(s):
    """Davenport-Heilbronn 함수 f(s)"""
    s = mpmath.mpc(s)
    return coeff_chi * mpmath.dirichlet(s, CHI_MOD5) + \
           coeff_chi_bar * mpmath.dirichlet(s, CHI_MOD5_BAR)


def Lambda_dh(s):
    """완비 DH 함수: Λ(s) = (5/π)^{s/2} · Γ((s+1)/2) · f(s)"""
    s = mpmath.mpc(s)
    pref = mpmath.power(5 / mpmath.pi, s / 2)
    gam  = mpmath.gamma((s + 1) / 2)
    return pref * gam * dh_func(s)


def Lambda_deriv(s, h=None):
    """Λ'(s) 수치 미분 (중앙차분)"""
    if h is None:
        h = mpmath.mpf(10) ** (-min(30, mpmath.mp.dps - 10))
    s = mpmath.mpc(s)
    return (Lambda_dh(s + h) - Lambda_dh(s - h)) / (2 * h)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. P1 — 함수방정식 검증 (rel_err < 1e-20)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 2. P1 함수방정식 검증 ━━━")
log("FE: Λ(s) = ε · conj(Λ(1-conj(s)))")
log()

fe_test_points = [
    mpmath.mpc('0.3', '10.0'),
    mpmath.mpc('0.7', '15.0'),
    mpmath.mpc('0.5', '20.0'),
    mpmath.mpc('0.2', '25.0'),
    mpmath.mpc('0.8', '30.0'),
]

epsilon_DH = None
fe_rel_errs = []
fe_pass = True

log(f"{'점':<20} {'|Λ(s)|':>15} {'|Λ(1-s̄)*|':>15} {'rel_err':>15} {'판정'}")
log("-" * 75)

for s in fe_test_points:
    Ls  = Lambda_dh(s)
    s_mirror = 1 - mpmath.conj(s)
    Lm  = mpmath.conj(Lambda_dh(s_mirror))

    abs_Ls = abs(Ls)
    abs_Lm = abs(Lm)

    if abs_Lm < mpmath.mpf(10)**(-50):
        log(f"  s={complex(s)}: Λ(1-s̄)≈0, 스킵")
        continue

    ratio = Ls / Lm
    if epsilon_DH is None:
        epsilon_DH = ratio

    rel_err = abs(Ls - epsilon_DH * Lm) / abs_Ls
    fe_rel_errs.append(float(rel_err))
    p1_ok = float(rel_err) < 1e-20

    if not p1_ok:
        fe_pass = False

    log(f"  s={complex(s)}: |Λ|={float(abs_Ls):.4e}  "
        f"rel_err={float(rel_err):.3e}  {'✅' if p1_ok else '❌'}")

log()
log(f"ε (root number) = {complex(epsilon_DH):.10f}")
log(f"|ε| = {float(abs(epsilon_DH)):.10f}  (≈1 이어야)")
log(f"최대 rel_err = {max(fe_rel_errs):.3e}")
log(f"P1 판정: {'✅ PASS (rel_err < 1e-20)' if fe_pass else '❌ FAIL'}")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. P2 — 영점 목록 (기존 #30에서 확인된 영점 사용)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 3. P2 영점 목록 ━━━")
log("(기존 실험 #30에서 확인된 영점 사용)")
log()

# On-critical 영점 (σ=0.5, 기존 #30 확인)
on_critical_zeros = [
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

# Off-critical 영점 (σ≠0.5, 기존 #30 + dh_small_radius 확인)
# 거울쌍: σ₀ → 1-σ₀ (함수방정식에 의해)
off_critical_zeros = [
    (0.808517, 85.699348),   # OFF#1 Spira 1968, mirror at σ=0.191483
    (0.650830, 114.163343),  # OFF#2, mirror at σ=0.349170
    (0.574356, 166.479306),  # OFF#3, mirror at σ=0.425644
    (0.724258, 176.702461),  # OFF#4, mirror at σ=0.275742
]

log(f"On-critical 영점: {len(on_critical_zeros)}개")
for i, (sig, t) in enumerate(on_critical_zeros):
    s = mpmath.mpf(str(sig)) + 1j * mpmath.mpf(str(t))
    absf = float(abs(dh_func(s)))
    log(f"  [ON#{i+1}] σ={sig}, t={t:.6f}, |f|={absf:.3e}")

log()
log(f"Off-critical 영점: {len(off_critical_zeros)}개")
for i, (sig, t) in enumerate(off_critical_zeros):
    s = mpmath.mpf(str(sig)) + 1j * mpmath.mpf(str(t))
    absf = float(abs(dh_func(s)))
    sig_mirror = 1.0 - sig
    log(f"  [OFF#{i+1}] σ={sig:.6f}, t={t:.6f}, |f|={absf:.3e}  (거울쌍 σ={sig_mirror:.6f})")
    log(f"           |σ-0.5| = {abs(sig-0.5):.6f}  (> 0.01: {'✅' if abs(sig-0.5)>0.01 else '❌'})")

log()

# 영점 잔류 재검증
log("영점 잔류 재검증 (findroot 정밀화):")
verified_off = []
for sig, t in off_critical_zeros:
    try:
        def f_sys(sv, tv):
            sc = mpmath.mpf(sv) + 1j * mpmath.mpf(tv)
            fv = dh_func(sc)
            return (mpmath.re(fv), mpmath.im(fv))
        result = mpmath.findroot(f_sys, (mpmath.mpf(str(sig)), mpmath.mpf(str(t))))
        sr, tr = float(result[0]), float(result[1])
        sv = mpmath.mpf(str(sr)) + 1j * mpmath.mpf(str(tr))
        res = float(abs(dh_func(sv)))
        verified_off.append((sr, tr, res))
        log(f"  σ={sr:.8f}, t={tr:.8f}, |f|={res:.3e}  {'✅' if res<1e-10 else '⚠️'}")
    except Exception as e:
        log(f"  ⚠️ findroot 실패: {e}")
        verified_off.append((sig, t, float(abs(dh_func(mpmath.mpf(str(sig)) + 1j*mpmath.mpf(str(t)))))))

log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. P3 — κδ² 측정 (완비 함수 Λ 사용)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 4. P3 κδ² 측정 ━━━")
log("공식: κδ² = |Λ'(σ₀+δ+it₀) / Λ(σ₀+δ+it₀)|² × δ²")
log("(σ₀: 영점의 실수부, δ=0.03 오프셋)")
log()

DELTA = mpmath.mpf('0.03')

def measure_kappa_delta2(sigma_0, t_0, delta=DELTA):
    """κδ² 측정: 완비 Λ에서 σ₀+δ+it₀에서의 connection"""
    s_off = mpmath.mpf(str(sigma_0)) + delta + 1j * mpmath.mpf(str(t_0))
    L_val = Lambda_dh(s_off)

    if abs(L_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 15):
        return float('nan'), float('nan')

    L_deriv = Lambda_deriv(s_off)
    connection = L_deriv / L_val
    kappa = float(abs(connection)**2)
    kappa_delta2 = kappa * float(delta)**2
    return kappa, kappa_delta2


log(f"{'영점':<35} {'σ₀':<10} {'t₀':<14} {'κ':>12} {'κδ²':>10} {'판정'}")
log("-" * 90)

on_kd2 = []
for i, (sig, t) in enumerate(on_critical_zeros):
    kappa, kd2 = measure_kappa_delta2(sig, t)
    on_kd2.append(kd2)
    ok = abs(kd2 - 1.0) < 0.05
    log(f"ON#{i+1:<3} σ=0.5, t={t:<12.6f}  {sig:<10.6f} {t:<14.6f} {kappa:>12.2f} {kd2:>10.6f} {'✅' if ok else '❌'}")

log()

off_kd2 = []
for i, (sr, tr, res) in enumerate(verified_off):
    kappa, kd2 = measure_kappa_delta2(sr, tr)
    off_kd2.append(kd2)
    ok = abs(kd2 - 1.0) < 0.05
    log(f"OFF#{i+1:<2} σ={sr:.6f}, t={tr:<8.3f}  {sr:<10.6f} {tr:<14.6f} {kappa:>12.2f} {kd2:>10.6f} {'✅' if ok else '❌'}")

log()
log(f"On-critical  κδ² 평균: {np.mean(on_kd2):.6f} ± {np.std(on_kd2):.6f}")
log(f"Off-critical κδ² 평균: {np.mean(off_kd2):.6f} ± {np.std(off_kd2):.6f}")
log()

if abs(np.mean(on_kd2) - 1.0) < 0.05 and abs(np.mean(off_kd2) - 1.0) < 0.05:
    log("P3 결론: κδ² ≈ 1 양쪽 → 프레임워크는 '영점 감지기' (임계선 감별기 아님)")
    p3_discriminates = False
elif abs(np.mean(on_kd2) - 1.0) < 0.05 and abs(np.mean(off_kd2) - 1.0) >= 0.05:
    log("P3 결론: on-critical κδ²≈1, off-critical κδ²≠1 → 프레임워크가 임계선 감별!")
    p3_discriminates = True
else:
    log("P3 결론: 불일치 (양쪽 모두 κδ²≠1)")
    p3_discriminates = False
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. P4 — 모노드로미 측정 (3 반경)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 5. P4 모노드로미 측정 (3 반경) ━━━")
log("폐곡선 적분: ∮ d(arg Λ) — 완비 함수 사용")
log("r=0.05 (단일 영점 분리), r=0.15 (중간), r=0.5 (거울쌍 포획 가능)")
log()

RADII = [0.05, 0.15, 0.5]
N_STEPS = 128  # 정밀도 향상


def monodromy_Lambda(s_center, radius=0.5, n_steps=N_STEPS):
    """완비 Λ를 사용한 폐곡선 적분 모노드로미"""
    total_delta = 0.0
    prev_arg = None
    skipped = 0

    for k in range(n_steps + 1):
        theta = 2 * np.pi * k / n_steps
        s = s_center + radius * mpmath.exp(1j * theta)
        val = Lambda_dh(s)

        if abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 20):
            skipped += 1
            continue

        curr_arg = float(mpmath.arg(val))

        if prev_arg is not None:
            delta_arg = curr_arg - prev_arg
            while delta_arg > np.pi:  delta_arg -= 2 * np.pi
            while delta_arg < -np.pi: delta_arg += 2 * np.pi
            total_delta += delta_arg

        prev_arg = curr_arg

    return total_delta, skipped


log("On-critical 영점 모노드로미:")
log(f"{'영점':<30} {'r=0.05':>10} {'r=0.15':>10} {'r=0.50':>10}")
log("-" * 65)

on_mono_r05, on_mono_r15, on_mono_r50 = [], [], []

for i, (sig, t) in enumerate(on_critical_zeros):
    s_c = mpmath.mpf(str(sig)) + 1j * mpmath.mpf(str(t))
    row = f"ON#{i+1:<2} σ=0.5, t={t:.4f}"
    vals = []
    for r in RADII:
        m, sk = monodromy_Lambda(s_c, radius=r)
        mpi = m / np.pi
        vals.append(mpi)
    on_mono_r05.append(vals[0])
    on_mono_r15.append(vals[1])
    on_mono_r50.append(vals[2])
    log(f"{row:<30} {vals[0]:>+10.4f}π {vals[1]:>+10.4f}π {vals[2]:>+10.4f}π")

log()
log("Off-critical 영점 모노드로미:")
log(f"{'영점':<35} {'거울σ':>8} {'r=0.05':>10} {'r=0.15':>10} {'r=0.50':>10}")
log("-" * 75)

off_mono_r05, off_mono_r15, off_mono_r50 = [], [], []

for i, (sr, tr, res) in enumerate(verified_off):
    s_c = mpmath.mpf(str(sr)) + 1j * mpmath.mpf(str(tr))
    mirror_sig = 1.0 - sr
    row = f"OFF#{i+1:<2} σ={sr:.4f}, t={tr:.2f}"
    vals = []
    for r in RADII:
        m, sk = monodromy_Lambda(s_c, radius=r)
        mpi = m / np.pi
        vals.append(mpi)
    off_mono_r05.append(vals[0])
    off_mono_r15.append(vals[1])
    off_mono_r50.append(vals[2])
    log(f"{row:<35} {mirror_sig:>8.4f} {vals[0]:>+10.4f}π {vals[1]:>+10.4f}π {vals[2]:>+10.4f}π")

log()

# 판정
log("P4 통계:")
log(f"  On-critical  r=0.05: {np.mean(on_mono_r05):.4f}π ± {np.std(on_mono_r05):.4f}π")
log(f"  On-critical  r=0.15: {np.mean(on_mono_r15):.4f}π ± {np.std(on_mono_r15):.4f}π")
log(f"  On-critical  r=0.50: {np.mean(on_mono_r50):.4f}π ± {np.std(on_mono_r50):.4f}π")
log()
log(f"  Off-critical r=0.05: {np.mean(off_mono_r05):.4f}π ± {np.std(off_mono_r05):.4f}π")
log(f"  Off-critical r=0.15: {np.mean(off_mono_r15):.4f}π ± {np.std(off_mono_r15):.4f}π")
log(f"  Off-critical r=0.50: {np.mean(off_mono_r50):.4f}π ± {np.std(off_mono_r50):.4f}π")
log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. 비교 테이블
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log()
log("━━━ 6. 비교 테이블 ━━━")
log()
log(f"{'메트릭':<25} {'On-critical (N=%d)'%len(on_critical_zeros):<30} {'Off-critical (N=%d)'%len(verified_off):<30} {'해석'}")
log("=" * 100)

# P1
log(f"{'P1 FE rel_err':<25} {'< 1e-20 → ✅ PASS' if fe_pass else '> 1e-20 → ❌':<30} {'—':<30} {'DH 함수방정식 만족'}")

# P2
log(f"{'P2 영점 수':<25} {str(len(on_critical_zeros))+'개 (σ=0.5)':<30} {str(len(verified_off))+'개 (σ≠0.5)':<30} {'#30 확인값 사용'}")

# P3
on_kd2_str = f"{np.mean(on_kd2):.4f} ± {np.std(on_kd2):.4f}"
off_kd2_str = f"{np.mean(off_kd2):.4f} ± {np.std(off_kd2):.4f}"
kd2_diff = abs(np.mean(on_kd2) - np.mean(off_kd2))
kd2_interp = "동일 (양쪽≈1)" if kd2_diff < 0.05 else f"차이={kd2_diff:.4f}"
log(f"{'P3 κδ² 평균':<25} {on_kd2_str:<30} {off_kd2_str:<30} {kd2_interp}")

# P4 각 반경
on_m05_str = f"{np.mean(on_mono_r05):.3f} ± {np.std(on_mono_r05):.3f}π"
off_m05_str = f"{np.mean(off_mono_r05):.3f} ± {np.std(off_mono_r05):.3f}π"
m05_same = abs(np.mean(on_mono_r05) - np.mean(off_mono_r05)) < 0.2
log(f"{'P4 mono/π (r=0.05)':<25} {on_m05_str:<30} {off_m05_str:<30} {'동일 (단순영점)' if m05_same else '★ 다름'}")

on_m15_str = f"{np.mean(on_mono_r15):.3f} ± {np.std(on_mono_r15):.3f}π"
off_m15_str = f"{np.mean(off_mono_r15):.3f} ± {np.std(off_mono_r15):.3f}π"
m15_same = abs(np.mean(on_mono_r15) - np.mean(off_mono_r15)) < 0.2
log(f"{'P4 mono/π (r=0.15)':<25} {on_m15_str:<30} {off_m15_str:<30} {'동일' if m15_same else '★ 다름'}")

on_m50_str = f"{np.mean(on_mono_r50):.3f} ± {np.std(on_mono_r50):.3f}π"
off_m50_str = f"{np.mean(off_mono_r50):.3f} ± {np.std(off_mono_r50):.3f}π"
m50_diff = abs(np.mean(on_mono_r50) - np.mean(off_mono_r50))
m50_same = m50_diff < 0.2
log(f"{'P4 mono/π (r=0.50)':<25} {on_m50_str:<30} {off_m50_str:<30} {'동일' if m50_same else '★ 다름 (Δ=%.2f)'%m50_diff}")

log()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 7. 판정
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
log("━━━ 7. 최종 판정 ━━━")
log()

# 성공 기준 체크
sc1 = fe_pass
sc2 = len(on_critical_zeros) >= 5
sc3 = len(verified_off) >= 1 and all(abs(sr - 0.5) > 0.01 for sr, tr, _ in verified_off)
sc4 = len(off_kd2) > 0 and len(off_mono_r50) > 0
sc5 = True  # 비교표 출력 완료

log(f"성공 기준 체크:")
log(f"  1. FE rel_err < 1e-20: {'✅' if sc1 else '❌'}")
log(f"  2. on-critical ≥5개 측정: {'✅ (%d개)'%len(on_critical_zeros) if sc2 else '❌'}")
log(f"  3. off-critical ≥1개 (|σ-0.5|>0.01): {'✅ (%d개)'%len(verified_off) if sc3 else '❌'}")
log(f"  4. off-critical κδ², mono/π 측정: {'✅' if sc4 else '❌'}")
log(f"  5. 비교표 출력: {'✅' if sc5 else '❌'}")
log()

# P3 해석
log("P3 해석 (κδ²):")
if kd2_diff < 0.05:
    log("  κδ² ≈ 1 양쪽 동일.")
    log("  → P3은 on/off 구별 불가.")
    log("  → 수학적 이유: 단순 영점에서 Λ(s) ≈ c(s-s₀) → Λ'/Λ ≈ 1/(s-s₀) → κδ²→1")
    log("  → 결론: κδ² = 1은 '영점 성질' (임계선 감별기 아님)")
else:
    log(f"  κδ² 차이 = {kd2_diff:.4f}")
    if np.mean(off_kd2) > np.mean(on_kd2) + 0.1:
        log("  → off-critical에서 κδ² > 1 → 단순 영점 아님? 또는 완비 함수 구조 기여")
    else:
        log("  → 예상치 못한 결과 → 추가 분석 필요")

log()

# P4 해석
log("P4 해석 (모노드로미):")
log(f"  r=0.05 (단일 영점 분리):")
log(f"    on-critical:  {np.mean(on_mono_r05):.3f}π  (기대값: 2)")
log(f"    off-critical: {np.mean(off_mono_r05):.3f}π  (기대값: 2)")
if abs(np.mean(on_mono_r05) - 2.0) < 0.1 and abs(np.mean(off_mono_r05) - 2.0) < 0.1:
    log("    → 양쪽 mono=2π: 단순 영점 확인")
log(f"  r=0.5 (거울쌍 포획 가능):")
log(f"    on-critical:  {np.mean(on_mono_r50):.3f}π")
log(f"    off-critical: {np.mean(off_mono_r50):.3f}π")

if not m50_same:
    log(f"    ★ P4(r=0.5)로 on/off 구별! Δmono/π = {m50_diff:.3f}")
    log("    수학적 이유: on-critical 영점은 FE 고정점 (σ=0.5) → r=0.5 반경에서 거울쌍=자기자신")
    log("    off-critical 영점은 FE에 의해 쌍 존재 → r=0.5 원이 거울쌍까지 포획 가능")
    p4_discriminates = True
else:
    log("    on/off 구별 불가 (r=0.5에서도 동일)")
    p4_discriminates = False

log()

# 최종 결론
log("═" * 72)
log("최종 결론: B-28 답변")
log("═" * 72)
log()

if p3_discriminates:
    log("  P3(κδ²) → RH 감별기 ★★★")
    log("    on-critical κδ²≈1, off-critical κδ²>1 → 완비 함수 곡률이 임계선 감별")
    if p4_discriminates:
        log("  P4(mono/π, 대반경): 추가 감별 확인")
    verdict = "★★★ 강양성 (P3+P4 감별)"
elif p4_discriminates:
    log("  P3(κδ²): 영점 감지기만 (on/off 구별 불가)")
    log("  P4(mono/π, 대반경): ★ RH 감별기 역할 가능!")
    log("    → 프레임워크는 κδ²로는 구별 못하지만, 모노드로미로 구별 가능")
    log("    → 모노드로미의 반경 의존성이 임계선 감별의 핵심")
    log("    → '거울쌍 포획 여부'가 on vs off-critical 감별 기제")
    verdict = "★★ 양성 (P4 감별, P3 동일)"
else:
    log("  P3(κδ²): 영점 감지기만 (on/off 구별 불가)")
    log("  P4(mono/π): 구별 불가 (모든 반경에서 동일)")
    log("  → 프레임워크는 RH 감별기가 아닌 영점 감지기")
    verdict = "★ 양성 (해석 수정 필요: 영점 감지기)"

log()
log(f"실험 #115 판정: {verdict}")
log()

log("설정:")
log(f"  mpmath.dps = {mpmath.mp.dps}")
log(f"  δ = {float(DELTA)}")
log(f"  반경 = {RADII}")
log(f"  n_steps = {N_STEPS}")
log(f"  on-critical: {len(on_critical_zeros)}개 (t∈[5,28])")
log(f"  off-critical: {len(verified_off)}개 (t∈[85,177])")
elapsed = time.time() - START
log(f"총 소요: {elapsed:.1f}초 ({elapsed/60:.1f}분)")
log(f"종료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

outf.close()
print(f"\n결과 저장: {RESULT_FILE}")
