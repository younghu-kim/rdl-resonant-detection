#!/usr/bin/env python3
"""
[Project RDL] #66 — 감마 인자 역보정(deconvolution) 실험
=============================================================================
가설: σ-유일성 FAIL은 감마 인자가 영점 신호를 "매끄럽게 퍼뜨린" 결과.
      감마 인자의 σ-의존 위상을 제거하면 σ=0.5의 특이성이 복원되는가?

방법:
  원래: Λ(s) = Q^s · γ(s) · L(s)  →  Re(Λ) 부호변환 카운트
  보정: Λ_corr(s) = Λ(s) · exp(-i·[arg γ(s) - arg γ(1/2+it)])
        = σ=0.5 기준으로 감마 인자의 위상 회전을 되돌린 것.

  만약 감마 인자의 위상이 "렌즈"라면, 역보정 후 σ=0.5가 극대로 복원되어야 함.

테스트: GL(3) sym²(11a1) + GL(2) 11a1 (비교용)

결과 파일: results/gl3_gamma_deconv_66.txt
"""
import sys, os, time
import numpy as np
from scipy.special import loggamma

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "results", "gl3_gamma_deconv_66.txt")
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))
def flush():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━ 공통 파라미터 ━━━━━━
C_SHIFT = 4.0
N_GH    = 50
GH_X, GH_W = np.polynomial.hermite.hermgauss(N_GH)
T_MIN, T_MAX = 2.0, 50.0
SIGMA_VALS = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # 세밀 스윕

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# GL(3) sym²(11a1): N=121, μ=(0,1,1), ε=+1
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def gl3_log_gamma(s):
    """log γ(s) for GL(3) sym²: Γ_ℝ(s)·Γ_ℝ(s+1)²"""
    return (-(3*s+2)/2 * np.log(np.pi) +
            loggamma(s/2) + 2*loggamma((s+1)/2))

def gl3_gamma_phase(sigma, t):
    """arg γ(σ+it) for GL(3)"""
    s = complex(sigma, t)
    lg = gl3_log_gamma(s)
    return lg.imag

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# GL(2) 11a1: N=11, Γ_ℂ(s) = Γ_ℝ(s)·Γ_ℝ(s+1), ε=+1
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def gl2_log_gamma(s):
    """log γ(s) for GL(2) 11a1: Γ_ℝ(s)·Γ_ℝ(s+1) = Γ_ℂ(s)"""
    return (-(2*s+1)/2 * np.log(np.pi) +
            loggamma(s/2) + loggamma((s+1)/2))

def gl2_gamma_phase(sigma, t):
    """arg γ(σ+it) for GL(2)"""
    s = complex(sigma, t)
    lg = gl2_log_gamma(s)
    return lg.imag

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# GL(1) ζ(s): N=1, Γ_ℝ(s), ε 해당없음 (참고용)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def gl1_log_gamma(s):
    """log γ(s) for GL(1) ζ: Γ_ℝ(s)"""
    return (-s/2 * np.log(np.pi) + loggamma(s/2))

def gl1_gamma_phase(sigma, t):
    s = complex(sigma, t)
    lg = gl1_log_gamma(s)
    return lg.imag

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 11a1 계수 (GL(2) + GL(3) sym² 공유)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_11a1_an(limit):
    known = {2:-2,3:-1,5:1,7:-2,11:1,13:4,17:-2,19:0,23:-1,
             29:0,31:7,37:3,41:-8,43:-6,47:8,53:-6,59:5,61:12}
    sieve = [True]*(limit+1); sieve[0]=sieve[1]=False
    for i in range(2, int(limit**0.5)+1):
        if sieve[i]:
            for j in range(i*i, limit+1, i): sieve[j] = False
    primes = [i for i in range(2, limit+1) if sieve[i]]
    ap = {}
    for p in primes:
        if p in known: ap[p] = known[p]; continue
        cnt = 0
        for x in range(p):
            d = (4*x*x*x - 4*x*x - 40*x - 79) % p
            if d == 0: cnt += 1
            elif pow(d, (p-1)//2, p) == 1: cnt += 2
        ap[p] = p - cnt
    apk = {}
    for p in primes:
        apk[(p,0)] = 1; apk[(p,1)] = ap[p]
        pk = p; k = 1
        while pk*p <= limit:
            pk *= p; k += 1
            if p == 11: apk[(p,k)] = ap[p]**k
            else: apk[(p,k)] = ap[p]*apk[(p,k-1)] - p*apk[(p,k-2)]
    an = [0]*(limit+1); an[1] = 1
    for n in range(2, limit+1):
        temp = n; result = 1
        for p in primes:
            if p*p > temp: break
            if temp % p == 0:
                k = 0
                while temp % p == 0: k += 1; temp //= p
                result *= apk.get((p,k), 0)
        if temp > 1: result *= ap.get(temp, 0)
        an[n] = result
    return an, ap, primes

def compute_sym2_cn(an, primes, limit):
    """sym²(11a1) Dirichlet 계수"""
    cpk = {}
    for p in primes:
        if p > limit: break
        at2 = an[p]**2 / p
        if p == 11:
            for k in range(20):
                cpk[(p,k)] = 11.0**(-k)
                if p**(k+1) > limit: break
        else:
            c1 = at2 - 1.0
            cpk[(p,0)] = 1.0; cpk[(p,1)] = c1
            pk = p; k = 1
            while pk*p <= limit:
                pk *= p; k += 1
                bkm1 = cpk[(p,k-1)]; bkm2 = cpk[(p,k-2)]
                bkm3 = cpk.get((p,k-3), 0.0)
                cpk[(p,k)] = c1*bkm1 - c1*bkm2 + bkm3 if k > 2 else c1*c1 - c1
    cn = np.zeros(limit+1)
    cn[1] = 1.0
    for n in range(2, limit+1):
        temp = n; result = 1.0
        for p in primes:
            if p*p > temp: break
            if temp % p == 0:
                k = 0
                while temp % p == 0: k += 1; temp //= p
                if (p,k) in cpk: result *= cpk[(p,k)]
                else: result = 0.0; break
        if result != 0.0 and temp > 1:
            if (temp,1) in cpk: result *= cpk[(temp,1)]
            elif temp in [pp for pp in primes if pp <= limit]:
                result *= (an[temp]**2/temp - 1.0)
            else: result = 0.0
        cn[n] = result
    return cn

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Λ 계산 — Mellin-Barnes AFE
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def Lambda_direct(s, log_gamma_fn, Q, cn_nz, log_n_nz, cn_vals_nz):
    s = complex(s)
    Qs = Q ** s
    gamma_s = np.exp(log_gamma_fn(s))
    powers = np.exp(-s * log_n_nz)
    L_s = np.dot(cn_vals_nz, powers)
    return Qs * gamma_s * L_s

def Lambda_AFE(s, log_gamma_fn, Q, eps, cn_nz, log_n_nz, cn_vals_nz,
               c=C_SHIFT):
    s = complex(s)
    s1 = 1.0 - s
    ec2 = np.exp(c**2)
    total = 0.0 + 0.0j
    for k in range(N_GH):
        y = GH_X[k]; w = GH_W[k]
        wk = c + 1j*y
        phase = np.exp(2j * c * y) / wk
        u_fwd = s + wk; u_bwd = s1 + wk
        lam_f = Lambda_direct(u_fwd, log_gamma_fn, Q, cn_nz, log_n_nz, cn_vals_nz)
        lam_b = Lambda_direct(u_bwd, log_gamma_fn, Q, cn_nz, log_n_nz, cn_vals_nz)
        total += w * (lam_f + eps * lam_b) * phase
    return ec2 / (2*np.pi) * total

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# σ-유일성 측정 (원본 + 보정)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def measure_sign_changes(sigma, t_arr, log_gamma_fn, Q, eps,
                         cn_nz, log_n_nz, cn_vals_nz, gamma_phase_fn,
                         correct=False):
    """
    correct=False: 원본 Re(Λ(σ+it)) 부호변환
    correct=True:  감마 위상 역보정 후 Re(Λ_corr) 부호변환
      Λ_corr(σ+it) = Λ(σ+it) · exp(-i·[arg γ(σ+it) - arg γ(0.5+it)])
    """
    vals = []
    for t in t_arr:
        lam = Lambda_AFE(complex(sigma, t), log_gamma_fn, Q, eps,
                         cn_nz, log_n_nz, cn_vals_nz)
        if correct:
            # σ-의존 감마 위상 제거
            phase_sigma = gamma_phase_fn(sigma, t)
            phase_half  = gamma_phase_fn(0.5, t)
            delta_phase  = phase_sigma - phase_half
            # Q^s 의 σ-의존 위상도 제거
            delta_phase += (sigma - 0.5) * np.log(Q) * t / abs(t) if abs(t) > 0.01 else 0
            # 전체 보정: 위상이 σ=0.5 기준과 같아지도록
            lam = lam * np.exp(-1j * delta_phase)
        vals.append(lam.real)
    vals = np.array(vals)
    return int(np.sum(np.diff(np.sign(vals)) != 0))

def run_sigma_sweep(name, log_gamma_fn, gamma_phase_fn, Q, eps,
                    cn_nz, log_n_nz, cn_vals_nz):
    """하나의 L-함수에 대해 원본/보정 σ-유일성 측정"""
    t_arr = np.arange(T_MIN, T_MAX, 1.0)

    log(f"\n  {'σ':>5s} | {'원본':>6s} | {'보정':>6s} | {'Δ(γ위상 범위)':>14s}")
    log(f"  {'─'*5}─┼─{'─'*6}─┼─{'─'*6}─┼─{'─'*14}")

    orig_sc = {}
    corr_sc = {}

    for sigma in SIGMA_VALS:
        # 원본
        sc_orig = measure_sign_changes(sigma, t_arr, log_gamma_fn, Q, eps,
                                       cn_nz, log_n_nz, cn_vals_nz,
                                       gamma_phase_fn, correct=False)
        # 보정
        sc_corr = measure_sign_changes(sigma, t_arr, log_gamma_fn, Q, eps,
                                       cn_nz, log_n_nz, cn_vals_nz,
                                       gamma_phase_fn, correct=True)

        # 감마 위상 범위 (σ vs 0.5의 위상 차이)
        phases = [gamma_phase_fn(sigma, t) - gamma_phase_fn(0.5, t) for t in t_arr]
        phase_range = max(phases) - min(phases)

        orig_sc[sigma] = sc_orig
        corr_sc[sigma] = sc_corr
        log(f"  {sigma:5.2f} | {sc_orig:6d} | {sc_corr:6d} | {phase_range:14.4f}")

    # 판정
    orig_max = max(orig_sc.values())
    corr_max = max(corr_sc.values())
    orig_half_is_max = orig_sc[0.5] == orig_max
    corr_half_is_max = corr_sc[0.5] == corr_max

    log(f"\n  원본: σ=0.5 부호변환={orig_sc[0.5]}, 최대={orig_max} → "
        f"{'PASS ✅' if orig_half_is_max else 'FAIL ❌'}")
    log(f"  보정: σ=0.5 부호변환={corr_sc[0.5]}, 최대={corr_max} → "
        f"{'PASS ✅' if corr_half_is_max else 'FAIL ❌'}")

    if not orig_half_is_max and corr_half_is_max:
        log(f"  ⭐⭐ 복원 성공! 감마 인자 역보정으로 σ-유일성 PASS 달성")
    elif not orig_half_is_max and not corr_half_is_max:
        log(f"  ⚠️ 복원 실패. 감마 위상만으로는 설명 불충분 — 진폭 효과도 있을 수 있음")
    else:
        log(f"  (원본이 이미 PASS)")

    return orig_sc, corr_sc


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 70)
log("결과 #66 — 감마 인자 역보정(Deconvolution) 실험")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"가설: σ-유일성 FAIL = 감마 인자의 σ-의존 위상이 만든 '매끄러움'")
log(f"테스트: 감마 위상을 σ=0.5 기준으로 역보정 → σ-유일성 복원 시도")
log(f"σ 스윕: {SIGMA_VALS}")
log(f"t 범위: [{T_MIN}, {T_MAX}], Δt=1.0")
log()

t_total = time.time()

# ── 계수 준비 ──
N_COEFF = 500
an, ap, primes = compute_11a1_an(N_COEFF)

# GL(3) sym²(11a1) 계수
cn3 = compute_sym2_cn(an, primes, N_COEFF)
nz3 = np.where(np.abs(cn3[1:]) > 1e-15)[0]
cn3_nz = cn3[nz3 + 1]
log3_n_nz = np.log((nz3 + 1).astype(float))

# GL(2) 11a1 계수
an2 = np.array([0] + [an[i] if i < len(an) else 0 for i in range(1, N_COEFF+1)], dtype=float)
nz2 = np.where(np.abs(an2[1:]) > 1e-15)[0]
cn2_nz = an2[nz2 + 1]
log2_n_nz = np.log((nz2 + 1).astype(float))

log(f"GL(3) sym²(11a1): {len(nz3)}개 비영 계수")
log(f"GL(2) 11a1: {len(nz2)}개 비영 계수")
flush()

# ━━━ Phase 1: 감마 위상 σ-의존성 분석 ━━━
log(f"\n{'='*70}")
log(f"Phase 1: 감마 인자 위상의 σ-의존성 분석")
log(f"{'='*70}")

log(f"\n[감마 위상 차이: arg γ(σ+it) − arg γ(0.5+it)]")
log(f"  t=10 기준:")
for label, gfn in [("GL(1)", gl1_gamma_phase),
                    ("GL(2)", gl2_gamma_phase),
                    ("GL(3)", gl3_gamma_phase)]:
    phases = {s: gfn(s, 10.0) - gfn(0.5, 10.0) for s in SIGMA_VALS}
    log(f"  {label}: " + ", ".join(f"σ={s:.1f}:{phases[s]:+.4f}" for s in [0.3, 0.5, 0.7, 0.9]))

log(f"\n  t=30 기준:")
for label, gfn in [("GL(1)", gl1_gamma_phase),
                    ("GL(2)", gl2_gamma_phase),
                    ("GL(3)", gl3_gamma_phase)]:
    phases = {s: gfn(s, 30.0) - gfn(0.5, 30.0) for s in SIGMA_VALS}
    log(f"  {label}: " + ", ".join(f"σ={s:.1f}:{phases[s]:+.4f}" for s in [0.3, 0.5, 0.7, 0.9]))
flush()

# ━━━ Phase 2: GL(3) σ-유일성 원본 vs 보정 ━━━
log(f"\n{'='*70}")
log(f"Phase 2: GL(3) sym²(11a1) — σ-유일성 원본 vs 감마 역보정")
log(f"{'='*70}")
log(f"  N=121, Q=11, ε=+1, γ=Γ_ℝ(s)·Γ_ℝ(s+1)²")

t2 = time.time()
gl3_orig, gl3_corr = run_sigma_sweep(
    "GL(3) sym²(11a1)", gl3_log_gamma, gl3_gamma_phase,
    11.0, 1, cn3_nz, log3_n_nz, cn3_nz)
log(f"  소요: {time.time()-t2:.1f}초")
flush()

# ━━━ Phase 3: GL(2) σ-유일성 원본 vs 보정 ━━━
# GL(2) 11a1: 함수방정식 Λ(s) = Λ(2-s) → s→2-s 대칭
# Mellin-Barnes AFE에서 s1 = 2-s 로 수정 필요
# 간단히: GL(2)는 s→2-s이므로 별도 AFE 래퍼 사용

def Lambda_AFE_gl2(s, cn_nz, log_n_nz, cn_vals_nz, c=C_SHIFT):
    """GL(2) AFE: s→2-s 대칭"""
    s = complex(s)
    s1 = 2.0 - s  # GL(2): 2-s
    ec2 = np.exp(c**2)
    total = 0.0 + 0.0j
    for k in range(N_GH):
        y = GH_X[k]; w = GH_W[k]
        wk = c + 1j*y
        phase = np.exp(2j * c * y) / wk
        u_fwd = s + wk; u_bwd = s1 + wk
        lam_f = Lambda_direct(u_fwd, gl2_log_gamma, np.sqrt(11.0)/(2*np.pi),
                              cn_nz, log_n_nz, cn_vals_nz)
        lam_b = Lambda_direct(u_bwd, gl2_log_gamma, np.sqrt(11.0)/(2*np.pi),
                              cn_nz, log_n_nz, cn_vals_nz)
        total += w * (lam_f + 1 * lam_b) * phase  # ε=+1 for 11a1
    return ec2 / (2*np.pi) * total

def measure_sc_gl2(sigma, t_arr, gamma_phase_fn, correct=False):
    vals = []
    for t in t_arr:
        # GL(2) 임계선: σ=1 (s→2-s 대칭)
        lam = Lambda_AFE_gl2(complex(sigma, t), cn2_nz, log2_n_nz, cn2_nz)
        if correct:
            phase_sigma = gamma_phase_fn(sigma, t)
            phase_one   = gamma_phase_fn(1.0, t)  # GL(2) 임계선 σ=1
            delta_phase  = phase_sigma - phase_one
            Q_gl2 = np.sqrt(11.0)/(2*np.pi)
            delta_phase += (sigma - 1.0) * np.log(Q_gl2) * t / abs(t) if abs(t) > 0.01 else 0
            lam = lam * np.exp(-1j * delta_phase)
        vals.append(lam.real)
    return int(np.sum(np.diff(np.sign(np.array(vals))) != 0))

log(f"\n{'='*70}")
log(f"Phase 3: GL(2) 11a1 — σ-유일성 원본 vs 감마 역보정")
log(f"{'='*70}")
log(f"  N=11, Q=√11/(2π), ε=+1, γ=Γ_ℝ(s)·Γ_ℝ(s+1), 임계선 σ=1.0")

t3 = time.time()
gl2_sigma_vals = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5]  # GL(2): σ=1 중심
t_arr = np.arange(T_MIN, T_MAX, 1.0)

log(f"\n  {'σ':>5s} | {'원본':>6s} | {'보정':>6s}")
log(f"  {'─'*5}─┼─{'─'*6}─┼─{'─'*6}")

gl2_orig = {}; gl2_corr = {}
for sigma in gl2_sigma_vals:
    sc_o = measure_sc_gl2(sigma, t_arr, gl2_gamma_phase, correct=False)
    sc_c = measure_sc_gl2(sigma, t_arr, gl2_gamma_phase, correct=True)
    gl2_orig[sigma] = sc_o; gl2_corr[sigma] = sc_c
    log(f"  {sigma:5.2f} | {sc_o:6d} | {sc_c:6d}")

gl2_o_max = max(gl2_orig.values())
gl2_c_max = max(gl2_corr.values())
gl2_o_pass = gl2_orig.get(1.0, 0) == gl2_o_max
gl2_c_pass = gl2_corr.get(1.0, 0) == gl2_c_max

log(f"\n  원본: σ=1.0 부호변환={gl2_orig.get(1.0,0)}, 최대={gl2_o_max} → "
    f"{'PASS ✅' if gl2_o_pass else 'FAIL ❌'}")
log(f"  보정: σ=1.0 부호변환={gl2_corr.get(1.0,0)}, 최대={gl2_c_max} → "
    f"{'PASS ✅' if gl2_c_pass else 'FAIL ❌'}")

if not gl2_o_pass and gl2_c_pass:
    log(f"  ⭐⭐ GL(2)에서도 복원 성공!")
elif not gl2_o_pass and not gl2_c_pass:
    log(f"  ⚠️ GL(2)에서도 복원 실패.")

log(f"  소요: {time.time()-t3:.1f}초")
flush()

# ━━━ Phase 4: 종합 비교 ━━━
log(f"\n{'='*70}")
log(f"종합 — 감마 역보정 실험 결과")
log(f"{'='*70}")
log()

# GL(3) 결과
g3_o = "FAIL" if not (gl3_orig[0.5] == max(gl3_orig.values())) else "PASS"
g3_c = "FAIL" if not (gl3_corr[0.5] == max(gl3_corr.values())) else "PASS"

# GL(2) 결과
g2_o = "PASS" if gl2_o_pass else "FAIL"
g2_c = "PASS" if gl2_c_pass else "FAIL"

log(f"  | {'L-함수':15s} | {'원본':6s} | {'보정':6s} | {'복원':6s} |")
log(f"  |{'─'*17}|{'─'*8}|{'─'*8}|{'─'*8}|")
log(f"  | {'GL(3) sym²11a1':15s} | {g3_o:6s} | {g3_c:6s} | {'✅ 성공' if g3_o=='FAIL' and g3_c=='PASS' else '❌ 실패' if g3_o=='FAIL' else '—':6s} |")
log(f"  | {'GL(2) 11a1':15s} | {g2_o:6s} | {g2_c:6s} | {'✅ 성공' if g2_o=='FAIL' and g2_c=='PASS' else '❌ 실패' if g2_o=='FAIL' else '—':6s} |")

log()
if g3_o == "FAIL" and g3_c == "PASS":
    log("  ⭐⭐⭐ 핵심 결론: 감마 인자의 σ-의존 위상이 σ-유일성 FAIL의 원인!")
    log("  → 감마 인자 = '렌즈'. 역보정 = '안경'. 영점 신호는 항상 σ=0.5에 있었다.")
    log("  → 프레임워크 한계가 아니라 측정 방법의 문제였음을 증명.")
elif g3_o == "FAIL" and g3_c == "FAIL":
    log("  결론: 감마 위상만으로는 σ-유일성 FAIL을 설명할 수 없음.")
    log("  → 감마 인자의 진폭(magnitude) 효과도 기여하거나,")
    log("  →   영점 구조 자체가 σ=0.5를 선호하지 않을 가능성.")
    log("  → 추가 실험: 진폭 보정 또는 정규화된 함수 Z(s) = Λ(s)/|Λ(s)| 위상 분석.")
else:
    log("  원본이 이미 PASS — 역보정 불필요.")

total_time = time.time() - t_total
log(f"\n  총 소요: {total_time:.1f}초 ({total_time/60:.1f}분)")
log(f"  시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
flush()
log(f"\n결과 저장: {OUTFILE}")
