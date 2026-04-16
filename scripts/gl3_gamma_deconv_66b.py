#!/usr/bin/env python3
"""
[Project RDL] #66b — 감마 인자 전체 역보정 (위상 + 진폭)
=============================================================================
#66 결과: 위상만 보정 → σ-유일성 복원 실패
#66b 시도: 감마 인자 전체(위상+진폭)를 제거하고 L(s) 수준에서 부호변환 측정

방법:
  Level 1: Λ(s) = Q^s · γ(s) · L(s)           ← 원본 (FAIL)
  Level 2: Λ(s) / |γ(s)|                       ← 진폭만 제거
  Level 3: Λ(s) / γ(s) = Q^s · L(s)           ← 감마 전체 제거
  Level 4: Λ(s) / (Q^s · γ(s)) = L(s)         ← 전부 제거 (순수 L-함수)
  Level 5: Z(s) = Λ(s)/|Λ(s)| (위상만 추출)   ← 크기 무관 위상 분석

핵심 질문: 어느 수준에서 σ=0.5의 특이성이 나타나는가?
"""
import sys, os, time
import numpy as np
from scipy.special import loggamma

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "results", "gl3_gamma_deconv_66b.txt")
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))
def flush():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━ 파라미터 ━━━━━━
C_SHIFT = 4.0; N_GH = 50
GH_X, GH_W = np.polynomial.hermite.hermgauss(N_GH)
T_MIN, T_MAX = 2.0, 50.0
SIGMA_VALS = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
N_COND = 121; Q_COND = 11.0; EPSILON = 1

# ━━━━━━ 감마 인자 ━━━━━━
def log_gamma_gl3(s):
    return (-(3*s+2)/2 * np.log(np.pi) +
            loggamma(s/2) + 2*loggamma((s+1)/2))

# ━━━━━━ 11a1 + sym² 계수 ━━━━━━
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
    cn = np.zeros(limit+1); cn[1] = 1.0
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

# ━━━━━━ AFE ━━━━━━
def Lambda_direct(s, cn_nz, log_n_nz, cn_vals_nz):
    s = complex(s)
    Qs = Q_COND ** s
    gamma_s = np.exp(log_gamma_gl3(s))
    powers = np.exp(-s * log_n_nz)
    L_s = np.dot(cn_vals_nz, powers)
    return Qs * gamma_s * L_s

def Lambda_AFE(s, cn_nz, log_n_nz, cn_vals_nz, c=C_SHIFT):
    s = complex(s); s1 = 1.0 - s
    ec2 = np.exp(c**2)
    total = 0.0 + 0.0j
    for k in range(N_GH):
        y = GH_X[k]; w = GH_W[k]
        wk = c + 1j*y
        phase = np.exp(2j * c * y) / wk
        u_fwd = s + wk; u_bwd = s1 + wk
        lam_f = Lambda_direct(u_fwd, cn_nz, log_n_nz, cn_vals_nz)
        lam_b = Lambda_direct(u_bwd, cn_nz, log_n_nz, cn_vals_nz)
        total += w * (lam_f + EPSILON * lam_b) * phase
    return ec2 / (2*np.pi) * total

# ━━━━━━ 다단계 보정 ━━━━━━
def compute_levels(sigma, t, cn_nz, log_n_nz, cn_vals_nz):
    """
    5단계 보정값 반환.
    Level 1: Λ(s) 원본
    Level 2: Λ(s)/|γ(s)| (진폭 제거)
    Level 3: Λ(s)/γ(s) = Q^s · L(s)
    Level 4: Λ(s)/(Q^s · γ(s)) = L(s)
    Level 5: Λ(s)/|Λ(s)| (위상만)
    """
    s = complex(sigma, t)
    lam = Lambda_AFE(s, cn_nz, log_n_nz, cn_vals_nz)
    log_g = log_gamma_gl3(s)
    gamma_s = np.exp(log_g)
    Qs = Q_COND ** s

    lv1 = lam                                          # 원본
    lv2 = lam / abs(gamma_s) if abs(gamma_s) > 1e-300 else lam  # 진폭만 제거
    lv3 = lam / gamma_s if abs(gamma_s) > 1e-300 else lam       # 감마 전체 제거
    lv4 = lam / (Qs * gamma_s) if abs(Qs * gamma_s) > 1e-300 else lam  # L(s)
    lv5 = lam / abs(lam) if abs(lam) > 1e-300 else lam          # 위상만

    return lv1, lv2, lv3, lv4, lv5

def count_sign_changes(vals):
    arr = np.array([v.real for v in vals])
    return int(np.sum(np.diff(np.sign(arr)) != 0))

# ━━━━━━ 메인 ━━━━━━
log("=" * 70)
log("결과 #66b — 감마 인자 전체 역보정 (5-Level 분석)")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"GL(3) sym²(11a1), N=121, Q=11, ε=+1")
log()

t_total = time.time()

# 계수
N_COEFF = 500
an, ap, primes = compute_11a1_an(N_COEFF)
cn = compute_sym2_cn(an, primes, N_COEFF)
nz_idx = np.where(np.abs(cn[1:]) > 1e-15)[0]
cn_nz = cn[nz_idx + 1]
log_n_nz = np.log((nz_idx + 1).astype(float))
log(f"비영 계수: {len(nz_idx)}개 / {N_COEFF}개")
flush()

# ━━━ 5-Level σ-유일성 스윕 ━━━
log(f"\n{'='*70}")
log(f"5-Level σ-유일성 스윕")
log(f"{'='*70}")
log(f"  Level 1: Λ(s) 원본")
log(f"  Level 2: Λ(s)/|γ(s)| — 감마 진폭 제거")
log(f"  Level 3: Λ(s)/γ(s) = Q^s·L(s) — 감마 전체 제거")
log(f"  Level 4: Λ(s)/(Q^s·γ(s)) = L(s) — 순수 L-함수")
log(f"  Level 5: Λ(s)/|Λ(s)| — 위상만 추출")
log()

t_arr = np.arange(T_MIN, T_MAX, 1.0)
results = {lv: {} for lv in range(1, 6)}
level_names = {1: "Λ(s)", 2: "Λ/|γ|", 3: "Λ/γ", 4: "L(s)", 5: "Phase"}

for sigma in SIGMA_VALS:
    ts = time.time()
    vals = {lv: [] for lv in range(1, 6)}
    for t in t_arr:
        lv1, lv2, lv3, lv4, lv5 = compute_levels(sigma, t, cn_nz, log_n_nz, cn_nz)
        vals[1].append(lv1)
        vals[2].append(lv2)
        vals[3].append(lv3)
        vals[4].append(lv4)
        vals[5].append(lv5)

    for lv in range(1, 6):
        results[lv][sigma] = count_sign_changes(vals[lv])

    elapsed = time.time() - ts
    log(f"  σ={sigma:.2f}: " + " | ".join(f"Lv{lv}={results[lv][sigma]:2d}" for lv in range(1, 6))
        + f"  ({elapsed:.1f}s)")
    flush()

# ━━━ 결과 테이블 ━━━
log(f"\n{'='*70}")
log(f"결과 테이블 — 부호변환 수")
log(f"{'='*70}")
header = f"  {'σ':>5s}"
for lv in range(1, 6):
    header += f" | {level_names[lv]:>8s}"
log(header)
log(f"  {'─'*5}" + "─┼─".join(["─"*8]*5))

for sigma in SIGMA_VALS:
    row = f"  {sigma:5.2f}"
    for lv in range(1, 6):
        sc = results[lv][sigma]
        # σ=0.5에서 최대인지 표시
        is_max_at_half = results[lv][0.5] >= max(results[lv].values())
        row += f" | {sc:8d}"
    log(row)

# ━━━ 각 Level 판정 ━━━
log(f"\n{'='*70}")
log(f"각 Level σ-유일성 판정")
log(f"{'='*70}")

for lv in range(1, 6):
    sc_half = results[lv][0.5]
    sc_max = max(results[lv].values())
    best_sigma = [s for s in SIGMA_VALS if results[lv][s] == sc_max]
    is_pass = sc_half == sc_max
    log(f"  Level {lv} ({level_names[lv]:>8s}): σ=0.5 → {sc_half}, "
        f"최대 → {sc_max} (at σ={best_sigma}) "
        f"→ {'✅ PASS' if is_pass else '❌ FAIL'}")

# ━━━ 핵심 분석 ━━━
log(f"\n{'='*70}")
log(f"분석")
log(f"{'='*70}")

# 어떤 Level에서 복원되는가?
restored = None
for lv in range(2, 6):
    if results[lv][0.5] == max(results[lv].values()):
        restored = lv
        break

if restored:
    log(f"\n  ⭐⭐ Level {restored} ({level_names[restored]})에서 σ-유일성 복원!")
    log(f"  → Level 1~{restored-1}은 FAIL, Level {restored}~5는 PASS")
    if restored == 2:
        log(f"  → 감마 진폭이 주범. 위상은 무관.")
    elif restored == 3:
        log(f"  → 감마 전체(진폭+위상)를 제거해야 복원.")
        log(f"  → 진폭과 위상 모두 기여. 감마 인자 = '렌즈' 가설 부분 확인.")
    elif restored == 4:
        log(f"  → Q^s 인자까지 제거해야 복원.")
        log(f"  → conductor 스케일링도 관여.")
    elif restored == 5:
        log(f"  → 위상 정보만으로 복원. 크기 정보가 방해였음.")
else:
    log(f"\n  ⚠️ 어떤 Level에서도 σ-유일성 복원되지 않음.")
    log(f"  → 감마 인자/Q/진폭이 원인이 아님.")
    log(f"  → 가능한 해석:")
    log(f"    (a) GL(3) 영점 구조 자체가 σ=0.5를 선호하지 않음")
    log(f"    (b) Mellin-Barnes AFE의 수치적 한계")
    log(f"    (c) 부호변환 수가 아닌 다른 측도(예: κ 피크 높이)가 필요")
    log(f"  → 추가 실험: 부호변환 수 대신 위상 점프 크기의 σ-의존성 분석")

total_time = time.time() - t_total
log(f"\n  총 소요: {total_time:.1f}초 ({total_time/60:.1f}분)")
log(f"  시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
flush()
log(f"\n결과 저장: {OUTFILE}")
