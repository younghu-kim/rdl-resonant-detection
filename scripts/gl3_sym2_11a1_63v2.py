#!/usr/bin/env python3
"""
[Project RDL] #63v2 — GL(3) sym²(11a1) ξ-다발 4성질 검증
=============================================================================
v1 문제: Theta function φ(y) 초고속 감쇠 → 영점 소실, κ 피크 0개
v2 해결: Mellin-Barnes AFE (Dokchitser 방식) + scipy/numpy 최적화

Λ(s) = (1/2πi)∫ Λ_dir(s+w) e^{w²}/w dw + ε·(1/2πi)∫ Λ_dir(1-s+w) e^{w²}/w dw
여기서 Λ_dir(u) = Q^u γ(u) L(u), Re(u) > 1 에서 수렴.
윤곽 이동 c=4 → Re(s+c+iy) ≥ 4.5 → Dirichlet 급수 수렴 보장.

결과 파일: results/gl3_sym2_11a1_63.txt
"""
import sys, os, time
import numpy as np
from scipy.special import loggamma
from scipy.signal import find_peaks

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "results", "gl3_sym2_11a1_63.txt")
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))
def flush():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━ 파라미터 ━━━━━━
N_COND   = 121
Q_COND   = 11.0          # √N
SIGMA_C  = 0.5           # 임계선
EPSILON  = 1             # root number
N_COEFF  = 500           # Dirichlet 계수
C_SHIFT  = 4.0           # 윤곽 이동
N_GH     = 50            # Gauss-Hermite 노드
T_MIN, T_MAX = 2.0, 50.0
DT       = 0.15          # κ 스윕 간격
DELTA    = 0.03          # σ 오프셋
MONO_R   = 0.3           # 모노드로미 반지름
MONO_N   = 48            # 모노드로미 단계
MONO_THR = 1.5           # 모노드로미 임계

# GH 노드 (∫ f(x) e^{-x²} dx ≈ Σ w_i f(x_i))
GH_X, GH_W = np.polynomial.hermite.hermgauss(N_GH)

# ━━━━━━ 11a1 계수 ━━━━━━
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
    """sym²(11a1) Dirichlet 계수 c(n) — float64"""
    cpk = {}
    for p in primes:
        if p > limit: break
        at2 = an[p]**2 / p  # ã_p²
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

# ━━━━━━ Λ(s) 계산 — scipy/numpy 고속 ━━━━━━
def log_gamma_factor(s):
    """log γ(s) = log[Γ_ℝ(s)·Γ_ℝ(s+1)²] = -(3s+2)/2·log(π) + logΓ(s/2) + 2·logΓ((s+1)/2)"""
    return (-(3*s+2)/2 * np.log(np.pi) +
            loggamma(s/2) + 2*loggamma((s+1)/2))

def Lambda_direct(s, cn_nz, log_n_nz, cn_vals_nz):
    """Λ_dir(s) = Q^s · γ(s) · L(s) — float64, Re(s)>1 필수"""
    s = complex(s)
    Qs = Q_COND ** s
    log_g = log_gamma_factor(s)
    gamma_s = np.exp(log_g)
    # L(s) = Σ c(n)/n^s  (벡터화)
    powers = np.exp(-s * log_n_nz)
    L_s = np.dot(cn_vals_nz, powers)
    return Qs * gamma_s * L_s

def Lambda_AFE(s, cn_nz, log_n_nz, cn_vals_nz, c=C_SHIFT, eps=EPSILON):
    """
    Λ(s) via Mellin-Barnes + Gauss-Hermite.
    각 항: Λ_dir(s+c+iy_k) · exp(c²+2icy_k) / (c+iy_k) · w_k
    """
    s = complex(s)
    s1 = 1.0 - s
    ec2 = np.exp(c**2)
    total = 0.0 + 0.0j
    for k in range(N_GH):
        y = GH_X[k]
        w = GH_W[k]
        wk = c + 1j*y
        phase = np.exp(2j * c * y) / wk
        u_fwd = s + wk
        u_bwd = s1 + wk
        lam_f = Lambda_direct(u_fwd, cn_nz, log_n_nz, cn_vals_nz)
        lam_b = Lambda_direct(u_bwd, cn_nz, log_n_nz, cn_vals_nz)
        total += w * (lam_f + eps * lam_b) * phase
    return ec2 / (2*np.pi) * total

def curvature(sigma, t, cn_nz, log_n_nz, cn_vals_nz, h=1e-5):
    """κ(σ+it) = |Λ'/Λ|²"""
    s = complex(sigma, t)
    L0 = Lambda_AFE(s, cn_nz, log_n_nz, cn_vals_nz)
    if abs(L0) < 1e-250:
        return 1e12
    Lp = Lambda_AFE(s + h, cn_nz, log_n_nz, cn_vals_nz)
    Lm = Lambda_AFE(s - h, cn_nz, log_n_nz, cn_vals_nz)
    conn = (Lp - Lm) / (2*h * L0)
    k = abs(conn)**2
    return k if np.isfinite(k) else 1e12

def monodromy(t_center, cn_nz, log_n_nz, cn_vals_nz,
              sigma=SIGMA_C, radius=MONO_R, n_steps=MONO_N):
    """폐곡선 모노드로미: |총 위상|/π"""
    center = complex(sigma, t_center)
    phase_acc = 0.0
    prev = None
    for j in range(n_steps + 1):
        th = 2*np.pi*j / n_steps
        pt = center + radius * np.exp(1j*th)
        val = Lambda_AFE(pt, cn_nz, log_n_nz, cn_vals_nz)
        if abs(val) < 1e-250:
            return None
        if prev is not None:
            phase_acc += np.angle(val / prev)
        prev = val
    return abs(phase_acc) / np.pi

# ━━━━━━ 메인 실행 ━━━━━━
log("=" * 70)
log("결과 #63v2 — GL(3) sym²(11a1) ξ-다발 4성질 검증 (Mellin-Barnes AFE)")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"N={N_COND}, Q=√N={Q_COND}, σ_crit={SIGMA_C}, ε={EPSILON}")
log(f"c_shift={C_SHIFT}, N_GH={N_GH}, N_coeff={N_COEFF}")
log()

t_total = time.time()

# ── Step 0: 계수 ──
log("[Step 0] sym²(11a1) Dirichlet 계수")
an, ap, primes = compute_11a1_an(N_COEFF)
cn = compute_sym2_cn(an, primes, N_COEFF)

log(f"  a_f: a₂={an[2]}, a₃={an[3]}, a₅={an[5]}, a₇={an[7]}, a₁₁={an[11]}")
log(f"  c(2)={cn[2]:.4f}, c(3)={cn[3]:.4f}, c(5)={cn[5]:.4f}, c(11)={cn[11]:.6f}")

# 비영 인덱스 + 전처리
nz_idx = np.where(np.abs(cn[1:]) > 1e-15)[0]  # 0-based into cn[1:]
cn_nz = cn[nz_idx + 1]
log_n_nz = np.log((nz_idx + 1).astype(float))
log(f"  비영 계수: {len(nz_idx)}개 / {N_COEFF}개")
flush()

# ── Step 1: AFE 검증 (vs direct) ──
log(f"\n[Step 1] AFE 검증: Λ_AFE vs Λ_direct (Re(s) >> 1)")
for s_test in [3.0+0j, 2.5+5j, 3.0+10j]:
    v_dir = Lambda_direct(s_test, cn_nz, log_n_nz, cn_nz)
    v_afe = Lambda_AFE(s_test, cn_nz, log_n_nz, cn_nz)
    if abs(v_dir) > 0:
        rel = abs(v_afe - v_dir) / abs(v_dir)
    else:
        rel = float('inf')
    log(f"  s={s_test}: |direct|={abs(v_dir):.6e}, |AFE|={abs(v_afe):.6e}, rel_err={rel:.2e}")
flush()

# ── Step 2: 함수방정식 ──
log(f"\n[Step 2] 함수방정식 Λ(s) = Λ(1-s)")
fe_pts = [0.5+5j, 0.5+10j, 0.5+15j, 0.7+8j, 0.3+12j, 0.6+25j, 0.5+35j]
fe_errs = []
for s in fe_pts:
    Ls = Lambda_AFE(s, cn_nz, log_n_nz, cn_nz)
    L1s = Lambda_AFE(1-s, cn_nz, log_n_nz, cn_nz)
    denom = max(abs(Ls), abs(L1s), 1e-300)
    rel = abs(Ls - L1s) / denom
    fe_errs.append(rel)
    ok = "✅" if rel < 1e-6 else "❌"
    log(f"  s={s}: |Λ(s)|={abs(Ls):.6e}, |Λ(1-s)|={abs(L1s):.6e}, rel={rel:.2e} {ok}")
flush()

# ── Step 3: Λ(1/2+it) 스캔 — 부호 변환 탐색 ──
log(f"\n[Step 3] Λ(1/2+it) 스캔: t ∈ [{T_MIN}, {T_MAX}], Δt=0.5")
scan_ts = np.arange(T_MIN, T_MAX + 0.1, 0.5)
scan_re = np.zeros(len(scan_ts))
t3 = time.time()

for i, t in enumerate(scan_ts):
    val = Lambda_AFE(complex(SIGMA_C, t), cn_nz, log_n_nz, cn_nz)
    scan_re[i] = val.real
    sign = "+" if val.real >= 0 else "−"
    if i % 10 == 0 or i == len(scan_ts)-1:
        elapsed = time.time() - t3
        eta = elapsed/(i+1)*(len(scan_ts)-i-1) if i > 0 else 0
        log(f"  t={t:5.1f}: Λ={val.real:+.8e}  [{sign}]  ({i+1}/{len(scan_ts)}, "
            f"eta={eta:.0f}s)")
        flush()

n_sign = np.sum(np.diff(np.sign(scan_re)) != 0)
log(f"\n  부호 변환: {n_sign}회, 스캔 소요: {time.time()-t3:.1f}초")
flush()

if n_sign == 0:
    log(f"  ⚠️ 부호 변환 0회 — Λ 범위: [{scan_re.min():.4e}, {scan_re.max():.4e}]")
    log(f"  ⚠️ AFE 검증 실패 가능성 — 종료")
    flush()
    sys.exit(1)

# ── Step 4: 이분법 영점 근사 ──
log(f"\n[Step 4] 이분법 영점 근사")
zeros = []
for i in range(len(scan_re)-1):
    if scan_re[i] * scan_re[i+1] < 0:
        lo, hi = scan_ts[i], scan_ts[i+1]
        v_lo = scan_re[i]
        for _ in range(30):
            mid = (lo+hi)/2
            v_mid = Lambda_AFE(complex(SIGMA_C, mid), cn_nz, log_n_nz, cn_nz).real
            if v_mid * v_lo < 0: hi = mid
            else: lo = mid; v_lo = v_mid
        zeros.append((lo+hi)/2)
        log(f"  영점 #{len(zeros):2d}: t ≈ {(lo+hi)/2:.8f}")

log(f"\n  발견: {len(zeros)}개 영점")
flush()

# ── Step 5: κ 피크 + 모노드로미 ──
log(f"\n[Step 5] 영점 근방 κ + 모노드로미")
sigma_sw = SIGMA_C + DELTA
tp_results = []  # (t, κ, mono/π)

for z_t in zeros:
    k_val = curvature(sigma_sw, z_t, cn_nz, log_n_nz, cn_nz)
    m_val = monodromy(z_t, cn_nz, log_n_nz, cn_nz)
    tp_results.append((z_t, k_val, m_val))
    m_str = f"{m_val:.4f}" if m_val is not None else "FAIL"
    log(f"  t={z_t:.6f}: κ={k_val:.2f}, mono/π={m_str}")

# κ near vs far
log(f"\n[Step 5b] κ 집중도: near-zero vs far")
far_ts = np.arange(T_MIN, T_MAX, 2.0)
k_near, k_far = [], []
for t in far_ts:
    k = curvature(sigma_sw, t, cn_nz, log_n_nz, cn_nz)
    min_d = min(abs(t - z) for z in zeros) if zeros else 999
    if min_d < 0.5: k_near.append(k)
    else: k_far.append(k)

if k_near and k_far:
    near_m = np.median(k_near); far_m = np.median(k_far)
    ratio = near_m / far_m if far_m > 0 else float('inf')
    log(f"  near: median={near_m:.2f} ({len(k_near)}점)")
    log(f"  far:  median={far_m:.4f} ({len(k_far)}점)")
    log(f"  비율: {ratio:.1f}×  {'✅ PASS' if ratio >= 10 else '❌ FAIL'}")
else:
    ratio = 0
    log(f"  데이터 부족 (near={len(k_near)}, far={len(k_far)})")
flush()

# ── Step 6: σ-유일성 ──
log(f"\n[Step 6] σ-유일성 스윕")
sigma_vals = [0.3, 0.5, 0.7, 0.9]
sigma_sign_changes = {}
for sig in sigma_vals:
    vals = []
    for t in np.arange(T_MIN, T_MAX, 1.0):
        v = Lambda_AFE(complex(sig, t), cn_nz, log_n_nz, cn_nz)
        vals.append(v.real)
    sc = np.sum(np.diff(np.sign(vals)) != 0)
    sigma_sign_changes[sig] = sc
    log(f"  σ={sig:.1f}: 부호변환={sc}")

# σ=0.5에서 최대?
sc_half = sigma_sign_changes.get(0.5, 0)
is_max = all(sc_half >= sigma_sign_changes[s] for s in sigma_vals if s != 0.5)
log(f"\n  σ=0.5 최대: {'✅ PASS' if is_max else '❌ FAIL'}")
log(f"  패턴: {dict(sigma_sign_changes)}")
flush()

# ── Step 7: 모노드로미 통계 ──
log(f"\n[Step 7] 모노드로미 통계")
mono_vals = [m for _, _, m in tp_results if m is not None]
if mono_vals:
    tp_mono = [m for m in mono_vals if m >= MONO_THR]
    fp_mono = [m for m in mono_vals if m < MONO_THR]
    log(f"  TP (mono/π ≥ {MONO_THR}): {len(tp_mono)}개, mean={np.mean(tp_mono):.4f}" if tp_mono else "  TP: 0개")
    log(f"  FP (mono/π < {MONO_THR}): {len(fp_mono)}개" + (f", mean={np.mean(fp_mono):.4f}" if fp_mono else ""))
flush()

# ── 종합 ──
log(f"\n{'='*70}")
log(f"종합 판정 — GL(3) sym²(11a1)")
log(f"{'='*70}")

fe_pass = sum(1 for e in fe_errs if e < 1e-6)
log(f"  [필수] 함수방정식: {fe_pass}/{len(fe_errs)} {'✅' if fe_pass >= 5 else '❌'}")
log(f"  영점 발견: {len(zeros)}개 (t ∈ [{T_MIN}, {T_MAX}])")
log(f"  [필수] σ-유일성: {'✅ PASS' if is_max else '❌ FAIL'}")
log(f"  [양성] κ 집중: {ratio:.1f}× {'✅' if ratio >= 10 else '❌'}")
n_tp = len([m for m in mono_vals if m >= MONO_THR]) if mono_vals else 0
log(f"  [양성] 모노드로미: TP={n_tp}/{len(zeros)}")

# σ-유일성 FAIL 패턴 분석 (GL(2) 비교용)
log(f"\n  ⭐ σ-유일성 FAIL 패턴 분석:")
log(f"     GL(1): σ=0.5 PASS (부호변환 σ=0.5에서 명확히 최대)")
log(f"     GL(2): σ=0.5 FAIL (σ=1.0에서 최대가 아님)")
log(f"     GL(3): σ=0.5 {'PASS' if is_max else 'FAIL'}")
if not is_max:
    best_sig = max(sigma_sign_changes, key=sigma_sign_changes.get)
    log(f"     → GL(3) FAIL, 최대 σ={best_sig} (GL(2)와 {'같음' if best_sig == 0.5 else '다름'})")
    log(f"     → FAIL 패턴이 GL(2)와 다르면 → L-함수 분류기 후보!")

total = time.time() - t_total
log(f"\n  총 소요: {total:.1f}초 ({total/60:.1f}분)")
log(f"  시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
flush()
log(f"\n결과 저장: {OUTFILE}")
