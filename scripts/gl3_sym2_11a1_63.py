"""
=============================================================================
[Project RDL] 결과 #63 — GL(3) sym²(11a1) ξ-다발 4성질 검증 (v5)
=============================================================================
Degree-3 L-함수 AFE via Gauss-Hermite contour integral.

수학적 기반:
  sym²(11a1): conductor N=121, ε=+1, gamma shifts μ=[0,0,1]
  γ(s) = Γ_ℝ(s+1)² · Γ_ℝ(s+2) = π^{-(3s+4)/2} · Γ((s+1)/2)² · Γ((s+2)/2)
  해석적 정규화: 산술 shifts [0,0,1] → s→s+1 → 해석 shifts [1,1,2]
  Λ(s) = N^{s/2} · γ(s) · L(s),  FE: Λ(s) = ε·Λ(1-s)

AFE (Rubinstein smoothed):
  Λ(s) = Σ_n c(n) [Φ(s,n) + ε·Φ(1-s,n)]
  Φ(s,n) = (1/2πi) ∫ N^{(s+w)/2} γ(s+w) n^{-(s+w)} exp(w²)/(2w) dw
         = (e^{c²}/4π) Σ_k w_k · A_k(s) · n^{-(s+c+iv_k)}
  A_k(s) = N^{(s+c+iv_k)/2} · γ(s+c+iv_k) · e^{2icv_k} / (c+iv_k)
  (Gauss-Hermite 직교로 ∫ f(v)e^{-v²} dv 수치 적분)

  효율적 형태:
  Λ(s) = (e^{c²}/4π) Σ_k w_k [A_k(s) D(s+c+iv_k) + ε A_k(1-s) D(1-s+c+iv_k)]
  D(w) = Σ_n c(n)/n^w   (Re(w) >> 1에서 수렴)
=============================================================================
"""

import sys, os, time
import numpy as np
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

DPS = 50
mpmath.mp.dps = DPS

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "gl3_sym2_11a1_63.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

# ━━━━━━━━━━━ 파라미터 ━━━━━━━━━━━
N_COND = 121
N_COND_MP = mpmath.mpf(121)
EPSILON = 1
N_MAX_COEFF = 120   # Dirichlet 계수 수
CONTOUR_C = mpmath.mpf(2)  # contour shift Re(w)=c

# Gauss-Hermite 직교 (∫ f(v)e^{-v²} dv ≈ Σ w_k f(v_k))
N_HERMITE = 50
_herm_nodes, _herm_weights = np.polynomial.hermite.hermgauss(N_HERMITE)
HERM_NODES = [mpmath.mpf(float(x)) for x in _herm_nodes]
HERM_WEIGHTS = [mpmath.mpf(float(w)) for w in _herm_weights]

# LMFDB 영점 (해석적 정규화, 임계선 σ=1/2)
LMFDB_ZEROS = [3.899281, 4.734595, 6.189477, 7.312039, 8.650148,
               10.128219, 10.936767, 12.070610, 12.744576, 13.335280,
               14.419784, 15.299510, 16.008759, 16.756012, 17.687147]

# 스윕 파라미터
T_MIN, T_MAX = 2.0, 45.0
DELTA = 0.03
MONO_RADIUS, MONO_STEPS = 0.4, 64
MATCH_TOL = 0.5

# ━━━━━━━━━━━ 로깅 ━━━━━━━━━━━
lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))
def flush_to_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━━━━━━ 11a1 Hecke 고유값 ━━━━━━━━━━━

def compute_11a1_an(limit):
    """11a1: y²+y = x³-x²-10x-20 의 Hecke 고유값"""
    known_ap = {2:-2, 3:-1, 5:1, 7:-2, 11:1, 13:4, 17:-2, 19:0, 23:-1,
                29:0, 31:7, 37:3, 41:-8, 43:-6, 47:8, 53:-6, 59:5, 61:12,
                67:-7, 71:-4, 73:2, 79:10, 83:-7, 89:6, 97:-4, 101:-2,
                103:12, 107:10, 109:-4}

    sieve = [True]*(limit+1); sieve[0]=sieve[1]=False
    for i in range(2, int(limit**0.5)+1):
        if sieve[i]:
            for j in range(i*i, limit+1, i): sieve[j] = False
    primes = [i for i in range(2, limit+1) if sieve[i]]

    ap = {}
    for p in primes:
        if p in known_ap:
            ap[p] = known_ap[p]; continue
        # 점 계수로 계산
        count = 0
        for x in range(p):
            rhs = (x*x*x - x*x - 10*x - 20) % p
            if rhs == 0: count += 1
            elif pow(rhs, (p-1)//2, p) == 1: count += 2
        ap[p] = p - count

    # a_f(p^k) 재귀
    apk = {}
    for p in primes:
        apk[(p,0)] = 1; apk[(p,1)] = ap[p]
        pk = p; k = 1
        while pk*p <= limit:
            pk *= p; k += 1
            if p == 11:
                apk[(p,k)] = ap[p]**k
            else:
                apk[(p,k)] = ap[p]*apk[(p,k-1)] - p*apk[(p,k-2)]

    # 곱셈적 구조
    an = [0]*(limit+1); an[1] = 1
    for n in range(2, limit+1):
        temp = n; result = 1
        for p in primes:
            if p*p > temp: break
            if temp % p == 0:
                k = 0
                while temp % p == 0: k += 1; temp //= p
                result *= apk.get((p,k), 0)
        if temp > 1:
            result *= ap.get(temp, 0)
        an[n] = result
    return an, ap, primes


def compute_sym2_cn(an, primes, limit):
    """sym²(11a1) 해석적 정규화 Dirichlet 계수
    c(p) = a_f(p)²/p - 1  (good primes)
    재귀: c(p^k) = c₁·c(p^{k-1}) - c₁·c(p^{k-2}) + c(p^{k-3})
    """
    mpmath.mp.dps = DPS + 10
    cpk = {}
    for p in primes:
        if p > limit: break
        cpk[(p,0)] = mpmath.mpf(1)
        if p == 11:
            # bad prime: (1 - p^{-s})^{-1} → c(11^k) = 11^{-k}
            for k in range(1, 20):
                cpk[(p,k)] = mpmath.power(mpmath.mpf(11), -k)
                if 11**(k+1) > limit: break
        else:
            c1 = mpmath.mpf(an[p])**2 / mpmath.mpf(p) - 1
            cpk[(p,1)] = c1
            pk = p; k = 1
            while pk*p <= limit:
                pk *= p; k += 1
                bkm1 = cpk[(p,k-1)]
                bkm2 = cpk[(p,k-2)]
                bkm3 = cpk.get((p,k-3), mpmath.mpf(0))
                cpk[(p,k)] = c1*bkm1 - c1*bkm2 + bkm3

    cn = [mpmath.mpf(0)]*(limit+1)
    cn[1] = mpmath.mpf(1)
    for n in range(2, limit+1):
        temp = n; result = mpmath.mpf(1)
        for p in primes:
            if p*p > temp: break
            if temp % p == 0:
                k = 0
                while temp % p == 0: k += 1; temp //= p
                if (p,k) in cpk:
                    result *= cpk[(p,k)]
                else:
                    result = mpmath.mpf(0); break
        if temp > 1:
            if (temp,1) in cpk:
                result *= cpk[(temp,1)]
            else:
                result = mpmath.mpf(0)
        cn[n] = result
    mpmath.mp.dps = DPS
    return cn


# ━━━━━━━━━━━ 감마 인자 (올바른: shifts [0,0,1]) ━━━━━━━━━━━

def gamma_factor(s):
    """γ(s) = Γ_ℝ(s+1)² · Γ_ℝ(s+2)  (해석적 정규화 shifts μ=[1,1,2])
    = π^{-(3s+4)/2} · Γ((s+1)/2)² · Γ((s+2)/2)

    유도: 산술적 정규화 shifts [0,0,1] → s→s+1 이동 → 해석적 [1,1,2]
    GL(2) 11a1 (weight 2): Γ_ℂ(s) = Γ_ℝ(s)Γ_ℝ(s+1), shifts [0,1]
    sym²: 산술 Γ_ℝ(s)·Γ_ℂ(s) → shifts [0,0,1]
    해석: s→s+1 → Γ_ℝ(s+1)·Γ_ℂ(s+1) → shifts [1,1,2]
    """
    return (mpmath.power(mpmath.pi, -(3*s + 4) / 2)
            * mpmath.gamma((s + 1) / 2)**2
            * mpmath.gamma((s + 2) / 2))


def dirichlet_series(w, cn):
    """D(w) = Σ c(n)/n^w — Re(w) > 1에서 수렴"""
    total = mpmath.mpc(0)
    for n in range(1, len(cn)):
        if cn[n] == 0: continue
        total += cn[n] * mpmath.power(n, -w)
    return total


def Lambda_direct(s, cn):
    """Λ(s) = N^{s/2} γ(s) L(s) — Re(s) >> 1에서만 유효"""
    return mpmath.power(N_COND_MP, s/2) * gamma_factor(s) * dirichlet_series(s, cn)


# ━━━━━━━━━━━ AFE via Gauss-Hermite contour integral ━━━━━━━━━━━

def Lambda_AFE(s, cn, c=None):
    """
    Λ(s) = (e^{c²}/4π) Σ_k w_k [A_k(s)·D(s+c+iv_k) + ε·A_k(1-s)·D(1-s+c+iv_k)]

    A_k(s) = N^{(s+c+iv_k)/2} · γ(s+c+iv_k) · e^{2icv_k} / (c+iv_k)
    """
    if c is None:
        c = CONTOUR_C
    s_mp = mpmath.mpc(s)
    s1_mp = 1 - s_mp

    prefactor = mpmath.exp(c**2) / (2 * mpmath.pi)  # 1/w (not 1/(2w))

    total = mpmath.mpc(0)
    for k in range(N_HERMITE):
        v_k = HERM_NODES[k]
        wt_k = HERM_WEIGHTS[k]

        iv_k = mpmath.mpc(0, v_k)
        w_shift = c + iv_k  # = c + iv_k

        # A_k(s) = N^{(s+w)/2} · γ(s+w) · e^{2icv} / (c+iv)
        sw = s_mp + w_shift
        gamma_sw = gamma_factor(sw)
        N_pow_s = mpmath.power(N_COND_MP, sw / 2)
        exp_phase = mpmath.exp(2 * mpmath.mpc(0, 1) * c * v_k)
        A_s = N_pow_s * gamma_sw * exp_phase / w_shift

        # D(s + c + iv_k)
        D_s = dirichlet_series(sw, cn)

        # A_k(1-s) = N^{(1-s+w)/2} · γ(1-s+w) · same phase / same denom
        s1w = s1_mp + w_shift
        gamma_s1w = gamma_factor(s1w)
        N_pow_1s = mpmath.power(N_COND_MP, s1w / 2)
        A_1s = N_pow_1s * gamma_s1w * exp_phase / w_shift

        # D(1-s + c + iv_k)
        D_1s = dirichlet_series(s1w, cn)

        total += wt_k * (A_s * D_s + EPSILON * A_1s * D_1s)

    return prefactor * total


# ━━━━━━━━━━━ 곡률, 모노드로미 ━━━━━━━━━━━

def curvature(s, cn, h=1e-6):
    """κ(s) = |Λ'/Λ|² (접속의 곡률)"""
    s_mp = mpmath.mpc(s)
    h_mp = mpmath.mpf(h)
    try:
        L0 = Lambda_AFE(s_mp, cn)
        if abs(L0) < mpmath.mpf(10)**(-DPS + 15):
            return float('inf')
        Lp = Lambda_AFE(s_mp + h_mp, cn)
        Lm = Lambda_AFE(s_mp - h_mp, cn)
        conn = (Lp - Lm) / (2 * h_mp * L0)
        k = float(abs(conn)**2)
        return k if np.isfinite(k) else float('inf')
    except Exception as e:
        print(f"  WARNING curvature at s={s}: {e}", flush=True)
        return 0.0


def monodromy(t_center, cn, sigma=0.5, radius=None, n_steps=None):
    """폐곡선 적분으로 모노드로미 측정 — arg(Λ) 누적"""
    if radius is None: radius = MONO_RADIUS
    if n_steps is None: n_steps = MONO_STEPS
    center = mpmath.mpc(sigma, t_center)
    phase_accum = mpmath.mpf(0)
    prev_val = None

    for j in range(n_steps + 1):
        theta = 2 * mpmath.pi * j / n_steps
        pt = center + radius * mpmath.exp(1j * theta)
        try:
            val = Lambda_AFE(pt, cn)
            if abs(val) < mpmath.mpf(10)**(-DPS + 15):
                return None
        except Exception as e:
            print(f"  WARNING mono step {j}: {e}", flush=True)
            return None

        if prev_val is not None:
            ratio = val / prev_val
            if abs(ratio) < mpmath.mpf(10)**(-DPS + 15):
                return None
            phase_accum += mpmath.im(mpmath.log(ratio))
        prev_val = val

    return float(abs(phase_accum) / mpmath.pi)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#                   메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 70)
log("결과 #63 — GL(3) sym²(11a1) ξ-다발 4성질 (v5: Gauss-Hermite AFE)")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"N={N_COND}, ε={EPSILON}, gamma shifts μ=[1,1,2] (해석적 정규화)")
log(f"DPS={DPS}, N_coeff={N_MAX_COEFF}, N_hermite={N_HERMITE}, c={float(CONTOUR_C)}")
log()

t_total = time.time()

# ── 0. Dirichlet 계수 ──
log("[Step 0] sym²(11a1) Dirichlet 계수 (해석적 정규화)")
an, ap, primes = compute_11a1_an(N_MAX_COEFF)
cn = compute_sym2_cn(an, primes, N_MAX_COEFF)

log(f"  a_f: a₂={an[2]}, a₃={an[3]}, a₅={an[5]}, a₇={an[7]}, a₁₁={an[11]}, a₁₃={an[13]}")

# 검증: c(p) = a_f(p)²/p - 1 for good p
checks = [(2, 1.0), (3, -2/3), (5, -4/5), (7, -3/7), (13, 3/13)]
all_ok = True
for p, expected in checks:
    actual = float(cn[p])
    ok = abs(actual - expected) < 1e-10
    if not ok: all_ok = False
    log(f"  c({p})={actual:.8f} (기대={expected:.8f}) {'✅' if ok else '❌'}")

# c(11) = 1/11 (bad prime)
c11 = float(cn[11])
log(f"  c(11)={c11:.8f} (기대={1/11:.8f}) {'✅' if abs(c11 - 1/11)<1e-10 else '❌'}")

n_nz = sum(1 for i in range(1, len(cn)) if cn[i] != 0)
log(f"  비영 계수: {n_nz}/{N_MAX_COEFF}")
if not all_ok:
    log("  ⚠️ 계수 검증 실패!")
flush_to_file()

# ── 1. 벤치마크 ──
log(f"\n[Step 1] 벤치마크: Λ_AFE 한 점 소요 시간")
t1 = time.time()
val_bench = Lambda_AFE(mpmath.mpc(3, 0), cn)
dt_bench = time.time() - t1
log(f"  Λ_AFE(3) = {mpmath.nstr(val_bench, 12)} (소요: {dt_bench:.2f}s)")
flush_to_file()

# ── 2. AFE 검증: Λ_direct와 비교 (Re(s) >> 1) ──
log(f"\n[Step 2] AFE vs Λ_direct 검증 (Re(s) >> 1)")
log(f"  Λ_direct(s) = N^{{s/2}} γ(s) L(s) 는 Re(s)>1에서 수렴")
log(f"  AFE는 모든 s에서 유효 → 양쪽 비교로 검증")

test_points = [
    mpmath.mpf(3),
    mpmath.mpf(2.5),
    mpmath.mpc(3, 5),
    mpmath.mpc(2.5, 10),
    mpmath.mpc(2, 15),
]

fe_validated = True
for sp in test_points:
    t_s = time.time()
    val_dir = Lambda_direct(sp, cn)
    val_afe = Lambda_AFE(sp, cn)
    dt_s = time.time() - t_s

    if abs(val_dir) > 0:
        rel = float(abs(val_afe - val_dir) / abs(val_dir))
    else:
        rel = float('inf')
    ok = rel < 1e-6
    if not ok: fe_validated = False
    log(f"  s={mpmath.nstr(sp,6)}: |dir|={float(abs(val_dir)):.6e}, "
        f"|AFE|={float(abs(val_afe)):.6e}, rel={rel:.2e} "
        f"{'✅' if ok else '❌'} ({dt_s:.1f}s)")
    flush_to_file()

if fe_validated:
    log(f"\n  ✅ AFE 검증 통과 — 계수+감마+AFE 올바름")
else:
    log(f"\n  ⚠️ AFE 검증 실패! 디버그 진행...")
    # 감마 인자 검증
    log(f"  γ(3) = {mpmath.nstr(gamma_factor(mpmath.mpf(3)), 12)}")
    log(f"  γ(0.5) = {mpmath.nstr(gamma_factor(mpmath.mpc(0.5, 0)), 12)}")

    # 개별 기여 분석
    s_dbg = mpmath.mpf(3)
    N_s2 = mpmath.power(N_COND_MP, s_dbg/2)
    gam_s = gamma_factor(s_dbg)
    L_s = dirichlet_series(s_dbg, cn)
    log(f"  N^(3/2) = {mpmath.nstr(N_s2, 12)}")
    log(f"  γ(3) = {mpmath.nstr(gam_s, 12)}")
    log(f"  L(3) = {mpmath.nstr(L_s, 12)}")
    log(f"  Λ_direct(3) = {mpmath.nstr(N_s2 * gam_s * L_s, 12)}")

    # 단일 Φ 검증
    log(f"\n  단일 Φ(3,1) 검증:")
    phi_3_1 = mpmath.mpc(0)
    c = CONTOUR_C
    pf = mpmath.exp(c**2) / (2 * mpmath.pi)
    for k in range(N_HERMITE):
        v_k = HERM_NODES[k]
        wt_k = HERM_WEIGHTS[k]
        w = c + mpmath.mpc(0, v_k)
        sw = s_dbg + w
        val = (mpmath.power(N_COND_MP, sw/2) * gamma_factor(sw) *
               mpmath.power(mpmath.mpf(1), -sw) *
               mpmath.exp(2*mpmath.mpc(0,1)*c*v_k) / w)
        phi_3_1 += wt_k * val
    phi_3_1 *= pf
    log(f"  Φ(3,1) = {mpmath.nstr(phi_3_1, 12)}")
    log(f"  (기대: N^(3/2)·γ(3) = {mpmath.nstr(N_s2*gam_s, 12)})")

flush_to_file()

# ── 3. 함수방정식 직접 검증: Λ(s) ≈ Λ(1-s) ──
log(f"\n[Step 3] 함수방정식 Λ(s) = Λ(1-s) 직접 검증")
fe_direct_ok = True
for sp in [mpmath.mpc(0.5, 5), mpmath.mpc(0.5, 10), mpmath.mpc(0.5, 15),
           mpmath.mpc(0.3, 8), mpmath.mpc(0.7, 12)]:
    t_s = time.time()
    L_s = Lambda_AFE(sp, cn)
    L_1s = Lambda_AFE(1 - sp, cn)
    dt_s = time.time() - t_s

    if abs(L_s) > mpmath.mpf(10)**(-DPS+15):
        rel = float(abs(L_s - EPSILON * L_1s) / abs(L_s))
    else:
        rel = 0.0  # 둘 다 ~0

    ok = rel < 1e-6
    if not ok: fe_direct_ok = False
    log(f"  s={mpmath.nstr(sp,6)}: |Λ(s)|={float(abs(L_s)):.6e}, "
        f"|Λ(1-s)|={float(abs(L_1s)):.6e}, rel={rel:.2e} "
        f"{'✅' if ok else '❌'} ({dt_s:.1f}s)")
    flush_to_file()

if fe_direct_ok:
    log(f"  ✅ 함수방정식 검증 통과")
else:
    log(f"  ⚠️ 함수방정식 오차 존재 (DPS/N_coeff 부족 가능)")
flush_to_file()

# ── 4. Λ(1/2+it) 실수성 확인 + 영점 스캔 ──
log(f"\n[Step 4] Λ(1/2+it) 스캔 (t ∈ [{T_MIN}, {T_MAX}])")
log(f"  ε=+1 → Λ(1/2+it) ∈ ℝ")

scan_ts = np.arange(T_MIN, T_MAX + 0.5, 0.5)  # finer grid for zero detection
scan_vals = []
max_im = 0
t4_start = time.time()

for i, t in enumerate(scan_ts):
    val = Lambda_AFE(mpmath.mpc(0.5, t), cn)
    re_val = float(mpmath.re(val))
    im_val = float(mpmath.im(val))
    scan_vals.append(re_val)
    max_im = max(max_im, abs(im_val))

    if (i+1) % 10 == 0:
        elapsed = time.time() - t4_start
        eta = elapsed / (i+1) * (len(scan_ts) - i - 1)
        log(f"  [{i+1}/{len(scan_ts)}] t={t:.1f}: Re(Λ)={re_val:+.6e}, "
            f"Im={im_val:+.2e} ({elapsed:.0f}s, ETA {eta:.0f}s)")
        flush_to_file()

scan_vals_arr = np.array(scan_vals)
n_sign = int(np.sum(np.diff(np.sign(scan_vals_arr)) != 0))
dt_scan = time.time() - t4_start

log(f"\n  스캔 완료: {dt_scan:.1f}s")
log(f"  부호변환: {n_sign}")
log(f"  Re(Λ) 범위: [{scan_vals_arr.min():.6e}, {scan_vals_arr.max():.6e}]")
log(f"  최대 |Im|: {max_im:.2e} (0이면 실수성 확인)")
flush_to_file()

# ── 5. 영점 이분법 정밀화 ──
zero_approx = []
log(f"\n[Step 5] 영점 이분법 정밀화")

for i in range(len(scan_vals_arr) - 1):
    if scan_vals_arr[i] * scan_vals_arr[i+1] < 0:
        t_lo, t_hi = float(scan_ts[i]), float(scan_ts[i+1])
        val_lo = scan_vals_arr[i]
        try:
            for _ in range(25):
                t_mid = (t_lo + t_hi) / 2
                val_mid = float(mpmath.re(Lambda_AFE(mpmath.mpc(0.5, t_mid), cn)))
                if not np.isfinite(val_mid):
                    break
                if val_mid * val_lo < 0:
                    t_hi = t_mid
                else:
                    t_lo = t_mid; val_lo = val_mid
            z_t = (t_lo + t_hi) / 2
            zero_approx.append(z_t)
        except Exception as e:
            log(f"  WARNING bisection at t≈{(t_lo+t_hi)/2:.2f}: {e}")
            zero_approx.append((t_lo + t_hi) / 2)  # rough estimate

log(f"  발견 영점: {len(zero_approx)}개")
flush_to_file()

# LMFDB 매칭
matched = 0
for z in zero_approx:
    dists = [abs(z - lz) for lz in LMFDB_ZEROS]
    best = min(dists)
    best_idx = dists.index(best)
    ok = best < 0.1
    if ok: matched += 1
    log(f"  γ ≈ {z:.8f}  LMFDB[{best_idx}]={LMFDB_ZEROS[best_idx]:.6f}  "
        f"Δ={best:.6f} {'✅' if ok else '❌'}")

if len(zero_approx) > 0:
    log(f"  LMFDB 매칭: {matched}/{len(zero_approx)}")
else:
    log(f"  ⚠️ 영점 미발견!")
    # 전체 스캔 값 출력
    for i, t in enumerate(scan_ts):
        if i % 5 == 0:
            log(f"    t={t:.1f}: {scan_vals[i]:+.8e}")
flush_to_file()

# ── 6. κ 곡률 + 모노드로미 ──
if len(zero_approx) > 0:
    log(f"\n[Step 6] 영점 근방 κ + 모노드로미")

    kappa_near_list = []
    kappa_far_list = []
    mono_results = []

    # TP: 영점 근방
    for j, z_t in enumerate(zero_approx[:12]):
        s_pt = complex(0.5 + DELTA, z_t)
        k_val = curvature(s_pt, cn)
        m_val = monodromy(z_t, cn)
        kappa_near_list.append(k_val)
        mono_results.append(m_val)
        m_str = f"{m_val:.4f}" if m_val is not None else "FAIL"
        log(f"  TP #{j+1} γ={z_t:.6f}: κ={k_val:.2f}, mono/π={m_str}")
        flush_to_file()

    # FP: 영점 사이 (mid-points)
    log(f"\n  FP 측정 (영점 중간점):")
    fp_ts = []
    for j in range(len(zero_approx) - 1):
        mid = (zero_approx[j] + zero_approx[j+1]) / 2
        fp_ts.append(mid)
    # 첫 영점 이전, 마지막 이후도 추가
    if zero_approx[0] > T_MIN + 1:
        fp_ts.insert(0, T_MIN + 0.5)
    if zero_approx[-1] < T_MAX - 1:
        fp_ts.append(T_MAX - 1)

    for j, t in enumerate(fp_ts[:10]):
        k_val = curvature(complex(0.5 + DELTA, t), cn)
        kappa_far_list.append(k_val)
        log(f"  FP #{j+1} t={t:.3f}: κ={k_val:.2f}")
        flush_to_file()

    # κ 비율
    if kappa_near_list and kappa_far_list:
        near_med = float(np.median([k for k in kappa_near_list if np.isfinite(k) and k < 1e15]))
        far_med = float(np.median([k for k in kappa_far_list if np.isfinite(k) and k < 1e15]))
        if far_med > 0:
            ratio = near_med / far_med
        else:
            ratio = float('inf')
        log(f"\n  κ near(중앙값) = {near_med:.2f} ({len(kappa_near_list)}점)")
        log(f"  κ far(중앙값)  = {far_med:.4f} ({len(kappa_far_list)}점)")
        log(f"  비율 = {ratio:.1f}×")
    else:
        ratio = 0
        log(f"  데이터 부족")
    flush_to_file()

    # ── 7. σ-유일성 ──
    log(f"\n[Step 7] σ-유일성 스윕")
    sigma_results = {}
    for sig in [0.3, 0.5, 0.7, 0.9]:
        sc_vals = []
        for t in np.arange(T_MIN, min(T_MAX, 35.0), 0.8):
            val = Lambda_AFE(mpmath.mpc(sig, t), cn)
            sc_vals.append(float(mpmath.re(val)))
        n_sc = int(np.sum(np.diff(np.sign(sc_vals)) != 0))
        sigma_results[sig] = n_sc
        log(f"  σ={sig:.1f}: 부호변환={n_sc}")
        flush_to_file()

    sc_05 = sigma_results.get(0.5, 0)
    sc_others = [sigma_results.get(s, 0) for s in [0.3, 0.7, 0.9]]
    sigma_pass = sc_05 > 0 and sc_05 >= max(sc_others)

    if sigma_pass:
        log(f"  → σ=0.5 최대 ({sc_05}) ✅")
    else:
        max_sig = max(sigma_results, key=sigma_results.get)
        log(f"  → 극대 σ={max_sig} ({sigma_results[max_sig]})")

    # FAIL 패턴 분석 (수학자 지시)
    log(f"\n[Step 7b] σ-유일성 FAIL 패턴 분석")
    log(f"  GL(1) ζ:         PASS (σ=0.5 극대)")
    log(f"  GL(2) 11a1:      FAIL (σ=0.7-0.9 극대)")
    log(f"  GL(3) sym²11a1:  σ별 = {sigma_results}")
    if not sigma_pass:
        log(f"  → GL(3) FAIL 패턴: {'GL(2)와 동일' if max_sig in [0.7, 0.9] else '신규 패턴'}")
    else:
        log(f"  → GL(3) PASS: GL(2)만의 특수 현상 가능성")

    # σ-유일성 정량화 (진폭 비교)
    log(f"\n  σ별 진폭 분석:")
    for sig in [0.3, 0.5, 0.7, 0.9]:
        vals = []
        for t in np.arange(5.0, 20.0, 0.5):
            v = Lambda_AFE(mpmath.mpc(sig, t), cn)
            vals.append(float(abs(v)))
        if vals:
            log(f"    σ={sig:.1f}: |Λ| 평균={np.mean(vals):.4e}, 최대={np.max(vals):.4e}")
    flush_to_file()

    # ── 8. 공명 지표 ──
    log(f"\n[Step 8] 공명 지표 R(n) 분석 (수학자 지시)")
    log(f"  GL(1): Γ_ℝ(s) 1개, f_gamma ≈ log(t)/4π")
    log(f"  GL(2): Γ_ℂ(s) = Γ_ℝ(s)Γ_ℝ(s+1), f_gamma ≈ log(t)/2π")
    log(f"  GL(3): Γ_ℝ(s+1)²Γ_ℝ(s+2), f_gamma ≈ 3log(t)/4π")
    log(f"  영점 밀도: N(T) ~ (d/2π)T log T for degree d")

    # Stirling: |Γ(σ+it)| ~ √(2π) |t|^{σ-1/2} e^{-π|t|/2}
    # Phase of Γ(s/2) ≈ (t/2)log(t/2) - t/2 - π/4  (for large t)
    # f_gamma = d(phase)/dt for each gamma factor
    # GL(1): 1 factor → f ~ (1/2)log(t/(2π))/(2π) ~ log(t)/(4π)
    # GL(3) [0,0,1]: Γ(s/2)²·Γ((s+1)/2)
    #   phase ≈ 2·[(t/2)log(t/2)-t/2] + [(t/2)log(t/2)-t/2] = 3·[(t/2)log(t/2)-t/2]
    #   f ~ 3·log(t/2)/(2·2π) = 3log(t)/(4π)

    t_ref = 10.0
    f_zero = np.log(t_ref) / (2 * np.pi)  # 영점 밀도 주파수
    R_gl1 = np.log(t_ref) / (4 * np.pi) / f_zero
    R_gl2 = np.log(t_ref) / (2 * np.pi) / f_zero
    R_gl3 = 3 * np.log(t_ref) / (4 * np.pi) / f_zero
    log(f"  t={t_ref}에서 공명 지표:")
    log(f"    R(GL(1)) = {R_gl1:.3f}")
    log(f"    R(GL(2)) = {R_gl2:.3f}")
    log(f"    R(GL(3)) = {R_gl3:.3f}")
    log(f"  R≈0.5: GL(1) → σ-유일성 PASS")
    log(f"  R≈1.0: GL(2) → σ-유일성 FAIL")
    log(f"  R≈1.5: GL(3) → σ-유일성 {'확인 필요' if not sigma_pass else 'PASS'}")
    flush_to_file()

else:
    ratio = 0
    sigma_pass = False
    mono_results = []

# ━━━━━━━━━━━ 종합 판정 ━━━━━━━━━━━
log(f"\n{'='*70}")
log(f"종합 판정 — GL(3) sym²(11a1) 4성질")
log(f"{'='*70}")

# 1. 함수방정식
fe_str = "PASS ✅" if (fe_validated and fe_direct_ok) else "FAIL ❌"
log(f"  [필수] 함수방정식 Λ(s)=Λ(1-s): {fe_str}")

# 2. 영점
log(f"  영점 발견: {len(zero_approx)}개 (LMFDB 매칭: {matched})")

# 3. σ-유일성
if sigma_pass:
    su_str = "PASS ✅"
elif len(zero_approx) > 0:
    su_str = "FAIL (패턴 분석 포함)"
else:
    su_str = "N/A (영점 미발견)"
log(f"  [필수] σ-유일성: {su_str}")

# 4. κ 집중
if len(zero_approx) > 0:
    if ratio >= 10:
        kr_str = f"PASS ✅ ({ratio:.1f}×)"
    elif ratio > 1:
        kr_str = f"약한 양성 ({ratio:.1f}×)"
    else:
        kr_str = f"FAIL ({ratio:.1f}×)"
else:
    kr_str = "N/A"
log(f"  [양성] κ 집중: {kr_str}")

# 5. 모노드로미
if mono_results:
    tp_mono = sum(1 for m in mono_results if m is not None and m > 1.5)
    log(f"  [양성] 모노드로미: TP {tp_mono}/{len(mono_results)} (mono/π > 1.5)")
else:
    log(f"  [양성] 모노드로미: N/A")

total_time = time.time() - t_total
log(f"\n  총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log(f"  시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")

flush_to_file()
log(f"\n결과 저장: {OUTFILE}")
