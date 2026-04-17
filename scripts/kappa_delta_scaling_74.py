#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #74 — κ-δ 스케일링 GL(2)+GL(3) degree 비교
=============================================================================
목적:
  κ(s₀+δ) = 1/δ² + A(t₀) + O(δ²) 정리를 GL(2) Maass + GL(3) sym²(11a1)에서
  수치 확인하고, A(t₀)의 degree 의존성 분석.

대상:
  GL(1): ζ(s)  — 결과 #73 hardcode
  GL(2): EVEN Maass form (R=13.78) — #69b 동일 L-함수, scipy/numpy (빠름)
  GL(3): sym²(11a1) — #63 동일 L-함수 (mpmath DPS=50, γ=[1,1,2])

핵심 설계:
  - GL(2): scipy/numpy, 영점 #69b에서 scipy.minimize로 정밀화
  - GL(3): mpmath DPS=50, #63 영점 hardcode, #63 동일 감마인자 사용
  - A(t₀) = κ - 1/δ²  at δ=0.01

주의:
  - κ·δ²≈1은 수학적 귀결 (자명) — 검토자 경고 반영
  - 참신성은 A(t₀)의 degree별 값

결과: results/kappa_delta_scaling_74.txt
=============================================================================
"""

import sys, os, time, re
import numpy as np
from scipy.special import gamma as sc_gamma
from scipy.optimize import minimize_scalar
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

DPS_GL3 = 50
mpmath.mp.dps = DPS_GL3

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "kappa_delta_scaling_74.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))
def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

DELTAS = [0.01, 0.02, 0.03, 0.05, 0.10]

log("=" * 72)
log("결과 #74 — κ-δ 스케일링 GL(2)+GL(3) degree 비교")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"δ 목록: {DELTAS}")
log()
log("주의: κ·δ²≈1은 κ=1/δ²+A(t₀)+O(δ²) 정리의 수학적 귀결 (자명).")
log("참신성은 A(t₀)의 degree별 값과 구조.")
log()

t_total = time.time()

# ════════════════════════════════════════════════════════
#  SECTION 0: GL(1) 참조
# ════════════════════════════════════════════════════════

log("=" * 72)
log("[GL(1) 참조] ζ(s) — 결과 #73 (13영점, δ=0.01)")
log("=" * 72)

gl1_kappa = np.array([10000.47, 10000.85, 10000.67, 10001.47, 10000.81,
                       10001.21, 10001.48, 10000.90, 10002.46, 10001.28,
                       10001.10, 10001.53, 10002.73])
gl1_A = gl1_kappa - 10000.0
log(f"  n(zeros)={len(gl1_A)}, mean(A)={float(np.mean(gl1_A)):.4f}, range=[{float(np.min(gl1_A)):.3f},{float(np.max(gl1_A)):.3f}]")
log()
flush_file()

# ════════════════════════════════════════════════════════
#  SECTION 1: GL(2) Maass — scipy/numpy
# ════════════════════════════════════════════════════════

log("=" * 72)
log("[GL(2)] EVEN Maass (R≈13.78, SL(2,Z)) — scipy/numpy + zero refinement")
log("=" * 72)
log()

R_GL2   = 13.7797513518907389
N_GH2   = 50
C_GL2   = 2.0
N_MAX2  = 200

GH_X2, GH_W2 = np.polynomial.hermite.hermgauss(N_GH2)

# #69b 영점 (초기값)
GL2_ZEROS_INIT = np.array([
    2.8977232873, 5.5912455142, 21.0903738320, 23.1575048745,
    25.4393050253, 29.1892052114, 31.0617492139, 32.4527095973,
    34.0272720516, 36.9312471688, 38.9871129692, 40.4655165613,
])

def load_maass_cn():
    cf = os.path.join(os.path.dirname(os.path.abspath(__file__)), "maass_coeff13_raw.txt")
    if not os.path.exists(cf): return None, None
    with open(cf, 'r') as f: text = f.read()
    ms = re.findall(r'C\[(\d+)\]=\s*([-+]?\d+\.\d+)', text)
    cd = {int(n_): float(v_) for n_, v_ in ms if 1 <= int(n_) <= N_MAX2}
    if not cd: return None, None
    ns = sorted(cd); return np.array(ns, dtype=float), np.array([cd[n] for n in ns])

log("[Step 1a] GL(2) Maass 계수 로드")
n_gl2, cn_gl2 = load_maass_cn()
gl2_ok = n_gl2 is not None
if gl2_ok:
    log(f"  계수: {len(n_gl2)}개")
    log_n_gl2 = np.log(n_gl2)
else:
    log("  ❌ 로드 실패")

def Lambda_gl2(s):
    """Λ(s) via AFE (scipy/numpy, N=200, N_GH=50)"""
    s = complex(s); s1 = 1.0 - s
    iR = 1j * R_GL2
    ec2 = np.exp(C_GL2**2)
    total = 0.0j
    for k in range(N_GH2):
        y = GH_X2[k]; w = GH_W2[k]
        wk = C_GL2 + 1j*y
        ph = np.exp(2j*C_GL2*y) / wk
        for sv, eps_ in [(s+wk, 1), (s1+wk, 1)]:
            gf = np.pi**(-sv) * sc_gamma((sv+iR)/2) * sc_gamma((sv-iR)/2)
            Ls = np.dot(cn_gl2, np.exp(-sv * log_n_gl2))
            total += w * gf * Ls * ph * eps_
    return ec2 / (2*np.pi) * total

def kappa_gl2(t, delta):
    h = min(delta/3.0, 1e-5)
    s = complex(0.5+delta, t)
    try:
        L0 = Lambda_gl2(s)
        if abs(L0) < 1e-250: return float('inf')
        Lp = Lambda_gl2(s+h); Lm = Lambda_gl2(s-h)
        k = float(abs((Lp-Lm)/(2*h*L0))**2)
        return k if np.isfinite(k) else float('inf')
    except Exception as e:
        log(f"    WARNING GL2 kappa: {e}")
        return float('inf')

# 영점 정밀화
log("[Step 1b] GL(2) 영점 정밀화 (scipy.minimize on |Λ(0.5+it)|²)")
if not gl2_ok:
    gl2_zeros = []
else:
    gl2_zeros = []
    for j, t0 in enumerate(GL2_ZEROS_INIT):
        try:
            r = minimize_scalar(lambda t: abs(Lambda_gl2(complex(0.5, t)))**2,
                                bounds=(t0-0.3, t0+0.3), method='bounded',
                                options={'xatol':1e-8, 'maxiter':100})
            t_r = r.x; v = r.fun
            log(f"  #{j+1:>2} {t0:.5f} → {t_r:.8f}  |Λ|²={v:.2e}  {'✅' if v<1e-2 else '⚠️'}")
            if v < 1e2:  # 상당히 작으면 Accept
                gl2_zeros.append(t_r)
        except Exception as e:
            log(f"  #{j+1:>2} 실패: {e} → 초기값 사용")
            gl2_zeros.append(t0)
    log(f"  최종 GL(2) 영점: {len(gl2_zeros)}개")
    if not gl2_zeros:
        gl2_ok = False

flush_file()

log()
log("[Step 1c] GL(2) κ-δ 스윕")
gl2_dk = {}  # delta → kappas
gl2_data = []  # (delta, t, kappa)

if gl2_ok and len(gl2_zeros) >= 2:
    for delta in DELTAS:
        log(f"  ── δ={delta} ──")
        kappas = []
        for j, z in enumerate(gl2_zeros):
            k = kappa_gl2(z, delta)
            kappas.append(k); gl2_data.append((delta, z, k))
            if np.isfinite(k) and k < 1e15:
                A_v = k - 1.0/delta**2
                log(f"    #{j+1:>2} t={z:.6f}: κ={k:>12.2f}  κ·δ²={k*delta**2:.5f}  A={A_v:.4f}")
        gl2_dk[delta] = kappas
        fin = [k for k in kappas if np.isfinite(k) and 0<k<1e15]
        if fin: log(f"  → median={np.median(fin):.2f}, κ·δ²={np.median(fin)*delta**2:.5f}")
        flush_file()

# ════════════════════════════════════════════════════════
#  SECTION 2: GL(3) sym²(11a1) — mpmath DPS=50 (#63과 동일)
# ════════════════════════════════════════════════════════

log()
log("=" * 72)
log("[GL(3)] sym²(11a1) — mpmath DPS=50, γ=[1,1,2] (#63과 동일)")
log("=" * 72)
log()

# #63 파라미터
N_COND3     = mpmath.mpf(121)
EPS_GL3     = 1
N_MAX3      = 120
N_HERMITE3  = 50
CONTOUR_C3  = mpmath.mpf(2)

_hx, _hw = np.polynomial.hermite.hermgauss(N_HERMITE3)
HERM_X3 = [mpmath.mpf(float(x)) for x in _hx]
HERM_W3 = [mpmath.mpf(float(w)) for w in _hw]

# #63 영점 (LMFDB 매칭, mpmath 정밀화 결과)
GL3_ZEROS_63 = [
    3.87646463, 4.72367882, 6.15862826, 7.30312439,
    8.65033064, 10.12916736, 11.01310887, 12.10705922,
    12.72067789, 13.33562852, 14.39334338, 15.32994895,
]

# #63 계수 (mpmath)
def compute_11a1_an_mp(limit):
    known = {2:-2,3:-1,5:1,7:-2,11:1,13:4,17:-2,19:0,23:-1,
             29:0,31:7,37:3,41:-8,43:-6,47:8,53:-6,59:5,61:12,
             67:-7,71:-4,73:2,79:10,83:-7,89:6,97:-4,101:-2,
             103:12,107:10,109:-4}
    sieve = [True]*(limit+1); sieve[0]=sieve[1]=False
    for i in range(2, int(limit**0.5)+1):
        if sieve[i]:
            for j in range(i*i, limit+1, i): sieve[j]=False
    primes = [i for i in range(2, limit+1) if sieve[i]]
    ap = {}
    for p in primes:
        if p in known: ap[p]=known[p]; continue
        cnt=0
        for x in range(p):
            rhs=(x*x*x-x*x-10*x-20)%p
            if rhs==0: cnt+=1
            elif pow(rhs,(p-1)//2,p)==1: cnt+=2
        ap[p]=p-cnt
    apk={}
    for p in primes:
        apk[(p,0)]=1; apk[(p,1)]=ap[p]
        pk=p; kk=1
        while pk*p<=limit:
            pk*=p; kk+=1
            if p==11: apk[(p,kk)]=ap[p]**kk
            else: apk[(p,kk)]=ap[p]*apk[(p,kk-1)]-p*apk[(p,kk-2)]
    an=[0]*(limit+1); an[1]=1
    for n in range(2,limit+1):
        temp=n; result=1
        for p in primes:
            if p*p>temp: break
            if temp%p==0:
                kk=0
                while temp%p==0: kk+=1; temp//=p
                result*=apk.get((p,kk),0)
        if temp>1: result*=ap.get(temp,0)
        an[n]=result
    return an, primes

def compute_sym2_cn_mp(an, primes, limit):
    mpmath.mp.dps = DPS_GL3+10
    cpk={}
    for p in primes:
        if p>limit: break
        cpk[(p,0)]=mpmath.mpf(1)
        if p==11:
            for kk in range(1,20):
                cpk[(p,kk)]=mpmath.power(mpmath.mpf(11),-kk)
                if 11**(kk+1)>limit: break
        else:
            c1=mpmath.mpf(an[p])**2/mpmath.mpf(p)-1
            cpk[(p,1)]=c1
            pk=p; kk=1
            while pk*p<=limit:
                pk*=p; kk+=1
                b1=cpk[(p,kk-1)]; b2=cpk[(p,kk-2)]
                b3=cpk.get((p,kk-3),mpmath.mpf(0))
                cpk[(p,kk)]=c1*b1-c1*b2+b3
    cn=[mpmath.mpf(0)]*(limit+1); cn[1]=mpmath.mpf(1)
    for n in range(2,limit+1):
        temp=n; result=mpmath.mpf(1)
        for p in primes:
            if p*p>temp: break
            if temp%p==0:
                kk=0
                while temp%p==0: kk+=1; temp//=p
                if (p,kk) in cpk: result*=cpk[(p,kk)]
                else: result=mpmath.mpf(0); break
        if temp>1:
            if (temp,1) in cpk: result*=cpk[(temp,1)]
            else: result=mpmath.mpf(0)
        cn[n]=result
    mpmath.mp.dps=DPS_GL3
    return cn

def gamma_factor_gl3_mp(s):
    """#63과 동일: γ(s) = π^{-(3s+4)/2} · Γ((s+1)/2)² · Γ((s+2)/2)"""
    return (mpmath.power(mpmath.pi, -(3*s+4)/2)
            * mpmath.gamma((s+1)/2)**2
            * mpmath.gamma((s+2)/2))

def dirichlet_gl3_mp(w, cn_gl3):
    total = mpmath.mpc(0)
    for n in range(1, len(cn_gl3)):
        if cn_gl3[n]==0: continue
        total += cn_gl3[n] * mpmath.power(n, -w)
    return total

def Lambda_AFE_gl3_mp(s, cn_gl3):
    """#63 Lambda_AFE 그대로"""
    c = CONTOUR_C3
    s_mp = mpmath.mpc(s)
    s1_mp = 1 - s_mp
    prefactor = mpmath.exp(c**2) / (2*mpmath.pi)
    total = mpmath.mpc(0)
    for k in range(N_HERMITE3):
        vk = HERM_X3[k]; wt = HERM_W3[k]
        ivk = mpmath.mpc(0, vk)
        wshift = c + ivk
        exp_ph = mpmath.exp(2*mpmath.mpc(0,1)*c*vk)
        # s-term
        sw = s_mp + wshift
        A_s = mpmath.power(N_COND3, sw/2) * gamma_factor_gl3_mp(sw) * exp_ph / wshift
        D_s = dirichlet_gl3_mp(sw, cn_gl3)
        # 1-s term
        s1w = s1_mp + wshift
        A_1s = mpmath.power(N_COND3, s1w/2) * gamma_factor_gl3_mp(s1w) * exp_ph / wshift
        D_1s = dirichlet_gl3_mp(s1w, cn_gl3)
        total += wt * (A_s*D_s + EPS_GL3*A_1s*D_1s)
    return prefactor * total

def kappa_gl3_mp(t, delta, cn_gl3):
    """κ = |Λ'/Λ|² at s=0.5+δ+it, h=min(δ/3,1e-5) — mpmath"""
    h_f = min(delta/3.0, 1e-5)
    h = mpmath.mpf(str(h_f))
    s = mpmath.mpc(mpmath.mpf('0.5')+mpmath.mpf(str(delta)), mpmath.mpf(str(t)))
    try:
        L0 = Lambda_AFE_gl3_mp(s, cn_gl3)
        if abs(L0) < mpmath.mpf(10)**(-DPS_GL3+15):
            log(f"    WARNING: Λ≈0 at t={t}, δ={delta}")
            return float('inf')
        Lp = Lambda_AFE_gl3_mp(s+h, cn_gl3)
        Lm = Lambda_AFE_gl3_mp(s-h, cn_gl3)
        conn = (Lp-Lm)/(2*h*L0)
        k = float(abs(conn)**2)
        return k if np.isfinite(k) else float('inf')
    except Exception as e:
        log(f"    WARNING GL3 kappa at t={t}, δ={delta}: {e}")
        return float('inf')

# ─── GL(3) 계수 준비 ───
log("[Step 2a] GL(3) sym²(11a1) 계수 준비 (mpmath, N=120)")
t_cn = time.time()
an_gl3, primes_gl3 = compute_11a1_an_mp(N_MAX3)
cn_gl3 = compute_sym2_cn_mp(an_gl3, primes_gl3, N_MAX3)
log(f"  계수 준비 완료: {time.time()-t_cn:.1f}s")

# ─── GL(3) 영점 정밀화 (#63 영점에서 mpmath minimize) ───
log("[Step 2b] GL(3) 영점 정밀화 (#63 영점 → mpmath |Λ|² minimize)")

def Lambda_norm_gl3_float(t):
    val = Lambda_AFE_gl3_mp(mpmath.mpc(0.5, t), cn_gl3)
    return float(abs(val)**2)

gl3_zeros = []
for j, t0 in enumerate(GL3_ZEROS_63):
    try:
        r = minimize_scalar(Lambda_norm_gl3_float,
                            bounds=(t0-0.2, t0+0.2), method='bounded',
                            options={'xatol':1e-6, 'maxiter':50})
        t_r = r.x; v = r.fun
        log(f"  #{j+1:>2} {t0:.8f} → {t_r:.10f}  |Λ|²={v:.2e}")
        gl3_zeros.append(t_r)
    except Exception as e:
        log(f"  #{j+1:>2} 실패: {e} → 초기값")
        gl3_zeros.append(t0)

gl3_ok = len(gl3_zeros) >= 8
log(f"  최종 GL(3) 영점: {len(gl3_zeros)}개")
if not gl3_ok:
    log("  ⚠️ 영점 <8개")

flush_file()

# ─── GL(3) δ 스윕 ───
log()
log("[Step 2c] GL(3) κ-δ 스윕 (mpmath)")
gl3_dk = {}
gl3_data = []

if gl3_ok and len(gl3_zeros) >= 2:
    for delta in DELTAS:
        log(f"  ── δ={delta} ──")
        t_d = time.time()
        kappas = []
        for j, z in enumerate(gl3_zeros):
            k = kappa_gl3_mp(z, delta, cn_gl3)
            kappas.append(k); gl3_data.append((delta, z, k))
            if np.isfinite(k) and k < 1e15:
                A_v = k - 1.0/delta**2
                log(f"    #{j+1:>2} t={z:.6f}: κ={k:>12.2f}  κ·δ²={k*delta**2:.5f}  A={A_v:.4f}")
            else:
                log(f"    #{j+1:>2} t={z:.6f}: κ=FAIL")
        gl3_dk[delta] = kappas
        fin = [k for k in kappas if np.isfinite(k) and 0<k<1e15]
        if fin: log(f"  → n={len(fin)}, median={np.median(fin):.2f}, κ·δ²={np.median(fin)*delta**2:.5f}  ({time.time()-t_d:.1f}s)")
        flush_file()

# ════════════════════════════════════════════════════════
#  SECTION 3: A(t₀) 추출
# ════════════════════════════════════════════════════════

log()
log("=" * 72)
log("[Step 3] A(t₀) 추출 (A = κ - 1/δ²  at δ=0.01)")
log("=" * 72)
log()
log("  이론: κ(s₀+δ) = 1/δ² + A(t₀) + O(δ²)")
log()

D_REF = 0.01; INV2 = 1.0/D_REF**2

gl2_A = []
if gl2_ok and D_REF in gl2_dk:
    log("  GL(2) A 목록:")
    for j,(k,z) in enumerate(zip(gl2_dk[D_REF], gl2_zeros)):
        if np.isfinite(k) and 0<k<1e15:
            A=k-INV2; gl2_A.append(A)
            log(f"    #{j+1} t={z:.6f}: κ={k:.4f}, A={A:.4f}")

gl3_A = []
if gl3_ok and D_REF in gl3_dk:
    log("  GL(3) A 목록:")
    for j,(k,z) in enumerate(zip(gl3_dk[D_REF], gl3_zeros)):
        if np.isfinite(k) and 0<k<1e15:
            A=k-INV2; gl3_A.append(A)
            log(f"    #{j+1} t={z:.6f}: κ={k:.4f}, A={A:.4f}")
log()
flush_file()

# ════════════════════════════════════════════════════════
#  SECTION 4: κ·δ² 스케일링 요약
# ════════════════════════════════════════════════════════

log("=" * 72)
log("[Step 4] κ-δ 스케일링 요약")
log("=" * 72)
log()

def summarize(name, dk):
    log(f"  ── {name} ──")
    log(f"  {'δ':>5} | {'1/δ²':>8} | {'median':>13} | {'CV%':>6} | {'κ·δ²':>8} | 판정")
    log(f"  {'-'*60}")
    prods = []
    for d in DELTAS:
        fin = [k for k in dk.get(d,[]) if np.isfinite(k) and 0<k<1e15]
        if not fin: log(f"  {d:>5.2f} | {'N/A':>8} | {'N/A':>13} | {'N/A':>6} | {'N/A':>8} | ❌"); continue
        med = float(np.median(fin))
        cv  = float(np.std(fin)/np.mean(fin)*100) if np.mean(fin)>0 else 999
        prod = med*d**2
        prods.append(prod)
        v = ("★★★" if 0.95<=prod<=1.15 else "★★" if 0.90<=prod<=1.20 else "★" if 0.80<=prod<=1.30 else "❌")
        log(f"  {d:>5.2f} | {1/d**2:>8.0f} | {med:>13.2f} | {cv:>6.2f} | {prod:>8.5f} | {v}")
    if prods:
        log(f"  κ·δ² ∈ [{min(prods):.5f}, {max(prods):.5f}]  all_pass={all(0.95<=p<=1.15 for p in prods)}")
    log()

if gl2_ok and gl2_dk: summarize("GL(2) Maass EVEN (R≈13.78)", gl2_dk)
if gl3_ok and gl3_dk: summarize("GL(3) sym²(11a1)", gl3_dk)
flush_file()

# ════════════════════════════════════════════════════════
#  SECTION 5: 3-degree 비교표
# ════════════════════════════════════════════════════════

log("=" * 72)
log("[Step 5] ★ 3-degree A(t₀) 비교표 ★")
log("=" * 72)
log()
log(f"  {'degree':>7} | {'L-함수':>22} | {'mean(A)':>9} | {'range(A)':>20} | {'n':>4}")
log(f"  {'-'*72}")

rows = [("GL(1)", "ζ(s)",
         float(np.mean(gl1_A)), float(np.min(gl1_A)), float(np.max(gl1_A)), len(gl1_A))]

if gl2_A:
    rows.append(("GL(2)", "Maass R≈13.78",
                 float(np.mean(gl2_A)), float(min(gl2_A)), float(max(gl2_A)), len(gl2_A)))
else:
    rows.append(("GL(2)", "Maass R≈13.78", float('nan'), float('nan'), float('nan'), 0))

if gl3_A:
    rows.append(("GL(3)", "sym²(11a1)",
                 float(np.mean(gl3_A)), float(min(gl3_A)), float(max(gl3_A)), len(gl3_A)))
else:
    rows.append(("GL(3)", "sym²(11a1)", float('nan'), float('nan'), float('nan'), 0))

for deg, lname, m, mn, mx, n in rows:
    if np.isfinite(m):
        log(f"  {deg:>7} | {lname:>22} | {m:>9.4f} | [{mn:.3f}, {mx:.3f}]{' ':>7} | {n:>4}")
    else:
        log(f"  {deg:>7} | {lname:>22} | {'N/A':>9} | {'N/A':>20} | {n:>4}")

log()
valid = [(d, m) for d, _, m, *_ in rows if np.isfinite(m)]
if len(valid) >= 2:
    means_ = [m for _, m in valid]
    is_inc = all(means_[i]<=means_[i+1] for i in range(len(means_)-1))
    is_dec = all(means_[i]>=means_[i+1] for i in range(len(means_)-1))
    mono = "단조증가" if is_inc else ("단조감소" if is_dec else "비단조")
    log(f"  mean(A): {' → '.join(f'{d}:{m:.4f}' for d,m in valid)}")
    log(f"  단조성: {mono}")
    if is_inc:   log("  → A가 degree와 함께 증가: 고차 Γ-인자 구조가 A를 키움")
    elif is_dec: log("  → A가 degree와 함께 감소: 고차 Γ-인자에서 상쇄 발생")
    else:        log("  → 비단조: A의 degree 의존성은 단순하지 않음")
log()
flush_file()

# ════════════════════════════════════════════════════════
#  SECTION 6: 성공 기준
# ════════════════════════════════════════════════════════

total_time = time.time() - t_total
log("=" * 72)
log("[Step 6] 성공 기준 체크")
log("=" * 72)
log()

gl2_n  = len(gl2_zeros) if gl2_ok else 0
gl3_n  = len(gl3_zeros) if gl3_ok else 0
gl2_p  = sum(1 for _,_,k in gl2_data if np.isfinite(k) and k<1e15)
gl3_p  = sum(1 for _,_,k in gl3_data if np.isfinite(k) and k<1e15)

def chk(dk):
    if not dk: return False
    for d in DELTAS:
        fin = [k for k in dk.get(d,[]) if np.isfinite(k) and 0<k<1e15]
        if not fin: return False
        med = float(np.median(fin))
        if not (0.95 <= med*d**2 <= 1.15): return False
    return True

ok = [gl2_n>=8, gl3_n>=8, gl2_p>=40, gl3_p>=40,
      chk(gl2_dk), chk(gl3_dk), len(valid)>=2, total_time<120]
lbl = ["GL(2) 영점≥8", "GL(3) 영점≥8", "GL(2) ≥40pts", "GL(3) ≥40pts",
       "κ·δ² GL(2)", "κ·δ² GL(3)", "3-deg표", "<120s"]
vals = [gl2_n, gl3_n, gl2_p, gl3_p, '', '', len(valid), f"{total_time:.1f}s"]

for v, l, val in zip(ok, lbl, vals):
    log(f"  {'✅' if v else '❌'} {l}: {val}")

n_pass = sum(ok)
verdict = ("★★★ 강양성" if n_pass==8 else "★★ 양성" if n_pass>=6 else
           "★ 조건부" if n_pass>=4 else "❌ 음성")
log()
log(f"  최종 판정: {verdict} ({n_pass}/8)")
log()

log("─"*72)
log("논문 반영 메모 (수학자용):")
log("  B-18: A(t₀) degree 의존성 → 결과 #61")
for deg, lname, m, mn, mx, n in rows:
    if np.isfinite(m):
        log(f"    {deg}: A∈[{mn:.3f},{mx:.3f}], mean={m:.4f} (n={n})")
log("─"*72)
log()
log(f"총 소요: {total_time:.1f}s")
log("[완료]")

flush_file()
print(f"\n결과 저장: {OUTFILE}", flush=True)
