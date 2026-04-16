"""
Quick test: which gamma shifts give zeros matching LMFDB for sym²(11a1)?
Tests μ = [0,0,1], [0,1,1], [0,1,2], [1,1,2]
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))
import mpmath

DPS = 50
mpmath.mp.dps = DPS

N_COND = mpmath.mpf(121)
EPSILON = 1
N_HERMITE = 50
CONTOUR_C = mpmath.mpf(2)

_hn, _hw = np.polynomial.hermite.hermgauss(N_HERMITE)
HN = [mpmath.mpf(float(x)) for x in _hn]
HW = [mpmath.mpf(float(w)) for w in _hw]

LMFDB_ZEROS = [3.899281, 4.734595, 6.189477, 7.312039, 8.650148,
               10.128219, 10.936767, 12.070610, 12.744576, 13.335280]

# ── 계수 ──
def compute_cn(limit):
    known_ap = {2:-2, 3:-1, 5:1, 7:-2, 11:1, 13:4, 17:-2, 19:0, 23:-1,
                29:0, 31:7, 37:3, 41:-8, 43:-6, 47:8, 53:-6, 59:5, 61:12,
                67:-7, 71:-4, 73:2, 79:10, 83:-7, 89:6, 97:-4, 101:-2, 103:12, 107:10, 109:-4}
    sieve = [True]*(limit+1); sieve[0]=sieve[1]=False
    for i in range(2, int(limit**0.5)+1):
        if sieve[i]:
            for j in range(i*i, limit+1, i): sieve[j] = False
    primes = [i for i in range(2, limit+1) if sieve[i]]
    ap = {}
    for p in primes:
        if p in known_ap: ap[p] = known_ap[p]; continue
        count = 0
        for x in range(p):
            rhs = (x*x*x - x*x - 10*x - 20) % p
            if rhs == 0: count += 1
            elif pow(rhs, (p-1)//2, p) == 1: count += 2
        ap[p] = p - count

    # sym² coefficients (analytic normalization)
    cpk = {}
    for p in primes:
        if p > limit: break
        cpk[(p,0)] = mpmath.mpf(1)
        if p == 11:
            for k in range(1, 20):
                cpk[(p,k)] = mpmath.power(mpmath.mpf(11), -k)
                if 11**(k+1) > limit: break
        else:
            c1 = mpmath.mpf(ap[p])**2 / mpmath.mpf(p) - 1
            cpk[(p,1)] = c1
            pk = p; k = 1
            while pk*p <= limit:
                pk *= p; k += 1
                cpk[(p,k)] = c1*cpk[(p,k-1)] - c1*cpk[(p,k-2)] + cpk.get((p,k-3), mpmath.mpf(0))

    cn = [mpmath.mpf(0)]*(limit+1); cn[1] = mpmath.mpf(1)
    for n in range(2, limit+1):
        temp = n; result = mpmath.mpf(1)
        for p in primes:
            if p*p > temp: break
            if temp % p == 0:
                k = 0
                while temp % p == 0: k += 1; temp //= p
                if (p,k) in cpk: result *= cpk[(p,k)]
                else: result = mpmath.mpf(0); break
        if temp > 1:
            if (temp,1) in cpk: result *= cpk[(temp,1)]
            else: result = mpmath.mpf(0)
        cn[n] = result
    return cn

cn = compute_cn(150)
print(f"c(2)={float(cn[2]):.6f}, c(3)={float(cn[3]):.6f}, c(5)={float(cn[5]):.6f}", flush=True)

# ── 감마 인자 옵션 ──
def make_gamma(shifts):
    """shifts = list of μ_j values"""
    def G(s):
        result = mpmath.mpf(1)
        for mu in shifts:
            # Γ_ℝ(s+μ) = π^{-(s+μ)/2} Γ((s+μ)/2)
            result *= mpmath.power(mpmath.pi, -(s + mu) / 2) * mpmath.gamma((s + mu) / 2)
        return result
    return G

gamma_options = {
    "[0,0,1]": make_gamma([0, 0, 1]),
    "[0,1,1]": make_gamma([0, 1, 1]),
    "[0,1,2]": make_gamma([0, 1, 2]),
    "[1,1,2]": make_gamma([1, 1, 2]),
}

# ── AFE ──
def Lambda_AFE(s, cn, G, c=CONTOUR_C):
    s_mp = mpmath.mpc(s)
    s1_mp = 1 - s_mp
    pf = mpmath.exp(c**2) / (2 * mpmath.pi)

    total = mpmath.mpc(0)
    for k in range(N_HERMITE):
        v_k = HN[k]; wt_k = HW[k]
        iv_k = mpmath.mpc(0, v_k)
        w_shift = c + iv_k
        exp_phase = mpmath.exp(2 * mpmath.mpc(0, 1) * c * v_k)

        # s part
        sw = s_mp + w_shift
        A_s = mpmath.power(N_COND, sw/2) * G(sw) * exp_phase / w_shift
        D_s = sum(cn[n] * mpmath.power(n, -sw) for n in range(1, len(cn)) if cn[n] != 0)

        # 1-s part
        s1w = s1_mp + w_shift
        A_1s = mpmath.power(N_COND, s1w/2) * G(s1w) * exp_phase / w_shift
        D_1s = sum(cn[n] * mpmath.power(n, -s1w) for n in range(1, len(cn)) if cn[n] != 0)

        total += wt_k * (A_s * D_s + EPSILON * A_1s * D_1s)

    return pf * total

# ── 각 gamma로 LMFDB 첫 영점 근방 스캔 ──
print("\n=== 감마 인자 비교: LMFDB 첫 3영점 근방 Λ(1/2+it) ===", flush=True)

for name, G in gamma_options.items():
    print(f"\n--- {name} ---", flush=True)
    t0 = time.time()

    # 첫 영점 γ₁ ≈ 3.899 근방 스캔
    for t in [3.0, 3.5, 3.9, 4.0, 4.5, 5.0, 5.5, 6.0, 6.2, 6.5, 7.0, 7.3, 7.5, 8.0]:
        val = Lambda_AFE(mpmath.mpc(0.5, t), cn, G)
        re = float(mpmath.re(val))
        print(f"  t={t:5.1f}: Re={re:+.6e}", flush=True)

    dt = time.time() - t0
    print(f"  ({dt:.1f}s)", flush=True)

    # FE check
    s_test = mpmath.mpc(0.5, 10)
    L_s = Lambda_AFE(s_test, cn, G)
    L_1s = Lambda_AFE(1 - s_test, cn, G)
    if abs(L_s) > 1e-40:
        fe_err = float(abs(L_s - L_1s) / abs(L_s))
    else:
        fe_err = 0
    print(f"  FE check: |Λ(s)-Λ(1-s)|/|Λ(s)| = {fe_err:.2e}", flush=True)

print("\n=== 완료 ===", flush=True)
