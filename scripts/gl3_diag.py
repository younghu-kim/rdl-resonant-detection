"""
GL(3) sym²(11a1) 진단: AFE 불일치 원인 분석
"""
import mpmath, numpy as np, time
mpmath.mp.dps = 50

Q = mpmath.mpf(11)   # √121
PI = mpmath.pi

# ── 계수 (검증됨) ──
known_ap = {2:-2,3:-1,5:1,7:-2,11:1,13:4,17:-2,19:0,23:-1,
            29:0,31:7,37:3,41:-8,43:-6,47:8,53:-6,59:5,61:12,67:-7,71:-4,73:2,79:10}

def gamma_011(s):
    """μ=(0,1,1): Γ_ℝ(s)Γ_ℝ(s+1)²"""
    return (mpmath.power(PI, -(3*s+2)/2) *
            mpmath.gamma(s/2) * mpmath.gamma((s+1)/2)**2)

def gamma_012(s):
    """μ=(0,1,2): Γ_ℝ(s)Γ_ℝ(s+1)Γ_ℝ(s+2)"""
    return (mpmath.power(PI, -(3*s+3)/2) *
            mpmath.gamma(s/2) * mpmath.gamma((s+1)/2) * mpmath.gamma((s+2)/2))

# ── Test 1: Mellin transform of K₃ ──
# K₃(ξ) = 4∫ exp(-2πξ/v - πv²) dv
# M[K₃](s) should equal γ(s) = Γ_ℝ(s)Γ_ℝ(s+1)²
print("=== Test 1: Mellin transform of K₃ → γ(s) ===")

def K3(xi):
    """K₃(ξ) = 4∫₀^∞ exp(-2πξ/v - πv²) dv"""
    xi = mpmath.mpf(xi)
    if xi <= 0:
        return 4 * mpmath.quad(lambda v: mpmath.exp(-PI*v**2), [mpmath.mpf('1e-10'), 10])
    def f(v):
        if v < mpmath.mpf('1e-30'): return mpmath.mpf(0)
        return mpmath.exp(-2*PI*xi/v - PI*v**2)
    return 4 * mpmath.quad(f, [mpmath.mpf('1e-10'), 20], maxdegree=8)

for s_test in [2.0, 3.0, 4.0]:
    s = mpmath.mpf(s_test)
    # ∫₀^∞ K₃(t) t^{s-1} dt = γ(s)?
    def mellin_int(t):
        return K3(t) * mpmath.power(t, s - 1)
    t1 = time.time()
    m_val = mpmath.quad(mellin_int, [mpmath.mpf('1e-10'), 30], maxdegree=8)
    dt = time.time() - t1
    g011 = gamma_011(s)
    g012 = gamma_012(s)
    r011 = float(abs(m_val - g011) / abs(g011))
    r012 = float(abs(m_val - g012) / abs(g012))
    print(f"  s={s_test}: Mellin={float(m_val):.10e}")
    print(f"    γ₀₁₁={float(g011):.10e} rel={r011:.2e} {'✅' if r011<1e-6 else '❌'}")
    print(f"    γ₀₁₂={float(g012):.10e} rel={r012:.2e} {'✅' if r012<1e-6 else '❌'}")
    print(f"    ({dt:.1f}s)", flush=True)

# ── Test 2: F_s(n) via mpmath.quad vs Gauss-Laguerre ──
print("\n=== Test 2: F_s(1, 3) — quad vs Laguerre ===")

alpha_1 = 2 * mpmath.power(PI, mpmath.mpf('1.5')) / Q  # 2π^{3/2}/Q
s = mpmath.mpf(3)
coeff = 2 * mpmath.power(Q / (2*PI), s) / mpmath.power(PI, (s+1)/2)

# Method A: mpmath.quad (reference)
def integ_quad(t):
    if t < mpmath.mpf('1e-30'): return mpmath.mpf(0)
    a = alpha_1 / mpmath.sqrt(t)
    gi = mpmath.gammainc(s, a)
    return mpmath.power(t, (s-1)/2) * gi

t2a = time.time()
I_quad = mpmath.quad(integ_quad, [mpmath.mpf('1e-10'), 200], maxdegree=8)
F_quad = coeff * I_quad
dt2a = time.time() - t2a
print(f"  mpmath.quad: F_3(1) = {float(F_quad):.12e} ({dt2a:.1f}s)")

# Method B: Gauss-Laguerre 50
N_GL = 50
nodes, weights = np.polynomial.laguerre.laggauss(N_GL)
I_gl = mpmath.mpf(0)
for k in range(N_GL):
    t_k = mpmath.mpf(float(nodes[k]))
    w_k = mpmath.mpf(float(weights[k]))
    if t_k < mpmath.mpf('1e-30'): continue
    a_k = alpha_1 / mpmath.sqrt(t_k)
    gi = mpmath.gammainc(s, a_k)
    t_pow = mpmath.power(t_k, (s-1)/2)
    I_gl += w_k * t_pow * gi
F_gl = coeff * I_gl
print(f"  GL(50):     F_3(1) = {float(F_gl):.12e}")
print(f"  rel diff = {float(abs(F_gl-F_quad)/abs(F_quad)):.2e}")

# Method C: direct integral ∫₁^∞ K₃(x/Q) x^{s-1} dx
def F_direct_int(x):
    return K3(x / Q) * mpmath.power(x, s - 1)

t2c = time.time()
F_direct = mpmath.quad(F_direct_int, [1, 30], maxdegree=8)
dt2c = time.time() - t2c
print(f"  Direct:     F_3(1) = {float(F_direct):.12e} ({dt2c:.1f}s)")
print(f"  quad vs direct = {float(abs(F_quad-F_direct)/abs(F_direct)):.2e}", flush=True)

# ── Test 3: Check Λ_direct(3) decomposition ──
print("\n=== Test 3: Λ_direct(3) vs Σ c(n) (Q/n)^s γ(s) ===")

# Compute a few sym² coefficients
cn = {1: mpmath.mpf(1), 2: mpmath.mpf(1), 3: mpmath.mpf(-2)/3,
      4: mpmath.mpf(1), 5: mpmath.mpf(-4)/5, 7: mpmath.mpf(-3)/7,
      8: mpmath.mpf(1), 9: mpmath.mpf(1)/9, 11: mpmath.mpf(1)/11,
      13: mpmath.mpf(3)/13}

lam_direct = mpmath.mpf(0)
for n, c in cn.items():
    term = c * mpmath.power(Q, s) * gamma_011(s) / mpmath.power(n, s)
    lam_direct += term

# Full direct
cn_full_sum = mpmath.mpf(0)
for n, c in cn.items():
    cn_full_sum += c / mpmath.power(n, s)

lam_full = mpmath.power(Q, s) * gamma_011(s) * cn_full_sum
print(f"  Q^s γ(s) L_partial(s): {float(lam_full):.10e}")
print(f"  (using {len(cn)} terms)")

# ── Test 4: Full cycle: Λ_theta = Σ c(n)[F_s(n)+F_{1-s}(n)] vs Λ_direct ──
print("\n=== Test 4: Λ_theta(3) vs Λ_direct(3), n≤5 ===")

# F_s via quad for n=1..5
lam_theta = mpmath.mpf(0)
s1 = 1 - s  # = -2

cn_list = [(1, mpmath.mpf(1)), (2, mpmath.mpf(1)), (3, mpmath.mpf(-2)/3),
           (4, mpmath.mpf(1)), (5, mpmath.mpf(-4)/5)]

for n, c in cn_list:
    alpha_n = 2 * mpmath.power(PI, mpmath.mpf('1.5')) * n / Q
    coeff_s = 2 * mpmath.power(Q / (2*PI*n), s) / mpmath.power(PI, (s+1)/2)
    coeff_1s = 2 * mpmath.power(Q / (2*PI*n), s1) / mpmath.power(PI, (s1+1)/2)

    def make_int(s_val, al):
        def f(t):
            if t < mpmath.mpf('1e-30'): return mpmath.mpf(0)
            a = al / mpmath.sqrt(t)
            return mpmath.power(t, (s_val-1)/2) * mpmath.gammainc(s_val, a)
        return f

    I_s = mpmath.quad(make_int(s, alpha_n), [mpmath.mpf('1e-10'), 200], maxdegree=8)
    I_1s = mpmath.quad(make_int(s1, alpha_n), [mpmath.mpf('1e-10'), 200], maxdegree=8)

    Fs = coeff_s * I_s
    F1s = coeff_1s * I_1s
    contrib = c * (Fs + F1s)
    lam_theta += contrib
    print(f"  n={n}: c={float(c):.4f}, F_s={float(Fs):.8e}, F_{{1-s}}={float(F1s):.8e}, "
          f"contrib={float(contrib):.8e}")

lam_dir_partial = mpmath.power(Q, s) * gamma_011(s) * sum(c/mpmath.power(n, s) for n,c in cn_list)
print(f"\n  Λ_theta(3) ≈ {float(lam_theta):.10e}")
print(f"  Λ_direct(3) ≈ {float(lam_dir_partial):.10e} (same {len(cn_list)} terms)")
rel = float(abs(lam_theta - lam_dir_partial) / abs(lam_dir_partial))
print(f"  rel diff = {rel:.4e}")
print(flush=True)

# ── Test 5: Alternative — try Q=√11 instead of Q=11 ──
print("\n=== Test 5: Try Q=√11 (N=11 instead of N=121) ===")
Q_alt = mpmath.sqrt(11)
alpha_1_alt = 2 * mpmath.power(PI, mpmath.mpf('1.5')) / Q_alt

coeff_alt = 2 * mpmath.power(Q_alt / (2*PI), s) / mpmath.power(PI, (s+1)/2)
I_alt = mpmath.quad(lambda t: (mpmath.power(t, 1) * mpmath.gammainc(s, alpha_1_alt/mpmath.sqrt(t))
                               if t > mpmath.mpf('1e-30') else mpmath.mpf(0)),
                    [mpmath.mpf('1e-10'), 200], maxdegree=8)
F_alt = coeff_alt * I_alt
print(f"  F_3(1, Q=√11) = {float(F_alt):.10e}")
print(f"  Q_alt^3 γ(3) = {float(mpmath.power(Q_alt, 3) * gamma_011(3)):.10e}")
print(flush=True)

print("\n=== 진단 완료 ===")
