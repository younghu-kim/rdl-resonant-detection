#!/usr/bin/env python3
"""
[Project RDL] #72v2 — GL(4) sym³(Δ) 4성질 검증 + κ_near(d=4) 측정
=============================================================================
#72 기각 이유: L(s,Δ×Δ) = L(s,sym²Δ)·ζ(s) → 비원시적 (degree 3+1 분해)
#72v2: sym³(Δ) — Kim-Shahidi 정리로 확립된 진짜 원시적 degree 4 L-함수

L-함수: L(s, sym³Δ), degree=4, conductor N=1 (level 1 Δ)
  - Ramanujan Δ = Σ τ(n) q^n, weight 12
  - Satake: α_p + α_p^{-1} = t_p = τ(p)/p^{11/2} (Deligne: |α_p|=1)
  - sym³ Euler 인자 roots: {α³, α, α^{-1}, α^{-3}}
  - c_p = t_p³ - 2t_p  (소수에서)
  - prime power 재귀: c_{p^k} = e₁·c_{p^{k-1}} - e₂·c_{p^{k-2}} + e₁·c_{p^{k-3}} - c_{p^{k-4}}
    where e₁ = t³-2t, e₂ = t⁴-3t²+2
  - 완전 곱셈적: c_n = ∏_{p^a||n} c_{p^a}
  - Ramanujan bound: |c_p| ≤ 4 (degree 4, |α_p|=1이면 |α³+α+α^{-1}+α^{-3}| ≤ 4)

감마 인자: Hodge 구조 (33,0),(22,11),(11,22),(0,33) → μ=(0,1,11,12)
  Λ(s) = π^{-2s}·Γ(s/2)·Γ((s+1)/2)·Γ((s+11)/2)·Γ((s+12)/2)·L(s)
  log γ(s) = -2s·log(π) + logΓ(s/2) + logΓ((s+1)/2) + logΓ((s+11)/2) + logΓ((s+12)/2)

함수방정식: Λ(s) = Λ(1-s)  (ε=1, N=1, self-dual)

방법: Mellin-Barnes AFE (GL(3) v2 구조 → degree 4 확장)
      + 비자명 FE 검증 (Direct Dirichlet vs AFE 비교)
"""
import sys, os, time
import numpy as np
from scipy.special import loggamma
from scipy.signal import find_peaks

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "results", "gl4_sym3_delta_72v2.txt")
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))
def flush():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━ 파라미터 ━━━━━━
EPSILON  = 1             # root number (sym³Δ는 self-dual, ε=1)
N_COEFF  = 800           # Dirichlet 계수 수
C_SHIFT  = 4.0           # Mellin-Barnes 윤곽 이동
N_GH     = 60            # Gauss-Hermite 노드
T_MIN, T_MAX = 2.0, 50.0
DELTA    = 0.03          # σ 오프셋 (κ_near 측정: σ = 0.5 + DELTA)
MONO_R   = 0.3           # 모노드로미 반지름
MONO_N   = 64            # 모노드로미 단계
MONO_THR = 1.5           # 모노드로미 임계값

# GH 노드
GH_X, GH_W = np.polynomial.hermite.hermgauss(N_GH)


# ━━━━━━ 소수 체 ━━━━━━
def sieve_primes(limit):
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit+1, i):
                is_prime[j] = False
    return [i for i in range(2, limit+1) if is_prime[i]]


# ━━━━━━ Ramanujan τ(n) 계산 ━━━━━━
def compute_tau_fast(limit):
    """Ramanujan τ(n): Δ(q) = q · ∏_{n≥1}(1-q^n)^24.
    Python 정수 연산 사용 — float64는 p>150 에서 정밀도 손실.
    τ(m) = p[m-1] (Δ = q·p 이므로 shift).
    """
    # 정수 배열로 계산 (임의 정밀도)
    p = [0] * (limit + 1)
    p[0] = 1
    for n in range(1, limit + 1):
        # (1 - q^n)을 24번 곱함 (역방향 순서로 안전한 in-place)
        for _ in range(24):
            for k in range(limit, n - 1, -1):
                p[k] -= p[k - n]
    # τ(m) = p[m-1], float64로 변환
    tau = np.zeros(limit + 1, dtype=np.float64)
    for m in range(1, limit + 1):
        tau[m] = float(p[m - 1])
    return tau


# ━━━━━━ sym³(Δ) Dirichlet 계수 ━━━━━━
def compute_sym3_cn(tau, primes, limit):
    """sym³(Δ) Dirichlet 계수 c(n).

    Euler product at prime p:
      roots = {α³, α, α^{-1}, α^{-3}}  where α+α^{-1} = t_p = τ(p)/p^{11/2}

    Coefficients at prime powers (재귀):
      e₁ = t³ - 2t
      e₂ = t⁴ - 3t² + 2
      c_{p^0} = 1, c_{p^1} = e₁
      c_{p^k} = e₁·c_{p^{k-1}} - e₂·c_{p^{k-2}} + e₁·c_{p^{k-3}} - c_{p^{k-4}}

    Multiplicative: c_n = ∏_{p^a || n} c_{p^a}
    """
    # cpk[(p, k)] = c_{p^k}
    cpk = {}

    for p in primes:
        if p > limit:
            break
        tp = tau[p] / (p ** 5.5)  # t_p = τ(p) / p^{11/2}
        e1 = tp**3 - 2*tp
        e2 = tp**4 - 3*tp**2 + 2

        cpk[(p, 0)] = 1.0
        cpk[(p, 1)] = e1

        pk = p
        k = 1
        while pk * p <= limit:
            pk *= p
            k += 1
            # c_{p^k} = e₁·c_{p^{k-1}} - e₂·c_{p^{k-2}} + e₁·c_{p^{k-3}} - c_{p^{k-4}}
            c_km1 = cpk[(p, k-1)]
            c_km2 = cpk.get((p, k-2), 0.0)
            c_km3 = cpk.get((p, k-3), 0.0)
            c_km4 = cpk.get((p, k-4), 0.0)
            cpk[(p, k)] = e1*c_km1 - e2*c_km2 + e1*c_km3 - c_km4

    # 일반 n: 완전 곱셈성
    cn = np.zeros(limit + 1, dtype=np.float64)
    cn[0] = 0.0
    cn[1] = 1.0

    for n in range(2, limit + 1):
        temp = n
        result = 1.0
        for p in primes:
            if p * p > temp:
                break
            if temp % p == 0:
                k = 0
                while temp % p == 0:
                    k += 1
                    temp //= p
                result *= cpk.get((p, k), 0.0)
        if temp > 1:
            # temp는 소수
            result *= cpk.get((temp, 1), 0.0)
        cn[n] = result

    return cn


# ━━━━━━ sym³(Δ) 감마 인자 ━━━━━━
# μ = (0, 1, 11, 12) — Hodge 구조 (33,0),(22,11),(11,22),(0,33) 에서 유도
# Λ(s) = π^{-2s}·Γ(s/2)·Γ((s+1)/2)·Γ((s+11)/2)·Γ((s+12)/2)·L(s)
# log γ(s) = -2s·log(π) + logΓ(s/2) + logΓ((s+1)/2) + logΓ((s+11)/2) + logΓ((s+12)/2)

MU = [0.0, 1.0, 11.0, 12.0]  # μ 값

def log_gamma_factor_sym3(s):
    """sym³(Δ) 감마 인자 (로그)
    log γ(s) = -2s·log(π) + Σ_{j} logΓ((s+μ_j)/2)
    """
    s = complex(s)
    log_pi = np.log(np.pi)
    result = -2 * s * log_pi
    for mu in MU:
        result += loggamma((s + mu) / 2)
    return result


def Lambda_direct(s, cn_nz, log_n_nz, cn_vals_nz):
    """Λ_dir(s) = γ(s) · L(s), Re(s) > 1 필요"""
    s = complex(s)
    log_g = log_gamma_factor_sym3(s)
    gamma_s = np.exp(log_g)
    powers = np.exp(-s * log_n_nz)
    L_s = np.dot(cn_vals_nz, powers)
    return gamma_s * L_s


def L_direct_series(s, cn_nz, log_n_nz, cn_vals_nz):
    """L(s) = Σ c_n/n^s (직접 급수, Re(s) > 1 필요)"""
    s = complex(s)
    powers = np.exp(-s * log_n_nz)
    return np.dot(cn_vals_nz, powers)


def Lambda_AFE(s, cn_nz, log_n_nz, cn_vals_nz, c=C_SHIFT, eps=EPSILON):
    """
    Mellin-Barnes AFE:
    Λ(s) ≈ (e^{c²}/2π) Σ_k w_k [Λ_dir(s+c+iy_k) + ε·Λ_dir(1-s+c+iy_k)]
             × e^{2icy_k}/(c+iy_k)
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


def L_from_AFE(s, cn_nz, log_n_nz, cn_vals_nz):
    """Λ_AFE / γ(s) = L(s) — AFE 기반 해석적 연속"""
    s = complex(s)
    lam = Lambda_AFE(s, cn_nz, log_n_nz, cn_vals_nz)
    g = np.exp(log_gamma_factor_sym3(s))
    if abs(g) < 1e-300:
        return complex('nan')
    return lam / g


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
    return float(k) if np.isfinite(k) else 1e12


def monodromy(t_center, cn_nz, log_n_nz, cn_vals_nz,
              sigma=0.5, radius=MONO_R, n_steps=MONO_N):
    """폐곡선 모노드로미: 총 위상/π"""
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


# ━━━━━━ 메인 ━━━━━━
log("=" * 70)
log("결과 #72v2 — GL(4) sym³(Δ) 4성질 검증 + κ_near(d=4)")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"degree=4, conductor N=1, root number ε={EPSILON}")
log(f"L-함수: L(s, sym³Δ), Δ=Ramanujan delta (weight 12)")
log(f"감마: μ=(0,1,11,12), log γ(s) = -2s·log(π) + logΓ(s/2) + logΓ((s+1)/2) + logΓ((s+11)/2) + logΓ((s+12)/2)")
log(f"c_shift={C_SHIFT}, N_GH={N_GH}, N_coeff={N_COEFF}")
log(f"t 범위: [{T_MIN}, {T_MAX}], δ={DELTA}, mono_r={MONO_R}")
log()

t_total = time.time()

# ── Step 0: τ(n) 계산 ──
log("[Step 0] Ramanujan τ(n) 계산 (∏(1-q^n)^24 전개)")
t0 = time.time()
tau = compute_tau_fast(N_COEFF)
log(f"  τ(1)={tau[1]:.0f}, τ(2)={tau[2]:.0f}, τ(3)={tau[3]:.0f}, τ(4)={tau[4]:.0f}")
log(f"  τ(5)={tau[5]:.0f}, τ(6)={tau[6]:.0f}, τ(7)={tau[7]:.0f}")
known_tau = {1:1, 2:-24, 3:252, 4:-1472, 5:4830, 6:-6048, 7:-16744}
ok_count = sum(1 for n, v in known_tau.items() if abs(tau[n] - v) < 0.5)
log(f"  알려진 τ(n) 검증: {ok_count}/{len(known_tau)} 일치")
if ok_count < len(known_tau):
    for n, v in known_tau.items():
        if abs(tau[n] - v) >= 0.5:
            log(f"  ⚠️ τ({n}): 계산={tau[n]:.0f}, 기대={v}")
    log("  ⚠️ τ(n) 계산 오류 — 종료")
    flush()
    sys.exit(1)
log(f"  ✅ τ(n) 계산 정확 (소요: {time.time()-t0:.1f}s)")
flush()

# ── Step 0b: sym³ 계수 계산 ──
log(f"\n[Step 0b] sym³(Δ) Dirichlet 계수 c(n) 계산")
t0b = time.time()
primes = sieve_primes(N_COEFF)
cn = compute_sym3_cn(tau, primes, N_COEFF)

# 예상값 검증: c_2
t2 = tau[2] / (2.0 ** 5.5)
e1_2 = t2**3 - 2*t2
log(f"  t_2 = τ(2)/2^{{11/2}} = {t2:.6f}")
log(f"  c_2 = t_2³ - 2t_2 = {e1_2:.6f}  (기대 ≈ 0.9115)")
log(f"  계산된 c(2) = {cn[2]:.6f}")
if abs(cn[2] - e1_2) > 1e-6:
    log(f"  ⚠️ c(2) 불일치 — 계수 계산 오류")
    flush()
    sys.exit(1)

log(f"  c(1)={cn[1]:.4f} (기대=1)")
log(f"  c(2)={cn[2]:.6f}")
log(f"  c(3)={cn[3]:.6f}")
log(f"  c(4)={cn[4]:.6f}")
log(f"  c(5)={cn[5]:.6f}")

# Ramanujan bound: |c_p| ≤ 4
raman_viol = [(p, cn[p]) for p in primes[:50] if abs(cn[p]) > 4.0 + 1e-10]
if raman_viol:
    log(f"  ⚠️ Ramanujan bound 위반: {raman_viol[:5]}")
else:
    log(f"  ✅ Ramanujan bound |c_p| ≤ 4: 소수 p ≤ {primes[min(49,len(primes)-1)]} 모두 통과")
    max_cp = max(abs(cn[p]) for p in primes if p <= N_COEFF)
    log(f"     max|c_p| = {max_cp:.6f} (≤ 4 확인)")

log(f"  소요: {time.time()-t0b:.1f}s")

# 비영 인덱스
nz_idx = np.where(np.abs(cn[1:]) > 1e-12)[0]
cn_nz = cn[nz_idx + 1]
log_n_nz = np.log((nz_idx + 1).astype(float))
log(f"  비영 계수: {len(nz_idx)}개 / {N_COEFF}개")
flush()

# ── Step 1a: 비자명 FE 검증 (핵심!) ──
log(f"\n[Step 1a] ★ 비자명 FE 검증: L_direct vs L_AFE (Re(s)=2)")
log(f"  목표: |L_direct(s) - L_AFE(s)| / |L_direct(s)| < 1e-6")
log(f"  (#72 기각의 핵심 — Dirichlet 계수 + 감마 인자 동시 검증)")
nontrivial_pts = [2.0+5j, 2.0+10j, 2.0+15j, 2.0+20j, 2.0+25j]
nontrivial_errs = []
t1a = time.time()

for s in nontrivial_pts:
    try:
        # 직접 급수 (Re(s)=2에서 수렴)
        L_dir = L_direct_series(s, cn_nz, log_n_nz, cn_nz)
        # AFE 기반 L(s)
        L_afe = L_from_AFE(s, cn_nz, log_n_nz, cn_nz)
        denom = max(abs(L_dir), 1e-300)
        rel = abs(L_dir - L_afe) / denom
        nontrivial_errs.append(rel)
        ok = "✅" if rel < 1e-6 else "❌"
        log(f"  s={s}: |L_dir|={abs(L_dir):.4e}, |L_afe|={abs(L_afe):.4e}, rel={rel:.2e} {ok}")
    except Exception as e:
        log(f"  s={s}: ERROR {e}")
        nontrivial_errs.append(1.0)

nt_pass = sum(1 for e in nontrivial_errs if e < 1e-6)
max_nt = max(nontrivial_errs) if nontrivial_errs else 999
log(f"  소요: {time.time()-t1a:.1f}s")
log(f"\n  비자명 FE pass: {nt_pass}/{len(nontrivial_errs)}, max_rel={max_nt:.2e}")

if nt_pass < 3:
    log(f"  ❌ 비자명 FE 검증 실패 — 감마 μ 또는 계수 오류")
    log(f"  → 대안 μ=(1,2,11,12) 시도 필요")
    # 계속 진행 (AFE 자체는 동작하므로)
    log(f"  ⚠️ μ=(0,1,11,12) 실패. 대안 μ 시도...")
    # 감마 인자 변경 시도
    MU_ALT = [1.0, 2.0, 11.0, 12.0]
    def log_gamma_factor_alt(s):
        s = complex(s)
        log_pi = np.log(np.pi)
        result = -2 * s * log_pi
        for mu in MU_ALT:
            result += loggamma((s + mu) / 2)
        return result

    # 잠시 교체해서 테스트
    orig_lgf = log_gamma_factor_sym3

    def Lambda_direct_alt(s, cn_nz, log_n_nz, cn_vals_nz):
        s = complex(s)
        log_g = log_gamma_factor_alt(s)
        gamma_s = np.exp(log_g)
        powers = np.exp(-s * log_n_nz)
        L_s = np.dot(cn_vals_nz, powers)
        return gamma_s * L_s

    def Lambda_AFE_alt(s, cn_nz, log_n_nz, cn_vals_nz, c=C_SHIFT, eps=EPSILON):
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
            lam_f = Lambda_direct_alt(u_fwd, cn_nz, log_n_nz, cn_vals_nz)
            lam_b = Lambda_direct_alt(u_bwd, cn_nz, log_n_nz, cn_vals_nz)
            total += w * (lam_f + eps * lam_b) * phase
        return ec2 / (2*np.pi) * total

    log(f"  대안 μ=(1,2,11,12) 시도 중...")
    alt_errs = []
    for s in nontrivial_pts[:3]:
        try:
            L_dir = L_direct_series(s, cn_nz, log_n_nz, cn_nz)
            lam_alt = Lambda_AFE_alt(s, cn_nz, log_n_nz, cn_nz)
            g_alt = np.exp(log_gamma_factor_alt(s))
            L_afe_alt = lam_alt / g_alt if abs(g_alt) > 1e-300 else complex('nan')
            rel = abs(L_dir - L_afe_alt) / max(abs(L_dir), 1e-300)
            alt_errs.append(rel)
            log(f"  alt s={s}: rel={rel:.2e}")
        except Exception as e:
            log(f"  alt s={s}: ERROR {e}")
            alt_errs.append(1.0)

    alt_pass = sum(1 for e in alt_errs if e < 1e-6)
    if alt_pass >= 3:
        log(f"  ✅ 대안 μ=(1,2,11,12) 통과! 감마 인자 교체")
        MU = MU_ALT
        log_gamma_factor_sym3 = log_gamma_factor_alt
        Lambda_direct = Lambda_direct_alt
        Lambda_AFE_func = Lambda_AFE_alt
        log(f"  이후 모든 계산에 μ=(1,2,11,12) 사용")
        # 비자명 FE 통과로 갱신
        nt_pass = alt_pass
        nontrivial_errs = alt_errs
    else:
        log(f"  ❌ 대안도 실패. μ 오류 또는 급수 수렴 부족 (N_COEFF 증가 필요?)")
        log(f"  → 계속 진행 (N_COEFF=800이 부족할 수 있음. 내부 일관성 확인)")
else:
    log(f"  ✅ 비자명 FE 검증 통과 — 계수+감마 인자 정합 확인")
flush()

# ── Step 1b: AFE 함수방정식 검증 (참고용) ──
log(f"\n[Step 1b] AFE 함수방정식 검증 Λ_AFE(s) = Λ_AFE(1-s) (참고용)")
fe_pts = [0.5+5j, 0.5+10j, 0.5+15j, 0.7+8j, 0.3+12j]
fe_errs = []
t1b = time.time()
for s in fe_pts:
    try:
        Ls  = Lambda_AFE(s,   cn_nz, log_n_nz, cn_nz)
        L1s = Lambda_AFE(1-s, cn_nz, log_n_nz, cn_nz)
        denom = max(abs(Ls), abs(L1s), 1e-300)
        rel = abs(Ls - EPSILON*L1s) / denom
        fe_errs.append(rel)
        log(f"  s={s}: |Λ(s)|={abs(Ls):.4e}, |Λ(1-s)|={abs(L1s):.4e}, rel={rel:.2e}")
    except Exception as e:
        log(f"  s={s}: ERROR {e}")
        fe_errs.append(1.0)
log(f"  (주의: AFE 구조상 rel≈0은 자명할 수 있음 — Step 1a가 실질 검증)")
log(f"  소요: {time.time()-t1b:.1f}s")
flush()

# ── Step 2: 영점 탐색 ──
log(f"\n[Step 2] Λ(1/2+it) 부호 변환 탐색 (t ∈ [{T_MIN},{T_MAX}], Δt=0.2)")
scan_ts = np.arange(T_MIN, T_MAX + 0.05, 0.2)
scan_re = np.zeros(len(scan_ts))
t2_start = time.time()

for i, t in enumerate(scan_ts):
    try:
        val = Lambda_AFE(complex(0.5, t), cn_nz, log_n_nz, cn_nz)
        scan_re[i] = val.real
    except Exception as e:
        log(f"  WARNING t={t:.1f}: {e}")
        scan_re[i] = 0.0
    if i % 40 == 0 or i == len(scan_ts)-1:
        elapsed = time.time() - t2_start
        eta = elapsed/(i+1)*(len(scan_ts)-i-1) if i > 0 else 0
        log(f"  t={t:5.1f}: Λ={scan_re[i]:+.6e}  ({i+1}/{len(scan_ts)}, eta={eta:.0f}s)")
        flush()

n_sign = int(np.sum(np.diff(np.sign(scan_re)) != 0))
log(f"\n  부호 변환: {n_sign}회, 소요: {time.time()-t2_start:.1f}s")
if n_sign == 0:
    log(f"  ⚠️ 부호 변환 0회! 범위: [{scan_re.min():.4e}, {scan_re.max():.4e}]")
    flush()
    sys.exit(1)
flush()

# 이분법 영점 근사
log(f"\n[Step 2b] 이분법 영점 근사")
zeros = []
for i in range(len(scan_re)-1):
    if scan_re[i] * scan_re[i+1] < 0:
        lo, hi = scan_ts[i], scan_ts[i+1]
        v_lo = scan_re[i]
        for _ in range(40):
            mid = (lo+hi)/2
            try:
                v_mid = Lambda_AFE(complex(0.5, mid), cn_nz, log_n_nz, cn_nz).real
            except Exception:
                break
            if v_mid * v_lo < 0: hi = mid
            else: lo = mid; v_lo = v_mid
        zeros.append((lo+hi)/2)
        log(f"  영점 #{len(zeros):2d}: t ≈ {(lo+hi)/2:.8f}")

if len(zeros) == 0:
    log("  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
    flush()
    sys.exit(1)

log(f"\n  발견: {len(zeros)}개 영점")
flush()

# ── Step 3: κ_near 측정 + 모노드로미 ──
log(f"\n[Step 3] κ_near(d=4) 측정 + 모노드로미")
sigma_near = 0.5 + DELTA
tp_results = []

for idx, z_t in enumerate(zeros):
    t3 = time.time()
    try:
        k_val = curvature(sigma_near, z_t, cn_nz, log_n_nz, cn_nz)
    except Exception as e:
        k_val = float('nan')
        log(f"  WARNING κ 계산 오류 t={z_t:.4f}: {e}")

    try:
        m_val = monodromy(z_t, cn_nz, log_n_nz, cn_nz)
    except Exception as e:
        m_val = None
        log(f"  WARNING mono 계산 오류 t={z_t:.4f}: {e}")

    tp_results.append((z_t, k_val, m_val))
    m_str = f"{m_val:.4f}" if m_val is not None else "FAIL"
    log(f"  TP #{idx+1:2d} t={z_t:.6f}: κ={k_val:.2f}, mono/π={m_str}  ({time.time()-t3:.1f}s)")
    flush()

# κ_near 통계
kappa_vals = [k for _, k, _ in tp_results if np.isfinite(k) and k < 1e11]
if kappa_vals:
    k_mean = np.mean(kappa_vals)
    k_std  = np.std(kappa_vals)
    k_cv   = k_std / k_mean * 100 if k_mean > 0 else 999
    log(f"\n  κ_near(d=4) = {k_mean:.2f} ± {k_std:.2f} (CV={k_cv:.1f}%, n={len(kappa_vals)})")
else:
    k_mean = float('nan')
    k_cv   = 999
    log(f"\n  ⚠️ κ_near 계산 실패")

# ── Step 4: κ near/far 비율 ──
log(f"\n[Step 4] κ near/far 비율")
far_ts = np.arange(T_MIN, T_MAX, 1.5)
k_near_list, k_far_list = [], []
for t in far_ts:
    try:
        k = curvature(sigma_near, t, cn_nz, log_n_nz, cn_nz)
    except Exception as e:
        log(f"  WARNING far κ t={t:.1f}: {e}")
        continue
    if not np.isfinite(k) or k >= 1e11:
        continue
    min_d = min(abs(t - z) for z in zeros) if zeros else 999
    if min_d < 0.5:
        k_near_list.append(k)
    else:
        k_far_list.append(k)

if k_near_list and k_far_list:
    near_med = np.median(k_near_list)
    far_med  = np.median(k_far_list)
    ratio    = near_med / far_med if far_med > 0 else float('inf')
    log(f"  near: median={near_med:.2f} (n={len(k_near_list)})")
    log(f"  far:  median={far_med:.4f} (n={len(k_far_list)})")
    log(f"  ratio= {ratio:.1f}×  {'✅' if ratio >= 10 else '⚠️ 낮음 (영점 밀도 포화?)'}")
else:
    ratio = 0.0
    log(f"  데이터 부족 (near={len(k_near_list)}, far={len(k_far_list)})")
flush()

# ── Step 5: 모노드로미 통계 ──
log(f"\n[Step 5] 모노드로미 통계")
mono_vals = [m for _, _, m in tp_results if m is not None]
if mono_vals:
    tp_mono = [m for m in mono_vals if m >= MONO_THR]
    fp_mono = [m for m in mono_vals if m < MONO_THR]
    n_tp_mono = len(tp_mono)
    log(f"  TP (mono/π ≥ {MONO_THR}): {n_tp_mono}/{len(zeros)}, "
        + (f"mean={np.mean(tp_mono):.4f}" if tp_mono else ""))
    log(f"  FP: {len(fp_mono)}개")
    if tp_mono:
        log(f"  mono/π 범위: [{min(tp_mono):.4f}, {max(tp_mono):.4f}]")
else:
    n_tp_mono = 0
    log(f"  ⚠️ mono 계산 전부 실패")
flush()

# ── Step 6: σ-유일성 ──
log(f"\n[Step 6] σ-유일성 검사 (N=1 → PASS 예측)")
sigma_vals = [0.3, 0.5, 0.7, 0.9]
sigma_sc = {}
for sig in sigma_vals:
    vals = []
    for t in np.arange(T_MIN, T_MAX, 0.8):
        try:
            v = Lambda_AFE(complex(sig, t), cn_nz, log_n_nz, cn_nz)
            vals.append(v.real)
        except Exception as e:
            log(f"  WARNING σ={sig} t={t:.1f}: {e}")
            vals.append(0.0)
    sc = int(np.sum(np.diff(np.sign(vals)) != 0))
    sigma_sc[sig] = sc
    log(f"  σ={sig:.1f}: 부호변환={sc}")

sc_half = sigma_sc.get(0.5, 0)
is_uniq = sc_half > 0 and all(sc_half >= sigma_sc[s] for s in sigma_vals if s != 0.5)
log(f"\n  σ=0.5 최대: {'✅ PASS (N=1 패턴 확인)' if is_uniq else '❌ FAIL'}")
flush()

# ── 종합 ──
log(f"\n{'='*70}")
log(f"종합 판정 — GL(4) sym³(Δ)")
log(f"{'='*70}")
nt_pass_final = sum(1 for e in nontrivial_errs if e < 1e-6)
fe_pass_final = sum(1 for e in fe_errs if e < 1e-6) if fe_errs else 0
log(f"  [P0] ★비자명 FE 검증: {nt_pass_final}/{len(nontrivial_errs)}, max_rel={max_nt:.2e}")
log(f"       {'✅ PASS (계수+감마 정합)' if nt_pass_final >= 3 else '❌ FAIL'}")
log(f"  [P1] AFE 함수방정식: {fe_pass_final}/{len(fe_errs)}")
log(f"       (참고용 — 자명할 수 있음)")
log(f"  [P2] 영점 발견: {len(zeros)}개 (t∈[{T_MIN},{T_MAX}])")
log(f"       {'✅ PASS' if len(zeros) >= 5 else '⚠️ 부족'}")
if np.isfinite(k_mean):
    log(f"  [P3] κ_near(d=4): {k_mean:.2f} ± {k_std:.2f} (CV={k_cv:.1f}%, n={len(kappa_vals)})")
    log(f"       {'✅ PASS' if k_cv < 2 else '⚠️ 산포 큼'}")
    log(f"       수학자 예측: 1125 < κ < 1150")
log(f"  [P4] 모노드로미: {n_tp_mono}/{len(zeros)} TP")
log(f"       {'✅ PASS' if n_tp_mono == len(zeros) and len(zeros) > 0 else '⚠️ 일부 실패'}")
log(f"  [B-10] σ-유일성(N=1): {'✅ PASS' if is_uniq else '❌ FAIL (예상 밖)'}")
log(f"  κ ratio: {ratio:.1f}×")

log(f"\n  ━━ κ_near(d) 비교 (B-12 핵심) ━━")
log(f"  d=1 (ζ):           κ_near = 1112.32")
log(f"  d=2 (GL2 avg):     κ_near = 1114.60")
log(f"  d=3 (GL3 sym²):    κ_near = 1125.16")
if np.isfinite(k_mean):
    log(f"  d=4 (sym³Δ):       κ_near = {k_mean:.2f}  ← 이번 측정 (★)")
    gap_34 = k_mean - 1125.16
    gap_23 = 1125.16 - 1114.60
    gap_12 = 1114.60 - 1112.32
    log(f"  gap(1→2) = {gap_12:.2f}")
    log(f"  gap(2→3) = {gap_23:.2f}")
    log(f"  gap(3→4) = {gap_34:.2f}")
    if gap_34 > 0:
        log(f"  → 단조증가 4점 확인: ✅")
        if gap_23 > 0:
            accel = gap_34 / gap_23
            log(f"  → 가속비 (3→4)/(2→3) = {accel:.2f}×")
    elif gap_34 < -0.5:
        log(f"  → ⚠️ 단조증가 위반 (gap<0) — sym³ 영점 특성 검토 필요")
    else:
        log(f"  → 거의 동일 (gap≈0) — κ_near(d) 수렴 가능성")

log(f"\n  ━━ κ ratio degree 비교 (B-13) ━━")
log(f"  d=1: 2200.7×")
log(f"  d=2: ~600×")
log(f"  d=3: ~320×")
log(f"  d=4: {ratio:.1f}×")
if ratio > 0:
    log(f"  → degree 반비례 패턴: {'✅ 일관' if ratio < 320 else '⚠️ 예상 밖 (영점 밀도 포화?)'}")

total = time.time() - t_total
log(f"\n  총 소요: {total:.1f}초 ({total/60:.1f}분)")
log(f"  시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
flush()
log(f"\n결과 저장: {OUTFILE}")
flush()
