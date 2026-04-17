#!/usr/bin/env python3
"""
[Project RDL] #72 — GL(4) Rankin-Selberg L(s, Δ×Δ) 4성질 검증 + κ_near(d=4) 측정
=============================================================================
L-함수: Rankin-Selberg convolution L(s, Δ×Δ), degree=4, conductor N=1
  - Ramanujan Δ = Σ τ(n) q^n, weight 12
  - Dirichlet 계수: a_n(Δ×Δ) = Σ_{d|n} τ(d)·τ(n/d)
  - 감마 인자: Γ_ℝ(s+11)² · Γ_ℝ(s+12)²  (μ = 11,11,12,12)
    즉 log γ(s) = -2s·log(π) + 2·logΓ((s+11)/2) + 2·logΓ((s+12)/2)
  - root number ε = 1, Q = 1/π² (N=1, degree=4)
  - Λ(s) = π^{-2s} · Γ((s+11)/2)² · Γ((s+12)/2)² · L(s, Δ×Δ)
  - 함수방정식: Λ(s) = Λ(1-s)

방법: Mellin-Barnes AFE (GL(3) v2 구조 그대로 degree=4로 교체)
"""
import sys, os, time
import numpy as np
from scipy.special import loggamma
from scipy.signal import find_peaks

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "results", "gl4_rankin_selberg_72.txt")
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))
def flush():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━ 파라미터 ━━━━━━
EPSILON  = 1             # root number (Δ×Δ는 self-dual, ε=1)
N_COEFF  = 800           # Dirichlet 계수 수 (degree 4 → 더 많이 필요)
C_SHIFT  = 4.0           # Mellin-Barnes 윤곽 이동
N_GH     = 60            # Gauss-Hermite 노드
T_MIN, T_MAX = 2.0, 50.0
DELTA    = 0.03          # σ 오프셋 (κ_near 측정)
MONO_R   = 0.3           # 모노드로미 반지름
MONO_N   = 64            # 모노드로미 단계
MONO_THR = 1.5           # 모노드로미 임계값

# GH 노드
GH_X, GH_W = np.polynomial.hermite.hermgauss(N_GH)

# ━━━━━━ Ramanujan τ(n) 계수 계산 ━━━━━━
def compute_tau(limit):
    """Ramanujan τ(n): Δ = Σ τ(n) q^n, weight 12.
    Ramanujan tau 공식 (소수에서 직접 계산):
    τ(p) = 알려진 값 사용 + 재귀.
    여기서는 Δ(q) = q·∏_{n≥1}(1-q^n)^24 의 계수를 직접 계산.
    """
    # q·∏(1-q^n)^24 계수 계산 (eta^24)
    # eta(q) = q^{1/24}·∏(1-q^n) → eta^24 = q·∏(1-q^n)^24 = Σ τ(n) q^n
    tau = [0] * (limit + 1)
    # 먼저 ∏(1-q^n)^24 를 0부터 limit-1 차수까지 계산
    # Δ(q) = q · prod_{n>=1} (1-q^n)^24
    # 계수: tau[n] = n번째 항 (1-indexed: Δ = Σ_{n≥1} τ(n) q^n)

    # 방법: log(∏(1-q^n)^24) = 24 · Σ log(1-q^n)
    # 하지만 정수 계산 위해 직접 곱

    # p = ∏(1-q^n)^24, p[k] = k차 계수
    p = np.zeros(limit, dtype=np.float64)
    p[0] = 1.0
    for n in range(1, limit):
        # (1-q^n)^24 를 현재 p에 곱함
        # (1-x)^24 = Σ_{k=0}^{24} C(24,k)(-1)^k x^k
        from math import comb
        # 고차 q-급수에 (1-q^n)^24 곱하기
        for k in range(24, 0, -1):
            coeff = ((-1)**k) * comb(24, k)
            # p[j] -= coeff * p[j - k*n]  (단 j - k*n >= 0)
            step = k * n
            if step < limit:
                p[step:] += coeff * p[:limit - step]

    # tau[n] = p[n-1] (Δ = q · p(q) 이므로 tau[n] = 계수 of q^{n-1} in p)
    tau_arr = np.zeros(limit + 1, dtype=np.float64)
    for n in range(1, limit + 1):
        if n - 1 < len(p):
            tau_arr[n] = p[n - 1]
    return tau_arr


def compute_tau_fast(limit):
    """Ramanujan τ(n): Δ(q) = q · ∏_{n≥1}(1-q^n)^24.
    올바른 알고리즘: (1-q^n)을 24번 반복 적용.
    각 곱: p[n:] -= p[:limit+1-n].copy()  (안전한 in-place).
    τ(m) = p[m-1] (Δ = q·p 이므로 shift).
    """
    p = np.zeros(limit + 1, dtype=np.float64)
    p[0] = 1.0

    for n in range(1, limit + 1):
        if n > limit:
            break
        # (1 - q^n)을 24번 곱함
        for _ in range(24):
            p[n:] -= p[:limit + 1 - n].copy()

    # τ(m) = p[m-1]
    tau = np.zeros(limit + 1, dtype=np.float64)
    for m in range(1, limit + 1):
        tau[m] = p[m - 1]
    return tau


def compute_rankin_selberg_an(tau, limit):
    """L(s, Δ×Δ) 계수: a_n = Σ_{d|n} τ(d)·τ(n/d).
    이것은 단순히 τ와 자기 자신의 Dirichlet convolution.
    주의: L(s, Δ×Δ) = Σ a_n / n^s 에서 a_n = (τ*τ)(n)
    """
    an = np.zeros(limit + 1, dtype=np.float64)
    for n in range(1, limit + 1):
        # Σ_{d|n} τ(d)·τ(n/d)
        s = 0.0
        d = 1
        while d * d <= n:
            if n % d == 0:
                s += tau[d] * tau[n // d]
                if d != n // d:
                    s += tau[n // d] * tau[d]
            d += 1
        an[n] = s
    return an


# ━━━━━━ GL(4) 감마 인자 ━━━━━━
# Λ(s) = π^{-2s} · Γ((s+11)/2)^2 · Γ((s+12)/2)^2 · L(s)
# log Λ(s) = -2s·log(π) + 2·logΓ((s+11)/2) + 2·logΓ((s+12)/2) + log L(s)
# 함수방정식: Λ(s) = ε · Λ(1-s), ε=1 (Δ×Δ self-dual)
#
# 참고: Rankin-Selberg Δ×Δ의 정확한 감마 인자는
# π^{-2s} · Γ((s+11)/2)^2 · Γ((s+12)/2)^2
# (Δ weight=12 → k=12, αj = (k-1)/2 = 11/2, (k+1)/2 - 1/2 = k/2 = 6)
# Shahidi 공식: μ_1=μ_2=11, μ_3=μ_4=12 (or 0,0,11,11 depending on convention)
# 우리는 LMFDB 관례 따름: Γ_ℝ(s+μ) = π^{-(s+μ)/2} Γ((s+μ)/2)
# L(s,Δ×Δ)의 completed form: (2π)^{-2s} Γ(s+11) Γ(s+12) L(s) 형태도 있음
# 여기서는 μ_j = 0,1,11,12 (일부 reference) 또는 11,11,12,12 (다른 reference)
# 함수방정식으로 확인할 것.

def log_gamma_factor_gl4(s):
    """
    GL(4) Rankin-Selberg Δ×Δ 감마 인자.
    Λ(s) = π^{-2s} · Γ((s+11)/2)² · Γ((s+12)/2)² · L(s)
    log γ(s) = -2s·log(π) + 2·logΓ((s+11)/2) + 2·logΓ((s+12)/2)

    함수방정식 체크: s ↔ 1-s 대칭
    s=1/2: log γ(1/2) = -log(π) + 2·logΓ(23/4) + 2·logΓ(25/4)
    1-s=1/2: 동일 → ε=1 이면 Λ(1/2)=Λ(1/2) 자동
    """
    return (-2*s * np.log(np.pi) +
            2*loggamma((s + 11) / 2) +
            2*loggamma((s + 12) / 2))


def Lambda_direct(s, an_nz, log_n_nz, an_vals_nz):
    """Λ_dir(s) = γ(s) · L(s), Re(s) > 1 필요"""
    s = complex(s)
    log_g = log_gamma_factor_gl4(s)
    gamma_s = np.exp(log_g)
    powers = np.exp(-s * log_n_nz)
    L_s = np.dot(an_vals_nz, powers)
    return gamma_s * L_s


def Lambda_AFE(s, an_nz, log_n_nz, an_vals_nz, c=C_SHIFT, eps=EPSILON):
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
        lam_f = Lambda_direct(u_fwd, an_nz, log_n_nz, an_vals_nz)
        lam_b = Lambda_direct(u_bwd, an_nz, log_n_nz, an_vals_nz)
        total += w * (lam_f + eps * lam_b) * phase
    return ec2 / (2*np.pi) * total


def curvature(sigma, t, an_nz, log_n_nz, an_vals_nz, h=1e-5):
    """κ(σ+it) = |Λ'/Λ|²"""
    s = complex(sigma, t)
    L0 = Lambda_AFE(s, an_nz, log_n_nz, an_vals_nz)
    if abs(L0) < 1e-250:
        return 1e12
    Lp = Lambda_AFE(s + h, an_nz, log_n_nz, an_vals_nz)
    Lm = Lambda_AFE(s - h, an_nz, log_n_nz, an_vals_nz)
    conn = (Lp - Lm) / (2*h * L0)
    k = abs(conn)**2
    return float(k) if np.isfinite(k) else 1e12


def monodromy(t_center, an_nz, log_n_nz, an_vals_nz,
              sigma=0.5, radius=MONO_R, n_steps=MONO_N):
    """폐곡선 모노드로미: 총 위상/π"""
    center = complex(sigma, t_center)
    phase_acc = 0.0
    prev = None
    for j in range(n_steps + 1):
        th = 2*np.pi*j / n_steps
        pt = center + radius * np.exp(1j*th)
        val = Lambda_AFE(pt, an_nz, log_n_nz, an_vals_nz)
        if abs(val) < 1e-250:
            return None
        if prev is not None:
            phase_acc += np.angle(val / prev)
        prev = val
    return abs(phase_acc) / np.pi


# ━━━━━━ 메인 ━━━━━━
log("=" * 70)
log("결과 #72 — GL(4) Rankin-Selberg L(s, Δ×Δ) 4성질 검증 + κ_near(d=4)")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"degree=4, conductor N=1, root number ε={EPSILON}")
log(f"L-함수: L(s, Δ×Δ), Δ=Ramanujan delta (weight 12)")
log(f"감마: π^{{-2s}} · Γ((s+11)/2)² · Γ((s+12)/2)²")
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
# 알려진 값 체크: τ(1)=1, τ(2)=-24, τ(3)=252, τ(4)=-1472, τ(5)=4830
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

# ── Rankin-Selberg 계수 ──
log(f"\n[Step 0b] L(s, Δ×Δ) Dirichlet 계수 a_n = (τ*τ)(n)")
t0b = time.time()
an = compute_rankin_selberg_an(tau, N_COEFF)
log(f"  a(1)={an[1]:.0f} (기대=1)")
log(f"  a(2)={an[2]:.0f} (기대=τ(1)τ(2)+τ(2)τ(1)={2*tau[2]:.0f})")
log(f"  a(3)={an[3]:.0f}")
log(f"  a(4)={an[4]:.0f}")
# a(2) = 2τ(2) = 2*(-24) = -48
log(f"  소요: {time.time()-t0b:.1f}s")

# 비영 인덱스
nz_idx = np.where(np.abs(an[1:]) > 1e-10)[0]
an_nz = an[nz_idx + 1]
log_n_nz = np.log((nz_idx + 1).astype(float))
log(f"  비영 계수: {len(nz_idx)}개 / {N_COEFF}개")
flush()

# ── Step 1: 함수방정식 검증 ──
log(f"\n[Step 1] 함수방정식 Λ(s) = ε·Λ(1-s)")
fe_pts = [0.5+5j, 0.5+10j, 0.5+15j, 0.7+8j, 0.3+12j, 0.6+25j, 0.5+35j]
fe_errs = []
t1 = time.time()
for s in fe_pts:
    try:
        Ls  = Lambda_AFE(s,   an_nz, log_n_nz, an_nz)
        L1s = Lambda_AFE(1-s, an_nz, log_n_nz, an_nz)
        denom = max(abs(Ls), abs(L1s), 1e-300)
        rel = abs(Ls - EPSILON*L1s) / denom
        fe_errs.append(rel)
        ok = "✅" if rel < 1e-6 else "❌"
        log(f"  s={s}: |Λ(s)|={abs(Ls):.4e}, |Λ(1-s)|={abs(L1s):.4e}, rel={rel:.2e} {ok}")
    except Exception as e:
        log(f"  s={s}: ERROR {e}")
        fe_errs.append(1.0)
log(f"  소요: {time.time()-t1:.1f}s")

if not fe_errs or min(fe_errs) > 1e-3:
    log("  ⚠️ FE 검증 실패 — 감마 인자 점검 필요")
    log("  → 대안1 (μ=0,1,11,12) 시도 — 종료하고 v2 필요")
    flush()
    sys.exit(1)

fe_pass = sum(1 for e in fe_errs if e < 1e-6)
max_rel = max(fe_errs) if fe_errs else 999
log(f"\n  FE pass: {fe_pass}/{len(fe_errs)}, max_rel={max_rel:.2e}")
log(f"  {'✅ FE 통과' if fe_pass >= 5 else '❌ FE 실패'}")
flush()

if fe_pass < 3:
    log("  ⚠️ FE 3개 미만 통과 — 결과 신뢰도 낮음. 계속 진행하나 주의.")

# ── Step 2: 영점 탐색 ──
log(f"\n[Step 2] Λ(1/2+it) 부호 변환 탐색 (t ∈ [{T_MIN},{T_MAX}], Δt=0.3)")
scan_ts = np.arange(T_MIN, T_MAX + 0.1, 0.3)
scan_re = np.zeros(len(scan_ts))
t2 = time.time()

for i, t in enumerate(scan_ts):
    try:
        val = Lambda_AFE(complex(0.5, t), an_nz, log_n_nz, an_nz)
        scan_re[i] = val.real
    except Exception as e:
        scan_re[i] = 0.0
    if i % 30 == 0 or i == len(scan_ts)-1:
        elapsed = time.time() - t2
        eta = elapsed/(i+1)*(len(scan_ts)-i-1) if i > 0 else 0
        log(f"  t={t:5.1f}: Λ={scan_re[i]:+.6e}  ({i+1}/{len(scan_ts)}, eta={eta:.0f}s)")
        flush()

n_sign = int(np.sum(np.diff(np.sign(scan_re)) != 0))
log(f"\n  부호 변환: {n_sign}회, 소요: {time.time()-t2:.1f}s")
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
                v_mid = Lambda_AFE(complex(0.5, mid), an_nz, log_n_nz, an_nz).real
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
        k_val = curvature(sigma_near, z_t, an_nz, log_n_nz, an_nz)
    except Exception as e:
        k_val = float('nan')
        log(f"  WARNING κ 계산 오류 t={z_t:.4f}: {e}")

    try:
        m_val = monodromy(z_t, an_nz, log_n_nz, an_nz)
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
    log(f"\n  ⚠️ κ_near 계산 실패")

# ── Step 4: κ near/far 비율 ──
log(f"\n[Step 4] κ near/far 비율 (ratio = κ 집중도)")
far_ts = np.arange(T_MIN, T_MAX, 1.5)
k_near_list, k_far_list = [], []
for t in far_ts:
    try:
        k = curvature(sigma_near, t, an_nz, log_n_nz, an_nz)
    except Exception:
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
    log(f"  ratio= {ratio:.1f}×  {'✅' if ratio >= 10 else '❌'}")
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
            v = Lambda_AFE(complex(sig, t), an_nz, log_n_nz, an_nz)
            vals.append(v.real)
        except Exception:
            vals.append(0.0)
    sc = int(np.sum(np.diff(np.sign(vals)) != 0))
    sigma_sc[sig] = sc
    log(f"  σ={sig:.1f}: 부호변환={sc}")

sc_half = sigma_sc.get(0.5, 0)
is_uniq = all(sc_half >= sigma_sc[s] for s in sigma_vals if s != 0.5)
log(f"\n  σ=0.5 최대: {'✅ PASS (N=1 패턴 확인)' if is_uniq else '❌ FAIL'}")
flush()

# ── 종합 ──
log(f"\n{'='*70}")
log(f"종합 판정 — GL(4) Rankin-Selberg L(s, Δ×Δ)")
log(f"{'='*70}")
fe_pass_final = sum(1 for e in fe_errs if e < 1e-6)
log(f"  [P1] 함수방정식: {fe_pass_final}/{len(fe_errs)}, max_rel={max_rel:.2e}")
log(f"       {'✅ PASS' if fe_pass_final >= 5 else '⚠️ 부분'}")
log(f"  [P2] 영점 발견: {len(zeros)}개 (t∈[{T_MIN},{T_MAX}])")
log(f"       {'✅ PASS' if len(zeros) >= 5 else '⚠️ 부족'}")
if np.isfinite(k_mean):
    log(f"  [P3] κ_near(d=4): {k_mean:.2f} (CV={k_cv:.1f}%, n={len(kappa_vals)})")
    log(f"       {'✅ PASS' if k_cv < 2 else '⚠️ 산포 큼'}")
log(f"  [P4] 모노드로미: {n_tp_mono}/{len(zeros)} TP")
log(f"       {'✅ PASS' if n_tp_mono == len(zeros) else '⚠️ 일부 실패'}")
log(f"  [B-10] σ-유일성(N=1): {'✅ PASS' if is_uniq else '❌ FAIL (예상 밖)'}")
log(f"  κ ratio: {ratio:.1f}×")

log(f"\n  ━━ κ_near(d) 비교 ━━")
log(f"  d=1 (ζ):         κ_near = 1112.32")
log(f"  d=2 (GL2 avg):   κ_near = 1114.6")
log(f"  d=3 (GL3 sym²):  κ_near = 1125.16")
if np.isfinite(k_mean):
    log(f"  d=4 (Δ×Δ):      κ_near = {k_mean:.2f}  ← 이번 측정")
    gap_34 = k_mean - 1125.16
    gap_23 = 1125.16 - 1114.6
    gap_12 = 1114.6 - 1112.32
    log(f"  gap(3→4) = {gap_34:.2f}")
    log(f"  gap(2→3) = {gap_23:.2f}")
    log(f"  gap(1→2) = {gap_12:.2f}")
    if gap_34 > 0:
        log(f"  → 단조증가 4점 확인: ✅")
        accel = gap_34 / gap_23 if gap_23 > 0 else 0
        log(f"  → 가속비 (3→4)/(2→3) = {accel:.2f}×")
    else:
        log(f"  → ⚠️ 단조증가 위반 또는 노이즈")

log(f"\n  ━━ κ ratio degree 비교 ━━")
log(f"  d=1: 2200.7×")
log(f"  d=2: ~600×")
log(f"  d=3: ~320×")
log(f"  d=4: {ratio:.1f}×")

total = time.time() - t_total
log(f"\n  총 소요: {total:.1f}초 ({total/60:.1f}분)")
log(f"  시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
flush()
log(f"\n결과 저장: {OUTFILE}")
flush()
