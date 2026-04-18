#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #90 — A(t₀) Hadamard 분해: 정리 도출 + 수치 검증
=============================================================================
이론 (Proposition):
  ξ(s) = ξ(0)·Π_{n≥1}[(1-s/ρ_n)(1-s/ρ̄_n)]  (D=0: ξ'/ξ(1/2)=0 by symmetry)

  ξ'/ξ(s) = Σ_{n≥1}[1/(s-ρ_n) + 1/(s-ρ̄_n)]

  Near ρ₀=1/2+it₀: isolate 1/(s-ρ₀) term:
    g(s) = 1/(s-ρ̄₀) + Σ_{n≠0}[1/(s-ρ_n)+1/(s-ρ̄_n)]
    c₀ = g(ρ₀): all terms purely imaginary → Re(c₀)=0 ✓

  Laurent coefficients (ρ̄₀=1/2-it₀ contributes explicitly):
    B  = Im(c₀) = -1/(2t₀) + Σ_{n≠0}[-(1/(t₀-t_n)+1/(t₀+t_n))]
    H₁ = Re(c₁) = +1/(4t₀²) + Σ_{n≠0}[1/(t₀-t_n)²+1/(t₀+t_n)²]

  A(t₀) = B² + 2H₁

검증:
  방법 A: A_direct = κ(ρ₀+δ) - 1/δ²   at δ=0.01
  방법 B: Richardson 외삽 → B_rich, H₁_rich
  방법 C: B_had² + 2H₁_had  (Hadamard paired sum, N=100..5000)

결과: results/hadamard_A_decomposition_90.txt
=============================================================================
"""

import sys, os, time
import numpy as np

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))
import mpmath

DPS = 60
mpmath.mp.dps = DPS

from bundle_utils import xi_func, find_zeros_zeta

CACHE_PATH = os.path.expanduser(
    '~/Desktop/gdl_unified/outputs/cache/zeta_zeros_t_N5000.npy')
OUTFILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/hadamard_A_decomposition_90.txt')
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))

def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━━━━━━ Euler-Maclaurin 꼬리 보정 ━━━━━━━━━━━
try:
    from scipy import integrate as _sci_int

    def tail_correction(t0_val, t_last, mode='B'):
        """
        N번째 영점 이후의 꼬리 합을 적분으로 근사.
        zero density rho(t) = log(t/(2pi)) / (2pi)

        B_tail  = integral_{t_last}^{inf} 2*t0/(t^2-t0^2) * rho(t) dt
        H1_tail = integral_{t_last}^{inf} 2/t^2 * rho(t) dt
        """
        if mode == 'B':
            fn = lambda t: 2*t0_val/(t**2-t0_val**2) * np.log(t/(2*np.pi))/(2*np.pi)
        else:
            fn = lambda t: 2.0/t**2 * np.log(t/(2*np.pi))/(2*np.pi)
        val, _ = _sci_int.quad(fn, t_last, 1e7, limit=300)
        return val

    HAS_TAIL = True
except ImportError:
    HAS_TAIL = False
    def tail_correction(t0_val, t_last, mode='B'): return 0.0

# ━━━━━━━━━━━ 파라미터 ━━━━━━━━━━━
T_MIN, T_MAX = 10.0, 60.0
DELTA_DIRECT = 0.01
N_HAD_LIST = [100, 500, 1000, 5000]

log("=" * 72)
log("결과 #90 — A(t₀) Hadamard 분해: 수치 검증")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"DPS={DPS}, δ_direct={DELTA_DIRECT}, N_Had={N_HAD_LIST}")
log()

# ━━━━━━━━━━━ xi'/xi 수치 함수 ━━━━━━━━━━━
def conn_at(s):
    """ξ'/ξ(s) via numerical differentiation"""
    h = mpmath.mpf(10)**(-20)
    xi_v = xi_func(s)
    if abs(xi_v) < mpmath.mpf(10)**(-DPS+10):
        return mpmath.mpc(1e10, 0)
    xi_d = (xi_func(s+h) - xi_func(s-h)) / (2*h)
    return xi_d / xi_v

# ━━━━━━━━━━━ Step 1: 영점 탐색 ━━━━━━━━━━━
log("[Step 1] ζ(s) 영점 탐색 (t ∈ [10, 60])")
t0 = time.time()
zeros_local = find_zeros_zeta(T_MIN, T_MAX)
log(f"  {len(zeros_local)}개 발견 ({time.time()-t0:.1f}s)")
if len(zeros_local) == 0:
    log("⚠️ 영점 0개 — 탐색 로직 점검 필요"); flush_file(); sys.exit(1)
for j, z in enumerate(zeros_local):
    log(f"    #{j+1:>2}: t = {z:.10f}")
log()
flush_file()

# ━━━━━━━━━━━ Step 2: zetazero 로드 (캐시 or 즉시) ━━━━━━━━━━━
log(f"[Step 2] 영점 데이터 로드 (N≤5000)")
N_MAX_AVAIL = 0
if os.path.exists(CACHE_PATH):
    all_t = np.load(CACHE_PATH)
    N_MAX_AVAIL = len(all_t)
    log(f"  캐시 로드: {N_MAX_AVAIL}개 ({CACHE_PATH})")
else:
    log(f"  캐시 없음 — 즉시 계산 (N=2000, DPS=30)")
    mpmath.mp.dps = 30
    t_start = time.time()
    t_vals = []
    for n in range(1, 2001):
        t_vals.append(float(mpmath.zetazero(n).imag))
        if n % 500 == 0:
            log(f"  {n}/2000 ({time.time()-t_start:.0f}s)")
            flush_file()
    all_t = np.array(t_vals)
    N_MAX_AVAIL = len(all_t)
    np.save(CACHE_PATH, all_t)
    mpmath.mp.dps = DPS
    log(f"  계산 완료: {N_MAX_AVAIL}개 ({time.time()-t_start:.0f}s)")

# N_HAD_LIST를 가용 범위로 제한
N_HAD_LIST = [N for N in N_HAD_LIST if N <= N_MAX_AVAIL]
if N_MAX_AVAIL >= 5000:
    if 5000 not in N_HAD_LIST: N_HAD_LIST.append(5000)
else:
    N_HAD_LIST.append(N_MAX_AVAIL)
N_HAD_LIST = sorted(set(N_HAD_LIST))
log(f"  실제 사용 N 목록: {N_HAD_LIST}")
log()
flush_file()

# ━━━━━━━━━━━ Hadamard sum 함수 ━━━━━━━━━━━
def hadamard_B_H1(t0_val, t_vals_N, include_rhobar0=True):
    """
    B = Im(c₀), H₁ = Re(c₁) from Hadamard paired sum.

    ρ₀ = 1/2+it₀,  ρ̄₀ = 1/2-it₀
    Other zeros: ρ_n=1/2+it_n, ρ̄_n=1/2-it_n  (n=1..N, t_n≠t₀)

    c₀ contributions:
      From ρ̄₀: 1/(ρ₀-ρ̄₀) = 1/(2it₀) = -i/(2t₀)  →  Im = -1/(2t₀)
      From (ρ_n,ρ̄_n): Im = -(1/(t₀-t_n)+1/(t₀+t_n))

    c₁ contributions:
      From ρ̄₀: -1/(ρ₀-ρ̄₀)² = 1/(4t₀²)  [real]
      From (ρ_n,ρ̄_n): 1/(t₀-t_n)²+1/(t₀+t_n)²  [real]
    """
    mask = np.abs(t_vals_N - t0_val) > 1e-6
    tv = t_vals_N[mask]

    B_others = -np.sum(1.0/(t0_val - tv) + 1.0/(t0_val + tv))
    H1_others = np.sum(1.0/(t0_val - tv)**2 + 1.0/(t0_val + tv)**2)

    if include_rhobar0:
        B_rho0bar = -1.0 / (2 * t0_val)
        H1_rho0bar = 1.0 / (4 * t0_val**2)
    else:
        B_rho0bar = 0.0
        H1_rho0bar = 0.0

    B = B_rho0bar + B_others
    H1 = H1_rho0bar + H1_others
    A = B**2 + 2*H1
    return B, H1, A

# ━━━━━━━━━━━ Step 3: 각 영점별 계산 ━━━━━━━━━━━
log("[Step 3] 영점별 A(t₀) 계산")
log()

results = []

for idx, t0_val in enumerate(zeros_local):
    log(f"  ── 영점 #{idx+1}: t₀ = {t0_val:.8f} ──")
    rho0 = mpmath.mpc(mpmath.mpf('0.5'), mpmath.mpf(str(t0_val)))

    # ─── 방법 A: A_direct ───
    delta = mpmath.mpf(str(DELTA_DIRECT))
    s_near = rho0 + delta
    h_num = mpmath.mpf(10)**(-20)
    xi_near = xi_func(s_near)
    xi_d = (xi_func(s_near+h_num) - xi_func(s_near-h_num))/(2*h_num)
    conn_near = xi_d / xi_near
    kappa_near = float(abs(conn_near)**2)
    A_direct = kappa_near - 1.0/DELTA_DIRECT**2
    log(f"    A_direct (δ={DELTA_DIRECT}) = {A_direct:.6f}")

    # ─── 방법 B: Richardson c₀, c₁ ───
    eps_list = [1e-3, 5e-4, 2.5e-4, 1.25e-4]
    c0_raw = []
    for eps_f in eps_list:
        eps = mpmath.mpf(str(eps_f))
        c0_raw.append(conn_at(rho0 + eps) - 1/eps)
    # Richardson O(ε²) 소거 2회
    r1 = [(4*c0_raw[i+1]-c0_raw[i])/3 for i in range(3)]
    r2 = [(4*r1[i+1]-r1[i])/3 for i in range(2)]
    c0_rich = (4*r2[1]-r2[0])/3

    # c₁: [conn(ρ₀+ε) - 1/ε - c₀] / ε at smallest ε
    eps_min = mpmath.mpf(str(eps_list[-1]))
    c1_rich = (conn_at(rho0+eps_min) - 1/eps_min - c0_rich) / eps_min

    B_rich = float(c0_rich.imag)
    H1_rich = float(c1_rich.real)
    A_rich = B_rich**2 + 2*H1_rich
    log(f"    Richardson: B={B_rich:.6f}, H₁={H1_rich:.6f}, A_rich={A_rich:.6f}")

    # ─── 방법 C: Hadamard paired sum ───
    A_had = {}
    for N in N_HAD_LIST:
        t_N = all_t[:N]
        B_n, H1_n, A_n = hadamard_B_H1(t0_val, t_N, include_rhobar0=True)
        A_had[N] = (B_n, H1_n, A_n)

    Nmax = N_HAD_LIST[-1]
    log(f"    Hadamard: " + ", ".join([f"N={N}→{A_had[N][2]:.4f}" for N in N_HAD_LIST]))
    rel_err = abs(A_had[Nmax][2] - A_direct) / (abs(A_direct) + 1e-12)
    log(f"    상대오차 (N={Nmax}) = {rel_err*100:.3f}%  ({'✓ <1%' if rel_err<0.01 else '✗ ≥1%'})")

    # ─── Euler-Maclaurin 꼬리 보정 ───
    t_last = all_t[Nmax-1]  # last zero used
    B_Nmax, H1_Nmax, _ = A_had[Nmax]
    B_tail = tail_correction(t0_val, t_last, 'B')
    H1_tail = tail_correction(t0_val, t_last, 'H1')
    B_corr = B_Nmax + B_tail
    H1_corr = H1_Nmax + H1_tail
    A_corr = B_corr**2 + 2*H1_corr
    rel_err_corr = abs(A_corr - A_direct) / (abs(A_direct) + 1e-12)
    log(f"    EM보정: B_tail={B_tail:.5f}, H₁_tail={H1_tail:.5f}")
    log(f"    A_corrected={A_corr:.6f}  상대오차={rel_err_corr*100:.3f}%  ({'✓' if rel_err_corr<0.01 else '✗'})")

    results.append({
        't': t0_val, 'A_direct': A_direct,
        'B_rich': B_rich, 'H1_rich': H1_rich, 'A_rich': A_rich,
        'A_had': A_had, 'rel_err': rel_err, 'Nmax': Nmax,
        'A_corr': A_corr, 'rel_err_corr': rel_err_corr,
    })
    log()
    flush_file()

# ━━━━━━━━━━━ Step 4: 요약 표 ━━━━━━━━━━━
log("=" * 90)
log(f"요약 표: A(t₀) 직접 vs Hadamard (N={N_HAD_LIST[-1]})")
log("=" * 90)
Nmax = N_HAD_LIST[-1]
log(f"{'#':>3} {'t₀':>12} {'A_direct':>10} {'A_Had(N)':>12} {'A_corrected':>13} {'err_Had%':>10} {'err_corr%':>10} pass")
log("-" * 100)
n_pass_had = 0
n_pass_corr = 0
for i, r in enumerate(results):
    A_Hn = r['A_had'][Nmax][2]
    ep = r['rel_err']*100
    ec = r['rel_err_corr']*100
    ok_had = ep < 1.0
    ok_corr = ec < 1.0
    if ok_had: n_pass_had += 1
    if ok_corr: n_pass_corr += 1
    log(f"  {i+1:>2} {r['t']:>12.6f} {r['A_direct']:>10.4f} {A_Hn:>12.4f} "
        f"{r['A_corr']:>13.4f} {ep:>10.3f} {ec:>10.3f}  {'✓' if ok_had else '✗'}/{'✓' if ok_corr else '✗'}")
log()
log(f"  Hadamard(N={Nmax}) < 1%: {n_pass_had}/{len(results)}")
log(f"  EM보정 후 < 1%:         {n_pass_corr}/{len(results)}  (성공 기준: 10/13 이상)")
log()
flush_file()

# ━━━━━━━━━━━ Step 5: N 수렴 분석 ━━━━━━━━━━━
log("=" * 72)
log("N 수렴 분석 (t₁, t₉, t₁₃ 대표 영점)")
log("=" * 72)
for idx in [0, 8, 12]:
    if idx >= len(results): continue
    r = results[idx]
    log(f"  t₀ = {r['t']:.6f}  (A_direct={r['A_direct']:.4f})")
    log(f"  {'N':>6} {'B':>10} {'H₁':>12} {'A_Had':>10} {'Δ/A%':>8}")
    for N in N_HAD_LIST:
        if N not in r['A_had']: continue
        B_n, H1_n, A_n = r['A_had'][N]
        dp = abs(A_n - r['A_direct'])/(abs(r['A_direct'])+1e-12)*100
        log(f"  {N:>6} {B_n:>10.5f} {H1_n:>12.5f} {A_n:>10.4f} {dp:>7.3f}%")
    log()

# ━━━━━━━━━━━ Step 6: B-20 해석 ━━━━━━━━━━━
log("=" * 72)
log("B-20 해석: t₉=48.005 높은 A — 근접 영점 기여")
log("=" * 72)
t9_idx = 8
if t9_idx < len(results):
    r9 = results[t9_idx]
    t9 = r9['t']
    log(f"  t₉={t9:.6f}, A_direct={r9['A_direct']:.4f}")

    tv5000 = all_t[:5000] if N_MAX_AVAIL>=5000 else all_t
    mask_all = np.abs(tv5000 - t9) > 1e-6
    tv = tv5000[mask_all]
    mask_near = np.abs(tv - t9) < 5.0

    # 근접 영점 나열
    tv_near = tv[mask_near]
    log(f"  근접 영점 (|t₀-t_n|<5):")
    for tn in tv_near:
        H1_contrib = 1.0/(t9-tn)**2 + 1.0/(t9+tn)**2
        B_contrib  = -(1.0/(t9-tn) + 1.0/(t9+tn))
        log(f"    t_n={tn:.4f}, Δt={t9-tn:.4f}, B_c={B_contrib:.4f}, H₁_c={H1_contrib:.4f}")

    # 근접 vs 원격 H₁ 기여
    H1_near = np.sum(1.0/(t9-tv[mask_near])**2 + 1.0/(t9+tv[mask_near])**2)
    H1_far  = np.sum(1.0/(t9-tv[~mask_near])**2 + 1.0/(t9+tv[~mask_near])**2)
    B_near  = -np.sum(1.0/(t9-tv[mask_near]) + 1.0/(t9+tv[mask_near]))
    log(f"  H₁(근접<5): {H1_near:.4f}, H₁(원격): {H1_far:.4f}")

    # 비교: t₁ (낮은 A)
    r1 = results[0]
    t1 = r1['t']
    tv1 = tv5000[np.abs(tv5000-t1)>1e-6]
    m1_near = np.abs(tv1 - t1) < 5.0
    H1_1near = np.sum(1.0/(t1-tv1[m1_near])**2 + 1.0/(t1+tv1[m1_near])**2)
    log(f"  비교 t₁={t1:.4f}: H₁_near={H1_1near:.4f}, vs t₉ H₁_near={H1_near:.4f}")
    if H1_near > H1_1near:
        log(f"  → t₉의 높은 A는 근접 영점 H₁ 기여({H1_near:.3f}) > t₁({H1_1near:.3f})로 설명 가능")
    else:
        log(f"  → 근접 영점만으로 설명 불충분 (B² 항도 확인 필요)")

log()
log("=" * 72)
log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"결과 파일: {OUTFILE}")
log("=" * 72)
flush_file()
