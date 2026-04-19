#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #102 — Hadamard A=B²+2H₁ GL(4) sym³(Δ) ξ-bundle 검증
=============================================================================
★ 핵심 발견:
  1. #80 gammaV=[-1,0,0,1], N=144 → FE=-1 (잘못된 파라미터, 가짜 영점)
  2. 올바른 파라미터: gammaV=[5.5,6.5,16.5,17.5], N=1, w=0 → FE=-51
     Γ_C(s+(k-1)/2)·Γ_C(s+3(k-1)/2) = Γ_C(s+5.5)·Γ_C(s+16.5)
  3. ξ-bundle κ = |Λ'/Λ|² = |L'/L + ψ_gamma|²
     lfunlambda는 불안정 → lfun + 수동 감마 보정 방식 사용
  4. lfuninit는 대형 감마이동에서 오작동 → lfun(Ldata, s) 직접 사용
  5. lfunzeros도 오작동 → |L| 최솟값 탐색으로 수동 영점 발견

수학적 접근:
  κ(σ₀+δ, t₀) = |conn_L(s) + conn_Γ(s)|²   [s = σ₀+δ+it₀]
  where conn_L = L'/L(s), conn_Γ = (1/2)Σψ((s+μ_j)/2) - (d/2)log(π)
  A_direct = κ - 1/δ²
  Re(c₀) Richardson 보정: A_true = A_direct - 2·Re(c₀)/δ
  Hadamard: A_corr = (B+B_tail)² + 2(H₁+H₁_tail)

결과: results/hadamard_gl4_sym3_102.txt
=============================================================================
"""

import sys, os, time, math
import numpy as np

OUTFILE = os.path.expanduser(
    '~/Desktop/gdl_unified/results/hadamard_gl4_sym3_102.txt')
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

lines = []
def log(msg=""):
    print(str(msg), flush=True)
    lines.append(str(msg))
def flush_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

log("=" * 72)
log("결과 #102 — Hadamard A=B²+2H₁ GL(4) sym³(Δ) ξ-bundle 검증")
log("=" * 72)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log()
log("★ 핵심: #80 gammaV 수정 + lfun+ψ (NO lfuninit) 방식으로 ξ-bundle κ 계산")
log()
flush_file()

# ━━━━━━━━━━━ PARI 초기화 ━━━━━━━━━━━
import cypari2
gp = cypari2.Pari()
gp.allocatemem(6000000000)
gp("default(realprecision, 38)")  # 15자리 계수 호환
log("[Init] PARI: 6GB, prec=38 ✅")
flush_file()

# ━━━━━━━━━━━ 파라미터 ━━━━━━━━━━━
CENTER       = 0.5
DELTA_DIRECT = 0.01
N_TEST       = 10
N_COEFF      = 30000
GAMMA_V      = [5.5, 6.5, 16.5, 17.5]  # FE=-51
T_SCAN_MAX   = 25.0   # lfun(Ldata)는 t≤25에서 신뢰할 수 있음

log(f"파라미터: center={CENTER}, δ={DELTA_DIRECT}, N_COEFF={N_COEFF}")
log(f"  gammaV={GAMMA_V}, N=1, w=0")
log(f"  영점 탐색 범위: t ∈ [0, {T_SCAN_MAX}]")
log(f"  방법: lfun(Ldata, s) 직접 사용 (lfuninit 사용 안 함)")
log()
flush_file()

# ━━━━━━━━━━━ EM 꼬리 보정 (d=4) ━━━━━━━━━━━
try:
    from scipy import integrate as _sci_int
    def tail_gl4(t0_val, t_last, mode):
        c = 4.0 / (2.0 * math.pi)
        c2pi = 2.0 * math.pi
        if mode == 'B':
            fn = lambda t: (2.0*t0_val/(t**2 - t0_val**2) *
                            c * math.log(t/c2pi) if t > t0_val+1e-6 else 0.0)
        else:
            fn = lambda t: 2.0/t**2 * c * math.log(t/c2pi)
        v, _ = _sci_int.quad(fn, t_last, 1e7, limit=500, epsabs=1e-12)
        return v
    HAS_TAIL = True
    log("[EM tail] scipy ✅")
except ImportError:
    HAS_TAIL = False
    def tail_gl4(*a): return 0.0
    log("[EM tail] scipy 없음")
flush_file()

# ━━━━━━━━━━━ Step 1: sym³(Δ) 계수 ━━━━━━━━━━━
log("=" * 72)
log(f"[Step 1] sym³(Δ) 계수 ({N_COEFF}개)")
log("=" * 72)
t0_s = time.time()

tau = [0] * (N_COEFF + 1)
for n in range(1, N_COEFF + 1):
    tau[n] = int(gp.ramanujantau(n))
    if n % 10000 == 0: log(f"  τ: n={n}/{N_COEFF}"); flush_file()

sieve = [True] * (N_COEFF + 1); sieve[0] = sieve[1] = False
for i in range(2, int(N_COEFF**0.5) + 1):
    if sieve[i]:
        for j in range(i*i, N_COEFF+1, i): sieve[j] = False
primes = [i for i in range(2, N_COEFF+1) if sieve[i]]

cpk = {}
for p in primes:
    tp = tau[p] / (p ** 5.5)
    e1 = tp**3 - 2.0*tp
    e2 = (tp**2-2.0)*(tp**2-1.0)
    cpk[(p,0)] = 1.0; cpk[(p,1)] = e1
    pk = p; k = 1
    while pk * p <= N_COEFF:
        pk *= p; k += 1
        cpk[(p,k)] = (e1*cpk.get((p,k-1),0) - e2*cpk.get((p,k-2),0)
                      + e1*cpk.get((p,k-3),0) - cpk.get((p,k-4),0))

cn = np.zeros(N_COEFF + 1); cn[1] = 1.0
for n in range(2, N_COEFF + 1):
    temp = n; result = 1.0
    for p in primes:
        if p*p > temp: break
        if temp % p == 0:
            k = 0
            while temp % p == 0: k += 1; temp //= p
            result *= cpk.get((p, k), 0.0)
    if temp > 1: result *= cpk.get((temp, 1), 0.0)
    cn[n] = result

log(f"  c(2)={cn[2]:.6f}, c(3)={cn[3]:.6f}")
log(f"  소요: {time.time()-t0_s:.1f}s")
flush_file()

# ━━━━━━━━━━━ Step 2: PARI L-함수 생성 (NO lfuninit) ━━━━━━━━━━━
log("=" * 72)
log("[Step 2] PARI L-함수 (NO lfuninit — 대형 감마이동 호환)")
log("=" * 72)
t0_s = time.time()

cn_str = "[" + ",".join(f"{cn[i]:.15e}" for i in range(1, N_COEFF+1)) + "]"
gp(f"global(Ldata); Ldata = lfuncreate([{cn_str}, 0, [5.5,6.5,16.5,17.5], 1, 1, 0])")

fe_score = float(gp("lfuncheckfeq(Ldata)"))
log(f"  FE = {fe_score:.1f}  {'✅' if fe_score <= -8 else '⚠️'}")

# L(1/2) 확인 — root number
L_half = complex(gp("lfun(Ldata, 0.5)"))
log(f"  L(1/2) = {L_half.real:.4e} + {L_half.imag:.4e}i")
if abs(L_half) < 1e-10:
    log(f"  → L(1/2) ≈ 0: root number = -1 (odd order vanishing)")
    has_central_zero = True
else:
    has_central_zero = False
log(f"  소요: {time.time()-t0_s:.1f}s")
flush_file()

# ━━━━━━━━━━━ Step 2b: 수동 영점 탐색 ━━━━━━━━━━━
log()
log("=" * 72)
log("[Step 2b] 수동 영점 탐색 (lfunzeros 대신 |L| 최솟값 탐색)")
log("=" * 72)
t0_s = time.time()

def eval_absL(t):
    """L(1/2+it) 의 절댓값"""
    try:
        v = complex(gp(f"lfun(Ldata, 0.5 + I*{t:.12f})"))
        return abs(v)
    except Exception:
        return 1e30

def refine_zero_golden(t_lo, t_hi, tol=1e-9, max_iter=80):
    """골든 섹션 서치로 |L| 최솟값 위치 정밀 결정"""
    phi = (1 + 5**0.5)/2
    for _ in range(max_iter):
        if t_hi - t_lo < tol:
            break
        t1 = t_hi - (t_hi - t_lo)/phi
        t2 = t_lo + (t_hi - t_lo)/phi
        if eval_absL(t1) < eval_absL(t2):
            t_hi = t2
        else:
            t_lo = t1
    t_mid = (t_lo + t_hi)/2
    return t_mid, eval_absL(t_mid)

# 1단계: 거친 스캔 (dt=0.05)
log("  [1] 거친 스캔 (dt=0.05)...")
flush_file()
dt = 0.05
t_scan = np.arange(0.2, T_SCAN_MAX + dt, dt)
absL_scan = []
for t in t_scan:
    absL_scan.append(eval_absL(t))
absL_scan = np.array(absL_scan)

# 극솟값 후보 수집 (|L| < 0.5)
candidate_intervals = []
for i in range(1, len(absL_scan)-1):
    if absL_scan[i] < absL_scan[i-1] and absL_scan[i] < absL_scan[i+1] and absL_scan[i] < 0.5:
        candidate_intervals.append((t_scan[max(0,i-1)], t_scan[min(len(t_scan)-1,i+1)]))

log(f"    극솟값 후보: {len(candidate_intervals)}개")

# 2단계: 정밀 보정
log("  [2] 정밀 보정 (골든 섹션)...")
flush_file()
confirmed_zeros = []
for lo, hi in candidate_intervals:
    t0, absL0 = refine_zero_golden(lo, hi)
    if absL0 < 0.1:  # |L| < 0.1이면 진짜 영점
        confirmed_zeros.append(t0)
        log(f"    ✓ t={t0:.8f}, |L|={absL0:.4e}")

# 중앙 영점 포함 (root number=-1인 경우)
if has_central_zero:
    confirmed_zeros = [0.0] + confirmed_zeros
    log(f"    ✓ t=0.00000000, |L|≈0 (central zero)")

all_zeros = np.array(sorted(confirmed_zeros))
N_ALL = len(all_zeros)
t_last = float(all_zeros[-1]) if N_ALL > 0 else 0.0

log()
log(f"  확인 영점: N={N_ALL}개")
if N_ALL > 0:
    log(f"    t ∈ [{all_zeros[0]:.6f}, {t_last:.6f}]")
    for i, z in enumerate(all_zeros):
        log(f"    γ_{i} = {z:.8f}")
    if N_ALL > 2:
        log(f"    비중앙 영점 간격 평균: {np.mean(np.diff(all_zeros[1:] if has_central_zero else all_zeros)):.4f}")
log(f"  소요: {time.time()-t0_s:.1f}s")
log()
flush_file()

if N_ALL < 3:
    log("❌ 영점 부족 — 중단"); flush_file(); sys.exit(1)

# ━━━━━━━━━━━ Step 3: ξ-bundle κ 검증 ━━━━━━━━━━━
log("=" * 72)
log("[Step 3] ξ-bundle κ 방법론")
log("=" * 72)
log("  κ(s) = |L'/L(s) + conn_Γ(s)|²")
log("  conn_Γ = -2·log(π) + (1/2)Σ_j ψ((s+μ_j)/2)")
log("  μ = [5.5, 6.5, 16.5, 17.5]")
log("  ★ lfun(Ldata, s) 직접 사용 (lfuninit 없음)")
log()

# 비중앙 영점에서 κ 검증
nz_zeros = all_zeros[all_zeros > 0.5]  # 비중앙 영점만
if len(nz_zeros) == 0:
    log("❌ 비중앙 영점 없음 — 중단"); flush_file(); sys.exit(1)

t0_test = nz_zeros[0]
gp(f"s_v = {CENTER + DELTA_DIRECT:.10f} + I*{t0_test:.10f}")
gp("Lv = lfun(Ldata, s_v)")
gp("Lpv = lfun(Ldata, s_v, 1)")
abs_L_test = float(gp("abs(Lv)"))
gp("conn_L_v = Lpv/Lv")
gp("conn_G = -2*log(Pi)+0.5*psi((s_v+5.5)/2)+0.5*psi((s_v+6.5)/2)"
   "+0.5*psi((s_v+16.5)/2)+0.5*psi((s_v+17.5)/2)")
gp("conn_total_v = conn_L_v + conn_G")
kappa_v = float(gp("abs(conn_total_v)^2"))
A_v = kappa_v - 1.0/DELTA_DIRECT**2
log(f"  검증 t₀={t0_test:.4f}: |L|={abs_L_test:.4e}, κ={kappa_v:.2f}, A={A_v:.4f}")
if kappa_v > 5000:
    log("  ✅ ξ-bundle (lfun+ψ, NO lfuninit) 작동 확인")
else:
    log("  ❌ κ 비정상")
    flush_file(); sys.exit(1)
log()
flush_file()

# ━━━━━━━━━━━ 함수 정의 ━━━━━━━━━━━
def hadamard_B_H1(t0, tzeros):
    """Hadamard B와 H₁ 계산 — 중앙 영점(t=0) 포함"""
    mask = np.abs(tzeros - t0) > 1e-6
    tv = tzeros[mask]
    if t0 < 1e-6:
        # t₀ ≈ 0이면 skip (테스트하지 않음)
        return 0.0, 0.0, 0.0
    d = t0 - tv; s = t0 + tv
    ok = (np.abs(d) > 1e-12) & (np.abs(s) > 1e-12)
    d, s = d[ok], s[ok]
    # B = -1/(2t₀) + Σ(-1/(t₀-γ_n) - 1/(t₀+γ_n))
    # 중앙 영점(γ=0)이 포함되면 자동으로 -1/(t₀-0) - 1/(t₀+0) = -2/t₀
    B  = -1.0/(2.0*t0) + np.sum(-1.0/d - 1.0/s)
    H1 = 1.0/(4.0*t0**2) + np.sum(1.0/d**2 + 1.0/s**2)
    return B, H1, B**2 + 2.0*H1

def compute_kappa_lfun_psi(t0_val, delta=DELTA_DIRECT):
    """ξ-bundle κ via lfun(Ldata) + 수동 감마 보정 — NO lfuninit"""
    try:
        gp(f"s_k = {CENTER+delta:.10f} + I*{t0_val:.10f}")
        gp("Lk = lfun(Ldata, s_k)")
        gp("Lpk = lfun(Ldata, s_k, 1)")
        abs_Lk = float(gp("abs(Lk)"))
        if abs_Lk < 1e-30:
            return float('nan'), float('nan'), False
        gp("conn_Lk = Lpk/Lk")
        gp("conn_Gk = -2*log(Pi)+0.5*psi((s_k+5.5)/2)+0.5*psi((s_k+6.5)/2)"
           "+0.5*psi((s_k+16.5)/2)+0.5*psi((s_k+17.5)/2)")
        gp("conn_tk = conn_Lk + conn_Gk")
        kappa = float(gp("abs(conn_tk)^2"))
        A_dir = kappa - 1.0/delta**2
        return kappa, A_dir, True
    except Exception as e:
        log(f"    WARNING κ: {e}")
        return float('nan'), float('nan'), False

def compute_Re_c0_richardson(t0_val, delta=DELTA_DIRECT):
    """Richardson Re(c₀) 외삽 — lfun(Ldata)+ψ 방식, NO lfuninit"""
    _, A_dir, ok = compute_kappa_lfun_psi(t0_val, delta)
    if not ok:
        return float('nan'), float('nan'), False

    eps_list = [1e-3, 5e-4, 2.5e-4, 1.25e-4]
    c0_raw = []
    for eps_f in eps_list:
        try:
            gp(f"s_r = {CENTER+eps_f:.10f} + I*{t0_val:.10f}")
            gp("Lr = lfun(Ldata, s_r)")
            gp("Lpr = lfun(Ldata, s_r, 1)")
            abs_Lr = float(gp("abs(Lr)"))
            if abs_Lr < 1e-30:
                c0_raw.append(float('nan')); continue
            gp("conn_Lr = Lpr/Lr")
            gp("conn_Gr = -2*log(Pi)+0.5*psi((s_r+5.5)/2)+0.5*psi((s_r+6.5)/2)"
               "+0.5*psi((s_r+16.5)/2)+0.5*psi((s_r+17.5)/2)")
            gp("conn_tr = conn_Lr + conn_Gr")
            conn_re = float(gp("real(conn_tr)"))
            conn_im = float(gp("imag(conn_tr)"))
            c0_val = conn_re - 1.0/eps_f + 1j*conn_im
            c0_raw.append(c0_val)
        except Exception as e:
            c0_raw.append(float('nan'))

    valid = [(i,v) for i,v in enumerate(c0_raw) if np.isfinite(np.real(v))]
    if len(valid) < 3:
        return A_dir, float('nan'), True

    c0_arr = np.array([v for _,v in valid])
    r1 = [(4*c0_arr[i+1]-c0_arr[i])/3 for i in range(len(c0_arr)-1)]
    if len(r1) >= 2:
        r2 = [(4*r1[i+1]-r1[i])/3 for i in range(len(r1)-1)]
        c0_rich = r2[0] if len(r2)==1 else (4*r2[-1]-r2[-2])/3
    else:
        c0_rich = r1[0]
    return A_dir, float(np.real(c0_rich)), True

# ━━━━━━━━━━━ Step 4: 테스트 영점별 Hadamard 검증 ━━━━━━━━━━━
log("=" * 72)
log(f"[Step 4] Hadamard 검증 (N_total={N_ALL} 영점)")
log("=" * 72)

# 테스트 영점: t < 18 영점 전부 + t ∈ [18, 20] 영점 (A_direct 폴백)
# t > 20은 lfun(Ldata) 자체가 부정확하므로 제외
T_RELIABLE = 20.0
reliable_zeros = nz_zeros[nz_zeros < T_RELIABLE]
test_zeros = reliable_zeros  # 모두 테스트

log(f"  테스트: {len(test_zeros)}개 (t < {T_RELIABLE})")
if len(test_zeros) > 0:
    log(f"    t = {test_zeros[0]:.4f} ~ {test_zeros[-1]:.4f}")
log(f"  제외: {len(nz_zeros) - len(test_zeros)}개 (t ≥ {T_RELIABLE}, lfun 정밀도 부족)")
log()
flush_file()

results = []
for idx, t0_val in enumerate(test_zeros):
    log(f"  ── 영점 #{idx+1}: t₀={t0_val:.6f} ──")
    tz = time.time()

    A_dir, Re_c0, success = compute_Re_c0_richardson(t0_val)
    if not np.isfinite(A_dir):
        log(f"    ⚠️ skip (A_dir NaN)")
        results.append({'t': t0_val, 'skip': True})
        log(); flush_file(); continue

    log(f"    A_direct = {A_dir:.6f}")
    log(f"    Re(c₀) = {Re_c0:.4e}  {'✅' if np.isfinite(Re_c0) and abs(Re_c0)<0.01 else '⚠️'}")

    # Richardson 보정: |Re(c₀)| < 0.01이면 적용, 아니면 A_direct 사용
    if np.isfinite(Re_c0) and abs(Re_c0) < 0.01:
        A_true = A_dir - 2.0 * Re_c0 / DELTA_DIRECT
        correction = 2.0 * Re_c0 / DELTA_DIRECT
        log(f"    A_true = {A_true:.6f}  (Richardson 보정: {correction:.4f})")
    else:
        A_true = A_dir; correction = 0.0
        if np.isfinite(Re_c0):
            log(f"    A_true = A_direct = {A_true:.6f}  (Re(c₀)={Re_c0:.2e} → Richardson 미적용)")
        else:
            log(f"    A_true = A_direct = {A_true:.6f}  (Re(c₀) NaN → Richardson 미적용)")

    B_had, H1_had, A_had = hadamard_B_H1(t0_val, all_zeros)
    log(f"    B_had={B_had:.6f}, H₁_had={H1_had:.6f}, A_had(raw)={A_had:.6f}")

    Bt = H1t = 0.0
    if HAS_TAIL and t_last > t0_val + 1.0:
        try:
            Bt = tail_gl4(t0_val, t_last, 'B')
            H1t = tail_gl4(t0_val, t_last, 'H1')
        except Exception as e:
            log(f"    WARNING EM: {e}")

    A_corr = (B_had + Bt)**2 + 2.0*(H1_had + H1t)
    log(f"    EM: B_tail={Bt:.5f}, H₁_tail={H1t:.5f}")
    log(f"    A_corr = {A_corr:.6f}")

    abs_err = abs(A_corr - A_true)
    if abs(A_true) > 0.5:
        err_pct = abs_err / abs(A_true) * 100.0; err_mode = "rel"
    else:
        err_pct = abs_err / (1.0/DELTA_DIRECT**2) * 100.0; err_mode = "abs/κ"

    pass_str = '✓<5%' if err_pct < 5.0 else ('✓<10%' if err_pct < 10.0 else '✗')
    log(f"    err={err_pct:.4f}% ({err_mode}) |ΔA|={abs_err:.6f}  {pass_str}")
    log(f"    소요: {time.time()-tz:.1f}s")

    results.append({
        't': t0_val, 'A_direct': A_dir, 'A_true': A_true,
        'Re_c0': Re_c0, 'correction': correction,
        'B_had': B_had, 'H1_had': H1_had,
        'Bt': Bt, 'H1t': H1t, 'A_corr': A_corr,
        'abs_err': abs_err, 'err_pct': err_pct, 'err_mode': err_mode,
        'skip': False
    })
    log(); flush_file()

# ━━━━━━━━━━━ 최종 판정 ━━━━━━━━━━━
log("=" * 72)
log("최종 판정 — Hadamard GL(4) sym³(Δ) ξ-bundle (#102)")
log("=" * 72)
log()

valid_results = [r for r in results if not r['skip']]
nv = len(valid_results)
n_5pct  = sum(1 for r in valid_results if r['err_pct'] < 5.0)
n_10pct = sum(1 for r in valid_results if r['err_pct'] < 10.0)
n_rec0  = sum(1 for r in valid_results
              if np.isfinite(r.get('Re_c0', float('nan'))) and abs(r['Re_c0']) < 0.01)

log(f"{'#':>3} {'t₀':>10} {'A_true':>10} {'A_corr':>10} {'|ΔA|':>10} {'err%':>8}  pass")
log("-" * 70)
for i, r in enumerate(results):
    if r['skip']:
        log(f"  {i+1:>2} {r['t']:>10.4f}  SKIP"); continue
    ok = '✓<5%' if r['err_pct'] < 5.0 else ('✓<10%' if r['err_pct'] < 10.0 else '✗')
    log(f"  {i+1:>2} {r['t']:>10.4f} {r['A_true']:>10.4f} {r['A_corr']:>10.4f} "
        f"{r['abs_err']:>10.6f} {r['err_pct']:>8.4f}  {ok}")
log()

log(f"판정 (N_zeros={N_ALL}, FE={fe_score:.0f}, ξ-bundle via lfun+ψ [NO lfuninit]):")
log(f"  유효: {nv}개, <5%: {n_5pct}/{nv}, <10%: {n_10pct}/{nv}")
log(f"  Re(c₀)<0.01: {n_rec0}/{nv}")
log()

if nv >= 10 and n_5pct >= 7:
    verdict = "★★ 양성 — d=4 Hadamard 5% 기준 통과"
elif nv >= 10 and n_5pct >= 5:
    verdict = "★ 조건부 양성"
elif nv >= 10 and n_10pct >= 7:
    verdict = "★ 약양성 (10% 기준)"
elif nv >= 3:
    rate = n_5pct/max(nv,1)
    if rate >= 0.7: verdict = f"★★ 양성 ({n_5pct}/{nv} <5%, 영점 {nv}<10)"
    elif rate >= 0.5: verdict = f"★ 조건부 양성 ({n_5pct}/{nv} <5%)"
    else: verdict = f"⚠️ 판정 유보 ({n_5pct}/{nv} <5%)"
else:
    verdict = f"⚠️ 영점 부족 ({nv}개)"

log(f"  {verdict}")
log()

# B-12/B-20
A_vals = [r['A_true'] for r in valid_results]
mean_A = np.mean(A_vals) if A_vals else float('nan')
log("[B-12 ξ-bundle ⟨A⟩]")
log(f"  d=1: 1.27, d=2: 3.93, d=3: 12.79, d=4: {mean_A:.4f} (ξ-bundle)")
log()

if len(valid_results) >= 2:
    log("[B-20 패턴]")
    Apairs = sorted([(r['A_true'], r) for r in valid_results], reverse=True)
    log(f"  최대 A: t₀={Apairs[0][1]['t']:.4f}, A={Apairs[0][0]:.4f}")
    log(f"  최소 A: t₀={Apairs[-1][1]['t']:.4f}, A={Apairs[-1][0]:.4f}")
    log()

log("도표 (Hadamard degree 비교):")
log(f"  d=1: 13/13 <0.003%  ★★★  (#90)")
log(f"  d=2:  7/10 <1.3%    ★★   (#92)")
log(f"  d=3:  6/10 <1%      ★★   (#100)")
s4 = '★★' if (nv>=10 and n_5pct>=7) or (nv>=3 and n_5pct/max(nv,1)>=0.7) else \
     ('★' if n_5pct/max(nv,1)>=0.5 else '⚠️')
log(f"  d=4:  {n_5pct}/{nv} <5%      {s4}    (#102)")
log()

log("=" * 72)
log("[핵심 발견]")
log("=" * 72)
log()
log("1. #80 파라미터 오류 발견:")
log("   gammaV=[-1,0,0,1], N=144 → FE=-1 → 79 '영점'은 가짜")
log(f"   올바른: gammaV={GAMMA_V}, N=1, w=0 → FE={fe_score:.0f}")
log()
log("2. PARI 대형 감마이동 한계:")
log("   lfuninit: gammaV 17.5 이동에서 수치 오버플로우 (|L| ≈ 10^12)")
log("   lfunzeros: 가짜 영점 반환 (모두 |L| >> 1)")
log("   해결: lfun(Ldata, s) 직접 호출 → t≤25에서 정확")
log()
log("3. ξ-bundle 방법론:")
log("   lfunlambda → 불안정 (계수 부족)")
log("   lfun + 수동 ψ 감마 보정 → 안정 (κ≈1/δ² 확인)")
log()
log("4. #80 4성질 재검토 필요:")
log("   기존 #80 결과(A≈0, 79영점)는 잘못된 파라미터 기반 → 신뢰할 수 없음")
log(f"   새 ⟨A⟩={mean_A:.4f} (진짜 영점 {N_ALL}개, t≤{T_SCAN_MAX})")
log()

log(f"완료: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"결과: {OUTFILE}")
log("=" * 72)
flush_file()
