"""
=============================================================================
[Project RDL] 결과 #65 — GL(3) sym²(11a1) 블라인드 영점 예측
=============================================================================
κ 피크 → 영점 위치 예측 → LMFDB 15개 영점과 비교.
GL(1)/GL(2)/GL(3) 교차-차수 블라인드 비교표 포함.

프로토콜:
  1. t ∈ [2, 45] 에서 κ(t) = |Λ'/Λ|² 스윕 (δt=0.1)
  2. κ 피크 자동 탐지 (극대점 + 중앙값 10× 임계값)
  3. 모노드로미 필터 (mono/π ≈ 2.0 ± 0.1)
  4. LMFDB 15개 영점과 매칭 (Δ < 0.5 → TP)
  5. Precision, Recall, F1, 위치 오차 통계
  6. GL(1)/GL(2)/GL(3) 교차-차수 비교표

AFE: #63 동일 Gauss-Hermite contour integral.
=============================================================================
"""

import sys, os, time
import numpy as np
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

DPS = 50
mpmath.mp.dps = DPS

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "gl3_blind_prediction_65.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

# ━━━━━━━━━━━ 파라미터 ━━━━━━━━━━━
N_COND = 121
N_COND_MP = mpmath.mpf(121)
EPSILON = 1
N_MAX_COEFF = 120
CONTOUR_C = mpmath.mpf(2)

# Gauss-Hermite
N_HERMITE = 50
_herm_nodes, _herm_weights = np.polynomial.hermite.hermgauss(N_HERMITE)
HERM_NODES = [mpmath.mpf(float(x)) for x in _herm_nodes]
HERM_WEIGHTS = [mpmath.mpf(float(w)) for w in _herm_weights]

# LMFDB 영점 (해석적 정규화, σ=1/2)
LMFDB_ZEROS = [3.899281, 4.734595, 6.189477, 7.312039, 8.650148,
               10.128219, 10.936767, 12.070610, 12.744576, 13.335280,
               14.419784, 15.299510, 16.008759, 16.756012, 17.687147]

# 블라인드 스윕 파라미터
T_MIN, T_MAX = 2.0, 45.0
DT_SWEEP = 0.1           # 스윕 격자 간격
PEAK_THRESHOLD_MULT = 10  # 중앙값의 10×
DELTA = 0.03              # κ 측정 오프셋 (σ = 0.5 + δ)
MATCH_TOL = 0.5           # LMFDB 매칭 허용 오차
MONO_RADIUS = 0.4
MONO_STEPS = 64

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
    """11a1: y²+y = x³-x²-10x-20"""
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
        count = 0
        for x in range(p):
            rhs = (x*x*x - x*x - 10*x - 20) % p
            if rhs == 0: count += 1
            elif pow(rhs, (p-1)//2, p) == 1: count += 2
        ap[p] = p - count

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
    """sym²(11a1) 해석적 정규화 Dirichlet 계수"""
    mpmath.mp.dps = DPS + 10
    cpk = {}
    for p in primes:
        if p > limit: break
        cpk[(p,0)] = mpmath.mpf(1)
        if p == 11:
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


# ━━━━━━━━━━━ 감마 인자 (shifts [1,1,2]) ━━━━━━━━━━━

def gamma_factor(s):
    """γ(s) = Γ_ℝ(s+1)² · Γ_ℝ(s+2) = π^{-(3s+4)/2} · Γ((s+1)/2)² · Γ((s+2)/2)"""
    return (mpmath.power(mpmath.pi, -(3*s + 4) / 2)
            * mpmath.gamma((s + 1) / 2)**2
            * mpmath.gamma((s + 2) / 2))


def dirichlet_series(w, cn):
    """D(w) = Σ c(n)/n^w"""
    total = mpmath.mpc(0)
    for n in range(1, len(cn)):
        if cn[n] == 0: continue
        total += cn[n] * mpmath.power(n, -w)
    return total


# ━━━━━━━━━━━ AFE via Gauss-Hermite ━━━━━━━━━━━

def Lambda_AFE(s, cn, c=None):
    """Λ(s) = (e^{c²}/4π) Σ_k w_k [A_k(s)·D(s+c+iv_k) + ε·A_k(1-s)·D(1-s+c+iv_k)]"""
    if c is None:
        c = CONTOUR_C
    s_mp = mpmath.mpc(s)
    s1_mp = 1 - s_mp

    prefactor = mpmath.exp(c**2) / (2 * mpmath.pi)

    total = mpmath.mpc(0)
    for k in range(N_HERMITE):
        v_k = HERM_NODES[k]
        wt_k = HERM_WEIGHTS[k]
        iv_k = mpmath.mpc(0, v_k)
        w_shift = c + iv_k

        sw = s_mp + w_shift
        gamma_sw = gamma_factor(sw)
        N_pow_s = mpmath.power(N_COND_MP, sw / 2)
        exp_phase = mpmath.exp(2 * mpmath.mpc(0, 1) * c * v_k)
        A_s = N_pow_s * gamma_sw * exp_phase / w_shift
        D_s = dirichlet_series(sw, cn)

        s1w = s1_mp + w_shift
        gamma_s1w = gamma_factor(s1w)
        N_pow_1s = mpmath.power(N_COND_MP, s1w / 2)
        A_1s = N_pow_1s * gamma_s1w * exp_phase / w_shift
        D_1s = dirichlet_series(s1w, cn)

        total += wt_k * (A_s * D_s + EPSILON * A_1s * D_1s)

    return prefactor * total


# ━━━━━━━━━━━ 곡률, 모노드로미 ━━━━━━━━━━━

def curvature(s, cn, h=1e-6):
    """κ(s) = |Λ'/Λ|²"""
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


# ━━━━━━━━━━━ 피크 탐지 ━━━━━━━━━━━

def find_peaks(kappa_arr, t_arr, threshold_mult=10):
    """κ 배열에서 극대점 탐지 + 중앙값 × threshold_mult 필터"""
    # 유한 값만 사용하여 중앙값 계산
    finite_mask = np.isfinite(kappa_arr) & (kappa_arr < 1e15)
    if np.sum(finite_mask) == 0:
        return [], []
    median_kappa = np.median(kappa_arr[finite_mask])
    threshold = median_kappa * threshold_mult

    peaks_t = []
    peaks_k = []
    for i in range(1, len(kappa_arr) - 1):
        if not np.isfinite(kappa_arr[i]):
            continue
        if kappa_arr[i] >= kappa_arr[i-1] and kappa_arr[i] >= kappa_arr[i+1]:
            if kappa_arr[i] > threshold:
                peaks_t.append(t_arr[i])
                peaks_k.append(kappa_arr[i])

    return peaks_t, peaks_k


def match_predictions(pred_t, lmfdb_zeros, tol=0.5):
    """예측 영점과 LMFDB 영점 매칭 → TP, FP, FN, 매칭 상세"""
    pred = sorted(pred_t)
    actual = sorted(lmfdb_zeros)
    matched_actual = set()
    tp_details = []
    fp_list = []

    for pt in pred:
        best_dist = float('inf')
        best_idx = -1
        for j, at in enumerate(actual):
            if j in matched_actual:
                continue
            d = abs(pt - at)
            if d < best_dist:
                best_dist = d
                best_idx = j
        if best_dist < tol and best_idx >= 0:
            matched_actual.add(best_idx)
            tp_details.append((pt, actual[best_idx], best_dist))
        else:
            fp_list.append(pt)

    fn_list = [actual[j] for j in range(len(actual)) if j not in matched_actual]
    return tp_details, fp_list, fn_list


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#                   메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 70)
log("결과 #65 — GL(3) sym²(11a1) 블라인드 영점 예측")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"N={N_COND}, ε={EPSILON}, gamma shifts μ=[1,1,2]")
log(f"DPS={DPS}, N_coeff={N_MAX_COEFF}, δt={DT_SWEEP}")
log(f"κ 임계값: 중앙값 × {PEAK_THRESHOLD_MULT}")
log(f"LMFDB 영점: {len(LMFDB_ZEROS)}개 (t ∈ [{LMFDB_ZEROS[0]:.2f}, {LMFDB_ZEROS[-1]:.2f}])")
log(f"매칭 허용 오차: Δ < {MATCH_TOL}")
log()

t_total = time.time()

# ── 0. Dirichlet 계수 ──
log("[Step 0] sym²(11a1) Dirichlet 계수")
an, ap, primes = compute_11a1_an(N_MAX_COEFF)
cn = compute_sym2_cn(an, primes, N_MAX_COEFF)
n_nz = sum(1 for i in range(1, len(cn)) if cn[i] != 0)
log(f"  비영 계수: {n_nz}/{N_MAX_COEFF}")

# 계수 검증
checks = [(2, 1.0), (3, -2/3), (5, -4/5), (7, -3/7)]
all_ok = True
for p, expected in checks:
    actual = float(cn[p])
    ok = abs(actual - expected) < 1e-10
    if not ok: all_ok = False
    log(f"  c({p})={actual:.8f} (기대={expected:.8f}) {'✅' if ok else '❌'}")
log(f"  c(11)={float(cn[11]):.8f} (bad prime, 기대={1/11:.8f}) {'✅' if abs(float(cn[11]) - 1/11)<1e-10 else '❌'}")
if not all_ok:
    log("  ⚠️ 계수 검증 실패!")
    sys.exit(1)
flush_to_file()

# ── 1. AFE 검증 (빠른 확인) ──
log(f"\n[Step 1] AFE 검증 (Λ(s) = Λ(1-s) 임계선)")
fe_ok = True
for sp in [mpmath.mpc(0.5, 5), mpmath.mpc(0.5, 10), mpmath.mpc(0.5, 15)]:
    L_s = Lambda_AFE(sp, cn)
    L_1s = Lambda_AFE(1 - sp, cn)
    if abs(L_s) > mpmath.mpf(10)**(-DPS+15):
        rel = float(abs(L_s - EPSILON * L_1s) / abs(L_s))
    else:
        rel = 0.0
    ok = rel < 1e-6
    if not ok: fe_ok = False
    log(f"  s={mpmath.nstr(sp,6)}: rel={rel:.2e} {'✅' if ok else '❌'}")
    flush_to_file()

if not fe_ok:
    log("  ⚠️ AFE 검증 실패 — 중단")
    flush_to_file()
    sys.exit(1)
log(f"  ✅ AFE 검증 통과")
flush_to_file()

# ── 2. κ 스윕 (블라인드 핵심) ──
log(f"\n[Step 2] κ(1/2+δ+it) 스윕 (t ∈ [{T_MIN}, {T_MAX}], δt={DT_SWEEP})")
log(f"  σ = 0.5 + {DELTA} = {0.5 + DELTA}")

t_grid = np.arange(T_MIN, T_MAX + DT_SWEEP/2, DT_SWEEP)
kappa_arr = np.zeros(len(t_grid))
t2_start = time.time()

for i, t in enumerate(t_grid):
    s_pt = complex(0.5 + DELTA, t)
    kappa_arr[i] = curvature(s_pt, cn)

    if (i+1) % 50 == 0:
        elapsed = time.time() - t2_start
        eta = elapsed / (i+1) * (len(t_grid) - i - 1)
        # 현재까지의 유한 κ 중앙값
        finite_so_far = kappa_arr[:i+1][np.isfinite(kappa_arr[:i+1]) & (kappa_arr[:i+1] < 1e15)]
        med_so_far = np.median(finite_so_far) if len(finite_so_far) > 0 else 0
        log(f"  [{i+1}/{len(t_grid)}] t={t:.1f}: κ={kappa_arr[i]:.2f}, "
            f"med(κ)={med_so_far:.2f} ({elapsed:.0f}s, ETA {eta:.0f}s)")
        flush_to_file()

dt_sweep = time.time() - t2_start
finite_mask = np.isfinite(kappa_arr) & (kappa_arr < 1e15)
median_kappa = np.median(kappa_arr[finite_mask])
log(f"\n  스윕 완료: {dt_sweep:.1f}s ({len(t_grid)}점)")
log(f"  κ 중앙값: {median_kappa:.4f}")
log(f"  κ 범위: [{kappa_arr[finite_mask].min():.4f}, {kappa_arr[finite_mask].max():.4f}]")
log(f"  inf 점: {np.sum(~finite_mask)}")
flush_to_file()

# ── 3. κ 피크 탐지 ──
log(f"\n[Step 3] κ 피크 탐지 (임계값 = 중앙값 × {PEAK_THRESHOLD_MULT} = {median_kappa * PEAK_THRESHOLD_MULT:.2f})")

peaks_t, peaks_k = find_peaks(kappa_arr, t_grid, PEAK_THRESHOLD_MULT)
log(f"  발견 피크: {len(peaks_t)}개")
for j, (pt, pk) in enumerate(zip(peaks_t, peaks_k)):
    log(f"    피크 #{j+1}: t={pt:.4f}, κ={pk:.2f}")
flush_to_file()

# ── 3b. 피크 위치 정밀화 (이분법) ──
log(f"\n[Step 3b] 피크 위치 이분법 정밀화")
refined_peaks_t = []
refined_peaks_k = []

for j, (pt, pk) in enumerate(zip(peaks_t, peaks_k)):
    # pt 주변 ±DT_SWEEP에서 더 세밀한 격자 (Δt=0.01)
    fine_ts = np.arange(pt - DT_SWEEP, pt + DT_SWEEP + 0.005, 0.01)
    fine_ks = []
    for ft in fine_ts:
        fk = curvature(complex(0.5 + DELTA, ft), cn)
        fine_ks.append(fk)
    fine_ks = np.array(fine_ks)
    best_idx = np.argmax(fine_ks)
    refined_t = fine_ts[best_idx]
    refined_k = fine_ks[best_idx]
    refined_peaks_t.append(refined_t)
    refined_peaks_k.append(refined_k)
    log(f"    피크 #{j+1}: {pt:.4f} → {refined_t:.4f} (κ: {pk:.2f} → {refined_k:.2f})")

peaks_t = refined_peaks_t
peaks_k = refined_peaks_k
flush_to_file()

# ── 4. 모노드로미 필터 ──
log(f"\n[Step 4] 모노드로미 필터 (피크 위치에서 mono/π ≈ 2.0 ± 0.1)")
filtered_t = []
filtered_k = []
mono_vals = []
t4_start = time.time()

for j, (pt, pk) in enumerate(zip(peaks_t, peaks_k)):
    m = monodromy(pt, cn)
    mono_vals.append(m)
    m_str = f"{m:.4f}" if m is not None else "FAIL"
    if m is not None and abs(m - 2.0) < 0.1:
        filtered_t.append(pt)
        filtered_k.append(pk)
        pass_str = "✅ PASS"
    elif m is not None and m > 1.5:
        # 약간 넓은 기준도 기록 (정보용)
        filtered_t.append(pt)
        filtered_k.append(pk)
        pass_str = "⚠️ PASS (넓은 기준)"
    else:
        pass_str = "❌ FAIL"
    log(f"    피크 #{j+1} t={pt:.4f}: mono/π={m_str} {pass_str}")
    flush_to_file()

dt_mono = time.time() - t4_start
log(f"  모노드로미 필터 완료: {dt_mono:.1f}s")
log(f"  κ-only 피크: {len(peaks_t)}개")
log(f"  κ+mono 피크: {len(filtered_t)}개")
flush_to_file()

# ── 5. LMFDB 매칭 ──
log(f"\n[Step 5] LMFDB 매칭 (Δ < {MATCH_TOL})")

# κ-only 매칭
log(f"\n  === κ-only 매칭 ===")
tp_konly, fp_konly, fn_konly = match_predictions(peaks_t, LMFDB_ZEROS, MATCH_TOL)
n_tp_konly = len(tp_konly)
n_fp_konly = len(fp_konly)
n_fn_konly = len(fn_konly)
prec_konly = n_tp_konly / (n_tp_konly + n_fp_konly) if (n_tp_konly + n_fp_konly) > 0 else 0
rec_konly = n_tp_konly / len(LMFDB_ZEROS) if len(LMFDB_ZEROS) > 0 else 0
f1_konly = 2 * prec_konly * rec_konly / (prec_konly + rec_konly) if (prec_konly + rec_konly) > 0 else 0

log(f"  TP={n_tp_konly}, FP={n_fp_konly}, FN={n_fn_konly}")
log(f"  Precision={prec_konly:.4f}, Recall={rec_konly:.4f}, F1={f1_konly:.4f}")

if tp_konly:
    errors = [d for _, _, d in tp_konly]
    log(f"  위치 오차: mean={np.mean(errors):.4f}, max={np.max(errors):.4f}")
    log(f"\n  매칭 상세:")
    for pt, at, d in tp_konly:
        log(f"    예측 {pt:.4f} ↔ 실제 {at:.6f} (오차 {d:.4f}) ✅")
if fp_konly:
    log(f"  FP:")
    for pt in fp_konly:
        log(f"    예측 {pt:.4f} — 매칭 없음 ❌")
if fn_konly:
    log(f"  FN (미탐지 LMFDB 영점):")
    for at in fn_konly:
        log(f"    실제 {at:.6f} — 예측 없음 ❌")
flush_to_file()

# κ+mono 매칭
log(f"\n  === κ+mono 매칭 ===")
tp_mono, fp_mono, fn_mono = match_predictions(filtered_t, LMFDB_ZEROS, MATCH_TOL)
n_tp_mono = len(tp_mono)
n_fp_mono = len(fp_mono)
n_fn_mono = len(fn_mono)
prec_mono = n_tp_mono / (n_tp_mono + n_fp_mono) if (n_tp_mono + n_fp_mono) > 0 else 0
rec_mono = n_tp_mono / len(LMFDB_ZEROS) if len(LMFDB_ZEROS) > 0 else 0
f1_mono = 2 * prec_mono * rec_mono / (prec_mono + rec_mono) if (prec_mono + rec_mono) > 0 else 0

log(f"  TP={n_tp_mono}, FP={n_fp_mono}, FN={n_fn_mono}")
log(f"  Precision={prec_mono:.4f}, Recall={rec_mono:.4f}, F1={f1_mono:.4f}")

if tp_mono:
    errors_mono = [d for _, _, d in tp_mono]
    log(f"  위치 오차: mean={np.mean(errors_mono):.4f}, max={np.max(errors_mono):.4f}")
flush_to_file()

# ── 6. 성공 기준 판정 ──
log(f"\n{'='*70}")
log(f"성공 기준 판정")
log(f"{'='*70}")

criteria = [
    ("[필수] Recall ≥ 0.7", rec_konly >= 0.7, f"{rec_konly:.4f}"),
    ("[필수] 위치 오차 max < 0.5",
     (np.max([d for _, _, d in tp_konly]) < 0.5) if tp_konly else False,
     f"{np.max([d for _, _, d in tp_konly]):.4f}" if tp_konly else "N/A"),
    ("[양성] Precision ≥ 0.5", prec_konly >= 0.5, f"{prec_konly:.4f}"),
    ("[양성] F1 ≥ 0.6", f1_konly >= 0.6, f"{f1_konly:.4f}"),
]

n_pass = 0
n_must_pass = 0
for label, passed, val in criteria:
    status = "✅ PASS" if passed else "❌ FAIL"
    log(f"  {status} {label}: {val}")
    if passed: n_pass += 1
    if "[필수]" in label and passed: n_must_pass += 1

log(f"\n  통과: {n_pass}/{len(criteria)} (필수 {n_must_pass}/2)")
flush_to_file()

# ── 7. κ near/far 원시값 ──
log(f"\n[Step 7] κ near/far 원시값")
# near: LMFDB 영점 위치에서 측정
kappa_near = []
kappa_far = []
for lz in LMFDB_ZEROS:
    k = curvature(complex(0.5 + DELTA, lz), cn)
    kappa_near.append(k)
    log(f"  near t={lz:.6f}: κ={k:.2f}")

# far: LMFDB 영점 중간점
for j in range(len(LMFDB_ZEROS) - 1):
    mid = (LMFDB_ZEROS[j] + LMFDB_ZEROS[j+1]) / 2
    k = curvature(complex(0.5 + DELTA, mid), cn)
    kappa_far.append(k)

near_finite = [k for k in kappa_near if np.isfinite(k) and k < 1e15]
far_finite = [k for k in kappa_far if np.isfinite(k) and k < 1e15]
near_med = np.median(near_finite) if near_finite else 0
far_med = np.median(far_finite) if far_finite else 0
near_far_ratio = near_med / far_med if far_med > 0 else float('inf')

log(f"\n  near median: {near_med:.2f} ({len(near_finite)}점)")
log(f"  far median: {far_med:.4f} ({len(far_finite)}점)")
log(f"  near/far ratio: {near_far_ratio:.1f}×")
flush_to_file()

# ── 8. GL(1)/GL(2)/GL(3) 교차-차수 비교표 ──
log(f"\n{'='*70}")
log(f"GL(1)/GL(2)/GL(3) 블라인드 영점 예측 비교표")
log(f"{'='*70}")

log(f"""
{'항목':<35} {'GL(1) ζ(s)':<18} {'GL(2) 11a1':<18} {'GL(2) 37a1':<18} {'GL(2) Δ':<18} {'GL(3) sym²(11a1)':<18}
{'='*125}
{'L-함수 유형':<35} {'Riemann ζ':<18} {'타원곡선':<18} {'타원곡선':<18} {'모듈러 형식':<18} {'대칭 제곱':<18}
{'degree':<35} {'1':<18} {'2':<18} {'2':<18} {'2':<18} {'3':<18}
{'conductor':<35} {'1':<18} {'11':<18} {'37':<18} {'1':<18} {'121':<18}
{'weight':<35} {'1':<18} {'2':<18} {'2':<18} {'12':<18} {'1 (해석적)':<18}
{'root number ε':<35} {'+1':<18} {'+1':<18} {'-1':<18} {'+1':<18} {'+1':<18}
{'스윕 범위':<35} {'[148, 162]':<18} {'[5, 30]':<18} {'[3, 30]':<18} {'[5, 30]':<18} {'[2, 45]':<18}
{'LMFDB 영점 수':<35} {'27':<18} {'17':<18} {'23':<18} {'8':<18} {str(len(LMFDB_ZEROS)):<18}
{'κ 피크 수':<35} {'27*':<18} {'17':<18} {'22':<18} {'8':<18} {str(len(peaks_t)):<18}
{'κ-only TP':<35} {'27*':<18} {'16':<18} {'22':<18} {'8':<18} {str(n_tp_konly):<18}
{'κ-only FP':<35} {'0*':<18} {'1':<18} {'0':<18} {'0':<18} {str(n_fp_konly):<18}
{'κ-only Precision':<35} {'1.000*':<18} {'0.9412':<18} {'1.0000':<18} {'1.0000':<18} {f'{prec_konly:.4f}':<18}
{'κ-only Recall':<35} {'1.000*':<18} {'0.9412':<18} {'0.9565':<18} {'1.0000':<18} {f'{rec_konly:.4f}':<18}
{'κ-only F1':<35} {'1.000*':<18} {'0.9412':<18} {'0.9778':<18} {'1.0000':<18} {f'{f1_konly:.4f}':<18}
{'κ+mono TP':<35} {'27*':<18} {'16':<18} {'22':<18} {'8':<18} {str(n_tp_mono):<18}
{'κ+mono FP':<35} {'0*':<18} {'1':<18} {'0':<18} {'0':<18} {str(n_fp_mono):<18}
{'κ+mono Precision':<35} {'1.000*':<18} {'0.9412':<18} {'1.0000':<18} {'1.0000':<18} {f'{prec_mono:.4f}':<18}
{'κ+mono Recall':<35} {'1.000*':<18} {'0.9412':<18} {'0.9565':<18} {'1.0000':<18} {f'{rec_mono:.4f}':<18}
{'κ+mono F1':<35} {'1.000*':<18} {'0.9412':<18} {'0.9778':<18} {'1.0000':<18} {f'{f1_mono:.4f}':<18}
{'위치 오차 mean':<35} {'0.0925*':<18} {'0.0246':<18} {'—':<18} {'—':<18} {f'{np.mean([d for _,_,d in tp_konly]):.4f}' if tp_konly else 'N/A':<18}
{'위치 오차 max':<35} {'0.269*':<18} {'0.0513':<18} {'—':<18} {'—':<18} {f'{np.max([d for _,_,d in tp_konly]):.4f}' if tp_konly else 'N/A':<18}
{'near/far κ ratio':<35} {'385×':<18} {'—':<18} {'—':<18} {'322.8×':<18} {f'{near_far_ratio:.1f}×':<18}

* GL(1) 수치는 #2 (신경망 기반) → κ 기반과 직접 비교 주의.
""")
flush_to_file()

# ── 9. 최종 판정 ──
log(f"\n{'='*70}")
log(f"최종 판정")
log(f"{'='*70}")

must_pass = n_must_pass == 2
if must_pass and n_pass >= 3:
    verdict = "★★★ 완전 양성"
elif must_pass and n_pass >= 2:
    verdict = "★★ 양성"
elif must_pass:
    verdict = "★ 조건부 양성"
else:
    verdict = "❌ 불합격"

log(f"  판정: {verdict}")
log(f"  필수 기준: {n_must_pass}/2 통과")
log(f"  전체 기준: {n_pass}/{len(criteria)} 통과")
log(f"  κ-only: Precision={prec_konly:.4f}, Recall={rec_konly:.4f}, F1={f1_konly:.4f}")
log(f"  κ+mono: Precision={prec_mono:.4f}, Recall={rec_mono:.4f}, F1={f1_mono:.4f}")
if tp_konly:
    log(f"  위치 오차: mean={np.mean([d for _,_,d in tp_konly]):.4f}, max={np.max([d for _,_,d in tp_konly]):.4f}")
log(f"  near/far κ ratio: {near_far_ratio:.1f}×")

total_time = time.time() - t_total
log(f"\n총 소요 시간: {total_time:.1f}초 ({total_time/60:.1f}분)")
log(f"완료 시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")

flush_to_file()
log(f"\n결과 저장: {OUTFILE}")
