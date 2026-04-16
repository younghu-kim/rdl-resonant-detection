"""
=============================================================================
[Project RDL] 결과 #67 — GL(3) sym²(11a1) FP 21개 영점 검증
=============================================================================
#65에서 발견된 21개 "FP" (t > 17.7, LMFDB 범위 밖)가
실제 영점인지 고정밀 계산으로 검증.

방법:
  1. dps=80에서 Λ(1/2+it) 직접 계산 → |Λ| 크기 비교
  2. TP 15개 (LMFDB 확인된 영점)에서 동일 계산 → baseline
  3. t ∈ [18, 36] 에서 Δt=0.01 부호변환 탐색 → 독립 영점 위치
  4. 이분법으로 |Λ| 극소 위치 정밀 탐색
  5. FP 21개 중 실제 영점 수 정량 보고 + 교정된 P/R/F1

수학자 성공 기준:
  [필수] 21개 중 L값 계산 완료 ≥ 18개
  [양성] ≥ 15개가 |Λ| < 0.01 → 실제 영점
  [강양성] ≥ 19개 → 교정 P ≥ 0.95
  [음성] ≤ 5개만 영점
=============================================================================
"""

import sys, os, time
import numpy as np
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

# ━━━━━━━━━━━ 정밀도 설정 ━━━━━━━━━━━
DPS = 80
mpmath.mp.dps = DPS

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "gl3_fp_verification_67.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

# ━━━━━━━━━━━ 파라미터 ━━━━━━━━━━━
N_COND = 121
N_COND_MP = mpmath.mpf(121)
EPSILON = 1
N_MAX_COEFF = 200   # dps 80에서 더 많은 계수 사용
CONTOUR_C = mpmath.mpf(3)  # 고 t를 위해 윤곽 약간 넓힘

# Gauss-Hermite
N_HERMITE = 60  # 고정밀을 위해 점 수 증가
_herm_nodes, _herm_weights = np.polynomial.hermite.hermgauss(N_HERMITE)
HERM_NODES = [mpmath.mpf(float(x)) for x in _herm_nodes]
HERM_WEIGHTS = [mpmath.mpf(float(w)) for w in _herm_weights]

# LMFDB 영점 (TP baseline)
LMFDB_ZEROS = [3.899281, 4.734595, 6.189477, 7.312039, 8.650148,
               10.128219, 10.936767, 12.070610, 12.744576, 13.335280,
               14.419784, 15.299510, 16.008759, 16.756012, 17.687147]

# #65에서 발견된 FP 21개 (정밀화 후 위치)
FP_LOCATIONS = [18.7500, 19.2100, 20.1300, 21.1000, 23.0800,
                23.5700, 24.1600, 24.7800, 26.0400, 26.5600,
                27.2400, 29.6500, 30.0100, 30.5600, 31.4600,
                31.8900, 32.3900, 33.1200, 33.6800, 34.4900,
                35.0600]

# ━━━━━━━━━━━ 로깅 ━━━━━━━━━━━
lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))
def flush_to_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")


# ━━━━━━━━━━━ 11a1 Hecke 고유값 (65와 동일) ━━━━━━━━━━━

def compute_11a1_an(limit):
    """11a1: y²+y = x³-x²-10x-20"""
    known_ap = {2:-2, 3:-1, 5:1, 7:-2, 11:1, 13:4, 17:-2, 19:0, 23:-1,
                29:0, 31:7, 37:3, 41:-8, 43:-6, 47:8, 53:-6, 59:5, 61:12,
                67:-7, 71:-4, 73:2, 79:10, 83:-7, 89:6, 97:-4, 101:-2,
                103:12, 107:10, 109:-4, 113:2, 127:-10, 131:-4, 137:8,
                139:-8, 149:14, 151:-6, 157:4, 163:-4, 167:-2, 173:-10,
                179:-6, 181:16, 191:4, 193:-8, 197:-12}

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
    mpmath.mp.dps = DPS + 20
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
    """Λ(s) via Gauss-Hermite contour integral"""
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


# ━━━━━━━━━━━ 부호변환 탐색 ━━━━━━━━━━━

def sign_change_search(t_start, t_end, dt, cn):
    """Re(Λ(1/2+it)) 부호변환 탐색 → 영점 후보 위치 반환"""
    sign_changes = []
    t_vals = np.arange(t_start, t_end + dt/2, dt)
    prev_val = None
    prev_t = None

    n_total = len(t_vals)
    for idx, t in enumerate(t_vals):
        if idx % 50 == 0 and idx > 0:
            elapsed = time.time() - _sc_start
            eta = elapsed / idx * (n_total - idx)
            print(f"  부호변환 탐색: [{idx}/{n_total}] t={t:.2f} ({elapsed:.0f}s, ETA {eta:.0f}s)", flush=True)

        try:
            val = Lambda_AFE(mpmath.mpc(0.5, t), cn)
            re_val = float(mpmath.re(val))
        except Exception as e:
            print(f"  WARNING sign_change at t={t:.4f}: {e}", flush=True)
            prev_val = None
            prev_t = None
            continue

        if prev_val is not None and prev_t is not None:
            if (prev_val > 0 and re_val < 0) or (prev_val < 0 and re_val > 0):
                sign_changes.append((prev_t, t))

        prev_val = re_val
        prev_t = t

    return sign_changes


def bisect_zero(t_lo, t_hi, cn, max_iter=30, tol=1e-10):
    """이분법으로 Re(Λ(1/2+it))=0 위치 정밀 탐색"""
    val_lo = float(mpmath.re(Lambda_AFE(mpmath.mpc(0.5, t_lo), cn)))
    val_hi = float(mpmath.re(Lambda_AFE(mpmath.mpc(0.5, t_hi), cn)))

    if val_lo * val_hi > 0:
        return None, None

    for _ in range(max_iter):
        t_mid = (t_lo + t_hi) / 2
        if t_hi - t_lo < tol:
            break
        val_mid = float(mpmath.re(Lambda_AFE(mpmath.mpc(0.5, t_mid), cn)))
        if val_lo * val_mid <= 0:
            t_hi = t_mid
            val_hi = val_mid
        else:
            t_lo = t_mid
            val_lo = val_mid

    t_zero = (t_lo + t_hi) / 2
    abs_Lambda = float(abs(Lambda_AFE(mpmath.mpc(0.5, t_zero), cn)))
    return t_zero, abs_Lambda


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#                   메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 70)
log("결과 #67 — GL(3) sym²(11a1) FP 21개 영점 검증")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"N={N_COND}, ε={EPSILON}, dps={DPS}, N_coeff={N_MAX_COEFF}")
log(f"Gauss-Hermite 점수: {N_HERMITE}, 윤곽 c={float(CONTOUR_C)}")
log(f"검증 대상: TP {len(LMFDB_ZEROS)}개 + FP {len(FP_LOCATIONS)}개 = {len(LMFDB_ZEROS)+len(FP_LOCATIONS)}개")
log()

t_total_start = time.time()

# ── 0. Dirichlet 계수 (N_coeff=200) ──
log("[Step 0] sym²(11a1) Dirichlet 계수 (N_coeff=200, dps=80)")
an, ap, primes = compute_11a1_an(N_MAX_COEFF)
cn = compute_sym2_cn(an, primes, N_MAX_COEFF)
n_nz = sum(1 for i in range(1, len(cn)) if cn[i] != 0)
log(f"  비영 계수: {n_nz}/{N_MAX_COEFF}")

# 계수 검증
checks = [(2, 1.0), (3, -2/3), (5, -4/5), (7, -3/7), (11, 1/11)]
for n, expected in checks:
    actual = float(cn[n])
    ok = "✅" if abs(actual - expected) < 1e-6 else "❌"
    log(f"  c({n})={actual:.8f} (기대={expected:.8f}) {ok}")
log()

# ── 1. AFE 검증 (함수방정식) ──
log("[Step 1] AFE 검증 (Λ(s) = Λ(1-s), dps=80)")
for t_test in [5.0, 15.0, 25.0, 35.0]:
    s_test = mpmath.mpc(0.5, t_test)
    L1 = Lambda_AFE(s_test, cn)
    L2 = Lambda_AFE(1 - s_test, cn)
    if abs(L1) > 0:
        rel = float(abs(L1 - L2) / abs(L1))
    else:
        rel = float(abs(L1 - L2))
    log(f"  s=(0.5 + {t_test}j): |Λ(s)-Λ(1-s)|/|Λ(s)| = {rel:.2e} {'✅' if rel < 1e-10 else '❌'}")
log()

# ── 2. TP baseline: |Λ(1/2+it)| at LMFDB 영점 ──
log("[Step 2] TP baseline: LMFDB 15개 영점에서 |Λ(1/2+it)| 계산")
tp_abs_values = []
for i, t_zero in enumerate(LMFDB_ZEROS):
    try:
        val = Lambda_AFE(mpmath.mpc(0.5, t_zero), cn)
        abs_val = float(abs(val))
        tp_abs_values.append(abs_val)
        log(f"  TP #{i+1:2d} t={t_zero:.6f}: |Λ| = {abs_val:.6e}")
    except Exception as e:
        log(f"  TP #{i+1:2d} t={t_zero:.6f}: FAIL ({e})")
        tp_abs_values.append(None)

tp_valid = [v for v in tp_abs_values if v is not None]
if tp_valid:
    tp_median = np.median(tp_valid)
    tp_max = max(tp_valid)
    tp_min = min(tp_valid)
    log(f"\n  TP |Λ| 통계: median={tp_median:.6e}, min={tp_min:.6e}, max={tp_max:.6e}")
    log(f"  TP 유효: {len(tp_valid)}/{len(LMFDB_ZEROS)}")
log()
flush_to_file()

# ── 3. FP 검증: |Λ(1/2+it)| at 21개 "FP" 위치 ──
log("[Step 3] FP 검증: 21개 'FP' 위치에서 |Λ(1/2+it)| 계산")
fp_abs_values = []
fp_results = []  # (t, abs_val, is_zero)
for i, t_fp in enumerate(FP_LOCATIONS):
    t0 = time.time()
    try:
        val = Lambda_AFE(mpmath.mpc(0.5, t_fp), cn)
        abs_val = float(abs(val))
        fp_abs_values.append(abs_val)

        # TP baseline과 비교하여 영점 판정
        # TP max 의 10배 이내이면 영점으로 판정 (넉넉한 기준)
        # 더 엄격한 기준: abs_val < 0.01
        is_zero_strict = abs_val < 0.01
        is_zero_loose = abs_val < 0.1

        elapsed = time.time() - t0
        status = "✅ 영점" if is_zero_strict else ("⚠️ 근접" if is_zero_loose else "❌ 비영점")
        fp_results.append((t_fp, abs_val, is_zero_strict, is_zero_loose))
        log(f"  FP #{i+1:2d} t={t_fp:.4f}: |Λ| = {abs_val:.6e} {status} ({elapsed:.1f}s)")
    except Exception as e:
        fp_abs_values.append(None)
        fp_results.append((t_fp, None, False, False))
        log(f"  FP #{i+1:2d} t={t_fp:.4f}: FAIL ({e})")

fp_valid = [v for v in fp_abs_values if v is not None]
if fp_valid:
    fp_median = np.median(fp_valid)
    fp_max = max(fp_valid)
    fp_min = min(fp_valid)
    log(f"\n  FP |Λ| 통계: median={fp_median:.6e}, min={fp_min:.6e}, max={fp_max:.6e}")
    log(f"  FP 유효: {len(fp_valid)}/{len(FP_LOCATIONS)}")
log()

# FP 영점 수 카운트
n_fp_computed = sum(1 for r in fp_results if r[1] is not None)
n_fp_zero_strict = sum(1 for r in fp_results if r[2])  # |Λ| < 0.01
n_fp_zero_loose = sum(1 for r in fp_results if r[3])   # |Λ| < 0.1
log(f"  FP 영점 판정 요약:")
log(f"    계산 완료: {n_fp_computed}/21")
log(f"    엄격 기준 (|Λ|<0.01): {n_fp_zero_strict}/21")
log(f"    완화 기준 (|Λ|<0.1):  {n_fp_zero_loose}/21")
log()
flush_to_file()

# ── 4. 독립 부호변환 탐색 (t ∈ [18, 36], Δt=0.01) ──
log("[Step 4] 독립 부호변환 탐색 (t ∈ [18, 36], Δt=0.05)")
log(f"  Re(Λ(1/2+it)) 부호변환으로 영점 위치 독립 확인")
_sc_start = time.time()
sc_results = sign_change_search(18.0, 36.0, 0.05, cn)
sc_elapsed = time.time() - _sc_start
log(f"  부호변환 탐색 완료: {sc_elapsed:.1f}s")
log(f"  발견된 부호변환: {len(sc_results)}개")
log()
flush_to_file()

# ── 5. 이분법 정밀화 ──
if len(sc_results) > 0:
    log("[Step 5] 이분법 정밀화 (부호변환 → 정확한 영점 위치)")
    independent_zeros = []
    for i, (t_lo, t_hi) in enumerate(sc_results):
        t_zero, abs_at_zero = bisect_zero(t_lo, t_hi, cn, max_iter=40, tol=1e-12)
        if t_zero is not None:
            independent_zeros.append((t_zero, abs_at_zero))
            log(f"  영점 #{i+1:2d}: t={t_zero:.8f}, |Λ|={abs_at_zero:.6e}")
        else:
            log(f"  영점 #{i+1:2d}: 이분법 실패 (t_lo={t_lo:.4f}, t_hi={t_hi:.4f})")
    log(f"\n  독립 확인 영점: {len(independent_zeros)}개 (t ∈ [18, 36])")
    log()
    flush_to_file()

    # ── 6. FP vs 독립 영점 매칭 ──
    log("[Step 6] FP 21개 ↔ 독립 영점 매칭 (Δ < 0.1)")
    matched_fp = 0
    unmatched_fp = []
    for fp_idx, t_fp in enumerate(FP_LOCATIONS):
        if t_fp < 18.0 or t_fp > 36.0:
            continue
        best_dist = float('inf')
        best_iz = None
        for iz_t, iz_abs in independent_zeros:
            d = abs(t_fp - iz_t)
            if d < best_dist:
                best_dist = d
                best_iz = iz_t

        if best_dist < 0.1:
            matched_fp += 1
            log(f"  FP #{fp_idx+1:2d} t={t_fp:.4f} ↔ 독립 영점 t={best_iz:.8f} (Δ={best_dist:.6f}) ✅")
        else:
            unmatched_fp.append(t_fp)
            log(f"  FP #{fp_idx+1:2d} t={t_fp:.4f} — 매칭 없음 (최소 거리={best_dist:.4f}) ❌")

    log(f"\n  FP→독립 영점 매칭: {matched_fp}/{len([t for t in FP_LOCATIONS if 18.0 <= t <= 36.0])}")
    log()
else:
    log("  ⚠️ 부호변환 없음 — 독립 영점 확인 불가")
    independent_zeros = []
    matched_fp = 0
    log()

flush_to_file()

# ── 7. TP 비교 (정규화): |Λ| at TP vs FP ──
log("[Step 7] TP vs FP 비교 — |Λ| 분포 비교")
log()
log(f"  {'위치':>12s} {'t':>10s} {'|Λ|':>14s} {'판정':>8s}")
log(f"  {'─'*12} {'─'*10} {'─'*14} {'─'*8}")
for i, t_zero in enumerate(LMFDB_ZEROS):
    v = tp_abs_values[i]
    if v is not None:
        log(f"  {'TP #'+str(i+1):>12s} {t_zero:10.6f} {v:14.6e} {'✅':>8s}")
log()
for i, (t_fp, abs_val, strict, loose) in enumerate(fp_results):
    if abs_val is not None:
        status = "✅" if strict else ("⚠️" if loose else "❌")
        log(f"  {'FP #'+str(i+1):>12s} {t_fp:10.4f} {abs_val:14.6e} {status:>8s}")
    else:
        log(f"  {'FP #'+str(i+1):>12s} {t_fp:10.4f} {'FAIL':>14s} {'❌':>8s}")
log()

# ── 8. 교정된 메트릭 ──
log("=" * 70)
log("교정된 성능 메트릭")
log("=" * 70)
log()

# 원래 메트릭
log(f"[원래 (#65)] TP=15, FP=21, FN=0")
log(f"  P=0.4167, R=1.0000, F1=0.5882")
log()

# 교정 방법 1: |Λ| 직접 계산 기반
log(f"[교정 방법 1: |Λ| 직접 계산]")
log(f"  엄격 기준 (|Λ|<0.01): {n_fp_zero_strict}/21 FP가 실제 영점")
log(f"  완화 기준 (|Λ|<0.1):  {n_fp_zero_loose}/21 FP가 실제 영점")

new_tp_strict = 15 + n_fp_zero_strict
new_fp_strict = 21 - n_fp_zero_strict
p_strict = new_tp_strict / (new_tp_strict + new_fp_strict) if (new_tp_strict + new_fp_strict) > 0 else 0
r_strict = new_tp_strict / (new_tp_strict + 0)  # FN=0 유지
f1_strict = 2 * p_strict * r_strict / (p_strict + r_strict) if (p_strict + r_strict) > 0 else 0
log(f"  교정 P={p_strict:.4f}, R={r_strict:.4f}, F1={f1_strict:.4f}")
log()

new_tp_loose = 15 + n_fp_zero_loose
new_fp_loose = 21 - n_fp_zero_loose
p_loose = new_tp_loose / (new_tp_loose + new_fp_loose) if (new_tp_loose + new_fp_loose) > 0 else 0
r_loose = new_tp_loose / (new_tp_loose + 0)
f1_loose = 2 * p_loose * r_loose / (p_loose + r_loose) if (p_loose + r_loose) > 0 else 0
log(f"  (완화) 교정 P={p_loose:.4f}, R={r_loose:.4f}, F1={f1_loose:.4f}")
log()

# 교정 방법 2: 독립 부호변환 기반
if len(independent_zeros) > 0:
    log(f"[교정 방법 2: 독립 부호변환 탐색]")
    log(f"  t ∈ [18, 36]에서 독립 영점 {len(independent_zeros)}개 발견")
    log(f"  FP 21개 중 {matched_fp}개가 독립 영점과 매칭 (Δ<0.1)")

    new_tp_ind = 15 + matched_fp
    new_fp_ind = 21 - matched_fp
    p_ind = new_tp_ind / (new_tp_ind + new_fp_ind) if (new_tp_ind + new_fp_ind) > 0 else 0
    f1_ind = 2 * p_ind * 1.0 / (p_ind + 1.0) if (p_ind + 1.0) > 0 else 0
    log(f"  교정 P={p_ind:.4f}, R=1.0000, F1={f1_ind:.4f}")
    log()

# ── 9. Weyl 법칙 비교 ──
log("[참고] Weyl 밀도 비교")
# GL(3) N(T) ~ (T/(2π))^3 * ... 대략적 밀도 추정
# LMFDB 범위: t ∈ [3.9, 17.7] → 15개 → 밀도 = 15/13.8 = 1.087/unit
# FP 범위: t ∈ [18.75, 35.06] → 21개 → 밀도 = 21/16.3 = 1.288/unit
lmfdb_density = 15 / (17.687147 - 3.899281)
fp_density = 21 / (35.06 - 18.75)
log(f"  LMFDB 범위 밀도: {lmfdb_density:.3f}/unit (t ∈ [3.9, 17.7])")
log(f"  FP 범위 밀도: {fp_density:.3f}/unit (t ∈ [18.75, 35.06])")
log(f"  밀도 비: {fp_density/lmfdb_density:.3f} (Weyl 법칙 — 고 t 밀도 증가 예상)")
log()

# ── 최종 판정 ──
log("=" * 70)
log("최종 판정")
log("=" * 70)

# 성공 기준
log(f"  [필수] 계산 완료 ≥ 18/21: {n_fp_computed}/21 {'✅ PASS' if n_fp_computed >= 18 else '❌ FAIL'}")
log(f"  [양성] ≥ 15/21 실제 영점 (|Λ|<0.01): {n_fp_zero_strict}/21 {'✅ PASS' if n_fp_zero_strict >= 15 else '❌ FAIL'}")
log(f"  [강양성] ≥ 19/21 실제 영점: {n_fp_zero_strict}/21 {'✅ PASS' if n_fp_zero_strict >= 19 else '❌ FAIL'}")
log(f"  [음성] ≤ 5/21 영점: {n_fp_zero_strict}/21 {'⚠️ 음성' if n_fp_zero_strict <= 5 else 'N/A'}")

if n_fp_zero_strict >= 19:
    verdict = "★★★ 강양성 — FP 대부분 실제 영점, 교정 P ≥ 0.95"
elif n_fp_zero_strict >= 15:
    verdict = "★★ 양성 — FP 대부분 실제 영점"
elif n_fp_zero_strict >= 6:
    verdict = "★ 약양성 — FP 일부만 영점"
else:
    verdict = "❌ 음성 — FP는 실제 FP, κ 임계값 재조정 필요"

log(f"\n  판정: {verdict}")
log(f"\n  교정 메트릭 (엄격): P={p_strict:.4f}, R=1.0000, F1={f1_strict:.4f}")
if len(independent_zeros) > 0:
    log(f"  교정 메트릭 (독립): P={p_ind:.4f}, R=1.0000, F1={f1_ind:.4f}")

total_time = time.time() - t_total_start
log(f"\n총 소요 시간: {total_time:.1f}초 ({total_time/60:.1f}분)")
log(f"완료 시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")

flush_to_file()
print(f"\n결과 저장: {OUTFILE}")
