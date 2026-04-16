"""
=============================================================================
[Project RDL] 결과 #68 — GL(3) κ-conductor 스케일링 법칙 정량화
=============================================================================
기존 2점(sym²(11a1) N=121 κ=283×, sym²(37a1) N=1369 κ=222×)에
4점 추가: sym²(19a1) N=361, sym²(43a1) N=1849, sym²(53a1) N=2809, sym²(67a1) N=4489

목표: log(κ) vs log(N) 피팅으로 α 추정 (κ ~ N^{-α})

방법론: #63/#64와 동일한 AFE + Gauss-Hermite.
μ=[1,1,2] 감마 인자, dps=50, 각 곡선별 FE 검증 + 영점 탐색 + κ near/far 비율.

주의:
- 23a1 비존재 (dim S_2(Γ_0(23))=2, Hecke 고유값 비유리수) → 53a1로 대체
- sym²(E) conductor = N_E²/gcd(N_E,f)². 소수 conductor는 N=p².
- 각 곡선의 bad prime 처리: 포인트 카운팅으로 자동 결정
=============================================================================
"""

import sys, os, time
import numpy as np
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

DPS = 50
mpmath.mp.dps = DPS

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "gl3_conductor_scaling_68.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

# ━━━━━━━━━━━ Gauss-Hermite 직교 ━━━━━━━━━━━
N_HERMITE = 50
_herm_nodes, _herm_weights = np.polynomial.hermite.hermgauss(N_HERMITE)
HERM_NODES = [mpmath.mpf(float(x)) for x in _herm_nodes]
HERM_WEIGHTS = [mpmath.mpf(float(w)) for w in _herm_weights]

CONTOUR_C = mpmath.mpf(2)
EPSILON = 1  # sym² always ε=+1

# ━━━━━━━━━━━ 곡선 정의 ━━━━━━━━━━━
# [a1,a2,a3,a4,a6], conductor, N_coeff
CURVES = {
    '19a1': {
        'ainvs': [0, 1, 1, -9, -15],
        'conductor_E': 19,
        'N_coeff': 200,
        'lmfdb_ap': {2: 0, 3: -2, 5: 3, 7: -1},  # 검증용
        'bad_type': 'split',  # LMFDB: split multiplicative at 19
    },
    '43a1': {
        'ainvs': [0, 1, 1, 0, 0],
        'conductor_E': 43,
        'N_coeff': 200,
        'lmfdb_ap': {2: -2, 3: -2, 5: -4, 7: None},
        'bad_type': 'nonsplit',  # LMFDB: nonsplit multiplicative at 43
    },
    '53a1': {
        'ainvs': [1, -1, 1, 0, 0],
        'conductor_E': 53,
        'N_coeff': 220,
        'lmfdb_ap': {},  # 검증 데이터 없음 — FE로 검증
        'bad_type': 'auto',  # 포인트 카운팅으로 자동 결정
    },
    '67a1': {
        'ainvs': [0, 1, 1, -12, -21],
        'conductor_E': 67,
        'N_coeff': 250,
        'lmfdb_ap': {2: 2, 3: -2, 5: 2, 7: -2},
        'bad_type': 'split',  # LMFDB: split multiplicative at 67
    },
}

# ━━━━━━━━━━━ 로깅 ━━━━━━━━━━━
lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))
def flush_to_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")


# ━━━━━━━━━━━ 범용 포인트 카운팅 ━━━━━━━━━━━

def count_points_Fp(ainvs, p):
    """F_p 위의 아핀 점 수 계산 (사영점 +1 별도)
    Weierstrass: y² + a1*x*y + a3*y = x³ + a2*x² + a4*x + a6
    """
    a1, a2, a3, a4, a6 = ainvs
    count = 0
    for x in range(p):
        # y² + (a1*x + a3)*y - (x³ + a2*x² + a4*x + a6) = 0 mod p
        # Let A = a1*x + a3, B = -(x³ + a2*x² + a4*x + a6)
        # y² + A*y + B = 0
        A = (a1 * x + a3) % p
        B = (-(x * x * x + a2 * x * x + a4 * x + a6)) % p

        # Discriminant: A² - 4B
        disc = (A * A - 4 * B) % p
        if disc == 0:
            count += 1  # double root
        elif p == 2:
            # For p=2, check y=0 and y=1
            for y in range(2):
                val = (y * y + A * y + B) % 2
                if val == 0:
                    count += 1
        else:
            # Legendre symbol
            leg = pow(disc, (p - 1) // 2, p)
            if leg == 1:
                count += 2
    return count


def compute_ap(ainvs, p):
    """aₚ = p + 1 - #E(F_p) (projective)"""
    affine = count_points_Fp(ainvs, p)
    return p + 1 - (affine + 1)  # +1 for point at infinity = p - affine


def compute_an_general(ainvs, limit, bad_prime):
    """범용 Hecke 고유값 aₙ (n ≤ limit)"""
    # Sieve primes
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False
    primes = [i for i in range(2, limit + 1) if sieve[i]]

    # Compute a_p for all primes
    ap = {}
    for p in primes:
        ap[p] = compute_ap(ainvs, p)

    # a_f(p^k) 재귀
    apk = {}
    for p in primes:
        apk[(p, 0)] = 1
        apk[(p, 1)] = ap[p]
        pk = p
        k = 1
        while pk * p <= limit:
            pk *= p
            k += 1
            if p == bad_prime:
                apk[(p, k)] = ap[p] ** k  # bad prime
            else:
                apk[(p, k)] = ap[p] * apk[(p, k - 1)] - p * apk[(p, k - 2)]

    # 곱셈적 구조
    an = [0] * (limit + 1)
    an[1] = 1
    for n in range(2, limit + 1):
        temp = n
        result = 1
        for p in primes:
            if p * p > temp:
                break
            if temp % p == 0:
                k = 0
                while temp % p == 0:
                    k += 1
                    temp //= p
                result *= apk.get((p, k), 0)
        if temp > 1:
            result *= ap.get(temp, 0)
        an[n] = result

    return an, ap, primes


def compute_sym2_cn(an, ap, primes, limit, bad_prime):
    """sym²(E) 해석적 정규화 Dirichlet 계수
    Good prime p: c(p) = a_f(p)²/p - 1
    Bad prime (split, a_p=1): c(p^k) = p^{-k}
    Bad prime (nonsplit, a_p=-1): c(p^k) = (-1)^k · p^{-k}
    재귀 (good): c(p^k) = c₁·c(p^{k-1}) - c₁·c(p^{k-2}) + c(p^{k-3})
    """
    mpmath.mp.dps = DPS + 10
    a_bad = ap[bad_prime]

    cpk = {}
    for p in primes:
        if p > limit:
            break
        cpk[(p, 0)] = mpmath.mpf(1)
        if p == bad_prime:
            for k in range(1, 20):
                cpk[(p, k)] = mpmath.power(mpmath.mpf(a_bad), k) * mpmath.power(mpmath.mpf(p), -k)
                if p ** (k + 1) > limit:
                    break
        else:
            c1 = mpmath.mpf(an[p]) ** 2 / mpmath.mpf(p) - 1
            cpk[(p, 1)] = c1
            pk = p
            k = 1
            while pk * p <= limit:
                pk *= p
                k += 1
                bkm1 = cpk[(p, k - 1)]
                bkm2 = cpk[(p, k - 2)]
                bkm3 = cpk.get((p, k - 3), mpmath.mpf(0))
                cpk[(p, k)] = c1 * bkm1 - c1 * bkm2 + bkm3

    cn = [mpmath.mpf(0)] * (limit + 1)
    cn[1] = mpmath.mpf(1)
    for n in range(2, limit + 1):
        temp = n
        result = mpmath.mpf(1)
        for p in primes:
            if p * p > temp:
                break
            if temp % p == 0:
                k = 0
                while temp % p == 0:
                    k += 1
                    temp //= p
                if (p, k) in cpk:
                    result *= cpk[(p, k)]
                else:
                    result = mpmath.mpf(0)
                    break
        if temp > 1:
            if (temp, 1) in cpk:
                result *= cpk[(temp, 1)]
            else:
                result = mpmath.mpf(0)
        cn[n] = result
    mpmath.mp.dps = DPS
    return cn


# ━━━━━━━━━━━ 감마 인자 (μ=[1,1,2], #63/#64 동일) ━━━━━━━━━━━

def gamma_factor(s):
    """γ(s) = Γ_ℝ(s+1)² · Γ_ℝ(s+2)
    = π^{-(3s+4)/2} · Γ((s+1)/2)² · Γ((s+2)/2)
    """
    return (mpmath.power(mpmath.pi, -(3 * s + 4) / 2)
            * mpmath.gamma((s + 1) / 2) ** 2
            * mpmath.gamma((s + 2) / 2))


def dirichlet_series(w, cn):
    """D(w) = Σ c(n)/n^w"""
    total = mpmath.mpc(0)
    for n in range(1, len(cn)):
        if cn[n] == 0:
            continue
        total += cn[n] * mpmath.power(n, -w)
    return total


def Lambda_AFE(s, cn, N_cond):
    """AFE via Gauss-Hermite contour integral"""
    N_MP = mpmath.mpf(N_cond)
    c = CONTOUR_C
    s_mp = mpmath.mpc(s)
    s1_mp = 1 - s_mp

    prefactor = mpmath.exp(c ** 2) / (2 * mpmath.pi)

    total = mpmath.mpc(0)
    for k in range(N_HERMITE):
        v_k = HERM_NODES[k]
        wt_k = HERM_WEIGHTS[k]

        iv_k = mpmath.mpc(0, v_k)
        w_shift = c + iv_k

        # A_k(s)
        sw = s_mp + w_shift
        gamma_sw = gamma_factor(sw)
        N_pow_s = mpmath.power(N_MP, sw / 2)
        exp_phase = mpmath.exp(2 * mpmath.mpc(0, 1) * c * v_k)
        A_s = N_pow_s * gamma_sw * exp_phase / w_shift

        D_s = dirichlet_series(sw, cn)

        # A_k(1-s)
        s1w = s1_mp + w_shift
        gamma_s1w = gamma_factor(s1w)
        N_pow_1s = mpmath.power(N_MP, s1w / 2)
        A_1s = N_pow_1s * gamma_s1w * exp_phase / w_shift

        D_1s = dirichlet_series(s1w, cn)

        total += wt_k * (A_s * D_s + EPSILON * A_1s * D_1s)

    return prefactor * total


def curvature(t, cn, N_cond, delta=0.03, h=1e-6):
    """κ(s) = |Λ'/Λ|² at s = 0.5+δ+it"""
    s = mpmath.mpc(0.5 + delta, t)
    h_mp = mpmath.mpf(h)
    try:
        L0 = Lambda_AFE(s, cn, N_cond)
        if abs(L0) < mpmath.mpf(10) ** (-DPS + 15):
            return float('inf')
        Lp = Lambda_AFE(s + h_mp, cn, N_cond)
        Lm = Lambda_AFE(s - h_mp, cn, N_cond)
        conn = (Lp - Lm) / (2 * h_mp * L0)
        k = float(abs(conn) ** 2)
        return k if np.isfinite(k) else float('inf')
    except Exception as e:
        print(f"  WARNING curvature at t={t}: {e}", flush=True)
        return 0.0


def monodromy(t_center, cn, N_cond, radius=0.4, n_steps=32):
    """폐곡선 적분으로 모노드로미 측정"""
    center = mpmath.mpc(0.5, t_center)
    phase_accum = mpmath.mpf(0)
    prev_val = None

    for j in range(n_steps + 1):
        theta = 2 * mpmath.pi * j / n_steps
        pt = center + radius * mpmath.exp(1j * theta)
        try:
            val = Lambda_AFE(pt, cn, N_cond)
            if abs(val) < mpmath.mpf(10) ** (-DPS + 15):
                return None
        except Exception as e:
            return None

        if prev_val is not None:
            ratio = val / prev_val
            if abs(ratio) < mpmath.mpf(10) ** (-DPS + 15):
                return None
            phase_accum += mpmath.im(mpmath.log(ratio))
        prev_val = val

    return float(abs(phase_accum) / mpmath.pi)


# ━━━━━━━━━━━ 단일 곡선 분석 함수 ━━━━━━━━━━━

def analyze_curve(label, info):
    """한 곡선에 대한 전체 분석: FE, 영점, κ, 모노드로미"""
    ainvs = info['ainvs']
    p_bad = info['conductor_E']
    N_cond = p_bad ** 2  # sym² conductor
    N_coeff = info['N_coeff']

    log(f"\n{'='*70}")
    log(f"곡선: {label} — sym²({label}) conductor N={N_cond}")
    log(f"{'='*70}")
    log(f"  Weierstrass: {ainvs}")
    log(f"  Bad prime: {p_bad}, N_coeff={N_coeff}")

    t_curve = time.time()

    # 0. Hecke 고유값 + sym² 계수
    log(f"\n  [Step 0] Hecke 고유값 계산")
    an, ap, primes = compute_an_general(ainvs, N_coeff, p_bad)

    # Bad prime a_p 확인
    a_bad = ap[p_bad]
    bad_type = "split" if a_bad == 1 else "nonsplit" if a_bad == -1 else f"unusual({a_bad})"
    log(f"  a_{p_bad} = {a_bad} ({bad_type})")

    if abs(a_bad) != 1:
        log(f"  ⚠️ a_{p_bad}={a_bad} ≠ ±1 — 곡선 모델 오류 가능")
        return None

    # LMFDB a_p 검증
    lmfdb = info.get('lmfdb_ap', {})
    for p_check, expected in lmfdb.items():
        if expected is not None:
            actual = ap.get(p_check, None)
            ok = actual == expected
            log(f"  a_{p_check}={actual} (LMFDB: {expected}) {'✅' if ok else '❌'}")
            if not ok:
                log(f"  ⚠️ a_p 불일치 — 모델 오류!")

    # sym² 계수
    cn = compute_sym2_cn(an, ap, primes, N_coeff, p_bad)
    n_nz = sum(1 for i in range(1, len(cn)) if cn[i] != 0)
    log(f"  sym² 비영 계수: {n_nz}/{N_coeff}")

    # sym² 계수 검증 (good primes)
    for p_check in [2, 3, 5, 7]:
        if p_check != p_bad:
            expected = ap[p_check] ** 2 / p_check - 1
            actual = float(cn[p_check])
            ok = abs(actual - expected) < 1e-10
            log(f"  c({p_check})={actual:.6f} (a²/p-1={expected:.6f}) {'✅' if ok else '❌'}")

    flush_to_file()

    # 1. 함수방정식 검증 (임계선 3점)
    log(f"\n  [Step 1] 함수방정식 Λ(s) = Λ(1-s)")
    fe_pts = [mpmath.mpc(0.5, 5), mpmath.mpc(0.5, 10), mpmath.mpc(0.5, 15)]
    fe_ok = True
    fe_max_rel = 0.0

    for sp in fe_pts:
        t_s = time.time()
        L_s = Lambda_AFE(sp, cn, N_cond)
        L_1s = Lambda_AFE(1 - sp, cn, N_cond)
        dt_s = time.time() - t_s

        if abs(L_s) > mpmath.mpf(10) ** (-DPS + 15):
            rel = float(abs(L_s - EPSILON * L_1s) / abs(L_s))
        else:
            rel = 0.0

        fe_max_rel = max(fe_max_rel, rel)
        ok = rel < 1e-6
        if not ok:
            fe_ok = False
        log(f"    s={mpmath.nstr(sp, 6)}: rel={rel:.2e} {'✅' if ok else '❌'} ({dt_s:.1f}s)")
        flush_to_file()

    if not fe_ok:
        log(f"  ❌ 함수방정식 실패 (max_rel={fe_max_rel:.2e}) — 이 곡선 건너뜀")
        flush_to_file()
        return None

    log(f"  ✅ 함수방정식 통과 (max_rel={fe_max_rel:.2e})")
    flush_to_file()

    # 2. 영점 스캔 (t ∈ [2, 25] — 축소 범위, 스케일링에 충분)
    T_MIN, T_MAX = 2.0, 25.0
    log(f"\n  [Step 2] Λ(1/2+it) 스캔 (t ∈ [{T_MIN}, {T_MAX}])")

    scan_ts = np.arange(T_MIN, T_MAX + 0.5, 0.5)
    scan_vals = []
    t_scan = time.time()

    for i, t in enumerate(scan_ts):
        val = Lambda_AFE(mpmath.mpc(0.5, t), cn, N_cond)
        re_val = float(mpmath.re(val))
        scan_vals.append(re_val)

        if (i + 1) % 15 == 0:
            elapsed = time.time() - t_scan
            eta = elapsed / (i + 1) * (len(scan_ts) - i - 1)
            log(f"    [{i+1}/{len(scan_ts)}] t={t:.1f}: Re(Λ)={re_val:+.4e} ({elapsed:.0f}s, ETA {eta:.0f}s)")
            flush_to_file()

    scan_arr = np.array(scan_vals)
    dt_scan = time.time() - t_scan
    log(f"  스캔 완료: {dt_scan:.1f}s")

    # 3. 영점 이분법 정밀화
    zero_approx = []
    for i in range(len(scan_arr) - 1):
        if scan_arr[i] * scan_arr[i + 1] < 0:
            t_lo, t_hi = float(scan_ts[i]), float(scan_ts[i + 1])
            val_lo = scan_arr[i]
            try:
                for _ in range(15):  # 15회 이분법 (정밀도 ~1e-5 충분)
                    t_mid = (t_lo + t_hi) / 2
                    val_mid = float(mpmath.re(Lambda_AFE(mpmath.mpc(0.5, t_mid), cn, N_cond)))
                    if not np.isfinite(val_mid):
                        break
                    if val_mid * val_lo < 0:
                        t_hi = t_mid
                    else:
                        t_lo = t_mid
                        val_lo = val_mid
                z_t = (t_lo + t_hi) / 2
                zero_approx.append(z_t)
            except Exception as e:
                zero_approx.append((t_lo + t_hi) / 2)

    n_zeros = len(zero_approx)
    log(f"  발견 영점: {n_zeros}개")
    if n_zeros == 0:
        log(f"  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
        flush_to_file()
        return None

    for j, z in enumerate(zero_approx[:10]):
        log(f"    #{j+1}: γ = {z:.8f}")
    if n_zeros > 10:
        log(f"    ... ({n_zeros - 10}개 추가)")
    flush_to_file()

    # 4. κ near/far 측정
    log(f"\n  [Step 3] κ near/far 비율")
    DELTA = 0.03
    n_tp = min(6, n_zeros)  # 6 TP (충분한 median 추정)
    kappa_near = []
    mono_results = []

    for j, z_t in enumerate(zero_approx[:n_tp]):
        k_val = curvature(z_t, cn, N_cond, delta=DELTA)
        m_val = monodromy(z_t, cn, N_cond)
        kappa_near.append(k_val)
        mono_results.append(m_val)
        m_str = f"{m_val:.4f}" if m_val is not None else "FAIL"
        log(f"    TP #{j+1} γ={z_t:.6f}: κ={k_val:.2f}, mono/π={m_str}")
        flush_to_file()

    # FP (중간점)
    kappa_far = []
    fp_ts = []
    for j in range(len(zero_approx) - 1):
        mid = (zero_approx[j] + zero_approx[j + 1]) / 2
        fp_ts.append(mid)
    if zero_approx[0] > T_MIN + 1:
        fp_ts.insert(0, T_MIN + 0.5)

    n_fp = min(5, len(fp_ts))
    for j, t in enumerate(fp_ts[:n_fp]):
        k_val = curvature(t, cn, N_cond, delta=DELTA)
        kappa_far.append(k_val)
        log(f"    FP #{j+1} t={t:.3f}: κ={k_val:.2f}")
        flush_to_file()

    # κ 비율 계산
    near_fin = [k for k in kappa_near if np.isfinite(k) and k < 1e15]
    far_fin = [k for k in kappa_far if np.isfinite(k) and k < 1e15]

    if near_fin and far_fin:
        near_med = float(np.median(near_fin))
        far_med = float(np.median(far_fin))
        if far_med > 0:
            ratio = near_med / far_med
        else:
            ratio = float('inf')
        near_cv = float(np.std(near_fin) / np.mean(near_fin) * 100) if np.mean(near_fin) > 0 else 999
    else:
        near_med = 0
        far_med = 0
        ratio = 0
        near_cv = 999

    log(f"\n    κ near median = {near_med:.2f} (n={len(near_fin)}, CV={near_cv:.1f}%)")
    log(f"    κ far median  = {far_med:.4f} (n={len(far_fin)})")
    log(f"    κ ratio = {ratio:.1f}×")

    # 모노드로미 요약
    tp_mono = sum(1 for m in mono_results if m is not None and 1.5 < m < 2.5)
    total_mono = len([m for m in mono_results if m is not None])
    log(f"    모노드로미: {tp_mono}/{total_mono} TP (mono/π≈2.0)")

    total_time = time.time() - t_curve
    log(f"\n  소요: {total_time:.1f}s ({total_time/60:.1f}분)")
    flush_to_file()

    return {
        'label': label,
        'conductor_E': p_bad,
        'N_cond': N_cond,
        'n_zeros': n_zeros,
        'kappa_ratio': ratio,
        'kappa_near_med': near_med,
        'kappa_far_med': far_med,
        'near_cv': near_cv,
        'fe_max_rel': fe_max_rel,
        'mono_tp': tp_mono,
        'mono_total': total_mono,
        'a_bad': a_bad,
        'bad_type': bad_type,
        'time_s': total_time,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#                   메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 70)
log("결과 #68 — GL(3) κ-conductor 스케일링 법칙 정량화")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"DPS={DPS}, N_hermite={N_HERMITE}, c={float(CONTOUR_C)}")
log(f"gamma shifts μ=[1,1,2] (해석적 정규화)")
log(f"대상 곡선: {list(CURVES.keys())}")
log(f"기존 데이터: sym²(11a1) N=121 κ=283×, sym²(37a1) N=1369 κ=222×")
log()

t_total = time.time()
results = []

# 기존 데이터 포인트 (하드코딩)
existing_data = [
    {'label': '11a1', 'N_cond': 121, 'kappa_ratio': 283.3, 'kappa_near_med': 1123.53,
     'n_zeros': 46, 'fe_max_rel': 0.0, 'mono_tp': 12, 'mono_total': 12},
    {'label': '37a1', 'N_cond': 1369, 'kappa_ratio': 222.1, 'kappa_near_med': 1127.21,
     'n_zeros': 44, 'fe_max_rel': 0.0, 'mono_tp': 15, 'mono_total': 15},
]

# 4곡선 분석
for label in CURVES:
    result = analyze_curve(label, CURVES[label])
    if result is not None:
        results.append(result)
    else:
        log(f"\n  ⚠️ {label} 실패 — 스케일링 분석에서 제외")
    flush_to_file()

# ━━━━━━━━━━━ 스케일링 분석 ━━━━━━━━━━━
log(f"\n{'='*70}")
log(f"κ-conductor 스케일링 분석")
log(f"{'='*70}")

# 전체 데이터 통합
all_data = existing_data + results
n_total = len(all_data)

log(f"\n전체 데이터: {n_total}점")
log(f"\n{'#':3s} | {'곡선':12s} | {'N':6s} | {'κ ratio':10s} | {'κ near':10s} | {'영점':5s} | {'FE':10s} | {'mono':8s}")
log(f"{'─'*3} | {'─'*12} | {'─'*6} | {'─'*10} | {'─'*10} | {'─'*5} | {'─'*10} | {'─'*8}")

for i, d in enumerate(all_data):
    fe_str = f"{d.get('fe_max_rel', 0):.1e}"
    mono_str = f"{d.get('mono_tp', '?')}/{d.get('mono_total', '?')}"
    log(f"{i+1:3d} | {d['label']:12s} | {d['N_cond']:6d} | {d['kappa_ratio']:10.1f} | {d.get('kappa_near_med', 0):10.2f} | {d.get('n_zeros', 0):5d} | {fe_str:10s} | {mono_str:8s}")

flush_to_file()

# log-log 회귀
if n_total >= 3:
    log_N = np.array([np.log(d['N_cond']) for d in all_data])
    log_kappa = np.array([np.log(d['kappa_ratio']) for d in all_data if d['kappa_ratio'] > 0])

    if len(log_kappa) >= 3:
        # 유효 데이터만
        valid = [(d['N_cond'], d['kappa_ratio']) for d in all_data if d['kappa_ratio'] > 0]
        log_N_v = np.array([np.log(n) for n, k in valid])
        log_k_v = np.array([np.log(k) for n, k in valid])

        # 선형 회귀: log(κ) = β₀ + β₁·log(N) → κ = C·N^β₁ → α = -β₁
        n_v = len(valid)
        mean_x = np.mean(log_N_v)
        mean_y = np.mean(log_k_v)
        ss_xx = np.sum((log_N_v - mean_x) ** 2)
        ss_xy = np.sum((log_N_v - mean_x) * (log_k_v - mean_y))
        ss_yy = np.sum((log_k_v - mean_y) ** 2)

        if ss_xx > 0:
            beta1 = ss_xy / ss_xx
            beta0 = mean_y - beta1 * mean_x
            alpha = -beta1

            # R²
            ss_res = np.sum((log_k_v - (beta0 + beta1 * log_N_v)) ** 2)
            R2 = 1 - ss_res / ss_yy if ss_yy > 0 else 0

            # 표준 오차 (beta1)
            if n_v > 2:
                se_beta1 = np.sqrt(ss_res / (n_v - 2) / ss_xx)
                t_stat = beta1 / se_beta1 if se_beta1 > 0 else float('inf')
            else:
                se_beta1 = float('inf')
                t_stat = 0

            C_fit = np.exp(beta0)

            log(f"\n{'─'*50}")
            log(f"log-log 회귀: log(κ) = {beta0:.4f} + {beta1:.4f} · log(N)")
            log(f"  → κ = {C_fit:.2f} · N^({beta1:.4f})")
            log(f"  → α = -β₁ = {alpha:.4f}")
            log(f"  R² = {R2:.4f}")
            log(f"  SE(β₁) = {se_beta1:.4f}")
            log(f"  α = {alpha:.4f} ± {se_beta1:.4f} (1σ)")
            log(f"  데이터 수: {n_v}점")

            # Spearman 순위 상관
            from scipy.stats import spearmanr
            try:
                r_s, p_val = spearmanr([n for n, k in valid], [k for n, k in valid])
                log(f"  Spearman r_s = {r_s:.4f} (p = {p_val:.4f})")
            except Exception:
                log(f"  (scipy 미설치 — Spearman 계산 생략)")

            # 잔차
            log(f"\n  잔차 분석:")
            for i, (n, k) in enumerate(valid):
                pred = C_fit * n ** beta1
                resid = np.log(k) - (beta0 + beta1 * np.log(n))
                log(f"    {all_data[i]['label']:12s}: N={n:6d}, κ={k:7.1f}, 예측={pred:7.1f}, 잔차={resid:+.4f}")

            # 모델 비교: 멱법칙 vs 상수 vs 로그
            log(f"\n  모델 비교:")

            # 상수 모델: κ = mean(κ)
            mean_k = np.mean([k for n, k in valid])
            ss_const = np.sum((np.array([k for n, k in valid]) - mean_k) ** 2)
            log(f"    [A] 상수: κ = {mean_k:.1f}, SS = {ss_const:.1f}")

            # 멱법칙: κ = C·N^α
            pred_power = np.array([C_fit * n ** beta1 for n, k in valid])
            ss_power = np.sum((np.array([k for n, k in valid]) - pred_power) ** 2)
            log(f"    [B] 멱법칙: κ = {C_fit:.1f}·N^({beta1:.4f}), SS = {ss_power:.1f}, R²={R2:.4f}")

            # 로그 모델: κ = a + b·log(N)
            mean_logN = np.mean([np.log(n) for n, k in valid])
            k_arr = np.array([k for n, k in valid])
            logN_arr = np.array([np.log(n) for n, k in valid])
            ss_logN_xx = np.sum((logN_arr - mean_logN) ** 2)
            ss_logN_xy = np.sum((logN_arr - mean_logN) * (k_arr - mean_k))
            if ss_logN_xx > 0:
                b_log = ss_logN_xy / ss_logN_xx
                a_log = mean_k - b_log * mean_logN
                pred_log = a_log + b_log * logN_arr
                ss_log = np.sum((k_arr - pred_log) ** 2)
                R2_log = 1 - ss_log / ss_const if ss_const > 0 else 0
                log(f"    [C] 로그: κ = {a_log:.1f} + {b_log:.1f}·log(N), SS = {ss_log:.1f}, R²={R2_log:.4f}")
            else:
                log(f"    [C] 로그: 계산 불가")

            log(f"\n  최적 모델: {'멱법칙' if R2 > R2_log else '로그'} (R² 기준)")
        else:
            log(f"\n  ⚠️ 분산 0 — 회귀 불가")
    else:
        log(f"\n  ⚠️ 유효 데이터 부족 (κ>0: {len(log_kappa)})")
else:
    log(f"\n  ⚠️ 데이터 부족 ({n_total}점)")

flush_to_file()

# κ 보편상수 분석
log(f"\n{'─'*50}")
log(f"κ near median 보편상수 분석")
all_near = [d.get('kappa_near_med', 0) for d in all_data if d.get('kappa_near_med', 0) > 100]
if len(all_near) >= 2:
    med_all = float(np.median(all_near))
    cv_all = float(np.std(all_near) / np.mean(all_near) * 100)
    log(f"  전체 κ near median: {med_all:.2f} (n={len(all_near)}, CV={cv_all:.1f}%)")
    for d in all_data:
        if d.get('kappa_near_med', 0) > 100:
            dev = (d['kappa_near_med'] - med_all) / med_all * 100
            log(f"    {d['label']:12s}: κ_near={d['kappa_near_med']:.2f} (편차 {dev:+.2f}%)")

# ━━━━━━━━━━━ 종합 판정 ━━━━━━━━━━━
log(f"\n{'='*70}")
log(f"종합 판정")
log(f"{'='*70}")

n_success = len(results)
n_total_pts = len(all_data)
log(f"  새 곡선 성공: {n_success}/4")
log(f"  총 데이터점: {n_total_pts}")

if n_total_pts >= 4:
    log(f"  ✅ 최소 4점 기준 충족")
else:
    log(f"  ❌ 최소 4점 미달")

if n_total_pts >= 6:
    log(f"  ✅ 6점 이상 기준 충족")
else:
    log(f"  ⚠️ 6점 미달 ({n_total_pts}점)")

# FE 전체 통과 여부
fe_all_pass = all(d.get('fe_max_rel', 1) < 1e-6 for d in all_data)
log(f"  FE 전체 통과: {'✅' if fe_all_pass else '❌'}")

total_time = time.time() - t_total
log(f"\n  총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log(f"  시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")

flush_to_file()
log(f"\n결과 저장: {OUTFILE}")
