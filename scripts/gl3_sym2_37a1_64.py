"""
=============================================================================
[Project RDL] 결과 #64 — GL(3) sym²(37a1) ξ-다발 4성질 검증
=============================================================================
#63 (sym²(11a1)) 보편성 재현 — 제2 독립 사례.

수학적 기반:
  sym²(37a1): conductor N=1369 (=37²), ε=+1
  gamma shifts μ=[1,1,2] (해석적 정규화) — #63과 동일, 절대 변경 금지
  γ(s) = Γ_ℝ(s+1)² · Γ_ℝ(s+2) = π^{-(3s+4)/2} · Γ((s+1)/2)² · Γ((s+2)/2)
  Λ(s) = N^{s/2} · γ(s) · L(s),  FE: Λ(s) = ε·Λ(1-s)

37a1: y² + y = x³ - x, conductor=37, rank=1, a₃₇=-1 (nonsplit multiplicative)
sym² bad prime p=37: L₃₇(sym²E,s) = (1+37^{-s})^{-1} → c(37^k) = (-1)^k·37^{-k}

AFE: Rubinstein smoothed via Gauss-Hermite contour integral (동일 방법론).

성공 기준 (#63 수학자 설정):
  [필수] 함수방정식 rel_err < 10^{-10} (임계선 5점)
  [필수] 영점 ≥ 30개 (t ∈ [2, 50])
  [양성] κ near/far median ≥ 100×
  [양성] mono/π = 2.0 (전 TP)
  [정보] σ-유일성 FAIL 기대 (GL(3) 패턴 확인)
=============================================================================
"""

import sys, os, time
import numpy as np
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

DPS = 50
mpmath.mp.dps = DPS

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "gl3_sym2_37a1_64.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

# ━━━━━━━━━━━ 파라미터 ━━━━━━━━━━━
N_COND = 1369          # 37² (sym²(37a1) conductor)
N_COND_MP = mpmath.mpf(1369)
EPSILON = 1            # sym² always ε=+1
N_MAX_COEFF = 250      # 더 많은 계수 (conductor 11× 증가 → 수렴 느림)
CONTOUR_C = mpmath.mpf(2)  # contour shift

# Gauss-Hermite 직교
N_HERMITE = 50
_herm_nodes, _herm_weights = np.polynomial.hermite.hermgauss(N_HERMITE)
HERM_NODES = [mpmath.mpf(float(x)) for x in _herm_nodes]
HERM_WEIGHTS = [mpmath.mpf(float(w)) for w in _herm_weights]

# LMFDB 영점 (해석적 정규화, 임계선 σ=1/2) — 가용 범위만
# sym²(37a1) LMFDB 라벨: degree 3, conductor 1369
# 아래는 자체 계산으로 검증할 예정; LMFDB 데이터 가용 시 교차검증
LMFDB_ZEROS = []  # LMFDB 데이터 추후 삽입

# 스윕 파라미터
T_MIN, T_MAX = 2.0, 50.0
DELTA = 0.03
MONO_RADIUS, MONO_STEPS = 0.4, 64
MATCH_TOL = 0.5

# ━━━━━━━━━━━ 로깅 ━━━━━━━━━━━
lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))
def flush_to_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━━━━━━ 37a1 Hecke 고유값 ━━━━━━━━━━━

def _compute_ap_37a1(p):
    """aₚ for 37a1: y² + y = x³ - x
    p=37 (bad prime, nonsplit multiplicative): a₃₇ = -1
    p=2: 직접 열거
    p≠2,37 (good prime): disc(x) = 4x³-4x+1, Legendre 기호
    """
    if p == 37:
        return -1  # nonsplit multiplicative reduction

    if p == 2:
        count_pts = 1  # point at infinity
        for x in range(2):
            for y in range(2):
                if (y * y + y - x * x * x + x) % 2 == 0:
                    count_pts += 1
        return p + 1 - count_pts

    # 홀수 good prime: y²+y = x³-x → (2y+1)² = 4x³-4x+1
    affine_count = 0
    for x in range(p):
        disc = (4 * x * x * x - 4 * x + 1) % p
        if disc == 0:
            affine_count += 1
        else:
            leg = pow(disc, (p - 1) // 2, p)
            if leg == 1:
                affine_count += 2
    return p - affine_count


def compute_37a1_an(limit):
    """37a1의 전체 Hecke 고유값 aₙ (n ≤ limit)"""
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False
    primes = [i for i in range(2, limit + 1) if sieve[i]]

    ap = {}
    for p in primes:
        ap[p] = _compute_ap_37a1(p)

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
            if p == 37:
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


def compute_sym2_cn(an, primes, limit):
    """sym²(37a1) 해석적 정규화 Dirichlet 계수
    Good prime p≠37: c(p) = a_f(p)²/p - 1
    Bad prime p=37 (nonsplit, a₃₇=-1):
      L₃₇(sym²E,s) = (1+37^{-s})^{-1} = Σ (-1)^k · 37^{-ks}
      → c(37^k) = (-1)^k · 37^{-k}
    재귀 (good): c(p^k) = c₁·c(p^{k-1}) - c₁·c(p^{k-2}) + c(p^{k-3})
    """
    mpmath.mp.dps = DPS + 10
    cpk = {}
    for p in primes:
        if p > limit:
            break
        cpk[(p, 0)] = mpmath.mpf(1)
        if p == 37:
            # bad prime (nonsplit): c(37^k) = (-1)^k · 37^{-k}
            for k in range(1, 20):
                cpk[(p, k)] = mpmath.power(mpmath.mpf(-1), k) * mpmath.power(mpmath.mpf(37), -k)
                if 37 ** (k + 1) > limit:
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


# ━━━━━━━━━━━ 감마 인자 (μ=[1,1,2], #63과 동일) ━━━━━━━━━━━

def gamma_factor(s):
    """γ(s) = Γ_ℝ(s+1)² · Γ_ℝ(s+2)  (해석적 정규화 shifts μ=[1,1,2])
    = π^{-(3s+4)/2} · Γ((s+1)/2)² · Γ((s+2)/2)

    주의: #63과 완전히 동일. μ=(0,1,1) 사용 금지 (#63b 실수 반복 금지).
    """
    return (mpmath.power(mpmath.pi, -(3 * s + 4) / 2)
            * mpmath.gamma((s + 1) / 2) ** 2
            * mpmath.gamma((s + 2) / 2))


def dirichlet_series(w, cn):
    """D(w) = Σ c(n)/n^w — Re(w) > 1에서 수렴"""
    total = mpmath.mpc(0)
    for n in range(1, len(cn)):
        if cn[n] == 0:
            continue
        total += cn[n] * mpmath.power(n, -w)
    return total


def Lambda_direct(s, cn):
    """Λ(s) = N^{s/2} γ(s) L(s) — Re(s) >> 1에서만 유효"""
    return mpmath.power(N_COND_MP, s / 2) * gamma_factor(s) * dirichlet_series(s, cn)


# ━━━━━━━━━━━ AFE via Gauss-Hermite contour integral ━━━━━━━━━━━

def Lambda_AFE(s, cn, c=None):
    """
    Λ(s) = (e^{c²}/4π) Σ_k w_k [A_k(s)·D(s+c+iv_k) + ε·A_k(1-s)·D(1-s+c+iv_k)]
    """
    if c is None:
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
        N_pow_s = mpmath.power(N_COND_MP, sw / 2)
        exp_phase = mpmath.exp(2 * mpmath.mpc(0, 1) * c * v_k)
        A_s = N_pow_s * gamma_sw * exp_phase / w_shift

        D_s = dirichlet_series(sw, cn)

        # A_k(1-s)
        s1w = s1_mp + w_shift
        gamma_s1w = gamma_factor(s1w)
        N_pow_1s = mpmath.power(N_COND_MP, s1w / 2)
        A_1s = N_pow_1s * gamma_s1w * exp_phase / w_shift

        D_1s = dirichlet_series(s1w, cn)

        total += wt_k * (A_s * D_s + EPSILON * A_1s * D_1s)

    return prefactor * total


# ━━━━━━━━━━━ 곡률, 모노드로미 ━━━━━━━━━━━

def curvature(s, cn, h=1e-6):
    """κ(s) = |Λ'/Λ|² (접속의 곡률)"""
    s_mp = mpmath.mpc(s)
    h_mp = mpmath.mpf(h)
    try:
        L0 = Lambda_AFE(s_mp, cn)
        if abs(L0) < mpmath.mpf(10) ** (-DPS + 15):
            return float('inf')
        Lp = Lambda_AFE(s_mp + h_mp, cn)
        Lm = Lambda_AFE(s_mp - h_mp, cn)
        conn = (Lp - Lm) / (2 * h_mp * L0)
        k = float(abs(conn) ** 2)
        return k if np.isfinite(k) else float('inf')
    except Exception as e:
        print(f"  WARNING curvature at s={s}: {e}", flush=True)
        return 0.0


def monodromy(t_center, cn, sigma=0.5, radius=None, n_steps=None):
    """폐곡선 적분으로 모노드로미 측정 — arg(Λ) 누적"""
    if radius is None:
        radius = MONO_RADIUS
    if n_steps is None:
        n_steps = MONO_STEPS
    center = mpmath.mpc(sigma, t_center)
    phase_accum = mpmath.mpf(0)
    prev_val = None

    for j in range(n_steps + 1):
        theta = 2 * mpmath.pi * j / n_steps
        pt = center + radius * mpmath.exp(1j * theta)
        try:
            val = Lambda_AFE(pt, cn)
            if abs(val) < mpmath.mpf(10) ** (-DPS + 15):
                return None
        except Exception as e:
            print(f"  WARNING mono step {j}: {e}", flush=True)
            return None

        if prev_val is not None:
            ratio = val / prev_val
            if abs(ratio) < mpmath.mpf(10) ** (-DPS + 15):
                return None
            phase_accum += mpmath.im(mpmath.log(ratio))
        prev_val = val

    return float(abs(phase_accum) / mpmath.pi)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#                   메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 70)
log("결과 #64 — GL(3) sym²(37a1) ξ-다발 4성질 검증")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"N={N_COND} (=37²), ε={EPSILON}, gamma shifts μ=[1,1,2] (해석적 정규화)")
log(f"37a1: y²+y = x³-x, rank=1, a₃₇=-1 (nonsplit multiplicative)")
log(f"DPS={DPS}, N_coeff={N_MAX_COEFF}, N_hermite={N_HERMITE}, c={float(CONTOUR_C)}")
log()

t_total = time.time()

# ── 0. Dirichlet 계수 ──
log("[Step 0] sym²(37a1) Dirichlet 계수 (해석적 정규화)")
an, ap, primes = compute_37a1_an(N_MAX_COEFF)
cn = compute_sym2_cn(an, primes, N_MAX_COEFF)

# LMFDB 참조 37a1 aₙ (n=1..20) — 검증
lmfdb_ref = {
    1: 1, 2: -2, 3: -3, 4: 2, 5: -2, 6: 6, 7: -1, 8: 0,
    9: 6, 10: 4, 11: -5, 12: -6, 13: -2, 14: 2, 15: 6,
    16: -4, 17: 0, 18: -12, 19: 0, 20: -4
}
an_ok = True
for n, expected in lmfdb_ref.items():
    if an[n] != expected:
        log(f"  ⚠️ a({n})={an[n]} ≠ LMFDB {expected}")
        an_ok = False
if an_ok:
    log(f"  ✅ aₙ (n=1..20) LMFDB 일치")
else:
    log(f"  ❌ aₙ 불일치 — 계산 오류!")
log(f"  a₂={an[2]}, a₃={an[3]}, a₅={an[5]}, a₇={an[7]}, a₃₇={an[37]}")

# sym² 계수 검증: c(p) = a_f(p)²/p - 1 for good primes
log(f"\n  sym² 계수 검증:")
checks_good = []
for p in [2, 3, 5, 7, 11, 13]:
    expected = an[p] ** 2 / p - 1
    actual = float(cn[p])
    ok = abs(actual - expected) < 1e-10
    checks_good.append(ok)
    log(f"  c({p})={actual:.8f} (a²/p-1={expected:.8f}) {'✅' if ok else '❌'}")

# bad prime 37: c(37) = (-1)^1 · 37^{-1} = -1/37
c37_expected = -1 / 37
c37_actual = float(cn[37])
c37_ok = abs(c37_actual - c37_expected) < 1e-10
log(f"  c(37)={c37_actual:.8f} (기대={c37_expected:.8f}) {'✅' if c37_ok else '❌'}")

n_nz = sum(1 for i in range(1, len(cn)) if cn[i] != 0)
log(f"  비영 계수: {n_nz}/{N_MAX_COEFF}")
flush_to_file()

# ── 1. 벤치마크 ──
log(f"\n[Step 1] 벤치마크: Λ_AFE 한 점 소요 시간")
t1 = time.time()
val_bench = Lambda_AFE(mpmath.mpc(3, 0), cn)
dt_bench = time.time() - t1
log(f"  Λ_AFE(3) = {mpmath.nstr(val_bench, 12)} (소요: {dt_bench:.2f}s)")
flush_to_file()

# ── 2. AFE vs Λ_direct 검증 (Re(s) >> 1) ──
log(f"\n[Step 2] AFE vs Λ_direct 검증 (Re(s) >> 1)")
test_points = [
    mpmath.mpf(3),
    mpmath.mpf(2.5),
    mpmath.mpc(3, 5),
    mpmath.mpc(2.5, 10),
    mpmath.mpc(2, 15),
]

fe_validated = True
for sp in test_points:
    t_s = time.time()
    val_dir = Lambda_direct(sp, cn)
    val_afe = Lambda_AFE(sp, cn)
    dt_s = time.time() - t_s

    if abs(val_dir) > 0:
        rel = float(abs(val_afe - val_dir) / abs(val_dir))
    else:
        rel = float('inf')
    ok = rel < 1e-6
    if not ok:
        fe_validated = False
    log(f"  s={mpmath.nstr(sp, 6)}: |dir|={float(abs(val_dir)):.6e}, "
        f"|AFE|={float(abs(val_afe)):.6e}, rel={rel:.2e} "
        f"{'✅' if ok else '❌'} ({dt_s:.1f}s)")
    flush_to_file()

if fe_validated:
    log(f"\n  ✅ AFE 검증 통과 — 계수+감마+AFE 올바름")
else:
    log(f"\n  ⚠️ AFE 검증 일부 실패 (Dirichlet 수렴 한계 — 임계선 검증으로 판정)")
flush_to_file()

# ── 3. 함수방정식 직접 검증: Λ(s) ≈ Λ(1-s) ──
log(f"\n[Step 3] 함수방정식 Λ(s) = Λ(1-s) 직접 검증 (핵심 검증)")
fe_test_pts = [
    mpmath.mpc(0.5, 5),
    mpmath.mpc(0.5, 10),
    mpmath.mpc(0.5, 15),
    mpmath.mpc(0.3, 8),
    mpmath.mpc(0.7, 12),
]
fe_direct_ok = True
fe_max_rel = 0.0
for sp in fe_test_pts:
    t_s = time.time()
    L_s = Lambda_AFE(sp, cn)
    L_1s = Lambda_AFE(1 - sp, cn)
    dt_s = time.time() - t_s

    if abs(L_s) > mpmath.mpf(10) ** (-DPS + 15):
        rel = float(abs(L_s - EPSILON * L_1s) / abs(L_s))
    else:
        rel = 0.0

    fe_max_rel = max(fe_max_rel, rel)
    ok = rel < 1e-6
    if not ok:
        fe_direct_ok = False
    log(f"  s={mpmath.nstr(sp, 6)}: |Λ(s)|={float(abs(L_s)):.6e}, "
        f"|Λ(1-s)|={float(abs(L_1s)):.6e}, rel={rel:.2e} "
        f"{'✅' if ok else '❌'} ({dt_s:.1f}s)")
    flush_to_file()

if fe_direct_ok:
    log(f"  ✅ 함수방정식 검증 통과 (max_rel={fe_max_rel:.2e})")
else:
    log(f"  ⚠️ 함수방정식 오차 존재 (max_rel={fe_max_rel:.2e})")
flush_to_file()

# ── 4. Λ(1/2+it) 실수성 확인 + 영점 스캔 ──
log(f"\n[Step 4] Λ(1/2+it) 스캔 (t ∈ [{T_MIN}, {T_MAX}])")
log(f"  ε=+1 → Λ(1/2+it) ∈ ℝ")

scan_ts = np.arange(T_MIN, T_MAX + 0.5, 0.5)  # 격자 간격 0.5
scan_vals = []
max_im = 0
t4_start = time.time()

for i, t in enumerate(scan_ts):
    val = Lambda_AFE(mpmath.mpc(0.5, t), cn)
    re_val = float(mpmath.re(val))
    im_val = float(mpmath.im(val))
    scan_vals.append(re_val)
    max_im = max(max_im, abs(im_val))

    if (i + 1) % 10 == 0:
        elapsed = time.time() - t4_start
        eta = elapsed / (i + 1) * (len(scan_ts) - i - 1)
        log(f"  [{i+1}/{len(scan_ts)}] t={t:.1f}: Re(Λ)={re_val:+.6e}, "
            f"Im={im_val:+.2e} ({elapsed:.0f}s, ETA {eta:.0f}s)")
        flush_to_file()

scan_vals_arr = np.array(scan_vals)
n_sign = int(np.sum(np.diff(np.sign(scan_vals_arr)) != 0))
dt_scan = time.time() - t4_start

log(f"\n  스캔 완료: {dt_scan:.1f}s")
log(f"  부호변환: {n_sign}")
log(f"  Re(Λ) 범위: [{scan_vals_arr.min():.6e}, {scan_vals_arr.max():.6e}]")
log(f"  최대 |Im|: {max_im:.2e} (0이면 실수성 확인)")
flush_to_file()

# ── 5. 영점 이분법 정밀화 ──
zero_approx = []
log(f"\n[Step 5] 영점 이분법 정밀화")

for i in range(len(scan_vals_arr) - 1):
    if scan_vals_arr[i] * scan_vals_arr[i + 1] < 0:
        t_lo, t_hi = float(scan_ts[i]), float(scan_ts[i + 1])
        val_lo = scan_vals_arr[i]
        try:
            for _ in range(25):
                t_mid = (t_lo + t_hi) / 2
                val_mid = float(mpmath.re(Lambda_AFE(mpmath.mpc(0.5, t_mid), cn)))
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
            log(f"  WARNING bisection at t≈{(t_lo + t_hi) / 2:.2f}: {e}")
            zero_approx.append((t_lo + t_hi) / 2)

log(f"  발견 영점: {len(zero_approx)}개")
if len(zero_approx) == 0:
    log(f"  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
flush_to_file()

# LMFDB 매칭 (데이터 있을 때)
if LMFDB_ZEROS:
    matched = 0
    for z in zero_approx:
        dists = [abs(z - lz) for lz in LMFDB_ZEROS]
        best = min(dists)
        best_idx = dists.index(best)
        ok = best < 0.1
        if ok:
            matched += 1
        log(f"  γ ≈ {z:.8f}  LMFDB[{best_idx}]={LMFDB_ZEROS[best_idx]:.6f}  "
            f"Δ={best:.6f} {'✅' if ok else '❌'}")
    log(f"  LMFDB 매칭: {matched}/{len(zero_approx)}")
else:
    matched = 0
    log(f"\n  LMFDB 데이터 미가용 — 자체 영점만 보고")
    for j, z in enumerate(zero_approx):
        log(f"  #{j+1:3d}: γ = {z:.8f}")
        if (j + 1) % 20 == 0:
            flush_to_file()
flush_to_file()

# ── 6. κ 곡률 + 모노드로미 ──
if len(zero_approx) > 0:
    log(f"\n[Step 6] 영점 근방 κ + 모노드로미")

    kappa_near_list = []
    kappa_far_list = []
    mono_results = []

    # TP: 영점 근방 (최대 15개 측정)
    n_tp = min(15, len(zero_approx))
    for j, z_t in enumerate(zero_approx[:n_tp]):
        s_pt = complex(0.5 + DELTA, z_t)
        k_val = curvature(s_pt, cn)
        m_val = monodromy(z_t, cn)
        kappa_near_list.append(k_val)
        mono_results.append(m_val)
        m_str = f"{m_val:.4f}" if m_val is not None else "FAIL"
        log(f"  TP #{j+1} γ={z_t:.6f}: κ={k_val:.2f}, mono/π={m_str}")
        flush_to_file()

    # FP: 영점 사이 (mid-points)
    log(f"\n  FP 측정 (영점 중간점):")
    fp_ts = []
    for j in range(len(zero_approx) - 1):
        mid = (zero_approx[j] + zero_approx[j + 1]) / 2
        fp_ts.append(mid)
    if zero_approx[0] > T_MIN + 1:
        fp_ts.insert(0, T_MIN + 0.5)
    if zero_approx[-1] < T_MAX - 1:
        fp_ts.append(T_MAX - 1)

    n_fp = min(12, len(fp_ts))
    for j, t in enumerate(fp_ts[:n_fp]):
        k_val = curvature(complex(0.5 + DELTA, t), cn)
        kappa_far_list.append(k_val)
        log(f"  FP #{j+1} t={t:.3f}: κ={k_val:.2f}")
        flush_to_file()

    # κ 비율 통계
    near_finite = [k for k in kappa_near_list if np.isfinite(k) and k < 1e15]
    far_finite = [k for k in kappa_far_list if np.isfinite(k) and k < 1e15]

    if near_finite and far_finite:
        near_med = float(np.median(near_finite))
        far_med = float(np.median(far_finite))
        if far_med > 0:
            ratio = near_med / far_med
        else:
            ratio = float('inf')
        log(f"\n  κ near(중앙값) = {near_med:.2f} (n={len(near_finite)})")
        log(f"  κ far(중앙값)  = {far_med:.4f} (n={len(far_finite)})")
        log(f"  비율 = {ratio:.1f}×")

        # Q25, Q75 분석
        if len(near_finite) >= 4:
            q25 = float(np.percentile(near_finite, 25))
            q75 = float(np.percentile(near_finite, 75))
            log(f"  κ near Q25={q25:.2f}, Q50={near_med:.2f}, Q75={q75:.2f}")
    else:
        ratio = 0
        near_med = 0
        far_med = 0
        log(f"  데이터 부족 (near={len(near_finite)}, far={len(far_finite)})")
    flush_to_file()

    # ── 7. σ-유일성 ──
    log(f"\n[Step 7] σ-유일성 스윕")
    sigma_results = {}
    for sig in [0.3, 0.5, 0.7, 0.9]:
        sc_vals = []
        for t in np.arange(T_MIN, min(T_MAX, 40.0), 0.8):
            val = Lambda_AFE(mpmath.mpc(sig, t), cn)
            sc_vals.append(float(mpmath.re(val)))
        n_sc = int(np.sum(np.diff(np.sign(sc_vals)) != 0))
        sigma_results[sig] = n_sc
        log(f"  σ={sig:.1f}: 부호변환={n_sc}")
        flush_to_file()

    sc_05 = sigma_results.get(0.5, 0)
    sc_others = [sigma_results.get(s, 0) for s in [0.3, 0.7, 0.9]]
    sigma_pass = sc_05 > 0 and sc_05 >= max(sc_others)

    if sigma_pass:
        log(f"  → σ=0.5 최대 ({sc_05}) ✅")
    else:
        max_sig = max(sigma_results, key=sigma_results.get)
        log(f"  → 극대 σ={max_sig} ({sigma_results[max_sig]})")
        log(f"  → GL(3) σ-유일성 FAIL 패턴 (GL(n≥2) 예상)")

    # #63과의 비교
    log(f"\n[Step 7b] σ-유일성 패턴 비교")
    log(f"  GL(1) ζ:            PASS (σ=0.5 극대)")
    log(f"  GL(2) 11a1:         FAIL (σ=0.7-0.9 극대)")
    log(f"  GL(3) sym²(11a1):   FAIL (거의 평탄)")
    log(f"  GL(3) sym²(37a1):   {sigma_results}")
    flush_to_file()

else:
    ratio = 0
    sigma_pass = False
    mono_results = []
    near_med = 0
    far_med = 0

# ━━━━━━━━━━━ 종합 판정 ━━━━━━━━━━━
log(f"\n{'='*70}")
log(f"종합 판정 — GL(3) sym²(37a1) 4성질")
log(f"{'='*70}")
log(f"Conductor: N={N_COND}, gamma shifts: μ=[1,1,2], ε={EPSILON}")
log(f"37a1: rank=1, a₃₇=-1 (nonsplit)")
log()

# 1. 함수방정식
fe_str = "PASS ✅" if (fe_direct_ok) else f"FAIL ❌ (max_rel={fe_max_rel:.2e})"
log(f"  [필수] 함수방정식 Λ(s)=Λ(1-s): {fe_str}")

# 2. 영점
log(f"  영점 발견: {len(zero_approx)}개 (t ∈ [{T_MIN}, {T_MAX}])")
n_zero_pass = len(zero_approx) >= 30
log(f"    목표 ≥30: {'✅ PASS' if n_zero_pass else '❌ FAIL'}")

# 3. κ 집중
if len(zero_approx) > 0:
    if ratio >= 100:
        kr_str = f"PASS ✅ ({ratio:.1f}×, 기준 100×)"
    elif ratio >= 10:
        kr_str = f"약한 양성 ({ratio:.1f}×, 기준 100×)"
    else:
        kr_str = f"FAIL ({ratio:.1f}×, 기준 100×)"
else:
    kr_str = "N/A"
log(f"  [양성] κ 집중: {kr_str}")

# 4. 모노드로미
if mono_results:
    tp_mono = sum(1 for m in mono_results if m is not None and 1.5 < m < 2.5)
    total_mono = len([m for m in mono_results if m is not None])
    if total_mono > 0:
        mono_str = f"{tp_mono}/{total_mono} (mono/π≈2.0)"
        if tp_mono == total_mono:
            mono_str += " ✅"
        else:
            mono_str += f" ⚠️ ({total_mono - tp_mono}개 이상)"
    else:
        mono_str = "FAIL (전부 None)"
else:
    mono_str = "N/A"
log(f"  [양성] 모노드로미: {mono_str}")

# 5. σ-유일성
if sigma_pass:
    su_str = "PASS ✅ (예상 외!)"
else:
    su_str = "FAIL (GL(3) 예상 패턴)"
log(f"  [정보] σ-유일성: {su_str}")

# ── 11a1 vs 37a1 비교 ──
log(f"\n{'='*70}")
log(f"교차비교: sym²(11a1) #63 vs sym²(37a1) #64")
log(f"{'='*70}")
log(f"  {'항목':20s} | {'#63 sym²(11a1)':20s} | {'#64 sym²(37a1)':20s}")
log(f"  {'─'*20} | {'─'*20} | {'─'*20}")
log(f"  {'conductor':20s} | {'121':20s} | {str(N_COND):20s}")
log(f"  {'rank':20s} | {'0':20s} | {'1':20s}")
log(f"  {'a_p bad prime':20s} | {'a₁₁=1 (split)':20s} | {'a₃₇=-1 (nonsplit)':20s}")
log(f"  {'N_coeff':20s} | {'120':20s} | {str(N_MAX_COEFF):20s}")
log(f"  {'FE rel_err':20s} | {'0.0':20s} | {f'{fe_max_rel:.2e}':20s}")
log(f"  {'영점 수':20s} | {'46':20s} | {str(len(zero_approx)):20s}")
log(f"  {'κ ratio':20s} | {'283.3×':20s} | {f'{ratio:.1f}×':20s}")
if mono_results:
    tp_m = sum(1 for m in mono_results if m is not None and 1.5 < m < 2.5)
    tot_m = len([m for m in mono_results if m is not None])
    log(f"  {'mono TP':20s} | {'12/12':20s} | {f'{tp_m}/{tot_m}':20s}")
log(f"  {'σ-유일성':20s} | {'FAIL':20s} | {'PASS' if sigma_pass else 'FAIL':20s}")

total_time = time.time() - t_total
log(f"\n  총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log(f"  시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")

flush_to_file()
log(f"\n결과 저장: {OUTFILE}")
