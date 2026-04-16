#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #69 — Maass form L-함수 4성질 검증
=============================================================================
대상: 첫 번째 Maass cusp form on SL(2,Z)
  - 스펙트럼 파라미터 R = 9.53369526135355755...
  - 대칭: ODD (symmetryclass=1)  ← 수학자 "even" 주장 정정
  - Conductor N=1
  - Root number ε = +1 (SL(2,Z) 레벨 1)

감마 인자 (ODD Maass form):
  γ(s) = Γ_R(s+1+iR) · Γ_R(s+1-iR)
       = π^{-(s+1)} Γ((s+1+iR)/2) Γ((s+1-iR)/2)

Λ(s) = γ(s) · L(s)    [N=1 이므로 Q^s = 1]
함수방정식: Λ(s) = Λ(1-s)

Fourier 계수: Booker-Strömbergsson-Venkatesh 데이터 (1000자리 정밀도, 455개)
출처: https://www2.math.uu.se/~astrombe/emaass/psl2z/coeff9

4성질:
  1. σ=1/2 유일성
  2. κ-피크 (κ_near 측정 → GL(3) 1125와 비교)
  3. 모노드로미 (mono/π = ±2)
  4. FE 잔차

결과 파일: results/maass_4property_69.txt
=============================================================================
"""

import sys, os, time, re
import numpy as np
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

DPS = 50
mpmath.mp.dps = DPS

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "maass_4property_69.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

# ━━━━━━━━━━━ 로깅 ━━━━━━━━━━━
lines = []
def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))
def flush():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━━━━━━ 파라미터 ━━━━━━━━━━━
R_SPECTRAL = mpmath.mpf("9.53369526135355755434423523592877032382125639510725198237579046")
EPSILON = 1  # root number (SL(2,Z) 레벨 1)
N_COND = 1   # conductor

# Gauss-Hermite 직교
N_HERMITE = 50
_herm_nodes, _herm_weights = np.polynomial.hermite.hermgauss(N_HERMITE)
HERM_NODES = [mpmath.mpf(float(x)) for x in _herm_nodes]
HERM_WEIGHTS = [mpmath.mpf(float(w)) for w in _herm_weights]

CONTOUR_C = mpmath.mpf(2)  # AFE contour shift

# κ 스윕 파라미터
T_MIN, T_MAX = 2.0, 50.0
DT_SCAN = 0.5      # 영점 스캔 간격
DELTA = 0.03        # σ 오프셋
MONO_R = 0.4        # 모노드로미 반지름
MONO_N = 48         # 모노드로미 단계

# ━━━━━━━━━━━ Fourier 계수 로드 ━━━━━━━━━━━
def load_coefficients():
    """Strömbergsson 데이터 파일에서 Maass form 계수 로드"""
    coeff_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "maass_coeff9_raw.txt")

    if not os.path.exists(coeff_file):
        log(f"⚠️ 계수 파일 없음: {coeff_file}")
        log("  Booker-Strömbergsson-Venkatesh 데이터 다운로드 시도...")
        try:
            import requests
            r = requests.get(
                'http://www2.math.uu.se/~astrombe/emaass/psl2z/coeff9',
                timeout=30, verify=False
            )
            with open(coeff_file, 'w') as f:
                f.write(r.text)
            log(f"  ✅ 다운로드 성공 ({len(r.text)} chars)")
        except Exception as e:
            log(f"  ❌ 다운로드 실패: {e}")
            return None

    with open(coeff_file, 'r') as f:
        text = f.read()

    # Parse R value
    R_match = re.search(r'R=(\d+\.\d+)', text)
    R_parsed = R_match.group(1)[:60] if R_match else "UNKNOWN"

    # Parse symmetry class
    sym_match = re.search(r'symmetryclass=(\d+)', text)
    sym_class = int(sym_match.group(1)) if sym_match else -1

    # Parse coefficients C[n]
    coeff_matches = re.findall(r'C\[(\d+)\]=\s*([-+]?\d+\.\d+)', text)
    cn_dict = {}
    for n_str, v_str in coeff_matches:
        n = int(n_str)
        # 50자리까지만 사용 (DPS=50)
        cn_dict[n] = mpmath.mpf(v_str[:55])

    return R_parsed, sym_class, cn_dict


# ━━━━━━━━━━━ 감마 인자 (ODD Maass form) ━━━━━━━━━━━

def gamma_factor_maass(s):
    """
    ODD Maass form 감마 인자:
    γ(s) = π^{-(s+1)} · Γ((s+1+iR)/2) · Γ((s+1-iR)/2)
    """
    iR = mpmath.mpc(0, R_SPECTRAL)
    return (mpmath.power(mpmath.pi, -(s + 1))
            * mpmath.gamma((s + 1 + iR) / 2)
            * mpmath.gamma((s + 1 - iR) / 2))


def dirichlet_series(w, cn_dict, N_max):
    """D(w) = Σ_{n=1}^{N_max} a(n)/n^w"""
    total = mpmath.mpc(0)
    for n in range(1, N_max + 1):
        if n not in cn_dict or cn_dict[n] == 0:
            continue
        total += cn_dict[n] * mpmath.power(n, -w)
    return total


def Lambda_AFE(s, cn_dict, N_max):
    """
    Λ(s) via Mellin-Barnes AFE + Gauss-Hermite.

    Λ(s) = (exp(c²)/(2π)) Σ_k w_k [A(s+c+iy_k) D(s+c+iy_k) + ε A(1-s+c+iy_k) D(1-s+c+iy_k)] exp(2icy_k) / (c+iy_k)

    여기서 A(u) = γ(u), D(u) = L(u) = Σ a(n)/n^u
    N=1이므로 N^{u/2} = 1 (conductor 기여 없음)
    """
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

        exp_phase = mpmath.exp(2 * mpmath.mpc(0, 1) * c * v_k)

        # A_k(s): gamma at s+w_shift
        sw = s_mp + w_shift
        gamma_sw = gamma_factor_maass(sw)
        # N=1 이므로 N^{sw/2} = 1
        A_s = gamma_sw * exp_phase / w_shift
        D_s = dirichlet_series(sw, cn_dict, N_max)

        # A_k(1-s)
        s1w = s1_mp + w_shift
        gamma_s1w = gamma_factor_maass(s1w)
        A_1s = gamma_s1w * exp_phase / w_shift
        D_1s = dirichlet_series(s1w, cn_dict, N_max)

        total += wt_k * (A_s * D_s + EPSILON * A_1s * D_1s)

    return prefactor * total


def curvature(t, cn_dict, N_max, delta=DELTA, h=1e-6):
    """κ(s) = |Λ'/Λ|² at s = 0.5+δ+it"""
    s = mpmath.mpc(0.5 + delta, t)
    h_mp = mpmath.mpf(h)
    try:
        L0 = Lambda_AFE(s, cn_dict, N_max)
        if abs(L0) < mpmath.mpf(10) ** (-DPS + 15):
            return float('inf')
        Lp = Lambda_AFE(s + h_mp, cn_dict, N_max)
        Lm = Lambda_AFE(s - h_mp, cn_dict, N_max)
        conn = (Lp - Lm) / (2 * h_mp * L0)
        k = float(abs(conn) ** 2)
        return k if np.isfinite(k) else float('inf')
    except Exception as e:
        print(f"  WARNING curvature at t={t}: {e}", flush=True)
        return 0.0


def monodromy(t_center, cn_dict, N_max, radius=MONO_R, n_steps=MONO_N):
    """폐곡선 적분으로 모노드로미 측정"""
    center = mpmath.mpc(0.5, t_center)
    phase_accum = mpmath.mpf(0)
    prev_val = None

    for j in range(n_steps + 1):
        theta = 2 * mpmath.pi * j / n_steps
        pt = center + radius * mpmath.exp(1j * theta)
        try:
            val = Lambda_AFE(pt, cn_dict, N_max)
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


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#                   메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

log("=" * 70)
log("결과 #69 — Maass form L-함수 4성질 검증")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"DPS={DPS}, N_hermite={N_HERMITE}, c={float(CONTOUR_C)}")
log(f"R = {mpmath.nstr(R_SPECTRAL, 30)}")
log(f"대칭: ODD (symmetryclass=1)")
log(f"감마: π^{{-(s+1)}} Γ((s+1+iR)/2) Γ((s+1-iR)/2)")
log(f"ε = {EPSILON}, N = {N_COND}")
log()

t_total = time.time()

# ━━━━━━ Step 0: 계수 로드 ━━━━━━

log("[Step 0] Fourier 계수 로드 (Booker-Strömbergsson-Venkatesh)")
result = load_coefficients()
if result is None:
    log("❌ 계수 로드 실패 — 중단")
    flush()
    sys.exit(1)

R_parsed, sym_class, cn_dict = result
N_max = max(cn_dict.keys())

log(f"  R (파일) = {R_parsed}...")
log(f"  symmetryclass = {sym_class} ({'ODD' if sym_class == 1 else 'EVEN' if sym_class == 0 else '???'})")
log(f"  계수 수: {len(cn_dict)} (n=1..{N_max})")
log(f"  a(1) = {float(cn_dict[1]):.6f}")
log(f"  a(2) = {float(cn_dict[2]):.15f}")
log(f"  a(3) = {float(cn_dict[3]):.15f}")
log(f"  a(5) = {float(cn_dict[5]):.15f}")
log(f"  a(7) = {float(cn_dict[7]):.15f}")
log(f"  a(11) = {float(cn_dict[11]):.15f}")

# Hecke 곱셈성 검증: a(6) = a(2)*a(3)
a6_check = cn_dict[2] * cn_dict[3]
a6_diff = abs(float(cn_dict[6]) - float(a6_check))
log(f"  Hecke 곱셈성: a(6)={float(cn_dict[6]):.10f}, a(2)*a(3)={float(a6_check):.10f}, diff={a6_diff:.2e}")

# a(4) = a(2)^2 - 1 (Hecke 재귀)
a4_check = cn_dict[2]**2 - 1
a4_diff = abs(float(cn_dict[4]) - float(a4_check))
log(f"  Hecke 재귀: a(4)={float(cn_dict[4]):.10f}, a(2)²-1={float(a4_check):.10f}, diff={a4_diff:.2e}")

# a(9) = a(3)^2 - 1
a9_check = cn_dict[3]**2 - 1
a9_diff = abs(float(cn_dict[9]) - float(a9_check))
log(f"  Hecke 재귀: a(9)={float(cn_dict[9]):.10f}, a(3)²-1={float(a9_check):.10f}, diff={a9_diff:.2e}")

if sym_class != 1:
    log(f"  ⚠️ symmetryclass={sym_class}, 예상: 1 (ODD)")

log()
flush()

# ━━━━━━ Step 1: 함수방정식 검증 ━━━━━━

log("[Step 1] 함수방정식 Λ(s) = εΛ(1-s) on critical line")
fe_pts = [
    mpmath.mpc(0.5, 5),
    mpmath.mpc(0.5, 10),
    mpmath.mpc(0.5, 15),
    mpmath.mpc(0.5, 25),
]
fe_ok = True
fe_max_rel = 0.0

for sp in fe_pts:
    t_s = time.time()
    L_s = Lambda_AFE(sp, cn_dict, N_max)
    L_1s = Lambda_AFE(1 - sp, cn_dict, N_max)
    dt_s = time.time() - t_s

    if abs(L_s) > mpmath.mpf(10) ** (-DPS + 15):
        rel = float(abs(L_s - EPSILON * L_1s) / abs(L_s))
    else:
        rel = 0.0  # 영점 근방

    fe_max_rel = max(fe_max_rel, rel)
    ok = rel < 1e-6
    if not ok:
        fe_ok = False
    log(f"  s={mpmath.nstr(sp, 6)}: |Λ|={float(abs(L_s)):.4e}, rel={rel:.2e} {'✅' if ok else '❌'} ({dt_s:.1f}s)")
    flush()

if fe_ok:
    log(f"  ✅ 함수방정식 통과 (max_rel={fe_max_rel:.2e})")
else:
    log(f"  ❌ 함수방정식 실패 (max_rel={fe_max_rel:.2e})")
    log(f"  경고: ODD 감마 인자 사용 중. 재확인 필요.")
log()
flush()

# FE 실패 시에도 계속 진행 (진단 목적)
# 만약 실패하면 EVEN 감마 인자로 재시도가 필요할 수 있음

# ━━━━━━ Step 2: Λ(1/2+it) 스캔 — 영점 탐색 ━━━━━━

log(f"[Step 2] Λ(1/2+it) 스캔 (t ∈ [{T_MIN}, {T_MAX}])")
scan_ts = np.arange(T_MIN, T_MAX + DT_SCAN, DT_SCAN)
scan_vals = []
t_scan = time.time()

for i, t in enumerate(scan_ts):
    val = Lambda_AFE(mpmath.mpc(0.5, t), cn_dict, N_max)
    re_val = float(mpmath.re(val))
    scan_vals.append(re_val)

    if (i + 1) % 20 == 0:
        elapsed = time.time() - t_scan
        eta = elapsed / (i + 1) * (len(scan_ts) - i - 1)
        log(f"  [{i+1}/{len(scan_ts)}] t={t:.1f}: Re(Λ)={re_val:+.4e} ({elapsed:.0f}s, ETA {eta:.0f}s)")
        flush()

scan_arr = np.array(scan_vals)
dt_scan = time.time() - t_scan
log(f"  스캔 완료: {dt_scan:.1f}s ({dt_scan/60:.1f}분)")

# ━━━━━━ Step 3: 영점 이분법 정밀화 ━━━━━━

zero_approx = []
for i in range(len(scan_arr) - 1):
    if scan_arr[i] * scan_arr[i + 1] < 0:
        t_lo, t_hi = float(scan_ts[i]), float(scan_ts[i + 1])
        val_lo = scan_arr[i]
        try:
            for _ in range(20):  # 20회 이분법
                t_mid = (t_lo + t_hi) / 2
                val_mid = float(mpmath.re(Lambda_AFE(mpmath.mpc(0.5, t_mid), cn_dict, N_max)))
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
            log(f"  WARNING bisection at t≈{(t_lo+t_hi)/2:.3f}: {e}")
            zero_approx.append((t_lo + t_hi) / 2)

n_zeros = len(zero_approx)
log(f"\n  발견 영점: {n_zeros}개")
if n_zeros == 0:
    log(f"  ⚠️ 영점 0개 — 탐색 로직 점검 필요")
    flush()
    # 영점이 없어도 진단을 위해 계속
else:
    for j, z in enumerate(zero_approx[:15]):
        log(f"    #{j+1}: γ = {z:.10f}")
    if n_zeros > 15:
        log(f"    ... ({n_zeros - 15}개 추가)")

log()
flush()

# ━━━━━━ Step 4: κ near/far 측정 + 모노드로미 ━━━━━━

if n_zeros >= 2:
    log(f"[Step 3] κ near/far 비율 + 모노드로미")
    n_tp = min(8, n_zeros)
    kappa_near = []
    mono_results = []

    for j, z_t in enumerate(zero_approx[:n_tp]):
        t_k = time.time()
        k_val = curvature(z_t, cn_dict, N_max, delta=DELTA)
        m_val = monodromy(z_t, cn_dict, N_max)
        dt_k = time.time() - t_k
        kappa_near.append(k_val)
        mono_results.append(m_val)
        m_str = f"{m_val:.4f}" if m_val is not None else "FAIL"
        log(f"    TP #{j+1} γ={z_t:.8f}: κ={k_val:.2f}, mono/π={m_str} ({dt_k:.1f}s)")
        flush()

    # FP (영점 사이 중간점)
    kappa_far = []
    fp_ts = []
    for j in range(len(zero_approx) - 1):
        mid = (zero_approx[j] + zero_approx[j + 1]) / 2
        fp_ts.append(mid)
    if zero_approx[0] > T_MIN + 1:
        fp_ts.insert(0, T_MIN + 0.5)

    n_fp = min(6, len(fp_ts))
    for j, t in enumerate(fp_ts[:n_fp]):
        t_k = time.time()
        k_val = curvature(t, cn_dict, N_max, delta=DELTA)
        dt_k = time.time() - t_k
        kappa_far.append(k_val)
        log(f"    FP #{j+1} t={t:.6f}: κ={k_val:.4f} ({dt_k:.1f}s)")
        flush()

    # κ 통계
    near_fin = [k for k in kappa_near if np.isfinite(k) and k < 1e15 and k > 0]
    far_fin = [k for k in kappa_far if np.isfinite(k) and k < 1e15 and k > 0]

    if near_fin and far_fin:
        near_med = float(np.median(near_fin))
        far_med = float(np.median(far_fin))
        ratio = near_med / far_med if far_med > 0 else float('inf')
        near_cv = float(np.std(near_fin) / np.mean(near_fin) * 100) if np.mean(near_fin) > 0 else 999
    else:
        near_med = 0
        far_med = 0
        ratio = 0
        near_cv = 999

    log(f"\n  κ near median = {near_med:.2f} (n={len(near_fin)}, CV={near_cv:.1f}%)")
    log(f"  κ far median  = {far_med:.4f} (n={len(far_fin)})")
    log(f"  κ ratio = {ratio:.1f}×")

    # 모노드로미 요약
    mono_pass = sum(1 for m in mono_results if m is not None and 1.5 < m < 2.5)
    total_mono = len([m for m in mono_results if m is not None])
    log(f"  모노드로미: {mono_pass}/{total_mono} (mono/π≈2.0)")

    # 상세 모노드로미
    for j, m in enumerate(mono_results):
        m_str = f"{m:.4f}" if m is not None else "FAIL"
        log(f"    TP #{j+1}: mono/π = {m_str}")

    log()
    flush()

# ━━━━━━ Step 5: σ-유일성 테스트 ━━━━━━

log("[Step 4] σ-유일성 테스트 (7 σ × 영점)")
if n_zeros >= 2:
    sigma_vals = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
    test_zeros = zero_approx[:min(3, n_zeros)]

    for z_t in test_zeros:
        log(f"\n  영점 γ={z_t:.6f}:")
        kappas_by_sigma = []
        for sigma in sigma_vals:
            k_val = curvature(z_t, cn_dict, N_max, delta=sigma - 0.5)
            kappas_by_sigma.append(k_val)
            log(f"    σ={sigma:.2f}: κ={k_val:.2f}")
            flush()

        # σ=0.5에서의 κ가 최대인지 확인
        idx_05 = sigma_vals.index(0.5)
        k_05 = kappas_by_sigma[idx_05]
        is_max = all(k_05 >= k for k in kappas_by_sigma if np.isfinite(k))
        log(f"    σ=0.5 최대? {'✅ YES' if is_max else '❌ NO'}")

    log()
    flush()
else:
    log("  ⚠️ 영점 부족 — σ-유일성 스킵")
    flush()

# ━━━━━━ 최종 요약 ━━━━━━

total_time = time.time() - t_total
log("=" * 70)
log("최종 요약")
log("=" * 70)
log(f"대상: Maass cusp form on SL(2,Z), R={mpmath.nstr(R_SPECTRAL, 15)}")
log(f"대칭: ODD (symmetryclass=1)")
log(f"감마: π^{{-(s+1)}} Γ((s+1+iR)/2) Γ((s+1-iR)/2)")
log(f"ε = {EPSILON}, N = {N_COND}")
log(f"계수: {len(cn_dict)}개 (Booker-Strömbergsson-Venkatesh)")
log(f"")
log(f"1. 함수방정식: {'✅ PASS' if fe_ok else '❌ FAIL'} (max_rel={fe_max_rel:.2e})")
log(f"2. 영점 발견: {n_zeros}개 (t ∈ [{T_MIN}, {T_MAX}])")

if n_zeros >= 2 and near_fin:
    log(f"3. κ_near median = {near_med:.2f} (CV={near_cv:.1f}%)")
    log(f"   κ_far median  = {far_med:.4f}")
    log(f"   κ ratio = {ratio:.1f}×")
    log(f"   GL(3) κ_near = 1125 → {'유사' if abs(near_med - 1125) / 1125 < 0.1 else '상이'} (차이 {abs(near_med-1125)/1125*100:.1f}%)")
    log(f"4. 모노드로미: {mono_pass}/{total_mono} (mono/π≈2.0)")
else:
    log(f"3. κ_near: 영점 부족으로 측정 불가")
    log(f"4. 모노드로미: 영점 부족")

log(f"")
log(f"총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log(f"")

# GL(3) 비교
if n_zeros >= 2 and near_fin:
    log("GL(2) Maass vs GL(3) sym² 비교:")
    log(f"  GL(3) sym² κ_near ≈ 1125 (6곡선, CV=0.15%)")
    log(f"  GL(2) Maass κ_near = {near_med:.2f} (CV={near_cv:.1f}%)")
    if abs(near_med - 1125) / 1125 < 0.05:
        log(f"  → κ_near ≈ 1125 보편상수 가설 지지 ✅")
    elif abs(near_med - 1125) / 1125 < 0.15:
        log(f"  → 유사하지만 구별 가능 — rank 의존성 가능")
    else:
        log(f"  → 확실히 다름 — L-함수 유형별 상이한 보편상수")

flush()
log("\n[완료]")
flush()
