#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] 결과 #69b — 2차 Maass form (EVEN, R=13.78) 4성질 검증
=============================================================================
대상: 첫 번째 EVEN Maass cusp form on SL(2,Z)
  - 스펙트럼 파라미터 R = 13.7797513518907389...
  - 대칭: EVEN (symmetryclass=0)
  - Conductor N=1
  - Root number ε = +1

감마 인자 (EVEN Maass form):
  γ(s) = π^{-s} Γ((s+iR)/2) Γ((s-iR)/2)

Λ(s) = γ(s) · L(s)
함수방정식: Λ(s) = Λ(1-s)

Fourier 계수: Booker-Strömbergsson-Venkatesh 데이터 (1000자리, 455개)
출처: https://www2.math.uu.se/~astrombe/emaass/psl2z/coeff13

목적:
  1. κ_near ≈ 1114 재현 → degree-의존 상수 확정
  2. σ-유일성 PASS 재현 → Maass 보편 특성 확인
  3. even parity → 감마 인자 다름 → 코드 유연성 + 물리적 의미 검증

결과 파일: results/maass_4property_69b.txt
=============================================================================
"""

import sys, os, time, re
import numpy as np
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

DPS = 50
mpmath.mp.dps = DPS

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "maass_4property_69b.txt"
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
R_SPECTRAL = mpmath.mpf("13.7797513518907389442436732815177125971551325687934870692523882216")
EPSILON = 1  # root number (SL(2,Z) 레벨 1)
N_COND = 1   # conductor

# Gauss-Hermite 직교
N_HERMITE = 50
_herm_nodes, _herm_weights = np.polynomial.hermite.hermgauss(N_HERMITE)
HERM_NODES = [mpmath.mpf(float(x)) for x in _herm_nodes]
HERM_WEIGHTS = [mpmath.mpf(float(w)) for w in _herm_weights]

CONTOUR_C = mpmath.mpf(2)  # AFE contour shift

# κ 스윕 파라미터
T_MIN, T_MAX = 2.0, 60.0   # 넓은 범위 (≥8 영점 확보)
DT_SCAN = 0.4               # 영점 스캔 간격 (촘촘하게)
DELTA = 0.03                 # σ 오프셋
MONO_R = 0.4                 # 모노드로미 반지름
MONO_N = 48                  # 모노드로미 단계

# ━━━━━━━━━━━ Fourier 계수 로드 ━━━━━━━━━━━
def load_coefficients():
    """Strömbergsson 데이터 파일에서 Maass form 계수 로드"""
    coeff_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "maass_coeff13_raw.txt")

    if not os.path.exists(coeff_file):
        log(f"⚠️ 계수 파일 없음: {coeff_file}")
        log("  Booker-Strömbergsson-Venkatesh 데이터 다운로드 시도...")
        try:
            import requests
            r = requests.get(
                'http://www2.math.uu.se/~astrombe/emaass/psl2z/coeff13',
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
        cn_dict[n] = mpmath.mpf(v_str[:55])

    return R_parsed, sym_class, cn_dict


# ━━━━━━━━━━━ 감마 인자 (EVEN Maass form) ━━━━━━━━━━━

def gamma_factor_maass(s):
    """
    EVEN Maass form 감마 인자:
    γ(s) = π^{-s} · Γ((s+iR)/2) · Γ((s-iR)/2)

    주의: ODD는 π^{-(s+1)} · Γ((s+1+iR)/2) · Γ((s+1-iR)/2)
    EVEN은 shift 없음!
    """
    iR = mpmath.mpc(0, R_SPECTRAL)
    return (mpmath.power(mpmath.pi, -s)
            * mpmath.gamma((s + iR) / 2)
            * mpmath.gamma((s - iR) / 2))


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
log("결과 #69b — 2차 Maass form (EVEN) L-함수 4성질 검증")
log("=" * 70)
log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log(f"DPS={DPS}, N_hermite={N_HERMITE}, c={float(CONTOUR_C)}")
log(f"R = {mpmath.nstr(R_SPECTRAL, 30)}")
log(f"대칭: EVEN (symmetryclass=0)")
log(f"감마: π^{{-s}} Γ((s+iR)/2) Γ((s-iR)/2)")
log(f"ε = {EPSILON}, N = {N_COND}")
log()

t_total = time.time()

# ━━━━━━ Step 0: 계수 로드 ━━━━━━

log("[Step 0] Fourier 계수 로드 (Booker-Strömbergsson-Venkatesh, R≈13.78)")
result = load_coefficients()
if result is None:
    log("❌ 계수 로드 실패 — 중단")
    flush()
    sys.exit(1)

R_parsed, sym_class, cn_dict = result
N_max = max(cn_dict.keys())

log(f"  R (파일) = {R_parsed}...")
log(f"  symmetryclass = {sym_class} ({'EVEN' if sym_class == 0 else 'ODD' if sym_class == 1 else '???'})")
log(f"  계수 수: {len(cn_dict)} (n=1..{N_max})")
log(f"  a(1) = {float(cn_dict[1]):.6f}")
log(f"  a(2) = {float(cn_dict[2]):.15f}")
log(f"  a(3) = {float(cn_dict[3]):.15f}")
log(f"  a(5) = {float(cn_dict[5]):.15f}")
log(f"  a(7) = {float(cn_dict[7]):.15f}")
log(f"  a(11) = {float(cn_dict[11]):.15f}")

if sym_class != 0:
    log(f"  ❌ symmetryclass={sym_class}, 예상: 0 (EVEN) — 중단")
    flush()
    sys.exit(1)

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

# a(8) = a(2)^3 - 2*a(2) (Hecke p^3 재귀)
a8_check = cn_dict[2]**3 - 2*cn_dict[2]
a8_diff = abs(float(cn_dict[8]) - float(a8_check))
log(f"  Hecke 재귀: a(8)={float(cn_dict[8]):.10f}, a(2)³-2a(2)={float(a8_check):.10f}, diff={a8_diff:.2e}")

log()
flush()

# ━━━━━━ Step 1: 함수방정식 검증 ━━━━━━

log("[Step 1] 함수방정식 Λ(s) = εΛ(1-s) on critical line")
fe_pts = [
    mpmath.mpc(0.5, 5),
    mpmath.mpc(0.5, 10),
    mpmath.mpc(0.5, 15),
    mpmath.mpc(0.5, 25),
    mpmath.mpc(0.5, 40),
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
    log(f"  경고: EVEN 감마 인자 사용 중. 재확인 필요.")
log()
flush()

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
dt_scan_elapsed = time.time() - t_scan
log(f"  스캔 완료: {dt_scan_elapsed:.1f}s ({dt_scan_elapsed/60:.1f}분)")

# ━━━━━━ Step 3: 영점 이분법 정밀화 ━━━━━━

zero_approx = []
for i in range(len(scan_arr) - 1):
    if scan_arr[i] * scan_arr[i + 1] < 0:
        t_lo, t_hi = float(scan_ts[i]), float(scan_ts[i + 1])
        val_lo = scan_arr[i]
        try:
            for _ in range(25):  # 25회 이분법 (정밀도 향상)
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
else:
    for j, z in enumerate(zero_approx[:20]):
        log(f"    #{j+1}: γ = {z:.10f}")
    if n_zeros > 20:
        log(f"    ... ({n_zeros - 20}개 추가)")

log()
flush()

# ━━━━━━ Step 4: κ near/far 측정 + 모노드로미 ━━━━━━

near_med = 0
far_med = 0
ratio = 0
near_cv = 999
mono_pass = 0
total_mono = 0
near_fin = []
kappa_near = []
mono_results = []

if n_zeros >= 2:
    log(f"[Step 3] κ near/far 비율 + 모노드로미")
    n_tp = min(12, n_zeros)  # 최대 12 TP (≥8 목표)

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

    n_fp = min(8, len(fp_ts))
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

# ━━━━━━ Step 5: σ-유일성 테스트 (≥8 영점) ━━━━━━

log("[Step 4] σ-유일성 테스트 (7σ × 영점, 목표 ≥8개)")
sigma_pass_count = 0
sigma_test_count = 0

if n_zeros >= 2:
    sigma_vals = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
    n_sigma_test = min(max(8, n_zeros), n_zeros)  # 최소 8개, 최대 전체
    test_zeros = zero_approx[:n_sigma_test]

    for z_idx, z_t in enumerate(test_zeros):
        log(f"\n  영점 #{z_idx+1} γ={z_t:.6f}:")
        kappas_by_sigma = []
        for sigma in sigma_vals:
            k_val = curvature(z_t, cn_dict, N_max, delta=sigma - 0.5)
            kappas_by_sigma.append(k_val)
            log(f"    σ={sigma:.2f}: κ={k_val:.4e}")
            flush()

        # σ=0.5에서의 κ가 최대인지 확인
        idx_05 = sigma_vals.index(0.5)
        k_05 = kappas_by_sigma[idx_05]
        is_max = all(k_05 >= k for k in kappas_by_sigma if np.isfinite(k))
        sigma_test_count += 1
        if is_max:
            sigma_pass_count += 1
        log(f"    σ=0.5 최대? {'✅ YES' if is_max else '❌ NO'} (κ(0.5)={k_05:.4e})")

    log(f"\n  σ-유일성 결과: {sigma_pass_count}/{sigma_test_count}")
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
log(f"대상: Maass cusp form on SL(2,Z), R={mpmath.nstr(R_SPECTRAL, 20)}")
log(f"대칭: EVEN (symmetryclass=0)")
log(f"감마: π^{{-s}} Γ((s+iR)/2) Γ((s-iR)/2)")
log(f"ε = {EPSILON}, N = {N_COND}")
log(f"계수: {len(cn_dict)}개 (Booker-Strömbergsson-Venkatesh)")
log(f"")
log(f"1. 함수방정식: {'✅ PASS' if fe_ok else '❌ FAIL'} (max_rel={fe_max_rel:.2e})")
log(f"2. 영점 발견: {n_zeros}개 (t ∈ [{T_MIN}, {T_MAX}])")

if n_zeros >= 2 and near_fin:
    log(f"3. κ_near median = {near_med:.2f} (n={len(near_fin)}, CV={near_cv:.1f}%)")
    log(f"   κ_far median  = {far_med:.4f}")
    log(f"   κ ratio = {ratio:.1f}×")

    # GL(3) vs GL(2) Maass 비교
    diff_gl3 = abs(near_med - 1125.16) / 1125.16 * 100
    diff_69 = abs(near_med - 1114.38) / 1114.38 * 100
    log(f"   GL(3) sym² κ_near = 1125.16 → 차이 {diff_gl3:.2f}%")
    log(f"   #69 (ODD, R=9.53) κ_near = 1114.38 → 차이 {diff_69:.2f}%")

    if diff_69 < 1.0:
        log(f"   → κ_near ≈ 1114: GL(2) Maass degree-의존 상수 확정")
    elif diff_gl3 < 1.0:
        log(f"   → κ_near ≈ 1125: 보편상수 유지, R=9.53 특이성")
    else:
        log(f"   → 제3의 값: 더 복잡한 구조 (R-의존? parity-의존?)")

    log(f"4. 모노드로미: {mono_pass}/{total_mono} (mono/π≈2.0)")
else:
    log(f"3. κ_near: 영점 부족으로 측정 불가")
    log(f"4. 모노드로미: 영점 부족")

log(f"5. σ-유일성: {sigma_pass_count}/{sigma_test_count}")
log(f"")
log(f"총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")
log(f"")

# TP별 상세 테이블
if n_zeros >= 2 and near_fin:
    log("TP별 상세:")
    log(f"{'TP':>4} | {'γ':>12} | {'κ_near':>10} | {'mono/π':>8}")
    log("-" * 50)
    for j in range(len(kappa_near)):
        m_str = f"{mono_results[j]:.4f}" if mono_results[j] is not None else "FAIL"
        log(f"  #{j+1:>2} | {zero_approx[j]:>12.6f} | {kappa_near[j]:>10.2f} | {m_str:>8}")

log()

# GL(2) Maass 비교 테이블
if n_zeros >= 2 and near_fin:
    log("GL(2) Maass 비교:")
    log(f"  #69 (ODD,  R=9.53):  κ_near = 1114.38, CV=0.1%, σ-유일 3/3")
    log(f"  #69b (EVEN, R=13.78): κ_near = {near_med:.2f}, CV={near_cv:.1f}%, σ-유일 {sigma_pass_count}/{sigma_test_count}")
    log(f"  GL(3) sym² 평균:      κ_near = 1125.16, CV=0.15%")
    log()
    log("핵심 판별:")
    if diff_69 < 1.0 and diff_gl3 > 0.5:
        log("  → ✅ κ_near ≈ 1114: degree-의존 상수 확정 (GL(2) Maass ≠ GL(3) sym²)")
    elif diff_gl3 < 1.0 and diff_69 > 0.5:
        log("  → ✅ κ_near ≈ 1125: 보편상수 유지, #69 특이성 가능")
    elif diff_69 < 1.0 and diff_gl3 < 1.0:
        log("  → ⚠️ 두 값 중간: 통계적 판별 어려움, 추가 Maass form 필요")
    else:
        log("  → ⚠️ 제3의 값: 더 복잡한 구조 (R-의존? parity-의존?)")

flush()
log("\n[완료]")
flush()
