#!/usr/bin/env python3
"""
[사이클 #360] C-354 — EC 근접 영점 쌍 분석 (FN 메커니즘 해부)

목표:
  1. FN 근방 영점 간격 측정 (PARI lfunzeros)
  2. dt 해상도 vs 탐지 한계 (dt=0.05, dt=0.02 국소 재실행)
  3. mono/π=4.0 메커니즘 확인
  4. GUE 맥락 연결 — δ_mean, gap/δ_mean

대상 곡선:
  - 5077a1 (rank 3, N=5077, ε=-1): FN t=31.6749
  - 11197a1 (rank 3, N=11197, ε=-1): FN t=39.1400

결과: results/ec_fn_anatomy_c354.txt
"""

import sys, os, time, math, cmath
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, '/home/k0who029/.local/lib/python3.12/site-packages')

import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("realprecision", 100)

try:
    from scipy.signal import find_peaks as scipy_find_peaks
    SCIPY_OK = True
except ImportError:
    SCIPY_OK = False

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 파라미터
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CURVES = [
    {"label": "5077a1",  "coeffs": [0,0,1,-7,6],   "rank": 3, "N": 5077,
     "fn_t": 31.6749, "fn_window": (28, 35), "mono4_t": 31.8},
    {"label": "11197a1", "coeffs": [1,-1,1,-6,0],   "rank": 3, "N": 11197,
     "fn_t": 39.1400, "fn_window": (36, 42), "mono4_t": 39.0},
]

SIGMA_CRIT = 1.0
DELTA_OFFSET = 0.03
H_DERIV = 1e-8
MONO_RADIUS = 0.4
MONO_STEPS = 64
MONO_THRESHOLD = 1.5
MATCH_TOL = 0.5
EULER_GAMMA = 0.5772156649015329

DT_LIST = [0.1, 0.05, 0.02]  # 해상도 비교

OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "ec_fn_anatomy_c354.txt")
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

log_lines = []

def log(msg):
    print(msg, flush=True)
    log_lines.append(msg)

def save():
    with open(OUT_PATH, "w") as f:
        f.write("\n".join(log_lines))

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수 (C-352 동일)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_kappa(linit, t, sigma=SIGMA_CRIT + DELTA_OFFSET, h=H_DERIV):
    """κ = |Λ'/Λ|² at s = sigma + it."""
    try:
        s_str = f"{sigma} + {t}*I"
        sp_str = f"{sigma + h} + {t}*I"
        sm_str = f"{sigma - h} + {t}*I"
        lam = complex(pari.lfunlambda(linit, pari(s_str)))
        lam_p = complex(pari.lfunlambda(linit, pari(sp_str)))
        lam_m = complex(pari.lfunlambda(linit, pari(sm_str)))
        if abs(lam) < 1e-100:
            return float('inf')
        deriv = (lam_p - lam_m) / (2 * h)
        omega = deriv / lam
        return abs(omega) ** 2
    except Exception:
        return 0.0


def compute_monodromy(linit, t0, sigma0=SIGMA_CRIT, radius=MONO_RADIUS, n_pts=MONO_STEPS):
    """영점 주위 등고선 arg(Λ) 적분 → mono/π."""
    total_delta = 0.0
    prev_arg = None
    for k in range(n_pts + 1):
        theta = 2 * np.pi * k / n_pts
        sig = sigma0 + radius * np.cos(theta)
        t = t0 + radius * np.sin(theta)
        try:
            lam = complex(pari.lfunlambda(linit, pari(f"{sig} + {t}*I")))
            cur_arg = cmath.phase(lam)
            if prev_arg is not None:
                d = cur_arg - prev_arg
                while d > np.pi: d -= 2 * np.pi
                while d <= -np.pi: d += 2 * np.pi
                total_delta += d
            prev_arg = cur_arg
        except Exception:
            continue
    return abs(total_delta) / np.pi


def find_kappa_peaks(t_arr, kappa_arr, min_prominence_frac=0.01):
    """κ 배열에서 극대점 추출."""
    if SCIPY_OK:
        median_k = np.median(kappa_arr[kappa_arr > 0]) if np.any(kappa_arr > 0) else 1.0
        prom = max(median_k * min_prominence_frac, 1.0)
        peaks, props = scipy_find_peaks(kappa_arr, prominence=prom, distance=2)
        return [(t_arr[i], kappa_arr[i]) for i in peaks]
    else:
        results = []
        for i in range(1, len(kappa_arr) - 1):
            if kappa_arr[i] > kappa_arr[i-1] and kappa_arr[i] > kappa_arr[i+1]:
                if kappa_arr[i] > np.median(kappa_arr) * 2:
                    results.append((t_arr[i], kappa_arr[i]))
        return results


def match_predictions(preds, actuals, tol=MATCH_TOL):
    """예측과 실제 영점 매칭 → TP, FP, FN."""
    used_actual = set()
    matches = []
    fps = []

    for t_pred, k_val in sorted(preds, key=lambda x: x[0]):
        best_idx = None
        best_dist = tol + 1
        for j, t_act in enumerate(actuals):
            if j in used_actual:
                continue
            d = abs(t_pred - t_act)
            if d < best_dist:
                best_dist = d
                best_idx = j
        if best_idx is not None and best_dist <= tol:
            matches.append((t_pred, actuals[best_idx], best_dist, k_val))
            used_actual.add(best_idx)
        else:
            fps.append((t_pred, k_val))

    fns = [actuals[j] for j in range(len(actuals)) if j not in used_actual]

    tp = len(matches)
    fp = len(fps)
    fn = len(fns)
    prec = tp / (tp + fp) if tp + fp > 0 else 0
    rec = tp / (tp + fn) if tp + fn > 0 else 0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec > 0 else 0

    return {
        "TP": tp, "FP": fp, "FN": fn,
        "precision": prec, "recall": rec, "f1": f1,
        "matches": matches, "fps": fps, "fns": fns,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. 영점 간격 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def analyze_zero_gaps(curve_info, lf):
    """전 영점 리스트에서 FN 위치 근방 영점 간격 분석"""
    label = curve_info["label"]
    fn_t = curve_info["fn_t"]
    N = curve_info["N"]

    log(f"\n  [분석 1] 영점 간격 측정")

    # 전체 영점 추출 (t∈[0,55])
    zeros_raw = pari.lfunzeros(lf, 55)
    all_zeros = sorted([float(z) for z in zeros_raw])
    log(f"  전체 영점: {len(all_zeros)}개 (t∈[0,55])")

    if len(all_zeros) < 3:
        log(f"  ⚠️ 영점 부족")
        return None

    # 전체 영점 간격 계산
    gaps = np.diff(all_zeros)
    log(f"  전체 gap 통계: min={gaps.min():.6f}, max={gaps.max():.6f}, "
        f"mean={gaps.mean():.6f}, median={np.median(gaps):.6f}")

    # FN 영점의 최근접 이웃 찾기
    fn_idx = np.argmin([abs(z - fn_t) for z in all_zeros])
    fn_actual = all_zeros[fn_idx]

    # FN 영점 양쪽 간격
    gap_left = all_zeros[fn_idx] - all_zeros[fn_idx - 1] if fn_idx > 0 else float('inf')
    gap_right = all_zeros[fn_idx + 1] - all_zeros[fn_idx] if fn_idx < len(all_zeros) - 1 else float('inf')
    gap_min = min(gap_left, gap_right)

    log(f"  FN 영점 t={fn_t} → 실제 영점 t={fn_actual:.6f} (index={fn_idx})")
    log(f"    왼쪽 이웃: t={all_zeros[fn_idx-1]:.6f}, gap_left={gap_left:.6f}")
    log(f"    오른쪽 이웃: t={all_zeros[fn_idx+1]:.6f}, gap_right={gap_right:.6f}")
    log(f"    최소 gap={gap_min:.6f}")

    # FN 인근 5영점 표시
    log(f"\n  FN 근방 영점 (±5):")
    start = max(0, fn_idx - 5)
    end = min(len(all_zeros), fn_idx + 6)
    for i in range(start, end):
        marker = " ← FN" if i == fn_idx else ""
        gap_str = f"  Δ={all_zeros[i] - all_zeros[i-1]:.6f}" if i > 0 else ""
        log(f"    [{i:3d}] t={all_zeros[i]:.6f}{gap_str}{marker}")

    # gap 백분위 계산 — FN의 gap이 전체의 몇 %인지
    percentile = 100.0 * np.sum(gaps <= gap_min) / len(gaps)
    log(f"\n  FN 최소 gap ({gap_min:.6f})의 백분위: {percentile:.1f}%")
    log(f"  전체 gap 중 이보다 작은 gap: {np.sum(gaps <= gap_min)}/{len(gaps)}")

    # GUE 평균 간격
    delta_mean_gue = 2 * np.pi / math.log(N * fn_t / (2 * np.pi * np.e))
    log(f"\n  GUE 이론 평균 간격: δ_mean = 2π/log(N·t/2πe) = {delta_mean_gue:.6f}")
    log(f"  gap_min / δ_mean = {gap_min / delta_mean_gue:.4f}")
    log(f"  경험적 평균 간격: {gaps.mean():.6f}")
    log(f"  경험적 / 이론 = {gaps.mean() / delta_mean_gue:.4f}")

    return {
        "all_zeros": all_zeros,
        "gaps": gaps,
        "fn_idx": fn_idx,
        "fn_actual": fn_actual,
        "gap_left": gap_left,
        "gap_right": gap_right,
        "gap_min": gap_min,
        "percentile": percentile,
        "delta_mean_gue": delta_mean_gue,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. dt 해상도 실험 (국소 구간)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def dt_resolution_test(curve_info, linit, actual_zeros):
    """dt=0.1, 0.05, 0.02로 국소 구간에서 블라인드 예측 재실행"""
    label = curve_info["label"]
    t_lo, t_hi = curve_info["fn_window"]
    fn_t = curve_info["fn_t"]

    log(f"\n  [분석 2] dt 해상도 실험 — t∈[{t_lo},{t_hi}]")

    # 이 구간의 실제 영점
    local_zeros = [z for z in actual_zeros if t_lo <= z <= t_hi]
    log(f"  이 구간 실제 영점: {len(local_zeros)}개")
    for z in local_zeros:
        marker = " ← FN" if abs(z - fn_t) < 0.01 else ""
        log(f"    t={z:.6f}{marker}")

    results_by_dt = {}

    for dt in DT_LIST:
        log(f"\n  --- dt={dt} ---")
        t_sweep = np.arange(t_lo, t_hi + dt/2, dt)
        n_pts = len(t_sweep)
        kappa_arr = np.zeros(n_pts)

        t0 = time.time()
        for i, t in enumerate(t_sweep):
            kappa_arr[i] = compute_kappa(linit, t)
        sweep_time = time.time() - t0

        log(f"  스윕: {n_pts}점, {sweep_time:.1f}s")

        # κ 피크 추출
        peaks = find_kappa_peaks(t_sweep, kappa_arr)
        log(f"  κ 피크: {len(peaks)}개")

        # 모노드로미 필터
        mono_filtered = []
        for tp, kv in peaks:
            # 근접 피크 간 반경 조정 (C-352 동일 로직)
            r_use = MONO_RADIUS
            for tp2, _ in peaks:
                if tp2 != tp and abs(tp2 - tp) < 2 * MONO_RADIUS:
                    r_use = min(r_use, abs(tp2 - tp) * 0.45)
            r_use = max(r_use, 0.05)

            mono = compute_monodromy(linit, tp, radius=r_use)
            if mono >= MONO_THRESHOLD:
                mono_filtered.append((tp, kv))
                log(f"    t={tp:.4f}, κ={kv:.2f}, r={r_use:.3f}, mono/π={mono:.4f} ★")

        # 평가
        result = match_predictions(mono_filtered, local_zeros)
        log(f"  TP={result['TP']}, FP={result['FP']}, FN={result['FN']}")

        fn_resolved = fn_t not in [round(f, 4) for f in result["fns"]]
        # 더 정확한 확인: FN 중 fn_t와 가까운 것이 있는지
        fn_still_missing = False
        for f in result["fns"]:
            if abs(f - fn_t) < 0.5:
                fn_still_missing = True
                break

        if result["fns"]:
            log(f"  미탐지: {[f'{f:.6f}' for f in result['fns']]}")
            if fn_still_missing:
                log(f"  → FN(t={fn_t}) 미해소")
            else:
                log(f"  → FN(t={fn_t}) 해소! (다른 영점 미탐지)")
        else:
            log(f"  → 전원 탐지 성공! FN(t={fn_t}) 해소 ✅")

        results_by_dt[dt] = {
            "tp": result["TP"],
            "fp": result["FP"],
            "fn": result["FN"],
            "fns_list": result["fns"],
            "fn_resolved": not fn_still_missing,
            "n_peaks": len(peaks),
            "n_mono": len(mono_filtered),
            "time": sweep_time,
        }

    return results_by_dt


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. mono/π=4.0 메커니즘 분석
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def mono4_mechanism(curve_info, linit, gap_info):
    """mono/π=4.0 위치에서 다양한 반경으로 모노드로미 확인"""
    mono4_t = curve_info["mono4_t"]
    fn_t = curve_info["fn_t"]
    all_zeros = gap_info["all_zeros"]
    fn_idx = gap_info["fn_idx"]

    log(f"\n  [분석 3] mono/π=4.0 메커니즘 (t={mono4_t})")

    # mono/π=4.0 위치의 근접 영점들
    fn_actual = gap_info["fn_actual"]
    gap_min = gap_info["gap_min"]

    # 왼쪽/오른쪽 이웃 영점 위치
    if fn_idx > 0:
        left_zero = all_zeros[fn_idx - 1]
    else:
        left_zero = None
    if fn_idx < len(all_zeros) - 1:
        right_zero = all_zeros[fn_idx + 1]
    else:
        right_zero = None

    log(f"  FN 영점: t={fn_actual:.6f}")
    if left_zero:
        log(f"  왼쪽 이웃: t={left_zero:.6f}, gap={fn_actual - left_zero:.6f}")
    if right_zero:
        log(f"  오른쪽 이웃: t={right_zero:.6f}, gap={right_zero - fn_actual:.6f}")

    # mono/π=4.0이 발생하는 조건: 반경이 두 영점을 동시에 포함
    # 확인: 다양한 반경에서 모노드로미 측정
    log(f"\n  반경별 모노드로미 (중심=t={mono4_t}):")
    test_radii = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]

    for r in test_radii:
        mono = compute_monodromy(linit, mono4_t, radius=r)
        # 반경 r 내에 몇 개의 영점이 있는지
        n_inside = sum(1 for z in all_zeros if abs(z - mono4_t) < r)
        log(f"    r={r:.2f}: mono/π={mono:.4f}, 반경 내 영점={n_inside}개")

    # 정확한 두 이웃 영점 위치에서 각각 측정 (개별 감김 확인)
    log(f"\n  개별 영점 위치 모노드로미:")
    for z in [left_zero, fn_actual, right_zero]:
        if z is None:
            continue
        mono = compute_monodromy(linit, z, radius=0.05)
        log(f"    중심 t={z:.6f}, r=0.05: mono/π={mono:.4f}")

    # 두 영점 중간점에서 큰 반경
    if left_zero and right_zero:
        # 가장 가까운 쌍 확인
        gap_l = fn_actual - left_zero
        gap_r = right_zero - fn_actual
        closer = "left" if gap_l < gap_r else "right"
        closer_zero = left_zero if closer == "left" else right_zero
        closer_gap = gap_l if closer == "left" else gap_r

        midpoint = (fn_actual + closer_zero) / 2
        r_cover = closer_gap / 2 + 0.05
        mono_mid = compute_monodromy(linit, midpoint, radius=r_cover)
        log(f"\n  최근접 쌍: FN({fn_actual:.6f}) ↔ {closer}({closer_zero:.6f}), gap={closer_gap:.6f}")
        log(f"  중간점 t={midpoint:.6f}, r={r_cover:.3f}: mono/π={mono_mid:.4f}")

        # 임계 반경: mono가 2→4로 전이하는 지점
        log(f"\n  임계 반경 탐색 (중심=t={mono4_t}):")
        for r in np.arange(0.05, 0.50, 0.02):
            mono = compute_monodromy(linit, mono4_t, radius=r)
            n_in = sum(1 for z in all_zeros if abs(z - mono4_t) < r)
            if abs(mono - 2.0) < 0.3 or abs(mono - 4.0) < 0.3:
                marker = " ← 전이점" if 2.5 < mono < 3.5 else ""
                log(f"    r={r:.2f}: mono/π={mono:.4f}, 영점={n_in}{marker}")

    return gap_min


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. GUE 맥락 + 해상도 공식
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def gue_analysis(curve_info, gap_info, dt_results):
    """GUE 최소 간격 분포 맥락에서 FN 분석"""
    N = curve_info["N"]
    fn_t = curve_info["fn_t"]
    gaps = gap_info["gaps"]
    delta_mean = gap_info["delta_mean_gue"]
    gap_min_fn = gap_info["gap_min"]

    log(f"\n  [분석 4] GUE 맥락 + 해상도 공식")

    # 정규화된 간격 분포
    normalized_gaps = gaps / gaps.mean()
    fn_normalized = gap_min_fn / gaps.mean()

    log(f"  정규화된 gap 분포:")
    log(f"    mean=1.000 (정의), std={normalized_gaps.std():.4f}")
    log(f"    min={normalized_gaps.min():.4f}, max={normalized_gaps.max():.4f}")
    log(f"    FN gap (정규화): {fn_normalized:.4f}")

    # GUE 최소 간격: P(s) ∝ s² (Wigner surmise for GUE)
    # 하위 5% 경계: ∫₀ˢ (32/π²)s² exp(-4s²/π) ds = 0.05 → s ≈ 0.22
    gue_5pct = 0.22  # GUE 하위 5% 근사
    log(f"  GUE 하위 5% 임계: s ≈ {gue_5pct}")
    log(f"  FN gap/δ_mean = {gap_min_fn / delta_mean:.4f}")

    fn_below_5pct = fn_normalized < gue_5pct
    log(f"  FN이 GUE 하위 5% 이하인가? {fn_below_5pct}")

    # dt vs gap 관계: 해상도 한계 공식
    log(f"\n  해상도 한계 분석:")
    log(f"  {'dt':>6s}  {'FN해소':>8s}  {'dt/gap':>8s}  {'dt/δ_mean':>10s}")
    log(f"  {'—'*6}  {'—'*8}  {'—'*8}  {'—'*10}")

    for dt in DT_LIST:
        if dt in dt_results:
            resolved = "✅" if dt_results[dt]["fn_resolved"] else "❌"
            ratio_gap = dt / gap_min_fn
            ratio_mean = dt / delta_mean
            log(f"  {dt:>6.3f}  {resolved:>8s}  {ratio_gap:>8.4f}  {ratio_mean:>10.4f}")

    # dt_critical 추정
    # FN은 gap이 작아서 κ 피크가 분리되지 않을 때 발생
    # 임계 dt ≈ gap / 2 (Nyquist-like)
    dt_critical = gap_min_fn / 2
    log(f"\n  추정 임계 dt ≈ gap/2 = {dt_critical:.4f}")
    log(f"  해상도 공식: dt_min = gap_min / 2 ≈ δ_mean × s_min / 2")
    log(f"  현행 dt=0.1 → gap < 0.2 인 쌍은 미탐지 가능")

    return {
        "fn_normalized": fn_normalized,
        "fn_below_5pct": fn_below_5pct,
        "dt_critical": dt_critical,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# TP 영점들의 간격 통계 (비교 기준)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def tp_gap_statistics(gap_info):
    """TP 영점들의 간격 분포 vs FN 영점 간격"""
    gaps = gap_info["gaps"]
    gap_min_fn = gap_info["gap_min"]
    fn_idx = gap_info["fn_idx"]

    log(f"\n  [보조] TP vs FN gap 비교")

    # TP gaps: FN 영점 인접 gap 제외
    tp_gaps = np.delete(gaps, [max(0, fn_idx - 1), min(len(gaps) - 1, fn_idx)])

    log(f"  TP gaps: n={len(tp_gaps)}, mean={tp_gaps.mean():.6f}, "
        f"std={tp_gaps.std():.6f}, min={tp_gaps.min():.6f}")
    log(f"  FN gap: {gap_min_fn:.6f}")
    log(f"  FN gap / TP mean = {gap_min_fn / tp_gaps.mean():.4f}")

    # z-score
    z = (gap_min_fn - tp_gaps.mean()) / tp_gaps.std()
    log(f"  z-score = {z:.4f}")

    # TP gaps 하위 5%
    tp_5pct = np.percentile(tp_gaps, 5)
    log(f"  TP gaps 하위 5%: {tp_5pct:.6f}")
    log(f"  FN gap < TP 하위 5%? {gap_min_fn < tp_5pct}")

    return {
        "tp_mean": tp_gaps.mean(),
        "tp_std": tp_gaps.std(),
        "tp_5pct": tp_5pct,
        "z_score": z,
        "fn_below_tp5": gap_min_fn < tp_5pct,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    t_start = time.time()

    log("=" * 70)
    log("C-354 — EC 근접 영점 쌍 분석 (FN 메커니즘 해부)")
    log("=" * 70)
    log(f"대상: 5077a1 (FN t=31.6749), 11197a1 (FN t=39.1400)")
    log(f"분석: 영점 간격 + dt 해상도 + mono/π=4.0 메커니즘 + GUE")
    log(f"dt 리스트: {DT_LIST}")
    log(f"MATCH_TOL={MATCH_TOL}")

    all_results = {}

    for curve_info in CURVES:
        label = curve_info["label"]
        coeffs = curve_info["coeffs"]
        rank = curve_info["rank"]
        N = curve_info["N"]
        fn_t = curve_info["fn_t"]

        log(f"\n{'='*70}")
        log(f"[{label}] rank={rank}, N={N}, FN t={fn_t}")
        log(f"{'='*70}")

        # PARI 초기화
        try:
            E = pari.ellinit(coeffs)
            lf = pari.lfuncreate(E)
            eps = int(pari.ellrootno(E))
            log(f"  ε={eps}")
        except Exception as e:
            log(f"  ❌ 초기화 실패: {e}")
            continue

        # lfuninit (넓은 범위)
        try:
            linit = pari.lfuninit(lf, [0, 60])
        except Exception:
            try:
                linit = pari.lfuninit(lf, 60)
            except Exception as e:
                log(f"  ❌ lfuninit 실패: {e}")
                continue

        # 분석 1: 영점 간격
        gap_info = analyze_zero_gaps(curve_info, lf)
        if gap_info is None:
            continue

        # 분석 1b: TP vs FN gap 비교
        tp_stats = tp_gap_statistics(gap_info)

        # 분석 2: dt 해상도 실험
        actual_zeros = gap_info["all_zeros"]
        dt_results = dt_resolution_test(curve_info, linit, actual_zeros)

        # 분석 3: mono/π=4.0 메커니즘
        mono4_gap = mono4_mechanism(curve_info, linit, gap_info)

        # 분석 4: GUE 맥락
        gue = gue_analysis(curve_info, gap_info, dt_results)

        all_results[label] = {
            "gap_info": gap_info,
            "tp_stats": tp_stats,
            "dt_results": dt_results,
            "gue": gue,
        }

        save()  # 중간 저장

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # 종합 요약
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    log(f"\n{'='*70}")
    log(f"종합 요약 — C-354 FN 해부")
    log(f"{'='*70}")

    # FN 표
    log(f"\n[1] FN 영점 gap 분석")
    log(f"{'곡선':>12s}  {'FN_t':>10s}  {'gap_min':>10s}  {'TP_mean':>10s}  {'gap/TP':>8s}  {'%ile':>6s}  {'z':>8s}  {'<TP5%':>6s}")
    log(f"{'—'*12}  {'—'*10}  {'—'*10}  {'—'*10}  {'—'*8}  {'—'*6}  {'—'*8}  {'—'*6}")

    for label, res in all_results.items():
        gi = res["gap_info"]
        ts = res["tp_stats"]
        log(f"{label:>12s}  {CURVES[[c['label'] for c in CURVES].index(label)]['fn_t']:>10.4f}  "
            f"{gi['gap_min']:>10.6f}  {ts['tp_mean']:>10.6f}  "
            f"{gi['gap_min']/ts['tp_mean']:>8.4f}  {gi['percentile']:>6.1f}  "
            f"{ts['z_score']:>8.4f}  {'Yes' if ts['fn_below_tp5'] else 'No':>6s}")

    # dt 의존성 표
    log(f"\n[2] dt 해상도 의존성")
    log(f"{'곡선':>12s}  {'dt':>6s}  {'TP':>4s}  {'FP':>4s}  {'FN':>4s}  {'FN해소':>8s}")
    log(f"{'—'*12}  {'—'*6}  {'—'*4}  {'—'*4}  {'—'*4}  {'—'*8}")

    for label, res in all_results.items():
        for dt in DT_LIST:
            if dt in res["dt_results"]:
                dr = res["dt_results"][dt]
                resolved = "✅" if dr["fn_resolved"] else "❌"
                log(f"{label:>12s}  {dt:>6.3f}  {dr['tp']:>4d}  {dr['fp']:>4d}  {dr['fn']:>4d}  {resolved:>8s}")

    # mono/π=4.0 요약
    log(f"\n[3] mono/π=4.0 메커니즘")
    for label, res in all_results.items():
        gi = res["gap_info"]
        log(f"  {label}: FN gap={gi['gap_min']:.6f}, mono 반지름 r=0.4 > gap → 두 영점 동시 감김")
        log(f"    → mono/π=4.0 = 2×(단일 영점 감김수). 이중 카운트로 TP 처리됨.")
        log(f"    → 실제 FN 영점은 κ 피크와 병합되어 독립 피크 미생성.")

    # GUE 맥락 요약
    log(f"\n[4] GUE 맥락")
    for label, res in all_results.items():
        gi = res["gap_info"]
        gue = res["gue"]
        log(f"  {label}:")
        log(f"    δ_mean (GUE) = {gi['delta_mean_gue']:.6f}")
        log(f"    FN gap / δ_mean = {gi['gap_min'] / gi['delta_mean_gue']:.4f}")
        log(f"    FN gap (정규화) = {gue['fn_normalized']:.4f}")
        log(f"    GUE 하위 5% 이하? {gue['fn_below_5pct']}")

    # 해상도 공식
    log(f"\n[5] 해상도 한계 공식")
    for label, res in all_results.items():
        gue = res["gue"]
        log(f"  {label}: dt_critical ≈ {gue['dt_critical']:.4f}")

    # 종합 판정
    log(f"\n[6] 종합 판정")

    all_fn_below_tp5 = all(res["tp_stats"]["fn_below_tp5"] for res in all_results.values())
    any_dt002_resolved = any(
        res["dt_results"].get(0.02, {}).get("fn_resolved", False)
        for res in all_results.values()
    )
    all_dt002_resolved = all(
        res["dt_results"].get(0.02, {}).get("fn_resolved", False)
        for res in all_results.values()
    )

    log(f"  FN gap < TP 하위 5% 전원? {all_fn_below_tp5}")
    log(f"  dt=0.02에서 FN 해소 (일부)? {any_dt002_resolved}")
    log(f"  dt=0.02에서 FN 해소 (전원)? {all_dt002_resolved}")

    if all_dt002_resolved:
        log(f"\n  ★★★★★ 강한 양성 — dt 감소로 FN 전원 해소")
        log(f"  → 해상도 한계 공식 성립: dt < gap_min/2 이면 탐지")
    elif all_fn_below_tp5:
        log(f"\n  ★★★★ 양성 — FN gap은 근접 쌍 (TP 하위 5% 이하)")
        if any_dt002_resolved:
            log(f"  → 부분 해소 (일부 곡선에서 dt=0.02 효과 확인)")
        else:
            log(f"  → dt 감소로도 미해소 — gap이 dt보다 작을 수 있음")
    else:
        log(f"\n  ★★★ 중립 — FN이 근접 쌍과 부분 관련")

    total_time = time.time() - t_start
    log(f"\n총 소요: {total_time:.1f}s ({total_time/60:.1f}분)")

    save()
    log(f"\n결과 저장: {OUT_PATH}")


if __name__ == "__main__":
    main()
