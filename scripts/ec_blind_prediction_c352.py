#!/usr/bin/env python3
"""
[사이클 #358] C-352 — EC 블라인드 영점 예측 확장 (8곡선) + κ_bg 이론 비교

목표:
  - C-351의 3곡선 → 8곡선 확장 (rank당 2곡선)
  - 기존: 11a1(r0), 389a1(r2), 5077a1(r3)
  - 추가: 43a1(r0), 37a1(r1), 79a1(r1), 571a1(r2), 11197a1(r3)
  - κ_bg 이론값 vs 경험적 κ_median 비교

프로토콜: C-351 동일 (σ=1.03, t∈[2,50], dt=0.1, MONO r=0.4, steps=64, threshold=1.5)

성공 기준:
  - 8곡선 전원 F1≥0.90 → ★★★★★ 양성 (EC 블라인드 확립)
  - κ_bg ratio CV<30% → 이론 양성

결과: results/ec_blind_prediction_c352.txt
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

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 실험 파라미터
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 기존 3곡선 (C-351) + 추가 5곡선 = 총 8곡선
# rank당 2곡선: rank 0×2, 1×2, 2×2, 3×2
CURVES = [
    # rank 0
    {"label": "11a1",    "coeffs": [0,-1,1,-10,-20],    "rank": 0, "N": 11},
    {"label": "43a1",    "coeffs": [0,1,1,0,0],         "rank": 0, "N": 43},
    # rank 1
    {"label": "37a1",    "coeffs": [0,0,1,-1,0],        "rank": 1, "N": 37},
    {"label": "79a1",    "coeffs": [1,1,1,-2,0],        "rank": 1, "N": 79},
    # rank 2
    {"label": "389a1",   "coeffs": [0,1,1,-2,0],        "rank": 2, "N": 389},
    {"label": "571a1",   "coeffs": [0,1,1,-4,2],        "rank": 2, "N": 571},
    # rank 3
    {"label": "5077a1",  "coeffs": [0,0,1,-7,6],        "rank": 3, "N": 5077},
    {"label": "11197a1", "coeffs": [1,-1,1,-6,0],       "rank": 3, "N": 11197},
]

T_MIN = 2.0
T_MAX_SWEEP = 50.0
DT = 0.1
SIGMA_CRIT = 1.0
DELTA_OFFSET = 0.03
H_DERIV = 1e-8
MONO_RADIUS = 0.4
MONO_STEPS = 64
MONO_THRESHOLD = 1.5
MATCH_TOL = 0.5

# Euler-Mascheroni constant
EULER_GAMMA = 0.5772156649015329

OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "ec_blind_prediction_c352.txt")
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

log_lines = []

def log(msg):
    print(msg, flush=True)
    log_lines.append(msg)

def save():
    with open(OUT_PATH, "w") as f:
        f.write("\n".join(log_lines))


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수 (C-351 동일)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

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
        median_k = np.median(kappa_arr[kappa_arr > 0])
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


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# κ_bg 이론 비교
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def kappa_bg_theory(N, degree=2):
    """
    이론적 κ_bg 예측 (degree 2 EC):
    κ_bg = (log(N_eff / 4π²) + γ)² + π²/4
    여기서 N_eff = N (conductor), γ = Euler-Mascheroni
    """
    log_term = math.log(N / (4 * math.pi**2)) + EULER_GAMMA
    return log_term**2 + (math.pi**2) / 4


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 곡선별 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def run_curve(curve_info):
    label = curve_info["label"]
    coeffs = curve_info["coeffs"]
    rank = curve_info["rank"]
    N = curve_info["N"]

    log(f"\n{'='*70}")
    log(f"[{label}] rank={rank}, N={N}")
    log(f"{'='*70}")

    t0 = time.time()

    # PARI 초기화
    try:
        E = pari.ellinit(coeffs)
        lf = pari.lfuncreate(E)
        eps = int(pari.ellrootno(E))
    except Exception as e:
        log(f"  ❌ 초기화 실패: {e}")
        return None

    log(f"  ε={eps}")

    # Step 0: 실제 영점 추출 (검증용으로만!)
    log(f"\n  [Step 0] 실제 영점 추출 (검증용, t∈[{T_MIN},{T_MAX_SWEEP}])")
    try:
        zeros_raw = pari.lfunzeros(lf, T_MAX_SWEEP + 5)
        all_zeros = sorted([float(z) for z in zeros_raw])
        actual_zeros = [z for z in all_zeros if T_MIN <= z <= T_MAX_SWEEP]
        log(f"  전체 영점: {len(all_zeros)}개, 범위 내: {len(actual_zeros)}개")
    except Exception as e:
        log(f"  ❌ lfunzeros 실패: {e}")
        return None

    if len(actual_zeros) < 3:
        log(f"  ⚠️ 영점 {len(actual_zeros)}개 — 부족, 스킵")
        return None

    # lfuninit
    try:
        linit = pari.lfuninit(lf, [0, T_MAX_SWEEP + 10])
    except Exception:
        try:
            linit = pari.lfuninit(lf, T_MAX_SWEEP + 10)
        except Exception as e:
            log(f"  ❌ lfuninit 실패: {e}")
            return None

    # Step 1: κ 밀집 스윕 (블라인드)
    t_sweep = np.arange(T_MIN, T_MAX_SWEEP + DT/2, DT)
    n_pts = len(t_sweep)
    kappa_arr = np.zeros(n_pts)

    log(f"\n  [Step 1] κ 스윕: σ={SIGMA_CRIT + DELTA_OFFSET:.3f}, "
        f"t∈[{T_MIN},{T_MAX_SWEEP}], dt={DT}")
    log(f"  총 {n_pts}점")

    t1 = time.time()
    for i, t in enumerate(t_sweep):
        kappa_arr[i] = compute_kappa(linit, t)
        if (i + 1) % 100 == 0:
            elapsed = time.time() - t1
            rate = (i + 1) / elapsed if elapsed > 0 else float('inf')
            remaining = (n_pts - i - 1) / rate if rate > 0 else 0
            log(f"  ... {i+1}/{n_pts} ({100*(i+1)/n_pts:.0f}%), "
                f"경과={elapsed:.0f}s, 잔여≈{remaining:.0f}s")

    sweep_time = time.time() - t1
    valid = kappa_arr[kappa_arr > 0]
    if len(valid) == 0:
        log(f"  ⚠️ 유효 κ 없음 — 스킵")
        return None
    log(f"  스윕 완료: {sweep_time:.1f}s")
    log(f"  κ 통계: min={valid.min():.4f}, max={valid.max():.2f}, "
        f"median={np.median(valid):.4f}, mean={valid.mean():.4f}")

    kappa_median = float(np.median(valid))

    # Step 2: κ 피크 추출
    log(f"\n  [Step 2] κ 피크 추출")
    peaks = find_kappa_peaks(t_sweep, kappa_arr)
    log(f"  → {len(peaks)}개 피크 발견")

    for i, (tp, kv) in enumerate(peaks[:10]):
        log(f"    후보 {i+1:2d}: t={tp:.4f}, κ={kv:.4f}")
    if len(peaks) > 10:
        log(f"    ... (이하 {len(peaks)-10}개 생략)")

    # Step 3: 모노드로미 필터
    log(f"\n  [Step 3] 모노드로미 필터 (r={MONO_RADIUS}, n={MONO_STEPS})")

    mono_filtered = []
    for i, (tp, kv) in enumerate(peaks):
        r_use = MONO_RADIUS
        for tp2, _ in peaks:
            if tp2 != tp and abs(tp2 - tp) < 2 * MONO_RADIUS:
                r_use = min(r_use, abs(tp2 - tp) * 0.45)
        r_use = max(r_use, 0.05)

        mono = compute_monodromy(linit, tp, radius=r_use)
        marker = "★ 영점" if mono > MONO_THRESHOLD else "  배제"
        if i < 5 or mono > MONO_THRESHOLD:
            log(f"    [{i+1:2d}] t={tp:.4f}, κ={kv:.4f}, r={r_use:.3f}, "
                f"mono/π={mono:.4f} ← {marker}")

        if mono > MONO_THRESHOLD:
            mono_filtered.append((tp, kv, mono))

    log(f"\n  κ+mono 필터 결과: {len(mono_filtered)}개 후보 → 영점 예측")

    # Step 4: 평가
    log(f"\n  [Step 4] 평가 — 실제 영점과 비교")
    log(f"  실제 영점 {len(actual_zeros)}개 (t∈[{T_MIN},{T_MAX_SWEEP}])")

    # κ+mono 평가
    mono_preds = [(tp, kv) for tp, kv, _ in mono_filtered]
    r_kmono = match_predictions(mono_preds, actual_zeros)
    log(f"  예측 {len(mono_preds)}개 vs 실제 {len(actual_zeros)}개 (tol={MATCH_TOL})")
    log(f"  TP={r_kmono['TP']}, FP={r_kmono['FP']}, FN={r_kmono['FN']}")
    log(f"  Precision={r_kmono['precision']:.4f}, "
        f"Recall={r_kmono['recall']:.4f}, F1={r_kmono['f1']:.4f}")

    if r_kmono['fns']:
        log(f"  미탐지: {[f'{z:.4f}' for z in r_kmono['fns']]}")
    if r_kmono['fps']:
        log(f"  오탐지: {[f't={tp:.4f}' for tp, _ in r_kmono['fps']]}")

    elapsed_total = time.time() - t0
    log(f"\n  총 소요: {elapsed_total:.0f}s")

    return {
        "label": label, "rank": rank, "N": N, "eps": eps,
        "n_actual": len(actual_zeros),
        "kmono": r_kmono,
        "kappa_median": kappa_median,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

if __name__ == "__main__":
    t_start = time.time()

    log("=" * 70)
    log("[Project RDL] C-352 — EC 블라인드 영점 예측 확장 (8곡선) + κ_bg 이론 비교")
    log("=" * 70)
    log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"곡선: {[c['label'] for c in CURVES]}")
    log(f"t∈[{T_MIN},{T_MAX_SWEEP}], dt={DT}, δ={DELTA_OFFSET}")
    log(f"MONO: r={MONO_RADIUS}, steps={MONO_STEPS}, threshold={MONO_THRESHOLD}")
    log(f"MATCH_TOL={MATCH_TOL}")
    log(f"scipy: {SCIPY_OK}")

    all_results = []
    for curve in CURVES:
        result = run_curve(curve)
        if result:
            all_results.append(result)
        save()  # 중간 저장

    # ━━━━━ 최종 요약 ━━━━━
    log(f"\n\n{'='*70}")
    log(f"최종 요약 — C-352 블라인드 예측 (8곡선)")
    log(f"{'='*70}")

    log(f"\n  [1] 블라인드 예측 성능")
    log(f"  {'곡선':>10s}  {'rank':>4s}  {'N':>6s}  {'ε':>2s}  {'영점':>4s}  "
        f"{'P':>6s}  {'R':>6s}  {'F1':>7s}")
    log(f"  {'-'*55}")

    for r in all_results:
        km = r["kmono"]
        log(f"  {r['label']:>10s}  {r['rank']:>4d}  {r['N']:>6d}  {r['eps']:>2d}  {r['n_actual']:>4d}  "
            f"{km['precision']:>6.3f}  {km['recall']:>6.3f}  {km['f1']:>7.3f}")

    # 기존 결과 (#52, #54) + C-351 + C-352 통합
    log(f"\n  [2] 전체 통합표 (기존 + C-352)")
    log(f"  {'곡선':>10s}  {'rank':>4s}  {'N':>6s}  {'method':>10s}  "
        f"{'P':>6s}  {'R':>6s}  {'F1':>7s}")
    log(f"  {'-'*60}")

    # #52: 11a1 mpmath
    log(f"  {'11a1':>10s}  {'0':>4s}  {'11':>6s}  {'mpmath#52':>10s}  "
        f"{'0.941':>6s}  {'0.941':>6s}  {'0.941':>7s}")
    # #54: 37a1 mpmath
    log(f"  {'37a1':>10s}  {'1':>4s}  {'37':>6s}  {'mpmath#54':>10s}  "
        f"{'1.000':>6s}  {'0.957':>6s}  {'0.978':>7s}")

    for r in all_results:
        km = r["kmono"]
        log(f"  {r['label']:>10s}  {r['rank']:>4d}  {r['N']:>6d}  {'PARI-C352':>10s}  "
            f"{km['precision']:>6.3f}  {km['recall']:>6.3f}  {km['f1']:>7.3f}")

    # ━━━━━ κ_bg 이론 비교 ━━━━━
    log(f"\n  [3] κ_bg 이론 vs 경험 비교")
    log(f"  이론: κ_bg = (log(N/4π²) + γ)² + π²/4  (degree 2 EC)")
    log(f"  경험: κ_median from σ={SIGMA_CRIT + DELTA_OFFSET} sweep")
    log(f"")
    log(f"  {'곡선':>10s}  {'N':>6s}  {'rank':>4s}  {'κ_theory':>10s}  {'κ_empirical':>12s}  {'ratio':>8s}")
    log(f"  {'-'*60}")

    ratios = []
    for r in all_results:
        k_th = kappa_bg_theory(r["N"])
        k_emp = r["kappa_median"]
        ratio = k_emp / k_th if k_th > 0 else float('nan')
        ratios.append(ratio)
        log(f"  {r['label']:>10s}  {r['N']:>6d}  {r['rank']:>4d}  {k_th:>10.4f}  {k_emp:>12.4f}  {ratio:>8.4f}")

    if ratios:
        ratios_arr = np.array(ratios)
        ratio_mean = np.mean(ratios_arr)
        ratio_std = np.std(ratios_arr, ddof=1) if len(ratios_arr) > 1 else 0
        ratio_cv = (ratio_std / ratio_mean * 100) if ratio_mean > 0 else float('nan')
        log(f"")
        log(f"  ratio 통계: mean={ratio_mean:.4f}, std={ratio_std:.4f}, CV={ratio_cv:.1f}%")

        if ratio_cv < 30:
            log(f"  → CV<30%: ★ 이론 양성 — κ_bg ∝ (log(N/4π²)+γ)²+π²/4 일관")
        else:
            log(f"  → CV≥30%: 이론 중립 — 곡선 간 ratio 산포 큼, 보정 필요")

    # ━━━━━ 종합 판정 ━━━━━
    log(f"\n  [4] 종합 판정")

    all_f1 = [r["kmono"]["f1"] for r in all_results]
    all_rec = [r["kmono"]["recall"] for r in all_results]
    n_total = len(all_results)
    n_f1_pass = sum(1 for f in all_f1 if f >= 0.90)
    n_f1_fail = n_total - n_f1_pass

    # 블라인드 예측 판정
    if n_f1_pass == n_total and all(f >= 0.90 for f in all_f1):
        log(f"  블라인드: ★★★★★ 양성 — {n_total}곡선 전원 F1≥0.90")
    elif n_f1_fail <= 2:
        log(f"  블라인드: ★★★ 조건부 양성 — {n_f1_fail}곡선 F1<0.90")
    else:
        log(f"  블라인드: ❌ 음성 — {n_f1_fail}곡선 F1<0.90")

    # 통합 통계 (기존 #52/#54 포함)
    all_tp = sum(r["kmono"]["TP"] for r in all_results)
    all_fn = sum(r["kmono"]["FN"] for r in all_results)
    all_fp = sum(r["kmono"]["FP"] for r in all_results)
    # 기존 #52: 16 TP (mpmath, 0 FP, 1 FN), #54: 22 TP (0 FP, 1 FN)
    legacy_tp = 16 + 22
    legacy_fn = 1 + 1
    total_tp = all_tp + legacy_tp
    total_fn = all_fn + legacy_fn
    total_fp = all_fp + 0  # 기존 FP=0
    total_det = total_tp + total_fn
    total_recall = total_tp / total_det if total_det > 0 else 0

    log(f"\n  통합 (기존+C-352): {total_tp}/{total_det} = {total_recall*100:.1f}% 검출")
    log(f"  TP={total_tp}, FP={total_fp}, FN={total_fn}")

    total_time = time.time() - t_start
    log(f"\n총 소요: {total_time:.0f}s ({total_time/60:.1f}분)")
    log(f"종료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    save()
    log(f"\n결과 저장: {OUT_PATH}")
