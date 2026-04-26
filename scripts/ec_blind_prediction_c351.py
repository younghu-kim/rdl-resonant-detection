#!/usr/bin/env python3
"""
[사이클 #357] C-351 — EC 블라인드 영점 예측 확장 (rank 0, 2, 3)

목표:
  - #52(11a1, rank 0)와 #54(37a1, rank 1)의 확장
  - rank 2 (389a1)와 rank 3 (5077a1)에서 블라인드 예측 검증
  - 11a1은 PARI 기반 독립 교차검증 (#52는 mpmath AFE)

프로토콜:
  1. PARI로 L-함수 초기화 + 실제 영점 추출 (검증용으로만)
  2. σ = 1.03에서 t ∈ [2, 50], dt=0.1 → ��(t) 밀집 스윕
  3. κ = |Λ'/Λ|² (lfunlambda 수치미분)
  4. scipy.find_peaks로 κ 극대점 추출 → 후보 목록
  5. 각 후보에 모노드로미 측정 (등고선 적분)
  6. 비교: κ-only vs κ+mono, P/R/F1

성공 기준:
  - Recall ≥ 0.85 (각 곡선)
  - F1 ≥ 0.80 (각 곡선)
  - rank 2, 3에서 rank 0, 1과 동등 수준

결과: results/ec_blind_prediction_c351.txt
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

CURVES = [
    {"label": "11a1",   "coeffs": [0,-1,1,-10,-20],   "rank": 0, "N": 11},
    {"label": "389a1",  "coeffs": [0,1,1,-2,0],       "rank": 2, "N": 389},
    {"label": "5077a1", "coeffs": [0,0,1,-7,6],       "rank": 3, "N": 5077},
]

T_MIN = 2.0
T_MAX_SWEEP = 50.0
DT = 0.1
SIGMA_CRIT = 1.0
DELTA_OFFSET = 0.03
H_DERIV = 1e-8           # κ 수치미분 간격
MONO_RADIUS = 0.4         # 모노드로미 반지름 (근접 영점 시 자동 축소)
MONO_STEPS = 64           # 등고선 분할 수
MONO_THRESHOLD = 1.5      # mono/π > threshold → 영점 판정
MATCH_TOL = 0.5           # 예측-실제 매칭 허용 오차

OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "results", "ec_blind_prediction_c351.txt")
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

log_lines = []

def log(msg):
    print(msg, flush=True)
    log_lines.append(msg)

def save():
    with open(OUT_PATH, "w") as f:
        f.write("\n".join(log_lines))


# ━━━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━���━━━━━━

def compute_kappa(linit, t, sigma=SIGMA_CRIT + DELTA_OFFSET, h=H_DERIV):
    """κ = |Λ'/Λ|² at s = sigma + it. 수치미분 사용."""
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
        # 수동 극대점 추출
        results = []
        for i in range(1, len(kappa_arr) - 1):
            if kappa_arr[i] > kappa_arr[i-1] and kappa_arr[i] > kappa_arr[i+1]:
                if kappa_arr[i] > np.median(kappa_arr) * 2:
                    results.append((t_arr[i], kappa_arr[i]))
        return results


def match_predictions(preds, actuals, tol=MATCH_TOL):
    """예측과 실제 영점 매칭 → TP, FP, FN, 매칭 상세."""
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


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━��━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인 루프
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━���━━━━━━━━━━━━━━━━━━���━━━━━━━

def run_curve(curve_info):
    label = curve_info["label"]
    coeffs = curve_info["coeffs"]
    rank = curve_info["rank"]
    N = curve_info["N"]

    log(f"\n{'='*70}")
    log(f"[{label}] rank={rank}, N={N}")
    log(f"{'='*70}")

    t0 = time.time()

    # PARI 초기���
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

    # Step 1: κ 밀집 스윕 (블라인드 — 영점 위치 미사용)
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
            rate = (i + 1) / elapsed
            remaining = (n_pts - i - 1) / rate
            log(f"  ... {i+1}/{n_pts} ({100*(i+1)/n_pts:.0f}%), "
                f"경과={elapsed:.0f}s, 잔여≈{remaining:.0f}s, "
                f"t={t:.1f}, κ={kappa_arr[i]:.2f}")

    sweep_time = time.time() - t1
    valid = kappa_arr[kappa_arr > 0]
    log(f"  스윕 ���료: {sweep_time:.1f}s")
    log(f"  κ 통계: min={valid.min():.4f}, max={valid.max():.2f}, "
        f"median={np.median(valid):.4f}, mean={valid.mean():.4f}")

    # Step 2: κ 피크 추출
    log(f"\n  [Step 2] κ 피크 추출")
    peaks = find_kappa_peaks(t_sweep, kappa_arr)
    log(f"  → {len(peaks)}개 피크 발��")

    for i, (tp, kv) in enumerate(peaks[:30]):
        log(f"    후보 {i+1:2d}: t={tp:.4f}, κ={kv:.4f}")
    if len(peaks) > 30:
        log(f"    ... (이하 {len(peaks)-30}개 생략)")

    # Step 3: 모���드로미 필터
    log(f"\n  [Step 3] 모노드로미 필터 (r={MONO_RADIUS}, n={MONO_STEPS})")
    log(f"  기준: mono/π > {MONO_THRESHOLD}")

    mono_filtered = []
    for i, (tp, kv) in enumerate(peaks):
        # 적응형 반지름: 근접 피크가 있으면 축소
        r_use = MONO_RADIUS
        for tp2, _ in peaks:
            if tp2 != tp and abs(tp2 - tp) < 2 * MONO_RADIUS:
                r_use = min(r_use, abs(tp2 - tp) * 0.45)
        r_use = max(r_use, 0.05)

        mono = compute_monodromy(linit, tp, radius=r_use)
        marker = "★ 영점" if mono > MONO_THRESHOLD else "  배제"
        log(f"    [{i+1:2d}] t={tp:.4f}, κ={kv:.4f}, r={r_use:.3f}, "
            f"mono/π={mono:.4f} ← {marker}")

        if mono > MONO_THRESHOLD:
            mono_filtered.append((tp, kv, mono))

    log(f"\n  κ+mono 필터 결과: {len(mono_filtered)}개 후보 → 영점 예측")

    # Step 4: 평가
    log(f"\n  [Step 4] 평가 — 실제 영점과 비교")
    log(f"  실제 영점 {len(actual_zeros)}개 (t∈[{T_MIN},{T_MAX_SWEEP}]):")
    for j, z in enumerate(actual_zeros[:10]):
        log(f"    γ_{j+1:2d} = {z:.8f}")
    if len(actual_zeros) > 10:
        log(f"    ... (이하 {len(actual_zeros)-10}개)")

    # κ-only 평가
    log(f"\n  [평가: κ-only (피크 전체)]")
    r_konly = match_predictions(peaks, actual_zeros)
    log(f"  예��� {len(peaks)}개 vs 실제 {len(actual_zeros)}개 (tol={MATCH_TOL})")
    log(f"  TP={r_konly['TP']}, FP={r_konly['FP']}, FN={r_konly['FN']}")
    log(f"  Precision={r_konly['precision']:.4f}, "
        f"Recall={r_konly['recall']:.4f}, F1={r_konly['f1']:.4f}")
    if r_konly['matches']:
        errs = [m[2] for m in r_konly['matches']]
        log(f"  위치 오차: mean={np.mean(errs):.4f}, max={max(errs):.4f}")

    # ��+mono ���가
    log(f"\n  [평가: ��+mono (모노드���미 필터 후)]")
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
        "label": label, "rank": rank, "N": N,
        "n_actual": len(actual_zeros),
        "konly": r_konly, "kmono": r_kmono,
    }


# ━━━━━━━━━━━━��━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 실행
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

if __name__ == "__main__":
    log("=" * 70)
    log("[Project RDL] C-351 — EC 블라인드 영점 예측 확장")
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

    # 최종 요약
    log(f"\n\n{'='*70}")
    log(f"최종 요약")
    log(f"{'='*70}")

    log(f"\n  {'곡선':>10s}  {'rank':>4s}  {'N':>6s}  {'영점':>4s}  "
        f"{'κ-P':>5s}  {'κ-R':>5s}  {'κ-F1':>5s}  "
        f"{'κ+m-P':>6s}  {'κ+m-R':>6s}  {'κ+m-F1':>7s}")
    log(f"  {'-'*75}")

    for r in all_results:
        ko = r["konly"]
        km = r["kmono"]
        log(f"  {r['label']:>10s}  {r['rank']:>4d}  {r['N']:>6d}  {r['n_actual']:>4d}  "
            f"{ko['precision']:>5.3f}  {ko['recall']:>5.3f}  {ko['f1']:>5.3f}  "
            f"{km['precision']:>6.3f}  {km['recall']:>6.3f}  {km['f1']:>7.3f}")

    # 기존 결과 (#52, #54)와 통합 비��
    log(f"\n  [기존 결과 포함 통합표]")
    log(f"  {'곡선':>10s}  {'rank':>4s}  {'N':>6s}  {'method':>10s}  "
        f"{'κ+m-P':>6s}  {'κ+m-R':>6s}  {'κ+m-F1':>7s}")
    log(f"  {'-'*60}")

    # #52: 11a1 mpmath
    log(f"  {'11a1':>10s}  {'0':>4s}  {'11':>6s}  {'mpmath#52':>10s}  "
        f"{'0.941':>6s}  {'0.941':>6s}  {'0.941':>7s}")
    # #54: 37a1 mpmath
    log(f"  {'37a1':>10s}  {'1':>4s}  {'37':>6s}  {'mpmath#54':>10s}  "
        f"{'1.000':>6s}  {'0.957':>6s}  {'0.978':>7s}")

    for r in all_results:
        km = r["kmono"]
        log(f"  {r['label']:>10s}  {r['rank']:>4d}  {r['N']:>6d}  {'PARI-C351':>10s}  "
            f"{km['precision']:>6.3f}  {km['recall']:>6.3f}  {km['f1']:>7.3f}")

    # 판정
    log(f"\n  [판정]")
    all_f1 = [r["kmono"]["f1"] for r in all_results]
    all_rec = [r["kmono"]["recall"] for r in all_results]
    if all(f >= 0.80 for f in all_f1) and all(r >= 0.85 for r in all_rec):
        log(f"  ★★★★★ 양성 — 전곡선 F1≥0.80, Recall≥0.85")
    elif all(f >= 0.60 for f in all_f1):
        log(f"  ★★★ 조건부 양성 — F1≥0.60이나 기준 미달 곡선 존재")
    else:
        log(f"  ❌ 음성 — F1<0.60 곡선 존재")

    total_time = sum(1 for _ in [])  # placeholder
    log(f"\n총 소요: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    save()
    log(f"\n결과 저���: {OUT_PATH}")
