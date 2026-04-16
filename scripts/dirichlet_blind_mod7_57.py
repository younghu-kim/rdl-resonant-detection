"""
=============================================================================
[Project RDL] 결과 #57 — 디리클레 L(s,χ) 블라인드 영점 예측 (χ mod 7)
=============================================================================

목표:
  - GL(1) 디리클레 L-함수에서 블라인드 영점 예측
  - #40/#40b에서 AFE 정확성 검증 완료된 χ mod 7 사용
  - GL(2) 블라인드(#52/#54/#56)와 동일한 파이프라인: κ 스윕 → find_peaks → 모노드로미 필터

프로토콜:
  1. σ = 1/2 + 0.03 따라 t ∈ [5, 40], Δt=0.2 스윕
  2. κ(0.53 + it) = |Λ'/Λ|² 계산
  3. scipy.signal.find_peaks 극대점 추출 → 후보 목록
  4. 각 후보에 모노드로미 측정 (radius=0.4, n_steps=64)
  5. 비교 기준:
     - κ-only: 피크 후보 전체
     - κ+mono: mono/π > 1.5 필터 후
  6. LMFDB 영점과 사후 비교 (검증)
  7. Precision, Recall, F1 측정

성공 기준:
  - Recall ≥ 0.7  [필수]
  - 위치 오차 < 0.5  [필수]
  - Precision ≥ 0.5  [양성 근거]
  - F1 ≥ 0.6  [양성 근거]
  - κ+mono가 κ-only보다 Precision +10%p  [양성 근거, bonus]

주의:
  ★    σ_crit = 1/2 (디리클레 GL(1), GL(2) σ=1과 다름)
  ★★   χ mod 7: 복소 지표 (order 6), a=1 (χ(-1)=-1, 홀수)
  ★★★  mpmath.dirichlet(s, chi) 직접 사용
  ★★★★ log-space monodromy (#40b 수정 방식, underflow 방지)
  ★★★★★ 블라인드 = LMFDB 영점 예측 단계에서 사용 금지

결과 파일: results/dirichlet_blind_mod7_57.txt
=============================================================================
"""

import sys
import os
import time
import numpy as np
import cmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))
sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

import mpmath

try:
    from scipy.signal import find_peaks as scipy_find_peaks
    SCIPY_OK = True
except ImportError:
    SCIPY_OK = False

OUTFILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "dirichlet_blind_mod7_57.txt"
)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 실험 파라미터
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

SIGMA_CRIT     = 0.5        # ⚠️ 디리클레 GL(1) 임계선 σ=1/2
DELTA_OFFSET   = 0.03       # κ 측정 오프셋
T_MIN          = 5.0        # 수학자 지시: t∈[5,40]
T_MAX          = 40.0
DT             = 0.2        # 수학자 지시: Δt=0.2
MONO_RADIUS    = 0.4        # 수학자 지시: 표준 파라미터
MONO_STEPS     = 64         # 수학자 지시: 64단계
MONO_THRESHOLD = 1.5        # mono/π > 1.5 → 영점 판정
MATCH_TOL      = 0.5        # 수학자 지시: 표준 파라미터

# DPS
DPS = 80                    # #40b 표준

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# χ mod 7 지표 정의 (#40b에서 정확히 복사)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
_w6 = cmath.exp(2j * cmath.pi / 6)
_chi7_raw = [0, 1, _w6**2, _w6**1, _w6**4, _w6**5, _w6**3]
_chi7 = [mpmath.mpc(c.real, c.imag) for c in _chi7_raw]

CHAR_MOD7 = {
    'chi': _chi7,
    'q': 7,
    'a': 1,          # χ(-1)=-1 → a=1 (홀수 지표)
    'label': 'χ₇ (mod 7, 복소, order 6)',
}

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# LMFDB 영점 (검증용으로만! 예측 단계에서 사용 금지)
# #40에서 하드코딩된 136개 영점 중 t∈[5,40] 범위 추출
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
LMFDB_ZEROS_CHI7 = np.array([
    17.16141654, 19.65122423, 21.65252507, 24.15466454, 25.68439459,
    27.13547138, 28.64452754, 31.84774664, 32.74739063, 34.35044122,
    36.44798399, 37.37953168, 39.57467224,
])
# → 이 변수는 evaluate_predictions()에서만 사용

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 로깅
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

lines = []

def log(msg=""):
    print(msg, flush=True)
    lines.append(str(msg))

def flush_to_file():
    with open(OUTFILE, "w") as f:
        f.write("\n".join(lines) + "\n")

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 디리클레 L-함수: completed_L (#40b에서 정확히 복사)
# Λ(s, χ) = (q/π)^{s/2} Γ((s+a)/2) L(s, χ)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def completed_L(s, ci):
    """완비 L-함수: Λ(s, χ) = (q/π)^{s/2} Γ((s+a)/2) L(s, χ)"""
    q = mpmath.mpf(ci['q'])
    a = mpmath.mpf(ci['a'])
    L_val = mpmath.dirichlet(s, ci['chi'])
    gamma_val = mpmath.gamma((s + a) / 2)
    prefactor = mpmath.power(q / mpmath.pi, s / 2)
    return prefactor * gamma_val * L_val


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# log-space arg 계산 (#40b 수정 방식, underflow 방지)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def log_arg_lambda(s, ci):
    """
    arg(Λ(s,χ)) = Im[log(Λ)]
    log(Λ) = (s/2)*log(q/π) + loggamma((s+a)/2) + log(L(s,χ))
    반환: (Im(log_Λ), L_abs) — L_abs는 0 체크용
    """
    q = mpmath.mpf(ci['q'])
    a = mpmath.mpf(ci['a'])

    # 항 1: (s/2)*log(q/π)
    log_prefactor = (s / 2) * mpmath.log(q / mpmath.pi)

    # 항 2: loggamma((s+a)/2)
    log_gamma = mpmath.loggamma((s + a) / 2)

    # 항 3: log(L(s,χ))
    L_val = mpmath.dirichlet(s, ci['chi'])
    L_abs = abs(L_val)
    if L_abs < mpmath.mpf(10)**(-mpmath.mp.dps + 15):
        return None, float(L_abs)

    log_L = mpmath.log(L_val)

    log_Lambda = log_prefactor + log_gamma + log_L
    return float(mpmath.im(log_Lambda)), float(L_abs)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# κ 계산 (log-space 미분)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def curvature_dirichlet(sigma, t, ci, h=1e-6):
    """κ(σ+it) = |Λ'/Λ|² via 중앙차분"""
    saved_dps = mpmath.mp.dps
    mpmath.mp.dps = DPS
    try:
        s_mp = mpmath.mpc(sigma, t)
        h_mp = mpmath.mpf(str(h))
        L0 = completed_L(s_mp, ci)
        if abs(L0) < mpmath.mpf(10)**(-DPS + 15):
            return 1e12  # 영점 근처
        Lp = completed_L(s_mp + h_mp, ci)
        Lm = completed_L(s_mp - h_mp, ci)
        conn = (Lp - Lm) / (2 * h_mp * L0)
        k = float(abs(conn)**2)
        return k if np.isfinite(k) else 1e12
    except Exception as e:
        print(f"    WARNING curvature t={t:.4f}: {e}", flush=True)
        return 0.0
    finally:
        mpmath.mp.dps = saved_dps


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 모노드로미 (log-space, #40b 수정 방식)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def monodromy_dirichlet(t_center, ci, radius=MONO_RADIUS, n_steps=MONO_STEPS):
    """
    Λ(s,χ) 주위 폐곡선 모노드로미 (log-space arg 누적).
    중심: σ_crit + it_center.
    Returns: |mono|/π (영점이면 ≈2.0, 비영점이면 ≈0.0), None이면 실패
    """
    saved_dps = mpmath.mp.dps
    mpmath.mp.dps = DPS
    try:
        s_center = mpmath.mpc(SIGMA_CRIT, t_center)
        total_delta = 0.0
        prev_arg = None
        skip_count = 0

        for k in range(n_steps + 1):
            theta = 2 * mpmath.pi * k / n_steps
            s = s_center + radius * mpmath.exp(1j * theta)
            try:
                curr_arg, L_abs = log_arg_lambda(s, ci)
                if curr_arg is None:
                    skip_count += 1
                    continue
            except Exception as e:
                skip_count += 1
                continue

            if prev_arg is not None:
                diff = curr_arg - prev_arg
                while diff > np.pi:
                    diff -= 2 * np.pi
                while diff < -np.pi:
                    diff += 2 * np.pi
                total_delta += diff
            prev_arg = curr_arg

        if skip_count > n_steps // 4:
            return None  # 너무 많은 skip
        return abs(total_delta) / np.pi
    finally:
        mpmath.mp.dps = saved_dps


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 0: 함수방정식 사전 검증
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def verify_functional_equation(ci):
    """함수방정식 Λ(s) ≈ ε̄·Λ(1-s̄) 검증 (2개 테스트 점)"""
    log("\n[Step 0] 함수방정식 사전 검증")
    saved_dps = mpmath.mp.dps
    mpmath.mp.dps = DPS

    test_points = [
        mpmath.mpc('0.7', '5'),
        mpmath.mpc('0.4', '12'),
    ]

    all_ok = True
    for s in test_points:
        Ls = completed_L(s, ci)
        s_conj_bar = 1 - mpmath.conj(s)
        Ls_conj = completed_L(s_conj_bar, ci)
        # Λ(s) / Λ(1-s̄) should have |ratio| ≈ 1
        if abs(Ls_conj) > 0:
            ratio = Ls / Ls_conj
            rel_err = abs(abs(ratio) - 1)
            log(f"  Λ({mpmath.nstr(s, 4)}) / Λ(1-s̄): |ratio|={float(abs(ratio)):.10f}, "
                f"rel_err={float(rel_err):.2e}")
            if rel_err > 1e-5:
                log(f"  ⚠️ 함수방정식 오차 크다!")
                all_ok = False
        else:
            log(f"  ⚠️ Λ(1-s̄) ≈ 0")
            all_ok = False

    mpmath.mp.dps = saved_dps
    log(f"  판정: {'✅ PASS' if all_ok else '❌ FAIL'}")
    return all_ok


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 1: 곡률 스윕 (블라인드 — LMFDB 영점 미사용)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def sweep_curvature(ci):
    """σ = SIGMA_CRIT + DELTA_OFFSET 따라 t 스윕. κ 배열 + t 배열 반환."""
    sigma_sweep = SIGMA_CRIT + DELTA_OFFSET
    ts = np.arange(T_MIN, T_MAX + DT/2, DT)
    kappas = []

    log(f"\n[Step 1] κ 스윕: σ={sigma_sweep:.3f}, t ∈ [{T_MIN}, {T_MAX}], Δt={DT}")
    log(f"  총 {len(ts)}점 계산 예정")

    t0 = time.time()
    for i, t in enumerate(ts):
        k = curvature_dirichlet(sigma_sweep, t, ci)
        kappas.append(k)
        if (i+1) % 30 == 0:
            elapsed = time.time() - t0
            eta = elapsed / (i+1) * (len(ts) - i - 1)
            log(f"  ... {i+1}/{len(ts)} ({(i+1)/len(ts)*100:.0f}%), "
                f"경과={elapsed:.0f}초, 잔여≈{eta:.0f}초, "
                f"t={t:.2f}, κ={k:.2f}")
            flush_to_file()

    kappas = np.array(kappas)
    log(f"\n  스윕 완료: {time.time()-t0:.1f}초")
    log(f"  κ 통계: min={kappas.min():.4f}, max={kappas.max():.2f}, "
        f"median={np.median(kappas):.4f}, mean={kappas.mean():.4f}")

    return ts, kappas


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 2: κ 피크 추출 (블라인드)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def extract_kappa_peaks(ts, kappas):
    """κ 극대점 추출. scipy.find_peaks 사용."""
    log(f"\n[Step 2] κ 피크 추출")

    log_k = np.log1p(kappas)

    if SCIPY_OK:
        median_logk = np.median(log_k)
        prominence_thresh = max(0.3, median_logk * 0.3)
        peaks, props = scipy_find_peaks(log_k, prominence=prominence_thresh, distance=2)
        log(f"  scipy.find_peaks: prominence≥{prominence_thresh:.3f}, distance≥2")
        log(f"  → {len(peaks)}개 피크 발견")

        # 피크 부족 시 threshold 낮추기
        if len(peaks) < 5:
            prominence_thresh2 = max(0.1, median_logk * 0.1)
            peaks2, _ = scipy_find_peaks(log_k, prominence=prominence_thresh2, distance=2)
            log(f"  threshold 낮춤 (prominence≥{prominence_thresh2:.3f}): {len(peaks2)}개 피크")
            if len(peaks2) > len(peaks):
                peaks = peaks2

        peak_ts = ts[peaks]
        peak_kappas = kappas[peaks]
    else:
        # 수동 극대점 탐색
        log("  scipy 없음 — 수동 극대점 탐색")
        peaks_list = []
        for i in range(1, len(kappas)-1):
            if kappas[i] > kappas[i-1] and kappas[i] > kappas[i+1]:
                peaks_list.append(i)
        peaks = np.array(peaks_list, dtype=int)
        peak_ts = ts[peaks]
        peak_kappas = kappas[peaks]
        log(f"  → {len(peaks)}개 극대점 발견")

    if len(peaks) == 0:
        log("  ⚠️ 피크 0개 — κ 분포 이상!")
        return np.array([]), np.array([])

    log(f"\n  피크 목록 (κ-only 후보):")
    for i, (tp, kp) in enumerate(zip(peak_ts, peak_kappas)):
        log(f"    후보 {i+1:2d}: t={tp:.4f}, κ={kp:.4f}")

    return peak_ts, peak_kappas


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 3: 모노드로미 필터 (블라인드)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def apply_monodromy_filter(peak_ts, peak_kappas, ci):
    """각 κ 피크 후보에 모노드로미 측정. mono/π > MONO_THRESHOLD → 영점 후보."""
    log(f"\n[Step 3] 모노드로미 필터 (radius={MONO_RADIUS}, n_steps={MONO_STEPS})")
    log(f"  기준: mono/π > {MONO_THRESHOLD}")
    log(f"  대상: {len(peak_ts)}개 후보")

    mono_values = []
    filtered_ts = []
    filtered_kappas = []

    for i, (t_cand, k_cand) in enumerate(zip(peak_ts, peak_kappas)):
        # radius: 다른 후보와의 거리 고려
        other_ts = [pt for j, pt in enumerate(peak_ts) if j != i]
        if other_ts:
            nearest = min(abs(t_cand - ot) for ot in other_ts)
            radius = min(MONO_RADIUS, nearest * 0.45)
            radius = max(radius, 0.1)
        else:
            radius = MONO_RADIUS

        mono_pi = monodromy_dirichlet(t_cand, ci, radius=radius)
        if mono_pi is None:
            log(f"    [{i+1:2d}] t={t_cand:.4f}, κ={k_cand:.4f}, mono=FAIL (skip)")
            mono_values.append(np.nan)
            continue

        mono_values.append(mono_pi)
        marker = " ← ★ 영점" if mono_pi > MONO_THRESHOLD else " (비영점)"
        log(f"    [{i+1:2d}] t={t_cand:.4f}, κ={k_cand:.4f}, "
            f"r={radius:.3f}, mono/π={mono_pi:.4f}{marker}")

        if mono_pi > MONO_THRESHOLD:
            filtered_ts.append(t_cand)
            filtered_kappas.append(k_cand)

    mono_values = np.array(mono_values)
    filtered_ts = np.array(filtered_ts)
    filtered_kappas = np.array(filtered_kappas)

    log(f"\n  κ+mono 필터 결과: {len(filtered_ts)}개 후보 → 영점 예측")
    if len(filtered_ts) == 0:
        log("  ⚠️ 모노드로미 필터 후 0개 — threshold 조정 필요")

    return mono_values, filtered_ts, filtered_kappas


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 4: 평가 (LMFDB와 비교)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def evaluate_predictions(pred_ts, label, true_zeros=None):
    """예측 영점 목록 vs LMFDB 영점 비교. |t_pred - t_true| < MATCH_TOL → 매칭."""
    if true_zeros is None:
        true_zeros = LMFDB_ZEROS_CHI7

    n_true = len(true_zeros)
    pred_ts = np.array(pred_ts) if len(pred_ts) > 0 else np.array([])

    log(f"\n[평가: {label}]")
    log(f"  예측 {len(pred_ts)}개 vs 실제 {n_true}개 (매칭 허용 오차: {MATCH_TOL})")

    if len(pred_ts) == 0:
        log(f"  예측 0개 → P=0, R=0, F1=0")
        return 0.0, 0.0, 0.0, []

    # 매칭 (greedy: 가장 가까운 쌍)
    true_matched = [False] * n_true
    pred_matched = [False] * len(pred_ts)
    matched_pairs = []
    errors = []

    dists = []
    for pi, pt in enumerate(pred_ts):
        for ti, tt in enumerate(true_zeros):
            dists.append((abs(pt - tt), pi, ti, pt, tt))
    dists.sort()

    for dist, pi, ti, pt, tt in dists:
        if pred_matched[pi] or true_matched[ti]:
            continue
        if dist < MATCH_TOL:
            pred_matched[pi] = True
            true_matched[ti] = True
            matched_pairs.append((pt, tt, dist))
            errors.append(dist)

    tp = len(matched_pairs)
    fp = len(pred_ts) - tp
    fn = n_true - tp

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    log(f"  TP={tp}, FP={fp}, FN={fn}")
    log(f"  Precision={precision:.4f}, Recall={recall:.4f}, F1={f1:.4f}")

    if matched_pairs:
        log(f"  위치 오차: mean={np.mean(errors):.4f}, max={np.max(errors):.4f}")
        log(f"\n  매칭 상세:")
        for pt, tt, d in sorted(matched_pairs, key=lambda x: x[1]):
            log(f"    예측 {pt:.4f} ↔ 실제 {tt:.8f} (오차 {d:.4f})")
    else:
        log("  매칭 없음")

    # 미탐지 영점
    undetected = [true_zeros[i] for i in range(n_true) if not true_matched[i]]
    if undetected:
        log(f"\n  미탐지 실제 영점: {[f'{t:.4f}' for t in undetected]}")

    # 오탐 예측
    false_positives = [pred_ts[i] for i in range(len(pred_ts)) if not pred_matched[i]]
    if false_positives:
        log(f"  오탐 예측: {[f'{t:.4f}' for t in false_positives]}")

    return precision, recall, f1, matched_pairs


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# near/far κ 비율 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_near_far_ratio(ts, kappas, true_zeros):
    """영점 근처 (|t-γ|<0.3) vs 먼 곳 (|t-γ|>1.0) κ 중앙값 비율"""
    near_mask = np.zeros(len(ts), dtype=bool)
    far_mask = np.ones(len(ts), dtype=bool)
    for gz in true_zeros:
        dists = np.abs(ts - gz)
        near_mask |= (dists < 0.3)
        far_mask &= (dists > 1.0)

    near_vals = kappas[near_mask]
    far_vals = kappas[far_mask]

    if len(near_vals) > 0 and len(far_vals) > 0:
        near_med = np.median(near_vals)
        far_med = np.median(far_vals)
        if far_med > 0:
            ratio = near_med / far_med
        else:
            ratio = float('inf')
        return ratio, near_med, far_med, len(near_vals), len(far_vals)
    return 0.0, 0.0, 0.0, 0, 0


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    ci = CHAR_MOD7

    log("=" * 70)
    log("결과 #57 — 디리클레 L(s,χ) 블라인드 영점 예측 (χ mod 7)")
    log("=" * 70)
    log(f"시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"지표: {ci['label']}")
    log(f"conductor q={ci['q']}, a={ci['a']} (χ(-1)=-1, 홀수 지표)")
    log(f"SIGMA_CRIT={SIGMA_CRIT} (GL(1) 디리클레), DELTA_OFFSET={DELTA_OFFSET}")
    log(f"t ∈ [{T_MIN}, {T_MAX}], Δt={DT} → {int((T_MAX - T_MIN)/DT + 1)}점")
    log(f"LMFDB 알려진 영점: {len(LMFDB_ZEROS_CHI7)}개 (검증용으로만)")
    log(f"MONO_RADIUS={MONO_RADIUS}, MONO_STEPS={MONO_STEPS}, MONO_THRESHOLD={MONO_THRESHOLD}")
    log(f"MATCH_TOL={MATCH_TOL}")
    log(f"dps={DPS}")
    log(f"scipy 사용: {SCIPY_OK}")
    log()

    t_total_start = time.time()

    # ── Step 0: 함수방정식 검증 ──
    mpmath.mp.dps = DPS
    fe_ok = verify_functional_equation(ci)
    flush_to_file()
    if not fe_ok:
        log("⚠️ 함수방정식 검증 실패 — 계속 진행하지만 주의 필요")

    # ── Step 1: κ 스윕 (블라인드) ──
    ts, kappas = sweep_curvature(ci)
    flush_to_file()

    # ── Step 2: κ 피크 추출 (블라인드) ──
    peak_ts, peak_kappas = extract_kappa_peaks(ts, kappas)
    flush_to_file()

    if len(peak_ts) == 0:
        log("\n⚠️⚠️ κ 피크 0개 — 실험 중단!")
        flush_to_file()
        sys.exit(1)

    # ── Step 3: 모노드로미 필터 (블라인드) ──
    mono_values, filtered_ts, filtered_kappas = apply_monodromy_filter(peak_ts, peak_kappas, ci)
    flush_to_file()

    # ── Step 4: 평가 (LMFDB 사용) ──
    log("\n" + "=" * 70)
    log("[Step 4] 평가 — LMFDB 영점과 비교 (검증 단계)")
    log("=" * 70)
    log(f"\nLMFDB χ mod 7 영점 ({len(LMFDB_ZEROS_CHI7)}개, t∈[{T_MIN},{T_MAX}]):")
    for i, tz in enumerate(LMFDB_ZEROS_CHI7):
        log(f"  γ_{i+1:2d} = {tz:.8f}")

    # κ-only 평가
    P_konly, R_konly, F1_konly, pairs_konly = evaluate_predictions(
        peak_ts, "κ-only (피크 전체)", LMFDB_ZEROS_CHI7
    )
    flush_to_file()

    # κ+mono 평가
    P_kmono, R_kmono, F1_kmono, pairs_kmono = evaluate_predictions(
        filtered_ts, "κ+mono (모노드로미 필터 후)", LMFDB_ZEROS_CHI7
    )
    flush_to_file()

    # near/far κ 비율
    nf_ratio, near_med, far_med, n_near, n_far = compute_near_far_ratio(
        ts, kappas, LMFDB_ZEROS_CHI7
    )

    # ── 최종 보고 ──
    log("\n" + "=" * 70)
    log("최종 결과 요약")
    log("=" * 70)
    log()

    log("[ κ 스윕 ]")
    log(f"  t 범위: [{T_MIN}, {T_MAX}], Δt={DT}, 총 {len(ts)}점")
    log(f"  σ 오프셋: {SIGMA_CRIT} + {DELTA_OFFSET} = {SIGMA_CRIT + DELTA_OFFSET}")
    log(f"  κ 중앙값: {np.median(kappas):.4f}, 최대: {kappas.max():.2f}")
    log()

    log("[ near/far κ 비율 ]")
    log(f"  near (|t-γ|<0.3): median={near_med:.2f} ({n_near}점)")
    log(f"  far  (|t-γ|>1.0): median={far_med:.2f} ({n_far}점)")
    log(f"  ratio: {nf_ratio:.1f}×")
    log()

    log("[ 후보 비교 ]")
    log(f"  κ-only  후보: {len(peak_ts)}개")
    log(f"  κ+mono  후보: {len(filtered_ts)}개")
    log()

    log("[ 성능 비교 ]")
    log(f"{'기준':<15} {'Precision':>12} {'Recall':>10} {'F1':>8} {'TP':>5} {'FP':>5}")
    log("-" * 55)
    n_true = len(LMFDB_ZEROS_CHI7)
    tp_konly = len(pairs_konly)
    fp_konly = len(peak_ts) - tp_konly
    tp_kmono = len(pairs_kmono)
    fp_kmono = len(filtered_ts) - tp_kmono
    log(f"{'κ-only':<15} {P_konly:>12.4f} {R_konly:>10.4f} {F1_konly:>8.4f} "
        f"{tp_konly:>5} {fp_konly:>5}")
    log(f"{'κ+mono':<15} {P_kmono:>12.4f} {R_kmono:>10.4f} {F1_kmono:>8.4f} "
        f"{tp_kmono:>5} {fp_kmono:>5}")
    log()

    # Precision 개선량
    delta_P = P_kmono - P_konly
    log(f"  Precision 개선 (κ+mono - κ-only): {delta_P:+.4f} ({delta_P*100:+.1f}%p)")
    log()

    log("[ 성공 기준 판정 ]")
    criterion_results = []

    # 필수 기준
    r_recall = R_kmono >= 0.7
    log(f"  {'✅' if r_recall else '❌'} [필수] Recall ≥ 0.7: "
        f"{R_kmono:.4f} {'PASS' if r_recall else 'FAIL'}")
    criterion_results.append(('recall_필수', r_recall))

    # 위치 오차
    if pairs_kmono:
        max_err = max(d for _, _, d in pairs_kmono)
        r_pos = max_err < MATCH_TOL
    else:
        max_err = float('nan')
        r_pos = False
    log(f"  {'✅' if r_pos else '❌'} [필수] 위치 오차 < {MATCH_TOL}: "
        f"max={max_err:.4f} {'PASS' if r_pos else 'FAIL'}")
    criterion_results.append(('position_필수', r_pos))

    # 양성 근거 기준
    r_precision = P_kmono >= 0.5
    log(f"  {'✅' if r_precision else '❌'} [양성] Precision ≥ 0.5: "
        f"{P_kmono:.4f} {'PASS' if r_precision else 'FAIL'}")
    criterion_results.append(('precision_양성', r_precision))

    r_f1 = F1_kmono >= 0.6
    log(f"  {'✅' if r_f1 else '❌'} [양성] F1 ≥ 0.6: "
        f"{F1_kmono:.4f} {'PASS' if r_f1 else 'FAIL'}")
    criterion_results.append(('f1_양성', r_f1))

    r_improvement = delta_P > 0.10
    log(f"  {'✅' if r_improvement else '❌'} [양성] κ+mono Precision 개선 > +10%p: "
        f"{delta_P*100:+.1f}%p {'PASS' if r_improvement else 'FAIL'}")
    criterion_results.append(('improvement_양성', r_improvement))

    log()
    n_pass = sum(v for _, v in criterion_results)
    n_total = len(criterion_results)
    n_mandatory = sum(1 for k, v in criterion_results if '필수' in k)
    n_mandatory_pass = sum(v for k, v in criterion_results if '필수' in k)
    log(f"  통과: {n_pass}/{n_total} (필수 {n_mandatory_pass}/{n_mandatory})")
    log()

    # 최종 판정
    mandatory_ok = n_mandatory_pass == n_mandatory
    if mandatory_ok and n_pass >= 4:
        verdict = "★★★ 완전 양성"
    elif mandatory_ok and n_pass >= 3:
        verdict = "★★ 양성"
    elif mandatory_ok:
        verdict = "★ 조건부 양성 (필수 기준 충족, 양성 기준 부족)"
    else:
        verdict = "❌ 음성 (필수 기준 미충족)"

    log(f"최종 판정: {verdict}")
    log()

    # 블라인드 보편성 종합표
    log("=" * 70)
    log("블라인드 예측 보편성 종합표")
    log("=" * 70)
    log()
    log(f"{'L-함수':<20} {'유형':<15} {'σ_crit':<8} {'영점수':<8} {'TP':<5} {'FP':<5} {'P':<8} {'R':<8} {'F1':<8}")
    log("-" * 90)
    log(f"{'ζ(s) (#2)':<20} {'GL(1) Riemann':<15} {'0.5':<8} {'7':<8} {'7':<5} {'0':<5} {'1.000':<8} {'1.000':<8} {'1.000':<8}")
    log(f"{'11a1 (#52)':<20} {'GL(2) EC':<15} {'1.0':<8} {'17':<8} {'16':<5} {'1':<5} {'0.941':<8} {'0.941':<8} {'0.941':<8}")
    log(f"{'37a1 (#54)':<20} {'GL(2) EC':<15} {'1.0':<8} {'23':<8} {'22':<5} {'0':<5} {'1.000':<8} {'0.957':<8} {'0.978':<8}")
    log(f"{'Δ (#56)':<20} {'GL(2) Modular':<15} {'6.0':<8} {'8':<8} {'8':<5} {'0':<5} {'1.000':<8} {'1.000':<8} {'1.000':<8}")
    log(f"{'χ mod 7 (#57)':<20} {'GL(1) Dirichlet':<15} {'0.5':<8} {f'{len(LMFDB_ZEROS_CHI7)}':<8} {f'{tp_kmono}':<5} {f'{fp_kmono}':<5} {f'{P_kmono:.3f}':<8} {f'{R_kmono:.3f}':<8} {f'{F1_kmono:.3f}':<8}")
    log()

    elapsed = time.time() - t_total_start
    log(f"총 소요 시간: {elapsed:.1f}초 ({elapsed/60:.1f}분)")
    log(f"완료 시각: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    flush_to_file()
    print(f"\n결과 저장: {OUTFILE}", flush=True)


if __name__ == '__main__':
    main()
