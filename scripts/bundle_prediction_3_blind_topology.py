"""
=============================================================================
[Project RDL] 다발 예측 실험 #3: 위상적 영점 예측 (블라인드)
=============================================================================
핵심 질문: 다발 기하학적 정보만으로 미지 영점을 예측할 수 있는가?

방법:
  1. 훈련 구간 [100, 150]에서 모노드로미 패턴 학습
  2. 테스트 구간 [150, 200]에서 3가지 기준으로 영점 예측:
     (a) 기존: |F₂| 극소점 (baseline)
     (b) 곡률: κ(1/2+it) 극대점
     (c) 이중 기준: κ 극대 + 모노드로미 |Δarg| > π/2
  3. 세 방법의 precision/recall 비교

예측: 이중 기준(c)의 정밀도가 (a), (b)보다 유의미하게 높음
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))

mpmath.mp.dps = 80  # 높은 t에서 ξ underflow 방지

# ─── 설정 ───
T_TRAIN = (14.0, 30.0)
T_TEST = (30.0, 50.0)
N_POINTS = 1000


# ─── 기본 함수 ───

def xi_func(s):
    half = mpmath.mpf('0.5')
    return half * s * (s - 1) * mpmath.power(mpmath.pi, -s/2) * mpmath.gamma(s/2) * mpmath.zeta(s)

def L_func(s):
    """L(s) = ξ'/ξ via numerical differentiation with adaptive step"""
    h = mpmath.mpf(1) / mpmath.mpf(10**20)
    xi_val = xi_func(s)
    xi_abs = abs(xi_val)
    # 영점 근방: 매우 높은 곡률
    if xi_abs < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e10, 0)
    xi_d = (xi_func(s + h) - xi_func(s - h)) / (2 * h)
    return xi_d / xi_val

def compute_curvature_profile(t_min, t_max, n_points):
    """t 구간에서 임계선 위 곡률 프로파일 계산"""
    ts = np.linspace(t_min, t_max, n_points)
    kappas = np.zeros(n_points)
    args = np.zeros(n_points)

    for i, t in enumerate(ts):
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        xi_val = xi_func(s)
        xi_abs = abs(xi_val)

        if xi_abs < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
            kappas[i] = 1e10
            args[i] = 0
            continue

        L = L_func(s)
        kappas[i] = float(abs(L)**2)
        args[i] = float(mpmath.arg(xi_val))

    return ts, kappas, args


def compute_monodromy_profile(ts, args):
    """위상 배열에서 모노드로미 프로파일 계산"""
    monos = np.zeros(len(ts))
    for i in range(1, len(ts)):
        delta = args[i] - args[i-1]
        while delta > np.pi:
            delta -= 2 * np.pi
        while delta < -np.pi:
            delta += 2 * np.pi
        monos[i] = delta
    return monos


def find_curvature_peaks(ts, kappas, min_prominence=10.0):
    """곡률 극대점 찾기"""
    peaks = []
    for i in range(1, len(kappas) - 1):
        if kappas[i] > kappas[i-1] and kappas[i] > kappas[i+1]:
            if kappas[i] > min_prominence * np.median(kappas):
                peaks.append(i)
    return np.array(peaks, dtype=int)


def evaluate_predictions(pred_t, true_t, tolerance=0.5):
    """precision, recall, F1 계산"""
    if len(pred_t) == 0:
        return 0.0, 0.0, 0.0

    tp = 0
    matched = set()
    for pt in pred_t:
        dists = np.abs(true_t - pt)
        best = np.argmin(dists)
        if dists[best] < tolerance and best not in matched:
            tp += 1
            matched.add(best)

    precision = tp / len(pred_t) if len(pred_t) > 0 else 0
    recall = tp / len(true_t) if len(true_t) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return precision, recall, f1


def get_zeros_in_range(t_min, t_max):
    """mpmath로 구간 내 리만 영점 구하기"""
    zeros = []
    n = 1
    while True:
        t = float(mpmath.zetazero(n).imag)
        if t > t_max:
            break
        if t >= t_min:
            zeros.append(t)
        n += 1
    return np.array(zeros)


# ─── 메인 ───

def main():
    out_path = os.path.expanduser(
        '~/Desktop/gdl_unified/outputs/analysis/bundle_prediction_blind_topology.txt'
    )

    print("=" * 70)
    print("다발 예측 실험 #3: 위상적 영점 예측 (블라인드)")
    print("=" * 70)

    # 테스트 구간 영점 (정답)
    print(f"\n정답 영점 수집: t∈[{T_TEST[0]}, {T_TEST[1]}]")
    true_zeros = get_zeros_in_range(T_TEST[0], T_TEST[1])
    print(f"  영점 수: {len(true_zeros)}")

    # 곡률 + 위상 프로파일 계산 (테스트 구간)
    print(f"\n곡률/위상 프로파일 계산: {N_POINTS}점...")
    ts, kappas, args = compute_curvature_profile(T_TEST[0], T_TEST[1], N_POINTS)
    monos = compute_monodromy_profile(ts, args)

    # ── 방법 (b): 곡률 극대점 ──
    print("\n[방법 b] 곡률 극대점")
    peak_idx = find_curvature_peaks(ts, kappas, min_prominence=5.0)
    pred_b = ts[peak_idx]
    p_b, r_b, f1_b = evaluate_predictions(pred_b, true_zeros)
    print(f"  예측: {len(pred_b)}개, P={p_b:.1%}, R={r_b:.1%}, F1={f1_b:.3f}")

    # ── 방법 (c): 이중 기준 (곡률 + 모노드로미) ──
    print("\n[방법 c] 이중 기준: 곡률 극대 + |mono| > π/2")
    mono_threshold = np.pi / 2
    dual_idx = []
    for idx in peak_idx:
        # 극대점 주변 ±5 index에서 모노드로미 검사
        window = monos[max(0, idx-5):min(len(monos), idx+6)]
        if np.max(np.abs(window)) > mono_threshold:
            dual_idx.append(idx)

    pred_c = ts[np.array(dual_idx)] if dual_idx else np.array([])
    p_c, r_c, f1_c = evaluate_predictions(pred_c, true_zeros)
    print(f"  예측: {len(pred_c)}개, P={p_c:.1%}, R={r_c:.1%}, F1={f1_c:.3f}")

    # ── 방법 (d): 모노드로미만 (점프 위치) ──
    print("\n[방법 d] 모노드로미 점프 위치만 (|Δarg| > π/2)")
    jump_idx = np.where(np.abs(monos) > mono_threshold)[0]
    pred_d = ts[jump_idx]
    p_d, r_d, f1_d = evaluate_predictions(pred_d, true_zeros)
    print(f"  예측: {len(pred_d)}개, P={p_d:.1%}, R={r_d:.1%}, F1={f1_d:.3f}")

    # 비교 표
    print("\n" + "=" * 70)
    print("비교 요약")
    print("=" * 70)
    print(f"{'방법':<30} {'N_pred':>8} {'P':>8} {'R':>8} {'F1':>8}")
    print("-" * 70)
    print(f"{'(b) 곡률 극대':<30} {len(pred_b):>8} {p_b:>8.1%} {r_b:>8.1%} {f1_b:>8.3f}")
    print(f"{'(c) 곡률+모노드로미':<30} {len(pred_c):>8} {p_c:>8.1%} {r_c:>8.1%} {f1_c:>8.3f}")
    print(f"{'(d) 모노드로미만':<30} {len(pred_d):>8} {p_d:>8.1%} {r_d:>8.1%} {f1_d:>8.3f}")

    # 곡률 통계
    print(f"\n곡률 통계:")
    print(f"  전체 median κ: {np.median(kappas):.2e}")
    print(f"  극대점 median κ: {np.median(kappas[peak_idx]):.2e}")
    print(f"  영점 근방 median κ: ", end="")
    zero_kappas = []
    for tz in true_zeros[:10]:
        idx = np.argmin(np.abs(ts - tz))
        zero_kappas.append(kappas[idx])
    print(f"{np.median(zero_kappas):.2e}")

    # 결과 저장
    with open(out_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("다발 예측 실험 #3: 위상적 영점 예측 (블라인드)\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"구간: t∈[{T_TEST[0]}, {T_TEST[1]}], {N_POINTS}점\n")
        f.write(f"정답 영점: {len(true_zeros)}개\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"{'방법':<30} {'N_pred':>8} {'P':>8} {'R':>8} {'F1':>8}\n")
        f.write("-" * 70 + "\n")
        f.write(f"{'(b) 곡률 극대':<30} {len(pred_b):>8} {p_b:>8.1%} {r_b:>8.1%} {f1_b:>8.3f}\n")
        f.write(f"{'(c) 곡률+모노드로미':<30} {len(pred_c):>8} {p_c:>8.1%} {r_c:>8.1%} {f1_c:>8.3f}\n")
        f.write(f"{'(d) 모노드로미만':<30} {len(pred_d):>8} {p_d:>8.1%} {r_d:>8.1%} {f1_d:>8.3f}\n")

        f.write(f"\n곡률 극대점 상세:\n")
        f.write(f"{'t':>12} {'κ':>14} {'|mono|/π':>12} {'true_zero?':>12}\n")
        f.write("-" * 55 + "\n")
        for idx in peak_idx[:30]:
            t = ts[idx]
            k = kappas[idx]
            m = np.max(np.abs(monos[max(0,idx-5):min(len(monos),idx+6)]))
            dist = np.min(np.abs(true_zeros - t))
            is_true = "YES" if dist < 0.5 else "no"
            f.write(f"{t:>12.4f} {k:>14.2e} {m/np.pi:>12.4f} {is_true:>12}\n")

    print(f"\n결과 저장: {out_path}")
    print("완료.")


if __name__ == '__main__':
    main()
