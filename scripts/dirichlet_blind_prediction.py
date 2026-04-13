"""
=============================================================================
[Project RDL] 디리클레 L-함수 블라인드 영점 예측
=============================================================================
다발 예측 실험 #3의 디리클레 확장:
  훈련 구간에서 학습한 곡률/모노드로미 패턴으로 테스트 구간의 영점을 예측.

대상 지표:
  - χ mod 3, χ mod 4, χ mod 5

방법:
  1. 훈련 구간 t∈[10, 25]에서 곡률/모노드로미 프로파일 분석
  2. 테스트 구간 t∈[25, 40]에서 3가지 기준으로 영점 예측:
     (b) 곡률: κ(1/2+it) 극대점
     (c) 이중 기준: κ 극대 + 모노드로미 |Δarg| > π/2
     (d) 모노드로미만: |Δarg| > π/2
  3. Precision / Recall / F1 비교 (지표별, 방법별)

예측: 이중 기준(c)의 정밀도가 단일 기준보다 유의미하게 높음
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified'))

mpmath.mp.dps = 80

# ─── 설정 ───
T_TRAIN = (10.0, 25.0)
T_TEST = (25.0, 40.0)
N_POINTS = 800  # 테스트 구간 스캔 밀도

# ─── 지표 정의 (dirichlet_bundle_verification.py와 동일) ───
CHARACTERS = {
    'χ mod 3': {
        'chi': [0, 1, -1],
        'q': 3, 'a': 1,
        'label': 'χ₃ (mod 3)',
    },
    'χ mod 4': {
        'chi': [0, 1, 0, -1],
        'q': 4, 'a': 1,
        'label': 'χ₄ (mod 4)',
    },
    'χ mod 5': {
        'chi': [0, 1, 1j, -1j, -1],
        'q': 5, 'a': 1,
        'label': 'χ₅ (mod 5)',
    },
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 핵심 함수
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def completed_L(s, char_info):
    """완비 L-함수: Λ(s, χ) = (q/π)^{s/2} Γ((s+a)/2) L(s, χ)"""
    q = mpmath.mpf(char_info['q'])
    a = mpmath.mpf(char_info['a'])
    chi = char_info['chi']
    L_val = mpmath.dirichlet(s, chi)
    gamma_val = mpmath.gamma((s + a) / 2)
    prefactor = mpmath.power(q / mpmath.pi, s / 2)
    return prefactor * gamma_val * L_val


def connection_L(s, char_info):
    """접속 L(s) = Λ'/Λ (수치 미분)"""
    h = mpmath.mpf(1) / mpmath.mpf(10**20)
    Lambda_val = completed_L(s, char_info)
    if abs(Lambda_val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
        return mpmath.mpc(1e10, 0)
    Lambda_d = (completed_L(s + h, char_info) - completed_L(s - h, char_info)) / (2 * h)
    return Lambda_d / Lambda_val


def curvature_at(s, char_info):
    """곡률 κ = |Λ'/Λ|²"""
    L = connection_L(s, char_info)
    return float(abs(L)**2)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 영점 탐색
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def find_zeros_dirichlet(char_info, t_min, t_max, n_scan=2000):
    """임계선 위 Re(Λ) / Im(Λ) 부호 변화로 영점 탐색 + findroot 정밀화"""
    ts = np.linspace(t_min, t_max, n_scan)
    zeros = []

    # Re(Λ) 부호 변화
    prev_re = None
    prev_t = None
    for t in ts:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        val = completed_L(s, char_info)
        curr_re = mpmath.re(val)

        if prev_re is not None and prev_re * curr_re < 0:
            try:
                def f_real(t_var):
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                    return mpmath.re(completed_L(sv, char_info))
                t_zero = float(mpmath.findroot(f_real, (prev_t, t)))
                if not zeros or abs(t_zero - zeros[-1]) > 0.1:
                    zeros.append(t_zero)
            except Exception:
                pass
        prev_re = curr_re
        prev_t = t

    # Im(Λ) 부호 변화 (보충)
    prev_im = None
    prev_t = None
    for t in ts:
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        val = completed_L(s, char_info)
        curr_im = mpmath.im(val)

        if prev_im is not None and prev_im * curr_im < 0:
            try:
                def f_imag(t_var):
                    sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(t_var)
                    return mpmath.im(completed_L(sv, char_info))
                t_zero = float(mpmath.findroot(f_imag, (prev_t, t)))
                sv = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t_zero))
                if abs(completed_L(sv, char_info)) < mpmath.mpf('1e-10'):
                    if not any(abs(t_zero - z) < 0.1 for z in zeros):
                        zeros.append(t_zero)
            except Exception:
                pass
        prev_im = curr_im
        prev_t = t

    zeros.sort()
    return np.array(zeros)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 프로파일 계산
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def compute_profiles(char_info, t_min, t_max, n_points):
    """구간 내 곡률 + 위상 프로파일 계산"""
    ts = np.linspace(t_min, t_max, n_points)
    kappas = np.zeros(n_points)
    args = np.zeros(n_points)

    for i, t in enumerate(ts):
        s = mpmath.mpf('0.5') + 1j * mpmath.mpf(str(t))
        val = completed_L(s, char_info)

        if abs(val) < mpmath.mpf(10)**(-mpmath.mp.dps + 10):
            kappas[i] = 1e10
            args[i] = 0
            continue

        L = connection_L(s, char_info)
        kappas[i] = float(abs(L)**2)
        args[i] = float(mpmath.arg(val))

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


def find_curvature_peaks(ts, kappas, min_prominence=5.0):
    """곡률 극대점 찾기 (중앙값 대비 min_prominence배 이상)"""
    peaks = []
    med = np.median(kappas)
    for i in range(1, len(kappas) - 1):
        if kappas[i] > kappas[i-1] and kappas[i] > kappas[i+1]:
            if kappas[i] > min_prominence * med:
                peaks.append(i)
    return np.array(peaks, dtype=int)


def evaluate_predictions(pred_t, true_t, tolerance=0.5):
    """Precision, Recall, F1 계산"""
    if len(pred_t) == 0:
        return 0.0, 0.0, 0.0
    if len(true_t) == 0:
        return 0.0, 0.0, 0.0

    tp = 0
    matched = set()
    for pt in pred_t:
        dists = np.abs(true_t - pt)
        best = np.argmin(dists)
        if dists[best] < tolerance and best not in matched:
            tp += 1
            matched.add(best)

    precision = tp / len(pred_t)
    recall = tp / len(true_t)
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    return precision, recall, f1


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 단일 지표 블라인드 예측
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def blind_predict_single(char_name, char_info, out_file):
    """단일 지표에 대해 블라인드 영점 예측 수행"""
    label = char_info['label']

    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"{'='*70}")

    out_file.write(f"\n{'='*70}\n")
    out_file.write(f"  {label}\n")
    out_file.write(f"{'='*70}\n")

    # ── 정답 영점 (테스트 구간) ──
    print(f"\n  정답 영점 수집: t∈[{T_TEST[0]}, {T_TEST[1]}]")
    true_zeros = find_zeros_dirichlet(char_info, T_TEST[0], T_TEST[1])
    print(f"    영점 수: {len(true_zeros)}")
    for i, tz in enumerate(true_zeros):
        print(f"      #{i+1}: t = {tz:.6f}")

    out_file.write(f"\n정답 영점 ({len(true_zeros)}개):\n")
    for i, tz in enumerate(true_zeros):
        out_file.write(f"  #{i+1}: t = {tz:.6f}\n")

    if len(true_zeros) == 0:
        print("  (영점 미발견, 건너뜀)")
        out_file.write("  영점 미발견, 건너뜀\n")
        return None

    # ── 곡률/위상 프로파일 (테스트 구간) ──
    print(f"\n  곡률/위상 프로파일 계산: {N_POINTS}점...")
    ts, kappas, args = compute_profiles(char_info, T_TEST[0], T_TEST[1], N_POINTS)
    monos = compute_monodromy_profile(ts, args)

    # ── 방법 (b): 곡률 극대점 ──
    print(f"\n  [방법 b] 곡률 극대점")
    peak_idx = find_curvature_peaks(ts, kappas, min_prominence=5.0)
    pred_b = ts[peak_idx] if len(peak_idx) > 0 else np.array([])
    p_b, r_b, f1_b = evaluate_predictions(pred_b, true_zeros)
    print(f"    예측: {len(pred_b)}개, P={p_b:.1%}, R={r_b:.1%}, F1={f1_b:.3f}")

    # ── 방법 (c): 이중 기준 (곡률 + 모노드로미) ──
    print(f"\n  [방법 c] 이중 기준: 곡률 극대 + |mono| > π/2")
    mono_threshold = np.pi / 2
    dual_idx = []
    for idx in peak_idx:
        window = monos[max(0, idx-5):min(len(monos), idx+6)]
        if np.max(np.abs(window)) > mono_threshold:
            dual_idx.append(idx)

    pred_c = ts[np.array(dual_idx, dtype=int)] if dual_idx else np.array([])
    p_c, r_c, f1_c = evaluate_predictions(pred_c, true_zeros)
    print(f"    예측: {len(pred_c)}개, P={p_c:.1%}, R={r_c:.1%}, F1={f1_c:.3f}")

    # ── 방법 (d): 모노드로미만 ──
    print(f"\n  [방법 d] 모노드로미 점프 위치만 (|Δarg| > π/2)")
    jump_idx = np.where(np.abs(monos) > mono_threshold)[0]
    pred_d = ts[jump_idx]
    p_d, r_d, f1_d = evaluate_predictions(pred_d, true_zeros)
    print(f"    예측: {len(pred_d)}개, P={p_d:.1%}, R={r_d:.1%}, F1={f1_d:.3f}")

    # ── 비교 표 ──
    print(f"\n  {'방법':<30} {'N_pred':>8} {'P':>8} {'R':>8} {'F1':>8}")
    print(f"  {'-'*62}")
    print(f"  {'(b) 곡률 극대':<30} {len(pred_b):>8} {p_b:>8.1%} {r_b:>8.1%} {f1_b:>8.3f}")
    print(f"  {'(c) 곡률+모노드로미':<30} {len(pred_c):>8} {p_c:>8.1%} {r_c:>8.1%} {f1_c:>8.3f}")
    print(f"  {'(d) 모노드로미만':<30} {len(pred_d):>8} {p_d:>8.1%} {r_d:>8.1%} {f1_d:>8.3f}")

    out_file.write(f"\n{'방법':<30} {'N_pred':>8} {'P':>8} {'R':>8} {'F1':>8}\n")
    out_file.write("-" * 62 + "\n")
    out_file.write(f"{'(b) 곡률 극대':<30} {len(pred_b):>8} {p_b:>8.1%} {r_b:>8.1%} {f1_b:>8.3f}\n")
    out_file.write(f"{'(c) 곡률+모노드로미':<30} {len(pred_c):>8} {p_c:>8.1%} {r_c:>8.1%} {f1_c:>8.3f}\n")
    out_file.write(f"{'(d) 모노드로미만':<30} {len(pred_d):>8} {p_d:>8.1%} {r_d:>8.1%} {f1_d:>8.3f}\n")

    # ── 곡률 극대점 상세 ──
    out_file.write(f"\n곡률 극대점 상세:\n")
    out_file.write(f"{'t':>12} {'κ':>14} {'|mono|/π':>12} {'true_zero?':>12}\n")
    out_file.write("-" * 55 + "\n")
    for idx in peak_idx[:30]:
        t = ts[idx]
        k = kappas[idx]
        m = np.max(np.abs(monos[max(0,idx-5):min(len(monos),idx+6)]))
        dist = np.min(np.abs(true_zeros - t)) if len(true_zeros) > 0 else 999
        is_true = "YES" if dist < 0.5 else "no"
        out_file.write(f"{t:>12.4f} {k:>14.2e} {m/np.pi:>12.4f} {is_true:>12}\n")

    return {
        'n_zeros': len(true_zeros),
        'methods': {
            'curvature': (p_b, r_b, f1_b, len(pred_b)),
            'dual': (p_c, r_c, f1_c, len(pred_c)),
            'monodromy': (p_d, r_d, f1_d, len(pred_d)),
        }
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 메인
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    out_path = os.path.expanduser(
        '~/Desktop/gdl_unified/results/dirichlet_blind_prediction.txt'
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    print("=" * 70)
    print("디리클레 L-함수 블라인드 영점 예측")
    print(f"훈련: t∈[{T_TRAIN[0]}, {T_TRAIN[1]}], "
          f"테스트: t∈[{T_TEST[0]}, {T_TEST[1]}]")
    print(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}")
    print("=" * 70)

    with open(out_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("디리클레 L-함수 블라인드 영점 예측\n")
        f.write(f"훈련: t∈[{T_TRAIN[0]}, {T_TRAIN[1]}], "
                f"테스트: t∈[{T_TEST[0]}, {T_TEST[1]}]\n")
        f.write(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"정밀도: {mpmath.mp.dps} 자릿수, 스캔 밀도: {N_POINTS}점\n")
        f.write("=" * 70 + "\n")

        all_results = {}
        for char_name, char_info in CHARACTERS.items():
            t0 = time.time()
            result = blind_predict_single(char_name, char_info, f)
            elapsed = time.time() - t0
            all_results[char_name] = result
            f.write(f"\n  (소요 시간: {elapsed:.0f}초)\n")
            print(f"\n  [소요 시간: {elapsed:.0f}초]")

        # ── 교차 비교 표 ──
        f.write(f"\n\n{'='*70}\n")
        f.write("지표간 비교 (이중 기준 c 기준)\n")
        f.write(f"{'='*70}\n\n")

        print(f"\n\n{'='*70}")
        print("지표간 비교 (이중 기준 c 기준)")
        print(f"{'='*70}")

        header = f"{'지표':<15} {'영점수':>8} {'N_pred':>8} {'P':>8} {'R':>8} {'F1':>8}"
        print(header)
        print("-" * 60)
        f.write(header + "\n")
        f.write("-" * 60 + "\n")

        for char_name, r in all_results.items():
            if r is None:
                line = f"{char_name:<15} {'N/A':>8}"
            else:
                p, rec, f1, n_pred = r['methods']['dual']
                line = (f"{char_name:<15} {r['n_zeros']:>8} {n_pred:>8} "
                        f"{p:>8.1%} {rec:>8.1%} {f1:>8.3f}")
            print(line)
            f.write(line + "\n")

    print(f"\n결과 저장: {out_path}")
    print("완료.")


if __name__ == '__main__':
    main()
