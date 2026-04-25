"""
=============================================================================
[Project RDL] 다발 예측 실험 #3: 위상적 영점 블라인드 예측 (C-300)
=============================================================================
핵심 질문: 곡률(κ) + 모노드로미(±π) 이중 기준이 FP를 제거하면서
          recall을 유지하는가?

방법:
  테스트 구간 [150, 200]에서 3가지 기준으로 영점 예측:
    (a) 곡률 극대점만 (baseline)
    (b) 모노드로미 점프만 (|Δarg| > π/2)
    (c) 이중 기준: 곡률 극대 + 모노드로미 폐곡선 |mono| > π/2
  세 방법의 precision/recall/F1 비교 + 기존 NN 결과와 대조

성공 기준 (수학자 지정):
  - Recall ≥ 95%, Precision ≥ 80%, F1 ≥ 87%
  - 강양성: P ≥ 90% + R ≥ 98%

참조: blind_zero_prediction.txt (NN |F₂| 극소 — P≈19.6%, R=100%)
=============================================================================
"""

import sys, os, time
import numpy as np
import mpmath

sys.path.insert(0, os.path.expanduser('~/Desktop/gdl_unified/scripts'))

mpmath.mp.dps = 80  # t>100에서 필수

from bundle_utils import (
    xi_func, connection_zeta, curvature_zeta,
    monodromy_contour, find_zeros_zeta,
    compute_curvature_profile, find_curvature_peaks,
    compute_monodromy_profile, evaluate_predictions,
)

# ─── 설정 ───
T_TEST = (150.0, 200.0)
N_POINTS = 2000       # Δt ≈ 0.025 (tolerance=0.5보다 훨씬 세밀)
TOLERANCE = 0.5       # 적중 판정 (기존과 동일)
MONO_THRESHOLD = np.pi / 2  # 모노드로미 점프 기준
KAPPA_PROMINENCE = 5.0      # 곡률 극대 기준 (중앙값 대비)
CONTOUR_RADIUS = 0.5        # 모노드로미 폐곡선 반지름


def main():
    results_path = os.path.expanduser(
        '~/Desktop/gdl_unified/results/blind_topology_c301.txt'
    )
    analysis_path = os.path.expanduser(
        '~/Desktop/gdl_unified/outputs/analysis/bundle_prediction_blind_topology.txt'
    )

    print("=" * 70)
    print("C-300: 위상적 영점 블라인드 예측 (이중 기준)")
    print("=" * 70)
    print(f"구간: t∈[{T_TEST[0]}, {T_TEST[1]}], {N_POINTS}점")
    print(f"dps={mpmath.mp.dps}, tolerance={TOLERANCE}")
    t_start = time.time()

    # ─── 1. 정답 영점 수집 ───
    print(f"\n[1] 정답 영점 수집...")
    true_zeros = find_zeros_zeta(T_TEST[0], T_TEST[1])
    print(f"  영점 수: {len(true_zeros)}")
    if len(true_zeros) == 0:
        print("⚠️ 영점 0개 — 탐색 로직 점검 필요")
        return
    if len(true_zeros) < 20:
        print(f"⚠️ 영점 {len(true_zeros)}개 — 기준 20개 미만, 주의")
    for i, z in enumerate(true_zeros[:5]):
        print(f"  γ_{i+1} = {z:.6f}")
    if len(true_zeros) > 5:
        print(f"  ... (총 {len(true_zeros)}개)")

    # ─── 2. 곡률 + 위상 프로파일 계산 ───
    print(f"\n[2] 곡률/위상 프로파일 계산 ({N_POINTS}점)...")
    ts, kappas, args = compute_curvature_profile(T_TEST[0], T_TEST[1], N_POINTS)
    monos = compute_monodromy_profile(ts, args)
    elapsed_profile = time.time() - t_start
    print(f"  완료 ({elapsed_profile:.1f}s)")
    print(f"  κ 통계: median={np.median(kappas[kappas < 1e9]):.2e}, "
          f"max={np.max(kappas):.2e}")

    # ─── 3. 방법 (a): 곡률 극대점만 (baseline) ───
    print(f"\n[3] 방법 (a): 곡률 극대점 (prominence>{KAPPA_PROMINENCE}×median)")
    peak_idx = find_curvature_peaks(ts, kappas, min_prominence=KAPPA_PROMINENCE)
    pred_a = ts[peak_idx]
    p_a, r_a, f1_a = evaluate_predictions(pred_a, true_zeros, TOLERANCE)
    print(f"  예측: {len(pred_a)}개, P={p_a:.1%}, R={r_a:.1%}, F1={f1_a:.3f}")

    # ─── 4. 방법 (b): 모노드로미 점프만 ───
    print(f"\n[4] 방법 (b): 모노드로미 점프 (|Δarg| > π/2)")
    jump_idx = np.where(np.abs(monos) > MONO_THRESHOLD)[0]
    # 클러스터링: 연속된 점프는 하나로 (최대 κ 점 선택)
    jump_clusters = []
    if len(jump_idx) > 0:
        cluster = [jump_idx[0]]
        for j in range(1, len(jump_idx)):
            if jump_idx[j] - jump_idx[j-1] <= 3:  # 3 index 이내 = 같은 영점
                cluster.append(jump_idx[j])
            else:
                # 클러스터 대표: 최대 |mono| 점
                best = cluster[np.argmax(np.abs(monos[cluster]))]
                jump_clusters.append(best)
                cluster = [jump_idx[j]]
        best = cluster[np.argmax(np.abs(monos[cluster]))]
        jump_clusters.append(best)

    pred_b = ts[np.array(jump_clusters, dtype=int)] if jump_clusters else np.array([])
    p_b, r_b, f1_b = evaluate_predictions(pred_b, true_zeros, TOLERANCE)
    print(f"  원시 점프: {len(jump_idx)}개, 클러스터: {len(jump_clusters)}개")
    print(f"  예측: {len(pred_b)}개, P={p_b:.1%}, R={r_b:.1%}, F1={f1_b:.3f}")

    # ─── 5. 방법 (c): 이중 기준 (곡률 극대 + 폐곡선 모노드로미) ───
    print(f"\n[5] 방법 (c): 이중 기준 (κ 극대 + 폐곡선 모노드로미)")
    print(f"  곡률 극대점 {len(peak_idx)}개에 대해 폐곡선 모노드로미 계산...")

    dual_results = []  # (t, κ, mono, is_dual)
    dual_idx = []
    n_mono_pass = 0

    for i, idx in enumerate(peak_idx):
        t_val = ts[idx]
        k_val = kappas[idx]

        # 폐곡선 모노드로미 계산 (핵심: eps 차분이 아닌 적분)
        mono_val = monodromy_contour(t_val, radius=CONTOUR_RADIUS, n_steps=64)

        is_pass = abs(mono_val) > MONO_THRESHOLD
        dual_results.append((t_val, k_val, mono_val, is_pass))

        if is_pass:
            dual_idx.append(idx)
            n_mono_pass += 1

        if (i + 1) % 10 == 0 or i == len(peak_idx) - 1:
            print(f"    {i+1}/{len(peak_idx)}: t={t_val:.4f}, "
                  f"κ={k_val:.2e}, mono/π={mono_val/np.pi:.4f}, "
                  f"pass={'✅' if is_pass else '❌'}")

    pred_c = ts[np.array(dual_idx, dtype=int)] if dual_idx else np.array([])
    p_c, r_c, f1_c = evaluate_predictions(pred_c, true_zeros, TOLERANCE)
    print(f"  모노드로미 통과: {n_mono_pass}/{len(peak_idx)}")
    print(f"  예측: {len(pred_c)}개, P={p_c:.1%}, R={r_c:.1%}, F1={f1_c:.3f}")

    # ─── 6. 비교 표 ───
    elapsed_total = time.time() - t_start
    print(f"\n{'=' * 70}")
    print(f"비교 요약 (총 {elapsed_total:.1f}s)")
    print(f"{'=' * 70}")
    print(f"정답 영점: {len(true_zeros)}개, 구간: t∈[{T_TEST[0]}, {T_TEST[1]}]")
    print(f"\n{'방법':<35} {'N_pred':>7} {'P':>8} {'R':>8} {'F1':>8}")
    print("-" * 70)
    print(f"{'(a) 곡률 극대만':<35} {len(pred_a):>7} {p_a:>8.1%} {r_a:>8.1%} {f1_a:>8.3f}")
    print(f"{'(b) 모노드로미 점프만':<35} {len(pred_b):>7} {p_b:>8.1%} {r_b:>8.1%} {f1_b:>8.3f}")
    print(f"{'(c) 이중기준 (κ+mono)':<35} {len(pred_c):>7} {p_c:>8.1%} {r_c:>8.1%} {f1_c:>8.3f}")
    print(f"{'(ref) NN |F₂| (기존)':<35} {'138':>7} {'19.6%':>8} {'100%':>8} {'0.327':>8}")

    # ─── 7. 성공 기준 판정 ───
    print(f"\n성공 기준 판정 (이중 기준 c):")
    criteria = {
        'Recall ≥ 95%': r_c >= 0.95,
        'Precision ≥ 80%': p_c >= 0.80,
        'F1 ≥ 87%': f1_c >= 0.87,
    }
    strong = p_c >= 0.90 and r_c >= 0.98

    for name, met in criteria.items():
        print(f"  {'✅' if met else '❌'} {name}: {met}")
    print(f"  {'🏆' if strong else '⚠️'} 강양성 (P≥90% + R≥98%): {strong}")

    all_met = all(criteria.values())
    if strong:
        verdict = "★★★★★ 강양성"
    elif all_met:
        verdict = "★★★★ 양성"
    elif sum(criteria.values()) >= 2:
        verdict = "★★★ 부분 양성"
    else:
        verdict = "★★ 음성"
    print(f"\n종합 판정: {verdict}")

    # ─── 8. 상세 분석 ───
    print(f"\n곡률 극대점 상세:")
    print(f"{'t':>10} {'κ':>12} {'mono/π':>10} {'dual?':>6} {'dist':>8} {'TP?':>5}")
    print("-" * 55)
    for t_val, k_val, mono_val, is_pass in dual_results:
        dist = np.min(np.abs(true_zeros - t_val))
        is_tp = "YES" if dist < TOLERANCE else "no"
        tag = '✅' if is_pass else '❌'
        print(f"{t_val:>10.4f} {k_val:>12.2e} {mono_val/np.pi:>10.4f} "
              f"{tag:>6} {dist:>8.4f} {is_tp:>5}")

    # FP 분석
    print(f"\n[FP 분석 — 이중 기준]")
    fp_count = 0
    fp_mono_count = 0
    for t_val, k_val, mono_val, is_pass in dual_results:
        dist = np.min(np.abs(true_zeros - t_val))
        if dist >= TOLERANCE:  # FP in 곡률
            fp_count += 1
            if is_pass:  # FP가 이중기준도 통과
                fp_mono_count += 1
                print(f"  ⚠️ FP 잔존: t={t_val:.4f}, κ={k_val:.2e}, mono/π={mono_val/np.pi:.4f}")
            else:
                print(f"  ✅ FP 제거: t={t_val:.4f}, κ={k_val:.2e}, mono/π={mono_val/np.pi:.4f}")
    print(f"  곡률 FP: {fp_count}개, 이중기준 잔존 FP: {fp_mono_count}개")
    if fp_count > 0:
        print(f"  FP 제거율: {(fp_count - fp_mono_count) / fp_count:.1%}")

    # ─── 9. 결과 저장 ───
    os.makedirs(os.path.dirname(results_path), exist_ok=True)
    os.makedirs(os.path.dirname(analysis_path), exist_ok=True)

    content = []
    content.append("=" * 70)
    content.append("C-300: 위상적 영점 블라인드 예측 (이중 기준)")
    content.append(f"날짜: {time.strftime('%Y-%m-%d %H:%M')}")
    content.append(f"구간: t∈[{T_TEST[0]}, {T_TEST[1]}], {N_POINTS}점, dps={mpmath.mp.dps}")
    content.append(f"정답 영점: {len(true_zeros)}개")
    content.append(f"곡률 prominence: {KAPPA_PROMINENCE}×median, 모노 threshold: π/2")
    content.append(f"폐곡선 반지름: {CONTOUR_RADIUS}, tolerance: {TOLERANCE}")
    content.append(f"총 시간: {elapsed_total:.1f}s")
    content.append("=" * 70)
    content.append("")
    content.append(f"{'방법':<35} {'N_pred':>7} {'P':>8} {'R':>8} {'F1':>8}")
    content.append("-" * 70)
    content.append(f"{'(a) 곡률 극대만':<35} {len(pred_a):>7} {p_a:>8.1%} {r_a:>8.1%} {f1_a:>8.3f}")
    content.append(f"{'(b) 모노드로미 점프만':<35} {len(pred_b):>7} {p_b:>8.1%} {r_b:>8.1%} {f1_b:>8.3f}")
    content.append(f"{'(c) 이중기준 (κ+mono)':<35} {len(pred_c):>7} {p_c:>8.1%} {r_c:>8.1%} {f1_c:>8.3f}")
    content.append(f"{'(ref) NN |F₂| (기존)':<35} {'138':>7} {'19.6%':>8} {'100%':>8} {'0.327':>8}")
    content.append("")
    content.append(f"종합 판정: {verdict}")
    content.append(f"  Recall ≥ 95%: {'✅' if r_c >= 0.95 else '❌'} ({r_c:.1%})")
    content.append(f"  Precision ≥ 80%: {'✅' if p_c >= 0.80 else '❌'} ({p_c:.1%})")
    content.append(f"  F1 ≥ 87%: {'✅' if f1_c >= 0.87 else '❌'} ({f1_c:.3f})")
    content.append(f"  강양성 (P≥90%+R≥98%): {'🏆' if strong else '⚠️'}")
    content.append("")

    # FP 분석
    content.append("FP 분석:")
    content.append(f"  곡률 FP: {fp_count}개")
    content.append(f"  이중기준 잔존 FP: {fp_mono_count}개")
    if fp_count > 0:
        content.append(f"  FP 제거율: {(fp_count - fp_mono_count) / fp_count:.1%}")
    content.append("")

    # 상세 테이블
    content.append("곡률 극대점 상세:")
    content.append(f"{'t':>10} {'κ':>12} {'mono/π':>10} {'dual?':>6} {'dist':>8} {'TP?':>5}")
    content.append("-" * 55)
    for t_val, k_val, mono_val, is_pass in dual_results:
        dist = np.min(np.abs(true_zeros - t_val))
        is_tp = "YES" if dist < TOLERANCE else "no"
        tag = "✅" if is_pass else "❌"
        content.append(f"{t_val:>10.4f} {k_val:>12.2e} {mono_val/np.pi:>10.4f} "
                       f"{tag:>6} {dist:>8.4f} {is_tp:>5}")
    content.append("")

    # 정답 영점 목록
    content.append("정답 영점:")
    for i, z in enumerate(true_zeros):
        content.append(f"  γ_{i+1} = {z:.6f}")

    text = "\n".join(content)
    for path in [results_path, analysis_path]:
        with open(path, 'w') as f:
            f.write(text)
        print(f"\n저장: {path}")

    print("\n완료.")


if __name__ == '__main__':
    main()
