#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GUE/Poisson 중간점 비국소 대조 실험 — 결과 #20 후보
RMT 기원 검증: GUE vs Poisson vs ζ/디리클레 (결과 #18-#19 기준)

수학자 지시 사이클 #28
"""

import numpy as np
from scipy.stats import spearmanr, ttest_1samp
import datetime
import os
import sys

# ── 설정 ─────────────────────────────────────────────────────────────────
np.random.seed(42)

N_GUE    = 500    # GUE 행렬 크기
N_RUNS   = 20     # 반복 수
N_POISSON = 234   # Poisson 점 수 (ζ 영점 수와 동일)
CAP_VAL  = 1e18   # κ cap
BULK_FRAC = 0.8   # bulk 비율 (중앙 80%)
WINDOW   = 25     # unfolding 국소 창 폭 (±WINDOW)

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, '..', 'results')
OUTPUT_FILE = os.path.join(RESULTS_DIR, 'gue_poisson_midpoint_comparison.txt')

# ── 이전 결과 기준값 ─────────────────────────────────────────────────────
ZR_RHO_A  = -0.031   # ζ (결과 #18) ρ(a)
ZR_RHO_B  = -0.654   # ζ (결과 #18) ρ(b)
ZR_PCORR  = -0.666   # ζ (결과 #18) 편상관
DIR_RHO_A = -0.025   # χ₃~χ₅ 평균 (결과 #19) ρ(a)
DIR_RHO_B = -0.594   # χ₃~χ₅ 평균 (결과 #19) ρ(b)
DIR_PCORR = -0.606   # χ₃~χ₅ 평균 (결과 #19) 편상관


# ── 유틸 함수 ────────────────────────────────────────────────────────────
def partial_corr_spearman(x, y, z):
    """Spearman 편상관 ρ(x,y|z) — 표준 1차 편상관 공식"""
    rxy, _ = spearmanr(x, y)
    rxz, _ = spearmanr(x, z)
    ryz, _ = spearmanr(y, z)
    denom = np.sqrt(max((1 - rxz**2) * (1 - ryz**2), 1e-30))
    return (rxy - rxz * ryz) / denom


def compute_kappa(eigs_all, idx_n, idx_np1, m):
    """
    κ = R(m)² — NN 제외 resolvent 제곱
    NN 항(idx_n, idx_np1)은 중간점에서 정확히 상쇄되므로 명시적 제외.
    """
    mask = np.ones(len(eigs_all), dtype=bool)
    mask[idx_n]   = False
    mask[idx_np1] = False
    diffs = m - eigs_all[mask]
    # 수치 안전 필터 (극히 작은 분모 방지)
    safe = np.abs(diffs) > 1e-15
    if not np.any(safe):
        return 0.0, False
    R_m = np.sum(1.0 / diffs[safe])
    kappa_raw = R_m ** 2
    capped = kappa_raw > CAP_VAL
    return min(kappa_raw, CAP_VAL), capped


# ── GUE 실험 ─────────────────────────────────────────────────────────────
def run_gue_experiment():
    print("=" * 60)
    print(f"GUE 실험: {N_GUE}×{N_GUE} 행렬, {N_RUNS}회 반복")
    print("=" * 60)
    sys.stdout.flush()

    all_runs = []

    for run in range(N_RUNS):
        rng = np.random.RandomState(run)
        A   = rng.randn(N_GUE, N_GUE) + 1j * rng.randn(N_GUE, N_GUE)
        H   = (A + A.conj().T) / 2.0
        eigs = np.sort(np.linalg.eigvalsh(H).real)  # 정렬된 실수 고유값

        n       = len(eigs)
        i_start = int((1 - BULK_FRAC) / 2 * n)   # =50
        i_end   = int((1 + BULK_FRAC) / 2 * n)   # =450

        gaps_list, gap_next_list, kappa_list = [], [], []
        n_capped = 0

        for i in range(i_start + 1, i_end - 1):
            # 쌍: (eigs[i-1], eigs[i])  →  midpoint m
            # gap_next: (eigs[i], eigs[i+1])
            raw_gap      = eigs[i]   - eigs[i-1]
            raw_gap_next = eigs[i+1] - eigs[i]
            m            = (eigs[i-1] + eigs[i]) / 2.0

            # 국소 평균 간격 (unfolding)
            w_s = max(0, i - WINDOW)
            w_e = min(n - 1, i + WINDOW)
            local_mean = np.mean(np.diff(eigs[w_s:w_e + 1]))
            if local_mean <= 0:
                continue

            gap      = raw_gap / local_mean
            gap_next = raw_gap_next / local_mean

            # κ 계산
            kappa, capped = compute_kappa(eigs, i - 1, i, m)
            if capped:
                n_capped += 1

            gaps_list.append(gap)
            gap_next_list.append(gap_next)
            kappa_list.append(kappa)

        if len(kappa_list) < 10:
            print(f"  Run {run+1:2d}: WARNING — 데이터 부족 ({len(kappa_list)}개)")
            sys.stdout.flush()
            continue

        gaps_arr     = np.array(gaps_list)
        gap_next_arr = np.array(gap_next_list)
        kappa_arr    = np.array(kappa_list)

        n_pairs  = len(kappa_arr)
        cap_pct  = 100.0 * n_capped / n_pairs

        rho_a, _   = spearmanr(kappa_arr, gaps_arr)
        rho_b, p_b = spearmanr(kappa_arr, gap_next_arr)
        pcorr = partial_corr_spearman(kappa_arr, gap_next_arr, gaps_arr)

        all_runs.append({
            'n': n_pairs, 'cap_pct': cap_pct,
            'rho_a': rho_a, 'rho_b': rho_b, 'p_b': p_b, 'pcorr': pcorr,
            'kappa_mean': float(np.mean(kappa_arr)),
            'kappa_std':  float(np.std(kappa_arr)),
        })

        print(f"  Run {run+1:2d}: n={n_pairs}, cap={cap_pct:.1f}%, "
              f"ρ(a)={rho_a:+.3f}, ρ(b)={rho_b:+.3f}, "
              f"p(b)={p_b:.2e}, 편상관={pcorr:+.3f}")
        sys.stdout.flush()

    return all_runs


# ── Poisson 실험 ──────────────────────────────────────────────────────────
def run_poisson_experiment():
    print("\n" + "=" * 60)
    print(f"Poisson 실험: {N_POISSON}점 균등분포, {N_RUNS}회 반복")
    print("=" * 60)
    sys.stdout.flush()

    all_runs = []

    for run in range(N_RUNS):
        rng = np.random.RandomState(run + 1000)
        pts = np.sort(rng.uniform(0, N_POISSON, N_POISSON))

        mean_gap_global = np.mean(np.diff(pts))  # ≈ 1

        gaps_list, gap_next_list, kappa_list = [], [], []
        n_capped = 0

        for i in range(1, N_POISSON - 1):
            m            = (pts[i-1] + pts[i]) / 2.0
            raw_gap      = pts[i]   - pts[i-1]
            raw_gap_next = pts[i+1] - pts[i]

            gap      = raw_gap / mean_gap_global
            gap_next = raw_gap_next / mean_gap_global

            # κ: NN 제외 resolvent^2
            mask = np.ones(N_POISSON, dtype=bool)
            mask[i-1] = False
            mask[i]   = False
            diffs = m - pts[mask]
            safe  = np.abs(diffs) > 1e-15
            if not np.any(safe):
                continue
            R_m       = np.sum(1.0 / diffs[safe])
            kappa_raw = R_m ** 2
            capped    = kappa_raw > CAP_VAL
            if capped:
                n_capped += 1
            kappa = min(kappa_raw, CAP_VAL)

            gaps_list.append(gap)
            gap_next_list.append(gap_next)
            kappa_list.append(kappa)

        gaps_arr     = np.array(gaps_list)
        gap_next_arr = np.array(gap_next_list)
        kappa_arr    = np.array(kappa_list)

        n_pairs = len(kappa_arr)
        cap_pct = 100.0 * n_capped / n_pairs if n_pairs > 0 else 0.0

        rho_a, _   = spearmanr(kappa_arr, gaps_arr)
        rho_b, p_b = spearmanr(kappa_arr, gap_next_arr)
        pcorr = partial_corr_spearman(kappa_arr, gap_next_arr, gaps_arr)

        all_runs.append({
            'n': n_pairs, 'cap_pct': cap_pct,
            'rho_a': rho_a, 'rho_b': rho_b, 'p_b': p_b, 'pcorr': pcorr,
        })

        print(f"  Run {run+1:2d}: n={n_pairs}, cap={cap_pct:.1f}%, "
              f"ρ(a)={rho_a:+.3f}, ρ(b)={rho_b:+.3f}, "
              f"p(b)={p_b:.2e}, 편상관={pcorr:+.3f}")
        sys.stdout.flush()

    return all_runs


# ── 요약 & 판정 ───────────────────────────────────────────────────────────
def summarize_runs(all_runs, label):
    rho_a_arr = np.array([r['rho_a'] for r in all_runs])
    rho_b_arr = np.array([r['rho_b'] for r in all_runs])
    pcorr_arr = np.array([r['pcorr'] for r in all_runs])
    cap_arr   = np.array([r['cap_pct'] for r in all_runs])

    _, p_ttest = ttest_1samp(rho_b_arr, 0.0)

    print(f"\n{label} 종합 ({len(all_runs)}회):")
    print(f"  ρ(a) gap:       {np.mean(rho_a_arr):+.3f} ± {np.std(rho_a_arr):.3f}")
    print(f"  ρ(b) gap_next:  {np.mean(rho_b_arr):+.3f} ± {np.std(rho_b_arr):.3f}  "
          f"(t-test vs 0: p={p_ttest:.2e})")
    print(f"  편상관:          {np.mean(pcorr_arr):+.3f} ± {np.std(pcorr_arr):.3f}")
    print(f"  cap%:            {np.mean(cap_arr):.1f}%")
    sys.stdout.flush()

    return {
        'label':    label,
        'n_runs':   len(all_runs),
        'rho_a':    (float(np.mean(rho_a_arr)), float(np.std(rho_a_arr))),
        'rho_b':    (float(np.mean(rho_b_arr)), float(np.std(rho_b_arr))),
        'pcorr':    (float(np.mean(pcorr_arr)), float(np.std(pcorr_arr))),
        'cap':      float(np.mean(cap_arr)),
        'p_b_ttest': float(p_ttest),
        'all_runs': all_runs,
    }


def determine_verdict(gue_s, poi_s):
    g = abs(gue_s['rho_b'][0])
    p = abs(poi_s['rho_b'][0])
    if g > 0.3 and p < 0.15:
        return "양성-RMT: 반발 통계(RMT)의 보편적 성질 (L-함수 고유 아님)"
    elif g < 0.15 and p < 0.15:
        return "양성-고유: L-함수 산술 구조 고유의 현상 (RMT와 무관)"
    elif g > 0.3 and p > 0.3:
        return "음성: Hadamard 급수의 수학적 필연 (점 과정 무관)"
    else:
        return f"미결: GUE|ρ|={g:.3f}, Poisson|ρ|={p:.3f} (기준 경계)"


# ── 메인 실행 ─────────────────────────────────────────────────────────────
if __name__ == '__main__':
    t0 = datetime.datetime.now()
    print(f"[{t0.strftime('%H:%M:%S')}] GUE/Poisson 중간점 비국소 실험 시작")
    print(f"  GUE: {N_GUE}×{N_GUE}, {N_RUNS}회 / Poisson: {N_POISSON}점, {N_RUNS}회")
    print()
    sys.stdout.flush()

    gue_runs = run_gue_experiment()
    poi_runs = run_poisson_experiment()

    gue_s = summarize_runs(gue_runs, "GUE")
    poi_s = summarize_runs(poi_runs, "Poisson")

    verdict = determine_verdict(gue_s, poi_s)

    t1 = datetime.datetime.now()
    elapsed = (t1 - t0).total_seconds()

    # ── 결과 파일 작성 ──────────────────────────────────────────────────
    lines = []
    lines.append("=" * 72)
    lines.append("GUE/Poisson 중간점 비국소 대조 실험 — 결과 #20 후보")
    lines.append(f"실행 시각: {t0.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"소요 시간: {elapsed:.1f}초")
    lines.append("=" * 72)
    lines.append("")
    lines.append("실험 설계:")
    lines.append(f"  GUE    : {N_GUE}×{N_GUE} Hermitian 행렬, bulk {int(BULK_FRAC*100)}%, {N_RUNS}회")
    lines.append(f"  Poisson: {N_POISSON}점 균등분포 [0, {N_POISSON}), {N_RUNS}회")
    lines.append(f"  κ      : R(m)² (NN 명시적 제외, resolvent 제곱)")
    lines.append(f"  unfolding: 국소 평균 간격 sliding window ±{WINDOW}")
    lines.append(f"  cap    : {CAP_VAL:.0e} (초과 비율 확인)")
    lines.append("")
    lines.append("-" * 72)
    lines.append("GUE 개별 실행 결과:")
    for i, r in enumerate(gue_runs):
        lines.append(f"  Run {i+1:2d}: n={r['n']}, cap={r['cap_pct']:.1f}%, "
                     f"ρ(a)={r['rho_a']:+.3f}, ρ(b)={r['rho_b']:+.3f}, "
                     f"p(b)={r['p_b']:.2e}, 편상관={r['pcorr']:+.3f}")
    lines.append("")
    lines.append("Poisson 개별 실행 결과:")
    for i, r in enumerate(poi_runs):
        lines.append(f"  Run {i+1:2d}: n={r['n']}, cap={r['cap_pct']:.1f}%, "
                     f"ρ(a)={r['rho_a']:+.3f}, ρ(b)={r['rho_b']:+.3f}, "
                     f"p(b)={r['p_b']:.2e}, 편상관={r['pcorr']:+.3f}")
    lines.append("")
    lines.append("=" * 72)
    lines.append("비교표:")
    lines.append("")
    header = f"{'점 과정':<22} | {'ρ(a) gap':>10} | {'ρ(b) gap_next':>15} | {'편상관':>10} | 판정"
    lines.append(header)
    lines.append("-" * 76)

    # ζ (결과 #18)
    lines.append(f"{'ζ (결과 #18)':<22} | {ZR_RHO_A:+10.3f} | {ZR_RHO_B:+15.3f} | {ZR_PCORR:+10.3f} | ✅ 양성")
    # χ₃~χ₅ 평균 (결과 #19)
    lines.append(f"{'χ₃~χ₅ 평균 (결과 #19)':<22} | {DIR_RHO_A:+10.3f} | {DIR_RHO_B:+15.3f} | {DIR_PCORR:+10.3f} | ✅ 양성")

    # GUE
    ga_m, ga_s = gue_s['rho_a']
    gb_m, gb_s = gue_s['rho_b']
    gc_m, gc_s = gue_s['pcorr']
    gue_judge = "★ RMT" if abs(gb_m) > 0.3 else ("★ 없음" if abs(gb_m) < 0.15 else "? 경계")
    lines.append(f"{'GUE (20회 평균)':<22} | {ga_m:+.3f}±{ga_s:.3f} | "
                 f"{gb_m:+.3f}±{gb_s:.3f}       | "
                 f"{gc_m:+.3f}±{gc_s:.3f} | {gue_judge}")

    # Poisson
    pa_m, pa_s = poi_s['rho_a']
    pb_m, pb_s = poi_s['rho_b']
    pc_m, pc_s = poi_s['pcorr']
    poi_judge = "★ 없음" if abs(pb_m) < 0.15 else ("★ 존재" if abs(pb_m) > 0.3 else "? 경계")
    lines.append(f"{'Poisson (20회)':<22} | {pa_m:+.3f}±{pa_s:.3f} | "
                 f"{pb_m:+.3f}±{pb_s:.3f}       | "
                 f"{pc_m:+.3f}±{pc_s:.3f} | {poi_judge}")

    lines.append("")
    lines.append("-" * 72)
    lines.append("통계 검정 (t-test, H₀: 20회 평균 ρ(b)=0):")
    lines.append(f"  GUE     : t-test p = {gue_s['p_b_ttest']:.3e}")
    lines.append(f"  Poisson : t-test p = {poi_s['p_b_ttest']:.3e}")
    lines.append("")
    lines.append("인프라 확인:")
    lines.append(f"  GUE cap%     : {gue_s['cap']:.1f}%  {'✅ OK' if gue_s['cap'] < 5 else '⚠️ 초과'}")
    lines.append(f"  Poisson cap% : {poi_s['cap']:.1f}%  {'✅ OK' if poi_s['cap'] < 5 else '⚠️ 초과'}")
    nn_ok = abs(ga_m) < 0.3
    lines.append(f"  GUE ρ(a)|nn 상쇄 : |{ga_m:.3f}| < 0.3  {'✅ OK' if nn_ok else '⚠️ FAIL'}")
    lines.append("")
    lines.append("성공 기준 분류:")
    lines.append(f"  GUE  |ρ(b)| = {abs(gb_m):.3f}  (기준: >0.3 = RMT, <0.15 = 고유)")
    lines.append(f"  Poisson |ρ(b)| = {abs(pb_m):.3f}  (기준: <0.15 = RMT 또는 고유 확인)")
    lines.append("")
    lines.append(f"★ 최종 판정: {verdict}")
    lines.append("")
    lines.append("=" * 72)

    result_text = '\n'.join(lines)
    print()
    print(result_text)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.write(result_text + '\n')
    print(f"\n결과 파일 저장: {OUTPUT_FILE}")
    sys.stdout.flush()
