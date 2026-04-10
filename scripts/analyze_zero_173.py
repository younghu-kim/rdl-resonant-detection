#!/usr/bin/env python3
"""
analyze_zero_173.py — t=173.4115 근방 확대 관찰 (N1).

세 개 지표에서 "최악점"으로 교차한 이 영점에 국소 구조가 있는지 검사.
- F₂ 난이도: 0.1149 (t∈[100,200] 최대)
- arg Z_out 시드 비일관성 1위 (cstd=1.62)
- 고립실험 per-영점 MSE

비교 대조군:
- t=174.75 (직후 영점, 시드-가장-일관 1위, cstd=0.13)
- t=167.18 (세 번째 일관, cstd=0.21)
- 평균 영점 (나머지 47개)

출력: outputs/analysis/zero_173_analysis.{json,txt}
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import torch
from scipy.signal import hilbert

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "outputs" / "analysis"
OVERNIGHT = ROOT / "outputs" / "overnight"
OUT_JSON = ANALYSIS / "zero_173_analysis.json"
OUT_TXT = ANALYSIS / "zero_173_analysis.txt"

T_FOCUS = 173.4115  # 중심
T_CMP_EASY = 174.7542  # 대조군 (직후 영점, 가장 일관)
WINDOW = 1.5  # 중심 기준 ±1.5 확대


def main():
    lines = []
    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    log("=" * 72)
    log("  analyze_zero_173 — t=173.4115 근방 확대 관찰 (N1)")
    log("=" * 72)

    # ------------------------------------------------------------------
    # 데이터 로드
    # ------------------------------------------------------------------
    cache = torch.load(OVERNIGHT / "xi_cache_t100-200_n500.pt",
                       weights_only=True)
    t = cache["t"].numpy()
    xr = cache["xi_real"].numpy()
    xr_norm = xr / np.abs(xr).max()
    z = hilbert(xr_norm)

    with open(ANALYSIS / "phase_isolation_experiment.json") as f:
        iso = json.load(f)
    zeros = np.array(iso["zeros"])
    seeds = sorted(iso["results"].keys())

    with open(ANALYSIS / "angle_at_zeros_inventory.json") as f:
        inv = json.load(f)
    cstd_per_zero = np.array(inv["D1_seed_consistency"]["circular_std_per_zero"])
    xi_prime = np.array(inv["D2_xi_prime_sign"]["xi_prime_at_zero"])

    # ------------------------------------------------------------------
    # 1. 기본 정보: 영점 인덱스 및 이웃 간격
    # ------------------------------------------------------------------
    idx_focus = int(np.argmin(np.abs(zeros - T_FOCUS)))
    log(f"\n[1] 기본 정보")
    log(f"    target t={T_FOCUS}, 로드된 zero_idx={idx_focus}, "
        f"실제값={zeros[idx_focus]:.6f}")

    # 이웃 간격
    gaps_before = []
    gaps_after = []
    for k in [1, 2, 3]:
        if idx_focus - k >= 0:
            gaps_before.append(zeros[idx_focus - k + 1] - zeros[idx_focus - k])
        if idx_focus + k < len(zeros):
            gaps_after.append(zeros[idx_focus + k] - zeros[idx_focus + k - 1])
    log(f"    이전 3개 간격: {[f'{g:.4f}' for g in gaps_before[::-1]]}")
    log(f"    이후 3개 간격: {[f'{g:.4f}' for g in gaps_after]}")
    all_gaps = np.diff(zeros)
    log(f"    [100,200] 전체 간격 통계: "
        f"평균={all_gaps.mean():.4f}, 중앙={np.median(all_gaps):.4f}, "
        f"최소={all_gaps.min():.4f}, 최대={all_gaps.max():.4f}")

    # 173.41의 전/후 간격이 전체 분포에서 어디에 있나
    gap_prev = zeros[idx_focus] - zeros[idx_focus - 1]
    gap_next = zeros[idx_focus + 1] - zeros[idx_focus]
    pct_prev = float((all_gaps < gap_prev).mean() * 100)
    pct_next = float((all_gaps < gap_next).mean() * 100)
    log(f"\n    t=173.41 근접 간격:")
    log(f"      이전 gap = {gap_prev:.4f}  (전체에서 하위 {pct_prev:.0f}%)")
    log(f"      이후 gap = {gap_next:.4f}  (전체에서 하위 {pct_next:.0f}%)")

    # ------------------------------------------------------------------
    # 2. 국소 ξ 형상: t=[T_FOCUS-1.5, T_FOCUS+1.5] 확대 덤프
    # ------------------------------------------------------------------
    log(f"\n[2] 국소 ξ(½+it) 확대 (t ∈ [{T_FOCUS-WINDOW:.2f}, {T_FOCUS+WINDOW:.2f}])")
    mask = (t > T_FOCUS - WINDOW) & (t < T_FOCUS + WINDOW)
    t_loc = t[mask]
    xr_loc = xr[mask]
    arg_z_loc = np.angle(z[mask])

    # 영점 근방 샘플
    log(f"    샘플 {mask.sum()}개 (dt≈{(t_loc[1]-t_loc[0]):.4f})")
    log(f"    ξ_real 국소 극값: min={xr_loc.min():.3e}, max={xr_loc.max():.3e}")

    # 부호변화 위치 (이 창 안의 영점들)
    sc = np.where(np.diff(np.sign(xr_loc)) != 0)[0]
    log(f"    이 창의 영점 {len(sc)}개:")
    for i in sc:
        # 선형 보간으로 영점 위치 정밀화
        t0 = t_loc[i] - xr_loc[i] * (t_loc[i+1]-t_loc[i]) / (xr_loc[i+1]-xr_loc[i])
        log(f"        t={t0:.6f}  (샘플 {i} 와 {i+1} 사이)")

    # ------------------------------------------------------------------
    # 3. |ξ| 국소 거동: 영점들 사이의 peak 높이 비교
    # ------------------------------------------------------------------
    log(f"\n[3] |ξ_real| 국소 peak 분석")
    log(f"    인접 영점 쌍 사이의 peak 높이가 작으면 Lehmer-유사 — 두 영점이")
    log(f"    근접해서 그 사이에서 ξ가 충분히 오르지 못함.")

    # t=173.41 전후 영점 사이의 peak
    # zeros 배열에서 [idx_focus-1, idx_focus, idx_focus+1, idx_focus+2] 구간
    def peak_between(t_lo, t_hi):
        m = (t > t_lo) & (t < t_hi)
        if not m.any():
            return 0.0, None
        sub = np.abs(xr[m])
        t_sub = t[m]
        j = int(np.argmax(sub))
        return float(sub[j]), float(t_sub[j])

    for i in range(idx_focus - 2, idx_focus + 3):
        if i < 0 or i + 1 >= len(zeros):
            continue
        peak, tpeak = peak_between(zeros[i], zeros[i+1])
        mark = " <-- focus pair" if i == idx_focus or i == idx_focus - 1 else ""
        log(f"    [{zeros[i]:.4f}, {zeros[i+1]:.4f}] peak |ξ|={peak:.3e} "
            f"@t={tpeak:.3f}{mark}")

    # 전체 창의 모든 peak 와 비교
    peaks_all = []
    for i in range(len(zeros) - 1):
        p, _ = peak_between(zeros[i], zeros[i+1])
        peaks_all.append(p)
    peaks_all = np.array(peaks_all)
    # 스케일 정규화: ξ는 t 따라 빠르게 감쇠하므로 국소 평균으로 나눈 상대 peak
    log(f"\n    전체 [100,200]의 인접 영점 사이 peak 통계:")
    log(f"        mean={peaks_all.mean():.3e}, median={np.median(peaks_all):.3e}")
    log(f"        min={peaks_all.min():.3e}, max={peaks_all.max():.3e}")

    # t=173.41 전후 peak 의 상대 크기 (그 국소 이웃 평균 대비)
    def local_relative(i):
        # i번째 간격의 peak vs 주변 5개 간격 평균
        p = peaks_all[i]
        lo = max(0, i - 2)
        hi = min(len(peaks_all), i + 3)
        neighbors = np.delete(peaks_all[lo:hi], i - lo)
        return p / neighbors.mean() if len(neighbors) else float("nan")

    rel_prev = local_relative(idx_focus - 1)
    rel_next = local_relative(idx_focus)
    log(f"\n    t=173.41 전후 peak 상대크기 (주변 이웃 평균 대비):")
    log(f"      [172.50, 173.41] peak 상대 = {rel_prev:.3f}")
    log(f"      [173.41, 174.75] peak 상대 = {rel_next:.3f}")
    log(f"    (1.0보다 작으면 Lehmer-유사 압축)")

    # ------------------------------------------------------------------
    # 4. t=173.41 vs 174.75 대조 — 네트워크 학습 각
    # ------------------------------------------------------------------
    log(f"\n[4] 네트워크 arg Z_out — t=173.41 (어려운) vs t=174.75 (쉬운)")
    idx_easy = int(np.argmin(np.abs(zeros - T_CMP_EASY)))
    log(f"                          시드={seeds}")
    def angles_at(idx):
        return [iso["results"][s]["phi_z_circ"][idx] for s in seeds]
    a_hard = angles_at(idx_focus)
    a_easy = angles_at(idx_easy)
    log(f"    t=173.41 arg Z_out: {[f'{a:+.3f}' for a in a_hard]}  "
        f"cstd={cstd_per_zero[idx_focus]:.3f}")
    log(f"    t=174.75 arg Z_out: {[f'{a:+.3f}' for a in a_easy]}  "
        f"cstd={cstd_per_zero[idx_easy]:.3f}")
    log(f"    ξ'(173.41)={xi_prime[idx_focus]:+.3e}  sign={int(np.sign(xi_prime[idx_focus]))}")
    log(f"    ξ'(174.75)={xi_prime[idx_easy]:+.3e}  sign={int(np.sign(xi_prime[idx_easy]))}")

    # ------------------------------------------------------------------
    # 5. 해석신호 각 근방 거동
    # ------------------------------------------------------------------
    log(f"\n[5] arg z(t) 국소 거동 (±π/2 점프)")
    # t=173.41 근방 샘플
    mask_f = (t > T_FOCUS - 0.3) & (t < T_FOCUS + 0.3)
    t_f = t[mask_f]
    a_f = np.angle(z[mask_f])
    log(f"    t=173.41 ±0.3 구간 arg z (샘플 {mask_f.sum()}개):")
    for tt, aa in zip(t_f[::max(1,len(t_f)//8)], a_f[::max(1,len(t_f)//8)]):
        log(f"        t={tt:.4f}  arg z={aa:+.4f}")

    mask_e = (t > T_CMP_EASY - 0.3) & (t < T_CMP_EASY + 0.3)
    t_e = t[mask_e]
    a_e = np.angle(z[mask_e])
    log(f"\n    t=174.75 ±0.3 구간 arg z (샘플 {mask_e.sum()}개):")
    for tt, aa in zip(t_e[::max(1,len(t_e)//8)], a_e[::max(1,len(t_e)//8)]):
        log(f"        t={tt:.4f}  arg z={aa:+.4f}")

    # ------------------------------------------------------------------
    # 6. 두 영점 사이 다른 "어려운" 영점과의 교차
    # ------------------------------------------------------------------
    log(f"\n[6] result_t100-200 에서 hardest/easiest 교차")
    with open(OVERNIGHT / "result_t100-200.json") as f:
        r100 = json.load(f)
    log(f"    hardest_zeros: {r100.get('hardest_zeros')}")
    log(f"    easiest_zeros: {r100.get('easiest_zeros')}")
    # t=173.41 과 네트워크 cstd 최악 5 영점이 어디까지 hardest 와 겹치는지
    hardest_t = [x[0] for x in r100.get("hardest_zeros", []) if isinstance(x, list)]
    worst_net_idx = np.argsort(-cstd_per_zero)[:5]
    worst_net_t = zeros[worst_net_idx].tolist()
    log(f"    F₂ hardest t: {hardest_t}")
    log(f"    cstd worst t: {[round(x,4) for x in worst_net_t]}")
    overlap = set(round(x, 2) for x in hardest_t) & set(round(x, 2) for x in worst_net_t)
    log(f"    교집합 (2자리): {overlap}")

    # ------------------------------------------------------------------
    # 저장
    # ------------------------------------------------------------------
    out = {
        "target_t": T_FOCUS,
        "target_idx": idx_focus,
        "gap_prev": float(gap_prev),
        "gap_next": float(gap_next),
        "gap_percentile_prev": pct_prev,
        "gap_percentile_next": pct_next,
        "all_gaps_stats": {
            "mean": float(all_gaps.mean()),
            "median": float(np.median(all_gaps)),
            "min": float(all_gaps.min()),
            "max": float(all_gaps.max()),
        },
        "peaks_all": peaks_all.tolist(),
        "peak_before_target": float(peaks_all[idx_focus - 1]),
        "peak_after_target": float(peaks_all[idx_focus]),
        "peak_rel_before": float(rel_prev),
        "peak_rel_after": float(rel_next),
        "angles_seed_at_hard": a_hard,
        "angles_seed_at_easy": a_easy,
        "cstd_hard": float(cstd_per_zero[idx_focus]),
        "cstd_easy": float(cstd_per_zero[idx_easy]),
        "xi_prime_hard": float(xi_prime[idx_focus]),
        "xi_prime_easy": float(xi_prime[idx_easy]),
        "local_t": t_f.tolist(),
        "local_arg_z": a_f.tolist(),
        "local_xi_real": xr[mask_f].tolist(),
        "F2_hardest_r100": r100.get("hardest_zeros"),
        "F2_easiest_r100": r100.get("easiest_zeros"),
        "cstd_worst_5_t": worst_net_t,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")
    log(f"저장: {OUT_TXT}")


if __name__ == "__main__":
    main()
