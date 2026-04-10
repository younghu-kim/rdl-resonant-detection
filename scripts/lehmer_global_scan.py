#!/usr/bin/env python3
"""
lehmer_global_scan.py — 전 xi_cache 대역에 Lehmer 압축 지표 적용 (N1c).

N1/N1a 결과: peak_rel = (인접 영점 쌍 사이 |ξ| peak) / (주변 이웃 peak 평균)
이 ≲ 0.1 인 쌍은 Lehmer-유사. t=173.41(0.080) 과 t=107.17(0.092) 가 양쪽
모두 F₂ 탐지기가 어려워한 영점이며, 두 독립 지표의 교차 최악점이었다.

이 스크립트는 52개 xi_cache 파일 전체에 이 지표를 적용하여:
  1. 각 대역의 "Lehmer 후보 쌍" 추출 (peak_rel < τ, τ ∈ {0.1, 0.2, 0.3})
  2. 각 대역의 result_t*.json 에서 hardest_zeros 목록 로드
  3. Lehmer-flagged 영점과 F₂-hardest 영점의 교차율 통계
  4. 전역 집계: 14,441개 영점 중 Lehmer 몇 %, F₂-hardest 몇 %, 교차 몇 %
  5. 교차 영점에서 어느 쪽 멤버(공백쪽/쌍-안쪽)인지 N1a 비대칭 가설 검증

출력: outputs/analysis/lehmer_global_scan.{json,txt}
"""
import json
import re
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "outputs" / "analysis"
OVERNIGHT = ROOT / "outputs" / "overnight"
OUT_JSON = ANALYSIS / "lehmer_global_scan.json"
OUT_TXT = ANALYSIS / "lehmer_global_scan.txt"

TAUS = [0.10, 0.20, 0.30]  # peak_rel 임계


def find_zeros_and_peaks(t, xr):
    """sign change 위치를 선형 보간으로 정밀화. 인접 쌍 사이 |ξ| peak 반환."""
    sgn = np.sign(xr)
    sc = np.where(np.diff(sgn) != 0)[0]
    zeros = []
    for i in sc:
        # 선형 보간
        x0, x1 = xr[i], xr[i+1]
        if x1 == x0:
            zeros.append(float(t[i]))
        else:
            z = t[i] - x0 * (t[i+1] - t[i]) / (x1 - x0)
            zeros.append(float(z))
    zeros = np.array(zeros)
    if len(zeros) < 2:
        return zeros, np.empty(0)

    peaks = np.empty(len(zeros) - 1)
    for k in range(len(zeros) - 1):
        m = (t > zeros[k]) & (t < zeros[k+1])
        if not m.any():
            peaks[k] = 0.0
        else:
            peaks[k] = float(np.abs(xr[m]).max())
    return zeros, peaks


def peak_rel(peaks, k, radius=2):
    """k 번째 peak 를 좌우 radius 개 이웃 peak 평균으로 정규화."""
    lo = max(0, k - radius)
    hi = min(len(peaks), k + radius + 1)
    nbr_idx = [j for j in range(lo, hi) if j != k]
    if not nbr_idx:
        return float('nan')
    nbr_mean = float(np.mean([peaks[j] for j in nbr_idx]))
    if nbr_mean <= 0:
        return float('nan')
    return peaks[k] / nbr_mean


def band_key_from_cache(name):
    """ 'xi_cache_t100-200_n500.pt' -> 't100-200' """
    m = re.match(r"xi_cache_(t[\d\-]+)_n\d+\.pt", name)
    return m.group(1) if m else None


def main():
    lines = []
    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    log("=" * 72)
    log("  lehmer_global_scan — 전 대역 Lehmer 압축 통계 (N1c)")
    log("=" * 72)

    cache_files = sorted(OVERNIGHT.glob("xi_cache_*.pt"))
    result_files = {rf.stem.replace("result_", ""): rf
                    for rf in OVERNIGHT.glob("result_t*.json")}
    log(f"  xi_cache: {len(cache_files)} 개")
    log(f"  result  : {len(result_files)} 개")

    per_band = []
    for cf in cache_files:
        key = band_key_from_cache(cf.name)
        if key is None:
            continue
        cache = torch.load(cf, weights_only=True)
        t = cache["t"].numpy()
        xr = cache["xi_real"].numpy()

        zeros, peaks = find_zeros_and_peaks(t, xr)
        if len(peaks) < 3:
            continue

        # peak_rel per pair
        rel = np.array([peak_rel(peaks, k) for k in range(len(peaks))])

        # Lehmer flagged zeros: 쌍 (k, k+1) 의 한 쪽 영점을 flagged 처리
        # 여기서는 "pair-level" flagging: rel[k] < τ 이면 인덱스 k 와 k+1 두 영점
        # 모두 Lehmer-멤버.
        flagged = {tau: set() for tau in TAUS}
        for k, r in enumerate(rel):
            if np.isnan(r):
                continue
            for tau in TAUS:
                if r < tau:
                    flagged[tau].add(k)      # 왼쪽 영점
                    flagged[tau].add(k+1)    # 오른쪽 영점

        # F₂ hardest 로드 (해당 대역)
        rf = result_files.get(key)
        hardest_t = []
        n_zeros_in_band = None
        if rf:
            with open(rf) as f:
                rd = json.load(f)
            hz = rd.get("hardest_zeros", [])
            if isinstance(hz, list):
                hardest_t = [float(item[0]) for item in hz
                             if isinstance(item, list) and len(item) >= 1]
            try:
                n_zeros_in_band = int(rd.get("n_zeros"))
            except (TypeError, ValueError):
                n_zeros_in_band = None

        # F₂-hardest 영점의 인덱스 (가장 가까운 zero)
        hardest_idx = []
        for ht in hardest_t:
            if len(zeros) == 0:
                continue
            j = int(np.argmin(np.abs(zeros - ht)))
            if abs(zeros[j] - ht) < 0.01:  # 1/100 매칭
                hardest_idx.append(j)

        # 교차: hardest 가 Lehmer-flagged 인가?
        band_stats = {
            "key": key,
            "t_range": [float(t.min()), float(t.max())],
            "n_zeros": len(zeros),
            "n_pairs": len(peaks),
            "n_hardest_F2": len(hardest_idx),
            "peak_rel_min": float(np.nanmin(rel)) if len(rel) else None,
            "peak_rel_median": float(np.nanmedian(rel)) if len(rel) else None,
        }
        for tau in TAUS:
            n_flag = len(flagged[tau])
            n_cross = len([j for j in hardest_idx if j in flagged[tau]])
            band_stats[f"lehmer_{tau}_n"] = n_flag
            band_stats[f"lehmer_{tau}_frac"] = n_flag / max(1, len(zeros))
            band_stats[f"cross_{tau}_hardest_in_lehmer"] = n_cross
            band_stats[f"cross_{tau}_frac"] = (n_cross / max(1, len(hardest_idx))
                                                if hardest_idx else None)

        # N1a 비대칭 가설: flagged 쌍에서 "공백쪽" 멤버가 더 F₂-hard 인가?
        # 쌍 k (영점 k, k+1) 이 Lehmer 압축이라 할 때, 왼쪽 영점 k 의 왼쪽 gap
        # 이 오른쪽 영점 k+1 의 오른쪽 gap 과 비교. 큰 gap 쪽 멤버가 어떤 쪽인지.
        gaps = np.diff(zeros)
        asym_stats = {"lehmer_pairs": 0, "hardest_on_gap_side": 0,
                      "hardest_on_tight_side": 0, "hardest_ambiguous": 0}
        for k, r in enumerate(rel):
            if np.isnan(r) or r >= 0.2:
                continue
            # k 번째 pair = 영점 k 와 k+1. gap_prev(k) 와 gap_next(k+1)
            gap_left = gaps[k-1] if k-1 >= 0 else None
            gap_right = gaps[k+1] if k+1 < len(gaps) else None
            if gap_left is None or gap_right is None:
                continue
            asym_stats["lehmer_pairs"] += 1
            # 공백 쪽 = 더 큰 gap 쪽 영점
            if gap_left > gap_right:
                gap_side = k       # 왼쪽 영점이 공백 맞닿음
                tight_side = k + 1
            else:
                gap_side = k + 1
                tight_side = k
            # hardest_idx 에 어느 쪽이 들었나
            in_gap = gap_side in hardest_idx
            in_tight = tight_side in hardest_idx
            if in_gap and not in_tight:
                asym_stats["hardest_on_gap_side"] += 1
            elif in_tight and not in_gap:
                asym_stats["hardest_on_tight_side"] += 1
            elif in_gap and in_tight:
                asym_stats["hardest_ambiguous"] += 1
            # 아무쪽도 hardest 가 아니면 집계 안 함 (sample size 감안)
        band_stats["asymmetry"] = asym_stats

        # 교차 영점 상세 (tau=0.2)
        cross_details = []
        for j in hardest_idx:
            if j in flagged[0.20]:
                cross_details.append({
                    "t": float(zeros[j]),
                    "peak_rel_left": float(rel[j-1]) if j-1 >= 0 else None,
                    "peak_rel_right": float(rel[j]) if j < len(rel) else None,
                })
        band_stats["cross_details_tau0.2"] = cross_details

        per_band.append(band_stats)

    # ------------------------------------------------------------------
    # 전역 집계
    # ------------------------------------------------------------------
    log("\n" + "=" * 72)
    log("  대역별 요약 (Lehmer τ=0.2 기준)")
    log("=" * 72)
    log(f"  {'key':>16s} {'nZ':>5s} {'L_n':>4s} {'L%':>5s} "
        f"{'hrd':>4s} {'cross':>6s} {'cross%':>7s}")
    for b in per_band[:10]:
        c = b['cross_0.2_frac']
        cs = f"{c*100:>6.1f}%" if c is not None else "   --"
        log(f"  {b['key']:>16s} {b['n_zeros']:>5d} "
            f"{b['lehmer_0.2_n']:>4d} {b['lehmer_0.2_frac']*100:>4.1f}% "
            f"{b['n_hardest_F2']:>4d} "
            f"{b['cross_0.2_hardest_in_lehmer']:>6d} {cs}")
    log("    ...")
    for b in per_band[-10:]:
        c = b['cross_0.2_frac']
        cs = f"{c*100:>6.1f}%" if c is not None else "   --"
        log(f"  {b['key']:>16s} {b['n_zeros']:>5d} "
            f"{b['lehmer_0.2_n']:>4d} {b['lehmer_0.2_frac']*100:>4.1f}% "
            f"{b['n_hardest_F2']:>4d} "
            f"{b['cross_0.2_hardest_in_lehmer']:>6d} {cs}")

    # 전역 숫자
    log("\n" + "=" * 72)
    log("  전역 집계")
    log("=" * 72)
    total_zeros = sum(b["n_zeros"] for b in per_band)
    total_hardest = sum(b["n_hardest_F2"] for b in per_band)
    log(f"  총 영점: {total_zeros}")
    log(f"  총 F₂-hardest (밴드별 상위 5): {total_hardest}")

    for tau in TAUS:
        total_flag = sum(b[f"lehmer_{tau}_n"] for b in per_band)
        total_cross = sum(b[f"cross_{tau}_hardest_in_lehmer"] for b in per_band)
        log(f"\n  τ = {tau}")
        log(f"    Lehmer flagged 영점: {total_flag} "
            f"({total_flag/max(1,total_zeros)*100:.2f}%)")
        log(f"    F₂-hardest 중 Lehmer-flagged: {total_cross}/{total_hardest} "
            f"= {total_cross/max(1,total_hardest)*100:.1f}%")

        # 랜덤 기대값과 비교
        p = total_flag / max(1, total_zeros)
        expected = total_hardest * p
        observed = total_cross
        if expected > 0:
            enrichment = observed / expected
            log(f"    랜덤 기대값: {expected:.2f}  "
                f"→ 관찰/기대 = {enrichment:.2f}×")
            # 이항 검정 (대략적 z 점수)
            if p > 0 and p < 1:
                var = total_hardest * p * (1 - p)
                z = (observed - expected) / max(1e-9, np.sqrt(var))
                log(f"    이항 z 점수 ≈ {z:.2f}  "
                    f"({'강한 농축' if z > 3 else '농축' if z > 2 else '약함'})")

    # 비대칭 집계
    log("\n" + "=" * 72)
    log("  N1a 비대칭 가설 검증 — '공백쪽 멤버가 어려운가'")
    log("=" * 72)
    total_pairs = sum(b["asymmetry"]["lehmer_pairs"] for b in per_band)
    gap_hard = sum(b["asymmetry"]["hardest_on_gap_side"] for b in per_band)
    tight_hard = sum(b["asymmetry"]["hardest_on_tight_side"] for b in per_band)
    amb = sum(b["asymmetry"]["hardest_ambiguous"] for b in per_band)
    log(f"  검사된 Lehmer 쌍 (τ=0.2, 좌우 gap 정의 가능): {total_pairs}")
    log(f"  hardest 가 공백쪽에만: {gap_hard}")
    log(f"  hardest 가 타이트쪽에만: {tight_hard}")
    log(f"  hardest 가 양쪽: {amb}")
    total_single = gap_hard + tight_hard
    if total_single > 0:
        p_gap = gap_hard / total_single
        log(f"  단측 사례 중 공백쪽 비율: {p_gap*100:.1f}% "
            f"(N1a 가설 예측: >50%)")

    # 저장
    out = {
        "taus": TAUS,
        "per_band": per_band,
        "global": {
            "total_zeros": total_zeros,
            "total_hardest": total_hardest,
            "asymmetry": {
                "total_pairs": total_pairs,
                "gap_side_only": gap_hard,
                "tight_side_only": tight_hard,
                "both_sides": amb,
            },
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")
    log(f"저장: {OUT_TXT}")


if __name__ == "__main__":
    main()
