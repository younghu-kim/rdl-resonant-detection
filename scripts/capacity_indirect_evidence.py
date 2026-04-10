#!/usr/bin/env python3
"""
capacity_indirect_evidence.py — N2d + N2e.

표현 예산 가설의 간접 증거 수집 (학습 없음).

N2d: 파라미터 밀도 (n_params / n_zeros) vs 전역 mean|F₂|
  가설: params/zero 가 낮은 밴드가 전반적으로 어려움

N2e: 군집 t∈[165,173] 의 물리 특성
  - Hardy Z 진폭 (|xi_real| 국소 최대)
  - 인접 영점 간격 분포
  - ξ' 크기 분포 (접선 기울기)
  - 밴드의 다른 영점과 비교
"""
import json
import math
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "outputs" / "analysis"
OVERNIGHT = ROOT / "outputs" / "overnight"
OUT_JSON = ANALYSIS / "capacity_indirect_evidence.json"
OUT_TXT = ANALYSIS / "capacity_indirect_evidence.txt"


def main():
    lines = []
    def log(s=""):
        print(s, flush=True); lines.append(s)

    log("=" * 72)
    log("  capacity_indirect_evidence — N2d + N2e")
    log("=" * 72)

    # ================================================================
    # N2d — 파라미터 밀도 vs 난이도 (51 밴드)
    # ================================================================
    log("\n[N2d] 파라미터 밀도 vs 전역 난이도")
    log("-" * 72)

    result_files = sorted(OVERNIGHT.glob("result_t*.json"))
    rows = []
    for rf in result_files:
        with open(rf) as f:
            d = json.load(f)
        n_params = d.get("n_params")
        n_zeros = d.get("n_zeros")
        f2_mean = d.get("F2_mean")
        f2_max = d.get("F2_max")
        hidden = d.get("hidden")
        epochs = d.get("epochs")
        t_min = d.get("t_min")
        t_max = d.get("t_max")
        if None in (n_params, n_zeros, f2_mean) or n_zeros == 0:
            continue
        # 탐지율 파싱
        dr = d.get("detection_rate", "0/0")
        if isinstance(dr, str) and "/" in dr:
            a, b = dr.split("/")
            det_rate = float(a) / float(b) if float(b) > 0 else 0.0
        else:
            det_rate = 0.0
        rows.append({
            "band": rf.stem.replace("result_", ""),
            "t_min": t_min, "t_max": t_max,
            "n_params": n_params, "n_zeros": n_zeros,
            "params_per_zero": n_params / n_zeros,
            "F2_mean": f2_mean, "F2_max": f2_max,
            "hidden": hidden, "epochs": epochs,
            "detect_rate": det_rate,
        })

    log(f"    밴드 수: {len(rows)}")
    # 밴드 요약 (params/zero 오름차순 = 용량 부족 순)
    rows_sorted = sorted(rows, key=lambda r: r["params_per_zero"])
    log(f"\n    {'band':>16s} {'n_z':>5s} {'hid':>4s} {'p/z':>7s} "
        f"{'F2_mean':>10s} {'F2_max':>10s} {'det':>6s}")
    log("    " + "-" * 66)
    for r in rows_sorted[:10]:
        log(f"    {r['band']:>16s} {r['n_zeros']:>5d} {r['hidden']:>4d} "
            f"{r['params_per_zero']:>7.0f} {r['F2_mean']:>10.4f} "
            f"{r['F2_max']:>10.4f} {r['detect_rate']:>5.0%}")
    log("    ...")
    for r in rows_sorted[-5:]:
        log(f"    {r['band']:>16s} {r['n_zeros']:>5d} {r['hidden']:>4d} "
            f"{r['params_per_zero']:>7.0f} {r['F2_mean']:>10.4f} "
            f"{r['F2_max']:>10.4f} {r['detect_rate']:>5.0%}")

    # 상관: log(params/zero) vs F2_mean
    ppz = np.array([r["params_per_zero"] for r in rows])
    f2m = np.array([r["F2_mean"] for r in rows])
    f2x = np.array([r["F2_max"] for r in rows])
    det = np.array([r["detect_rate"] for r in rows])
    log_ppz = np.log(ppz)

    log("\n    상관 (전역, N={}):".format(len(rows)))
    log(f"      corr(log(p/z), F2_mean)  = {np.corrcoef(log_ppz, f2m)[0,1]:+.4f}")
    log(f"      corr(log(p/z), F2_max)   = {np.corrcoef(log_ppz, f2x)[0,1]:+.4f}")
    log(f"      corr(log(p/z), det_rate) = {np.corrcoef(log_ppz, det)[0,1]:+.4f}")
    log(f"    (용량 가설: log(p/z) 와 F₂ 는 음의 상관, det_rate 와는 양의 상관)")

    # 주의: overnight_exploration.py 는 밴드마다 hidden 을 바꿨으므로 단일 hidden 으로
    # 묶어서 비교해야 순수 "밴드 크기" 효과가 나옴
    log("\n    hidden 고정 후 밴드 비교:")
    from collections import defaultdict
    by_h = defaultdict(list)
    for r in rows:
        by_h[r["hidden"]].append(r)
    for h, group in sorted(by_h.items()):
        if len(group) < 3:
            continue
        ppzh = np.array([x["params_per_zero"] for x in group])
        f2mh = np.array([x["F2_mean"] for x in group])
        corr = float(np.corrcoef(np.log(ppzh), f2mh)[0, 1])
        log(f"      hidden={h:3d}  N={len(group):2d}  "
            f"corr(log(p/z), F2_mean) = {corr:+.4f}")

    n2d_result = {
        "rows": rows,
        "corr_log_ppz_F2_mean": float(np.corrcoef(log_ppz, f2m)[0, 1]),
        "corr_log_ppz_F2_max": float(np.corrcoef(log_ppz, f2x)[0, 1]),
        "corr_log_ppz_det_rate": float(np.corrcoef(log_ppz, det)[0, 1]),
        "by_hidden": {
            str(h): {
                "n": len(g),
                "corr": float(np.corrcoef(
                    np.log([x["params_per_zero"] for x in g]),
                    [x["F2_mean"] for x in g])[0, 1])
                if len(g) >= 3 else None,
            } for h, g in by_h.items()
        },
    }

    # ================================================================
    # N2e — 난이도 군집 t∈[165,173] 물리 특성
    # ================================================================
    log("\n\n[N2e] 난이도 군집 물리 특성 (t∈[165,173])")
    log("-" * 72)

    cache = torch.load(OVERNIGHT / "xi_cache_t100-200_n500.pt",
                       weights_only=True)
    t = cache["t"].numpy()
    xr = cache["xi_real"].numpy()

    # 영점 찾기 (부호변화 + 선형 보간)
    sgn = np.sign(xr)
    sc = np.where(np.diff(sgn) != 0)[0]
    zeros = []
    for i in sc:
        x0, x1 = xr[i], xr[i+1]
        if x1 == x0:
            zeros.append(float(t[i]))
        else:
            zeros.append(float(t[i] - x0 * (t[i+1] - t[i]) / (x1 - x0)))
    zeros = np.array(zeros)
    log(f"    밴드 내 영점 수 (부호변화): {len(zeros)}")

    # 인접 영점 간 |ξ| peak
    peaks = []
    for k in range(len(zeros) - 1):
        m = (t > zeros[k]) & (t < zeros[k+1])
        peaks.append(float(np.abs(xr[m]).max()) if m.any() else 0.0)
    peaks = np.array(peaks)
    gaps = np.diff(zeros)

    # ξ' 수치 미분
    xp = np.gradient(xr, t)
    xi_prime_at_zero = np.array([
        float(xp[int(np.argmin(np.abs(t - z)))]) for z in zeros
    ])

    # 국소 곡률 (2차 미분)
    xpp = np.gradient(xp, t)

    # 군집 영역 마스크
    cluster_mask = (zeros >= 165.0) & (zeros <= 173.5)
    other_mask = ~cluster_mask
    c_idx = np.where(cluster_mask)[0]
    o_idx = np.where(other_mask)[0]
    log(f"    군집 영점 수 (t∈[165,173.5]): {cluster_mask.sum()}")
    log(f"    그 외: {other_mask.sum()}")

    log("\n    특성 비교 (군집 vs 그 외):")
    log(f"    {'metric':>24s} {'군집 mean':>14s} {'그 외 mean':>14s} {'비율':>8s}")

    def cmp(name, arr_c, arr_o):
        mc, mo = float(np.mean(arr_c)), float(np.mean(arr_o))
        ratio = mc / mo if mo != 0 else float('nan')
        log(f"    {name:>24s} {mc:>14.4e} {mo:>14.4e} {ratio:>7.3f}×")
        return mc, mo, ratio

    # |ξ'| 군집 vs 나머지 (영점에서의 기울기 크기)
    stats = {}
    stats["abs_xi_prime"] = cmp("|ξ'(영점)|",
                                np.abs(xi_prime_at_zero[cluster_mask]),
                                np.abs(xi_prime_at_zero[other_mask]))

    # gap 주변 — 군집 포함하는 pair 추출
    c_pair_idx = [k for k in range(len(gaps))
                  if cluster_mask[k] or cluster_mask[k+1]]
    o_pair_idx = [k for k in range(len(gaps)) if k not in c_pair_idx]
    stats["gap"] = cmp("인접 간격",
                       gaps[c_pair_idx], gaps[o_pair_idx])
    stats["peak"] = cmp("인접 pair |ξ| peak",
                        peaks[c_pair_idx], peaks[o_pair_idx])

    # 국소 |ξ|, |ξ'|, |ξ''| — 군집 구간 내 t 샘플 전체
    m_t_c = (t >= 165.0) & (t <= 173.5)
    m_t_o = ~m_t_c
    stats["local_abs_xi"] = cmp("국소 |ξ|",
                                np.abs(xr[m_t_c]), np.abs(xr[m_t_o]))
    stats["local_abs_xi_prime"] = cmp("국소 |ξ'|",
                                       np.abs(xp[m_t_c]), np.abs(xp[m_t_o]))
    stats["local_abs_xi_pp"] = cmp("국소 |ξ''|",
                                    np.abs(xpp[m_t_c]), np.abs(xpp[m_t_o]))

    # 대역폭: 군집 구간 국소 주파수 = |xi'|/|xi| 중앙값 (영점 근처는 제외)
    eps = np.percentile(np.abs(xr), 5)
    mask_nontrivial = np.abs(xr) > eps
    if mask_nontrivial.any():
        local_freq = np.abs(xp) / (np.abs(xr) + 1e-30)
        lf_c = local_freq[m_t_c & mask_nontrivial]
        lf_o = local_freq[m_t_o & mask_nontrivial]
        stats["local_freq"] = cmp("국소 주파수 |ξ'/ξ|",
                                   lf_c, lf_o)

    log("\n    **해석**:")
    ratio_xi_pp = stats["local_abs_xi_pp"][2]
    ratio_freq = stats.get("local_freq", (0, 0, 1))[2]
    if ratio_xi_pp > 1.2:
        log(f"    → 군집 영역은 국소 곡률 |ξ''| 이 {ratio_xi_pp:.2f}× 더 큼. "
            f"표현 비용 증가 가설 지지.")
    elif ratio_xi_pp < 0.8:
        log(f"    → 군집은 오히려 곡률이 작음 ({ratio_xi_pp:.2f}×). 가설 기각.")
    else:
        log(f"    → 곡률 차이 미미 ({ratio_xi_pp:.2f}×). 중립.")

    if ratio_freq > 1.2:
        log(f"    → 국소 주파수 {ratio_freq:.2f}× — 군집은 고주파 구역, "
            f"유한 대역폭 모델이 어려워할 만함.")

    n2e_result = {
        "cluster_range": [165.0, 173.5],
        "n_cluster_zeros": int(cluster_mask.sum()),
        "n_other_zeros": int(other_mask.sum()),
        "comparisons": {k: {"cluster_mean": v[0], "other_mean": v[1],
                            "ratio": v[2]} for k, v in stats.items()},
    }

    # ================================================================
    # 저장 및 최종 요약
    # ================================================================
    log("\n" + "=" * 72)
    log("  종합")
    log("=" * 72)
    log(f"  N2d: log(params/zero) vs F₂_mean 전역 상관 = "
        f"{n2d_result['corr_log_ppz_F2_mean']:+.4f}")
    log(f"       (음수 & |r|>0.3 이면 용량 가설 지지)")
    log(f"  N2e: 군집 국소 |ξ''| 비율 = "
        f"{stats['local_abs_xi_pp'][2]:.3f}×")
    log(f"       (>1.2 이면 군집이 고곡률 구역 = 표현 비용 큼)")

    out = {"N2d": n2d_result, "N2e": n2e_result}
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")


if __name__ == "__main__":
    main()
