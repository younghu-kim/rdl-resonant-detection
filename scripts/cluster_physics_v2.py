#!/usr/bin/env python3
"""
cluster_physics_v2.py — N2e 재시도. 지수 감쇠 envelope 정규화.

이전 버전은 |ξ(½+it)| 의 exp(-πt/4) 감쇠 때문에 t=165 vs t=140 을 직접
비교하면 10 자리 차이가 났다. 여기서는:

  1. local envelope E(t) = rolling max |ξ_real| (창 폭 5)
  2. 정규화 ξ̂(t) = ξ_real(t) / E(t)  → 약 [-1, 1]
  3. 정규화된 영역에서 |ξ̂'|, |ξ̂''|, 국소 주파수 계산
  4. 군집 영역 vs 나머지 비교
"""
import json
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "outputs" / "analysis"
OVERNIGHT = ROOT / "outputs" / "overnight"
OUT_JSON = ANALYSIS / "cluster_physics_v2.json"
OUT_TXT = ANALYSIS / "cluster_physics_v2.txt"

CLUSTER_LO, CLUSTER_HI = 165.0, 173.5


def rolling_max(x, window):
    """센터된 rolling max."""
    out = np.empty_like(x)
    h = window // 2
    for i in range(len(x)):
        lo = max(0, i - h)
        hi = min(len(x), i + h + 1)
        out[i] = np.max(np.abs(x[lo:hi]))
    return out


def main():
    lines = []
    def log(s=""):
        print(s, flush=True); lines.append(s)

    log("=" * 72)
    log("  cluster_physics_v2 — envelope 정규화 군집 분석")
    log("=" * 72)

    cache = torch.load(OVERNIGHT / "xi_cache_t100-200_n500.pt",
                       weights_only=True)
    t = cache["t"].numpy()
    xr = cache["xi_real"].numpy()
    dt = float(t[1] - t[0])
    log(f"  샘플: N={len(t)}  dt={dt:.4f}")

    # 1. envelope
    # 창 폭: 인접 영점 평균 간격 ~ 2, 창 ±1 → 샘플 수
    window_t = 2.0
    window_samples = max(3, int(window_t / dt))
    if window_samples % 2 == 0:
        window_samples += 1
    E = rolling_max(xr, window_samples)
    E = np.where(E < 1e-300, 1e-300, E)
    xr_n = xr / E

    log(f"  envelope 창 = {window_samples} 샘플 ({window_t} in t)")
    log(f"  원본 |ξ_real|:   min={np.abs(xr).min():.2e}  max={np.abs(xr).max():.2e}")
    log(f"  정규화 |ξ̂|:     min={np.abs(xr_n).min():.2e}  max={np.abs(xr_n).max():.2e}")

    # 2. 미분
    xp_n = np.gradient(xr_n, t)
    xpp_n = np.gradient(xp_n, t)

    # 3. 마스크
    m_c = (t >= CLUSTER_LO) & (t <= CLUSTER_HI)
    m_o = ~m_c & (t > 110) & (t < 195)  # 경계 제외
    log(f"  군집 샘플: {m_c.sum()}   대조군 샘플: {m_o.sum()}")

    def cmp(name, arr_c, arr_o):
        mc = float(np.mean(arr_c))
        mo = float(np.mean(arr_o))
        mdc = float(np.median(arr_c))
        mdo = float(np.median(arr_o))
        ratio = mc / mo if mo != 0 else float('nan')
        log(f"    {name:>22s}  mean:{mc:>10.4f} vs {mo:>10.4f}  "
            f"median:{mdc:>8.4f} vs {mdo:>8.4f}  ratio:{ratio:>6.2f}×")
        return {"name": name, "cluster_mean": mc, "other_mean": mo,
                "cluster_median": mdc, "other_median": mdo, "ratio": ratio}

    stats = []
    log("\n  정규화 비교 (군집 vs 대조군):")
    stats.append(cmp("|ξ̂|", np.abs(xr_n[m_c]), np.abs(xr_n[m_o])))
    stats.append(cmp("|ξ̂'|", np.abs(xp_n[m_c]), np.abs(xp_n[m_o])))
    stats.append(cmp("|ξ̂''|", np.abs(xpp_n[m_c]), np.abs(xpp_n[m_o])))

    # 4. 국소 주파수 (영점 밀도 직접)
    sgn = np.sign(xr)
    # 군집 내 부호변화
    sc_c = (np.diff(sgn[m_c]) != 0).sum()
    sc_o = (np.diff(sgn[m_o]) != 0).sum()
    t_c_range = float(t[m_c].max() - t[m_c].min())
    t_o_range = float(t[m_o].max() - t[m_o].min())
    density_c = sc_c / t_c_range
    density_o = sc_o / t_o_range
    log(f"\n    영점 밀도 (부호변화/단위 t):")
    log(f"      군집: {density_c:.3f}   대조군: {density_o:.3f}   "
        f"비율: {density_c/density_o:.2f}×")

    # 5. 논리: 정규화 envelope 은 크기를 죽였으니 |ξ̂''| 이 크면 진짜 곡률 高
    log("\n  **해석**:")
    ratio_pp = stats[2]["ratio"]
    if ratio_pp > 1.3:
        log(f"    → |ξ̂''| 비율 {ratio_pp:.2f}× > 1.3: 군집은 고곡률 구역. "
            f"유한 대역폭 모델이 비싸게 표현해야 함.")
    elif ratio_pp < 0.77:
        log(f"    → |ξ̂''| 비율 {ratio_pp:.2f}× < 0.77: 군집은 저곡률. "
            f"곡률 가설 기각, 난이도의 다른 원인 필요.")
    else:
        log(f"    → |ξ̂''| 비율 {ratio_pp:.2f}×: 대조군과 유사. "
            f"곡률은 난이도와 무관.")

    # 6. 군집 내 개별 영점 근방 확대
    log("\n  군집 각 영점 ±0.5 창에서 정규화 곡률 mean:")
    cluster_zeros_expected = [165.5, 167.2, 169.1, 169.9, 173.4]
    cluster_zero_stats = []
    for z in cluster_zeros_expected:
        mm = (t > z - 0.5) & (t < z + 0.5)
        if mm.sum() > 0:
            xpp_mean = float(np.mean(np.abs(xpp_n[mm])))
            xp_mean = float(np.mean(np.abs(xp_n[mm])))
            cluster_zero_stats.append({
                "t": z, "local_abs_xpp_n": xpp_mean,
                "local_abs_xp_n": xp_mean,
            })
            log(f"    t={z:.1f}: |ξ̂''|={xpp_mean:.3f}  |ξ̂'|={xp_mean:.3f}")

    # 전체 밴드 평균 비교
    log(f"  전체 밴드 |ξ̂''| 평균: {float(np.mean(np.abs(xpp_n))):.3f}")
    log(f"  전체 밴드 |ξ̂'| 평균: {float(np.mean(np.abs(xp_n))):.3f}")

    out = {
        "cluster_range": [CLUSTER_LO, CLUSTER_HI],
        "window_samples": window_samples,
        "stats": stats,
        "zero_density_cluster": density_c,
        "zero_density_other": density_o,
        "cluster_zero_details": cluster_zero_stats,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")


if __name__ == "__main__":
    main()
