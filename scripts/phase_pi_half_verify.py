#!/usr/bin/env python3
"""
=============================================================================
|φ - π/2| 정밀 검증
=============================================================================
phase가 영점에서 π/2로 수렴한다는 가설 검증.
기존 f2_model_independence.json 재분석 (학습 0).
"""

import os
import json
import math
import itertools
import numpy as np
from scipy.stats import spearmanr

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
JSON_PATH = os.path.join(PROJECT_ROOT, "outputs", "analysis",
                         "f2_model_independence.json")
OUT_PATH = os.path.join(PROJECT_ROOT, "outputs", "analysis",
                        "phase_pi_half_verify.txt")

PI_HALF = math.pi / 2.0
HIDDEN_SIZES = [32, 64, 128]
SEEDS = [42, 123, 777]


def main():
    with open(JSON_PATH) as f:
        data = json.load(f)

    range_data = data["t100-200"]
    n_zeros = range_data["n_zeros"]
    raw = range_data["results"]

    # 키 변환 + dev 배열 계산 (φ - π/2 와 |φ - π/2|)
    results = {}
    for key_str, v in raw.items():
        parts = key_str.split("_")
        h = int(parts[0][1:])
        s = int(parts[1][1:])
        phi = np.array(v["phase"])
        results[(h, s)] = {
            "phi": phi,
            "dev": phi - PI_HALF,
            "abs_dev": np.abs(phi - PI_HALF),
        }

    lines = []
    lines.append("=" * 70)
    lines.append(f"  |φ - π/2| 정밀 검증 — t∈[100,200], {n_zeros}개 영점")
    lines.append(f"  π/2 = {PI_HALF:.10f}")
    lines.append("=" * 70)

    # ========================================================================
    # 1. 모델별 통계
    # ========================================================================
    lines.append("")
    lines.append("=" * 70)
    lines.append("  1. 모델별 |φ - π/2| 통계 (영점 50개 평균)")
    lines.append("=" * 70)
    lines.append(f"  {'모델':12s} {'φ 평균':>10s} {'φ-π/2 평균':>14s} {'|φ-π/2| 평균':>14s} {'|φ-π/2| 표준편차':>16s} {'|φ-π/2| 최대':>14s}")

    summary_by_h = {h: [] for h in HIDDEN_SIZES}
    for h in HIDDEN_SIZES:
        for s in SEEDS:
            r = results[(h, s)]
            phi_mean = r["phi"].mean()
            dev_mean = r["dev"].mean()
            abs_mean = r["abs_dev"].mean()
            abs_std = r["abs_dev"].std()
            abs_max = r["abs_dev"].max()
            label = f"h={h:>3},s={s}"
            lines.append(f"  {label:12s} {phi_mean:>10.6f} {dev_mean:>+14.6f} {abs_mean:>14.6f} {abs_std:>16.6f} {abs_max:>14.6f}")
            summary_by_h[h].append(abs_mean)

    # ========================================================================
    # 2. 모델 크기별 수렴 (|φ-π/2| 평균)
    # ========================================================================
    lines.append("")
    lines.append("=" * 70)
    lines.append("  2. 모델 크기별 |φ - π/2| 평균 (시드 평균/표준편차)")
    lines.append("=" * 70)
    lines.append(f"  {'hidden':>8s} {'평균':>14s} {'표준편차':>14s}")
    for h in HIDDEN_SIZES:
        arr = np.array(summary_by_h[h])
        lines.append(f"  {h:>8d} {arr.mean():>14.6f} {arr.std():>14.6f}")

    lines.append("")
    lines.append("  → 모델 클수록 |φ-π/2| 평균이 작아지면 = 진짜 수렴 (학습 한계)")
    lines.append("  → 변하지 않으면 = 시스템 노이즈 (이론적 잔차)")

    # ========================================================================
    # 3. 영점별 |φ-π/2| 의 모델 간 Spearman
    # ========================================================================
    lines.append("")
    lines.append("=" * 70)
    lines.append("  3. |φ - π/2| 영점별 패턴의 모델 간 Spearman ρ")
    lines.append("=" * 70)

    keys = [(h, s) for h in HIDDEN_SIZES for s in SEEDS]
    n = len(keys)

    # 행렬
    rho_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                rho_mat[i, j] = 1.0
            else:
                rho, _ = spearmanr(results[keys[i]]["abs_dev"],
                                   results[keys[j]]["abs_dev"])
                rho_mat[i, j] = rho

    header = "          " + "  ".join(f"h{k[0]:>3}s{k[1]}" for k in keys)
    lines.append(header)
    for i, ki in enumerate(keys):
        row = f"h{ki[0]:>3}s{ki[1]}  "
        row += "  ".join(f"{rho_mat[i,j]:>7.3f}" for j in range(n))
        lines.append(row)

    cross_model_rhos = []
    for h1, h2 in itertools.combinations(HIDDEN_SIZES, 2):
        for s in SEEDS:
            rho, p = spearmanr(results[(h1, s)]["abs_dev"],
                               results[(h2, s)]["abs_dev"])
            cross_model_rhos.append(rho)

    seed_rhos = []
    for h in HIDDEN_SIZES:
        for s1, s2 in itertools.combinations(SEEDS, 2):
            rho, _ = spearmanr(results[(h, s1)]["abs_dev"],
                               results[(h, s2)]["abs_dev"])
            seed_rhos.append(rho)

    lines.append("")
    lines.append("  요약:")
    lines.append(f"    모델 간 평균 ρ: {np.mean(cross_model_rhos):.4f} ± {np.std(cross_model_rhos):.4f}")
    lines.append(f"    시드 간 평균 ρ: {np.mean(seed_rhos):.4f} ± {np.std(seed_rhos):.4f}")

    # ========================================================================
    # 4. φ - π/2 의 부호 분포 (대칭성)
    # ========================================================================
    lines.append("")
    lines.append("=" * 70)
    lines.append("  4. φ - π/2 부호 분포 (대칭성 검사)")
    lines.append("=" * 70)
    for h in HIDDEN_SIZES:
        for s in SEEDS:
            dev = results[(h, s)]["dev"]
            n_pos = int((dev > 0).sum())
            n_neg = int((dev < 0).sum())
            mean_dev = dev.mean()
            lines.append(f"  h={h:>3},s={s}: +/- = {n_pos}/{n_neg}, 평균편차 = {mean_dev:+.6f}")

    # ========================================================================
    # 5. |F₂| vs |φ-π/2| 상관 (영점 어려움 = 위상 미수렴?)
    # ========================================================================
    lines.append("")
    lines.append("=" * 70)
    lines.append("  5. |F₂| vs |φ-π/2| 영점별 상관 (모델 내)")
    lines.append("=" * 70)
    f2_phi_rhos = []
    for h in HIDDEN_SIZES:
        for s in SEEDS:
            f2 = np.array(raw[f"h{h}_s{s}"]["F2"])
            abs_dev = results[(h, s)]["abs_dev"]
            rho, p = spearmanr(f2, abs_dev)
            f2_phi_rhos.append(rho)
            lines.append(f"  h={h:>3},s={s}: ρ(|F₂|, |φ-π/2|) = {rho:+.4f} (p={p:.2e})")
    lines.append(f"  → 평균 ρ = {np.mean(f2_phi_rhos):+.4f}")

    # ========================================================================
    # 6. 최종 판정
    # ========================================================================
    lines.append("")
    lines.append("=" * 70)
    lines.append("  6. 최종 판정")
    lines.append("=" * 70)

    all_abs_dev = np.concatenate([results[k]["abs_dev"] for k in keys])
    global_mean = all_abs_dev.mean()
    global_max = all_abs_dev.max()
    global_p99 = np.percentile(all_abs_dev, 99)

    lines.append(f"  전체 |φ - π/2|:")
    lines.append(f"    평균 = {global_mean:.6f}")
    lines.append(f"    99% 분위 = {global_p99:.6f}")
    lines.append(f"    최대 = {global_max:.6f}")
    lines.append("")

    h_means = {h: np.mean(summary_by_h[h]) for h in HIDDEN_SIZES}
    decreasing = h_means[32] > h_means[64] > h_means[128]
    lines.append(f"  모델 크기 따른 수렴: h32={h_means[32]:.4f} → h64={h_means[64]:.4f} → h128={h_means[128]:.4f}")
    lines.append(f"  단조 감소: {'예 (수렴)' if decreasing else '아니오 (시스템 노이즈)'}")
    lines.append("")

    cm_mean = np.mean(cross_model_rhos)
    lines.append(f"  |φ-π/2| 패턴 모델 간 ρ = {cm_mean:.4f}")
    if cm_mean > 0.5:
        verdict = "→ 진짜 영점별 구조 존재"
    elif cm_mean > 0.2:
        verdict = "→ 약한 보존, 부분 신호"
    else:
        verdict = "→ |φ-π/2| 자체는 노이즈, 평균값만 의미 있음"
    lines.append(f"  {verdict}")

    output = "\n".join(lines)
    print(output)
    with open(OUT_PATH, "w") as f:
        f.write(output)
    print(f"\n저장: {OUT_PATH}")


if __name__ == "__main__":
    main()
