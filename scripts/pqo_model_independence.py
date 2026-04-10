#!/usr/bin/env python3
"""
=============================================================================
cos² PQO / 위상(phase) 모델 독립성 검증
=============================================================================
f2_model_independence.json에 저장된 9개 모델의 cos2/phase 배열을 재분석.
|F₂|가 모델 의존(ρ=0.0327)으로 판정되었으니, PQO가 독립인지 따로 본다.

추가 학습 없음 — 기존 결과 재사용.
"""

import os
import json
import itertools
import numpy as np
from scipy.stats import spearmanr

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
JSON_PATH = os.path.join(PROJECT_ROOT, "outputs", "analysis",
                         "f2_model_independence.json")
OUT_PATH = os.path.join(PROJECT_ROOT, "outputs", "analysis",
                        "pqo_model_independence.txt")

HIDDEN_SIZES = [32, 64, 128]
SEEDS = [42, 123, 777]


def analyze_field(results, field_name, label):
    """주어진 필드(cos2 또는 phase)에 대해 모델/시드 간 Spearman."""
    lines = []
    lines.append("=" * 70)
    lines.append(f"  {label} 모델 독립성 분석")
    lines.append("=" * 70)

    keys = list(results.keys())
    n = len(keys)

    # 1. 전체 상관 행렬
    rho_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                rho_matrix[i, j] = 1.0
            else:
                rho, _ = spearmanr(results[keys[i]][field_name],
                                   results[keys[j]][field_name])
                rho_matrix[i, j] = rho

    lines.append("")
    lines.append("  전체 Spearman 행렬:")
    header = "          " + "  ".join(f"h{k[0]:>3}s{k[1]}" for k in keys)
    lines.append(header)
    for i, ki in enumerate(keys):
        row = f"h{ki[0]:>3}s{ki[1]}  "
        row += "  ".join(f"{rho_matrix[i,j]:>7.3f}" for j in range(n))
        lines.append(row)

    # 2. 모델 간 상관 (시드 평균)
    lines.append("")
    lines.append("  모델 크기 간 상관 (모든 시드):")
    cross_model_rhos = []
    for h1, h2 in itertools.combinations(HIDDEN_SIZES, 2):
        for s in SEEDS:
            k1, k2 = (h1, s), (h2, s)
            if k1 in results and k2 in results:
                rho, p = spearmanr(results[k1][field_name],
                                   results[k2][field_name])
                cross_model_rhos.append(rho)
                lines.append(f"    h{h1} vs h{h2} (s={s}): ρ={rho:>7.4f} (p={p:.2e})")

    # 3. 시드 간 상관
    lines.append("")
    lines.append("  시드 간 상관 (같은 모델):")
    seed_rhos = []
    for h in HIDDEN_SIZES:
        for s1, s2 in itertools.combinations(SEEDS, 2):
            k1, k2 = (h, s1), (h, s2)
            if k1 in results and k2 in results:
                rho, p = spearmanr(results[k1][field_name],
                                   results[k2][field_name])
                seed_rhos.append(rho)
                lines.append(f"    h={h:>3}, s{s1} vs s{s2}: ρ={rho:>7.4f} (p={p:.2e})")

    # 4. 요약 + 판정
    lines.append("")
    lines.append("  요약:")
    if cross_model_rhos:
        mc = np.mean(cross_model_rhos)
        sc = np.std(cross_model_rhos)
        lines.append(f"    모델 간 평균 ρ: {mc:.4f} ± {sc:.4f}")
    if seed_rhos:
        ms = np.mean(seed_rhos)
        ss = np.std(seed_rhos)
        lines.append(f"    시드 간 평균 ρ: {ms:.4f} ± {ss:.4f}")

    lines.append("")
    lines.append("  판정:")
    if cross_model_rhos:
        mr = np.mean(cross_model_rhos)
        if mr > 0.7:
            verdict = f"강한 모델 독립성 — {label}는 임계선의 내재적 성질"
        elif mr > 0.3:
            verdict = f"부분 독립 — {label}는 일부 보존, 일부 모델 의존"
        else:
            verdict = f"모델 의존 — {label}도 모델 아티팩트"
        lines.append(f"    모델 간 평균 ρ = {mr:.4f}")
        lines.append(f"    → {verdict}")

    # 5. 분포 통계 (필드 평균/표준편차)
    lines.append("")
    lines.append(f"  {label} 분포 통계 (모델 크기별):")
    for h in HIDDEN_SIZES:
        means = []
        stds = []
        for s in SEEDS:
            k = (h, s)
            if k in results:
                arr = np.array(results[k][field_name])
                means.append(arr.mean())
                stds.append(arr.std())
        if means:
            lines.append(f"    h={h:>3}: mean={np.mean(means):.4f} ± {np.std(means):.4f}, "
                         f"std={np.mean(stds):.4f}")

    return "\n".join(lines), cross_model_rhos, seed_rhos


def main():
    with open(JSON_PATH) as f:
        data = json.load(f)

    range_data = data["t100-200"]
    n_zeros = range_data["n_zeros"]
    raw_results = range_data["results"]

    # "h32_s42" → (32, 42) 키 변환
    results = {}
    for key_str, v in raw_results.items():
        parts = key_str.split("_")
        h = int(parts[0][1:])
        s = int(parts[1][1:])
        results[(h, s)] = v

    print(f"로드: {n_zeros}개 영점, {len(results)}개 모델 결과")
    print()

    all_lines = []
    all_lines.append("=" * 70)
    all_lines.append(f"  cos² PQO / phase 모델 독립성 — t∈[100,200], {n_zeros}개 영점")
    all_lines.append("=" * 70)

    # cos² PQO 분석
    text_cos2, cm_cos2, sm_cos2 = analyze_field(results, "cos2", "cos² PQO")
    all_lines.append("")
    all_lines.append(text_cos2)

    all_lines.append("")
    all_lines.append("")

    # phase 분석
    text_phase, cm_phase, sm_phase = analyze_field(results, "phase", "위상 (phase)")
    all_lines.append(text_phase)

    # 비교 요약
    all_lines.append("")
    all_lines.append("=" * 70)
    all_lines.append("  최종 비교: |F₂| vs cos² PQO vs phase")
    all_lines.append("=" * 70)
    all_lines.append(f"  |F₂|     모델 간 ρ = 0.0327 (이전 분석)")
    if cm_cos2:
        all_lines.append(f"  cos² PQO 모델 간 ρ = {np.mean(cm_cos2):.4f}")
    if cm_phase:
        all_lines.append(f"  phase    모델 간 ρ = {np.mean(cm_phase):.4f}")

    output = "\n".join(all_lines)
    print(output)

    with open(OUT_PATH, "w") as f:
        f.write(output)
    print(f"\n저장: {OUT_PATH}")


if __name__ == "__main__":
    main()
