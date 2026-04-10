#!/usr/bin/env python3
"""
=============================================================================
|F₂| 모델 독립성 검증 실험
=============================================================================
같은 데이터 구간에서 다른 크기의 모델을 훈련하고,
영점별 |F₂| 패턴이 보존되는지 검증한다.

- 보존 (Spearman ρ > 0.7) → |F₂|는 임계선의 내재적 성질 (수학적 발견)
- 비보존 (ρ < 0.3) → |F₂|는 모델 아티팩트

모델 크기: hidden = 32, 64, 128, 256 (각 3회 시드 반복)
구간: t ∈ [10000, 10100] (주), t ∈ [10100, 10200] (검증)

실행: ~/qrop_env/bin/python3 -u scripts/f2_model_independence.py
"""

import os
import sys
import json
import time
import copy
import itertools

import torch
from scipy.stats import spearmanr
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_ROOT)

from gdl.rdl.constants import R_CONST, PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.pipeline.xi_feature_dataset import (
    get_xi_feature_dataloaders, compute_zeros_in_range,
)

OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs", "analysis")
OVERNIGHT_DIR = os.path.join(PROJECT_ROOT, "outputs", "overnight")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_zeros_from_overnight(t_min, t_max):
    """기존 야간 탐사 JSON에서 영점 목록 로드 (mpmath 재계산 회피)."""
    tag = f"result_t{int(t_min)}-{int(t_max)}.json"
    path = os.path.join(OVERNIGHT_DIR, tag)
    if os.path.exists(path):
        with open(path) as f:
            data = json.load(f)
        zeros = data.get("zeros_list", [])
        if zeros:
            print(f"  [캐시] {path} 에서 {len(zeros)}개 영점 로드")
            return zeros
    # 캐시 없으면 직접 계산
    print(f"  [계산] mpmath로 영점 계산 (느림)...")
    return compute_zeros_in_range(t_min, t_max)

# =============================================================================
# 실험 설정
# =============================================================================
HIDDEN_SIZES = [32, 64, 128]
SEEDS = [42, 123, 777]
EPOCHS = 200  # 야간 탐사와 동일 — 충분한 수렴 보장
LR = 0.001
IN_FEATURES = 64
BATCH_SIZE = 32

# 실험 구간: 중위 50영점 (통계 검정력 확보)
RANGES = [
    (100, 200, 500, "중위 50 영점"),
]


# =============================================================================
# 단일 모델 훈련 + F₂ 추출
# =============================================================================
def train_and_extract_f2(train_loader, val_loader, dataset, zeros_list,
                         hidden, seed, epochs, device, label=""):
    """모델 훈련 후 영점별 |F₂|, cos², phase 추출."""

    torch.manual_seed(seed)
    np.random.seed(seed)

    # 모델 생성
    model = MasterResonantNetwork(
        in_features=IN_FEATURES,
        hidden_features=hidden,
        out_features=2,
        num_layers=3,
        channel_type="paper3ch",
        damping_mode="paper",
    )
    n_params = sum(p.numel() for p in model.parameters())

    # 손실 + 옵티마이저
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1,
        lambda_tgt=1.0, lambda_pqo=0.5,
        pqo_mode="cos2",
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)

    # 영점 특징 벡터
    zero_features = dataset.get_features_at_t(zeros_list).to(device)

    best_val = float("inf")
    best_state = None

    t0 = time.time()

    for epoch in range(1, epochs + 1):
        # 훈련
        model.train()
        train_loss_acc = 0.0
        n_batch = 0
        for X_batch, _ in train_loader:
            X_in = X_batch.to(device, dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            optimizer.zero_grad(set_to_none=True)
            outputs = model(X_in)
            total_loss, metrics = loss_fn(**outputs)
            if torch.isnan(total_loss) or torch.isinf(total_loss):
                continue
            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            train_loss_acc += total_loss.item()
            n_batch += 1

        # 검증
        model.eval()
        val_loss_acc = 0.0
        n_val = 0
        with torch.enable_grad():
            for X_batch, _ in val_loader:
                X_in = X_batch.to(device, dtype=PrecisionManager.REAL_DTYPE)
                X_in.requires_grad_(True)
                outputs = model(X_in)
                loss, _ = loss_fn(**outputs)
                if not (torch.isnan(loss) or torch.isinf(loss)):
                    val_loss_acc += loss.item()
                    n_val += 1
        val_loss = val_loss_acc / max(1, n_val)

        if val_loss < best_val:
            best_val = val_loss
            best_state = copy.deepcopy(model.state_dict())

        # 진행 출력 (50 에폭마다)
        if epoch % 50 == 0 or epoch == epochs:
            train_loss = train_loss_acc / max(1, n_batch)
            print(f"    {label} ep {epoch:03d}: train={train_loss:.4f} val={val_loss:.4f}")

    train_time = time.time() - t0

    # 최적 모델로 F₂ 추출
    if best_state:
        model.load_state_dict(best_state)

    model.eval()
    with torch.enable_grad():
        X_z = zero_features.clone().requires_grad_(True)
        out_z = model(X_z)

    Z_out = out_z["Z_out"]
    phi = out_z["phi"]
    psi = out_z["psi"]
    L_G = out_z["L_G"]

    # F₂ 계산
    phi_real = phi.to(dtype=PrecisionManager.REAL_DTYPE)
    rotation = torch.complex(torch.cos(phi_real), -torch.sin(phi_real))
    psi_complex = psi.to(dtype=PrecisionManager.COMPLEX_DTYPE)
    residual = rotation * (L_G - psi_complex)
    F2 = residual.imag.mean(dim=-1).detach().cpu()

    # PQO
    phase_z = torch.angle(Z_out).mean(dim=-1).detach().cpu()
    L = R_CONST.PQO_L_DEFAULT
    cos2_pqo = (torch.cos(L * phase_z / 2.0) ** 2)

    detected = sum(1 for f in F2.abs().tolist() if f < 0.1)

    return {
        "hidden": hidden,
        "seed": seed,
        "n_params": n_params,
        "train_time": round(train_time, 1),
        "best_val_loss": round(best_val, 6),
        "F2": F2.abs().tolist(),
        "cos2": cos2_pqo.tolist(),
        "phase": phase_z.tolist(),
        "F2_mean": round(float(F2.abs().mean()), 6),
        "F2_max": round(float(F2.abs().max()), 6),
        "detected": detected,
        "n_zeros": len(zeros_list),
    }


# =============================================================================
# 상관 분석
# =============================================================================
def analyze_correlations(results, zeros_list):
    """모델 간, 시드 간 Spearman 상관 분석."""

    lines = []
    all_keys = list(results.keys())  # (hidden, seed) 튜플

    # 1. 전체 쌍 상관 행렬
    lines.append("=" * 70)
    lines.append("  1. 전체 Spearman 상관 행렬 (|F₂| 패턴)")
    lines.append("=" * 70)

    n = len(all_keys)
    rho_matrix = np.zeros((n, n))
    p_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i == j:
                rho_matrix[i, j] = 1.0
                p_matrix[i, j] = 0.0
            else:
                rho, p = spearmanr(results[all_keys[i]]["F2"],
                                   results[all_keys[j]]["F2"])
                rho_matrix[i, j] = rho
                p_matrix[i, j] = p

    # 행렬 출력
    header = "          " + "  ".join(f"h{k[0]:>3}s{k[1]}" for k in all_keys)
    lines.append(header)
    for i, ki in enumerate(all_keys):
        row = f"h{ki[0]:>3}s{ki[1]}  "
        row += "  ".join(f"{rho_matrix[i,j]:>7.3f}" for j in range(n))
        lines.append(row)

    # 2. 모델 간 상관 (다른 hidden, 같은 시드 42)
    lines.append("")
    lines.append("=" * 70)
    lines.append("  2. 모델 크기 간 상관 (seed=42 고정)")
    lines.append("=" * 70)

    model_pairs = []
    for h1, h2 in itertools.combinations(HIDDEN_SIZES, 2):
        k1, k2 = (h1, 42), (h2, 42)
        if k1 in results and k2 in results:
            rho, p = spearmanr(results[k1]["F2"], results[k2]["F2"])
            lines.append(f"  hidden {h1:>3} vs {h2:>3}: ρ = {rho:.4f} (p = {p:.2e})")
            model_pairs.append((h1, h2, rho, p))

    # 3. 시드 간 상관 (같은 hidden, 다른 시드)
    lines.append("")
    lines.append("=" * 70)
    lines.append("  3. 시드 간 상관 (같은 모델, 다른 초기화)")
    lines.append("=" * 70)

    seed_rhos = []
    for h in HIDDEN_SIZES:
        for s1, s2 in itertools.combinations(SEEDS, 2):
            k1, k2 = (h, s1), (h, s2)
            if k1 in results and k2 in results:
                rho, p = spearmanr(results[k1]["F2"], results[k2]["F2"])
                lines.append(f"  hidden={h:>3}, seed {s1} vs {s2}: ρ = {rho:.4f} (p = {p:.2e})")
                seed_rhos.append(rho)

    # 4. 모델 크기별 평균 상관
    lines.append("")
    lines.append("=" * 70)
    lines.append("  4. 요약 통계")
    lines.append("=" * 70)

    # 모델 간 (시드 평균)
    cross_model_rhos = []
    for h1, h2 in itertools.combinations(HIDDEN_SIZES, 2):
        for s in SEEDS:
            k1, k2 = (h1, s), (h2, s)
            if k1 in results and k2 in results:
                rho, _ = spearmanr(results[k1]["F2"], results[k2]["F2"])
                cross_model_rhos.append(rho)

    if cross_model_rhos:
        mean_cross = np.mean(cross_model_rhos)
        std_cross = np.std(cross_model_rhos)
        lines.append(f"  모델 간 평균 ρ: {mean_cross:.4f} ± {std_cross:.4f}")
    if seed_rhos:
        mean_seed = np.mean(seed_rhos)
        std_seed = np.std(seed_rhos)
        lines.append(f"  시드 간 평균 ρ: {mean_seed:.4f} ± {std_seed:.4f}")

    # 5. 판정
    lines.append("")
    lines.append("=" * 70)
    lines.append("  5. 판정")
    lines.append("=" * 70)

    if cross_model_rhos:
        mean_rho = np.mean(cross_model_rhos)
        if mean_rho > 0.7:
            verdict = "강한 모델 독립성 — |F₂|는 임계선의 내재적 성질 (수학적 발견)"
        elif mean_rho > 0.3:
            verdict = "부분 독립 — 일부 구조는 내재적, 일부는 모델 의존"
        else:
            verdict = "모델 의존 — |F₂|는 모델 아티팩트"
        lines.append(f"  모델 간 평균 ρ = {mean_rho:.4f}")
        lines.append(f"  → {verdict}")

    # 6. Rank 보존 (top-10 겹침)
    lines.append("")
    lines.append("=" * 70)
    lines.append("  6. Top-10 어려운 영점 Rank 보존")
    lines.append("=" * 70)

    n_zeros = len(zeros_list)
    top_k = min(10, n_zeros // 2)

    for h in HIDDEN_SIZES:
        top_sets = []
        for s in SEEDS:
            k = (h, s)
            if k in results:
                f2_arr = results[k]["F2"]
                ranked = sorted(range(len(f2_arr)), key=lambda i: -f2_arr[i])
                top_sets.append(set(ranked[:top_k]))

        if len(top_sets) >= 2:
            # 시드 간 top-k 겹침
            overlap = top_sets[0]
            for ts in top_sets[1:]:
                overlap = overlap & ts
            lines.append(f"  hidden={h:>3}: 시드 간 top-{top_k} 겹침 = {len(overlap)}/{top_k}")

    # 모델 간 top-k 겹침 (seed=42)
    model_top_sets = []
    for h in HIDDEN_SIZES:
        k = (h, 42)
        if k in results:
            f2_arr = results[k]["F2"]
            ranked = sorted(range(len(f2_arr)), key=lambda i: -f2_arr[i])
            model_top_sets.append((h, set(ranked[:top_k])))

    if len(model_top_sets) >= 2:
        for i in range(len(model_top_sets)):
            for j in range(i+1, len(model_top_sets)):
                h1, s1 = model_top_sets[i]
                h2, s2 = model_top_sets[j]
                overlap = len(s1 & s2)
                lines.append(f"  hidden {h1} vs {h2} (seed=42): top-{top_k} 겹침 = {overlap}/{top_k}")

    # 7. |F₂| 크기 의존성
    lines.append("")
    lines.append("=" * 70)
    lines.append("  7. |F₂| 크기별 평균 (체계적 의존성 확인)")
    lines.append("=" * 70)

    for h in HIDDEN_SIZES:
        f2_means = []
        for s in SEEDS:
            k = (h, s)
            if k in results:
                f2_means.append(results[k]["F2_mean"])
        if f2_means:
            lines.append(f"  hidden={h:>3}: |F₂| mean = {np.mean(f2_means):.6f} ± {np.std(f2_means):.6f}")

    return "\n".join(lines), {
        "cross_model_rhos": cross_model_rhos,
        "seed_rhos": seed_rhos,
        "rho_matrix": rho_matrix.tolist(),
        "keys": [list(k) for k in all_keys],
    }


# =============================================================================
# 메인
# =============================================================================
def main():
    PrecisionManager.setup_precision()
    device = torch.device("cpu")

    all_range_results = {}

    for t_min, t_max, num_points, range_desc in RANGES:
        range_tag = f"t{int(t_min)}-{int(t_max)}"
        print(f"\n{'='*70}")
        print(f"  구간: t∈[{t_min},{t_max}] — {range_desc}")
        print(f"{'='*70}")

        # 1. 영점 로드 (기존 야간 탐사 캐시 우선)
        print(f"  영점 로드 중...")
        t0 = time.time()
        zeros_list = load_zeros_from_overnight(t_min, t_max)
        n_zeros = len(zeros_list)
        print(f"  {n_zeros}개 영점 ({time.time()-t0:.1f}초)")

        if n_zeros < 5:
            print(f"  SKIP: 영점 부족")
            continue

        # 2. 데이터 로더 (1회 — 모든 모델이 공유)
        #    캐시 파일명을 야간 탐사와 동일한 정수 형식으로 맞춤
        print(f"  데이터 준비 중...")
        cache_dir = os.path.join(PROJECT_ROOT, "outputs", "overnight")
        cache_name = f"xi_cache_t{int(t_min)}-{int(t_max)}_n{num_points}.pt"
        cache_path = os.path.join(cache_dir, cache_name)

        if os.path.exists(cache_path):
            print(f"  [캐시] {cache_path} 직접 로드")
            cache_data = torch.load(cache_path, weights_only=True)
        else:
            print(f"  [경고] 캐시 없음, mpmath 계산 (느림)")
            from gdl.rdl.pipeline.xi_feature_dataset import build_xi_cache
            cache_data = build_xi_cache(t_min, t_max, num_points, zeros_list=zeros_list)
            torch.save(cache_data, cache_path)

        from gdl.rdl.pipeline.xi_feature_dataset import XiFeatureDataset
        dataset = XiFeatureDataset(cache_data, in_features=IN_FEATURES)

        val_size = int(len(dataset) * 0.2)
        train_size = len(dataset) - val_size
        train_ds, val_ds = torch.utils.data.random_split(
            dataset, [train_size, val_size],
            generator=torch.Generator().manual_seed(0)  # 데이터 분할 고정
        )
        train_loader = torch.utils.data.DataLoader(
            train_ds, batch_size=BATCH_SIZE, shuffle=True, drop_last=True)
        val_loader = torch.utils.data.DataLoader(
            val_ds, batch_size=BATCH_SIZE, shuffle=False, drop_last=False)

        # 3. 모델별 훈련
        results = {}
        total_runs = len(HIDDEN_SIZES) * len(SEEDS)
        run_count = 0

        for hidden in HIDDEN_SIZES:
            for seed in SEEDS:
                run_count += 1
                label = f"[{run_count}/{total_runs}] h={hidden},s={seed}"
                print(f"\n  {label}")

                result = train_and_extract_f2(
                    train_loader, val_loader, dataset, zeros_list,
                    hidden, seed, EPOCHS, device, label
                )
                results[(hidden, seed)] = result

                det = result["detected"]
                print(f"    → |F₂|={result['F2_mean']:.4f}, "
                      f"탐지={det}/{n_zeros}, "
                      f"시간={result['train_time']}초")

        # 4. 상관 분석
        print(f"\n\n{'='*70}")
        print(f"  상관 분석: t∈[{t_min},{t_max}]")
        print(f"{'='*70}\n")

        analysis_text, analysis_data = analyze_correlations(results, zeros_list)
        print(analysis_text)

        all_range_results[range_tag] = {
            "zeros_list": zeros_list,
            "n_zeros": n_zeros,
            "results": {f"h{k[0]}_s{k[1]}": v for k, v in results.items()},
            "analysis": analysis_data,
        }

    # 5. 결과 저장
    # JSON
    json_path = os.path.join(OUTPUT_DIR, "f2_model_independence.json")
    with open(json_path, "w") as f:
        json.dump(all_range_results, f, indent=2, ensure_ascii=False)
    print(f"\n저장: {json_path}")

    # 텍스트 요약
    txt_path = os.path.join(OUTPUT_DIR, "f2_model_independence.txt")
    summary_lines = []
    summary_lines.append("=" * 70)
    summary_lines.append("  |F₂| 모델 독립성 검증 — 최종 결과")
    summary_lines.append(f"  실행: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    summary_lines.append("=" * 70)

    for range_tag, data in all_range_results.items():
        summary_lines.append(f"\n\n{'='*70}")
        summary_lines.append(f"  구간: {range_tag} ({data['n_zeros']}개 영점)")
        summary_lines.append("=" * 70)

        # 재분석 (results dict 복원)
        results_restored = {}
        for key_str, v in data["results"].items():
            # "h32_s42" → (32, 42)
            parts = key_str.split("_")
            h = int(parts[0][1:])
            s = int(parts[1][1:])
            results_restored[(h, s)] = v

        text, _ = analyze_correlations(results_restored, data["zeros_list"])
        summary_lines.append(text)

    with open(txt_path, "w") as f:
        f.write("\n".join(summary_lines))
    print(f"저장: {txt_path}")

    print(f"\n\n완료!")


if __name__ == "__main__":
    main()
