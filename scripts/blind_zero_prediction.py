#!/usr/bin/env python3
"""
blind_zero_prediction.py — 블라인드 영점 예측 실험.

핵심 질문:
  학습 범위 밖의 영점을 |F₂| 극소값으로 예측할 수 있는가?

설계:
  (A) 순방향: t∈[100,150] 학습 → t∈[150,200] 예측
  (B) 역방향: t∈[150,200] 학습 → t∈[100,150] 예측

평가:
  - 밀집 격자(2000점) 위 |F₂| 극소값 탐색
  - 실제 영점과 최근접 극소값 거리 < 임계값(0.3) → 적중(hit)
  - 정밀도(precision), 재현율(recall), F1 스코어

구성:
  in_features ∈ {64, 128}, seeds = [42, 7, 123, 2026, 314], ep=200
"""
import copy
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "5")
os.environ.setdefault("MKL_NUM_THREADS", "5")
import numpy as np
import torch
torch.set_num_threads(5)

from scipy.signal import argrelmin

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.pipeline.xi_feature_dataset import get_xi_feature_dataloaders

OUTPUT_DIR = PROJECT_ROOT / "outputs" / "overnight"
ANALYSIS = PROJECT_ROOT / "outputs" / "analysis"
OUT_JSON = ANALYSIS / "blind_zero_prediction.json"
OUT_TXT = ANALYSIS / "blind_zero_prediction.txt"

# ---------- 실험 설정 ----------
IN_FEATURES_LIST = [64, 128]
SEEDS = [42, 7, 123, 2026, 314]
HIDDEN = 64
EPOCHS = 200
NUM_POINTS = 500
DENSE_GRID_N = 2000
HIT_THRESHOLD = 0.3  # 전형적 영점 간격 ~0.6 의 절반

# 실험 방향 정의: (이름, 학습범위, 예측범위)
DIRECTIONS = [
    ("forward",  100, 150, 150, 200),   # 순방향
    ("reverse",  150, 200, 100, 150),   # 역방향
]


def train_one(train_loader, val_loader, in_features, seed, label):
    """단일 모델 학습 — early stopping 포함."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = MasterResonantNetwork(
        in_features=in_features, hidden_features=HIDDEN, out_features=2,
        num_layers=3, channel_type="paper3ch", damping_mode="paper",
    )
    n_params = sum(p.numel() for p in model.parameters())
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1,
        lambda_tgt=1.0, lambda_pqo=0.5, pqo_mode="cos2",
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    best_val = float("inf")
    best_state = None
    best_ep = 0
    t0 = time.time()

    for epoch in range(1, EPOCHS + 1):
        model.train()
        for X_batch, _ in train_loader:
            X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)
            optimizer.zero_grad(set_to_none=True)
            outputs = model(X_in)
            total_loss, _ = loss_fn(**outputs)
            if torch.isnan(total_loss) or torch.isinf(total_loss):
                continue
            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()

        model.eval()
        va = 0.0
        nv = 0
        with torch.enable_grad():
            for X_batch, _ in val_loader:
                X_in = X_batch.to(dtype=PrecisionManager.REAL_DTYPE)
                X_in.requires_grad_(True)
                o = model(X_in)
                vl, _ = loss_fn(**o)
                if not (torch.isnan(vl) or torch.isinf(vl)):
                    va += vl.item()
                    nv += 1
        vl_m = va / max(1, nv)
        if vl_m < best_val:
            best_val = vl_m
            best_state = copy.deepcopy(model.state_dict())
            best_ep = epoch
        if epoch % 50 == 0 or epoch == EPOCHS:
            print(f"    {label} ep{epoch:03d}: val={vl_m:.5f} best={best_val:.5f}@{best_ep}",
                  flush=True)

    if best_state:
        model.load_state_dict(best_state)
    return model, best_val, best_ep, time.time() - t0, n_params


def eval_F2(model, features):
    """F₂ 잔차 계산 — features 는 [N, in_features] 텐서."""
    model.eval()
    with torch.enable_grad():
        X = features.clone().requires_grad_(True)
        out = model(X)
    phi = out["phi"].detach()
    psi = out["psi"].detach()
    L_G = out["L_G"].detach()
    phi_real = phi.to(dtype=PrecisionManager.REAL_DTYPE)
    rot = torch.complex(torch.cos(phi_real), -torch.sin(phi_real))
    psi_c = psi.to(dtype=PrecisionManager.COMPLEX_DTYPE)
    return (rot * (L_G - psi_c)).imag.mean(dim=-1).cpu().numpy()


def find_local_minima(t_grid, abs_f2, order=5):
    """
    |F₂(t)| 배열에서 극소값 위치 탐색.
    order: 양쪽 order 개의 이웃보다 작아야 극소값으로 인정.
    반환: 극소값에 해당하는 t 값 배열.
    """
    indices = argrelmin(abs_f2, order=order)[0]
    return t_grid[indices], abs_f2[indices]


def match_predictions(predicted_t, actual_zeros, threshold):
    """
    예측 영점과 실제 영점 매칭.

    반환: (hits, precision, recall, f1, details)
      - hits: 적중 개수 (실제 영점 기준)
      - precision: 적중 예측 / 전체 예측
      - recall: 적중 실제 / 전체 실제
      - details: 각 실제 영점에 대한 최근접 예측과 거리
    """
    actual = np.array(actual_zeros)
    pred = np.array(predicted_t)
    details = []
    hit_count = 0
    matched_pred = set()  # 이미 매칭된 예측 인덱스

    for z in actual:
        if len(pred) == 0:
            details.append({"actual": float(z), "nearest_pred": None,
                            "distance": float("inf"), "hit": False})
            continue
        dists = np.abs(pred - z)
        nearest_idx = int(np.argmin(dists))
        dist = float(dists[nearest_idx])
        hit = dist < threshold
        if hit:
            hit_count += 1
            matched_pred.add(nearest_idx)
        details.append({
            "actual": float(z),
            "nearest_pred": float(pred[nearest_idx]),
            "distance": dist,
            "hit": hit,
        })

    n_pred = len(pred)
    n_actual = len(actual)
    # 정밀도: 실제 영점과 매칭된 예측의 비율
    precision = len(matched_pred) / n_pred if n_pred > 0 else 0.0
    recall = hit_count / n_actual if n_actual > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)
          if (precision + recall) > 0 else 0.0)

    return hit_count, precision, recall, f1, details


def run_direction(direction_name, train_min, train_max, eval_min, eval_max,
                  train_zeros, test_zeros, all_zeros, log_fn):
    """하나의 방향(순방향/역방향) 실험 전체를 수행."""
    log = log_fn
    log(f"\n{'=' * 72}")
    log(f"  방향: {direction_name}")
    log(f"  학습: t∈[{train_min},{train_max}] ({len(train_zeros)}개 영점)")
    log(f"  예측: t∈[{eval_min},{eval_max}] ({len(test_zeros)}개 영점)")
    log(f"{'=' * 72}")

    dense_grid = np.linspace(eval_min, eval_max, DENSE_GRID_N)
    dense_grid_list = dense_grid.tolist()

    direction_results = {}

    for ifx in IN_FEATURES_LIST:
        log(f"\n  --- in_features={ifx} ---")

        # 학습 범위로 데이터셋 구축
        train_loader, val_loader, dataset = get_xi_feature_dataloaders(
            in_features=ifx, batch_size=32,
            t_min=train_min, t_max=train_max, num_points=NUM_POINTS,
            cache_dir=str(OUTPUT_DIR), zeros_list=train_zeros,
        )

        # 밀집 격자에서의 특징 벡터 (학습 범위의 푸리에 기저로 외삽)
        dense_features = dataset.get_features_at_t(dense_grid_list)

        # 실제 영점 위치의 특징 벡터 (참고용: 학습 범위 내 영점)
        train_zero_features = dataset.get_features_at_t(train_zeros)

        seed_results = {}

        for seed in SEEDS:
            log(f"\n    seed={seed}:")
            model, bv, be, dt, n_params = train_one(
                train_loader, val_loader, ifx, seed,
                f"{direction_name}/if={ifx}/s={seed}",
            )

            # --- 학습 범위 내 |F₂| (건전성 확인) ---
            f2_train = eval_F2(model, train_zero_features)
            train_abs_mean = float(np.abs(f2_train).mean())
            log(f"      학습 범위 |F₂| mean={train_abs_mean:.5f} (건전성)")

            # --- 블라인드 영역 밀집 격자 평가 ---
            f2_dense = eval_F2(model, dense_features)
            abs_f2_dense = np.abs(f2_dense)

            # 극소값 탐색
            pred_t, pred_vals = find_local_minima(dense_grid, abs_f2_dense, order=5)
            log(f"      밀집 격자 극소값: {len(pred_t)}개 발견")

            # 매칭
            hits, prec, rec, f1, details = match_predictions(
                pred_t, test_zeros, HIT_THRESHOLD,
            )
            log(f"      적중: {hits}/{len(test_zeros)}  "
                f"정밀도={prec:.3f}  재현율={rec:.3f}  F1={f1:.3f}")

            # 개별 영점 상세
            for d in details:
                marker = "HIT" if d["hit"] else "MISS"
                np_str = (f"{d['nearest_pred']:.4f}" if d["nearest_pred"] is not None
                          else "N/A")
                log(f"        t={d['actual']:.4f} → 최근접={np_str} "
                    f"(거리={d['distance']:.4f}) [{marker}]")

            seed_results[seed] = {
                "seed": seed, "best_val": bv, "best_epoch": be,
                "train_time": dt, "n_params": n_params,
                "train_F2_mean": train_abs_mean,
                "n_minima": len(pred_t),
                "predicted_zeros": pred_t.tolist(),
                "predicted_vals": pred_vals.tolist(),
                "hits": hits, "precision": prec, "recall": rec, "f1": f1,
                "details": details,
                "dense_F2": f2_dense.tolist(),
            }

        # in_features 요약
        recalls = [seed_results[s]["recall"] for s in SEEDS]
        f1s = [seed_results[s]["f1"] for s in SEEDS]
        precs = [seed_results[s]["precision"] for s in SEEDS]
        log(f"\n    if={ifx} 앙상블 요약:")
        log(f"      재현율: {np.mean(recalls):.3f} ± {np.std(recalls):.3f}")
        log(f"      정밀도: {np.mean(precs):.3f} ± {np.std(precs):.3f}")
        log(f"      F1:     {np.mean(f1s):.3f} ± {np.std(f1s):.3f}")

        direction_results[ifx] = {
            "in_features": ifx,
            "per_seed": {str(s): seed_results[s] for s in SEEDS},
            "ensemble": {
                "recall_mean": float(np.mean(recalls)),
                "recall_std": float(np.std(recalls)),
                "precision_mean": float(np.mean(precs)),
                "precision_std": float(np.std(precs)),
                "f1_mean": float(np.mean(f1s)),
                "f1_std": float(np.std(f1s)),
            },
        }

    return direction_results


def main():
    PrecisionManager.setup_precision()
    ANALYSIS.mkdir(parents=True, exist_ok=True)

    lines = []
    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    log("=" * 72)
    log("  blind_zero_prediction — 블라인드 영점 예측 실험")
    log("=" * 72)
    log(f"  in_features: {IN_FEATURES_LIST}")
    log(f"  seeds: {SEEDS}")
    log(f"  epochs: {EPOCHS}, hidden: {HIDDEN}")
    log(f"  밀집 격자: {DENSE_GRID_N}점, 적중 임계값: {HIT_THRESHOLD}")

    # 영점 로드
    with open(OUTPUT_DIR / "result_t100-200.json") as f:
        rd = json.load(f)
    all_zeros = rd["zeros_list"]
    all_zeros_arr = np.array(all_zeros)
    log(f"  전체 영점 수 (t∈[100,200]): {len(all_zeros)}")

    # 영점 분할
    train_fwd = [z for z in all_zeros if 100 <= z <= 150]
    test_fwd = [z for z in all_zeros if 150 < z <= 200]
    train_rev = [z for z in all_zeros if 150 <= z <= 200]
    test_rev = [z for z in all_zeros if 100 <= z < 150]

    log(f"  순방향: 학습 {len(train_fwd)}개, 예측 {len(test_fwd)}개")
    log(f"  역방향: 학습 {len(train_rev)}개, 예측 {len(test_rev)}개")

    t_start = time.time()
    all_results = {}

    # 순방향: [100,150] → [150,200]
    all_results["forward"] = run_direction(
        "forward", 100, 150, 150, 200,
        train_fwd, test_fwd, all_zeros, log,
    )

    # 역방향: [150,200] → [100,150]
    all_results["reverse"] = run_direction(
        "reverse", 150, 200, 100, 150,
        train_rev, test_rev, all_zeros, log,
    )

    total_time = time.time() - t_start

    # ==================== 최종 요약 ====================
    log(f"\n{'=' * 72}")
    log("  최종 요약")
    log(f"{'=' * 72}")
    log(f"  총 실험 시간: {total_time:.1f}s")

    log(f"\n  {'방향':>8s} {'if':>4s} {'recall':>10s} {'precision':>10s} {'F1':>10s}")
    log(f"  {'-'*8:>8s} {'-'*4:>4s} {'-'*10:>10s} {'-'*10:>10s} {'-'*10:>10s}")
    for direction in ["forward", "reverse"]:
        for ifx in IN_FEATURES_LIST:
            ens = all_results[direction][ifx]["ensemble"]
            log(f"  {direction:>8s} {ifx:>4d} "
                f"{ens['recall_mean']:.3f}±{ens['recall_std']:.3f} "
                f"{ens['precision_mean']:.3f}±{ens['precision_std']:.3f} "
                f"{ens['f1_mean']:.3f}±{ens['f1_std']:.3f}")

    # 판정
    log(f"\n  판정 기준: 재현율 ≥ 0.5 → 유의미한 외삽 예측력")
    best_recall = 0.0
    best_config = ""
    for direction in ["forward", "reverse"]:
        for ifx in IN_FEATURES_LIST:
            r = all_results[direction][ifx]["ensemble"]["recall_mean"]
            if r > best_recall:
                best_recall = r
                best_config = f"{direction}/if={ifx}"

    if best_recall >= 0.5:
        log(f"  결과: 블라인드 예측 성공 ({best_config}, recall={best_recall:.3f})")
    elif best_recall >= 0.3:
        log(f"  결과: 부분적 예측력 ({best_config}, recall={best_recall:.3f})")
    else:
        log(f"  결과: 외삽 예측 실패 ({best_config}, recall={best_recall:.3f})")
        log("  해석: 푸리에 기저가 학습 범위에 종속 — 외삽 시 기저 표현력 상실")

    # JSON 저장 (dense_F2 는 용량이 크므로 별도 키로 분리)
    save_data = {
        "config": {
            "in_features_list": IN_FEATURES_LIST,
            "seeds": SEEDS,
            "hidden": HIDDEN,
            "epochs": EPOCHS,
            "dense_grid_n": DENSE_GRID_N,
            "hit_threshold": HIT_THRESHOLD,
        },
        "zeros": {
            "all": all_zeros,
            "forward_train": train_fwd,
            "forward_test": test_fwd,
            "reverse_train": train_rev,
            "reverse_test": test_rev,
        },
        "results": {},
    }
    # 결과 저장 (dense_F2 제외하여 JSON 크기 절약)
    for direction in ["forward", "reverse"]:
        save_data["results"][direction] = {}
        for ifx in IN_FEATURES_LIST:
            dr = all_results[direction][ifx]
            compact_seeds = {}
            for s_key, s_val in dr["per_seed"].items():
                compact = {k: v for k, v in s_val.items() if k != "dense_F2"}
                compact_seeds[s_key] = compact
            save_data["results"][direction][ifx] = {
                "in_features": ifx,
                "per_seed": compact_seeds,
                "ensemble": dr["ensemble"],
            }

    with open(OUT_JSON, "w") as f:
        json.dump(save_data, f, indent=2, ensure_ascii=False)
    with open(OUT_TXT, "w") as f:
        f.write("\n".join(lines))
    log(f"\n저장: {OUT_JSON}")
    log(f"      {OUT_TXT}")


if __name__ == "__main__":
    main()
