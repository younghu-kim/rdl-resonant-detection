#!/usr/bin/env python3
"""
=============================================================================
야간 자동 탐사: ξ-함수 영점의 수학적 풍경 매핑
=============================================================================
t∈[10,50]부터 t∈[5000,10000]까지 여러 구간을 순회하며
RDL 모델을 훈련하고 수학적 특징을 기록한다.

기록 항목:
- 각 구간의 영점 개수, 탐지율, cos² PQO 품질
- 게이지 임계 전이 시점 (φ std 급감 에폭)
- 영점별 |F₂| 난이도 분포 (Lehmer 현상 추적)
- 모델 크기 대비 탐지 성능 스케일링

예상 실행 시간: ~5-6시간

실행: ~/qrop_env/bin/python3 scripts/overnight_exploration.py
"""

import os
import sys
import json
import time
import copy
import math
from collections import defaultdict

import torch
import torch.nn as nn
from torch.utils.data import DataLoader

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_ROOT)

from gdl.rdl.constants import R_CONST, PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.pipeline.xi_feature_dataset import (
    get_xi_feature_dataloaders, compute_zeros_in_range, XiFeatureDataset,
)

# 출력 디렉토리
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs", "overnight")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# 탐사 구간 정의
# =============================================================================
EXPLORATION_RANGES = [
    # (t_min, t_max, hidden_dim, epochs, description)
    (10,     50,     32,  200,  "기본 범위 (검증 기준선)"),
    (100,    200,    64,  200,  "고위 영점 (기존 검증)"),
    (200,    500,    64,  300,  "미탐사 1구간"),
    (500,    1000,   96,  300,  "미탐사 2구간"),
    (1000,   2000,   96,  400,  "미탐사 3구간"),
    (2000,   5000,   128, 400,  "미탐사 4구간"),
    (5000,   10000,  128, 500,  "미탐사 5구간 (최고위)"),
]


# =============================================================================
# 단일 구간 훈련 + 수학적 특징 추출
# =============================================================================
def run_single_range(t_min, t_max, hidden, epochs, description, run_idx):
    """한 구간을 훈련하고 수학적 특징을 반환."""

    range_tag = f"t{int(t_min)}-{int(t_max)}"
    print(f"\n{'='*70}")
    print(f"  [{run_idx}] {description}: t∈[{t_min},{t_max}], hidden={hidden}, ep={epochs}")
    print(f"{'='*70}")

    # 1. 영점 계산
    t0 = time.time()
    print(f"  영점 계산 중...")
    try:
        zeros_list = compute_zeros_in_range(t_min, t_max)
    except Exception as e:
        print(f"  ERROR: 영점 계산 실패 — {e}")
        return None
    n_zeros = len(zeros_list)
    dt_zeros = time.time() - t0
    print(f"  {n_zeros}개 영점 발견 ({dt_zeros:.1f}초)")

    if n_zeros == 0:
        print(f"  SKIP: 영점 없음")
        return None

    # 2. 데이터 준비
    print(f"  데이터 준비 중...")
    in_features = 64
    num_points = min(500, max(200, n_zeros * 10))

    try:
        train_loader, val_loader, dataset = get_xi_feature_dataloaders(
            in_features=in_features,
            batch_size=32,
            t_min=t_min, t_max=t_max,
            num_points=num_points,
            cache_dir=OUTPUT_DIR,
            zeros_list=zeros_list,
        )
    except Exception as e:
        print(f"  ERROR: 데이터 준비 실패 — {e}")
        return None

    # 3. 모델
    model = MasterResonantNetwork(
        in_features=in_features,
        hidden_features=hidden,
        out_features=2,
        num_layers=3,
        channel_type="paper3ch",
        damping_mode="paper",
    )
    n_params = sum(p.numel() for p in model.parameters())
    print(f"  모델: {n_params:,} params")

    # 4. 손실 + 옵티마이저 (항상 Adam + cos²)
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=0.1,
        lambda_tgt=1.0, lambda_pqo=0.5,
        pqo_mode="cos2",
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    device = torch.device("cpu")

    # 5. 영점 특징 벡터 사전 생성 (평가용)
    zero_features = dataset.get_features_at_t(zeros_list).to(device)

    # 6. 훈련 + 에폭별 기록
    epoch_records = []
    best_val = float("inf")
    best_state = None

    print(f"  훈련 시작...")
    t_train_start = time.time()

    for epoch in range(1, epochs + 1):
        # --- Train ---
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

        train_loss = train_loss_acc / max(1, n_batch)

        # --- Validate ---
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

        # --- 영점 평가 (매 10 에폭 + 처음 5 + 마지막) ---
        if epoch % 10 == 0 or epoch <= 5 or epoch == epochs:
            with torch.enable_grad():
                X_z = zero_features.clone().requires_grad_(True)
                out_z = model(X_z)

            Z_out = out_z["Z_out"]
            phi = out_z["phi"]
            psi = out_z["psi"]
            L_G = out_z["L_G"]

            # F₂
            phi_real = phi.to(dtype=PrecisionManager.REAL_DTYPE)
            rotation = torch.complex(torch.cos(phi_real), -torch.sin(phi_real))
            psi_complex = psi.to(dtype=PrecisionManager.COMPLEX_DTYPE)
            residual = rotation * (L_G - psi_complex)
            F2 = residual.imag.mean(dim=-1).detach()

            # PQO
            phase_z = torch.angle(Z_out).mean(dim=-1).detach()
            L = R_CONST.PQO_L_DEFAULT
            pqo_cos2 = (torch.cos(L * phase_z / 2.0) ** 2).detach()

            # 게이지 통계
            phi_std = phi.std().item()

            rec = {
                "epoch": epoch,
                "train_loss": train_loss,
                "val_loss": val_loss,
                "F2_abs": [abs(f) for f in F2.cpu().tolist()],
                "F2_mean": float(F2.abs().mean()),
                "F2_max": float(F2.abs().max()),
                "cos2_pqo": pqo_cos2.cpu().tolist(),
                "cos2_max": float(pqo_cos2.max()),
                "phase_mean": float(phase_z.mean()),
                "phi_std": phi_std,
            }
            epoch_records.append(rec)

            if epoch % 50 == 0 or epoch == epochs:
                detected = sum(1 for f in F2.abs().cpu().tolist() if f < 0.1)
                print(f"    [ep {epoch:03d}] loss={train_loss:.4f} val={val_loss:.4f} "
                      f"|F2|={float(F2.abs().mean()):.4f} cos2_max={float(pqo_cos2.max()):.4f} "
                      f"phi_std={phi_std:.4f} detect={detected}/{n_zeros}")

    t_train = time.time() - t_train_start

    # 7. 최적 모델로 최종 평가
    if best_state:
        model.load_state_dict(best_state)

    model.eval()
    with torch.enable_grad():
        X_z = zero_features.clone().requires_grad_(True)
        out_z = model(X_z)

    Z_out = out_z["Z_out"]
    phi = out_z["phi"]
    L_G = out_z["L_G"]
    psi = out_z["psi"]

    phi_real = phi.to(dtype=PrecisionManager.REAL_DTYPE)
    rotation = torch.complex(torch.cos(phi_real), -torch.sin(phi_real))
    psi_complex = psi.to(dtype=PrecisionManager.COMPLEX_DTYPE)
    residual = rotation * (L_G - psi_complex)
    final_F2 = residual.imag.mean(dim=-1).detach().cpu()

    phase_z = torch.angle(Z_out).mean(dim=-1).detach().cpu()
    L = R_CONST.PQO_L_DEFAULT
    final_cos2 = (torch.cos(L * phase_z / 2.0) ** 2)

    detected = sum(1 for f in final_F2.abs().tolist() if f < 0.1)

    # 영점별 난이도 정렬
    zero_difficulty = sorted(
        zip(zeros_list, final_F2.abs().tolist(), final_cos2.tolist()),
        key=lambda x: -x[1]  # |F₂| 내림차순 = 어려운 순
    )

    # 게이지 전이 시점 탐지 (phi_std가 처음으로 급감하는 에폭)
    phi_stds = [(r["epoch"], r["phi_std"]) for r in epoch_records]
    transition_epoch = None
    if len(phi_stds) >= 3:
        for i in range(1, len(phi_stds)):
            if phi_stds[i][1] < phi_stds[0][1] * 0.5:
                transition_epoch = phi_stds[i][0]
                break

    # 결과 요약
    result = {
        "range": f"[{t_min},{t_max}]",
        "t_min": t_min,
        "t_max": t_max,
        "description": description,
        "n_zeros": n_zeros,
        "hidden": hidden,
        "epochs": epochs,
        "n_params": n_params,
        "train_time_sec": round(t_train, 1),
        "zeros_compute_sec": round(dt_zeros, 1),
        "detection_rate": f"{detected}/{n_zeros}",
        "detected": detected,
        "F2_mean": round(float(final_F2.abs().mean()), 6),
        "F2_max": round(float(final_F2.abs().max()), 6),
        "cos2_max": round(float(final_cos2.max()), 6),
        "cos2_mean": round(float(final_cos2.mean()), 6),
        "best_val_loss": round(best_val, 6),
        "gauge_transition_epoch": transition_epoch,
        "hardest_zeros": [(round(t, 4), round(f2, 6)) for t, f2, _ in zero_difficulty[:5]],
        "easiest_zeros": [(round(t, 4), round(f2, 6)) for t, f2, _ in zero_difficulty[-5:]],
        "epoch_records": epoch_records,
        "zeros_list": zeros_list,
        "final_F2": final_F2.tolist(),
        "final_cos2": final_cos2.tolist(),
        "final_phase": phase_z.tolist(),
    }

    # 개별 결과 저장
    result_path = os.path.join(OUTPUT_DIR, f"result_{range_tag}.json")
    with open(result_path, "w") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    print(f"  저장: {result_path}")

    return result


# =============================================================================
# 전체 탐사 요약 생성
# =============================================================================
def generate_summary(all_results):
    """전체 결과를 하나의 요약 파일로."""

    summary_path = os.path.join(OUTPUT_DIR, "exploration_summary.json")
    summary_txt_path = os.path.join(OUTPUT_DIR, "exploration_summary.txt")

    # JSON
    summary_data = []
    for r in all_results:
        summary_data.append({
            "range": r["range"],
            "n_zeros": r["n_zeros"],
            "detected": r["detected"],
            "detection_rate": r["detection_rate"],
            "F2_mean": r["F2_mean"],
            "F2_max": r["F2_max"],
            "cos2_max": r["cos2_max"],
            "n_params": r["n_params"],
            "train_time_sec": r["train_time_sec"],
            "gauge_transition_epoch": r["gauge_transition_epoch"],
            "hardest_zeros": r["hardest_zeros"][:3],
        })

    with open(summary_path, "w") as f:
        json.dump(summary_data, f, indent=2, ensure_ascii=False)

    # 텍스트 요약
    lines = []
    lines.append("=" * 80)
    lines.append("  RDL Xi-Function Overnight Exploration — 수학적 풍경 요약")
    lines.append("=" * 80)
    lines.append("")
    lines.append(f"{'Range':>20} {'Zeros':>6} {'Detect':>8} {'|F2| mean':>10} "
                 f"{'cos2 max':>10} {'Params':>8} {'Time':>8} {'Transition':>10}")
    lines.append("-" * 80)

    for r in all_results:
        lines.append(
            f"{r['range']:>20} {r['n_zeros']:>6} {r['detection_rate']:>8} "
            f"{r['F2_mean']:>10.4f} {r['cos2_max']:>10.4f} "
            f"{r['n_params']:>8,} {r['train_time_sec']:>7.0f}s "
            f"{str(r['gauge_transition_epoch']):>10}"
        )

    lines.append("")
    lines.append("=" * 80)
    lines.append("  수학적 관찰")
    lines.append("=" * 80)
    lines.append("")

    # 스케일링 분석
    if len(all_results) >= 2:
        f2_means = [(r["t_max"], r["F2_mean"]) for r in all_results]
        lines.append("  [1] |F₂| vs 높이 (영점 탐지 난이도 스케일링)")
        for t_max, f2m in f2_means:
            bar = "#" * int(min(50, f2m * 500))
            lines.append(f"      t<{t_max:>6}: |F₂|={f2m:.4f} {bar}")

        lines.append("")
        lines.append("  [2] cos² PQO 품질 vs 높이 (Maslov 보정 유효 범위)")
        for r in all_results:
            lines.append(f"      t<{r['t_max']:>6}: cos²_max={r['cos2_max']:.4f} "
                         f"cos²_mean={r['cos2_mean']:.6f}")

        lines.append("")
        lines.append("  [3] 게이지 임계 전이 (Kuramoto 동기화)")
        for r in all_results:
            ep = r["gauge_transition_epoch"]
            lines.append(f"      t<{r['t_max']:>6}: 전이 에폭 = {ep}")

        lines.append("")
        lines.append("  [4] 가장 어려운 영점 (Lehmer 현상 후보)")
        for r in all_results:
            if r["hardest_zeros"]:
                hardest = r["hardest_zeros"][0]
                lines.append(f"      t<{r['t_max']:>6}: t={hardest[0]:.2f}, |F₂|={hardest[1]:.4f}")

    txt = "\n".join(lines)
    with open(summary_txt_path, "w") as f:
        f.write(txt)

    print(f"\n{txt}")
    print(f"\n저장: {summary_txt_path}")
    print(f"저장: {summary_path}")


# =============================================================================
# 메인
# =============================================================================
def generate_next_ranges(last_t_max):
    """마지막 구간 이후 자동으로 다음 구간들을 생성. 멈출 때까지 계속."""
    # 10000 이후부터는 100 단위로 균등 저장
    t_min = last_t_max
    if t_min < 10000:
        width = 5000
    else:
        width = 100

    t_max = t_min + width

    # 높이에 따라 모델 크기 조정
    if t_max <= 10000:
        hidden, epochs = 128, 500
    else:
        # 100 단위 구간: 영점 ~15-20개이므로 작은 모델로 빠르게
        hidden, epochs = 64, 300

    desc = f"자동 탐사 t∈[{int(t_min)},{int(t_max)}]"
    return [(t_min, t_max, hidden, epochs, desc)]


def main():
    PrecisionManager.setup_precision()
    torch.manual_seed(42)

    print("=" * 80)
    print("  RDL Xi-Function Continuous Exploration")
    print(f"  초기 구간: {len(EXPLORATION_RANGES)}개 → 이후 무한 확장")
    print(f"  시작: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  중지: kill {os.getpid()} 또는 Ctrl+C")
    print("=" * 80)

    all_results = []
    t_total_start = time.time()
    run_idx = 0

    # Phase 1: 사전 정의 구간
    ranges_queue = list(EXPLORATION_RANGES)

    while True:
        if not ranges_queue:
            # Phase 2: 자동 확장 — 마지막 구간 이후부터 계속
            last_t_max = all_results[-1]["t_max"] if all_results else 10000
            ranges_queue = generate_next_ranges(last_t_max)
            print(f"\n>>> 자동 확장: 다음 구간 {ranges_queue[0][0]}-{ranges_queue[0][1]}")

        t_min, t_max, hidden, epochs, desc = ranges_queue.pop(0)
        run_idx += 1

        # 이미 완료된 구간 스킵
        result_path = os.path.join(OUTPUT_DIR, f"result_t{int(t_min)}-{int(t_max)}.json")
        if os.path.exists(result_path):
            print(f"\n  [{run_idx}] {desc} — 이미 완료, 스킵")
            # 기존 결과 로드
            with open(result_path) as f:
                existing = json.load(f)
            all_results.append(existing)
            continue

        try:
            result = run_single_range(t_min, t_max, hidden, epochs, desc, run_idx)
            if result:
                all_results.append(result)
        except KeyboardInterrupt:
            print(f"\n\n사용자 중지 (Ctrl+C)")
            break
        except Exception as e:
            print(f"\n  FATAL ERROR in [{t_min},{t_max}]: {e}")
            import traceback
            traceback.print_exc()
            continue

        # 매 구간 완료 시 요약 갱신
        if all_results:
            generate_summary(all_results)

        elapsed = (time.time() - t_total_start) / 3600
        print(f"\n  [경과 {elapsed:.1f}시간, 완료 {len(all_results)}개 구간]")

    t_total = time.time() - t_total_start
    print(f"\n\n총 탐사 시간: {t_total/3600:.1f}시간")
    print(f"완료 구간: {len(all_results)}개")
    print(f"종료: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    if all_results:
        generate_summary(all_results)


if __name__ == "__main__":
    main()
