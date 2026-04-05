#!/usr/bin/env python3
"""
=============================================================================
[Project RDL] Xi-Function Deep Learning Pipeline — 실전 훈련 & 수학적 특징 추출
=============================================================================
ξ(1/2+it) 데이터로 공명 딥러닝 모델을 훈련하고,
훈련 과정에서 수학적 특징(위상 궤적, F₂ 판별식, PQO 격자, 게이지 진화)을 추출한다.

실행: python scripts/run_xi_pipeline.py --optimizer adam --pqo-mode cos2 --epochs 200
"""

import os
import sys
import math
import json
import time
import copy
import argparse
from collections import defaultdict

import torch
import torch.nn as nn
from torch.utils.data import DataLoader

# gdl_unified 루트를 경로에 추가
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_ROOT)

from gdl.rdl.constants import R_CONST, PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.optim.pggd import PGGD
from gdl.rdl.pipeline.xi_feature_dataset import (
    get_xi_feature_dataloaders, KNOWN_ZEROS, XiFeatureDataset,
    compute_zeros_in_range,
)

# 출력 디렉토리
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs")
PLOT_DIR = os.path.join(OUTPUT_DIR, "plots")
LOG_DIR = os.path.join(OUTPUT_DIR, "logs")


# =============================================================================
# 수학적 특징 추출기
# =============================================================================
class MathFeatureExtractor:
    """훈련 루프에서 epoch마다 수학적 상태를 캡처"""

    def __init__(self, dataset: XiFeatureDataset, device: torch.device, zeros_list=None):
        self.dataset = dataset
        self.device = device
        self.zeros_list = zeros_list if zeros_list is not None else KNOWN_ZEROS

        self.zero_features = dataset.get_features_at_t(self.zeros_list).to(device)
        self.zero_t = torch.tensor(self.zeros_list, dtype=PrecisionManager.REAL_DTYPE)

        t_min, t_max = dataset.t.min().item(), dataset.t.max().item()
        self.dense_t = torch.linspace(t_min, t_max, 200, dtype=PrecisionManager.REAL_DTYPE)
        self.dense_features = dataset.get_features_at_t(self.dense_t.tolist()).to(device)

        self.records = {
            "epochs": [],
            "loss_total": [], "loss_res": [], "loss_curv": [],
            "loss_tgt": [], "loss_pqo": [], "val_loss": [],
            "F2_at_zeros": [],
            "xi_pred_amp_at_zeros": [],
            "pqo_at_zeros": [],
            "pqo_cos2_at_zeros": [],
            "phase_at_zeros": [],
            "phase_curve": [],
            "amp_curve": [],
            "phi_mean": [], "phi_std": [],
            "psi_mean": [], "psi_std": [],
            "optim_phi_hist": [],
        }

    @torch.no_grad()
    def capture(self, model, optimizer, epoch, train_metrics, val_loss):
        """한 epoch 종료 후 수학적 특징 스냅샷"""
        model.eval()

        self.records["epochs"].append(epoch)
        self.records["loss_total"].append(train_metrics["Total"])
        self.records["loss_res"].append(train_metrics["L_res"])
        self.records["loss_curv"].append(train_metrics["L_curv"])
        self.records["loss_tgt"].append(train_metrics["L_tgt"])
        self.records["loss_pqo"].append(train_metrics.get("L_pqo", 0.0))
        self.records["val_loss"].append(val_loss)

        with torch.enable_grad():
            X_z = self.zero_features.clone().requires_grad_(True)
            out_z = model(X_z)

        Z_out = out_z["Z_out"]
        phi = out_z["phi"]
        psi = out_z["psi"]
        L_G = out_z["L_G"]

        phi_real = phi.to(dtype=PrecisionManager.REAL_DTYPE)
        rotation = torch.complex(torch.cos(phi_real), -torch.sin(phi_real))
        psi_complex = psi.to(dtype=PrecisionManager.COMPLEX_DTYPE)
        residual = rotation * (L_G - psi_complex)
        F2 = residual.imag.mean(dim=-1)

        z_amp = Z_out.abs().mean(dim=-1)
        phase_z = torch.angle(Z_out).mean(dim=-1)
        L = R_CONST.PQO_L_DEFAULT
        pqo_sin2 = torch.sin(L * phase_z / 2.0) ** 2
        pqo_cos2 = torch.cos(L * phase_z / 2.0) ** 2

        self.records["F2_at_zeros"].append(F2.cpu().tolist())
        self.records["xi_pred_amp_at_zeros"].append(z_amp.cpu().tolist())
        self.records["pqo_at_zeros"].append(pqo_sin2.cpu().tolist())
        self.records["pqo_cos2_at_zeros"].append(pqo_cos2.cpu().tolist())
        self.records["phase_at_zeros"].append(phase_z.cpu().tolist())

        if epoch % 10 == 0 or epoch <= 5:
            with torch.enable_grad():
                X_d = self.dense_features.clone().requires_grad_(True)
                out_d = model(X_d)
            phase_curve = torch.angle(out_d["Z_out"]).mean(dim=-1).cpu().tolist()
            amp_curve = out_d["Z_out"].abs().mean(dim=-1).cpu().tolist()
            self.records["phase_curve"].append((epoch, phase_curve))
            self.records["amp_curve"].append((epoch, amp_curve))

        self.records["phi_mean"].append(phi.mean().item())
        self.records["phi_std"].append(phi.std().item())
        self.records["psi_mean"].append(psi.mean().item())
        self.records["psi_std"].append(psi.std().item())

        if epoch % 20 == 0 and isinstance(optimizer, PGGD):
            phi_opts = []
            for group in optimizer.param_groups:
                for p in group["params"]:
                    if p in optimizer.state and "phi_opt" in optimizer.state[p]:
                        phi_opts.append(
                            optimizer.state[p]["phi_opt"].flatten().cpu().tolist()[:50]
                        )
                        break
                if phi_opts:
                    break
            if phi_opts:
                self.records["optim_phi_hist"].append((epoch, phi_opts[0]))

        model.train()

    def save(self, path):
        with open(path, "w") as f:
            json.dump(self.records, f, indent=2)
        print(f"[Feature] 수학적 특징 저장: {path}")


# =============================================================================
# 계측 트레이너
# =============================================================================
class InstrumentedTrainer:
    def __init__(self, model, train_loader, val_loader, optimizer,
                 loss_fn, device, extractor, max_grad_norm=5.0):
        self.model = model.to(device)
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.optimizer = optimizer
        self.loss_fn = loss_fn.to(device)
        self.device = device
        self.extractor = extractor
        self.max_grad_norm = max_grad_norm

    def train_epoch(self):
        self.model.train()
        metrics_acc = defaultdict(float)
        num_batches = 0
        nan_count = 0

        for X_batch, _ in self.train_loader:
            X_in = X_batch.to(self.device, dtype=PrecisionManager.REAL_DTYPE)
            X_in.requires_grad_(True)

            self.optimizer.zero_grad(set_to_none=True)
            outputs = self.model(X_in)
            total_loss, metrics = self.loss_fn(**outputs)

            if torch.isnan(total_loss) or torch.isinf(total_loss):
                nan_count += 1
                if nan_count > 5:
                    print("[WARN] NaN/Inf 5회 초과, epoch 중단")
                    break
                continue

            total_loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), self.max_grad_norm)
            self.optimizer.step()

            for k, v in metrics.items():
                metrics_acc[k] += v
            num_batches += 1

        if num_batches == 0:
            return {"Total": float("nan"), "L_res": 0, "L_curv": 0, "L_tgt": 0, "L_pqo": 0}
        return {k: v / num_batches for k, v in metrics_acc.items()}

    @torch.no_grad()
    def validate_epoch(self):
        self.model.eval()
        total_val = 0.0
        num_batches = 0

        with torch.enable_grad():
            for X_batch, _ in self.val_loader:
                X_in = X_batch.to(self.device, dtype=PrecisionManager.REAL_DTYPE)
                X_in.requires_grad_(True)
                outputs = self.model(X_in)
                loss, _ = self.loss_fn(**outputs)
                if not (torch.isnan(loss) or torch.isinf(loss)):
                    total_val += loss.item()
                    num_batches += 1

        return total_val / max(1, num_batches)

    def fit(self, epochs, print_freq=5):
        print(f"\n{'='*80}")
        print(f"[RDL] Xi-함수 공명 딥러닝 훈련 시작 (Epochs: {epochs})")
        print(f"{'='*80}\n")

        best_val = float("inf")
        best_model_state = None

        for epoch in range(1, epochs + 1):
            t0 = time.time()
            train_metrics = self.train_epoch()
            val_loss = self.validate_epoch()
            dt = time.time() - t0

            self.extractor.capture(self.model, self.optimizer, epoch, train_metrics, val_loss)

            if val_loss < best_val:
                best_val = val_loss
                best_model_state = copy.deepcopy(self.model.state_dict())

            if epoch % print_freq == 0 or epoch <= 3 or epoch == epochs:
                F2_vals = self.extractor.records["F2_at_zeros"][-1]
                F2_mean = sum(abs(v) for v in F2_vals) / len(F2_vals)
                print(
                    f"[Epoch {epoch:03d}] "
                    f"Loss: {train_metrics['Total']:.4f} "
                    f"(Res:{train_metrics['L_res']:.3f} Curv:{train_metrics['L_curv']:.3f} "
                    f"Tgt:{train_metrics['L_tgt']:.3f} PQO:{train_metrics.get('L_pqo',0):.3f}) "
                    f"Val:{val_loss:.4f} "
                    f"|F2|_zeros:{F2_mean:.4f} "
                    f"[{dt:.1f}s]"
                )

        print(f"\n{'='*80}")
        print(f"[RDL] 훈련 완료. 최적 Val Loss: {best_val:.4f}")
        print(f"{'='*80}")

        if best_model_state:
            self.model.load_state_dict(best_model_state)

        return self.extractor.records


# =============================================================================
# 메인 실행
# =============================================================================
def main():
    parser = argparse.ArgumentParser(description="RDL Xi-Function Pipeline")
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--num-points", type=int, default=500)
    parser.add_argument("--in-features", type=int, default=64)
    parser.add_argument("--hidden", type=int, default=32)
    parser.add_argument("--layers", type=int, default=3)
    parser.add_argument("--lr", type=float, default=0.005)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--optimizer", choices=["pggd", "adam"], default="pggd")
    parser.add_argument("--lambda-curv", type=float, default=0.1)
    parser.add_argument("--lambda-pqo", type=float, default=0.5)
    parser.add_argument("--pqo-mode", choices=["sin2", "cos2"], default="sin2",
                        help="sin2: 원래 Gate A, cos2: Maslov 보정")
    parser.add_argument("--eta-phi", type=float, default=0.0,
                        help="PGGD 복원 포텐셜 강도")
    parser.add_argument("--t-min", type=float, default=10.0)
    parser.add_argument("--t-max", type=float, default=50.0)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    PrecisionManager.setup_precision()
    torch.manual_seed(args.seed)
    device = torch.device("cpu")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(PLOT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    # 동적 영점 계산
    if args.t_min >= 50.0 or args.t_max > 50.0:
        print(f"\n[Phase 0] t in [{args.t_min},{args.t_max}] ���위 영점 계산 중...")
        active_zeros = compute_zeros_in_range(args.t_min, args.t_max)
        print(f"  {len(active_zeros)}개 영점 발견")
    else:
        active_zeros = [z for z in KNOWN_ZEROS if args.t_min <= z <= args.t_max]

    # 모듈 레벨 영점 교체
    import gdl.rdl.pipeline.xi_feature_dataset as xi_mod
    xi_mod.KNOWN_ZEROS = active_zeros

    print(f"\n[RDL] Xi-Function Deep Learning Pipeline")
    print(f"  Optimizer: {args.optimizer.upper()}")
    print(f"  Epochs: {args.epochs}, Points: {args.num_points}")
    print(f"  t range: [{args.t_min}, {args.t_max}], Known zeros: {len(active_zeros)}")
    print(f"  Architecture: {args.in_features} -> {args.hidden} -> 2, Layers: {args.layers}")
    print(f"  PQO mode: {args.pqo_mode}")

    # 1. 데이터
    print(f"\n[Phase 1] 데이터 준비...")
    train_loader, val_loader, dataset = get_xi_feature_dataloaders(
        in_features=args.in_features,
        batch_size=args.batch_size,
        t_min=args.t_min, t_max=args.t_max,
        num_points=args.num_points,
        cache_dir=OUTPUT_DIR,
        zeros_list=active_zeros,
    )
    print(f"  데이터셋: {len(dataset)} samples")

    # 2. 모델
    print(f"\n[Phase 2] 모델 조립...")
    model = MasterResonantNetwork(
        in_features=args.in_features,
        hidden_features=args.hidden,
        out_features=2,
        num_layers=args.layers,
        channel_type="paper3ch",
        damping_mode="paper",
    )
    param_count = sum(p.numel() for p in model.parameters())
    print(f"  파라미터: {param_count:,}")

    # 3. 손실 + 옵티마이저
    loss_fn = TotalResonanceLoss(
        lambda_res=1.0, lambda_curv=args.lambda_curv,
        lambda_tgt=1.0, lambda_pqo=args.lambda_pqo,
        pqo_mode=args.pqo_mode,
    )

    if args.optimizer == "pggd":
        optimizer = PGGD(model.parameters(), lr=args.lr, tau_g=0.04, h=0.04,
                         eta_phi=args.eta_phi)
    else:
        optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    # 4. 훈련
    extractor = MathFeatureExtractor(dataset, device, zeros_list=active_zeros)
    trainer = InstrumentedTrainer(
        model, train_loader, val_loader, optimizer,
        loss_fn, device, extractor, max_grad_norm=5.0,
    )

    print(f"\n[Phase 3] 훈련 시작...")
    t_start = time.time()
    records = trainer.fit(args.epochs, print_freq=5)
    t_total = time.time() - t_start
    print(f"\n총 훈련 시간: {t_total:.1f}초 ({t_total/60:.1f}분)")

    # 5. 저장
    log_path = os.path.join(LOG_DIR, f"xi_{args.optimizer}_ep{args.epochs}.json")
    extractor.save(log_path)

    # 6. 영점 성적표
    final_F2 = records["F2_at_zeros"][-1]
    final_pqo_cos2 = records.get("pqo_cos2_at_zeros", [None])[-1]
    final_phase = records["phase_at_zeros"][-1]

    print(f"\n{'='*60}")
    print(f"  영점 탐지 최종 성적표")
    print(f"{'='*60}")
    print(f"{'t':>10} {'|F2|':>10} {'cos2':>8} {'phase':>8}")
    print(f"{'-'*40}")
    for i, zt in enumerate(active_zeros):
        if i >= len(final_F2):
            break
        f2 = abs(final_F2[i])
        pc2 = final_pqo_cos2[i] if final_pqo_cos2 else 0
        ph = final_phase[i]
        marker = " <<" if f2 < 0.1 else ""
        print(f"{zt:10.4f} {f2:10.4f} {pc2:8.4f} {ph:8.4f}{marker}")

    detected = sum(1 for f in final_F2 if abs(f) < 0.1)
    print(f"\nF2 < 0.1 영점: {detected}/{len(active_zeros)}")
    print(f"\n[완료] 출력: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
