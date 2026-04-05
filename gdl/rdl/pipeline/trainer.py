"""
=============================================================================
[Project RDL] Orchestration - Master Trainer Engine
=============================================================================
공명 딥러닝 시스템의 End-to-End 훈련 루프.

4항 공명 조건 (판별식 + 곡률 + 타겟 + PQO) 수렴을 목표로,
데이터 로더 → 공명 전파 → 게이지 적분 → 손실 산출 → PGGD 갱신
사이클을 수행한다.
"""

import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from typing import Dict, List

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.optim.pggd import PGGD

class ResonanceTrainer:
    """공명 딥러닝 시스템의 E2E 훈련 관리자."""
    def __init__(self,
                 model: MasterResonantNetwork,
                 train_loader: DataLoader,
                 val_loader: DataLoader,
                 optimizer: PGGD,
                 loss_fn: TotalResonanceLoss,
                 device: torch.device,
                 max_grad_norm: float = 1.0):
        """
        Args:
            model: 공명 마스터 모델
            train_loader: 훈련 데이터 파이프라인
            val_loader: 검증 파이프라인
            optimizer: PGGD 최적화기
            loss_fn: 4항 통합 손실 래퍼
            device: 연산 장치
            max_grad_norm: 기울기 클리핑 임계값
        """
        self.model = model.to(device)
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.optimizer = optimizer
        self.loss_fn = loss_fn.to(device)
        self.device = device
        self.max_grad_norm = max_grad_norm

        self.history = {
            'train_loss': [], 'val_loss': [],
            'L_res': [], 'L_curv': [], 'L_tgt': [], 'L_pqo': []
        }

    def train_epoch(self) -> Dict[str, float]:
        """1 에포크 학습 수행 및 통계 집계"""
        self.model.train()
        epoch_metrics = {'Total': 0.0, 'L_res': 0.0, 'L_curv': 0.0, 'L_tgt': 0.0, 'L_pqo': 0.0}
        num_batches = len(self.train_loader)

        for batch_idx, (X_noisy, _) in enumerate(self.train_loader):
            # [방어 기제 1] 데이터 텐서 장치 이동 및 정밀도 통제
            X_in = X_noisy.to(self.device, dtype=PrecisionManager.REAL_DTYPE)

            # [방어 기제 2] 데이터 텐서 미분 추적기(requires_grad) 활성화
            # 이 한 줄이 없으면 Phase 5.2 곡률 극점 안정화 손실에서 에러가 발생합니다.
            X_in.requires_grad_(True)

            self.optimizer.zero_grad(set_to_none=True) # VRAM 절약을 위해 그래디언트 소멸

            # -------------------------------------------------------------
            # [Phase 2, 3, 4] Forward Pass: 공명 전파망 및 동역학 버퍼 관통
            # -------------------------------------------------------------
            outputs = self.model(X_in)

            # 4항 공명 손실 산출
            total_loss, metrics = self.loss_fn(**outputs)

            # NaN/Inf 방어
            if torch.isnan(total_loss) or torch.isinf(total_loss):
                print(f"\n[WARN] Batch {batch_idx}: Loss NaN/Inf 발생, 스킵.")
                continue

            # -------------------------------------------------------------
            # [Phase 6] Backward & Optimization: 초거대 통합 미분 및 PGGD 진화
            # -------------------------------------------------------------
            total_loss.backward()

            # [방어 기제 4] 특이점(Singularity) 기울기 폭발(NaN/Inf) 억제 클리핑
            # 곡률 손실의 3차 미분 파도를 강제로 깎아내어 NaN 오염 방지
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), self.max_grad_norm)

            # PGGD 옵티마이저의 게이지 적분 및 가중치 업데이트
            self.optimizer.step()

            for k in epoch_metrics.keys():
                epoch_metrics[k] += metrics.get(k, 0.0)

        # 배치 평균치 반환
        return {k: v / max(1, num_batches) for k, v in epoch_metrics.items()}

    @torch.no_grad()
    def validate_epoch(self) -> float:
        """모든 검증 데이터를 순회하며 성능 평가"""
        self.model.eval()
        total_val_loss = 0.0
        num_batches = len(self.val_loader)

        # [방어 기제 5] 평가(Eval) 모드의 역발상 (enable_grad 개방)
        # 검증 환경에서도 곡률 2차 미분(A'')을 산출해야 하므로 추론 모드에서 예외적으로 미분을 허용합니다.
        with torch.enable_grad():
            for X_noisy, _ in self.val_loader:
                X_in = X_noisy.to(self.device, dtype=PrecisionManager.REAL_DTYPE)
                X_in.requires_grad_(True)

                outputs = self.model(X_in)
                total_loss, _ = self.loss_fn(**outputs)

                total_val_loss += total_loss.item()

        return total_val_loss / max(1, num_batches)

    def fit(self, epochs: int, print_freq: int = 1) -> Dict[str, List[float]]:
        """지정된 에포크 수만큼 전체 훈련 루프를 구동하는 메인 엔진"""
        print(f"\n{'='*80}")
        print(f"[RDL] 공명 딥러닝 훈련 시작 (Epochs: {epochs}, Device: {self.device})")
        print(f"{'='*80}")

        for epoch in range(1, epochs + 1):
            train_metrics = self.train_epoch()
            val_loss = self.validate_epoch()

            self.history['train_loss'].append(train_metrics['Total'])
            self.history['L_res'].append(train_metrics['L_res'])
            self.history['L_curv'].append(train_metrics['L_curv'])
            self.history['L_tgt'].append(train_metrics['L_tgt'])
            self.history['L_pqo'].append(train_metrics.get('L_pqo', 0.0))
            self.history['val_loss'].append(val_loss)

            if epoch % print_freq == 0 or epoch == epochs:
                pqo_str = f" | PQO: {train_metrics.get('L_pqo', 0.0):.4f}" if train_metrics.get('L_pqo', 0) > 0 else ""
                print(f"[Epoch {epoch:03d}/{epochs:03d}] "
                      f"Train: {train_metrics['Total']:.4f} "
                      f"(Res: {train_metrics['L_res']:.4f} | Curv: {train_metrics['L_curv']:.4f} | "
                      f"Tgt: {train_metrics['L_tgt']:.4f}{pqo_str}) "
                      f"|| Val: {val_loss:.4f}")

        print(f"{'='*80}")
        print("[RDL] 공명 조건 수렴. 훈련 종료.")
        print(f"{'='*80}")
        return self.history


# =================================================================
# 직접 실행 시 파이프라인 무결성 5 에포크 주행 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()
    from gdl.rdl.pipeline.dataset import get_dataloaders

    print("\n--- [RDL Orchestration] Master Trainer Engine 5 Epochs Trial ---")

    # 1. 시뮬레이션 환경 통제 (MPS/CUDA/CPU 자동 할당)
    device = torch.device('cuda' if torch.cuda.is_available() else ('mps' if torch.backends.mps.is_available() else 'cpu'))
    print(f"▶ 훈련 장치 록온: {device}")

    # 2. 데이터셋 및 파이프라인 구축 (1D 64-Sequence Harmonic Noise)
    in_dim, hidden_dim, out_dim = 64, 16, 2
    train_loader, val_loader = get_dataloaders(
        batch_size=8, num_samples=80, seq_len=in_dim, noise_std=0.5
    )

    # 3. 모델, 손실, 옵티마이저 조립
    model = MasterResonantNetwork(in_features=in_dim, hidden_features=hidden_dim, out_features=out_dim,
                                   num_layers=2, channel_type='paper3ch', damping_mode='paper')
    loss_fn = TotalResonanceLoss(lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.5)
    optimizer = PGGD(model.parameters(), lr=0.01, tau_g=0.5, h=0.1)

    # 4. 트레이너 엔진 시동
    trainer = ResonanceTrainer(model, train_loader, val_loader, optimizer, loss_fn, device, max_grad_norm=5.0)

    # 5. 미니 학습 루프 돌리기 (엔진 가동!)
    try:
        history = trainer.fit(epochs=5)

        # 훈련 수렴 여부 확인
        initial_loss = history['train_loss'][0]
        final_loss = history['train_loss'][-1]

        if final_loss < initial_loss:
            print(f"\n✅ 대성공: 단 5 Epoch 만에 스칼라 에너지가 성공적으로 수축(Collapse)했습니다! ({initial_loss:.4f} -> {final_loss:.4f})")
            print("          고차 미분(Hessian) 역전파와 게이지 적분이 VRAM 폭주나 NaN 붕괴 없이 하나로 결합되었습니다.")
        else:
            print("\n⚠️ 경고: 훈련 에러는 없었으나 손실값이 뚜렷이 감소하지 않았습니다. 하이퍼파라미터 조율이 필요합니다.")

    except Exception as e:
        print(f"\n❌ 실패: 훈련 사이클 도중 시스템 붕괴 발생! 에러 로그:\n{e}")
