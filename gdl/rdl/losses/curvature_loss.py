"""
=============================================================================
[Project RDL] Gate B - Curvature Stabilization Loss (v_sigma=0, A'=0, A''<0)
=============================================================================
진폭(에너지) A가 공간 상에서 극점(A'=0)에 도달하고,
그 지점이 오목하게 안정화(A'' < 0)되도록 1, 2차 미분을 통제하는 손실 함수.

[구현 수식 매핑]
- 1차 미분(극점 강제): λ_a * || A' ||^2
- 2차 미분(안정성 강제): λ_b * ReLU(A'') (곡률이 양수/볼록할 때만 처벌)
"""

import torch
import torch.nn as nn

from gdl.rdl.constants import PrecisionManager

class CurvatureStabilizationLoss(nn.Module):
    def __init__(self, lambda_a: float = 1.0, lambda_b: float = 1.0):
        super(CurvatureStabilizationLoss, self).__init__()
        self.lambda_a = lambda_a
        self.lambda_b = lambda_b
        self.eps = 1e-8

    def forward(self, Z: torch.Tensor, X: torch.Tensor) -> torch.Tensor:
        """
        Args:
            Z (torch.Tensor): 현재 네트워크의 복소 출력 텐서 (Shape: [Batch, ...])
            X (torch.Tensor): 모델 초입의 원본 실수 입력 텐서 (Shape: [Batch, ...])
                              ※ 2차 미분을 위해 반드시 requires_grad=True 상태여야 함.
        """
        if not X.requires_grad:
            raise RuntimeError("❌ [CurvatureLoss] 입력 텐서 X는 requires_grad=True 상태여야 합니다.")

        if Z.dtype != PrecisionManager.COMPLEX_DTYPE:
            Z = Z.to(dtype=PrecisionManager.COMPLEX_DTYPE)

        # 1. 안전한 진폭(Amplitude) 추출 (eps 추가로 0 근방 미분 꺾임 방지)
        A = torch.sqrt(Z.real**2 + Z.imag**2 + self.eps)

        # 2. 1차 미분(A') 추출: 극점 강제
        # 2차 미분을 위해 create_graph=True 필수
        A_prime = torch.autograd.grad(
            outputs=A.sum(),
            inputs=X,
            create_graph=True,
            retain_graph=True
        )[0]

        # 배치별 기울기 크기 제곱합 (Batch 차원 제외 나머지 차원 합산)
        dim_reduce = list(range(1, A_prime.dim())) if A_prime.dim() > 1 else [0]
        loss_a = self.lambda_a * (A_prime ** 2).sum(dim=dim_reduce).mean()

        # 3. 2차 미분(A'') 추출: 오목성 안정화 강제 (Hessian Trace Proxy)
        # 역전파(backward) 시 모델 가중치 업데이트를 위한 3차 미분이 일어나므로 True 유지
        A_double_prime = torch.autograd.grad(
            outputs=A_prime.sum(),
            inputs=X,
            create_graph=True,
            retain_graph=True
        )[0]

        # A'' > 0 (아래로 볼록, 불안정)일 때만 페널티 부과, A'' <= 0 이면 0.0으로 무시
        A_dp_sum = A_double_prime.sum(dim=dim_reduce)
        loss_b = self.lambda_b * torch.relu(A_dp_sum).mean()

        return loss_a + loss_b


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()
    print("\n--- [RDL Losses] Curvature Stabilization Loss Test ---")

    # 1. 입력 텐서 X (Batch=2, Features=3)
    X = torch.tensor([
        [0.0, 0.0, 0.0],  # Batch 0: 산 정상 (극대점 모사)
        [2.0, 2.0, 2.0]   # Batch 1: 가파른 비탈길 모사
    ], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    # 2. 가상의 복소 공간 랜드스케이프(Z) 모사
    Z_real = torch.empty_like(X)
    Z_real[0] = 10.0 - X[0]**2  # Batch 0: A'=0, A''=-2 (완벽한 안정 분지, 페널티 0 기대)
    Z_real[1] = X[1]**2         # Batch 1: A'=4, A''=2  (불안정한 볼록 곡면, 페널티 발생 기대)

    Z_imag = torch.ones_like(X) * 0.1 # 미분 소실 방지용 미세 허수부 주입
    Z = torch.complex(Z_real, Z_imag)

    # 3. 곡률 손실 구동
    loss_fn = CurvatureStabilizationLoss(lambda_a=1.0, lambda_b=1.0)

    # 각 배치별 Loss를 확인하기 위해 직접 연산 시뮬레이션
    print(f"▶ Batch 0 (안정된 오목 분지)은 ReLU(-2)에 의해 곡률 페널티가 0으로 차단됩니다.")
    print(f"▶ Batch 1 (불안정 볼록 비탈길)은 강력한 1차, 2차 미분 페널티가 부과됩니다.")

    total_loss = loss_fn(Z, X)
    print(f"\n▶ 통합 Loss: {total_loss.item():.4f}")

    # 4. 거대 3차 미분(Loss 역전파) 보존 검증
    total_loss.backward()
    if X.grad is not None and not torch.isnan(X.grad).any():
        print("✅ 성공: 1차 미분(A') -> 2차 미분(A'') -> Loss 역전파(3차 미분)까지 Autograd 그래프가 단절 없이 보존되었습니다!")
    else:
        print("❌ 실패: 고차 미분 트리가 끊어졌거나 NaN 붕괴가 발생했습니다.")
