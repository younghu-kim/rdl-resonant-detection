"""
=============================================================================
[Project RDL] Target Phase Matching Loss
=============================================================================
네트워크의 최종 복소 출력 Z^{(L)}의 위상각(arg Z)이
사전에 정의된 이상 위상 타겟 곡선 Ψ(X)와 일치하도록 강제하는 목적 함수.

[구현 수식 매핑]
- 단순 수식: L_match = || arg(Z^{(L)}) - Ψ(X) ||^2
- 방어 수식: L_match = || atan2(sin(Δφ), cos(Δφ)) ||^2
  (2π 주기를 넘나들 때 오차가 폭발하는 현상을 막기 위해,
   삼각함수 투영을 통한 '최단 각도 거리'로 정밀한 MSE를 계산합니다.)
"""

import math
import torch
import torch.nn as nn

from gdl.rdl.constants import PrecisionManager

class TargetPhaseMatchingLoss(nn.Module):
    def __init__(self, reduction: str = 'mean'):
        super(TargetPhaseMatchingLoss, self).__init__()
        if reduction not in ['mean', 'sum', 'none']:
            raise ValueError(f"❌ [TargetLoss Error] 지원하지 않는 reduction 방식입니다: {reduction}")
        self.reduction = reduction

    def forward(self, Z_out: torch.Tensor, Psi_target: torch.Tensor) -> torch.Tensor:
        """
        Args:
            Z_out (torch.Tensor): 네트워크 최종 출력 복소 텐서 (Shape: [Batch, ...])
            Psi_target (torch.Tensor): [Phase 2]에서 생성된 이상 위상 타겟 (Shape: [Batch, ...])

        Returns:
            torch.Tensor: 위상 래핑(Wrapping) 모순이 제거된 최소 각도 거리 제곱 오차(MSE)
        """
        # [방어 기제 1] 정밀도 무결성 확인 및 타입 승격
        if Z_out.dtype != PrecisionManager.COMPLEX_DTYPE:
            Z_out = Z_out.to(dtype=PrecisionManager.COMPLEX_DTYPE)
        if Psi_target.dtype != PrecisionManager.REAL_DTYPE:
            Psi_target = Psi_target.to(dtype=PrecisionManager.REAL_DTYPE)

        # 1. 네트워크 출력의 현재 위상각 추출: arg(Z) ∈ [-π, π]
        # torch.angle() 은 내부적으로 atan2(imag, real)를 호출하며 C^∞ 미분을 보장합니다.
        current_phase = torch.angle(Z_out)

        # Broadcasting 호환성을 위한 차원 정렬
        while Psi_target.dim() < current_phase.dim():
            Psi_target = Psi_target.unsqueeze(-1)

        # 2. 위상각 단순 차이 도출
        phase_diff = current_phase - Psi_target

        # 3. [방어 기제 2] 위상 래핑(Phase Wrapping) 모순 해결 (-π ~ π 구간 락킹)
        # 예: current=3.14 (π), target=-3.14 (-π) 일 때 단순 차이는 6.28 이지만,
        # 원 위에서는 같은 위치이므로 실제 차이는 0.0 이어야 합니다.
        shortest_diff = torch.atan2(torch.sin(phase_diff), torch.cos(phase_diff))

        # 4. 스케일링 된 최소 각도 거리 제곱 (MSE)
        loss_tensor = shortest_diff ** 2

        if self.reduction == 'mean':
            return torch.mean(loss_tensor)
        elif self.reduction == 'sum':
            return torch.sum(loss_tensor)
        else:
            return loss_tensor

    def extra_repr(self) -> str:
        return f"Equation: L_match = || min_dist(arg(Z), Ψ) ||^2, reduction='{self.reduction}'"


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()
    print("\n--- [RDL Losses] Target Phase Matching Loss Test ---")

    batch_size, out_dim = 2, 3

    # 1. 테스트 시나리오 구성
    # Batch 0: 일반적인 오차 상황 (현재: 0.5 rad, 타겟: 1.0 rad) -> 차이 0.5, 제곱 0.25
    # Batch 1: 극단적 위상 래핑 모순 상황 (현재: π, 타겟: -π)
    # -> 물리적으로 동일 위치이므로 단순 뺄셈(6.28^2)이 아닌 기하학적 최소 거리(0.0)가 도출되어야 함.

    Z_out_b0 = torch.complex(
        torch.cos(torch.tensor([0.5]*out_dim, dtype=PrecisionManager.REAL_DTYPE)),
        torch.sin(torch.tensor([0.5]*out_dim, dtype=PrecisionManager.REAL_DTYPE))
    )
    Z_out_b1 = torch.complex(
        torch.cos(torch.tensor([math.pi]*out_dim, dtype=PrecisionManager.REAL_DTYPE)),
        torch.sin(torch.tensor([math.pi]*out_dim, dtype=PrecisionManager.REAL_DTYPE))
    )
    Z_out = torch.stack([Z_out_b0, Z_out_b1], dim=0).requires_grad_(True)

    Psi_target = torch.tensor([
        [1.0]*out_dim,
        [-math.pi]*out_dim
    ], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    # 2. 매칭 손실 함수 구동
    loss_fn = TargetPhaseMatchingLoss(reduction='none')
    loss_tensor = loss_fn(Z_out, Psi_target)

    print("▶ 배치별 목표 위상 매칭 손실 (L_match) 값:")
    print(f"  └─ Batch 0 (일반 오차 0.5): {loss_tensor[0].mean().item():.4f} (기대값: 0.2500)")
    print(f"  └─ Batch 1 (π와 -π 래핑):   {loss_tensor[1].mean().item():.8f} (기대값: 0.0000, 래핑 방어 성공)\n")

    if loss_tensor[1].mean().item() < 1e-10:
        print("✅ 성공: 원 위에서 동일한 위치인 π와 -π의 차이를 0.0으로 완벽히 인식하여 래핑(Wrapping) 모순을 100% 방어했습니다.")
    else:
        print(f"❌ 실패: π와 -π를 다른 위치로 오판하여 {loss_tensor[1].mean().item():.2f} 의 거짓(False) 페널티를 부여했습니다.")

    # 3. 역전파 그래프 보존 검증
    loss_fn_mean = TargetPhaseMatchingLoss(reduction='mean')
    total_loss = loss_fn_mean(Z_out, Psi_target)
    total_loss.backward()

    if Z_out.grad is not None and not torch.isnan(Z_out.grad).any():
        print("✅ 성공: atan2(sin, cos) 기반의 비선형 최소 거리 계산이 단절 없이 C^∞ 역전파 트리를 완벽히 보존했습니다!\n")
    else:
        print("❌ 실패: 위상 래핑 우회 중 미분 꺾임(NaN)이 발생했습니다.")
