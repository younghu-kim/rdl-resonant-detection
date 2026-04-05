"""
=============================================================================
[Project RDL] Phase Space Embedding - Ideal Phase Target Generator
=============================================================================
이 모듈은 실수 입력 데이터 X의 스케일(||X||)을 바탕으로,
네트워크가 도달해야 할 궁극적인 '이상 위상(Ideal Phase)' 타겟 Ψ(X)를 생성합니다.
고전 수론의 정수 'n'을 텐서 데이터의 Norm(||X||)으로 연속화하여 대체합니다.

[구현 수학 식]
Ψ(X) = α * log(||X||) + β * ( sin(γ * log(||X||)) / ||X||^δ )
"""

import torch
import torch.nn as nn

from gdl.rdl.constants import R_CONST, PrecisionManager

class IdealPhaseTargetGenerator(nn.Module):
    """
    입력 텐서 X로부터 이상 위상 목표값 Ψ(X)를 계산하는 레이어.
    학습 파라미터 없이 R_CONST의 상수만으로 궤도를 생성합니다 (eq:psi_n).
    """
    def __init__(self):
        super(IdealPhaseTargetGenerator, self).__init__()

        # 상수 매핑 (수학적 불변성 유지)
        self.alpha = R_CONST.ALPHA
        self.beta  = R_CONST.BETA
        self.gamma = R_CONST.GAMMA
        self.delta = R_CONST.DELTA
        self.eps   = R_CONST.EPSILON

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x (torch.Tensor): 원본 실수 입력 텐서 (Shape: [Batch, ...])

        Returns:
            torch.Tensor: 계산된 이상 위상 타겟 Ψ(X) (Shape: [Batch, 1], dtype: Float64)
        """
        # [방어 기제 1] Float64 정밀도 강제
        if x.dtype != PrecisionManager.REAL_DTYPE:
            x = x.to(dtype=PrecisionManager.REAL_DTYPE)

        # 1. 배치 단위 데이터 크기(Norm) ||X|| 추출
        # Batch 차원(dim=0)을 유지하고 나머지 모든 차원을 묶어서 L2 Norm 계산을 위한 제곱합 도출
        x_flat = x.view(x.size(0), -1)
        x_sq_sum = torch.sum(x_flat**2, dim=1, keepdim=True)

        # [방어 기제 2] 영벡터(0)에서 미분 시 NaN 붕괴를 막기 위한 Epsilon-safe Norm (C^∞ 보장)
        # torch.linalg.norm이나 단순 clamp를 사용하면 x=0일 때 미분 그래프가 파괴되므로,
        # 제곱합에 eps를 직접 더한 뒤 루트를 씌워 완벽한 미분 가능 곡면을 만듭니다.
        safe_norm = torch.sqrt(x_sq_sum + self.eps)

        # 2. 공통 연산: log(||X||)
        log_norm = torch.log(safe_norm)

        # 3. 주성분 (단조 상승 추세): α * log(||X||)
        term1 = self.alpha * log_norm

        # 4. 진동 및 감쇠 항: β * sin(γ * log(||X||)) / ||X||^δ
        oscillation = torch.sin(self.gamma * log_norm)
        damping = torch.pow(safe_norm, self.delta)
        term2 = self.beta * (oscillation / damping)

        # 5. 최종 이상 위상 조립
        psi_x = term1 + term2

        return psi_x

    def extra_repr(self) -> str:
        return f'α={self.alpha}, β={self.beta}, γ={self.gamma}, δ={self.delta}, eps={self.eps:.1e}'

# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Layers] Ideal Phase Target Generator Test ---")

    # 1. 0벡터(특이점)를 포함하여 크기가 각기 다른 데이터 텐서 생성 (Batch=3)
    x_input = torch.tensor([
        [0.0, 0.0, 0.0],       # Singularity (특이점)
        [1.0, 2.0, 3.0],       # 일반 데이터
        [10.0, 20.0, 30.0]     # 대규모 스케일
    ], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    print(f"입력 데이터 X:\n{x_input.detach().cpu().numpy()}\n")

    # 2. Generator 인스턴스화 및 Forward Pass
    psi_gen = IdealPhaseTargetGenerator()
    psi_out = psi_gen(x_input)

    print(f"생성된 이상 위상 타겟 Ψ(X) (Shape: {psi_out.shape}):\n{psi_out.detach().cpu().numpy()}\n")

    # 3. NaN / Inf 발생 여부 검증
    if torch.isnan(psi_out).any() or torch.isinf(psi_out).any():
        print("❌ 실패: Ψ(X) 계산 중 NaN 또는 Inf가 발생했습니다 (특이점 붕괴).")
    else:
        print("✅ 성공: 영벡터(Zero-vector) 입력 시에도 특이점을 안전하게 우회하여 Ψ(X)가 계산되었습니다.")

    # 4. 미분 그래프 보존(Autograd) 확인
    loss = psi_out.sum()
    loss.backward()

    if x_input.grad is not None and not torch.isnan(x_input.grad).any():
        print("✅ 성공: 이상 위상 역전파(Backpropagation) 시 미분 그래프가 완벽히 보존되었습니다.\n")
    else:
        print("❌ 실패: 미분 중 단절 또는 기울기 폭발(NaN)이 발생했습니다.\n")
