"""
=============================================================================
[Project RDL] Gauge State Dynamics - Signal Separator
=============================================================================
이 모듈은 Phase 3의 다원군 공명 전파 네트워크에서 도출된
복소 유효 로그기울기(L_G)를 실수부(Real)와 허수부(Imaginary)로 분리합니다.

[구현 수학 식 매핑]
- E28: v_t(λ_j) = Re(L_G(λ_j)) (진행 방향의 평균 위상 속도, Signal)
- E29: v_σ(λ_j) = Im(L_G(λ_j)) (측면 방향의 순간적 위상 노이즈, Noise)

분리된 v_t는 게이지 엔진의 내부 상태(Phase)를 업데이트하는 주 동력원으로 사용되고,
분리된 v_σ는 위상 노이즈의 강도를 측정하여 반응 속도를 조절하는 데 쓰입니다.
"""

import torch
import torch.nn as nn
from typing import Tuple

from gdl.rdl.constants import PrecisionManager

class SignalSeparator(nn.Module):
    """
    복소 유효 로그기울기 텐서 L_G를 두 개의 독립적인 실수(Real) 텐서(v_t, v_sigma)로
    해체(Decoupling)하는 레이어. 내부 학습 파라미터는 없으며, Autograd 미분 그래프를 온전히 유지합니다.
    """
    def __init__(self):
        super(SignalSeparator, self).__init__()

    def forward(self, L_G: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Args:
            L_G (torch.Tensor): Phase 3에서 추출된 유효 복소 궤도 텐서
                                (Shape: [Batch, ..., Features], dtype: Complex128)

        Returns:
            Tuple[torch.Tensor, torch.Tensor]:
                - v_t: 진행 방향 위상 속도 신호 (Re(L_G), dtype: Float64)
                - v_sigma: 노이즈 궤도 이탈률 텐서 (Im(L_G), dtype: Float64)
        """
        # [수학적 방어 기제 1] 복소수 타입 무결성 검증
        # 입력이 실수(Real)일 경우 허수부 노이즈가 강제로 증발하여 게이지 동역학이 정지하므로
        # 반드시 Complex128 타입이어야 합니다.
        if not L_G.is_complex():
            # 실수 텐서가 잘못 들어온 경우, 강제로 허수부를 0으로 갖는 복소 텐서로 승격
            L_G = torch.complex(L_G, torch.zeros_like(L_G))

        if L_G.dtype != PrecisionManager.COMPLEX_DTYPE:
            L_G = L_G.to(dtype=PrecisionManager.COMPLEX_DTYPE)

        # 1. 실수부 추출: v_t (E28)
        # PyTorch의 .real은 원본 텐서의 View를 반환하여 O(1)의 속도를 보장하며
        # 미분 연결(Autograd Graph)을 완벽히 유지합니다.
        v_t = L_G.real

        # 2. 허수부 추출: v_sigma (E29)
        # 이 노이즈는 [Phase 4.2]의 적응형 감쇠율(α_eff)을 결정하는 핵심 지표로 사용됩니다.
        v_sigma = L_G.imag

        return v_t, v_sigma

    def extra_repr(self) -> str:
        return 'Performs: v_t = Re(L_G), v_σ = Im(L_G)'


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Dynamics] Signal Separator Test ---")

    # 1. 테스트 환경 설정 (Batch=2, Seq=3, Features=4)
    batch_size, seq_len, out_dim = 2, 3, 4

    # 2. 가상의 복소 유효 로그기울기 L_G 생성 (Requires_grad=True)
    real_signal = torch.randn((batch_size, seq_len, out_dim), dtype=PrecisionManager.REAL_DTYPE)
    imag_noise  = torch.randn((batch_size, seq_len, out_dim), dtype=PrecisionManager.REAL_DTYPE) * 0.5
    L_G_mock = torch.complex(real_signal, imag_noise).requires_grad_(True)

    print(f"입력 복소 텐서 L_G Shape: {L_G_mock.shape} | Dtype: {L_G_mock.dtype}\n")

    # 3. 분리기 모듈 인스턴스화 및 Forward Pass
    separator = SignalSeparator()
    v_t_out, v_sigma_out = separator(L_G_mock)

    print(f"추출된 전진 속도(v_t) Shape:       {v_t_out.shape} | Dtype: {v_t_out.dtype}")
    print(f"추출된 위상 노이즈(v_sigma) Shape: {v_sigma_out.shape} | Dtype: {v_sigma_out.dtype}\n")

    # 4. 차원 및 타입 무결성 검증
    expected_shape = (batch_size, seq_len, out_dim)
    if (v_t_out.shape == expected_shape) and (v_sigma_out.shape == expected_shape):
        if (v_t_out.dtype == PrecisionManager.REAL_DTYPE) and (v_sigma_out.dtype == PrecisionManager.REAL_DTYPE):
            print("✅ 성공: 복소 텐서가 두 개의 독립적인 실수 물리량 텐서(v_t, v_σ)로 완벽히 해체(Decoupling)되었습니다.")
        else:
            print("❌ 실패: 분리된 텐서의 정밀도가 Float64가 아닙니다.")
    else:
        print("❌ 실패: 차원 분리 중 손실이 발생했습니다.")

    # 5. 역전파(Autograd) 누수 방어(No Leakage) 및 분리 독립성 검증
    # v_t(신호)에만 오차를 발생시켜 역전파했을 때, 원본 L_G의 실수부에만 Gradient가 흐르는지 확인
    loss_t = v_t_out.sum()
    loss_t.backward()

    grad_L_G = L_G_mock.grad

    if grad_L_G is not None:
        if torch.allclose(grad_L_G.real, torch.ones_like(v_t_out)) and torch.allclose(grad_L_G.imag, torch.zeros_like(v_sigma_out)):
            print("✅ 성공: 실수부(v_t) 오차 역전파가 허수부(v_σ)로 누수(Leakage)되지 않고 완벽히 독립 보존됩니다!\n")
        else:
            print("❌ 실패: 역전파 과정에서 실수부와 허수부의 오차가 섞여버렸습니다(Leakage).\n")
    else:
        print("❌ 실패: L_G로의 미분 연결이 끊어졌습니다.\n")
