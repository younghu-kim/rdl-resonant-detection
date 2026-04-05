"""
=============================================================================
[Project RDL] Multigroup Resonant Forward - Effective Log-Derivative Assembler
=============================================================================
이 모듈은 다원군 프레임에서 K개의 방향으로 분할되어 추출된 '방향성 로그 미분(L_k)'을
최소자승법으로 구해진 '기하학적 가중치(w_k)'와 결합(Assembly)하여,
최종적인 단일 '유효 로그기울기 텐서(L_G)'로 수축(Collapse)시킵니다.

[구현 수학 식 매핑]
- E27: L_G(λ_j) = Σ_k w_{j,k} * L_k(λ_j)
- E28: v(λ_j) = Re(L_G(λ_j)) (Phase 4 게이지 동역학으로 넘길 실수부 추출의 기반)
"""

import torch
import torch.nn as nn

from gdl.rdl.constants import PrecisionManager

class EffectiveLogDerivativeAssembler(nn.Module):
    """
    K개의 방향으로 추출된 로그 미분 텐서 L_k와 최적화된 가중치 w_k를
    내적(Dot Product)하여 최종 유효 로그기울기 L_G를 산출하는 모듈.
    학습 가능한 파라미터(nn.Parameter)가 없는 순수 수학적 결합(Reduction) 연산자입니다.
    """
    def __init__(self):
        super(EffectiveLogDerivativeAssembler, self).__init__()

    def forward(self, L_k: torch.Tensor, w_k: torch.Tensor) -> torch.Tensor:
        """
        Args:
            L_k (torch.Tensor): [Block 3.1]에서 추출된 방향성 로그 미분 텐서.
                                (Shape: [Batch, ..., K, Features], dtype: Complex128)
            w_k (torch.Tensor): [Block 1.4]에서 계산된 다원군 기하학적 가중치.
                                (Shape: [Batch, ..., K], dtype: Float64)

        Returns:
            torch.Tensor: 조립된 유효 로그기울기 텐서 L_G
                          (Shape: [Batch, ..., Features], dtype: Complex128)
        """
        # [수학적 방어 기제 1] 데이터 타입 강제 확인 (무결성 유지)
        if L_k.dtype != PrecisionManager.COMPLEX_DTYPE:
            L_k = L_k.to(dtype=PrecisionManager.COMPLEX_DTYPE)
        if w_k.dtype != PrecisionManager.REAL_DTYPE:
            w_k = w_k.to(dtype=PrecisionManager.REAL_DTYPE)

        # [수학적 방어 기제 2] 채널 수(K) 차원 정합성 검증
        # w_k의 마지막 차원(K)과 L_k의 뒤에서 두 번째 차원(K)이 일치해야 합니다.
        K_L = L_k.size(-2)
        K_w = w_k.size(-1)
        if K_L != K_w:
            raise ValueError(f"❌ [Assembler Error] 채널 수(K) 불일치: L_k의 K={K_L}, w_k의 K={K_w}")

        # 1. 가중치 텐서 Broadcasting 차원 확장
        # w_k: [Batch, ..., K] -> [Batch, ..., K, 1]
        # L_k의 Feature 차원에 스칼라 곱셈이 Broadcasting 될 수 있도록 축을 추가합니다.
        w_expanded = w_k.unsqueeze(-1)

        # 2. 유효 로그기울기 결합 (수식 E27: Σ w_k * L_k)
        # L_k (Complex) 와 w_expanded (Real) 의 요소별 곱셈(Element-wise)을 수행한 뒤,
        # K 차원(dim=-2)에 대해 합산(Sum)하여 다원군 채널 차원을 소멸(Collapse)시킵니다.
        L_G = torch.sum(w_expanded * L_k, dim=-2)

        # Output Shape: [Batch, ..., Features]
        return L_G

    def extra_repr(self) -> str:
        return 'Performs: L_G = Σ (w_k * L_k)'


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()
    from gdl.rdl.core.multigroup import MultigroupFrame
    from gdl.rdl.core.ridge_solver import RidgeWeightSolver

    print("\n--- [RDL Layers] Effective Log-Derivative Assembler Test ---")

    # 1. 환경 설정 (Batch=2, Seq=3, K=5(5ch), Features=8)
    batch_size, seq_len, K_channels, out_dim = 2, 3, 5, 8

    # 2. Block 3.1의 결과물을 모사하는 L_k (복소수)
    real_part = torch.randn((batch_size, seq_len, K_channels, out_dim), dtype=PrecisionManager.REAL_DTYPE)
    imag_part = torch.randn((batch_size, seq_len, K_channels, out_dim), dtype=PrecisionManager.REAL_DTYPE)
    L_k_mock = torch.complex(real_part, imag_part).requires_grad_(True)

    # 3. Block 1.4 릿지 솔버를 이용한 '진짜' 가중치 w_k 생성
    frame_5ch = MultigroupFrame(channel_type='5ch')
    solver = RidgeWeightSolver(frame_5ch)

    # 임의의 경로 접선 방향 (σ', τ') 생성 및 가중치 도출
    sigma_p = torch.randn((batch_size, seq_len), dtype=PrecisionManager.REAL_DTYPE) * 0.1
    tau_p   = torch.ones((batch_size, seq_len), dtype=PrecisionManager.REAL_DTYPE)
    w_k_mock = solver.solve(sigma_p, tau_p)  # Shape: [2, 3, 5]
    w_k_mock.requires_grad_(True)

    print(f"입력 L_k Shape: {L_k_mock.shape} | Dtype: {L_k_mock.dtype}")
    print(f"입력 w_k Shape:  {w_k_mock.shape}    | Dtype: {w_k_mock.dtype}\n")

    # 4. 어셈블러 인스턴스화 및 Forward Pass
    assembler = EffectiveLogDerivativeAssembler()
    L_G_out = assembler(L_k_mock, w_k_mock)

    print(f"조립된 유효 로그기울기 L_G Shape: {L_G_out.shape} (기대값: [{batch_size}, {seq_len}, {out_dim}])")

    # 5. 차원 붕괴(Collapse) 검증
    expected_shape = (batch_size, seq_len, out_dim)
    if L_G_out.shape == expected_shape:
        print("✅ 성공: K개의 다원군 채널이 완벽하게 하나의 유효 궤도(L_G)로 융합(Collapse)되었습니다.")
    else:
        print(f"❌ 실패: 출력 차원에 오류가 있습니다. 실제값: {L_G_out.shape}")

    # 6. 미분 역전파(Autograd) 다중 경로 보존 검증
    # L_k와 w_k 두 갈래로 나뉘었던 미분 트리가 L_G에서 다시 만나 역방향으로 흐르는지 확인
    loss = L_G_out.abs().sum()
    loss.backward()

    if L_k_mock.grad is not None and w_k_mock.grad is not None:
        if not torch.isnan(L_k_mock.grad).any() and not torch.isnan(w_k_mock.grad).any():
            print("✅ 성공: 복소 텐서(L_k)와 실수 텐서(w_k)의 혼합 내적 후에도 미분 그래프가 단절 없이 100% 보존되었습니다!\n")
        else:
            print("❌ 실패: 역전파 중 NaN이 발생했습니다.\n")
    else:
        print("❌ 실패: L_k 또는 w_k로의 역전파 연결이 끊어졌습니다.\n")
