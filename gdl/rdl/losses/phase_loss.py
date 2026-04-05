"""
=============================================================================
[Project RDL] Discriminant Loss - F₂(t) = 0 Residual
=============================================================================
판별식 F₂(t) = Im{e^{-iφ}(L̂ - φ̇)} = 0의 제곱 손실 형태.

[구현 수식 매핑 - unified_en.tex 참조]
- eq:F2: F₂(t) = Im{e^{-iφ}(L̂ - φ̇)} = 0 (공명 조건의 허수부)
- L_res = || Im{e^{-iφ}(L_G - ψ·τ')} ||²

기하학적 의미: 잔차 벡터를 위상 -φ만큼 회전시킨 후
궤도에 직교하는 허수부(노이즈) 성분만을 0으로 강제한다.
"""

import torch
import torch.nn as nn
from typing import Optional

from gdl.rdl.constants import PrecisionManager

class DiscriminantLoss(nn.Module):
    """
    판별식 손실: L_res = ||Im{e^{-iφ}(L̂ - φ̇)}||² (eq:F2)

    L̂(= L_G)과 게이지 예측 궤도(ψ·τ')의 직교 어긋남을 측정한다.
    """
    def __init__(self, reduction: str = 'mean'):
        """
        Args:
            reduction (str): 'mean', 'sum', 'none' 중 배치/피처 단위 손실 축소 방식.
        """
        super(DiscriminantLoss, self).__init__()
        if reduction not in ['mean', 'sum', 'none']:
            raise ValueError(f"❌ [PhaseLoss Error] 지원하지 않는 reduction 방식입니다: {reduction}")
        self.reduction = reduction

    def forward(self,
                L_G: torch.Tensor,
                phi: torch.Tensor,
                psi: torch.Tensor,
                tau_prime: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Args:
            L_G (torch.Tensor): [Phase 3] 추출 유효 복소 흐름 속도 (Shape: [Batch, ..., Features], Complex128)
            phi (torch.Tensor): [Phase 4] 게이지 버퍼의 내부 누적 위상 상태 (Shape: [Batch, ..., Features], Float64)
            psi (torch.Tensor): [Phase 4] 게이지 버퍼의 위상 속도 관성 (Shape: [Batch, ..., Features], Float64)
            tau_prime (torch.Tensor, optional): 현재 경로의 진행 궤도 속도 (Shape: [Batch, ...], Float64)

        Returns:
            torch.Tensor: 계산된 위상 정렬 잔차 손실(L_res) 값
        """
        # [방어 기제 1] 정밀도 무결성 강제 점검 (Loss 계산 시 미세 오차 폭주 차단)
        if L_G.dtype != PrecisionManager.COMPLEX_DTYPE:
            L_G = L_G.to(dtype=PrecisionManager.COMPLEX_DTYPE)

        if tau_prime is None:
            tau_prime = torch.ones_like(psi)

        real_tensors = [phi, psi, tau_prime]
        for i, t in enumerate(real_tensors):
            if t.dtype != PrecisionManager.REAL_DTYPE:
                real_tensors[i] = t.to(dtype=PrecisionManager.REAL_DTYPE)
        phi, psi, tau_prime = real_tensors

        # [방어 기제 2] Broadcasting 조율 (tau_prime 차원 맞춤)
        while tau_prime.dim() < psi.dim():
            tau_prime = tau_prime.unsqueeze(-1)

        # -------------------------------------------------------------
        # [Step 1] 예상 잔차(Residual) 벡터 도출
        # -------------------------------------------------------------
        # Diff = L_G - (psi * tau')
        # L_G는 복소수이고 예측치(psi*tau')는 실수이므로,
        # 복소수 뺄셈 시 L_G의 실수부(Real)에서만 값이 차감됩니다.
        expected_flow = psi * tau_prime
        diff_real = L_G.real - expected_flow
        diff_imag = L_G.imag  # 허수부는 원본 그대로 유지됨

        # -------------------------------------------------------------
        # [Step 2] 복소 위상 회전 및 직교 노이즈(Im) O(1) 초고속 추출
        # -------------------------------------------------------------
        # 목표 수식: Im { e^{-iφ} * (diff_real + i * diff_imag) }
        # [극강 최적화] 복소 지수 텐서(e^{-iφ})를 생성하여 복소 행렬곱을 수행하면
        # GPU VRAM 할당 폭주 및 부동소수점 오차가 누적됩니다.
        # 오일러 공식을 통해 대수적으로 전개하여 실수(Real) 연산만으로 허수부를 즉시 빼냅니다.
        # Im( (cos(-φ) + i*sin(-φ)) * (diff_real + i*diff_imag) )
        # = Im( (cos(φ) - i*sin(φ)) * (diff_real + i*diff_imag) )
        # = cos(φ) * diff_imag - sin(φ) * diff_real

        cos_phi = torch.cos(phi)
        sin_phi = torch.sin(phi)

        # 위상 공간을 비틀어(Rotation) 궤도 측면으로 튀어나간 횡방향 직교 노이즈 에너지 도출
        orthogonal_noise = cos_phi * diff_imag - sin_phi * diff_real

        # -------------------------------------------------------------
        # [Step 3] 손실 에너지(L2 Norm Squared) 계산 및 집계
        # -------------------------------------------------------------
        # L_res = || 직교 노이즈 ||^2
        loss_tensor = orthogonal_noise ** 2

        if self.reduction == 'mean':
            return torch.mean(loss_tensor)
        elif self.reduction == 'sum':
            return torch.sum(loss_tensor)
        else:
            return loss_tensor

    def extra_repr(self) -> str:
        return f"eq:F2: L_res = ||Im{{e^(-iφ)(L̂ - φ̇)}}||², reduction='{self.reduction}'"


# 하위 호환 별칭
PhaseAlignmentLoss = DiscriminantLoss


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Losses] Phase Alignment Residual Loss Test ---")

    # 1. 텐서 설정 (Batch=2, Features=3)
    batch_size, out_dim = 2, 3

    # 2. 극한의 정렬/어긋남 시나리오 모사
    # Batch 0: 완벽한 정렬 (L_G = ψτ', 허수부 0) -> 궤도를 정확히 타고 흐름 (Loss 0.0)
    # Batch 1: 극심한 직교 어긋남 (L_G 실수부는 맞으나, 측면 허수부 노이즈가 10.0으로 폭발) -> 거대한 Loss

    phi_mock = torch.zeros((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)
    psi_mock = torch.ones((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)
    tau_p_mock = torch.ones((batch_size, 1), dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    # L_G_aligned: 실수부 1.0 (ψτ'와 일치), 허수부 0.0
    L_G_aligned = torch.complex(
        torch.ones((1, out_dim), dtype=PrecisionManager.REAL_DTYPE),
        torch.zeros((1, out_dim), dtype=PrecisionManager.REAL_DTYPE)
    )

    # L_G_misaligned: 실수부 1.0 (ψτ'와 일치), 허수부 10.0 (강력한 측면 노이즈)
    L_G_misaligned = torch.complex(
        torch.ones((1, out_dim), dtype=PrecisionManager.REAL_DTYPE),
        torch.full((1, out_dim), 10.0, dtype=PrecisionManager.REAL_DTYPE)
    )
    L_G_mock = torch.cat([L_G_aligned, L_G_misaligned], dim=0).requires_grad_(True)

    # 3. 손실 함수 구동 (reduction='none'으로 배치별 오차 독립 확인)
    loss_fn_none = PhaseAlignmentLoss(reduction='none')
    loss_tensor = loss_fn_none(L_G_mock, phi_mock, psi_mock, tau_p_mock)

    print("배치별 위상 정렬 잔차 손실 (L_res) 값:")
    print(f"  └─ Batch 0 (완벽 정렬, 측면 노이즈 0): {loss_tensor[0].mean().item():.8f} (기대값: 0.0)")
    print(f"  └─ Batch 1 (극심한 횡방향 노이즈 10): {loss_tensor[1].mean().item():.4f} (기대값: 100.0 / 10^2)\n")

    is_aligned = torch.allclose(loss_tensor[0], torch.zeros_like(loss_tensor[0]), atol=1e-12)
    is_punished = torch.allclose(loss_tensor[1], torch.full_like(loss_tensor[1], 100.0), atol=1e-5)

    if is_aligned and is_punished:
        print("✅ 성공: 위상 궤도(Real)를 얼마나 빨리/느리게 가는지는 처벌하지 않고, 오직 궤도를 이탈하는 측면 직교 노이즈(Imaginary)만을 완벽히 추출하여 응징합니다.")
    else:
        print("❌ 실패: 기하학적 직교 노이즈 추출에 오류가 있습니다.")

    # 4. 거대 4-Way 역전파(BPTT) 보존 및 미분 독립성 검증
    loss_fn_mean = PhaseAlignmentLoss(reduction='mean')
    total_loss = loss_fn_mean(L_G_mock, phi_mock, psi_mock, tau_p_mock)
    total_loss.backward()

    if all(t.grad is not None for t in [L_G_mock, phi_mock, psi_mock, tau_p_mock]):
        if not any(torch.isnan(t.grad).any() for t in [L_G_mock, phi_mock, psi_mock, tau_p_mock]):
            print("✅ 성공: Loss의 붕괴 압력이 [L_G(전파망), φ(위상), ψ(관성), τ'(진행)] 라는 4개의 다른 차원으로 찢어지며 완벽하게 역전파(Backprop)됩니다!\n")

            # Batch 0의 완벽히 정렬된 데이터로 흐르는 기울기는 정확히 0으로 차단되어야 함
            if torch.allclose(L_G_mock.grad[0], torch.zeros_like(L_G_mock.grad[0])):
                print("✅ 성공: 완벽한 정렬에 도달한 텐서(Batch 0)에는 미분값 0이 흘러 불필요한 가중치 요동(Oscillation)을 기하학적으로 원천 차단합니다!\n")
            else:
                print("❌ 실패: 정답에 도달했음에도 미분 누수(Leakage)가 발생하여 궤도를 어지럽힙니다.\n")
        else:
            print("❌ 실패: 위상 손실 역전파 중 기울기 폭발(NaN) 발생\n")
    else:
        print("❌ 실패: 4-Way 미분 그래프가 끊어졌습니다.\n")
