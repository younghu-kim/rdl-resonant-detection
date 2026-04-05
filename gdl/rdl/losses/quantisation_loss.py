"""
=============================================================================
[Project RDL] Phase Quantisation Loss - sin²(L·arg(Z)/2) 격자 정렬 손실
=============================================================================
위상 양자화 연산자 P(φ;L)를 손실 함수로 적용하여
네트워크 출력 Z의 위상이 격자 {2πk/L}에 정렬되도록 유도한다.

[구현 수식 매핑 - unified_en.tex 참조]
- eq:pqo: P(φ;L) = sin²(Lφ/2)
- 비고 4.1: QROP-Net sin² 손실과의 구조적 대응
- L_pqo = mean(P(arg(Z); L))

격자점에서 L_pqo = 0 (완벽 양자화), 중점에서 L_pqo = 1 (최대 이탈).
"""

import torch
import torch.nn as nn

from gdl.rdl.constants import R_CONST, PrecisionManager
from gdl.rdl.layers.phase_quantisation import PhaseQuantisationOperator, MaslovCorrectedPQO


class PhaseQuantisationLoss(nn.Module):
    """
    L_pqo = mean(P(arg(Z); L))

    mode='sin2' (기본): P = sin²(L·arg(Z)/2) — 격자점 {2πk/L}에서 최소
    mode='cos2' (Maslov 보정): P = cos²(L·arg(Z)/2) — 전환점 {π/2+mπ}에서 최소

    cos2 모드는 실수값 함수(Z(t))의 영점 탐지에 사용한다.
    영점에서 arg(Z) -> π/2 (부호변화 정칙화)이므로 cos²(π/2) = 0,
    F₂ = 0과 PQO = 0이 동시에 달성된다.
    """

    def __init__(self, L: int = None, reduction: str = 'mean', mode: str = 'sin2'):
        """
        Args:
            L: 양자화 레벨. None이면 R_CONST.PQO_L_DEFAULT (=2)
            reduction: 'mean', 'sum', 'none'
            mode: 'sin2' (원래 Gate A) 또는 'cos2' (Maslov 보정)
        """
        super(PhaseQuantisationLoss, self).__init__()
        if mode == 'cos2':
            self.pqo = MaslovCorrectedPQO(L=L)
        else:
            self.pqo = PhaseQuantisationOperator(L=L)
        self.reduction = reduction
        self.mode = mode

    def forward(self, Z_out: torch.Tensor) -> torch.Tensor:
        """
        Args:
            Z_out: 네트워크 출력 복소 텐서 (임의 shape)

        Returns:
            위상 양자화 손실 스칼라 (reduction='mean'/'sum') 또는 텐서 ('none')
        """
        phase = torch.angle(Z_out)
        penalty = self.pqo(phase)

        if self.reduction == 'mean':
            return penalty.mean()
        elif self.reduction == 'sum':
            return penalty.sum()
        return penalty

    def extra_repr(self) -> str:
        return f"L={self.pqo.L}, mode='{self.mode}', reduction='{self.reduction}'"


# =================================================================
# 직접 실행 시 무결성 테스트
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Losses] Phase Quantisation Loss Test ---")

    loss_fn = PhaseQuantisationLoss(L=2)
    print(f"손실 함수: {loss_fn}")

    # 테스트 1: 격자점 위상 -> 손실 0
    Z_lattice = torch.complex(
        torch.tensor([1.0, -1.0, 1.0], dtype=torch.float64),
        torch.tensor([0.0, 0.0, 0.0], dtype=torch.float64)
    )  # arg = 0, pi, 0
    loss_lattice = loss_fn(Z_lattice)
    print(f"\n[Test 1] 격자점 위상(0, pi, 0): L_pqo = {loss_lattice.item():.8f} (기대값: 0.0)")

    # 테스트 2: 중점 위상 -> 손실 1
    Z_mid = torch.complex(
        torch.tensor([0.0, 0.0], dtype=torch.float64),
        torch.tensor([1.0, -1.0], dtype=torch.float64)
    )  # arg = pi/2, -pi/2
    loss_mid = loss_fn(Z_mid)
    print(f"[Test 2] 중점 위상(pi/2, -pi/2): L_pqo = {loss_mid.item():.4f} (기대값: 1.0)")

    # 테스트 3: autograd
    Z_auto = torch.complex(
        torch.tensor([1.0, 0.5], dtype=torch.float64, requires_grad=True),
        torch.tensor([0.3, 0.8], dtype=torch.float64, requires_grad=True)
    )
    loss_auto = loss_fn(Z_auto)
    loss_auto.backward()
    ok3 = Z_auto.real.grad is not None and not torch.isnan(Z_auto.real.grad).any()
    print(f"[Test 3] Autograd: {'OK' if ok3 else 'FAIL'}")

    if abs(loss_lattice.item()) < 1e-12 and abs(loss_mid.item() - 1.0) < 1e-12 and ok3:
        print("\n[OK] PhaseQuantisationLoss 전체 테스트 통과")
    else:
        print("\n[FAIL] 일부 테스트 실패")
