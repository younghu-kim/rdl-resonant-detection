"""
=============================================================================
[Project RDL] Phase Quantisation Operator - P(phi; L) = sin^2(L * phi / 2)
=============================================================================
QROP-Net과 공명소수법칙을 잇는 다리 방정식(Bridge Equation).

[구현 수식 매핑 - unified_en.tex 참조]
- eq:pqo: P(phi; L) = sin^2(L * phi / 2)
- 정의 3.1: 위상 양자화 연산자
- 비고 4.1: QROP-Net sin^2 손실과의 구조적 대응

[양자화 격자]
- P = 0  ⟺  phi = 2*pi*k/L  (k 정수) — 격자점 위
- P = 1  ⟺  phi = 2*pi*(k+0.5)/L — 격자 중점 (최대 불확실성)

[L 값별 의미]
- L=2: Gate A (Δv = m*pi) — 공명소수법칙 임계선 조건
- L=2048: QROP-Net u-계수 양자화 (Kyber 격자)
- L=8: QROP-Net v-계수 양자화
"""

import torch
import torch.nn as nn

from gdl.rdl.constants import R_CONST, PrecisionManager


class PhaseQuantisationOperator(nn.Module):
    """
    위상 양자화 연산자: P(phi; L) = sin^2(L * phi / 2)

    논문 참조: unified_en.tex, 정의 3.1, eq:pqo

    이 연산자는 위상 값이 격자 {2*pi*k/L}에 얼마나 가까운지를 측정한다.
    격자점에서 P=0 (완벽 양자화), 격자 중점에서 P=1 (최대 이탈).

    기울기 ∇P = (L/2) * sin(L * phi)는 가장 가까운 격자점을 향해 조향한다.
    """

    def __init__(self, L: int = None):
        """
        Args:
            L (int): 양자화 레벨. None이면 R_CONST.PQO_L_DEFAULT (=2) 사용.
                     L=2: Gate A (공명 조건)
                     L=2048: QROP-Net u-계수
                     L=8: QROP-Net v-계수
        """
        super().__init__()
        self.L = L if L is not None else R_CONST.PQO_L_DEFAULT

    def forward(self, phi: torch.Tensor) -> torch.Tensor:
        """
        P(phi; L) = sin^2(L * phi / 2)

        Args:
            phi: 위상 텐서 (임의 shape)

        Returns:
            동일 shape의 양자화 페널티 텐서, 범위 [0, 1]
        """
        return torch.sin(self.L * phi / 2.0) ** 2

    def gradient(self, phi: torch.Tensor) -> torch.Tensor:
        """
        ∇P = (L/2) * sin(L * phi)

        가장 가까운 격자점 방향의 기울기를 명시적으로 반환한다.
        autograd와 독립적으로 사용 가능 (예: PGGD 위상 조향).

        Args:
            phi: 위상 텐서 (임의 shape)

        Returns:
            동일 shape의 기울기 텐서
        """
        return (self.L / 2.0) * torch.sin(self.L * phi)

    def lattice_distance(self, phi: torch.Tensor) -> torch.Tensor:
        """
        가장 가까운 격자점까지의 위상 거리를 반환한다.

        격자점: phi_k = 2*pi*k/L
        거리: min_k |phi - phi_k|

        Args:
            phi: 위상 텐서

        Returns:
            격자 거리 텐서, 범위 [0, pi/L]
        """
        spacing = 2.0 * torch.pi / self.L
        # phi를 [0, spacing) 구간으로 래핑
        remainder = torch.remainder(phi, spacing)
        # 가장 가까운 격자점까지의 거리 (0 또는 spacing)
        return torch.minimum(remainder, spacing - remainder)

    def __repr__(self):
        return f"<PhaseQuantisationOperator: L={self.L}, lattice_spacing={2*3.141592653589793/self.L:.4f}>"


class MaslovCorrectedPQO(nn.Module):
    """
    Maslov 보정 위상 양자화 연산자: P_cos(phi; L) = cos^2(L * phi / 2)

    실수값 함수(Hardy Z-함수 등)의 영점에서 arg(Z) -> pi/2 (부호변화 정칙화).
    sin^2(pi/2) = 1 이므로 원래 PQO는 영점에서 최대가 되는 반면,
    cos^2(pi/2) = 0 이므로 이 연산자는 영점에서 최소가 된다.

    이는 준고전 근사의 Maslov 지수 보정 (mu=2, mu*pi/4 = pi/2)에 대응한다.
    Gate A와 판별식 F2 = 0 조건이 동시에 만족 가능해진다.

    격자 구조:
    - P_cos = 0 at phi = pi/2 + m*pi (영점/전환점)
    - P_cos = 1 at phi = m*pi (비영점)
    """

    def __init__(self, L: int = None):
        super().__init__()
        self.L = L if L is not None else R_CONST.PQO_L_DEFAULT

    def forward(self, phi: torch.Tensor) -> torch.Tensor:
        """P_cos(phi; L) = cos^2(L * phi / 2)"""
        return torch.cos(self.L * phi / 2.0) ** 2

    def gradient(self, phi: torch.Tensor) -> torch.Tensor:
        """dP_cos/dphi = -(L/2) * sin(L * phi)"""
        return -(self.L / 2.0) * torch.sin(self.L * phi)

    def lattice_distance(self, phi: torch.Tensor) -> torch.Tensor:
        """Maslov 격자점 (pi/2 + m*pi)까지의 거리"""
        spacing = torch.pi / self.L  # cos^2의 영점 간격은 pi/L
        shifted = phi - torch.pi / (2.0 * self.L)
        remainder = torch.remainder(shifted, spacing)
        return torch.minimum(remainder, spacing - remainder)

    def __repr__(self):
        return f"<MaslovCorrectedPQO: L={self.L}, cos^2 variant>"


# =================================================================
# 직접 실행 시 무결성 테스트
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL] Phase Quantisation Operator Test ---")

    pqo = PhaseQuantisationOperator(L=2)
    print(f"연산자: {pqo}")

    # 테스트 1: 격자점에서 P=0 확인 (phi = k*pi for L=2)
    lattice_points = torch.tensor([0.0, torch.pi, 2*torch.pi, -torch.pi],
                                   dtype=PrecisionManager.REAL_DTYPE)
    p_at_lattice = pqo(lattice_points)
    print(f"\n[Test 1] 격자점 phi = k*pi 에서의 P 값:")
    print(f"  phi: {lattice_points.numpy()}")
    print(f"  P:   {p_at_lattice.numpy()}")
    lattice_ok = torch.allclose(p_at_lattice, torch.zeros_like(p_at_lattice), atol=1e-14)
    print(f"  P=0 확인: {'OK' if lattice_ok else 'FAIL'}")

    # 테스트 2: 격자 중점에서 P=1 확인 (phi = (k+0.5)*pi for L=2)
    midpoints = torch.tensor([torch.pi/2, 3*torch.pi/2, -torch.pi/2],
                              dtype=PrecisionManager.REAL_DTYPE)
    p_at_mid = pqo(midpoints)
    print(f"\n[Test 2] 격자 중점 phi = (k+0.5)*pi 에서의 P 값:")
    print(f"  phi: {midpoints.numpy()}")
    print(f"  P:   {p_at_mid.numpy()}")
    mid_ok = torch.allclose(p_at_mid, torch.ones_like(p_at_mid), atol=1e-14)
    print(f"  P=1 확인: {'OK' if mid_ok else 'FAIL'}")

    # 테스트 3: 기울기가 격자점에서 0인지 확인
    grad_at_lattice = pqo.gradient(lattice_points)
    print(f"\n[Test 3] 격자점에서 기울기:")
    print(f"  grad: {grad_at_lattice.numpy()}")
    grad_ok = torch.allclose(grad_at_lattice, torch.zeros_like(grad_at_lattice), atol=1e-14)
    print(f"  grad=0 확인: {'OK' if grad_ok else 'FAIL'}")

    # 테스트 4: autograd 호환성
    phi_auto = torch.tensor([0.3, 1.2, 2.7], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)
    p_auto = pqo(phi_auto)
    p_auto.sum().backward()
    print(f"\n[Test 4] Autograd 기울기:")
    print(f"  autograd: {phi_auto.grad.numpy()}")
    print(f"  명시적:   {pqo.gradient(phi_auto.detach()).numpy()}")
    autograd_ok = torch.allclose(phi_auto.grad, pqo.gradient(phi_auto.detach()), atol=1e-12)
    print(f"  일치 확인: {'OK' if autograd_ok else 'FAIL'}")

    # 테스트 5: QROP-Net 대응 (L=2048)
    pqo_qrop = PhaseQuantisationOperator(L=2048)
    phi_qrop = torch.tensor([0.0, 2*torch.pi/2048, torch.pi/2048],
                             dtype=PrecisionManager.REAL_DTYPE)
    p_qrop = pqo_qrop(phi_qrop)
    print(f"\n[Test 5] QROP-Net L=2048:")
    print(f"  phi: {phi_qrop.numpy()}")
    print(f"  P:   {p_qrop.numpy()}")
    print(f"  P(0)=0 확인: {'OK' if abs(p_qrop[0].item()) < 1e-14 else 'FAIL'}")
    print(f"  P(2pi/L)=0 확인: {'OK' if abs(p_qrop[1].item()) < 1e-14 else 'FAIL'}")
    print(f"  P(pi/L)=1 확인: {'OK' if abs(p_qrop[2].item() - 1.0) < 1e-10 else 'FAIL'}")

    if all([lattice_ok, mid_ok, grad_ok, autograd_ok]):
        print("\n[OK] 위상 양자화 연산자 전체 테스트 통과")
    else:
        print("\n[FAIL] 일부 테스트 실패")
