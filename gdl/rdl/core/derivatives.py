"""
=============================================================================
[Project RDL] Core Math Primitives - Numerical Derivatives & Curvature Extractor
=============================================================================
이 모듈은 진폭(A) 또는 임의의 텐서에 대해 1차 도함수(Gradient, A')와
2차 도함수(Curvature/Hessian Trace, A'')를 추출합니다.

PyTorch의 autograd 연산 그래프를 파괴하지 않고(create_graph=True)
고차 미분을 수행함으로써, 추출된 곡률 자체가 다시 손실 함수(Loss)에 포함되어
최종 가중치 역전파(Backpropagation)로 이어질 수 있도록 설계되었습니다.
"""

import torch
from typing import Tuple

from gdl.rdl.constants import PrecisionManager

def get_first_derivative(y: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
    """
    텐서 y의 x에 대한 1차 도함수(Gradient, y')를 추출합니다.
    공명 조건 중 '진폭 극점(A' = 0)'을 확인하고 손실(Loss)로 강제하기 위해 사용됩니다.

    Args:
        y (torch.Tensor): 목적 함수 텐서 (예: 진폭 A)
        x (torch.Tensor): 입력 변수 텐서 (예: 입력 데이터 또는 위상 흐름 t)

    Returns:
        torch.Tensor: x와 동일한 차원을 갖는 1차 도함수 텐서
    """
    if not x.requires_grad:
        raise ValueError("❌ [Curvature Error] 입력 텐서 x는 반드시 requires_grad=True 상태여야 합니다.")

    # y가 다차원 텐서(Batch)일 경우를 대비해 ones_like 텐서 생성 (Jacobian-Vector Product의 역할)
    # 이는 y.sum()을 한 뒤 미분하는 것과 수학적으로 동일하며, 각 요소의 독립적인 미분값을 안전하게 추출합니다.
    grad_outputs = torch.ones_like(y, requires_grad=False)

    # [수학적 방어 기제: create_graph=True]
    # 이 옵션이 없으면 1차 미분을 구하는 즉시 파이토치가 메모리를 비워버립니다.
    # 우리는 A'를 구한 뒤 A''까지 구하고, 이 값들을 Loss에 넣어 또 미분해야 하므로
    # 연산 그래프를 영구적으로 보존(retain_graph)하고 확장(create_graph)해야 합니다.
    first_deriv = torch.autograd.grad(
        outputs=y,
        inputs=x,
        grad_outputs=grad_outputs,
        create_graph=True,
        retain_graph=True,
        only_inputs=True,
        allow_unused=True
    )[0]

    if first_deriv is None:
        # x와 y 사이에 미분적 연결 고리가 없는 경우 안전하게 0 텐서 반환
        return torch.zeros_like(x)

    return first_deriv


def get_resonance_derivatives(y: torch.Tensor, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    텐서 y의 x에 대한 1차 도함수(A')와 2차 도함수(Curvature, A'')를 동시에 추출합니다.
    요소별(Element-wise) 매핑에 대한 Hessian의 대각합(Trace / Laplacian)과 동치입니다.

    공명 조건 중 '곡률 안정 및 에너지 응축(A'' < 0)'을 강제하기 위해 사용됩니다.

    Args:
        y (torch.Tensor): 목적 함수 텐서 (예: 진폭 A)
        x (torch.Tensor): 입력 변수 텐서

    Returns:
        Tuple[torch.Tensor, torch.Tensor]: (1차 도함수 텐서, 2차 도함수 텐서)
    """
    # 1단계: 1차 미분(Gradient) 추출
    first_deriv = get_first_derivative(y, x)

    # 2단계: 1차 미분에 대해 다시 x로 미분하여 2차 미분(Curvature) 추출
    # (1차 도함수를 목적 함수로 삼아 한 번 더 미분)
    second_deriv = get_first_derivative(first_deriv, x)

    return first_deriv, second_deriv


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Core Math] Derivatives & Curvature Test ---")

    # 1. 2차 함수 생성: y = -x^2 + 4x
    # (이 함수는 x=2에서 극대값 y=4를 가지며, 곡률은 항상 -2로 오목(안정)한 함수입니다.)
    # 공명 수론에서 에너지가 응축되는 진폭 A(t)의 이상적인 모습을 모사합니다.
    x = torch.tensor([1.0, 2.0, 3.0], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)
    y = -(x**2) + 4*x

    print(f"입력 x (t):       {x.detach().cpu().numpy()}")
    print(f"진폭 y (A):       {y.detach().cpu().numpy()}")

    # 2. 도함수 및 곡률 추출
    dy_dx, d2y_dx2 = get_resonance_derivatives(y, x)

    print(f"1차 미분 (A'):    {dy_dx.detach().cpu().numpy()}")
    print(f"2차 곡률 (A''):   {d2y_dx2.detach().cpu().numpy()}\n")

    # 3. 공명 조건 (Gate B) 판별
    # x=2.0 위치(인덱스 1)에서 A'=0, A''<0 인지 확인
    is_critical_point = torch.abs(dy_dx[1]) < 1e-10
    is_stable_curvature = d2y_dx2[1] < 0.0

    if is_critical_point and is_stable_curvature:
        print("✅ 성공: 극점(A'=0)과 곡률 안정성(A''<0)이 성공적으로 포착되었습니다. (Gate B 통과)")
    else:
        print("❌ 실패: 곡률 추출에 문제가 있습니다.")

    # 4. 3차 미분(Loss 역전파) 가능성 테스트
    try:
        # 가상의 손실(Loss) 함수: A'를 0으로 만들고 A''가 양수(불안정)일 때 페널티 부여
        dummy_loss = (dy_dx**2).sum() + (torch.relu(d2y_dx2)).sum()
        dummy_loss.backward()
        print("✅ 성공: Loss에 대한 역전파(Backprop)가 완료되었습니다. 3차 미분 그래프가 완벽히 보존되었습니다!\n")
    except Exception as e:
        print(f"❌ 실패: 미분 그래프가 끊어졌습니다. 에러: {e}\n")
