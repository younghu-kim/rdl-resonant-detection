"""
=============================================================================
[Project RDL] Core Math Primitives - Complex Tensor Safe Operations
=============================================================================
이 모듈은 복소 텐서의 위상(Phase)과 진폭(Amplitude)을 추출하고,
나눗셈 및 로그 연산 시 발생하는 영분모(Zero-Division) 및 특이점 붕괴를 수학적으로 방어합니다.
모든 연산은 PyTorch의 Autograd 그래프를 완벽하게 유지(Differentiable)하도록 설계되었습니다.
"""

import torch
from typing import Tuple

from gdl.rdl.constants import R_CONST

def get_phase_amp(z: torch.Tensor, eps: float = R_CONST.EPSILON) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    복소 텐서 Z에서 미분 연결(Autograd Graph)을 유지하며 위상(Φ)과 진폭(A)을 추출합니다.

    [수학적 방어 기제]
    일반적인 z.abs()는 z=0일 때 V자 형태로 꺾여 1차/2차 미분 시 NaN 붕괴를 일으킵니다.
    이를 막기 위해 A = sqrt(Re(Z)^2 + Im(Z)^2 + ε) 형태의 Epsilon-safe
    방식을 사용하여 무한 번 미분 가능(C^∞)한 상태를 보장합니다.

    Args:
        z (torch.Tensor): 복소 텐서 (Complex128 강제)
        eps (float): 미분 붕괴 방지용 초미세 상수

    Returns:
        Tuple[torch.Tensor, torch.Tensor]: 위상(Phase) 텐서, 진폭(Amplitude) 텐서 (Float64)
    """
    real = z.real
    imag = z.imag

    # 1. 진폭(A): 특이점(0)에서의 미분 발산을 막기 위한 Epsilon-safe Euclidean Norm
    amp = torch.sqrt(real**2 + imag**2 + eps)

    # 2. 위상(Φ): arctan2(Im, Re)를 사용하여 -π ~ π 값을 안전하게 반환
    # real과 imag가 완벽한 0일 때 발생할 수 있는 극단적 그래디언트 요동 및 NaN을 막기 위해 미세 마진 삽입
    mask = (real == 0.0) & (imag == 0.0)
    safe_real = torch.where(mask, real + eps, real)
    phase = torch.atan2(imag, safe_real)

    return phase, amp


def safe_complex_div(num: torch.Tensor, den: torch.Tensor, eps: float = R_CONST.EPSILON) -> torch.Tensor:
    """
    영분모 방지(Epsilon-safe) 복소수 나눗셈 연산을 수행합니다.

    [수학적 방어 기제]
    복소 나눗셈 N / D 에서 단순히 분모에 eps를 더하면 (D + eps), 복소수의 '위상(Phase) 각도'가 틀어집니다.
    공명 수론에서 위상의 왜곡은 '정렬 실패'를 의미하므로 치명적입니다.
    따라서 N / D = (N * D^*) / |D|^2 로 변환한 뒤, 실수 스칼라인 |D|^2 에만 eps를 더하여
    위상은 100% 보존하고 진폭만 붕괴하지 않도록 방어합니다.

    Args:
        num (torch.Tensor): 분자 복소 텐서
        den (torch.Tensor): 분모 복소 텐서
        eps (float): 영분모 방지용 Epsilon

    Returns:
        torch.Tensor: 위상 왜곡 없이 안전하게 계산된 나눗셈 결과 복소 텐서
    """
    # 1. |D|^2 + ε (안전한 분모 실수 스칼라)
    den_mag_sq = den.real**2 + den.imag**2
    safe_den = den_mag_sq + eps

    # 2. N * D^* (분자에 분모의 켤레 복소수 곱셈)
    # torch.conj()는 내부적으로 미분 시 켤레 미분을 올바르게 추적함
    num_div = num * torch.conj(den)

    # 3. 최종 복소수 나눗셈 (safe_den은 실수이므로 Broadcasting 나눗셈 적용)
    return num_div / safe_den


def safe_complex_log(z: torch.Tensor, eps: float = R_CONST.EPSILON) -> torch.Tensor:
    """
    (Phase 3의 방향성 로그 미분 D_k log Z 계산 시 특이점 붕괴를 방지하기 위해 선제적으로 추가)
    복소 텐서의 자연로그 log(Z)를 안전하게 계산합니다.

    log(Z) = log(|Z|) + i * arg(Z)

    Args:
        z (torch.Tensor): 입력 복소 텐서
        eps (float): log(0) 폭주 (-inf) 방지용 마진

    Returns:
        torch.Tensor: 복소 로그 텐서
    """
    phase, amp = get_phase_amp(z, eps=eps)

    # amp는 get_phase_amp 내부에서 이미 eps가 더해져 있으므로 log(0)이 발생하지 않음
    log_amp = torch.log(amp)

    return torch.complex(log_amp, phase)


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    from gdl.rdl.constants import PrecisionManager
    PrecisionManager.setup_precision()

    print("\n--- [RDL Core Math] Complex Ops Test ---")

    # 1. 0을 포함한 복소 텐서 생성
    z = torch.tensor([3.0 + 4.0j, 0.0 + 0.0j], dtype=PrecisionManager.COMPLEX_DTYPE, requires_grad=True)
    print(f"입력 Z:\n{z.detach().cpu().numpy()}")

    # 2. Phase & Amp 추출 및 미분 안정성 테스트
    phase, amp = get_phase_amp(z)

    # 2차 미분(Hessian) 연산 시뮬레이션
    amp_sum = amp.sum()

    # 1차 미분
    grad_z = torch.autograd.grad(amp_sum, z, create_graph=True)[0]

    # 2차 미분 (Hessian Trace 근사)
    grad_z_sum = grad_z.real.sum() + grad_z.imag.sum()
    grad2_z = torch.autograd.grad(grad_z_sum, z)[0]

    if torch.isnan(grad2_z).any() or torch.isinf(grad2_z).any():
        print("❌ 실패: 진폭 2차 미분(Hessian) 그래프에 NaN/Inf가 발생했습니다.")
    else:
        print("✅ 성공: 0(Singularity)에서도 2차 미분이 폭주하지 않고 매끄러운(C-infinity) 상태를 유지합니다.")
        print(f"   └─ 1차 미분값: {grad_z.detach().cpu().numpy()}")
        print(f"   └─ 2차 미분값: {grad2_z.detach().cpu().numpy()}")

    # 3. 안전한 나눗셈 & 로그 테스트
    num = torch.tensor([1.0 + 1.0j, 1.0 + 1.0j], dtype=PrecisionManager.COMPLEX_DTYPE)
    den = torch.tensor([1.0 + 0.0j, 0.0 + 0.0j], dtype=PrecisionManager.COMPLEX_DTYPE) # 두 번째 요소가 0

    div_res = safe_complex_div(num, den)
    log_res = safe_complex_log(den)

    if torch.isnan(div_res).any() or torch.isnan(log_res).any():
        print("❌ 실패: 나눗셈 또는 로그 연산에서 NaN이 발생했습니다.")
    else:
        print("✅ 성공: 영분모(Zero-Denominator) 및 Log(0)가 위상 왜곡 없이 완벽하게 방어되었습니다.")
        print(f"   └─ 안전 나눗셈 결과: {div_res.detach().cpu().numpy()}")
        print(f"   └─ 안전 로그 결과:   {log_res.detach().cpu().numpy()}\n")
