"""
=============================================================================
[Project RDL] Core Math Primitives - Multigroup Frame Generator
=============================================================================
다원군(Multigroup) 구조를 위한 방향 벡터와 미소 복소 증분 텐서를 생성/캐싱한다.

[구현 수식 매핑 - unified_en.tex 참조]
- eq:channel_vector: e_k = (cos(theta_k), lambda_k * sin(theta_k))
- eq:channel_set: (-30, 1.10, 0.40), (15, 0.85, 0.35), (60, 1.00, 0.25)
- eq:aniso_scale: M(t) = diag(1/mu(t), 1/nu(t))
- eq:displacement: delta_k(t) = h * (v_{k,x} + i * v_{k,y})

레거시 등방 채널:
- iso3: {-60, 0, +60} (대칭, 논문에 없음)
- 5ch: {-72, -36, 0, +36, +72} (대칭, 논문에 없음)
"""

import torch
import math
from typing import Tuple

from gdl.rdl.constants import R_CONST, PrecisionManager

class MultigroupFrame:
    """
    다원군 프레임을 생성하고 텐서 연산을 위한 상수 벡터들을 유지하는 캐시 클래스.
    이 클래스의 모든 텐서는 requires_grad=False로 설정되어,
    오류 역전파 시 연산 그래프에 불필요한 노드가 생기는 것을 방지합니다.
    """

    # 지원하는 채널별 각도 배열 (단위: Degree)
    # 0도를 기준으로 완벽한 좌우 대칭을 이루어 미러(Mirror) 불변성을 보장함.
    SUPPORTED_CHANNELS = {
        'iso3': [-60.0, 0.0, 60.0],
        '5ch': [-72.0, -36.0, 0.0, 36.0, 72.0]
    }

    def __init__(self, channel_type: str = 'iso3', h: float = None, device: torch.device = None):
        """
        Args:
            channel_type (str): 'iso3' (3채널 기본) 또는 '5ch' (5채널 확장)
            h (float): 중앙차분 미소 스텝. None일 경우 전역 상수 R_CONST.H 사용
            device (torch.device): 텐서가 상주할 디바이스 (CPU/GPU)
        """
        if channel_type not in self.SUPPORTED_CHANNELS:
            raise ValueError(f"❌ [Multigroup Error] 지원하지 않는 채널 타입입니다: {channel_type}. 'iso3' 또는 '5ch'를 사용하세요.")

        self.channel_type = channel_type
        self.h = h if h is not None else R_CONST.H
        self.device = device if device is not None else torch.device('cpu')

        self._build_frame()

    def _build_frame(self):
        """내부적으로 각도와 증분 텐서를 생성하여 멤버 변수로 고정(Cache)합니다."""

        # 1. 각도(Degree) -> 라디안(Radian) 변환 텐서 생성 (Float64 강제)
        degrees = self.SUPPORTED_CHANNELS[self.channel_type]
        self.K = len(degrees)

        radians = [math.radians(d) for d in degrees]
        self.theta = torch.tensor(radians, dtype=PrecisionManager.REAL_DTYPE, device=self.device)

        # [방어 기제] 프레임 텐서는 방향을 찌를 '기준 좌표계'이므로 미분 그래프에서 제외
        self.theta.requires_grad_(False)

        # 2. 방향 벡터 성분 추출 v_k = (cos(Θ_k), sin(Θ_k)) (E14, E15)
        self._v_x = torch.cos(self.theta)
        self._v_y = torch.sin(self.theta)

        # [K, 2] 형태의 실수 행렬로 조립 (Block 1.4 Ridge Solver의 M 행렬 생성 시 사용)
        self._v_k = torch.stack([self._v_x, self._v_y], dim=1)

        # 3. 복소 미소 증분 텐서 δ_k = h * (v_x + i * v_y) 조립 (E16)
        # 딥러닝 Forward 패스에서 Z(x + δ_k) 로 데이터를 입체적으로 스캔할 때 쓰입니다.
        real_part = self.h * self._v_x
        imag_part = self.h * self._v_y
        self._delta_k = torch.complex(real_part, imag_part)
        self._delta_k.requires_grad_(False)

    @property
    def delta_k(self) -> torch.Tensor:
        """사전 캐싱된 미소 증분 복소 텐서 δ_k 반환 (Shape: [K])"""
        return self._delta_k

    @property
    def v_k_components(self) -> Tuple[torch.Tensor, torch.Tensor]:
        """사전 캐싱된 방향 벡터의 (x성분, y성분) 반환 (Block 1.4 릿지 솔버 조립용)"""
        return self._v_x, self._v_y

    @property
    def v_k(self) -> torch.Tensor:
        """방향 벡터 행렬을 반환합니다. (Shape: [K, 2])"""
        return self._v_k

    def to(self, device: torch.device):
        """
        프레임 내부의 모든 텐서를 지정된 디바이스(예: 'cuda')로 안전하게 이동시킵니다.
        (PyTorch의 nn.Module.to() 와 동일한 디자인 패턴)
        """
        self.device = device
        self.theta = self.theta.to(device)
        self._v_x = self._v_x.to(device)
        self._v_y = self._v_y.to(device)
        self._v_k = self._v_k.to(device)
        self._delta_k = self._delta_k.to(device)
        return self

    def __repr__(self):
        return f"<MultigroupFrame: {self.channel_type} (K={self.K}), h={self.h:.2e}, device={self.device}>"


class AnisotropicMultigroupFrame:
    """
    논문 정합 비등방 다원군 프레임 (unified_en.tex).

    등방 프레임(MultigroupFrame)과 달리:
    1. 채널별 비등방 축척 lambda_k 지원 (eq:channel_vector)
    2. 논문 권장 가중치 w_k 사전 지정 (eq:channel_set)
    3. 비등방 축척 행렬 M(t) 지원 (eq:aniso_scale)

    수식:
        e_k = (cos(theta_k), lambda_k * sin(theta_k))
        delta_k = h * (e_{k,x} + i * e_{k,y})
    """

    def __init__(self, h: float = None, device: torch.device = None):
        """
        Args:
            h: 중심차분 미소단계. None이면 R_CONST.H_DIFF 사용.
            device: 텐서 디바이스
        """
        self.h = h if h is not None else R_CONST.H_DIFF
        self.device = device if device is not None else torch.device('cpu')

        self.channel_type = 'paper3ch'
        self.K = len(R_CONST.PAPER_CHANNEL_ANGLES)

        # 논문 eq:channel_set에서 직접 가져온 구성
        self._angles = R_CONST.PAPER_CHANNEL_ANGLES       # (-30, 15, 60)
        self._lambdas = R_CONST.PAPER_CHANNEL_LAMBDAS     # (1.10, 0.85, 1.00)
        self._weights = R_CONST.PAPER_CHANNEL_WEIGHTS      # (0.40, 0.35, 0.25)

        self._build_frame()

    def _build_frame(self):
        """비등방 방향 벡터와 복소 증분을 생성한다."""
        dtype = PrecisionManager.REAL_DTYPE
        dev = self.device

        radians = torch.tensor(
            [math.radians(d) for d in self._angles], dtype=dtype, device=dev
        )
        lambdas = torch.tensor(self._lambdas, dtype=dtype, device=dev)

        # eq:channel_vector: e_k = (cos(theta_k), lambda_k * sin(theta_k))
        self._v_x = torch.cos(radians)
        self._v_y = lambdas * torch.sin(radians)

        self._v_k = torch.stack([self._v_x, self._v_y], dim=1)  # [K, 2]

        # eq:displacement: delta_k = h * (v_x + i * v_y)
        self._delta_k = torch.complex(self.h * self._v_x, self.h * self._v_y)
        self._delta_k.requires_grad_(False)

        # 논문 가중치 (eq:channel_set)
        self._w_k = torch.tensor(self._weights, dtype=dtype, device=dev)
        self._w_k.requires_grad_(False)

    @property
    def delta_k(self) -> torch.Tensor:
        """복소 미소 증분 delta_k (Shape: [K])"""
        return self._delta_k

    @property
    def v_k_components(self) -> Tuple[torch.Tensor, torch.Tensor]:
        """방향 벡터 (v_x, v_y) 성분"""
        return self._v_x, self._v_y

    @property
    def v_k(self) -> torch.Tensor:
        """방향 벡터 행렬 [K, 2]"""
        return self._v_k

    @property
    def weights(self) -> torch.Tensor:
        """논문 사전 지정 가중치 w_k (Shape: [K]), sum=1"""
        return self._w_k

    def to(self, device: torch.device):
        """모든 텐서를 지정 디바이스로 이동"""
        self.device = device
        self._v_x = self._v_x.to(device)
        self._v_y = self._v_y.to(device)
        self._v_k = self._v_k.to(device)
        self._delta_k = self._delta_k.to(device)
        self._w_k = self._w_k.to(device)
        return self

    def __repr__(self):
        return (f"<AnisotropicMultigroupFrame: paper3ch (K={self.K}), "
                f"h={self.h:.2e}, lambdas={self._lambdas}, device={self.device}>")


# =================================================================
# 직접 실행 시 모듈 무결성 테스트
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Core Math] Multigroup Frame Generator Test ---")

    # 1. iso3 (3채널) 생성 및 테스트
    frame_iso3 = MultigroupFrame(channel_type='iso3')
    delta_iso3 = frame_iso3.delta_k
    vk_x3, vk_y3 = frame_iso3.v_k_components

    print("\n[iso3 Channel (-60°, 0°, 60°)]")
    print(f"  └─ v_k.x: {vk_x3.numpy()}")
    print(f"  └─ v_k.y: {vk_y3.numpy()}")
    print(f"  └─ delta_k: \n{delta_iso3.numpy()}")

    # 2. 5ch (5채널) 생성 및 테스트
    frame_5ch = MultigroupFrame(channel_type='5ch')
    delta_5ch = frame_5ch.delta_k
    vk_x5, vk_y5 = frame_5ch.v_k_components

    print("\n[5ch Channel (-72°, -36°, 0°, 36°, 72°)]")
    print(f"  └─ v_k.x: {vk_x5.numpy()}")
    print(f"  └─ v_k.y: {vk_y5.numpy()}")
    print(f"  └─ delta_k: \n{delta_5ch.numpy()}")

    # 3. 무결성 및 대칭성 검증
    # 좌우 대칭 합이 0인지 확인 (예: sin(-60) + sin(60) = 0)
    sym_check_iso3 = torch.abs(vk_y3.sum()) < 1e-15
    sym_check_5ch = torch.abs(vk_y5.sum()) < 1e-15

    # 증분의 절대값(Norm)이 정확히 상수 h와 일치하는지 확인
    amp_iso3 = torch.abs(delta_iso3)
    amp_check = torch.allclose(amp_iso3, torch.tensor(R_CONST.H, dtype=PrecisionManager.REAL_DTYPE))

    if sym_check_iso3 and sym_check_5ch and amp_check:
        print("\n✅ 성공: 다원군 프레임이 정밀한 미러 대칭성(Mirror Symmetry)을 유지하며 생성되었습니다.")
        print(f"   (미소 증분의 크기는 |δ_k| = {R_CONST.H} 로 완벽히 캐싱되어 병목 없이 호출 가능합니다.)\n")
    else:
        print("\n❌ 실패: 채널의 좌우 대칭성 또는 진폭이 깨졌습니다.")
