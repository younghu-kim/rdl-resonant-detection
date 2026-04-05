"""
=============================================================================
[Project RDL] Phase Space Embedding - Complex Linear Mapping Layer
=============================================================================
이 모듈은 현실의 실수(Real) 형태 데이터를 복소 위상 흐름의 세계(Z^(0))로
강제 승격(Ascension)시키는 게이트웨이(Gateway) 신경망 레이어입니다.

일반적인 딥러닝 임베딩과 달리, 입력 데이터 X에 대해 실수부(Real)와
허수부(Imaginary)를 담당하는 두 개의 독립적인 가중치 네트워크를 통과시켜,
초기 진폭(Amplitude)과 위상(Phase)의 기하학적 자유도를 부여합니다.

[구현 수학 식]
Z^{(0)} = \sigma(X W_{Re} + b_{Re}) + i \cdot \sigma(X W_{Im} + b_{Im})
"""

import torch
import torch.nn as nn

from gdl.rdl.constants import PrecisionManager

class ComplexLinearEmbedding(nn.Module):
    """
    실수 텐서 X를 복소 위상 공간 Z^(0)로 사상(Mapping)하는 임베딩 레이어.
    공명 딥러닝 시스템의 가장 첫 번째 관문입니다.
    """
    def __init__(self, in_features: int, out_features: int, activation: nn.Module = None):
        """
        Args:
            in_features (int): 입력 실수 데이터의 차원 크기
            out_features (int): 승격될 복소 위상 공간의 차원 크기
            activation (nn.Module, optional): 비선형 활성화 함수 (기본값: SiLU).
                                              2차 미분(곡률 A'')을 계산해야 하므로
                                              반드시 C^∞(무한 번 미분 가능) 함수를 권장합니다.
        """
        super(ComplexLinearEmbedding, self).__init__()

        self.in_features = in_features
        self.out_features = out_features

        # [수학적 방어 기제 1] 가중치를 Float64 초고정밀도로 강제 생성
        # nn.Linear는 기본적으로 float32로 가중치를 만들기 때문에 이를 통제합니다.
        factory_kwargs = {'dtype': PrecisionManager.REAL_DTYPE}

        # 1. 실수부(Real) 투영을 위한 독립적 가중치 (에너지/진폭 기저)
        self.fc_re = nn.Linear(in_features, out_features, **factory_kwargs)

        # 2. 허수부(Imaginary) 투영을 위한 독립적 가중치 (흐름/위상 기저)
        self.fc_im = nn.Linear(in_features, out_features, **factory_kwargs)

        # 3. 비선형 활성화 함수 (주어지지 않으면 SiLU 사용)
        # SiLU(x * sigmoid(x))나 Tanh는 0 근방에서 매끄럽게 미분 가능하여 곡률 붕괴를 막습니다.
        self.activation = activation if activation is not None else nn.SiLU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        실수 텐서 X를 복소 텐서 Z로 사상합니다.

        Args:
            x (torch.Tensor): 입력 실수 텐서 (Shape: [Batch, ..., in_features])

        Returns:
            torch.Tensor: 승격된 복소 텐서 Z (Shape: [Batch, ..., out_features], dtype: Complex128)
        """
        # [수학적 방어 기제 2] 데이터 로더에서 들어온 텐서가 Float64가 아닐 경우 강제 캐스팅
        # 외부 데이터의 정밀도 손실을 원천 차단합니다.
        if x.dtype != PrecisionManager.REAL_DTYPE:
            x = x.to(dtype=PrecisionManager.REAL_DTYPE)

        # 독립적인 공간으로의 투영 및 비선형 활성화
        real_part = self.activation(self.fc_re(x))
        imag_part = self.activation(self.fc_im(x))

        # 실수부와 허수부를 결합하여 초고정밀 복소 텐서 Z 조립
        z_complex = torch.complex(real_part, imag_part)

        return z_complex

    def extra_repr(self) -> str:
        return f'in_features={self.in_features}, out_features={self.out_features}, activation={self.activation}'


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Layers] Complex Linear Embedding Test ---")

    # 1. 가상의 입력 데이터 생성 (Batch=4, Features=10)
    batch_size = 4
    in_dim = 10
    out_dim = 16

    # 의도적으로 float32 타입의 데이터를 생성하여 방어 기제가 작동하는지 확인
    x_real = torch.randn((batch_size, in_dim), dtype=torch.float32, requires_grad=True)
    print(f"입력 실수 데이터 X (Shape: {x_real.shape}) | Dtype: {x_real.dtype}")

    # 2. 임베딩 레이어 인스턴스화
    embedder = ComplexLinearEmbedding(in_features=in_dim, out_features=out_dim)

    # 3. Forward Pass: 복소 위상 공간으로 승격
    z_out = embedder(x_real)

    print(f"승격된 복소 데이터 Z (Shape: {z_out.shape}) | Dtype: {z_out.dtype}")

    # 4. 검증: Z의 데이터 타입 확인
    if z_out.dtype == PrecisionManager.COMPLEX_DTYPE:
        print("✅ 성공: 실수 데이터가 128-bit 초고정밀도 복소 위상 공간으로 완벽히 승격되었습니다.")
    else:
        print("❌ 실패: 차원 또는 데이터 타입이 올바르지 않습니다.")

    # 5. 미분 그래프(Autograd) 연결 및 2차 미분 가능성 확인
    amp_sq = z_out.real**2 + z_out.imag**2
    loss = amp_sq.sum()

    # 1차 미분
    grad1 = torch.autograd.grad(loss, x_real, create_graph=True)[0]

    # 2차 미분
    grad2 = torch.autograd.grad(grad1.sum(), x_real)[0]

    if grad2 is not None and not torch.isnan(grad2).any():
        print("✅ 성공: C^∞ 활성화 함수(SiLU)를 통해 2차 미분(곡률 A'') 그래프가 단절이나 NaN 없이 보존됩니다.\n")
    else:
        print("❌ 실패: 2차 미분 중 단절이 발생했습니다.\n")
