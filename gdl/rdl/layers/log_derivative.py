"""
=============================================================================
[Project RDL] Multigroup Resonant Forward - Directional Log-Derivative Extractor
=============================================================================
이 모듈은 기존 딥러닝의 CNN 필터나 Attention 메커니즘을 완전히 대체하는
'방향성 로그 미분(Directional Log-Derivative)' 추출기입니다.

입력된 복소 텐서 Z를 위상 공간의 좌표(s)로 취급하고, 내부의 학습 가능한 복소 함수 f(Z)에 대해
다원군 채널(K개) 방향으로 미소 증분(δ_k)을 더하고 뺀 뒤,
로그 도메인 중앙차분을 수행하여 폭주가 억제된 초정밀 미분 텐서 L_k를 산출합니다.

[구현 수학 식 매핑]
- f(Z) = Z * W^T + b (복소 선형 공간 변환, 신경망의 가중치 역할)
- E20: D_k log_xi(s) ≈ [log_xi(s + δ_k) - log_xi(s - δ_k)] / (2h)
- E22: L_k(Z) = D_k log f(Z)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.core.complex_ops import safe_complex_log
from gdl.rdl.core.multigroup import MultigroupFrame

class DirectionalLogDerivative(nn.Module):
    """
    입력 복소 텐서 Z의 K개 방향(δ_k)에 대한 로그 미분 텐서 L_k를 병렬 추출하는 레이어.
    """
    def __init__(self, in_features: int, out_features: int, frame: MultigroupFrame):
        """
        Args:
            in_features (int): 입력 복소 텐서의 Feature 차원 크기
            out_features (int): 출력 복소 텐서의 Feature 차원 크기 (다음 레이어로 전달될 크기)
            frame (MultigroupFrame): 사전 캐싱된 다원군 프레임 (δ_k 텐서 보유)
        """
        super(DirectionalLogDerivative, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.frame = frame
        self.h = frame.h

        # [수학적 방어 기제 1] 복소 가중치 초기화 분리
        # PyTorch의 Complex Parameter 최적화를 Adam 등 옵티마이저에서 안전하게 수행하기 위해,
        # 실수부와 허수부를 독립적인 Float64 실수 파라미터로 선언하고 연산 시 결합합니다.
        # Xavier(Glorot) 초기화 기법 적용
        std = (1.0 / in_features) ** 0.5
        factory_kwargs = {'dtype': PrecisionManager.REAL_DTYPE}

        self.weight_re = nn.Parameter(torch.randn(out_features, in_features, **factory_kwargs) * std)
        self.weight_im = nn.Parameter(torch.randn(out_features, in_features, **factory_kwargs) * std)

        self.bias_re = nn.Parameter(torch.zeros(out_features, **factory_kwargs))
        self.bias_im = nn.Parameter(torch.zeros(out_features, **factory_kwargs))

    def _complex_transform(self, z: torch.Tensor) -> torch.Tensor:
        """
        내부 복소 선형 변환 함수 f(Z) = Z * W^T + b
        복소수 행렬곱: (Re(Z)Re(W) - Im(Z)Im(W)) + i(Re(Z)Im(W) + Im(Z)Re(W))
        """
        real_out = F.linear(z.real, self.weight_re, self.bias_re) - \
                   F.linear(z.imag, self.weight_im, None)
        imag_out = F.linear(z.real, self.weight_im, self.bias_im) + \
                   F.linear(z.imag, self.weight_re, None)

        return torch.complex(real_out, imag_out)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        """
        Args:
            z (torch.Tensor): 승격된 입력 복소 텐서 (Shape: [Batch, ..., in_features])

        Returns:
            torch.Tensor: 방향성 로그 미분 텐서 L_k
                          (Shape: [Batch, ..., K, out_features])
        """
        # [방어 기제 2] Float64 / Complex128 강제 점검
        if z.dtype != PrecisionManager.COMPLEX_DTYPE:
            z = z.to(dtype=PrecisionManager.COMPLEX_DTYPE)

        # 프레임의 미소 증분 텐서 로드
        self.frame.to(z.device)
        delta_k = self.frame.delta_k  # Shape: [K]
        K = self.frame.K

        # [수학적 방어 기제 3] 브로드캐스팅(Broadcasting) 차원 정렬
        # for-loop를 사용하면 GPU 메모리 병목이 발생하므로,
        # Z를 [Batch, ..., 1, Features]로, δ_k를 [1, ..., K, 1] 로 맞추어
        # O(1) 배치 연산으로 K개 방향의 평행우주를 동시에 스캔합니다.
        z_expanded = z.unsqueeze(-2)

        shape_diff = [1] * z.dim()
        shape_diff.insert(-1, K)
        delta_k_expanded = delta_k.view(*shape_diff)

        # 1. Z 주변의 K개 방향 미소 증분 텐서 생성
        # Shape: [Batch, ..., K, in_features]
        z_plus = z_expanded + delta_k_expanded
        z_minus = z_expanded - delta_k_expanded

        # 2. 신경망 변환 f(Z_plus), f(Z_minus) 통과
        # F.linear는 마지막 차원(in_features -> out_features)만 변환하므로 K차원은 완벽히 보존됩니다.
        f_plus = self._complex_transform(z_plus)
        f_minus = self._complex_transform(z_minus)

        # 3. [방어 기제 4] 로그 도메인 안전 변환 (Singularity Evasion)
        # Phase 1에서 만든 Epsilon-safe 로그 함수를 통해 0으로 수렴하는 특이점 구간에서의 NaN 폭주를 원천 차단합니다.
        log_f_plus = safe_complex_log(f_plus)
        log_f_minus = safe_complex_log(f_minus)

        # 4. 중앙차분(Central Difference)을 통한 방향성 로그 미분 추출 (수식 E20)
        # L_k = (log f(Z+δ_k) - log f(Z-δ_k)) / 2h
        L_k = (log_f_plus - log_f_minus) / (2.0 * self.h)

        # Output Shape: [Batch, ..., K, out_features]
        return L_k

    def extra_repr(self) -> str:
        return f'in_features={self.in_features}, out_features={self.out_features}, K={self.frame.K}, h={self.h:.2e}'


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Layers] Directional Log-Derivative Extractor Test ---")

    # 1. 5채널(5ch) 다원군 프레임 생성
    frame_5ch = MultigroupFrame(channel_type='5ch')
    print(f"로드된 프레임: {frame_5ch}")

    # 2. 로그 미분 추출기 인스턴스화
    in_dim = 4
    out_dim = 8
    extractor = DirectionalLogDerivative(in_features=in_dim, out_features=out_dim, frame=frame_5ch)

    # 3. 가상의 복소 데이터 생성 (Batch=2, Seq=3, Features=4)
    # 시계열/그리드 데이터의 복합 구조를 모사 (0벡터에 가까운 노이즈 포함)
    z_real = torch.randn((2, 3, in_dim), dtype=PrecisionManager.REAL_DTYPE) * 0.01
    z_imag = torch.randn((2, 3, in_dim), dtype=PrecisionManager.REAL_DTYPE) * 0.01
    z_input = torch.complex(z_real, z_imag).requires_grad_(True)

    print(f"입력 복소 텐서 Z Shape: {z_input.shape} | Dtype: {z_input.dtype}")

    # 4. Forward Pass: K개 방향의 로그 미분 추출
    L_k_out = extractor(z_input)

    print(f"\n추출된 방향성 로그 미분 텐서 L_k Shape: {L_k_out.shape} (기대값: [2, 3, {frame_5ch.K}, {out_dim}])")

    # 5. 차원 및 NaN 무결성 검증
    expected_shape = (2, 3, frame_5ch.K, out_dim)
    if L_k_out.shape == expected_shape:
        print(f"✅ 성공: 브로드캐스팅(Broadcasting)을 통한 다원군 차원(K={frame_5ch.K}) 확장이 완벽합니다.")
    else:
        print(f"❌ 실패: L_k의 차원에 오류가 있습니다. 실제값: {L_k_out.shape}")

    if not torch.isnan(L_k_out).any() and not torch.isinf(L_k_out).any():
        print(f"✅ 성공: 극한의 차분(h={frame_5ch.h})과 로그 변환 중 특이점 붕괴(NaN/Inf)가 100% 방어되었습니다.")
    else:
        print("❌ 실패: L_k 텐서 내부에 NaN/Inf가 발생했습니다.")

    # 6. 미분 역전파(Autograd) 보존 검증
    loss = L_k_out.abs().sum()
    loss.backward()

    if z_input.grad is not None and not torch.isnan(z_input.grad).any():
        print("✅ 성공: 복소 공간 미소 이동 및 로그 차분에 대한 역전파 그래프가 100% 보존되었습니다!\n")
    else:
        print("❌ 실패: 미분 그래프가 단절되었거나 역전파 중 NaN이 발생했습니다.\n")
