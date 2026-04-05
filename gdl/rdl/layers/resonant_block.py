"""
=============================================================================
[Project RDL] Multigroup Resonant Forward - Integrated Resonant Block
=============================================================================
다원군 프레임의 K개 채널에서 방향성 로그 미분(L_k)을 추출하고,
가중합하여 유효 로그기울기(L_G = Σ w_k · L_k)를 출력하는 핵심 은닉층.

[구현 수식 매핑 - unified_en.tex 참조]
- eq:central_diff: D_k log f ≈ [f(z+δ_k) - f(z-δ_k)] / (2δ_k)
- eq:lhat: L̂ = Σ w_k · D_k log f
- eq:channel_set: 논문 3채널 (-30°, 15°, 60°) + 가중치 (0.40, 0.35, 0.25)
- eq:channel_vector: e_k = (cos θ_k, λ_k sin θ_k) 비등방 방향 벡터

채널 모드:
- 'paper3ch': 논문 정합 비등방 3채널 (AnisotropicMultigroupFrame + FixedWeightSolver)
- 'iso3': 레거시 등방 3채널 (MultigroupFrame + RidgeWeightSolver)
- '5ch': 레거시 등방 5채널 (MultigroupFrame + RidgeWeightSolver)
"""

import torch
import torch.nn as nn
from typing import Optional

from gdl.rdl.constants import R_CONST, PrecisionManager
from gdl.rdl.core.multigroup import MultigroupFrame, AnisotropicMultigroupFrame
from gdl.rdl.core.ridge_solver import RidgeWeightSolver, FixedWeightSolver
from gdl.rdl.layers.log_derivative import DirectionalLogDerivative
from gdl.rdl.layers.assembler import EffectiveLogDerivativeAssembler

class ResonantBlock(nn.Module):
    """
    다원군 공명 전파 은닉층. L̂ = Σ w_k · D_k log f (eq:lhat)

    채널 모드별 프레임-솔버 조합:
    - 'paper3ch': AnisotropicMultigroupFrame + FixedWeightSolver (논문 정합)
    - 'iso3'/'5ch': MultigroupFrame + RidgeWeightSolver (레거시)
    """
    def __init__(self, in_features: int, out_features: int,
                 channel_type: str = 'paper3ch', rho: Optional[float] = None):
        """
        Args:
            in_features (int): 입력 복소 텐서의 차원 수
            out_features (int): 출력 복소 텐서의 차원 수
            channel_type (str): 'paper3ch' (논문 기본), 'iso3', '5ch' (레거시)
            rho (float, optional): 릿지 정규화 계수 (레거시 채널 전용)
        """
        super(ResonantBlock, self).__init__()

        self.in_features = in_features
        self.out_features = out_features
        self.channel_type = channel_type
        self.rho = rho if rho is not None else R_CONST.RHO

        # 채널 모드에 따라 프레임과 솔버를 분기
        if channel_type == 'paper3ch':
            # 논문 정합: 비등방 프레임 + 고정 가중치 (eq:channel_set)
            self.frame = AnisotropicMultigroupFrame()
            self.solver = FixedWeightSolver(frame=self.frame)
        else:
            # 레거시: 등방 프레임 + 릿지 최적화
            self.frame = MultigroupFrame(channel_type=self.channel_type)
            self.solver = RidgeWeightSolver(frame=self.frame, rho=self.rho)

        # K개 방향성 로그 미분 추출기 (eq:central_diff)
        self.extractor = DirectionalLogDerivative(
            in_features=in_features,
            out_features=out_features,
            frame=self.frame
        )

        # 유효 로그기울기 어셈블러 L̂ = Σ w_k · L_k (eq:lhat)
        self.assembler = EffectiveLogDerivativeAssembler()

    def forward(self,
                z: torch.Tensor,
                sigma_prime: Optional[torch.Tensor] = None,
                tau_prime: Optional[torch.Tensor] = None,
                dynamic_rho: Optional[float] = None) -> torch.Tensor:
        """
        Args:
            z (torch.Tensor): 입력 복소 텐서 Z^{(l)} (Shape: [Batch, ..., in_features])
            sigma_prime (torch.Tensor, optional): 현재 상태의 실수부 위상 변화율 (노이즈 궤도)
            tau_prime (torch.Tensor, optional): 현재 상태의 허수부 위상 변화율 (진행 궤도)
            dynamic_rho (float, optional): 스윕(Sweep) 시 동적으로 변경할 릿지 계수

        Returns:
            torch.Tensor: 노이즈가 상쇄된 유효 로그기울기 텐서 L_G
                          (Shape: [Batch, ..., out_features], dtype: Complex128)
        """
        # [방어 기제 1] 입력 데이터 Z의 초고정밀도 강제
        if z.dtype != PrecisionManager.COMPLEX_DTYPE:
            z = z.to(dtype=PrecisionManager.COMPLEX_DTYPE)

        device = z.device

        # 디바이스 동기화 (GPU 이동 시 프레임/솔버 텐서도 함께 이동)
        if self.frame.device != device:
            self.frame.to(device)
            if self.channel_type == 'paper3ch':
                self.solver.to(device)
            else:
                self.solver = RidgeWeightSolver(frame=self.frame, rho=self.rho)

        # 궤도 정보(σ', τ') 누락 시 기본 정방향 궤도 자동 주입 (σ'=0, τ'=1)
        batch_shape = z.shape[:-1]

        if sigma_prime is None:
            sigma_prime = torch.zeros(batch_shape, dtype=PrecisionManager.REAL_DTYPE, device=device)
        elif sigma_prime.dtype != PrecisionManager.REAL_DTYPE:
            sigma_prime = sigma_prime.to(dtype=PrecisionManager.REAL_DTYPE, device=device)

        if tau_prime is None:
            tau_prime = torch.ones(batch_shape, dtype=PrecisionManager.REAL_DTYPE, device=device)
        elif tau_prime.dtype != PrecisionManager.REAL_DTYPE:
            tau_prime = tau_prime.to(dtype=PrecisionManager.REAL_DTYPE, device=device)

        # [Step 1] 방향성 로그 미분 D_k log f (eq:central_diff)
        # Shape: [Batch, ..., K, out_features]
        L_k = self.extractor(z)

        # [Step 2] 채널 가중치 w_k 산출 (paper3ch: 고정, 레거시: 릿지 최적화)
        # Shape: [Batch, ..., K]
        w_k = self.solver.solve(sigma_prime, tau_prime, dynamic_rho=dynamic_rho)

        # [Step 3] 유효 로그기울기 L̂ = Σ w_k · L_k (eq:lhat)
        # Shape: [Batch, ..., out_features]
        L_G = self.assembler(L_k, w_k)

        return L_G

    def extra_repr(self) -> str:
        return f'in_features={self.in_features}, out_features={self.out_features}, channel={self.channel_type}, ρ={self.rho:.1e}'


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Layers] Integrated Resonant Block Test ---")

    # 1. 하이퍼파라미터 설정 (Batch=2, Seq=4, in=8, out=16)
    batch_size, seq_len, in_dim, out_dim = 2, 4, 8, 16

    # 2. ResonantBlock 인스턴스화 (논문 정합 paper3ch 모드)
    res_block = ResonantBlock(in_features=in_dim, out_features=out_dim, channel_type='paper3ch')
    print(f"공명 블록 (paper3ch): \n{res_block}\n")

    # 3. 데이터 준비 (z: 복소 텐서, σ'/τ': 경로 텐서)
    z_real = torch.randn((batch_size, seq_len, in_dim), dtype=PrecisionManager.REAL_DTYPE)
    z_imag = torch.randn((batch_size, seq_len, in_dim), dtype=PrecisionManager.REAL_DTYPE)
    z_input = torch.complex(z_real, z_imag).requires_grad_(True)

    # 외부 게이지 상태 (Phase 4)를 모사
    sigma_p = torch.randn((batch_size, seq_len), dtype=PrecisionManager.REAL_DTYPE, requires_grad=True) * 0.1
    tau_p   = torch.ones((batch_size, seq_len), dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    print(f"입력 데이터 Z_in Shape:  {z_input.shape}")
    print(f"진행 방향 궤도 σ' Shape: {sigma_p.shape}, τ' Shape: {tau_p.shape}\n")

    # 4. Forward Pass: 데이터가 3개의 모듈을 관통하며 노이즈가 필터링 됨
    L_G_out = res_block(z_input, sigma_p, tau_p)

    print(f"출력 유효 궤도 L_G Shape: {L_G_out.shape} (기대값: [{batch_size}, {seq_len}, {out_dim}])")

    # 5. 차원 무결성 검증
    expected_shape = (batch_size, seq_len, out_dim)
    if L_G_out.shape == expected_shape:
        print("[OK] L_G 출력 차원 정합 확인")
    else:
        print(f"❌ 실패: 블록 출력 차원에 오류가 있습니다. 실제값: {L_G_out.shape}")

    # 6. 미분 역전파(Autograd) 다중 경로 보존 검증 (3-Way Gradient Flow)
    loss = L_G_out.abs().sum()
    loss.backward()

    if (z_input.grad is not None) and (sigma_p.grad is not None) and (tau_p.grad is not None):
        if not torch.isnan(z_input.grad).any() and not torch.isnan(sigma_p.grad).any():
            print("✅ 성공: 오차(Loss)가 Z(데이터), σ'(노이즈 궤도), τ'(진행 속도)의 세 갈래로 완벽히 역전파됩니다!\n")
        else:
            print("❌ 실패: 역전파 중 NaN 붕괴가 발생했습니다.\n")
    else:
        print("❌ 실패: 미분 역전파 연결이 어딘가 끊어졌습니다.\n")
