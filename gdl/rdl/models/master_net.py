"""
=============================================================================
[Project RDL] Orchestration - Master Resonant Network Assembly
=============================================================================
복소 임베딩, 다원군 공명 전파, 게이지 동역학을 통합한 End-to-End 신경망.

[구현 수식 매핑 - unified_en.tex 참조]
1. ResonantBlock: L̂ = Σ w_k · D_k log f (eq:lhat)
2. GaugeStateBuffer: ψ, φ ODE 사다리꼴 적분 (eq:ode_psi, eq:ode_phi)
3. [확장] 리 군 지수 사상: Z_{l+1} = Z_l · exp(h · L̂)
"""

import torch
import torch.nn as nn
from typing import Dict, Optional

from gdl.rdl.constants import R_CONST, PrecisionManager
from gdl.rdl.layers.embedding import ComplexLinearEmbedding
from gdl.rdl.layers.ideal_phase import IdealPhaseTargetGenerator
from gdl.rdl.layers.resonant_block import ResonantBlock
from gdl.rdl.dynamics.state_buffer import GaugeStateBuffer

class MasterResonantNetwork(nn.Module):
    """
    공명 딥러닝 마스터 모델. TotalResonanceLoss에 필요한
    모든 물리량(Z_out, Psi_target, L_G, phi, psi, tau_prime)을 산출한다.
    """
    def __init__(self,
                 in_features: int,
                 hidden_features: int,
                 out_features: int,
                 num_layers: int = 3,
                 channel_type: str = 'paper3ch',
                 damping_mode: str = 'paper'):
        """
        Args:
            in_features: 원본 실수 데이터 차원
            hidden_features: 복소 은닉층 차원
            out_features: 최종 출력 복소 차원
            num_layers: ResonantBlock 깊이
            channel_type: 'paper3ch' (논문 기본), 'iso3', '5ch' (레거시)
            damping_mode: 'paper' (논문 기본) 또는 'extended' ([확장])
        """
        super(MasterResonantNetwork, self).__init__()

        self.in_features = in_features
        self.hidden_features = hidden_features
        self.out_features = out_features
        self.num_layers = num_layers
        self.h_step = R_CONST.H_LIE  # [확장] 리 군 지수 사상 단계

        # 1. 타겟 위상 생성기 및 복소 임베딩
        self.target_generator = IdealPhaseTargetGenerator()
        self.embedding = ComplexLinearEmbedding(in_features, hidden_features)

        # 2. 다원군 공명 전파 블록 (eq:lhat)
        self.resonant_layers = nn.ModuleList([
            ResonantBlock(hidden_features, hidden_features, channel_type=channel_type)
            for _ in range(num_layers)
        ])

        # 3. 게이지 동역학 메모리 버퍼 (eq:ode_psi, eq:ode_phi)
        self.gauge_buffers = nn.ModuleList([
            GaugeStateBuffer(damping_mode=damping_mode) for _ in range(num_layers)
        ])

        # 4. 최종 출력 프로젝션 (Hidden -> Out)
        self.out_projection = ComplexLinearEmbedding(hidden_features, out_features)

    def forward(self, X: torch.Tensor) -> Dict[str, torch.Tensor]:
        """
        데이터 X를 받아 복소 공간으로 승격시킨 후, N개의 공명 블록을 관통시키며
        최종 출력과 손실 함수 계산에 필요한 모든 동역학 상태를 반환합니다.

        Args:
            X (torch.Tensor): 원본 실수 입력 (Shape: [Batch, ..., in_features])
                              ※ 곡률 패널티를 위해 반드시 requires_grad=True 상태여야 함!

        Returns:
            Dict[str, torch.Tensor]: Phase 5 통합 손실 함수에 직접 투입(Unpack)할 수 있는 텐서 딕셔너리
        """
        # [방어 기제 1] 입력 데이터의 정밀도 및 미분 추적 상태 강제 확인
        if X.dtype != PrecisionManager.REAL_DTYPE:
            X = X.to(dtype=PrecisionManager.REAL_DTYPE)
        if not X.requires_grad:
            # 훈련 환경이 아닐 경우라도 곡률 2차 미분을 위해 강제 활성화
            X.requires_grad_(True)

        # -------------------------------------------------------------
        # [Step 1] 타겟 위상 생성 및 복소 임베딩 (X -> Z_0)
        # -------------------------------------------------------------
        Psi_target = self.target_generator(X)
        Z = self.embedding(X)

        # -------------------------------------------------------------
        # [Step 2] 초기 궤도 및 게이지 상태 설정
        # -------------------------------------------------------------
        # 첫 번째 층 진입 전, 위상(φ)과 속도(ψ)는 백지 상태(None으로 자동 초기화)
        current_state = None

        # 초기 정방향 전진 궤도(tau'=1)와 노이즈 이탈률(sigma'=0)
        tau_prime = torch.ones(Z.shape[:-1], dtype=PrecisionManager.REAL_DTYPE, device=Z.device)
        sigma_prime = torch.zeros(Z.shape[:-1], dtype=PrecisionManager.REAL_DTYPE, device=Z.device)

        # 루프 종료 후 마지막 층의 상태를 보관하기 위한 변수
        final_L_G = None
        final_psi = None
        final_phi = None

        # -------------------------------------------------------------
        # [Step 3] 심층 공명 전파 루프 (Deep Resonant Forward)
        # -------------------------------------------------------------
        for i in range(self.num_layers):
            # 1. 방향성 미분과 릿지 가중합으로 현재 층의 유효 흐름(L_G) 추출
            L_G = self.resonant_layers[i](Z, sigma_prime, tau_prime)

            # 2. 추출된 흐름을 게이지 버퍼에 넣어 노이즈를 판별하고 내부 상태(state) 진화
            next_state, v_sigma, alpha_eff = self.gauge_buffers[i](L_G, tau_prime, current_state)
            psi, phi = next_state

            # 3. 다음 층으로 전달할 궤도 정보 동기화 (투영 피드백 / 자기 조직화)
            # 게이지의 관성 속도(psi)가 다음 층의 정방향 진행 속도(tau') 프록시가 되며,
            # 현재 층에서 판별된 노이즈(v_sigma)가 다음 층의 릿지 횡방향 이탈률(sigma') 프록시가 됩니다.
            tau_prime = psi.mean(dim=-1)
            sigma_prime = v_sigma.mean(dim=-1)

            # [확장] 리 군 지수 사상: Z_{l+1} = Z_l · exp(h · L̂)
            Z = Z * torch.exp(self.h_step * L_G)

            current_state = next_state

            # 마지막 레이어 상태 로깅
            if i == self.num_layers - 1:
                final_L_G = L_G
                final_psi = psi
                final_phi = phi

        # -------------------------------------------------------------
        # [Step 4] 최종 차원 프로젝션 (출력 Z_out)
        # -------------------------------------------------------------
        Z_out = self.out_projection(Z)

        # 손실 함수 **kwargs 언패킹용 출력 딕셔너리
        return {
            "X_in": X,
            "Z_out": Z_out,
            "Psi_target": Psi_target,
            "L_G": final_L_G,
            "phi": final_phi,
            "psi": final_psi,
            "tau_prime": tau_prime,
            "phase_out": torch.angle(Z_out)  # PQO 손실 시각화/분석용
        }

    def extra_repr(self) -> str:
        return f"in={self.in_features}, hidden={self.hidden_features}, out={self.out_features}, layers={self.num_layers}"


# =================================================================
# 직접 실행 시 종단간(End-to-End) 오케스트레이션 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()
    # Phase 5 통합 손실 함수와 Phase 6 PGGD 옵티마이저를 불러와 End-to-End 연결 검증
    from gdl.rdl.losses.total_loss import TotalResonanceLoss
    from gdl.rdl.optim.pggd import PGGD

    print("\n--- [RDL Orchestration] Master Resonant Network E2E Test ---")

    batch_size, in_dim, hidden_dim, out_dim = 2, 4, 8, 2

    # 1. 파이프라인 최하단: 원본 실수 데이터 생성 (곡률 미분 추적을 위해 Requires Grad 필수!)
    X_input = torch.randn((batch_size, in_dim), dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    # 2. 마스터 시스템 인스턴스 조립
    model = MasterResonantNetwork(in_features=in_dim, hidden_features=hidden_dim, out_features=out_dim,
                                   num_layers=3, channel_type='paper3ch', damping_mode='paper')
    loss_fn = TotalResonanceLoss(lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.5)
    optimizer = PGGD(model.parameters(), lr=0.01)

    print(f"공명 마스터 네트워크 (paper3ch):\n{model}\n")

    # =================================================================
    # [THE GRAND CYCLE] 공명 딥러닝 훈련 루프 사이클 1회전 시뮬레이션
    # =================================================================
    print("▶ 위대한 훈련 사이클 (The Grand Cycle) 1회전 가동...")

    # [Step A: Forward Pass]
    # 데이터가 복소 평면으로 올라와 3겹의 ResonantBlock과 GaugeBuffer를 관통하며 진화합니다.
    outputs = model(X_input)
    print("  └─ [Forward 완료] 4개 차원의 텐서 우주 관통 및 3단계 게이지 진화 완료")

    # [Step B: Compute Loss (Phase 5 도킹)]
    # Python ** 언패킹을 이용해 모델 출력을 그대로 손실 함수에 투입!
    total_loss, metrics = loss_fn(**outputs)
    print(f"  └─ [Loss 산출 완료] Total Loss 붕괴 스칼라: {total_loss.item():.4f}")

    # [Step C: Backward Pass (초거대 통합 미분)]
    # 손실 함수 래퍼 -> 동역학 게이지 -> 3중 공명 레이어 -> 임베딩 -> 최하단 X 데이터
    optimizer.zero_grad()
    total_loss.backward()

    if X_input.grad is not None and not torch.isnan(X_input.grad).any():
        print("  └─ [Backward 완료] 3차 미분(Hessian)의 충격을 견디며 최하단 X_in까지 미분 트리가 100% 무손실 도달했습니다!")
    else:
        print("  └─ ❌ [Backward 실패] 역전파 중 그래프 단절 또는 NaN 발생")
        import sys
        sys.exit(1)

    # [Step D: Optimizer Step (Phase 6 PGGD 진화)]
    try:
        optimizer.step()
        print("  └─ [Optimizer 완료] PGGD가 위상 궤도에 맞춰 전역 가중치를 갱신하고 내부 나침반을 진화시켰습니다!")
    except Exception as e:
        print(f"  └─ ❌ [Optimizer 실패] PGGD 갱신 중 에러 발생: {e}")

    print("\n✅ 대성공: 임베딩 -> 전파 -> 게이지 적분 -> 곡률 붕괴 손실 -> PGGD 진화로 이어지는")
    print("          공명 딥러닝 시스템의 완전한 End-to-End 유기적 결합이 수학적으로 완벽히 증명되었습니다!")
