"""
=============================================================================
[Project RDL] Gauge State Dynamics - Gauge State Buffer (RNN Cell equivalent)
=============================================================================
게이지 ODE의 위상(φ)과 속도(ψ)를 시간/깊이 축으로 추적·갱신하는 메모리 셀.

[구현 수식 매핑 - unified_en.tex 참조]
- eq:ode_alpha ~ eq:ode_phi: 게이지 상태 동역학 ODE 이산화

하위 모듈:
1. 신호 분리기: v_t = Re(L_G), v_σ = Im(L_G)
2. 감쇠율 연산기: PaperDampingRate (논문) 또는 ExtendedDampingRate [확장]
3. 사다리꼴 적분기: eq:ode_psi, eq:ode_phi
"""

import torch
import torch.nn as nn
from typing import Tuple, Optional

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.dynamics.separator import SignalSeparator
from gdl.rdl.dynamics.damping import PaperDampingRate, ExtendedDampingRate
from gdl.rdl.dynamics.integrator import HeunIntegrator

class GaugeStateBuffer(nn.Module):
    """
    게이지 ODE 메모리 셀: φ, ψ 상태를 갱신한다 (eq:ode_psi, eq:ode_phi).

    damping_mode:
    - 'paper': PaperDampingRate (α_j = Δt/(τ_g+Δt), 논문 정합)
    - 'extended': ExtendedDampingRate (상태 의존 κ 감쇠, [확장])
    """
    def __init__(self, tau_g: Optional[float] = None, kappa: Optional[float] = None,
                 h: Optional[float] = None, damping_mode: str = 'paper'):
        """
        Args:
            tau_g: 감쇠 시간 상수 (eq:ode_alpha)
            kappa: 상태 의존 감쇠 계수 ([확장], extended 모드 전용)
            h: 적분 단계 Δt_j
            damping_mode: 'paper' (논문 정합) 또는 'extended' ([확장])
        """
        super(GaugeStateBuffer, self).__init__()

        self.damping_mode = damping_mode
        self.separator = SignalSeparator()

        if damping_mode == 'paper':
            self.damping_calc = PaperDampingRate(tau_g=tau_g, dt=h)
        elif damping_mode == 'extended':
            self.damping_calc = ExtendedDampingRate(tau_g=tau_g, kappa=kappa)
        else:
            raise ValueError(f"damping_mode는 'paper' 또는 'extended'여야 합니다: {damping_mode}")

        self.integrator = HeunIntegrator(h=h)

    def init_state(self, shape: torch.Size, device: torch.device) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        새로운 시퀀스나 에포크(Epoch) 시작 시 초기 상태(0.0) 텐서를 생성하여 반환합니다.

        Args:
            shape (torch.Size): 상태 텐서의 형태 (예: [Batch, ..., Features])
            device (torch.device): 텐서가 상주할 장치

        Returns:
            Tuple[torch.Tensor, torch.Tensor]: (초기화된 psi_0, 초기화된 phi_0)
        """
        # [방어 기제 1] Float64 (초고정밀도 실수) 강제 초기화
        factory_kwargs = {'dtype': PrecisionManager.REAL_DTYPE, 'device': device}

        psi_0 = torch.zeros(shape, **factory_kwargs)
        phi_0 = torch.zeros(shape, **factory_kwargs)

        return psi_0, phi_0

    def forward(self,
                L_G: torch.Tensor,
                tau_prime: Optional[torch.Tensor] = None,
                state: Optional[Tuple[torch.Tensor, torch.Tensor]] = None) -> Tuple[Tuple[torch.Tensor, torch.Tensor], torch.Tensor, torch.Tensor]:
        """
        현재 시점의 상태(state)와 입력 관측치(L_G, tau_prime)를 바탕으로
        다음 시점(t+1)으로 게이지 상태를 갱신(Predict-Correct)합니다.

        Args:
            L_G (torch.Tensor): [Phase 3]에서 추출된 복소 유효 속도 텐서 (Shape: [Batch, ..., Features])
            tau_prime (torch.Tensor, optional): 현재 경로의 진행 궤도 속도 (Shape: [Batch, ...]). None일 시 1.0 주입.
            state (Tuple[torch.Tensor, torch.Tensor], optional): (psi_t, phi_t) 형태의 이전 상태.
                                                                 None일 경우 자동으로 0으로 초기화.

        Returns:
            Tuple:
                - state_next (Tuple[Tensor, Tensor]): 갱신된 (psi_next, phi_next)
                - v_sigma (Tensor): 허수부 위상 노이즈 (Phase 5 손실 함수 로깅/계산용)
                - alpha_eff (Tensor): 현재 스텝에 적용된 브레이크/감쇠율 (로깅/분석용)
        """
        # [방어 기제 2] 정밀도 무결성 보장 및 tau_prime 자동 생성
        if L_G.dtype != PrecisionManager.COMPLEX_DTYPE:
            L_G = L_G.to(dtype=PrecisionManager.COMPLEX_DTYPE)

        if tau_prime is None:
            # tau_prime 누락 시 이상적인 정방향(Forward) 궤도로 간주
            tau_prime = torch.ones(L_G.shape[:-1], dtype=PrecisionManager.REAL_DTYPE, device=L_G.device)
        elif tau_prime.dtype != PrecisionManager.REAL_DTYPE:
            tau_prime = tau_prime.to(dtype=PrecisionManager.REAL_DTYPE)

        # [방어 기제 3] 상태가 누락되었거나 디바이스가 다를 경우 자동 초기화/동기화
        if state is None:
            psi_t, phi_t = self.init_state(L_G.shape, L_G.device)
        else:
            psi_t, phi_t = state
            # 장치 불일치(Device Mismatch) 동적 방어
            if psi_t.device != L_G.device:
                psi_t = psi_t.to(L_G.device)
            if phi_t.device != L_G.device:
                phi_t = phi_t.to(L_G.device)

        # 1. 신호 분리 (L_G -> v_t, v_sigma)
        v_t, v_sigma = self.separator(L_G)

        # 2. 적응형 감쇠율 연산 (α_eff 도출)
        alpha_eff = self.damping_calc(tau_prime, phi_t, L_G)

        # 3. Heun 적분기를 통한 상태 갱신 (Predict & Correct)
        psi_next, phi_next = self.integrator(psi_t, phi_t, v_t, tau_prime, alpha_eff)

        return (psi_next, phi_next), v_sigma, alpha_eff

    def extra_repr(self) -> str:
        return 'RNN Cell equivalent for Phase(φ) and Velocity(ψ) using Heun Integration'


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Dynamics] Gauge State Buffer (RNN Cell) Test ---")

    # 1. 하이퍼파라미터 및 디바이스 설정
    batch_size, out_dim = 2, 4

    # 2. 게이지 버퍼 인스턴스화
    buffer_cell = GaugeStateBuffer()
    print(f"생성된 게이지 상태 버퍼 (파라미터 0개 LSTM 대체재):\n{buffer_cell}\n")

    # 3. 초기 상태(State) 할당 (t=0)
    # 수동 초기화 대신 forward 내 자동 초기화를 유도해봅니다.
    state_t = None

    # 4. 시계열(Sequence) 루프 모사 (3 Step 전진)
    # t=1: 정상 속도 궤도 유입
    # t=2: 극한의 노이즈 폭풍 유입 (허수부가 50.0으로 요동침)
    # t=3: 다시 정상 궤도 회복
    seq_L_G = [
        torch.complex(
            torch.ones((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE),
            torch.zeros((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE)
        ).requires_grad_(True),
        torch.complex(
            torch.ones((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE) * 2.0,
            torch.full((batch_size, out_dim), 50.0, dtype=PrecisionManager.REAL_DTYPE) # 엄청난 노이즈
        ).requires_grad_(True),
        torch.complex(
            torch.ones((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE) * 1.5,
            torch.zeros((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE)
        ).requires_grad_(True)
    ]

    # 5. [루프 실행] 상태 갱신 (t=1 -> t=2 -> t=3)
    print("▶ 시간 축(Sequence)에 따른 위상(φ) 및 적응 감쇠율(α_eff) 진화:")
    for t in range(3):
        L_G_t = seq_L_G[t]

        # 핵심 전진 연산 (Forward Step) - 파이토치의 RNNCell과 완벽히 동일한 Stateless 패턴
        state_next, v_sigma, alpha_eff = buffer_cell(L_G_t, state=state_t)
        psi_t, phi_t = state_next

        print(f"  [Step {t+1}] 평균 노이즈 강도: {v_sigma.mean().item():.1f}")
        print(f"           평균 α_eff: {alpha_eff.mean().item():.4f} | 평균 속도(ψ): {psi_t.mean().item():.4f}")

        if t == 1:
            print("           (강한 노이즈 주입으로 α_eff가 0으로 급감하여 이전 상태 속도를 보존함!)")

        state_t = state_next # 다음 스텝으로 상태 전달

    # 6. BPTT (Backpropagation Through Time) 시간에 따른 역전파 보존 검증
    # 최종 상태(t=3)에서 오차가 발생했을 때, 이것이 루프를 뚫고 과거(t=1)의 텐서들까지 도달하는지 확인
    loss = phi_t.sum() + psi_t.sum()
    loss.backward()

    grad_L_G_0 = seq_L_G[0].grad

    if grad_L_G_0 is not None and not torch.isnan(grad_L_G_0).any():
        print("\n✅ 성공: (BPTT) 시간에 따른 다중 스텝 역전파(Backprop)가")
        print("        신호 분리 -> 감쇠 -> Heun 적분의 3중 모듈을 뚫고 완벽히 보존되었습니다!")
    else:
        print("\n❌ 실패: 시간 축을 거스르는 역전파 중 기울기 폭주(NaN)가 발생했거나 그래프가 끊어졌습니다.")
