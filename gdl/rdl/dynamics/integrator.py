"""
=============================================================================
[Project RDL] Gauge State Dynamics - Trapezoidal Discrete Integrator
=============================================================================
게이지 ODE의 위상(φ)과 위상 속도(ψ)를 사다리꼴 이산화로 갱신하는 적분기.

[구현 수식 매핑 - unified_en.tex 참조]
- eq:ode_psi: ψ_{j+1} = (1 - α_j) · ψ_j + α_j · v_t  (지수 이동 평균 속도 평활화)
- eq:ode_phi: φ_{j+1} = φ_j + (Δt/2) · (ψ_j + ψ_{j+1})  (사다리꼴 적분)

적분 단계 h = R_CONST.H_GAUGE (eq:ode_alpha의 Δt_j)

* Out-of-place 텐서 연산으로 Autograd 그래프를 보존합니다.
"""

import math
import torch
import torch.nn as nn
from typing import Tuple, Optional

from gdl.rdl.constants import R_CONST, PrecisionManager

class HeunIntegrator(nn.Module):
    """
    게이지 ODE 사다리꼴 적분기 (eq:ode_psi, eq:ode_phi).

    무상태(Stateless) 수학적 적분기로, 내부에 학습 파라미터가 없다.
    """
    def __init__(self, h: Optional[float] = None):
        """
        Args:
            h (float, optional): 적분 단계 Δt_j. None이면 R_CONST.H_GAUGE 사용.
        """
        super(HeunIntegrator, self).__init__()
        self.h = h if h is not None else R_CONST.H_GAUGE

    def forward(self,
                psi_t: torch.Tensor,
                phi_t: torch.Tensor,
                v_t: torch.Tensor,
                tau_prime: torch.Tensor,
                alpha_eff: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Args:
            psi_t (torch.Tensor): 현재 시점의 내부 위상 속도 상태 (Shape: [Batch, ..., Features])
            phi_t (torch.Tensor): 현재 시점의 내부 위상 상태 (Shape: [Batch, ..., Features])
            v_t (torch.Tensor): 현재 들어온 새로운 신호 속도 (Re(L_G)) (Shape: [Batch, ..., Features])
            tau_prime (torch.Tensor): 현재 경로의 진행 속도 궤도 (Shape: [Batch, ...])
            alpha_eff (torch.Tensor): [Block 4.2]에서 도출된 적응형 감쇠율 (Shape: [Batch, ..., Features])

        Returns:
            Tuple[torch.Tensor, torch.Tensor]:
                - psi_next: 다음 시점으로 갱신된 위상 속도 상태
                - phi_next: 다음 시점으로 갱신된 위상 상태 (Heun 적분 보정 완료 및 -π~π 래핑)
        """
        # [방어 기제 1] 모든 연산이 Float64(실수부)에서만 수행되는지 강제 점검
        # 게이지 엔진의 상태 메모리는 적분이 누적되므로 절대 단 1비트의 정밀도 손실도 허용해선 안 됩니다.
        tensors = [psi_t, phi_t, v_t, tau_prime, alpha_eff]
        for i, t in enumerate(tensors):
            if t.dtype != PrecisionManager.REAL_DTYPE:
                tensors[i] = t.to(dtype=PrecisionManager.REAL_DTYPE)
        psi_t, phi_t, v_t, tau_prime, alpha_eff = tensors

        # [방어 기제 2] Broadcasting 차원 확장
        # tau_prime은 Feature 차원이 없으므로, (Batch, 1) 형태로 맞춰주어
        # Feature 차원(phi, psi)과 충돌 없이 요소별 곱셈(Element-wise)이 되도록 조율합니다.
        while tau_prime.dim() < psi_t.dim():
            tau_prime = tau_prime.unsqueeze(-1)

        # [Step 1] 속도 평활화 (eq:ode_psi)
        # ψ_{j+1} = (1 - α_j) · ψ_j + α_j · v_t
        psi_next = (1.0 - alpha_eff) * psi_t + alpha_eff * v_t

        # [Step 2] 사다리꼴 위상 적분 (eq:ode_phi)
        # φ_{j+1} = φ_j + (Δt/2) · (ψ_j + ψ_{j+1})
        phi_next = phi_t + 0.5 * self.h * (psi_t + psi_next) * tau_prime

        # -------------------------------------------------------------
        # [Step 3] 위상(Phase) 텐서의 무한 팽창(Overflow) 방지 래핑
        # -------------------------------------------------------------
        # 위상각은 2π를 주기로 회전하므로, 무한히 커지면 Float64의 한계에 부딪혀 수치 불안정을 야기합니다.
        # atan2(sin, cos)를 활용하여 역전파(Autograd) C^∞ 미분 그래프를 온전히 보존한 채로
        # 각도를 -π ~ +π 구간으로 안전하게 매핑(Wrapping)합니다.
        phi_next_wrapped = torch.atan2(torch.sin(phi_next), torch.cos(phi_next))

        return psi_next, phi_next_wrapped

    def extra_repr(self) -> str:
        return f"Method=Heun(Predictor-Corrector), h={self.h:.2e}"


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Dynamics] Heun's Discrete Integrator Test ---")

    # 1. 텐서 설정 (Batch=1, Features=2)
    batch_size, out_dim = 1, 2

    # 2. 초기 상태 (t=0)
    # 현재 위상: 0 라디안, 현재 위상 속도: 1.0 rad/s
    phi_t = torch.tensor([[0.0, 0.0]], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)
    psi_t = torch.tensor([[1.0, 1.0]], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    # 3. 입력 관측치 (t=1)
    # 새로운 위상 속도 관측치: 3.0 rad/s
    # 경로 진행 속도(τ'): 1.0
    v_t = torch.tensor([[3.0, 3.0]], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)
    tau_p = torch.tensor([[1.0]], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    # 4. 적응형 감쇠율 테스트 (α_eff = 0.5 가정)
    # 기존 속도(1.0)와 관측 속도(3.0)의 정확히 중간인 2.0으로 갱신되어야 함
    alpha_eff = torch.tensor([[0.5, 0.5]], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    # 5. 적분기 구동 (테스트 가시성을 위해 극단적인 h=1.0 강제 주입)
    test_h = 1.0
    integrator = HeunIntegrator(h=test_h)
    psi_next, phi_next = integrator(psi_t, phi_t, v_t, tau_p, alpha_eff)

    print(f"초기 상태 (t=0):   φ={phi_t[0,0].item():.1f}, ψ={psi_t[0,0].item():.1f}")
    print(f"관측치    (t=1):   v_t={v_t[0,0].item():.1f}, α_eff={alpha_eff[0,0].item():.1f}, τ'={tau_p[0,0].item():.1f}")

    print(f"\n갱신된 위상 속도 (ψ_next): {psi_next[0,0].item():.4f} (기대값: 2.0000)")
    # 오일러 적분 기대값: φ = 0 + 1.0(h) * 1.0(현재속도) * 1(τ') = 1.0
    # Heun 적분 기대값: φ = 0 + 0.5(h) * (1.0 + 2.0) * 1(τ') = 1.5
    print(f"갱신된 누적 위상 (φ_next): {phi_next[0,0].item():.4f} (기대값: 1.5000 / Heun 사다리꼴 적분)\n")

    # 6. 주기성 래핑(Phase Wrapping) C^∞ 미분 보존 검증
    # phi_t를 3.14 (약 π)로 두고 업데이트하면 -π 쪽으로 부드럽게 순환해야 함
    phi_t_wrap = torch.tensor([[math.pi - 0.1, 0.0]], dtype=PrecisionManager.REAL_DTYPE)
    _, phi_next_wrap = integrator(psi_t, phi_t_wrap, v_t, tau_p, alpha_eff)

    if phi_next_wrap[0, 0].item() < 0:
        print(f"✅ 성공: 위상각(Phase)이 π를 초과하자 즉시 -π 쪽으로 순환(Wrapping)하여 메모리 팽창을 방어합니다. (결과: {phi_next_wrap[0, 0].item():.4f})")
    else:
        print("❌ 실패: 위상각 래핑 기제가 작동하지 않아 Overflow 위험이 있습니다.")

    # 7. 5-Way 역전파(Autograd) 보존 검증 (다중 입력 특이점 붕괴 방어)
    loss = psi_next.sum() + phi_next.sum()
    loss.backward()

    if all(t.grad is not None for t in [psi_t, phi_t, v_t, tau_p, alpha_eff]):
        if not any(torch.isnan(t.grad).any() for t in [psi_t, phi_t, v_t, tau_p, alpha_eff]):
            print("✅ 성공: 5개의 텐서(과거 2, 현재 3) 모두로 찢어지는 5-Way 역전파가 단절 없이 100% 보존되었습니다!\n")
        else:
            print("❌ 실패: 역전파 중 NaN이 발생했습니다.\n")
    else:
        print("❌ 실패: 미분 그래프가 어딘가 끊어졌습니다.\n")
