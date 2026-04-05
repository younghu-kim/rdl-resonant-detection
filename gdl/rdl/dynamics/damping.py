"""
=============================================================================
[Project RDL] Gauge State Dynamics - Damping Rate Calculators
=============================================================================
게이지 상태 갱신에 사용할 감쇠율(alpha)을 계산하는 모듈.

[구현 수식 매핑 - unified_en.tex 참조]

1. PaperDampingRate (논문 정합):
   - eq:ode_alpha: α_j = Δt_j / (τ_g + Δt_j)
   - 고정 적분 단계 H_GAUGE를 사용하는 단순 감쇠

2. ExtendedDampingRate [확장: 논문 외 독자적 기여]:
   - α_eff = [|τ'|/(τ_g+|τ'|)] × [1/(1+κ|Re(e^{-iφ}L_G)|)]
   - 1항: 경로 속도 비례 기본 감쇠
   - 2항: 위상 어긋남에 따른 상태 의존 감쇠 (κ 계수)
"""

import torch
import torch.nn as nn
from typing import Optional

from gdl.rdl.constants import R_CONST, PrecisionManager


class PaperDampingRate(nn.Module):
    """
    논문 정합 감쇠율: α_j = Δt_j / (τ_g + Δt_j)

    논문 참조: unified_en.tex, eq:ode_alpha

    고정 적분 단계(H_GAUGE)와 감쇠 상수(τ_g)만으로 결정되는
    상수 감쇠율. 경로 속도나 내부 위상에 의존하지 않는다.
    """

    def __init__(self, tau_g: Optional[float] = None, dt: Optional[float] = None):
        """
        Args:
            tau_g: 감쇠 시간 상수. None이면 R_CONST.TAU_G (=0.04)
            dt: 적분 단계 Δt_j. None이면 R_CONST.H_GAUGE (=0.04)
        """
        super(PaperDampingRate, self).__init__()
        self.tau_g = tau_g if tau_g is not None else R_CONST.TAU_G
        self.dt = dt if dt is not None else R_CONST.H_GAUGE
        # 사전 계산: α = Δt / (τ_g + Δt)
        self._alpha = self.dt / (self.tau_g + self.dt)

    def forward(self, tau_prime: torch.Tensor, phi: torch.Tensor,
                L_G: torch.Tensor) -> torch.Tensor:
        """
        인터페이스는 ExtendedDampingRate과 동일하나, tau_prime/phi/L_G를 사용하지 않고
        상수 α_j를 phi와 동일한 shape으로 브로드캐스트하여 반환한다.

        Returns:
            torch.Tensor: 상수 감쇠율 α_j (phi와 동일 shape)
        """
        return torch.full_like(phi, self._alpha)

    def extra_repr(self) -> str:
        return f'τ_g={self.tau_g:.3f}, Δt={self.dt:.3f}, α={self._alpha:.4f} (논문 eq:ode_alpha)'


class ExtendedDampingRate(nn.Module):
    """
    [확장: 논문 외 독자적 기여] 상태 의존 적응형 감쇠율.

    α_eff = [|τ'|/(τ_g+|τ'|)] × [1/(1+κ|Re(e^{-iφ}L_G)|)]

    1항: 경로 속도 비례 기본 감쇠
    2항: 위상-신호 어긋남에 따른 방어적 감쇠 (κ 계수)
    """
    def __init__(self, tau_g: Optional[float] = None, kappa: Optional[float] = None):
        super(ExtendedDampingRate, self).__init__()

        # 상수 매핑 (수학적 불변성 유지)
        self.tau_g = tau_g if tau_g is not None else R_CONST.TAU_G
        self.kappa = kappa if kappa is not None else R_CONST.KAPPA
        self.eps   = R_CONST.EPSILON

    def _safe_abs(self, x: torch.Tensor) -> torch.Tensor:
        """
        절댓값 |x|의 x=0 에서 발생하는 미분 꺾임(NaN 폭주)을 막기 위해
        C^∞(무한 번 미분 가능) 형태의 평활화된(Euclidean norm) 절댓값을 반환합니다.
        """
        return torch.sqrt(x**2 + self.eps)

    def forward(self, tau_prime: torch.Tensor, phi: torch.Tensor, L_G: torch.Tensor) -> torch.Tensor:
        """
        Args:
            tau_prime (torch.Tensor): 현재 배치의 경로 진행 속도 (Shape: [Batch, ...], Real)
            phi (torch.Tensor): 게이지 버퍼에 저장된 현재 내부 순수 위상 상태 (Shape: [Batch, ..., Features], Real)
            L_G (torch.Tensor): [Phase 3]에서 넘어온 유효 로그기울기 (Shape: [Batch, ..., Features], Complex)

        Returns:
            torch.Tensor: 계산된 적응형 감쇠율 α_eff (Shape: [Batch, ..., Features], Real)
                          값의 범위는 항상 0 <= α_eff < 1 을 완벽히 보장합니다.
        """
        # [방어 기제 1] 입력 텐서 정밀도 확인 및 승격
        if tau_prime.dtype != PrecisionManager.REAL_DTYPE:
            tau_prime = tau_prime.to(dtype=PrecisionManager.REAL_DTYPE)
        if phi.dtype != PrecisionManager.REAL_DTYPE:
            phi = phi.to(dtype=PrecisionManager.REAL_DTYPE)
        if L_G.dtype != PrecisionManager.COMPLEX_DTYPE:
            L_G = L_G.to(dtype=PrecisionManager.COMPLEX_DTYPE)

        # Broadcasting을 위해 차원이 다를 경우 자동으로 확장 (예: [Batch, Seq] -> [Batch, Seq, 1])
        while tau_prime.dim() < phi.dim():
            tau_prime = tau_prime.unsqueeze(-1)

        # -------------------------------------------------------------------
        # [Step 1] Base Damping Factor 계산: |τ'| / (τ_g + |τ'|)
        # 진행 속도에 비례하는 기본 관성 필터 (항상 0 ~ 1 미만의 값을 가짐)
        # -------------------------------------------------------------------
        abs_tau = self._safe_abs(tau_prime)
        # tau_g가 전역 상수로 양수이므로 분모가 0이 될 위험은 없으나, 절대 붕괴 방지를 위해 eps가 내장된 safe_abs 사용
        base_damping = abs_tau / (self.tau_g + abs_tau)

        # -------------------------------------------------------------------
        # [Step 2] State Factor (어긋남 페널티) 계산: 1 / (1 + κ * |Re(e^{-iφ} L_G)|)
        # -------------------------------------------------------------------
        # [극강의 최적화] e^{-iφ} L_G 연산을 복소수 텐서 생성 오버헤드 없이,
        # 오일러 공식을 통해 대수적으로 전개하여 O(1) 실수 연산으로 치환합니다.
        # Re( (cos(φ) - i*sin(φ)) * (Re(L_G) + i*Im(L_G)) )
        # = Re(L_G)*cos(φ) + Im(L_G)*sin(φ)
        cos_phi = torch.cos(phi)
        sin_phi = torch.sin(phi)

        # O(1) 초고속 기하 내적 연산 (복소수 곱셈 메모리 할당 제거)
        real_projection = L_G.real * cos_phi + L_G.imag * sin_phi

        # 절대값(abs) 취득 시 특이점 붕괴 방지를 위해 Epsilon-safe abs 사용
        safe_abs_projection = self._safe_abs(real_projection)

        state_factor = 1.0 / (1.0 + self.kappa * safe_abs_projection)

        # -------------------------------------------------------------------
        # [Step 3] 최종 적응형 감쇠율 조합 (α_eff)
        # -------------------------------------------------------------------
        alpha_eff = base_damping * state_factor

        return alpha_eff

    def extra_repr(self) -> str:
        return f'τ_g={self.tau_g:.3f}, κ={self.kappa:.2f}, eps={self.eps:.1e} [확장]'


# 하위 호환 별칭
AdaptiveDampingRate = ExtendedDampingRate


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Dynamics] Adaptive Damping Rate Test ---")

    # 1. 텐서 설정 (Batch=3, Features=4)
    batch_size, out_dim = 3, 4

    # 2. 극한의 상황을 모사한 데이터 입력
    # Batch 0: 정상 진행 (τ'=1.0), 노이즈 없음 (L_G가 거의 0에 수렴) -> 엑셀(높은 α_eff)
    # Batch 1: 정상 진행 (τ'=1.0), 노이즈 극심 (L_G가 거대하게 요동침) -> 브레이크(낮은 α_eff)
    # Batch 2: 정지 궤도 (τ'=0.0) -> 정보 차단 (α_eff = 0)
    tau_p = torch.tensor([[1.0], [1.0], [0.0]], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    # 내부 기억 위상(phi): 모두 0 라디안으로 고정
    phi = torch.zeros((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    # 외부 입력 L_G
    L_G_real = torch.tensor([[0.01]*out_dim, [100.0]*out_dim, [10.0]*out_dim], dtype=PrecisionManager.REAL_DTYPE)
    L_G_imag = torch.tensor([[0.01]*out_dim, [100.0]*out_dim, [10.0]*out_dim], dtype=PrecisionManager.REAL_DTYPE)
    L_G = torch.complex(L_G_real, L_G_imag).requires_grad_(True)

    # 3. 감쇠율 연산기 구동
    damping_calc = AdaptiveDampingRate()
    alpha_eff = damping_calc(tau_p, phi, L_G)

    print("계산된 적응형 감쇠율 (α_eff):")
    print(f"  └─ Batch 0 (τ'=1, 안정된 궤도)  : {alpha_eff[0, 0].item():.4f} (기대값: 최대치 근접 / 적극 수용)")
    print(f"  └─ Batch 1 (τ'=1, 극심한 노이즈): {alpha_eff[1, 0].item():.4f} (기대값: 0.0 근접 / 방어 기제 발동)")
    print(f"  └─ Batch 2 (τ'=0, 궤도 정지)    : {alpha_eff[2, 0].item():.4f} (기대값: 0.0 수렴 / 업데이트 차단)\n")

    # 4. 물리적 제약조건 (0 <= α < 1) 검증
    is_safe = torch.all((alpha_eff >= 0.0) & (alpha_eff < 1.0))
    if is_safe:
        print("✅ 성공: 모든 감쇠율이 [0, 1)의 절대 물리적 안전 지대 안에 완벽히 갇혀 발산을 억제합니다.")
    else:
        print("❌ 실패: 감쇠율이 1을 초과하여 발산했거나 0 미만으로 역전되었습니다.")

    # 5. 역전파 (3-Way Autograd) 검증 (특이점 미분 붕괴 방어)
    loss = alpha_eff.sum()
    loss.backward()

    if tau_p.grad is not None and phi.grad is not None and L_G.grad is not None:
        if not torch.isnan(tau_p.grad).any() and not torch.isnan(phi.grad).any() and not torch.isnan(L_G.grad).any():
            print("✅ 성공: L_G, phi, tau' 로 찢어지는 3-Way 역전파 그래프가 절댓값 미분 특이점에서도 꺾임 없이 100% 보존되었습니다!\n")
        else:
            print("❌ 실패: 역전파 중 절댓값 미분 꺾임에 의해 NaN이 발생했습니다.\n")
    else:
        print("❌ 실패: 미분 그래프가 끊어졌습니다.\n")
