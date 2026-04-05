"""
=============================================================================
[Project RDL] Total Resonance Loss - 4항 통합 손실 래퍼
=============================================================================
공명 딥러닝의 4가지 손실 조건을 가중합하는 마스터 래퍼.

[구현 수식 매핑 - unified_en.tex 참조]
1. Discriminant Loss (L_res):  eq:F2, Im{e^{-iφ}(L̂ - φ̇)} = 0
2. Curvature Loss (L_curv):    Gate B, v_σ=0 (A'=0, A''<0)
3. Target Loss (L_tgt):        출력 위상 arg(Z)과 Ψ(n) 매칭
4. PQO Loss (L_pqo):           eq:pqo, sin²(L·arg(Z)/2) → 격자 정렬

L_total = λ_res·L_res + λ_curv·L_curv + λ_tgt·L_tgt + λ_pqo·L_pqo
"""

import torch
import torch.nn as nn
from typing import Tuple, Dict, Optional

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.losses.phase_loss import DiscriminantLoss
from gdl.rdl.losses.curvature_loss import CurvatureStabilizationLoss
from gdl.rdl.losses.target_loss import TargetPhaseMatchingLoss
from gdl.rdl.losses.quantisation_loss import PhaseQuantisationLoss

class TotalResonanceLoss(nn.Module):
    """
    4항 통합 손실: L_total = λ_res·L_res + λ_curv·L_curv + λ_tgt·L_tgt + λ_pqo·L_pqo
    """
    def __init__(self,
                 lambda_res: float = 1.0,
                 lambda_curv: float = 0.5,
                 lambda_tgt: float = 1.0,
                 lambda_pqo: float = 0.0,
                 lambda_a: float = 1.0,
                 lambda_b: float = 1.0,
                 pqo_L: int = None,
                 pqo_mode: str = 'sin2',
                 reduction: str = 'mean'):
        """
        Args:
            lambda_res: 판별식 손실 가중치 (eq:F2)
            lambda_curv: 곡률 손실 가중치 (Gate B)
            lambda_tgt: 타겟 위상 매칭 손실 가중치
            lambda_pqo: 위상 양자화 손실 가중치 (eq:pqo). 0.0이면 비활성 (하위 호환)
            lambda_a: 곡률 1차 미분 페널티 (내부)
            lambda_b: 곡률 2차 미분 페널티 (내부)
            pqo_L: PQO 양자화 레벨. None이면 R_CONST.PQO_L_DEFAULT
            pqo_mode: 'sin2' (원래 Gate A) 또는 'cos2' (Maslov 보정).
                      cos2는 실수값 함수 영점 탐지에 사용 — F₂=0과 PQO=0 동시 달성.
            reduction: 배치 축소 방식
        """
        super(TotalResonanceLoss, self).__init__()

        self.lambda_res = lambda_res
        self.lambda_curv = lambda_curv
        self.lambda_tgt = lambda_tgt
        self.lambda_pqo = lambda_pqo

        self.loss_phase = DiscriminantLoss(reduction=reduction)
        self.loss_curvature = CurvatureStabilizationLoss(lambda_a=lambda_a, lambda_b=lambda_b)
        self.loss_target = TargetPhaseMatchingLoss(reduction=reduction)
        self.loss_pqo = PhaseQuantisationLoss(L=pqo_L, reduction=reduction, mode=pqo_mode)

    def forward(self,
                Z_out: torch.Tensor,
                X_in: torch.Tensor,
                Psi_target: torch.Tensor,
                L_G: torch.Tensor,
                phi: torch.Tensor,
                psi: torch.Tensor,
                tau_prime: Optional[torch.Tensor] = None,
                **kwargs) -> Tuple[torch.Tensor, Dict[str, float]]:
        """
        Args:
            Z_out (torch.Tensor): 네트워크 최종 출력 복소 텐서 (Shape: [Batch, ...])
            X_in (torch.Tensor): 모델 초입 원본 실수 입력 (requires_grad=True 필수, 곡률 2차 미분용)
            Psi_target (torch.Tensor): 이상 위상 타겟 (Shape: [Batch, ...])
            L_G (torch.Tensor): 추출된 유효 복소 흐름 속도 (Shape: [Batch, ...])
            phi (torch.Tensor): 게이지 버퍼의 내부 누적 위상 상태 (Shape: [Batch, ...])
            psi (torch.Tensor): 게이지 버퍼의 위상 속도 관성 (Shape: [Batch, ...])
            tau_prime (torch.Tensor, optional): 현재 경로의 진행 궤도 속도 (Shape: [Batch, ...])

        Returns:
            Tuple:
                - total_loss (torch.Tensor): 역전파(backward)를 위한 스칼라 텐서
                - loss_dict (Dict): 터미널 로깅 및 모니터링을 위한 항목별 Loss (스칼라 값)
        """
        # [방어 기제 1] X_in 곡률 미분(Hessian) 그래프 생성 확인 (곡률 활성 시에만)
        if self.lambda_curv > 0 and not X_in.requires_grad:
            raise RuntimeError("❌ [TotalLoss Error] 곡률 2차 미분을 위해 X_in 텐서는 반드시 requires_grad=True 여야 합니다.")

        # -------------------------------------------------------------
        # 1. 위상 정렬 잔차 손실 (L_res)
        # -------------------------------------------------------------
        l_res = self.loss_phase(L_G, phi, psi, tau_prime)

        # -------------------------------------------------------------
        # 2. 곡률 극점 안정화 손실 (L_curv)
        # -------------------------------------------------------------
        l_curv = self.loss_curvature(Z_out, X_in) if self.lambda_curv > 0 else torch.tensor(0.0, device=Z_out.device)

        # -------------------------------------------------------------
        # 3. 목표 위상 매칭 손실 (L_tgt)
        # -------------------------------------------------------------
        l_tgt = self.loss_target(Z_out, Psi_target)

        # -------------------------------------------------------------
        # 4. 위상 양자화 손실 (eq:pqo)
        # -------------------------------------------------------------
        l_pqo = self.loss_pqo(Z_out) if self.lambda_pqo > 0 else torch.tensor(0.0, device=Z_out.device)

        # -------------------------------------------------------------
        # 5. 통합 가중합
        # -------------------------------------------------------------
        total_loss = (self.lambda_res * l_res + self.lambda_curv * l_curv
                      + self.lambda_tgt * l_tgt + self.lambda_pqo * l_pqo)

        loss_dict = {
            'L_res': l_res.item(),
            'L_curv': l_curv.item(),
            'L_tgt': l_tgt.item(),
            'L_pqo': l_pqo.item() if isinstance(l_pqo, torch.Tensor) else 0.0,
            'Total': total_loss.item()
        }

        return total_loss, loss_dict

    def extra_repr(self) -> str:
        return (f"λ_res={self.lambda_res}, λ_curv={self.lambda_curv}, "
                f"λ_tgt={self.lambda_tgt}, λ_pqo={self.lambda_pqo}")


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()
    print("\n--- [RDL Losses] Total Resonance Loss Wrapper Test ---")

    batch_size, out_dim = 2, 4

    # 1. 전역 시스템 입출력 모사 (requires_grad=True 필수 텐서들)
    # 모델 최하단 입력 데이터
    X_in = torch.randn((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    # 모델 최상단 출력 Z
    # 곡률 미분을 위해 X_in과 수학적 연관성을 맺어줌 (Z = -X_in^2 + 0.1i)
    Z_real = 10.0 - X_in**2
    Z_imag = torch.ones_like(X_in) * 0.1
    Z_out = torch.complex(Z_real, Z_imag)

    # 목표 및 동역학 상태
    Psi_target = torch.randn((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)
    L_G = torch.complex(
        torch.randn((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE),
        torch.randn((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE)
    ).requires_grad_(True)
    phi = torch.randn((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)
    psi = torch.randn((batch_size, out_dim), dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)
    tau_prime = torch.ones((batch_size, 1), dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    # 2. 통합 손실 함수 구동
    total_loss_fn = TotalResonanceLoss(lambda_res=1.0, lambda_curv=0.5, lambda_tgt=1.0)
    print(f"생성된 마스터 통합 로스: \n{total_loss_fn}\n")

    total_loss, metrics = total_loss_fn(Z_out, X_in, Psi_target, L_G, phi, psi, tau_prime)

    print("▶ 통합 공명 붕괴 손실 딕셔너리:")
    for k, v in metrics.items():
        print(f"  └─ {k:<7}: {v:.4f}")

    if metrics['Total'] > 0:
        print("\n✅ 성공: 세 가지 물리적/기하학적 페널티가 단일 스칼라 에너지로 완벽히 융합(Wrap)되었습니다.")
    else:
        print("\n❌ 실패: 가중합 로직에 오류가 발생했습니다.")

    # 3. 4항 통합 역전파 보존 검증
    # 이 단 한 번의 backward()로 3차 미분 곡률 트리를 포함해 7방향으로 미분이 뻗어나갑니다.
    total_loss.backward()

    all_tensors = [X_in, Psi_target, phi, psi, tau_prime, L_G]
    names = ["X_in(데이터, 곡률)", "Psi_tgt(타겟)", "phi(위상)", "psi(속도)", "tau'(궤도)", "L_G(흐름)"]

    all_safe = True
    print("\n▶ 통합 7-Way 초거대 역전파 생존 검증:")
    for name, t in zip(names, all_tensors):
        if t.grad is None:
            print(f"  ❌ 실패: {name} 텐서로 향하는 미분 그래프 단절")
            all_safe = False
        elif torch.isnan(t.grad).any():
            print(f"  ❌ 실패: {name} 미분 중 NaN 붕괴 발생")
            all_safe = False
        else:
            print(f"  ✅ {name:<17}: 정상 도달")

    if all_safe:
        print("\n✅ 대성공: 위상 회전, 타겟 매칭, 그리고 가장 파괴되기 쉬운 3차 미분(Hessian) 연산이")
        print("          동시에 맞물렸음에도 역전파 트리가 한 방울의 누수 없이 100% 생존했습니다!")
