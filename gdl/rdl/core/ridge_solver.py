"""
=============================================================================
[Project RDL] Core Math Primitives - Weight Solvers for Multigroup Channels
=============================================================================
다원군(Multigroup) 채널의 방향성 미분 정보를 결합하는 가중치(w_k) 솔버.

[구현 수식 매핑 - unified_en.tex 참조]
- eq:channel_set: 논문 사전 지정 가중치 w_k = (0.40, 0.35, 0.25)
- eq:lhat: L̂ = Σ w_k · D_k log f (가중 결합)

[릿지 솔버 - 레거시/과결정 채널용]
- b_j = [σ'(λ_j), τ'(λ_j), 1]^T (경로 접선 + 합=1 제약)
- M_j(:,k) = [v_{k,x}, v_{k,y}, 1]^T (채널 방향 설계 행렬)
- w = (M^T M + ρI)^{-1} M^T b (릿지 정규화 최소자승해)

[최적화 기법]
- M 행렬은 방향 각도(상수)로 고정되어 있으므로, 딥러닝 훈련 루프 내 병목을 없애기 위해
  초기화 시점에 솔버 행렬 S = (M^T M + ρI)^{-1} M^T 를 사전 연산(Pre-compute)하여 캐싱합니다.
- 배치(Batch) 연산 시 역행렬을 매번 구하지 않고 O(1) 텐서 행렬곱만 수행하여 GPU 효율을 극대화합니다.
"""

import torch
from typing import Optional

from gdl.rdl.constants import R_CONST, PrecisionManager
from gdl.rdl.core.multigroup import MultigroupFrame, AnisotropicMultigroupFrame


class FixedWeightSolver:
    """
    논문 사전 지정 가중치를 반환하는 솔버 (unified_en.tex, eq:channel_set).

    논문의 3-채널 구성은 가중치가 사전 결정되어 있으므로
    릿지 최적화가 불필요하다. 이 솔버는 AnisotropicMultigroupFrame과
    함께 사용하여 w_k = (0.40, 0.35, 0.25)를 즉시 반환한다.

    RidgeWeightSolver와 동일한 인터페이스를 제공하여 상위 모듈에서
    솔버를 교체하더라도 코드 변경이 최소화된다.
    """

    def __init__(self, frame: AnisotropicMultigroupFrame):
        """
        Args:
            frame (AnisotropicMultigroupFrame): 논문 정합 비등방 프레임
        """
        self.frame = frame
        self.device = frame.device
        self.K = frame.K
        self._w_k = frame.weights  # [K] 텐서, sum=1

    def solve(self, sigma_prime: torch.Tensor, tau_prime: torch.Tensor,
              dynamic_rho=None) -> torch.Tensor:
        """
        입력 경로 접선과 무관하게 논문 고정 가중치를 배치 크기에 맞춰 반환한다.

        Args:
            sigma_prime: [Batch, ...] (인터페이스 호환용, 실제 사용하지 않음)
            tau_prime: [Batch, ...] (인터페이스 호환용, 실제 사용하지 않음)
            dynamic_rho: 무시됨 (인터페이스 호환용)

        Returns:
            torch.Tensor: [Batch, ..., K] 논문 고정 가중치
        """
        original_shape = sigma_prime.shape
        # w_k를 입력 텐서의 배치 차원에 맞게 브로드캐스트
        expand_dims = list(original_shape) + [self.K]
        w = self._w_k.expand(expand_dims)
        return w

    def to(self, device: torch.device):
        """디바이스 이동"""
        self.device = device
        self._w_k = self._w_k.to(device)
        return self

    def __repr__(self):
        return f"<FixedWeightSolver: w_k={self._w_k.tolist()}, device={self.device}>"


class RidgeWeightSolver:
    """
    다원군 프레임의 방향 미분 가중치를 배치(Batched) 단위로 초고속 계산하는 클래스.
    """
    def __init__(self, frame: MultigroupFrame, rho: float = None):
        """
        Args:
            frame (MultigroupFrame): 다원군 프레임 인스턴스 (방향 벡터 v_k 포함)
            rho (float): 릿지 정규화 계수. 기본값은 R_CONST.RHO 사용.
        """
        self.frame = frame
        self.device = frame.device
        self.K = frame.K

        # [방어 기제] 3채널(iso3)은 미지수와 제약식이 3개로 완벽히 일치하여
        # 완벽한 닫힌형(Closed-form) 해가 존재하므로 ρ=0.0으로 강제합니다.
        if frame.channel_type == 'iso3':
            self.rho = 0.0
        else:
            self.rho = rho if rho is not None else R_CONST.RHO

        self._build_matrices()

    def _build_matrices(self):
        """
        최소자승법을 위한 설계 행렬 M을 조립하고, O(1) 호출을 위한 솔버 행렬을 사전 연산합니다.
        """
        # 1. 채널 행렬 M 조립 (Shape: [3, K])
        # 1행: v_x, 2행: v_y, 3행: 1 (가중치 총합=1 제약조건)
        v_x, v_y = self.frame.v_k_components
        ones = torch.ones_like(v_x)
        self.M = torch.stack([v_x, v_y, ones], dim=0)

        # 2. M^T 구성 (Shape: [K, 3])
        self.M_T = self.M.t()

        # 3. M^T M 구성 (Shape: [K, K])
        self.M_T_M = torch.matmul(self.M_T, self.M)

        # 4. 식별 행렬 I (Shape: [K, K])
        self.I = torch.eye(self.K, dtype=PrecisionManager.REAL_DTYPE, device=self.device)

        # 5. 기본 A 행렬 (M^T M + ρI)
        self.A = self.M_T_M + self.rho * self.I

        # [핵심 최적화] S = A^{-1} M^T 를 미리 구하여 캐싱 (Shape: [K, 3])
        # 이후 w = S * b 만으로 가중치를 즉시 산출할 수 있습니다.
        # 수치적 안정성을 위해 .inverse() 대신 torch.linalg.solve 활용
        try:
            self.solver_matrix = torch.linalg.solve(self.A, self.M_T)
        except torch._C._LinAlgError:
            # 완벽한 0의 rho로 인해 특이행렬(Singular) 붕괴 시 안전 마진 추가 후 복구
            fallback_rho = R_CONST.EPSILON
            A_fallback = self.M_T_M + fallback_rho * self.I
            self.solver_matrix = torch.linalg.solve(A_fallback, self.M_T)

        # 역전파(Autograd)에서 이 행렬들이 가중치 업데이트 대상이 되지 않도록 고정
        self.M.requires_grad_(False)
        self.M_T.requires_grad_(False)
        self.A.requires_grad_(False)
        self.solver_matrix.requires_grad_(False)

    def solve(self, sigma_prime: torch.Tensor, tau_prime: torch.Tensor, dynamic_rho: Optional[float] = None) -> torch.Tensor:
        """
        주어진 배치(Batch)의 경로 접선 방향(sigma', tau')에 대해 최적의 채널 가중치 w를 계산합니다.

        Args:
            sigma_prime (torch.Tensor): 경로의 실수부 위상 변화율 (Shape: [Batch, ...])
            tau_prime (torch.Tensor): 경로의 허수부 위상 변화율 (Shape: [Batch, ...])
            dynamic_rho (float, optional): 기본 ρ 외에 스윕(Sweep) 등으로 동적 릿지 계수를 부여할 때 사용

        Returns:
            torch.Tensor: 최적화된 다원군 가중치 w (Shape: [Batch, ..., K])
        """
        original_shape = sigma_prime.shape

        # 텐서 평탄화 (Batch 단위 병렬 처리를 위해 1D로 통일)
        sig_p = sigma_prime.contiguous().view(-1)
        tau_p = tau_prime.contiguous().view(-1)
        B = sig_p.shape[0]

        # 1. 타겟 벡터 b 구성 (Shape: [B, 3, 1])
        # b_j = [sigma'_j, tau'_j, 1]^T
        ones_b = torch.ones_like(sig_p)
        b = torch.stack([sig_p, tau_p, ones_b], dim=1).unsqueeze(-1)

        if dynamic_rho is None or dynamic_rho == self.rho:
            # [O(1) 캐싱 연산] solver_matrix [K, 3] @ b [B, 3, 1] -> 브로드캐스팅 -> [B, K, 1]
            w_unsqueeze = torch.matmul(self.solver_matrix, b)
        else:
            # 동적 ρ(rho)가 주어졌을 때는 실시간으로 시스템 방정식을 풉니다.
            A_dynamic = self.M_T_M + dynamic_rho * self.I
            A_batch = A_dynamic.unsqueeze(0).expand(B, -1, -1)   # [B, K, K]

            M_T_batch = self.M_T.unsqueeze(0).expand(B, -1, -1)  # [B, K, 3]
            M_T_b = torch.matmul(M_T_batch, b)                   # [B, K, 1]

            # A * w = M^T b -> w = A^{-1} M^T b
            w_unsqueeze = torch.linalg.solve(A_batch, M_T_b)

        # 2. [B, K, 1] -> [B, K] -> 원본 Shape 복원 (예: [Batch, Seq, K])
        w = w_unsqueeze.squeeze(-1).view(*original_shape, self.K)

        return w

# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    PrecisionManager.setup_precision()

    print("\n--- [RDL Core Math] Ridge Weight Solver Test ---")

    # 1. 테스트용 입력 데이터 생성 (Batch=2)
    # Batch 1: σ'=0.0 (이상적 정렬), τ'=1.0 (정방향 진행)
    # Batch 2: σ'=0.5 (위상 곡률 발생), τ'=-1.0 (미러/역방향 진행)
    sigma_prime = torch.tensor([0.0, 0.5], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)
    tau_prime   = torch.tensor([1.0, -1.0], dtype=PrecisionManager.REAL_DTYPE, requires_grad=True)

    # ---------------------------------------------------------
    # [Test A] iso3 채널 (K=3, ρ=0) - 완벽 해 도출
    # ---------------------------------------------------------
    frame_iso3 = MultigroupFrame(channel_type='iso3')
    solver_iso3 = RidgeWeightSolver(frame_iso3)
    w_iso3 = solver_iso3.solve(sigma_prime, tau_prime)

    print("\n[Test A] iso3 Channel (ρ=0.0) 가중치 w 결과:")
    print(w_iso3.detach().cpu().numpy())

    sum_w_iso3 = w_iso3.sum(dim=-1)
    print(f"  └─ 가중치 총합 (기대값 1.0): {sum_w_iso3.detach().cpu().numpy()}")

    # 방향성 복원 검증 (M * w ≈ b 확인)
    # w_iso3: [2, 3] -> [2, 1, 3], M^T: [3, 3]
    b_recon_iso3 = torch.matmul(w_iso3.unsqueeze(1), solver_iso3.M.t()).squeeze(1)
    print(f"  └─ 복원된 σ' (기대값 [0.0, 0.5]): {b_recon_iso3[:, 0].detach().cpu().numpy()}")
    print(f"  └─ 복원된 τ' (기대값 [1.0, -1.0]): {b_recon_iso3[:, 1].detach().cpu().numpy()}")

    # ---------------------------------------------------------
    # [Test B] 5ch 채널 (K=5, ρ=3e-3) - 릿지 정규화
    # ---------------------------------------------------------
    frame_5ch = MultigroupFrame(channel_type='5ch')
    solver_5ch = RidgeWeightSolver(frame_5ch)
    w_5ch = solver_5ch.solve(sigma_prime, tau_prime)

    print(f"\n[Test B] 5ch Channel (ρ={R_CONST.RHO}) 가중치 w 결과:")
    print(w_5ch.detach().cpu().numpy())

    sum_w_5ch = w_5ch.sum(dim=-1)
    print(f"  └─ 가중치 총합 (기대값 ~1.0): {sum_w_5ch.detach().cpu().numpy()}")

    # ---------------------------------------------------------
    # [Test C] 미분(Autograd) 연결 테스트
    # ---------------------------------------------------------
    loss = w_5ch.sum()
    loss.backward()

    if sigma_prime.grad is not None and tau_prime.grad is not None:
        print("\n✅ 성공: 사전 연산(Pre-compute) 최적화를 적용했음에도 Autograd 미분 그래프가 완벽히 보존되었습니다!")
    else:
        print("\n❌ 실패: 미분 그래프가 단절되었습니다.")
