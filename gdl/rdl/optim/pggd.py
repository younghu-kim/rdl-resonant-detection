"""
=============================================================================
[Project RDL] PGGD (Phase-Gauge Gradient Descent) Optimizer
=============================================================================
[확장: 논문 외 독자적 기여] 게이지 ODE (eq:gauge_ode) 기반 최적화기.

판별식 회전을 통한 위상 정렬 경사 필터링:
1. 상태 초기화: 가중치별 Float64 위상(φ_opt) 및 속도(ψ_opt) 할당
2. 경사 필터링: Im{e^{-iφ}(g - ψ)} → 직교 노이즈만 추출하여 유효 경사 도출
3. 게이지 적분: 사다리꼴 이산화(eq:ode_psi, eq:ode_phi)로 내부 상태 진화

적분 단계: H_GAUGE (eq:ode_alpha의 Δt_j)
"""

import torch
from torch.optim import Optimizer
from typing import Iterable, Optional

from gdl.rdl.constants import R_CONST, PrecisionManager

class PGGD(Optimizer):
    """
    [확장] Phase-Gauge Gradient Descent (위상-게이지 경사 하강법).

    게이지 ODE 사다리꼴 적분으로 내부 위상 상태를 진화시키며,
    판별식 회전(Im{e^{-iφ}·Δg})으로 노이즈를 필터링한 유효 경사를 적용한다.
    """
    def __init__(self,
                 params: Iterable[torch.Tensor],
                 lr: float = 1e-3,
                 tau_g: Optional[float] = None,
                 kappa: Optional[float] = None,
                 h: Optional[float] = None,
                 eta_phi: float = 0.0):
        """
        Args:
            params: 최적화 대상 파라미터
            lr: 학습률
            tau_g: 속도 관성 (기본값 R_CONST.TAU_G)
            kappa: 확장 감쇠 민감도 (기본값 R_CONST.KAPPA)
            h: 게이지 적분 스텝 (기본값 R_CONST.H_GAUGE)
            eta_phi: 복원 포텐셜 강도. V(φ) = -eta_phi·cos(φ)의 기울기
                     eta_phi·sin(φ)를 위상 갱신에 빼서 에르고딕 확산을 방지한다.
                     0.0이면 비활성 (기존 동작과 동일).
        """
        if lr < 0.0:
            raise ValueError(f"[PGGD Error] 학습률(lr)은 0.0 이상이어야 합니다: {lr}")

        defaults = dict(
            lr=lr,
            tau_g=tau_g if tau_g is not None else R_CONST.TAU_G,
            kappa=kappa if kappa is not None else R_CONST.KAPPA,
            h=h if h is not None else R_CONST.H_GAUGE,
            eta_phi=eta_phi
        )
        super(PGGD, self).__init__(params, defaults)

    @torch.no_grad()
    def step(self, closure=None) -> Optional[float]:
        """ 단일 최적화 스텝 수행 (가중치 갱신 및 내부 게이지 적분 진화) """
        loss = None
        if closure is not None:
            with torch.enable_grad():
                loss = closure()

        for group in self.param_groups:
            lr = group['lr']
            tau_g = group['tau_g']  # 속도 관성 갱신율 (Adam의 1-beta1 역할 대체)
            h_step = group['h']     # 이산 적분 스텝 크기
            eta_phi = group['eta_phi']  # 복원 포텐셜 강도 (Kuramoto형 결합)

            for p in group['params']:
                if p.grad is None:
                    continue

                grad = p.grad
                if grad.is_sparse:
                    raise RuntimeError("❌ [PGGD Error] 희소(Sparse) 기울기는 지원하지 않습니다.")

                state = self.state[p]

                # =========================================================
                # [Phase 6.1] 상태 초기화 (Step 0)
                # =========================================================
                if len(state) == 0:
                    state['step'] = 0
                    # 옵티마이저 내부 메모리는 누적 오차 방지를 위해 무조건 Float64로 강제 승격
                    state['phi_opt'] = torch.zeros_like(p, memory_format=torch.preserve_format, dtype=PrecisionManager.REAL_DTYPE)
                    state['psi_opt'] = torch.zeros_like(p, memory_format=torch.preserve_format, dtype=PrecisionManager.REAL_DTYPE)

                state['step'] += 1
                phi_opt = state['phi_opt']
                psi_opt = state['psi_opt']

                # =========================================================
                # [Phase 6.2] 판별식 회전을 통한 위상 정렬 경사 필터링
                # =========================================================
                # Re{e^{-iφ}·(g - ψ + i·g_imag)} = cos(φ)·(g-ψ) + sin(φ)·g_imag
                # 정렬된 유효 경사 성분을 추출한다 (직교 노이즈 제거).
                if grad.is_complex():
                    g_real = grad.real.to(dtype=PrecisionManager.REAL_DTYPE)
                    g_imag = grad.imag.to(dtype=PrecisionManager.REAL_DTYPE)
                else:
                    g_real = grad.to(dtype=PrecisionManager.REAL_DTYPE)
                    g_imag = torch.zeros_like(g_real)

                diff_real = g_real - psi_opt
                diff_imag = g_imag

                # Re{e^{-iφ}·(diff_real + i·diff_imag)} 오일러 전개
                sin_phi = torch.sin(phi_opt)
                cos_phi = torch.cos(phi_opt)
                tilde_g_t = cos_phi * diff_real + sin_phi * diff_imag

                # =========================================================
                # [Phase 6.3] 모델 가중치 갱신 및 Heun 게이지 적분
                # =========================================================

                # 1. 모델 가중치(θ) 갱신 (θ_{t+1} = θ_t - η * tilde_g_t)
                # 모델 파라미터(p)의 원래 정밀도(예: Float32)로 다시 캐스팅하여 적용 (정밀도 충돌 방어)
                p.add_(tilde_g_t.to(dtype=p.dtype), alpha=-lr)

                # 2. 속도 평활화 (eq:ode_psi 적용)
                v_t = g_real
                psi_next = (1.0 - tau_g) * psi_opt + tau_g * v_t

                # 3. 사다리꼴 위상 적분 (eq:ode_phi 적용)
                phi_next = phi_opt + 0.5 * h_step * (psi_opt + psi_next)

                # 3-b. 복원 포텐셜: -∇V(φ) = -eta_phi·sin(φ)
                # V(φ) = -eta_phi·cos(φ) → φ=0 안정 평형, φ=±π 불안정 평형
                # Weyl 등분포에 의한 경사 소멸을 방지 (Kuramoto형 결합)
                if eta_phi > 0.0:
                    phi_next = phi_next - eta_phi * torch.sin(phi_opt)

                # 4. 위상 래핑 (-π ~ π)
                phi_next_wrapped = torch.atan2(torch.sin(phi_next), torch.cos(phi_next))

                # 5. 내부 메모리(State Dict) 안전 덮어쓰기 (In-place copy_)
                # 메모리 주소를 유지하여 파이토치 옵티마이저의 VRAM 파편화 방지
                state['psi_opt'].copy_(psi_next)
                state['phi_opt'].copy_(phi_next_wrapped)

        return loss

    def extra_repr(self) -> str:
        return "Phase-Gauge Gradient Descent (PGGD) Optimizer - Final"


# =================================================================
# 직접 실행 시 모듈 무결성 테스트 (Sanity Check)
# =================================================================
if __name__ == "__main__":
    import math
    PrecisionManager.setup_precision()
    import torch.nn as nn

    print("\n--- [RDL Optim] PGGD Final Integration Test (Block 6.3) ---")

    # 1. 테스트 설정 (초기 가중치 10.0)
    dummy_weight = nn.Parameter(torch.tensor([10.0], dtype=torch.float32))
    lr_test = 0.1
    tau_g_test = 0.5
    h_test = 1.0

    optimizer = PGGD([dummy_weight], lr=lr_test, tau_g=tau_g_test, h=h_test)

    # 2. 상태 수동 세팅을 위한 헛스텝 (Step 0)
    dummy_weight.grad = torch.tensor([10.0], dtype=torch.float32)
    optimizer.step()
    state = optimizer.state[dummy_weight]

    # [시나리오 D: Heun 적분을 통한 옵티마이저의 능동적 진화 증명]
    # 위상(φ) = 90도 (π/2), 속도(ψ) = 0.0 으로 강제 세팅 (궤도 직교 이탈 상태)
    # 일반 기울기(g_t) = 10.0 입력

    state['phi_opt'].fill_(math.pi / 2)
    state['psi_opt'].fill_(0.0)
    dummy_weight.data.fill_(10.0)
    dummy_weight.grad = torch.tensor([10.0], dtype=torch.float32)

    # ==========================
    # 수동 수학적 예측 (Hand-calculation)
    # diff_real = 10.0 - 0.0 = 10.0
    # tilde_g_t = cos(90)*0 - sin(90)*10.0 = -10.0
    # 파라미터 갱신 = 10.0 - (0.1 * -10.0) = 11.0
    # ψ_next = (1 - 0.5)*0.0 + 0.5*(10.0) = 5.0
    # φ_next = (π/2) + 0.5 * 1.0 * (0.0 + 5.0) = 1.5708 + 2.5 = 4.0708 라디안
    # 위상 래핑: 4.0708 > π 이므로 -2.2124 라디안으로 래핑 됨
    # ==========================

    # Step 구동!
    optimizer.step()

    w_out = dummy_weight.data.item()
    psi_out = state['psi_opt'].item()
    phi_out = state['phi_opt'].item()

    print("▶ PGGD 옵티마이저 진화 스텝 결과:")
    print(f"  └─ 모델 파라미터(θ) 갱신: {w_out:.4f} (기대값: 11.0000)")
    print(f"  └─ 내부 속도(ψ) 갱신:    {psi_out:.4f} (기대값: 5.0000)")
    print(f"  └─ 내부 위상(φ) 갱신:    {phi_out:.4f} (기대값: -2.2124, Heun 적분 & 래핑 완료)")

    success = math.isclose(w_out, 11.0, abs_tol=1e-4) and \
              math.isclose(psi_out, 5.0, abs_tol=1e-4) and \
              math.isclose(phi_out, -2.2124, abs_tol=1e-3)

    if success:
        print("\n✅ 대성공: 파라미터가 갱신됨과 동시에, 옵티마이저 내부의 '위상 나침반(Gauge)'이 ")
        print("          Heun 사다리꼴 적분을 통해 완벽하게 미래 시점으로 진화(Evolve)했습니다!")
    else:
        print("\n❌ 실패: 가중치 업데이트 또는 게이지 적분 계산에 오차가 발생했습니다.")

    # 3. VRAM 메모리 누수 방벽 점검
    ptr_phi = state['phi_opt'].data_ptr()
    dummy_weight.grad = torch.tensor([5.0], dtype=torch.float32)
    optimizer.step() # Step 2

    if state['phi_opt'].data_ptr() == ptr_phi:
        print("✅ 성공: State Dictionary의 텐서 메모리 주소가 완벽히 보존되었습니다. (In-place Copy 성공, VRAM 누수 없음)")
    else:
        print("❌ 실패: 매 스텝마다 새로운 텐서가 할당되어 메모리 누수(OOM) 위험이 있습니다.")
