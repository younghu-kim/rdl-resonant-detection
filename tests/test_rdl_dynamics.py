"""
=============================================================================
[Test] RDL Gauge Dynamics
=============================================================================
SignalSeparator, PaperDampingRate, ExtendedDampingRate,
HeunIntegrator, GaugeStateBuffer 모듈 테스트.
"""

import pytest
import torch

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.dynamics.separator import SignalSeparator
from gdl.rdl.dynamics.damping import PaperDampingRate, ExtendedDampingRate
from gdl.rdl.dynamics.integrator import HeunIntegrator
from gdl.rdl.dynamics.state_buffer import GaugeStateBuffer

BATCH = 3
FEATURES = 4
DTYPE_R = PrecisionManager.REAL_DTYPE
DTYPE_C = PrecisionManager.COMPLEX_DTYPE


def _make_complex(shape):
    """테스트용 복소 텐서 생성."""
    real = torch.randn(shape, dtype=DTYPE_R)
    imag = torch.randn(shape, dtype=DTYPE_R)
    return torch.complex(real, imag)


class TestSignalSeparator:
    """SignalSeparator: L_G -> (v_t, v_sigma) 분리."""

    def test_signal_separator(self):
        """복소 L_G에서 실수 v_t, v_sigma를 올바르게 추출해야 한다."""
        separator = SignalSeparator()
        L_G = _make_complex((BATCH, FEATURES))

        v_t, v_sigma = separator(L_G)

        # 실수 타입 검증
        assert not v_t.is_complex(), f"v_t가 복소수: {v_t.dtype}"
        assert not v_sigma.is_complex(), f"v_sigma가 복소수: {v_sigma.dtype}"

        # shape 보존 검증
        assert v_t.shape == (BATCH, FEATURES)
        assert v_sigma.shape == (BATCH, FEATURES)

        # 값 정합 검증
        assert torch.allclose(v_t, L_G.real), "v_t != Re(L_G)"
        assert torch.allclose(v_sigma, L_G.imag), "v_sigma != Im(L_G)"


class TestPaperDampingRate:
    """PaperDampingRate: 상수 감쇠율."""

    def test_paper_damping_rate(self):
        """출력이 (0, 1) 범위의 상수여야 한다."""
        damping = PaperDampingRate()

        tau_prime = torch.ones((BATCH,), dtype=DTYPE_R)
        phi = torch.zeros((BATCH, FEATURES), dtype=DTYPE_R)
        L_G = _make_complex((BATCH, FEATURES))

        alpha = damping(tau_prime, phi, L_G)

        assert alpha.shape == phi.shape, f"shape 불일치: {alpha.shape}"
        assert torch.all(alpha > 0.0), f"감쇠율이 0 이하: {alpha.min()}"
        assert torch.all(alpha < 1.0), f"감쇠율이 1 이상: {alpha.max()}"
        assert not torch.isnan(alpha).any(), "NaN 발생"


class TestExtendedDampingRate:
    """ExtendedDampingRate: 상태 의존 적응형 감쇠율."""

    def test_extended_damping_rate(self):
        """출력이 유효 범위 [0, 1) 안에 있어야 한다."""
        damping = ExtendedDampingRate()

        tau_prime = torch.tensor([[1.0], [0.5], [0.0]], dtype=DTYPE_R)
        phi = torch.randn((BATCH, FEATURES), dtype=DTYPE_R)
        L_G = _make_complex((BATCH, FEATURES))

        alpha = damping(tau_prime, phi, L_G)

        assert alpha.shape == phi.shape, f"shape 불일치: {alpha.shape}"
        assert torch.all(alpha >= 0.0), f"감쇠율이 0 미만: {alpha.min()}"
        assert torch.all(alpha < 1.0), f"감쇠율이 1 이상: {alpha.max()}"
        assert not torch.isnan(alpha).any(), "NaN 발생"


class TestHeunIntegrator:
    """HeunIntegrator: 사다리꼴 적분기."""

    def test_heun_integrator_shape(self):
        """출력 (psi_next, phi_next)의 shape이 입력과 동일해야 한다."""
        integrator = HeunIntegrator()

        psi_t = torch.randn((BATCH, FEATURES), dtype=DTYPE_R)
        phi_t = torch.randn((BATCH, FEATURES), dtype=DTYPE_R)
        v_t = torch.randn((BATCH, FEATURES), dtype=DTYPE_R)
        tau_prime = torch.ones((BATCH,), dtype=DTYPE_R)
        alpha_eff = torch.full((BATCH, FEATURES), 0.5, dtype=DTYPE_R)

        psi_next, phi_next = integrator(psi_t, phi_t, v_t, tau_prime, alpha_eff)

        assert psi_next.shape == psi_t.shape, \
            f"psi shape 불일치: {psi_next.shape} != {psi_t.shape}"
        assert phi_next.shape == phi_t.shape, \
            f"phi shape 불일치: {phi_next.shape} != {phi_t.shape}"

        # NaN 없음
        assert not torch.isnan(psi_next).any(), "psi_next에 NaN"
        assert not torch.isnan(phi_next).any(), "phi_next에 NaN"


class TestGaugeStateBuffer:
    """GaugeStateBuffer: 게이지 ODE 메모리 셀."""

    def test_gauge_buffer_initial_state(self):
        """state=None일 때 자동으로 0으로 초기화되어야 한다."""
        buffer = GaugeStateBuffer(damping_mode='paper')
        L_G = _make_complex((BATCH, FEATURES))

        # state=None으로 호출
        (psi_next, phi_next), v_sigma, alpha_eff = buffer(L_G, state=None)

        # 출력이 유효한 텐서여야 함
        assert psi_next.shape == (BATCH, FEATURES)
        assert phi_next.shape == (BATCH, FEATURES)
        assert v_sigma.shape == (BATCH, FEATURES)

        # NaN 없음
        assert not torch.isnan(psi_next).any(), "psi_next에 NaN"
        assert not torch.isnan(phi_next).any(), "phi_next에 NaN"
        assert not torch.isnan(v_sigma).any(), "v_sigma에 NaN"
        assert not torch.isnan(alpha_eff).any(), "alpha_eff에 NaN"
