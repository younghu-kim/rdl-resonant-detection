"""
=============================================================================
[Test] RDL Core Math Modules
=============================================================================
complex_ops, multigroup, ridge_solver 모듈의 단위 테스트.
"""

import math
import pytest
import torch

from gdl.rdl.constants import R_CONST, PrecisionManager
from gdl.rdl.core.complex_ops import get_phase_amp, safe_complex_div, safe_complex_log
from gdl.rdl.core.multigroup import MultigroupFrame, AnisotropicMultigroupFrame
from gdl.rdl.core.ridge_solver import FixedWeightSolver, RidgeWeightSolver


# ------------------------------------------------------------------
# complex_ops 테스트
# ------------------------------------------------------------------

class TestGetPhaseAmp:
    """get_phase_amp: 복소 텐서에서 위상/진폭 추출."""

    def test_get_phase_amp(self):
        """위상은 [-pi, pi], 진폭은 >= 0 이어야 한다."""
        z = torch.tensor(
            [3.0 + 4.0j, -1.0 + 0.0j, 0.0 - 2.0j, 1.0 + 1.0j],
            dtype=PrecisionManager.COMPLEX_DTYPE,
        )
        phase, amp = get_phase_amp(z)

        # 위상 범위 검증
        assert torch.all(phase >= -math.pi), f"위상이 -pi 미만: {phase}"
        assert torch.all(phase <= math.pi), f"위상이 pi 초과: {phase}"

        # 진폭 비음수 검증
        assert torch.all(amp >= 0.0), f"진폭이 음수: {amp}"

        # NaN 없음
        assert not torch.isnan(phase).any()
        assert not torch.isnan(amp).any()


class TestSafeComplexDiv:
    """safe_complex_div: 영분모 방지 복소 나눗셈."""

    def test_safe_complex_div_no_nan(self):
        """분모가 0에 가까워도 NaN이 발생하지 않아야 한다."""
        num = torch.tensor(
            [1.0 + 1.0j, 2.0 + 3.0j],
            dtype=PrecisionManager.COMPLEX_DTYPE,
        )
        den = torch.tensor(
            [1e-300 + 1e-300j, 0.0 + 0.0j],
            dtype=PrecisionManager.COMPLEX_DTYPE,
        )
        result = safe_complex_div(num, den)

        assert not torch.isnan(result).any(), f"NaN 발생: {result}"
        assert not torch.isinf(result.real).any(), f"Inf 발생 (실수부): {result}"
        assert not torch.isinf(result.imag).any(), f"Inf 발생 (허수부): {result}"


class TestSafeComplexLog:
    """safe_complex_log: 영점 근처에서도 안전한 복소 로그."""

    def test_safe_complex_log_no_nan(self):
        """입력이 0에 가까워도 NaN이 발생하지 않아야 한다."""
        z = torch.tensor(
            [1e-300 + 0.0j, 0.0 + 0.0j, 1.0 + 0.0j],
            dtype=PrecisionManager.COMPLEX_DTYPE,
        )
        result = safe_complex_log(z)

        assert not torch.isnan(result).any(), f"NaN 발생: {result}"
        assert not torch.isinf(result.real).any(), f"Inf 발생 (실수부): {result}"


# ------------------------------------------------------------------
# multigroup 테스트
# ------------------------------------------------------------------

class TestMultigroupFrame:
    """MultigroupFrame: 등방 다원군 프레임 생성."""

    def test_multigroup_frame_paper3ch(self):
        """AnisotropicMultigroupFrame('paper3ch')은 3채널, 논문 각도(-30, 15, 60)를 사용한다."""
        frame = AnisotropicMultigroupFrame()

        assert frame.K == 3, f"채널 수가 3이 아님: {frame.K}"
        assert frame.channel_type == 'paper3ch'

        # 논문 각도 검증 (도 -> 라디안)
        expected_angles_deg = (-30.0, 15.0, 60.0)
        expected_radians = [math.radians(d) for d in expected_angles_deg]

        # delta_k에서 각도를 역추출하는 대신, 내부 _angles를 확인
        for i, (expected_deg, expected_rad) in enumerate(zip(expected_angles_deg, expected_radians)):
            assert frame._angles[i] == expected_deg, \
                f"채널 {i} 각도 불일치: {frame._angles[i]} != {expected_deg}"


class TestAnisotropicFrameWeights:
    """AnisotropicMultigroupFrame: 가중치 합이 1이어야 한다."""

    def test_anisotropic_frame_weights(self):
        frame = AnisotropicMultigroupFrame()
        w_sum = frame.weights.sum().item()
        assert abs(w_sum - 1.0) < 1e-10, f"가중치 합이 1이 아님: {w_sum}"


# ------------------------------------------------------------------
# ridge_solver 테스트
# ------------------------------------------------------------------

class TestFixedWeightSolver:
    """FixedWeightSolver: 논문 고정 가중치 반환."""

    def test_fixed_weight_solver(self):
        """가중치가 [0.40, 0.35, 0.25]이어야 한다."""
        frame = AnisotropicMultigroupFrame()
        solver = FixedWeightSolver(frame)

        expected = [0.40, 0.35, 0.25]
        actual = solver._w_k.tolist()

        for exp, act in zip(expected, actual):
            assert abs(exp - act) < 1e-10, f"가중치 불일치: {actual} != {expected}"


class TestRidgeWeightSolver:
    """RidgeWeightSolver: 출력 shape 검증."""

    def test_ridge_solver_output_shape(self):
        """솔버 출력 가중치의 길이가 채널 수와 동일해야 한다."""
        frame = MultigroupFrame(channel_type='iso3')
        solver = RidgeWeightSolver(frame)

        batch_size = 4
        sigma_prime = torch.randn(batch_size, dtype=PrecisionManager.REAL_DTYPE)
        tau_prime = torch.randn(batch_size, dtype=PrecisionManager.REAL_DTYPE)

        weights = solver.solve(sigma_prime, tau_prime)

        assert weights.shape == (batch_size, frame.K), \
            f"출력 shape 불일치: {weights.shape} != ({batch_size}, {frame.K})"
