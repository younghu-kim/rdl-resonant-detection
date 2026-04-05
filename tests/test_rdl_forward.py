"""
=============================================================================
[Test] RDL Forward Pass End-to-End
=============================================================================
MasterResonantNetwork의 전방향 전파 및 역전파 테스트.
"""

import pytest
import torch

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.models.master_net import MasterResonantNetwork


# 테스트 공통 설정
IN_DIM = 4
HIDDEN_DIM = 8
OUT_DIM = 2
NUM_LAYERS = 2
BATCH_SIZE = 2

EXPECTED_KEYS = {
    "X_in", "Z_out", "Psi_target", "L_G",
    "phi", "psi", "tau_prime", "phase_out",
}


@pytest.fixture
def model():
    """MasterResonantNetwork 인스턴스."""
    return MasterResonantNetwork(
        in_features=IN_DIM,
        hidden_features=HIDDEN_DIM,
        out_features=OUT_DIM,
        num_layers=NUM_LAYERS,
    )


@pytest.fixture
def sample_input():
    """테스트용 입력 텐서 [Batch, in_features]."""
    return torch.randn(
        (BATCH_SIZE, IN_DIM),
        dtype=PrecisionManager.REAL_DTYPE,
        requires_grad=True,
    )


class TestMasterNetForward:
    """MasterResonantNetwork forward 패스 출력 검증."""

    def test_master_net_forward(self, model, sample_input):
        """출력 딕셔너리에 필요한 모든 키가 존재해야 한다."""
        outputs = model(sample_input)

        assert isinstance(outputs, dict), "출력이 딕셔너리가 아님"
        for key in EXPECTED_KEYS:
            assert key in outputs, f"출력 딕셔너리에 '{key}' 키가 없음"

    def test_master_net_output_shapes(self, model, sample_input):
        """Z_out은 복소수, phi/psi는 실수이며 shape이 배치와 일치해야 한다."""
        outputs = model(sample_input)

        Z_out = outputs["Z_out"]
        phi = outputs["phi"]
        psi = outputs["psi"]

        # Z_out은 복소수 타입
        assert Z_out.is_complex(), f"Z_out이 복소수가 아님: {Z_out.dtype}"

        # Z_out shape: [Batch, out_features]
        assert Z_out.shape[0] == BATCH_SIZE, \
            f"Z_out 배치 크기 불일치: {Z_out.shape[0]} != {BATCH_SIZE}"

        # phi, psi는 실수
        assert not phi.is_complex(), f"phi가 복소수임: {phi.dtype}"
        assert not psi.is_complex(), f"psi가 복소수임: {psi.dtype}"

        # phi, psi 배치 차원 일치
        assert phi.shape[0] == BATCH_SIZE, \
            f"phi 배치 크기 불일치: {phi.shape[0]} != {BATCH_SIZE}"
        assert psi.shape[0] == BATCH_SIZE, \
            f"psi 배치 크기 불일치: {psi.shape[0]} != {BATCH_SIZE}"

    def test_master_net_gradient_flow(self, model, sample_input):
        """역전파 시 입력 X까지 그래디언트가 도달해야 한다."""
        outputs = model(sample_input)

        # Z_out의 절댓값 합으로 스칼라 손실 생성
        loss = outputs["Z_out"].abs().sum()
        loss.backward()

        assert sample_input.grad is not None, "입력 텐서에 그래디언트가 없음"
        assert not torch.isnan(sample_input.grad).any(), \
            "입력 그래디언트에 NaN 존재"

    def test_master_net_no_nan(self, model, sample_input):
        """모든 출력 텐서에 NaN이 없어야 한다."""
        outputs = model(sample_input)

        for key, tensor in outputs.items():
            if isinstance(tensor, torch.Tensor):
                if tensor.is_complex():
                    assert not torch.isnan(tensor.real).any(), \
                        f"{key} 실수부에 NaN 존재"
                    assert not torch.isnan(tensor.imag).any(), \
                        f"{key} 허수부에 NaN 존재"
                else:
                    assert not torch.isnan(tensor).any(), \
                        f"{key}에 NaN 존재"
