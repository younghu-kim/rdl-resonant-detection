"""
=============================================================================
[Test] RDL Loss Functions
=============================================================================
DiscriminantLoss, CurvatureStabilizationLoss, TargetPhaseMatchingLoss,
PhaseQuantisationLoss, TotalResonanceLoss 모듈 테스트.
"""

import pytest
import torch

from gdl.rdl.constants import PrecisionManager
from gdl.rdl.losses.phase_loss import DiscriminantLoss
from gdl.rdl.losses.curvature_loss import CurvatureStabilizationLoss
from gdl.rdl.losses.target_loss import TargetPhaseMatchingLoss
from gdl.rdl.losses.quantisation_loss import PhaseQuantisationLoss
from gdl.rdl.losses.total_loss import TotalResonanceLoss
from gdl.rdl.models.master_net import MasterResonantNetwork

DTYPE_R = PrecisionManager.REAL_DTYPE
DTYPE_C = PrecisionManager.COMPLEX_DTYPE
BATCH = 2
FEATURES = 4


class TestDiscriminantLoss:
    """DiscriminantLoss: 판별식 잔차 손실."""

    def test_discriminant_loss_scalar(self):
        """출력이 스칼라이며 0 이상이어야 한다."""
        loss_fn = DiscriminantLoss(reduction='mean')

        L_G = torch.complex(
            torch.randn(BATCH, FEATURES, dtype=DTYPE_R),
            torch.randn(BATCH, FEATURES, dtype=DTYPE_R),
        )
        phi = torch.randn(BATCH, FEATURES, dtype=DTYPE_R)
        psi = torch.randn(BATCH, FEATURES, dtype=DTYPE_R)

        loss = loss_fn(L_G, phi, psi)

        assert loss.dim() == 0, f"스칼라가 아님: dim={loss.dim()}"
        assert loss.item() >= 0.0, f"손실이 음수: {loss.item()}"
        assert not torch.isnan(loss), f"NaN 발생"


class TestTotalResonanceLoss:
    """TotalResonanceLoss: 4항 통합 손실."""

    @pytest.fixture
    def model_and_outputs(self):
        """MasterResonantNetwork으로 forward pass 실행 후 outputs 반환."""
        in_dim, hidden_dim, out_dim = 4, 8, 2
        model = MasterResonantNetwork(
            in_features=in_dim,
            hidden_features=hidden_dim,
            out_features=out_dim,
            num_layers=2,
        )
        X = torch.randn(BATCH, in_dim, dtype=DTYPE_R, requires_grad=True)
        outputs = model(X)
        return model, X, outputs

    def test_total_loss_returns_tuple(self, model_and_outputs):
        """반환값이 (loss_tensor, metrics_dict) 튜플이어야 한다."""
        _, _, outputs = model_and_outputs
        loss_fn = TotalResonanceLoss(
            lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.5,
        )

        result = loss_fn(**outputs)

        assert isinstance(result, tuple), f"튜플이 아님: {type(result)}"
        assert len(result) == 2, f"튜플 길이가 2가 아님: {len(result)}"

        total_loss, metrics = result
        assert isinstance(total_loss, torch.Tensor), "첫 번째 요소가 텐서가 아님"
        assert total_loss.dim() == 0, f"손실이 스칼라가 아님: dim={total_loss.dim()}"
        assert isinstance(metrics, dict), "두 번째 요소가 딕셔너리가 아님"

    def test_total_loss_backward(self, model_and_outputs):
        """역전파 시 그래디언트가 흐르어야 한다."""
        model, X, outputs = model_and_outputs
        loss_fn = TotalResonanceLoss(
            lambda_res=1.0, lambda_curv=0.1, lambda_tgt=1.0, lambda_pqo=0.5,
        )

        total_loss, metrics = loss_fn(**outputs)
        total_loss.backward()

        # 입력 X까지 그래디언트 도달 확인
        assert X.grad is not None, "입력 X에 그래디언트가 없음"
        assert not torch.isnan(X.grad).any(), "입력 그래디언트에 NaN 존재"

        # 모델 파라미터에도 그래디언트 존재 확인
        has_grad = False
        for p in model.parameters():
            if p.grad is not None:
                has_grad = True
                break
        assert has_grad, "모델 파라미터에 그래디언트가 전혀 없음"
