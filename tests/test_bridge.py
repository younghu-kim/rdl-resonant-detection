"""
=============================================================================
[Test] RDL <-> GDL 5G Bridge
=============================================================================
rdl_gauge_domain, gauge_rdl_adapter 모듈 테스트.
RDL의 U(1) 게이지 대칭과 GDL의 SO(2) 구조군 간 변환을 검증한다.
"""

import pytest
import torch

from gdl.rdl.constants import PrecisionManager
from gdl.bridge.rdl_gauge_domain import (
    time_grid_to_gauge_domain,
    complex_to_gauge_signal,
    gauge_signal_to_complex,
)
from gdl.bridge.gauge_rdl_adapter import RDLGaugeAdapter
from gdl.domains.gauges.domain import GaugeDomain
from gdl.domains.gauges.signal import GaugeSignal

NUM_NODES = 10
IN_CH = 3
OUT_CH = 2
BATCH = 2


@pytest.fixture
def gauge_domain():
    """1D 시간 격자로 생성한 GaugeDomain."""
    return time_grid_to_gauge_domain(num_points=NUM_NODES, t_min=10.0, t_max=50.0)


@pytest.fixture
def sample_complex(gauge_domain):
    """테스트용 복소 텐서 [B, N, C]."""
    return torch.complex(
        torch.randn(BATCH, NUM_NODES, IN_CH, dtype=torch.float64),
        torch.randn(BATCH, NUM_NODES, IN_CH, dtype=torch.float64),
    )


@pytest.fixture
def sample_gauge_signal(gauge_domain, sample_complex):
    """테스트용 GaugeSignal [B, N, C, 2]."""
    return complex_to_gauge_signal(sample_complex, gauge_domain)


class TestTimeGridDomain:
    """time_grid_to_gauge_domain: 1D 격자 -> GaugeDomain 변환."""

    def test_time_grid_domain(self, gauge_domain):
        """올바른 노드 수와 양방향 엣지를 생성해야 한다."""
        # 노드 수 검증
        assert gauge_domain.num_nodes == NUM_NODES, \
            f"노드 수 불일치: {gauge_domain.num_nodes} != {NUM_NODES}"

        # 엣지 수: (N-1) 양방향 = 2*(N-1)
        expected_edges = 2 * (NUM_NODES - 1)
        actual_edges = gauge_domain.edge_index.size(1)
        assert actual_edges == expected_edges, \
            f"엣지 수 불일치: {actual_edges} != {expected_edges}"

        # 양방향 검증: 각 (i, i+1) 쌍에 대해 역방향 (i+1, i)도 존재
        src = gauge_domain.edge_index[0]
        dst = gauge_domain.edge_index[1]
        edge_set = set(zip(src.tolist(), dst.tolist()))
        for i in range(NUM_NODES - 1):
            assert (i, i + 1) in edge_set, f"정방향 엣지 ({i}, {i+1}) 없음"
            assert (i + 1, i) in edge_set, f"역방향 엣지 ({i+1}, {i}) 없음"


class TestComplexGaugeRoundtrip:
    """complex -> gauge -> complex 왕복 변환 정합성."""

    def test_complex_gauge_roundtrip(self, gauge_domain, sample_complex):
        """왕복 변환 시 값이 보존되어야 한다 (float32 중간 변환 허용오차)."""
        # complex -> gauge signal
        signal = complex_to_gauge_signal(sample_complex, gauge_domain)

        # gauge signal -> complex
        recovered = gauge_signal_to_complex(signal)

        # float32 중간 변환으로 인한 정밀도 손실 허용
        # recovered is complex128, sample_complex is complex128
        assert torch.allclose(
            sample_complex,
            recovered,
            atol=1e-6,
            rtol=1e-5,
        ), "왕복 변환 후 값이 보존되지 않음"


class TestRDLGaugeAdapter:
    """RDLGaugeAdapter: RDL -> GDL 5G 등변 레이어."""

    @pytest.fixture
    def adapter(self):
        """RDLGaugeAdapter 인스턴스."""
        return RDLGaugeAdapter(
            num_nodes=NUM_NODES,
            in_channels=IN_CH,
            out_channels=OUT_CH,
            hidden_features=16,
            num_layers=2,
        )

    def test_rdl_adapter_forward(self, adapter, sample_gauge_signal):
        """어댑터가 올바른 shape의 GaugeSignal을 출력해야 한다."""
        output = adapter(sample_gauge_signal)

        assert isinstance(output, GaugeSignal), f"출력 타입 불일치: {type(output)}"

        # 출력 shape: [B, N, out_ch, 2]
        expected_shape = (BATCH, NUM_NODES, OUT_CH, 2)
        assert output.features.shape == expected_shape, \
            f"출력 shape 불일치: {output.features.shape} != {expected_shape}"

    def test_rdl_adapter_gradient(self, adapter, gauge_domain):
        """어댑터를 통해 그래디언트가 흐르어야 한다."""
        # requires_grad가 가능한 입력 생성
        features = torch.randn(
            BATCH, NUM_NODES, IN_CH, 2,
            dtype=torch.float32,
            requires_grad=True,
        )
        signal = GaugeSignal(domain=gauge_domain, features=features)

        output = adapter(signal)
        loss = output.features.sum()
        loss.backward()

        assert features.grad is not None, "입력 features에 그래디언트가 없음"
        assert not torch.isnan(features.grad).any(), "그래디언트에 NaN 존재"

    def test_domain_preservation(self, adapter, sample_gauge_signal):
        """출력의 도메인이 입력의 도메인과 동일해야 한다."""
        output = adapter(sample_gauge_signal)

        in_domain = sample_gauge_signal.domain
        out_domain = output.domain

        assert out_domain.num_nodes == in_domain.num_nodes, \
            "도메인 노드 수 불일치"
        assert torch.equal(out_domain.edge_index, in_domain.edge_index), \
            "도메인 엣지 불일치"
        assert torch.equal(out_domain.transport_angles, in_domain.transport_angles), \
            "도메인 평행 이동 각도 불일치"
