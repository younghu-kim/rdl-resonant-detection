"""
등변 레이어(Equivariant Layer) 출력 형상 및 도메인 보존 테스트.

각 도메인의 등변 레이어가 올바른 시그널 타입과 형상을 유지하는지 검증한다.
"""

import pytest
import torch

# GDLFactory를 먼저 임포트하여 순환 임포트를 방지한다.
from gdl.core.factory import GDLFactory

from gdl.domains.grids import GridDomain, GridSignal
from gdl.domains.grids.layers import GridConvolution, GridNonLinearity

from gdl.domains.groups import GroupDomain, GroupSignal
from gdl.domains.groups.layers import GroupConvolution, GroupNonLinearity

from gdl.domains.graphs import SetDomain, GraphDomain, NodeSignal
from gdl.domains.graphs.layers import (
    DeepSetsLayer,
    GraphMessagePassingLayer,
    TransformerEquivariantLayer,
    GraphNonLinearity,
)

from gdl.domains.geodesics import GeodesicDomain, GeodesicSignal
from gdl.domains.geodesics.layers import IntrinsicMeshConvolution, GeodesicNonLinearity

from gdl.domains.gauges import GaugeDomain, GaugeSignal
from gdl.domains.gauges.layers import GaugeMessagePassingLayer, GaugeNonLinearity

from gdl.domains.temporal import TimeDomain, TimeSignal
from gdl.domains.temporal.layers import TimeDifferenceGatedRNN, TimeNonLinearity


# ---------------------------------------------------------------------------
# 헬퍼: 도메인별 시그널 + 레이어 생성
# ---------------------------------------------------------------------------

def _make_grid(batch, in_ch, out_ch):
    dom = GridDomain(spatial_shape=(8, 8))
    sig = GridSignal(dom, torch.randn(batch, in_ch, 8, 8))
    layer = GridConvolution(in_ch, out_ch, kernel_size=3, num_spatial_dims=2)
    return sig, layer, GridSignal, GridDomain

def _make_group(batch, in_ch, out_ch):
    dom = GroupDomain(l_max=1)
    # (l_max+1)^2 = 4 spectral components
    sig = GroupSignal(dom, torch.randn(batch, in_ch, 4))
    layer = GroupConvolution(in_ch, out_ch, l_max=1)
    return sig, layer, GroupSignal, GroupDomain

def _make_set(batch, in_ch, out_ch):
    dom = SetDomain(num_nodes=5)
    sig = NodeSignal(dom, torch.randn(batch, 5, in_ch))
    layer = DeepSetsLayer(in_ch, out_ch)
    return sig, layer, NodeSignal, SetDomain

def _make_graph(batch, in_ch, out_ch):
    ei = torch.tensor([[0,1,1,2,2,3,3,0],[1,0,2,1,3,2,0,3]], dtype=torch.long)
    ew = torch.ones(8)
    dom = GraphDomain(num_nodes=4, edge_index=ei, edge_weight=ew)
    sig = NodeSignal(dom, torch.randn(batch, 4, in_ch))
    layer = GraphMessagePassingLayer(in_ch, out_ch)
    return sig, layer, NodeSignal, GraphDomain

def _make_attention(batch, in_ch, out_ch):
    # Transformer requires embed_dim == in_channels and does not change channel count
    dom = SetDomain(num_nodes=4)
    sig = NodeSignal(dom, torch.randn(batch, 4, in_ch))
    layer = TransformerEquivariantLayer(embed_dim=in_ch, num_heads=1)
    # output channels == input channels for attention
    return sig, layer, NodeSignal, SetDomain

def _make_geodesic(batch, in_ch, out_ch):
    verts = torch.tensor([[1.,1.,1.],[1.,-1.,-1.],[-1.,1.,-1.],[-1.,-1.,1.]])
    faces = torch.tensor([[0,1,2],[0,2,3],[0,3,1],[1,3,2]], dtype=torch.long)
    dom = GeodesicDomain(vertices=verts, faces=faces)
    sig = GeodesicSignal(dom, torch.randn(batch, 4, in_ch))
    layer = IntrinsicMeshConvolution(in_ch, out_ch, K=2)
    return sig, layer, GeodesicSignal, GeodesicDomain

def _make_gauge(batch, in_ch, out_ch):
    ei = torch.tensor([[0,1,1,2,2,3,3,0],[1,0,2,1,3,2,0,3]], dtype=torch.long)
    angles = torch.randn(8)
    dom = GaugeDomain(num_nodes=4, edge_index=ei, transport_angles=angles)
    sig = GaugeSignal(dom, torch.randn(batch, 4, in_ch, 2))
    layer = GaugeMessagePassingLayer(in_ch, out_ch, aggr="mean")
    return sig, layer, GaugeSignal, GaugeDomain

def _make_temporal(batch, in_ch, out_ch):
    dom = TimeDomain(sequence_length=6)
    sig = TimeSignal(dom, torch.randn(batch, 6, in_ch))
    layer = TimeDifferenceGatedRNN(in_ch, out_ch)
    return sig, layer, TimeSignal, TimeDomain


DOMAIN_MAKERS = {
    "grid": _make_grid,
    "group": _make_group,
    "set": _make_set,
    "graph": _make_graph,
    "attention": _make_attention,
    "geodesic": _make_geodesic,
    "gauge": _make_gauge,
    "temporal": _make_temporal,
}


# ---------------------------------------------------------------------------
# 테스트
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("domain_name", list(DOMAIN_MAKERS.keys()))
def test_equivariant_layer_output_type(domain_name, batch_size):
    """등변 레이어의 출력이 동일한 시그널 타입인지 확인."""
    maker = DOMAIN_MAKERS[domain_name]
    in_ch, out_ch = 3, 8
    signal, layer, expected_signal_cls, _ = maker(batch_size, in_ch, out_ch)

    with torch.no_grad():
        output = layer(signal)

    assert isinstance(output, expected_signal_cls), (
        f"{domain_name}: 출력 타입이 {expected_signal_cls.__name__}이어야 하지만 "
        f"{type(output).__name__}을 받았습니다."
    )


@pytest.mark.parametrize("domain_name", list(DOMAIN_MAKERS.keys()))
def test_equivariant_layer_output_shape(domain_name, batch_size):
    """등변 레이어의 출력 형상이 올바른지 확인."""
    maker = DOMAIN_MAKERS[domain_name]
    in_ch, out_ch = 3, 8
    signal, layer, _, _ = maker(batch_size, in_ch, out_ch)

    with torch.no_grad():
        output = layer(signal)

    out_feat = output.features
    # 배치 차원은 항상 보존
    assert out_feat.shape[0] == batch_size, (
        f"{domain_name}: 배치 차원 불일치 (expected {batch_size}, got {out_feat.shape[0]})"
    )
    # 텐서가 비어있지 않은지
    assert out_feat.numel() > 0, f"{domain_name}: 출력 텐서가 비어있습니다."


@pytest.mark.parametrize("domain_name", list(DOMAIN_MAKERS.keys()))
def test_equivariant_layer_preserves_domain(domain_name, batch_size):
    """등변 레이어가 도메인 정보를 보존하는지 확인."""
    maker = DOMAIN_MAKERS[domain_name]
    in_ch, out_ch = 3, 8
    signal, layer, _, expected_domain_cls = maker(batch_size, in_ch, out_ch)

    with torch.no_grad():
        output = layer(signal)

    assert isinstance(output.domain, expected_domain_cls), (
        f"{domain_name}: 도메인 타입이 {expected_domain_cls.__name__}이어야 하지만 "
        f"{type(output.domain).__name__}을 받았습니다."
    )


# ---------------------------------------------------------------------------
# 비선형 레이어 테스트 (형상 보존)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("domain_name,nonlin_cls", [
    ("grid", GridNonLinearity),
    ("group", GroupNonLinearity),
    ("set", GraphNonLinearity),
    ("graph", GraphNonLinearity),
    ("geodesic", GeodesicNonLinearity),
    ("gauge", GaugeNonLinearity),
    ("temporal", TimeNonLinearity),
])
def test_nonlinearity_preserves_shape(domain_name, nonlin_cls, batch_size):
    """비선형 레이어가 시그널 형상을 변경하지 않는지 확인."""
    maker = DOMAIN_MAKERS[domain_name]
    in_ch = 4
    signal, _, _, _ = maker(batch_size, in_ch, in_ch)
    nonlin = nonlin_cls()

    with torch.no_grad():
        output = nonlin(signal)

    assert output.features.shape == signal.features.shape, (
        f"{domain_name}: 비선형 레이어가 형상을 변경했습니다 "
        f"({signal.features.shape} -> {output.features.shape})"
    )
