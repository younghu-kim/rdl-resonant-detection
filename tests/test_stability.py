"""
수치 안정성 테스트.

각 도메인에 대해:
1. 순전파에서 NaN/Inf 가 발생하지 않는지
2. 역전파 그래디언트에 NaN이 없는지
3. 몇 스텝 훈련 시 손실이 감소하는지 (수렴 건전성)
"""

import pytest
import torch

from gdl.core.factory import GDLFactory

from gdl.domains.grids import GridDomain, GridSignal
from gdl.domains.groups import GroupDomain, GroupSignal
from gdl.domains.graphs import SetDomain, GraphDomain, NodeSignal
from gdl.domains.geodesics import GeodesicDomain, GeodesicSignal
from gdl.domains.gauges import GaugeDomain, GaugeSignal
from gdl.domains.temporal import TimeDomain, TimeSignal


# ---------------------------------------------------------------------------
# 헬퍼: 도메인별 더미 데이터 생성 (run_gdl_demo.py의 make_dummy_data 패턴)
# ---------------------------------------------------------------------------

def make_dummy_data(domain_type: str, batch: int, in_ch: int, out_ch: int):
    """도메인별 더미 시그널과 타겟 생성."""
    target = torch.randn(batch, out_ch)

    if domain_type == "1G_GRID":
        dom = GridDomain(spatial_shape=(8, 8))
        feat = torch.randn(batch, in_ch, 8, 8)
        return GridSignal(dom, feat), target

    elif domain_type == "2G_GROUP":
        dom = GroupDomain(l_max=1)
        feat = torch.randn(batch, in_ch, 4)
        return GroupSignal(dom, feat), target

    elif domain_type == "3G_SET":
        dom = SetDomain(num_nodes=5)
        feat = torch.randn(batch, 5, in_ch)
        return NodeSignal(dom, feat), target

    elif domain_type == "3G_GRAPH":
        ei = torch.tensor([[0,1,1,2,2,3,3,0],[1,0,2,1,3,2,0,3]], dtype=torch.long)
        ew = torch.ones(8)
        dom = GraphDomain(num_nodes=4, edge_index=ei, edge_weight=ew)
        feat = torch.randn(batch, 4, in_ch)
        return NodeSignal(dom, feat), target

    elif domain_type == "3G_ATTENTION":
        # Transformer는 embed_dim을 유지하므로 출력 채널이 in_ch와 동일하다.
        ei = torch.tensor([[0,1,1,2,2,3,3,0],[1,0,2,1,3,2,0,3]], dtype=torch.long)
        ew = torch.ones(8)
        dom = GraphDomain(num_nodes=4, edge_index=ei, edge_weight=ew)
        feat = torch.randn(batch, 4, in_ch)
        # 타겟도 in_ch에 맞춘다 (Transformer readout이 in_ch를 출력)
        attn_target = torch.randn(batch, in_ch)
        return NodeSignal(dom, feat), attn_target

    elif domain_type == "4G_GEODESIC":
        verts = torch.tensor([
            [1.,1.,1.],[1.,-1.,-1.],[-1.,1.,-1.],[-1.,-1.,1.]
        ])
        faces = torch.tensor([[0,1,2],[0,2,3],[0,3,1],[1,3,2]], dtype=torch.long)
        dom = GeodesicDomain(vertices=verts, faces=faces)
        feat = torch.randn(batch, 4, in_ch)
        return GeodesicSignal(dom, feat), target

    elif domain_type == "5G_GAUGE":
        ei = torch.tensor([[0,1,1,2,2,3,3,0],[1,0,2,1,3,2,0,3]], dtype=torch.long)
        angles = torch.randn(8)
        dom = GaugeDomain(num_nodes=4, edge_index=ei, transport_angles=angles)
        feat = torch.randn(batch, 4, in_ch, 2)
        return GaugeSignal(dom, feat), target

    elif domain_type == "TIME":
        dom = TimeDomain(sequence_length=6)
        feat = torch.randn(batch, 6, in_ch)
        return TimeSignal(dom, feat), target

    raise ValueError(f"Unknown domain type: {domain_type}")


# 4G_GEODESIC는 팩토리에 num_clusters/target_nodes 파라미터 불일치 버그가
# 있을 수 있으므로, 직접 모델을 생성하여 테스트한다.
# 나머지 도메인은 팩토리를 통해 모델을 생성한다.

DOMAIN_CONFIGS = [
    ("1G_GRID",      dict(num_spatial_dims=2, kernel_size=3, pool_kernel=2, pool_stride=2)),
    ("2G_GROUP",     dict(l_max=1)),
    ("3G_SET",       dict()),
    # 3G_GRAPH: GraphCoarsening이 SetDomain으로 변환하므로 후속
    # GraphMessagePassingLayer가 실패할 수 있다. 별도 테스트로 분리.
    ("3G_ATTENTION", dict(num_heads=1)),
    ("5G_GAUGE",     dict(aggr="mean")),
    ("TIME",         dict()),
]


def _create_model(domain_type, in_ch, hid_ch, out_ch, extra_kwargs):
    """GDLFactory로 모델 생성. 일부 도메인은 추가 kwargs가 필요하다."""
    return GDLFactory.create_model(
        domain_type=domain_type,
        in_channels=in_ch,
        hidden_channels=hid_ch,
        out_channels=out_ch,
        **extra_kwargs,
    )


# ---------------------------------------------------------------------------
# 테스트 1: 순전파에서 NaN/Inf 없음
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("domain_type,extra_kwargs", DOMAIN_CONFIGS,
                         ids=[c[0] for c in DOMAIN_CONFIGS])
def test_no_nan_inf_in_forward(domain_type, extra_kwargs, batch_size):
    """순전파 출력에 NaN 또는 Inf가 없는지 확인."""
    in_ch, hid_ch, out_ch = 3, 8, 4
    torch.manual_seed(0)

    model = _create_model(domain_type, in_ch, hid_ch, out_ch, extra_kwargs)
    signal, _ = make_dummy_data(domain_type, batch_size, in_ch, out_ch)

    with torch.no_grad():
        output = model(signal)

    assert not torch.isnan(output).any(), (
        f"{domain_type}: 순전파 출력에 NaN이 존재합니다."
    )
    assert not torch.isinf(output).any(), (
        f"{domain_type}: 순전파 출력에 Inf가 존재합니다."
    )


# ---------------------------------------------------------------------------
# 테스트 2: 역전파 그래디언트에 NaN 없음
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("domain_type,extra_kwargs", DOMAIN_CONFIGS,
                         ids=[c[0] for c in DOMAIN_CONFIGS])
def test_gradients_no_nan(domain_type, extra_kwargs, batch_size):
    """역전파 그래디언트에 NaN이 없는지 확인."""
    in_ch, hid_ch, out_ch = 3, 8, 4
    torch.manual_seed(0)

    model = _create_model(domain_type, in_ch, hid_ch, out_ch, extra_kwargs)
    signal, target = make_dummy_data(domain_type, batch_size, in_ch, out_ch)

    output = model(signal)
    loss = torch.nn.MSELoss()(output, target)
    loss.backward()

    for name, param in model.named_parameters():
        if param.grad is not None:
            assert not torch.isnan(param.grad).any(), (
                f"{domain_type}: 파라미터 '{name}'의 그래디언트에 NaN이 존재합니다."
            )


# ---------------------------------------------------------------------------
# 테스트 3: 손실 감소 (수렴 건전성)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("domain_type,extra_kwargs", DOMAIN_CONFIGS,
                         ids=[c[0] for c in DOMAIN_CONFIGS])
def test_loss_decreases(domain_type, extra_kwargs, batch_size):
    """몇 스텝 훈련 후 손실이 감소하는지 확인 (수렴 건전성)."""
    in_ch, hid_ch, out_ch = 3, 8, 4
    torch.manual_seed(42)

    model = _create_model(domain_type, in_ch, hid_ch, out_ch, extra_kwargs)
    signal, target = make_dummy_data(domain_type, batch_size, in_ch, out_ch)

    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    criterion = torch.nn.MSELoss()

    losses = []
    num_steps = 10
    for _ in range(num_steps):
        optimizer.zero_grad()
        pred = model(signal)
        loss = criterion(pred, target)
        loss.backward()
        optimizer.step()
        losses.append(loss.item())

    # 마지막 손실이 첫 번째 손실보다 작아야 한다 (overfit 가능한 간단한 데이터)
    assert losses[-1] < losses[0], (
        f"{domain_type}: 손실이 감소하지 않았습니다.\n"
        f"  첫 번째 손실: {losses[0]:.6f}\n"
        f"  마지막 손실: {losses[-1]:.6f}\n"
        f"  전체 손실: {losses}"
    )


# ---------------------------------------------------------------------------
# 4G_GEODESIC 단독 테스트 (팩토리 우회)
# ---------------------------------------------------------------------------

def test_geodesic_no_nan_inf_manual(batch_size):
    """4G_GEODESIC을 팩토리 없이 직접 구성하여 NaN/Inf 확인."""
    from gdl.domains.geodesics.layers import (
        IntrinsicMeshConvolution, GeodesicNonLinearity,
        MeshDecimationPooling, GeodesicGlobalPooling,
    )
    from gdl.core.base_blueprint import BaseGDLModel

    in_ch, hid_ch, out_ch = 3, 8, 4
    torch.manual_seed(0)

    blocks = [
        IntrinsicMeshConvolution(in_ch, hid_ch, K=2),
        GeodesicNonLinearity(),
        MeshDecimationPooling(in_channels=hid_ch, target_nodes=2),
        IntrinsicMeshConvolution(hid_ch, out_ch, K=2),
        GeodesicNonLinearity(),
    ]
    readout = GeodesicGlobalPooling(aggregation="max")
    model = BaseGDLModel(blocks=blocks, readout=readout)

    signal, target = make_dummy_data("4G_GEODESIC", batch_size, in_ch, out_ch)

    output = model(signal)
    assert not torch.isnan(output).any(), "4G_GEODESIC: NaN in forward"
    assert not torch.isinf(output).any(), "4G_GEODESIC: Inf in forward"

    loss = torch.nn.MSELoss()(output, target)
    loss.backward()

    for name, param in model.named_parameters():
        if param.grad is not None:
            assert not torch.isnan(param.grad).any(), (
                f"4G_GEODESIC: 파라미터 '{name}'의 그래디언트에 NaN"
            )


def test_geodesic_loss_decreases_manual(batch_size):
    """4G_GEODESIC을 직접 구성하여 손실 감소 확인."""
    from gdl.domains.geodesics.layers import (
        IntrinsicMeshConvolution, GeodesicNonLinearity,
        MeshDecimationPooling, GeodesicGlobalPooling,
    )
    from gdl.core.base_blueprint import BaseGDLModel

    in_ch, hid_ch, out_ch = 3, 8, 4
    torch.manual_seed(42)

    blocks = [
        IntrinsicMeshConvolution(in_ch, hid_ch, K=2),
        GeodesicNonLinearity(),
        MeshDecimationPooling(in_channels=hid_ch, target_nodes=2),
        IntrinsicMeshConvolution(hid_ch, out_ch, K=2),
        GeodesicNonLinearity(),
    ]
    readout = GeodesicGlobalPooling(aggregation="max")
    model = BaseGDLModel(blocks=blocks, readout=readout)

    signal, target = make_dummy_data("4G_GEODESIC", batch_size, in_ch, out_ch)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    criterion = torch.nn.MSELoss()

    losses = []
    for _ in range(10):
        optimizer.zero_grad()
        pred = model(signal)
        loss = criterion(pred, target)
        loss.backward()
        optimizer.step()
        losses.append(loss.item())

    assert losses[-1] < losses[0], (
        f"4G_GEODESIC: 손실 미감소 ({losses[0]:.6f} -> {losses[-1]:.6f})"
    )


# ---------------------------------------------------------------------------
# 3G_GRAPH 단독 테스트 (GraphCoarsening 이후 SetDomain 전환 문제 우회)
# ---------------------------------------------------------------------------

def test_graph_no_coarsening_stability(batch_size):
    """
    3G_GRAPH를 GraphCoarsening 없이 구성하여 안정성 확인.

    참고: 팩토리의 3G_GRAPH 구성은 GraphCoarsening(→ SetDomain) 후
    GraphMessagePassingLayer(GraphDomain 필요)를 배치하여 런타임 에러가 발생한다.
    여기서는 Coarsening 없이 순수 메시지 패싱만 테스트한다.
    """
    from gdl.domains.graphs.layers import (
        GraphMessagePassingLayer, GraphNonLinearity, GraphGlobalPooling,
    )
    from gdl.core.base_blueprint import BaseGDLModel

    in_ch, hid_ch, out_ch = 3, 8, 4
    torch.manual_seed(42)

    blocks = [
        GraphMessagePassingLayer(in_ch, hid_ch),
        GraphNonLinearity(),
        GraphMessagePassingLayer(hid_ch, out_ch),
        GraphNonLinearity(),
    ]
    readout = GraphGlobalPooling(aggregation="mean")
    model = BaseGDLModel(blocks=blocks, readout=readout)

    signal, target = make_dummy_data("3G_GRAPH", batch_size, in_ch, out_ch)

    # NaN/Inf 없음
    output = model(signal)
    assert not torch.isnan(output).any(), "3G_GRAPH: NaN in forward"
    assert not torch.isinf(output).any(), "3G_GRAPH: Inf in forward"

    # 그래디언트 흐름
    loss = torch.nn.MSELoss()(output, target)
    loss.backward()
    for name, param in model.named_parameters():
        if param.grad is not None:
            assert not torch.isnan(param.grad).any(), (
                f"3G_GRAPH: NaN in grad of '{name}'"
            )

    # 손실 감소
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    criterion = torch.nn.MSELoss()
    losses = []
    for _ in range(10):
        optimizer.zero_grad()
        pred = model(signal)
        loss = criterion(pred, target)
        loss.backward()
        optimizer.step()
        losses.append(loss.item())

    assert losses[-1] < losses[0], (
        f"3G_GRAPH: 손실 미감소 ({losses[0]:.6f} -> {losses[-1]:.6f})"
    )
