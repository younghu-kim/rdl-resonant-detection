"""
전역 불변(Global Invariant) 레이어 테스트.

각 도메인의 GlobalInvariant(리드아웃)이 도메인 독립적인 텐서를 생성하는지,
그리고 출력 형상이 [batch, channels]인지 검증한다.
그래프 도메인에서는 노드 순열 불변성도 테스트한다.
"""

import pytest
import torch

# GDLFactory를 먼저 임포트하여 순환 임포트를 방지한다.
from gdl.core.factory import GDLFactory
from gdl.core.base_signal import BaseSignal

from gdl.domains.grids import GridDomain, GridSignal
from gdl.domains.grids.layers import GridGlobalAveragePooling

from gdl.domains.groups import GroupDomain, GroupSignal
from gdl.domains.groups.layers import GroupGlobalPooling

from gdl.domains.graphs import SetDomain, GraphDomain, NodeSignal
from gdl.domains.graphs.layers import GraphGlobalPooling

from gdl.domains.geodesics import GeodesicDomain, GeodesicSignal
from gdl.domains.geodesics.layers import GeodesicGlobalPooling

from gdl.domains.gauges import GaugeDomain, GaugeSignal
from gdl.domains.gauges.layers import GaugeGlobalPooling

from gdl.domains.temporal import TimeDomain, TimeSignal
from gdl.domains.temporal.layers import TimeGlobalPooling


# ---------------------------------------------------------------------------
# 헬퍼: 도메인별 시그널 + GlobalInvariant 생성
# ---------------------------------------------------------------------------

def _global_grid(batch, ch):
    dom = GridDomain(spatial_shape=(8, 8))
    sig = GridSignal(dom, torch.randn(batch, ch, 8, 8))
    pool = GridGlobalAveragePooling()
    expected_out_ch = ch
    return sig, pool, expected_out_ch

def _global_group(batch, ch):
    dom = GroupDomain(l_max=1)
    sig = GroupSignal(dom, torch.randn(batch, ch, 4))
    pool = GroupGlobalPooling(use_power_spectrum=True)
    # power spectrum: channels * (l_max + 1) = ch * 2
    expected_out_ch = ch * 2
    return sig, pool, expected_out_ch

def _global_set(batch, ch):
    dom = SetDomain(num_nodes=5)
    sig = NodeSignal(dom, torch.randn(batch, 5, ch))
    pool = GraphGlobalPooling(aggregation="sum")
    expected_out_ch = ch
    return sig, pool, expected_out_ch

def _global_graph(batch, ch):
    ei = torch.tensor([[0,1,1,2,2,3,3,0],[1,0,2,1,3,2,0,3]], dtype=torch.long)
    ew = torch.ones(8)
    dom = GraphDomain(num_nodes=4, edge_index=ei, edge_weight=ew)
    sig = NodeSignal(dom, torch.randn(batch, 4, ch))
    pool = GraphGlobalPooling(aggregation="mean")
    expected_out_ch = ch
    return sig, pool, expected_out_ch

def _global_geodesic(batch, ch):
    verts = torch.tensor([[1.,1.,1.],[1.,-1.,-1.],[-1.,1.,-1.],[-1.,-1.,1.]])
    faces = torch.tensor([[0,1,2],[0,2,3],[0,3,1],[1,3,2]], dtype=torch.long)
    dom = GeodesicDomain(vertices=verts, faces=faces)
    sig = GeodesicSignal(dom, torch.randn(batch, 4, ch))
    pool = GeodesicGlobalPooling(aggregation="max")
    expected_out_ch = ch
    return sig, pool, expected_out_ch

def _global_gauge(batch, ch):
    ei = torch.tensor([[0,1,1,2,2,3,3,0],[1,0,2,1,3,2,0,3]], dtype=torch.long)
    angles = torch.randn(8)
    dom = GaugeDomain(num_nodes=4, edge_index=ei, transport_angles=angles)
    sig = GaugeSignal(dom, torch.randn(batch, 4, ch, 2))
    pool = GaugeGlobalPooling(aggregation="mean")
    expected_out_ch = ch
    return sig, pool, expected_out_ch

def _global_temporal(batch, ch):
    dom = TimeDomain(sequence_length=6)
    sig = TimeSignal(dom, torch.randn(batch, 6, ch))
    pool = TimeGlobalPooling(aggregation="last")
    expected_out_ch = ch
    return sig, pool, expected_out_ch


GLOBAL_MAKERS = {
    "grid": _global_grid,
    "group": _global_group,
    "set": _global_set,
    "graph": _global_graph,
    "geodesic": _global_geodesic,
    "gauge": _global_gauge,
    "temporal": _global_temporal,
}


# ---------------------------------------------------------------------------
# 테스트: 출력 타입 (raw Tensor, not Signal)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("domain_name", list(GLOBAL_MAKERS.keys()))
def test_global_invariant_returns_raw_tensor(domain_name, batch_size):
    """GlobalInvariant의 출력이 BaseSignal이 아닌 순수 torch.Tensor인지 확인."""
    maker = GLOBAL_MAKERS[domain_name]
    ch = 4
    signal, pool, _ = maker(batch_size, ch)

    with torch.no_grad():
        output = pool(signal)

    assert isinstance(output, torch.Tensor), (
        f"{domain_name}: 출력이 torch.Tensor여야 합니다 (got {type(output).__name__})"
    )
    assert not isinstance(output, BaseSignal), (
        f"{domain_name}: 출력이 BaseSignal이면 안 됩니다 — 도메인 정보가 제거되어야 합니다."
    )


# ---------------------------------------------------------------------------
# 테스트: 출력 형상 [batch, channels]
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("domain_name", list(GLOBAL_MAKERS.keys()))
def test_global_invariant_output_shape(domain_name, batch_size):
    """GlobalInvariant 출력 형상이 [batch, channels]인지 확인."""
    maker = GLOBAL_MAKERS[domain_name]
    ch = 4
    signal, pool, expected_out_ch = maker(batch_size, ch)

    with torch.no_grad():
        output = pool(signal)

    assert output.ndim == 2, (
        f"{domain_name}: 출력 텐서가 2D여야 합니다 (got {output.ndim}D)"
    )
    assert output.shape[0] == batch_size, (
        f"{domain_name}: 배치 차원 불일치 (expected {batch_size}, got {output.shape[0]})"
    )
    assert output.shape[1] == expected_out_ch, (
        f"{domain_name}: 채널 차원 불일치 (expected {expected_out_ch}, got {output.shape[1]})"
    )


# ---------------------------------------------------------------------------
# 테스트: 그래프 순열 불변성
# ---------------------------------------------------------------------------

def test_graph_global_pooling_permutation_invariance(batch_size):
    """
    그래프 GlobalPooling이 노드 순열에 불변인지 검증.

    동일한 그래프에서 노드 순서만 바꾼 뒤 GlobalPooling 결과가 동일해야 한다.
    sum/mean 집계는 순서와 무관하므로 정확히 일치해야 한다.
    """
    num_nodes = 5
    ch = 4
    torch.manual_seed(42)

    # 원본 시그널
    dom = SetDomain(num_nodes=num_nodes)
    features = torch.randn(batch_size, num_nodes, ch)
    signal_orig = NodeSignal(dom, features)

    # 무작위 순열 적용
    perm = torch.randperm(num_nodes)
    features_perm = features[:, perm, :]
    signal_perm = NodeSignal(dom, features_perm)

    for agg in ["sum", "mean"]:
        pool = GraphGlobalPooling(aggregation=agg)

        with torch.no_grad():
            out_orig = pool(signal_orig)
            out_perm = pool(signal_perm)

        assert torch.allclose(out_orig, out_perm, atol=1e-6), (
            f"GraphGlobalPooling(aggregation='{agg}')이 순열 불변이 아닙니다.\n"
            f"  원본 출력: {out_orig}\n"
            f"  순열 출력: {out_perm}"
        )


def test_graph_global_pooling_max_permutation_invariance(batch_size):
    """
    max 집계 역시 순열 불변인지 검증.
    max 연산은 순서와 무관하다.
    """
    num_nodes = 5
    ch = 4
    torch.manual_seed(123)

    dom = SetDomain(num_nodes=num_nodes)
    features = torch.randn(batch_size, num_nodes, ch)
    signal_orig = NodeSignal(dom, features)

    perm = torch.randperm(num_nodes)
    features_perm = features[:, perm, :]
    signal_perm = NodeSignal(dom, features_perm)

    pool = GraphGlobalPooling(aggregation="max")

    with torch.no_grad():
        out_orig = pool(signal_orig)
        out_perm = pool(signal_perm)

    assert torch.allclose(out_orig, out_perm, atol=1e-6), (
        "GraphGlobalPooling(aggregation='max')이 순열 불변이 아닙니다."
    )
