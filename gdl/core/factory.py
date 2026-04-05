"""
GDL Factory: Auto-assembles BaseGDLModel pipelines for all 8 geometric domains.

Demonstrates Bronstein et al.'s "Common Blueprint" — regardless of geometry,
architecture is universally:
[EquivariantLayer -> NonLinearity -> (Optional) LocalPooling] x N -> GlobalInvariant
"""

from __future__ import annotations

import typing

from .base_blueprint import BaseGDLModel

# 1G: Grids
from gdl.domains.grids.layers import (
    GridConvolution,
    GridNonLinearity,
    GridPooling,
    GridGlobalAveragePooling,
)
# 2G: Groups
from gdl.domains.groups.layers import (
    GroupConvolution,
    GroupNonLinearity,
    GroupPooling,
    GroupGlobalPooling,
)
# 3G: Graphs and Sets
from gdl.domains.graphs.layers import (
    DeepSetsLayer,
    GraphMessagePassingLayer,
    TransformerEquivariantLayer,
    GraphNonLinearity,
    GraphCoarsening,
    GraphGlobalPooling,
)
# 4G: Geodesics
from gdl.domains.geodesics.layers import (
    IntrinsicMeshConvolution,
    GeodesicNonLinearity,
    MeshDecimationPooling,
    GeodesicGlobalPooling,
)
# 5G: Gauges
from gdl.domains.gauges.layers import (
    GaugeMessagePassingLayer,
    GaugeNonLinearity,
    GaugeGlobalPooling,
)
# Time
from gdl.domains.temporal.layers import (
    TimeDifferenceGatedRNN,
    TimeNonLinearity,
    TimeGlobalPooling,
)


class GDLFactory:
    """
    GDL 통합 프레임워크의 최종 시스템 팩토리.

    도메인 타입 문자열과 채널 차원만으로 완전한 신경망 파이프라인을 자동 조립한다.

    지원 도메인:
    - "1G_GRID": 평행이동 등변 (CNN)
    - "2G_GROUP": 회전 등변 (Spherical CNN)
    - "3G_SET": 순열 등변 (Deep Sets)
    - "3G_GRAPH": 순열 등변 + 엣지 (GNN)
    - "3G_ATTENTION": 순열 등변 + 밀집 엣지 (Transformer)
    - "4G_GEODESIC": 등거리 등변 (Intrinsic Mesh CNN)
    - "5G_GAUGE": 게이지 등변 (Directional Gauge CNN)
    - "TIME": 시간 워핑 등변 (Difference-Gated RNN)
    """

    @staticmethod
    def create_model(
        domain_type: str,
        in_channels: int,
        hidden_channels: int,
        out_channels: int,
        **kwargs: typing.Any,
    ) -> BaseGDLModel:
        """도메인 타입에 따라 BaseGDLModel 파이프라인을 자동 조립."""

        blocks: list[typing.Any] = []
        readout: typing.Any = None

        domain_type = domain_type.upper()

        if domain_type == "1G_GRID":
            num_spatial_dims = kwargs.get("num_spatial_dims", 2)
            kernel_size = kwargs.get("kernel_size", 3)
            pool_kernel = kwargs.get("pool_kernel", 2)
            pool_stride = kwargs.get("pool_stride", 2)
            blocks = [
                GridConvolution(in_channels, hidden_channels,
                                kernel_size=kernel_size,
                                num_spatial_dims=num_spatial_dims),
                GridNonLinearity(),
                GridPooling(kernel_size=pool_kernel, stride=pool_stride,
                            num_spatial_dims=num_spatial_dims),
                GridConvolution(hidden_channels, out_channels,
                                kernel_size=kernel_size,
                                num_spatial_dims=num_spatial_dims),
                GridNonLinearity(),
            ]
            readout = GridGlobalAveragePooling()

        elif domain_type == "2G_GROUP":
            l_max = kwargs.get("l_max", 1)
            blocks = [
                GroupConvolution(in_channels, hidden_channels, l_max=l_max),
                GroupNonLinearity(),
                GroupPooling(target_l_max=0),
                GroupConvolution(hidden_channels, out_channels, l_max=0),
                GroupNonLinearity(),
            ]
            readout = GroupGlobalPooling(use_power_spectrum=True)

        elif domain_type == "3G_SET":
            blocks = [
                DeepSetsLayer(in_channels, hidden_channels),
                GraphNonLinearity(),
                DeepSetsLayer(hidden_channels, out_channels),
                GraphNonLinearity(),
            ]
            readout = GraphGlobalPooling(aggregation="sum")

        elif domain_type == "3G_GRAPH":
            num_clusters = kwargs.get("num_clusters", 2)
            blocks = [
                GraphMessagePassingLayer(in_channels, hidden_channels),
                GraphNonLinearity(),
                GraphCoarsening(in_channels=hidden_channels, num_clusters=num_clusters),
                # GraphCoarsening 후 SetDomain이므로 DeepSetsLayer 사용
                DeepSetsLayer(hidden_channels, out_channels),
                GraphNonLinearity(),
            ]
            readout = GraphGlobalPooling(aggregation="mean")

        elif domain_type == "3G_ATTENTION":
            num_heads = kwargs.get("num_heads", 1)
            # Transformer는 embed_dim 보존 → DeepSets로 차원 변환
            blocks = [
                TransformerEquivariantLayer(embed_dim=in_channels, num_heads=num_heads),
                GraphNonLinearity(),
                DeepSetsLayer(in_channels, out_channels),
                GraphNonLinearity(),
            ]
            readout = GraphGlobalPooling(aggregation="mean")

        elif domain_type == "4G_GEODESIC":
            K = kwargs.get("K", 2)
            target_nodes = kwargs.get("target_nodes", 2)
            blocks = [
                IntrinsicMeshConvolution(in_channels, hidden_channels, K=K),
                GeodesicNonLinearity(),
                MeshDecimationPooling(in_channels=hidden_channels,
                                      target_nodes=target_nodes),
                IntrinsicMeshConvolution(hidden_channels, out_channels, K=K),
                GeodesicNonLinearity(),
            ]
            readout = GeodesicGlobalPooling(aggregation="max")

        elif domain_type == "5G_GAUGE":
            aggr = kwargs.get("aggr", "mean")
            blocks = [
                GaugeMessagePassingLayer(in_channels, hidden_channels, aggr=aggr),
                GaugeNonLinearity(),
                GaugeMessagePassingLayer(hidden_channels, out_channels, aggr=aggr),
                GaugeNonLinearity(),
            ]
            readout = GaugeGlobalPooling(aggregation="mean")

        elif domain_type == "TIME":
            blocks = [
                TimeDifferenceGatedRNN(in_channels, hidden_channels),
                TimeNonLinearity(),
                TimeDifferenceGatedRNN(hidden_channels, out_channels),
                TimeNonLinearity(),
            ]
            readout = TimeGlobalPooling(aggregation="last")

        else:
            raise ValueError(f"Unknown GDL domain type: {domain_type}")

        return BaseGDLModel(blocks=blocks, readout=readout)
