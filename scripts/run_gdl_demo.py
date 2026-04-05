#!/usr/bin/env python3
"""
=============================================================================
GDL Unified Framework Demo — 8개 기하학 도메인 전체 순회
=============================================================================
GDLFactory로 8가지 도메인을 자동 조립하고 각각 3 에폭 훈련하여
Bronstein et al.의 Common Blueprint가 모든 기하학에서 동작함을 검증한다.

실행: python scripts/run_gdl_demo.py
"""

import os
import sys
import time

import torch

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_ROOT)

from gdl.core.factory import GDLFactory
from gdl.core.base_blueprint import BaseGDLModel

# 도메인별 시그널 / 도메인 클래스
from gdl.domains.grids import GridDomain, GridSignal
from gdl.domains.groups import GroupDomain, GroupSignal
from gdl.domains.graphs import SetDomain, GraphDomain, NodeSignal
from gdl.domains.geodesics import GeodesicDomain, GeodesicSignal
from gdl.domains.gauges import GaugeDomain, GaugeSignal
from gdl.domains.temporal import TimeDomain, TimeSignal


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

    elif domain_type in ("3G_GRAPH", "3G_ATTENTION"):
        ei = torch.tensor([[0, 1, 1, 2, 2, 3, 3, 0],
                           [1, 0, 2, 1, 3, 2, 0, 3]], dtype=torch.long)
        ew = torch.ones(8)
        dom = GraphDomain(num_nodes=4, edge_index=ei, edge_weight=ew)
        feat = torch.randn(batch, 4, in_ch)
        return NodeSignal(dom, feat), target

    elif domain_type == "4G_GEODESIC":
        verts = torch.tensor([
            [1., 1., 1.], [1., -1., -1.], [-1., 1., -1.], [-1., -1., 1.]
        ])
        faces = torch.tensor([[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]], dtype=torch.long)
        dom = GeodesicDomain(vertices=verts, faces=faces)
        feat = torch.randn(batch, 4, in_ch)
        return GeodesicSignal(dom, feat), target

    elif domain_type == "5G_GAUGE":
        ei = torch.tensor([[0, 1, 1, 2, 2, 3, 3, 0],
                           [1, 0, 2, 1, 3, 2, 0, 3]], dtype=torch.long)
        angles = torch.randn(8)
        dom = GaugeDomain(num_nodes=4, edge_index=ei, transport_angles=angles)
        feat = torch.randn(batch, 4, in_ch, 2)
        return GaugeSignal(dom, feat), target

    elif domain_type == "TIME":
        dom = TimeDomain(sequence_length=6)
        feat = torch.randn(batch, 6, in_ch)
        return TimeSignal(dom, feat), target

    raise ValueError(f"Unknown: {domain_type}")


DOMAIN_TYPES = [
    "1G_GRID", "2G_GROUP", "3G_SET", "3G_GRAPH",
    "3G_ATTENTION", "4G_GEODESIC", "5G_GAUGE", "TIME",
]


def main():
    print("=" * 70)
    print("  GDL Unified Framework — 8 Domain End-to-End Demo")
    print("=" * 70)

    batch, in_ch, hid_ch, out_ch = 2, 3, 8, 4
    results = {}

    for dtype in DOMAIN_TYPES:
        print(f"\n--- {dtype} ---")
        t0 = time.time()

        try:
            model = GDLFactory.create_model(
                domain_type=dtype,
                in_channels=in_ch,
                hidden_channels=hid_ch,
                out_channels=out_ch,
                num_clusters=2, num_heads=1, K=2,
                num_spatial_dims=2, kernel_size=3,
                pool_kernel=2, pool_stride=2,
                l_max=1, aggr="mean",
            )

            signal, target = make_dummy_data(dtype, batch, in_ch, out_ch)
            optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
            criterion = torch.nn.MSELoss()

            losses = []
            for step in range(3):
                optimizer.zero_grad()
                pred = model(signal)
                loss = criterion(pred, target)
                loss.backward()
                optimizer.step()
                losses.append(loss.item())

            has_grad = any(p.grad is not None for p in model.parameters())
            dt = time.time() - t0

            n_params = sum(p.numel() for p in model.parameters())
            status = "OK" if has_grad and not any(torch.isnan(torch.tensor(l)) for l in losses) else "FAIL"
            results[dtype] = status
            print(f"  Params: {n_params:,} | Loss: {losses[-1]:.4f} | Grad: {has_grad} | {dt:.2f}s | {status}")

        except Exception as e:
            results[dtype] = f"ERROR: {e}"
            print(f"  ERROR: {e}")

    # 최종 요약
    print(f"\n{'=' * 70}")
    print("  Summary")
    print(f"{'=' * 70}")
    ok = sum(1 for v in results.values() if v == "OK")
    for d, s in results.items():
        print(f"  {d:20s} {s}")
    print(f"\n  {ok}/{len(DOMAIN_TYPES)} domains passed")


if __name__ == "__main__":
    main()
