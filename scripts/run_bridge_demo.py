#!/usr/bin/env python3
"""
=============================================================================
RDL ↔ GDL 5G Bridge Demo
=============================================================================
RDL의 MasterResonantNetwork을 GDL 5G Gauge 프레임워크 안에서 실행하여
두 시스템의 수학적 통합을 검증한다.

수학적 대응:
  GDL SO(2) 구조군 = RDL U(1) 위상 회전
  GDL 평행 이동 = RDL 게이지 ODE 적분
  GDL Schur 가중치 = RDL 복소 선형 연산

실행: python scripts/run_bridge_demo.py
"""

import os
import sys

import torch

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, PROJECT_ROOT)

from gdl.rdl.constants import PrecisionManager
from gdl.bridge.rdl_gauge_domain import (
    time_grid_to_gauge_domain,
    complex_to_gauge_signal,
    gauge_signal_to_complex,
)
from gdl.bridge.gauge_rdl_adapter import RDLGaugeAdapter
from gdl.domains.gauges.domain import GaugeDomain
from gdl.domains.gauges.signal import GaugeSignal


def main():
    PrecisionManager.setup_precision()

    print("=" * 70)
    print("  RDL <-> GDL 5G Bridge Demo")
    print("=" * 70)

    # 1. 시간 격자를 GaugeDomain으로 변환
    num_nodes = 8
    domain = time_grid_to_gauge_domain(num_nodes, t_min=10.0, t_max=50.0)
    print(f"\n[1] GaugeDomain 생성: {domain}")
    print(f"    노드: {domain.num_nodes}, 엣지: {domain.edge_index.size(1)}")
    print(f"    평행이동 각도 범위: [{domain.transport_angles.min():.4f}, {domain.transport_angles.max():.4f}]")

    # 2. 더미 복소 데이터 → GaugeSignal
    batch_size = 2
    in_channels = 3
    out_channels = 2
    z_input = torch.randn(batch_size, num_nodes, in_channels, dtype=torch.complex128)
    signal = complex_to_gauge_signal(z_input, domain)
    print(f"\n[2] GaugeSignal 생성: {signal.features.shape}")
    print(f"    복소 → 접선벡터 변환 확인: z[0,0,0] = {z_input[0,0,0]:.4f}")
    print(f"    features[0,0,0,:] = {signal.features[0,0,0,:]}")

    # 3. RDL 어댑터 생성 및 전방향 전파
    adapter = RDLGaugeAdapter(
        num_nodes=num_nodes,
        in_channels=in_channels,
        out_channels=out_channels,
        hidden_features=16,
        num_layers=2,
    )

    n_params = sum(p.numel() for p in adapter.parameters())
    print(f"\n[3] RDLGaugeAdapter 생성: {n_params:,} parameters")

    out_signal = adapter(signal)
    print(f"    입력: {signal.features.shape} → 출력: {out_signal.features.shape}")

    # 4. 출력을 다시 복소수로 변환
    z_output = gauge_signal_to_complex(out_signal)
    print(f"\n[4] 출력 복소 변환: {z_output.shape}, dtype={z_output.dtype}")

    # 5. 미분 가능성 검증
    print(f"\n[5] 미분 검증...")
    loss = out_signal.features.abs().mean()
    loss.backward()

    has_grad = any(
        p.grad is not None and not torch.isnan(p.grad).any()
        for p in adapter.parameters()
    )
    print(f"    Gradient flow: {'OK' if has_grad else 'BROKEN'}")

    # 6. 도메인 보존 검증
    same_domain = (out_signal.domain.num_nodes == signal.domain.num_nodes)
    print(f"    Domain preserved: {'OK' if same_domain else 'FAIL'}")

    print(f"\n{'=' * 70}")
    status = "SUCCESS" if has_grad and same_domain else "FAIL"
    print(f"  Bridge Demo: {status}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
