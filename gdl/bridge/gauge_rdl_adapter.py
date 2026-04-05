"""
RDL → GDL 5G 어댑터: MasterResonantNetwork을 EquivariantLayer[GaugeSignal]로 래핑.

수학적 의미:
GDL의 게이지 등변 레이어는 좌표 표현이 게이지 독립적이어야 한다 (Bronstein eq 23).
RDL의 U(1) 게이지 ODE는 정확히 이 조건을 만족한다:
- 위상 φ는 좌표 선택(게이지)에 의존하지만
- 판별식 F₂ = Im{e^{-iφ}(L̂-ψ̇)} = 0 조건은 게이지 불변

이 어댑터는 RDL의 복소 연산을 GDL의 [B,N,C,2] 접선 벡터 표현으로 변환하여
5G 블루프린트 파이프라인에 삽입할 수 있게 한다.
"""

import torch
import torch.nn as nn

from gdl.core.base_blueprint import EquivariantLayer
from gdl.domains.gauges.signal import GaugeSignal
from gdl.domains.gauges.domain import GaugeDomain
from gdl.rdl.models.master_net import MasterResonantNetwork

from .rdl_gauge_domain import gauge_signal_to_complex, complex_to_gauge_signal


class RDLGaugeAdapter(EquivariantLayer[GaugeSignal]):
    """
    MasterResonantNetwork을 GDL 5G EquivariantLayer로 래핑.

    변환 흐름:
    1. GaugeSignal [B, N, C, 2] → 복소 텐서 [B, N*C]
    2. MasterResonantNetwork 전방향 전파
    3. 복소 출력 → GaugeSignal [B, N, C', 2]

    이 레이어는 GDL 블루프린트의 B (EquivariantLayer) 위치에 삽입 가능하며,
    RDL의 9단계 파이프라인을 5G 게이지 등변 프레임워크 안에서 실행한다.
    """

    def __init__(
        self,
        num_nodes: int,
        in_channels: int,
        out_channels: int,
        hidden_features: int = 32,
        num_layers: int = 3,
        channel_type: str = "paper3ch",
    ) -> None:
        super().__init__()

        self.num_nodes = num_nodes
        self.in_channels = in_channels
        self.out_channels = out_channels

        # RDL 입력 차원 = 노드 수 × 채널 수 × 2 (실수/허수)
        rdl_in = num_nodes * in_channels * 2
        rdl_out = num_nodes * out_channels

        self.rdl_net = MasterResonantNetwork(
            in_features=rdl_in,
            hidden_features=hidden_features,
            out_features=rdl_out,
            num_layers=num_layers,
            channel_type=channel_type,
        )

        # RDL 출력 딕셔너리에서 추출할 키
        self._output_key = "Z_out"

    def forward(self, signal: GaugeSignal) -> GaugeSignal:
        """
        GaugeSignal → RDL 처리 → GaugeSignal.

        Args:
            signal: 입력 GaugeSignal [B, N, C, 2]

        Returns:
            처리된 GaugeSignal [B, N, C', 2]
        """
        domain: GaugeDomain = signal.domain
        features = signal.features  # [B, N, C, 2]
        batch_size = features.size(0)

        # 1. [B, N, C, 2] → [B, N*C*2] (실수 평탄화)
        x_flat = features.reshape(batch_size, -1).to(torch.float64)

        # 2. RDL 전방향 전파
        rdl_out = self.rdl_net(x_flat)
        z_out = rdl_out[self._output_key]  # [B, rdl_out] 복소

        # 3. 복소 출력 → [B, N, C'] → GaugeSignal [B, N, C', 2]
        z_reshaped = z_out.reshape(batch_size, self.num_nodes, self.out_channels)

        return complex_to_gauge_signal(z_reshaped, domain)

    @property
    def last_rdl_outputs(self) -> dict:
        """마지막 forward에서 RDL이 생성한 전체 출력 딕셔너리 접근용."""
        return getattr(self, "_last_rdl_outputs", {})

    def extra_repr(self) -> str:
        return (
            f"nodes={self.num_nodes}, in_ch={self.in_channels}, "
            f"out_ch={self.out_channels}, rdl={self.rdl_net.extra_repr()}"
        )
