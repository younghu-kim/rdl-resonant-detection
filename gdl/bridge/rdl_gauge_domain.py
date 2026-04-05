"""
RDL ↔ GDL 5G 도메인 매핑.

수학적 대응:
- RDL U(1) 위상 회전 = GDL SO(2) 구조군
- RDL 게이지 ODE 적분 = GDL 평행 이동 R(ρ_ij)
- RDL 복소 선형 연산 = GDL Schur 가중치 [[a,-b],[b,a]]
- RDL 진폭-위상 분리 = GDL 노름 게이팅

이 모듈은 RDL의 1D 시간축 격자를 GaugeDomain 그래프로 변환한다.
시간 격자의 인접 노드 간 위상차가 평행 이동 각도가 된다.
"""

import torch
import math

from gdl.domains.gauges.domain import GaugeDomain
from gdl.domains.gauges.signal import GaugeSignal


def time_grid_to_gauge_domain(
    num_points: int,
    t_min: float = 10.0,
    t_max: float = 50.0,
    device: torch.device | str = "cpu",
) -> GaugeDomain:
    """
    1D 시간 격자를 5G GaugeDomain으로 변환.

    시간 격자의 각 점이 노드, 인접 점 간 연결이 엣지,
    시간 간격에 비례하는 위상차가 평행 이동 각도가 된다.

    Args:
        num_points: 격자 점 수 (= 노드 수)
        t_min: 시간 범위 시작
        t_max: 시간 범위 끝
        device: 연산 장치

    Returns:
        GaugeDomain: 1D 체인 그래프 위의 게이지 도메인
    """
    # 1D 체인 그래프: 0-1, 1-2, ..., (N-2)-(N-1) 양방향
    sources = []
    targets = []
    for i in range(num_points - 1):
        sources.extend([i, i + 1])
        targets.extend([i + 1, i])

    edge_index = torch.tensor([sources, targets], dtype=torch.long, device=device)

    # 평행 이동 각도: dt에 비례하는 U(1) 위상차
    # RDL에서 게이지 ODE가 시간 간격을 따라 위상을 적분하는 것에 대응
    dt = (t_max - t_min) / (num_points - 1)
    base_angle = dt * math.pi / (t_max - t_min)  # 정규화된 위상 증분

    # 양방향: 정방향 +angle, 역방향 -angle
    angles = []
    for _ in range(num_points - 1):
        angles.extend([base_angle, -base_angle])

    transport_angles = torch.tensor(angles, dtype=torch.float64, device=device)

    return GaugeDomain(
        num_nodes=num_points,
        edge_index=edge_index,
        transport_angles=transport_angles,
    )


def complex_to_gauge_signal(
    z: torch.Tensor,
    domain: GaugeDomain,
) -> GaugeSignal:
    """
    복소 텐서 [B, N, C] → GaugeSignal [B, N, C, 2].

    RDL의 복소 특성을 GDL 5G의 2D 접선 벡터로 변환.
    실수부 = x 성분, 허수부 = y 성분.

    Args:
        z: 복소 텐서 [Batch, Nodes, Channels]
        domain: GaugeDomain

    Returns:
        GaugeSignal: [B, N, C, 2] 접선 벡터장
    """
    if z.is_complex():
        features = torch.stack([z.real, z.imag], dim=-1)
    else:
        # 이미 실수인 경우 y=0
        features = torch.stack([z, torch.zeros_like(z)], dim=-1)

    return GaugeSignal(domain=domain, features=features.to(torch.float32))


def gauge_signal_to_complex(signal: GaugeSignal) -> torch.Tensor:
    """
    GaugeSignal [B, N, C, 2] → 복소 텐서 [B, N, C].

    GDL 5G의 2D 접선 벡터를 RDL의 복소 특성으로 변환.

    Args:
        signal: GaugeSignal [B, N, C, 2]

    Returns:
        복소 텐서 [B, N, C]
    """
    x = signal.features[..., 0]  # 실수부
    y = signal.features[..., 1]  # 허수부
    return torch.complex(x.to(torch.float64), y.to(torch.float64))
