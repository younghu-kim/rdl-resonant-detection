import torch
import torch.nn as nn

from gdl.core import (
    EquivariantLayer,
    GeometricNonLinearity,
    GlobalInvariant,
)
from .domain import GaugeDomain
from .signal import GaugeSignal


class ParallelTransport(nn.Module):
    """
    5G Utility: Parallel Transport of Tangent Vectors along Edges.

    To compare or aggregate vectors from neighboring nodes (j) to a central node (i)
    on a curved manifold, we cannot simply add them because their local coordinate
    frames (Gauges) are completely independent and misaligned.

    We must mathematically "transport" the vector from j's tangent space to i's tangent
    space using the Levi-Civita connection (represented by the transport angles rho_{ij}).
    This active rotation rigorously aligns the neighbor's vector into the center node's frame.
    """

    def __init__(self) -> None:
        super().__init__()

    def forward(self, signal: GaugeSignal) -> torch.Tensor:
        """
        Transports all neighbor vectors X_j to the local coordinate frames of their targets X_i.

        Returns:
            transported_j: Tensor of shape [Batch, NumEdges, Channels, 2]
                           representing the aligned neighbor messages ready for aggregation.
        """
        domain: GaugeDomain = signal.domain

        # x shape: [Batch, Nodes, Channels, 2]
        x = signal.features

        # edge_index shape: [2, NumEdges]. row: source (j), col: target (i)
        row, col = domain.edge_index

        # 1. Gather the neighbor vectors (X_j) for all edges
        # Shape: [Batch, NumEdges, Channels, 2]
        x_j = x[:, row, :, :]

        # 2. Extract the transport angles rho_{ij}
        # Shape: [NumEdges]
        rho = domain.transport_angles

        # 3. Construct the 2D Rotation Matrix R(rho_{ij}) for each edge
        # Active rotation to physically align the tangent spaces across the curved manifold
        c = torch.cos(rho)
        s = torch.sin(rho)

        # R shape: [NumEdges, 2, 2]
        R = torch.stack(
            [
                torch.stack([c, -s], dim=-1),
                torch.stack([s, c], dim=-1),
            ],
            dim=-2,
        )

        # 4. Apply Parallel Transport (Rotation) to X_j
        # e: edges, b: batch, c: channels, i/j: 2D vector components
        transported_j = torch.einsum("eij, becj -> beci", R, x_j)

        return transported_j


class GaugeNonLinearity(GeometricNonLinearity[GaugeSignal]):
    """
    5G Non-Linearity: Norm-based gating for Tangent Vectors.

    Applying a standard ReLU directly to the (x, y) components of a tangent vector
    destroys the physical direction of the vector and instantly breaks Gauge Equivariance.

    Instead, we compute the magnitude (L2 Norm) of the vector, apply ReLU to the norm
    to create a scalar gate, and multiply it back. This strictly preserves the
    original 2D orientation while providing necessary deep learning non-linearity.
    """

    def __init__(self) -> None:
        super().__init__()
        self.relu = nn.ReLU()

    def forward(self, signal: GaugeSignal) -> GaugeSignal:
        # x shape: [Batch, Nodes, Channels, 2]
        x = signal.features

        # Calculate the magnitude (length) of the 2D vector along the last dimension
        # Shape: [Batch, Nodes, Channels, 1]
        norm = torch.norm(x, p=2, dim=-1, keepdim=True)

        # Avoid division by zero
        safe_norm = norm + 1e-8

        # Apply non-linearity strictly to the invariant scalar magnitude
        gate = self.relu(norm)

        # Scale the original unit vector (direction is perfectly preserved)
        out_features = x * (gate / safe_norm)

        return GaugeSignal(domain=signal.domain, features=out_features)


class GaugeMessagePassingLayer(EquivariantLayer[GaugeSignal]):
    """
    5G Equivariant Layer: Directional Message Passing on Gauge Bundles (Gauge CNN).

    Guarantees Gauge Equivariance (SO(2)^N) by:
    1. Using Parallel Transport to mathematically align neighbor vectors (j)
       into the local tangent space of the central node (i).
    2. Restricting the learnable weight matrices to commute with SO(2) rotations.
       According to Schur's Lemma for SO(2), the weight must act like complex
       multiplication, parameterized as a 2x2 block matrix [[a, -b], [b, a]].
    3. Strictly prohibiting any vector bias additions, as an arbitrary 2D bias
       does not rotate with the Gauge, destroying the physical symmetry.
    """

    def __init__(self, in_channels: int, out_channels: int, aggr: str = "mean"):
        super().__init__()

        valid_aggs = ["sum", "mean"]
        if aggr not in valid_aggs:
            raise ValueError(f"Unsupported aggregation: {aggr}")
        self.aggr = aggr
        self.out_channels = out_channels

        # Parallel Transport module to align coordinate frames across edges
        self.transport = ParallelTransport()

        # We parameterize the SO(2)-commuting weights as 'a' (scale/real) and 'b' (rotation/imag)
        std = in_channels**-0.5
        self.w_self_a = nn.Parameter(torch.randn(out_channels, in_channels) * std)
        self.w_self_b = nn.Parameter(torch.randn(out_channels, in_channels) * std)

        self.w_neigh_a = nn.Parameter(torch.randn(out_channels, in_channels) * std)
        self.w_neigh_b = nn.Parameter(torch.randn(out_channels, in_channels) * std)

        # CRUCIAL GEOMETRIC CONSTRAINT:
        # No bias term! Adding a fixed [x, y] vector breaks gauge symmetry.

    def _build_so2_weight(self, a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
        """
        Constructs a complex-like 2x2 block weight matrix that commutes with SO(2).
        Returns shape: [C_out, C_in, 2, 2]
        """
        # Row 1: [a, -b]
        row1 = torch.stack([a, -b], dim=-1)
        # Row 2: [b, a]
        row2 = torch.stack([b, a], dim=-1)
        # Combine: [[a, -b], [b, a]]
        return torch.stack([row1, row2], dim=-2)

    def forward(self, signal: GaugeSignal) -> GaugeSignal:
        domain: GaugeDomain = signal.domain

        # x shape: [Batch, Nodes, C_in, 2]
        x = signal.features
        batch_size, num_nodes, in_channels, _ = x.shape

        # 1. Self-Connection Update (Transforming node i's own vectors)
        W_self = self._build_so2_weight(self.w_self_a, self.w_self_b)
        # Einsum: o(C_out), c(C_in), b(Batch), n(Nodes), i/j(2D dims)
        out_self = torch.einsum("ocij, bncj -> bnoi", W_self, x)

        # 2. Message Function: Align neighbor vectors into target's tangent space
        row, col = domain.edge_index
        if row.size(0) == 0:
            return GaugeSignal(domain=domain, features=out_self)

        # transported_j shape: [Batch, NumEdges, C_in, 2]
        transported_j = self.transport(signal)

        # Transform the aligned neighbor messages using SO(2)-commuting weights
        W_neigh = self._build_so2_weight(self.w_neigh_a, self.w_neigh_b)
        msg_j = torch.einsum("ocij, becj -> beoi", W_neigh, transported_j)

        # 3. Aggregation Function: Scatter add the messages to the target nodes
        aggr_msg = torch.zeros(
            batch_size,
            num_nodes,
            self.out_channels,
            2,
            device=x.device,
            dtype=x.dtype,
        )

        # col shape: [NumEdges] -> expand to match [Batch, NumEdges, C_out, 2]
        scatter_idx = col.view(1, -1, 1, 1).expand(
            batch_size, -1, self.out_channels, 2
        )
        aggr_msg.scatter_add_(dim=1, index=scatter_idx, src=msg_j)

        if self.aggr == "mean":
            degree = torch.zeros(num_nodes, device=x.device, dtype=x.dtype)
            ones = torch.ones_like(col, dtype=x.dtype)
            degree.scatter_add_(dim=0, index=col, src=ones)
            degree = degree.clamp(min=1.0).view(1, -1, 1, 1)
            aggr_msg = aggr_msg / degree

        # 4. Final Update
        out_features = out_self + aggr_msg

        return GaugeSignal(domain=domain, features=out_features)


class GaugeGlobalPooling(GlobalInvariant[GaugeSignal]):
    """
    5G Global Invariant: Destroys the Tangent Bundle geometry to extract Gauge-Invariant features.

    The local (x, y) components of a tangent vector are entirely dependent on the arbitrary
    choice of Gauge, so they are physically meaningless on a global scale.
    However, the length (L2 Norm) of the vector is an absolute physical invariant
    (Gauge Invariant) because SO(2) rotations are orthogonal matrices.
    """

    def __init__(self, aggregation: str = "mean"):
        super().__init__()
        valid_aggs = ["sum", "mean", "max"]
        if aggregation not in valid_aggs:
            raise ValueError(f"Unsupported aggregation: {aggregation}")
        self.aggregation = aggregation

    def forward(self, signal: GaugeSignal) -> torch.Tensor:
        x = signal.features  # Shape: [Batch, Nodes, Channels, 2]

        # 1. Extract the Gauge-Invariant physical property: Vector Magnitude (L2 Norm)
        # Shape: [Batch, Nodes, Channels]
        invariant_norm = torch.norm(x, p=2, dim=-1)

        # 2. Pool across the spatial graph structure (dim=1)
        if self.aggregation == "sum":
            invariant_tensor = torch.sum(invariant_norm, dim=1)
        elif self.aggregation == "mean":
            invariant_tensor = torch.mean(invariant_norm, dim=1)
        elif self.aggregation == "max":
            invariant_tensor = torch.max(invariant_norm, dim=1)[0]
        else:
            raise RuntimeError("Unreachable aggregation state.")

        # Final Shape: [Batch, Channels]
        return invariant_tensor
