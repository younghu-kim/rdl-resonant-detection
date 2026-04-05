import math

import torch

from gdl.core import BaseTransform
from .domain import GaugeDomain
from .signal import GaugeSignal


class GaugeRotation(BaseTransform[GaugeDomain, GaugeSignal]):
    """
    5G Transform: Local Gauge Transformation on the Tangent Bundle (SO(2)^N).

    Since there is no canonical coordinate system (Gauge) on a curved manifold,
    every node 'i' can arbitrarily rotate its local 2D reference frame by an angle theta_i.

    A Gauge Equivariant model must produce the same physical vector predictions
    regardless of how these reference frames are chosen.

    Mathematically, when node 'i' rotates its frame by theta_i (active rotation of axes):
    1. The vector components (x, y) at node 'i' rotate by -theta_i (passive rotation).
    2. The parallel transport angle rho_{ij} (from j to i) shifts by (-theta_i + theta_j)
       because both the source and target frames have moved.
    """

    def __init__(self, gauge_angles: torch.Tensor):
        """
        Args:
            gauge_angles: FloatTensor of shape [num_nodes]. The arbitrary
                          angle theta_i by which node i rotates its local Gauge.
        """
        if gauge_angles.dim() != 1:
            raise ValueError(
                f"gauge_angles must be a 1D tensor, got {gauge_angles.dim()}D."
            )
        self.gauge_angles = gauge_angles

    def forward_domain(self, domain: GaugeDomain) -> GaugeDomain:
        if self.gauge_angles.size(0) != domain.num_nodes:
            raise ValueError(
                f"Number of gauge angles ({self.gauge_angles.size(0)}) must match "
                f"the number of nodes ({domain.num_nodes})."
            )

        device = domain.transport_angles.device
        dtype = domain.transport_angles.dtype
        theta = self.gauge_angles.to(device=device, dtype=dtype)

        # When frame i rotates by theta_i and frame j rotates by theta_j,
        # the transport angle rho_{ij} from j to i must compensate for both changes.
        # edge_index: row 0 is source (j), row 1 is target (i).
        row, col = domain.edge_index

        theta_j = theta[row]
        theta_i = theta[col]

        # New transport angle: rho'_{ij} = rho_{ij} - theta_i + theta_j
        new_transport_angles = domain.transport_angles - theta_i + theta_j

        # Keep angles normalized between [-pi, pi] for stability
        new_transport_angles = (new_transport_angles + math.pi) % (2 * math.pi) - math.pi

        return GaugeDomain(
            num_nodes=domain.num_nodes,
            edge_index=domain.edge_index,
            transport_angles=new_transport_angles,
        )

    def forward_signal(self, signal: GaugeSignal) -> GaugeSignal:
        domain = signal.domain
        device = signal.features.device
        dtype = signal.features.dtype

        theta = self.gauge_angles.to(device=device, dtype=dtype)

        # Passive Rotation: If the coordinate frame rotates counter-clockwise by theta,
        # the coordinates of the physical vector appear to rotate clockwise by -theta.
        c = torch.cos(-theta)
        s = torch.sin(-theta)

        # Construct a [Nodes, 2, 2] rotation matrix for each node
        # R_i = [[cos(-theta_i), -sin(-theta_i)],
        #        [sin(-theta_i),  cos(-theta_i)]]
        R = torch.stack(
            [
                torch.stack([c, -s], dim=-1),
                torch.stack([s, c], dim=-1),
            ],
            dim=-2,
        )

        # Apply the rotation to the vector features.
        # features shape: [Batch, Nodes, Channels, 2]
        # R shape: [Nodes, 2, 2]
        x = signal.features

        # Einsum: b: batch, n: nodes, c: channels, i/j: 2D vector components
        rotated_features = torch.einsum("nij, bncj -> bnci", R, x)

        # The signal must ride on the newly rotated domain geometry!
        deformed_domain = self.forward_domain(domain)

        return GaugeSignal(
            domain=deformed_domain,
            features=rotated_features,
        )
