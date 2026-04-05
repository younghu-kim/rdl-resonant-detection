import torch

from gdl.core import BaseDomain


class GaugeDomain(BaseDomain):
    """
    5G Domain: Gauge Bundles and Tangent Spaces on Manifolds.

    Unlike 4G which deals with scalar fields (like temperature), 5G deals with
    Vector Fields (like wind direction) living in the 2D tangent plane of each node.

    To represent vectors, each node requires an arbitrary local coordinate frame (a "Gauge").
    To exchange vectors between neighboring nodes on a curved surface, we must align
    their coordinate frames using a rotation called "Parallel Transport".
    """

    def __init__(
        self,
        num_nodes: int,
        edge_index: torch.Tensor,
        transport_angles: torch.Tensor,
    ):
        """
        Args:
            num_nodes: Total number of points on the manifold.
            edge_index: LongTensor of shape [2, num_edges]. row: source (j), col: target (i).
            transport_angles: FloatTensor of shape [num_edges]. The angle rho_{ij}
                              required to parallel transport a vector from node j's
                              tangent space to node i's tangent space.
        """
        if edge_index.dim() != 2 or edge_index.size(0) != 2:
            raise ValueError(
                f"edge_index must have shape [2, num_edges], got {edge_index.shape}"
            )

        num_edges = edge_index.size(1)
        if transport_angles.dim() != 1 or transport_angles.size(0) != num_edges:
            raise ValueError(
                f"transport_angles must have shape [{num_edges}], "
                f"got {transport_angles.shape}"
            )

        self.num_nodes = num_nodes
        self.edge_index = edge_index
        self.transport_angles = transport_angles

    def to(self, device: torch.device | str) -> "GaugeDomain":
        return GaugeDomain(
            num_nodes=self.num_nodes,
            edge_index=self.edge_index.to(device),
            transport_angles=self.transport_angles.to(device),
        )

    def __repr__(self) -> str:
        return f"GaugeDomain(nodes={self.num_nodes}, edges={self.edge_index.size(1)})"
