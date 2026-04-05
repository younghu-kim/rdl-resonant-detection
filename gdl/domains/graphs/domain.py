import torch

from gdl.core.base_domain import BaseDomain


class SetDomain(BaseDomain):
    """
    3G Domain (Sets): A collection of N discrete nodes with NO edges.
    Symmetry Group: Permutation Group S_n.
    """

    def __init__(self, num_nodes: int) -> None:
        if num_nodes <= 0:
            raise ValueError("num_nodes must be strictly positive.")
        self.num_nodes = num_nodes

    def to(self, device: torch.device | str) -> "SetDomain":
        # Sets have no heavy structural tensors.
        return self

    def __repr__(self) -> str:
        return f"SetDomain(num_nodes={self.num_nodes})"


class GraphDomain(SetDomain):
    """
    3G Domain (Graphs): A collection of N nodes connected by edges.
    Unlike Grids, the geometry is explicitly defined by an Adjacency Matrix (edge_index).
    """

    def __init__(
        self,
        num_nodes: int,
        edge_index: torch.Tensor,
        edge_weight: torch.Tensor | None = None,
    ) -> None:
        super().__init__(num_nodes)

        if edge_index.dim() != 2 or edge_index.size(0) != 2:
            raise ValueError(
                f"edge_index must have shape [2, num_edges], got {edge_index.shape}"
            )

        self.edge_index = edge_index
        self.edge_weight = edge_weight

    def to(self, device: torch.device | str) -> "GraphDomain":
        # CRUCIAL: Graphs have heavy structural tensors that MUST move to the GPU!
        # This solves the classic silent bug where features go to GPU but adjacency stays on CPU.
        return GraphDomain(
            num_nodes=self.num_nodes,
            edge_index=self.edge_index.to(device),
            edge_weight=self.edge_weight.to(device) if self.edge_weight is not None else None,
        )

    def __repr__(self) -> str:
        num_edges = self.edge_index.size(1)
        return f"GraphDomain(num_nodes={self.num_nodes}, num_edges={num_edges})"
