import torch

from gdl.core.base_transform import BaseTransform
from .domain import SetDomain, GraphDomain
from .signal import NodeSignal


class NodePermutation(BaseTransform[SetDomain, NodeSignal]):
    """
    Action of the Permutation Symmetry Group S_n on Sets and Graphs.

    If the order of nodes in memory changes, the geometric structure itself
    is perfectly preserved (Isomorphism). Deep Learning models on these domains
    must be Equivariant to this transformation.
    """

    def __init__(self, perm: torch.Tensor) -> None:
        """
        Args:
            perm: A 1D tensor containing a valid permutation of indices from 0 to N-1.
        """
        if perm.dim() != 1:
            raise ValueError(f"Permutation tensor must be 1D, got {perm.dim()}D.")
        self.perm = perm

    def forward_domain(self, domain: SetDomain) -> SetDomain:
        """
        Applies the permutation to the spatial structure.
        For Sets, nothing structurally changes.
        For Graphs, the adjacency matrix (edge_index) MUST be remapped to
        reflect the new node indices.
        """
        if domain.num_nodes != self.perm.size(0):
            raise ValueError(
                f"Permutation length ({self.perm.size(0)}) does not match "
                f"domain node count ({domain.num_nodes})."
            )

        if isinstance(domain, GraphDomain):
            # To remap edges, we need the INVERSE permutation mapping.
            inv_perm = torch.empty_like(self.perm)
            inv_perm[self.perm] = torch.arange(
                self.perm.size(0), device=self.perm.device, dtype=self.perm.dtype,
            )

            # Remap edge_index values to their new node indices
            new_edge_index = inv_perm[domain.edge_index]

            return GraphDomain(
                num_nodes=domain.num_nodes,
                edge_index=new_edge_index,
                edge_weight=domain.edge_weight,
            )

        # For pure Sets (No edges), the abstract domain is unaltered
        return SetDomain(num_nodes=domain.num_nodes)

    def forward_signal(self, signal: NodeSignal) -> NodeSignal:
        """
        Shuffles the node features according to the permutation matrix.
        Simulates rho(g)x(u) = x(g^{-1} u).
        """
        domain = signal.domain

        if domain.num_nodes != self.perm.size(0):
            raise ValueError(
                f"Permutation length ({self.perm.size(0)}) does not match "
                f"signal node count ({domain.num_nodes})."
            )

        # Reorder features along the Nodes dimension (dim=1)
        # Expected shape: [Batch, Nodes, Channels]
        permuted_features = signal.features[:, self.perm, :]

        # We must bind these shuffled features to the newly shuffled domain
        permuted_domain = self.forward_domain(domain)

        return NodeSignal(domain=permuted_domain, features=permuted_features)
