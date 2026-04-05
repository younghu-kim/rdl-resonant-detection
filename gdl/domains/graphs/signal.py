import torch

from gdl.core.base_signal import BaseSignal
from .domain import SetDomain


class NodeSignal(BaseSignal):
    """
    3G Signal: Features living on the nodes of a Set or Graph.

    Tensor Shape strictly enforced: [Batch, NumNodes, Channels]
    (Note: This differs from 1G GridSignal [B, C, Spatial] to perfectly align
    with standard Transformers and GNN node-as-sequence paradigms).
    """

    def __init__(self, domain: SetDomain, features: torch.Tensor) -> None:
        # SetDomain is the parent of GraphDomain, so this accepts both.
        super().__init__(domain, features)
        self._validate_shape()

    def _validate_shape(self) -> None:
        """Ensures the feature tensor mathematically matches the node count."""
        set_domain: SetDomain = self.domain

        if self.features.dim() != 3:
            raise ValueError(
                f"NodeSignal expects 3D tensor [Batch, Nodes, Channels], "
                f"got {self.features.dim()}D."
            )

        actual_nodes = self.features.size(1)
        if actual_nodes != set_domain.num_nodes:
            raise ValueError(
                f"Geometry mismatch: Domain declares {set_domain.num_nodes} nodes, "
                f"but signal features contain {actual_nodes} nodes."
            )

    @property
    def domain(self) -> SetDomain:
        return self._domain  # type: ignore[return-value]

    def to(self, device: torch.device | str) -> "NodeSignal":
        return NodeSignal(
            domain=self.domain.to(device),
            features=self.features.to(device),
        )
