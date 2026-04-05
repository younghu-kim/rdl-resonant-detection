import torch

from gdl.core import BaseSignal
from .domain import GaugeDomain


class GaugeSignal(BaseSignal):
    """
    5G Signal: Vector Fields living in local Tangent Spaces.

    Strictly enforces the tensor shape convention: [Batch, Nodes, Channels, 2]
    The final dimension of size 2 explicitly represents the (x, y) components
    of the 2D vector defined relative to each node's local Gauge frame.
    """

    def __init__(self, domain: GaugeDomain, features: torch.Tensor):
        super().__init__(domain, features)
        self._validate_shape()

    def _validate_shape(self) -> None:
        """Ensures the feature tensor strictly represents 2D tangent vectors."""
        gauge_domain: GaugeDomain = self.domain

        if self.features.dim() != 4:
            raise ValueError(
                f"GaugeSignal expects a 4D tensor [Batch, Nodes, Channels, 2], "
                f"got {self.features.dim()}D."
            )

        actual_nodes = self.features.size(1)
        if actual_nodes != gauge_domain.num_nodes:
            raise ValueError(
                f"Geometry mismatch: Domain declares {gauge_domain.num_nodes} nodes, "
                f"but signal features contain {actual_nodes} nodes."
            )

        vector_dim = self.features.size(3)
        if vector_dim != 2:
            raise ValueError(
                f"Tangent vectors must be strictly 2D (x, y components). "
                f"Got final dimension of size {vector_dim}."
            )

    @property
    def domain(self) -> GaugeDomain:
        return self._domain  # type: ignore[return-value]

    def to(self, device: torch.device | str) -> "GaugeSignal":
        return GaugeSignal(
            domain=self.domain.to(device),
            features=self.features.to(device),
        )
