import torch

from gdl.core import BaseSignal
from .domain import GeodesicDomain


class GeodesicSignal(BaseSignal):
    """
    4G Signal: Features living on the surface of a 3D Mesh.

    Tensor Shape strictly enforced: [Batch, Nodes, Channels]
    (Nodes corresponds to the vertices of the GeodesicDomain. This seamlessly
    aligns with the sequence paradigm of 3G Graph Neural Networks).
    """

    def __init__(self, domain: GeodesicDomain, features: torch.Tensor):
        super().__init__(domain, features)
        self._validate_shape()

    def _validate_shape(self) -> None:
        """Ensures the feature tensor mathematically matches the vertex count."""
        geodesic_domain: GeodesicDomain = self.domain

        if self.features.dim() != 3:
            raise ValueError(
                f"GeodesicSignal expects 3D tensor [Batch, Nodes, Channels], "
                f"got {self.features.dim()}D."
            )

        actual_nodes = self.features.size(1)
        if actual_nodes != geodesic_domain.num_nodes:
            raise ValueError(
                f"Geometry mismatch: Domain declares {geodesic_domain.num_nodes} vertices, "
                f"but signal features contain {actual_nodes} nodes."
            )

    @property
    def domain(self) -> GeodesicDomain:
        return self._domain  # type: ignore[return-value]

    def to(self, device: torch.device | str) -> "GeodesicSignal":
        return GeodesicSignal(
            domain=self.domain.to(device),
            features=self.features.to(device),
        )
