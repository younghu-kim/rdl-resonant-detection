import torch

from gdl.core.base_signal import BaseSignal
from .domain import GridDomain


class GridSignal(BaseSignal):
    """
    1G Signal: Data living on a regular GridDomain.

    Strictly enforces the tensor shape convention: [Batch, Channels, *Spatial_Dims]
    """

    def __init__(self, domain: GridDomain, features: torch.Tensor) -> None:
        super().__init__(domain, features)
        self._validate_shape()

    def _validate_shape(self) -> None:
        """Ensures the feature tensor mathematically matches the grid geometry."""
        grid_domain: GridDomain = self.domain
        expected_dims = grid_domain.num_spatial_dims + 2  # +2 for [Batch, Channels]

        if self.features.dim() != expected_dims:
            raise ValueError(
                f"Shape mismatch: GridDomain is {grid_domain.num_spatial_dims}D, "
                f"so features must be {expected_dims}D [B, C, *Spatial]. "
                f"Got {self.features.dim()}D tensor."
            )

        spatial_shape_of_features = tuple(self.features.shape[2:])
        if spatial_shape_of_features != grid_domain.spatial_shape:
            raise ValueError(
                f"Geometry mismatch: Feature spatial shape {spatial_shape_of_features} "
                f"does not match GridDomain shape {grid_domain.spatial_shape}."
            )

    @property
    def domain(self) -> GridDomain:
        return self._domain  # type: ignore[return-value]

    def to(self, device: torch.device | str) -> "GridSignal":
        return GridSignal(
            domain=self.domain.to(device),
            features=self.features.to(device),
        )
