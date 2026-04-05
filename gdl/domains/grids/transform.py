from typing import Tuple

import torch

from gdl.core.base_transform import BaseTransform
from .domain import GridDomain
from .signal import GridSignal


class GridTranslation(BaseTransform[GridDomain, GridSignal]):
    """
    Action of the Translation Symmetry Group on the Grid.
    Represents shifting the grid spatially.
    """

    def __init__(self, shifts: Tuple[int, ...]) -> None:
        """
        Args:
            shifts: Number of steps to shift along each spatial dimension.
                    e.g., (2, -1) means shift down by 2, left by 1 in a 2D grid.
        """
        self.shifts = shifts

    def forward_domain(self, domain: GridDomain) -> GridDomain:
        # Translating a boundless Euclidean grid (or assuming periodic boundaries)
        # does not change its topology or shape.
        return domain

    def forward_signal(self, signal: GridSignal) -> GridSignal:
        """
        Shifts the features along the spatial dimensions.
        Simulates rho(g)x(u) = x(g^{-1} u) using periodic boundary conditions (roll).
        """
        grid_domain: GridDomain = signal.domain

        if len(self.shifts) != grid_domain.num_spatial_dims:
            raise ValueError(
                f"Number of shifts ({len(self.shifts)}) must match "
                f"domain dimensions ({grid_domain.num_spatial_dims})."
            )

        # Spatial dimensions are at the end: [-num_spatial_dims, ..., -1]
        dims = tuple(range(-grid_domain.num_spatial_dims, 0))
        shifted_features = torch.roll(signal.features, shifts=self.shifts, dims=dims)

        return GridSignal(domain=grid_domain, features=shifted_features)
