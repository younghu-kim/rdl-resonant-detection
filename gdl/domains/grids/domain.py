from typing import Tuple

import torch

from gdl.core.base_domain import BaseDomain


class GridDomain(BaseDomain):
    """
    1G Domain: Regular Euclidean Grid (1D, 2D, or 3D).

    Unlike Graphs or Meshes which require explicit adjacency matrices or coordinates,
    a regular grid is implicit. The geometry is fully defined by its spatial dimensions.

    Examples:
        - 1D Audio: spatial_shape=(Length,)
        - 2D Image: spatial_shape=(Height, Width)
        - 3D Voxel: spatial_shape=(Depth, Height, Width)
    """

    def __init__(self, spatial_shape: Tuple[int, ...]) -> None:
        if not spatial_shape:
            raise ValueError("spatial_shape cannot be empty.")
        self.spatial_shape = spatial_shape
        self.num_spatial_dims = len(spatial_shape)

    def to(self, device: torch.device | str) -> "GridDomain":
        # Grids are implicit domains; there are no heavy structural tensors
        # to move to the GPU.
        return self

    def __repr__(self) -> str:
        return f"GridDomain({self.num_spatial_dims}D, shape={self.spatial_shape})"
