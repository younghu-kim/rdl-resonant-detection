import torch
import torch.nn as nn

from gdl.core.base_blueprint import (
    EquivariantLayer,
    GeometricNonLinearity,
    LocalPooling,
    GlobalInvariant,
)
from .domain import GridDomain
from .signal import GridSignal


class GridConvolution(EquivariantLayer[GridSignal]):
    """
    1G Equivariant Layer: Convolution on a Grid.

    Ensures Translation Equivariance: f(T_shift(x)) = T_shift(f(x))
    Uses 'circular' padding to perfectly match the periodic boundary
    condition (torch.roll) used in GridTranslation, preventing edge artifacts
    from breaking mathematical symmetries.
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        kernel_size: int,
        num_spatial_dims: int,
        bias: bool = True,
    ) -> None:
        super().__init__()
        self.num_spatial_dims = num_spatial_dims

        # Kernel must be odd to maintain symmetric center
        if kernel_size % 2 == 0:
            raise ValueError("kernel_size must be odd to maintain symmetric padding.")

        # To maintain the same spatial dimensions (stride=1)
        padding = kernel_size // 2

        self.conv: nn.Module
        # Dynamically select the correct PyTorch Conv module
        if num_spatial_dims == 1:
            self.conv = nn.Conv1d(
                in_channels, out_channels, kernel_size,
                padding=padding, padding_mode="circular", bias=bias,
            )
        elif num_spatial_dims == 2:
            self.conv = nn.Conv2d(
                in_channels, out_channels, kernel_size,
                padding=padding, padding_mode="circular", bias=bias,
            )
        elif num_spatial_dims == 3:
            self.conv = nn.Conv3d(
                in_channels, out_channels, kernel_size,
                padding=padding, padding_mode="circular", bias=bias,
            )
        else:
            raise ValueError(
                f"GridConvolution supports 1D, 2D, or 3D grids. Got {num_spatial_dims}D."
            )

    def forward(self, signal: GridSignal) -> GridSignal:
        """Applies translation-equivariant convolution to the GridSignal."""
        grid_domain: GridDomain = signal.domain

        if grid_domain.num_spatial_dims != self.num_spatial_dims:
            raise ValueError(
                f"Layer expects {self.num_spatial_dims}D grid, "
                f"but received signal on {grid_domain.num_spatial_dims}D grid."
            )

        # Apply the convolution over the features
        out_features = self.conv(signal.features)

        # Return a new GridSignal on the SAME domain (spatial shape is preserved)
        return GridSignal(domain=grid_domain, features=out_features)


class GridNonLinearity(GeometricNonLinearity[GridSignal]):
    """
    1G Non-Linearity: Point-wise activation on a Grid.

    Since 1G grid features (like RGB or Audio) are typically scalar fields,
    a standard point-wise ReLU perfectly preserves Translation Equivariance.
    """

    def __init__(self, inplace: bool = False) -> None:
        super().__init__()
        self.activation = nn.ReLU(inplace=inplace)

    def forward(self, signal: GridSignal) -> GridSignal:
        out_features = self.activation(signal.features)
        # Domain remains unchanged
        return GridSignal(domain=signal.domain, features=out_features)


class GridPooling(LocalPooling[GridSignal]):
    """
    1G Local Pooling: Coarse-graining the Grid.

    Reduces the spatial resolution, creating 'Scale Separation'.
    Crucially, it dynamically generates a NEW GridDomain representing
    the coarsened geometry to strictly follow the GDL Blueprint.
    """

    def __init__(self, kernel_size: int, stride: int, num_spatial_dims: int) -> None:
        super().__init__()
        self.num_spatial_dims = num_spatial_dims

        self.pool: nn.Module
        if num_spatial_dims == 1:
            self.pool = nn.AvgPool1d(kernel_size, stride=stride)
        elif num_spatial_dims == 2:
            self.pool = nn.AvgPool2d(kernel_size, stride=stride)
        elif num_spatial_dims == 3:
            self.pool = nn.AvgPool3d(kernel_size, stride=stride)
        else:
            raise ValueError(f"Unsupported spatial dims: {num_spatial_dims}")

    def forward(self, signal: GridSignal) -> GridSignal:
        grid_domain: GridDomain = signal.domain

        if grid_domain.num_spatial_dims != self.num_spatial_dims:
            raise ValueError("Dimension mismatch between signal domain and pooling layer.")

        # Pool features
        out_features = self.pool(signal.features)

        # DYNAMIC DOMAIN MUTATION: The geometry has shrunk!
        new_spatial_shape = tuple(out_features.shape[2:])
        coarsened_domain = GridDomain(spatial_shape=new_spatial_shape)

        return GridSignal(domain=coarsened_domain, features=out_features)


class GridGlobalAveragePooling(GlobalInvariant[GridSignal]):
    """
    1G Global Invariant: Destroys the grid geometry completely.

    Averages features across all spatial dimensions to produce a fixed-size vector.
    Ensures absolute Translation Invariance: f(T(x)) = f(x).
    """

    def __init__(self) -> None:
        super().__init__()

    def forward(self, signal: GridSignal) -> torch.Tensor:
        grid_domain: GridDomain = signal.domain

        # The spatial dims are at the end. E.g., for 2D Grid: dims = (-2, -1)
        dims_to_reduce = tuple(range(-grid_domain.num_spatial_dims, 0))

        # Average across all spatial dimensions
        invariant_tensor = torch.mean(signal.features, dim=dims_to_reduce)

        # The domain is now destroyed. Return a raw tensor.
        return invariant_tensor
