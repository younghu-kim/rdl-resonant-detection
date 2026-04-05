import torch
import torch.nn as nn

from gdl.core import (
    EquivariantLayer,
    GeometricNonLinearity,
    LocalPooling,
    GlobalInvariant,
)
from .domain import GroupDomain
from .signal import GroupSignal


class GroupConvolution(EquivariantLayer[GroupSignal]):
    """
    2G Equivariant Layer: Spectral Convolution on the Continuous Group (SO(3)).

    Guarantees SO(3) Rotation Equivariance by applying learnable weights
    strictly within each frequency degree (l) of the Spherical Harmonics.

    Mathematical GDL Formulation (Schur's Lemma):
    To commute with the Wigner-D matrices (W * D^l = D^l * W),
    a convolution in the spectral domain must apply the EXACT SAME linear
    transformation across the (2l+1) components of degree l.
    """

    def __init__(self, in_channels: int, out_channels: int, l_max: int, bias: bool = True):
        super().__init__()
        self.l_max = l_max
        self.in_channels = in_channels
        self.out_channels = out_channels

        # We need a separate learnable weight matrix for each frequency degree l.
        # Shape of each weight: [out_channels, in_channels]
        self.weights_l = nn.ParameterList([
            nn.Parameter(torch.randn(out_channels, in_channels) / (in_channels ** 0.5))
            for _ in range(l_max + 1)
        ])

        # CRUCIAL GEOMETRIC CONSTRAINT:
        # Bias can ONLY be applied to the rotation-invariant scalar component (l=0).
        # Adding a fixed non-zero vector bias to a 3D vector (l>=1) would immediately
        # shift the origin and break the rotational symmetry of the space!
        if bias:
            self.bias_l0 = nn.Parameter(torch.zeros(out_channels))
        else:
            self.register_parameter("bias_l0", None)

    def forward(self, signal: GroupSignal) -> GroupSignal:
        domain: GroupDomain = signal.domain

        if domain.l_max != self.l_max:
            raise ValueError(
                f"Dimension mismatch: Layer expects l_max={self.l_max}, "
                f"but signal domain has l_max={domain.l_max}."
            )

        # x shape: [Batch, C_in, (l_max+1)^2]
        x = signal.features

        out_blocks: list[torch.Tensor] = []

        # Apply convolution independently to each frequency degree l
        for l in range(self.l_max + 1):
            # The mathematical components for degree l are stored from index l^2 to (l+1)^2 - 1
            start_idx = l ** 2
            end_idx = (l + 1) ** 2

            # Extract the (2l+1) spatial components for this degree
            # Shape: [Batch, C_in, 2l+1]
            x_l = x[..., start_idx:end_idx]

            # Extract the learnable weight for degree l
            # Shape: [C_out, C_in]
            W_l = self.weights_l[l]

            # 1. Linear combination over the channel dimension.
            # We strictly DO NOT mix the spatial (m) components to preserve SO(3) symmetry.
            # o: out_channels, c: in_channels, b: batch, m: spatial components (2l+1)
            out_l = torch.einsum("oc, bcm -> bom", W_l, x_l)

            # 2. Add bias ONLY to the invariant scalar component (l=0)
            if l == 0 and self.bias_l0 is not None:
                # bias_l0 shape: [C_out] -> [1, C_out, 1] for broadcasting
                out_l = out_l + self.bias_l0.view(1, -1, 1)

            out_blocks.append(out_l)

        # Reconstruct the full spectral tensor by concatenating along the frequency dimension
        # Final shape: [Batch, C_out, (l_max+1)^2]
        out_features = torch.cat(out_blocks, dim=-1)

        return GroupSignal(domain=domain, features=out_features)


class GroupNonLinearity(GeometricNonLinearity[GroupSignal]):
    """
    2G Non-Linearity: SO(3)-Equivariant Activation (Norm-based Gating).

    Standard pointwise activations (like ReLU) break rotation equivariance if applied
    directly to vector/tensor components (l > 0).
    Solution: Apply ReLU only to the invariant scalar (l=0) or the L2-Norm of the vectors,
    preserving the original spatial direction perfectly intact.
    """

    def __init__(self) -> None:
        super().__init__()
        self.relu = nn.ReLU()

    def forward(self, signal: GroupSignal) -> GroupSignal:
        domain: GroupDomain = signal.domain
        x = signal.features

        out_blocks: list[torch.Tensor] = []

        for l in range(domain.l_max + 1):
            start_idx = l ** 2
            end_idx = (l + 1) ** 2

            # x_l shape: [Batch, Channels, 2l+1]
            x_l = x[..., start_idx:end_idx]

            if l == 0:
                # l=0 is an invariant scalar, safe to use standard ReLU
                out_l = self.relu(x_l)
            else:
                # l > 0 are equivariant vectors/tensors.
                # Calculate the L2 Norm along the component dimension (dim=-1)
                norm = torch.norm(x_l, p=2, dim=-1, keepdim=True)  # [Batch, Channels, 1]

                # Avoid division by zero
                safe_norm = norm + 1e-8

                # Apply ReLU to the invariant norm to create a gating scalar
                gate = self.relu(norm)

                # Scale the original vector (preserves its 3D direction perfectly)
                out_l = x_l * (gate / safe_norm)

            out_blocks.append(out_l)

        out_features = torch.cat(out_blocks, dim=-1)

        return GroupSignal(domain=domain, features=out_features)


class GroupPooling(LocalPooling[GroupSignal]):
    """
    2G Local Pooling: Spectral Truncation (Low-pass filtering).

    In the spherical spectral domain, reducing the spatial resolution is mathematically
    equivalent to truncating the higher frequency Spherical Harmonics coefficients.
    We drop any degree l > target_l_max.

    The geometry dynamically mutates: The sphere's frequency limit (l_max) is reduced.
    """

    def __init__(self, target_l_max: int):
        super().__init__()
        if target_l_max < 0:
            raise ValueError("target_l_max must be non-negative.")
        self.target_l_max = target_l_max

    def forward(self, signal: GroupSignal) -> GroupSignal:
        domain: GroupDomain = signal.domain

        if self.target_l_max >= domain.l_max:
            raise ValueError(
                f"target_l_max ({self.target_l_max}) must be strictly less than "
                f"current domain l_max ({domain.l_max}) to perform pooling."
            )

        # The number of coefficients to keep is (target_l_max + 1)^2
        target_coeffs = (self.target_l_max + 1) ** 2

        # Slicing the tensor removes the high-frequency bands completely
        # Shape: [Batch, Channels, target_coeffs]
        out_features = signal.features[..., :target_coeffs]

        # DYNAMIC DOMAIN MUTATION: The geometry's resolution has shrunk!
        coarsened_domain = GroupDomain(l_max=self.target_l_max)

        return GroupSignal(domain=coarsened_domain, features=out_features)


class GroupGlobalPooling(GlobalInvariant[GroupSignal]):
    """
    2G Global Invariant: Destroys the spherical geometry to extract SO(3)-Invariant features.

    Guarantees absolute Rotation Invariance: f(R * X) = f(X).
    We can extract invariants in two ways:
    1. Spherical Integral (l=0): The DC component is intrinsically invariant.
    2. Power Spectrum: The length (L2 norm) of the vector at any degree l is invariant
       under Wigner-D orthogonal rotations.
    """

    def __init__(self, use_power_spectrum: bool = True):
        super().__init__()
        self.use_power_spectrum = use_power_spectrum

    def forward(self, signal: GroupSignal) -> torch.Tensor:
        domain: GroupDomain = signal.domain
        x = signal.features

        if not self.use_power_spectrum:
            # Simplest invariant: Extract ONLY the l=0 scalar component (Spherical Average)
            # Output Shape: [Batch, Channels]
            invariant_tensor = x[..., 0]
        else:
            # Power Spectrum: Extract the invariant scalar for l=0,
            # and the L2 norm for EACH frequency degree l > 0.
            invariant_blocks: list[torch.Tensor] = []

            for l in range(domain.l_max + 1):
                start_idx = l ** 2
                end_idx = (l + 1) ** 2

                # Shape: [Batch, Channels, 2l+1]
                x_l = x[..., start_idx:end_idx]

                if l == 0:
                    # l=0 is already an invariant scalar
                    # Shape: [Batch, Channels]
                    invariant_blocks.append(x_l.squeeze(-1))
                else:
                    # The L2 norm of x_l is invariant because Wigner-D matrices are orthogonal!
                    # Shape: [Batch, Channels]
                    norm_l = torch.norm(x_l, p=2, dim=-1)
                    invariant_blocks.append(norm_l)

            # Stack along the last dimension
            # Shape: [Batch, Channels, l_max + 1]
            stacked_invariants = torch.stack(invariant_blocks, dim=-1)

            # Flatten the channels and spectral norm dimensions into a single vector per batch
            # Final Shape: [Batch, Channels * (l_max + 1)]
            batch_size = x.size(0)
            invariant_tensor = stacked_invariants.view(batch_size, -1)

        return invariant_tensor
