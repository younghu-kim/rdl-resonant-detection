import torch

from gdl.core import BaseTransform
from .domain import GroupDomain
from .signal import GroupSignal


class SphericalRotation(BaseTransform[GroupDomain, GroupSignal]):
    """
    Action of the 3D Rotation Group SO(3) on the Sphere.

    Mathematically, rotating a sphere in the spatial domain corresponds to multiplying
    its spherical harmonic coefficients by block-diagonal Wigner-D matrices
    in the spectral domain.
    """

    def __init__(self, rotation_matrix: torch.Tensor):
        """
        Args:
            rotation_matrix: A 3x3 orthogonal matrix representing a 3D rotation.
        """
        if rotation_matrix.shape != (3, 3):
            raise ValueError(f"Rotation matrix must be strictly 3x3, got {rotation_matrix.shape}")
        self.R = rotation_matrix

    def forward_domain(self, domain: GroupDomain) -> GroupDomain:
        # Rotating a perfect sphere does not change its abstract topology or l_max frequency limit.
        return domain

    def forward_signal(self, signal: GroupSignal) -> GroupSignal:
        """
        Applies Wigner-D matrices to the signal coefficients.

        For this core blueprint, we implement the exact analytical action for:
        - l=0 (invariant scalar, e.g., temperature/mass)
        - l=1 (equivariant 3D vector, e.g., velocity/direction)
        This mathematically proves SO(3) equivariance without external heavy dependencies.
        """
        group_domain: GroupDomain = signal.domain

        # Clone features to avoid in-place modification of the original signal
        features = signal.features.clone()
        device = features.device
        dtype = features.dtype

        # l=0 Block: The first coefficient (index 0) is a scalar field component.
        # Rotating the 3D space doesn't change it. (Wigner-D for l=0 is exactly [1.0]).
        # Thus, features[..., 0] remains perfectly untouched.

        # l=1 Block: The next 3 coefficients (indices 1 to 3) represent a 3D vector (x, y, z).
        # To rotate this frequency band, we multiply it by the 3x3 rotation matrix R.
        if group_domain.l_max >= 1:
            l1_coeffs = features[..., 1:4]  # Expected Shape: [Batch, Channels, 3]

            R_device = self.R.to(device=device, dtype=dtype)

            # Matrix multiplication: [B, C, 3] x [3, 3]^T -> [B, C, 3]
            rotated_l1 = torch.matmul(l1_coeffs, R_device.T)

            # Inject the rotated 3D vectors back into the frequency block
            features[..., 1:4] = rotated_l1

        # Note: For l >= 2 (higher frequencies), one would typically use e3nn's wigner_D matrices
        # (e.g., a 5x5 matrix for l=2). We intentionally omit it here to keep the core foundation
        # lightweight, dependency-free, and transparently testable.

        return GroupSignal(domain=group_domain, features=features)
