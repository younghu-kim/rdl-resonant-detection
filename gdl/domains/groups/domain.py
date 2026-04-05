import torch

from gdl.core import BaseDomain


class GroupDomain(BaseDomain):
    """
    2G Domain: Continuous Groups and Homogeneous Spaces (e.g., Sphere S^2).

    Unlike a flat grid, the sphere has no edges and continuous rotational symmetry SO(3).
    To process data symmetrically without pole distortions, we operate purely in the
    spectral domain using Spherical Harmonics up to a maximum degree (l_max).
    """

    def __init__(self, l_max: int):
        if l_max < 0:
            raise ValueError("l_max must be a non-negative integer.")
        self.l_max = l_max
        # Total number of spherical harmonic coefficients up to degree l_max is (l_max + 1)^2
        self.num_coefficients = (l_max + 1) ** 2

    def to(self, device: torch.device | str) -> "GroupDomain":
        # The abstract spectral domain has no heavy structural tensors to move.
        return self

    def __repr__(self) -> str:
        return f"GroupDomain(l_max={self.l_max}, coeffs={self.num_coefficients})"
