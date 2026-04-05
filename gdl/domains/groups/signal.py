import torch

from gdl.core import BaseSignal
from .domain import GroupDomain


class GroupSignal(BaseSignal):
    """
    2G Signal: Features on the GroupDomain represented via Spherical Harmonics coefficients.

    Strictly enforces the tensor shape convention: [Batch, Channels, (l_max + 1)^2]
    The spatial dimensions are replaced by the frequency (spectral) dimension.
    """

    def __init__(self, domain: GroupDomain, features: torch.Tensor):
        super().__init__(domain, features)
        self._validate_shape()

    def _validate_shape(self) -> None:
        """Ensures the feature tensor mathematically matches the spectral frequency limit."""
        group_domain: GroupDomain = self.domain

        if self.features.dim() != 3:
            raise ValueError(
                f"GroupSignal expects 3D tensor [Batch, Channels, Coefficients], "
                f"got {self.features.dim()}D."
            )

        expected_coeffs = group_domain.num_coefficients
        actual_coeffs = self.features.size(2)

        if actual_coeffs != expected_coeffs:
            raise ValueError(
                f"Geometry mismatch: Domain l_max={group_domain.l_max} requires "
                f"{expected_coeffs} coefficients, but signal has {actual_coeffs}."
            )

    @property
    def domain(self) -> GroupDomain:
        return self._domain  # type: ignore[return-value]

    def to(self, device: torch.device | str) -> "GroupSignal":
        return GroupSignal(
            domain=self.domain.to(device),
            features=self.features.to(device),
        )
