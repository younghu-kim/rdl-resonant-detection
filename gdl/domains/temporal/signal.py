import torch

from gdl.core import BaseSignal
from .domain import TimeDomain


class TimeSignal(BaseSignal):
    """
    Time Signal: Features flowing along the TimeDomain.

    Strictly enforces the tensor shape convention: [Batch, SequenceLength, Channels].
    This seamlessly aligns with PyTorch's `batch_first=True` RNN convention.
    """

    def __init__(self, domain: TimeDomain, features: torch.Tensor):
        super().__init__(domain, features)
        self._validate_shape()

    def _validate_shape(self) -> None:
        """Ensures the feature tensor mathematically matches the sequence length."""
        time_domain: TimeDomain = self.domain

        if self.features.dim() != 3:
            raise ValueError(
                f"TimeSignal expects 3D tensor [Batch, SequenceLength, Channels], "
                f"got {self.features.dim()}D."
            )

        actual_seq_len = self.features.size(1)
        if actual_seq_len != time_domain.sequence_length:
            raise ValueError(
                f"Geometry mismatch: Domain expects sequence length "
                f"{time_domain.sequence_length}, "
                f"but signal features have length {actual_seq_len}."
            )

    @property
    def domain(self) -> TimeDomain:
        return self._domain  # type: ignore[return-value]

    def to(self, device: torch.device | str) -> "TimeSignal":
        return TimeSignal(
            domain=self.domain.to(device),
            features=self.features.to(device),
        )
