import torch

from gdl.core import BaseTransform
from .domain import TimeDomain
from .signal import TimeSignal


class TimeWarping(BaseTransform[TimeDomain, TimeSignal]):
    """
    Time Transform: Dynamic Time-Warping (Non-linear 1D Diffeomorphism).

    According to GDL principles (Section 5.8), temporal data (like speech or handwriting)
    must be invariant to the speed at which it unfolds.

    In a discrete framework, this time-warping is simulated by repeating each frame
    in the sequence a variable number of times. For example:
    Input:  [A, B, C]
    Warped: [A, A, A, B, C, C]

    A strictly Time-Warping Invariant model (like a carefully designed Gated RNN)
    must produce the exact same final state regardless of these redundant repetitions.
    """

    def __init__(self, repeats: torch.Tensor):
        """
        Args:
            repeats: LongTensor of shape [sequence_length].
                     Specifies how many times each frame in the sequence should be repeated.
                     Must contain strictly positive integers (>= 1).
        """
        if repeats.dim() != 1:
            raise ValueError(f"repeats must be a 1D tensor, got {repeats.dim()}D.")
        if not torch.all(repeats >= 1):
            raise ValueError(
                "All repeat counts must be at least 1 (no frame dropping allowed)."
            )

        self.repeats = repeats

    def forward_domain(self, domain: TimeDomain) -> TimeDomain:
        if self.repeats.size(0) != domain.sequence_length:
            raise ValueError(
                f"The length of the repeats tensor ({self.repeats.size(0)}) must match "
                f"the domain sequence length ({domain.sequence_length})."
            )

        # The new sequence length is simply the sum of all frame repetitions
        new_sequence_length = int(self.repeats.sum().item())

        # DYNAMIC MUTATION: Time stretches elastically!
        return TimeDomain(sequence_length=new_sequence_length)

    def forward_signal(self, signal: TimeSignal) -> TimeSignal:
        domain = signal.domain
        device = signal.features.device

        # Ensure repeats is on the same device
        repeats_dev = self.repeats.to(device=device)

        # x shape: [Batch, SequenceLength, Channels]
        x = signal.features

        # Apply torch.repeat_interleave along the time dimension (dim=1)
        warped_features = torch.repeat_interleave(x, repeats_dev, dim=1)

        warped_domain = self.forward_domain(domain)

        return TimeSignal(
            domain=warped_domain,
            features=warped_features,
        )
