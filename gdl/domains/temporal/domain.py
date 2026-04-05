import torch

from gdl.core import BaseDomain


class TimeDomain(BaseDomain):
    """
    Time Domain: 1D Temporal Sequence Geometry.

    Unlike spatial domains (1G~5G), Time has a strict unidirectional causality
    (past -> future).

    According to GDL principles (Sections 5.7-5.8), the primary symmetry of time
    is Time-Warping Invariance (1D Diffeomorphism). The rate at which events unfold
    can stretch or compress dynamically. In our discrete representation, this
    corresponds to changes in the sequence_length via frame repetition or dropping.
    """

    def __init__(self, sequence_length: int):
        if sequence_length <= 0:
            raise ValueError("sequence_length must be strictly positive.")
        self.sequence_length = sequence_length

    def to(self, device: torch.device | str) -> "TimeDomain":
        # The abstract temporal domain has no structural tensors to move.
        return self

    def __repr__(self) -> str:
        return f"TimeDomain(seq_len={self.sequence_length})"
