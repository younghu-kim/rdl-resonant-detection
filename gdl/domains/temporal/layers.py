import torch
import torch.nn as nn

from gdl.core import (
    EquivariantLayer,
    GeometricNonLinearity,
    GlobalInvariant,
)
from .domain import TimeDomain
from .signal import TimeSignal


class TimeDifferenceGatedRNN(EquivariantLayer[TimeSignal]):
    """
    Time Domain Layer: Difference-Gated RNN for Time-Warping Equivariance.

    According to GDL Section 5.8, to achieve strict invariance/equivariance to
    Time-Warping (frame repetitions like [A, A, B]), the network MUST possess a
    Gating mechanism that halts state updates when no new information arrives.

    To mathematically guarantee this blueprint without prior training (for our TDD suite),
    we implement an explicit Difference-Gate. The update gate activates proportionally
    to the change between the current frame and the previous frame (x_t - x_{t-1}).
    When a frame is exactly repeated (Time-Warping), the difference is exactly zero,
    the gate closes entirely, and the hidden state is perfectly copied: h_t = h_{t-1}.
    """

    def __init__(self, in_channels: int, out_channels: int):
        super().__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels

        # Candidate hidden state generator
        self.W_c = nn.Linear(in_channels + out_channels, out_channels, bias=True)

        # Update gate logic based strictly on the feature difference.
        # CRUCIAL GEOMETRIC CONSTRAINT: No bias is used so that diff=0 strictly implies u_t=0.
        self.W_u = nn.Linear(in_channels, out_channels, bias=False)

    def forward(self, signal: TimeSignal) -> TimeSignal:
        domain: TimeDomain = signal.domain

        # x shape: [Batch, SequenceLength, InChannels]
        x = signal.features
        batch_size, seq_len, _ = x.shape
        device = x.device
        dtype = x.dtype

        h_t = torch.zeros(batch_size, self.out_channels, device=device, dtype=dtype)
        h_states: list[torch.Tensor] = []

        # Keep track of previous input to compute the causal difference
        x_prev = torch.zeros(batch_size, self.in_channels, device=device, dtype=dtype)

        for t in range(seq_len):
            x_t = x[:, t, :]

            # 1. Temporal Difference
            if t == 0:
                # First step: the entire input is "new" information
                diff = x_t
            else:
                diff = x_t - x_prev

            x_prev = x_t

            # 2. Candidate State
            combined = torch.cat([x_t, h_t], dim=-1)
            c_t = torch.tanh(self.W_c(combined))

            # 3. Difference Gate: Activated only if there is a change in the input.
            # Using torch.tanh(abs) ensures it acts as a gating scalar in [0, 1) per channel.
            u_t = torch.tanh(torch.abs(self.W_u(diff)))

            # 4. Gated Update
            # Time-Warping Defense: When u_t is 0, h_t is identical to h_{t-1}
            h_t = u_t * c_t + (1.0 - u_t) * h_t

            h_states.append(h_t)

        # Reconstruct the temporal sequence of hidden states.
        # Shape: [Batch, SequenceLength, OutChannels]
        out_features = torch.stack(h_states, dim=1)

        return TimeSignal(domain=domain, features=out_features)


class TimeNonLinearity(GeometricNonLinearity[TimeSignal]):
    """
    Time Domain Non-Linearity: Point-wise Activation.

    Applying a standard ReLU to each frame independently perfectly commutes
    with Time-Warping (frame repetition), maintaining 1D temporal equivariance.
    """

    def __init__(self, inplace: bool = False):
        super().__init__()
        self.activation = nn.ReLU(inplace=inplace)

    def forward(self, signal: TimeSignal) -> TimeSignal:
        out_features = self.activation(signal.features)
        return TimeSignal(domain=signal.domain, features=out_features)


class TimeGlobalPooling(GlobalInvariant[TimeSignal]):
    """
    Time Global Invariant: Extracts the Final Hidden State.

    Due to the causal (directional) nature of time, the entire history of the sequence
    is summarized in the final hidden state of the RNN.
    Extracting only the final frame (t = -1) destroys the time dimension completely
    and guarantees strict invariance to the sequence length and time-warping.
    """

    def __init__(self, aggregation: str = "last"):
        super().__init__()
        valid_aggs = ["last", "max", "mean"]
        if aggregation not in valid_aggs:
            raise ValueError(f"Unsupported aggregation: {aggregation}")
        self.aggregation = aggregation

    def forward(self, signal: TimeSignal) -> torch.Tensor:
        # x shape: [Batch, SequenceLength, Channels]
        x = signal.features

        if self.aggregation == "last":
            # Extract strictly the final time step
            invariant_tensor = x[:, -1, :]
        elif self.aggregation == "max":
            # Max pooling over time is intrinsically invariant to frame repetition
            invariant_tensor = torch.max(x, dim=1)[0]
        elif self.aggregation == "mean":
            invariant_tensor = torch.mean(x, dim=1)
        else:
            raise RuntimeError("Unreachable aggregation state.")

        return invariant_tensor
