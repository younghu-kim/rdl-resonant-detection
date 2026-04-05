from abc import ABC, abstractmethod
from typing import Generic, TypeVar, Sequence

import torch
import torch.nn as nn

from .base_signal import BaseSignal

# 블루프린트 블록 전체에서 사용하는 시그널 타입 변수
T_Signal = TypeVar("T_Signal", bound=BaseSignal)


# ─────────────────────────────────────────────────────────────
# Block 1: 등변 레이어 (Equivariant Layer)
# ─────────────────────────────────────────────────────────────
class EquivariantLayer(nn.Module, ABC, Generic[T_Signal]):
    """
    Blueprint Block 1: Local Equivariant Layer.

    This layer mixes information among neighbors in a way that respects the
    geometric symmetry of the domain. It guarantees Equivariance.

    Mathematical condition: f(ρ(g)x) = ρ'(g)f(x)
    """

    def __init__(self) -> None:
        super().__init__()

    @abstractmethod
    def forward(self, signal: T_Signal) -> T_Signal:
        """
        Applies a local equivariant transformation to the input signal.

        Args:
            signal (T_Signal): The input signal living on a specific domain.

        Returns:
            T_Signal: The transformed signal. The output domain is typically
                      the same as the input domain.
        """
        pass


# ─────────────────────────────────────────────────────────────
# Block 2: 기하학적 비선형성 (Geometric Non-Linearity)
# ─────────────────────────────────────────────────────────────
class GeometricNonLinearity(nn.Module, ABC, Generic[T_Signal]):
    """
    Blueprint Block 2: Geometric Non-Linearity.

    Applies a point-wise non-linear activation to the signal.
    Crucially, it must be designed carefully so as not to break the geometric
    properties (e.g., orientation) of the features.

    Mathematical condition: σ(ρ(g)x) = ρ'(g)σ(x)

    WARNING: For simple scalar features, this can be standard ReLU.
    However, for vector/tensor features (e.g., 3D velocity), applying ReLU channel-wise
    breaks geometric symmetries. In such cases, gated non-linearities or norm-based
    activations must be implemented in the derived classes.
    """

    def __init__(self) -> None:
        super().__init__()

    @abstractmethod
    def forward(self, signal: T_Signal) -> T_Signal:
        """
        Applies a geometry-preserving non-linear activation to the input signal.

        Args:
            signal (T_Signal): The input signal living on a specific domain.

        Returns:
            T_Signal: The non-linearly transformed signal. The domain remains unchanged.
        """
        pass


# ─────────────────────────────────────────────────────────────
# Block 3: 로컬 풀링 (Local Pooling / Coarse-graining)
# ─────────────────────────────────────────────────────────────
class LocalPooling(nn.Module, ABC, Generic[T_Signal]):
    """
    Blueprint Block 3: Local Pooling (Coarse-graining).

    This layer reduces the spatial resolution of the domain (coarse-graining)
    while aggregating the signal features within local neighborhoods.
    This creates 'Scale Separation', expanding the receptive field and providing
    stability against local deformations.

    Mathematical condition: P : 𝒳(Ω) -> 𝒳(Ω') where Ω' ⊆ Ω (or structurally coarser).
    """

    def __init__(self) -> None:
        super().__init__()

    @abstractmethod
    def forward(self, signal: T_Signal) -> T_Signal:
        """
        Applies a local pooling operation to coarsen the domain and aggregate signals.

        Args:
            signal (T_Signal): The input signal living on the fine domain Ω.

        Returns:
            T_Signal: A NEW signal object living on a NEW, coarser domain Ω'.
                      Both the features and the underlying domain structure MUST be modified.
        """
        pass


# ─────────────────────────────────────────────────────────────
# Block 4: 전역 불변 레이어 (Global Invariant / Readout)
# ─────────────────────────────────────────────────────────────
class GlobalInvariant(nn.Module, ABC, Generic[T_Signal]):
    """
    Blueprint Block 4: Global Invariant Layer (Readout).

    This layer destroys the spatial/geometric structure of the domain entirely,
    aggregating the signal into a fixed-size vector (or set of vectors).
    This guarantees Invariance to the symmetry group of the domain.

    Mathematical condition: f(ρ(g)x) = f(x)
    """

    def __init__(self) -> None:
        super().__init__()

    @abstractmethod
    def forward(self, signal: T_Signal) -> torch.Tensor:
        """
        Applies a global invariant aggregation to the input signal.

        Args:
            signal (T_Signal): The input signal living on the domain.

        Returns:
            torch.Tensor: The aggregated invariant feature tensor.
                          Notice that it returns a raw torch.Tensor, NOT a BaseSignal,
                          because the geometric domain information has been stripped away.
                          Typical shape: [Batch_size, Hidden_channels]
        """
        pass


# ─────────────────────────────────────────────────────────────
# GDL 파이프라인 모델 (Blueprint Pipeline Template)
# ─────────────────────────────────────────────────────────────
class BaseGDLModel(nn.Module, ABC, Generic[T_Signal]):
    """
    The Geometric Deep Learning Blueprint Pipeline.

    Combines the four fundamental blocks to create a G-invariant function f: 𝒳(Ω) -> 𝒴
    as defined in Bronstein et al. (2021):
    f = A ∘ σ_J ∘ B_J ∘ P_{J-1} ∘ ... ∘ P_1 ∘ σ_1 ∘ B_1

    Where:
    - B: EquivariantLayer
    - σ: GeometricNonLinearity
    - P: LocalPooling
    - A: GlobalInvariant
    """

    def __init__(
        self,
        blocks: Sequence[EquivariantLayer | GeometricNonLinearity | LocalPooling],
        readout: GlobalInvariant,
    ):
        """
        Args:
            blocks: A sequence of B, σ, and P layers that maintain Signal context.
            readout: The final A layer that destroys the domain and outputs a Tensor.
        """
        super().__init__()

        # 모든 중간 블록이 GDL ABC 인스턴스인지 검증
        for idx, block in enumerate(blocks):
            if not isinstance(block, (EquivariantLayer, GeometricNonLinearity, LocalPooling)):
                raise TypeError(
                    f"Block at index {idx} must be an EquivariantLayer, GeometricNonLinearity, "
                    f"or LocalPooling. Got {type(block)}"
                )

        if not isinstance(readout, GlobalInvariant):
            raise TypeError(f"Readout must be a GlobalInvariant layer. Got {type(readout)}")

        self.blocks = nn.ModuleList(blocks)  # type: ignore[arg-type]
        self.readout = readout

    def forward(self, signal: T_Signal) -> torch.Tensor:
        """
        Executes the GDL blueprint.

        Args:
            signal (T_Signal): The initial input signal.

        Returns:
            torch.Tensor: The final invariant feature representation.
        """
        current_signal = signal

        # 1. B, σ, P 시퀀스를 통과 (도메인 기하학이 보존/조립화됨)
        for block in self.blocks:
            current_signal = block(current_signal)

        # 2. 최종 전역 리드아웃 A (도메인 기하학이 완전히 제거됨)
        out_tensor = self.readout(current_signal)

        return out_tensor
