from abc import ABC, abstractmethod
import torch

from .base_domain import BaseDomain


class BaseSignal(ABC):
    """
    Abstract Base Class representing a signal X(Ω, C) living on a domain Ω.
    This strictly binds a feature tensor to a specific Geometric Domain.
    """

    def __init__(self, domain: BaseDomain, features: torch.Tensor):
        self._domain = domain
        self._features = features

    @property
    def domain(self) -> BaseDomain:
        """Returns the underlying geometric domain."""
        return self._domain

    @property
    def features(self) -> torch.Tensor:
        """Returns the raw feature tensor."""
        return self._features

    @features.setter
    def features(self, value: torch.Tensor) -> None:
        self._features = value

    @property
    def shape(self) -> torch.Size:
        """Returns the shape of the signal."""
        return self._features.shape

    @property
    def device(self) -> torch.device:
        """Returns the device of the feature tensor."""
        return self._features.device

    @abstractmethod
    def to(self, device: torch.device | str) -> 'BaseSignal':
        """Moves the signal and its domain to the specified device."""
        pass

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(shape={self.shape}, domain={self.domain.__class__.__name__})"
