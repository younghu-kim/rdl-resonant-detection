from abc import ABC, abstractmethod
import torch


class BaseDomain(ABC):
    """
    Abstract Base Class representing a geometric domain Ω.
    Examples: 2D Grid, 3D Mesh, Graph, Spherical Manifold, etc.
    """

    @abstractmethod
    def to(self, device: torch.device | str) -> 'BaseDomain':
        """Moves the domain data (if any structural tensors exist) to the specified device."""
        pass
