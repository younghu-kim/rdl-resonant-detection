from abc import ABC, abstractmethod
from typing import TypeVar, Generic

from .base_domain import BaseDomain
from .base_signal import BaseSignal

# 도메인과 시그널에 대한 유연한 타입 바인딩
T_Domain = TypeVar('T_Domain', bound=BaseDomain)
T_Signal = TypeVar('T_Signal', bound=BaseSignal)


class BaseTransform(ABC, Generic[T_Domain, T_Signal]):
    """
    Abstract Base Class representing an action of a symmetry group element g ∈ 𝔊.
    Transforms can be applied to either the Domain itself, or the Signal on the domain.
    """

    @abstractmethod
    def forward_domain(self, domain: T_Domain) -> T_Domain:
        """
        Applies the geometric transformation to the Domain.
        E.g., rotating the spatial coordinates of a mesh or applying a random permutation to graph nodes.
        """
        pass

    @abstractmethod
    def forward_signal(self, signal: T_Signal) -> T_Signal:
        """
        Applies the geometric transformation to the Signal.
        Crucially, this often involves applying the inverse transformation to the signal's coordinates
        or representation, following the representation theory rule: [ρ(g)x](u) = x(g⁻¹ u).
        """
        pass

    def __call__(self, obj: T_Signal | T_Domain) -> T_Signal | T_Domain:
        """
        Syntactic sugar to apply the transformation directly.
        Routes to the appropriate method based on the object type.
        """
        if isinstance(obj, BaseSignal):
            return self.forward_signal(obj)  # type: ignore[arg-type]
        elif isinstance(obj, BaseDomain):
            return self.forward_domain(obj)  # type: ignore[arg-type]
        else:
            raise TypeError(
                f"BaseTransform can only be applied to BaseSignal or BaseDomain, got {type(obj)}"
            )

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
