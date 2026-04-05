import torch

from gdl.core import BaseTransform
from .domain import GeodesicDomain
from .signal import GeodesicSignal


class MeshRigidIsometry(BaseTransform[GeodesicDomain, GeodesicSignal]):
    """
    4G Transform: Rigid Isometry (Rotation + Translation) of a 3D Mesh.

    In Intrinsic Geodesic Deep Learning, models must be invariant to Isometric
    Deformations (where intrinsic geodesic distances on the surface are preserved).

    A strict and mathematically computable subset of isometries is Rigid 3D Transformation.
    By applying a random 3D rotation and translation to the extrinsic XYZ coordinates
    of the mesh and forcing the domain to RECOMPUTE its Cotangent Laplacian, we can
    mathematically prove that our intrinsic representation is perfectly invariant
    to how the mesh is posed or placed in 3D space.
    """

    def __init__(self, rotation_matrix: torch.Tensor, translation_vector: torch.Tensor):
        if rotation_matrix.shape != (3, 3):
            raise ValueError("Rotation matrix must be 3x3.")
        if translation_vector.shape != (3,):
            raise ValueError("Translation vector must be 1D with 3 elements.")

        self.R = rotation_matrix
        self.t = translation_vector

    def forward_domain(self, domain: GeodesicDomain) -> GeodesicDomain:
        device = domain.vertices.device
        dtype = domain.vertices.dtype

        R = self.R.to(device=device, dtype=dtype)
        t = self.t.to(device=device, dtype=dtype)

        # Apply Extrinsic 3D Transform: V_new = V @ R^T + t
        new_vertices = torch.matmul(domain.vertices, R.T) + t

        # DYNAMIC RECOMPUTATION:
        # We explicitly pass laplacian=None to trigger _compute_normalized_cotangent_laplacian().
        # Because the transformation preserves all intrinsic edge lengths, angles, and areas,
        # the newly computed Laplacian matrix will be mathematically identical to the original!
        return GeodesicDomain(
            vertices=new_vertices,
            faces=domain.faces,
            laplacian=None,
        )

    def forward_signal(self, signal: GeodesicSignal) -> GeodesicSignal:
        # The intrinsic features (e.g., color, heat, curvature) ride along with the
        # deformed surface without altering their scalar values.
        deformed_domain = self.forward_domain(signal.domain)

        return GeodesicSignal(
            domain=deformed_domain,
            features=signal.features.clone(),
        )
