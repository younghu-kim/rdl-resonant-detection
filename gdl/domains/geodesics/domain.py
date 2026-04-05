import torch

from gdl.core import BaseDomain


class GeodesicDomain(BaseDomain):
    """
    4G Domain: Geodesic Manifolds (e.g., 3D Meshes).

    A 3D mesh is an irregular discretization of a continuous 2D surface.
    Unlike 3G Graphs where edges denote arbitrary relationships, edges here
    represent physical proximity and curvature on a continuous manifold.

    The fundamental geometric invariant of a manifold is the Laplace-Beltrami Operator (LBO).
    In discrete mesh form, it is rigorously approximated using the Cotangent Laplacian matrix,
    which perfectly captures the intrinsic curvature invariant to isometric bending.
    """

    def __init__(
        self,
        vertices: torch.Tensor,
        faces: torch.Tensor,
        laplacian: torch.Tensor | None = None,
    ):
        """
        Args:
            vertices: FloatTensor of shape [num_nodes, 3] representing 3D coordinates.
            faces: LongTensor of shape [num_faces, 3] representing triangle indices.
            laplacian: Optional precomputed Laplacian. If None, it computes dynamically.
        """
        if vertices.dim() != 2 or vertices.size(1) != 3:
            raise ValueError(f"vertices must have shape [num_nodes, 3], got {vertices.shape}")
        if faces.dim() != 2 or faces.size(1) != 3:
            raise ValueError(f"faces must have shape [num_faces, 3], got {faces.shape}")

        self.vertices = vertices
        self.faces = faces
        self.num_nodes = vertices.size(0)
        self.num_faces = faces.size(0)

        # The true geometry is intrinsic, captured by the Laplacian.
        if laplacian is not None:
            self.laplacian = laplacian
        else:
            self.laplacian = self._compute_normalized_cotangent_laplacian()

    def _compute_normalized_cotangent_laplacian(self) -> torch.Tensor:
        """
        Computes the symmetric normalized discrete Laplace-Beltrami Operator
        (Cotangent weights) for a triangle mesh without requiring external heavy C++ libraries.
        Returns a dense [N, N] tensor with eigenvalues bounded in [0, 2] for spectral filtering.
        """
        v = self.vertices
        f = self.faces
        N = self.num_nodes

        # Edge vectors of each triangle
        # Shape: [num_faces, 3]
        v0, v1, v2 = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
        e0 = v2 - v1
        e1 = v0 - v2
        e2 = v1 - v0

        # Compute dot products for angles (Numerator of cotangent: cos(theta))
        dot0 = torch.sum(e1 * (-e2), dim=-1)
        dot1 = torch.sum(e2 * (-e0), dim=-1)
        dot2 = torch.sum(e0 * (-e1), dim=-1)

        # Cross product magnitudes for areas (Denominator of cotangent: sin(theta))
        cross_mag0 = torch.norm(torch.cross(e1, -e2, dim=-1), dim=-1)
        cross_mag1 = torch.norm(torch.cross(e2, -e0, dim=-1), dim=-1)
        cross_mag2 = torch.norm(torch.cross(e0, -e1, dim=-1), dim=-1)

        # Cotangents: cos/sin = dot / cross_mag. Add epsilon to prevent division by zero.
        eps = 1e-6
        cot0 = dot0 / (cross_mag0 + eps)
        cot1 = dot1 / (cross_mag1 + eps)
        cot2 = dot2 / (cross_mag2 + eps)

        # Build the weight matrix W (Adjacency weighted by cotangents)
        W = torch.zeros((N, N), device=v.device, dtype=v.dtype)

        # Accumulate weights symmetrically for each edge in each triangle
        # Edge (1, 2) uses cot0 (angle at vertex 0)
        W.index_put_((f[:, 1], f[:, 2]), cot0 * 0.5, accumulate=True)
        W.index_put_((f[:, 2], f[:, 1]), cot0 * 0.5, accumulate=True)

        # Edge (2, 0) uses cot1 (angle at vertex 1)
        W.index_put_((f[:, 2], f[:, 0]), cot1 * 0.5, accumulate=True)
        W.index_put_((f[:, 0], f[:, 2]), cot1 * 0.5, accumulate=True)

        # Edge (0, 1) uses cot2 (angle at vertex 2)
        W.index_put_((f[:, 0], f[:, 1]), cot2 * 0.5, accumulate=True)
        W.index_put_((f[:, 1], f[:, 0]), cot2 * 0.5, accumulate=True)

        # To ensure stability and bounded spectrum [0, 2], we use the symmetric normalized Laplacian
        # L_sym = I - D^{-1/2} W D^{-1/2}
        D_vec = torch.sum(W, dim=1)
        D_vec_clamped = torch.clamp(D_vec, min=eps)
        D_inv_sqrt = torch.diag(torch.pow(D_vec_clamped, -0.5))

        I = torch.eye(N, device=v.device, dtype=v.dtype)
        L_sym = I - torch.matmul(D_inv_sqrt, torch.matmul(W, D_inv_sqrt))

        return L_sym

    def to(self, device: torch.device | str) -> "GeodesicDomain":
        return GeodesicDomain(
            vertices=self.vertices.to(device),
            faces=self.faces.to(device),
            laplacian=self.laplacian.to(device),
        )

    def __repr__(self) -> str:
        return f"GeodesicDomain(nodes={self.num_nodes}, faces={self.num_faces})"
