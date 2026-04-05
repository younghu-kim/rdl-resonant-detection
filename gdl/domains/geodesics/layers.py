import torch
import torch.nn as nn

from gdl.core import (
    EquivariantLayer,
    GeometricNonLinearity,
    LocalPooling,
    GlobalInvariant,
)
from .domain import GeodesicDomain
from .signal import GeodesicSignal


class IntrinsicMeshConvolution(EquivariantLayer[GeodesicSignal]):
    """
    4G Equivariant Layer: Intrinsic Convolution on a Geodesic Manifold (Mesh).

    Guarantees Isometry Equivariance (invariant to bending and rigid transformations)
    by operating strictly on the spectral domain of the Cotangent Laplacian.

    Mathematical GDL Formulation (ChebNet):
    Instead of costly O(N^3) eigendecomposition, we define spectral filters as
    polynomials of the Laplacian matrix. We use Chebyshev polynomials of the
    first kind T_k(L_tilde) for numerical stability and strict K-hop localization.
    """

    def __init__(self, in_channels: int, out_channels: int, K: int = 3, bias: bool = True):
        super().__init__()
        if K < 1:
            raise ValueError("Polynomial degree K must be at least 1.")

        self.in_channels = in_channels
        self.out_channels = out_channels
        self.K = K

        # Chebyshev Filter Coefficients
        # Shape: [K, in_channels, out_channels]
        self.weights = nn.Parameter(
            torch.randn(K, in_channels, out_channels) / (in_channels ** 0.5)
        )

        if bias:
            self.bias = nn.Parameter(torch.zeros(out_channels))
        else:
            self.register_parameter("bias", None)

    def forward(self, signal: GeodesicSignal) -> GeodesicSignal:
        domain: GeodesicDomain = signal.domain

        # x shape: [Batch, Nodes, Channels]
        x = signal.features

        # The Cotangent Laplacian from 8.1 is symmetrically normalized in [0, 2].
        # We must shift its eigenvalues to [-1, 1] for Chebyshev polynomial stability.
        L_sym = domain.laplacian
        N = L_sym.size(0)
        I = torch.eye(N, device=L_sym.device, dtype=L_sym.dtype)
        L_tilde = L_sym - I

        # T_0(L_tilde) x = x
        x0 = x

        # Filter output accumulation: x0 @ W_0
        out_features = torch.matmul(x0, self.weights[0])

        if self.K > 1:
            # T_1(L_tilde) x = L_tilde @ x
            # L_tilde is [N, N], x0 is [Batch, N, Channels]
            x1 = torch.einsum("nm, bmc -> bnc", L_tilde, x0)

            # Accumulate: x1 @ W_1
            out_features = out_features + torch.matmul(x1, self.weights[1])

            # Recurrence relation: T_k(x) = 2 * L_tilde @ T_{k-1}(x) - T_{k-2}(x)
            for k in range(2, self.K):
                x2 = 2 * torch.einsum("nm, bmc -> bnc", L_tilde, x1) - x0

                # Accumulate: x2 @ W_k
                out_features = out_features + torch.matmul(x2, self.weights[k])

                # Shift variables for the next degree
                x0, x1 = x1, x2

        if self.bias is not None:
            # Broadcast bias over [Batch, Nodes, OutChannels]
            out_features = out_features + self.bias.view(1, 1, -1)

        return GeodesicSignal(domain=domain, features=out_features)


class GeodesicNonLinearity(GeometricNonLinearity[GeodesicSignal]):
    """
    4G Non-Linearity: Point-wise activation on the Manifold.

    Since features on a geodesic manifold (like heat or color at a vertex)
    are scalar fields, a simple pointwise ReLU perfectly preserves Isometry Equivariance.
    """

    def __init__(self, inplace: bool = False):
        super().__init__()
        self.activation = nn.ReLU(inplace=inplace)

    def forward(self, signal: GeodesicSignal) -> GeodesicSignal:
        out_features = self.activation(signal.features)
        return GeodesicSignal(domain=signal.domain, features=out_features)


class MeshDecimationPooling(LocalPooling[GeodesicSignal]):
    """
    4G Local Pooling: Mesh Decimation (Coarsening).

    Reduces the spatial resolution of the manifold by clustering N vertices into
    K super-vertices (target_nodes) using a learned soft-assignment matrix.

    Crucially, to maintain the intrinsic blueprint, the physical 3D coordinates (vertices)
    are also pooled, and a NEW GeodesicDomain is dynamically constructed.
    For simplicity in this blueprint, the face topology (edges) of the coarsened mesh
    is abstracted away (set to empty), treating the pooled geometry as a continuous
    manifold without explicit triangle faces, but with a dynamically recomputed and
    stabilized Laplacian based on the soft-assignment of the original Laplacian.
    """

    def __init__(self, in_channels: int, target_nodes: int):
        super().__init__()
        if target_nodes <= 0:
            raise ValueError("target_nodes must be strictly positive.")

        self.target_nodes = target_nodes
        # Scorer calculates the probability of each vertex merging into one of the K clusters
        self.scorer = nn.Linear(in_channels, target_nodes)

    def forward(self, signal: GeodesicSignal) -> GeodesicSignal:
        domain: GeodesicDomain = signal.domain

        # x shape: [Batch, N, C]
        x = signal.features

        # 1. Compute Soft-Assignment Matrix S. Shape: [Batch, N, K]
        # S_b,i,k is the probability that vertex i belongs to super-vertex k
        S = torch.softmax(self.scorer(x), dim=-1)

        # 2. Coarsen Features: X_new = S^T @ X
        # Shape: [Batch, K, N] @ [Batch, N, C] -> [Batch, K, C]
        out_features = torch.bmm(S.transpose(1, 2), x)

        # 3. DYNAMIC DOMAIN MUTATION:
        # To strictly follow the GDL Blueprint, the geometric domain MUST shrink.

        # 3.1 Coarsen the physical 3D coordinates (Vertices)
        # Average S across the batch to get a deterministic geometric transformation.
        # Shape: [N, K]
        S_mean = S.mean(dim=0)
        # Normalize columns so they sum to 1 (Weighted Average of Coordinates)
        S_mean_norm = S_mean / (S_mean.sum(dim=0, keepdim=True) + 1e-8)

        # V_new = S_norm^T @ V
        # Shape: [K, N] @ [N, 3] -> [K, 3]
        new_vertices = torch.matmul(S_mean_norm.T, domain.vertices)

        # 3.2 Coarsen the Intrinsic Laplacian Matrix
        # L_coarse = S_norm^T @ L_old @ S_norm
        # Shape: [K, N] @ [N, N] @ [N, K] -> [K, K]
        L_old = domain.laplacian
        L_coarse = torch.matmul(S_mean_norm.T, torch.matmul(L_old, S_mean_norm))

        # To guarantee spectral stability for the next ChebNet layer, we must
        # ensure the new Laplacian is symmetrically normalized with spectrum in [0, 2].
        # We extract a pseudo-adjacency matrix from L_coarse (assuming L_old = I - A_old)
        I_K = torch.eye(self.target_nodes, device=L_coarse.device, dtype=L_coarse.dtype)
        A_coarse = I_K - L_coarse

        # Ensure non-negativity and remove self-loops
        A_coarse = torch.relu(A_coarse) * (1.0 - I_K)

        # Symmetrically normalize the new pseudo-adjacency
        eps = 1e-6
        D_vec = torch.sum(A_coarse, dim=1).clamp(min=eps)
        D_inv_sqrt = torch.diag(torch.pow(D_vec, -0.5))
        new_laplacian = I_K - torch.matmul(D_inv_sqrt, torch.matmul(A_coarse, D_inv_sqrt))

        # 3.3 Construct the new GeodesicDomain
        # The explicit faces (triangles) are abstracted during this soft-pooling.
        # We pass an empty tensor for faces, but crucially, we preserve the
        # INTRINSIC GEOMETRY by explicitly injecting the stabilized pooled laplacian!
        empty_faces = torch.empty((0, 3), dtype=torch.long, device=new_vertices.device)

        coarsened_domain = GeodesicDomain(
            vertices=new_vertices,
            faces=empty_faces,
            laplacian=new_laplacian,
        )

        return GeodesicSignal(domain=coarsened_domain, features=out_features)


class GeodesicGlobalPooling(GlobalInvariant[GeodesicSignal]):
    """
    4G Global Invariant: Destroys the intrinsic geometry completely.

    Aggregates all N vertex features into a single fixed-size vector per batch.
    Guarantees absolute Isometry Invariance: f(Iso(x)) = f(x).
    """

    def __init__(self, aggregation: str = "mean"):
        super().__init__()
        valid_aggs = ["sum", "mean", "max"]
        if aggregation not in valid_aggs:
            raise ValueError(f"Unsupported aggregation: {aggregation}")
        self.aggregation = aggregation

    def forward(self, signal: GeodesicSignal) -> torch.Tensor:
        x = signal.features  # Shape: [Batch, Nodes(N), Channels(C)]

        # Reduce the spatial (vertex) dimension (dim=1)
        if self.aggregation == "sum":
            invariant_tensor = torch.sum(x, dim=1)
        elif self.aggregation == "mean":
            invariant_tensor = torch.mean(x, dim=1)
        elif self.aggregation == "max":
            invariant_tensor = torch.max(x, dim=1)[0]
        else:
            raise RuntimeError("Unreachable aggregation state.")

        # The manifold structure is fully destroyed. Return a raw vector.
        # Final Shape: [Batch, Channels]
        return invariant_tensor
