import torch
import torch.nn as nn

from gdl.core.base_blueprint import (
    EquivariantLayer,
    GeometricNonLinearity,
    LocalPooling,
    GlobalInvariant,
)
from .domain import SetDomain, GraphDomain
from .signal import NodeSignal


class DeepSetsLayer(EquivariantLayer[NodeSignal]):
    """
    3G Equivariant Layer (Sets): Deep Sets formulation.

    Guarantees Permutation Equivariance for unordered collections of nodes without edges.
    Mathematically proven by Zaheer et al. (2017), the only linear permutation
    equivariant operation on a set is of the form:
    f(X)_i = W1 * X_i + W2 * aggregate(X)

    Since the aggregation function (sum/mean/max) is symmetric, shuffling
    the input order of nodes simply shuffles the output order identically.
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        aggregation: str = "mean",
        bias: bool = True,
    ) -> None:
        super().__init__()

        valid_aggs = ["sum", "mean", "max"]
        if aggregation not in valid_aggs:
            raise ValueError(f"Unsupported aggregation: {aggregation}")
        self.aggregation = aggregation

        # W1: operates on each node strictly individually
        self.lin_self = nn.Linear(in_channels, out_channels, bias=bias)

        # W2: operates on the global invariant context of the entire set
        self.lin_global = nn.Linear(in_channels, out_channels, bias=False)

    def forward(self, signal: NodeSignal) -> NodeSignal:
        domain: SetDomain = signal.domain

        # Shape is strictly [Batch, Nodes, Channels]
        x = signal.features

        # 1. Independent node transformation: W1 @ X_i
        node_features = self.lin_self(x)

        # 2. Symmetric global aggregation: Pool(X) -> [Batch, 1, Channels]
        if self.aggregation == "sum":
            global_context = torch.sum(x, dim=1, keepdim=True)
        elif self.aggregation == "mean":
            global_context = torch.mean(x, dim=1, keepdim=True)
        elif self.aggregation == "max":
            global_context = torch.max(x, dim=1, keepdim=True)[0]
        else:
            raise RuntimeError("Unreachable aggregation state.")

        # 3. Transform the context: W2 @ Pool(X) -> [Batch, 1, Channels]
        global_features = self.lin_global(global_context)

        # 4. Broadcast the global context back to all nodes (Broadcasting over dim=1)
        # Final shape: [Batch, Nodes, Channels]
        out_features = node_features + global_features

        return NodeSignal(domain=domain, features=out_features)


class GraphMessagePassingLayer(EquivariantLayer[NodeSignal]):
    """
    3G Equivariant Layer (Graphs): Local Message Passing over Edges.

    Guarantees Permutation Equivariance for irregular topologies (Graphs) defined by
    an Adjacency Matrix (edge_index).

    Mathematical GDL Formulation (e.g., simplified GCN / GraphSAGE):
    h_i^{(l+1)} = W_self * h_i^{(l)} + Aggregate_{j in N(i)}( W_neigh * h_j^{(l)} )

    This pure PyTorch implementation demonstrates the core mechanism of GNNs
    without relying on PyTorch Geometric's heavy MessagePassing base class.
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        aggr: str = "sum",
        bias: bool = True,
    ) -> None:
        super().__init__()

        valid_aggs = ["sum", "mean"]
        if aggr not in valid_aggs:
            raise ValueError(
                f"GraphMessagePassingLayer only supports {valid_aggs}, got '{aggr}'."
            )
        self.aggr = aggr

        # W_self: Transforms the central node's own features
        self.lin_self = nn.Linear(in_channels, out_channels, bias=bias)

        # W_neigh: Transforms the neighboring nodes' features (the messages)
        self.lin_neigh = nn.Linear(in_channels, out_channels, bias=False)

    def forward(self, signal: NodeSignal) -> NodeSignal:
        domain = signal.domain

        # This layer strictly requires edges to function correctly!
        if not isinstance(domain, GraphDomain):
            raise TypeError(
                f"GraphMessagePassingLayer requires a GraphDomain with edges. "
                f"Got {type(domain).__name__}. If you have no edges, use DeepSetsLayer."
            )

        # x shape: [Batch, Nodes, Channels]
        x = signal.features
        batch_size, num_nodes, _ = x.shape
        out_channels = self.lin_self.out_features

        # 1. Update Function (Self): Transform root node's own features
        out_self = self.lin_self(x)

        # edge_index shape: [2, NumEdges]. row: source (j), col: target (i)
        row, col = domain.edge_index

        # Edge case: Graph with no edges
        if row.size(0) == 0:
            return NodeSignal(domain=domain, features=out_self)

        # 2. Message Function: Extract and transform neighbor features
        # We index the node dimension (dim=1). msg_j shape: [Batch, NumEdges, Channels]
        x_j = x[:, row, :]
        msg_j = self.lin_neigh(x_j)

        # Apply edge weights if present
        if domain.edge_weight is not None:
            # Reshape for broadcasting: [NumEdges] -> [1, NumEdges, 1]
            weights = domain.edge_weight.view(1, -1, 1)
            msg_j = msg_j * weights

        # 3. Aggregation Function (Scatter Add)
        # Prepare an empty tensor for aggregated messages: [Batch, Nodes, Channels]
        aggr_msg = torch.zeros(
            batch_size, num_nodes, out_channels,
            device=x.device, dtype=x.dtype,
        )

        # Expand col (target indices) to match the message channel dimensions
        # col shape: [NumEdges] -> [1, NumEdges, 1] -> [Batch, NumEdges, Channels]
        scatter_idx = col.view(1, -1, 1).expand(batch_size, -1, out_channels)

        # Scatter the messages into the target nodes using PyTorch's built-in add
        aggr_msg.scatter_add_(dim=1, index=scatter_idx, src=msg_j)

        if self.aggr == "mean":
            # Compute degrees (number of incoming edges) to calculate the mean
            degree = torch.zeros(num_nodes, device=x.device, dtype=x.dtype)
            ones = torch.ones_like(col, dtype=x.dtype)
            degree.scatter_add_(dim=0, index=col, src=ones)

            # Avoid division by zero by clamping degree to min=1.0
            degree = degree.clamp(min=1.0).view(1, -1, 1)

            # Divide aggregated sum by degree
            aggr_msg = aggr_msg / degree

        # 4. Final Update: Combine self features and aggregated neighbor messages
        out_features = out_self + aggr_msg

        return NodeSignal(domain=domain, features=out_features)


class TransformerEquivariantLayer(EquivariantLayer[NodeSignal]):
    """
    3G Equivariant Layer (Sets/Complete Graphs): Multi-Head Self-Attention.

    A Transformer is mathematically equivalent to a Graph Neural Network (GNN)
    operating on a Complete Graph (where every node is connected to every other node),
    with dynamically computed edge weights (attention scores) based on node features.

    Without Positional Encodings, treating nodes purely as a set, this layer
    is perfectly Permutation Equivariant: f(P * X) = P * f(X).
    """

    def __init__(
        self,
        embed_dim: int,
        num_heads: int,
        dropout: float = 0.0,
        bias: bool = True,
    ) -> None:
        super().__init__()

        # We rely on PyTorch's highly optimized scaled dot-product attention
        self.mha = nn.MultiheadAttention(
            embed_dim=embed_dim,
            num_heads=num_heads,
            dropout=dropout,
            bias=bias,
            batch_first=True,
        )

        # Standard structural component of a Transformer block (Update Function)
        self.norm = nn.LayerNorm(embed_dim)

    def forward(self, signal: NodeSignal) -> NodeSignal:
        domain = signal.domain

        # x shape strictly guaranteed: [Batch, Nodes, Channels(embed_dim)]
        x = signal.features

        # 1. Message Passing on a Complete Graph (Self-Attention)
        # Each node (Query) gathers information from all other nodes (Keys/Values).
        attn_output, _ = self.mha(query=x, key=x, value=x, need_weights=False)

        # 2. Update Function: Residual connection + Layer Normalization
        # Out shape: [Batch, Nodes, Channels]
        out_features = self.norm(x + attn_output)

        # 3. Return a new Signal dynamically bound to the same abstract geometry
        return NodeSignal(domain=domain, features=out_features)


class GraphNonLinearity(GeometricNonLinearity[NodeSignal]):
    """
    3G Non-Linearity: Point-wise activation on Graph/Set node features.

    Since node features are scalar fields per channel, standard ReLU
    preserves Permutation Equivariance (point-wise ops commute with permutation).
    """

    def __init__(self, inplace: bool = False) -> None:
        super().__init__()
        self.activation = nn.ReLU(inplace=inplace)

    def forward(self, signal: NodeSignal) -> NodeSignal:
        out_features = self.activation(signal.features)
        return NodeSignal(domain=signal.domain, features=out_features)


class GraphCoarsening(LocalPooling[NodeSignal]):
    """
    3G Local Pooling: Graph/Set Coarsening (Node Clustering).

    Uses a learnable soft-assignment matrix to cluster N nodes into K super-nodes.
    This is mathematically inspired by differentiable pooling (DiffPool).
    Since the assignment matrix is computed via a permutation-equivariant linear layer,
    the resulting K super-nodes are perfectly invariant to the input node order.

    The geometry dynamically mutates: A Graph/Set of N nodes becomes a Set of K nodes.
    """

    def __init__(self, in_channels: int, num_clusters: int) -> None:
        super().__init__()
        if num_clusters <= 0:
            raise ValueError("num_clusters must be positive.")
        self.num_clusters = num_clusters

        # Scorer calculates the probability of each node belonging to each of the K clusters
        self.scorer = nn.Linear(in_channels, num_clusters)

    def forward(self, signal: NodeSignal) -> NodeSignal:
        x = signal.features  # Shape: [Batch, Nodes(N), Channels(C)]

        # S: Assignment matrix. Shape: [Batch, N, K]
        # S_b,i,k is the probability that node i belongs to cluster k
        S = torch.softmax(self.scorer(x), dim=-1)

        # X_coarsened = S^T @ X
        # Shape: [Batch, K, N] @ [Batch, N, C] -> [Batch, K, C]
        out_features = torch.bmm(S.transpose(1, 2), x)

        # DYNAMIC DOMAIN MUTATION:
        # The original graph or set topology is collapsed into K disconnected super-nodes.
        new_domain = SetDomain(num_nodes=self.num_clusters)

        return NodeSignal(domain=new_domain, features=out_features)


class GraphGlobalPooling(GlobalInvariant[NodeSignal]):
    """
    3G Global Invariant: Destroys the node geometry completely.

    Aggregates all N node features into a single fixed-size vector per batch.
    Guarantees absolute Permutation Invariance: f(P * X) = f(X).
    """

    def __init__(self, aggregation: str = "mean") -> None:
        super().__init__()
        valid_aggs = ["sum", "mean", "max"]
        if aggregation not in valid_aggs:
            raise ValueError(f"Unsupported aggregation: {aggregation}")
        self.aggregation = aggregation

    def forward(self, signal: NodeSignal) -> torch.Tensor:
        x = signal.features  # Shape: [Batch, Nodes(N), Channels(C)]

        # Reduce the spatial (node) dimension (dim=1)
        if self.aggregation == "sum":
            invariant_tensor = torch.sum(x, dim=1)
        elif self.aggregation == "mean":
            invariant_tensor = torch.mean(x, dim=1)
        elif self.aggregation == "max":
            # torch.max returns a tuple (values, indices)
            invariant_tensor = torch.max(x, dim=1)[0]
        else:
            raise RuntimeError("Unreachable aggregation state.")

        # The node structure (Set/Graph) is fully destroyed. Return a raw vector.
        # Final Shape: [Batch, Channels]
        return invariant_tensor
