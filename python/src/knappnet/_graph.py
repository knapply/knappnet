"""The backend graph and the table of implemented functions.

A NetworkX graph is converted once per edge attribute (NetworkX caches the
result in ``G.__networkx_cache__``). The converted graph keeps a reference
to the original, so node attributes are read from it at call time and
``convert_to_nx`` returns it unchanged.
"""

from __future__ import annotations

import inspect
import math

import networkx as nx

from . import _core

BACKEND_NAME = "knappnet"



class BackendGraph:
    """A NetworkX graph converted to an engine graph (internal)."""

    __networkx_backend__ = BACKEND_NAME
    __networkx_cache__ = None

    __slots__ = ("graph", "core", "weight", "default", "attrs", "_order", "_cores")

    def __init__(self, graph, core, weight, default, attrs=None):
        self.graph = graph
        self.core = core
        self.weight = weight
        self.default = default
        # Every edge attribute the call asked for ({name: default}); the
        # core holds the first when there is only one.
        self.attrs = attrs or ({weight: default} if weight is not None else {})
        self._order = None
        self._cores = None

    def weighted_core(self, attr, default):
        """A core whose edge weights are the attribute ``attr`` (``default``
        where missing); the core itself if it already holds them."""
        if attr is None or (attr == self.weight and self.core.weighted):
            return self.core
        if self._cores is None:
            self._cores = {}
        key = (attr, float(default))
        core = self._cores.get(key)
        if core is None:
            G = self.graph
            directed = self.core.directed
            try:
                core = _core.from_networkx(G._node, G._succ if directed else G._adj, directed,
                                           G.is_multigraph(), attr, float(default))
            except (TypeError, ValueError) as err:
                raise NotImplementedError("edge attribute is not a number") from err
            self._cores[key] = core
        return core

    def edge_values(self, attr, default):
        """One float per edge of the core, in its edge order: the edge
        attribute ``attr`` (``default`` where missing), read from the
        NetworkX graph unless the core already holds it."""
        if attr is None:
            return [float(default)] * self.core.m
        if attr == self.weight and self.core.weighted:
            return list(memoryview(_core.edge_weights(self.core)).cast("d"))
        adj = self.graph._succ if self.core.directed else self.graph._adj
        out = []
        multi = self.graph.is_multigraph()
        for u, v in _core.edge_keys(self.core):
            data = adj[u][v]
            if multi:
                raise NotImplementedError("edge attributes of multigraphs are not read here")
            x = data.get(attr, default)
            try:
                out.append(float(x))
            except (TypeError, ValueError) as err:
                raise NotImplementedError("edge attribute is not a number") from err
        return out

    @property
    def input_order(self):
        """True if the engine's input-order traversals visit neighbours as
        NetworkX does: a node's successors (or neighbours) in its adjacency
        order, and predecessors by increasing node position."""
        if self._order is None:
            self._order = _adjacency_in_input_order(self)
        return self._order

    # Read-only conveniences used by the implementations.
    @property
    def nodes(self):
        return self.core.nodes

    @property
    def index(self):
        return self.core.index

    def __len__(self):
        return self.core.n

    def is_directed(self):
        return self.core.directed

    def is_multigraph(self):
        return self.graph.is_multigraph()

    def __repr__(self):
        return (f"<knappnet graph: {self.core.n} nodes, {self.core.m} edges, "
                f"{'directed' if self.core.directed else 'undirected'}, "
                f"weight={self.weight!r}>")


def _adjacency_in_input_order(B):
    """Whether a graph rebuilt from its edges in G.edges() order has the same
    neighbour order as G (see knappnet's engine order.rs)."""
    G = B.graph
    index = B.core.index
    if G.is_multigraph():
        return False
    if not B.core.directed:
        for v, nbrs in G._adj.items():
            iv = index[v]
            last, high = -1, False
            for u in nbrs:
                iu = index[u]
                if iu < iv:
                    if high or iu < last:
                        return False
                    last = iu
                else:
                    high = True
        return True
    for nbrs in G._pred.values():
        last = -1
        for u in nbrs:
            iu = index[u]
            if iu < last:
                return False
            last = iu
    return True


def convert_from_nx(G, edge_attrs=None, node_attrs=None, preserve_edge_attrs=False,
                    preserve_node_attrs=False, preserve_graph_attrs=False, name=None,
                    graph_name=None):
    """Converts a NetworkX graph, reading at most one numeric edge attribute."""
    if isinstance(G, BackendGraph):
        return G
    from ._checks import plain_storage

    if not isinstance(G, nx.Graph) or not plain_storage(G):
        raise NotImplementedError(f"cannot convert {type(G).__name__}")
    weight, default = None, 1.0
    attrs = None
    if isinstance(edge_attrs, dict) and edge_attrs:
        attrs = {k: (math.nan if d is None else d) for k, d in edge_attrs.items()}
        if len(edge_attrs) > 1:
            # Read through BackendGraph.edge_values instead.
            weight = None
        else:
            ((weight, default),) = attrs.items()
    elif isinstance(preserve_edge_attrs, dict):
        attrs = preserve_edge_attrs.get(graph_name) or {}
        if len(attrs) == 1:
            ((weight, default),) = attrs.items()
    directed = G.is_directed()
    adj = G._succ if directed else G._adj
    core = _core.from_networkx(G._node, adj, directed, G.is_multigraph(), weight,
                               float(default))
    return BackendGraph(G, core, weight, default, attrs)


def convert_to_nx(obj, *, name=None):
    """Results are already NetworkX objects; backend graphs give back their source."""
    if isinstance(obj, BackendGraph):
        return obj.graph
    return obj


# ---- Implemented functions ------------------------------------------------

IMPLEMENTATIONS = {}  # NetworkX dispatch name -> implementation
CHECKS = {}  # NetworkX dispatch name -> can_run check (bound arguments -> reason)
NOTES = {}  # NetworkX dispatch name -> what is supported, for the docs


def implements(name, check=None, declined=""):
    """Registers an implementation of the NetworkX function ``name``.

    ``check`` receives the call's arguments bound to NetworkX's signature
    (defaults applied) and returns None to accept the call or a string
    saying why the backend declines it, in which case NetworkX runs its own
    implementation.
    """

    def wrap(fn):
        if name in IMPLEMENTATIONS:
            raise ValueError(f"{name} implemented twice")
        IMPLEMENTATIONS[name] = fn
        CHECKS[name] = check
        NOTES[name] = declined
        return fn

    return wrap


_SIGNATURES = {}


def bind(name, args, kwargs):
    """Binds a call to NetworkX's signature of ``name`` (defaults applied)."""
    sig = _SIGNATURES.get(name)
    if sig is None:
        func = nx.utils.backends._registered_algorithms[name]
        sig = _SIGNATURES[name] = inspect.signature(func.orig_func)
    bound = sig.bind(*args, **kwargs)
    bound.apply_defaults()
    return bound.arguments


def threads():
    """Thread count for engine calls (0 = every core)."""
    cfg = getattr(nx.config.backends, BACKEND_NAME, None)
    value = getattr(cfg, "threads", 0) if cfg is not None else 0
    return int(value or 0)
