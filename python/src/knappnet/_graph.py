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

    __slots__ = ("graph", "core", "weight", "default")

    def __init__(self, graph, core, weight, default):
        self.graph = graph
        self.core = core
        self.weight = weight
        self.default = default

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
    if isinstance(edge_attrs, dict) and edge_attrs:
        if len(edge_attrs) > 1:
            raise NotImplementedError("only one edge attribute is supported")
        ((weight, default),) = edge_attrs.items()
        if default is None:
            default = math.nan
    elif isinstance(preserve_edge_attrs, dict):
        attrs = preserve_edge_attrs.get(graph_name) or {}
        if len(attrs) == 1:
            ((weight, default),) = attrs.items()
    directed = G.is_directed()
    adj = G._succ if directed else G._adj
    core = _core.from_networkx(G._node, adj, directed, G.is_multigraph(), weight,
                               float(default))
    return BackendGraph(G, core, weight, default)


def convert_to_nx(obj, *, name=None):
    """Results are already NetworkX objects; backend graphs give back their source."""
    if isinstance(obj, BackendGraph):
        return obj.graph
    return obj


# ---- Implemented functions ------------------------------------------------

IMPLEMENTATIONS = {}  # NetworkX dispatch name -> implementation
CHECKS = {}  # NetworkX dispatch name -> can_run check (bound arguments -> reason)
NOTES = {}  # NetworkX dispatch name -> what is supported, for the docs


def implements(name, check=None, notes=""):
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
        NOTES[name] = notes
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
