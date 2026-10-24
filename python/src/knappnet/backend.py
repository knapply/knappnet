"""The NetworkX backend interface (the ``networkx.backends`` entry point).

NetworkX looks up each dispatchable function by name on this class and
calls ``can_run`` and ``should_run`` before converting a graph with
``convert_from_nx``. Conversions are cached by NetworkX in
``G.__networkx_cache__``, so repeated calls on an unchanged graph convert
it once.
"""

from __future__ import annotations

import networkx as nx

from . import _algorithms  # noqa: F401  (registers the implementations)
from ._checks import plain_storage
from ._graph import (
    BACKEND_NAME,
    CHECKS,
    IMPLEMENTATIONS,
    BackendGraph,
    bind,
    convert_from_nx,
    convert_to_nx,
)


def _graph_argument(name, arguments):
    func = nx.utils.backends._registered_algorithms[name]
    for gname in func.graphs:
        G = arguments.get(gname)
        if G is not None:
            return G
    return None


def can_run(name, args, kwargs):
    if name not in IMPLEMENTATIONS:
        return False
    try:
        arguments = bind(name, args, kwargs)
    except TypeError:
        return "the arguments do not match NetworkX's signature"
    G = _graph_argument(name, arguments)
    if isinstance(G, BackendGraph):
        G = G.graph
    if G is None or not isinstance(G, nx.Graph):
        return "the graph argument is not a NetworkX graph"
    if not plain_storage(G):
        return "graph classes that override NetworkX's storage are left to NetworkX"
    check = CHECKS.get(name)
    reason = check(arguments) if check is not None else None
    return reason if reason else True


# Linear-time functions, which NetworkX runs about as fast as one conversion
# costs: with automatic dispatch they run on the engine only once the graph
# has a cached conversion (a previous call converted it).
LINEAR_TIME = frozenset({
    "degree_centrality", "in_degree_centrality", "out_degree_centrality",
    "connected_components", "number_connected_components", "is_connected",
    "node_connected_component", "strongly_connected_components",
    "number_strongly_connected_components", "is_strongly_connected",
    "weakly_connected_components", "number_weakly_connected_components",
    "is_weakly_connected", "core_number", "single_source_shortest_path_length",
    "single_target_shortest_path_length", "single_source_dijkstra_path_length",
    "has_path", "descendants", "ancestors", "articulation_points", "bridges",
    "has_bridges", "biconnected_components", "is_biconnected", "overall_reciprocity",
    "reciprocity",
})


def _has_cached_conversion(G):
    cache = getattr(G, "__networkx_cache__", None)
    if not cache:
        return False
    return bool(cache.get("backends", {}).get(BACKEND_NAME))


def should_run(name, args, kwargs):
    try:
        arguments = bind(name, args, kwargs)
    except TypeError:
        return False
    G = _graph_argument(name, arguments)
    if isinstance(G, BackendGraph):
        return True
    cfg = getattr(nx.config.backends, BACKEND_NAME, None)
    min_nodes = getattr(cfg, "min_nodes", 0) if cfg is not None else 0
    if G is not None and len(G) < (min_nodes or 0):
        return f"graphs with fewer than {min_nodes} nodes run faster in NetworkX"
    if name in LINEAR_TIME and G is not None and not _has_cached_conversion(G):
        return "linear-time function on a graph not yet converted"
    return True


def on_start_tests(items):
    """Marks NetworkX tests whose expectations the backend does not share."""
    from ._testing import mark_known_differences

    mark_known_differences(items)


class BackendInterface:
    convert_from_nx = staticmethod(convert_from_nx)
    convert_to_nx = staticmethod(convert_to_nx)
    can_run = staticmethod(can_run)
    should_run = staticmethod(should_run)
    on_start_tests = staticmethod(on_start_tests)


for _name, _func in IMPLEMENTATIONS.items():
    setattr(BackendInterface, _name, staticmethod(_func))
del _name, _func
