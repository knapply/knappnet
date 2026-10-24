"""Distance indices on the engine: Schultz, Gutman, and hyper-Wiener
indices, harmonic diameter, local efficiency, and closeness vitality.

The engine's indices take NetworkX's definitions (netcapi.h, "Distance
measures": index 1 = Schultz, 2 = Gutman, 4 = hyper-Wiener as NetworkX
computes it, 5 = harmonic diameter).
"""

from __future__ import annotations

import math

from . import _checks as chk
from . import _core
from ._algorithms import (
    ALL, ANY, OUT, SIMPLE, _connected, _hashable, _node_dict, _positive_weights, _weighted,
)
from ._graph import implements, threads
from ._paths import doubles

DISTANCE = 0
SCHULTZ, GUTMAN, HYPER_WIENER_NX, HARMONIC_DIAMETER = 1, 2, 4, 5

_WEIGHT = chk.weight_attr()
_UNDIRECTED = chk.all_of(SIMPLE, chk.undirected, _WEIGHT)


def _index(G, index, weight):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    return _core.distance_index(G.core, index, OUT, w, DISTANCE, threads())


def _index_if_connected(G, index, weight):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    if not _connected(G):
        return math.inf
    return _core.distance_index(G.core, index, OUT, w, DISTANCE, threads())


@implements("schultz_index", _UNDIRECTED, declined="directed graphs; multigraphs; weight functions; weights <= 0.")
def schultz_index(G, weight=None):
    return _index_if_connected(G, SCHULTZ, weight)


@implements("gutman_index", _UNDIRECTED, declined="directed graphs; multigraphs; weight functions; weights <= 0.")
def gutman_index(G, weight=None):
    return _index_if_connected(G, GUTMAN, weight)


@implements("hyper_wiener_index", _UNDIRECTED, declined="directed graphs; multigraphs; weight functions; weights <= 0.")
def hyper_wiener_index(G, weight=None):
    return _index_if_connected(G, HYPER_WIENER_NX, weight)


@implements(
    "harmonic_diameter",
    chk.all_of(ANY, _WEIGHT, chk.is_none("sp")),
    declined="sp; weight functions; weights <= 0.",
)
def harmonic_diameter(G, sp=None, *, weight=None):
    n = len(G)
    if n < 2:
        return math.nan
    return _index(G, HARMONIC_DIAMETER, weight)


@implements("local_efficiency", chk.all_of(SIMPLE, chk.undirected), declined="directed graphs; multigraphs.")
def local_efficiency(G):
    if len(G) == 0:
        raise NotImplementedError("empty graph")
    return _core.local_efficiency(G.core, 1, threads())


def _check_vitality(a):
    G = chk.graph_of(a)
    v = a.get("node")
    if v is not None and not (_hashable(v) and v in G):
        return "node is not a node of the graph"
    return None


@implements(
    "closeness_vitality",
    chk.all_of(SIMPLE, _WEIGHT, chk.is_none("wiener_index"), _check_vitality, chk.nonempty),
    declined="multigraphs; weight functions; weights <= 0; wiener_index.",
)
def closeness_vitality(G, node=None, weight=None, wiener_index=None):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    if len(G) < 2:
        raise NotImplementedError("fewer than two nodes")
    vals = doubles(_core.closeness_vitality(G.core, OUT, w, DISTANCE, threads()))
    result = _node_dict(G, vals)
    if node is None:
        return result
    return result[node]
