"""Cores, assortativity and mixing, neighbour-degree measures, cut measures,
square clustering, generalized degree, matchings, and greedy colouring.

Graph-valued results (the core and truss subgraphs) are NetworkX graphs of
the input's class. The functions that return a matching or a colouring run
on the engine only where NetworkX's answer is determined by the graph (the
maximal matching scans the edges in order; the colourings that break ties
by node order need adjacency in input order).
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Collection

import networkx as nx

from . import _checks as chk
from . import _core
from ._algorithms import ALL, IN, OUT, SIMPLE, _hashable, _node_dict, _weighted
from ._graph import implements, threads
from ._paths import doubles

_WEIGHT = chk.weight_attr()
_SIMPLE_NO_LOOPS = chk.all_of(SIMPLE, chk.no_selfloops)
_UNDIRECTED_SIMPLE = chk.all_of(SIMPLE, chk.undirected)

_KINDS = {"out": OUT, "in": IN, "in+out": ALL}


# ---- Helpers --------------------------------------------------------------


def _indices(G, nodes):
    index = G.index
    return [index[v] for v in nodes]


def _integral_weights(G, weight):
    """True if every value of the edge attribute (1 where missing) is an
    int, in which case NetworkX's sums over it are ints."""
    if weight is None:
        return True
    adj = G.graph._succ if G.core.directed else G.graph._adj
    return all(type(d.get(weight, 1)) is int for nbrs in adj.values() for d in nbrs.values())


def _codes(values):
    """Distinct values as codes 0..k-1 (equal as Python dict keys), with the
    value of each code; unhashable or NaN values are left to NetworkX."""
    seen = {}
    labels = []
    out = []
    for x in values:
        try:
            if x != x:
                raise NotImplementedError("NaN attribute values")
            c = seen.get(x)
        except TypeError as err:
            raise NotImplementedError("unhashable attribute values") from err
        if c is None:
            c = seen[x] = len(labels)
            labels.append(x)
        out.append(float(c))
    return out, labels


def _degree_kinds(G, x, y):
    if G.is_directed():
        return _KINDS[x], _KINDS[y]
    return ALL, ALL


def _check_xy(a):
    G = chk.graph_of(a)
    if G.is_directed():
        for p in ("x", "y"):
            if a.get(p) not in ("in", "out"):
                return f"{p} is not 'in' or 'out'"
    return None


# ---- Cores ------------------------------------------------------------------

_CORE_CHECK = chk.all_of(_SIMPLE_NO_LOOPS, chk.is_none("core_number"), chk.optional_integer("k"))


def _core_graph(G, which, k):
    nodes = _core.core_subset(G.core, which, 0 if k is None else k, int(k is not None))
    names = G.nodes
    return G.graph.subgraph([names[i] for i in nodes]).copy()


@implements("k_core", _CORE_CHECK, declined="multigraphs; self-loops; core_number.")
def k_core(G, k=None, core_number=None):
    return _core_graph(G, 0, k)


@implements("k_shell", _CORE_CHECK, declined="multigraphs; self-loops; core_number.")
def k_shell(G, k=None, core_number=None):
    return _core_graph(G, 1, k)


@implements("k_crust", _CORE_CHECK, declined="multigraphs; self-loops; core_number.")
def k_crust(G, k=None, core_number=None):
    return _core_graph(G, 2, k)


@implements(
    "k_corona",
    chk.all_of(_SIMPLE_NO_LOOPS, chk.is_none("core_number"), chk.integer("k")),
    declined="multigraphs; self-loops; core_number.",
)
def k_corona(G, k, core_number=None):
    return _core_graph(G, 3, k)


@implements(
    "k_truss",
    chk.all_of(_UNDIRECTED_SIMPLE, chk.no_selfloops, chk.integer("k")),
    declined="directed graphs; multigraphs; self-loops.",
)
def k_truss(G, k):
    edges, nodes = _core.k_truss(G.core, k)
    keys = _core.edge_keys(G.core)
    keep = set(edges)
    H = G.graph.copy()
    H.remove_edges_from([keys[i] for i in range(len(keys)) if i not in keep])
    keep_nodes = {G.nodes[i] for i in nodes}
    H.remove_nodes_from([v for v in G.nodes if v not in keep_nodes])
    return H


@implements(
    "onion_layers",
    chk.all_of(_UNDIRECTED_SIMPLE, chk.no_selfloops),
    declined="directed graphs; multigraphs; self-loops.",
)
def onion_layers(G):
    return _node_dict(G, _core.onion_layers(G.core))


# ---- Assortativity and mixing ------------------------------------------------

_NODES = chk.optional_nodes("nodes")


def _subset_args(G, nodes):
    return None if nodes is None else _indices(G, nodes)


def _mixing_dict(values, counts, k, labels, normalized):
    c = memoryview(counts).cast("d").tolist()
    d = {labels[i]: {} for i in range(k)}
    total = 0
    for i in range(k):
        row = d[labels[i]]
        for j in range(k):
            n = c[i * k + j]
            if n:
                row[labels[j]] = int(n)
                total += int(n)
    if normalized:
        for row in d.values():
            for j in row:
                row[j] /= total
    return d


def _degree_mixing(G, x, y, weight, nodes):
    if weight is not None:
        raise NotImplementedError("weighted degrees")
    kx, ky = _degree_kinds(G, x, y)
    values, counts, k = _core.mixing(G.core, 0, kx, ky, 0, None, _subset_args(G, nodes))
    labels = [int(v) for v in memoryview(values).cast("d")]
    return counts, k, labels


@implements(
    "degree_mixing_dict",
    chk.all_of(SIMPLE, chk.is_none("weight"), _check_xy, _NODES),
    declined="multigraphs; weights.",
)
def degree_mixing_dict(G, x="out", y="in", weight=None, nodes=None, normalized=False):
    counts, k, labels = _degree_mixing(G, x, y, weight, nodes)
    return _mixing_dict(None, counts, k, labels, normalized)


def _matrix(G, counts, k, labels, mapping, normalized):
    import numpy as np

    n = len(mapping)
    a = np.zeros((n, n))
    c = memoryview(counts).cast("d").tolist()
    for i in range(k):
        mi = mapping.get(labels[i])
        if mi is None:
            continue
        for j in range(k):
            mj = mapping.get(labels[j])
            if mj is not None and c[i * k + j]:
                a[mi, mj] = c[i * k + j]
    if normalized:
        a = a / a.sum()
    return a


def _mapping_ok(a):
    m = a.get("mapping")
    if not isinstance(m, dict) or not m:
        return "without a mapping the matrix order is NetworkX's own"
    idx = sorted(m.values()) if all(isinstance(v, int) and not isinstance(v, bool) for v in m.values()) else None
    if idx is None or idx != list(range(len(m))):
        return "mapping must number its keys 0..n-1"
    return None


@implements(
    "degree_mixing_matrix",
    chk.all_of(SIMPLE, chk.is_none("weight"), _check_xy, _NODES, _mapping_ok),
    declined="multigraphs; weights; no mapping.",
)
def degree_mixing_matrix(G, x="out", y="in", weight=None, nodes=None, normalized=True, mapping=None):
    counts, k, labels = _degree_mixing(G, x, y, weight, nodes)
    return _matrix(G, counts, k, labels, mapping, normalized)


def _node_values_none(G, attribute):
    node = G.graph._node
    return [node[v].get(attribute, None) for v in node]


def _attribute_mixing_run(G, attribute, nodes):
    values = _node_values_none(G, attribute)
    codes, labels = _codes(values)
    vals, counts, k = _core.mixing(G.core, 1, ALL, ALL, 0, codes, _subset_args(G, nodes))
    keys = [labels[int(v)] for v in memoryview(vals).cast("d")]
    return counts, k, keys


@implements(
    "attribute_mixing_dict",
    chk.all_of(SIMPLE, _NODES),
    declined="multigraphs; unhashable attribute values.",
)
def attribute_mixing_dict(G, attribute, nodes=None, normalized=False):
    counts, k, keys = _attribute_mixing_run(G, attribute, nodes)
    return _mixing_dict(None, counts, k, keys, normalized)


@implements(
    "attribute_mixing_matrix",
    chk.all_of(SIMPLE, _NODES, _mapping_ok),
    declined="multigraphs; unhashable attribute values; no mapping.",
)
def attribute_mixing_matrix(G, attribute, nodes=None, mapping=None, normalized=True):
    counts, k, keys = _attribute_mixing_run(G, attribute, nodes)
    return _matrix(G, counts, k, keys, mapping, normalized)


def _check_neighbor_kinds(source_default):
    def check(a):
        G = chk.graph_of(a)
        for p in ("source", "target"):
            v = a.get(p)
            if G.is_directed():
                if v not in ("in", "out", "in+out"):
                    return f"{p} is not 'in', 'out' or 'in+out'"
            elif v != source_default:
                return f"{p} is only supported for directed graphs"
        return None

    return check


@implements(
    "average_neighbor_degree",
    chk.all_of(SIMPLE, _WEIGHT, _check_neighbor_kinds("out"), _NODES),
    declined="multigraphs; weight functions; a single node as nodes.",
)
def average_neighbor_degree(G, source="out", target="out", nodes=None, weight=None):
    w = _weighted(G, weight)
    ks, kt = _degree_kinds(G, source, target)
    sel = list(G.graph) if nodes is None else list(nodes)
    vals = doubles(_core.neighbor_degree(G.core, ks, kt, w, _subset_args(G, nodes)))
    return dict(zip(sel, vals))


@implements(
    "average_degree_connectivity",
    chk.all_of(SIMPLE, _WEIGHT, _check_neighbor_kinds("in+out"), _NODES),
    declined="multigraphs; weight functions; a single node as nodes.",
)
def average_degree_connectivity(G, source="in+out", target="in+out", nodes=None, weight=None):
    w = _weighted(G, weight)
    ks, kt = _degree_kinds(G, source, target)
    degrees, values = _core.degree_connectivity(G.core, ks, kt, w, _subset_args(G, nodes))
    return {int(k): v for k, v in zip(doubles(degrees), doubles(values))}


# ---- Cut measures -----------------------------------------------------------


def _check_cut_sets(with_t):
    def check(a):
        G = chk.graph_of(a)
        for p in ("S", "T") if with_t else ("S",):
            v = a.get(p)
            if v is None and p == "T":
                continue
            if not isinstance(v, Collection) or isinstance(v, (str, bytes)):
                return f"{p} must be a collection of nodes"
            items = list(v)
            if not items:
                return f"{p} is empty"
            try:
                if not all(x in G for x in items) or len(set(items)) != len(items):
                    return f"{p} has repeated items or items that are not nodes"
            except TypeError:
                return f"{p} has unhashable items"
        return None

    return check


def _measures(G, S, T, weight):
    w = _weighted(G, weight)
    s = _indices(G, S)
    t = None if T is None else _indices(G, T)
    return doubles(_core.cut_measures(G.core, s, t, w))


def _finite(x):
    if not math.isfinite(x):
        raise NotImplementedError("division by zero is left to NetworkX")
    return x


def _total(G, weight, x):
    return int(x) if _integral_weights(G, weight) else x


_CUT = chk.all_of(SIMPLE, _WEIGHT, _check_cut_sets(True))
_CUT_S = chk.all_of(SIMPLE, _WEIGHT, _check_cut_sets(False))
_NODE_S = chk.all_of(SIMPLE, _check_cut_sets(False))


@implements("cut_size", _CUT, declined="multigraphs; weight functions; unusual S or T.")
def cut_size(G, S, T=None, weight=None):
    if G.is_directed() and T is None:
        raise NotImplementedError("T is required for directed graphs")
    return _total(G, weight, _measures(G, S, T, weight)[0])


@implements("volume", _CUT_S, declined="multigraphs; weight functions; unusual S.")
def volume(G, S, weight=None):
    return _total(G, weight, _measures(G, S, None, weight)[1])


@implements("normalized_cut_size", _CUT, declined="multigraphs; weight functions; unusual S or T; zero volume.")
def normalized_cut_size(G, S, T=None, weight=None):
    return _finite(_measures(G, S, T, weight)[3])


@implements("conductance", _CUT, declined="multigraphs; weight functions; unusual S or T; zero volume.")
def conductance(G, S, T=None, weight=None):
    return _finite(_measures(G, S, T, weight)[4])


@implements("edge_expansion", _CUT, declined="multigraphs; weight functions; unusual S or T.")
def edge_expansion(G, S, T=None, weight=None):
    return _finite(_measures(G, S, T, weight)[5])


@implements("mixing_expansion", _CUT, declined="multigraphs; weight functions; unusual S or T.")
def mixing_expansion(G, S, T=None, weight=None):
    if G.is_directed() and T is None:
        raise NotImplementedError("T is required for directed graphs")
    return _finite(_measures(G, S, T, weight)[6])


@implements("node_expansion", _NODE_S, declined="multigraphs; unusual S.")
def node_expansion(G, S):
    return _finite(_measures(G, S, None, None)[7])


@implements("boundary_expansion", _NODE_S, declined="multigraphs; unusual S.")
def boundary_expansion(G, S):
    return _finite(_measures(G, S, None, None)[8])


# ---- Clustering variants ----------------------------------------------------


@implements(
    "square_clustering",
    chk.all_of(_SIMPLE_NO_LOOPS, _NODES),
    declined="multigraphs; self-loops.",
)
def square_clustering(G, nodes=None):
    result = _node_dict(G, doubles(_core.square_clustering(G.core, threads())))
    if nodes is None:
        return result
    if _hashable(nodes) and nodes in G.graph:
        return result[nodes]
    return {v: result[v] for v in nodes}


@implements(
    "generalized_degree",
    chk.all_of(_UNDIRECTED_SIMPLE, chk.no_selfloops, _NODES),
    declined="directed graphs; multigraphs; self-loops.",
)
def generalized_degree(G, nodes=None):
    offsets, ks, counts = _core.generalized_degree(G.core, threads())
    out = {}
    for i, v in enumerate(G.nodes):
        c = Counter()
        for j in range(offsets[i], offsets[i + 1]):
            c[ks[j]] = counts[j]
        out[v] = c
    if nodes is None:
        return out
    if _hashable(nodes) and nodes in G.graph:
        return out[nodes]
    return {v: out[v] for v in nodes}


# ---- Matching ----------------------------------------------------------------


@implements(
    "maximal_matching",
    chk.all_of(_UNDIRECTED_SIMPLE),
    declined="directed graphs; multigraphs.",
)
def maximal_matching(G):
    keys = _core.edge_keys(G.core)
    return {keys[i] for i in _core.matching(G.core, 2, 0, 0)}


def _matching_pairs(a):
    G = chk.graph_of(a)
    m = a.get("matching")
    if isinstance(m, dict):
        m = list(m.items())
    elif not isinstance(m, (set, frozenset, list, tuple)):
        return "matching must be a set, list, tuple or dict"
    for e in m:
        if not isinstance(e, tuple) or len(e) != 2:
            return "matching has an item that is not a pair"
        try:
            if e[0] not in G or e[1] not in G:
                return "matching has a pair with a node outside the graph"
        except TypeError:
            return "matching has unhashable nodes"
    return None


_MATCHING = chk.all_of(SIMPLE, chk.undirected, _matching_pairs)


def _check_matching(G, matching):
    if isinstance(matching, dict):
        from networkx.algorithms.matching import matching_dict_to_set

        matching = matching_dict_to_set(matching)
    pairs = list(matching)
    index = G.index
    for u, v in pairs:
        if u == v:
            return (False, False, False)
    us = [index[u] for u, _ in pairs]
    vs = [index[v] for _, v in pairs]
    return _core.check_matching(G.core, us, vs)


@implements("is_matching", _MATCHING, declined="directed graphs; multigraphs; malformed matchings.")
def is_matching(G, matching):
    return _check_matching(G, matching)[0]


@implements("is_maximal_matching", _MATCHING, declined="directed graphs; multigraphs; malformed matchings.")
def is_maximal_matching(G, matching):
    return _check_matching(G, matching)[1]


@implements("is_perfect_matching", _MATCHING, declined="directed graphs; multigraphs; malformed matchings.")
def is_perfect_matching(G, matching):
    return _check_matching(G, matching)[2]


# ---- Colouring ---------------------------------------------------------------

# Only the strategies whose order is fixed by the degrees (and the node
# order for ties): the others pick from sets, or at random, in NetworkX.
_STRATEGIES = {
    "largest_first": 0,
    "saturation_largest_first": 6,
    "DSATUR": 6,
}


def _check_strategy(a):
    s = a.get("strategy")
    if not isinstance(s, str) or s not in _STRATEGIES:
        return "strategy is not one of the deterministic named strategies"
    if a.get("interchange"):
        return "interchange"
    return None


@implements(
    "greedy_color",
    chk.all_of(_SIMPLE_NO_LOOPS, _check_strategy),
    declined="multigraphs; self-loops; interchange; strategies other than largest_first and saturation_largest_first.",
)
def greedy_color(G, strategy="largest_first", interchange=False):
    return _node_dict(G, _core.greedy_color(G.core, _STRATEGIES[strategy]))
