"""Shortest paths on the engine: lengths with positive or negative weights,
negative-cycle tests, dense all-pairs lengths, and the functions that return
paths or predecessor lists.

Distances are exact wherever NetworkX's are. Which of several equally short
paths NetworkX returns depends on the order in its adjacency dictionaries
(its breadth-first and Dijkstra searches keep the first predecessor they
meet), which the engine does not follow, so the functions that return paths
or predecessor lists run on the engine only for sources from which every
shortest path is unique (the engine reports every predecessor on a shortest
path, comparing lengths exactly as NetworkX does); otherwise they are left
to NetworkX.
"""

from __future__ import annotations

import math
from collections import defaultdict

import networkx as nx

from . import _checks as chk
from . import _core
from ._algorithms import (
    ALL, ANY, IN, OUT, SIMPLE, _check_cutoff, _hashable, _nan_to_zero, _node_dict,
    _nonnegative_weights, _positive_weights, _source_blocks, _weighted,
)
from ._graph import implements, threads

INF = math.inf
BELLMAN_FORD, GOLDBERG_RADZIK = 0, 1
EXACT = 1  # tie rule: lengths are equal only when equal as floating-point numbers


# ---- Helpers --------------------------------------------------------------


def doubles(b):
    return memoryview(b).cast("d").tolist()


def int64s(b):
    return memoryview(b).cast("q").tolist()


def _graph(G, weight, mode=OUT, collapse=False, lengths=None):
    """A path graph over the converted edge attribute (unit lengths when
    weight is None)."""
    w = _weighted(G, weight) if lengths is None else 0
    return _core.path_graph(G.core, mode, w, 0, lengths, int(collapse))


def _num(x, weighted):
    """NetworkX's unweighted lengths are integers; weighted ones are floats."""
    return x if weighted else int(x)


def _cut(cutoff):
    return INF if cutoff is None else float(cutoff)


def _index(G, v):
    return G.index[v]


def _tie_free(preds):
    for lst, _ in preds:
        if len(lst) > 1:
            raise NotImplementedError("tied shortest paths follow NetworkX's search order")


def _one_search(G, weight, sources, mode=OUT, cutoff=None, collapse=False):
    """Distances (list) and the single predecessor of every node (-1 for none)
    of a search from the node indices ``sources``; raises NotImplementedError if
    any node has two predecessors on shortest paths."""
    w = _weighted(G, weight)
    _positive_weights(G, w)
    pg = _graph(G, weight, mode, collapse)
    dist_b, preds = _core.sp_predecessors(pg, sources, _cut(cutoff), 0, EXACT)
    _tie_free(preds)
    dist = doubles(dist_b)
    pred = [lst[0] if lst else -1 for lst, _ in preds]
    return dist, pred, w


def _order(dist):
    """Reached nodes by increasing distance (ties by index)."""
    return sorted((i for i, d in enumerate(dist) if d < INF), key=lambda i: (dist[i], i))


def _paths_from(G, dist, pred, forward):
    """Paths to every reached node (forward) or from it (backward, for
    searches from a target), keyed by label."""
    nodes = G.nodes
    paths = {}
    index_path = {}
    for i in _order(dist):
        p = pred[i]
        if p < 0:
            ip = [i]
        elif forward:
            ip = index_path[p] + [i]
        else:
            ip = [i] + index_path[p]
        index_path[i] = ip
        paths[nodes[i]] = [nodes[j] for j in ip]
    return paths


def _distances(G, dist, weighted, cutoff=None):
    nodes = G.nodes
    return {nodes[i]: _num(dist[i], weighted) for i in _order(dist)
            if cutoff is None or dist[i] <= cutoff}


# ---- Dijkstra and Bellman-Ford: lengths ------------------------------------

_LEN = chk.all_of(ANY, chk.weight_attr(), chk.numeric_or_none("cutoff"))


@implements(
    "multi_source_dijkstra_path_length",
    chk.all_of(_LEN, chk.nodes_all_in_graph("sources")),
    declined="weight functions; weights <= 0.",
)
def multi_source_dijkstra_path_length(G, sources, cutoff=None, weight="weight"):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    idx = [_index(G, s) for s in sources]
    pg = _graph(G, weight, OUT, True)
    tree = _core.sp_tree(pg, idx, _cut(cutoff), -1, 0)
    return _distances(G, doubles(tree[0]), w, cutoff)


@implements(
    "dijkstra_path_length",
    chk.all_of(ANY, chk.weight_attr(), chk.nodes_in_graph("source", "target")),
    declined="weight functions; weights <= 0; unreachable targets.",
)
def dijkstra_path_length(G, source, target, weight="weight"):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    if source == target:
        return 0
    pg = _graph(G, weight, OUT, True)
    t = _index(G, target)
    tree = _core.sp_tree(pg, [_index(G, source)], INF, t, 0)
    d = doubles(tree[0])[t]
    if d == INF:
        raise NotImplementedError("no path; NetworkX raises")
    return _num(d, w)


def _bellman_ford_dist(G, weight, source):
    """Distances from one source with weights of any sign (NotImplementedError
    on a negative cycle, which NetworkX reports as an exception)."""
    w = _weighted(G, weight)
    pg = _graph(G, weight, OUT, True)
    dist_b, _, _, cycle = _core.sp_negative(pg, [_index(G, source)], BELLMAN_FORD, EXACT)
    if cycle:
        raise NotImplementedError("negative cycle; NetworkX raises")
    return doubles(dist_b), w


@implements(
    "single_source_bellman_ford_path_length",
    chk.all_of(ANY, chk.weight_attr(), chk.nodes_in_graph("source")),
    declined="weight functions; negative cycles.",
)
def single_source_bellman_ford_path_length(G, source, weight="weight"):
    dist, w = _bellman_ford_dist(G, weight, source)
    return _distances(G, dist, w)


@implements(
    "bellman_ford_path_length",
    chk.all_of(ANY, chk.weight_attr(), chk.nodes_in_graph("source", "target")),
    declined="weight functions; negative cycles; unreachable targets.",
)
def bellman_ford_path_length(G, source, target, weight="weight"):
    if source == target:
        return 0
    dist, w = _bellman_ford_dist(G, weight, source)
    d = dist[_index(G, target)]
    if d == INF:
        raise NotImplementedError("no path; NetworkX raises")
    return _num(d, w)


@implements(
    "all_pairs_bellman_ford_path_length",
    chk.all_of(ANY, chk.weight_attr()),
    declined="weight functions; negative cycles.",
)
def all_pairs_bellman_ford_path_length(G, weight="weight"):
    w = _weighted(G, weight)
    n = len(G)
    pg = _graph(G, weight, OUT, True)
    if not w or G.core.wmin >= 0:
        def rows():
            for block in _source_blocks(n):
                flat = doubles(_core.sp_all_pairs(pg, block, INF, threads(), 0)[0])
                yield block, flat
    else:
        dist_b, _, cycle = _core.sp_johnson(pg, threads())
        if cycle:
            raise NotImplementedError("negative cycle; NetworkX raises")
        flat_all = doubles(dist_b)

        def rows():
            yield list(range(n)), flat_all
    return _bf_rows(G, rows(), w)


def _bf_rows(G, blocks, w):
    nodes = G.nodes
    for block, flat in blocks:
        k = len(block)
        # Column-major: source row r, node v at v * k + r.
        for r, i in enumerate(block):
            yield nodes[i], _distances(G, flat[r::k], w)


@implements(
    "negative_edge_cycle",
    chk.all_of(ANY, chk.weight_attr()),
    declined="weight functions.",
)
def negative_edge_cycle(G, weight="weight", heuristic=True):
    w = _weighted(G, weight)
    if G.core.m == 0:
        return False
    pg = _graph(G, weight, OUT, True)
    return bool(_core.sp_negative_cycle(pg, -1))


# ---- Dense all-pairs ---------------------------------------------------------


def _floyd(G, weight):
    w = _weighted(G, weight)
    pg = _graph(G, weight, OUT, True)
    dist_b, pred_b, neg = _core.sp_floyd_warshall(pg, threads())
    if neg:
        raise NotImplementedError("negative cycle; NetworkX raises")
    n = len(G)
    return doubles(dist_b), int64s(pred_b), n


def _rows(flat, n):
    """Column-major n-by-n values as rows: element (u, v) is flat[v * n + u]."""
    return [flat[u::n] for u in range(n)]


@implements(
    "floyd_warshall",
    chk.all_of(ANY, chk.weight_attr()),
    declined="weight functions; negative cycles.",
)
def floyd_warshall(G, weight="weight"):
    dist, _, n = _floyd(G, weight)
    nodes = G.nodes
    out = {}
    for u, row in zip(nodes, _rows(dist, n)):
        d = defaultdict(lambda: INF)
        d.update(zip(nodes, row))
        out[u] = d
    return out


@implements(
    "floyd_warshall_numpy",
    chk.all_of(ANY, chk.weight_attr()),
    declined="weight functions; negative cycles; nodelist other than None.",
)
def floyd_warshall_numpy(G, nodelist=None, weight="weight"):
    import numpy as np

    if nodelist is not None:
        raise NotImplementedError("nodelist is left to NetworkX")
    dist, _, n = _floyd(G, weight)
    return np.array(dist, dtype=float).reshape(n, n).T.copy()


@implements(
    "floyd_warshall_predecessor_and_distance",
    chk.all_of(ANY, chk.weight_attr()),
    declined="weight functions; negative cycles.",
)
def floyd_warshall_predecessor_and_distance(G, weight="weight"):
    dist, pred, n = _floyd(G, weight)
    nodes = G.nodes
    dist_rows = _rows(dist, n)
    pred_rows = _rows(pred, n)
    selfloop = set(v for v, nbrs in G.graph._adj.items() if v in nbrs)
    out_pred = {}
    for i, u in enumerate(nodes):
        d = {}
        for j, p in enumerate(pred_rows[i]):
            if p >= 0 and j != i:
                d[nodes[j]] = nodes[p]
        if u in selfloop:
            d[u] = u  # NetworkX leaves the self-loop's predecessor in place
        if d:
            out_pred[u] = d
    out_dist = {}
    for u, row in zip(nodes, dist_rows):
        dd = defaultdict(lambda: INF)
        dd.update(zip(nodes, row))
        out_dist[u] = dd
    return out_pred, out_dist


# ---- Functions that return paths (unique shortest paths only) -----------------


@implements(
    "single_source_shortest_path",
    chk.all_of(SIMPLE, chk.nodes_in_graph("source"), chk.numeric_or_none("cutoff")),
    declined="multigraphs; sources with tied shortest paths.",
)
def single_source_shortest_path(G, source, cutoff=None):
    dist, pred, _ = _one_search(G, None, [_index(G, source)], OUT, cutoff)
    paths = _paths_from(G, dist, pred, True)
    return {v: p for v, p in paths.items() if cutoff is None or len(p) - 1 <= cutoff}


@implements(
    "single_target_shortest_path",
    chk.all_of(SIMPLE, chk.nodes_in_graph("target"), chk.numeric_or_none("cutoff")),
    declined="multigraphs; targets with tied shortest paths.",
)
def single_target_shortest_path(G, target, cutoff=None):
    dist, pred, _ = _one_search(G, None, [_index(G, target)], IN, cutoff)
    paths = _paths_from(G, dist, pred, False)
    return {v: p for v, p in paths.items() if cutoff is None or len(p) - 1 <= cutoff}


@implements(
    "all_pairs_shortest_path",
    chk.all_of(SIMPLE, chk.numeric_or_none("cutoff")),
    declined="multigraphs; graphs with tied shortest paths.",
)
def all_pairs_shortest_path(G, cutoff=None):
    # A tie anywhere declines the call before the first result is produced.
    pg = _graph(G, None, OUT)
    for i in range(len(G)):
        _, preds = _core.sp_predecessors(pg, [i], _cut(cutoff), 0, EXACT)
        _tie_free(preds)
    return ((v, single_source_shortest_path(G, v, cutoff)) for v in G.nodes)


@implements(
    "predecessor",
    chk.all_of(SIMPLE, chk.nodes_in_graph("source"), chk.numeric_or_none("cutoff"),
               chk.is_none("return_seen", "target")),
    declined="multigraphs; return_seen; target; sources with tied shortest paths.",
)
def predecessor(G, source, target=None, cutoff=None, return_seen=None):
    # NetworkX stops after level `cutoff` and treats 0 as no cutoff.
    lim = cutoff if cutoff else None
    dist, pred, _ = _one_search(G, None, [_index(G, source)], OUT, lim)
    nodes = G.nodes
    return {nodes[i]: ([nodes[pred[i]]] if pred[i] >= 0 else []) for i in _order(dist)}


def _bidirectional_unique(G, weight, source, target):
    """The unique shortest path from source to target, or NotImplementedError."""
    w = _weighted(G, weight)
    _positive_weights(G, w)
    s, t = _index(G, source), _index(G, target)
    pg = _graph(G, weight, OUT)
    _, preds = _core.sp_predecessors(pg, [s], INF, 0, EXACT)
    # Only the part of the tree on shortest paths to the target needs to be unique.
    path = [t]
    v = t
    while v != s:
        lst = preds[v][0]
        if not lst:
            raise NotImplementedError("no path; NetworkX raises")
        if len(lst) > 1:
            raise NotImplementedError("tied shortest paths follow NetworkX's search order")
        v = lst[0]
        path.append(v)
    path.reverse()
    return [G.nodes[i] for i in path], w


@implements(
    "bidirectional_shortest_path",
    chk.all_of(SIMPLE, chk.nodes_in_graph("source", "target")),
    declined="multigraphs; pairs with tied shortest paths or no path.",
)
def bidirectional_shortest_path(G, source, target):
    if source == target:
        return [source]
    return _bidirectional_unique(G, None, source, target)[0]


@implements(
    "bidirectional_dijkstra",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.nodes_in_graph("source", "target")),
    declined="multigraphs; weight functions; weights <= 0; pairs with tied shortest paths or no path.",
)
def bidirectional_dijkstra(G, source, target, weight="weight"):
    if source == target:
        return 0, [source]
    path, w = _bidirectional_unique(G, weight, source, target)
    pg = _graph(G, weight, OUT)
    t = _index(G, target)
    tree = _core.sp_tree(pg, [_index(G, source)], INF, t, 0)
    return _num(doubles(tree[0])[t], w), path


@implements(
    "shortest_path",
    chk.all_of(SIMPLE, chk.weight_attr(),
               chk.equals("method", "dijkstra"),
               lambda a: None if (a["source"] is None or (_hashable(a["source"]) and a["source"] in chk.graph_of(a)))
               and (a["target"] is None or (_hashable(a["target"]) and a["target"] in chk.graph_of(a)))
               else "source and target must be nodes of the graph"),
    declined="multigraphs; method='bellman-ford'; weight functions; weights <= 0; pairs with tied shortest paths or no path.",
)
def shortest_path(G, source=None, target=None, weight=None, method="dijkstra"):
    if source is None and target is None:
        return all_pairs_dijkstra_path(G, weight=weight) if weight is not None \
            else all_pairs_shortest_path(G)
    if source is None:
        return single_target_shortest_path(G, target) if weight is None \
            else _single_target_dijkstra_path(G, target, weight)
    if target is None:
        return single_source_shortest_path(G, source) if weight is None \
            else single_source_dijkstra_path(G, source, weight=weight)
    if source == target:
        return [source]
    return _bidirectional_unique(G, weight, source, target)[0]


def _single_target_dijkstra_path(G, target, weight):
    dist, pred, _ = _one_search(G, weight, [_index(G, target)], IN)
    return _paths_from(G, dist, pred, False)


@implements(
    "single_source_dijkstra_path",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.nodes_in_graph("source"),
               chk.numeric_or_none("cutoff")),
    declined="multigraphs; weight functions; weights <= 0; sources with tied shortest paths.",
)
def single_source_dijkstra_path(G, source, cutoff=None, weight="weight"):
    dist, pred, _ = _one_search(G, weight, [_index(G, source)], OUT, cutoff)
    paths = _paths_from(G, dist, pred, True)
    return {v: p for v, p in paths.items() if cutoff is None or dist[G.index[v]] <= cutoff}


@implements(
    "single_source_dijkstra",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.nodes_in_graph("source"),
               chk.numeric_or_none("cutoff"),
               lambda a: None if a["target"] is None or (_hashable(a["target"]) and a["target"] in chk.graph_of(a))
               else "target is not a node of the graph"),
    declined="multigraphs; weight functions; weights <= 0; sources with tied shortest paths; unreachable targets.",
)
def single_source_dijkstra(G, source, target=None, cutoff=None, weight="weight"):
    w = _weighted(G, weight)
    if target is not None and source == target:
        return 0, [source]
    dist, pred, _ = _one_search(G, weight, [_index(G, source)], OUT, cutoff)
    paths = _paths_from(G, dist, pred, True)
    if target is None:
        d = _distances(G, dist, w, cutoff)
        return d, {v: paths[v] for v in d}
    t = _index(G, target)
    if dist[t] == INF or (cutoff is not None and dist[t] > cutoff):
        raise NotImplementedError("no path; NetworkX raises")
    return _num(dist[t], w), paths[target]


@implements(
    "dijkstra_path",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.nodes_in_graph("source", "target")),
    declined="multigraphs; weight functions; weights <= 0; pairs with tied shortest paths or no path.",
)
def dijkstra_path(G, source, target, weight="weight"):
    if source == target:
        return [source]
    return _bidirectional_unique(G, weight, source, target)[0]


@implements(
    "dijkstra_predecessor_and_distance",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.nodes_in_graph("source"),
               chk.numeric_or_none("cutoff")),
    declined="multigraphs; weight functions; weights <= 0; sources with tied shortest paths.",
)
def dijkstra_predecessor_and_distance(G, source, cutoff=None, weight="weight"):
    dist, pred, w = _one_search(G, weight, [_index(G, source)], OUT, cutoff)
    nodes = G.nodes
    order = [i for i in _order(dist) if cutoff is None or dist[i] <= cutoff]
    return ({nodes[i]: ([nodes[pred[i]]] if pred[i] >= 0 else []) for i in order},
            {nodes[i]: _num(dist[i], w) for i in order})


def _negative_tree(G, weight, source):
    """Distances and unique predecessors from one source with weights of any
    sign; NotImplementedError on a negative cycle or a tied predecessor."""
    w = _weighted(G, weight)
    if w and G.core.wzeros:
        raise NotImplementedError("zero-weight edges can tie with the source")
    pg = _graph(G, weight, OUT)
    dist_b, _, preds, cycle = _core.sp_negative(pg, [_index(G, source)], BELLMAN_FORD, EXACT)
    if cycle:
        raise NotImplementedError("negative cycle; NetworkX raises")
    _tie_free(preds)
    dist = doubles(dist_b)
    return dist, [lst[0] if lst else -1 for lst, _ in preds], w


@implements(
    "bellman_ford_predecessor_and_distance",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.nodes_in_graph("source"), chk.is_none("target")),
    declined="multigraphs; weight functions; target; negative cycles; sources with tied shortest paths.",
)
def bellman_ford_predecessor_and_distance(G, source, target=None, weight="weight",
                                          heuristic=False):
    dist, pred, w = _negative_tree(G, weight, source)
    nodes = G.nodes
    order = _order(dist)
    return ({nodes[i]: ([nodes[pred[i]]] if pred[i] >= 0 else []) for i in order},
            {nodes[i]: _num(dist[i], w) for i in order})


def _bf_paths(G, weight, source):
    dist, pred, w = _negative_tree(G, weight, source)
    index_path = {}
    paths = {}
    # Without ties, distances can be negative, so order by tree depth instead.
    children = defaultdict(list)
    roots = []
    for i, p in enumerate(pred):
        if dist[i] == INF:
            continue
        (children[p] if p >= 0 else roots).append(i)
    stack = list(roots)
    while stack:
        i = stack.pop()
        p = pred[i]
        index_path[i] = [i] if p < 0 else index_path[p] + [i]
        paths[G.nodes[i]] = [G.nodes[j] for j in index_path[i]]
        stack.extend(children[i])
    return dist, paths, w


@implements(
    "single_source_bellman_ford_path",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.nodes_in_graph("source")),
    declined="multigraphs; weight functions; negative cycles; sources with tied shortest paths.",
)
def single_source_bellman_ford_path(G, source, weight="weight"):
    return _bf_paths(G, weight, source)[1]


@implements(
    "bellman_ford_path",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.nodes_in_graph("source", "target")),
    declined="multigraphs; weight functions; negative cycles; pairs with tied shortest paths or no path.",
)
def bellman_ford_path(G, source, target, weight="weight"):
    if source == target:
        return [source]
    paths = _bf_paths(G, weight, source)[1]
    if target not in paths:
        raise NotImplementedError("no path; NetworkX raises")
    return paths[target]


@implements(
    "single_source_bellman_ford",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.nodes_in_graph("source"),
               lambda a: None if a["target"] is None or (_hashable(a["target"]) and a["target"] in chk.graph_of(a))
               else "target is not a node of the graph"),
    declined="multigraphs; weight functions; negative cycles; sources with tied shortest paths; unreachable targets.",
)
def single_source_bellman_ford(G, source, target=None, weight="weight"):
    if target is not None and source == target:
        return 0, [source]
    dist, paths, w = _bf_paths(G, weight, source)
    if target is None:
        d = {G.nodes[i]: _num(dist[i], w) for i in _order(dist)}
        return d, paths
    if target not in paths:
        raise NotImplementedError("no path; NetworkX raises")
    return _num(dist[G.index[target]], w), paths[target]


@implements(
    "all_pairs_dijkstra_path",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.numeric_or_none("cutoff")),
    declined="multigraphs; weight functions; weights <= 0; graphs with tied shortest paths.",
)
def all_pairs_dijkstra_path(G, cutoff=None, weight="weight"):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    pg = _graph(G, weight, OUT)
    for i in range(len(G)):
        _, preds = _core.sp_predecessors(pg, [i], _cut(cutoff), 0, EXACT)
        _tie_free(preds)
    return ((v, single_source_dijkstra_path(G, v, cutoff, weight)) for v in G.nodes)


@implements(
    "all_pairs_dijkstra",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.numeric_or_none("cutoff")),
    declined="multigraphs; weight functions; weights <= 0; graphs with tied shortest paths.",
)
def all_pairs_dijkstra(G, cutoff=None, weight="weight"):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    pg = _graph(G, weight, OUT)
    for i in range(len(G)):
        _, preds = _core.sp_predecessors(pg, [i], _cut(cutoff), 0, EXACT)
        _tie_free(preds)
    return ((v, single_source_dijkstra(G, v, None, cutoff, weight)) for v in G.nodes)


@implements(
    "all_pairs_bellman_ford_path",
    chk.all_of(SIMPLE, chk.weight_attr()),
    declined="multigraphs; weight functions; negative cycles; graphs with tied shortest paths.",
)
def all_pairs_bellman_ford_path(G, weight="weight"):
    _weighted(G, weight)
    # Computed up front, so that a tie or a negative cycle declines the call.
    results = [(v, single_source_bellman_ford_path(G, v, weight)) for v in G.nodes]
    return iter(results)
