"""Spanning trees, directed acyclic graphs, dominators, and cycles.

Where NetworkX's answer depends on the order of its adjacency dictionaries
(topological orders, the cycle basis, the depth-first search for a cycle),
the engine reproduces it for a graph whose adjacency order is the one a
rebuild from the edge list would give (``BackendGraph.input_order``); other
graphs are left to NetworkX. Where it depends on the iteration order of a
set (all simple cycles, the lowest common ancestor in a DAG with several of
them), the calls are left to NetworkX.
"""

from __future__ import annotations

import math
from itertools import combinations_with_replacement

import networkx as nx

from . import _checks as chk
from . import _core
from ._algorithms import ALL, ANY, SIMPLE, _hashable
from ._flow import _number, _values
from ._graph import implements, threads
from ._paths import doubles, int64s

KRUSKAL = 0
INPUT_ORDER = 1

_UNDIRECTED_SIMPLE = chk.all_of(SIMPLE, chk.undirected)
_DIRECTED_SIMPLE = chk.all_of(SIMPLE, chk.directed)


# ---- Spanning trees ---------------------------------------------------------------


def _check_mst(a):
    G = chk.graph_of(a)
    if G.is_directed():
        return "NetworkX does not define this for directed graphs"
    if G.is_multigraph():
        return "multigraphs are not supported"
    if a.get("algorithm") != "kruskal":
        return "only Kruskal's algorithm orders ties as the engine does"
    w = a.get("weight")
    if w is not None and not isinstance(w, str):
        return "weight must be an attribute name"
    return None


def _forest(G, weight, maximum):
    """The edges (u, v, data) of NetworkX's Kruskal forest, in its order."""
    if len(G) == 0:
        return []
    core = G.weighted_core(weight, 1)
    if weight is not None and not core.wfinite:
        raise NotImplementedError("weights must be finite numbers")
    edges, _ = _core.spanning_tree(core, int(weight is not None), int(maximum), KRUSKAL)
    keys = _core.edge_keys(G.core)
    adj = G.graph._adj
    return [(keys[i][0], keys[i][1], adj[keys[i][0]][keys[i][1]]) for i in edges]


def _edge_tuples(edges, data):
    if data:
        return iter(edges)
    return iter([(u, v) for u, v, _ in edges])


@implements("minimum_spanning_edges", _check_mst, declined="directed graphs; multigraphs; algorithms other than Kruskal; non-finite weights.")
def minimum_spanning_edges(G, algorithm="kruskal", weight="weight", keys=True, data=True, ignore_nan=False):
    return _edge_tuples(_forest(G, weight, False), data)


@implements("maximum_spanning_edges", _check_mst, declined="directed graphs; multigraphs; algorithms other than Kruskal; non-finite weights.")
def maximum_spanning_edges(G, algorithm="kruskal", weight="weight", keys=True, data=True, ignore_nan=False):
    return _edge_tuples(_forest(G, weight, True), data)


def _tree(G, edges):
    T = G.graph.__class__()
    T.graph.update(G.graph.graph)
    T.add_nodes_from(G.graph.nodes.items())
    T.add_edges_from(edges)
    return T


@implements("minimum_spanning_tree", _check_mst, declined="directed graphs; multigraphs; algorithms other than Kruskal; non-finite weights.")
def minimum_spanning_tree(G, weight="weight", algorithm="kruskal", ignore_nan=False):
    return _tree(G, _forest(G, weight, False))


@implements("maximum_spanning_tree", _check_mst, declined="directed graphs; multigraphs; algorithms other than Kruskal; non-finite weights.")
def maximum_spanning_tree(G, weight="weight", algorithm="kruskal", ignore_nan=False):
    return _tree(G, _forest(G, weight, True))


def _check_kruskal(a):
    reason = _check_mst({**a, "algorithm": "kruskal"})
    if reason:
        return reason
    if a.get("partition") is not None:
        return "partition"
    if not isinstance(a.get("minimum"), bool):
        return "minimum must be a bool"
    return None


@implements("kruskal_mst_edges", _check_kruskal, declined="directed graphs; multigraphs; partition; non-finite weights.")
def kruskal_mst_edges(G, minimum, weight="weight", keys=True, data=True, ignore_nan=False, partition=None):
    return _edge_tuples(_forest(G, weight, not minimum), data)


# ---- Topological orders -------------------------------------------------------------------


def _ordered(G):
    """Engine order code for a graph whose adjacency order the engine follows."""
    if not G.input_order:
        raise NotImplementedError("the order depends on adjacency order")
    return INPUT_ORDER


def _acyclic(G):
    if not _core.is_dag(G.core):
        raise NotImplementedError("a graph with a cycle is left to NetworkX")


def _generations(G):
    _acyclic(G)
    nodes, gen = _core.topological_sort(G.core, _ordered(G))
    names = G.nodes
    out = []
    last = -1
    for v in nodes:
        g = gen[v]
        if g != last:
            out.append([])
            last = g
        out[-1].append(names[v])
    return out


@implements("topological_generations", _DIRECTED_SIMPLE, declined="undirected graphs; multigraphs; graphs with a cycle; graphs whose adjacency is not in input order.")
def topological_generations(G):
    return iter(_generations(G))


@implements("topological_sort", _DIRECTED_SIMPLE, declined="undirected graphs; multigraphs; graphs with a cycle; graphs whose adjacency is not in input order.")
def topological_sort(G):
    return iter([v for gen in _generations(G) for v in gen])


def _ascending(G):
    names = list(G.graph)
    try:
        return all(a < b for a, b in zip(names, names[1:]))
    except TypeError:
        return False


@implements("lexicographical_topological_sort", chk.all_of(_DIRECTED_SIMPLE, chk.is_none("key")),
            declined="undirected graphs; multigraphs; key; graphs with a cycle; nodes not in ascending order.")
def lexicographical_topological_sort(G, key=None):
    if not _ascending(G):
        raise NotImplementedError("the engine breaks ties by node position")
    _acyclic(G)
    names = G.nodes
    return iter([names[i] for i in _core.lex_topological_sort(G.core)])


_MAX_SORTS = 100_000


@implements("all_topological_sorts", _DIRECTED_SIMPLE,
            declined="undirected graphs; multigraphs; graphs with a cycle; graphs whose adjacency is not in input order; more than 100000 orders.")
def all_topological_sorts(G):
    _acyclic(G)
    sorts = _core.all_topological_sorts(G.core, _ordered(G), _MAX_SORTS + 1)
    if len(sorts) > _MAX_SORTS:
        raise NotImplementedError("too many orders to hold at once")
    names = G.nodes
    return iter([[names[i] for i in s] for s, _ in sorts])


@implements("is_directed_acyclic_graph", chk.all_of(chk.no_multigraph), declined="multigraphs.")
def is_directed_acyclic_graph(G):
    return G.is_directed() and _core.is_dag(G.core)


@implements("has_cycle", _DIRECTED_SIMPLE, declined="undirected graphs; multigraphs.")
def has_cycle(G):
    return not _core.is_dag(G.core)


# ---- Closure and reduction -------------------------------------------------------------------


@implements("transitive_closure", chk.all_of(_DIRECTED_SIMPLE, chk.equals("reflexive", None, True, False)),
            declined="undirected graphs; multigraphs; reflexive other than None, True and False.")
def transitive_closure(G, reflexive=False):
    code = 1 if reflexive is None else 2 if reflexive is True else 0
    lists = _core.transitive_closure(G.core, code)
    names = G.nodes
    TC = G.graph.copy()
    adj = TC._adj
    for u, (targets, _) in zip(names, lists):
        have = adj[u]
        TC.add_edges_from((u, names[v]) for v in targets if names[v] not in have)
    return TC


@implements("transitive_reduction", _DIRECTED_SIMPLE, declined="undirected graphs; multigraphs; graphs with a cycle.")
def transitive_reduction(G):
    _acyclic(G)
    names = G.nodes
    TR = nx.DiGraph()
    TR.add_nodes_from(names)
    for u, (kept, _) in zip(names, _core.transitive_reduction(G.core)):
        TR.add_edges_from((u, names[v]) for v in kept)
    return TR


def _check_dag_weight(a):
    w = a.get("weight")
    if w is None:
        if a.get("default_weight") != 1 or isinstance(a.get("default_weight"), bool):
            return "without a weight attribute every edge has default_weight"
        return None
    if not isinstance(w, str):
        return "weight must be an attribute name"
    d = a.get("default_weight")
    if not isinstance(d, (int, float)) or isinstance(d, bool) or not math.isfinite(d):
        return "default_weight must be a finite number"
    return None


@implements("dag_longest_path_length", chk.all_of(_DIRECTED_SIMPLE, _check_dag_weight),
            declined="undirected graphs; multigraphs; graphs with a cycle; weight functions; non-finite weights.")
def dag_longest_path_length(G, weight="weight", default_weight=1):
    _acyclic(G)
    if weight is None:
        _, length = _core.dag_longest_path(G.core, 0, _ordered(G))
        return int(length)
    values, integral = _values(G, weight, default_weight)
    if not all(math.isfinite(x) for x in values):
        raise NotImplementedError("weights must be finite numbers")
    core = G.weighted_core(weight, default_weight)
    _, length = _core.dag_longest_path(core, 1, 0)
    return _number(length, integral)


@implements("dag_longest_path", chk.all_of(_DIRECTED_SIMPLE, _check_dag_weight, chk.is_none("topo_order")),
            declined="undirected graphs; multigraphs; graphs with a cycle; topo_order; non-finite weights.")
def dag_longest_path(G, weight="weight", default_weight=1, topo_order=None):
    _acyclic(G)
    order = _ordered(G)
    values, _ = _values(G, weight, default_weight) if weight is not None else ([], True)
    if not all(math.isfinite(x) for x in values):
        raise NotImplementedError("weights must be finite numbers")
    core = G.weighted_core(weight, default_weight)
    path, _ = _core.dag_longest_path(core, int(weight is not None), order)
    names = G.nodes
    return [names[i] for i in path]


# ---- Dominators and common ancestors ------------------------------------------------------------


def _start_in_graph(a):
    G = chk.graph_of(a)
    try:
        if a.get("start") in G:
            return None
    except TypeError:
        pass
    return "start is not a node of the graph"


@implements("immediate_dominators", chk.all_of(_DIRECTED_SIMPLE, _start_in_graph),
            declined="undirected graphs; multigraphs; a start that is not a node.")
def immediate_dominators(G, start):
    idom = _core.dominators(G.core, G.index[start])
    names = G.nodes
    s = G.index[start]
    return {names[v]: names[p] for v, p in enumerate(idom) if p >= 0 and v != s}


@implements("dominance_frontiers", chk.all_of(_DIRECTED_SIMPLE, _start_in_graph),
            declined="undirected graphs; multigraphs; a start that is not a node.")
def dominance_frontiers(G, start):
    rows = _core.dominance_frontiers(G.core, G.index[start])
    names = G.nodes
    return {names[u]: {names[v] for v in items} for u, (items, reach) in enumerate(rows) if reach == 1.0}


def _in_forest(G):
    pred = G.graph._pred
    return all(len(p) <= 1 for p in pred.values())


def _check_lca_graph(a):
    reason = _DIRECTED_SIMPLE(a)
    if reason:
        return reason
    G = chk.graph_of(a)
    if not _in_forest(_Wrap(G)):
        return "a node with several parents may have several lowest common ancestors"
    return None


class _Wrap:
    def __init__(self, graph):
        self.graph = graph


@implements("lowest_common_ancestor", chk.all_of(_check_lca_graph, chk.nodes_in_graph("node1", "node2")),
            declined="undirected graphs; multigraphs; a node with several parents; graphs with a cycle.")
def lowest_common_ancestor(G, node1, node2, default=None):
    _acyclic(G)
    out = _core.lca(G.core, [G.index[node1]], [G.index[node2]], -2, 0)
    return G.nodes[out[0]] if out[0] >= 0 else default


def _check_pairs(a):
    pairs = a.get("pairs")
    if pairs is None:
        return None
    G = chk.graph_of(a)
    from collections.abc import Collection

    if not isinstance(pairs, Collection) or isinstance(pairs, (str, bytes)):
        return "pairs must be a collection of pairs"
    try:
        for p in pairs:
            if not isinstance(p, (tuple, list)) or len(p) != 2 or p[0] not in G or p[1] not in G:
                return "pairs must be pairs of nodes of the graph"
    except TypeError:
        return "pairs must be pairs of nodes of the graph"
    return None


@implements("all_pairs_lowest_common_ancestor", chk.all_of(_check_lca_graph, _check_pairs),
            declined="undirected graphs; multigraphs; a node with several parents; graphs with a cycle.")
def all_pairs_lowest_common_ancestor(G, pairs=None):
    _acyclic(G)
    names = G.nodes
    index = G.index
    if pairs is None:
        pair_list = list(combinations_with_replacement(names, 2))
    else:
        pair_list = list(dict.fromkeys(tuple(p) for p in pairs))
    out = []
    for lo in range(0, len(pair_list), 1 << 16):
        chunk = pair_list[lo:lo + (1 << 16)]
        res = _core.lca(G.core, [index[v] for v, _ in chunk], [index[w] for _, w in chunk], -2, 0)
        out.extend(((v, w), names[r]) for (v, w), r in zip(chunk, res) if r >= 0)
    return iter(out)


# ---- Cycles ---------------------------------------------------------------------------------------


def _root_ok(a):
    r = a.get("root")
    if r is None:
        return None
    G = chk.graph_of(a)
    try:
        if r in G:
            return None
    except TypeError:
        pass
    return "root is not a node of the graph"


@implements("cycle_basis", chk.all_of(_UNDIRECTED_SIMPLE, _root_ok),
            declined="directed graphs; multigraphs; graphs whose adjacency is not in input order.")
def cycle_basis(G, root=None):
    order = _ordered(G)
    lists = _core.cycle_basis(G.core, order, -1 if root is None else G.index[root])
    names = G.nodes
    return [[names[i] for i in c] for c, _ in lists]


@implements("girth", chk.all_of(_UNDIRECTED_SIMPLE, chk.no_selfloops), declined="directed graphs; multigraphs; self-loops.")
def girth(G):
    g = _core.girth(G.core)
    return g if math.isinf(g) else int(g)


# ---- Rich club -------------------------------------------------------------------------------------


@implements("rich_club_coefficient", chk.all_of(_UNDIRECTED_SIMPLE, chk.no_selfloops, chk.equals("normalized", False)),
            declined="directed graphs; multigraphs; self-loops; normalized=True (it rewires with NetworkX's random stream).")
def rich_club_coefficient(G, normalized=True, Q=100, seed=None, *, n_samples=1):
    if G.core.m == 0:
        raise NotImplementedError("no edges")
    return dict(enumerate(doubles(_core.rich_club(G.core, threads()))))


def _check_find_cycle(a):
    reason = chk.no_multigraph(a)
    if reason:
        return reason
    if a.get("orientation") not in (None, "original", "reverse", "ignore"):
        return "orientation is not one of the four NetworkX accepts"
    s = a.get("source")
    if s is not None:
        G = chk.graph_of(a)
        try:
            if s not in G:
                return "source is not a single node of the graph"
        except TypeError:
            return "source is not a single node of the graph"
    return None


_ORIENTATIONS = {None: 0, "original": 0, "reverse": 1, "ignore": 2}


@implements("find_cycle", _check_find_cycle,
            declined="multigraphs; a source that is not a single node; graphs whose adjacency is not in input order.")
def find_cycle(G, source=None, orientation=None):
    order = _ordered(G)
    directed = G.core.directed
    code = _ORIENTATIONS[orientation] if directed else 0
    edges, flipped = _core.find_cycle(G.core, -1 if source is None else G.index[source], code, order)
    if not edges:
        raise NotImplementedError("no cycle is left to NetworkX to report")
    keys = _core.edge_keys(G.core)
    out = []
    if directed:
        for i, rev in zip(edges, flipped):
            u, v = keys[i]
            tup = (u, v)
            flag = "reverse" if rev else "forward"
            out.append(tup if orientation is None else (*tup, flag))
        return out
    # An undirected edge has no direction in the engine: the walk's own
    # orientation follows from the chain, each edge leaving where the last
    # one arrived and the first one leaving the end that the second does not
    # touch.
    pairs = [keys[i] for i in edges]
    head = None
    if len(pairs) > 1:
        head = next(x for x in pairs[0] if x in pairs[1])
    for j, (u, v) in enumerate(pairs):
        if j == 0:
            tail, head = (u, u) if len(pairs) == 1 else ((u, v) if head == v else (v, u))
        else:
            tail, head = head, (v if head == u else u)
        out.append((tail, head) if orientation is None else (tail, head, "forward"))
    return out
