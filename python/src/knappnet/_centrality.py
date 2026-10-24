"""Centrality variants on the engine: betweenness for sources and targets,
percolation, load, current-flow, communicability, resistance, dispersion,
VoteRank, Laplacian centrality, trophic levels, and reaching centrality.

The engine runs each with the NetworkX convention flags that
crates/capi/include/netcapi.h documents ("Centrality: path variants,
current flow, walks"); validation/bindings/centrality_extra.py and the
manifest in validation/manifest.d/centrality_extra.yaml show which call of
NetworkX each engine call reproduces. Calls on which NetworkX raises an
exception (disconnected graphs for the current-flow measures, for
example) are left to NetworkX.
"""

from __future__ import annotations

import math

import networkx as nx

from . import _checks as chk
from . import _core
from ._algorithms import (
    ALL, ANY, IN, OUT, SIMPLE, _connected, _hashable, _node_dict, _nonnegative_weights,
    _positive_weights, _weighted,
)
from ._graph import implements, threads
from ._paths import doubles, int64s

DISTANCE, STRENGTH = 0, 1
UNORDERED, NETWORKX, NETWORKX_SUBSET = 0, 1, 2
DISTINCT, COLLAPSE = 0, 1
PROPORTIONAL, NETWORKX_SPLIT = 0, 1
EXACT = 1


def _rows(flat, n):
    """Column-major n-by-n values as a list of rows."""
    return [flat[u::n] for u in range(n)]


def _edge_dict(G, values):
    return dict(zip(_core.edge_keys(G.core), values))


def _members(G, items, what, dedupe):
    """Indices of a collection of nodes; NotImplementedError unless all are nodes."""
    try:
        items = list(items)
    except TypeError as err:
        raise NotImplementedError(f"{what} must be a collection") from err
    index = G.index
    try:
        idx = [index[v] for v in items]
    except (KeyError, TypeError) as err:
        raise NotImplementedError(f"{what} must be nodes of the graph") from err
    if dedupe:
        idx = sorted(set(idx))
    elif len(set(idx)) != len(idx):
        raise NotImplementedError(f"repeated {what} are left to NetworkX")
    return idx


def _hashable_items(a, *params):
    for p in params:
        try:
            items = list(a.get(p))
            for v in items:
                hash(v)
        except TypeError:
            return f"{p} must be a collection of nodes"
    return None


_WEIGHT = chk.weight_attr()


def _orderable(G):
    """NetworkX sorts (length, node) pairs in load centrality, which raises
    for labels that cannot be compared."""
    try:
        sorted(G.nodes)
    except TypeError as err:
        raise NotImplementedError("node labels cannot be ordered; NetworkX raises") from err


def _rcm_edge_dict(G, values):
    """Edge values keyed as NetworkX's current-flow functions key them: each
    edge as (u, v) with u first in the reverse Cuthill-McKee ordering."""
    order = list(nx.utils.reverse_cuthill_mckee_ordering(G.graph))
    pos = {v: i for i, v in enumerate(order)}
    out = {}
    for (u, v), x in zip(_core.edge_keys(G.core), values):
        out[(u, v) if pos[u] < pos[v] else (v, u)] = x
    return out


# ---- Betweenness for sources and targets ------------------------------------


@implements(
    "betweenness_centrality_subset",
    chk.all_of(ANY, _WEIGHT, lambda a: _hashable_items(a, "sources", "targets")),
    declined="weight functions; weights <= 0; repeated sources; sources or targets that are not nodes.",
)
def betweenness_centrality_subset(G, sources, targets, normalized=False, weight=None):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    s = _members(G, sources, "sources", False)
    t = _members(G, targets, "targets", True)
    nodes, _ = _core.betweenness_subset(G.core, s, t, w, DISTANCE, 0, int(bool(normalized)),
                                        NETWORKX_SUBSET, COLLAPSE, PROPORTIONAL, EXACT,
                                        threads(), 0)
    return _node_dict(G, doubles(nodes))


@implements(
    "edge_betweenness_centrality_subset",
    chk.all_of(SIMPLE, _WEIGHT, lambda a: _hashable_items(a, "sources", "targets")),
    declined="multigraphs; weight functions; weights <= 0; repeated sources; sources or targets that are not nodes.",
)
def edge_betweenness_centrality_subset(G, sources, targets, normalized=False, weight=None):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    s = _members(G, sources, "sources", False)
    t = _members(G, targets, "targets", True)
    _, edges = _core.betweenness_subset(G.core, s, t, w, DISTANCE, 0, int(bool(normalized)),
                                        NETWORKX_SUBSET, COLLAPSE, NETWORKX_SPLIT, EXACT,
                                        threads(), 1)
    return _edge_dict(G, doubles(edges))


def _check_percolation(a):
    G = chk.graph_of(a)
    if len(G) < 3:
        return "graphs with fewer than three nodes are left to NetworkX"
    states = a.get("states")
    if states is not None:
        if not isinstance(states, dict) or any(v not in states for v in G):
            return "states must give a value for every node"
    return None


@implements(
    "percolation_centrality",
    chk.all_of(ANY, _WEIGHT, _check_percolation),
    declined="weight functions; weights <= 0; graphs with fewer than three nodes; states that are not complete.",
)
def percolation_centrality(G, attribute="percolation", states=None, weight=None):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    if states is None:
        data = G.graph._node
        states = {v: data[v].get(attribute, 1) for v in G.nodes}
    try:
        x = [float(states[v]) for v in G.nodes]
    except (TypeError, ValueError) as err:
        raise NotImplementedError("states must be numbers") from err
    vals = doubles(_core.percolation(G.core, x, w, DISTANCE, COLLAPSE, EXACT, threads()))
    if not all(math.isfinite(v) for v in vals):
        raise NotImplementedError("undefined values; NetworkX raises or returns them as is")
    n = len(G)
    return _node_dict(G, vals)


# ---- Load centrality ----------------------------------------------------------


def _check_cutoff_num(a):
    c = a.get("cutoff")
    if c is None or c is False:
        return None
    if isinstance(c, (int, float)) and not isinstance(c, bool):
        return None
    return "cutoff must be a number"


@implements(
    "newman_betweenness_centrality",
    chk.all_of(SIMPLE, _WEIGHT, _check_cutoff_num),
    declined="multigraphs; weight functions; weights <= 0; cutoff of 0.",
)
def newman_betweenness_centrality(G, v=None, cutoff=None, normalized=True, weight=None):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    _orderable(G)
    if cutoff is not None and cutoff is not False and not cutoff > 0:
        raise NotImplementedError("cutoff of 0 is left to NetworkX")
    c = -1.0 if not cutoff else float(cutoff)
    vals = doubles(_core.load(G.core, w, DISTANCE, int(bool(normalized)),
                              math.inf if c < 0 else c, 1, EXACT, threads()))
    result = _node_dict(G, vals)
    if v is None:
        return result
    if not _hashable(v) or v not in result:
        raise NotImplementedError("v is not a node")
    return result[v]


@implements(
    "edge_load_centrality",
    chk.all_of(SIMPLE, _check_cutoff_num),
    declined="multigraphs; cutoff of 0.",
)
def edge_load_centrality(G, cutoff=False):
    _orderable(G)
    if cutoff is not None and cutoff is not False and not cutoff > 0:
        raise NotImplementedError("cutoff of 0 is left to NetworkX")
    c = math.inf if not cutoff else float(cutoff)
    vals = doubles(_core.edge_load(G.core, c, threads()))
    out = {}
    keys = _core.edge_keys(G.core)
    # Two values per edge: (source, target) and (target, source).
    for i, (u, v) in enumerate(keys):
        out[(u, v)] = vals[2 * i]
        out[(v, u)] = vals[2 * i + 1]
    return out


# ---- Current flow ---------------------------------------------------------------

_CF = chk.all_of(SIMPLE, chk.undirected, _WEIGHT,
                 chk.equals("solver", "full", "lu", "cg"))


def _current_flow_ready(G, w):
    _positive_weights(G, w)
    if len(G) < 3 or not _connected(G) or G.core.loops:
        raise NotImplementedError("left to NetworkX")


@implements(
    "current_flow_betweenness_centrality", _CF,
    declined="directed graphs; multigraphs; weight functions; weights <= 0; disconnected graphs; graphs of fewer than three nodes; self-loops.",
)
def current_flow_betweenness_centrality(G, normalized=True, weight=None, dtype=float,
                                        solver="full"):
    w = _weighted(G, weight)
    _current_flow_ready(G, w)
    nodes, _ = _core.cf_betweenness(G.core, w, STRENGTH, int(bool(normalized)), threads(), 0)
    return _node_dict(G, doubles(nodes))


@implements(
    "edge_current_flow_betweenness_centrality", _CF,
    declined="directed graphs; multigraphs; weight functions; weights <= 0; disconnected graphs; graphs of fewer than three nodes; self-loops.",
)
def edge_current_flow_betweenness_centrality(G, normalized=True, weight=None, dtype=float,
                                             solver="full"):
    w = _weighted(G, weight)
    _current_flow_ready(G, w)
    _, edges = _core.cf_betweenness(G.core, w, STRENGTH, int(bool(normalized)), threads(), 1)
    return _rcm_edge_dict(G, doubles(edges))


@implements(
    "current_flow_betweenness_centrality_subset",
    chk.all_of(_CF, lambda a: _hashable_items(a, "sources", "targets")),
    declined="directed graphs; multigraphs; weight functions; weights <= 0; disconnected graphs; graphs of fewer than three nodes; self-loops.",
)
def current_flow_betweenness_centrality_subset(G, sources, targets, normalized=True,
                                               weight=None, dtype=float, solver="lu"):
    w = _weighted(G, weight)
    _current_flow_ready(G, w)
    s = _members(G, sources, "sources", False)
    t = _members(G, targets, "targets", False)
    nodes, _ = _core.cf_betweenness_subset(G.core, s, t, 1, w, STRENGTH, int(bool(normalized)),
                                           threads(), 0)
    return _node_dict(G, doubles(nodes))


@implements(
    "edge_current_flow_betweenness_centrality_subset",
    chk.all_of(_CF, lambda a: _hashable_items(a, "sources", "targets")),
    declined="directed graphs; multigraphs; weight functions; weights <= 0; disconnected graphs; graphs of fewer than three nodes; self-loops.",
)
def edge_current_flow_betweenness_centrality_subset(G, sources, targets, normalized=True,
                                                    weight=None, dtype=float, solver="lu"):
    w = _weighted(G, weight)
    _current_flow_ready(G, w)
    s = _members(G, sources, "sources", False)
    t = _members(G, targets, "targets", False)
    _, edges = _core.cf_betweenness_subset(G.core, s, t, 1, w, STRENGTH, int(bool(normalized)),
                                           threads(), 1)
    return _rcm_edge_dict(G, doubles(edges))


# ---- Walks and communicability ---------------------------------------------------

_BINARY = chk.all_of(SIMPLE, chk.undirected)


@implements(
    "communicability_betweenness_centrality", _BINARY,
    declined="directed graphs; multigraphs.",
)
def communicability_betweenness_centrality(G):
    if len(G) == 0 or not _connected(G):
        # NetworkX divides by zero between components and returns NaN.
        raise NotImplementedError("disconnected graphs are left to NetworkX")
    return _node_dict(G, doubles(_core.communicability_betweenness(G.core, threads())))


def _communicability_matrix(G):
    n = len(G)
    return _rows(doubles(_core.communicability(G.core)), n)


@implements("communicability", _BINARY, declined="directed graphs; multigraphs.")
def communicability(G):
    rows = _communicability_matrix(G)
    nodes = G.nodes
    return {u: dict(zip(nodes, row)) for u, row in zip(nodes, rows)}


@implements("communicability_exp", _BINARY, declined="directed graphs; multigraphs.")
def communicability_exp(G):
    return communicability(G)


@implements("subgraph_centrality", _BINARY, declined="directed graphs; multigraphs.")
def subgraph_centrality(G, *, normalized=False):
    if normalized:
        raise NotImplementedError("normalized scores need the largest eigenvalue")
    rows = _communicability_matrix(G)
    vals = [rows[i][i] for i in range(len(G))]
    if normalized:
        top = max(vals) if vals else 1.0
        vals = [v / top for v in vals]
    return _node_dict(G, vals)


@implements("subgraph_centrality_exp", _BINARY, declined="directed graphs; multigraphs.")
def subgraph_centrality_exp(G, *, normalized=False):
    return subgraph_centrality(G, normalized=normalized)


@implements("estrada_index", _BINARY, declined="directed graphs; multigraphs.")
def estrada_index(G):
    rows = _communicability_matrix(G)
    return sum(rows[i][i] for i in range(len(G)))


# ---- Resistance and the random walk -----------------------------------------------


def _check_resistance(a):
    G = chk.graph_of(a)
    for p in ("nodeA", "nodeB"):
        v = a.get(p)
        if v is not None and not (_hashable(v) and v in G):
            return f"{p} is not a node of the graph"
    return None


@implements(
    "resistance_distance",
    chk.all_of(SIMPLE, chk.undirected, chk.nonempty, _WEIGHT, _check_resistance),
    declined="directed graphs; multigraphs; weight functions; weights <= 0; disconnected graphs.",
)
def resistance_distance(G, nodeA=None, nodeB=None, weight=None, invert_weight=True):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    if not _connected(G):
        raise NotImplementedError("disconnected graphs are left to NetworkX")
    n = len(G)
    as_ = DISTANCE if invert_weight else STRENGTH
    rows = _rows(doubles(_core.resistance_distance(G.core, w, as_, threads())), n)
    nodes = G.nodes
    index = G.index
    if nodeA is not None and nodeB is not None:
        return rows[index[nodeA]][index[nodeB]]
    if nodeA is not None:
        return dict(zip(nodes, rows[index[nodeA]]))
    if nodeB is not None:
        j = index[nodeB]
        return {u: rows[i][j] for i, u in enumerate(nodes)}
    return {u: dict(zip(nodes, rows[i])) for i, u in enumerate(nodes)}


@implements(
    "effective_graph_resistance",
    chk.all_of(SIMPLE, chk.undirected, chk.nonempty, _WEIGHT),
    declined="directed graphs; multigraphs; weight functions; weights <= 0.",
)
def effective_graph_resistance(G, weight=None, invert_weight=True):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    if not _connected(G):
        return math.inf
    as_ = DISTANCE if invert_weight else STRENGTH
    return _core.effective_resistance(G.core, w, as_, threads())


@implements(
    "kemeny_constant",
    chk.all_of(SIMPLE, chk.undirected, chk.nonempty, chk.equals("weight", "weight")),
    declined="directed graphs; multigraphs; any weight other than the attribute 'weight' "
             "(NetworkX reads that attribute for the transition matrix whatever the argument); "
             "disconnected graphs; negative weights.",
)
def kemeny_constant(G, *, weight=None):
    w = _weighted(G, weight)
    _nonnegative_weights(G, w)
    if not _connected(G):
        raise NotImplementedError("disconnected graphs are left to NetworkX")
    return _core.kemeny(G.core, w)


@implements(
    "second_order_centrality",
    chk.all_of(SIMPLE, chk.undirected, chk.nonempty, _WEIGHT),
    declined="directed graphs; multigraphs; weight functions; disconnected graphs; negative weights.",
)
def second_order_centrality(G, weight="weight"):
    if len(G) < 2 or weight != "weight":
        raise NotImplementedError("left to NetworkX")
    w = _weighted(G, weight)
    _nonnegative_weights(G, w)
    if not _connected(G):
        raise NotImplementedError("disconnected graphs are left to NetworkX")
    return _node_dict(G, doubles(_core.second_order(G.core, w)))


# ---- Dispersion, VoteRank, Laplacian centrality ---------------------------------------


def _check_dispersion(a):
    G = chk.graph_of(a)
    for p in ("u", "v"):
        x = a.get(p)
        if x is not None and not (_hashable(x) and x in G):
            return f"{p} is not a node of the graph"
    return None


@implements(
    "dispersion",
    chk.all_of(SIMPLE, chk.undirected, _check_dispersion, chk.no_selfloops),
    declined="directed graphs; multigraphs; self-loops.",
)
def dispersion(G, u=None, v=None, normalized=True, alpha=1.0, b=0.0, c=0.0):
    index = G.index
    nodes = G.nodes
    adj = G.graph._succ if G.core.directed else G.graph._adj

    def ordered(x):
        return [y for y in adj[x]]

    if u is not None and v is not None:
        if v not in adj[u]:
            raise NotImplementedError("u and v are not adjacent")
        pairs = [(u, v)]
    elif u is None and v is None:
        pairs = [(x, y) for x in G.nodes for y in ordered(x)]
    elif u is not None:
        pairs = [(u, y) for y in ordered(u)]
    else:
        pairs = [(v, y) for y in ordered(v)]
    flat = []
    for x, y in pairs:
        flat += (index[x], index[y])
    vals = doubles(_core.dispersion(G.core, flat, int(bool(normalized)), float(alpha), float(b),
                                    float(c)))
    if u is not None and v is not None:
        return vals[0]
    if u is None and v is None:
        out = {x: {} for x in nodes}
        for (x, y), val in zip(pairs, vals):
            out[x][y] = val
        return out
    if u is not None:
        return {y: val for (_, y), val in zip(pairs, vals)}
    return {y: val for (_, y), val in zip(pairs, vals)}


@implements("voterank", ANY, declined="")
def voterank(G, number_of_nodes=None):
    n = len(G)
    if number_of_nodes is not None and not isinstance(number_of_nodes, int):
        raise NotImplementedError("number_of_nodes must be an integer")
    count = -1 if number_of_nodes is None or number_of_nodes > n else number_of_nodes
    if count == 0:
        return []
    return [G.nodes[i] for i in _core.voterank(G.core, count)]


@implements(
    "laplacian_centrality",
    chk.all_of(SIMPLE, chk.undirected, _WEIGHT, chk.is_none("nodelist")),
    declined="directed graphs; multigraphs; weight functions; nodelist; graphs without edges.",
)
def laplacian_centrality(G, normalized=True, nodelist=None, weight="weight", walk_type=None,
                         alpha=0.95):
    w = _weighted(G, weight)
    if len(G) == 0 or G.core.m == 0:
        raise NotImplementedError("left to NetworkX")
    vals = doubles(_core.laplacian_centrality(G.core, w, int(bool(normalized))))
    return _node_dict(G, vals)


# ---- Trophic structure and reaching centrality ------------------------------------------


def _trophic_ready(G, w):
    _positive_weights(G, w)
    # NetworkX needs a basal node and a path from one to every node.
    indeg = [d for _, d in G.graph.in_degree()]
    if 0 not in indeg:
        raise NotImplementedError("no basal node; NetworkX raises")


@implements(
    "trophic_levels",
    chk.all_of(SIMPLE, chk.directed, chk.nonempty, _WEIGHT),
    declined="undirected graphs; multigraphs; weight functions; graphs where some node cannot be reached from a basal node.",
)
def trophic_levels(G, weight="weight"):
    w = _weighted(G, weight)
    _trophic_ready(G, w)
    vals = doubles(_core.trophic(G.core, 0, w, 0))
    if not all(math.isfinite(x) for x in vals):
        raise NotImplementedError("trophic levels are undefined; NetworkX raises")
    out = _node_dict(G, vals)
    indeg = dict(G.graph.in_degree())
    for v in out:
        if indeg[v] == 0:
            out[v] = 1
    return out


@implements(
    "trophic_differences",
    chk.all_of(SIMPLE, chk.directed, chk.nonempty, _WEIGHT),
    declined="undirected graphs; multigraphs; weight functions; graphs where some node cannot be reached from a basal node.",
)
def trophic_differences(G, weight="weight"):
    w = _weighted(G, weight)
    _trophic_ready(G, w)
    vals = doubles(_core.trophic(G.core, 1, w, 0))
    if not all(math.isfinite(x) for x in vals):
        raise NotImplementedError("trophic levels are undefined; NetworkX raises")
    return _edge_dict(G, vals)


@implements(
    "trophic_incoherence_parameter",
    chk.all_of(SIMPLE, chk.directed, chk.nonempty, _WEIGHT),
    declined="undirected graphs; multigraphs; weight functions; graphs where some node cannot be reached from a basal node.",
)
def trophic_incoherence_parameter(G, weight="weight", cannibalism=False):
    w = _weighted(G, weight)
    _trophic_ready(G, w)
    val = _core.trophic(G.core, 2, w, int(bool(cannibalism)))
    x = doubles(val)[0]
    if not math.isfinite(x):
        raise NotImplementedError("undefined; NetworkX raises")
    return x


def _reaching_ready(G, w):
    if w:
        # Among tied shortest paths NetworkX averages the weights along the
        # path its search meets first; the engine's can differ.
        raise NotImplementedError("weighted reaching centrality depends on tied paths")
    _positive_weights(G, w)
    if G.core.m == 0:
        raise NotImplementedError("graphs without edges are left to NetworkX")


@implements(
    "global_reaching_centrality",
    chk.all_of(SIMPLE, chk.nonempty, _WEIGHT),
    declined="multigraphs; weights; graphs without edges.",
)
def global_reaching_centrality(G, weight=None, normalized=True):
    w = _weighted(G, weight)
    _reaching_ready(G, w)
    if len(G) < 2:
        raise NotImplementedError("left to NetworkX")
    _, g = _core.reaching(G.core, w, int(bool(normalized)), threads())
    if math.isnan(g):
        raise NotImplementedError("undefined")
    return g


@implements(
    "local_reaching_centrality",
    chk.all_of(SIMPLE, chk.nonempty, _WEIGHT, chk.node_in_graph("v"), chk.is_none("paths")),
    declined="multigraphs; weights; paths; graphs without edges.",
)
def local_reaching_centrality(G, v, paths=None, weight=None, normalized=True):
    w = _weighted(G, weight)
    _reaching_ready(G, w)
    if len(G) < 2:
        raise NotImplementedError("left to NetworkX")
    local, _ = _core.reaching(G.core, w, int(bool(normalized)), threads())
    x = doubles(local)[G.index[v]]
    if math.isnan(x):
        raise NotImplementedError("undefined")
    return x
