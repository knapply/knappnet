"""Implementations of NetworkX dispatchable functions on the engine.

Each function has NetworkX's signature and returns what NetworkX returns:
dicts keyed by the original node labels, generators where NetworkX yields,
and the same orientation and normalization. The engine runs with the
NetworkX convention wherever conventions differ between tools (the
convention tables in crates/engine/src/conventions.rs and the NetworkX
oracles in validation/oracles/*networkx*.py record which settings match
NetworkX). Calls that the engine cannot reproduce exactly are declined in
``can_run`` (see _checks.py) or raise NotImplementedError here, so that
NetworkX runs its own implementation.

Codes passed to the engine follow crates/capi/include/netcapi.h: mode
1 = out, 2 = in, 3 = all; loops 1 = an undirected self-loop counts once
(NetworkX's adjacency matrix); unreachable 0 = ignore, 1 = Wasserman-Faust.
"""

from __future__ import annotations

import math
import random

import networkx as nx

from . import _checks as chk
from . import _core
from ._graph import implements, threads

OUT, IN, ALL = 1, 2, 3
LOOPS_ONCE = 1
UNREACHABLE_IGNORE, UNREACHABLE_WF = 0, 1

# Dense distance blocks are kept below this many entries.
_BLOCK_ENTRIES = 1 << 22


# ---- Helpers -------------------------------------------------------------


def _weighted(G, weight):
    """1 if the call uses the converted edge attribute, else 0."""
    if weight is None:
        return 0
    if G.weight != weight:
        raise NotImplementedError("graph was converted without this edge attribute")
    if not G.core.wfinite or math.isnan(G.core.wmin):
        raise NotImplementedError("weights must be finite numbers")
    return 1


def _positive_weights(G, w):
    if w and not G.core.wmin > 0:
        raise NotImplementedError("non-positive weights are left to NetworkX")


def _nonnegative_weights(G, w):
    if w and not G.core.wmin >= 0:
        raise NotImplementedError("negative weights are left to NetworkX")


def _node_dict(G, values):
    return dict(zip(G.nodes, values))


def _subset(G, result, nodes):
    """NetworkX's `nodes` argument: None, one node (a scalar result), or an nbunch."""
    if nodes is None:
        return result
    if nodes in G.graph:
        return result[nodes]
    return {v: result[v] for v in G.graph.nbunch_iter(nodes)}


def _nan_to_zero(values):
    return [0.0 if math.isnan(x) else x for x in values]


def _ints(values):
    return [int(x) for x in values]


def _seed64(seed):
    """A 64-bit engine seed from NetworkX's seed argument."""
    if seed is None:
        return random.getrandbits(64)
    if isinstance(seed, int) and not isinstance(seed, bool):
        return seed % (1 << 64)
    return nx.utils.create_py_random_state(seed).getrandbits(64)


def _membership(G, communities):
    """Labels of a partition, or None if `communities` is not a partition of G."""
    index = G.index
    memb = [-1] * len(index)
    for label, comm in enumerate(communities):
        for v in comm:
            i = index.get(v) if _hashable(v) else None
            if i is None or memb[i] != -1:
                return None
            memb[i] = label
    if -1 in memb:
        return None
    return memb


def _nx_version():
    parts = []
    for x in nx.__version__.split(".")[:2]:
        digits = "".join(c for c in x if c.isdigit())
        parts.append(int(digits or 0))
    return tuple(parts)


def _hashable(v):
    try:
        hash(v)
    except TypeError:
        return False
    return True


def _groups(G, membership):
    groups = {}
    for v, label in zip(G.nodes, membership):
        groups.setdefault(label, set()).add(v)
    return list(groups.values())


def _connected(G):
    """True if G is connected (strongly, for directed graphs)."""
    kind = 1 if G.core.directed else 0
    return _core.component_count(G.core, kind) == 1


def _source_blocks(n):
    size = max(1, _BLOCK_ENTRIES // max(n, 1))
    for start in range(0, n, size):
        yield list(range(start, min(n, start + size)))


SIMPLE = chk.all_of(chk.nonempty, chk.no_multigraph)
ANY = chk.nonempty


# ---- Centrality: degree --------------------------------------------------


def _degree_centrality(G, mode):
    n = len(G)
    if n <= 1:
        return {v: 1 for v in G.nodes}
    s = 1.0 / (n - 1.0)
    return _node_dict(G, [d * s for d in _core.degree(G.core, mode, 1, 0)])


@implements("degree_centrality", ANY)
def degree_centrality(G):
    return _degree_centrality(G, ALL)


@implements("in_degree_centrality", chk.all_of(ANY, chk.directed))
def in_degree_centrality(G):
    return _degree_centrality(G, IN)


@implements("out_degree_centrality", chk.all_of(ANY, chk.directed))
def out_degree_centrality(G):
    return _degree_centrality(G, OUT)


# ---- Centrality: shortest paths ------------------------------------------


@implements(
    "betweenness_centrality",
    chk.all_of(SIMPLE, chk.is_none("k"), chk.equals("endpoints", False), chk.weight_attr()),
    declined="k (sampling); endpoints=True; multigraphs; weight functions; weights <= 0.",
)
def betweenness_centrality(G, k=None, normalized=True, weight=None, endpoints=False,
                           seed=None):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    vals = _core.betweenness(G.core, w, int(bool(normalized)), threads())
    return _node_dict(G, vals)


@implements(
    "edge_betweenness_centrality",
    chk.all_of(SIMPLE, chk.is_none("k"), chk.weight_attr()),
    declined="k (sampling); multigraphs; weight functions; weights <= 0.",
)
def edge_betweenness_centrality(G, k=None, normalized=True, weight=None, seed=None):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    vals = _core.edge_betweenness(G.core, w, int(bool(normalized)), threads())
    return dict(zip(_core.edge_keys(G.core), vals))


@implements(
    "closeness_centrality",
    chk.all_of(SIMPLE, chk.is_none("u", "sp"), chk.weight_attr("distance")),
    declined="u (one node); sp; multigraphs; distance functions; weights <= 0.",
)
def closeness_centrality(G, u=None, distance=None, wf_improved=True, *, sp=None):
    w = _weighted(G, distance)
    _positive_weights(G, w)
    if len(G) == 1:
        return {G.nodes[0]: 0.0}
    rule = UNREACHABLE_WF if wf_improved else UNREACHABLE_IGNORE
    vals = _core.closeness(G.core, 0, IN, w, 1, rule, threads())
    return _node_dict(G, _nan_to_zero(vals))


@implements(
    "harmonic_centrality",
    chk.all_of(SIMPLE, chk.is_none("sources"), chk.weight_attr("distance")),
    declined="sources; multigraphs; distance functions; weights <= 0.",
)
def harmonic_centrality(G, nbunch=None, distance=None, sources=None):
    w = _weighted(G, distance)
    _positive_weights(G, w)
    vals = _core.closeness(G.core, 1, IN, w, 0, 0, threads())
    result = _node_dict(G, vals)
    if nbunch is None:
        return result
    return {v: result[v] for v in set(G.graph.nbunch_iter(nbunch))}


# ---- Centrality: spectral and walk-based ----------------------------------


def _all_zero(values):
    return all(v == 0 for v in values)


@implements(
    "eigenvector_centrality",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.at_least("max_iter", 100)),
    declined="graphs that are not (strongly) connected; nstart with negative values; max_iter below 100; weights <= 0.",
)
def eigenvector_centrality(G, max_iter=100, tol=1e-06, nstart=None, weight=None):
    if nstart is not None:
        vals = list(nstart.values())
        if _all_zero(vals) or any(v < 0 for v in vals):
            raise NotImplementedError("nstart is left to NetworkX")
    return _eigenvector(G, weight)


@implements(
    "eigenvector_centrality_numpy",
    chk.all_of(SIMPLE, chk.weight_attr()),
    declined="graphs of at most two nodes; weights <= 0.",
)
def eigenvector_centrality_numpy(G, weight=None, max_iter=50, tol=0):
    if len(G) <= 2:
        # NetworkX's sparse eigensolver refuses graphs this small.
        raise NotImplementedError("left to NetworkX")
    return _eigenvector(G, weight)


def _eigenvector(G, weight):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    vals, _, connected = _core.eigenvector(G.core, IN, w, 1, LOOPS_ONCE, 0.0, 0, threads())
    if not connected:
        raise NotImplementedError("disconnected graphs are left to NetworkX")
    return _node_dict(G, vals)


def _beta_vector(G, beta):
    if isinstance(beta, dict):
        if set(beta) != set(G.nodes):
            raise NotImplementedError("beta must have a value for every node")
        return [float(beta[v]) for v in G.nodes]
    try:
        return [float(beta)]
    except (TypeError, ValueError) as err:
        raise NotImplementedError("beta must be a number") from err


@implements(
    "katz_centrality",
    chk.all_of(SIMPLE, chk.weight_attr(), chk.at_least("max_iter", 1000)),
    declined="graphs of one node; alpha >= 1 / (largest eigenvalue), where NetworkX does not converge; max_iter below 1000; negative weights.",
)
def katz_centrality(G, alpha=0.1, beta=1.0, max_iter=1000, tol=1e-06, nstart=None,
                    normalized=True, weight=None):
    return _katz(G, alpha, beta, normalized, weight)


@implements(
    "katz_centrality_numpy",
    chk.all_of(SIMPLE, chk.weight_attr()),
    declined="graphs of one node; alpha >= 1 / (largest eigenvalue); negative weights.",
)
def katz_centrality_numpy(G, alpha=0.1, beta=1.0, normalized=True, weight=None):
    return _katz(G, alpha, beta, normalized, weight)


def _katz(G, alpha, beta, normalized, weight):
    if len(G) < 2:
        raise NotImplementedError("left to NetworkX")
    w = _weighted(G, weight)
    _nonnegative_weights(G, w)
    if not alpha > 0:
        raise NotImplementedError("alpha must be positive")
    b = _beta_vector(G, beta)
    vals = _core.katz(G.core, float(alpha), b, IN, w, LOOPS_ONCE, int(bool(normalized)),
                      0.0, 0, threads())
    return _node_dict(G, vals)


def _dangling_code(G, dangling, personalization):
    """0 = jump by the personalization vector, 1 = uniformly; None otherwise."""
    if dangling is None:
        return 0
    d = [float(dangling.get(v, 0)) for v in G.nodes]
    total = sum(d)
    if not total > 0:
        return None
    d = [x / total for x in d]
    n = len(G)
    if all(math.isclose(x, 1.0 / n, rel_tol=1e-15, abs_tol=0) for x in d):
        return 1
    if personalization is not None:
        p = [float(personalization.get(v, 0)) for v in G.nodes]
        ps = sum(p)
        if ps > 0 and all(math.isclose(x, y / ps, rel_tol=1e-15) for x, y in zip(d, p)):
            return 0
    return None


@implements(
    "pagerank",
    chk.all_of(ANY, chk.weight_attr(), chk.at_least("max_iter", 100)),
    declined="dangling vectors other than uniform or the personalization; alpha outside (0, 1); max_iter below 100; negative weights.",
)
def pagerank(G, alpha=0.85, personalization=None, max_iter=100, tol=1e-06, nstart=None,
             weight="weight", dangling=None):
    w = _weighted(G, weight)
    _nonnegative_weights(G, w)
    if not 0 < alpha < 1:
        raise NotImplementedError("alpha outside (0, 1) is left to NetworkX")
    pers = None
    if personalization is not None:
        pers = [float(personalization.get(v, 0)) for v in G.nodes]
        if not sum(pers) > 0 or min(pers) < 0:
            raise NotImplementedError("personalization is left to NetworkX")
    code = _dangling_code(G, dangling, personalization)
    if code is None:
        raise NotImplementedError("this dangling vector is left to NetworkX")
    vals = _core.pagerank(G.core, float(alpha), int(G.core.directed), w, pers, code,
                          LOOPS_ONCE, 0.0, 0, threads())
    return _node_dict(G, vals)


def _check_hits_version(a):
    # Before NetworkX 3.6 (no `method` argument), hits used SciPy's sparse
    # SVD and returned unscaled singular vectors with arbitrary signs when
    # normalized=False.
    if "method" not in a and not a.get("normalized", True):
        return "normalized=False is left to NetworkX before 3.6"
    return None


@implements(
    "hits",
    chk.all_of(SIMPLE, chk.is_none("nstart"), chk.at_least("max_iter", 100),
             chk.equals("method", None, "power_iteration"), _check_hits_version),
    declined="nstart; max_iter below 100; graphs whose leading singular value is repeated (scores not unique); multigraphs; negative weights.",
)
def hits(G, max_iter=100, tol=1e-08, nstart=None, normalized=True, *,
         method="power_iteration"):
    w = 1 if G.weight is not None else 0
    if w:
        _nonnegative_weights(G, w)
        if not G.core.wfinite:
            raise NotImplementedError("weights must be finite")
    scale = 1 if normalized else 0
    h, a, _, unique = _core.hits(G.core, w, LOOPS_ONCE, scale, 0.0, 0, threads())
    if not unique:
        raise NotImplementedError("the scores are not unique; left to NetworkX")
    return _node_dict(G, h), _node_dict(G, a)


@implements(
    "current_flow_closeness_centrality",
    chk.all_of(SIMPLE, chk.undirected, chk.weight_attr(), chk.equals("solver", "lu", "full", "cg")),
    declined="graphs of one node; multigraphs; weights <= 0.",
)
def current_flow_closeness_centrality(G, weight=None, dtype=float, solver="lu"):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    n = len(G)
    if n < 2 or not _connected(G):
        raise NotImplementedError("disconnected graphs are left to NetworkX")
    vals = _core.information_centrality(G.core, w)
    return _node_dict(G, [x / n for x in vals])


# ---- Centrality: group ---------------------------------------------------


def _group_indices(G, S):
    index = G.index
    try:
        idx = [index[v] for v in S]
    except (KeyError, TypeError) as err:
        raise NotImplementedError("group members must be nodes") from err
    if len(set(idx)) != len(idx) or not idx or len(idx) >= len(G):
        raise NotImplementedError("this group is left to NetworkX")
    return idx


@implements(
    "group_degree_centrality", SIMPLE,
    declined="multigraphs; groups with repeated members.",
)
def group_degree_centrality(G, S):
    idx = _group_indices(G, S)
    return _core.group_centrality(G.core, idx, 0, OUT, 0, 1, threads())


@implements(
    "group_in_degree_centrality", chk.all_of(SIMPLE, chk.directed),
    declined="multigraphs; groups with repeated members.",
)
def group_in_degree_centrality(G, S):
    idx = _group_indices(G, S)
    return _core.group_centrality(G.core, idx, 0, IN, 0, 1, threads())


@implements(
    "group_out_degree_centrality", chk.all_of(SIMPLE, chk.directed),
    declined="multigraphs; groups with repeated members.",
)
def group_out_degree_centrality(G, S):
    idx = _group_indices(G, S)
    return _core.group_centrality(G.core, idx, 0, OUT, 0, 1, threads())


@implements(
    "group_closeness_centrality", chk.all_of(SIMPLE, chk.weight_attr()),
    declined="multigraphs; groups with repeated members; weights <= 0.",
)
def group_closeness_centrality(G, S, weight=None):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    idx = _group_indices(G, S)
    return _core.group_centrality(G.core, idx, 1, IN, w, 1, threads())


@implements(
    "group_betweenness_centrality",
    chk.all_of(SIMPLE, chk.equals("endpoints", False), chk.weight_attr()),
    declined="endpoints=True; multigraphs; groups with repeated members; weights <= 0.",
)
def group_betweenness_centrality(G, C, normalized=True, weight=None, endpoints=False):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    try:
        first = next(iter(C))
    except (StopIteration, TypeError) as err:
        raise NotImplementedError("empty group") from err
    single = _hashable(first) and first in G.index
    groups = [C] if single else C
    out = []
    for grp in groups:
        idx = _group_indices(G, grp)
        out.append(_core.group_centrality(G.core, idx, 2, ALL, w, int(bool(normalized)),
                                          threads()))
    return out[0] if single else out


# ---- Structural holes ----------------------------------------------------


def _structural_holes(G, measure, nodes, weight):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    vals = _core.structural_holes(G.core, measure, w, threads())
    result = _node_dict(G, vals)
    if nodes is None:
        return result
    return {v: result[v] for v in nodes}


@implements(
    "constraint",
    chk.all_of(SIMPLE, chk.no_selfloops, chk.weight_attr()),
    declined="self-loops; multigraphs; weights <= 0.",
)
def constraint(G, nodes=None, weight=None):
    result = _structural_holes(G, 0, nodes, weight)
    if G.core.directed and _nx_version() < (3, 5):
        # NetworkX 3.4 gives NaN to nodes without successors.
        succ = G.graph._succ
        for v in result:
            if len(succ[v]) == 0:
                result[v] = math.nan
    return result


def _check_effective_size(a):
    # With nodes=None (and SciPy installed) NetworkX 3.7 takes a vectorized
    # path that scales each mutual tie z_jq by q's strongest tie instead of
    # j's, so for weighted or directed graphs it disagrees with its own
    # per-node path and with Burt's definition. Those calls are left to it.
    G = chk.graph_of(a)
    if a.get("nodes") is None and (G.is_directed() or a.get("weight") is not None):
        return "NetworkX's vectorized path is left to NetworkX"
    return None


@implements(
    "effective_size",
    chk.all_of(SIMPLE, chk.no_selfloops, chk.weight_attr(), _check_effective_size),
    declined="self-loops; multigraphs; weighted or directed graphs with nodes=None; weights <= 0.",
)
def effective_size(G, nodes=None, weight=None):
    result = _structural_holes(G, 1, nodes, weight)
    if G.core.directed:
        # NetworkX's per-node path gives NaN to nodes without successors.
        succ = G.graph._succ
        for v in result:
            if all(u == v for u in succ[v]):
                result[v] = math.nan
    return result


# ---- Clustering, reciprocity, triads ---------------------------------------


@implements(
    "transitivity", chk.all_of(SIMPLE, chk.undirected),
    declined="multigraphs; directed graphs.",
)
def transitivity(G):
    t = _core.transitivity(G.core, 1, threads())
    return 0 if t == 0 else t


@implements(
    "clustering",
    chk.all_of(SIMPLE, chk.undirected, chk.weight_attr()),
    declined="multigraphs; directed graphs; weights <= 0.",
)
def clustering(G, nodes=None, weight=None):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    vals = _core.clustering(G.core, 0, 2 if w else 0, 1, threads())
    result = _node_dict(G, [0 if x == 0 else x for x in vals])
    return _subset(G, result, nodes)


@implements(
    "average_clustering",
    chk.all_of(SIMPLE, chk.undirected, chk.weight_attr()),
    declined="multigraphs; directed graphs; weights <= 0.",
)
def average_clustering(G, nodes=None, weight=None, count_zeros=True):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    if nodes is None and count_zeros:
        return _core.clustering(G.core, 1, 2 if w else 0, 1, threads())
    vals = _node_dict(G, _core.clustering(G.core, 0, 2 if w else 0, 1, threads()))
    c = [vals[v] for v in G.graph.nbunch_iter(nodes)] if nodes is not None else list(
        vals.values())
    if not count_zeros:
        c = [v for v in c if abs(v) > 0]
    if not c:
        raise NotImplementedError("no nodes to average; left to NetworkX")
    return sum(c) / len(c)


@implements("overall_reciprocity",
    chk.all_of(SIMPLE, chk.directed))
def overall_reciprocity(G):
    m = G.core.m
    if m == 0:
        raise NotImplementedError("empty graphs are left to NetworkX")
    loops = G.core.loops
    if m == loops:
        return 0.0
    r = _core.reciprocity(G.core, 0)
    # The engine ignores self-loops; NetworkX counts them as unreciprocated arcs.
    return r * (m - loops) / m if loops else r


@implements(
    "reciprocity",
    chk.all_of(SIMPLE, chk.directed, chk.is_none("nodes")),
    declined="nodes (per-node reciprocity).",
)
def reciprocity(G, nodes=None):
    return overall_reciprocity(G)


_TRIADS = ("003", "012", "102", "021D", "021U", "021C", "111D", "111U", "030T",
           "030C", "201", "120D", "120U", "120C", "210", "300")


@implements(
    "triadic_census",
    chk.all_of(SIMPLE, chk.directed, chk.is_none("nodelist"), chk.no_selfloops),
    declined="nodelist; self-loops.",
)
def triadic_census(G, nodelist=None):
    counts = _core.triad_census(G.core, threads())
    return dict(zip(_TRIADS, _ints(counts)))


# ---- Cores and cliques ----------------------------------------------------


@implements("core_number", chk.all_of(SIMPLE, chk.no_selfloops))
def core_number(G):
    return _node_dict(G, _ints(_core.coreness(G.core, ALL)))


@implements(
    "find_cliques",
    chk.all_of(SIMPLE, chk.undirected, chk.is_none("nodes")),
    declined="nodes; multigraphs.",
)
def find_cliques(G, nodes=None):
    return iter(_core.max_cliques(G.core, 0, threads()))


# ---- Components -------------------------------------------------------------


def _components(G, kind):
    return (set(c) for c in _core.components(G.core, kind))


@implements("connected_components", chk.all_of(ANY, chk.undirected))
def connected_components(G):
    return _components(G, 0)


@implements("number_connected_components", chk.all_of(ANY, chk.undirected))
def number_connected_components(G):
    return _core.component_count(G.core, 0)


@implements("is_connected", chk.all_of(ANY, chk.undirected))
def is_connected(G):
    return _core.component_count(G.core, 0) == 1


@implements("node_connected_component",
            chk.all_of(ANY, chk.undirected, chk.node_in_graph("n")))
def node_connected_component(G, n):
    return set(_one_source(G, n, OUT, 0))


@implements("strongly_connected_components", chk.all_of(ANY, chk.directed))
def strongly_connected_components(G):
    return (set(c) for c in _core.strong_components(G.core))


@implements("number_strongly_connected_components", chk.all_of(ANY, chk.directed))
def number_strongly_connected_components(G):
    return _core.component_count(G.core, 1)


@implements("is_strongly_connected", chk.all_of(ANY, chk.directed))
def is_strongly_connected(G):
    return _core.component_count(G.core, 1) == 1


@implements("weakly_connected_components", chk.all_of(ANY, chk.directed))
def weakly_connected_components(G):
    return _components(G, 0)


@implements("number_weakly_connected_components", chk.all_of(ANY, chk.directed))
def number_weakly_connected_components(G):
    return _core.component_count(G.core, 0)


@implements("is_weakly_connected", chk.all_of(ANY, chk.directed))
def is_weakly_connected(G):
    return _core.component_count(G.core, 0) == 1


@implements(
    "articulation_points", chk.all_of(SIMPLE, chk.undirected),
    declined="multigraphs.",
)
def articulation_points(G):
    cut, _ = _core.cut_structure(G.core)
    return (v for v, flag in zip(G.nodes, cut) if flag)


@implements(
    "bridges", chk.all_of(SIMPLE, chk.undirected, chk.is_none("root")),
    declined="root; multigraphs.",
)
def bridges(G, root=None):
    _, flags = _core.cut_structure(G.core)
    return (e for e, flag in zip(_core.edge_keys(G.core), flags) if flag)


@implements(
    "has_bridges", chk.all_of(SIMPLE, chk.undirected, chk.is_none("root")),
    declined="root; multigraphs.",
)
def has_bridges(G, root=None):
    _, flags = _core.cut_structure(G.core)
    return any(flags)


@implements(
    "biconnected_components", chk.all_of(SIMPLE, chk.undirected, chk.no_selfloops),
    declined="self-loops; multigraphs.",
)
def biconnected_components(G):
    return (set(c) for c in _core.bicomponents(G.core))


@implements(
    "is_biconnected", chk.all_of(SIMPLE, chk.undirected, chk.no_selfloops),
    declined="self-loops; multigraphs.",
)
def is_biconnected(G):
    comps = _core.bicomponents(G.core)
    return len(comps) == 1 and len(comps[0]) == len(G)


# ---- Shortest path lengths ------------------------------------------------


def _distance_dicts(G, sources, mode, w, cutoff=None):
    c = -1.0 if cutoff is None else float(cutoff)
    return _core.distances(G.core, sources, mode, w, threads(), 0 if w else 1, c)


def _one_source(G, source, mode, w, cutoff=None):
    return _distance_dicts(G, [G.index[source]], mode, w, cutoff)[0]


def _all_pairs(G, w, cutoff=None):
    nodes = G.nodes
    for block in _source_blocks(len(G)):
        for i, d in zip(block, _distance_dicts(G, block, OUT, w, cutoff)):
            yield nodes[i], d


def _check_cutoff(cutoff):
    if cutoff is not None and not isinstance(cutoff, (int, float)):
        raise NotImplementedError("cutoff must be a number")


@implements("single_source_shortest_path_length",
            chk.all_of(ANY, chk.node_in_graph("source")))
def single_source_shortest_path_length(G, source, cutoff=None):
    _check_cutoff(cutoff)
    return _one_source(G, source, OUT, 0, cutoff)


@implements("single_target_shortest_path_length",
            chk.all_of(ANY, chk.node_in_graph("target")))
def single_target_shortest_path_length(G, target, cutoff=None):
    _check_cutoff(cutoff)
    d = _one_source(G, target, IN, 0, cutoff)
    if _nx_version() < (3, 5):
        # NetworkX 3.4 returns an iterator of (node, length) pairs.
        return iter(sorted(d.items(), key=lambda kv: kv[1]))
    return d


@implements("all_pairs_shortest_path_length", ANY)
def all_pairs_shortest_path_length(G, cutoff=None):
    _check_cutoff(cutoff)
    return _all_pairs(G, 0, cutoff)


@implements(
    "single_source_dijkstra_path_length",
            chk.all_of(ANY, chk.node_in_graph("source"), chk.weight_attr()),
    declined="weight functions; weights <= 0.",
)
def single_source_dijkstra_path_length(G, source, cutoff=None, weight="weight"):
    _check_cutoff(cutoff)
    w = _weighted(G, weight)
    _positive_weights(G, w)
    return _one_source(G, source, OUT, w, cutoff)


@implements(
    "all_pairs_dijkstra_path_length", chk.all_of(ANY, chk.weight_attr()),
    declined="weight functions; weights <= 0.",
)
def all_pairs_dijkstra_path_length(G, cutoff=None, weight="weight"):
    _check_cutoff(cutoff)
    w = _weighted(G, weight)
    _positive_weights(G, w)
    return _all_pairs(G, w, cutoff)


def _check_shortest_path_length(a):
    G = chk.graph_of(a)
    for p in ("source", "target"):
        v = a.get(p)
        if v is not None and not (_hashable(v) and v in G):
            return f"{p} is not a node of the graph"
    if a.get("method") not in ("dijkstra", "bellman-ford"):
        return "method is not supported"
    return None


@implements(
    "shortest_path_length",
            chk.all_of(ANY, chk.weight_attr(), _check_shortest_path_length),
    declined="weight functions; weights <= 0.",
)
def shortest_path_length(G, source=None, target=None, weight=None, method="dijkstra"):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    if source is None and target is None:
        return _all_pairs(G, w)
    if target is None:
        return _one_source(G, source, OUT, w)
    if source is None:
        return _one_source(G, target, IN, w)
    d = _one_source(G, source, OUT, w)
    if target not in d:
        raise NotImplementedError("no path; NetworkX raises")
    return d[target]


@implements("has_path", chk.all_of(ANY, chk.node_in_graph("source"), chk.node_in_graph("target")))
def has_path(G, source, target):
    return target in _one_source(G, source, OUT, 0)


@implements("descendants", chk.all_of(ANY, chk.node_in_graph("source")))
def descendants(G, source):
    return set(_one_source(G, source, OUT, 0)) - {source}


@implements("ancestors", chk.all_of(ANY, chk.node_in_graph("source")))
def ancestors(G, source):
    return set(_one_source(G, source, IN, 0)) - {source}


def _check_aspl_method(a):
    if a.get("method") not in (None, "unweighted", "dijkstra", "bellman-ford"):
        return "method is left to NetworkX"
    return None


@implements(
    "average_shortest_path_length",
            chk.all_of(ANY, chk.weight_attr(), _check_aspl_method),
    declined="Floyd-Warshall methods; weight functions; weights <= 0.",
)
def average_shortest_path_length(G, weight=None, method=None):
    if method is None:
        method = "unweighted" if weight is None else "dijkstra"
    w = 0 if method == "unweighted" else _weighted(G, weight)
    _positive_weights(G, w)
    n = len(G)
    if n == 1:
        return 0
    if not _connected(G):
        raise NotImplementedError("disconnected graphs are left to NetworkX")
    summary = _core.path_summary(G.core, OUT, w, 1, threads())
    return summary[2] / (n * (n - 1))


@implements("global_efficiency", chk.all_of(ANY, chk.undirected))
def global_efficiency(G):
    if len(G) < 2:
        return 0
    return _core.compactness(G.core, OUT, 0, threads())


@implements(
    "wiener_index",
    chk.all_of(ANY, chk.weight_attr()),
    declined="weight functions; weights <= 0.",
)
def wiener_index(G, weight=None):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    if not _connected(G):
        return float("inf")
    total = _core.path_summary(G.core, OUT, w, 1, threads())[2]
    if not w:
        total = int(total)
    return total if G.core.directed else total / 2


# ---- Distance measures ----------------------------------------------------


def _eccentricities(G, weight):
    w = _weighted(G, weight)
    _positive_weights(G, w)
    if not _connected(G):
        raise NotImplementedError("disconnected graphs are left to NetworkX")
    vals = _core.eccentricity(G.core, OUT, w, threads())
    return _node_dict(G, vals if w else _ints(vals))


_DM = chk.all_of(ANY, chk.is_none("e"), chk.weight_attr())


@implements(
    "eccentricity", chk.all_of(ANY, chk.is_none("v", "sp"), chk.weight_attr()),
    declined="v; sp; weight functions; weights <= 0.",
)
def eccentricity(G, v=None, sp=None, weight=None):
    return _eccentricities(G, weight)


@implements(
    "diameter", _DM,
    declined="e; weight functions; weights <= 0.",
)
def diameter(G, e=None, usebounds=False, weight=None):
    return max(_eccentricities(G, weight).values())


@implements(
    "radius", _DM,
    declined="e; weight functions; weights <= 0.",
)
def radius(G, e=None, usebounds=False, weight=None):
    return min(_eccentricities(G, weight).values())


@implements(
    "center", _DM,
    declined="e; weight functions; weights <= 0.",
)
def center(G, e=None, usebounds=False, weight=None):
    ecc = _eccentricities(G, weight)
    r = min(ecc.values())
    return [v for v in ecc if ecc[v] == r]


@implements(
    "periphery", _DM,
    declined="e; weight functions; weights <= 0.",
)
def periphery(G, e=None, usebounds=False, weight=None):
    ecc = _eccentricities(G, weight)
    d = max(ecc.values())
    return [v for v in ecc if ecc[v] == d]


# ---- Assortativity ----------------------------------------------------------


def _check_degree_xy(a):
    G = chk.graph_of(a)
    if G.is_directed() and (a.get("x"), a.get("y")) != ("out", "in"):
        return "only x='out', y='in' is supported for directed graphs"
    if a.get("nodes") is not None:
        return "nodes is left to NetworkX"
    return None


_DEG_ASSORT = chk.all_of(SIMPLE, chk.weight_attr(), _check_degree_xy)


@implements(
    "degree_assortativity_coefficient", _DEG_ASSORT,
    declined="x, y other than 'out', 'in' on directed graphs; nodes; multigraphs; negative weights.",
)
def degree_assortativity_coefficient(G, x="out", y="in", weight=None, nodes=None):
    w = _weighted(G, weight)
    _nonnegative_weights(G, w)
    r = _core.assortativity_degree(G.core, w, int(G.core.directed), LOOPS_ONCE, 1)
    if math.isnan(r):
        raise NotImplementedError("undefined; left to NetworkX")
    return r


@implements(
    "degree_pearson_correlation_coefficient", _DEG_ASSORT,
    declined="x, y other than 'out', 'in' on directed graphs; nodes; multigraphs; negative weights.",
)
def degree_pearson_correlation_coefficient(G, x="out", y="in", weight=None, nodes=None):
    return degree_assortativity_coefficient(G, x, y, weight, nodes)


def _node_values(G, attribute):
    data = G.graph._node
    try:
        return [data[v][attribute] for v in G.nodes]
    except (KeyError, TypeError) as err:
        raise NotImplementedError("every node needs the attribute") from err


@implements(
    "attribute_assortativity_coefficient",
            chk.all_of(SIMPLE, chk.is_none("nodes")),
    declined="nodes; multigraphs.",
)
def attribute_assortativity_coefficient(G, attribute, nodes=None):
    values = _node_values(G, attribute)
    codes = {}
    try:
        types = [codes.setdefault(v, len(codes)) for v in values]
    except TypeError as err:
        raise NotImplementedError("attribute values must be hashable") from err
    r = _core.assortativity_nominal(G.core, types, len(codes), int(G.core.directed),
                                    LOOPS_ONCE, 1)
    if math.isnan(r):
        raise NotImplementedError("undefined; left to NetworkX")
    return r


@implements(
    "numeric_assortativity_coefficient",
            chk.all_of(SIMPLE, chk.is_none("nodes")),
    declined="nodes; multigraphs.",
)
def numeric_assortativity_coefficient(G, attribute, nodes=None):
    values = _node_values(G, attribute)
    try:
        x = [float(v) for v in values]
    except (TypeError, ValueError) as err:
        raise NotImplementedError("attribute values must be numbers") from err
    r = _core.assortativity_values(G.core, x, None, int(G.core.directed), LOOPS_ONCE, 1)
    if math.isnan(r):
        raise NotImplementedError("undefined; left to NetworkX")
    return r


# ---- Communities -------------------------------------------------------------


@implements(
    "modularity", chk.all_of(SIMPLE, chk.weight_attr()),
    declined="multigraphs; negative weights.",
)
def modularity(G, communities, weight="weight", resolution=1):
    if not isinstance(communities, list):
        communities = list(communities)
    memb = _membership(G, communities)
    if memb is None:
        raise NotImplementedError("not a partition; NetworkX raises")
    if G.core.m == 0:
        # NetworkX 3.7 returns 0 here and earlier versions raise.
        raise NotImplementedError("graphs without edges are left to NetworkX")
    w = _weighted(G, weight)
    _nonnegative_weights(G, w)
    q = _core.modularity(G.core, memb, float(resolution), w)
    return 0 if math.isnan(q) else q


@implements(
    "louvain_communities",
    chk.all_of(SIMPLE, chk.undirected, chk.weight_attr(), chk.equals("threshold", 1e-07),
             chk.is_none("max_level")),
    declined="directed graphs; multigraphs; threshold; max_level; negative weights.",
)
def louvain_communities(G, weight="weight", resolution=1, threshold=1e-07, max_level=None,
                        seed=None):
    w = _weighted(G, weight)
    _nonnegative_weights(G, w)
    memb = _core.louvain(G.core, float(resolution), w, _seed64(seed))
    return _groups(G, memb)



# ---- Bipartite --------------------------------------------------------------


def _modes(G, nodes):
    top = set(nodes)
    if not top <= set(G.index):
        raise NotImplementedError("nodes must be nodes of the graph")
    return [0 if v in top else 1 for v in G.nodes], top


_BIP = chk.all_of(SIMPLE, chk.undirected)


@implements(
    "bipartite_degree_centrality", _BIP,
    declined="directed graphs; multigraphs.",
)
def bipartite_degree_centrality(G, nodes):
    types, top = _modes(G, nodes)
    if not top or len(top) == len(G):
        raise NotImplementedError("both modes must be non-empty")
    vals = _core.twomode_centrality(G.core, types, 0, 1, 0, threads())
    return _node_dict(G, vals)


@implements(
    "bipartite_betweenness_centrality", _BIP,
    declined="directed graphs; multigraphs.",
)
def bipartite_betweenness_centrality(G, nodes):
    types, top = _modes(G, nodes)
    if not top or len(top) == len(G):
        raise NotImplementedError("both modes must be non-empty")
    vals = _core.twomode_centrality(G.core, types, 2, 1, 0, threads())
    return _node_dict(G, vals)


@implements(
    "bipartite_closeness_centrality", _BIP,
    declined="directed graphs; multigraphs.",
)
def bipartite_closeness_centrality(G, nodes, normalized=True):
    types, top = _modes(G, nodes)
    if not top or len(top) == len(G):
        raise NotImplementedError("both modes must be non-empty")
    rule = UNREACHABLE_WF if normalized else UNREACHABLE_IGNORE
    vals = _core.twomode_centrality(G.core, types, 1, 1, rule, threads())
    return _node_dict(G, _nan_to_zero(vals))
