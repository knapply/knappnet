"""Link prediction, Weisfeiler-Lehman hashes, SimRank, graph edit distance,
and isomorphism tests.

Link prediction scores come from the engine in batches while the candidate
pairs are enumerated by NetworkX's own ``non_edges`` (whose order depends on
set iteration), so the generators yield the same triples in the same order.
Functions that return a mapping between isomorphic graphs, and the ones that
draw random paths (Panther), are left to NetworkX: which mapping, and which
random stream, is not determined by the graphs.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Collection

import networkx as nx

from . import _checks as chk
from . import _core
from ._algorithms import ALL, SIMPLE, _hashable, _node_dict, _weighted
from ._graph import implements, threads
from ._paths import doubles

_BATCH = 1 << 14

COMMON_NEIGHBORS, JACCARD, ADAMIC_ADAR, RESOURCE_ALLOCATION = 0, 1, 2, 3
PREFERENTIAL, CN_CENTRALITY, CN_SH, RA_SH, WITHIN_INTER = 4, 5, 6, 7, 8

_UNDIRECTED_SIMPLE = chk.all_of(SIMPLE, chk.undirected)


# ---- Link prediction --------------------------------------------------------


def _check_ebunch(a):
    e = a.get("ebunch")
    if e is None:
        return None
    G = chk.graph_of(a)
    if not isinstance(e, Collection) or isinstance(e, (str, bytes)):
        return "ebunch must be a collection of pairs"
    try:
        for p in e:
            if not isinstance(p, (tuple, list)) or len(p) != 2:
                return "ebunch has an item that is not a pair"
            if p[0] not in G or p[1] not in G:
                return "ebunch has a node outside the graph"
            if p[0] == p[1]:
                return "ebunch has a pair of one node"
    except TypeError:
        return "ebunch has unhashable nodes"
    return None


def _check_community(a):
    G = chk.graph_of(a)
    attr = a.get("community")
    try:
        if any(attr not in d for d in G._node.values()):
            return "a node has no community"
    except TypeError:
        return "community is not an attribute name"
    return None


def _community_codes(G, attr):
    codes, out = {}, []
    for v in G.graph._node.values():
        x = v[attr]
        try:
            if x != x:
                raise NotImplementedError("NaN community")
            out.append(codes.setdefault(x, len(codes)))
        except TypeError as err:
            raise NotImplementedError("unhashable community") from err
    return out


def _scores(G, method, ebunch, alpha=0.0, delta=0.0, communities=None, integral=False):
    index = G.index
    if ebunch is not None:
        pairs = [(u, v) for u, v in ebunch]
    else:
        pairs = None

    def convert(chunk):
        us = [index[u] for u, _ in chunk]
        vs = [index[v] for _, v in chunk]
        vals = doubles(_core.link_prediction(G.core, method, alpha, delta, communities, us, vs, threads()))
        return [int(x) for x in vals] if integral else vals

    if pairs is not None:
        scores = convert(pairs) if pairs else []
        return ((u, v, s) for (u, v), s in zip(pairs, scores))

    def run():
        it = nx.non_edges(G.graph)
        while True:
            chunk = list(itertools.islice(it, _BATCH))
            if not chunk:
                return
            for (u, v), s in zip(chunk, convert(chunk)):
                yield u, v, s

    return run()


def _link(name, method, integral=False):
    def make(fn):
        return implements(
            name,
            chk.all_of(_UNDIRECTED_SIMPLE, chk.no_selfloops, _check_ebunch),
            declined="directed graphs; multigraphs; self-loops; an ebunch that is not a collection of pairs.",
        )(fn)

    return make


@implements(
    "common_neighbor_centrality",
    chk.all_of(_UNDIRECTED_SIMPLE, chk.no_selfloops, _check_ebunch, chk.numeric_or_none("alpha")),
    declined="directed graphs; multigraphs; self-loops; an ebunch that is not a collection of pairs.",
)
def common_neighbor_centrality(G, ebunch=None, alpha=0.8):
    if alpha is None:
        raise NotImplementedError("alpha must be a number")
    return _scores(G, CN_CENTRALITY, ebunch, alpha=float(alpha), integral=alpha == 1)


@_link("resource_allocation_index", RESOURCE_ALLOCATION)
def resource_allocation_index(G, ebunch=None):
    return _scores(G, RESOURCE_ALLOCATION, ebunch)


@_link("jaccard_coefficient", JACCARD)
def jaccard_coefficient(G, ebunch=None):
    return _scores(G, JACCARD, ebunch)


@_link("adamic_adar_index", ADAMIC_ADAR)
def adamic_adar_index(G, ebunch=None):
    return _scores(G, ADAMIC_ADAR, ebunch)


@_link("preferential_attachment", PREFERENTIAL)
def preferential_attachment(G, ebunch=None):
    return _scores(G, PREFERENTIAL, ebunch, integral=True)


_COMMUNITY = chk.all_of(
    _UNDIRECTED_SIMPLE, chk.no_selfloops, _check_ebunch, _check_community
)


@implements(
    "cn_soundarajan_hopcroft",
    _COMMUNITY,
    declined="directed graphs; multigraphs; self-loops; nodes without a community; unhashable communities.",
)
def cn_soundarajan_hopcroft(G, ebunch=None, community="community"):
    return _scores(G, CN_SH, ebunch, communities=_community_codes(G, community), integral=True)


@implements(
    "ra_index_soundarajan_hopcroft",
    _COMMUNITY,
    declined="directed graphs; multigraphs; self-loops; nodes without a community; unhashable communities.",
)
def ra_index_soundarajan_hopcroft(G, ebunch=None, community="community"):
    return _scores(G, RA_SH, ebunch, communities=_community_codes(G, community))


@implements(
    "within_inter_cluster",
    chk.all_of(_COMMUNITY, chk.numeric_or_none("delta")),
    declined="directed graphs; multigraphs; self-loops; nodes without a community; delta <= 0.",
)
def within_inter_cluster(G, ebunch=None, delta=0.001, community="community"):
    if delta is None or delta <= 0:
        raise NotImplementedError("delta must be positive")
    return _scores(G, WITHIN_INTER, ebunch, delta=float(delta), communities=_community_codes(G, community))


# ---- Weisfeiler-Lehman hashes -----------------------------------------------


def _check_wl(a):
    for p in ("iterations", "digest_size"):
        v = a.get(p)
        if not isinstance(v, int) or isinstance(v, bool):
            return f"{p} must be an integer"
    if a["iterations"] <= 0:
        return "iterations must be positive"
    if not 1 <= a["digest_size"] <= 64:
        return "digest_size must be between 1 and 64"
    G = chk.graph_of(a)
    na, ea = a.get("node_attr"), a.get("edge_attr")
    try:
        if na is not None and any(na not in d for d in G._node.values()):
            return "a node lacks node_attr"
        if ea is not None:
            adj = G._succ if G.is_directed() else G._adj
            if any(ea not in d for nbrs in adj.values() for d in nbrs.values()):
                return "an edge lacks edge_attr"
    except TypeError:
        return "attribute names must be hashable"
    return None


def _labels(G, node_attr, edge_attr):
    nl = el = None
    if node_attr is not None:
        nl = [str(d[node_attr]) for d in G.graph._node.values()]
        if not all(s.isascii() for s in nl):
            raise NotImplementedError("non-ASCII labels")
    if edge_attr is not None:
        el = [str(d[edge_attr]) for d in _edge_data(G)]
        if not all(s.isascii() for s in el):
            raise NotImplementedError("non-ASCII labels")
    return nl, el


def _edge_data(G):
    adj = G.graph._succ if G.core.directed else G.graph._adj
    return [adj[u][v] for u, v in _core.edge_keys(G.core)]


_WL = chk.all_of(SIMPLE, _check_wl)


@implements(
    "weisfeiler_lehman_graph_hash",
    _WL,
    declined="multigraphs; missing or non-ASCII attributes; digest_size outside 1..64.",
)
def weisfeiler_lehman_graph_hash(G, edge_attr=None, node_attr=None, iterations=3, digest_size=16):
    nl, el = _labels(G, node_attr, edge_attr)
    if G.core.directed:
        import warnings

        warnings.warn(
            "The hashes produced for directed graphs changed in version v3.5 due to a bugfix to "
            "track in and out edges separately (see documentation).",
            UserWarning,
            stacklevel=2,
        )
    return _core.wl_hash(G.core, nl, el, iterations, digest_size, 0, 0, threads())


@implements(
    "weisfeiler_lehman_subgraph_hashes",
    chk.all_of(_WL, lambda a: None if isinstance(a.get("include_initial_labels"), bool) else "include_initial_labels"),
    declined="multigraphs; missing or non-ASCII attributes; digest_size outside 1..64.",
)
def weisfeiler_lehman_subgraph_hashes(G, edge_attr=None, node_attr=None, iterations=3, digest_size=16,
                                      include_initial_labels=False):
    nl, el = _labels(G, node_attr, edge_attr)
    if G.core.directed:
        import warnings

        warnings.warn(
            "The hashes produced for directed graphs changed in v3.5 due to a bugfix (see documentation).",
            UserWarning,
            stacklevel=2,
        )
    rows = _core.wl_hash(G.core, nl, el, iterations, digest_size, 1, int(include_initial_labels), threads())
    return dict(zip(G.nodes, rows))


# ---- SimRank -----------------------------------------------------------------


def _check_simrank(a):
    G = chk.graph_of(a)
    for p in ("source", "target"):
        v = a.get(p)
        if v is not None:
            try:
                if v not in G:
                    return f"{p} is not a node of the graph"
            except TypeError:
                return f"{p} is not a node of the graph"
    for p in ("importance_factor", "tolerance"):
        v = a.get(p)
        if not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v):
            return f"{p} must be a number"
    m = a.get("max_iterations")
    if not isinstance(m, int) or isinstance(m, bool) or m < 1:
        return "max_iterations must be a positive integer"
    if any("weight" in d for nbrs in G._adj.values() for d in nbrs.values()):
        return "edge weights"
    return None


@implements(
    "simrank_similarity",
    chk.all_of(SIMPLE, _check_simrank),
    declined="multigraphs; graphs with a weight edge attribute; a failure to converge.",
)
def simrank_similarity(G, source=None, target=None, importance_factor=0.9, max_iterations=1000,
                       tolerance=0.0001):
    n = len(G)
    vals = doubles(_core.simrank(G.core, float(importance_factor), max_iterations, float(tolerance),
                                 0, 0, threads()))
    names = G.nodes
    index = G.index
    if source is not None and target is not None:
        return vals[index[source] * n + index[target]]
    if source is not None:
        row = vals[index[source] * n:(index[source] + 1) * n]
        return dict(zip(names, row))
    return {u: dict(zip(names, vals[i * n:(i + 1) * n])) for i, u in enumerate(names)}


# ---- Graph edit distance -----------------------------------------------------


def _two(a):
    G1, G2 = chk.graph_of(a, "G1"), chk.graph_of(a, "G2")
    if G1.is_multigraph() or G2.is_multigraph():
        return "multigraphs are not supported"
    return None


def _check_edit_distance(a):
    reason = _two(a)
    if reason:
        return reason
    for p in ("node_match", "edge_match", "node_subst_cost", "node_del_cost", "node_ins_cost",
              "edge_subst_cost", "edge_del_cost", "edge_ins_cost", "timeout"):
        if a.get(p) is not None:
            return f"{p} is left to NetworkX"
    G1, G2 = chk.graph_of(a, "G1"), chk.graph_of(a, "G2")
    if G1.is_directed() != G2.is_directed():
        return "graphs of different types"
    if nx.number_of_selfloops(G1) or nx.number_of_selfloops(G2):
        return "self-loops are priced differently by the engine"
    r = a.get("roots")
    if r is not None:
        try:
            u, v = r
            if u not in G1 or v not in G2:
                return "roots are not nodes"
        except (TypeError, ValueError):
            return "roots must be a pair of nodes"
    ub = a.get("upper_bound")
    if ub is not None and (not isinstance(ub, (int, float)) or isinstance(ub, bool) or math.isnan(ub)):
        return "upper_bound must be a number"
    return None


@implements(
    "graph_edit_distance",
    _check_edit_distance,
    declined="multigraphs; self-loops; match functions and cost functions; timeout.",
)
def graph_edit_distance(G1, G2, node_match=None, edge_match=None, node_subst_cost=None,
                        node_del_cost=None, node_ins_cost=None, edge_subst_cost=None,
                        edge_del_cost=None, edge_ins_cost=None, roots=None, upper_bound=None,
                        timeout=None):
    if roots is None:
        r1 = r2 = -1
    else:
        r1, r2 = G1.index[roots[0]], G2.index[roots[1]]
    ub = math.inf if upper_bound is None else float(upper_bound)
    return _core.edit_distance(G1.core, G2.core, r1, r2, ub, 0)


# ---- Isomorphism -------------------------------------------------------------

ISO, SUBGRAPH, MONO = 0, 1, 2


def _same_type_simple(a, first="G1", second="G2"):
    reason = _two({"G1": a.get(first), "G2": a.get(second)})
    if reason:
        return reason
    if chk.graph_of(a, first).is_directed() != chk.graph_of(a, second).is_directed():
        return "graphs of different types"
    return None


def _check_is_isomorphic(a):
    if a.get("node_match") is not None or a.get("edge_match") is not None:
        return "match functions are left to NetworkX"
    return _same_type_simple(a)


@implements(
    "is_isomorphic",
    _check_is_isomorphic,
    declined="multigraphs; node_match and edge_match; graphs of different types.",
)
def is_isomorphic(G1, G2, node_match=None, edge_match=None):
    if G1.core.n != G2.core.n or G1.core.m != G2.core.m:
        return False
    if G1.core.n == 0:
        return True
    return _core.vf2pp_count(G1.core, G2.core, None, None, ISO, 1, 1) > 0


def _label_codes(FG, SG, node_label, default_label):
    if node_label is None:
        return None, None
    codes = {}

    def of(G):
        out = []
        for d in G.graph._node.values():
            x = d.get(node_label, default_label)
            try:
                if x != x:
                    raise NotImplementedError("NaN label")
                out.append(codes.setdefault(x, len(codes)))
            except TypeError as err:
                raise NotImplementedError("unhashable label") from err
        return out

    return of(FG), of(SG)


def _check_vf2pp(a):
    reason = _same_type_simple(a, "FG", "SG")
    if reason:
        return reason
    nl = a.get("node_label")
    if nl is not None:
        try:
            hash(nl)
        except TypeError:
            return "node_label must be hashable"
    return None


def _morphism_exists(FG, SG, node_label, default_label, problem):
    if FG.core.n == 0 or SG.core.n == 0:
        return False
    lf, ls = _label_codes(FG, SG, node_label, default_label)
    return _core.vf2pp_count(FG.core, SG.core, lf, ls, problem, 0, 1) > 0


@implements("vf2pp_is_isomorphic", _check_vf2pp, declined="multigraphs; graphs of different types.")
def vf2pp_is_isomorphic(FG, SG, node_label=None, default_label=None):
    return _morphism_exists(FG, SG, node_label, default_label, ISO)


@implements("vf2pp_subgraph_is_isomorphic", _check_vf2pp, declined="multigraphs; graphs of different types.")
def vf2pp_subgraph_is_isomorphic(FG, SG, node_label=None, default_label=None):
    return _morphism_exists(FG, SG, node_label, default_label, SUBGRAPH)


@implements("vf2pp_is_monomorphic", _check_vf2pp, declined="multigraphs; graphs of different types.")
def vf2pp_is_monomorphic(FG, SG, node_label=None, default_label=None):
    return _morphism_exists(FG, SG, node_label, default_label, MONO)


def _properties_check(default):
    def check(a):
        reason = _two(a)
        if reason:
            return reason
        props = a.get("properties", default)
        if not isinstance(props, str):
            return "properties must be a string"
        G1, G2 = chk.graph_of(a, "G1"), chk.graph_of(a, "G2")
        if set(props) & {"t", "c"} and (G1.is_directed() or G2.is_directed()):
            return "triangles and cliques are undirected measures"
        return None

    return check


def _bits(properties):
    p = set(properties)
    return (1 if "d" in p else 0) | (2 if "t" in p else 0) | (4 if "c" in p else 0)


@implements("could_be_isomorphic", _properties_check("dtc"), declined="multigraphs; directed graphs with t or c.")
def could_be_isomorphic(G1, G2, *, properties="dtc"):
    if G1.core.n != G2.core.n:
        return False
    return _core.could_be_isomorphic(G1.core, G2.core, _bits(properties), threads())


@implements("fast_could_be_isomorphic", _properties_check("dt"), declined="multigraphs; directed graphs.")
def fast_could_be_isomorphic(G1, G2):
    return could_be_isomorphic(G1, G2, properties="dt")


@implements("faster_could_be_isomorphic", _properties_check("d"), declined="multigraphs.")
def faster_could_be_isomorphic(G1, G2):
    return could_be_isomorphic(G1, G2, properties="d")
