"""Maximum flow and minimum cut values, minimum-cost flow cost, Gomory-Hu
trees, connectivity, minimum st-cuts, and the component structure of
directed graphs.

Functions that return a flow (a dict of dicts) or a set of paths are left to
NetworkX: many different flows and path sets have the same value. The cuts
returned here are the ones NetworkX picks: the source side is every node
that cannot reach the sink in the residual network, a choice that does not
depend on the order of the search.
"""

from __future__ import annotations

import math

import networkx as nx

from . import _checks as chk
from . import _core
from ._algorithms import ALL, ANY, SIMPLE, _hashable
from ._graph import implements, threads
from ._paths import doubles, int64s

INF = math.inf
NODE, EDGE = 0, 1
SMALLEST, LARGEST = 0, 1
NETWORKX_CUT = 1


# ---- Helpers ----------------------------------------------------------------


def _values(G, attr, default):
    """The values of an edge attribute in the core's edge order, and whether
    every finite value is an int (so that NetworkX's sums are ints)."""
    adj = G.graph._succ if G.core.directed else G.graph._adj
    out = []
    integral = True
    for u, v in _core.edge_keys(G.core):
        x = adj[u][v].get(attr, default)
        if not isinstance(x, (int, float)) or isinstance(x, bool) or x != x:
            raise NotImplementedError("edge attribute is not a number")
        if type(x) is not int and math.isfinite(x):
            integral = False
        out.append(float(x))
    return out, integral


def _capacities(G, capacity):
    caps, integral = _values(G, capacity, INF)
    if any(c < 0 for c in caps):
        raise NotImplementedError("negative capacities are left to NetworkX")
    return caps, integral


def _number(x, integral):
    if not math.isfinite(x):
        raise NotImplementedError("an unbounded flow is left to NetworkX")
    return int(x) if integral else x


def _node_pair(a, s="_s", t="_t"):
    G = chk.graph_of(a)
    u, v = a.get(s), a.get(t)
    try:
        if u not in G or v not in G:
            return "terminals must be nodes of the graph"
    except TypeError:
        return "terminals must be nodes of the graph"
    if u == v:
        return "the terminals are the same node"
    return None


def _no_kwargs(a):
    if a.get("kwargs"):
        return "extra keyword arguments"
    return None


def _capacity_name(a):
    c = a.get("capacity")
    if c is None or isinstance(c, str):
        return None
    return "capacity must be an attribute name"


_FLOW = chk.all_of(SIMPLE, chk.is_none("flow_func"), _no_kwargs, _capacity_name, _node_pair)


def _terminals(G, s, t):
    return [G.index[s]], [G.index[t]]


def _flow_value(G, s, t, capacity):
    caps, integral = _capacities(G, capacity)
    src, snk = _terminals(G, s, t)
    value, _ = _core.max_flow(G.core, src, snk, caps, 0)
    return _number(value, integral)


@implements("maximum_flow_value", _FLOW, declined="multigraphs; flow_func; extra keyword arguments; negative capacities; unbounded flows.")
def maximum_flow_value(flowG, _s, _t, capacity="capacity", flow_func=None, **kwargs):
    return _flow_value(flowG, _s, _t, capacity)


@implements("minimum_cut_value", _FLOW, declined="multigraphs; flow_func; extra keyword arguments; negative capacities; unbounded flows.")
def minimum_cut_value(flowG, _s, _t, capacity="capacity", flow_func=None, **kwargs):
    return _flow_value(flowG, _s, _t, capacity)


@implements("minimum_cut", _FLOW, declined="multigraphs; flow_func; extra keyword arguments; negative capacities; unbounded flows.")
def minimum_cut(flowG, _s, _t, capacity="capacity", flow_func=None, **kwargs):
    caps, integral = _capacities(flowG, capacity)
    src, snk = _terminals(flowG, _s, _t)
    value, side, _ = _core.min_cut(flowG.core, src, snk, caps, LARGEST)
    value = _number(value, integral)
    names = flowG.nodes
    reachable = {v for v, f in zip(names, side) if f}
    return value, (reachable, set(names) - reachable)


# ---- Minimum-cost flow ---------------------------------------------------------


def _check_cost_flow(a):
    G = chk.graph_of(a)
    if not G.is_directed():
        return "NetworkX does not define this for undirected graphs"
    for p in ("demand", "capacity", "weight"):
        if not isinstance(a.get(p), str):
            return f"{p} must be an attribute name"
    if nx.number_of_selfloops(G):
        return "self-loops are handled by NetworkX"
    return None


@implements("min_cost_flow_cost", chk.all_of(SIMPLE, _check_cost_flow),
            declined="undirected graphs; multigraphs; self-loops; infeasible, unbounded or non-numeric inputs.")
def min_cost_flow_cost(G, demand="demand", capacity="capacity", weight="weight"):
    caps, i1 = _capacities(G, capacity)
    costs, i2 = _values(G, weight, 0)
    dem = []
    i3 = True
    for d in G.graph._node.values():
        x = d.get(demand, 0)
        if not isinstance(x, (int, float)) or isinstance(x, bool) or not math.isfinite(x):
            raise NotImplementedError("demand is not a finite number")
        if type(x) is not int:
            i3 = False
        dem.append(float(x))
    if not i3:
        raise NotImplementedError("non-integer demands are left to NetworkX")
    total, _ = _core.min_cost_flow(G.core, dem, caps, costs, 0)
    return _number(total, i1 and i2 and i3)


# ---- Gomory-Hu tree --------------------------------------------------------------


@implements("gomory_hu_tree", chk.all_of(SIMPLE, chk.undirected, chk.is_none("flow_func"), _capacity_name),
            declined="directed graphs; multigraphs; flow_func; negative or infinite capacities.")
def gomory_hu_tree(G, capacity="capacity", flow_func=None):
    caps, integral = _capacities(G, capacity)
    parent, weight = _core.gomory_hu(G.core, caps)
    names = G.nodes
    T = nx.Graph()
    T.add_nodes_from(G.graph)
    w = doubles(weight)
    T.add_weighted_edges_from((names[i], names[parent[i]], _number(w[i], integral))
                              for i in range(1, len(names)))
    return T


# ---- Connectivity -------------------------------------------------------------------


def _check_terminals(a):
    s, t = a.get("s"), a.get("t")
    if s is None and t is None:
        return None
    if s is None or t is None:
        return "both terminals are needed"
    return _node_pair(a, "s", "t")


_CONN = chk.all_of(SIMPLE, chk.is_none("flow_func"), _check_terminals)


def _connectivity(G, s, t, kind):
    if s is None:
        return _core.connectivity(G.core, kind, threads())
    return _core.local_connectivity(G.core, G.index[s], G.index[t], kind)


@implements("node_connectivity", _CONN, declined="multigraphs; flow_func; incomplete or equal terminals.")
def node_connectivity(G, s=None, t=None, flow_func=None):
    return _connectivity(G, s, t, NODE)


@implements("edge_connectivity", chk.all_of(_CONN, chk.is_none("cutoff"), chk.no_selfloops),
            declined="multigraphs; self-loops; flow_func; cutoff; incomplete or equal terminals.")
def edge_connectivity(G, s=None, t=None, flow_func=None, cutoff=None):
    if s is None and G.core.directed and len(G) < 2:
        raise NotImplementedError("NetworkX raises for a directed graph of one node")
    return _connectivity(G, s, t, EDGE)


def _local_check(a):
    for p in ("flow_func", "auxiliary", "residual", "cutoff"):
        if a.get(p) is not None:
            return f"{p} is left to NetworkX"
    return _node_pair(a, "s", "t")


@implements("local_node_connectivity", chk.all_of(SIMPLE, _local_check),
            declined="multigraphs; flow_func, auxiliary, residual and cutoff; equal terminals.")
def local_node_connectivity(G, s, t, flow_func=None, auxiliary=None, residual=None, cutoff=None):
    return _connectivity(G, s, t, NODE)


@implements("local_edge_connectivity", chk.all_of(SIMPLE, _local_check),
            declined="multigraphs; flow_func, auxiliary, residual and cutoff; equal terminals.")
def local_edge_connectivity(G, s, t, flow_func=None, auxiliary=None, residual=None, cutoff=None):
    return _connectivity(G, s, t, EDGE)


@implements("average_node_connectivity", chk.all_of(SIMPLE, chk.is_none("flow_func")),
            declined="multigraphs; flow_func.")
def average_node_connectivity(G, flow_func=None):
    _, avg = _core.all_pairs_connectivity(G.core, threads())
    return 0 if len(G) < 2 else avg


def _check_nbunch(a):
    nb = a.get("nbunch")
    if nb is None:
        return None
    return chk.nodes_all_in_graph("nbunch", allow_empty=True)(a)


@implements("all_pairs_node_connectivity", chk.all_of(SIMPLE, chk.is_none("flow_func"), _check_nbunch),
            declined="multigraphs; flow_func; nbunch with items outside the graph.")
def all_pairs_node_connectivity(G, nbunch=None, flow_func=None):
    mat, _ = _core.all_pairs_connectivity(G.core, threads())
    n = len(G)
    vals = doubles(mat)
    names = G.nodes
    index = G.index
    nodes = list(G.graph) if nbunch is None else list(set(nbunch))
    out = {u: {} for u in nodes}
    directed = G.core.directed
    for i, u in enumerate(nodes):
        for v in nodes[i + 1:] if not directed else nodes:
            if u == v:
                continue
            out[u][v] = int(vals[index[v] * n + index[u]])
            if not directed:
                out[v][u] = out[u][v]
    return out


def _st_cut_check(a):
    for p in ("flow_func", "auxiliary", "residual"):
        if a.get(p) is not None:
            return f"{p} is left to NetworkX"
    return _node_pair(a, "s", "t")


@implements("minimum_st_edge_cut", chk.all_of(SIMPLE, _st_cut_check),
            declined="multigraphs; flow_func, auxiliary and residual; equal terminals.")
def minimum_st_edge_cut(G, s, t, flow_func=None, auxiliary=None, residual=None):
    a, b = _core.min_edge_cut(G.core, G.index[s], G.index[t], threads())
    names = G.nodes
    return {(names[u], names[v]) for u, v in zip(a, b)}


@implements("minimum_st_node_cut", chk.all_of(SIMPLE, _st_cut_check),
            declined="multigraphs; flow_func, auxiliary and residual; equal terminals.")
def minimum_st_node_cut(G, s, t, flow_func=None, auxiliary=None, residual=None):
    names = G.nodes
    return {names[i] for i in _core.min_node_cut(G.core, G.index[s], G.index[t], NETWORKX_CUT, threads())}


def _terminals_given(a):
    s, t = a.get("s"), a.get("t")
    if s is None or t is None:
        return "a global cut depends on NetworkX's iteration order"
    return _node_pair(a, "s", "t")


@implements("minimum_node_cut", chk.all_of(SIMPLE, chk.is_none("flow_func"), _terminals_given),
            declined="multigraphs; flow_func; global cuts (no s and t).")
def minimum_node_cut(G, s=None, t=None, flow_func=None):
    return minimum_st_node_cut(G, s, t)


@implements("minimum_edge_cut", chk.all_of(SIMPLE, chk.is_none("flow_func"), _terminals_given),
            declined="multigraphs; flow_func; global cuts (no s and t).")
def minimum_edge_cut(G, s=None, t=None, flow_func=None):
    return minimum_st_edge_cut(G, s, t)


# ---- Bridges and directed components ----------------------------------------------------


@implements("local_bridges", chk.all_of(SIMPLE, chk.undirected, chk.is_none("weight"), chk.no_selfloops),
            declined="directed graphs; multigraphs; self-loops; weight.")
def local_bridges(G, with_span=True, weight=None):
    flags, span = _core.local_bridges(G.core)
    keys = _core.edge_keys(G.core)
    spans = doubles(span)
    if with_span is not True:
        return ((keys[i][0], keys[i][1]) for i, f in enumerate(flags) if f)
    out = []
    for i, f in enumerate(flags):
        if f:
            s = spans[i]
            out.append((keys[i][0], keys[i][1], s if math.isinf(s) else int(s)))
    return iter(out)


def _scc_order(G):
    if not G.input_order:
        raise NotImplementedError("component order depends on adjacency order")
    return 1


@implements("condensation", chk.all_of(SIMPLE, chk.directed, chk.is_none("scc")),
            declined="undirected graphs; multigraphs; scc.")
def condensation(G, scc=None):
    membership, count = _core.condensation(G.core, _scc_order(G))
    names = G.nodes
    members = [set() for _ in range(count)]
    mapping = {}
    for v, c in zip(names, membership):
        members[c].add(v)
        mapping[v] = c
    C = nx.DiGraph()
    C.graph["mapping"] = mapping
    C.add_nodes_from(range(count))
    C.add_edges_from((mapping[u], mapping[v]) for u, v in G.graph.edges() if mapping[u] != mapping[v])
    nx.set_node_attributes(C, dict(enumerate(members)), "members")
    return C


def _attracting(G):
    order = _scc_order(G)
    membership, count = _core.condensation(G.core, order)
    attract = _core.attracting(G.core, order)
    comps = [set() for _ in range(count)]
    for v, c in zip(G.nodes, membership):
        comps[c].add(v)
    sinks = sorted({c for c in attract if c >= 0})
    return [comps[c] for c in sinks]


_DIRECTED = chk.all_of(SIMPLE, chk.directed)


@implements("attracting_components", _DIRECTED, declined="undirected graphs; multigraphs.")
def attracting_components(G):
    return iter(_attracting(G))


@implements("number_attracting_components", _DIRECTED, declined="undirected graphs; multigraphs.")
def number_attracting_components(G):
    return len(_attracting(G))


@implements("is_attracting_component", _DIRECTED, declined="undirected graphs; multigraphs.")
def is_attracting_component(G):
    ac = _attracting(G)
    return len(ac) == 1 and len(ac[0]) == len(G)


@implements("is_semiconnected", _DIRECTED, declined="undirected graphs; multigraphs.")
def is_semiconnected(G):
    return _core.semiconnected(G.core)


# ---- Components of undirected graphs -------------------------------------------------------------


@implements("k_components", chk.all_of(SIMPLE, chk.undirected, chk.is_none("flow_func")),
            declined="directed graphs; multigraphs; flow_func.")
def k_components(G, flow_func=None):
    lists = _core.k_components(G.core, threads())
    names = G.nodes
    by_k = {}
    for items, k in lists:
        by_k.setdefault(int(k), []).append({names[i] for i in items})
    top = max(by_k, default=0)
    return {k: by_k.get(k, []) for k in range(top, 0, -1)}


def _components_of(G, k):
    membership = _core.k_edge_components(G.core, k, threads())
    groups = {}
    for v, m in zip(G.nodes, membership):
        groups.setdefault(m, set()).add(v)
    return [groups[m] for m in sorted(groups)]


@implements("k_edge_components", chk.all_of(SIMPLE, chk.undirected, chk.integer("k", 1)),
            declined="directed graphs; multigraphs; k below 1. For k of 3 or more the components come in the engine's order.")
def k_edge_components(G, k):
    if k >= 3 and len(G) < 2:
        raise NotImplementedError("NetworkX yields nothing for a single node")
    return iter(_components_of(G, k))


@implements("bridge_components", chk.all_of(SIMPLE, chk.undirected), declined="directed graphs; multigraphs.")
def bridge_components(G):
    return iter(_components_of(G, 2))


@implements("is_k_edge_connected", chk.all_of(SIMPLE, chk.undirected, chk.no_selfloops, chk.integer("k", 1)),
            declined="directed graphs; multigraphs; self-loops; k below 1.")
def is_k_edge_connected(G, k):
    if len(G) < k + 1:
        return False
    if any(d < k for _, d in G.graph.degree()):
        return False
    return _core.connectivity(G.core, EDGE, threads()) >= k


def _terminals_in_graph(a):
    return _node_pair(a, "s", "t")


@implements("is_locally_k_edge_connected", chk.all_of(SIMPLE, chk.undirected, chk.no_selfloops, chk.integer("k", 1), _terminals_in_graph),
            declined="directed graphs; multigraphs; self-loops; k below 1; equal terminals.")
def is_locally_k_edge_connected(G, s, t, k):
    degree = G.graph.degree
    if degree(s) < k or degree(t) < k:
        return False
    return _core.local_connectivity(G.core, G.index[s], G.index[t], EDGE) >= k
