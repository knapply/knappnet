"""Engine runners for the centrality variants (crates/capi/src/centrality.rs):
betweenness for sources and targets, sampled betweenness, percolation and
load centrality, edge load, current-flow measures, resistance distance,
Kemeny's constant, second-order centrality, communicability, dispersion,
VoteRank, Laplacian centrality, trophic levels, and reaching centrality.

Markers: NaN is written as UNDEFINED = -999 and infinity as INFINITE =
-1000 (no value compared here can take either). Matrices are returned as
one tuple per row.

Manifest engine options:
  fn: cx_betweenness_subset  node (`output: nodes`) or edge (`output:
      edges`) values for `sources` and `targets`, each `all`, `even` (even
      node indices), or `thirds` (indices divisible by 3), or, with
      `sample: k`, the k sources NetworkX's betweenness_centrality draws
      with random.Random(`seed`); rescale (unordered | networkx |
      networkx_subset), parallel (distinct | collapse), edge_split
      (proportional | networkx), endpoints, normalized, weighted, tie
  fn: cx_percolation         with states x_v = ((37 v) mod 11) / 10
  fn: cx_load, cx_edge_load  (cutoff, networkx_ties)
  fn: cx_current_flow        node or edge values (`output`), normalized
  fn: cx_current_flow_subset S = even, T = thirds; networkx: true | false
  fn: cx_current_flow_approx the pairs NetworkX's approximation draws with
      random.Random(`seed`) on its reverse Cuthill-McKee order
  fn: cx_current_flow_closeness, cx_resistance (matrix),
      cx_effective_resistance, cx_kemeny, cx_second_order
  fn: cx_communicability (matrix), cx_communicability_betweenness
  fn: cx_dispersion          every ordered pair (u, v), v a distinct
      neighbor of u, in increasing order of (u, v)
  fn: cx_voterank            the chosen nodes, in order
  fn: cx_laplacian           (normalized)
  fn: cx_trophic             what: levels | differences | incoherence
  fn: cx_reaching            what: local | global
  weighted: use weights; weights_as: distance | strength.
"""

from __future__ import annotations

import ctypes
import math
import random

from engine import ERR_LEN, EngineError, load, register

UNDEFINED = -999.0
INFINITE = -1000.0
WEIGHTS_AS = {"distance": 0, "strength": 1}
RESCALE = {"unordered": 0, "networkx": 1, "networkx_subset": 2}
PARALLEL = {"distinct": 0, "collapse": 1}
SPLIT = {"proportional": 0, "networkx": 1}
TIES = {"relative": 0, "exact": 1}

_declared = False


def node_subset(n: int, kind: str) -> list[int]:
    """The node sets named in the manifest; the oracles repeat this."""
    if kind == "all":
        return list(range(n))
    if kind == "even":
        return [v for v in range(n) if v % 2 == 0]
    if kind == "thirds":
        return [v for v in range(n) if v % 3 == 0]
    raise ValueError(kind)


def states(n: int) -> list[float]:
    return [((37 * v) % 11) / 10 for v in range(n)]


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64, u64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double, ctypes.c_uint64
    f64p, i64p = ctypes.POINTER(f64), ctypes.POINTER(i64)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    sig = {
        "net_betweenness_subset": [h, i64p, i64, i64p, i64, i32, i32, i32, i32, i32, i32, i32,
                                   i32, i32, f64p, f64p, cp, sz],
        "net_percolation_centrality": [h, f64p, i32, i32, i32, i32, i32, f64p, cp, sz],
        "net_load_centrality": [h, i32, i32, i32, f64, i32, i32, i32, f64p, cp, sz],
        "net_edge_load": [h, f64, i32, f64p, cp, sz],
        "net_current_flow_betweenness": [h, i32, i32, i32, i32, f64p, f64p, cp, sz],
        "net_current_flow_betweenness_subset": [h, i64p, i64, i64p, i64, i32, i32, i32, i32,
                                                i32, f64p, f64p, cp, sz],
        "net_current_flow_betweenness_approx": [h, i64p, i64, f64, i64, f64, u64, i32, i32,
                                                i32, i32, f64p, cp, sz],
        "net_current_flow_closeness": [h, i32, i32, i32, f64p, cp, sz],
        "net_resistance_distance": [h, i32, i32, i32, f64p, cp, sz],
        "net_effective_graph_resistance": [h, i32, i32, i32, f64p, cp, sz],
        "net_kemeny_constant": [h, i32, f64p, cp, sz],
        "net_second_order_centrality": [h, i32, f64p, cp, sz],
        "net_communicability": [h, f64p, cp, sz],
        "net_communicability_betweenness": [h, i32, f64p, cp, sz],
        "net_dispersion": [h, i64p, i64, i32, f64, f64, f64, f64p, cp, sz],
        "net_voterank": [h, i64, i64p, i64p, cp, sz],
        "net_laplacian_centrality": [h, i32, i32, f64p, cp, sz],
        "net_trophic": [h, i32, i32, i32, f64p, cp, sz],
        "net_reaching_centrality": [h, i32, i32, f64p, f64p, cp, sz],
    }
    for name, args in sig.items():
        fn = getattr(lib, name)
        fn.argtypes = args
        fn.restype = i32
    _declared = True
    return lib


def _mark(x: float) -> float:
    if math.isnan(x):
        return UNDEFINED
    if math.isinf(x):
        return INFINITE
    return x


def _call(fn, handle, size, *args):
    out = (ctypes.c_double * max(size, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    if fn(handle, *args, out, err, ERR_LEN):
        raise EngineError(err.value.decode())
    return [_mark(x) for x in list(out)[:size]]


def _nodes(xs):
    arr = (ctypes.c_int64 * max(len(xs), 1))(*xs)
    return arr, len(xs)


def _weights(spec):
    return int(spec.get("weighted", False)), WEIGHTS_AS[spec.get("weights_as", "distance")]


def _matrix(flat, n):
    return [tuple(flat[j * n + i] for j in range(n)) for i in range(n)]


def sampled_sources(n: int, k: int, seed: int) -> list[int]:
    """The k sources NetworkX's betweenness_centrality(G, k, seed=seed)
    uses on a graph whose nodes are 0..n-1 in order."""
    return random.Random(seed).sample(list(range(n)), k)


@register("cx_betweenness_subset")
def _subset(g, spec, fx):
    lib = _declare()
    n, m = g.n, len(fx.src)
    if "sample" in spec:
        k = min(int(spec["sample"]), n)
        src = sampled_sources(n, k, int(spec.get("seed", 1)))
    else:
        src = node_subset(n, spec.get("sources", "all"))
    tgt = node_subset(n, spec.get("targets", "all"))
    s, ns = _nodes(src)
    t, nt = _nodes(tgt)
    w, wa = _weights(spec)
    nodes = (ctypes.c_double * max(n, 1))()
    edges = (ctypes.c_double * max(m, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    status = lib.net_betweenness_subset(
        g.handle, s, ns, t, nt, w, wa, int(spec.get("endpoints", False)),
        int(spec.get("normalized", False)), RESCALE[spec.get("rescale", "unordered")],
        PARALLEL[spec.get("parallel", "distinct")], SPLIT[spec.get("edge_split", "proportional")],
        TIES[spec.get("tie", "relative")], 0, nodes, edges, err, ERR_LEN)
    if status:
        raise EngineError(err.value.decode())
    if spec.get("output", "nodes") == "edges":
        return [_mark(x) for x in list(edges)[:m]]
    return [_mark(x) for x in list(nodes)[:n]]


@register("cx_percolation")
def _percolation(g, spec, fx):
    lib = _declare()
    x = (ctypes.c_double * max(g.n, 1))(*states(g.n))
    w, wa = _weights(spec)
    return _call(lib.net_percolation_centrality, g.handle, g.n, x, w, wa,
                 PARALLEL[spec.get("parallel", "collapse")], TIES[spec.get("tie", "exact")], 0)


@register("cx_load")
def _load(g, spec, fx):
    lib = _declare()
    w, wa = _weights(spec)
    return _call(lib.net_load_centrality, g.handle, g.n, w, wa,
                 int(spec.get("normalized", False)), float(spec.get("cutoff", math.inf)),
                 int(spec.get("networkx_ties", False)), TIES[spec.get("tie", "exact")], 0)


@register("cx_edge_load")
def _edge_load(g, spec, fx):
    lib = _declare()
    m = len(fx.src)
    flat = _call(lib.net_edge_load, g.handle, 2 * m, float(spec.get("cutoff", math.inf)), 0)
    return [(flat[2 * e], flat[2 * e + 1]) for e in range(m)]


def _flow(spec):
    w, wa = int(spec.get("weighted", False)), WEIGHTS_AS[spec.get("weights_as", "strength")]
    return w, wa


@register("cx_current_flow")
def _current_flow(g, spec, fx):
    lib = _declare()
    n, m = g.n, len(fx.src)
    nodes = (ctypes.c_double * max(n, 1))()
    edges = (ctypes.c_double * max(m, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    if lib.net_current_flow_betweenness(g.handle, *_flow(spec), int(spec.get("normalized", True)),
                                        0, nodes, edges, err, ERR_LEN):
        raise EngineError(err.value.decode())
    if spec.get("output", "nodes") == "edges":
        return [_mark(x) for x in list(edges)[:m]]
    return [_mark(x) for x in list(nodes)[:n]]


@register("cx_current_flow_subset")
def _current_flow_subset(g, spec, fx):
    lib = _declare()
    n, m = g.n, len(fx.src)
    s, ns = _nodes(node_subset(n, spec.get("sources", "even")))
    t, nt = _nodes(node_subset(n, spec.get("targets", "thirds")))
    nodes = (ctypes.c_double * max(n, 1))()
    edges = (ctypes.c_double * max(m, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    if lib.net_current_flow_betweenness_subset(
            g.handle, s, ns, t, nt, int(spec.get("networkx", False)), *_flow(spec),
            int(spec.get("normalized", True)), 0, nodes, edges, err, ERR_LEN):
        raise EngineError(err.value.decode())
    if spec.get("output", "nodes") == "edges":
        return [_mark(x) for x in list(edges)[:m]]
    return [_mark(x) for x in list(nodes)[:n]]


def approximation_pairs(fx, seed: int, epsilon: float, sample_weight: float = 1.0):
    """The pairs NetworkX's approximate_current_flow_betweenness_centrality
    draws (in original node labels): random.Random(seed).sample(range(n), 2)
    k times on its reverse Cuthill-McKee relabeling."""
    import networkx as nx
    from networkx.utils import reverse_cuthill_mckee_ordering

    g = nx.Graph()
    g.add_nodes_from(range(fx.n))
    g.add_edges_from((a, b) for a, b in zip(fx.src, fx.tgt) if a != b)
    n = fx.n
    if n < 3:
        return []
    ordering = list(reverse_cuthill_mckee_ordering(g))
    nb = (n - 1.0) * (n - 2.0)
    cstar = n * (n - 1) / nb
    k = int(sample_weight * math.ceil((cstar / epsilon) ** 2 * math.log(n)))
    rng = random.Random(seed)
    pairs = []
    for _ in range(k):
        s, t = rng.sample(range(n), 2)
        pairs.append((ordering[s], ordering[t]))
    return pairs


@register("cx_current_flow_approx")
def _current_flow_approx(g, spec, fx):
    lib = _declare()
    pairs = approximation_pairs(fx, int(spec.get("seed", 1)), float(spec.get("epsilon", 0.5)))
    flat = [x for p in pairs for x in p]
    arr = (ctypes.c_int64 * max(len(flat), 1))(*flat)
    return _call(lib.net_current_flow_betweenness_approx, g.handle, g.n, arr, len(pairs), 0.5,
                 10000, 1.0, 0, *_flow(spec), int(spec.get("normalized", True)), 0)


@register("cx_current_flow_closeness")
def _cfc(g, spec, fx):
    lib = _declare()
    return _call(lib.net_current_flow_closeness, g.handle, g.n, *_flow(spec), 0)


@register("cx_resistance")
def _resistance(g, spec, fx):
    lib = _declare()
    w, wa = int(spec.get("weighted", False)), WEIGHTS_AS[spec.get("weights_as", "distance")]
    return _matrix(_call(lib.net_resistance_distance, g.handle, g.n * g.n, w, wa, 0), g.n)


@register("cx_effective_resistance")
def _effective(g, spec, fx):
    lib = _declare()
    w, wa = int(spec.get("weighted", False)), WEIGHTS_AS[spec.get("weights_as", "distance")]
    return _call(lib.net_effective_graph_resistance, g.handle, 1, w, wa, 0)


@register("cx_kemeny")
def _kemeny(g, spec, fx):
    lib = _declare()
    return _call(lib.net_kemeny_constant, g.handle, 1, int(spec.get("weighted", False)))


@register("cx_second_order")
def _second_order(g, spec, fx):
    lib = _declare()
    return _call(lib.net_second_order_centrality, g.handle, g.n,
                 int(spec.get("weighted", False)))


@register("cx_communicability")
def _communicability(g, spec, fx):
    lib = _declare()
    return _matrix(_call(lib.net_communicability, g.handle, g.n * g.n), g.n)


@register("cx_communicability_betweenness")
def _cbc(g, spec, fx):
    lib = _declare()
    return _call(lib.net_communicability_betweenness, g.handle, g.n, 0)


def dispersion_pairs(fx) -> list[tuple[int, int]]:
    """Every ordered pair (u, v), v a distinct (out-)neighbor of u, sorted."""
    pairs = set()
    for a, b in zip(fx.src, fx.tgt):
        pairs.add((a, b))
        if not fx.directed:
            pairs.add((b, a))
    return sorted(pairs)


@register("cx_dispersion")
def _dispersion(g, spec, fx):
    lib = _declare()
    pairs = dispersion_pairs(fx)
    flat = [x for p in pairs for x in p]
    arr = (ctypes.c_int64 * max(len(flat), 1))(*flat)
    return _call(lib.net_dispersion, g.handle, len(pairs), arr, len(pairs),
                 int(spec.get("normalized", True)), float(spec.get("alpha", 1.0)),
                 float(spec.get("b", 0.0)), float(spec.get("c", 0.0)))


@register("cx_voterank")
def _voterank(g, spec, fx):
    lib = _declare()
    out = (ctypes.c_int64 * max(g.n, 1))()
    count = ctypes.c_int64()
    err = ctypes.create_string_buffer(ERR_LEN)
    if lib.net_voterank(g.handle, int(spec.get("count", -1)), out, ctypes.byref(count), err,
                        ERR_LEN):
        raise EngineError(err.value.decode())
    return [float(x) for x in list(out)[:count.value]]


@register("cx_laplacian")
def _laplacian(g, spec, fx):
    lib = _declare()
    return _call(lib.net_laplacian_centrality, g.handle, g.n, int(spec.get("weighted", False)),
                 int(spec.get("normalized", True)))


@register("cx_trophic")
def _trophic(g, spec, fx):
    lib = _declare()
    what = {"levels": 0, "differences": 1, "incoherence": 2}[spec["what"]]
    size = {0: g.n, 1: len(fx.src), 2: 1}[what]
    return _call(lib.net_trophic, g.handle, size, what, int(spec.get("weighted", False)),
                 int(spec.get("cannibalism", False)))


@register("cx_reaching")
def _reaching(g, spec, fx):
    lib = _declare()
    local = (ctypes.c_double * max(g.n, 1))()
    glob = ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    if lib.net_reaching_centrality(g.handle, int(spec.get("weighted", False)),
                                   int(spec.get("normalized", True)), local, ctypes.byref(glob),
                                   err, ERR_LEN):
        raise EngineError(err.value.decode())
    if spec.get("what", "local") == "global":
        return [_mark(glob.value)]
    return [_mark(x) for x in list(local)[:g.n]]
