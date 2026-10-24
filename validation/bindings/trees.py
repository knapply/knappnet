"""Engine runners for spanning trees, branchings and arborescences, and
directed acyclic graphs (crates/capi/src/trees.rs).

Derived networks, which the oracles build the same way:
  dag(fx)       the arcs a -> b with a < b of the fixture, in file order
                (an undirected edge {a, b} becomes a -> b), on all n nodes:
                a DAG for every fixture
  bfs_tree(fx)  the breadth-first tree from node 0 over the fixture read
                as undirected (neighbours in increasing order), arcs from
                parent to child
  nx_order(fx)  the fixture's edges in the order NetworkX's
                MultiGraph.edges() reports them for a graph built from the
                file (nodes 0..n - 1, then the edges in file order)
Weights are the fixture's (with `weighted`) or 1.

Encodings: a missing result (no spanning arborescence, no DAG) is -1;
edges are (u, v, weight) with u <= v for undirected networks; a node set
or order is one tuple per line.

Manifest engine options:
  fn: spanning_tree      the forest's edges; maximum, algorithm, weighted,
                         nx_order (feed the edges in NetworkX's order)
  fn: spanning_weight    [weight]; maximum, weighted
  fn: branching          [weight] (-1 if none); kind, weighted
  fn: topological_sort   the order; order: sorted | input; on dag(fx)
  fn: generations        each node's generation; on dag(fx)
  fn: lexicographic      the order; on dag(fx)
  fn: all_topological    the first 50 orders, one per line; on dag(fx)
  fn: is_dag             [0 or 1]; on the fixture
  fn: closure            the pairs (u, v); reflexive; on the fixture
  fn: reduction          the pairs (u, v); on dag(fx)
  fn: longest_path       the path's nodes; weighted; on dag(fx)
  fn: longest_length     [length]; weighted; on dag(fx)
  fn: idom               each node's immediate dominator from node 0
  fn: frontiers          one tuple per node (-1 if unreachable)
  fn: dag_lca            each pair u <= v's ancestor (-1 if none), pairs
                         in row order; on dag(fx), n <= 40
  fn: dag_lca_check      [pairs answered, 1 if each answer is a lowest
                         common ancestor]; on dag(fx)
  fn: tree_lca           each pair u <= v of tree nodes; on bfs_tree(fx)
"""

from __future__ import annotations

import ctypes
from collections import deque

from engine import ERR_LEN, EngineError, Graph, load, register

MISSING = -1.0
ALGORITHM = {"kruskal": 0, "prim": 1, "boruvka": 2}
KIND = {"max_branching": 0, "min_branching": 1, "minimal": 2, "max_arborescence": 3,
        "min_arborescence": 4}
ORDER = {"sorted": 0, "input": 1}
REFLEXIVE = {"cycles": 0, "never": 1, "always": 2}
LCA_MAX_NODES = 40

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    f64p, i64p, i32p = ctypes.POINTER(f64), ctypes.POINTER(i64), ctypes.POINTER(i32)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    hp = ctypes.POINTER(h)
    lib.net_spanning_tree.argtypes = [h, i32, i32, i32, i64p, i64p, f64p, cp, sz]
    lib.net_branching.argtypes = [h, i32, i32, i64, i64p, i64p, f64p, cp, sz]
    lib.net_topological_sort.argtypes = [h, i32, i64p, i64p, cp, sz]
    lib.net_lexicographical_topological_sort.argtypes = [h, i64p, cp, sz]
    lib.net_all_topological_sorts.argtypes = [h, i32, i64, hp, cp, sz]
    lib.net_is_dag.argtypes = [h, i32p, cp, sz]
    lib.net_transitive_closure.argtypes = [h, i32, hp, cp, sz]
    lib.net_transitive_reduction.argtypes = [h, hp, cp, sz]
    lib.net_dag_longest_path.argtypes = [h, i32, i32, i64p, i64p, f64p, cp, sz]
    lib.net_immediate_dominators.argtypes = [h, i64, i64p, cp, sz]
    lib.net_dominance_frontiers.argtypes = [h, i64, hp, cp, sz]
    lib.net_lowest_common_ancestors.argtypes = [h, i64p, i64p, i64, i64, i32, i64p, cp, sz]
    for fn in ("net_spanning_tree", "net_branching", "net_topological_sort",
               "net_lexicographical_topological_sort", "net_all_topological_sorts",
               "net_is_dag", "net_transitive_closure", "net_transitive_reduction",
               "net_dag_longest_path", "net_immediate_dominators", "net_dominance_frontiers",
               "net_lowest_common_ancestors"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


def _lists(lib, handle):
    from bindings_connectivity import read_lists
    return read_lists(lib, handle)


def weights_of(fx, weighted):
    return list(fx.weights) if weighted and fx.weighted else [1.0] * len(fx.src)


def dag_edges(fx):
    """(src, tgt, positions) of the derived DAG, in file order."""
    keep = [e for e, (a, b) in enumerate(zip(fx.src, fx.tgt)) if a < b]
    return [fx.src[e] for e in keep], [fx.tgt[e] for e in keep], keep


def dag_graph(fx, weighted=False):
    s, t, keep = dag_edges(fx)
    w = [fx.weights[e] for e in keep] if weighted and fx.weighted else None
    return Graph(fx.n, s, t, w, True)


def bfs_tree(fx):
    adj = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            adj[a].add(b)
            adj[b].add(a)
    parent = {0: None} if fx.n else {}
    q = deque([0] if fx.n else [])
    arcs = []
    while q:
        v = q.popleft()
        for w in sorted(adj[v]):
            if w not in parent:
                parent[w] = v
                arcs.append((v, w))
                q.append(w)
    return arcs, sorted(parent)


def nx_order(fx):
    """Positions of the fixture's edges in NetworkX MultiGraph.edges() order."""
    import networkx as nx
    g = (nx.MultiDiGraph if fx.directed else nx.MultiGraph)()
    g.add_nodes_from(range(fx.n))
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, key=e)
    return [k for _, _, k in g.edges(keys=True)]


def _edge_tuple(fx, e, w):
    a, b = fx.src[e], fx.tgt[e]
    if not fx.directed:
        a, b = min(a, b), max(a, b)
    return (a, b, w[e])


@register("spanning_tree")
def _run_spanning_tree(g, spec, fx):
    lib = _declare()
    weighted = spec.get("weighted", False) and fx.weighted
    w = weights_of(fx, weighted)
    order = nx_order(fx) if spec.get("nx_order") else list(range(len(fx.src)))
    h = Graph(fx.n, [fx.src[e] for e in order], [fx.tgt[e] for e in order],
              [w[e] for e in order] if weighted else None, fx.directed)
    edges = (ctypes.c_int64 * max(fx.n, 1))()
    size, weight = ctypes.c_int64(), ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_spanning_tree(h.handle, int(weighted), int(spec.get("maximum", False)),
                                 ALGORITHM[spec.get("algorithm", "kruskal")], edges,
                                 ctypes.byref(size), ctypes.byref(weight), err, ERR_LEN), err)
    return [_edge_tuple(fx, order[edges[i]], w) for i in range(size.value)]


@register("spanning_weight")
def _run_spanning_weight(g, spec, fx):
    lib = _declare()
    edges = (ctypes.c_int64 * max(fx.n, 1))()
    size, weight = ctypes.c_int64(), ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_spanning_tree(g.handle, int(spec.get("weighted", False) and fx.weighted),
                                 int(spec.get("maximum", False)), 0, edges, ctypes.byref(size),
                                 ctypes.byref(weight), err, ERR_LEN), err)
    return [weight.value]


@register("branching")
def _run_branching(g, spec, fx):
    lib = _declare()
    edges = (ctypes.c_int64 * max(fx.n, 1))()
    size, weight = ctypes.c_int64(), ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    weighted = spec.get("weighted", False) and fx.weighted
    status = lib.net_branching(g.handle, int(weighted), KIND[spec["kind"]], -1,
                               edges, ctypes.byref(size), ctypes.byref(weight), err, ERR_LEN)
    if status:
        if "no spanning arborescence" in err.value.decode():
            return [MISSING]
        raise EngineError(err.value.decode())
    return [weight.value]


def _topo(g, fx, order):
    lib = _declare()
    nodes = (ctypes.c_int64 * max(fx.n, 1))()
    gen = (ctypes.c_int64 * max(fx.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_topological_sort(g.handle, ORDER[order], nodes, gen, err, ERR_LEN), err)
    return list(nodes[:fx.n]), list(gen[:fx.n])


@register("topological_sort")
def _run_topological_sort(g, spec, fx):
    return [float(v) for v in _topo(dag_graph(fx), fx, spec["order"])[0]]


@register("generations")
def _run_generations(g, spec, fx):
    return [float(v) for v in _topo(dag_graph(fx), fx, spec.get("order", "input"))[1]]


@register("lexicographic")
def _run_lexicographic(g, spec, fx):
    lib = _declare()
    h = dag_graph(fx)
    nodes = (ctypes.c_int64 * max(fx.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_lexicographical_topological_sort(h.handle, nodes, err, ERR_LEN), err)
    return [float(v) for v in nodes[:fx.n]]


@register("all_topological")
def _run_all_topological(g, spec, fx):
    lib = _declare()
    h = dag_graph(fx)
    hd = ctypes.c_void_p()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_all_topological_sorts(h.handle, ORDER[spec.get("order", "input")], 50,
                                         ctypes.byref(hd), err, ERR_LEN), err)
    lists, _ = _lists(lib, hd)
    return [tuple(float(v) for v in o) for o in lists]


@register("is_dag")
def _run_is_dag(g, spec, fx):
    lib = _declare()
    out = ctypes.c_int32()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_is_dag(g.handle, ctypes.byref(out), err, ERR_LEN), err)
    return [float(out.value)]


def _pairs(lib, handle):
    lists, _ = _lists(lib, handle)
    return [(float(u), float(v)) for u, vs in enumerate(lists) for v in vs]


@register("closure")
def _run_closure(g, spec, fx):
    lib = _declare()
    hd = ctypes.c_void_p()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_transitive_closure(g.handle, REFLEXIVE[spec["reflexive"]], ctypes.byref(hd),
                                      err, ERR_LEN), err)
    return _pairs(lib, hd)


@register("reduction")
def _run_reduction(g, spec, fx):
    lib = _declare()
    h = dag_graph(fx)
    hd = ctypes.c_void_p()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_transitive_reduction(h.handle, ctypes.byref(hd), err, ERR_LEN), err)
    return _pairs(lib, hd)


def _longest(fx, weighted):
    weighted = weighted and fx.weighted
    lib = _declare()
    h = dag_graph(fx, weighted)
    path = (ctypes.c_int64 * max(fx.n, 1))()
    size, length = ctypes.c_int64(), ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_dag_longest_path(h.handle, int(weighted), 1, path, ctypes.byref(size),
                                    ctypes.byref(length), err, ERR_LEN), err)
    return [float(v) for v in path[:size.value]], length.value


@register("longest_path")
def _run_longest_path(g, spec, fx):
    return _longest(fx, spec.get("weighted", False))[0]


@register("longest_length")
def _run_longest_length(g, spec, fx):
    return [_longest(fx, spec.get("weighted", False))[1]]


@register("idom")
def _run_idom(g, spec, fx):
    if fx.n == 0:
        return []
    lib = _declare()
    out = (ctypes.c_int64 * fx.n)()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_immediate_dominators(g.handle, 0, out, err, ERR_LEN), err)
    return [float(v) for v in out]


@register("frontiers")
def _run_frontiers(g, spec, fx):
    lib = _declare()
    hd = ctypes.c_void_p()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_dominance_frontiers(g.handle, 0, ctypes.byref(hd), err, ERR_LEN), err)
    lists, values = _lists(lib, hd)
    return [tuple(float(v) for v in l) if r else (MISSING,) for l, r in zip(lists, values)]


def _lca(h, pairs, root, order="sorted"):
    lib = _declare()
    k = len(pairs)
    a = (ctypes.c_int64 * max(k, 1))(*[p[0] for p in pairs])
    b = (ctypes.c_int64 * max(k, 1))(*[p[1] for p in pairs])
    out = (ctypes.c_int64 * max(k, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_lowest_common_ancestors(h.handle, a, b, k, root, ORDER[order], out, err,
                                           ERR_LEN), err)
    return list(out[:k])


def all_pairs(nodes):
    return [(u, v) for i, u in enumerate(nodes) for v in nodes[i:]]


@register("dag_lca")
def _run_dag_lca(g, spec, fx):
    h = dag_graph(fx)
    return [float(x) for x in _lca(h, all_pairs(list(range(fx.n))), -2,
                                   spec.get("order", "input"))]


@register("dag_lca_check")
def _run_dag_lca_check(g, spec, fx):
    h = dag_graph(fx)
    pairs = all_pairs(list(range(fx.n)))
    answers = _lca(h, pairs, -2)
    s, t, _ = dag_edges(fx)
    succ = [set() for _ in range(fx.n)]
    for a, b in zip(s, t):
        succ[a].add(b)
    # Ancestors (inclusive) by repeated relaxation in topological order.
    anc = [{v} for v in range(fx.n)]
    for v in range(fx.n):  # arcs go from lower to higher index
        for w in succ[v]:
            anc[w] |= anc[v]
    ok, answered = True, 0
    for (u, v), x in zip(pairs, answers):
        common = anc[u] & anc[v]
        if x < 0:
            ok &= not common
            continue
        answered += 1
        ok &= x in common and not (succ[x] & common)
    return [float(answered), 1.0 if ok else 0.0]


@register("tree_lca")
def _run_tree_lca(g, spec, fx):
    arcs, nodes = bfs_tree(fx)
    h = Graph(fx.n, [a for a, _ in arcs], [b for _, b in arcs], None, True)
    if fx.n == 0:
        return []
    return [float(x) for x in _lca(h, all_pairs(nodes), 0)]
