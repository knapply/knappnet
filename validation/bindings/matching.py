"""Engine runners for matchings, greedy colouring, random maximal
independent sets, square clustering, and generalized degree
(crates/capi/src/matching.rs).

Derived inputs, which the oracles build the same way:
  nx_order(fx)   the fixture's distinct non-loop pairs in the order of
                 NetworkX's Graph.edges() for a graph built from the file
                 (nodes 0..n - 1, then the edges in file order); the
                 engine is given its edges in this order for
                 maximal_matching
  candidate      every other edge of the file (positions 0, 2, 4, ...),
                 as pairs, for the matching tests

Encodings: a matching summary is (size, total weight, 1 if valid), size
-1 for weighted maximum-weight matchings (optimal ones may differ in size); a
matching is one (u, v) line per edge with u < v; colourings and square
clustering are node-aligned; generalized degree is one line per node
holding k1 c1 k2 c2 ... (empty for no neighbours); checks are 0 or 1.

Manifest engine options:
  fn: max_weight_matching   weighted, max_cardinality
  fn: min_weight_matching   weighted
  fn: maximal_matching      edges given in nx_order
  fn: check_matching        pairs (candidate | maximal)
  fn: greedy_color          strategy
  fn: mis_check             (independent, maximal, contains node 0's start)
  fn: square_clustering
  fn: generalized_degree
"""

from __future__ import annotations

import ctypes

from engine import ERR_LEN, EngineError, Graph, load, register

STRATEGY = {"largest_first": 0, "random_sequential": 1, "smallest_last": 2,
            "independent_set": 3, "connected_sequential_bfs": 4,
            "connected_sequential_dfs": 5, "saturation_largest_first": 6}

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64, u64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double, ctypes.c_uint64
    f64p, i64p, i32p = ctypes.POINTER(f64), ctypes.POINTER(i64), ctypes.POINTER(i32)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    lib.net_max_weight_matching.argtypes = [h, i32, i32, i64p, i64p, cp, sz]
    lib.net_min_weight_matching.argtypes = [h, i32, i64p, i64p, cp, sz]
    lib.net_maximal_matching.argtypes = [h, i64p, i64p, cp, sz]
    lib.net_check_matching.argtypes = [h, i64p, i64p, i64, i32p, cp, sz]
    lib.net_greedy_color.argtypes = [h, i32, u64, i64p, cp, sz]
    lib.net_maximal_independent_set.argtypes = [h, i64p, i64, u64, i64p, i64p, cp, sz]
    lib.net_square_clustering.argtypes = [h, i32, f64p, cp, sz]
    lib.net_generalized_degree.argtypes = [h, i32, i64p, i64p, i64p, i64, i64p, cp, sz]
    for fn in ("net_max_weight_matching", "net_min_weight_matching", "net_maximal_matching",
               "net_check_matching", "net_greedy_color", "net_maximal_independent_set",
               "net_square_clustering", "net_generalized_degree"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


def nx_order(fx):
    adj = {v: {} for v in range(fx.n)}
    for a, b in zip(fx.src, fx.tgt):
        adj[a][b] = True
        adj[b][a] = True
    seen, out = set(), []
    for u in range(fx.n):
        for v in adj[u]:
            if v not in seen and u != v:
                out.append((u, v))
        seen.add(u)
    return out


def weights_of(fx):
    """The best (largest) weight of every pair, and the smallest."""
    hi, lo = {}, {}
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        k = (min(a, b), max(a, b))
        w = fx.weights[i] if fx.weighted else 1.0
        hi[k] = max(hi.get(k, w), w)
        lo[k] = min(lo.get(k, w), w)
    return hi, lo


def _edges(lib, g, fn, *args):
    err = ctypes.create_string_buffer(ERR_LEN)
    m = max(len(g.fx.src), 1)
    out = (ctypes.c_int64 * m)()
    found = ctypes.c_int64()
    _check(fn(g.handle, *args, out, ctypes.byref(found), err, ERR_LEN), err)
    return [out[i] for i in range(found.value)]


def _checks(lib, g, n, pairs):
    err = ctypes.create_string_buffer(ERR_LEN)
    us = (ctypes.c_int64 * max(len(pairs), 1))(*[p[0] for p in pairs])
    vs = (ctypes.c_int64 * max(len(pairs), 1))(*[p[1] for p in pairs])
    out = (ctypes.c_int32 * 3)()
    _check(lib.net_check_matching(g.handle, us, vs, len(pairs), out, err, ERR_LEN), err)
    return list(out)


def _summary(lib, g, fx, edges, weighted, low=False):
    hi, lo = weights_of(fx)
    tab = lo if low else hi
    pairs = [(min(fx.src[e], fx.tgt[e]), max(fx.src[e], fx.tgt[e])) for e in edges]
    total = sum(tab[p] if weighted else 1.0 for p in pairs)
    return [(len(pairs), total, _checks(lib, g, fx.n, pairs)[0])]


@register("max_weight_matching")
def _mwm(g, spec, fx):
    lib = _declare()
    g.fx = fx
    e = _edges(lib, g, lib.net_max_weight_matching, int(bool(spec.get("weighted"))),
               int(bool(spec.get("max_cardinality"))))
    out = _summary(lib, g, fx, e, bool(spec.get("weighted")))
    if spec.get("weighted") and not spec.get("max_cardinality"):
        # Heaviest matchings of different sizes can tie; compare weights.
        out = [(-1,) + out[0][1:]]
    return out


@register("min_weight_matching")
def _minwm(g, spec, fx):
    lib = _declare()
    g.fx = fx
    e = _edges(lib, g, lib.net_min_weight_matching, int(bool(spec.get("weighted"))))
    return _summary(lib, g, fx, e, bool(spec.get("weighted")), low=True)


def _maximal_pairs(lib, fx):
    order = nx_order(fx)
    h = Graph(fx.n, [a for a, _ in order], [b for _, b in order], None, False)
    h.fx = type("F", (), {"src": [a for a, _ in order]})
    e = _edges(lib, h, lib.net_maximal_matching)
    return [order[i] for i in e], h


@register("maximal_matching")
def _maximal(g, spec, fx):
    lib = _declare()
    pairs, _ = _maximal_pairs(lib, fx)
    return [(min(p), max(p)) for p in pairs]


@register("check_matching")
def _check_matching(g, spec, fx):
    lib = _declare()
    if spec["pairs"] == "candidate":
        pairs = [(fx.src[i], fx.tgt[i]) for i in range(0, len(fx.src), 2)]
    else:
        pairs, _ = _maximal_pairs(lib, fx)
    return _checks(lib, g, fx.n, pairs)


@register("greedy_color")
def _greedy_color(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    out = (ctypes.c_int64 * max(fx.n, 1))()
    _check(lib.net_greedy_color(g.handle, STRATEGY[spec["strategy"]], 1, out, err, ERR_LEN), err)
    return [out[i] for i in range(fx.n)]


@register("mis_check")
def _mis_check(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    if fx.n == 0:
        return [-1]
    looped = {a for a, b in zip(fx.src, fx.tgt) if a == b}
    start = [v for v in range(fx.n) if v not in looped][:1]
    sa = (ctypes.c_int64 * 1)(*start) if start else None
    nodes = (ctypes.c_int64 * fx.n)()
    found = ctypes.c_int64()
    _check(lib.net_maximal_independent_set(g.handle, sa, len(start), 11, nodes,
                                           ctypes.byref(found), err, ERR_LEN), err)
    got = {nodes[i] for i in range(found.value)}
    nbrs = {v: set() for v in range(fx.n)}
    for a, b in zip(fx.src, fx.tgt):
        nbrs[a].add(b)
        nbrs[b].add(a)
    independent = all(not (nbrs[v] & got - {v}) for v in got)
    maximal = all(v in got or nbrs[v] & got for v in range(fx.n))
    return [(int(independent), int(maximal), int(set(start) <= got))]


@register("square_clustering")
def _square_clustering(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    out = (ctypes.c_double * max(fx.n, 1))()
    _check(lib.net_square_clustering(g.handle, 0, out, err, ERR_LEN), err)
    return [out[i] for i in range(fx.n)]


@register("generalized_degree")
def _generalized_degree(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    cap = max(2 * len(fx.src), 1)
    off = (ctypes.c_int64 * (fx.n + 1))()
    ks = (ctypes.c_int64 * cap)()
    cs = (ctypes.c_int64 * cap)()
    total = ctypes.c_int64()
    _check(lib.net_generalized_degree(g.handle, 0, off, ks, cs, cap, ctypes.byref(total), err,
                                      ERR_LEN), err)
    out = []
    for v in range(fx.n):
        row = []
        for i in range(off[v], off[v + 1]):
            row += [ks[i], cs[i]]
        out.append(tuple(row))
    return out
