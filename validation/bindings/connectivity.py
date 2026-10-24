"""Engine runners for connectivity (crates/capi/src/connectivity.rs).

Pairs are source 0 and target n - 1. Results that do not exist are
written as -1 (fewer than two nodes; a node cut between adjacent nodes in
the engine convention). Spans of local bridges are written per edge as
the span, -2 for an infinite span, and -1 for an edge that is not a local
bridge. Cut and path results that are not unique are compared through
[size, valid], where valid is 1 when the engine's set really separates
(or its paths really are disjoint paths from s to t); this check is done
here, independently of the engine.

Manifest engine options:
  fn: local_connectivity   [value]; kind: node | edge
  fn: connectivity         [value]; kind: node | edge
  fn: all_pairs_node       n * n values, row-major, zero diagonal
  fn: average_node         [value]
  fn: st_edge_cut          the cut's pairs (u, v)
  fn: st_node_cut          the cut's nodes (sorted); convention
  fn: st_node_cut_check    [size, valid]; convention
  fn: node_cut_check       global, [size, valid]; convention
  fn: edge_cut_check       global, [size, valid]
  fn: disjoint_paths_check [count, valid]; kind
  fn: k_components         each component as k * 1000 + node
  fn: k_edge_components    a label per node; k
  fn: local_bridges        per edge, as above
  fn: condensation         a label per node; order: sorted | input
  fn: attracting           a label per node (-1 outside); order
  fn: semiconnected        [0 or 1]
  fn: attracting_flags     1 if the node is in an attracting component
  fn: strong_partition     strongly connected components (as a partition)
"""

from __future__ import annotations

import ctypes

from engine import ERR_LEN, EngineError, load, register

MISSING = -1.0
KIND = {"node": 0, "edge": 1}
CONVENTION = {"engine": 0, "networkx": 1}
ORDER = {"sorted": 0, "input": 1}

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
    lib.net_local_connectivity.argtypes = [h, i64, i64, i32, i64p, cp, sz]
    lib.net_connectivity.argtypes = [h, i32, i32, i64p, cp, sz]
    lib.net_all_pairs_node_connectivity.argtypes = [h, i32, f64p, f64p, cp, sz]
    lib.net_min_edge_cut.argtypes = [h, i64, i64, i32, i64p, i64p, i64p, cp, sz]
    lib.net_min_node_cut.argtypes = [h, i64, i64, i32, i32, i64p, i64p, cp, sz]
    lib.net_disjoint_paths.argtypes = [h, i64, i64, i32, hp, cp, sz]
    lib.net_k_components.argtypes = [h, i32, hp, cp, sz]
    lib.net_k_edge_components.argtypes = [h, i64, i32, i32p, cp, sz]
    lib.net_local_bridges.argtypes = [h, i32p, f64p, cp, sz]
    lib.net_condensation.argtypes = [h, i32, i32p, i64p, cp, sz]
    lib.net_attracting_components.argtypes = [h, i32, i32p, cp, sz]
    lib.net_is_semiconnected.argtypes = [h, i32p, cp, sz]
    lib.net_lists_count.argtypes = [h]
    lib.net_lists_count.restype = i64
    lib.net_lists_total.argtypes = [h]
    lib.net_lists_total.restype = i64
    lib.net_lists_copy.argtypes = [h, i64p, i64p, f64p, cp, sz]
    lib.net_lists_free.argtypes = [h]
    lib.net_lists_free.restype = None
    for fn in ("net_local_connectivity", "net_connectivity", "net_all_pairs_node_connectivity",
               "net_min_edge_cut", "net_min_node_cut", "net_disjoint_paths", "net_k_components",
               "net_k_edge_components", "net_local_bridges", "net_condensation",
               "net_attracting_components", "net_is_semiconnected", "net_lists_copy"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


def read_lists(lib, handle):
    """(lists, values) from a ListsHandle, which is then freed."""
    try:
        count = lib.net_lists_count(handle)
        total = lib.net_lists_total(handle)
        offsets = (ctypes.c_int64 * (count + 1))()
        items = (ctypes.c_int64 * max(total, 1))()
        values = (ctypes.c_double * max(count, 1))()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(lib.net_lists_copy(handle, offsets, items, values, err, ERR_LEN), err)
        lists = [list(items[offsets[i]:offsets[i + 1]]) for i in range(count)]
        return lists, list(values)[:count]
    finally:
        lib.net_lists_free(handle)


def adjacency(fx):
    """Simple-graph successor sets (neighbour sets if undirected)."""
    out = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a == b:
            continue
        out[a].add(b)
        if not fx.directed:
            out[b].add(a)
    return out


def reaches(adj, s, t, removed=(), removed_arcs=()):
    seen = {s} | set(removed)
    stack = [s]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if (v, w) in removed_arcs:
                continue
            if w == t:
                return True
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return False


def disconnects_nodes(fx, cut):
    """Whether removing `cut` leaves at most one node or a network that is
    not strongly connected (connected, if undirected)."""
    adj = adjacency(fx)
    rest = [v for v in range(fx.n) if v not in set(cut)]
    if len(rest) <= 1:
        return True
    r = rest[0]
    return any(not reaches(adj, r, v, cut) or not reaches(adj, v, r, cut) for v in rest[1:])


def disconnects_arcs(fx, pairs):
    adj = adjacency(fx)
    removed = set(pairs) | ({(b, a) for a, b in pairs} if not fx.directed else set())
    return any(not reaches(adj, 0, v, (), removed) or not reaches(adj, v, 0, (), removed)
               for v in range(1, fx.n))


def _min_edge_cut(g, fx, s, t):
    lib = _declare()
    m = max(len(fx.src), 1)
    a, b = (ctypes.c_int64 * m)(), (ctypes.c_int64 * m)()
    size = ctypes.c_int64()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_min_edge_cut(g.handle, s, t, 0, a, b, ctypes.byref(size), err, ERR_LEN), err)
    return [(a[i], b[i]) for i in range(size.value)]


def _min_node_cut(g, fx, s, t, convention):
    lib = _declare()
    out = (ctypes.c_int64 * max(fx.n, 1))()
    size = ctypes.c_int64()
    err = ctypes.create_string_buffer(ERR_LEN)
    status = lib.net_min_node_cut(g.handle, s, t, CONVENTION[convention], 0, out,
                                  ctypes.byref(size), err, ERR_LEN)
    if status:
        if "adjacent" in err.value.decode():
            return None
        raise EngineError(err.value.decode())
    return list(out[:size.value])


@register("local_connectivity")
def _local(g, spec, fx):
    if fx.n < 2:
        return [MISSING]
    lib = _declare()
    out = ctypes.c_int64()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_local_connectivity(g.handle, 0, fx.n - 1, KIND[spec["kind"]],
                                      ctypes.byref(out), err, ERR_LEN), err)
    return [float(out.value)]


@register("connectivity")
def _global(g, spec, fx):
    lib = _declare()
    out = ctypes.c_int64()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_connectivity(g.handle, KIND[spec["kind"]], 0, ctypes.byref(out), err,
                                ERR_LEN), err)
    return [float(out.value)]


def _all_pairs(g, fx):
    lib = _declare()
    out = (ctypes.c_double * max(fx.n * fx.n, 1))()
    avg = ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_all_pairs_node_connectivity(g.handle, 0, out, ctypes.byref(avg), err,
                                               ERR_LEN), err)
    n = fx.n
    # Column-major to row-major.
    return [out[j * n + i] for i in range(n) for j in range(n)], avg.value


@register("all_pairs_node")
def _run_all_pairs(g, spec, fx):
    return _all_pairs(g, fx)[0]


@register("average_node")
def _run_average(g, spec, fx):
    return [_all_pairs(g, fx)[1]]


@register("st_edge_cut")
def _run_st_edge_cut(g, spec, fx):
    if fx.n < 2:
        return []
    return _min_edge_cut(g, fx, 0, fx.n - 1)


@register("st_node_cut")
def _run_st_node_cut(g, spec, fx):
    if fx.n < 2:
        return [MISSING]
    cut = _min_node_cut(g, fx, 0, fx.n - 1, spec["convention"])
    return [MISSING] if cut is None else sorted(float(v) for v in cut)


@register("st_node_cut_check")
def _run_st_node_cut_check(g, spec, fx):
    if fx.n < 2:
        return [MISSING, MISSING]
    cut = _min_node_cut(g, fx, 0, fx.n - 1, spec["convention"])
    if cut is None:
        return [MISSING, MISSING]
    ok = not reaches(adjacency(fx), 0, fx.n - 1, cut) and 0 not in cut and fx.n - 1 not in cut
    return [float(len(cut)), 1.0 if ok else 0.0]


@register("node_cut_check")
def _run_node_cut_check(g, spec, fx):
    cut = _min_node_cut(g, fx, -1, -1, spec["convention"])
    return [float(len(cut)), 1.0 if disconnects_nodes(fx, cut) else 0.0]


@register("edge_cut_check")
def _run_edge_cut_check(g, spec, fx):
    cut = _min_edge_cut(g, fx, -1, -1)
    ok = fx.n <= 1 or disconnects_arcs(fx, cut)
    return [float(len(cut)), 1.0 if ok else 0.0]


@register("disjoint_paths_check")
def _run_disjoint(g, spec, fx):
    if fx.n < 2:
        return [MISSING, MISSING]
    lib = _declare()
    hd = ctypes.c_void_p()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_disjoint_paths(g.handle, 0, fx.n - 1, KIND[spec["kind"]], ctypes.byref(hd),
                                  err, ERR_LEN), err)
    paths, _ = read_lists(lib, hd)
    adj = adjacency(fx)
    s, t = 0, fx.n - 1
    ok = True
    used_arcs, used_nodes = set(), set()
    for p in paths:
        ok &= len(p) >= 2 and p[0] == s and p[-1] == t and len(set(p)) == len(p)
        for a, b in zip(p, p[1:]):
            ok &= b in adj[a]
            arc = (a, b) if fx.directed else (min(a, b), max(a, b))
            ok &= arc not in used_arcs
            used_arcs.add(arc)
        if spec["kind"] == "node":
            inner = set(p[1:-1])
            ok &= not (inner & used_nodes)
            used_nodes |= inner
    return [float(len(paths)), 1.0 if ok else 0.0]


@register("k_components")
def _run_k_components(g, spec, fx):
    lib = _declare()
    hd = ctypes.c_void_p()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_k_components(g.handle, 0, ctypes.byref(hd), err, ERR_LEN), err)
    lists, values = read_lists(lib, hd)
    return [tuple(float(k * 1000 + v) for v in comp) for comp, k in zip(lists, values)]


@register("k_edge_components")
def _run_k_edge(g, spec, fx):
    lib = _declare()
    out = (ctypes.c_int32 * max(fx.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_k_edge_components(g.handle, spec["k"], 0, out, err, ERR_LEN), err)
    return [float(x) for x in out[:fx.n]]


@register("local_bridges")
def _run_local_bridges(g, spec, fx):
    lib = _declare()
    m = len(fx.src)
    flags = (ctypes.c_int32 * max(m, 1))()
    span = (ctypes.c_double * max(m, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_local_bridges(g.handle, flags, span, err, ERR_LEN), err)
    out = []
    for e in range(m):
        if not flags[e]:
            out.append(-1.0)
        elif span[e] == float("inf"):
            out.append(-2.0)
        else:
            out.append(span[e])
    return out


@register("condensation")
def _run_condensation(g, spec, fx):
    lib = _declare()
    out = (ctypes.c_int32 * max(fx.n, 1))()
    count = ctypes.c_int64()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_condensation(g.handle, ORDER[spec["order"]], out, ctypes.byref(count), err,
                                ERR_LEN), err)
    return [float(x) for x in out[:fx.n]]


@register("attracting")
def _run_attracting(g, spec, fx):
    lib = _declare()
    out = (ctypes.c_int32 * max(fx.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_attracting_components(g.handle, ORDER[spec["order"]], out, err, ERR_LEN), err)
    return [float(x) for x in out[:fx.n]]


@register("semiconnected")
def _run_semiconnected(g, spec, fx):
    lib = _declare()
    out = ctypes.c_int32()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_is_semiconnected(g.handle, ctypes.byref(out), err, ERR_LEN), err)
    return [float(out.value)]


@register("attracting_flags")
def _run_attracting_flags(g, spec, fx):
    return [1.0 if x >= 0 else 0.0 for x in _run_attracting(g, {"order": "sorted"}, fx)]


@register("strong_partition")
def _run_strong_partition(g, spec, fx):
    return _run_condensation(g, {"order": "sorted"}, fx)
