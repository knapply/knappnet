"""Engine runners for Weisfeiler-Lehman hashes, SimRank, Panther, and
graph edit distance (crates/capi/src/similarity.rs).

Derived inputs, which the oracles build the same way:
  labels         node labels str(v % 3), edge labels str(e % 2) for edge e
                 in file order
  hub(fx)        the Panther source: the node of largest degree in the
                 simple undirected graph (lowest index on ties)
  ged_pair(fx)   G2 for the edit distance: the fixture without its first
                 edge, plus edges (0, n - 1) (if absent) and (n, min(1, n - 1)), on
                 n + 1 nodes; node labels v % 2, edge labels by position
                 % 2 in each graph's edge list

Encodings: a hash is its 32 (or 2 digest_size) hexadecimal digits cut into
8-digit integers, one tuple per hash; subgraph hashes are one tuple per
node holding all its hashes; a matrix is its entries in row-major order;
Panther scores are node-aligned (-1 for the source and removed nodes);
an engine error on a fixture an entry declares undefined is [-1].

Manifest engine options:
  fn: wl_hash            iterations, digest_size, node_labels, edge_labels
  fn: wl_subgraph        the same and include_initial
  fn: simrank            stop (allclose | relative), weighted, tolerance
  fn: panther            eps, seed (k = every other node)
  fn: ged                labels (none | nodes | both)
"""

from __future__ import annotations

import ctypes

from engine import ERR_LEN, EngineError, Graph, load, register

STOP = {"allclose": 0, "relative": 1}
GED_MAX_NODES = 10

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64, u64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double, ctypes.c_uint64
    f64p, i64p, i32p = ctypes.POINTER(f64), ctypes.POINTER(i64), ctypes.POINTER(i32)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    cpp = ctypes.POINTER(ctypes.c_char_p)
    lib.net_wl_graph_hash.argtypes = [h, cpp, cpp, i64, i64, i32, cp, sz, cp, sz]
    lib.net_wl_subgraph_hashes.argtypes = [h, cpp, cpp, i64, i64, i32, i32, cp, sz, cp, sz]
    lib.net_simrank.argtypes = [h, f64, i64, f64, i32, i32, i32, f64p, cp, sz]
    lib.net_panther_similarity.argtypes = [h, i64, i64, i64, f64, f64, f64, i32, u64, i32,
                                           i64p, f64p, i64, i64p, cp, sz]
    lib.net_graph_edit_distance.argtypes = [h, h, i64p, i64p, i64p, i64p, f64p, i64, i64,
                                            f64, i64, f64p, i32p, cp, sz]
    for fn in ("net_wl_graph_hash", "net_wl_subgraph_hashes", "net_simrank",
               "net_panther_similarity", "net_graph_edit_distance"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


def hex_ints(text: str) -> tuple:
    return tuple(int(text[i:i + 8], 16) for i in range(0, len(text), 8))


def _labels(values):
    arr = (ctypes.c_char_p * max(len(values), 1))(*[v.encode() for v in values])
    return arr


def _wl_args(g, spec, fx):
    nl = _labels([str(v % 3) for v in range(fx.n)]) if spec.get("node_labels") else None
    el = _labels([str(e % 2) for e in range(len(fx.src))]) if spec.get("edge_labels") else None
    return nl, el, int(spec.get("iterations", 3)), int(spec.get("digest_size", 16))


@register("wl_hash")
def _wl_hash(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    nl, el, it, ds = _wl_args(g, spec, fx)
    out = ctypes.create_string_buffer(2 * ds)
    status = lib.net_wl_graph_hash(g.handle, nl, el, it, ds, 0, out, 2 * ds, err, ERR_LEN)
    if status and "parallel edges" in err.value.decode():
        return [-1]
    _check(status, err)
    return [hex_ints(out.raw.decode())]


@register("wl_subgraph")
def _wl_subgraph(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    nl, el, it, ds = _wl_args(g, spec, fx)
    init = int(bool(spec.get("include_initial")))
    per = (it + init) * 2 * ds
    size = fx.n * per
    out = ctypes.create_string_buffer(max(size, 1))
    status = lib.net_wl_subgraph_hashes(g.handle, nl, el, it, ds, init, 0, out, size, err,
                                        ERR_LEN)
    if status and "parallel edges" in err.value.decode():
        return [-1]
    _check(status, err)
    text = out.raw[:size].decode()
    return [hex_ints(text[v * per:(v + 1) * per]) for v in range(fx.n)]


@register("simrank")
def _simrank(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    n = fx.n
    out = (ctypes.c_double * max(n * n, 1))()
    _check(lib.net_simrank(g.handle, 0.9, 1000, float(spec.get("tolerance", 1e-4)),
                           int(bool(spec.get("weighted"))), STOP[spec.get("stop", "allclose")],
                           0, out, err, ERR_LEN), err)
    # Column-major and symmetric; read row-major.
    return [out[j * n + i] for i in range(n) for j in range(n)]


def hub(fx) -> int:
    nbrs = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            nbrs[a].add(b)
            nbrs[b].add(a)
    return max(range(fx.n), key=lambda v: (len(nbrs[v]), -v))


@register("panther")
def _panther(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    n = fx.n
    if n < 2:
        return [-1]
    s = hub(fx)
    touched = set(fx.src) | set(fx.tgt)
    k = len(touched) - 1
    nodes = (ctypes.c_int64 * (k + 1))()
    scores = (ctypes.c_double * (k + 1))()
    found = ctypes.c_int64()
    status = lib.net_panther_similarity(g.handle, s, k, 5, 0.5, 0.1, float(spec["eps"]), 0,
                                        int(spec.get("seed", 1)), 0, nodes, scores, k + 1,
                                        ctypes.byref(found), err, ERR_LEN)
    if status and "isolated source" in err.value.decode():
        return [-1]
    _check(status, err)
    out = [-1.0] * n
    for i in range(found.value):
        out[nodes[i]] = scores[i]
    return out


def ged_pair(fx):
    n = fx.n
    edges = list(zip(fx.src, fx.tgt))[1:]
    present = set(edges) | (set() if fx.directed else {(b, a) for a, b in edges})
    if (0, n - 1) not in present:
        edges.append((0, n - 1))
    edges.append((n, min(1, n - 1)))
    return n + 1, [a for a, _ in edges], [b for _, b in edges]


@register("ged")
def _ged(g, spec, fx):
    if fx.n > GED_MAX_NODES:
        return [-1]
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    n2, s2, t2 = ged_pair(fx)
    g2 = Graph(n2, s2, t2, None, fx.directed)
    labels = spec.get("labels", "none")
    nl1 = nl2 = el1 = el2 = None
    if labels in ("nodes", "both"):
        nl1 = (ctypes.c_int64 * max(fx.n, 1))(*[v % 2 for v in range(fx.n)])
        nl2 = (ctypes.c_int64 * n2)(*[v % 2 for v in range(n2)])
    if labels == "both":
        el1 = (ctypes.c_int64 * max(len(fx.src), 1))(*[e % 2 for e in range(len(fx.src))])
        el2 = (ctypes.c_int64 * max(len(s2), 1))(*[e % 2 for e in range(len(s2))])
    out = ctypes.c_double()
    found = ctypes.c_int32()
    _check(lib.net_graph_edit_distance(g.handle, g2.handle, nl1, nl2, el1, el2, None, -1, -1,
                                       float("inf"), 0, ctypes.byref(out),
                                       ctypes.byref(found), err, ERR_LEN), err)
    return [out.value]
