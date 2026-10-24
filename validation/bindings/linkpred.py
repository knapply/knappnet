"""Engine runners for link prediction (crates/capi/src/similarity.rs).

Community labels, which the oracles assign the same way, are node % 3.

Encodings: every result is one line per pair, (u, v, score), pairs in
lexicographic order.

Manifest engine options:
  fn: link_prediction   method (common_neighbors, jaccard, adamic_adar,
                        resource_allocation, preferential_attachment,
                        cnc, cn_sh, ra_sh, wic), alpha, delta, and
                        pairs: non_edges (every non-edge u < v) or
                        all (every pair u < v, edges included)
"""

from __future__ import annotations

import ctypes

from engine import ERR_LEN, EngineError, load, register

METHOD = {"common_neighbors": 0, "jaccard": 1, "adamic_adar": 2, "resource_allocation": 3,
          "preferential_attachment": 4, "cnc": 5, "cn_sh": 6, "ra_sh": 7, "wic": 8}

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    f64p, i64p = ctypes.POINTER(f64), ctypes.POINTER(i64)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    lib.net_link_prediction.argtypes = [h, i32, f64, f64, i64p, i64p, i64p, i64, i32, f64p,
                                        cp, sz]
    lib.net_non_edge_count.argtypes = [h, i64p, cp, sz]
    lib.net_link_prediction_non_edges.argtypes = [h, i32, f64, f64, i64p, i32, i64p, i64p,
                                                  f64p, i64, cp, sz]
    for fn in ("net_link_prediction", "net_non_edge_count", "net_link_prediction_non_edges"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


@register("link_prediction")
def _link_prediction(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    method = METHOD[spec["method"]]
    alpha = float(spec.get("alpha", 0.8))
    delta = float(spec.get("delta", 0.001))
    comm = (ctypes.c_int64 * max(g.n, 1))(*[v % 3 for v in range(g.n)])
    if spec["pairs"] == "non_edges":
        count = ctypes.c_int64()
        _check(lib.net_non_edge_count(g.handle, ctypes.byref(count), err, ERR_LEN), err)
        k = count.value
        us = (ctypes.c_int64 * max(k, 1))()
        vs = (ctypes.c_int64 * max(k, 1))()
        sc = (ctypes.c_double * max(k, 1))()
        _check(lib.net_link_prediction_non_edges(g.handle, method, alpha, delta, comm, 0, us,
                                                 vs, sc, k, err, ERR_LEN), err)
        return [(us[i], vs[i], sc[i]) for i in range(k)]
    pairs = [(u, v) for u in range(g.n) for v in range(u + 1, g.n)]
    k = len(pairs)
    us = (ctypes.c_int64 * max(k, 1))(*[p[0] for p in pairs])
    vs = (ctypes.c_int64 * max(k, 1))(*[p[1] for p in pairs])
    sc = (ctypes.c_double * max(k, 1))()
    _check(lib.net_link_prediction(g.handle, method, alpha, delta, comm, us, vs, k, 0, sc, err,
                                   ERR_LEN), err)
    return [(u, v, sc[i]) for i, (u, v) in enumerate(pairs)]
