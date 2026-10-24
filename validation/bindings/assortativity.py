"""Engine runners for NetworkX-style mixing: mixing counts, assortativity
over node subsets, average neighbour degree, and average degree
connectivity (crates/capi/src/assortativity.rs).

Derived inputs, which the oracles build the same way:
  labels    node v carries label v % 3 (attributes) and value
            (v * 7) % 5 (numeric attributes)
  subset    the even-numbered nodes (`subset: true`)

Encodings: mixing counts are one line (x, y, count) per non-zero entry in
increasing (x, y) order; a coefficient is one value, NaN written as -999;
average neighbour degree is node-aligned (subset nodes only, in order);
average degree connectivity is one line (k, k_nn) per degree.

Manifest engine options:
  fn: mixing        kind (degree | attribute), x, y, weighted, subset
  fn: assort        coefficient (degree | attribute | numeric), x, y,
                    weighted, subset
  fn: avg_nbr_deg   source, target, weighted, subset
  fn: avg_deg_conn  source, target, weighted, subset
"""

from __future__ import annotations

import ctypes
import math

from engine import ERR_LEN, EngineError, load, register

KIND = {"out": 1, "in": 2, "in+out": 3}
COEF = {"degree": 0, "attribute": 1, "numeric": 2}
UNDEFINED = -999.0

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    f64p, i64p = ctypes.POINTER(f64), ctypes.POINTER(i64)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    lib.net_degree_mixing.argtypes = [h, i32, i32, i32, i64p, i64, i64, f64p, f64p, i64p, cp, sz]
    lib.net_attribute_mixing.argtypes = [h, f64p, i64p, i64, i64, f64p, f64p, i64p, cp, sz]
    lib.net_assortativity_subset.argtypes = [h, i32, i32, i32, i32, f64p, i64p, i64, f64p, cp,
                                             sz]
    lib.net_average_neighbor_degree.argtypes = [h, i32, i32, i32, i64p, i64, f64p, cp, sz]
    lib.net_average_degree_connectivity.argtypes = [h, i32, i32, i32, i64p, i64, f64p, f64p,
                                                    i64p, cp, sz]
    for fn in ("net_degree_mixing", "net_attribute_mixing", "net_assortativity_subset",
               "net_average_neighbor_degree", "net_average_degree_connectivity"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


def subset(fx, spec):
    return [v for v in range(fx.n) if v % 2 == 0] if spec.get("subset") else None


def _nodes(fx, spec):
    s = subset(fx, spec)
    if s is None:
        return None, 0
    return (ctypes.c_int64 * max(len(s), 1))(*s), len(s)


def labels(fx, kind):
    if kind == "numeric":
        return [float((v * 7) % 5) for v in range(fx.n)]
    return [float(v % 3) for v in range(fx.n)]


@register("mixing")
def _mixing(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    nodes, nn = _nodes(fx, spec)
    k = ctypes.c_int64()
    cap = 1
    while True:
        vals = (ctypes.c_double * cap)()
        counts = (ctypes.c_double * (cap * cap))()
        if spec["kind"] == "degree":
            _check(lib.net_degree_mixing(g.handle, KIND[spec.get("x", "out")],
                                         KIND[spec.get("y", "in")],
                                         int(bool(spec.get("weighted"))), nodes, nn, cap, vals,
                                         counts, ctypes.byref(k), err, ERR_LEN), err)
        else:
            lab = (ctypes.c_double * max(fx.n, 1))(*labels(fx, "attribute"))
            _check(lib.net_attribute_mixing(g.handle, lab, nodes, nn, cap, vals, counts,
                                            ctypes.byref(k), err, ERR_LEN), err)
        if k.value <= cap:
            break
        cap = k.value  # too small: ask again with the size reported
    kk = k.value
    return [(vals[i], vals[j], counts[i * kk + j]) for i in range(kk) for j in range(kk)
            if counts[i * kk + j] != 0]


@register("assort")
def _assort(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    nodes, nn = _nodes(fx, spec)
    coef = spec["coefficient"]
    vals = None
    if coef != "degree":
        vals = (ctypes.c_double * max(fx.n, 1))(*labels(fx, coef))
    out = ctypes.c_double()
    _check(lib.net_assortativity_subset(g.handle, COEF[coef], KIND[spec.get("x", "out")],
                                        KIND[spec.get("y", "in")], int(bool(spec.get("weighted"))),
                                        vals, nodes, nn, ctypes.byref(out), err, ERR_LEN), err)
    return [UNDEFINED if math.isnan(out.value) else out.value]


@register("avg_nbr_deg")
def _avg_nbr_deg(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    nodes, nn = _nodes(fx, spec)
    size = nn if nodes is not None else fx.n
    out = (ctypes.c_double * max(size, 1))()
    _check(lib.net_average_neighbor_degree(g.handle, KIND[spec.get("source", "out")],
                                           KIND[spec.get("target", "out")],
                                           int(bool(spec.get("weighted"))), nodes, nn, out, err,
                                           ERR_LEN), err)
    return [out[i] for i in range(size)]


@register("avg_deg_conn")
def _avg_deg_conn(g, spec, fx):
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    nodes, nn = _nodes(fx, spec)
    cap = max(fx.n, 1)
    ks = (ctypes.c_double * cap)()
    vs = (ctypes.c_double * cap)()
    c = ctypes.c_int64()
    _check(lib.net_average_degree_connectivity(g.handle, KIND[spec.get("source", "in+out")],
                                               KIND[spec.get("target", "in+out")],
                                               int(bool(spec.get("weighted"))), nodes, nn, ks, vs,
                                               ctypes.byref(c), err, ERR_LEN), err)
    return [(ks[i], vs[i]) for i in range(c.value)]
