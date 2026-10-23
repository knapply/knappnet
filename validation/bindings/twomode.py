"""Engine runners for ego networks and two-mode networks.

Undefined results (NaN, such as the density of an ego network with fewer
than two alters) are reported as UNDEFINED = -1. None of these measures can
be negative, and compare.py treats NaN as unequal to everything, so the
oracles write the same marker.

Two-mode measures need the node types, which two-mode fixtures record as
`type` in meta.json (0 = first mode, 1 = second mode). On fixtures without
types the runners return no values, and the manifest declares those
fixtures undefined for every oracle.

A projection is returned as one (from, to, weight) triple per edge, from
<= to, ordered by from and then to, with 0-based node indices of the
original network.
"""

from __future__ import annotations

import ctypes
import json
import math

from common import FIXTURES
from engine import ERR_LEN, EngineError, load, register

UNDEFINED = -1.0

MODES = {"out": 1, "in": 2, "all": 3}
MEASURES = {"degree": 0, "closeness": 1, "betweenness": 2}
UNREACHABLE = {"ignore": 0, "wasserman_faust": 1}
WEIGHTS = {"count": 0, "newman": 1, "product": 2}
SUMMARY = ["size", "ties", "pairs", "density", "components", "betweenness"]

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    i32p, i64p, f64p = (ctypes.POINTER(i32), ctypes.POINTER(i64), ctypes.POINTER(f64))
    h = ctypes.c_void_p
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    lib.net_ego_sizes.argtypes = [h, i64p, i64, i64, i32, i64, i32, f64p, cp, sz]
    lib.net_ego_members.argtypes = [h, i64p, i64, i64, i32, i64, i32, i64p, i64p, i64,
                                    cp, sz]
    lib.net_ego_summary.argtypes = [h, i64, i32, i32, f64p, f64p, f64p, f64p, f64p, f64p,
                                    cp, sz]
    lib.net_twomode_density.argtypes = [h, i32p, f64p, cp, sz]
    lib.net_twomode_centrality.argtypes = [h, i32p, i32, i32, i32, i32, f64p, cp, sz]
    lib.net_twomode_eigenvector.argtypes = [h, i32p, f64, i64, i32, f64p, f64p, i32p,
                                            cp, sz]
    lib.net_project.argtypes = [h, i32p, i32, i32, i32, i32, ctypes.POINTER(h), cp, sz]
    lib.net_projection_count.argtypes = [h]
    lib.net_projection_count.restype = i64
    lib.net_projection_copy.argtypes = [h, i64p, i64p, f64p, cp, sz]
    lib.net_projection_free.argtypes = [h]
    for fn in ("net_ego_sizes", "net_ego_members", "net_ego_summary", "net_twomode_density",
               "net_twomode_centrality", "net_twomode_eigenvector", "net_project",
               "net_projection_copy"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status: int, err) -> None:
    if status:
        raise EngineError(err.value.decode())


def _defined(values):
    return [UNDEFINED if math.isnan(v) else v for v in values]


def node_types(fx):
    """The fixture's node types, or None for a one-mode fixture."""
    meta = json.loads((FIXTURES / fx.name / "meta.json").read_text())
    return meta.get("type")


# ---- Ego networks -----------------------------------------------------------

@register("ego_size")
def _ego_size(g, spec, fx):
    lib = _declare()
    out = (ctypes.c_double * max(g.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_ego_sizes(g.handle, None, 0, spec.get("order", 1),
                             MODES[spec.get("mode", "all")], spec.get("mindist", 0), 0, out,
                             err, ERR_LEN), err)
    return list(out)[:g.n]


@register("ego_members")
def _ego_members(g, spec, fx):
    """Every node's neighborhood, as one group per node."""
    lib = _declare()
    n = g.n
    order, mode, mindist = spec.get("order", 1), MODES[spec.get("mode", "all")], spec.get(
        "mindist", 0)
    egos = (ctypes.c_int64 * max(n, 1))(*range(n))
    sizes = (ctypes.c_double * max(n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_ego_sizes(g.handle, egos, n, order, mode, mindist, 0, sizes, err,
                             ERR_LEN), err)
    total = int(sum(list(sizes)[:n]))
    offsets = (ctypes.c_int64 * (n + 1))()
    members = (ctypes.c_int64 * max(total, 1))()
    _check(lib.net_ego_members(g.handle, egos, n, order, mode, mindist, 0, offsets, members,
                               total, err, ERR_LEN), err)
    return [tuple(float(members[k]) for k in range(offsets[i], offsets[i + 1]))
            for i in range(n)]


@register("ego_summary")
def _ego_summary(g, spec, fx):
    lib = _declare()
    n = g.n
    bufs = [(ctypes.c_double * max(n, 1))() for _ in SUMMARY]
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_ego_summary(g.handle, spec.get("order", 1), MODES[spec.get("mode", "all")],
                               0, *bufs, err, ERR_LEN), err)
    return _defined(list(bufs[SUMMARY.index(spec["measure"])])[:n])


# ---- Two-mode networks ------------------------------------------------------

def _types_array(types):
    return (ctypes.c_int32 * max(len(types), 1))(*types)


@register("twomode")
def _twomode(g, spec, fx):
    types = node_types(fx)
    if types is None:
        return []
    lib = _declare()
    t = _types_array(types)
    err = ctypes.create_string_buffer(ERR_LEN)
    measure = spec["measure"]
    if measure == "density":
        out = ctypes.c_double()
        _check(lib.net_twomode_density(g.handle, t, ctypes.byref(out), err, ERR_LEN), err)
        return _defined([out.value])
    out = (ctypes.c_double * max(g.n, 1))()
    if measure == "eigenvector":
        value = ctypes.c_double()
        connected = ctypes.c_int32()
        _check(lib.net_twomode_eigenvector(g.handle, t, 0.0, 0, 0, out, ctypes.byref(value),
                                           ctypes.byref(connected), err, ERR_LEN), err)
    else:
        _check(lib.net_twomode_centrality(g.handle, t, MEASURES[measure],
                                          int(spec.get("normalized", True)),
                                          UNREACHABLE[spec.get("unreachable", "ignore")], 0,
                                          out, err, ERR_LEN), err)
    return _defined(list(out)[:g.n])


@register("project")
def _project(g, spec, fx):
    types = node_types(fx)
    if types is None:
        return []
    lib = _declare()
    err = ctypes.create_string_buffer(ERR_LEN)
    p = ctypes.c_void_p()
    _check(lib.net_project(g.handle, _types_array(types), int(spec["second"]),
                           WEIGHTS[spec.get("weight", "count")], int(spec.get("loops", False)),
                           0, ctypes.byref(p), err, ERR_LEN), err)
    try:
        m = lib.net_projection_count(p)
        a = (ctypes.c_int64 * max(m, 1))()
        b = (ctypes.c_int64 * max(m, 1))()
        w = (ctypes.c_double * max(m, 1))()
        _check(lib.net_projection_copy(p, a, b, w, err, ERR_LEN), err)
        return [(float(a[k]), float(b[k]), w[k]) for k in range(m)]
    finally:
        lib.net_projection_free(p)
