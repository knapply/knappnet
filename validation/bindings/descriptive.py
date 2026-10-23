"""Engine runners for reciprocity, transitivity and clustering, distances,
mean distance, diameter, eccentricity, radius, reachability, cut points,
and bridges.

Non-finite results are written with markers, since compare.py treats NaN
as unequal to everything: UNDEFINED = -999 for NaN and INFINITE = -1000
for infinity. No measure here can take either value (the reciprocity
correlation lies in [-1, 1]), and every oracle for these measures writes
the same markers.
"""

from __future__ import annotations

import ctypes
import math

from engine import ERR_LEN, EngineError, load, register

UNDEFINED = -999.0
INFINITE = -1000.0

MODES = {"out": 1, "in": 2, "all": 3}
WEIGHTS_AS = {"distance": 0, "strength": 1}
RECIPROCITY = {"default": 0, "ratio": 1, "dyadic": 2, "correlation": 3}
WEIGHTING = {"none": 0, "barrat": 1, "onnela": 2}
UNDEFINED_AS = {"nan": 0, "zero": 1}

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    f64p, i64p, i32p = ctypes.POINTER(f64), ctypes.POINTER(i64), ctypes.POINTER(i32)
    h = ctypes.c_void_p
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    lib.net_reciprocity.argtypes = [h, i32, f64p, cp, sz]
    lib.net_transitivity.argtypes = [h, i32, i32, f64p, cp, sz]
    lib.net_clustering.argtypes = [h, i32, i32, i32, i32, f64p, cp, sz]
    lib.net_distances.argtypes = [h, i64p, i64, i32, i32, i32, i32, f64p, cp, sz]
    lib.net_path_summary.argtypes = [h, i64p, i64, i32, i32, i32, i32, i32, f64p, cp, sz]
    lib.net_eccentricity.argtypes = [h, i32, i32, i32, i32, f64p, cp, sz]
    lib.net_reachability.argtypes = [h, i32, i32, f64p, cp, sz]
    lib.net_cut_structure.argtypes = [h, i32p, i32p, cp, sz]
    lib.net_edge_counts.argtypes = [h, f64p, cp, sz]
    for fn in ("net_reciprocity", "net_transitivity", "net_clustering", "net_distances",
               "net_path_summary", "net_eccentricity", "net_reachability",
               "net_cut_structure", "net_edge_counts"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status: int, err) -> None:
    if status:
        raise EngineError(err.value.decode())


def _mark(x: float) -> float:
    if math.isnan(x):
        return UNDEFINED
    if math.isinf(x):
        return INFINITE
    return x


def _doubles(fn, handle, size, *args):
    out = (ctypes.c_double * max(size, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(fn(handle, *args, out, err, ERR_LEN), err)
    return [_mark(x) for x in list(out)[:size]]


def _path_args(spec):
    return (MODES[spec.get("mode", "out")], int(spec.get("weighted", False)),
            WEIGHTS_AS[spec.get("weights_as", "distance")])


@register("reciprocity")
def _reciprocity(g, spec, fx):
    lib = _declare()
    return _doubles(lib.net_reciprocity, g.handle, 1, RECIPROCITY[spec["mode"]])


@register("transitivity")
def _transitivity(g, spec, fx):
    lib = _declare()
    return _doubles(lib.net_transitivity, g.handle, 1,
                    UNDEFINED_AS[spec.get("undefined", "nan")], 0)


def _clustering(g, spec, average: bool):
    lib = _declare()
    return _doubles(lib.net_clustering, g.handle, 1 if average else g.n, int(average),
                    WEIGHTING[spec.get("weighting", "none")],
                    UNDEFINED_AS[spec.get("undefined", "nan")], 0)


@register("clustering")
def _clustering_runner(g, spec, fx):
    return _clustering(g, spec, average=False)


@register("average_clustering")
def _average_clustering_runner(g, spec, fx):
    return _clustering(g, spec, average=True)


@register("distances")
def _distances(g, spec, fx):
    """One tuple per source (row), in node order."""
    lib = _declare()
    n = g.n
    flat = _doubles(lib.net_distances, g.handle, n * n, None, -1, *_path_args(spec), 0)
    # Column-major: entry (i, j) is at j * n + i.
    return [tuple(flat[j * n + i] for j in range(n)) for i in range(n)]


def _summary(g, spec):
    lib = _declare()
    return _doubles(lib.net_path_summary, g.handle, 7, None, -1, *_path_args(spec),
                    int(spec.get("unconnected", True)), 0)


@register("mean_distance")
def _mean_distance(g, spec, fx):
    return [_summary(g, spec)[5]]


@register("diameter")
def _diameter(g, spec, fx):
    return [_summary(g, spec)[6]]


@register("eccentricity")
def _eccentricity(g, spec, fx):
    lib = _declare()
    return _doubles(lib.net_eccentricity, g.handle, g.n, *_path_args(spec), 0)


@register("radius")
def _radius(g, spec, fx):
    ecc = _eccentricity(g, spec, fx)
    return [min(ecc) if ecc else UNDEFINED]


@register("reachability")
def _reachability(g, spec, fx):
    lib = _declare()
    return _doubles(lib.net_reachability, g.handle, g.n, MODES[spec["mode"]], 0)


def _cut_structure(g, fx):
    lib = _declare()
    m = len(fx.src)
    cut = (ctypes.c_int32 * max(g.n, 1))()
    bridge = (ctypes.c_int32 * max(m, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_cut_structure(g.handle, cut, bridge, err, ERR_LEN), err)
    return list(cut)[:g.n], list(bridge)[:m]


@register("cut_points")
def _cut_points(g, spec, fx):
    """Indices of the cut points."""
    cut, _ = _cut_structure(g, fx)
    return [float(v) for v in range(g.n) if cut[v]]


@register("bridges")
def _bridges(g, spec, fx):
    """Each bridge as its pair of endpoints, smaller first. A bridge is
    never a parallel edge, so the pair identifies it."""
    _, bridge = _cut_structure(g, fx)
    return [tuple(sorted((float(fx.src[e]), float(fx.tgt[e]))))
            for e in range(len(fx.src)) if bridge[e]]
