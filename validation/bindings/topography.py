"""Engine runners for density, components, fragmentation, centralization,
closeness, harmonic centrality, and eigenvector centrality.

Undefined results (NaN, such as the closeness of a node that reaches no
other under igraph's convention) are reported as UNDEFINED = -1. None of
these measures can be negative, and compare.py treats NaN as unequal to
everything, so the oracles for these measures write the same marker.
"""

from __future__ import annotations

import ctypes
import math

from engine import ERR_LEN, EngineError, Graph, load, register

UNDEFINED = -1.0

MODES = {"out": 1, "in": 2, "all": 3}
WEIGHTS_AS = {"distance": 0, "strength": 1}
KINDS = {"weak": 0, "strong": 1}
UNREACHABLE = {"ignore": 0, "wasserman_faust": 1, "zero": 2, "max_plus_one": 3}
MEASURES = {"degree": 0, "betweenness": 1, "closeness": 2}
CONVENTIONS = {"igraph": 0, "networkx": 1, "sna": 2, "ucinet": 3}
SCALES = {"max": 0, "unit": 1}
LOOPS = {"twice": 0, "once": 1, "ignore": 2}
KEEP = {"max": 0, "min": 1}

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    f64p = ctypes.POINTER(f64)
    h = ctypes.c_void_p
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    lib.net_graph_simplify.argtypes = [h, i32, i32, ctypes.POINTER(h), cp, sz]
    lib.net_density.argtypes = [h, i32, f64p, cp, sz]
    lib.net_components.argtypes = [h, i32, ctypes.POINTER(i32), ctypes.POINTER(i64), cp, sz]
    lib.net_fragmentation.argtypes = [h, i32, f64p, cp, sz]
    lib.net_centralization.argtypes = [h, i32, i32, i32, i32, i32, i32, f64p, cp, sz]
    lib.net_closeness.argtypes = [h, i32, i32, i32, i32, i32, i32, i32, f64p, cp, sz]
    lib.net_eigenvector.argtypes = [h, i32, i32, i32, i32, f64, i64, i32, f64p, f64p,
                                    ctypes.POINTER(i32), cp, sz]
    for fn in ("net_graph_simplify", "net_density", "net_components", "net_fragmentation",
               "net_centralization", "net_closeness", "net_eigenvector"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status: int, err) -> None:
    if status:
        raise EngineError(err.value.decode())


def _defined(values):
    return [UNDEFINED if math.isnan(v) else v for v in values]


class _Simplified:
    """A handle to the adjacency-matrix reading of a graph."""

    def __init__(self, g: Graph, undirected: bool = False, keep: str = "max"):
        self.lib = _declare()
        self.n = g.n
        self.handle = ctypes.c_void_p()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_graph_simplify(g.handle, int(undirected), KEEP[keep],
                                           ctypes.byref(self.handle), err, ERR_LEN), err)

    def __del__(self):
        if getattr(self, "handle", None):
            self.lib.net_graph_free(self.handle)
            self.handle = None


def _target(g: Graph, spec: dict):
    """The graph itself, or its adjacency-matrix reading when the spec says
    `simple: true`."""
    return _Simplified(g) if spec.get("simple") else g


def _scalar(fn, handle, *args) -> float:
    out = ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(fn(handle, *args, ctypes.byref(out), err, ERR_LEN), err)
    return out.value


@register("density")
def _density(g, spec, fx):
    lib = _declare()
    t = _target(g, spec)
    return _defined([_scalar(lib.net_density, t.handle, int(spec["loops"]))])


@register("components")
def _components(g, spec, fx):
    lib = _declare()
    n = g.n
    member = (ctypes.c_int32 * max(n, 1))()
    count = ctypes.c_int64()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_components(g.handle, KINDS[spec["kind"]], member, ctypes.byref(count),
                              err, ERR_LEN), err)
    return [float(x) for x in list(member)[:n]]


@register("fragmentation")
def _fragmentation(g, spec, fx):
    lib = _declare()
    return _defined([_scalar(lib.net_fragmentation, g.handle, KINDS[spec["kind"]])])


@register("centralization")
def _centralization(g, spec, fx):
    lib = _declare()
    value = _scalar(lib.net_centralization, g.handle, MEASURES[spec["measure"]],
                    MODES[spec.get("mode", "all")], int(spec.get("loops", True)),
                    int(spec.get("normalized", True)),
                    CONVENTIONS[spec.get("convention", "igraph")], 0)
    return _defined([value])


def _closeness(g, spec, harmonic: bool):
    lib = _declare()
    t = _target(g, spec)
    out = (ctypes.c_double * max(g.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_closeness(t.handle, int(harmonic), MODES[spec.get("mode", "out")],
                             int(spec.get("weighted", False)),
                             WEIGHTS_AS[spec.get("weights_as", "distance")],
                             int(spec.get("normalized", False)),
                             UNREACHABLE[spec.get("unreachable", "ignore")], 0, out,
                             err, ERR_LEN), err)
    return _defined(list(out)[:g.n])


@register("closeness")
def _closeness_runner(g, spec, fx):
    return _closeness(g, spec, harmonic=False)


@register("harmonic")
def _harmonic_runner(g, spec, fx):
    return _closeness(g, spec, harmonic=True)


@register("eigenvector")
def _eigenvector(g, spec, fx):
    lib = _declare()
    t = _Simplified(g, undirected=True) if spec.get("symmetrize") else _target(g, spec)
    out = (ctypes.c_double * max(g.n, 1))()
    value = ctypes.c_double()
    connected = ctypes.c_int32()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_eigenvector(t.handle, MODES[spec.get("mode", "in")],
                               int(spec.get("weighted", False)), SCALES[spec.get("scale", "max")],
                               LOOPS[spec.get("loops", "twice")], 0.0, 0, 0, out,
                               ctypes.byref(value), ctypes.byref(connected), err, ERR_LEN), err)
    return _defined(list(out)[:g.n])


def _vector(fn, handle, n, *args):
    out = (ctypes.c_double * max(n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(fn(handle, *args, out, err, ERR_LEN), err)
    return list(out)[:n]


@register("degree_simple")
def _degree_simple(g, spec, fx):
    """Degree of the adjacency-matrix reading (sna, UCINET)."""
    t = _Simplified(g, keep=spec.get("simple", "max"))
    return _vector(t.lib.net_degree, t.handle, g.n, MODES[spec["mode"]], 0, 0)


@register("betweenness_simple")
def _betweenness_simple(g, spec, fx):
    """Betweenness with parallel edges collapsed: keep "min" is NetworkX's
    reading (the shortest parallel edge), "max" the adjacency matrix."""
    t = _Simplified(g, keep=spec.get("simple", "max"))
    return _vector(t.lib.net_betweenness, t.handle, g.n, int(spec["weighted"]),
                   WEIGHTS_AS[spec.get("weights_as", "distance")], int(spec["normalized"]), 0)
