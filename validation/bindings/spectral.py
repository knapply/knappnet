"""Engine runners for PageRank, Katz centrality, Bonacich power centrality,
hubs and authorities, assortativity, mixing matrices, the E-I index, and
within-group tie shares.

Golden-file conventions shared with the oracles in
validation/oracles/spectral_*_oracle.py and validation/oracles/r_spectral.R:

- An undefined value (NaN) is written as an empty line and read back as an
  empty group, so undefined matches undefined and nothing else. (Bonacich
  scores and assortativity can be negative, so no numeric marker is safe.)
- When the engine refuses a fixture (for example, Bonacich power
  centrality with an exponent whose I - exponent A is singular), the runner
  returns no values and the oracles write no golden file; the manifest lists
  such fixtures under `undefined_on`.
- Fixtures have no node attributes, so the manifest names deterministic
  rules (see `node_rule`): `mod3` is node index modulo 3, `id` gives every
  node its own category, `mod3plus1` is (v mod 3) + 1, `num7` is
  (7 v mod 11), and `num3` is (3 v mod 5).
- Some measures take a parameter relative to the spectral radius of the
  adjacency matrix the measure uses (`alpha_fraction`,
  `exponent_fraction`): the value is that fraction of 1 / lambda_max, with
  lambda_max computed by each oracle independently (here, by the engine's
  eigenvector centrality).
- Mixing matrices are written row by row, one row per line. Within-group
  shares are written one group per line: size, internal, external,
  observed, expected.
"""

from __future__ import annotations

import ctypes
import math

from engine import ERR_LEN, EngineError, load, register

MODES = {"out": 1, "in": 2, "all": 3}
LOOPS = {"twice": 0, "once": 1, "ignore": 2}
DANGLING = {"personalized": 0, "uniform": 1}
HITS_SCALE = {"max": 0, "sum": 1, "unit": 2}
SCOPES = {"global": 0, "group": 1, "node": 2}
COUNTING = {"standard": 0, "homophily": 1}

_declared = False


def _lib():
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    f64p = ctypes.POINTER(f64)
    i64p = ctypes.POINTER(i64)
    i32p = ctypes.POINTER(i32)
    h = ctypes.c_void_p
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    lib.net_pagerank.argtypes = [h, f64, i32, i32, f64p, i32, i32, f64, i64, i32, f64p, cp, sz]
    lib.net_katz.argtypes = [h, f64, f64p, i64, i32, i32, i32, i32, f64, i64, i32, f64p, cp, sz]
    lib.net_bonacich_power.argtypes = [h, f64, i32, i32, i32, i32, f64, i64, i32, f64p, cp, sz]
    lib.net_hits.argtypes = [h, i32, i32, i32, f64, i64, i32, f64p, f64p, f64p, i32p, cp, sz]
    lib.net_assortativity_nominal.argtypes = [h, i64p, i64, i32, i32, i32, f64p, cp, sz]
    lib.net_assortativity_values.argtypes = [h, f64p, f64p, i32, i32, i32, f64p, cp, sz]
    lib.net_assortativity_degree.argtypes = [h, i32, i32, i32, i32, f64p, cp, sz]
    lib.net_mixing_matrix.argtypes = [h, i64p, i64, i64p, i64, i32, i32, f64p, i32p, cp, sz]
    lib.net_ei_index.argtypes = [h, i64p, i64, i32, i32, i32, f64p, f64p, f64p, cp, sz]
    lib.net_group_shares.argtypes = [h, i64p, i64, i32, f64p, cp, sz]
    lib.net_graph_simplify.argtypes = [h, i32, i32, ctypes.POINTER(h), cp, sz]
    lib.net_eigenvector.argtypes = [h, i32, i32, i32, i32, f64, i64, i32, f64p, f64p, i32p,
                                    cp, sz]
    for fn in ("net_pagerank", "net_katz", "net_bonacich_power", "net_hits",
               "net_assortativity_nominal", "net_assortativity_values",
               "net_assortativity_degree", "net_mixing_matrix", "net_ei_index",
               "net_group_shares", "net_graph_simplify", "net_eigenvector"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def node_rule(rule: str, n: int) -> list:
    """Deterministic node attributes for fixtures without any."""
    if rule == "mod3":
        return [v % 3 for v in range(n)]
    if rule == "id":
        return list(range(n))
    if rule == "mod3plus1":
        return [float(v % 3 + 1) for v in range(n)]
    if rule == "num7":
        return [float(7 * v % 11) for v in range(n)]
    if rule == "num3":
        return [float(3 * v % 5) for v in range(n)]
    raise ValueError(f"unknown node rule {rule!r}")


def codes(values: list) -> tuple[list[int], int]:
    """Category codes in order of first appearance, and their number."""
    seen: dict = {}
    out = [seen.setdefault(v, len(seen)) for v in values]
    return out, len(seen)


def undefined_as_empty(values):
    return [() if isinstance(v, float) and math.isnan(v) else v for v in values]


def _check(status: int, err) -> None:
    if status:
        raise EngineError(err.value.decode())


def _f64(values):
    return (ctypes.c_double * max(len(values), 1))(*values)


def _i64(values):
    return (ctypes.c_int64 * max(len(values), 1))(*values)


class _Simplified:
    """The adjacency-matrix reading of a graph (parallel edges collapsed,
    self-loops dropped), as sna reads it."""

    def __init__(self, g):
        self.lib = _lib()
        self.n = g.n
        self.handle = ctypes.c_void_p()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_graph_simplify(g.handle, 0, 0, ctypes.byref(self.handle), err,
                                           ERR_LEN), err)

    def __del__(self):
        if getattr(self, "handle", None):
            self.lib.net_graph_free(self.handle)
            self.handle = None


def _target(g, spec):
    return _Simplified(g) if spec.get("simple") else g


def spectral_radius(t, mode: str, weighted: bool, loops: str) -> float:
    """Largest eigenvalue of the adjacency matrix a measure uses, from the
    engine's eigenvector centrality."""
    lib = _lib()
    out = (ctypes.c_double * max(t.n, 1))()
    value = ctypes.c_double()
    connected = ctypes.c_int32()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_eigenvector(t.handle, MODES[mode], int(weighted), 0, LOOPS[loops], 0.0, 0,
                               0, out, ctypes.byref(value), ctypes.byref(connected), err,
                               ERR_LEN), err)
    return value.value


def relative(spec: dict, key: str, t, mode: str, weighted: bool, loops: str):
    """`key` as given, or `key_fraction` / lambda_max."""
    if key in spec:
        return float(spec[key])
    lam = spectral_radius(t, mode, weighted, loops)
    if lam <= 0.0:
        return float(spec[key + "_fraction"])
    return float(spec[key + "_fraction"]) / lam


def _vector(fn, t, n, *args):
    out = (ctypes.c_double * max(n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    status = fn(t.handle, *args, out, err, ERR_LEN)
    if status:
        return None
    return undefined_as_empty(list(out)[:n])


def _scalar(fn, t, *args):
    out = ctypes.c_double()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(fn(t.handle, *args, ctypes.byref(out), err, ERR_LEN), err)
    return undefined_as_empty([out.value])


@register("pagerank")
def _pagerank(g, spec, fx):
    lib = _lib()
    p = spec.get("personalization")
    pers = _f64(node_rule(p, g.n)) if p else None
    values = _vector(lib.net_pagerank, g, g.n, float(spec.get("damping", 0.85)),
                     int(spec.get("directed", True)), int(spec.get("weighted", False)), pers,
                     DANGLING[spec.get("dangling", "personalized")],
                     LOOPS[spec.get("loops", "twice")], 0.0, 0, 0)
    return values if values is not None else []


@register("katz")
def _katz(g, spec, fx):
    lib = _lib()
    mode, weighted, loops = spec.get("mode", "in"), spec.get("weighted", False), \
        spec.get("loops", "once")
    alpha = relative(spec, "alpha", g, mode, weighted, loops)
    b = spec.get("beta", 1.0)
    beta = node_rule(b, g.n) if isinstance(b, str) else [float(b)]
    values = _vector(lib.net_katz, g, g.n, alpha, _f64(beta), len(beta), MODES[mode],
                     int(weighted), LOOPS[loops], int(spec.get("normalized", True)), 0.0, 0, 0)
    return values if values is not None else []


@register("bonacich")
def _bonacich(g, spec, fx):
    lib = _lib()
    t = _target(g, spec)
    mode, weighted, loops = spec.get("mode", "out"), spec.get("weighted", False), \
        spec.get("loops", "ignore")
    exponent = relative(spec, "exponent", t, "in" if mode == "out" else mode, weighted, loops)
    values = _vector(lib.net_bonacich_power, t, g.n, exponent, MODES[mode], int(weighted),
                     LOOPS[loops], int(spec.get("rescale", False)), 0.0, 0, 0)
    return values if values is not None else []


@register("hits")
def _hits(g, spec, fx):
    lib = _lib()
    n = g.n
    hubs = (ctypes.c_double * max(n, 1))()
    auth = (ctypes.c_double * max(n, 1))()
    value = ctypes.c_double()
    unique = ctypes.c_int32()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_hits(g.handle, int(spec.get("weighted", False)),
                        LOOPS[spec.get("loops", "twice")], HITS_SCALE[spec.get("scale", "max")],
                        0.0, 0, 0, hubs, auth, ctypes.byref(value), ctypes.byref(unique), err,
                        ERR_LEN), err)
    out = hubs if spec["which"] == "hubs" else auth
    return undefined_as_empty(list(out)[:n])


def _assort(spec):
    return (int(spec.get("directed", True)), LOOPS[spec.get("loops", "twice")],
            int(spec.get("normalized", True)))


@register("assortativity_nominal")
def _assortativity_nominal(g, spec, fx):
    lib = _lib()
    c, k = codes(node_rule(spec["types"], g.n))
    return _scalar(lib.net_assortativity_nominal, g, _i64(c), k, *_assort(spec))


@register("assortativity_values")
def _assortativity_values(g, spec, fx):
    lib = _lib()
    vals = _f64(node_rule(spec["values"], g.n))
    vin = _f64(node_rule(spec["values_in"], g.n)) if spec.get("values_in") else None
    return _scalar(lib.net_assortativity_values, g, vals, vin, *_assort(spec))


@register("assortativity_degree")
def _assortativity_degree(g, spec, fx):
    lib = _lib()
    return _scalar(lib.net_assortativity_degree, g, int(spec.get("weighted", False)),
                   *_assort(spec))


@register("mixing_matrix")
def _mixing_matrix(g, spec, fx):
    lib = _lib()
    d1, k1 = codes(node_rule(spec["dim1"], g.n))
    same = "dim2" not in spec
    d2, k2 = (d1, k1) if same else codes(node_rule(spec["dim2"], g.n))
    out = (ctypes.c_double * max(k1 * k2, 1))()
    transposed = ctypes.c_int32()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_mixing_matrix(g.handle, _i64(d1), k1, _i64(d2), k2, int(same),
                                 MODES[spec.get("direction", "all")], out,
                                 ctypes.byref(transposed), err, ERR_LEN), err)
    rows, cols = (k2, k1) if transposed.value else (k1, k2)
    flat = list(out)[:k1 * k2]
    return [tuple(flat[r * cols:(r + 1) * cols]) for r in range(rows)]


@register("ei_index")
def _ei_index(g, spec, fx):
    lib = _lib()
    c, k = codes(node_rule(spec["groups"], g.n))
    scope = spec.get("scope", "global")
    size = {"global": 1, "group": k, "node": g.n}[scope]
    out = (ctypes.c_double * max(size, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_ei_index(g.handle, _i64(c), k, SCOPES[scope],
                            MODES[spec.get("direction", "all")],
                            COUNTING[spec.get("counting", "standard")], out, None, None, err,
                            ERR_LEN), err)
    return undefined_as_empty(list(out)[:size])


@register("group_shares")
def _group_shares(g, spec, fx):
    lib = _lib()
    c, k = codes(node_rule(spec["groups"], g.n))
    out = (ctypes.c_double * max(5 * k, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_group_shares(g.handle, _i64(c), k, MODES[spec.get("direction", "all")],
                                out, err, ERR_LEN), err)
    flat = list(out)[:5 * k]
    rows = []
    for r in range(k):
        row = flat[5 * r:5 * r + 5]
        rows.append(tuple(-1.0 if math.isnan(x) else x for x in row))
    return rows
