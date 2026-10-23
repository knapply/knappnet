"""Engine runners for reach, edge betweenness, information centrality,
bridging, Burt's hierarchy, eigenvector centralization, regular and
automorphic equivalence, blockmodels, classical scaling, and geospatial
networks (validation/manifest.d/roles.yaml).

Golden-file conventions shared with validation/oracles/roles_*_oracle.py
and validation/oracles/r_roles.R:

- An undefined value (NaN) is written as an empty line and read back as an
  empty group, so undefined matches undefined and nothing else.
- Fixtures carry no partitions or coordinates, so the inputs these
  measures need are derived deterministically from the fixture by the
  rules in this file (`partition_for`, `coordinates`, `mds_input`), which
  every oracle reimplements. They describe data, not the measures.
- Block densities are written row by row, one block per line (k values).
- Classical scaling writes one line per point with its coordinates, and
  edge lists one line per edge: from, to, distance.
"""

from __future__ import annotations

import ctypes
import math
from collections import deque

from engine import ERR_LEN, EngineError, load, register

MODES = {"out": 1, "in": 2, "all": 3}
WEIGHTS_AS = {"distance": 0, "strength": 1}
CRS = {"lonlat": 0, "planar": 1}
EARTH_RADIUS_M = 6371008.7714

_declared = False


def _lib():
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    f64p, i64p, i32p = (ctypes.POINTER(f64), ctypes.POINTER(i64), ctypes.POINTER(i32))
    h = ctypes.c_void_p
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    lib.net_reach.argtypes = [h, f64, i32, i32, i32, i32, i32, f64p, cp, sz]
    lib.net_edge_betweenness.argtypes = [h, i32, i32, i32, i32, f64p, cp, sz]
    lib.net_information_centrality.argtypes = [h, i32, f64p, cp, sz]
    lib.net_bridging_coefficient.argtypes = [h, f64p, cp, sz]
    lib.net_betweenness.argtypes = [h, i32, i32, i32, i32, f64p, cp, sz]
    lib.net_structural_holes.argtypes = [h, i32, i32, i32, f64p, cp, sz]
    lib.net_centralization.argtypes = [h, i32, i32, i32, i32, i32, i32, f64p, cp, sz]
    lib.net_regular_equivalence.argtypes = [h, i64p, i64, i32, i32p, i32p, cp, sz]
    lib.net_automorphic_equivalence.argtypes = [h, i64, i32p, cp, sz]
    lib.net_block_model.argtypes = [h, i64p, i64, i64, i32, i32, f64p, f64p, f64p, f64p,
                                    cp, sz]
    lib.net_classical_mds.argtypes = [f64p, i64, i64, f64p, i32p, f64p, f64p, cp, sz]
    lib.net_mds_stress.argtypes = [f64p, i64, f64p, i64, f64p, cp, sz]
    lib.net_edge_lengths.argtypes = [h, f64p, f64p, i64, i32, f64, f64p, cp, sz]
    lib.net_proximity.argtypes = [f64p, f64p, i64, i32, f64, f64, i64, i32,
                                  ctypes.POINTER(h), cp, sz]
    lib.net_edgelist_size.argtypes = [h]
    lib.net_edgelist_size.restype = i64
    lib.net_edgelist_copy.argtypes = [h, i64p, i64p, f64p, cp, sz]
    lib.net_edgelist_free.argtypes = [h]
    lib.net_edgelist_free.restype = None
    lib.net_group_geometry.argtypes = [f64p, f64p, i64p, i64, i32, f64, i64, f64p, f64p,
                                       f64p, f64p, f64p, cp, sz]
    for fn in ("net_reach", "net_edge_betweenness", "net_information_centrality",
               "net_bridging_coefficient", "net_betweenness", "net_structural_holes",
               "net_centralization", "net_regular_equivalence",
               "net_automorphic_equivalence", "net_block_model", "net_classical_mds",
               "net_mds_stress", "net_edge_lengths", "net_proximity", "net_edgelist_copy",
               "net_group_geometry"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status: int, err) -> None:
    if status:
        raise EngineError(err.value.decode())


def _undefined_as_empty(values):
    return [() if isinstance(v, float) and math.isnan(v) else v for v in values]


def _doubles(values):
    return (ctypes.c_double * max(len(values), 1))(*values)


def _vector(fn, handle, size, *args):
    out = (ctypes.c_double * max(size, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(fn(handle, *args, out, err, ERR_LEN), err)
    return list(out)[:size]


# ---- Inputs derived from a fixture ------------------------------------------

def partition_for(rule: str, n: int) -> list[int]:
    """Deterministic partitions: `mod3` puts node v in block v % 3, `mod2`
    in block v % 2."""
    if rule == "mod3":
        return [v % 3 for v in range(n)]
    if rule == "mod2":
        return [v % 2 for v in range(n)]
    raise ValueError(f"unknown partition rule {rule!r}")


def _hash_unit(v: int, salt: int) -> float:
    """A number in [0, 1) from node index v, by a multiplicative hash in
    exact integer arithmetic (identical on every platform)."""
    h = ((v + 1) * 2654435761 + salt * 40503 + 12345) % 4294967296
    h = (h ^ (h >> 13)) * 1274126177 % 4294967296
    return h / 4294967296.0


def coordinates(crs: str, n: int) -> tuple[list[float], list[float]]:
    """Node coordinates: longitude in [-180, 180) and latitude uniform on
    the sphere (asin of a uniform value), or planar x, y in [0, 100)."""
    u = [_hash_unit(v, 1) for v in range(n)]
    w = [_hash_unit(v, 2) for v in range(n)]
    if crs == "lonlat":
        return ([360.0 * a - 180.0 for a in u],
                [math.degrees(math.asin(2.0 * b - 1.0)) for b in w])
    return ([100.0 * a for a in u], [100.0 * b for b in w])


def mds_input(fx) -> list[list[float]]:
    """Hop distances on the undirected, simple reading of the fixture, with
    unreachable pairs set to the largest finite distance plus 1 (1 when no
    pair is connected)."""
    n = fx.n
    nbrs = [set() for _ in range(n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            nbrs[a].add(b)
            nbrs[b].add(a)
    d = [[math.inf] * n for _ in range(n)]
    for s in range(n):
        d[s][s] = 0.0
        q = deque([s])
        while q:
            v = q.popleft()
            for w in nbrs[v]:
                if d[s][w] == math.inf:
                    d[s][w] = d[s][v] + 1.0
                    q.append(w)
    finite = [x for row in d for x in row if x != math.inf]
    fill = max(finite) + 1.0 if finite else 1.0
    return [[fill if x == math.inf else x for x in row] for row in d]


# ---- Centrality and brokerage ------------------------------------------------

@register("reach")
def _reach(g, spec, fx):
    lib = _lib()
    k = math.inf if spec["k"] == "inf" else float(spec["k"])
    values = _vector(lib.net_reach, g.handle, g.n, k, MODES[spec.get("mode", "out")],
                     int(spec.get("weighted", False)),
                     WEIGHTS_AS[spec.get("weights_as", "distance")],
                     int(spec.get("normalized", False)), 0)
    return _undefined_as_empty(values)


@register("edge_betweenness")
def _edge_betweenness(g, spec, fx):
    lib = _lib()
    m = len(fx.src)
    return _vector(lib.net_edge_betweenness, g.handle, m, int(spec.get("weighted", False)),
                   WEIGHTS_AS[spec.get("weights_as", "distance")],
                   int(spec.get("normalized", False)), 0)


@register("information_centrality")
def _information(g, spec, fx):
    lib = _lib()
    return _undefined_as_empty(_vector(lib.net_information_centrality, g.handle, g.n,
                                       int(spec.get("weighted", False))))


@register("bridging_coefficient")
def _bridging_coefficient(g, spec, fx):
    lib = _lib()
    return _undefined_as_empty(_vector(lib.net_bridging_coefficient, g.handle, g.n))


@register("bridging_centrality")
def _bridging_centrality(g, spec, fx):
    """Betweenness (igraph convention, unweighted) times the bridging
    coefficient, as the R function forms it."""
    lib = _lib()
    b = _vector(lib.net_betweenness, g.handle, g.n, 0, 0, 0, 0)
    c = _vector(lib.net_bridging_coefficient, g.handle, g.n)
    return _undefined_as_empty([x * y for x, y in zip(b, c)])


@register("hierarchy")
def _hierarchy(g, spec, fx):
    lib = _lib()
    return _undefined_as_empty(_vector(lib.net_structural_holes, g.handle, g.n, 3,
                                       int(spec.get("weighted", False)), 0))


@register("centralization_eigenvector")
def _centralization_eigenvector(g, spec, fx):
    lib = _lib()
    value = _vector(lib.net_centralization, g.handle, 1, 3, MODES["all"], 1,
                    int(spec.get("normalized", True)), 0, 0)
    return _undefined_as_empty(value)


# ---- Positions and roles -----------------------------------------------------

@register("regular_equivalence")
def _regular(g, spec, fx):
    lib = _lib()
    n = g.n
    seed = spec.get("seed")
    seed_arr = (ctypes.c_int64 * max(n, 1))(*partition_for(seed, n)) if seed else None
    out = (ctypes.c_int32 * max(n, 1))()
    rounds = ctypes.c_int32()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_regular_equivalence(g.handle, seed_arr, n if seed else 0,
                                       int(spec.get("weighted", False)), out,
                                       ctypes.byref(rounds), err, ERR_LEN), err)
    return [float(x) for x in list(out)[:n]]


@register("automorphic_equivalence")
def _automorphic(g, spec, fx):
    lib = _lib()
    n = g.n
    out = (ctypes.c_int32 * max(n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_automorphic_equivalence(g.handle, 0, out, err, ERR_LEN), err)
    return [float(x) for x in list(out)[:n]]


def _block_model(g, spec):
    lib = _lib()
    n = g.n
    labels = partition_for(spec["partition"], n)
    k = len(set(labels))
    arr = (ctypes.c_int64 * max(n, 1))(*labels)
    dens = (ctypes.c_double * max(k * k, 1))()
    sums = (ctypes.c_double * max(k * k, 1))()
    pairs = (ctypes.c_double * max(k * k, 1))()
    summary = (ctypes.c_double * 2)()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_block_model(g.handle, arr, n, k, int(spec.get("weighted", False)),
                               int(spec.get("loops", False)), dens, sums, pairs, summary,
                               err, ERR_LEN), err)
    return k, list(dens)[:k * k], list(summary)


@register("block_density")
def _block_density(g, spec, fx):
    k, dens, _ = _block_model(g, spec)
    rows = []
    for r in range(k):
        row = [dens[r + k * c] for c in range(k)]
        # An undefined density (a one-node diagonal block) is written as -1
        # inside a row, since densities are never negative.
        rows.append(tuple(-1.0 if math.isnan(x) else x for x in row))
    return rows


@register("block_fit")
def _block_fit(g, spec, fx):
    _, _, summary = _block_model(g, spec)
    return _undefined_as_empty(summary)


# ---- Scaling -----------------------------------------------------------------

def _mds(fx, k):
    lib = _lib()
    d = mds_input(fx)
    n = fx.n
    flat = _doubles([d[i][j] for j in range(n) for i in range(n)])
    pts = (ctypes.c_double * max(n * k, 1))()
    dims = ctypes.c_int32()
    eig = (ctypes.c_double * max(n, 1))()
    gof = (ctypes.c_double * 2)()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_classical_mds(flat, n, k, pts, ctypes.byref(dims), eig, gof, err,
                                 ERR_LEN), err)
    cols = dims.value
    points = [tuple(pts[i + n * c] for c in range(cols)) for i in range(n)]
    return d, points, list(eig)[:n], list(gof), cols


@register("classical_mds")
def _classical_mds(g, spec, fx):
    k = spec["k"]
    if fx.n <= k:
        return []
    what = spec["output"]
    d, points, eig, gof, cols = _mds(fx, k)
    if what == "points":
        return points
    if what == "eigenvalues":
        return eig
    if what == "gof":
        return gof
    if what == "stress":
        lib = _lib()
        n = fx.n
        flat = _doubles([d[i][j] for j in range(n) for i in range(n)])
        p = _doubles([points[i][c] for c in range(cols) for i in range(n)])
        out = ctypes.c_double()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(lib.net_mds_stress(flat, n, p, cols, ctypes.byref(out), err, ERR_LEN), err)
        return _undefined_as_empty([out.value])
    raise ValueError(f"unknown classical_mds output {what!r}")


# ---- Geospatial --------------------------------------------------------------

def _crs(spec):
    return CRS[spec["crs"]], EARTH_RADIUS_M


@register("edge_lengths")
def _edge_lengths(g, spec, fx):
    lib = _lib()
    x, y = coordinates(spec["crs"], fx.n)
    code, radius = _crs(spec)
    return _vector(lib.net_edge_lengths, g.handle, len(fx.src), _doubles(x), _doubles(y),
                   fx.n, code, radius)


@register("proximity")
def _proximity(g, spec, fx):
    lib = _lib()
    x, y = coordinates(spec["crs"], fx.n)
    code, radius = _crs(spec)
    h = ctypes.c_void_p()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_proximity(_doubles(x), _doubles(y), fx.n, code, radius,
                             float(spec.get("r", 0.0)), int(spec.get("k", 0)), 0,
                             ctypes.byref(h), err, ERR_LEN), err)
    try:
        m = lib.net_edgelist_size(h)
        a = (ctypes.c_int64 * max(m, 1))()
        b = (ctypes.c_int64 * max(m, 1))()
        d = (ctypes.c_double * max(m, 1))()
        _check(lib.net_edgelist_copy(h, a, b, d, err, ERR_LEN), err)
        return [(float(a[i]), float(b[i]), d[i]) for i in range(m)]
    finally:
        lib.net_edgelist_free(h)


@register("group_geometry")
def _group_geometry(g, spec, fx):
    lib = _lib()
    n = fx.n
    x, y = coordinates(spec["crs"], n)
    labels = partition_for(spec["partition"], n)
    k = len(set(labels))
    code, radius = _crs(spec)
    outs = [(ctypes.c_double * max(k, 1))() for _ in range(5)]
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_group_geometry(_doubles(x), _doubles(y),
                                  (ctypes.c_int64 * max(n, 1))(*labels), n, code, radius, k,
                                  *outs, err, ERR_LEN), err)
    return [tuple(o[j] for o in outs) for j in range(k)]
