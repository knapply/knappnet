"""Engine runners for shortest paths (crates/capi/src/paths.rs).

Every runner checks the paths it receives before reducing them to numbers:
each step of a returned path must be an arc of the path graph, and the
path's length must equal the reported distance. A violation raises, so a
wrong path cannot pass as a right distance.

Markers: an unreachable node is written as INFINITE = -1e9 and a negative
cycle as NEGATIVE_CYCLE = -2e9 (distances here can be negative, so the -1
and -999 markers of other areas do not work). Collections of node lists
(predecessor lists, sets of paths) are encoded as groups of integers so
that the `sets` comparison, which sorts within groups, keeps their
meaning (see `_encode_*`).

Manifest engine options:
  fn: paths_sssp            distances from every source, one tuple per source
  fn: paths_multi_source    distances from the nearest of nodes 0 and n - 1
  fn: paths_predecessors    predecessor lists of a search from node 0
  fn: paths_all_pairs       distance matrix from the all-pairs call
  fn: paths_bidirectional   pair distances by bidirectional search
  fn: paths_astar           pair distances by A* with heuristic
                            `factor` * (distance to the target)
  fn: paths_negative        distances from node 0 with lengths of any sign
                            (`method`: bellman_ford or goldberg_radzik)
  fn: paths_negative_cycle  1 if there is a negative cycle anywhere, else 0
  fn: paths_johnson         distance matrix by Johnson's algorithm
  fn: paths_floyd_warshall  distance matrix by Floyd-Warshall
  fn: paths_count           number of shortest paths from node 0
  fn: paths_all_shortest    every shortest path from node 0 to node n - 1
  fn: paths_k_shortest      lengths of the k shortest simple paths from node
                            0 to node n - 1
  mode: out | in | all; lengths: unit | weights | strength | potential |
  shift (see edge_lengths); via_weights (weights and strengths through
  the graph's weights instead of per-edge lengths); cutoff; collapse
  (parallel arcs); tie: relative | exact.
"""

from __future__ import annotations

import ctypes
import math

from engine import ERR_LEN, EngineError, load, register

INFINITE = -1e9
NEGATIVE_CYCLE = -2e9
MODES = {"out": 1, "in": 2, "all": 3}
TIES = {"relative": 0, "exact": 1}

_declared = False


def edge_lengths(fx, kind: str) -> list[float]:
    """Edge lengths for a manifest `lengths` option. The oracles repeat
    these formulas.

    potential: w + p(u) - p(v) with p(v) = v mod 5 (directed networks only;
               shortest paths are those of w, and no cycle is negative)
    shift:     w - 1.5 on every third edge (index divisible by 3), w on the
               others, which makes some cycles negative
    with w the weight (1 on unweighted networks)."""
    m = len(fx.src)
    w = fx.weights if fx.weights is not None else [1.0] * m
    if kind in ("unit",):
        return [1.0] * m
    if kind == "weights":
        return list(w)
    if kind == "strength":
        return [1.0 / x for x in w]
    if kind == "potential":
        return [w[e] + fx.src[e] % 5 - fx.tgt[e] % 5 for e in range(m)]
    if kind == "shift":
        return [w[e] - (1.5 if e % 3 == 0 else 0.0) for e in range(m)]
    raise ValueError(f"unknown lengths {kind!r}")


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double
    f64p, i64p, i32p = ctypes.POINTER(f64), ctypes.POINTER(i64), ctypes.POINTER(i32)
    h, cp, sz = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_size_t
    hp = ctypes.POINTER(h)
    sig = {
        "net_path_graph_new": [h, i32, i32, i32, f64p, i64, i32, hp, cp, sz],
        "net_path_lists_copy": [h, i64p, i64p, f64p, cp, sz],
        "net_shortest_paths": [h, i64p, i64, f64, i64, i32, f64p, i64p, i64p, i64p, cp, sz],
        "net_shortest_path_predecessors": [h, i64p, i64, f64, i32, i32, f64p, hp, cp, sz],
        "net_all_pairs_shortest_paths": [h, i64p, i64, f64, i32, f64p, i64p, cp, sz],
        "net_bidirectional_shortest_path": [h, i64, i64, hp, cp, sz],
        "net_astar_path": [h, i64, i64, f64p, f64, hp, cp, sz],
        "net_negative_shortest_paths": [h, i64p, i64, i32, i32, f64p, i64p, hp, hp, cp, sz],
        "net_negative_cycle": [h, i64, hp, cp, sz],
        "net_johnson": [h, i64p, i64, i32, f64p, i64p, hp, cp, sz],
        "net_floyd_warshall": [h, i32, f64p, i64p, i32p, cp, sz],
        "net_count_shortest_paths": [h, i64, i32, f64p, cp, sz],
        "net_all_shortest_paths": [h, i64, i64, i32, i64, hp, cp, sz],
        "net_k_shortest_paths": [h, i64, i64, i64, hp, cp, sz],
    }
    for name, args in sig.items():
        fn = getattr(lib, name)
        fn.argtypes = args
        fn.restype = i32
    for name in ("net_path_lists_count", "net_path_lists_total", "net_path_graph_arc_count"):
        getattr(lib, name).argtypes = [h]
        getattr(lib, name).restype = i64
    lib.net_path_graph_free.argtypes = [h]
    lib.net_path_lists_free.argtypes = [h]
    _declared = True
    return lib


def _check(status: int, err) -> None:
    if status:
        raise EngineError(err.value.decode())


def _mark(x: float) -> float:
    return INFINITE if math.isinf(x) else x


class PathGraph:
    """A path graph for a fixture and manifest options, with the lengths
    and arcs it was built from (for checking returned paths)."""

    def __init__(self, g, spec, fx):
        self.lib = _declare()
        self.n = g.n
        kind = spec.get("lengths", "unit")
        self.lengths = edge_lengths(fx, kind)
        mode = spec.get("mode", "out")
        self.handle = ctypes.c_void_p()
        err = ctypes.create_string_buffer(ERR_LEN)
        if spec.get("via_weights") and kind in ("weights", "strength"):
            # Through the graph's weights and weights_as rather than lengths.
            args = (1, {"weights": 0, "strength": 1}[kind], None, -1)
        else:
            lens = (ctypes.c_double * max(len(self.lengths), 1))(*self.lengths)
            args = (0, 0, lens, len(self.lengths))
        _check(self.lib.net_path_graph_new(g.handle, MODES[mode], *args,
                                           int(spec.get("collapse", False)),
                                           ctypes.byref(self.handle), err, ERR_LEN), err)
        # Shortest arc between each ordered pair, as the search sees them.
        self.arc = {}
        for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
            if math.isinf(self.lengths[e]):
                continue
            pairs = []
            if not fx.directed or mode == "all":
                pairs = [(a, b), (b, a)]
            elif mode == "out":
                pairs = [(a, b)]
            else:
                pairs = [(b, a)]
            for p in pairs:
                if p not in self.arc or self.lengths[e] < self.arc[p]:
                    self.arc[p] = self.lengths[e]
        self.fx = fx
        self.mode = mode

    def __del__(self):
        if getattr(self, "handle", None):
            self.lib.net_path_graph_free(self.handle)
            self.handle = None

    def path_length(self, path) -> float:
        total = 0.0
        for u, v in zip(path, path[1:]):
            if (u, v) not in self.arc:
                raise EngineError(f"returned path {path} uses {u} -> {v}, which is not an arc")
            total += self.arc[(u, v)]
        return total

    def lists(self, h) -> list[tuple[list[int], float]]:
        lib = self.lib
        count = lib.net_path_lists_count(h)
        total = lib.net_path_lists_total(h)
        off = (ctypes.c_int64 * (count + 1))()
        nodes = (ctypes.c_int64 * max(total, 1))()
        vals = (ctypes.c_double * max(count, 1))()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(lib.net_path_lists_copy(h, off, nodes, vals, err, ERR_LEN), err)
        lib.net_path_lists_free(h)
        return [(list(nodes[off[i]:off[i + 1]]), vals[i]) for i in range(count)]

    def call_lists(self, fn, *args):
        out = ctypes.c_void_p()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(fn(self.handle, *args, ctypes.byref(out), err, ERR_LEN), err)
        return self.lists(out)


def _close(a: float, b: float) -> bool:
    return a == b or abs(a - b) <= 1e-9 * max(1.0, abs(a), abs(b))


def _check_tree(pg: PathGraph, dist, parent, sources) -> None:
    """Every reached non-source node's parent arc closes its distance."""
    for v in range(pg.n):
        if math.isinf(dist[v]) or v in sources:
            continue
        p = parent[v]
        if p < 0:
            raise EngineError(f"node {v} is reached but has no parent")
        if (p, v) not in pg.arc or not _close(dist[p] + pg.arc[(p, v)], dist[v]):
            raise EngineError(f"parent {p} of node {v} does not give its distance")


def _search(pg: PathGraph, sources, cutoff=math.inf, reverse=False):
    n = pg.n
    src = (ctypes.c_int64 * max(len(sources), 1))(*sources)
    dist = (ctypes.c_double * max(n, 1))()
    par = (ctypes.c_int64 * max(n, 1))()
    edge = (ctypes.c_int64 * max(n, 1))()
    origin = (ctypes.c_int64 * max(n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(pg.lib.net_shortest_paths(pg.handle, src, len(sources), cutoff, -1, int(reverse),
                                     dist, par, edge, origin, err, ERR_LEN), err)
    return list(dist)[:n], list(par)[:n], list(edge)[:n], list(origin)[:n]


@register("paths_sssp")
def _sssp(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    cutoff = float(spec.get("cutoff", math.inf))
    rows = []
    for s in range(g.n):
        dist, par, edge, _ = _search(pg, [s], cutoff)
        _check_tree(pg, dist, par, {s})
        for v in range(g.n):
            if par[v] >= 0 and not _close(pg.lengths[edge[v]], pg.arc[(par[v], v)]):
                raise EngineError(f"parent edge of {v} is not the shortest arc from {par[v]}")
        rows.append(tuple(_mark(x) for x in dist))
    return rows


@register("paths_multi_source")
def _multi_source(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    sources = sorted({0, g.n - 1})
    dist, par, _, origin = _search(pg, sources)
    _check_tree(pg, dist, par, set(sources))
    for v in range(g.n):
        if not math.isinf(dist[v]) and origin[v] not in sources:
            raise EngineError(f"node {v} has origin {origin[v]}")
    return [_mark(x) for x in dist]


def _encode_lists(n: int, lists) -> list[tuple]:
    """Node v's list as the group (v (n + 1) + u for u in the list), or
    (v (n + 1) + n) when empty."""
    return [tuple(v * (n + 1) + u for u in members) or (v * (n + 1) + n,)
            for v, members in enumerate(lists)]


def _encode_paths(n: int, paths) -> list[tuple]:
    """A path as the group (i n + node_i), which sorts back into path order."""
    return [tuple(i * n + v for i, v in enumerate(p)) for p in paths]


@register("paths_predecessors")
def _predecessors(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    n = g.n
    src = (ctypes.c_int64 * 1)(0)
    dist = (ctypes.c_double * max(n, 1))()
    lists = pg.call_lists(pg.lib.net_shortest_path_predecessors, src, 1,
                          float(spec.get("cutoff", math.inf)), 0,
                          TIES[spec.get("tie", "exact")], dist)
    return _encode_lists(n, [members for members, _ in lists])


@register("paths_all_pairs")
def _all_pairs(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    n = g.n
    dist = (ctypes.c_double * max(n * n, 1))()
    par = (ctypes.c_int64 * max(n * n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(pg.lib.net_all_pairs_shortest_paths(pg.handle, None, -1,
                                               float(spec.get("cutoff", math.inf)), 0,
                                               dist, par, err, ERR_LEN), err)
    rows = []
    for s in range(n):
        d = [dist[v * n + s] for v in range(n)]
        p = [par[v * n + s] for v in range(n)]
        _check_tree(pg, d, p, {s})
        rows.append(tuple(_mark(x) for x in d))
    return rows


@register("paths_bidirectional")
def _bidirectional(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    rows = []
    for s in range(g.n):
        row = []
        for t in range(g.n):
            found = pg.call_lists(pg.lib.net_bidirectional_shortest_path, s, t)
            if not found:
                row.append(INFINITE)
                continue
            path, length = found[0]
            if path[0] != s or path[-1] != t or not _close(pg.path_length(path), length):
                raise EngineError(f"bidirectional path {path} for ({s}, {t}) is inconsistent")
            row.append(length)
        rows.append(tuple(row))
    return rows


@register("paths_astar")
def _astar(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    factor = float(spec.get("factor", 0.5))
    n = g.n
    rows = [[INFINITE] * n for _ in range(n)]
    for t in range(n):
        to_t, *_ = _search(pg, [t], reverse=True)
        h = (ctypes.c_double * max(n, 1))(*[0.0 if math.isinf(d) else factor * d for d in to_t])
        for s in range(n):
            found = pg.call_lists(pg.lib.net_astar_path, s, t, h, math.inf)
            if found:
                path, length = found[0]
                if path[0] != s or path[-1] != t or not _close(pg.path_length(path), length):
                    raise EngineError(f"A* path {path} for ({s}, {t}) is inconsistent")
                rows[s][t] = length
    return [tuple(r) for r in rows]


def _cycle_ok(pg: PathGraph, cycle) -> None:
    if len(cycle) < 2 or cycle[0] != cycle[-1] or pg.path_length(cycle) >= 0:
        raise EngineError(f"{cycle} is not a negative cycle")


@register("paths_negative")
def _negative(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    n = g.n
    method = {"bellman_ford": 0, "goldberg_radzik": 1}[spec.get("method", "bellman_ford")]
    src = (ctypes.c_int64 * 1)(0)
    dist = (ctypes.c_double * max(n, 1))()
    par = (ctypes.c_int64 * max(n, 1))()
    cycle = pg.call_lists(pg.lib.net_negative_shortest_paths, src, 1, method,
                          TIES[spec.get("tie", "exact")], dist, par, None)
    if cycle:
        _cycle_ok(pg, cycle[0][0])
        return [NEGATIVE_CYCLE]
    d, p = list(dist)[:n], list(par)[:n]
    _check_tree(pg, d, p, {0})
    return [_mark(x) for x in d]


@register("paths_negative_cycle")
def _negative_cycle(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    found = pg.call_lists(pg.lib.net_negative_cycle, -1)
    if found:
        _cycle_ok(pg, found[0][0])
    return [1.0 if found else 0.0]


def _matrix_rows(pg, dist, par, n):
    rows = []
    for s in range(n):
        d = [dist[v * n + s] for v in range(n)]
        if par is not None:
            p = [par[v * n + s] for v in range(n)]
            _check_tree(pg, d, p, {s})
        rows.append(tuple(_mark(x) for x in d))
    return rows


@register("paths_johnson")
def _johnson(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    n = g.n
    dist = (ctypes.c_double * max(n * n, 1))()
    par = (ctypes.c_int64 * max(n * n, 1))()
    cycle = pg.call_lists(pg.lib.net_johnson, None, -1, 0, dist, par)
    if cycle:
        _cycle_ok(pg, cycle[0][0])
        return [NEGATIVE_CYCLE]
    return _matrix_rows(pg, dist, par, n)


@register("paths_floyd_warshall")
def _floyd_warshall(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    n = g.n
    dist = (ctypes.c_double * max(n * n, 1))()
    pred = (ctypes.c_int64 * max(n * n, 1))()
    neg = ctypes.c_int32()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(pg.lib.net_floyd_warshall(pg.handle, 0, dist, pred, ctypes.byref(neg), err, ERR_LEN),
           err)
    if neg.value:
        return [NEGATIVE_CYCLE]
    # pred[u][v] is the node before v on the path from u: check it closes
    # the distance like a parent.
    for u in range(n):
        d = [dist[v * n + u] for v in range(n)]
        p = [pred[v * n + u] for v in range(n)]
        _check_tree(pg, d, p, {u})
    return _matrix_rows(pg, dist, None, n)


@register("paths_count")
def _count(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    n = g.n
    out = (ctypes.c_double * max(n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(pg.lib.net_count_shortest_paths(pg.handle, 0, TIES[spec.get("tie", "exact")], out,
                                           err, ERR_LEN), err)
    return list(out)[:n]


@register("paths_all_shortest")
def _all_shortest(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    n = g.n
    found = pg.call_lists(pg.lib.net_all_shortest_paths, 0, n - 1,
                          TIES[spec.get("tie", "exact")], 100000)
    for path, length in found:
        if path[0] != 0 or path[-1] != n - 1 or not _close(pg.path_length(path), length):
            raise EngineError(f"shortest path {path} is inconsistent")
    return _encode_paths(n, [p for p, _ in found])


@register("paths_k_shortest")
def _k_shortest(g, spec, fx):
    pg = PathGraph(g, spec, fx)
    n = g.n
    k = int(spec.get("k", 10))
    found = pg.call_lists(pg.lib.net_k_shortest_paths, 0, n - 1, k)
    seen = set()
    for path, length in found:
        if (path[0] != 0 or path[-1] != n - 1 or len(set(path)) != len(path)
                or not _close(pg.path_length(path), length) or tuple(path) in seen):
            raise EngineError(f"simple path {path} is inconsistent")
        seen.add(tuple(path))
    return [length for _, length in found]
