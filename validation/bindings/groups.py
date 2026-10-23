"""Engine runners for cohesive subsets, community partitions, topography
(core-periphery, compactness, small worlds), and random graphs."""

from __future__ import annotations

import ctypes

from engine import ERR_LEN, EngineError, load, register

_declared = False


def _lib():
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64, u64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double, ctypes.c_uint64
    i64p, f64p = ctypes.POINTER(i64), ctypes.POINTER(f64)
    i32p = ctypes.POINTER(i32)
    handle = ctypes.c_void_p
    hp = ctypes.POINTER(handle)
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    sigs = {
        "net_node_sets_count": ([handle], i64),
        "net_node_sets_total": ([handle], i64),
        "net_node_sets_copy": ([handle, i64p, i64p, cp, sz], i32),
        "net_node_sets_free": ([handle], None),
        "net_bicomponents": ([handle, hp, cp, sz], i32),
        "net_n_cliques": ([handle, i64, i64, i32, i32, hp, cp, sz], i32),
        "net_k_plexes": ([handle, i64, i64, f64, i32, hp, cp, sz], i32),
        "net_lambda_connectivity": ([handle, f64p, cp, sz], i32),
        "net_factions": ([handle, i64, i32, i64, i64, i64, u64, i32, i64p, f64p, f64p, cp, sz],
                         i32),
        "net_girvan_newman": ([handle, i32, i64p, f64p, i64p, i64p, f64p, i64p, cp, sz], i32),
        "net_leading_eigen": ([handle, i32, i64p, f64p, i64p, cp, sz], i32),
        "net_fast_greedy": ([handle, f64, i32, i64p, f64p, cp, sz], i32),
        "net_core_periphery": ([handle, i32, i32, i64, i64, i64, u64, i32, i32p, f64p, cp, sz],
                               i32),
        "net_core_periphery_continuous": ([handle, i32, i64, i64, u64, i32, f64p, f64p, i64p, cp, sz],
                                          i32),
        "net_compactness": ([handle, i32, i32, i32, i32, f64p, cp, sz], i32),
        "net_small_world": ([handle, i32, i64, u64, i32, f64p, cp, sz], i32),
        "net_edge_lists_count": ([handle], i64),
        "net_edge_lists_total": ([handle], i64),
        "net_edge_lists_copy": ([handle, i64p, i64p, i64p, cp, sz], i32),
        "net_edge_lists_free": ([handle], None),
        "net_random_gnp": ([i64, f64, i32, i32, i64, u64, i32, hp, cp, sz], i32),
        "net_random_gnm": ([i64, i64, i32, i32, i64, u64, i32, hp, cp, sz], i32),
        "net_random_configuration": ([i64, i64p, i64p, i32, i64, i64, u64, i32, hp, cp, sz], i32),
        "net_random_uman": ([i64, i64, i64, i64, u64, i32, hp, cp, sz], i32),
        "net_random_small_world": ([i64, i64, f64, i64, u64, i32, hp, cp, sz], i32),
        "net_random_preferential": ([i64, i64, i64, u64, i32, hp, cp, sz], i32),
        "net_random_sbm": ([i64, i64p, f64p, i32, i32, i64, u64, i32, hp, cp, sz], i32),
    }
    for name, (args, res) in sigs.items():
        fn = getattr(lib, name)
        fn.argtypes = args
        fn.restype = res
    _declared = True
    return lib


def _check(status, err):
    if status:
        raise EngineError(err.value.decode())


def _err():
    return ctypes.create_string_buffer(ERR_LEN)


def _sets(call) -> list[tuple[int, ...]]:
    """Runs call(out_handle_ref, err) and copies the node sets out."""
    lib = _lib()
    err = _err()
    h = ctypes.c_void_p()
    _check(call(ctypes.byref(h), err), err)
    try:
        count = lib.net_node_sets_count(h)
        total = lib.net_node_sets_total(h)
        offsets = (ctypes.c_int64 * (count + 1))()
        members = (ctypes.c_int64 * max(total, 1))()
        _check(lib.net_node_sets_copy(h, offsets, members, err, ERR_LEN), err)
    finally:
        lib.net_node_sets_free(h)
    return [tuple(members[offsets[i]:offsets[i + 1]]) for i in range(count)]


def _edge_lists(call) -> list[list[tuple[int, int]]]:
    lib = _lib()
    err = _err()
    h = ctypes.c_void_p()
    _check(call(ctypes.byref(h), err), err)
    try:
        count = lib.net_edge_lists_count(h)
        total = lib.net_edge_lists_total(h)
        offsets = (ctypes.c_int64 * (count + 1))()
        a = (ctypes.c_int64 * max(total, 1))()
        b = (ctypes.c_int64 * max(total, 1))()
        _check(lib.net_edge_lists_copy(h, offsets, a, b, err, ERR_LEN), err)
    finally:
        lib.net_edge_lists_free(h)
    return [list(zip(a[offsets[i]:offsets[i + 1]], b[offsets[i]:offsets[i + 1]]))
            for i in range(count)]


def bicomponents(g):
    lib = _lib()
    return _sets(lambda out, err: lib.net_bicomponents(g.handle, out, err, ERR_LEN))


def n_cliques(g, n: int, min_size: int = 0, clans: bool = False, threads: int = 0):
    lib = _lib()
    return _sets(lambda out, err: lib.net_n_cliques(g.handle, n, min_size, int(clans), threads,
                                                    out, err, ERR_LEN))


def k_plexes(g, k: int, min_size: int, max_steps: float = 1e8, threads: int = 0):
    lib = _lib()
    return _sets(lambda out, err: lib.net_k_plexes(g.handle, k, min_size, max_steps, threads,
                                                   out, err, ERR_LEN))


def lambda_connectivity(g) -> list[list[float]]:
    lib = _lib()
    err = _err()
    out = (ctypes.c_double * max(g.n * g.n, 1))()
    _check(lib.net_lambda_connectivity(g.handle, out, err, ERR_LEN), err)
    return [[out[j * g.n + i] for j in range(g.n)] for i in range(g.n)]


def factions(g, k: int, fit: str = "errors", starts: int = 20, max_stall: int = 200,
             tenure: int = 0, seed: int = 1, threads: int = 0):
    lib = _lib()
    err = _err()
    memb = (ctypes.c_int64 * max(g.n, 1))()
    errors, corr = ctypes.c_double(), ctypes.c_double()
    _check(lib.net_factions(g.handle, k, {"errors": 0, "correlation": 1}[fit], starts,
                            max_stall, tenure, seed, threads, memb, ctypes.byref(errors),
                            ctypes.byref(corr), err, ERR_LEN), err)
    return list(memb)[:g.n], errors.value, corr.value


def girvan_newman(g, m: int, threads: int = 0):
    lib = _lib()
    err = _err()
    memb = (ctypes.c_int64 * max(g.n, 1))()
    q = ctypes.c_double()
    rf = (ctypes.c_int64 * max(m, 1))()
    rt = (ctypes.c_int64 * max(m, 1))()
    path = (ctypes.c_double * (m + 1))()
    k = ctypes.c_int64()
    _check(lib.net_girvan_newman(g.handle, threads, memb, ctypes.byref(q), rf, rt, path,
                                 ctypes.byref(k), err, ERR_LEN), err)
    removed = list(zip(rf[:k.value], rt[:k.value]))
    return list(memb)[:g.n], q.value, removed, list(path)[:k.value + 1]


def leading_eigen(g, weighted: bool = False):
    lib = _lib()
    err = _err()
    memb = (ctypes.c_int64 * max(g.n, 1))()
    q = ctypes.c_double()
    splits = ctypes.c_int64()
    _check(lib.net_leading_eigen(g.handle, int(weighted), memb, ctypes.byref(q),
                                 ctypes.byref(splits), err, ERR_LEN), err)
    return list(memb)[:g.n], q.value


def fast_greedy(g, resolution: float = 1.0, weighted: bool = False):
    lib = _lib()
    err = _err()
    memb = (ctypes.c_int64 * max(g.n, 1))()
    q = ctypes.c_double()
    _check(lib.net_fast_greedy(g.handle, resolution, int(weighted), memb, ctypes.byref(q),
                               err, ERR_LEN), err)
    return list(memb)[:g.n], q.value


LINKS = {"missing": 0, "ones": 1, "zeros": 2}


def core_periphery(g, links: str = "missing", weighted: bool = False, starts: int = 20,
                   max_stall: int = 200, tenure: int = 0, seed: int = 1, threads: int = 0):
    lib = _lib()
    err = _err()
    core = (ctypes.c_int32 * max(g.n, 1))()
    corr = ctypes.c_double()
    _check(lib.net_core_periphery(g.handle, LINKS[links], int(weighted), starts, max_stall,
                                  tenure, seed, threads, core, ctypes.byref(corr), err,
                                  ERR_LEN), err)
    return [bool(c) for c in list(core)[:g.n]], corr.value


def core_periphery_continuous(g, weighted: bool = False, max_iter: int = 1000,
                              starts: int = 20, seed: int = 1, threads: int = 1):
    lib = _lib()
    err = _err()
    c = (ctypes.c_double * max(g.n, 1))()
    corr = ctypes.c_double()
    iters = ctypes.c_int64()
    _check(lib.net_core_periphery_continuous(g.handle, int(weighted), max_iter, starts, seed,
                                             threads, c, ctypes.byref(corr),
                                             ctypes.byref(iters), err, ERR_LEN), err)
    return list(c)[:g.n], corr.value


def compactness(g, mode: str = "out", weighted: bool = False, weights_as: str = "distance",
                threads: int = 0) -> float:
    lib = _lib()
    err = _err()
    out = ctypes.c_double()
    _check(lib.net_compactness(g.handle, {"out": 1, "in": 2, "all": 3}[mode], int(weighted),
                               {"distance": 0, "strength": 1}[weights_as], threads,
                               ctypes.byref(out), err, ERR_LEN), err)
    return out.value


def small_world(g, local: bool = False, reps: int = 20, seed: int = 1, threads: int = 0):
    lib = _lib()
    err = _err()
    out = (ctypes.c_double * 5)()
    _check(lib.net_small_world(g.handle, int(local), reps, seed, threads, out, err, ERR_LEN),
           err)
    return list(out)


def random_gnp(n, p, directed=False, loops=False, reps=1, seed=1, threads=0):
    lib = _lib()
    return _edge_lists(lambda out, err: lib.net_random_gnp(n, p, int(directed), int(loops), reps,
                                                           seed, threads, out, err, ERR_LEN))


def random_gnm(n, m, directed=False, loops=False, reps=1, seed=1, threads=0):
    lib = _lib()
    return _edge_lists(lambda out, err: lib.net_random_gnm(n, m, int(directed), int(loops), reps,
                                                           seed, threads, out, err, ERR_LEN))


def random_configuration(degrees, in_degrees=None, simple=False, max_tries=1000, reps=1,
                         seed=1, threads=0):
    lib = _lib()
    n = len(degrees)
    d = (ctypes.c_int64 * max(n, 1))(*degrees)
    di = None if in_degrees is None else (ctypes.c_int64 * max(n, 1))(*in_degrees)
    return _edge_lists(lambda out, err: lib.net_random_configuration(
        n, d, di, int(simple), max_tries, reps, seed, threads, out, err, ERR_LEN))


def random_uman(n, mutual, asym, reps=1, seed=1, threads=0):
    lib = _lib()
    return _edge_lists(lambda out, err: lib.net_random_uman(n, mutual, asym, reps, seed, threads,
                                                            out, err, ERR_LEN))


def random_small_world(n, k, p, reps=1, seed=1, threads=0):
    lib = _lib()
    return _edge_lists(lambda out, err: lib.net_random_small_world(n, k, p, reps, seed, threads,
                                                                   out, err, ERR_LEN))


def random_preferential(n, m, reps=1, seed=1, threads=0):
    lib = _lib()
    return _edge_lists(lambda out, err: lib.net_random_preferential(n, m, reps, seed, threads,
                                                                    out, err, ERR_LEN))


def random_sbm(sizes, probs, directed=False, loops=False, reps=1, seed=1, threads=0):
    lib = _lib()
    k = len(sizes)
    s = (ctypes.c_int64 * max(k, 1))(*sizes)
    p = (ctypes.c_double * max(k * k, 1))(*[x for row in probs for x in row])
    return _edge_lists(lambda out, err: lib.net_random_sbm(k, s, p, int(directed), int(loops),
                                                           reps, seed, threads, out, err,
                                                           ERR_LEN))


# ---- Runners --------------------------------------------------------------
# Undefined results (an engine error that the manifest declares, such as a
# core-periphery fit of a network without variation) are written as -1 by
# the runners and the oracles.


def _q(q: float) -> float:
    """Modularity lies in [-1/2, 1], so -1 marks it undefined (NaN)."""
    return -1.0 if q != q else q


@register("bicomponents")
def _bicomponents(g, spec, fx):
    return bicomponents(g)


@register("n_cliques")
def _n_cliques(g, spec, fx):
    return n_cliques(g, spec["n"], spec.get("min_size", 0), clans=False)


@register("n_clans")
def _n_clans(g, spec, fx):
    return n_cliques(g, spec["n"], spec.get("min_size", 0), clans=True)


@register("k_plexes")
def _k_plexes(g, spec, fx):
    return k_plexes(g, spec["k"], spec["min_size"])


@register("lambda_connectivity")
def _lambda(g, spec, fx):
    """The upper triangle of the connectivity matrix, row by row."""
    m = lambda_connectivity(g)
    return [m[i][j] for i in range(g.n) for j in range(i + 1, g.n)]


@register("lambda_sets")
def _lambda_sets(g, spec, fx):
    """Every lambda set: the classes of lambda >= k with at least two
    members, for every k from 1 to the largest connectivity."""
    m = lambda_connectivity(g)
    top = int(max((max(r) for r in m), default=0))
    out = set()
    for k in range(1, top + 1):
        seen = [False] * g.n
        for v in range(g.n):
            if seen[v]:
                continue
            cls = [u for u in range(g.n) if u == v or m[v][u] >= k]
            for u in cls:
                seen[u] = True
            if len(cls) > 1:
                out.add(tuple(cls))
    return sorted(out)


@register("factions_errors")
def _factions_errors(g, spec, fx):
    try:
        _, errors, _ = factions(g, spec["k"], "errors", seed=spec.get("seed", 1))
    except EngineError:
        return [-1.0]
    return [errors]


@register("factions_correlation")
def _factions_correlation(g, spec, fx):
    try:
        _, _, corr = factions(g, spec["k"], "correlation", seed=spec.get("seed", 1))
    except EngineError:
        return [-1.0]
    return [corr]


@register("girvan_newman")
def _girvan_newman(g, spec, fx):
    memb, q, _, _ = girvan_newman(g, len(fx.src))
    return memb if spec.get("output") == "membership" else [_q(q)]


@register("leading_eigen")
def _leading_eigen(g, spec, fx):
    memb, q = leading_eigen(g, spec.get("weighted", False))
    return memb if spec.get("output") == "membership" else [_q(q)]


@register("fast_greedy")
def _fast_greedy(g, spec, fx):
    memb, q = fast_greedy(g, spec.get("resolution", 1.0), spec.get("weighted", False))
    return memb if spec.get("output") == "membership" else [_q(q)]


@register("core_periphery")
def _core_periphery(g, spec, fx):
    try:
        _, corr = core_periphery(g, spec.get("links", "missing"), spec.get("weighted", False))
    except EngineError:
        return [-1.0]
    return [-1.0 if corr != corr else corr]


@register("core_periphery_continuous")
def _core_periphery_continuous(g, spec, fx):
    try:
        c, corr = core_periphery_continuous(g, spec.get("weighted", False))
    except EngineError:
        return [-1.0]
    return [corr] if spec.get("output", "fit") == "fit" else c


@register("compactness")
def _compactness(g, spec, fx):
    v = compactness(g, spec.get("mode", "out"), spec.get("weighted", False),
                    spec.get("weights_as", "distance"))
    return [-1.0 if v != v else v]


@register("small_world_observed")
def _small_world_observed(g, spec, fx):
    """The observed clustering and mean distance that enter the quotient."""
    out = small_world(g, spec.get("local", False), reps=1)
    return [-1.0 if v != v else v for v in out[:2]]
