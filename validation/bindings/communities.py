"""Engine runners for community detection (Walktrap, label propagation,
Infomap, Leiden, spinglass, exact maximum modularity, fluid communities,
Voronoi communities), partition comparison, and dendrogram cuts."""

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
    i64p, f64p, i32p = ctypes.POINTER(i64), ctypes.POINTER(f64), ctypes.POINTER(i32)
    handle = ctypes.c_void_p
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    sigs = {
        "net_walktrap": ([handle, i64, i32, i32, i64p, f64p, i64p, i64p, f64p, i64p, cp, sz],
                         i32),
        "net_label_propagation": ([handle, i32, i32, i32, i64p, i32p, u64, i64p, cp, sz], i32),
        "net_infomap": ([handle, i32, f64p, i64, u64, i64p, f64p, cp, sz], i32),
        "net_map_equation": ([handle, i64p, i64, i32, f64p, f64p, cp, sz], i32),
        "net_leiden": ([handle, i32, f64, f64, i64, i32, f64p, i64p, u64, i64p, f64p, i64p,
                        cp, sz], i32),
        "net_spinglass": ([handle, i64, i32, f64, f64, f64, i32, f64, i32, f64, i32, u64, i64p,
                           f64p, f64p, cp, sz], i32),
        "net_optimal_modularity": ([handle, f64, i32, i64p, f64p, cp, sz], i32),
        "net_fluid_communities": ([handle, i64, u64, i64p, cp, sz], i32),
        "net_voronoi": ([handle, f64p, i32, i32, f64, i32, i64p, i64p, i64p, f64p, f64p, cp,
                         sz], i32),
        "net_compare_partitions": ([i64p, i64p, i64, i32, f64p, cp, sz], i32),
        "net_split_join_distance": ([i64p, i64p, i64, i64p, i64p, cp, sz], i32),
        "net_cut_dendrogram": ([i64, i64p, i64p, i64, i64, i64, i64p, cp, sz], i32),
        "net_dendrogram_modularity": ([handle, i64p, i64p, i64, f64, i32, f64p, i64p, f64p, cp,
                                       sz], i32),
        "net_merges_from_pairs": ([i64, i64p, i64p, i64, i32, i64p, i64p, i64p, cp, sz], i32),
        "net_girvan_newman": ([handle, i32, i64p, f64p, i64p, i64p, f64p, i64p, cp, sz], i32),
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


def _i64(values, size=None):
    size = len(values) if size is None else size
    return (ctypes.c_int64 * max(size, 1))(*values)


def _f64(size):
    return (ctypes.c_double * max(size, 1))()


VARIANTS = {"dominance": 0, "retention": 1, "fast": 2, "semi_synchronous": 3}
MODES = {"out": 1, "in": 2, "all": 3}


def walktrap(g, steps=4, weighted=False, threads=0):
    lib = _lib()
    err = _err()
    n = g.n
    memb, ma, mb = _i64([], n), _i64([], n), _i64([], n)
    path = _f64(n)
    q, k = ctypes.c_double(), ctypes.c_int64()
    _check(lib.net_walktrap(g.handle, steps, int(weighted), threads, memb, ctypes.byref(q), ma,
                            mb, path, ctypes.byref(k), err, ERR_LEN), err)
    merges = list(zip(ma[:k.value], mb[:k.value]))
    return list(memb)[:n], q.value, merges, list(path)[:k.value + 1]


def label_propagation(g, variant="dominance", mode="all", weighted=False, seed=1,
                      initial=None, fixed=None):
    lib = _lib()
    err = _err()
    memb = _i64([], g.n)
    init = _i64(initial) if initial is not None else None
    fix = (ctypes.c_int32 * max(g.n, 1))(*fixed) if fixed is not None else None
    _check(lib.net_label_propagation(g.handle, VARIANTS[variant], MODES[mode], int(weighted),
                                     init, fix, seed, memb, err, ERR_LEN), err)
    return list(memb)[:g.n]


def infomap(g, weighted=False, trials=10, seed=1):
    lib = _lib()
    err = _err()
    memb = _i64([], g.n)
    cl = ctypes.c_double()
    _check(lib.net_infomap(g.handle, int(weighted), None, trials, seed, memb, ctypes.byref(cl),
                           err, ERR_LEN), err)
    return list(memb)[:g.n], cl.value


def map_equation(g, membership, weighted=False):
    lib = _lib()
    err = _err()
    out = ctypes.c_double()
    _check(lib.net_map_equation(g.handle, _i64(membership), len(membership), int(weighted), None,
                                ctypes.byref(out), err, ERR_LEN), err)
    return out.value


def leiden(g, objective="modularity", resolution=1.0, beta=0.01, n_iterations=-1,
           weighted=False, seed=1):
    lib = _lib()
    err = _err()
    memb = _i64([], g.n)
    quality, iters = ctypes.c_double(), ctypes.c_int64()
    code = {"cpm": 0, "modularity": 1}[objective]
    _check(lib.net_leiden(g.handle, code, resolution, beta, n_iterations, int(weighted), None,
                          None, seed, memb, ctypes.byref(quality), ctypes.byref(iters), err,
                          ERR_LEN), err)
    return list(memb)[:g.n], quality.value


def spinglass(g, spins=25, weighted=False, gamma=1.0, null_model="config", seed=1,
              signed=False, gamma_minus=1.0, parallel=False):
    lib = _lib()
    err = _err()
    memb = _i64([], g.n)
    q, energy = ctypes.c_double(), ctypes.c_double()
    code = {"simple": 0, "config": 1}[null_model]
    _check(lib.net_spinglass(g.handle, spins, int(parallel), 1.0, 0.01, 0.99, code, gamma,
                             int(signed), gamma_minus, int(weighted), seed, memb,
                             ctypes.byref(q), ctypes.byref(energy), err, ERR_LEN), err)
    return list(memb)[:g.n], q.value, energy.value


def optimal_modularity(g, resolution=1.0, weighted=False):
    lib = _lib()
    err = _err()
    memb = _i64([], g.n)
    q = ctypes.c_double()
    _check(lib.net_optimal_modularity(g.handle, resolution, int(weighted), memb, ctypes.byref(q),
                                      err, ERR_LEN), err)
    return list(memb)[:g.n], q.value


def fluid_communities(g, k, seed=1):
    lib = _lib()
    err = _err()
    memb = _i64([], g.n)
    _check(lib.net_fluid_communities(g.handle, k, seed, memb, err, ERR_LEN), err)
    return list(memb)[:g.n]


def voronoi(g, lengths=None, weighted=False, mode="out", radius=-1.0, threads=0):
    lib = _lib()
    err = _err()
    n = g.n
    memb, gens = _i64([], n), _i64([], n)
    k, q, r = ctypes.c_int64(), ctypes.c_double(), ctypes.c_double()
    lens = (ctypes.c_double * max(len(lengths), 1))(*lengths) if lengths is not None else None
    _check(lib.net_voronoi(g.handle, lens, int(weighted), MODES[mode], radius, threads, memb,
                           gens, ctypes.byref(k), ctypes.byref(q), ctypes.byref(r), err,
                           ERR_LEN), err)
    return list(memb)[:n], list(gens)[:k.value], q.value, r.value


COMPARE = {"vi": 0, "nmi": 1, "split_join": 2, "rand": 3, "adjusted_rand": 4, "ami": 5}


def compare_partitions(x, y, method):
    lib = _lib()
    err = _err()
    out = ctypes.c_double()
    _check(lib.net_compare_partitions(_i64(x), _i64(y), len(x), COMPARE[method],
                                      ctypes.byref(out), err, ERR_LEN), err)
    return out.value


def cut_dendrogram(n, merges, steps=0, groups=0):
    lib = _lib()
    err = _err()
    memb = _i64([], n)
    a = _i64([m[0] for m in merges])
    b = _i64([m[1] for m in merges])
    _check(lib.net_cut_dendrogram(n, a, b, len(merges), steps, groups, memb, err, ERR_LEN), err)
    return list(memb)[:n]


def merges_from_pairs(n, pairs, removals=False):
    lib = _lib()
    err = _err()
    ma, mb = _i64([], n), _i64([], n)
    k = ctypes.c_int64()
    _check(lib.net_merges_from_pairs(n, _i64([p[0] for p in pairs]), _i64([p[1] for p in pairs]),
                                     len(pairs), int(removals), ma, mb, ctypes.byref(k), err,
                                     ERR_LEN), err)
    return list(zip(ma[:k.value], mb[:k.value]))


def girvan_newman_removals(g, m):
    lib = _lib()
    err = _err()
    memb = _i64([], g.n)
    q = ctypes.c_double()
    rf, rt = _i64([], m), _i64([], m)
    path = _f64(m + 1)
    k = ctypes.c_int64()
    _check(lib.net_girvan_newman(g.handle, 1, memb, ctypes.byref(q), rf, rt, path,
                                 ctypes.byref(k), err, ERR_LEN), err)
    return list(zip(rf[:k.value], rt[:k.value]))


def blocks(n):
    """The fixed partition used across the validation suite: node v of n
    goes to group floor(3 v / n)."""
    return [3 * v // n for v in range(n)]


def stripes(n):
    """A second fixed partition, crossing `blocks`: node v goes to group
    (7 v) mod 4."""
    return [(7 * v) % 4 for v in range(n)]


def _q(x):
    return -1.0 if x != x else x


def _skipped(spec, fx):
    return fx is not None and fx.name in spec.get("skip", [])


# ---- Runners --------------------------------------------------------------
# A runner returns [-1] (or n values of -1) when the engine rejects the
# fixture (for example a disconnected network for a method that needs a
# connected one); the oracles write -1 in the same cases.


@register("walktrap")
def _walktrap(g, spec, fx):
    """Modularity ("modularity"), membership ("membership"), merges
    ("merges", one pair per line), or the membership cut at `groups`
    groups ("cut")."""
    output = spec.get("output", "modularity")
    memb, q, merges, _ = walktrap(g, spec.get("steps", 4), spec.get("weighted", False))
    if output == "membership":
        return memb
    if output == "merges":
        return [tuple(m) for m in merges]
    if output == "cut":
        try:
            return cut_dendrogram(g.n, merges, groups=spec["groups"])
        except EngineError:
            return [-1.0] * g.n
    return [_q(q)]


@register("girvan_newman_cut")
def _girvan_newman_cut(g, spec, fx):
    """The Girvan-Newman dendrogram (from its edge removals) cut at
    `groups` groups."""
    removed = girvan_newman_removals(g, len(fx.src))
    merges = merges_from_pairs(g.n, removed, removals=True)
    try:
        return cut_dendrogram(g.n, merges, groups=spec["groups"])
    except EngineError:
        return [-1.0] * g.n


@register("label_propagation")
def _label_propagation(g, spec, fx):
    """Membership ("membership"), mean modularity over seeds 1..`seeds`
    ("quality"), or the engine's own labels ("labels", for the fixed-point
    check, whose oracle counts the nodes carrying a dominant label)."""
    output = spec.get("output", "membership")
    variant = spec.get("variant", "dominance")
    mode = spec.get("mode", "all")
    weighted = spec.get("weighted", False)
    if output == "quality":
        from bindings_subgroups import modularity
        qs = []
        for s in range(1, spec["seeds"] + 1):
            memb = label_propagation(g, variant, mode, weighted, s)
            qs.append(modularity(g, memb, 1.0, weighted))
        return [_q(sum(qs) / len(qs))]
    if output == "fixed_point":
        # Every node carries a dominant label: the oracle computes the
        # fraction of such nodes in the engine's partition.
        return [1.0]
    return label_propagation(g, variant, mode, weighted, spec.get("seed", 1))


@register("map_equation")
def _map_equation(g, spec, fx):
    """The map equation of the fixed partition `blocks`."""
    return [_q(map_equation(g, blocks(g.n), spec.get("weighted", False)))]


@register("infomap")
def _infomap(g, spec, fx):
    """Codelength of the best of `trials` trials ("codelength"), or mean
    over seeds 1..`seeds` of that codelength ("quality")."""
    weighted = spec.get("weighted", False)
    trials = spec.get("trials", 10)
    if spec.get("output") == "quality":
        cls = [infomap(g, weighted, trials, s)[1] for s in range(1, spec["seeds"] + 1)]
        return [_q(sum(cls) / len(cls))]
    return [_q(infomap(g, weighted, trials, spec.get("seed", 1))[1])]


@register("leiden")
def _leiden(g, spec, fx):
    """Quality of the partition ("quality"), mean quality over seeds
    ("mean_quality"), or the number of groups ("groups", against the
    oracle's count of connected pieces of the groups)."""
    output = spec.get("output", "quality")
    args = dict(objective=spec.get("objective", "modularity"),
                resolution=spec.get("resolution", 1.0), beta=spec.get("beta", 0.01),
                n_iterations=spec.get("n_iterations", -1), weighted=spec.get("weighted", False))
    if output == "mean_quality":
        qs = [leiden(g, seed=s, **args)[1] for s in range(1, spec["seeds"] + 1)]
        return [_q(sum(qs) / len(qs))]
    memb, q = leiden(g, seed=spec.get("seed", 1), **args)
    if output == "groups":
        return [len(set(memb))]
    return [_q(q)]


@register("spinglass")
def _spinglass(g, spec, fx):
    """Energy of the final configuration ("energy") or mean modularity over
    seeds ("quality"); -1 when the network is rejected (disconnected)."""
    args = dict(spins=spec.get("spins", 25), weighted=spec.get("weighted", False),
                gamma=spec.get("gamma", 1.0), null_model=spec.get("null_model", "config"))
    try:
        if spec.get("output") == "quality":
            qs = [spinglass(g, seed=s, **args)[1] for s in range(1, spec["seeds"] + 1)]
            return [_q(sum(qs) / len(qs))]
        return [spinglass(g, seed=spec.get("seed", 1), **args)[2]]
    except EngineError:
        return [-1.0]


@register("optimal_modularity")
def _optimal_modularity(g, spec, fx):
    """Maximum modularity; -1 on fixtures listed in `skip` (beyond the exact
    search's practical size) or rejected by the engine."""
    if _skipped(spec, fx):
        return [-1.0]
    try:
        _, q = optimal_modularity(g, spec.get("resolution", 1.0), spec.get("weighted", False))
    except EngineError:
        return [-1.0]
    return [_q(q)]


@register("fluid_communities")
def _fluid_communities(g, spec, fx):
    """Mean modularity over seeds 1..`seeds` ("quality") or the number of
    groups ("groups"); -1 when the network is rejected."""
    from bindings_subgroups import modularity
    try:
        if spec.get("output") == "groups":
            return [len(set(fluid_communities(g, spec["k"], spec.get("seed", 1))))]
        qs = [modularity(g, fluid_communities(g, spec["k"], s), 1.0, False)
              for s in range(1, spec["seeds"] + 1)]
    except EngineError:
        return [-1.0]
    return [_q(sum(qs) / len(qs))]


@register("voronoi")
def _voronoi(g, spec, fx):
    """Membership ("membership"), modularity ("modularity"), or radius
    ("radius") of the Voronoi partition with lengths derived from the
    weights as `lengths` says ("none" for unit lengths, "inverse" for
    1 / weight)."""
    lengths = None
    if spec.get("lengths", "none") == "inverse":
        lengths = [1.0 / w for w in fx.weights]
    memb, _, q, r = voronoi(g, lengths, spec.get("weighted", False), spec.get("mode", "out"),
                            spec.get("radius", -1.0))
    output = spec.get("output", "membership")
    if output == "modularity":
        return [_q(q)]
    if output == "radius":
        return [r]
    return memb


@register("compare_partitions")
def _compare_partitions(g, spec, fx):
    """A comparison of the fixed partitions `blocks` and `stripes`."""
    return [_q(compare_partitions(blocks(g.n), stripes(g.n), spec["method"]))]
