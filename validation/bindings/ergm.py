"""Engine runners for exponential random graph models.

The fixtures carry no node attributes, so every ERGM oracle and this runner
derive the same deterministic attributes from the 0-based node index i:

    x[i] = i mod 5 + 1            numeric (nodecov, absdiff)
    a[i] = i mod 3                categorical with levels 0, 1, 2
    w[i, j] = ((7 i + 3 j) mod 5) / 2   dyadic covariate; for undirected
             networks the (min, max) ordering is used, so w is symmetric

Self-loops are dropped (ERGMs are defined on loopless networks) and
parallel edges collapse to one tie, which is the sociomatrix reading that
the R oracle uses to build `network` objects.

Term sets named in the manifest (`terms:`) are defined in TERM_SETS.
Results are flat lists: statistics in term order, then (for estimates)
coefficients, and for the exact variants the log-likelihood, expectations,
or coefficients as the manifest entry says.
"""

from __future__ import annotations

import ctypes
import math

from engine import ERR_LEN, EngineError, load, register

CODES = {
    "edges": 0, "mutual": 1, "nodecov": 2, "nodeicov": 3, "nodeocov": 4,
    "nodefactor": 5, "nodeifactor": 6, "nodeofactor": 7, "nodematch": 8,
    "nodematch_diff": 9, "absdiff": 10, "edgecov": 11, "kstar": 12, "istar": 13,
    "ostar": 14, "triangle": 15, "ttriple": 16, "ctriple": 17, "gwesp": 18,
    "gwdsp": 19, "gwnsp": 20, "gwdegree": 21, "gwidegree": 22, "gwodegree": 23,
    "degree": 24, "idegree": 25, "odegree": 26, "isolates": 27,
}

# (name, param, attribute index, levels). Attribute index 0 refers to x
# (numeric terms), a (categorical terms), or w (edgecov).
TERM_SETS = {
    "dyadic": [("edges", 0, -1, []), ("nodecov", 0, 0, []),
               ("nodefactor", 0, 0, [1, 2]), ("nodematch", 0, 0, []),
               ("nodematch_diff", 0, 0, [0, 1, 2]), ("absdiff", 1, 0, []),
               ("absdiff", 2, 0, []), ("edgecov", 0, 0, [])],
    "undirected": [("kstar", 2, -1, []), ("kstar", 3, -1, []), ("triangle", 0, -1, []),
                   ("gwesp", 0.5, -1, []), ("gwdsp", 0.3, -1, []), ("gwnsp", 0.3, -1, []),
                   ("gwdegree", 0.5, -1, []), ("degree", 1, -1, []), ("degree", 2, -1, []),
                   ("isolates", 0, -1, [])],
    "directed": [("mutual", 0, -1, []), ("istar", 2, -1, []), ("ostar", 2, -1, []),
                 ("ttriple", 0, -1, []), ("ctriple", 0, -1, []), ("triangle", 0, -1, []),
                 ("gwesp", 0.5, -1, []), ("gwdsp", 0.3, -1, []), ("gwnsp", 0.3, -1, []),
                 ("gwidegree", 0.5, -1, []), ("gwodegree", 0.5, -1, []),
                 ("idegree", 1, -1, []), ("odegree", 2, -1, []), ("isolates", 0, -1, []),
                 ("nodeicov", 0, 0, []), ("nodeocov", 0, 0, []),
                 ("nodeifactor", 0, 0, [1, 2]), ("nodeofactor", 0, 0, [1, 2])],
    "mple_dyadic": [("edges", 0, -1, []), ("nodecov", 0, 0, []), ("nodematch", 0, 0, []),
                    ("absdiff", 1, 0, [])],
    "mple_gwesp": [("edges", 0, -1, []), ("gwesp", 0.5, -1, [])],
    "mple_directed": [("edges", 0, -1, []), ("mutual", 0, -1, []), ("gwesp", 0.5, -1, [])],
    "edges": [("edges", 0, -1, [])],
    "edges_triangle": [("edges", 0, -1, []), ("triangle", 0, -1, [])],
    "edges_gwesp": [("edges", 0, -1, []), ("gwesp", 0.5, -1, [])],
    "edges_kstar2": [("edges", 0, -1, []), ("kstar", 2, -1, [])],
    # Temporal models: dyadic covariate 1 is the memory (stability) matrix
    # 2 y_prev - 1 and covariate 2 the transpose of y_prev (delayed
    # reciprocity).
    "tergm_undirected": [("edges", 0, -1, []), ("nodecov", 0, 0, []), ("edgecov", 0, 1, [])],
    "tergm_directed": [("edges", 0, -1, []), ("mutual", 0, -1, []), ("nodecov", 0, 0, []),
                       ("edgecov", 0, 1, []), ("edgecov", 0, 2, [])],
}


def second_wave(fx):
    """The fixture's ties after a deterministic rewiring: dyad (i, j), in
    canonical order, is toggled when (7 i + 3 j) mod 11 == 0. Returns the
    tie sets of both waves (loops dropped, parallel edges collapsed)."""
    first = set()
    for a, b in zip(fx.src, fx.tgt):
        if a == b:
            continue
        first.add((a, b) if fx.directed else (min(a, b), max(a, b)))
    second = set(first)
    for i in range(fx.n):
        for j in range(fx.n):
            if i == j or (not fx.directed and j < i):
                continue
            if (7 * i + 3 * j) % 11 == 0:
                second ^= {(i, j)}
    return sorted(first), sorted(second)


def memory_matrices(n, directed, first):
    """Column-major 2 y_prev - 1 and t(y_prev) for the tie set `first`."""
    prev = [[0.0] * n for _ in range(n)]
    for a, b in first:
        prev[a][b] = 1.0
        if not directed:
            prev[b][a] = 1.0
    stability = [0.0] * (n * n)
    reciprocal = [0.0] * (n * n)
    for i in range(n):
        for j in range(n):
            stability[j * n + i] = 2 * prev[i][j] - 1
            reciprocal[j * n + i] = prev[j][i]
    return stability, reciprocal

_declared = False


def _declare() -> ctypes.CDLL:
    global _declared
    lib = load()
    if _declared:
        return lib
    i32, i64, f64, u64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_double, ctypes.c_uint64
    i32p, i64p, f64p = ctypes.POINTER(i32), ctypes.POINTER(i64), ctypes.POINTER(f64)
    h = ctypes.c_void_p
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    lib.net_ergm_new.argtypes = [i64, i64p, i64p, i64, i32, i32p, f64p, i32p, i32p, i32p, i32p,
                                 i32, f64p, i32, i32p, i32, f64p, i32, i32p,
                                 ctypes.POINTER(h), cp, sz]
    lib.net_ergm_free.argtypes = [h]
    lib.net_ergm_width.argtypes = [h]
    lib.net_ergm_width.restype = i64
    lib.net_ergm_dyads.argtypes = [h]
    lib.net_ergm_dyads.restype = i64
    lib.net_ergm_stats.argtypes = [h, f64p, cp, sz]
    lib.net_ergm_mple.argtypes = [h, i32, f64, f64p, f64p, f64p, cp, sz]
    lib.net_ergm_exact.argtypes = [h, f64p, f64p, f64p, f64p, cp, sz]
    lib.net_ergm_exact_mle.argtypes = [h, f64p, i32, f64p, f64p, f64p, cp, sz]
    lib.net_ergm_mcmle.argtypes = [h, f64p, i64, i64, i64, i32, i32, u64, i32, i32, f64, f64,
                                   f64p, f64p, f64p, f64p, f64p, cp, sz]
    lib.net_ergm_loglik.argtypes = [h, f64p, i32, i64, i64, i64, i32, i32, u64, i32, f64p, cp, sz]
    lib.net_ergm_simulate.argtypes = [h, f64p, i64, i64, i64, i32, i32, u64, i32, i32, i32, f64p,
                                      f64p, ctypes.POINTER(h), cp, sz]
    lib.net_ergm_changestats.argtypes = [h, f64p, f64p, i64p, i64p, cp, sz]
    lib.net_tergm_mple.argtypes = [f64p, f64p, i32p, i64, i64, i32, i64, u64, i32, i32, f64,
                                   f64p, f64p, f64p, f64p, cp, sz]
    for fn in ("net_ergm_new", "net_ergm_stats", "net_ergm_mple", "net_ergm_exact",
               "net_ergm_exact_mle", "net_ergm_mcmle", "net_ergm_loglik", "net_ergm_simulate",
               "net_ergm_changestats", "net_tergm_mple"):
        getattr(lib, fn).restype = i32
    _declared = True
    return lib


def _check(status: int, err) -> None:
    if status:
        raise EngineError(err.value.decode())


def attributes(n: int, directed: bool):
    x = [float(i % 5 + 1) for i in range(n)]
    a = [i % 3 for i in range(n)]
    w = [0.0] * (n * n)
    for i in range(n):
        for j in range(n):
            lo, hi = (i, j) if directed else (min(i, j), max(i, j))
            w[j * n + i] = ((7 * lo + 3 * hi) % 5) / 2.0  # column-major
    return x, a, w


class _Ergm:
    """A model handle over a fixture with the shared attribute scheme."""

    def __init__(self, fx, terms: str, ties=None, extra_dyadic=()):
        self.lib = _declare()
        n = fx.n
        if ties is None:
            pairs = [(a, b) for a, b in zip(fx.src, fx.tgt) if a != b]
        else:
            pairs = list(ties)
        m = len(pairs)
        src = (ctypes.c_int64 * max(m, 1))(*[a for a, _ in pairs])
        tgt = (ctypes.c_int64 * max(m, 1))(*[b for _, b in pairs])
        spec = TERM_SETS[terms]
        k = len(spec)
        codes = (ctypes.c_int32 * k)(*[CODES[t[0]] for t in spec])
        params = (ctypes.c_double * k)(*[float(t[1]) for t in spec])
        attrs = (ctypes.c_int32 * k)(*[t[2] for t in spec])
        levels_flat = [l for t in spec for l in t[3]]
        offs, acc = [], 0
        for t in spec:
            offs.append(acc)
            acc += len(t[3])
        off = (ctypes.c_int32 * k)(*offs)
        lens = (ctypes.c_int32 * k)(*[len(t[3]) for t in spec])
        levels = (ctypes.c_int32 * max(len(levels_flat), 1))(*levels_flat)
        x, a, w = attributes(n, fx.directed)
        num = (ctypes.c_double * max(n, 1))(*x)
        cat = (ctypes.c_int32 * max(n, 1))(*a)
        dyadic = [w] + list(extra_dyadic)
        flat = [v for mat in dyadic for v in mat]
        dy = (ctypes.c_double * max(len(flat), 1))(*flat)
        self.handle = ctypes.c_void_p()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_ergm_new(n, src, tgt, m, int(fx.directed), codes, params, attrs, off,
                                     lens, levels, k, num, 1, cat, 1, dy, len(dyadic), None,
                                     ctypes.byref(self.handle), err, ERR_LEN), err)
        self.p = self.lib.net_ergm_width(self.handle)
        self.n = n

    def __del__(self):
        if getattr(self, "handle", None):
            self.lib.net_ergm_free(self.handle)
            self.handle = None

    def stats(self):
        out = (ctypes.c_double * self.p)()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_ergm_stats(self.handle, out, err, ERR_LEN), err)
        return list(out)

    def changestats(self):
        d = self.lib.net_ergm_dyads(self.handle)
        p = self.p
        x = (ctypes.c_double * max(d * p, 1))()
        y = (ctypes.c_double * max(d, 1))()
        tail = (ctypes.c_int64 * max(d, 1))()
        head = (ctypes.c_int64 * max(d, 1))()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_ergm_changestats(self.handle, x, y, tail, head, err, ERR_LEN), err)
        return list(x)[:d * p], list(y)[:d], d

    def mple(self):
        p = self.p
        coef = (ctypes.c_double * p)()
        cov = (ctypes.c_double * (p * p))()
        fit = (ctypes.c_double * 6)()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_ergm_mple(self.handle, 0, 0.0, coef, cov, fit, err, ERR_LEN), err)
        se = [math.sqrt(cov[k * p + k]) for k in range(p)]
        return list(coef), se, list(fit)

    def exact(self, theta):
        p = self.p
        th = (ctypes.c_double * p)(*theta)
        ll = ctypes.c_double()
        mean = (ctypes.c_double * p)()
        cov = (ctypes.c_double * (p * p))()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_ergm_exact(self.handle, th, ctypes.byref(ll), mean, cov, err, ERR_LEN), err)
        return ll.value, list(mean)

    def exact_mle(self):
        p = self.p
        start = (ctypes.c_double * p)()
        theta = (ctypes.c_double * p)()
        cov = (ctypes.c_double * (p * p))()
        info = (ctypes.c_double * 2)()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_ergm_exact_mle(self.handle, start, 200, theta, cov, info, err, ERR_LEN), err)
        return list(theta), info[0]

    def mcmle(self, burnin, interval, nsim, chains, seed, critical):
        p = self.p
        # Start from the pseudolikelihood estimate when it exists; when the
        # pseudolikelihood is separated (as for edges + triangle on two
        # disjoint triangles, where the exact MLE still exists) start from
        # zero, which is what the R interface does as well.
        try:
            coef0, _, fit = self.mple()
            if fit[4] != 1.0 or any(abs(c) > 10 for c in coef0):
                coef0 = [0.0] * p
        except EngineError:
            coef0 = [0.0] * p
        theta0 = (ctypes.c_double * p)(*coef0)
        theta = (ctypes.c_double * p)()
        cov = (ctypes.c_double * (p * p))()
        mc = (ctypes.c_double * (p * p))()
        info = (ctypes.c_double * 5)()
        stats = (ctypes.c_double * (chains * nsim * p))()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_ergm_mcmle(self.handle, theta0, burnin, interval, nsim, chains, 0, seed,
                                       0, 60, 0.05, critical, theta, cov, mc, info, stats, err,
                                       ERR_LEN), err)
        return list(theta), info[1] == 1.0

    def loglik(self, theta, steps, burnin, interval, nsim, chains, seed):
        p = self.p
        th = (ctypes.c_double * p)(*theta)
        out = (ctypes.c_double * 2)()
        err = ctypes.create_string_buffer(ERR_LEN)
        _check(self.lib.net_ergm_loglik(self.handle, th, steps, burnin, interval, nsim, chains, 0,
                                        seed, 0, out, err, ERR_LEN), err)
        dyads = self.lib.net_ergm_dyads(self.handle)
        return out[0] - dyads * math.log(2.0), out[1]


@register("ergm_stats")
def _ergm_stats(g, spec, fx):
    return _Ergm(fx, spec["terms"]).stats()


# Fixtures on which an estimate does not exist (separation, no dyads, too
# many dyads to enumerate) make the engine report an error; the manifest
# declares those fixtures undefined, so the runner returns nothing there.
def _or_nothing(fn):
    try:
        return fn()
    except EngineError:
        return []


@register("ergm_mple")
def _ergm_mple(g, spec, fx):
    def run():
        coef, se, _ = _Ergm(fx, spec["terms"]).mple()
        return se if spec.get("output") == "se" else coef
    return _or_nothing(run)


@register("ergm_exact_loglik")
def _ergm_exact_loglik(g, spec, fx):
    def run():
        ll, mean = _Ergm(fx, spec["terms"]).exact(spec["theta"])
        return [ll] if spec.get("output", "loglik") == "loglik" else mean
    return _or_nothing(run)


@register("ergm_exact_mle")
def _ergm_exact_mle(g, spec, fx):
    def run():
        theta, ll = _Ergm(fx, spec["terms"]).exact_mle()
        return [ll] if spec.get("output") == "loglik" else theta
    return _or_nothing(run)


# Chi-squared critical values for the convergence test, by degrees of
# freedom (0.95 quantiles).
CHISQ95 = {1: 3.841459, 2: 5.991465, 3: 7.814728, 4: 9.487729, 5: 11.07050}


@register("ergm_mcmle")
def _ergm_mcmle(g, spec, fx):
    if fx.name not in spec.get("fixtures", []):
        return []
    m = _Ergm(fx, spec["terms"])
    theta, converged = m.mcmle(spec.get("burnin", 20000), spec.get("interval", 200),
                               spec.get("nsim", 10000), spec.get("chains", 2),
                               spec.get("seed", 1), CHISQ95[m.p])
    if not converged:
        raise EngineError("MCMC MLE did not converge")
    if spec.get("output") == "loglik":
        ll, _ = m.loglik(theta, 16, spec.get("burnin", 20000), spec.get("interval", 200),
                         spec.get("nsim", 10000), spec.get("chains", 2), spec.get("seed", 1) + 1)
        return [ll]
    return theta


@register("tergm_mple")
def _tergm_mple(g, spec, fx):
    """Pooled pseudolikelihood over the one transition from the fixture to
    its deterministically rewired second wave."""
    if fx.name not in spec.get("fixtures", []):
        return []
    first, second = second_wave(fx)
    stability, reciprocal = memory_matrices(fx.n, fx.directed, first)
    m = _Ergm(fx, spec["terms"], ties=second, extra_dyadic=[stability, reciprocal])
    x, y, d = m.changestats()
    p = m.p
    lib = m.lib
    xa = (ctypes.c_double * max(len(x), 1))(*x)
    ya = (ctypes.c_double * max(d, 1))(*y)
    period = (ctypes.c_int32 * max(d, 1))()
    coef = (ctypes.c_double * p)()
    cov = (ctypes.c_double * (p * p))()
    fit = (ctypes.c_double * 5)()
    boot = (ctypes.c_double * 1)()
    err = ctypes.create_string_buffer(ERR_LEN)
    _check(lib.net_tergm_mple(xa, ya, period, d, p, 1, 0, 1, 0, 0, 0.0, coef, cov, fit, boot,
                              err, ERR_LEN), err)
    return list(coef)
