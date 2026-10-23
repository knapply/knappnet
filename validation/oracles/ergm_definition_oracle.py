#!/usr/bin/env python3
"""Golden values for exponential random graph models, computed directly
from the definitions of the statistics.

This oracle shares no code or strategy with the engine. Statistics are
matrix formulas on the sociomatrix (triangles as trace(A^3)/6, shared
partners as entries of A^2, stars as binomial coefficients of degrees),
where the engine counts incrementally from neighbor lists. Change
statistics are differences of whole-network statistics with a dyad set to
1 and to 0. The pseudolikelihood is maximized with scipy's BFGS on the
logistic log-likelihood (the engine uses iteratively reweighted least
squares). Exact likelihoods enumerate every network on the free dyads as
bit patterns with numpy, and the exact MLE maximizes the enumerated
log-likelihood with scipy (the engine uses Gray-code enumeration and
Newton's method).

Node attributes follow the scheme documented in validation/bindings/ergm.py:
x[i] = i mod 5 + 1, a[i] = i mod 3, w[i, j] = ((7 i + 3 j) mod 5) / 2 with
(min, max) ordering for undirected networks. Self-loops are dropped and
parallel edges collapse, as in a sociomatrix.

    python3 validation/oracles/ergm_definition_oracle.py
"""

from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.optimize import minimize
from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import GOLDEN, MANIFEST_DIR, applies, fixture_names, load_fixture  # noqa: E402
from common import write_values  # noqa: E402

ORACLE = "definition"
MAX_FREE_DYADS = 22

TERM_SETS = {
    "dyadic": [("edges", 0, []), ("nodecov", 0, []), ("nodefactor", 0, [1, 2]),
               ("nodematch", 0, []), ("nodematch_diff", 0, [0, 1, 2]), ("absdiff", 1, []),
               ("absdiff", 2, []), ("edgecov", 0, [])],
    "undirected": [("kstar", 2, []), ("kstar", 3, []), ("triangle", 0, []), ("gwesp", 0.5, []),
                   ("gwdsp", 0.3, []), ("gwnsp", 0.3, []), ("gwdegree", 0.5, []),
                   ("degree", 1, []), ("degree", 2, []), ("isolates", 0, [])],
    "directed": [("mutual", 0, []), ("istar", 2, []), ("ostar", 2, []), ("ttriple", 0, []),
                 ("ctriple", 0, []), ("triangle", 0, []), ("gwesp", 0.5, []), ("gwdsp", 0.3, []),
                 ("gwnsp", 0.3, []), ("gwidegree", 0.5, []), ("gwodegree", 0.5, []),
                 ("idegree", 1, []), ("odegree", 2, []), ("isolates", 0, []),
                 ("nodeicov", 0, []), ("nodeocov", 0, []), ("nodeifactor", 0, [1, 2]),
                 ("nodeofactor", 0, [1, 2])],
    "mple_dyadic": [("edges", 0, []), ("nodecov", 0, []), ("nodematch", 0, []),
                    ("absdiff", 1, [])],
    "mple_gwesp": [("edges", 0, []), ("gwesp", 0.5, [])],
    "mple_directed": [("edges", 0, []), ("mutual", 0, []), ("gwesp", 0.5, [])],
    "edges": [("edges", 0, [])],
    "edges_triangle": [("edges", 0, []), ("triangle", 0, [])],
    "edges_gwesp": [("edges", 0, []), ("gwesp", 0.5, [])],
    "edges_kstar2": [("edges", 0, []), ("kstar", 2, [])],
}


def manifest() -> dict:
    return yaml.safe_load((MANIFEST_DIR / "ergm.yaml").read_text())


def attributes(n: int, directed: bool):
    x = np.array([i % 5 + 1 for i in range(n)], dtype=float)
    a = np.array([i % 3 for i in range(n)])
    w = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            lo, hi = (i, j) if directed else (min(i, j), max(i, j))
            w[i, j] = ((7 * lo + 3 * hi) % 5) / 2.0
    return x, a, w


def sociomatrix(fx) -> np.ndarray:
    n = fx.n
    A = np.zeros((n, n))
    for i, j in zip(fx.src, fx.tgt):
        if i == j:
            continue
        A[i, j] = 1
        if not fx.directed:
            A[j, i] = 1
    return A


def gw(alpha: float, k):
    k = np.asarray(k, dtype=float)
    return np.where(k > 0, math.exp(alpha) * (1 - (1 - math.exp(-alpha)) ** k), 0.0)


def choose(d, k):
    d = np.asarray(d, dtype=float)
    out = np.ones_like(d)
    for i in range(k):
        out = out * (d - i) / (i + 1)
    return np.where(d >= k, out, 0.0)


def statistics(A: np.ndarray, directed: bool, terms, cov) -> list[float]:
    """Every statistic of the term list, by matrix formulas."""
    x, a, w = cov
    n = A.shape[0]
    od = A.sum(axis=1)
    idg = A.sum(axis=0)
    deg = od + idg if directed else od  # total degree (undirected: degree)
    half = 1.0 if directed else 0.5    # undirected sums count each tie twice
    A2 = A @ A
    if directed:
        pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
    else:
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    sp = np.array([A2[i, j] for i, j in pairs])
    tie = np.array([A[i, j] for i, j in pairs])
    out = []
    for name, param, levels in terms:
        if name == "edges":
            out.append(half * A.sum())
        elif name == "mutual":
            out.append((A * A.T).sum() / 2)
        elif name == "nodecov":
            out.append((deg * x).sum())
        elif name == "nodeocov":
            out.append((od * x).sum())
        elif name == "nodeicov":
            out.append((idg * x).sum())
        elif name in ("nodefactor", "nodeofactor", "nodeifactor"):
            d = {"nodefactor": deg, "nodeofactor": od, "nodeifactor": idg}[name]
            out.extend((d * (a == l)).sum() for l in levels)
        elif name == "nodematch":
            M = (a[:, None] == a[None, :]).astype(float)
            out.append(half * (A * M).sum())
        elif name == "nodematch_diff":
            for l in levels:
                M = ((a[:, None] == l) & (a[None, :] == l)).astype(float)
                out.append(half * (A * M).sum())
        elif name == "absdiff":
            D = np.abs(x[:, None] - x[None, :]) ** param
            out.append(half * (A * D).sum())
        elif name == "edgecov":
            out.append(half * (A * w).sum())
        elif name == "kstar":
            out.append(choose(deg, int(param)).sum())
        elif name == "istar":
            out.append(choose(idg, int(param)).sum())
        elif name == "ostar":
            out.append(choose(od, int(param)).sum())
        elif name == "triangle":
            A3 = np.trace(A2 @ A)
            if directed:
                out.append((A * A2).sum() + A3 / 3)  # transitive plus cyclic
            else:
                out.append(A3 / 6)
        elif name == "ttriple":
            out.append((A * A2).sum())
        elif name == "ctriple":
            out.append(np.trace(A2 @ A) / 3)
        elif name == "gwesp":
            out.append((tie * gw(param, sp)).sum())
        elif name == "gwdsp":
            out.append(gw(param, sp).sum())
        elif name == "gwnsp":
            out.append(((1 - tie) * gw(param, sp)).sum())
        elif name == "gwdegree":
            out.append(gw(param, deg).sum())
        elif name == "gwidegree":
            out.append(gw(param, idg).sum())
        elif name == "gwodegree":
            out.append(gw(param, od).sum())
        elif name == "degree":
            out.append(float((deg == param).sum()))
        elif name == "idegree":
            out.append(float((idg == param).sum()))
        elif name == "odegree":
            out.append(float((od == param).sum()))
        elif name == "isolates":
            out.append(float((deg == 0).sum()))
        else:
            raise ValueError(name)
    return [float(v) for v in out]


def dyads(n: int, directed: bool):
    if directed:
        return [(i, j) for i in range(n) for j in range(n) if i != j]
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def set_dyad(A, i, j, v, directed):
    B = A.copy()
    B[i, j] = v
    if not directed:
        B[j, i] = v
    return B


def change_matrix(A, directed, terms, cov):
    """Change statistics by full recomputation: g(y + ij) - g(y - ij)."""
    rows, y = [], []
    for i, j in dyads(A.shape[0], directed):
        with_ = statistics(set_dyad(A, i, j, 1, directed), directed, terms, cov)
        without = statistics(set_dyad(A, i, j, 0, directed), directed, terms, cov)
        rows.append([u - v for u, v in zip(with_, without)])
        y.append(A[i, j])
    return np.array(rows), np.array(y)


def logistic_mle(X, y):
    """Maximum likelihood logistic regression without intercept, by BFGS."""
    def nll(b):
        eta = X @ b
        return float(np.sum(np.logaddexp(0, eta) - y * eta))

    def grad(b):
        p = 1 / (1 + np.exp(-(X @ b)))
        return X.T @ (p - y)

    res = minimize(nll, np.zeros(X.shape[1]), jac=grad, method="BFGS",
                   options={"gtol": 1e-10, "maxiter": 5000})
    # Observed information for standard errors.
    p = 1 / (1 + np.exp(-(X @ res.x)))
    info = (X * (p * (1 - p))[:, None]).T @ X
    return res.x, np.sqrt(np.diag(np.linalg.inv(info))), res


def separated(X, y, coef):
    """A fit whose coefficients keep growing has no finite maximizer."""
    return bool(np.any(np.abs(coef) > 15))


# ---- Exact enumeration -------------------------------------------------------

def enumerate_stats(A, directed, terms, cov):
    """Statistics of every network on the free dyads, one row per network
    (bit d of the row index is dyad d), by vectorized definitions."""
    n = A.shape[0]
    x, a, w = cov
    ds = dyads(n, directed)
    D = len(ds)
    count = 1 << D
    idx = np.arange(count, dtype=np.int64)
    Y = {d: ((idx >> k) & 1).astype(np.float64) for k, d in enumerate(ds)}

    def y(i, j):
        if directed:
            return Y[(i, j)]
        return Y[(min(i, j), max(i, j))]

    def two_paths(i, j):
        s = np.zeros(count)
        for k in range(n):
            if k != i and k != j:
                s += y(i, k) * y(k, j)
        return s

    deg = [sum(y(i, j) + (y(j, i) if directed else 0) for j in range(n) if j != i) for i in range(n)]
    out = []
    for name, param, levels in terms:
        if name == "edges":
            out.append(sum(Y.values()))
        elif name == "triangle":
            if directed:
                tt = sum(y(i, j) * two_paths(i, j) for i, j in ds)
                cy = sum(y(i, j) * two_paths(j, i) for i, j in ds) / 3
                out.append(tt + cy)
            else:
                t = np.zeros(count)
                for i, j, k in itertools.combinations(range(n), 3):
                    t += y(i, j) * y(j, k) * y(i, k)
                out.append(t)
        elif name == "gwesp":
            s = np.zeros(count)
            for i, j in ds:
                s += y(i, j) * gw(param, two_paths(i, j))
            out.append(s)
        elif name == "kstar":
            out.append(sum(choose(d, int(param)) for d in deg))
        else:
            raise ValueError(f"exact enumeration does not cover {name}")
    return np.column_stack(out)


def exact_loglik(G, g_obs, theta):
    lin = G @ theta
    return float(g_obs @ theta - logsumexp(lin))


def exact_mean(G, theta):
    lin = G @ theta
    w = np.exp(lin - logsumexp(lin))
    return w @ G


def exact_mle(G, g_obs):
    def nll(t):
        return -exact_loglik(G, g_obs, t)

    def grad(t):
        return exact_mean(G, t) - g_obs

    res = minimize(nll, np.zeros(G.shape[1]), jac=grad, method="BFGS",
                   options={"gtol": 1e-11, "maxiter": 5000})
    return res.x, -res.fun


def main() -> None:
    spec_all = manifest()
    written = 0
    for name in fixture_names():
        fx = load_fixture(name)
        A = sociomatrix(fx)
        cov = attributes(fx.n, fx.directed)
        enumerations = {}
        for variant, spec in spec_all.items():
            if ORACLE not in spec["oracles"] or not applies(spec, fx):
                continue
            if name in spec["oracles"][ORACLE].get("undefined_on", []):
                continue
            eng = spec["engine"]
            terms = TERM_SETS[eng["terms"]]
            fn = eng["fn"]
            if fn == "ergm_stats":
                values = statistics(A, fx.directed, terms, cov)
            elif fn == "ergm_mple":
                X, y = change_matrix(A, fx.directed, terms, cov)
                coef, se, _ = logistic_mle(X, y)
                if separated(X, y, coef):
                    print(f"   {variant} on {name}: no finite pseudolikelihood maximizer; skipped")
                    continue
                values = list(se) if eng.get("output") == "se" else list(coef)
            elif fn in ("ergm_exact_loglik", "ergm_exact_mle", "ergm_mcmle"):
                # The Monte Carlo MLE is checked against the exact MLE where
                # enumeration is possible.
                if fn == "ergm_mcmle" and name not in eng.get("fixtures", []):
                    continue
                if len(dyads(fx.n, fx.directed)) > MAX_FREE_DYADS:
                    continue
                key = eng["terms"]
                if key not in enumerations:
                    enumerations[key] = enumerate_stats(A, fx.directed, terms, cov)
                G = enumerations[key]
                g_obs = np.array(statistics(A, fx.directed, terms, cov))
                if fn == "ergm_exact_loglik":
                    theta = np.array(eng["theta"], dtype=float)
                    if eng.get("output", "loglik") == "loglik":
                        values = [exact_loglik(G, g_obs, theta)]
                    else:
                        values = list(exact_mean(G, theta))
                else:
                    theta, ll = exact_mle(G, g_obs)
                    values = [ll] if eng.get("output") == "loglik" else list(theta)
            else:
                continue
            write_values(ORACLE, name, variant, values)
            written += 1
    (GOLDEN / ORACLE).mkdir(parents=True, exist_ok=True)
    (GOLDEN / ORACLE / "VERSION").write_text(
        "definition (direct evaluation, exact arithmetic)\n", newline="")
    print(f"Wrote {written} definition golden files (ergm).")


if __name__ == "__main__":
    main()
