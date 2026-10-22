#!/usr/bin/env python3
"""Validate the quadratic assignment procedure (QAP) against independent
implementations.

QAP works on n-by-n matrices rather than the graph fixtures, so it is
checked here instead of through the manifest. Three kinds of checks:

1. Least squares. Coefficients, standard errors, t-statistics, R-squared,
   adjusted R-squared, and F from the engine match numpy's least squares on
   the vectorized dyads (with listwise deletion of missing dyads), for
   directed and undirected data, with and without the diagonal and the
   intercept.

2. Exact p-values. A definition oracle in this file (numpy only, sharing no
   code with the engine) enumerates all n! relabelings for small n and
   computes the exact permutation distribution of the QAP correlation and
   of the MRQAP statistics, for Y-permutation and for double
   semi-partialling (DSP). The engine's exhaustive mode must reproduce the
   oracle's counts exactly, and its Monte Carlo p-values with many draws
   must lie within a binomial tolerance of the exact ones.

3. Reproducibility. Results are identical for every thread count.

The oracle follows the published procedures: the matrix is relabeled as
M[perm, perm] (rows and columns permuted together); Y-permutation
(Krackhardt 1988) permutes the response; DSP (Dekker, Krackhardt, and
Snijders 2007) replaces predictor k by its residuals on the other
predictors, permutes that residual matrix, and refits with y and the other
predictors, using the t-statistic of the permuted column. Statistics are
counted as tied with the observed one when they differ by at most 1e-10
times max(1, |observed|), the engine's documented rule.

    python3 validation/qap_check.py

Exits non-zero on any failure.
"""

from __future__ import annotations

import ctypes
import itertools
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from engine import ERR_LEN, EngineError, load  # noqa: E402

TIE_EPS = 1e-10
# Monte Carlo p-values must be within this many binomial standard errors of
# the exact value (plus a small absolute floor for p near 0 or 1). With
# about 60 comparisons, a 4.5-sigma band has a false-alarm chance near 1e-4.
Z = 4.5

failures: list[str] = []


def check(ok: bool, what: str) -> None:
    if not ok:
        failures.append(what)
        print(f"FAIL  {what}")


# ---- Engine binding -------------------------------------------------------

def _lib():
    lib = load()
    f64p = ctypes.POINTER(ctypes.c_double)
    i32, i64, u64 = ctypes.c_int32, ctypes.c_int64, ctypes.c_uint64
    lib.net_qap_cor.argtypes = [f64p, f64p, i64, i32, i32, i32, i64, u64, i32,
                                f64p, f64p, ctypes.c_char_p, ctypes.c_size_t]
    lib.net_qap_cor.restype = i32
    lib.net_mrqap.argtypes = [f64p, f64p, i64, i64, i32, i32, i32, i32, i32,
                              ctypes.c_double, i32, i64, u64, i32, f64p, f64p, f64p,
                              ctypes.c_char_p, ctypes.c_size_t]
    lib.net_mrqap.restype = i32
    return lib


LIB = None


def _colmajor(m: np.ndarray):
    flat = np.asarray(m, dtype=float).flatten(order="F")
    return (ctypes.c_double * max(len(flat), 1))(*flat)


def engine_qap_cor(y, x, mode="digraph", diag=False, reps=1000, seed=1,
                   exhaustive=False, threads=0, keep=False):
    n = y.shape[0]
    draws = math.factorial(n) if exhaustive else reps
    summary = (ctypes.c_double * 7)()
    dist = (ctypes.c_double * max(draws, 1))() if keep else None
    err = ctypes.create_string_buffer(ERR_LEN)
    status = LIB.net_qap_cor(_colmajor(y), _colmajor(x), n, int(mode == "graph"), int(diag),
                             int(exhaustive), reps, seed, threads, summary, dist, err, ERR_LEN)
    if status:
        raise EngineError(err.value.decode())
    s = list(summary)
    out = {"r": s[0], "dyads": int(s[1]), "draws": int(s[2]),
           "le": int(s[3]), "ge": int(s[4]), "abs": int(s[5]), "valid": int(s[6])}
    if keep:
        out["dist"] = list(dist)[:draws]
    return out


def engine_mrqap(y, xs, mode="digraph", diag=False, intercept=True, method="dsp",
                 statistic="t", reps=1000, seed=1, exhaustive=False, threads=0,
                 keep=False):
    n = y.shape[0]
    k = len(xs)
    p = k + int(intercept)
    draws = math.factorial(n) if exhaustive else reps
    xflat = np.concatenate([np.asarray(m, float).flatten(order="F") for m in xs]) if k else []
    xbuf = (ctypes.c_double * max(len(xflat), 1))(*xflat)
    coef = (ctypes.c_double * (7 * p))()
    fit = (ctypes.c_double * 10)()
    dist = (ctypes.c_double * max(draws * p, 1))() if keep else None
    err = ctypes.create_string_buffer(ERR_LEN)
    status = LIB.net_mrqap(_colmajor(y), xbuf, n, k, int(intercept), int(mode == "graph"),
                           int(diag), {"dsp": 0, "y": 1}[method], {"t": 0, "beta": 1}[statistic],
                           1e-7, int(exhaustive), reps, seed, threads, coef, fit, dist,
                           err, ERR_LEN)
    if status:
        raise EngineError(err.value.decode())
    table = np.array(list(coef)).reshape((7, p)).T
    f = list(fit)
    out = {"coef": table[:, 0], "se": table[:, 1], "t": table[:, 2],
           "le": table[:, 3].astype(int), "ge": table[:, 4].astype(int),
           "abs": table[:, 5].astype(int), "valid": table[:, 6].astype(int),
           "dyads": int(f[0]), "df": int(f[1]), "rss": f[2], "mss": f[3], "r2": f[4],
           "adj_r2": f[5], "sigma": f[6], "F": f[7], "method": int(f[8]),
           "draws": int(f[9])}
    if keep:
        out["dist"] = np.array(list(dist)[:draws * p]).reshape((p, draws)).T
    return out


# ---- Definition oracle ----------------------------------------------------

def dyad_mask(n: int, mode: str, diag: bool) -> np.ndarray:
    r, c = np.indices((n, n))
    if mode == "graph":
        mask = r > c
    else:
        mask = r != c
    if diag:
        mask = mask | (r == c)
    return mask


def symmetrize_lower(m: np.ndarray) -> np.ndarray:
    """Reads an undirected matrix from its lower triangle."""
    low = np.tril(m)
    return low + np.tril(m, -1).T


def vec(m: np.ndarray, mask: np.ndarray) -> np.ndarray:
    # Column-major order of the selected cells, as R's m[mask].
    return m.T[mask.T]


def oracle_ols(cols: list[np.ndarray], y: np.ndarray):
    X = np.column_stack(cols)
    ok = ~np.isnan(y) & ~np.isnan(X).any(axis=1)
    X, yy = X[ok], y[ok]
    N, p = X.shape
    if N <= p or np.linalg.matrix_rank(X) < p:
        return None
    beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
    fitted = X @ beta
    resid = yy - fitted
    rss = float(resid @ resid)
    sigma2 = rss / (N - p)
    se = np.sqrt(np.diag(np.linalg.inv(X.T @ X)) * sigma2)
    full_resid = np.full(len(y), np.nan)
    full_resid[ok] = resid
    return {"coef": beta, "se": se, "t": beta / se, "N": N, "rss": rss,
            "fitted": fitted, "resid": full_resid}


def oracle_summary(fit, intercept: bool, p: int):
    f = fit["fitted"]
    mss = float(((f - f.mean()) ** 2).sum()) if intercept else float((f ** 2).sum())
    rss, N = fit["rss"], fit["N"]
    r2 = mss / (mss + rss)
    df_int = 1 if intercept else 0
    adj = 1 - (1 - r2) * ((N - df_int) / (N - p))
    F = (mss / (p - df_int)) / (rss / (N - p)) if p > df_int else float("nan")
    return {"mss": mss, "r2": r2, "adj_r2": adj, "F": F}


def prep(m, mode):
    m = np.asarray(m, float)
    return symmetrize_lower(m) if mode == "graph" else m


def design(y, xs, mode, diag, intercept):
    n = y.shape[0]
    mask = dyad_mask(n, mode, diag)
    mats = ([np.ones((n, n))] if intercept else []) + [prep(x, mode) for x in xs]
    return mask, prep(y, mode), mats


def count(stats, observed):
    eps = TIE_EPS * max(1.0, abs(observed))
    s = np.asarray([v for v in stats if not np.isnan(v)])
    return {"le": int((s <= observed + eps).sum()), "ge": int((s >= observed - eps).sum()),
            "abs": int((np.abs(s) >= abs(observed) - eps).sum()), "valid": len(s)}


def oracle_cor(y, x, mode="digraph", diag=False):
    n = y.shape[0]
    mask = dyad_mask(n, mode, diag)
    yy, xx = prep(y, mode), prep(x, mode)
    yv = vec(yy, mask)

    def r(a, b):
        ok = ~np.isnan(a) & ~np.isnan(b)
        return float(np.corrcoef(a[ok], b[ok])[0, 1])

    obs = r(yv, vec(xx, mask))
    stats = []
    for perm in itertools.permutations(range(n)):
        p = np.array(perm)
        stats.append(r(yv, vec(xx[np.ix_(p, p)], mask)))
    return obs, count(stats, obs)


def oracle_mrqap_exact(y, xs, mode="digraph", diag=False, intercept=True, method="dsp",
                       statistic="t"):
    n = y.shape[0]
    mask, yy, mats = design(y, xs, mode, diag, intercept)
    yv = vec(yy, mask)
    cols = [vec(m, mask) for m in mats]
    p = len(cols)
    base = oracle_ols(cols, yv)
    key = "t" if statistic == "t" else "coef"
    observed = base[key]
    if p == 1:
        method = "y"
    perms = [np.array(q) for q in itertools.permutations(range(n))]
    stats = np.full((len(perms), p), np.nan)
    if method == "y":
        for i, q in enumerate(perms):
            f = oracle_ols(cols, vec(yy[np.ix_(q, q)], mask))
            if f is not None:
                stats[i] = f[key]
    else:
        for k in range(p):
            others = [cols[j] for j in range(p) if j != k]
            part = oracle_ols(others, cols[k])
            res = np.full((n, n), np.nan)
            res.T[mask.T] = part["resid"]
            if mode == "graph":
                # Residuals live on the lower triangle; read symmetrically.
                res = np.where(np.tril(np.ones((n, n), bool), -1 if not diag else 0), res, res.T)
            for i, q in enumerate(perms):
                f = oracle_ols(others + [vec(res[np.ix_(q, q)], mask)], yv)
                if f is not None:
                    stats[i, k] = f[key][-1]
    return base, [count(stats[:, k], observed[k]) for k in range(p)]


# ---- Data -----------------------------------------------------------------

def random_case(rng, n, k, binary=False, missing=0.0, mode="digraph"):
    if binary:
        y = (rng.random((n, n)) < 0.4).astype(float)
    else:
        y = rng.normal(size=(n, n))
    xs = []
    for i in range(k):
        x = rng.normal(size=(n, n)) + (0.6 * y if i == 0 else 0)
        if binary and i == 1:
            x = (rng.random((n, n)) < 0.5).astype(float)
        xs.append(x)
    if mode == "graph":
        y = symmetrize_lower(y)
        xs = [symmetrize_lower(x) for x in xs]
    if missing:
        y = y.copy()
        cells = rng.random((n, n)) < missing
        y[cells] = np.nan
        if mode == "graph":
            y = symmetrize_lower(y)
    return y, xs


# ---- Checks ---------------------------------------------------------------

def check_ols(rng):
    print("== least squares against numpy")
    for n in (5, 12, 30):
        for mode in ("digraph", "graph"):
            for diag in (False, True):
                for intercept in (True, False):
                    for missing in (0.0, 0.1):
                        y, xs = random_case(rng, n, 3, missing=missing, mode=mode)
                        tag = f"n={n} {mode} diag={diag} intercept={intercept} missing={missing}"
                        mask, yy, mats = design(y, xs, mode, diag, intercept)
                        fit = oracle_ols([vec(m, mask) for m in mats], vec(yy, mask))
                        p = len(mats)
                        summ = oracle_summary(fit, intercept, p)
                        eng = engine_mrqap(y, xs, mode, diag, intercept, reps=0)
                        for name, a, b in (("coef", eng["coef"], fit["coef"]),
                                           ("se", eng["se"], fit["se"]),
                                           ("t", eng["t"], fit["t"])):
                            check(np.allclose(a, b, rtol=1e-9, atol=1e-11),
                                  f"ols {name} {tag}: {a} vs {b}")
                        check(eng["dyads"] == fit["N"], f"ols dyads {tag}")
                        check(eng["df"] == fit["N"] - p, f"ols df {tag}")
                        for name in ("mss", "r2", "adj_r2", "F"):
                            check(math.isclose(eng[name], summ[name], rel_tol=1e-9, abs_tol=1e-12),
                                  f"ols {name} {tag}: {eng[name]} vs {summ[name]}")
                        check(math.isclose(eng["rss"], fit["rss"], rel_tol=1e-9),
                              f"ols rss {tag}")


def within_binomial(p_mc, p_exact, reps):
    se = math.sqrt(max(p_exact * (1 - p_exact), 0.0) / reps)
    return abs(p_mc - p_exact) <= Z * se + 2.0 / reps


def check_cor_exact(rng):
    print("== QAP correlation: exhaustive enumeration")
    reps = 40000
    cases = [
        ("digraph", False, False), ("digraph", True, False), ("graph", False, False),
        ("graph", True, False), ("digraph", False, True), ("graph", False, True),
    ]
    for n in (5, 6):
        for mode, diag, binary in cases:
            y, xs = random_case(rng, n, 2, binary=binary, mode=mode)
            x = xs[1] if binary else xs[0]
            tag = f"cor n={n} {mode} diag={diag} binary={binary}"
            obs, exact = oracle_cor(y, x, mode, diag)
            eng = engine_qap_cor(y, x, mode, diag, exhaustive=True)
            check(math.isclose(eng["r"], obs, rel_tol=1e-12, abs_tol=1e-14),
                  f"{tag} r {eng['r']} vs {obs}")
            for key in ("le", "ge", "abs", "valid"):
                check(eng[key] == exact[key], f"{tag} exhaustive {key} {eng[key]} vs {exact[key]}")
            mc = engine_qap_cor(y, x, mode, diag, reps=reps, seed=12345)
            for key in ("le", "ge", "abs"):
                pe = exact[key] / exact["valid"]
                pm = mc[key] / mc["valid"]
                check(within_binomial(pm, pe, mc["valid"]),
                      f"{tag} monte carlo {key} {pm:.4f} vs exact {pe:.4f}")


def check_mrqap_exact(rng):
    print("== MRQAP: exhaustive enumeration (DSP and Y-permutation)")
    reps = 20000
    cases = [
        (5, "digraph", False, True, "dsp", False),
        (5, "digraph", False, True, "y", False),
        (5, "graph", False, True, "dsp", False),
        (6, "graph", False, True, "y", False),
        (5, "digraph", True, False, "dsp", False),
        (5, "digraph", False, True, "dsp", True),
        (6, "digraph", False, True, "dsp", False),
    ]
    for n, mode, diag, intercept, method, binary in cases:
        y, xs = random_case(rng, n, 2, binary=binary, mode=mode)
        for statistic in ("t", "beta"):
            tag = (f"mrqap n={n} {mode} diag={diag} intercept={intercept} {method} "
                   f"{statistic} binary={binary}")
            base, exact = oracle_mrqap_exact(y, xs, mode, diag, intercept, method, statistic)
            eng = engine_mrqap(y, xs, mode, diag, intercept, method, statistic,
                               exhaustive=True)
            check(np.allclose(eng["coef"], base["coef"], rtol=1e-9), f"{tag} coef")
            for k, ex in enumerate(exact):
                for key in ("le", "ge", "abs", "valid"):
                    check(int(eng[key][k]) == ex[key],
                          f"{tag} col {k} exhaustive {key} {eng[key][k]} vs {ex[key]}")
            if statistic == "beta":
                continue
            mc = engine_mrqap(y, xs, mode, diag, intercept, method, statistic, reps=reps,
                              seed=2024)
            for k, ex in enumerate(exact):
                for key in ("le", "ge", "abs"):
                    pe = ex[key] / ex["valid"]
                    pm = mc[key][k] / mc["valid"][k]
                    check(within_binomial(pm, pe, int(mc["valid"][k])),
                          f"{tag} col {k} monte carlo {key} {pm:.4f} vs exact {pe:.4f}")


def check_missing_exact(rng):
    print("== missing dyads under permutation")
    y, xs = random_case(rng, 5, 2, missing=0.15)
    for method in ("dsp", "y"):
        base, exact = oracle_mrqap_exact(y, xs, method=method)
        eng = engine_mrqap(y, xs, method=method, exhaustive=True)
        check(np.allclose(eng["coef"], base["coef"], rtol=1e-9), f"missing {method} coef")
        for k, ex in enumerate(exact):
            for key in ("le", "ge", "abs", "valid"):
                check(int(eng[key][k]) == ex[key],
                      f"missing {method} col {k} {key} {eng[key][k]} vs {ex[key]}")
    x = xs[0].copy()
    x[1, 3] = np.nan
    obs, exact = oracle_cor(y, x)
    eng = engine_qap_cor(y, x, exhaustive=True)
    check(math.isclose(eng["r"], obs, rel_tol=1e-12), "missing cor r")
    for key in ("le", "ge", "abs", "valid"):
        check(eng[key] == exact[key], f"missing cor {key} {eng[key]} vs {exact[key]}")


def check_threads(rng):
    print("== reproducibility across thread counts")
    y, xs = random_case(rng, 15, 3)
    for method in ("dsp", "y"):
        runs = [engine_mrqap(y, xs, method=method, reps=500, seed=99, threads=t, keep=True)
                for t in (1, 2, 7)]
        for r in runs[1:]:
            check(np.array_equal(r["dist"], runs[0]["dist"], equal_nan=True),
                  f"threads {method} distribution")
            check(np.array_equal(r["ge"], runs[0]["ge"]), f"threads {method} counts")
    runs = [engine_qap_cor(y, xs[0], reps=700, seed=5, threads=t, keep=True) for t in (1, 3)]
    check(runs[0]["dist"] == runs[1]["dist"], "threads cor distribution")


def main() -> int:
    global LIB
    LIB = _lib()
    rng = np.random.default_rng(20261022)
    check_ols(rng)
    check_cor_exact(rng)
    check_mrqap_exact(rng)
    check_missing_exact(rng)
    check_threads(rng)
    if failures:
        print(f"{len(failures)} QAP check(s) failed")
        return 1
    print("All QAP checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
