#!/usr/bin/env python3
"""Validate the random graph generators behind random_*() and the null
models of cug_test() against their definitions.

The generators produce whole networks rather than node-level values, so
they are checked here instead of through the manifest. Three kinds of
checks, all independent of the engine (numpy only):

1. Exact constraints, on every replicate of every batch:
   - G(n, m) has exactly m distinct edges (no loops unless asked), and
     G(n, p) with p = 1 or p = 0 has every or no edge;
   - the configuration model reproduces the degree sequence exactly (the
     multigraph form, a loop counting 2), and its simple form has no loop or
     parallel edge while keeping the degrees;
   - U|MAN has exactly the requested numbers of mutual and asymmetric
     dyads;
   - Watts-Strogatz networks are simple with n k / 2 edges and p = 0 gives
     the ring lattice; preferential attachment gives a simple network with
     m (n - m) edges in which every node added after the initial star has
     degree at least m; block models tie only within blocks when the
     off-diagonal probabilities are 0.

2. Distributions, by exact chi-square tests on many replicates:
   - under G(n, m) every pair is tied with probability m / N over the N
     possible pairs, and under G(n, p) with probability p;
   - under U|MAN every dyad is mutual with probability M / D and asymmetric
     with probability A / D over the D dyads, and each direction of an
     asymmetric dyad is equally likely;
   - under the simple configuration model with regular degrees every pair
     is tied with probability d / (n - 1) (the uniform distribution over
     simple graphs with that degree sequence is symmetric in the nodes);
   - in the Watts-Strogatz model the fraction of lattice edges that were
     rewired is p.

3. Reproducibility: a seed gives the same replicates on 1 and 3 threads,
   and replicate r of a batch is the single draw of stream (seed, r).

    python3 validation/cug_check.py

Exits non-zero on any failure.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bindings.groups as b  # noqa: E402

# Chi-square tests are rejected at this p-value. With about 15 tests, a
# 1e-5 threshold has a false-alarm chance near 1.5e-4.
ALPHA = 1e-5

failures: list[str] = []


def check(ok: bool, what: str) -> None:
    if not ok:
        failures.append(what)
        print(f"FAIL  {what}")


def pairs_of(n: int, directed: bool, loops: bool = False):
    for a in range(n):
        for c in range(n):
            if a == c and not loops:
                continue
            if not directed and c < a:
                continue
            yield (a, c)


def is_simple(edges, directed: bool) -> bool:
    seen = set()
    for a, c in edges:
        if a == c:
            return False
        key = (a, c) if directed else (min(a, c), max(a, c))
        if key in seen:
            return False
        seen.add(key)
    return True


def degrees(edges, n: int):
    out = np.zeros(n, dtype=int)
    inn = np.zeros(n, dtype=int)
    for a, c in edges:
        out[a] += 1
        inn[c] += 1
    return out, inn


def chi_square_uniform(counts: np.ndarray, expected: np.ndarray, what: str) -> None:
    """Pearson chi-square of observed counts against expected counts."""
    mask = expected > 0
    stat = float((((counts - expected) ** 2) / np.where(mask, expected, 1))[mask].sum())
    df = int(mask.sum()) - 1
    if df <= 0:
        return
    p = float(chi2.sf(stat, df))
    check(p > ALPHA, f"{what}: chi-square {stat:.1f} on {df} df, p = {p:.2e}")


# ---- Exact constraints -------------------------------------------------------

def check_exact(rng: np.random.Generator) -> None:
    for trial in range(40):
        n = int(rng.integers(2, 25))
        directed = bool(rng.integers(0, 2))
        loops = bool(rng.integers(0, 2))
        seed = int(rng.integers(0, 2**40))
        total = n * (n - 1) // (1 if directed else 2) + (n if loops else 0)
        m = int(rng.integers(0, total + 1))
        for edges in b.random_gnm(n, m, directed, loops, reps=4, seed=seed):
            check(len(edges) == m, f"G(n, m): {len(edges)} edges, not {m} (n={n}, seed={seed})")
            keys = {(a, c) if directed else (min(a, c), max(a, c)) for a, c in edges}
            check(len(keys) == m, f"G(n, m): repeated edge (n={n}, m={m}, seed={seed})")
            check(loops or all(a != c for a, c in edges), f"G(n, m): loop without loops=True")
        full = b.random_gnp(n, 1.0, directed, loops, reps=1, seed=seed)[0]
        check(len(full) == total, f"G(n, 1): {len(full)} edges, not {total}")
        check(not b.random_gnp(n, 0.0, directed, loops, reps=1, seed=seed)[0],
              "G(n, 0) has edges")

        # Configuration model from the degrees of a random simple graph.
        base = b.random_gnm(n, min(m, n * (n - 1) // 2), False, False, reps=1, seed=seed)[0]
        d_out, d_in = degrees(base, n)
        und = d_out + d_in
        for edges in b.random_configuration(list(und), reps=3, seed=seed):
            o, i = degrees(edges, n)
            check(np.array_equal(o + i, und), f"configuration: degrees differ (n={n}, seed={seed})")
        try:
            simple = b.random_configuration(list(und), simple=True, max_tries=500, reps=2,
                                            seed=seed)
        except b.EngineError:
            simple = []
        for edges in simple:
            o, i = degrees(edges, n)
            check(np.array_equal(o + i, und), "simple configuration: degrees differ")
            check(is_simple(edges, False), "simple configuration: not simple")
        dbase = b.random_gnm(n, min(m, n * (n - 1)), True, False, reps=1, seed=seed)[0]
        d_out, d_in = degrees(dbase, n)
        for edges in b.random_configuration(list(d_out), list(d_in), reps=2, seed=seed):
            o, i = degrees(edges, n)
            check(np.array_equal(o, d_out) and np.array_equal(i, d_in),
                  f"directed configuration: degrees differ (n={n}, seed={seed})")

        # U|MAN.
        dyads = n * (n - 1) // 2
        mutual = int(rng.integers(0, dyads + 1))
        asym = int(rng.integers(0, dyads - mutual + 1))
        for edges in b.random_uman(n, mutual, asym, reps=3, seed=seed):
            arcs = set(edges)
            check(len(arcs) == len(edges) and all(a != c for a, c in edges), "U|MAN: not simple")
            mm = sum(1 for a, c in arcs if a < c and (c, a) in arcs)
            aa = sum(1 for a, c in arcs if (c, a) not in arcs)
            check((mm, aa) == (mutual, asym),
                  f"U|MAN: census ({mm}, {aa}) not ({mutual}, {asym}) (n={n}, seed={seed})")

        # Watts-Strogatz and preferential attachment.
        if n >= 4:
            k = 2 * int(rng.integers(1, (n - 1) // 2 + 1))
            p = float(rng.uniform(0, 1))
            for edges in b.random_small_world(n, k, p, reps=2, seed=seed):
                check(is_simple(edges, False) and len(edges) == n * k // 2,
                      f"Watts-Strogatz: not simple with n k / 2 edges (n={n}, k={k})")
            lattice = b.random_small_world(n, k, 0.0, reps=1, seed=seed)[0]
            o, i = degrees(lattice, n)
            check(bool(np.all(o + i == k)), "Watts-Strogatz p = 0: not a ring lattice")
            ring = {(min(v, (v + j) % n), max(v, (v + j) % n))
                    for v in range(n) for j in range(1, k // 2 + 1)}
            check(set(lattice) == ring, "Watts-Strogatz p = 0: not the ring lattice")
        ma = int(rng.integers(1, n))
        for edges in b.random_preferential(n, ma, reps=2, seed=seed):
            o, i = degrees(edges, n)
            check(is_simple(edges, False) and len(edges) == ma * (n - ma),
                  f"preferential attachment: not simple with m (n - m) edges (n={n}, m={ma})")
            check(bool(np.all((o + i)[ma + 1:] >= ma)), "preferential attachment: late node below m")

        # Block models.
        sizes = [n // 2, n - n // 2]
        probs = [[1.0, 0.0], [0.0, 1.0]]
        edges = b.random_sbm(sizes, probs, directed, loops, reps=1, seed=seed)[0]
        within = lambda v: v < sizes[0]  # noqa: E731
        check(all(within(a) == within(c) for a, c in edges), "SBM: tie across blocks with p = 0")
        expect = sum(s * (s - 1) // (1 if directed else 2) + (s if loops else 0) for s in sizes)
        check(len(edges) == expect, f"SBM: {len(edges)} edges, not {expect} (n={n})")


# ---- Distributions ----------------------------------------------------------

def check_distributions() -> None:
    reps = 4000
    n, m = 7, 9
    pairs = list(pairs_of(n, False))
    counts = np.zeros(len(pairs))
    index = {p: i for i, p in enumerate(pairs)}
    for edges in b.random_gnm(n, m, reps=reps, seed=11, threads=3):
        for e in edges:
            counts[index[e]] += 1
    chi_square_uniform(counts, np.full(len(pairs), reps * m / len(pairs)), "G(n, m) pair counts")

    p = 0.3
    counts = np.zeros(len(pairs))
    sizes = []
    for edges in b.random_gnp(n, p, reps=reps, seed=12, threads=3):
        sizes.append(len(edges))
        for e in edges:
            counts[index[e]] += 1
    chi_square_uniform(counts, np.full(len(pairs), reps * p), "G(n, p) pair counts")
    # The number of edges is binomial(N, p).
    N = len(pairs)
    hist = np.bincount(sizes, minlength=N + 1).astype(float)
    from scipy.stats import binom
    expected = reps * binom.pmf(np.arange(N + 1), N, p)
    # Pool the tails so every cell expects at least 5.
    lo = int(np.argmax(np.cumsum(expected) >= 5))
    hi = N - int(np.argmax(np.cumsum(expected[::-1]) >= 5))
    pooled_obs = np.array([hist[:lo + 1].sum()] + list(hist[lo + 1:hi]) + [hist[hi:].sum()])
    pooled_exp = np.array([expected[:lo + 1].sum()] + list(expected[lo + 1:hi])
                          + [expected[hi:].sum()])
    chi_square_uniform(pooled_obs, pooled_exp, "G(n, p) edge count")

    # U|MAN: each dyad mutual with probability M / D, asymmetric with A / D,
    # and the direction of an asymmetric dyad is a fair coin.
    n, mutual, asym = 6, 4, 5
    dyads = list(pairs_of(n, False))
    D = len(dyads)
    mut = np.zeros(D)
    asy = np.zeros(D)
    forward = 0
    for edges in b.random_uman(n, mutual, asym, reps=reps, seed=13, threads=3):
        arcs = set(edges)
        for i, (a, c) in enumerate(dyads):
            f, r = (a, c) in arcs, (c, a) in arcs
            if f and r:
                mut[i] += 1
            elif f or r:
                asy[i] += 1
                forward += f
    chi_square_uniform(mut, np.full(D, reps * mutual / D), "U|MAN mutual dyad positions")
    chi_square_uniform(asy, np.full(D, reps * asym / D), "U|MAN asymmetric dyad positions")
    total_asym = reps * asym
    chi_square_uniform(np.array([forward, total_asym - forward]),
                       np.array([total_asym / 2, total_asym / 2]), "U|MAN asymmetric direction")

    # Simple configuration model with a regular degree sequence: every
    # pair tied with probability d / (n - 1).
    n, d = 8, 3
    pairs = list(pairs_of(n, False))
    index = {p: i for i, p in enumerate(pairs)}
    counts = np.zeros(len(pairs))
    for edges in b.random_configuration([d] * n, simple=True, max_tries=10000, reps=reps,
                                        seed=14, threads=3):
        for a, c in edges:
            counts[index[(min(a, c), max(a, c))]] += 1
    chi_square_uniform(counts, np.full(len(pairs), reps * d / (n - 1)),
                       "simple configuration model pair counts (3-regular)")

    # Watts-Strogatz: the fraction of lattice edges rewired is p.
    n, k, p = 40, 4, 0.25
    ring = {(min(v, (v + j) % n), max(v, (v + j) % n)) for v in range(n) for j in range(1, 3)}
    kept = 0
    total = 0
    for edges in b.random_small_world(n, k, p, reps=500, seed=15, threads=3):
        kept += len(ring & set(edges))
        total += len(ring)
    rewired = total - kept
    chi_square_uniform(np.array([rewired, kept], dtype=float),
                       np.array([total * p, total * (1 - p)]), "Watts-Strogatz rewiring rate")


# ---- Reproducibility --------------------------------------------------------

def check_reproducibility() -> None:
    a = b.random_gnm(30, 60, reps=7, seed=99, threads=1)
    c = b.random_gnm(30, 60, reps=7, seed=99, threads=3)
    check(a == c, "G(n, m): thread count changed the replicates")
    for r in range(3):
        single = b.random_gnm(30, 60, reps=r + 1, seed=99, threads=1)[r]
        check(single == a[r], f"G(n, m): replicate {r} differs from the shorter batch")
    u1 = b.random_uman(12, 10, 20, reps=5, seed=5, threads=1)
    u3 = b.random_uman(12, 10, 20, reps=5, seed=5, threads=3)
    check(u1 == u3, "U|MAN: thread count changed the replicates")
    s1 = b.random_sbm([5, 6], [[0.5, 0.1], [0.1, 0.4]], reps=5, seed=6, threads=1)
    s3 = b.random_sbm([5, 6], [[0.5, 0.1], [0.1, 0.4]], reps=5, seed=6, threads=3)
    check(s1 == s3, "SBM: thread count changed the replicates")
    check(b.random_gnm(30, 60, reps=1, seed=98)[0] != a[0],
          "G(n, m): different seeds gave the same network")


def main() -> int:
    rng = np.random.default_rng(20261023)
    check_exact(rng)
    check_distributions()
    check_reproducibility()
    if failures:
        print(f"{len(failures)} random graph check(s) failed")
        return 1
    print("All random graph and CUG null model checks passed")
    return 0


if __name__ == "__main__":
    if {"-h", "--help"} & set(sys.argv[1:]):
        print(__doc__.strip())
        sys.exit(0)
    sys.exit(main())
