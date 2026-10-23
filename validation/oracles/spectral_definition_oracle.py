#!/usr/bin/env python3
"""Golden values for PageRank, Katz centrality, Bonacich power centrality,
hubs and authorities, assortativity, mixing matrices, the E-I index, and
within-group tie shares, computed directly from their definitions.

This oracle shares no code or algorithmic strategy with the engine. The
engine iterates (power iteration, fixed-point series, the Lanczos method);
here every walk-based measure is the solution of its defining linear
system or eigenproblem on the dense adjacency matrix, found with numpy
(LU solves and symmetric eigendecomposition). Assortativity, mixing
matrices, the E-I index, and tie shares are counted pair by pair in exact
rational arithmetic.

Undefined values are written as empty lines. Node attributes come from the
rules in validation/bindings/spectral.py (restated here, not imported).

    python3 validation/oracles/spectral_definition_oracle.py
"""

from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import GOLDEN, MANIFEST_DIR, applies, fixture_names, load_fixture  # noqa: E402
from common import write_values  # noqa: E402

ORACLE = "definition"
MAX_NODES = 80
UNDEFINED = ()


def rule(name, n):
    return {
        "mod3": [v % 3 for v in range(n)],
        "id": list(range(n)),
        "mod3plus1": [v % 3 + 1 for v in range(n)],
        "num7": [7 * v % 11 for v in range(n)],
        "num3": [3 * v % 5 for v in range(n)],
    }[name]


def categories(values):
    """Codes by first appearance."""
    seen = {}
    return [seen.setdefault(v, len(seen)) for v in values], len(seen)


def defined(values):
    return [UNDEFINED if (isinstance(v, float) and math.isnan(v)) else v for v in values]


# ---- Adjacency matrices -------------------------------------------------------

def adjacency(fx, weighted=False, loops="twice", undirected=False, simple=False):
    """A[i, j] = total weight of arcs i -> j. Undirected edges (and, with
    `undirected`, directed arcs) count both ways; an undirected self-loop
    contributes 2, 1, or 0 to the diagonal; a directed loop 1 or 0.
    `simple`: the binary adjacency matrix without its diagonal."""
    n = fx.n
    a = np.zeros((n, n))
    sym = undirected or not fx.directed
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        w = fx.weights[i] if (weighted and fx.weighted) else 1.0
        if s == t:
            if simple:
                continue
            factor = {"twice": 2.0 if sym else 1.0, "once": 1.0, "ignore": 0.0}[loops]
            a[s, s] += factor * w
        elif simple:
            a[s, t] = 1.0
            if sym:
                a[t, s] = 1.0
        else:
            a[s, t] += w
            if sym:
                a[t, s] += w
    return a


def spectral_radius(a):
    if a.shape[0] == 0:
        return 0.0
    return float(max(abs(np.linalg.eigvals(a))))


# ---- Walk-based measures ----------------------------------------------------------

def pagerank(fx, e):
    n = fx.n
    a = adjacency(fx, e.get("weighted", False), e.get("loops", "twice"),
                  undirected=not e.get("directed", True))
    d = float(e.get("damping", 0.85))
    p = np.ones(n) if not e.get("personalization") else np.array(
        rule(e["personalization"], n), dtype=float)
    p = p / p.sum()
    out = a.sum(axis=1)
    dangling = out == 0
    jump = p if e.get("dangling", "personalized") == "personalized" else np.ones(n) / n
    # Column-stochastic transition matrix M: M[i, j] = probability of
    # moving from j to i.
    m = np.zeros((n, n))
    for j in range(n):
        if dangling[j]:
            m[:, j] = jump
        else:
            m[:, j] = a[j, :] / out[j]
    x = np.linalg.solve(np.eye(n) - d * m, (1 - d) * p)
    return list(x / x.sum())


def relative(e, key, a):
    if key in e:
        return float(e[key])
    lam = spectral_radius(a)
    return float(e[key + "_fraction"]) / lam if lam > 0 else float(e[key + "_fraction"])


def katz(fx, e):
    n = fx.n
    mode = e.get("mode", "in")
    a = adjacency(fx, e.get("weighted", False), e.get("loops", "once"),
                  undirected=(mode == "all"))
    alpha = relative(e, "alpha", a)
    b = e.get("beta", 1.0)
    beta = np.array(rule(b, n), dtype=float) if isinstance(b, str) else np.full(n, float(b))
    m = a.T if mode == "in" else a
    x = np.linalg.solve(np.eye(n) - alpha * m, beta)
    if e.get("normalized", True):
        x = x / np.linalg.norm(x)
    return list(x)


def bonacich(fx, e):
    n = fx.n
    mode = e.get("mode", "out")
    a = adjacency(fx, e.get("weighted", False), e.get("loops", "ignore"),
                  undirected=(mode == "all"), simple=e.get("simple", False))
    if mode == "in":
        a = a.T
    beta = relative(e, "exponent", a)
    b = a @ np.ones(n)
    if not b.any():
        return [float("nan")] * n
    system = np.eye(n) - beta * a
    if np.linalg.cond(system) > 1e12:
        return None
    y = np.linalg.solve(system, b)
    if e.get("rescale", False):
        return list(y / y.sum())
    return list(y * math.sqrt(n / float(y @ y)))


def hits(fx, e):
    n = fx.n
    a = adjacency(fx, e.get("weighted", False), e.get("loops", "twice"))
    if not (a > 0).any():
        return [float("nan")] * n
    vals, vecs = np.linalg.eigh(a.T @ a)
    order = np.argsort(vals)[::-1]
    lead = vals[order[0]]
    if n > 1 and abs(vals[order[1]] - lead) <= 1e-8 * lead:
        return None
    auth = vecs[:, order[0]]
    if auth.sum() < 0:
        auth = -auth
    auth[np.abs(auth) < 1e-14] = 0.0
    hubs = a @ auth
    x = hubs if e["which"] == "hubs" else auth
    scale = e.get("scale", "max")
    x = x / (x.max() if scale == "max" else x.sum() if scale == "sum" else np.linalg.norm(x))
    return list(x)


# ---- Assortativity in exact arithmetic --------------------------------------------

def arcs(fx, directed, loops):
    """(source, target) for every arc: undirected edges (or, ignoring
    direction, directed arcs) both ways; an undirected loop 2, 1, or 0
    times; a directed loop once."""
    sym = not (fx.directed and directed)
    out = []
    for s, t in zip(fx.src, fx.tgt):
        if s == t:
            copies = {"twice": 2, "once": 1, "ignore": 0}[loops] if sym else 1
            out += [(s, s)] * copies
        elif sym:
            out += [(s, t), (t, s)]
        else:
            out.append((s, t))
    return out


def degrees(fx, weighted):
    """Out-, in-, and total degree (or strength); an undirected loop adds 2."""
    n = fx.n
    out, inc = [Fraction(0)] * n, [Fraction(0)] * n
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        w = Fraction(fx.weights[i]) if weighted else Fraction(1)
        out[s] += w
        inc[t] += w
        if not fx.directed:
            out[t] += w
            inc[s] += w
    total = [o + i for o, i in zip(out, inc)] if fx.directed else out
    return out, inc, total


def pearson(pairs, normalized=True):
    m = len(pairs)
    if m == 0:
        return float("nan")
    mx = sum(x for x, _ in pairs) / m
    my = sum(y for _, y in pairs) / m
    cxy = sum((x - mx) * (y - my) for x, y in pairs)
    if not normalized:
        return float(cxy / m)
    cxx = sum((x - mx) ** 2 for x, _ in pairs)
    cyy = sum((y - my) ** 2 for _, y in pairs)
    if cxx == 0 or cyy == 0:
        return float("nan")
    return float(cxy) / math.sqrt(float(cxx) * float(cyy))


def assortativity(fx, e):
    fn = e["fn"]
    directed, loops = e.get("directed", True), e.get("loops", "twice")
    normalized = e.get("normalized", True)
    a = arcs(fx, directed, loops)
    if fn == "assortativity_nominal":
        types = rule(e["types"], fx.n)
        m = len(a)
        if m == 0:
            return float("nan")
        cats = sorted(set(types))
        trace = Fraction(sum(1 for s, t in a if types[s] == types[t]), m)
        ab = sum(Fraction(sum(1 for s, _ in a if types[s] == c), m)
                 * Fraction(sum(1 for _, t in a if types[t] == c), m) for c in cats)
        if not normalized:
            return float(trace - ab)
        return float((trace - ab) / (1 - ab)) if ab != 1 else float("nan")
    if fn == "assortativity_values":
        x = [Fraction(v) for v in rule(e["values"], fx.n)]
        y = [Fraction(v) for v in rule(e["values_in"], fx.n)] if e.get("values_in") else x
        return pearson([(x[s], y[t]) for s, t in a], normalized)
    out, inc, total = degrees(fx, e.get("weighted", False))
    if fx.directed and directed:
        return pearson([(out[s], inc[t]) for s, t in a], normalized)
    return pearson([(total[s], total[t]) for s, t in a], normalized)


# ---- Mixing matrices, E-I, shares -------------------------------------------------

def mixing(fx, e):
    """Counted the way homophily documents its mixing matrices."""
    d1, k1 = categories(rule(e["dim1"], fx.n))
    same = "dim2" not in e
    d2, k2 = (d1, k1) if same else categories(rule(e["dim2"], fx.n))
    direction = e.get("direction", "all")
    m = [[0] * k2 for _ in range(k1)]
    pairs = []
    for s, t in zip(fx.src, fx.tgt):
        if not fx.directed:
            pairs += [(s, t), (t, s)]
        elif not same and direction == "all":
            pairs += [(s, t), (t, s)]
        elif not same and direction == "in":
            pairs.append((t, s))
        else:
            pairs.append((s, t))
    for a, b in pairs:
        m[d1[a]][d2[b]] += 1
    if not fx.directed and same:
        m = [[(m[i][j] + m[j][i]) if j > i else (m[i][i] if i == j else 0)
              for j in range(k1)] for i in range(k1)]
    if fx.directed and direction == "in":
        m = [list(col) for col in zip(*m)] if m else []
    return [tuple(float(x) for x in row) for row in m]


def ei(fx, e):
    groups = rule(e["groups"], fx.n)
    scope = e.get("scope", "global")
    direction = e.get("direction", "all")
    counting = e.get("counting", "standard")
    directed = fx.directed
    sent = not directed or direction != "in"
    received = not directed or direction != "out"
    edges = list(zip(fx.src, fx.tgt))

    def idx(i, ex):
        return float("nan") if i + ex == 0 else (ex - i) / (ex + i)

    if scope == "global":
        i = sum(1 for s, t in edges if groups[s] == groups[t])
        return [idx(i, len(edges) - i)]
    labels = sorted(set(groups))
    if scope == "group" and counting == "homophily":
        rows = mixing(fx, {"dim1": e["groups"], "direction": direction})
        return [idx(rows[c][c], sum(rows[c][c + 1:])) for c in range(len(rows))]
    if scope == "group":
        out = []
        for c in labels:
            i = sum(1 for s, t in edges if groups[s] == c and groups[t] == c)
            ex = sum(1 for s, t in edges if groups[s] != groups[t]
                     and ((sent and groups[s] == c) or (received and groups[t] == c)))
            out.append(idx(i, ex))
        return out
    out = []
    for v in range(fx.n):
        i = ex = 0
        for s, t in edges:
            if s == t == v:
                same_side = sent and received and counting == "homophily"
                i += 2 if same_side else 1
                continue
            if (sent and s == v) or (received and t == v):
                other = t if s == v else s
                if groups[other] == groups[v]:
                    i += 1
                else:
                    ex += 1
        out.append(idx(i, ex))
    return out


def shares(fx, e):
    groups = rule(e["groups"], fx.n)
    direction = e.get("direction", "all")
    sent = not fx.directed or direction != "in"
    received = not fx.directed or direction != "out"
    out = []
    for c in sorted(set(groups)):
        size = sum(1 for g in groups if g == c)
        internal = external = 0
        for s, t in zip(fx.src, fx.tgt):
            if s == t:
                continue
            for member, other, counts in ((s, t, sent), (t, s, received)):
                if counts and groups[member] == c:
                    if groups[other] == c:
                        internal += 1
                    else:
                        external += 1
        total = internal + external
        observed = internal / total if total else -1.0
        expected = (size - 1) / (fx.n - 1) if fx.n >= 2 else -1.0
        out.append((float(size), float(internal), float(external), observed, expected))
    return out


def compute(fx, e):
    fn = e["fn"]
    if fn == "pagerank":
        return pagerank(fx, e)
    if fn == "katz":
        return katz(fx, e)
    if fn == "bonacich":
        return bonacich(fx, e)
    if fn == "hits":
        return hits(fx, e)
    if fn.startswith("assortativity"):
        return [assortativity(fx, e)]
    if fn == "mixing_matrix":
        return mixing(fx, e)
    if fn == "ei_index":
        return ei(fx, e)
    if fn == "group_shares":
        return shares(fx, e)
    return None


def main() -> None:
    spec = yaml.safe_load((MANIFEST_DIR / "spectral.yaml").read_text())
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        if fx.n > MAX_NODES or fx.n == 0:
            continue
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(fx, v["engine"])
            if values is None:
                continue
            write_values(ORACLE, name, variant, defined(values))
            count += 1
    (GOLDEN / ORACLE).mkdir(parents=True, exist_ok=True)
    (GOLDEN / ORACLE / "VERSION").write_text(
        "definition (direct evaluation, exact arithmetic)\n", newline="")
    print(f"definition (spectral): wrote {count} golden files")


if __name__ == "__main__":
    main()
