#!/usr/bin/env python3
"""Golden values for density, components, fragmentation, closeness,
harmonic centrality, centralization, and eigenvector centrality, computed
directly from their definitions.

This oracle shares no code or algorithmic strategy with the engine. All
distances come from the Floyd-Warshall recurrence in exact rational
arithmetic (the engine uses breadth-first search, Dial's buckets, or
Dijkstra's algorithm in floating point); components and fragmentation come
from the boolean transitive closure (the engine uses Tarjan's algorithm and
breadth-first search); betweenness counts shortest paths pair by pair with
exact integers; eigenvectors come from a dense eigendecomposition with
numpy (the engine uses power iteration).

Undefined values (NaN) are written as -1, as the engine runner does.

    python3 validation/oracles/topography_definition_oracle.py
"""

from __future__ import annotations

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
UNDEFINED = -1.0
INF = None  # no path


def manifest() -> dict:
    return yaml.safe_load((MANIFEST_DIR / "topography.yaml").read_text())


# ---- Reading a fixture -------------------------------------------------------

def arcs(fx, simple=False, symmetrize=False):
    """(source, target, weight) for every arc. Undirected edges give two arcs
    (a loop gives one). simple: the adjacency-matrix reading, with parallel
    edges merged (largest weight) and loops dropped. symmetrize: as simple,
    and every arc mirrored."""
    out = []
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        w = Fraction(fx.weights[i]) if fx.weighted else Fraction(1)
        out.append((a, b, w))
        if (not fx.directed or symmetrize) and a != b:
            out.append((b, a, w))
    if simple or symmetrize:
        merged = {}
        for a, b, w in out:
            if a != b:
                merged[(a, b)] = max(merged.get((a, b), w), w)
        out = [(a, b, w) for (a, b), w in sorted(merged.items())]
    return out


def directed_arcs(fx, mode, simple=False):
    """Arcs whose direction is the direction of the paths `mode` asks for:
    out as given, in reversed, all both ways."""
    base = arcs(fx, simple)
    if not fx.directed or mode == "out":
        return base
    if mode == "in":
        return [(b, a, w) for a, b, w in base]
    return base + [(b, a, w) for a, b, w in base if a != b]


# ---- Distances by Floyd-Warshall in exact arithmetic ------------------------

def distances(n, arc_list, weighted, strength=False):
    d = [[INF] * n for _ in range(n)]
    for v in range(n):
        d[v][v] = Fraction(0)
    for a, b, w in arc_list:
        if a == b:
            continue
        length = (1 / w if strength else w) if weighted else Fraction(1)
        if d[a][b] is INF or length < d[a][b]:
            d[a][b] = length
    for k in range(n):
        dk = d[k]
        for i in range(n):
            dik = d[i][k]
            if dik is INF:
                continue
            di = d[i]
            for j in range(n):
                if dk[j] is INF:
                    continue
                alt = dik + dk[j]
                if di[j] is INF or alt < di[j]:
                    di[j] = alt
    return d


def closeness(n, d, rule, normalized):
    longest = max((x for row in d for x in row if x is not INF), default=Fraction(0))
    out = []
    for i in range(n):
        reach = [d[i][j] for j in range(n) if j != i and d[i][j] is not INF]
        r, total = len(reach), sum(reach, Fraction(0))
        if rule == "ignore":
            if r == 0:
                out.append(UNDEFINED)
            else:
                out.append(float(Fraction(r) / total if normalized else 1 / total))
        elif rule == "wasserman_faust":
            out.append(0.0 if r == 0 else float((Fraction(r) / total) * Fraction(r, n - 1)))
        elif rule == "zero":
            if n < 2:
                out.append(UNDEFINED)
            elif r < n - 1:
                out.append(0.0)
            else:
                out.append(float(Fraction(n - 1) / total if normalized else 1 / total))
        elif rule == "max_plus_one":
            if n < 2:
                out.append(UNDEFINED)
                continue
            total += (n - 1 - r) * (longest + 1)
            out.append(float(Fraction(n - 1) / total if normalized else 1 / total))
        else:
            raise ValueError(rule)
    return out


def harmonic(n, d, normalized):
    out = []
    for i in range(n):
        h = sum((1 / d[i][j] for j in range(n) if j != i and d[i][j] is not INF), Fraction(0))
        if normalized and n > 1:
            h /= n - 1
        out.append(float(h))
    return out


# ---- Components from the transitive closure ---------------------------------

def closure(n, arc_list):
    reach = [[i == j for j in range(n)] for i in range(n)]
    for a, b, _ in arc_list:
        reach[a][b] = True
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                ri, rk = reach[i], reach[k]
                for j in range(n):
                    if rk[j]:
                        ri[j] = True
    return reach


def same_component(fx, kind):
    n = fx.n
    if kind == "weak" or not fx.directed:
        r = closure(n, arcs(fx, symmetrize=True))
        return [[r[i][j] for j in range(n)] for i in range(n)]
    r = closure(n, arcs(fx))
    return [[r[i][j] and r[j][i] for j in range(n)] for i in range(n)]


def membership(fx, kind):
    same = same_component(fx, kind)
    # Label each node by the lowest-numbered node in its component.
    return [float(min(j for j in range(fx.n) if same[i][j])) for i in range(fx.n)]


def fragmentation(fx, kind):
    n = fx.n
    if n < 2:
        return UNDEFINED
    same = same_component(fx, kind)
    apart = sum(1 for i in range(n) for j in range(n) if i != j and not same[i][j])
    return float(Fraction(apart, n * (n - 1)))


# ---- Density -----------------------------------------------------------------

def density(fx, loops, simple):
    n = fx.n
    if simple:
        pairs = {(a, b) if fx.directed else tuple(sorted((a, b)))
                 for a, b in zip(fx.src, fx.tgt) if a != b}
        m = len(pairs)
    else:
        m = len(fx.src)
    if n == 0 or (n == 1 and not loops):
        return UNDEFINED
    if fx.directed:
        possible = n * n if loops else n * (n - 1)
    else:
        possible = n * (n + 1) // 2 if loops else n * (n - 1) // 2
    return float(Fraction(m, possible))


# ---- Centralization ------------------------------------------------------------

def degrees(fx, mode, loops, simple):
    out = [0] * fx.n
    for a, b, _ in (arcs(fx, simple=True) if simple else
                    [(a, b, 1) for a, b in zip(fx.src, fx.tgt)]):
        if simple and not fx.directed and a > b:
            continue  # each undirected edge appears twice in arcs()
        if a == b and not loops:
            continue
        if not fx.directed or mode != "in":
            out[a] += 1
        if not fx.directed or mode != "out":
            out[b] += 1
    return out


def betweenness(fx, simple):
    """Freeman betweenness: for every ordered pair (s, t), the share of
    shortest s-t paths through each other node, with parallel edges as
    distinct paths unless simple. Halved for undirected graphs."""
    n = fx.n
    a = arcs(fx, simple)
    d = distances(n, a, weighted=False)
    mult = {}
    for s, t, _ in a:
        if s != t:
            mult[(s, t)] = mult.get((s, t), 0) + 1
    sigma = [[0] * n for _ in range(n)]
    for s in range(n):
        sigma[s][s] = 1
        order = sorted((v for v in range(n) if d[s][v] is not INF), key=lambda v: d[s][v])
        for t in order:
            if t == s:
                continue
            sigma[s][t] = sum(sigma[s][u] * k for (u, w), k in mult.items()
                              if w == t and d[s][u] is not INF and d[s][u] + 1 == d[s][t])
    bc = [Fraction(0)] * n
    for s in range(n):
        for t in range(n):
            if s == t or d[s][t] is INF:
                continue
            for v in range(n):
                if v in (s, t) or d[s][v] is INF or d[v][t] is INF:
                    continue
                if d[s][v] + d[v][t] == d[s][t]:
                    bc[v] += Fraction(sigma[s][v] * sigma[v][t], sigma[s][t])
    if not fx.directed:
        bc = [b / 2 for b in bc]
    return bc


def freeman_max(n, directed, measure, mode, loops, convention):
    """Largest possible sum of differences from the maximum, from each
    tool's documentation (see crates/engine/src/conventions.rs)."""
    mode = mode if directed else "all"
    loops = loops and convention == "igraph"
    if measure == "degree":
        if not directed:
            return (n - 1) * (n if loops else n - 2)
        if mode == "all":
            return 2 * (n - 1) * ((n - 1) if loops else (n - 2))
        return (n - 1) * (n if loops else n - 1)
    if measure == "betweenness":
        t = Fraction((n - 1) ** 2 * (n - 2))
        return t if directed else t / 2
    if convention == "ucinet" or mode == "all":
        return Fraction((n - 1) * (n - 2), 2 * n - 3)
    return Fraction((n - 1) ** 2, n)


def centralization(fx, e):
    conv = e.get("convention", "igraph")
    simple = conv in ("sna", "ucinet")
    measure = e["measure"]
    mode = e.get("mode", "all") if fx.directed else "all"
    n = fx.n
    if n == 0:
        return UNDEFINED
    if measure == "degree":
        loops = e.get("loops", True) and conv == "igraph"
        scores = [Fraction(x) for x in degrees(fx, mode, loops, simple)]
    elif measure == "betweenness":
        scores = betweenness(fx, simple)
    else:
        rule = "zero" if conv == "sna" else "ignore"
        d = distances(n, directed_arcs(fx, mode, simple), weighted=False)
        if conv == "ucinet" and any(x is INF for row in d for x in row):
            return UNDEFINED
        scores = []
        for i in range(n):
            reach = [d[i][j] for j in range(n) if j != i and d[i][j] is not INF]
            r, total = len(reach), sum(reach, Fraction(0))
            if rule == "ignore":
                if r == 0:
                    return UNDEFINED
                scores.append(Fraction(r) / total)
            else:
                scores.append(Fraction(0) if r < n - 1 or n < 2 else Fraction(n - 1) / total)
    tmax = freeman_max(n, fx.directed, measure, mode, e.get("loops", True), conv)
    if tmax == 0:
        return UNDEFINED
    top = max(scores)
    return float(sum(top - c for c in scores) / tmax)


# ---- Eigenvector centrality by dense eigendecomposition ---------------------

def eigenvector(fx, e):
    """The principal eigenvector of the matrix whose rows sum the scores a
    node receives. Returns None when the principal eigenvalue is repeated,
    since the eigenvector is then not unique and no definition applies."""
    n = fx.n
    if n == 0:
        return []
    mode = e.get("mode", "in") if fx.directed else "all"
    if e.get("symmetrize"):
        mode = "all"
    loops = e.get("loops", "twice")
    weighted = e.get("weighted", False) and fx.weighted
    a = np.zeros((n, n))
    if e.get("symmetrize"):
        edges = [(s, t, w) for s, t, w in arcs(fx, symmetrize=True) if s < t]
        undirected = True
    elif e.get("simple"):
        edges = arcs(fx, simple=True)
        if not fx.directed:
            edges = [(s, t, w) for s, t, w in edges if s < t]
        undirected = not fx.directed
    else:
        edges = [(s, t, Fraction(fx.weights[i]) if fx.weighted else Fraction(1))
                 for i, (s, t) in enumerate(zip(fx.src, fx.tgt))]
        undirected = not fx.directed
    if mode == "all" and not undirected:
        undirected = True  # direction ignored: A + A^T
    # An undirected loop is two edge endpoints (igraph) or one (NetworkX); a
    # directed loop is one arc.
    loop_factor = {"twice": 2.0 if undirected else 1.0, "once": 1.0, "ignore": 0.0}[loops]
    for s, t, w in edges:
        x = float(w) if weighted else 1.0
        if s == t:
            a[s, s] += loop_factor * x
        elif undirected:
            a[s, t] += x
            a[t, s] += x
        else:
            a[s, t] += x
    # Row i of m sums the scores node i receives.
    m = a.T if mode == "in" else a
    if undirected:
        vals, vecs = np.linalg.eigh(m)
    else:
        vals, vecs = np.linalg.eig(m)
    order = np.argsort(-vals.real)
    lead = vals[order[0]].real
    if lead <= 1e-12:
        # Every eigenvalue is 0 (no cycles): following igraph, nodes with no
        # outgoing ties in the direction of the mode score 1, others 0.
        has_out = (np.abs(m) > 0).any(axis=0)
        x = np.where(has_out, 0.0, 1.0)
    else:
        if n > 1 and abs(vals[order[1]].real - lead) <= 1e-8 * max(1.0, lead):
            return None
        x = vecs[:, order[0]].real
        if x.sum() < 0:
            x = -x
        x[np.abs(x) < 1e-13] = 0.0
    if e.get("scale", "max") == "max":
        x = x / x.max()
    else:
        x = x / np.linalg.norm(x)
    return [float(v) for v in x]


# ---- Driver --------------------------------------------------------------------

def compute(fx, e):
    fn = e["fn"]
    if fn == "density":
        return [density(fx, e["loops"], e.get("simple", False))]
    if fn == "components":
        return membership(fx, e["kind"])
    if fn == "fragmentation":
        return [fragmentation(fx, e["kind"])]
    if fn in ("closeness", "harmonic"):
        mode = e.get("mode", "out")
        d = distances(fx.n, directed_arcs(fx, mode, e.get("simple", False)),
                      e.get("weighted", False), e.get("weights_as") == "strength")
        if fn == "harmonic":
            return harmonic(fx.n, d, e.get("normalized", False))
        return closeness(fx.n, d, e.get("unreachable", "ignore"), e.get("normalized", False))
    if fn == "centralization":
        return [centralization(fx, e)]
    if fn == "eigenvector":
        return eigenvector(fx, e)
    if fn == "degree_simple":
        return [float(x) for x in degrees(fx, e["mode"], False, simple=True)]
    if fn == "betweenness_simple" and not e["weighted"]:
        return [float(b) for b in betweenness(fx, simple=True)]
    return None


def main() -> None:
    spec = manifest()
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        if fx.n > MAX_NODES:
            continue
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(fx, v["engine"])
            if values is None:
                continue
            write_values(ORACLE, name, variant, values)
            count += 1
    (GOLDEN / ORACLE).mkdir(parents=True, exist_ok=True)
    (GOLDEN / ORACLE / "VERSION").write_text(
        "definition (direct evaluation, exact arithmetic)\n", newline="")
    print(f"definition (topography): wrote {count} golden files")


if __name__ == "__main__":
    main()
