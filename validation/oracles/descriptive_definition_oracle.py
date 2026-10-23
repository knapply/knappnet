#!/usr/bin/env python3
"""Golden values for reciprocity, transitivity and clustering, distances,
mean distance, diameter, eccentricity, radius, reachability, cut points,
and bridges, computed directly from their definitions.

This oracle shares no code or algorithmic strategy with the engine:
- reciprocity counts dyads in an adjacency set, and its correlation form
  is a Pearson correlation computed by numpy (the engine uses the closed
  form of Garlaschelli and Loffredo);
- transitivity counts ordered two-paths and the closed ones by brute
  force, and local coefficients enumerate every pair of neighbors (the
  engine lists each triangle once over a degree ordering);
- weighted coefficients evaluate the published sums over ordered neighbor
  pairs, in exact rational arithmetic for Barrat's;
- distances come from the Floyd-Warshall recurrence in exact rational
  arithmetic (the engine uses breadth-first search, Dial's buckets, or
  Dijkstra's algorithm in floating point);
- reachability comes from the boolean transitive closure;
- cut points and bridges are found by deleting each node or edge and
  counting components (the engine uses one depth-first search).

NaN is written as -999 and infinity as -1000, as the engine runner does.

    python3 validation/oracles/descriptive_definition_oracle.py
"""

from __future__ import annotations

import itertools
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,
                    load_fixture, write_values, write_version)  # noqa: E402

ORACLE = "definition"
MAX_NODES = 80
UNDEFINED = -999.0
INFINITE = -1000.0
INF = None  # no path


def manifest() -> dict:
    return load_area("descriptive")


def mark(x):
    if x is INF:
        return INFINITE
    x = float(x)
    if math.isnan(x):
        return UNDEFINED
    if math.isinf(x):
        return INFINITE
    return x


# ---- Readings of a fixture ---------------------------------------------------

def arc_set(fx):
    """Distinct arcs between distinct nodes (both directions if undirected)."""
    out = set()
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            out.add((a, b))
            if not fx.directed:
                out.add((b, a))
    return out


def neighbor_sets(fx):
    """The simple undirected reading: distinct neighbors other than self."""
    nb = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            nb[a].add(b)
            nb[b].add(a)
    return nb


def weight_matrix(fx):
    """Total weight between two distinct nodes over every edge joining them,
    in either direction, as exact fractions."""
    w = {}
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        x = Fraction(fx.weights[i])
        key = (min(a, b), max(a, b))
        w[key] = w.get(key, Fraction(0)) + x
    return w


def path_arcs(fx, mode, weighted, strength):
    """(source, target, length) along which paths in `mode` run."""
    out = []
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if weighted:
            w = Fraction(fx.weights[i])
            length = 1 / w if strength else w
        else:
            length = Fraction(1)
        if not fx.directed or mode == "all":
            out += [(a, b, length), (b, a, length)]
        elif mode == "out":
            out.append((a, b, length))
        else:
            out.append((b, a, length))
    return out


# ---- Reciprocity ------------------------------------------------------------------

def reciprocity(fx, mode):
    n = fx.n
    arcs = arc_set(fx)
    dyads = [(i, j) for i in range(n) for j in range(i + 1, n)]
    mutual = sum(1 for i, j in dyads if (i, j) in arcs and (j, i) in arcs)
    asym = sum(1 for i, j in dyads if ((i, j) in arcs) != ((j, i) in arcs))
    null = len(dyads) - mutual - asym
    if mode == "default":
        total = 2 * mutual + asym
        return UNDEFINED if total == 0 else float(Fraction(2 * mutual, total))
    if mode == "ratio":
        total = mutual + asym
        return UNDEFINED if total == 0 else float(Fraction(mutual, total))
    if mode == "dyadic":
        return UNDEFINED if not dyads else float(Fraction(mutual + null, len(dyads)))
    if mode == "correlation":
        pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
        if not pairs:
            return UNDEFINED
        x = np.array([float((i, j) in arcs) for i, j in pairs])
        y = np.array([float((j, i) in arcs) for i, j in pairs])
        if x.std() == 0 or y.std() == 0:
            return UNDEFINED
        return float(np.corrcoef(x, y)[0, 1])
    raise ValueError(mode)


# ---- Transitivity and clustering --------------------------------------------

def transitivity(fx, undefined):
    nb = neighbor_sets(fx)
    paths = closed = 0
    for j in range(fx.n):
        for i in nb[j]:
            for k in nb[j]:
                if i != k:
                    paths += 1
                    closed += k in nb[i]
    if paths == 0:
        return 0.0 if undefined == "zero" else UNDEFINED
    return float(Fraction(closed, paths))


def local(fx, weighting):
    """Local coefficients as exact values or floats; None below two neighbors."""
    nb = neighbor_sets(fx)
    w = weight_matrix(fx) if weighting != "none" else None

    def wt(a, b):
        return w[(min(a, b), max(a, b))]

    top = max(w.values()) if w else None
    out = []
    for i in range(fx.n):
        k = len(nb[i])
        if k < 2:
            out.append(None)
            continue
        ordered = [(j, h) for j in nb[i] for h in nb[i] if j != h]
        if weighting == "none":
            adjacent = sum(1 for j, h in ordered if h in nb[j])
            out.append(Fraction(adjacent, k * (k - 1)))
        elif weighting == "barrat":
            s = sum(wt(i, j) for j in nb[i])
            total = sum((wt(i, j) + wt(i, h)) / 2 for j, h in ordered if h in nb[j])
            out.append(total / (s * (k - 1)))
        elif weighting == "onnela":
            total = sum(math.cbrt(float(wt(i, j) * wt(i, h) * wt(j, h) / top ** 3))
                        for j, h in ordered if h in nb[j])
            out.append(total / (k * (k - 1)))
        else:
            raise ValueError(weighting)
    return out


def clustering(fx, weighting, undefined):
    return [(0.0 if undefined == "zero" else UNDEFINED) if c is None else float(c)
            for c in local(fx, weighting)]


def average_clustering(fx, weighting, undefined):
    c = local(fx, weighting)
    if undefined == "zero":
        values = [Fraction(0) if x is None else x for x in c]
    else:
        values = [x for x in c if x is not None]
    if not values:
        return UNDEFINED
    if weighting == "onnela":
        return float(sum(float(v) for v in values) / len(values))
    return float(sum(values, Fraction(0)) / len(values))


# ---- Distances ----------------------------------------------------------------

def floyd_warshall(n, arcs):
    d = [[INF] * n for _ in range(n)]
    for v in range(n):
        d[v][v] = Fraction(0)
    for a, b, length in arcs:
        if a != b and (d[a][b] is INF or length < d[a][b]):
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


def distance_matrix(fx, e):
    return floyd_warshall(fx.n, path_arcs(fx, e.get("mode", "out"), e.get("weighted", False),
                                          e.get("weights_as") == "strength"))


def finite_pairs(d):
    return [d[i][j] for i in range(len(d)) for j in range(len(d)) if i != j and d[i][j] is not INF]


def mean_distance(d, unconnected):
    n = len(d)
    pairs = finite_pairs(d)
    if not unconnected and len(pairs) < n * (n - 1):
        return INFINITE
    if not pairs:
        return UNDEFINED
    return float(sum(pairs, Fraction(0)) / len(pairs))


def diameter(d, unconnected):
    n = len(d)
    if n == 0:
        return UNDEFINED
    pairs = finite_pairs(d)
    if not unconnected and len(pairs) < n * (n - 1):
        return INFINITE
    return float(max(pairs, default=Fraction(0)))


def eccentricity(d):
    return [max((x for x in row if x is not INF), default=Fraction(0)) for row in d]


# ---- Reachability ---------------------------------------------------------------

def closure(n, arcs):
    reach = [[i == j for j in range(n)] for i in range(n)]
    for a, b in arcs:
        reach[a][b] = True
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                ri, rk = reach[i], reach[k]
                for j in range(n):
                    if rk[j]:
                        ri[j] = True
    return reach


def reachability(fx, mode):
    arcs = list(zip(fx.src, fx.tgt))
    if not fx.directed or mode == "all":
        arcs += [(b, a) for a, b in arcs]
    r = closure(fx.n, arcs)
    if fx.directed and mode == "in":
        return [float(sum(r[i][v] for i in range(fx.n)) - 1) for v in range(fx.n)]
    return [float(sum(r[v]) - 1) for v in range(fx.n)]


# ---- Cut points and bridges by deletion ----------------------------------------

def count_components(n, edges, removed_node=None):
    """Components of the undirected graph on the nodes other than
    removed_node, by depth-first search over an adjacency list."""
    adj = [[] for _ in range(n)]
    for a, b in edges:
        if removed_node in (a, b):
            continue
        adj[a].append(b)
        adj[b].append(a)
    seen = [False] * n
    count = 0
    for s in range(n):
        if s == removed_node or seen[s]:
            continue
        count += 1
        stack = [s]
        seen[s] = True
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if not seen[w]:
                    seen[w] = True
                    stack.append(w)
    return count


def cut_points(fx):
    edges = [(a, b) for a, b in zip(fx.src, fx.tgt) if a != b]
    base = count_components(fx.n, edges)
    out = []
    for v in range(fx.n):
        # Removing a node removes one component if it was isolated.
        isolated = all(v not in e for e in edges)
        if count_components(fx.n, edges, v) > base - (1 if isolated else 0):
            out.append(float(v))
    return out


def bridges(fx):
    edges = list(zip(fx.src, fx.tgt))
    base = count_components(fx.n, [e for e in edges if e[0] != e[1]])
    out = []
    for k, (a, b) in enumerate(edges):
        if a == b:
            continue
        rest = [e for i, e in enumerate(edges) if i != k and e[0] != e[1]]
        if count_components(fx.n, rest) > base:
            out.append((float(min(a, b)), float(max(a, b))))
    return out


# ---- Driver --------------------------------------------------------------------

def compute(fx, e):
    fn = e["fn"]
    if fn == "reciprocity":
        return [reciprocity(fx, e["mode"])]
    if fn == "transitivity":
        return [transitivity(fx, e.get("undefined", "nan"))]
    if fn == "clustering":
        return clustering(fx, e.get("weighting", "none"), e.get("undefined", "nan"))
    if fn == "average_clustering":
        return [average_clustering(fx, e.get("weighting", "none"), e.get("undefined", "nan"))]
    if fn == "distances":
        return [tuple(mark(x) for x in row) for row in distance_matrix(fx, e)]
    if fn == "mean_distance":
        return [mean_distance(distance_matrix(fx, e), e.get("unconnected", True))]
    if fn == "diameter":
        return [diameter(distance_matrix(fx, e), e.get("unconnected", True))]
    if fn == "eccentricity":
        return [float(x) for x in eccentricity(distance_matrix(fx, e))]
    if fn == "radius":
        ecc = eccentricity(distance_matrix(fx, e))
        return [float(min(ecc)) if ecc else UNDEFINED]
    if fn == "reachability":
        return reachability(fx, e["mode"])
    if fn == "cut_points":
        return cut_points(fx)
    if fn == "bridges":
        return bridges(fx)
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
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (descriptive): wrote {count} golden files")


if __name__ == "__main__":
    main()
