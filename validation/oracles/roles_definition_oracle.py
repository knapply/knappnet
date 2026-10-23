#!/usr/bin/env python3
"""Golden values for reach, proximity prestige, edge betweenness,
information centrality, bridging, Burt's hierarchy, eigenvector
centralization, regular and automorphic equivalence, blockmodels,
classical scaling, and geospatial measures, computed from their
definitions (validation/manifest.d/roles.yaml).

This oracle shares no code with the engine and uses different methods:

- Distances come from the Floyd-Warshall recurrence in exact rational
  arithmetic, and shortest-path counts from a dynamic program over nodes in
  order of exact distance (the engine uses breadth-first search, Dial's
  buckets, or Dijkstra's algorithm in floating point).
- Edge and node betweenness sum sigma_su sigma_vt / sigma_st pair by pair
  over every arc (u, v) on a shortest s-t path, in exact arithmetic (the
  engine uses Brandes' accumulation).
- Information centrality inverts B = L + J with numpy (the engine uses a
  Cholesky factorization of its own).
- Hierarchy follows Burt's constraint terms in exact arithmetic before the
  logarithms.
- Eigenvector centralization uses numpy's dense symmetric eigensolver (the
  engine uses power iteration).
- Regular equivalence is the greatest fixed point of the defining
  condition on pairs of nodes, found by deleting violating pairs (the
  engine refines a partition, CATREGE). Automorphic equivalence enumerates
  every permutation of networks with at most 8 nodes.
- Classical scaling uses numpy's eigh (the engine uses Householder
  tridiagonalization and implicit QL iterations).
- Geospatial distances use Vincenty's atan2 formula for a sphere (the
  engine uses the haversine formula), and proximity networks compare every
  pair of points (the engine uses a k-d tree).

Undefined values (NaN) are written as empty lines.

    python3 validation/oracles/roles_definition_oracle.py
"""

from __future__ import annotations

import itertools
import math
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "bindings"))
from common import GOLDEN, MANIFEST_DIR, applies, fixture_names, load_fixture  # noqa: E402
from common import write_values  # noqa: E402

ORACLE = "definition"
MAX_NODES = 80
EARTH_RADIUS_M = 6371008.7714
INF = None


def manifest() -> dict:
    return yaml.safe_load((MANIFEST_DIR / "roles.yaml").read_text())


# ---- Inputs derived from a fixture (data rules shared with the bindings) ------

def partition_for(rule: str, n: int) -> list[int]:
    if rule == "mod3":
        return [v % 3 for v in range(n)]
    if rule == "mod2":
        return [v % 2 for v in range(n)]
    raise ValueError(rule)


def hash_unit(v: int, salt: int) -> float:
    h = ((v + 1) * 2654435761 + salt * 40503 + 12345) % 4294967296
    h = (h ^ (h >> 13)) * 1274126177 % 4294967296
    return h / 4294967296.0


def coordinates(crs: str, n: int):
    u = [hash_unit(v, 1) for v in range(n)]
    w = [hash_unit(v, 2) for v in range(n)]
    if crs == "lonlat":
        return ([360.0 * a - 180.0 for a in u],
                [math.degrees(math.asin(2.0 * b - 1.0)) for b in w])
    return ([100.0 * a for a in u], [100.0 * b for b in w])


def mds_input(fx):
    """Hop distances on the undirected simple reading, by repeated
    relaxation; unreachable pairs get the largest finite distance plus 1."""
    n = fx.n
    d = [[math.inf] * n for _ in range(n)]
    for v in range(n):
        d[v][v] = 0.0
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            d[a][b] = d[b][a] = 1.0
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if d[i][k] + d[k][j] < d[i][j]:
                    d[i][j] = d[i][k] + d[k][j]
    finite = [x for row in d for x in row if x != math.inf]
    fill = max(finite) + 1.0 if finite else 1.0
    return [[fill if x == math.inf else x for x in row] for row in d]


# ---- Arcs, distances, and path counts --------------------------------------------

def arcs(fx, mode="out", weighted=False, strength=False):
    """(u, v, edge, length) for every arc in the direction of `mode`.
    Undirected edges give two arcs; loops are dropped (they are never on a
    shortest path)."""
    out = []
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        if weighted:
            w = Fraction(fx.weights[e])
            length = 1 / w if strength else w
        else:
            length = Fraction(1)
        if not fx.directed or mode == "all":
            out.append((a, b, e, length))
            out.append((b, a, e, length))
        elif mode == "out":
            out.append((a, b, e, length))
        else:
            out.append((b, a, e, length))
    return out


def floyd(n, arc_list):
    d = [[INF] * n for _ in range(n)]
    for v in range(n):
        d[v][v] = Fraction(0)
    for u, v, _, length in arc_list:
        if d[u][v] is INF or length < d[u][v]:
            d[u][v] = length
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


def path_counts(n, arc_list, d):
    """sigma[s][t]: number of shortest s-t paths, each arc distinct."""
    sigma = [[0] * n for _ in range(n)]
    into = defaultdict(list)
    for u, v, _, length in arc_list:
        into[v].append((u, length))
    for s in range(n):
        order = sorted((t for t in range(n) if d[s][t] is not INF), key=lambda t: d[s][t])
        sigma[s][s] = 1
        for t in order:
            if t == s:
                continue
            sigma[s][t] = sum(sigma[s][u] for u, length in into[t]
                              if d[s][u] is not INF and d[s][u] + length == d[s][t])
    return sigma


def edge_betweenness(fx, weighted, strength, normalized):
    n, m = fx.n, len(fx.src)
    al = arcs(fx, "out", weighted, strength)
    d = floyd(n, al)
    sigma = path_counts(n, al, d)
    eb = [Fraction(0)] * m
    for u, v, e, length in al:
        for s in range(n):
            if d[s][u] is INF:
                continue
            for t in range(n):
                if t == s or d[v][t] is INF or d[s][t] is INF:
                    continue
                if d[s][u] + length + d[v][t] == d[s][t]:
                    eb[e] += Fraction(sigma[s][u] * sigma[v][t], sigma[s][t])
    scale = Fraction(1) if fx.directed else Fraction(1, 2)
    if normalized:
        if n < 2:
            return [math.nan] * m
        pairs = n * (n - 1) if fx.directed else Fraction(n * (n - 1), 2)
        scale = scale / pairs
    return [float(x * scale) for x in eb]


def node_betweenness(fx):
    n = fx.n
    al = arcs(fx, "out")
    d = floyd(n, al)
    sigma = path_counts(n, al, d)
    b = [Fraction(0)] * n
    for v in range(n):
        for s in range(n):
            if s == v or d[s][v] is INF:
                continue
            for t in range(n):
                if t in (s, v) or d[v][t] is INF or d[s][t] is INF:
                    continue
                if d[s][v] + d[v][t] == d[s][t]:
                    b[v] += Fraction(sigma[s][v] * sigma[v][t], sigma[s][t])
    if not fx.directed:
        b = [x / 2 for x in b]
    return b


def reach(fx, spec):
    k = spec["k"]
    weighted = spec.get("weighted", False)
    strength = spec.get("weights_as") == "strength"
    mode = spec.get("mode", "out")
    d = floyd(fx.n, arcs(fx, mode, weighted, strength))
    limit = None if k == "inf" else Fraction(k)
    out = []
    for i in range(fx.n):
        c = sum(1 for j in range(fx.n) if j != i and d[i][j] is not INF
                and (limit is None or d[i][j] <= limit))
        if spec.get("normalized"):
            out.append(math.nan if fx.n < 2 else c / (fx.n - 1))
        else:
            out.append(float(c))
    return out


def proximity_prestige(fx, spec):
    """Lin (1976): (|I_i| / (n - 1)) / (mean distance from the influence
    domain I_i, the nodes that reach i); 0 for an empty domain."""
    weighted = spec.get("weighted", False)
    d = floyd(fx.n, arcs(fx, "out", weighted))
    n = fx.n
    out = []
    for i in range(n):
        dom = [d[j][i] for j in range(n) if j != i and d[j][i] is not INF]
        if not dom:
            out.append(0.0)
            continue
        size = len(dom)
        out.append(float(Fraction(size, n - 1) / (sum(dom, Fraction(0)) / size)))
    return out


# ---- Information centrality and bridging ------------------------------------------

def information_centrality(fx, weighted):
    n = fx.n
    total = defaultdict(float)
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        key = (a, b) if fx.directed or a < b else (b, a)
        if weighted:
            total[key] += fx.weights[e]
        else:
            total[key] = 1.0
    ties = {}
    for (a, b), x in total.items():
        key = (min(a, b), max(a, b))
        ties[key] = max(ties.get(key, 0.0), x)
    kept = sorted({v for pair in ties for v in pair})
    out = [0.0] * n
    if not kept:
        return out
    idx = {v: i for i, v in enumerate(kept)}
    # Connectivity by depth-first search.
    nbrs = defaultdict(set)
    for a, b in ties:
        nbrs[a].add(b)
        nbrs[b].add(a)
    seen, stack = {kept[0]}, [kept[0]]
    while stack:
        v = stack.pop()
        for w in nbrs[v]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    if len(seen) != len(kept):
        for v in kept:
            out[v] = math.nan
        return out
    k = len(kept)
    a = np.zeros((k, k))
    for (u, v), x in ties.items():
        a[idx[u], idx[v]] = a[idx[v], idx[u]] = x
    b = np.diag(a.sum(axis=1)) - a + np.ones((k, k))
    c = np.linalg.inv(b)
    trace = np.trace(c)
    row = c.sum(axis=1).mean()
    for v in kept:
        i = idx[v]
        out[v] = 1.0 / (c[i, i] + (trace - 2.0 * row) / k)
    return out


def distinct_neighbors(fx):
    nb = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            nb[a].add(b)
            nb[b].add(a)
    return nb


def bridging_coefficient(fx):
    nb = distinct_neighbors(fx)
    out = []
    for v in range(fx.n):
        if not nb[v]:
            out.append(None)
            continue
        s = sum(Fraction(1, len(nb[u])) for u in nb[v])
        out.append(Fraction(1, len(nb[v])) / s)
    return out


def bridging(fx, centrality):
    coef = bridging_coefficient(fx)
    if not centrality:
        return [math.nan if c is None else float(c) for c in coef]
    b = node_betweenness(fx)
    return [math.nan if c is None else float(c * x) for c, x in zip(coef, b)]


# ---- Hierarchy -----------------------------------------------------------------

def hierarchy(fx, weighted):
    n = fx.n
    a = [[Fraction(0)] * n for _ in range(n)]
    for e, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if s == t:
            continue
        w = Fraction(fx.weights[e]) if weighted else Fraction(1)
        a[s][t] += w
        if not fx.directed:
            a[t][s] += w
    z = [[a[i][j] + a[j][i] for j in range(n)] for i in range(n)]
    tot = [sum(row, Fraction(0)) for row in z]
    p = [[z[i][j] / tot[i] if tot[i] else Fraction(0) for j in range(n)] for i in range(n)]
    out = []
    for i in range(n):
        contacts = [j for j in range(n) if j != i and z[i][j] != 0]
        if len(contacts) < 2:
            out.append(math.nan)
            continue
        c = []
        for j in contacts:
            indirect = sum((p[i][q] * p[q][j] for q in contacts if q != j), Fraction(0))
            c.append((p[i][j] + indirect) ** 2)
        big_c = sum(c, Fraction(0))
        nn = len(contacts)
        h = 0.0
        for cj in c:
            r = float(cj * nn / big_c)
            h += r * math.log(r)
        out.append(h / (nn * math.log(nn)))
    return out


# ---- Eigenvector centralization ------------------------------------------------

def eigenvector_centralization(fx, normalized):
    """igraph convention: direction ignored, loops counted twice, scores
    scaled to a maximum of 1, theoretical maximum n - 2. Returns None when
    the leading eigenvalue is not simple (the scores are then not unique)."""
    n = fx.n
    if n == 0:
        return math.nan
    a = np.zeros((n, n))
    for s, t in zip(fx.src, fx.tgt):
        if s == t:
            a[s, s] += 2.0
        else:
            a[s, t] += 1.0
            a[t, s] += 1.0
    if not a.any():
        scores = np.ones(n)
    else:
        vals, vecs = np.linalg.eigh(a)
        if vals[-1] - vals[-2] < 1e-9 * max(1.0, abs(vals[-1])):
            return None
        v = vecs[:, -1]
        v = v if v.sum() > 0 else -v
        scores = v / np.abs(v).max()
        scores[np.abs(scores) < 1e-13] = 0.0
    raw = float((scores.max() - scores).sum())
    if not normalized:
        return raw
    return math.nan if n < 3 else raw / (n - 2)


# ---- Regular and automorphic equivalence ------------------------------------------

def tie_values(fx, weighted):
    val = defaultdict(float)
    for e, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if s == t:
            continue
        x = fx.weights[e] if weighted else 1.0
        for a, b in ([(s, t)] if fx.directed else [(s, t), (t, s)]):
            val[(a, b)] = val[(a, b)] + x if weighted else 1.0
    return val


def regular_equivalence(fx, seed, weighted):
    """Greatest relation R with: (i, j) in R implies equal seed labels and,
    for every k with tie type t = (a_ik, a_ki) != (0, 0), some m with the
    same tie type from j and (k, m) in R, and symmetrically."""
    n = fx.n
    val = tie_values(fx, weighted)
    typ = {}
    for i in range(n):
        for k in range(n):
            if i != k:
                t = (val.get((i, k), 0.0), val.get((k, i), 0.0))
                if t != (0.0, 0.0):
                    typ[(i, k)] = t
    alters = defaultdict(list)
    for (i, k), t in typ.items():
        alters[i].append((k, t))
    labels = partition_for(seed, n) if seed else [0] * n
    rel = [[labels[i] == labels[j] for j in range(n)] for i in range(n)]

    def matched(i, j):
        for k, t in alters[i]:
            if not any(t2 == t and rel[k][m] for m, t2 in alters[j]):
                return False
        return True

    changed = True
    while changed:
        changed = False
        for i in range(n):
            for j in range(n):
                if rel[i][j] and not (matched(i, j) and matched(j, i)):
                    rel[i][j] = rel[j][i] = False
                    changed = True
    out, cls = [], {}
    for i in range(n):
        rep = next(j for j in range(n) if rel[i][j])
        out.append(float(cls.setdefault(rep, len(cls))))
    return out


def automorphic_equivalence(fx):
    n = fx.n
    if n > 8:
        return None
    arcset = {(s, t) for s, t in zip(fx.src, fx.tgt)}
    if not fx.directed:
        arcset |= {(t, s) for s, t in arcset}
    orbit = list(range(n))
    for perm in itertools.permutations(range(n)):
        if all(((perm[s], perm[t]) in arcset) for s, t in arcset):
            for v in range(n):
                a, b = orbit[v], orbit[perm[v]]
                if a != b:
                    lo, hi = min(a, b), max(a, b)
                    orbit = [lo if x == hi else x for x in orbit]
    cls = {}
    return [float(cls.setdefault(o, len(cls))) for o in orbit]


# ---- Blockmodels ----------------------------------------------------------------

def block_model(fx, spec):
    n = fx.n
    weighted, loops = spec.get("weighted", False), spec.get("loops", False)
    labels = partition_for(spec["partition"], n)
    blocks = {}
    b = [blocks.setdefault(x, len(blocks)) for x in labels]
    k = len(blocks)
    a = [[Fraction(0)] * n for _ in range(n)]
    for e, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if s == t and not loops:
            continue
        x = Fraction(fx.weights[e]) if weighted else Fraction(1)
        pairs_ = [(s, t)] if fx.directed or s == t else [(s, t), (t, s)]
        for i, j in pairs_:
            a[i][j] = a[i][j] + x if weighted else Fraction(1)
    cells = [(i, j) for i in range(n) for j in range(n) if loops or i != j]
    sums = [[Fraction(0)] * k for _ in range(k)]
    cnt = [[0] * k for _ in range(k)]
    for i, j in cells:
        sums[b[i]][b[j]] += a[i][j]
        cnt[b[i]][b[j]] += 1
    dens = [[sums[r][s] / cnt[r][s] if cnt[r][s] else None for s in range(k)]
            for r in range(k)]
    if spec["fn"] == "block_density":
        return [tuple(-1.0 if x is None else float(x) for x in row) for row in dens]
    total = sum(cnt[r][s] for r in range(k) for s in range(k))
    overall = math.nan if total == 0 else float(
        sum(sums[r][s] for r in range(k) for s in range(k)) / total)
    xs = [a[i][j] for i, j in cells]
    ys = [dens[b[i]][b[j]] for i, j in cells]
    if not cells:
        return [overall, math.nan]
    mx = sum(xs, Fraction(0)) / len(xs)
    my = sum(ys, Fraction(0)) / len(ys)
    cov = sum(((x - mx) * (y - my) for x, y in zip(xs, ys)), Fraction(0))
    vx = sum(((x - mx) ** 2 for x in xs), Fraction(0))
    vy = sum(((y - my) ** 2 for y in ys), Fraction(0))
    r2 = math.nan if vx == 0 or vy == 0 else float(cov * cov / (vx * vy))
    return [overall, r2]


# ---- Classical scaling ----------------------------------------------------------

DEGENERATE = 1e-6


def classical_mds(fx, spec):
    """Returns the requested output, or None where the coordinates are not
    unique (equal eigenvalues among the first k + 1)."""
    k, what = spec["k"], spec["output"]
    n = fx.n
    if n <= k:
        return []
    d = np.array(mds_input(fx))
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ (d ** 2) @ j
    vals, vecs = np.linalg.eigh(b)
    vals, vecs = vals[::-1], vecs[:, ::-1]
    if what == "eigenvalues":
        return list(vals)
    if what == "gof":
        top = vals[:k].sum()
        return [top / np.abs(vals).sum(), top / np.maximum(vals, 0).sum()]
    scale = max(abs(vals[0]), 1e-300)
    floor = 1e-12 * abs(vals[0])
    dims = 0
    while dims < k and vals[dims] > floor:
        dims += 1
    # A dimension is unique only when its eigenvalue is separated from the
    # next one (an eigenvalue at the floor is zero up to rounding, so a
    # dimension above the floor is separated from it).
    gaps = [abs(vals[i] - vals[i + 1]) / scale for i in range(dims)
            if i + 1 < n and vals[i + 1] > floor]
    if any(g < DEGENERATE for g in gaps):
        return None
    pts = vecs[:, :dims] * np.sqrt(vals[:dims])
    for c in range(dims):
        col = pts[:, c]
        big = np.abs(col).max()
        first = next((x for x in col if abs(x) > 1e-8 * big), 0.0)
        if first < 0:
            pts[:, c] = -col
    if what == "points":
        return [tuple(row) for row in pts]
    if what == "stress":
        num = den = 0.0
        for i in range(n):
            for jj in range(i + 1, n):
                dist = float(np.sqrt(((pts[i] - pts[jj]) ** 2).sum()))
                num += (dist - d[i, jj]) ** 2
                den += dist ** 2
        return [math.sqrt(num / den) if den > 0 else math.nan]
    raise ValueError(what)


# ---- Geospatial -----------------------------------------------------------------

def great_circle(lon1, lat1, lon2, lat2, r=EARTH_RADIUS_M):
    """Vincenty's formula for a sphere (the atan2 form, accurate at every
    distance), a different formula from the engine's haversine."""
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dl = math.radians(lon2 - lon1)
    num = math.hypot(math.cos(p2) * math.sin(dl),
                     math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl))
    den = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    return r * math.atan2(num, den)


def metric(crs):
    if crs == "lonlat":
        return lambda x, y, i, j: great_circle(x[i], y[i], x[j], y[j])
    return lambda x, y, i, j: math.hypot(x[i] - x[j], y[i] - y[j])


def spatial(fx, spec):
    fn, crs = spec["fn"], spec["crs"]
    x, y = coordinates(crs, fx.n)
    dist = metric(crs)
    n = fx.n
    if fn == "edge_lengths":
        return [dist(x, y, s, t) for s, t in zip(fx.src, fx.tgt)]
    if fn == "proximity" and spec.get("k", 0) > 0:
        out = []
        for i in range(n):
            cand = sorted((dist(x, y, i, j), j) for j in range(n) if j != i)
            out += [(float(i), float(j), dd) for dd, j in cand[:spec["k"]]]
        return out
    if fn == "proximity":
        r = float(spec["r"])
        return [(float(i), float(j), dist(x, y, i, j)) for i in range(n)
                for j in range(i + 1, n) if dist(x, y, i, j) <= r]
    if fn == "group_geometry":
        labels = partition_for(spec["partition"], n)
        groups = list(dict.fromkeys(labels))
        out = []
        for gl in groups:
            mem = [i for i in range(n) if labels[i] == gl]
            if crs == "planar":
                cx = sum(x[i] for i in mem) / len(mem)
                cy = sum(y[i] for i in mem) / len(mem)
            else:
                vs = [(math.cos(math.radians(y[i])) * math.cos(math.radians(x[i])),
                       math.cos(math.radians(y[i])) * math.sin(math.radians(x[i])),
                       math.sin(math.radians(y[i]))) for i in mem]
                m = [sum(v[c] for v in vs) / len(vs) for c in range(3)]
                norm = math.sqrt(sum(c * c for c in m))
                cx = math.degrees(math.atan2(m[1], m[0]))
                cy = math.degrees(math.asin(m[2] / norm))
            xs, ys = x + [cx], y + [cy]
            ds = [dist(xs, ys, i, n) for i in mem]
            out.append((float(len(mem)), cx, cy, sum(ds) / len(ds),
                        math.sqrt(sum(t * t for t in ds) / len(ds))))
        return out
    raise ValueError(fn)


# ---- Driver ---------------------------------------------------------------------

def evaluate(fx, spec):
    fn = spec["fn"]
    if fn == "reach":
        return reach(fx, spec)
    if fn == "closeness":
        return proximity_prestige(fx, spec)
    if fn == "edge_betweenness":
        return edge_betweenness(fx, spec.get("weighted", False),
                                spec.get("weights_as") == "strength",
                                spec.get("normalized", False))
    if fn == "information_centrality":
        return information_centrality(fx, spec.get("weighted", False))
    if fn == "bridging_coefficient":
        return bridging(fx, False)
    if fn == "bridging_centrality":
        return bridging(fx, True)
    if fn == "hierarchy":
        return hierarchy(fx, spec.get("weighted", False))
    if fn == "centralization_eigenvector":
        v = eigenvector_centralization(fx, spec.get("normalized", True))
        return None if v is None else [v]
    if fn == "regular_equivalence":
        return regular_equivalence(fx, spec.get("seed"), spec.get("weighted", False))
    if fn == "automorphic_equivalence":
        return automorphic_equivalence(fx)
    if fn in ("block_density", "block_fit"):
        return block_model(fx, spec)
    if fn == "classical_mds":
        return classical_mds(fx, spec)
    if fn in ("edge_lengths", "proximity", "group_geometry"):
        return spatial(fx, spec)
    raise ValueError(f"no definition for {fn!r}")


def as_lines(values):
    out = []
    for v in values:
        if isinstance(v, tuple):
            out.append(v)
        elif v is None or (isinstance(v, float) and math.isnan(v)):
            out.append(())
        else:
            out.append(float(v))
    return out


def main() -> int:
    spec_all = manifest()
    written = skipped = 0
    for name in fixture_names():
        fx = load_fixture(name)
        if fx.n > MAX_NODES:
            continue
        for variant, spec in spec_all.items():
            if ORACLE not in spec["oracles"] or not applies(spec, fx):
                continue
            values = evaluate(fx, spec["engine"])
            if values is None:
                skipped += 1
                continue
            write_values(ORACLE, name, variant, as_lines(values))
            written += 1
    (GOLDEN / ORACLE / "VERSION").write_text("definition (direct evaluation, exact arithmetic)\n", newline="")
    print(f"Wrote {written} definition golden files for roles ({skipped} undefined).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
