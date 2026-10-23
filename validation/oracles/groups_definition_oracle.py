#!/usr/bin/env python3
"""Golden values for cohesive subsets, community partitions, and
topography, computed from the definitions.

Shares no code with the engine and, where possible, no algorithmic
strategy:

- n-cliques and n-clans (Luce 1950; Alba 1973; Mokken 1979): every set of
  nodes whose pairwise distances (breadth-first search in the whole
  network) are at most n is grown in increasing node order, and kept when
  no node outside is within n of every member; n-clans additionally need
  every pair within n inside the set.
- k-plexes (Seidman and Foster 1978): every k-plex is grown in increasing
  node order and kept when maximal and large enough.
- Edge connectivity: Menger's theorem, counting edge-disjoint paths with
  depth-first augmenting paths for every pair separately; lambda sets
  (Borgatti, Everett, and Shirey 1990) from those values.
- Factions: exhaustive search over every partition into exactly k
  non-empty groups (small networks only), with errors and the correlation
  counted directly.
- Girvan-Newman: edge betweenness from its definition,
  sum over pairs s, t of sigma_st(e) / sigma_st, with sigma counted from
  distance tables, removing the largest (ties within a relative 1e-9 go to
  the lexicographically first edge, the engine's documented rule) and
  keeping the component partition of highest exact modularity.
- Leading eigenvector (Newman 2006): numpy's symmetric eigensolver on each
  generalized modularity matrix, with the engine's documented rules. When
  an eigenvalue used for a decision is degenerate, or a component of the
  eigenvector is too close to zero for its sign to be meaningful, the
  partition is not unique and no golden file is written.
- Greedy modularity (Clauset, Newman, and Moore 2004): every connected pair
  of groups is scored in exact rational arithmetic at every step, with the
  engine's documented tie rule.
- Core-periphery (Borgatti and Everett 1999): exhaustive search over every
  core (small networks only) for the discrete model; for the continuous
  model, scipy's BFGS from many random starts maximizing the correlation
  with c_i c_j over c = t^2, and scipy's bounded L-BFGS-B from many starts
  minimizing the MINRES residual sum_{i != j} (a_ij - c_i c_j)^2 over
  c >= 0.
- Compactness (Borgatti 2006): all-pairs distances by Floyd-Warshall in
  exact arithmetic.
- Small-world inputs: transitivity by counting triangles and connected
  triples directly, mean local clustering from the same counts, and mean
  distance over connected pairs.

Results that are undefined are written as -1.

    python3 validation/oracles/groups_definition_oracle.py
"""

from __future__ import annotations

import itertools
import sys
from collections import deque
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_fixture,
                    load_manifest, write_values, write_version)  # noqa: E402

ORACLE = "definition"
MEASURES = {"n_cliques", "n_clans", "k_plexes", "lambda_connectivity", "lambda_sets",
            "factions_errors", "factions_correlation", "girvan_newman", "leading_eigen",
            "fast_greedy", "core_periphery", "core_periphery_continuous", "compactness",
            "small_world_observed"}
# Enumerations stop beyond this many sets; the fixture is then skipped.
MAX_SETS = 300_000
# Exhaustive searches are limited to this many partitions or cores.
MAX_PARTITIONS = 200_000
# Girvan-Newman from the definition costs about m^2 n^2 steps.
GN_MAX_WORK = 30_000_000


class TooMany(Exception):
    pass


def simple_adjacency(fx) -> list[set[int]]:
    adj = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            adj[a].add(b)
            adj[b].add(a)
    return adj


def bfs(adj, s, allowed=None):
    dist = {s: 0}
    q = deque([s])
    while q:
        v = q.popleft()
        for w in adj[v]:
            if w not in dist and (allowed is None or w in allowed):
                dist[w] = dist[v] + 1
                q.append(w)
    return dist


# ---- n-cliques, n-clans, k-plexes ------------------------------------------

def grow_maximal(n, compatible, extendable, min_size):
    """Grows every set S with compatible(S + [v]) in increasing order and
    reports S when no node outside can be added (extendable)."""
    out = []
    count = 0

    def grow(s):
        nonlocal count
        count += 1
        if count > MAX_SETS:
            raise TooMany
        if len(s) >= min_size and not any(extendable(s, v) for v in range(n) if v not in s):
            out.append(tuple(s))
        for v in range(s[-1] + 1, n):
            if compatible(s, v):
                grow(s + [v])

    for v in range(n):
        grow([v])
    return out


def n_cliques(fx, order, min_size, clans):
    adj = simple_adjacency(fx)
    dist = [bfs(adj, v) for v in range(fx.n)]

    def near(a, b):
        return dist[a].get(b, order + 1) <= order

    def compatible(s, v):
        return all(near(u, v) for u in s)

    found = grow_maximal(fx.n, compatible, compatible, min_size)
    if clans:
        def clan(s):
            members = set(s)
            return all(len([u for u, d in bfs(adj, v, members).items() if d <= order]) == len(s)
                       for v in s)
        found = [s for s in found if clan(s)]
    return found


def k_plexes(fx, k, min_size):
    adj = simple_adjacency(fx)

    def is_plex(s):
        return all(sum(1 for u in s if u != v and u not in adj[v]) <= k - 1 for v in s)

    def compatible(s, v):
        return is_plex(s + [v])

    return grow_maximal(fx.n, compatible, compatible, min_size)


# ---- Edge connectivity -----------------------------------------------------

def edge_disjoint_paths(adj, s, t):
    """Largest number of edge-disjoint s-t paths, by augmenting paths on the
    network with unit capacity in each direction of every edge."""
    cap = {}
    for v in range(len(adj)):
        for w in adj[v]:
            cap[(v, w)] = 1
    flow = 0
    while True:
        parent = {s: None}
        stack = [s]
        while stack and t not in parent:
            v = stack.pop()
            for w in adj[v]:
                if w not in parent and cap[(v, w)] > 0:
                    parent[w] = v
                    stack.append(w)
        if t not in parent:
            return flow
        w = t
        while parent[w] is not None:
            v = parent[w]
            cap[(v, w)] -= 1
            cap[(w, v)] += 1
            w = v
        flow += 1


def connectivity(fx):
    adj = simple_adjacency(fx)
    m = [[0] * fx.n for _ in range(fx.n)]
    for a in range(fx.n):
        for b in range(a + 1, fx.n):
            m[a][b] = m[b][a] = edge_disjoint_paths(adj, a, b)
    return m


def lambda_sets(m, n):
    top = max((max(r) for r in m), default=0)
    out = set()
    for k in range(1, top + 1):
        for v in range(n):
            cls = tuple(u for u in range(n) if u == v or m[v][u] >= k)
            # Check the class is a lambda set directly.
            inside = [m[a][b] for a in cls for b in cls if a < b]
            outside = [m[a][c] for a in cls for c in range(n) if c not in cls]
            if len(cls) > 1 and min(inside) > max(outside, default=-1):
                out.add(cls)
    return sorted(out)


# ---- Factions --------------------------------------------------------------

def tie_matrix(fx):
    a = [[0] * fx.n for _ in range(fx.n)]
    for s, t in zip(fx.src, fx.tgt):
        if s != t:
            a[s][t] = 1
            if not fx.directed:
                a[t][s] = 1
    return a


def pairs_of(fx):
    return [(i, j) for i in range(fx.n) for j in range(fx.n)
            if i != j and (fx.directed or i < j)]


def partitions(n, k):
    """Every partition of 0..n-1 into exactly k non-empty unlabelled groups,
    as label lists with labels in order of first appearance."""
    def rec(i, labels, used):
        if i == n:
            if used == k:
                yield list(labels)
            return
        if k - used > n - i:
            return
        for c in range(min(used + 1, k)):
            labels.append(c)
            yield from rec(i + 1, labels, max(used, c + 1))
            labels.pop()
    yield from rec(0, [], 0)


def correlation(x, y):
    n = len(x)
    mx, my = Fraction(sum(x), n), Fraction(sum(y), n)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx == 0 or syy == 0:
        return None
    return float(sxy) / (float(sxx) * float(syy)) ** 0.5


def stirling2(n, k):
    row = [1] + [0] * k
    for i in range(1, n + 1):
        row = [0] + [row[j - 1] + j * row[j] for j in range(1, k + 1)]
    return row[k]


def factions_best(fx, k, fit):
    if k > fx.n or stirling2(fx.n, k) > MAX_PARTITIONS:
        return None
    a = tie_matrix(fx)
    pairs = pairs_of(fx)
    x = [a[i][j] for i, j in pairs]
    best = None
    for labels in partitions(fx.n, k):
        y = [int(labels[i] == labels[j]) for i, j in pairs]
        if fit == "errors":
            v = -sum(1 for p, q in zip(x, y) if p != q)
        else:
            v = correlation(x, y)
            if v is None:
                continue
        if best is None or v > best:
            best = v
    if best is None:
        return -1.0
    return -best if fit == "errors" else best


# ---- Modularity and partitions ---------------------------------------------

def modularity_matrix(fx, weighted, simple):
    """Symmetric weighted adjacency with the modularity conventions (a loop
    of weight w adds 2w to the diagonal); `simple` collapses parallel edges
    and drops loops first."""
    n = fx.n
    a = [[Fraction(0)] * n for _ in range(n)]
    seen = set()
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if simple:
            key = (min(s, t), max(s, t))
            if s == t or key in seen:
                continue
            seen.add(key)
        w = Fraction(fx.weights[i]) if weighted else Fraction(1)
        if s == t:
            a[s][s] += 2 * w
        else:
            a[s][t] += w
            a[t][s] += w
    return a


def modularity(a, labels, gamma=Fraction(1)):
    n = len(a)
    k = [sum(r) for r in a]
    two_m = sum(k)
    if two_m == 0:
        return None
    q = Fraction(0)
    for i in range(n):
        for j in range(n):
            if labels[i] == labels[j]:
                q += a[i][j] - gamma * k[i] * k[j] / two_m
    return q / two_m


def component_labels(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    label = [-1] * n
    c = 0
    for v in range(n):
        if label[v] < 0:
            for u in bfs(adj, v):
                label[u] = c
            c += 1
    return label


def girvan_newman(fx, output):
    a = modularity_matrix(fx, False, True)
    n = fx.n
    edges = sorted({(min(s, t), max(s, t)) for s, t in zip(fx.src, fx.tgt) if s != t})
    if len(edges) ** 2 * n * n > GN_MAX_WORK:
        return None
    labels = component_labels(n, edges)
    best_q, best = modularity(a, labels), labels
    remaining = list(edges)
    while remaining:
        adj = [set() for _ in range(n)]
        for u, v in remaining:
            adj[u].add(v)
            adj[v].add(u)
        dist = [bfs(adj, s) for s in range(n)]
        sigma = [[0] * n for _ in range(n)]
        for s in range(n):
            order = sorted(dist[s], key=dist[s].get)
            sigma[s][s] = 1
            for v in order:
                for w in adj[v]:
                    if dist[s].get(w) == dist[s][v] + 1:
                        sigma[s][w] += sigma[s][v]
        eb = []
        for u, v in remaining:
            total = 0.0
            for s in range(n):
                for t in range(n):
                    if s == t or t not in dist[s]:
                        continue
                    d = dist[s][t]
                    through = 0
                    for x, y in ((u, v), (v, u)):
                        if x in dist[s] and y in dist[t] and dist[s][x] + 1 + dist[t][y] == d:
                            through += sigma[s][x] * sigma[y][t]
                    total += through / sigma[s][t]
            eb.append(total / 2)
        top = max(eb)
        cut = next(i for i, x in enumerate(eb) if x >= top - 1e-9 * abs(top))
        remaining.pop(cut)
        labels = component_labels(n, remaining)
        q = modularity(a, labels)
        if q is not None and (best_q is None or q > best_q):
            best_q, best = q, labels
    if output == "membership":
        return best
    return [-1.0 if best_q is None else float(best_q)]


class Degenerate(Exception):
    pass


def leading_eigen(fx, weighted, output):
    a = np.array([[float(x) for x in row] for row in modularity_matrix(fx, weighted, False)])
    n = fx.n
    k = a.sum(axis=1)
    two_m = k.sum()
    if two_m == 0:
        labels = list(range(n))
    else:
        labels = [0] * n
        nxt = 1
        queue = deque([list(range(n))])
        while queue:
            g = queue.popleft()
            if len(g) < 2:
                continue
            b = a[np.ix_(g, g)] - np.outer(k[g], k[g]) / two_m
            b -= np.diag(b.sum(axis=1))
            scale = np.abs(b).sum() / len(g)
            vals, vecs = np.linalg.eigh(b)
            lam = vals[-1]
            gap = vals[-1] - vals[-2]
            # B(g) always has the eigenvalue 0 (its rows sum to 0), which
            # the solver returns as roundoff of about 1e-15 * scale; the
            # engine's threshold is well clear of that, and only an
            # eigenvalue within roundoff of the threshold is ambiguous.
            if lam <= 1e-10 * scale:
                if abs(lam - 1e-10 * scale) < 1e-12 * scale:
                    raise Degenerate
                continue
            if gap < 1e-7 * scale:
                raise Degenerate
            u = vecs[:, -1]
            big = int(np.argmax(np.abs(u)))
            top = abs(u[big])
            u = u if u[big] > 0 else -u
            # Components between 1e-12 and 1e-7 of the largest could fall on
            # either side of the engine's threshold for zero.
            if np.any((np.abs(u) >= 1e-12 * top) & (np.abs(u) < 1e-7 * top)):
                raise Degenerate
            # Nodes with zero components go to the non-positive side, which
            # depends on the sign normalization; it is ambiguous when two
            # components of opposite sign have nearly the largest size.
            zero = np.abs(u) < 1e-12 * top
            rivals = (np.abs(np.abs(u) - top) < 1e-9 * top) & (u < 0)
            if zero.any() and rivals.any():
                raise Degenerate
            s = np.where(u > 1e-12 * abs(u[big]), 1.0, -1.0)
            gain = s @ b @ s
            if gain <= 1e-10 * scale * len(g):
                if abs(gain - 1e-10 * scale * len(g)) < 1e-12 * scale * len(g):
                    raise Degenerate
                continue
            plus = [g[i] for i in range(len(g)) if s[i] > 0]
            minus = [g[i] for i in range(len(g)) if s[i] < 0]
            if not plus or not minus:
                continue
            for v in minus:
                labels[v] = nxt
            nxt += 1
            queue.append(plus)
            queue.append(minus)
    if output == "membership":
        return labels
    q = modularity(modularity_matrix(fx, weighted, False), labels)
    return [-1.0 if q is None else float(q)]


def fast_greedy(fx, weighted, gamma, output):
    a = modularity_matrix(fx, weighted, False)
    n = fx.n
    k = [sum(r) for r in a]
    two_m = sum(k)
    if two_m == 0:
        return list(range(n)) if output == "membership" else [-1.0]
    groups = {v: {v} for v in range(n)}
    merges = []
    total, best_total, best_step = Fraction(0), Fraction(0), 0
    while True:
        best = None
        names = sorted(groups)
        for i, x in enumerate(names):
            for y in names[i + 1:]:
                w = sum(a[u][v] for u in groups[x] for v in groups[y])
                if w == 0:
                    continue
                kx = sum(k[u] for u in groups[x])
                ky = sum(k[v] for v in groups[y])
                score = two_m * w - gamma * kx * ky
                # Largest score; ties to the smallest (x, y).
                if best is None or score > best[0]:
                    best = (score, x, y)
        if best is None:
            break
        score, x, y = best
        total += score
        groups[x] |= groups.pop(y)
        merges.append((x, y))
        if total > best_total:
            best_total, best_step = total, len(merges)
    labels = list(range(n))
    for x, y in merges[:best_step]:
        old = labels[y]
        new = labels[x]
        labels = [new if c == old else c for c in labels]
    if output == "membership":
        return labels
    return [float(modularity(a, labels, gamma))]


# ---- Core-periphery --------------------------------------------------------

def value_matrix(fx, weighted):
    a = np.zeros((fx.n, fx.n))
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if s == t:
            continue
        w = fx.weights[i] if weighted else 1.0
        if weighted:
            a[s, t] += w
            if not fx.directed:
                a[t, s] += w
        else:
            a[s, t] = 1.0
            if not fx.directed:
                a[t, s] = 1.0
    if fx.directed:
        mask = ~np.eye(fx.n, dtype=bool)
    else:
        mask = np.triu(np.ones((fx.n, fx.n), dtype=bool), 1)
    return a, mask


def discrete_cp(fx, links, weighted):
    a, mask = value_matrix(fx, weighted)
    x_all = a[mask]
    if fx.n < 2 or np.var(x_all) == 0:
        return [-1.0]
    if 2 ** fx.n > MAX_PARTITIONS:
        return None
    best = -np.inf
    for bits in itertools.product([False, True], repeat=fx.n):
        c = np.array(bits)
        both = np.outer(c, c)
        neither = np.outer(~c, ~c)
        if links == "missing":
            keep = mask & (both | neither)
            d = both[keep].astype(float)
        elif links == "ones":
            keep = mask
            d = (~neither)[keep].astype(float)
        else:
            keep = mask
            d = both[keep].astype(float)
        x = a[keep]
        if len(x) < 2 or np.var(x) == 0 or np.var(d) == 0:
            continue
        r = np.corrcoef(x, d)[0, 1]
        best = max(best, r)
    # No core gives a pattern with both tied and untied pairs: undefined.
    return [-1.0 if best == -np.inf else float(best)]


def minres_cp(fx, weighted, output):
    """MINRES by bounded quasi-Newton from several starts; the residual, or
    the unit-norm coreness."""
    from scipy.optimize import Bounds
    a, mask = value_matrix(fx, weighted)
    x = a[mask]
    if fx.n < 2 or np.var(x) == 0:
        return [-1.0] if output == "residual" else [-1.0] * fx.n

    def objective(c):
        y = np.outer(c, c)[mask]
        r = x - y
        gm = np.zeros((fx.n, fx.n))
        gm[mask] = -2.0 * r
        return float(r @ r), (gm + gm.T) @ c

    rng = np.random.default_rng(2024)
    deg = (a.sum(axis=0) + a.sum(axis=1)) / 2.0
    starts = [deg / max(deg.max(), 1e-12)] + [rng.uniform(0.0, 1.0, fx.n) for _ in range(20)]
    best = None
    for c0 in starts:
        res = minimize(objective, c0, jac=True, method="L-BFGS-B",
                       bounds=Bounds(0.0, np.inf),
                       options={"ftol": 1e-15, "gtol": 1e-12, "maxiter": 20000})
        if best is None or res.fun < best.fun:
            best = res
    # Polish the best point with exact coordinate minimization, c_i =
    # max(0, sum_j (a_ij + a_ji) c_j / (k sum_j c_j^2)) over j != i, with k = 2
    # for directed and 1 for undirected networks, which makes the result independent of where the
    # quasi-Newton iterations stopped (they differ across platforms on flat
    # objectives such as path10's).
    c = best.x.copy()
    full = np.zeros((fx.n, fx.n))
    full[mask] = x
    for _ in range(200000):
        change = 0.0
        for i in range(fx.n):
            others = np.arange(fx.n) != i
            num = float(full[i, others] @ c[others] + full[others, i] @ c[others])
            # Each ordered pair is a term when directed, each unordered
            # pair when undirected.
            den = (2.0 if fx.directed else 1.0) * float(c[others] @ c[others])
            new = max(0.0, num / den) if den > 0 else 0.0
            change = max(change, abs(new - c[i]))
            c[i] = new
        if change < 1e-15:
            break
    if output == "residual":
        return [float(objective(c)[0])]
    norm = np.linalg.norm(c)
    return list(c / norm) if norm > 0 else [-1.0] * fx.n


def continuous_cp(fx, weighted):
    a, mask = value_matrix(fx, weighted)
    x = a[mask]
    if fx.n < 2 or np.var(x) == 0:
        return [-1.0]
    xc = x - x.mean()

    def neg_r(t):
        c = t * t
        y = np.outer(c, c)[mask]
        yc = y - y.mean()
        den = np.sqrt((xc @ xc) * (yc @ yc))
        if den == 0:
            return 1.0, np.zeros_like(t)
        r = (xc @ yc) / den
        # dr/dy over the pairs, then dr/dc = (G + G') c and dr/dt = 2 t dr/dc.
        gy = xc / den - r * yc / (yc @ yc)
        gm = np.zeros((fx.n, fx.n))
        gm[mask] = gy
        dc = (gm + gm.T) @ c
        return -r, -2.0 * t * dc

    rng = np.random.default_rng(12345)
    deg = a.sum(axis=0) + a.sum(axis=1) + 1.0
    starts = [np.sqrt(deg / np.linalg.norm(deg))]
    starts += [rng.uniform(0.05, 1.0, fx.n) for _ in range(30)]
    best = -np.inf
    for t0 in starts:
        # BFGS stops when its line search can make no progress, which a
        # stale Hessian approximation can cause well before the maximum;
        # restarting from the stopping point (a fresh approximation) until
        # nothing changes finds the maximum to full precision.
        value = np.inf
        while True:
            res = minimize(neg_r, t0, jac=True, method="BFGS",
                           options={"gtol": 1e-11, "maxiter": 5000})
            if res.fun >= value - 1e-13:
                break
            value, t0 = res.fun, res.x
        best = max(best, -value)
    return [float(best)]


# ---- Compactness and small worlds ------------------------------------------

def compactness(fx, weighted):
    n = fx.n
    if n < 2:
        return [-1.0]
    inf = None
    d = [[Fraction(0) if i == j else inf for j in range(n)] for i in range(n)]
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if s == t:
            continue
        w = Fraction(fx.weights[i]) if weighted else Fraction(1)
        for x, y in ((s, t),) if fx.directed else ((s, t), (t, s)):
            if d[x][y] is None or w < d[x][y]:
                d[x][y] = w
    for m in range(n):
        for i in range(n):
            if d[i][m] is None:
                continue
            for j in range(n):
                if d[m][j] is not None:
                    alt = d[i][m] + d[m][j]
                    if d[i][j] is None or alt < d[i][j]:
                        d[i][j] = alt
    total = sum(1 / d[i][j] for i in range(n) for j in range(n)
                if i != j and d[i][j] is not None)
    return [float(Fraction(total) / (n * (n - 1)))]


def small_world_observed(fx, local):
    adj = simple_adjacency(fx)
    n = fx.n
    tri = [sum(1 for u, v in itertools.combinations(sorted(adj[x]), 2) if v in adj[u])
           for x in range(n)]
    triples = [len(adj[x]) * (len(adj[x]) - 1) // 2 for x in range(n)]
    if local:
        c = sum(Fraction(t, p) if p else Fraction(0) for t, p in zip(tri, triples)) / n \
            if n else None
    else:
        c = Fraction(sum(tri), sum(triples)) if sum(triples) else None
    total, pairs = 0, 0
    for v in range(n):
        for u, dd in bfs(adj, v).items():
            if u != v:
                total += dd
                pairs += 1
    return [-1.0 if c is None else float(c), total / pairs if pairs else -1.0]


def compute(spec: dict, fx):
    e = spec["engine"]
    fn = e["fn"]
    if fn in ("n_cliques", "n_clans"):
        return n_cliques(fx, e["n"], e.get("min_size", 0), fn == "n_clans")
    if fn == "k_plexes":
        return k_plexes(fx, e["k"], e["min_size"])
    if fn == "lambda_connectivity":
        m = connectivity(fx)
        return [m[i][j] for i in range(fx.n) for j in range(i + 1, fx.n)]
    if fn == "lambda_sets":
        return lambda_sets(connectivity(fx), fx.n)
    if fn in ("factions_errors", "factions_correlation"):
        v = factions_best(fx, e["k"], fn.split("_")[1])
        return None if v is None else [v]
    if fn == "girvan_newman":
        return girvan_newman(fx, e.get("output", "modularity"))
    if fn == "leading_eigen":
        try:
            return leading_eigen(fx, e.get("weighted", False), e.get("output", "modularity"))
        except Degenerate:
            return None
    if fn == "fast_greedy":
        return fast_greedy(fx, e.get("weighted", False), Fraction(e.get("resolution", 1.0)),
                           e.get("output", "modularity"))
    if fn == "core_periphery":
        return discrete_cp(fx, e.get("links", "missing"), e.get("weighted", False))
    if fn == "core_periphery_continuous":
        if e.get("fit", "minres") == "minres":
            return minres_cp(fx, e.get("weighted", False), e.get("output", "residual"))
        return continuous_cp(fx, e.get("weighted", False))
    if fn == "compactness":
        return compactness(fx, e.get("weighted", False))
    if fn == "small_world_observed":
        return small_world_observed(fx, e.get("local", False))
    return None


def main() -> None:
    manifest = load_manifest()
    count = 0
    skipped = []
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, spec in manifest.items():
            if spec["engine"]["fn"] not in MEASURES:
                continue
            if ORACLE not in spec["oracles"] or not applies(spec, fx):
                continue
            if name in spec["oracles"][ORACLE].get("undefined_on", []):
                continue
            try:
                values = compute(spec, fx)
            except TooMany:
                values = None
            if values is None:
                skipped.append((variant, name))
                continue
            write_values(ORACLE, name, variant, values)
            count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (groups): wrote {count} golden files")
    for variant, name in skipped:
        print(f"  no answer: {variant} on {name}")


if __name__ == "__main__":
    main()
