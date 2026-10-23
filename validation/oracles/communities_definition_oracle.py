#!/usr/bin/env python3
"""Golden values for community detection and partition comparison,
computed from the definitions.

Shares no code with the engine and, where possible, no algorithmic
strategy:

- Walktrap (Pons and Latapy 2006): the walk distributions P^t_i. are
  computed in exact rational arithmetic, r^2 = sum_k (P_ik - P_jk)^2 / d_k
  (which needs no square roots), and at every step every adjacent pair of
  groups is scored by Ward's delta sigma exactly. Exact ties go to the
  pair with the smallest group names, the engine's documented rule; the
  fixtures with ties are reported. Modularity of each level is exact.
- Voronoi communities (Deritei et al. 2014; Molnar et al. 2024): edge
  clustering coefficients by counting common neighbors, local relative
  densities by counting edges, all-pairs distances by Floyd-Warshall in
  exact arithmetic, and every candidate radius evaluated with exact
  comparisons and the engine's documented tie rules.
- Map equation (Rosvall and Bergstrom 2008): visit rates in closed form
  (undirected) or by solving the PageRank linear system with numpy
  (directed), then the two-level codelength term by term.
- Exact maximum modularity: every partition (restricted growth strings)
  for n <= 10; otherwise the integer linear program of Brandes et al.
  (2008), solved by scipy's MILP interface (HiGHS) with a zero optimality
  gap, and the modularity of the partition it returns evaluated exactly.
- Leiden quality, spinglass energy, the map equation of the engine's
  Infomap partition, the connectedness of Leiden groups, and the label
  propagation fixed-point property are evaluated on the engine's
  partitions (seed 1), as for Louvain elsewhere in the suite.
- Partition comparison (Meila 2007; Danon et al. 2005; van Dongen 2000;
  Rand 1971; Hubert and Arabie 1985; Vinh, Epps, and Bailey 2010): from
  the contingency table, with the expected mutual information summed
  term by term with log-gamma functions.

Results that are undefined are written as -1.

    python3 validation/oracles/communities_definition_oracle.py
"""

from __future__ import annotations

import heapq
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_fixture,  # noqa: E402
                    load_manifest, write_values, write_version)

ORACLE = "definition"
MEASURES = {"walktrap", "voronoi", "map_equation", "infomap", "optimal_modularity", "leiden",
            "spinglass", "label_propagation", "compare_partitions", "fluid_communities"}
# Partitions enumerated exhaustively up to this many nodes; the integer
# program beyond, up to MILP_MAX_NODES.
ENUM_MAX_NODES = 10
MILP_MAX_NODES = 36


# Fixtures on which Walktrap meets exact ties (resolved by the engine's
# documented rule, smallest group names first), reported at the end.
TIED: set[str] = set()


def blocks(n):
    return [3 * v // n for v in range(n)]


def stripes(n):
    return [(7 * v) % 4 for v in range(n)]


def first_appearance(labels):
    seen = {}
    return [seen.setdefault(c, len(seen)) for c in labels]


def weight_of(fx, i, weighted):
    return Fraction(fx.weights[i]) if weighted else Fraction(1)


# ---- Modularity ------------------------------------------------------------

def modularity(fx, membership, weighted, gamma=Fraction(1)):
    """Exact modularity (Newman and Girvan 2004; directed form of Leicht
    and Newman 2008): a loop adds 2w to A_ii in an undirected network and
    w in a directed one. None without edges."""
    n = fx.n
    a = {}
    k_out = [Fraction(0)] * n
    k_in = [Fraction(0)] * n
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        w = weight_of(fx, i, weighted)
        if fx.directed:
            a[(s, t)] = a.get((s, t), 0) + w
            k_out[s] += w
            k_in[t] += w
        else:
            if s == t:
                a[(s, s)] = a.get((s, s), 0) + 2 * w
            else:
                a[(s, t)] = a.get((s, t), 0) + w
                a[(t, s)] = a.get((t, s), 0) + w
            k_out[s] += w
            k_out[t] += w
    if not fx.directed:
        k_in = k_out
    total = sum(k_out)
    if total == 0:
        return None
    q = sum(w for (i, j), w in a.items() if membership[i] == membership[j])
    groups = {}
    for v in range(n):
        o, i_ = groups.get(membership[v], (Fraction(0), Fraction(0)))
        groups[membership[v]] = (o + k_out[v], i_ + k_in[v])
    q -= gamma * sum(o * i_ for o, i_ in groups.values()) / total
    return q / total


def flt(q):
    return -1.0 if q is None else float(q)


# ---- Walktrap ----------------------------------------------------------------

def walktrap(fx, steps, weighted):
    """Returns (merges, path of exact modularities). Raises Tie."""
    n = fx.n
    rows = [dict() for _ in range(n)]
    other = [[Fraction(0), 0] for _ in range(n)]
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        w = weight_of(fx, i, weighted)
        rows[s][t] = rows[s].get(t, 0) + w
        if s != t:
            rows[t][s] = rows[t].get(s, 0) + w
            for v in (s, t):
                other[v][0] += w
                other[v][1] += 1
    for v in range(n):
        total, count = other[v]
        rows[v][v] = rows[v].get(v, 0) + (total / count if count else Fraction(1))
    deg = [sum(r.values()) for r in rows]

    def walk(i):
        x = {i: Fraction(1)}
        for _ in range(steps):
            y = {}
            for j, xj in x.items():
                for k, a in rows[j].items():
                    y[k] = y.get(k, 0) + xj * a / deg[j]
            x = y
        return x

    vec = {i: walk(i) for i in range(n)}
    size = {i: 1 for i in range(n)}
    nbrs = {i: {j for j in rows[i] if j != i} for i in range(n)}

    def delta(a, b):
        va, vb = vec[a], vec[b]
        r2 = sum((va.get(k, 0) - vb.get(k, 0)) ** 2 / deg[k] for k in set(va) | set(vb))
        return Fraction(size[a] * size[b], size[a] + size[b]) * r2 / n

    cache = {}
    for a in range(n):
        for b in nbrs[a]:
            if a < b:
                cache[(a, b)] = delta(a, b)
    merges = []
    while cache:
        best = min(cache.values())
        pairs = sorted(p for p, d in cache.items() if d == best)
        if len(pairs) > 1:
            TIED.add(fx.name)
        a, b = pairs[0]
        new = n + len(merges)
        merges.append((a, b))
        sa, sb = size[a], size[b]
        keys = set(vec[a]) | set(vec[b])
        vec[new] = {k: (sa * vec[a].get(k, 0) + sb * vec[b].get(k, 0)) / (sa + sb) for k in keys}
        size[new] = sa + sb
        around = (nbrs[a] | nbrs[b]) - {a, b}
        for p in [p for p in cache if a in p or b in p]:
            del cache[p]
        for x in around:
            nbrs[x] = (nbrs[x] - {a, b}) | {new}
            cache[(x, new)] = delta(x, new)
        nbrs[new] = around
        del vec[a], vec[b]
    return merges


def cut(n, merges, k):
    """Membership after the first k merges, labelled by first appearance."""
    members = {i: [i] for i in range(n)}
    for i, (a, b) in enumerate(merges[:k]):
        members[n + i] = members.pop(a) + members.pop(b)
    lab = [0] * n
    for ms in members.values():
        for v in ms:
            lab[v] = min(ms)
    return first_appearance(lab)


def walktrap_values(fx, e):
    merges = walktrap(fx, e.get("steps", 4), e.get("weighted", False))
    output = e.get("output", "modularity")
    if output == "merges":
        return [tuple(m) for m in merges]
    if output == "cut":
        k = e["groups"]
        if not 1 <= k <= fx.n or fx.n - k > len(merges):
            return [-1.0] * fx.n
        return cut(fx.n, merges, fx.n - k)
    path = [modularity(fx, cut(fx.n, merges, s), e.get("weighted", False))
            for s in range(len(merges) + 1)]
    if path[0] is None:
        best = 0
    else:
        best = max(range(len(path)), key=lambda s: (path[s], -s))
    if output == "membership":
        return cut(fx.n, merges, best)
    return [flt(path[best])]


# ---- Voronoi -----------------------------------------------------------------

def simple_neighbors(fx):
    nb = [set() for _ in range(fx.n)]
    for s, t in zip(fx.src, fx.tgt):
        if s != t:
            nb[s].add(t)
            nb[t].add(s)
    return nb


def voronoi(fx, e):
    n = fx.n
    weighted = e.get("weighted", False)
    mode = e.get("mode", "out")
    nb = simple_neighbors(fx)
    deg = [len(x) for x in nb]
    strength = [Fraction(0)] * n
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if s != t:
            w = weight_of(fx, i, weighted)
            strength[s] += w
            strength[t] += w
    density = []
    for v in range(n):
        closed = nb[v] | {v}
        inside = sum(1 for a in closed for b in nb[a] if b in closed and a < b)
        boundary = sum(1 for a in closed for b in nb[a] if b not in closed)
        density.append(Fraction(0) if inside + boundary == 0
                       else strength[v] * Fraction(inside, inside + boundary))
    order = sorted(range(n), key=lambda v: (-density[v], v))
    inf = None
    d = [[inf] * n for _ in range(n)]
    for v in range(n):
        d[v][v] = Fraction(0)

    def relax(a, b, ln):
        if d[a][b] is None or ln < d[a][b]:
            d[a][b] = ln

    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if s == t:
            continue
        base = Fraction(1.0 / fx.weights[i]) if e.get("lengths") == "inverse" else Fraction(1)
        lo = min(deg[s], deg[t])
        ln = Fraction(0) if lo <= 1 else base * Fraction(lo - 1, len(nb[s] & nb[t]) + 1)
        if not fx.directed or mode == "all":
            relax(s, t, ln)
            relax(t, s, ln)
        elif mode == "out":
            relax(s, t, ln)
        else:
            relax(t, s, ln)
    for k in range(n):
        dk = d[k]
        for i in range(n):
            dik = d[i][k]
            if dik is None:
                continue
            di = d[i]
            for j in range(n):
                if dk[j] is not None and (di[j] is None or dik + dk[j] < di[j]):
                    di[j] = dik + dk[j]

    def evaluate(radius):
        excluded = [False] * n
        gens = []
        for v in order:
            if excluded[v]:
                continue
            gens.append(v)
            excluded[v] = True
            for u in range(n):
                if d[v][u] is not None and (radius is None or d[v][u] < radius):
                    excluded[u] = True
        cells = []
        for u in range(n):
            best = None
            for r, gnode in enumerate(gens):
                du = d[gnode][u]
                if du is not None and (best is None or du < best[0]):
                    best = (du, r)
            cells.append(best[1])
        memb = first_appearance(cells)
        return memb, modularity(fx, memb, weighted)

    if e.get("radius", -1.0) >= 0:
        radius = Fraction(e["radius"])
        memb, q = evaluate(radius)
        return memb, q, radius
    values = sorted({d[i][j] for i in range(n) for j in range(n)
                     if i != j and d[i][j] is not None})
    if len(values) > 1000:
        k = len(values)
        values = [values[i * (k - 1) // 999] for i in range(1000)]
    values.append(None)  # infinity
    best = None
    for r in values:
        memb, q = evaluate(r)
        if best is None or (q is not None and best[1] is not None and q > best[1]):
            best = (memb, q, r)
    return best


def voronoi_values(fx, e):
    memb, q, r = voronoi(fx, e)
    if e.get("output", "membership") == "modularity":
        return [flt(q)]
    return memb


# ---- Map equation ------------------------------------------------------------

def plogp(x):
    return x * math.log2(x) if x > 0 else 0.0


def map_equation(fx, membership, weighted, tau=0.15):
    n = fx.n
    w = [fx.weights[i] if weighted else 1.0 for i in range(len(fx.src))]
    if not fx.directed:
        s = [0.0] * n
        exit_w = {}
        total = 0.0
        for (a, b), x in zip(zip(fx.src, fx.tgt), w):
            s[a] += x
            s[b] += x
            total += 2 * x
            if membership[a] != membership[b]:
                exit_w[membership[a]] = exit_w.get(membership[a], 0.0) + x
                exit_w[membership[b]] = exit_w.get(membership[b], 0.0) + x
        if total == 0:
            return None
        p = [x / total for x in s]
        q = {c: exit_w.get(c, 0.0) / total for c in set(membership)}
    else:
        out_w = [0.0] * n
        for (a, _), x in zip(zip(fx.src, fx.tgt), w):
            out_w[a] += x
        # Transition matrix T[a][b] (rows sum to 1) and stationary p = p T.
        t = np.zeros((n, n))
        for (a, b), x in zip(zip(fx.src, fx.tgt), w):
            if out_w[a] > 0:
                t[a, b] += (1 - tau) * x / out_w[a]
        for a in range(n):
            ta = tau if out_w[a] > 0 else 1.0
            t[a, :] += ta / n
        m = t.T - np.eye(n)
        m[-1, :] = 1.0
        rhs = np.zeros(n)
        rhs[-1] = 1.0
        p = list(np.linalg.solve(m, rhs))
        size = {}
        for v in range(n):
            size[membership[v]] = size.get(membership[v], 0) + 1
        q = {c: 0.0 for c in size}
        leave = [dict() for _ in range(n)]
        for (a, b), x in zip(zip(fx.src, fx.tgt), w):
            if membership[a] != membership[b]:
                leave[a][membership[b]] = leave[a].get(membership[b], 0.0) + x
        for a in range(n):
            c = membership[a]
            ta = tau if out_w[a] > 0 else 1.0
            outside = sum(leave[a].values())
            link = (1 - ta) * outside / out_w[a] if out_w[a] > 0 else 0.0
            q[c] += p[a] * (ta * (1 - size[c] / n) + link)
    flow = {}
    for v in range(n):
        flow[membership[v]] = flow.get(membership[v], 0.0) + p[v]
    qt = sum(q.values())
    return (plogp(qt) - 2 * sum(plogp(x) for x in q.values()) - sum(plogp(x) for x in p)
            + sum(plogp(q[c] + flow[c]) for c in q))


# ---- Engine partitions ------------------------------------------------------

def engine_graph(fx):
    from engine import Graph, _load_bindings
    _load_bindings()
    return Graph(fx.n, fx.src, fx.tgt, fx.weights, fx.directed), sys.modules["bindings_communities"]


def leiden_quality(fx, membership, objective, gamma, weighted):
    n = fx.n
    a_in = {}
    k = [Fraction(0)] * n
    two_m = Fraction(0)
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        w = weight_of(fx, i, weighted)
        k[s] += w
        k[t] += w
        two_m += 2 * w
        if membership[s] == membership[t]:
            a_in[membership[s]] = a_in.get(membership[s], 0) + 2 * w
    if two_m == 0:
        return None
    node = k if objective == "modularity" else [Fraction(1)] * n
    g = gamma / two_m if objective == "modularity" else gamma
    tot = {}
    for v in range(n):
        tot[membership[v]] = tot.get(membership[v], 0) + node[v]
    return (sum(a_in.values()) - g * sum(x * x for x in tot.values())) / two_m


def pieces(fx, membership):
    """Total number of connected pieces of the groups' induced subgraphs."""
    nb = simple_neighbors(fx)
    seen = [False] * fx.n
    count = 0
    for s in range(fx.n):
        if seen[s]:
            continue
        count += 1
        seen[s] = True
        stack = [s]
        while stack:
            v = stack.pop()
            for u in nb[v]:
                if not seen[u] and membership[u] == membership[v]:
                    seen[u] = True
                    stack.append(u)
    return count


def spinglass_energy(fx, membership, weighted, gamma, null_model):
    n = fx.n
    a = {}
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if s != t:
            key = (min(s, t), max(s, t))
            a[key] = a.get(key, 0) + weight_of(fx, i, weighted)
    k = [Fraction(0)] * n
    for (s, t), w in a.items():
        k[s] += w
        k[t] += w
    m = sum(a.values())
    h = Fraction(0)
    for i in range(n):
        for j in range(i + 1, n):
            if membership[i] != membership[j]:
                continue
            if m == 0:
                p = Fraction(0)
            elif null_model == "config":
                p = k[i] * k[j] / (2 * m)
            else:
                p = 2 * m / (n * (n - 1))
            h -= a.get((i, j), 0) - gamma * p
    return h


def dominant_fraction(fx, membership, weighted):
    n = fx.n
    tot = [dict() for _ in range(n)]
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        if s == t:
            continue
        w = weight_of(fx, i, weighted)
        tot[s][membership[t]] = tot[s].get(membership[t], 0) + w
        tot[t][membership[s]] = tot[t].get(membership[s], 0) + w
    good = sum(1 for v in range(n)
               if not tot[v] or tot[v].get(membership[v], 0) == max(tot[v].values()))
    return good / n if n else -1.0


# ---- Exact maximum modularity ------------------------------------------------

def restricted_growth(n):
    a = [0] * n

    def rec(i, top):
        if i == n:
            yield a
            return
        for c in range(top + 2):
            a[i] = c
            yield from rec(i + 1, max(top, c))

    if n == 0:
        yield []
    else:
        yield from rec(1, 0)


def optimal_modularity(fx, weighted, gamma):
    if fx.n <= ENUM_MAX_NODES:
        best = None
        for p in restricted_growth(fx.n):
            q = modularity(fx, p, weighted, gamma)
            if q is None:
                return None
            if best is None or q > best:
                best = q
        return best
    if fx.n > MILP_MAX_NODES:
        return "skip"
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import lil_matrix
    n = fx.n
    a = np.zeros((n, n))
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        w = fx.weights[i] if weighted else 1.0
        if fx.directed:
            a[s, t] += w
        elif s == t:
            a[s, s] += 2 * w
        else:
            a[s, t] += w
            a[t, s] += w
    kout, kin = a.sum(axis=1), a.sum(axis=0)
    total = kout.sum()
    if total == 0:
        return None
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    index = {p: x for x, p in enumerate(pairs)}
    c = np.array([-(a[i, j] + a[j, i] - float(gamma) * (kout[i] * kin[j] + kout[j] * kin[i])
                    / total) for i, j in pairs])
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                ij, jk, ik = index[(i, j)], index[(j, k)], index[(i, k)]
                rows += [((ij, 1), (jk, 1), (ik, -1)), ((ij, 1), (jk, -1), (ik, 1)),
                         ((ij, -1), (jk, 1), (ik, 1))]
    mat = lil_matrix((len(rows), len(pairs)))
    for r, terms in enumerate(rows):
        for col, v in terms:
            mat[r, col] = v
    res = milp(c, integrality=np.ones(len(pairs)), bounds=Bounds(0, 1),
               constraints=LinearConstraint(mat.tocsr(), -np.inf, 1),
               options={"mip_rel_gap": 0.0, "time_limit": 600})
    if not res.success:
        return "skip"
    x = np.round(res.x)
    label = list(range(n))

    def find(v):
        while label[v] != v:
            label[v] = label[label[v]]
            v = label[v]
        return v

    for (i, j), xv in zip(pairs, x):
        if xv > 0.5:
            label[find(i)] = find(j)
    return modularity(fx, [find(v) for v in range(n)], weighted, gamma)


# ---- Partition comparison ------------------------------------------------------

def compare(x, y, method):
    n = len(x)
    if n == 0:
        return None
    cells, a, b = {}, {}, {}
    for u, v in zip(x, y):
        cells[(u, v)] = cells.get((u, v), 0) + 1
        a[u] = a.get(u, 0) + 1
        b[v] = b.get(v, 0) + 1

    def h(sizes):
        return -sum(s / n * math.log(s / n) for s in sizes.values())

    mi = sum(c / n * math.log(n * c / (a[u] * b[v])) for (u, v), c in cells.items())
    if method == "vi":
        return max(h(a) + h(b) - 2 * mi, 0.0)
    if method == "nmi":
        return 1.0 if h(a) + h(b) == 0 else 2 * mi / (h(a) + h(b))
    if method == "split_join":
        rows, cols = {}, {}
        for (u, v), c in cells.items():
            rows[u] = max(rows.get(u, 0), c)
            cols[v] = max(cols.get(v, 0), c)
        return 2 * n - sum(rows.values()) - sum(cols.values())

    def c2(k):
        return k * (k - 1) // 2

    if method in ("rand", "adjusted_rand"):
        if n < 2:
            return None
        s = sum(c2(c) for c in cells.values())
        sa, sb, total = sum(c2(k) for k in a.values()), sum(c2(k) for k in b.values()), c2(n)
        if method == "rand":
            return (total + 2 * s - sa - sb) / total
        expected = Fraction(sa * sb, total)
        denom = Fraction(sa + sb, 2) - expected
        return 1.0 if denom == 0 else float((s - expected) / denom)
    if method == "ami":
        lg = math.lgamma
        emi = 0.0
        for ai in a.values():
            for bj in b.values():
                for nij in range(max(1, ai + bj - n), min(ai, bj) + 1):
                    logp = (lg(ai + 1) + lg(bj + 1) + lg(n - ai + 1) + lg(n - bj + 1)
                            - lg(n + 1) - lg(nij + 1) - lg(ai - nij + 1) - lg(bj - nij + 1)
                            - lg(n - ai - bj + nij + 1))
                    emi += nij / n * math.log(n * nij / (ai * bj)) * math.exp(logp)
        mean_h = (h(a) + h(b)) / 2
        both_trivial = (len(a) == 1 and len(b) == 1) or (len(a) == n and len(b) == n)
        if both_trivial:
            return 1.0
        denom = mean_h - emi
        if denom == 0:
            return None
        return (mi - emi) / denom
    raise ValueError(method)


# ---- Dispatch ----------------------------------------------------------------

def connected(fx):
    nb = simple_neighbors(fx)
    if fx.n == 0:
        return False
    seen = {0}
    stack = [0]
    while stack:
        v = stack.pop()
        for u in nb[v] - seen:
            seen.add(u)
            stack.append(u)
    return len(seen) == fx.n


def compute(spec: dict, fx):
    e = spec["engine"]
    fn = e["fn"]
    weighted = e.get("weighted", False)
    if fn == "walktrap":
        return walktrap_values(fx, e)
    if fn == "voronoi":
        return voronoi_values(fx, e)
    if fn == "map_equation":
        v = map_equation(fx, blocks(fx.n), weighted)
        return [-1.0 if v is None else v]
    if fn == "compare_partitions":
        v = compare(blocks(fx.n), stripes(fx.n), e["method"])
        return [-1.0 if v is None else float(v)]
    if fn == "optimal_modularity":
        if fx.name in e.get("skip", []):
            return [-1.0]
        v = optimal_modularity(fx, weighted, Fraction(e.get("resolution", 1.0)))
        if v == "skip":
            return None
        return [flt(v)]
    if fn == "fluid_communities":
        ok = connected(fx) and fx.n >= 2 and e["k"] <= fx.n
        return [e["k"] if ok else -1.0]
    g, binding = engine_graph(fx)
    if fn == "infomap":
        memb, _ = binding.infomap(g, weighted, e.get("trials", 10), e.get("seed", 1))
        v = map_equation(fx, memb, weighted)
        return [-1.0 if v is None else v]
    if fn == "leiden":
        objective = e.get("objective", "modularity")
        memb, _ = binding.leiden(g, objective, e.get("resolution", 1.0), e.get("beta", 0.01),
                                 e.get("n_iterations", -1), weighted, e.get("seed", 1))
        if e.get("output") == "groups":
            return [pieces(fx, memb)]
        return [flt(leiden_quality(fx, memb, objective, Fraction(e.get("resolution", 1.0)),
                                   weighted))]
    if fn == "spinglass":
        if not connected(fx):
            return [-1.0]
        memb, _, _ = binding.spinglass(g, e.get("spins", 25), weighted, e.get("gamma", 1.0),
                                       e.get("null_model", "config"), e.get("seed", 1))
        return [float(spinglass_energy(fx, memb, weighted, Fraction(e.get("gamma", 1.0)),
                                       e.get("null_model", "config")))]
    if fn == "label_propagation" and e.get("output") == "fixed_point":
        memb = binding.label_propagation(g, e.get("variant", "dominance"), e.get("mode", "all"),
                                         weighted, e.get("seed", 1))
        return [dominant_fraction(fx, memb, weighted)]
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
            values = compute(spec, fx)
            if values is None:
                skipped.append((variant, name))
                continue
            write_values(ORACLE, name, variant, values)
            count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (communities): wrote {count} golden files")
    for variant, name in skipped:
        print(f"  no answer: {variant} on {name}")
    if TIED:
        print("  Walktrap merge ties (engine rule applied): " + ", ".join(sorted(TIED)))


if __name__ == "__main__":
    main()
