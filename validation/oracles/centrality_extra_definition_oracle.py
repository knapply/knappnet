#!/usr/bin/env python3
"""Golden values for the centrality variants computed directly from their
definitions.

This oracle shares no code or algorithmic strategy with the engine:
- shortest-path betweenness (for sources and targets, sampled, with
  endpoints, of edges), percolation centrality: exact distances by the
  Floyd-Warshall recurrence in rational arithmetic, path counts sigma by
  the recurrence over predecessors, and for every pair (s, t) the share
  sigma_sv sigma_vt / sigma_st of each node (sigma_su sigma_vt / sigma_st
  of each arc) on a shortest path, summed pair by pair (the engine uses
  Brandes' accumulation per source);
- load centrality: for every pair (s, t) separately, a unit of flow
  pushed back from t to s and split equally among predecessors, in
  rational arithmetic (the engine accumulates over all targets of a
  source at once);
- current-flow measures, resistance distance: the Moore-Penrose inverse of
  the Laplacian (numpy), the current of every edge for every pair, and the
  throughput summed pair by pair (the engine sorts flow-matrix rows of a
  Cholesky-based grounded inverse);
- Kemeny's constant: the trace of the fundamental matrix minus 1 (the
  engine uses eigenvalues); second-order centrality: one linear solve per
  node for the mean first passage times, as published (the engine uses the
  fundamental matrix);
- communicability: scipy's matrix exponential (the engine uses an
  eigendecomposition per component);
- Laplacian centrality: the Laplacian energy before and after deleting each
  node, from explicit matrices (the engine uses a closed form);
- trophic levels: Gaussian elimination in rational arithmetic;
- reaching centrality (unweighted): reachability and exact distances;
- dispersion (undirected simple networks): the definition over pairs of
  common neighbors.

Node sets, states, and sampled sources and pairs repeat
validation/bindings/centrality_extra.py. NaN is written as -999 and
infinity as -1000.

    python3 validation/oracles/centrality_extra_definition_oracle.py
"""

from __future__ import annotations

import itertools
import math
import random
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy.linalg

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
MAX_NODES = 80
UNDEFINED = -999.0
INFINITE = -1000.0


def mark(x):
    if x is None:
        return UNDEFINED
    x = float(x)
    if math.isnan(x):
        return UNDEFINED
    if math.isinf(x):
        return INFINITE
    return x


def node_subset(n, kind):
    if kind == "all":
        return list(range(n))
    if kind == "even":
        return [v for v in range(n) if v % 2 == 0]
    if kind == "thirds":
        return [v for v in range(n) if v % 3 == 0]
    raise ValueError(kind)


def states(n):
    return [((37 * v) % 11) / 10 for v in range(n)]


def sampled_sources(n, k, seed):
    return random.Random(seed).sample(list(range(n)), k)


def approximation_pairs(fx, seed, epsilon):
    """The pairs NetworkX's approximation draws (see the runner)."""
    import networkx as nx
    from networkx.utils import reverse_cuthill_mckee_ordering

    n = fx.n
    if n < 3:
        return []
    g = nx.Graph()
    g.add_nodes_from(range(n))
    g.add_edges_from((a, b) for a, b in zip(fx.src, fx.tgt) if a != b)
    ordering = list(reverse_cuthill_mckee_ordering(g))
    cstar = n * (n - 1) / ((n - 1.0) * (n - 2.0))
    k = int(math.ceil((cstar / epsilon) ** 2 * math.log(n)))
    rng = random.Random(seed)
    out = []
    for _ in range(k):
        s, t = rng.sample(range(n), 2)
        out.append((ordering[s], ordering[t]))
    return out


# ---- Shortest paths, counted pair by pair -------------------------------------


class Paths:
    """Exact all-pairs distances (Floyd-Warshall in rationals), the shortest
    length and number of shortest parallel arcs of every ordered pair, and
    path counts."""

    def __init__(self, fx, weighted, distinct):
        n = fx.n
        self.n = n
        lens = [Fraction(fx.weights[e]) if weighted else Fraction(1) for e in range(len(fx.src))]
        self.edge_len = lens
        arcs = {}
        for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
            pairs = [(a, b)] if fx.directed else [(a, b), (b, a)]
            if a == b:
                continue
            for p in set(pairs):
                L = lens[e]
                if p not in arcs or L < arcs[p][0]:
                    arcs[p] = [L, 1]
                elif L == arcs[p][0]:
                    arcs[p][1] += 1
        if not distinct:
            for p in arcs:
                arcs[p][1] = 1
        self.arcs = arcs
        d = [[None] * n for _ in range(n)]
        for v in range(n):
            d[v][v] = Fraction(0)
        for (a, b), (L, _) in arcs.items():
            if d[a][b] is None or L < d[a][b]:
                d[a][b] = L
        for k in range(n):
            for i in range(n):
                if d[i][k] is None:
                    continue
                for j in range(n):
                    if d[k][j] is None:
                        continue
                    alt = d[i][k] + d[k][j]
                    if d[i][j] is None or alt < d[i][j]:
                        d[i][j] = alt
        self.d = d
        self.sigma = []
        for s in range(n):
            order = sorted((v for v in range(n) if d[s][v] is not None), key=lambda v: d[s][v])
            sig = [0] * n
            sig[s] = 1
            for v in order:
                if v == s:
                    continue
                sig[v] = sum(sig[u] * m for (u, w), (L, m) in arcs.items()
                             if w == v and d[s][u] is not None and d[s][u] + L == d[s][v])
            self.sigma.append(sig)

    def on_path(self, s, v, t):
        d = self.d
        return (d[s][v] is not None and d[v][t] is not None and d[s][t] is not None
                and d[s][v] + d[v][t] == d[s][t])

    def node_share(self, s, v, t):
        if v in (s, t) or not self.on_path(s, v, t):
            return Fraction(0)
        return Fraction(self.sigma[s][v] * self.sigma[v][t], self.sigma[s][t])

    def arc_share(self, s, u, v, t):
        """Share of one shortest arc u -> v in the shortest s-t paths."""
        d = self.d
        L = self.arcs[(u, v)][0]
        if d[s][u] is None or d[v][t] is None or d[s][t] is None or d[s][u] + L + d[v][t] != d[s][t]:
            return Fraction(0)
        return Fraction(self.sigma[s][u] * self.sigma[v][t], self.sigma[s][t])


def rescale_nodes(values, n, directed, normalized, endpoints, rescale, sampled):
    """The documented scaling conventions (see path_centrality.rs)."""
    if rescale == "unordered":
        half = Fraction(1) if directed else Fraction(1, 2)
        if normalized:
            pairs = n * (n - 1) if endpoints else (n - 1) * (n - 2)
            f = half / pairs if pairs > 0 else half
        else:
            f = half
        return [x * f for x in values]
    big = n if endpoints else n - 1
    if big < 2:
        return values
    corr = 1 if directed else 2
    if rescale == "networkx_subset" or sampled is None or endpoints:
        k = big if (rescale == "networkx_subset" or sampled is None) else len(sampled)
        f = Fraction(1, k * (big - 1)) if normalized else Fraction(big, k * corr)
        return [x * f for x in values]
    k = len(sampled)
    src = set(sampled)
    out = []
    for v, x in enumerate(values):
        if v in src:
            if k <= 1:
                out.append(None)
                continue
            f = Fraction(1, (k - 1) * (big - 1)) if normalized else Fraction(big, (k - 1) * corr)
        else:
            f = Fraction(1, k * (big - 1)) if normalized else Fraction(big, k * corr)
        out.append(x * f)
    return out


def rescale_edges(values, n, directed, normalized, rescale, sampled):
    if rescale == "unordered":
        half = Fraction(1) if directed else Fraction(1, 2)
        f = half / (n * (n - 1)) if normalized and n > 1 else half
        return [x * f for x in values]
    big = n
    if big < 2:
        return values
    corr = 1 if directed else 2
    k = len(sampled) if (rescale == "networkx" and sampled is not None) else big
    f = Fraction(1, k * (big - 1)) if normalized else Fraction(big, k * corr)
    return [x * f for x in values]


def betweenness(fx, e):
    n = fx.n
    weighted = bool(e.get("weighted"))
    distinct = e.get("parallel", "distinct") == "distinct"
    p = Paths(fx, weighted, distinct)
    sampled = None
    if "sample" in e:
        k = min(int(e["sample"]), n)
        sources = sampled_sources(n, k, int(e.get("seed", 1)))
        sampled = sources
    else:
        sources = node_subset(n, e.get("sources", "all"))
    targets = node_subset(n, e.get("targets", "all"))
    rescale = e.get("rescale", "unordered")
    normalized = bool(e.get("normalized", False))
    if e.get("output", "nodes") == "edges":
        m = len(fx.src)
        pair_value = {}
        for (u, v) in p.arcs:
            pair_value[(u, v)] = sum((p.arc_share(s, u, v, t) for s in sources for t in targets
                                      if s != t), Fraction(0))
        vals = [Fraction(0)] * m
        groups = {}
        for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
            key = (a, b) if fx.directed else (min(a, b), max(a, b))
            groups.setdefault(key, []).append(i)
        for key, edges in groups.items():
            a, b = key
            if a == b:
                continue
            arcs = [(a, b)] if fx.directed else [(a, b), (b, a)]
            shortest = min(p.edge_len[i] for i in edges)
            tied = [i for i in edges if p.edge_len[i] == shortest]
            for (u, v) in arcs:
                val = pair_value.get((u, v), Fraction(0))
                if distinct:
                    for i in tied:
                        vals[i] += val
                else:
                    for i in tied:
                        vals[i] += val / len(tied)
        return [mark(x) for x in rescale_edges(vals, n, fx.directed, normalized, rescale,
                                               sampled)]
    endpoints = bool(e.get("endpoints", False))
    tset = set(targets)
    vals = [Fraction(0)] * n
    for s in sources:
        for t in targets:
            if s == t:
                continue
            for v in range(n):
                vals[v] += p.node_share(s, v, t)
            if endpoints and p.d[s][t] is not None:
                vals[s] += 1
                vals[t] += 1
    del tset
    return [mark(x) for x in rescale_nodes(vals, n, fx.directed, normalized, endpoints, rescale,
                                           sampled)]


def percolation(fx, e):
    n = fx.n
    if n <= 2:
        return [UNDEFINED] * n
    p = Paths(fx, bool(e.get("weighted")), e.get("parallel", "collapse") == "distinct")
    x = [Fraction(v).limit_denominator(1000) for v in states(n)]
    total = sum(x)
    out = []
    for v in range(n):
        acc = Fraction(0)
        undefined = False
        for s in range(n):
            for t in range(n):
                if s == t or v in (s, t):
                    continue
                share = p.node_share(s, v, t)
                if total - x[v] == 0:
                    if p.d[s][v] is not None:
                        undefined = True
                    continue
                acc += share * x[s] / (total - x[v])
        out.append(UNDEFINED if undefined else float(acc / (n - 2)))
    return out


def load(fx, e):
    n = fx.n
    p = Paths(fx, bool(e.get("weighted")), False)
    cutoff = e.get("cutoff")
    lim = None if cutoff is None else Fraction(cutoff)
    vals = [Fraction(0)] * n
    for s in range(n):
        reach = [v for v in range(n) if p.d[s][v] is not None and (lim is None or p.d[s][v] <= lim)]
        preds = {v: [u for (u, w), (L, _) in p.arcs.items()
                     if w == v and u in reach and p.d[s][u] + L == p.d[s][v]] for v in reach}
        for t in reach:
            if t == s:
                continue
            flow = {t: Fraction(1)}
            for w in sorted((v for v in reach if p.on_path(s, v, t)), key=lambda v: -p.d[s][v]):
                f = flow.get(w, Fraction(0))
                if w == s or f == 0:
                    continue
                for u in preds[w]:
                    flow[u] = flow.get(u, Fraction(0)) + f / len(preds[w])
            for v, f in flow.items():
                if v not in (s, t):
                    vals[v] += f
    if e.get("normalized") and n > 2:
        vals = [x / ((n - 1) * (n - 2)) for x in vals]
    return [float(x) for x in vals]


# ---- Electrical measures --------------------------------------------------------


def conductances(fx, weighted, as_resistance):
    """Total conductance of every unordered pair, and per-edge conductances."""
    n = fx.n
    c = np.zeros((n, n))
    per_edge = []
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        w = fx.weights[i] if (weighted and fx.weighted) else 1.0
        g = 1.0 / w if as_resistance else w
        per_edge.append(0.0 if a == b else g)
        if a != b:
            c[a, b] += g
            c[b, a] += g
    return c, per_edge


def laplacian_pinv(c):
    lap = np.diag(c.sum(axis=1)) - c
    return np.linalg.pinv(lap, hermitian=True)


def components(c):
    n = c.shape[0]
    comp = [-1] * n
    k = 0
    for s in range(n):
        if comp[s] >= 0:
            continue
        stack = [s]
        comp[s] = k
        while stack:
            v = stack.pop()
            for u in np.nonzero(c[v])[0]:
                if comp[u] < 0:
                    comp[u] = k
                    stack.append(u)
        k += 1
    return comp, k


def throughputs(c, lp, s, t):
    """Node throughputs and edge currents (absolute, per unordered pair)."""
    b = np.zeros(c.shape[0])
    b[s], b[t] = 1.0, -1.0
    pot = lp @ b
    cur = np.abs(c * (pot[:, None] - pot[None, :]))
    return 0.5 * cur.sum(axis=1), cur


def current_flow(fx, e):
    n = fx.n
    weighted = bool(e.get("weighted"))
    as_res = e.get("weights_as", "strength") == "distance"
    c, per_edge = conductances(fx, weighted, as_res)
    comp, k = components(c)
    edges = e.get("output", "nodes") == "edges"
    m = len(fx.src)
    if k > 1:
        return [UNDEFINED] * (m if edges else n)
    lp = laplacian_pinv(c)
    normalized = e.get("normalized", True)
    fn = e["fn"]
    node = np.zeros(n)
    pair_cur = np.zeros((n, n))
    if fn == "cx_current_flow":
        if n < 3:
            return [0.0] * (m if edges else n)
        for s, t in itertools.combinations(range(n), 2):
            tau, cur = throughputs(c, lp, s, t)
            tau[s] = tau[t] = 0.0
            node += tau
            pair_cur += cur
        scale = 2.0 / ((n - 1) * (n - 2)) if normalized else 1.0
        node *= scale
        pair_cur *= 0.5 * scale
    elif fn == "cx_current_flow_subset":
        src = node_subset(n, e.get("sources", "even"))
        tgt = node_subset(n, e.get("targets", "thirds"))
        for s in src:
            for t in tgt:
                if s == t:
                    continue
                tau, cur = throughputs(c, lp, s, t)
                tau[s] = tau[t] = 0.0
                node += tau
                pair_cur += 0.5 * cur
        scale = 1.0 / ((n - 1) * (n - 2)) if normalized and n > 2 else 1.0
        node *= scale
        pair_cur *= scale
    elif fn == "cx_current_flow_approx":
        pairs = approximation_pairs(fx, int(e.get("seed", 1)), float(e.get("epsilon", 0.5)))
        if n < 3 or not pairs:
            return [0.0] * n
        cstar = n * (n - 1) / ((n - 1.0) * (n - 2.0))
        for s, t in pairs:
            tau, _ = throughputs(c, lp, s, t)
            tau[s] = tau[t] = 0.0
            node += tau
        node *= cstar / len(pairs)
        if not normalized:
            node *= (n - 1.0) * (n - 2.0) / 2.0
    if edges:
        out = []
        for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
            out.append(0.0 if a == b else float(pair_cur[a, b] * per_edge[i] / c[a, b]))
        return out
    return [float(x) for x in node]


def resistance_matrix(fx, weighted, as_res):
    n = fx.n
    c, _ = conductances(fx, weighted, as_res)
    comp, k = components(c)
    r = np.full((n, n), np.inf)
    for q in range(k):
        nodes = [v for v in range(n) if comp[v] == q]
        sub = c[np.ix_(nodes, nodes)]
        lp = laplacian_pinv(sub)
        for i, a in enumerate(nodes):
            for j, b in enumerate(nodes):
                r[a, b] = lp[i, i] + lp[j, j] - 2 * lp[i, j]
    np.fill_diagonal(r, 0.0)
    return r


def adjacency(fx, weighted):
    """Parallel edges summed, a self-loop once on the diagonal."""
    n = fx.n
    a = np.zeros((n, n))
    for i, (u, v) in enumerate(zip(fx.src, fx.tgt)):
        w = fx.weights[i] if (weighted and fx.weighted) else 1.0
        a[u, v] += w
        if u != v:
            a[v, u] += w
    return a


def kemeny(fx, weighted):
    a = adjacency(fx, weighted)
    _, k = components(a - np.diag(np.diag(a)))
    if k > 1:
        return UNDEFINED
    d = a.sum(axis=1)
    if fx.n == 1:
        return 0.0
    p = a / d[:, None]
    pi = d / d.sum()
    z = np.linalg.inv(np.eye(fx.n) - p + np.outer(np.ones(fx.n), pi))
    return float(np.trace(z) - 1.0)


def second_order(fx, weighted):
    n = fx.n
    a = adjacency(fx, weighted)
    _, k = components(a - np.diag(np.diag(a)))
    if k > 1:
        return [UNDEFINED] * n
    deg = a.sum(axis=1)
    dmax = deg.max()
    for i in range(n):
        if deg[i] < dmax:
            a[i, i] = dmax - deg[i]
    p = a / a.sum(axis=1)[:, None]
    out = []
    for i in range(n):
        q = p.copy()
        q[:, i] = 0.0
        mfpt = np.linalg.solve(np.eye(n) - q, np.ones(n))
        out.append(float(np.sqrt(2 * mfpt.sum() - n * (n + 1))))
    return out


def binary(fx):
    n = fx.n
    a = np.zeros((n, n))
    for u, v in zip(fx.src, fx.tgt):
        a[u, v] = a[v, u] = 1.0
    return a


def communicability_betweenness(fx):
    n = fx.n
    a = binary(fx)
    g = scipy.linalg.expm(a)
    out = []
    for r in range(n):
        b = a.copy()
        b[r, :] = 0
        b[:, r] = 0
        gr = scipy.linalg.expm(b)
        total = 0.0
        for p_ in range(n):
            for q in range(n):
                if p_ == q or r in (p_, q) or g[p_, q] == 0.0:
                    continue
                total += (g[p_, q] - gr[p_, q]) / g[p_, q]
        out.append(total)
    if n > 2:
        out = [x / ((n - 1) ** 2 - (n - 1)) for x in out]
    return out


def laplacian_centrality(fx, weighted, normalized):
    n = fx.n
    a = adjacency(fx, weighted)
    np.fill_diagonal(a, 0.0)
    lap = np.diag(a.sum(axis=1)) - a
    energy = float((lap ** 2).sum())
    out = []
    for i in range(n):
        keep = [j for j in range(n) if j != i]
        sub = lap[np.ix_(keep, keep)].copy()
        for k, j in enumerate(keep):
            sub[k, k] = lap[j, j] - abs(lap[j, i])
        drop = energy - float((sub ** 2).sum())
        if normalized:
            out.append(UNDEFINED if energy == 0 else drop / energy)
        else:
            out.append(drop)
    return out


def trophic_levels(fx, weighted, drop_loops=False):
    n = fx.n
    edges = [(a, b, Fraction(fx.weights[i]) if (weighted and fx.weighted) else Fraction(1))
             for i, (a, b) in enumerate(zip(fx.src, fx.tgt)) if not (drop_loops and a == b)]
    indeg = [0] * n
    for a, b, _ in edges:
        indeg[b] += 1
    basal = [v for v in range(n) if indeg[v] == 0]
    if not basal:
        return None, edges
    reach = set(basal)
    frontier = list(basal)
    while frontier:
        v = frontier.pop()
        for a, b, _ in edges:
            if a == v and b not in reach:
                reach.add(b)
                frontier.append(b)
    if len(reach) < n:
        return None, edges
    others = [v for v in range(n) if indeg[v] > 0]
    pos = {v: i for i, v in enumerate(others)}
    k = len(others)
    strength = [Fraction(0)] * n
    for a, b, w in edges:
        strength[b] += w
    # s_i - sum_j (w_ji / k_i) s_j = 1 over non-basal i, basal levels 1.
    m = [[Fraction(0)] * (k + 1) for _ in range(k)]
    for i in range(k):
        m[i][i] = Fraction(1)
        m[i][k] = Fraction(1)
    for a, b, w in edges:
        i = pos[b]
        if a in pos:
            m[i][pos[a]] -= w / strength[b]
        else:
            m[i][k] += w / strength[b]
    for col in range(k):
        piv = next((r for r in range(col, k) if m[r][col] != 0), None)
        if piv is None:
            return None, edges
        m[col], m[piv] = m[piv], m[col]
        for r in range(k):
            if r != col and m[r][col] != 0:
                f = m[r][col] / m[col][col]
                m[r] = [x - f * y for x, y in zip(m[r], m[col])]
    levels = [Fraction(1)] * n
    for i, v in enumerate(others):
        levels[v] = m[i][k] / m[i][i]
    return levels, edges


def trophic(fx, e):
    weighted = bool(e.get("weighted"))
    what = e["what"]
    if what == "levels":
        levels, _ = trophic_levels(fx, weighted)
        return [UNDEFINED] * fx.n if levels is None else [float(x) for x in levels]
    if what == "differences":
        levels, _ = trophic_levels(fx, weighted)
        if levels is None:
            return [UNDEFINED] * len(fx.src)
        return [float(levels[b] - levels[a]) for a, b in zip(fx.src, fx.tgt)]
    levels, edges = trophic_levels(fx, weighted, drop_loops=not e.get("cannibalism", False))
    if levels is None or not edges:
        return [UNDEFINED]
    diffs = {(a, b): levels[b] - levels[a] for a, b, _ in edges}
    vals = list(diffs.values())
    mean = sum(vals) / len(vals)
    var = sum((x - mean) ** 2 for x in vals) / len(vals)
    return [math.sqrt(var)]


def reaching(fx, e):
    n = fx.n
    if n <= 1 or not fx.src:
        return [UNDEFINED] * (1 if e.get("what") == "global" else n)
    p = Paths(fx, False, False)
    local = []
    for v in range(n):
        if fx.directed:
            local.append(Fraction(sum(1 for u in range(n) if u != v and p.d[v][u] is not None),
                                  n - 1))
        else:
            local.append(sum((1 / p.d[v][u] for u in range(n)
                              if u != v and p.d[v][u] is not None), Fraction(0)) / (n - 1))
    if e.get("what") == "global":
        top = max(local)
        return [float(sum(top - c for c in local) / (n - 1))]
    return [float(x) for x in local]


def dispersion(fx, e):
    nbrs = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            nbrs[a].add(b)
            nbrs[b].add(a)
    pairs = sorted({(a, b) for a, b in zip(fx.src, fx.tgt)} | {(b, a) for a, b in zip(fx.src, fx.tgt)})
    out = []
    for u, v in pairs:
        common = sorted(nbrs[u] & nbrs[v])
        total = 0
        for s, t in itertools.combinations(common, 2):
            if t in nbrs[s]:
                continue
            if (nbrs[s] & nbrs[t] & nbrs[u]) - {u, v}:
                continue
            total += 1
        emb = len(common)
        if e.get("normalized", True):
            alpha, b, c = float(e.get("alpha", 1.0)), float(e.get("b", 0.0)), float(e.get("c", 0.0))
            val = (total + b) ** alpha
            if emb + c != 0:
                val /= emb + c
            out.append(val)
        else:
            out.append(float(total))
    return out


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    n = fx.n
    if fn == "cx_betweenness_subset":
        return betweenness(fx, e)
    if fn == "cx_percolation":
        return percolation(fx, e)
    if fn == "cx_load":
        return load(fx, e)
    if fn in ("cx_current_flow", "cx_current_flow_subset", "cx_current_flow_approx"):
        return current_flow(fx, e)
    if fn == "cx_current_flow_closeness":
        if n <= 1:
            return [1.0] * n
        r = resistance_matrix(fx, bool(e.get("weighted")), e.get("weights_as", "strength") == "distance")
        return [mark(1.0 / r[v].sum()) if np.isfinite(r[v].sum()) else UNDEFINED for v in range(n)]
    if fn == "cx_resistance":
        r = resistance_matrix(fx, bool(e.get("weighted")), e.get("weights_as", "distance") == "distance")
        return [tuple(mark(r[i, j]) for j in range(n)) for i in range(n)]
    if fn == "cx_effective_resistance":
        if n == 0:
            return [UNDEFINED]
        r = resistance_matrix(fx, bool(e.get("weighted")), e.get("weights_as", "distance") == "distance")
        return [mark(sum(r[i, j] for i in range(n) for j in range(i + 1, n)))]
    if fn == "cx_kemeny":
        return [mark(kemeny(fx, bool(e.get("weighted"))))]
    if fn == "cx_second_order":
        return second_order(fx, bool(e.get("weighted")))
    if fn == "cx_communicability":
        g = scipy.linalg.expm(binary(fx))
        return [tuple(float(g[i, j]) for j in range(n)) for i in range(n)]
    if fn == "cx_communicability_betweenness":
        return communicability_betweenness(fx)
    if fn == "cx_laplacian":
        return [mark(x) for x in laplacian_centrality(fx, bool(e.get("weighted")),
                                                      e.get("normalized", True))]
    if fn == "cx_trophic":
        return trophic(fx, e)
    if fn == "cx_reaching":
        return reaching(fx, e)
    if fn == "cx_dispersion":
        return dispersion(fx, e)
    raise ValueError(f"no definition oracle for {fn}")


def main() -> int:
    manifest = load_area("centrality_extra")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        if fx.n > MAX_NODES:
            continue
        tags = set(fx.tags) | {"directed" if fx.directed else "undirected"}
        for variant, spec in manifest.items():
            ospec = spec["oracles"].get(ORACLE)
            if ospec is None or not applies(spec, fx):
                continue
            if name in ospec.get("undefined_on", []) or tags & set(ospec.get("undefined_when",
                                                                              [])):
                continue
            write_values(ORACLE, name, variant, compute(spec, fx))
            count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (centrality_extra): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
