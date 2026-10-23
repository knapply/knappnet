#!/usr/bin/env python3
"""Golden values for cohesive subsets, community partitions, and
topography from NetworkX (Hagberg, Schult, and Swart 2008).

Every measure here reads the network as an undirected simple graph
(`nx.Graph` built from the edges, self-loops removed), except compactness,
which uses directed distances on directed fixtures.

- Bicomponents: `nx.biconnected_components`.
- n-cliques: `nx.find_cliques` on the graph joining nodes at distance at
  most n (`nx.single_source_shortest_path_length` with a cutoff); n-clans:
  the n-cliques whose induced subgraph is connected with diameter at most
  n (`nx.diameter`).
- Edge connectivity: `nx.edge_connectivity(G, u, v)` for every pair (a
  separate maximum flow per pair); lambda sets: the classes of
  connectivity >= k with at least two members, for every k >= 1.
- Girvan-Newman: the highest `nx.community.modularity` over the
  partitions yielded by `nx.community.girvan_newman`.
- Greedy modularity: the modularity of
  `nx.community.greedy_modularity_communities`. NetworkX breaks ties
  between equally good mergers by the order of its floating-point heap,
  the engine by node order, so the two can end in different partitions.
  NetworkX's merge sequence is therefore replayed in exact arithmetic:
  when every merger it made was a best merger at the time (a pure
  tie-break difference) and its modularity differs from the one the
  engine's rule gives, no golden file is written and the fixture is
  reported, so that the manifest can list it under `undefined_on`. A
  merger that was not a best merger would be a NetworkX error, and the
  golden file is written so that the comparison fails loudly.
- Compactness: `nx.global_efficiency` on undirected unweighted fixtures;
  otherwise the mean of 1 / d over ordered pairs from
  `nx.shortest_path_length` (directed distances, weights as distances).
- Small-world inputs: `nx.transitivity` or `nx.average_clustering`, and
  the mean of `nx.shortest_path_length` over connected ordered pairs.

    python3 validation/oracles/groups_networkx_oracle.py
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

import networkx as nx
from networkx.algorithms.community.modularity_max import (
    _greedy_modularity_communities_generator as greedy_generator,
)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import GOLDEN, applies, fixture_names, load_fixture, write_values  # noqa: E402
from common import load_manifest  # noqa: E402

ORACLE = "networkx"
MEASURES = {"bicomponents", "n_cliques", "n_clans", "lambda_connectivity", "lambda_sets",
            "girvan_newman", "fast_greedy", "compactness", "small_world_observed"}


def simple(fx) -> nx.Graph:
    g = nx.Graph()
    g.add_nodes_from(range(fx.n))
    g.add_edges_from((a, b) for a, b in zip(fx.src, fx.tgt) if a != b)
    return g


def multigraph(fx, weighted: bool):
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, weight=fx.weights[i] if weighted else 1.0)
    return g


def power_cliques(g: nx.Graph, n: int) -> list[tuple[int, ...]]:
    p = nx.Graph()
    p.add_nodes_from(g)
    for v in g:
        for u, d in nx.single_source_shortest_path_length(g, v, cutoff=n).items():
            if u != v:
                p.add_edge(u, v)
    return [tuple(sorted(c)) for c in nx.find_cliques(p)]


def is_clan(g: nx.Graph, members, n: int) -> bool:
    if len(members) == 1:
        return True
    sub = g.subgraph(members)
    return nx.is_connected(sub) and nx.diameter(sub) <= n


def connectivity(g: nx.Graph):
    n = g.number_of_nodes()
    m = [[0] * n for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            m[a][b] = m[b][a] = nx.edge_connectivity(g, a, b) if g.degree(a) and g.degree(b) \
                else 0
    return m


def lambda_sets(m, n):
    top = max((max(r) for r in m), default=0)
    out = set()
    for k in range(1, top + 1):
        seen = [False] * n
        for v in range(n):
            if seen[v]:
                continue
            cls = [u for u in range(n) if u == v or m[v][u] >= k]
            for u in cls:
                seen[u] = True
            if len(cls) > 1:
                out.add(tuple(cls))
    return sorted(out)


class TieBreak(Exception):
    """NetworkX's greedy path differs from the engine's only by tie-breaking."""


def greedy_with_tie_check(fx, weighted: bool, gamma: float):
    """NetworkX's greedy modularity partition, checked in exact arithmetic:
    every merger NetworkX makes must be a best merger at that step. Returns
    (partition, consistent), where partition is the one
    `greedy_modularity_communities` returns (the first at which no merger
    would raise modularity)."""
    mg = multigraph(fx, weighted)
    n = fx.n
    gamma = Fraction(gamma)
    a = [[Fraction(0)] * n for _ in range(n)]
    k = [Fraction(0)] * n
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        w = Fraction(fx.weights[i]) if weighted else Fraction(1)
        k[s] += w
        k[t] += w
        if s != t:
            a[s][t] += w
            a[t][s] += w
    two_m = sum(k)
    gen = greedy_generator(mg, weight="weight" if weighted else None,
                           resolution=float(gamma))
    groups = [frozenset(c) for c in next(gen)]
    consistent = True
    while len(groups) > 1:
        try:
            dq = next(gen)
        except StopIteration:
            # No two groups are connected any more.
            break
        if dq < 0:
            break
        new = [frozenset(c) for c in next(gen)]
        merged = next(c for c in new if c not in groups)
        parts = [c for c in groups if c <= merged]
        assert len(parts) == 2 and parts[0] | parts[1] == merged
        gains = {}
        for i, x in enumerate(groups):
            for y in groups[i + 1:]:
                w = sum(a[u][v] for u in x for v in y)
                if w:
                    gains[(x, y)] = two_m * w - gamma * sum(k[u] for u in x) * sum(k[v] for v in y)
        best = max(gains.values())
        made = gains.get((parts[0], parts[1]), gains.get((parts[1], parts[0])))
        if made != best:
            consistent = False
        groups = new
    return groups, consistent


def membership_of(groups, n):
    label = [0] * n
    for c, grp in enumerate(groups):
        for v in grp:
            label[v] = c
    return label


def mean_reciprocal(fx, weighted: bool) -> float:
    if fx.n < 2:
        return -1.0
    if fx.directed:
        g = nx.DiGraph()
    else:
        g = nx.Graph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        w = fx.weights[i] if weighted else 1.0
        # Parallel edges: the shortest counts.
        if g.has_edge(a, b):
            g[a][b]["weight"] = min(g[a][b]["weight"], w)
        else:
            g.add_edge(a, b, weight=w)
    if not fx.directed and not weighted:
        return nx.global_efficiency(g)
    total = 0.0
    for v, dist in nx.shortest_path_length(g, weight="weight"):
        total += sum(1.0 / d for u, d in dist.items() if u != v)
    return total / (fx.n * (fx.n - 1))


def engine_rule_modularity(fx, weighted: bool, gamma: float) -> float:
    """Modularity of the greedy partition under the engine's tie rule
    (lowest-numbered groups merge first), in exact arithmetic, as the
    definition oracle computes it."""
    n = fx.n
    gamma = Fraction(gamma)
    a = [[Fraction(0)] * n for _ in range(n)]
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        w = Fraction(fx.weights[i]) if weighted else Fraction(1)
        if s == t:
            a[s][s] += 2 * w
        else:
            a[s][t] += w
            a[t][s] += w
    k = [sum(r) for r in a]
    two_m = sum(k)
    groups = {v: {v} for v in range(n)}
    total, best_total, best_groups = Fraction(0), Fraction(0), {v: {v} for v in range(n)}
    while True:
        best = None
        names = sorted(groups)
        for i, x in enumerate(names):
            for y in names[i + 1:]:
                w = sum(a[u][v] for u in groups[x] for v in groups[y])
                if w == 0:
                    continue
                score = two_m * w - gamma * sum(k[u] for u in groups[x]) * \
                    sum(k[v] for v in groups[y])
                if best is None or score > best[0]:
                    best = (score, x, y)
        if best is None:
            break
        total += best[0]
        groups[best[1]] |= groups.pop(best[2])
        if total > best_total:
            best_total = total
            best_groups = {x: set(m) for x, m in groups.items()}
    q = Fraction(0)
    for members in best_groups.values():
        for u in members:
            for v in members:
                q += a[u][v] - gamma * k[u] * k[v] / two_m
    return float(q / two_m)


def compute(spec: dict, fx):
    e = spec["engine"]
    fn = e["fn"]
    g = simple(fx)
    if fn == "bicomponents":
        return [tuple(sorted(c)) for c in nx.biconnected_components(g)]
    if fn in ("n_cliques", "n_clans"):
        cl = [c for c in power_cliques(g, e["n"]) if len(c) >= e.get("min_size", 0)]
        if fn == "n_clans":
            cl = [c for c in cl if is_clan(g, c, e["n"])]
        return cl
    if fn == "lambda_connectivity":
        m = connectivity(g)
        return [m[i][j] for i in range(fx.n) for j in range(i + 1, fx.n)]
    if fn == "lambda_sets":
        return lambda_sets(connectivity(g), fx.n)
    if fn == "girvan_newman":
        if g.number_of_edges() == 0:
            return [-1.0]
        best = max([nx.community.modularity(g, c) for c in nx.community.girvan_newman(g)]
                   + [nx.community.modularity(g, list(nx.connected_components(g)))])
        return [best]
    if fn == "fast_greedy":
        weight = "weight" if e.get("weighted", False) else None
        mg = multigraph(fx, e.get("weighted", False))
        if mg.number_of_edges() == 0:
            return list(range(fx.n)) if e.get("output") == "membership" else [-1.0]
        gamma = e.get("resolution", 1.0)
        groups, consistent = greedy_with_tie_check(fx, e.get("weighted", False), gamma)
        if e.get("output") == "membership":
            return membership_of(groups, fx.n)
        q = nx.community.modularity(mg, groups, weight=weight, resolution=gamma)
        if consistent and abs(q - engine_rule_modularity(fx, e.get("weighted", False),
                                                         gamma)) > 1e-9:
            raise TieBreak
        return [q]
    if fn == "compactness":
        assert e.get("mode", "out") == "out"
        return [mean_reciprocal(fx, e.get("weighted", False))]
    if fn == "small_world_observed":
        if e.get("local", False):
            c = nx.average_clustering(g) if fx.n else -1.0
        else:
            c = nx.transitivity(g)
            # NetworkX reports 0 when there are no connected triples; the
            # engine reports the value as undefined.
            if sum(d * (d - 1) for _, d in g.degree()) == 0:
                c = -1.0
        total, pairs = 0, 0
        for v, dist in nx.shortest_path_length(g):
            for u, d in dist.items():
                if u != v:
                    total += d
                    pairs += 1
        return [c, total / pairs if pairs else -1.0]
    return None


def main() -> None:
    manifest = load_manifest()
    count = 0
    tie_breaks = []
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
            except TieBreak:
                tie_breaks.append((variant, name))
                continue
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    (GOLDEN / ORACLE).mkdir(parents=True, exist_ok=True)
    (GOLDEN / ORACLE / "VERSION").write_text(f"networkx {nx.__version__}\n", newline="")
    print(f"networkx {nx.__version__} (groups): wrote {count} golden files")
    for variant, name in tie_breaks:
        print(f"  tie-break difference (list under undefined_on): {variant} on {name}")


if __name__ == "__main__":
    main()
