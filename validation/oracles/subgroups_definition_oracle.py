#!/usr/bin/env python3
"""Golden values for cohesive subgroups, computed from the definitions.

Shares no code or algorithmic strategy with the engine:

- Coreness (Seidman 1983): for k = 0, 1, 2, ..., repeatedly delete every
  node whose degree in the remaining network is below k; the nodes left
  form the k-core. A node's coreness is the largest k whose core contains
  it. Degrees count parallel edges and loops as igraph does (an undirected
  loop adds 2; a directed loop adds 1 out and 1 in).
- Maximal cliques: every clique of the undirected simple graph underlying
  the network is grown in increasing node order, and a clique is kept when
  no node at all is adjacent to every member.
- Modularity: the double sum over all ordered node pairs of
  [A_ij - gamma k_i k_j / 2m] delta(c_i, c_j) / 2m (undirected; Newman and
  Girvan 2004) or [A_ij - gamma k_i^out k_j^in / m] delta(c_i, c_j) / m
  (directed; Leicht and Newman 2008), in exact rational arithmetic.

For Louvain, partitions legitimately differ between implementations, so
this oracle evaluates the modularity of the partition the engine returns;
a match shows the engine reports the modularity of its own partition.

    python3 validation/oracles/subgroups_definition_oracle.py
"""

from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_fixture,
                    load_manifest, write_values, write_version)  # noqa: E402

ORACLE = "definition"
MAX_NODES = 80
MEASURES = {"coreness", "max_cliques", "clique_number", "modularity", "louvain_modularity"}


def coreness(fx, mode: str) -> list[float]:
    def degree(alive, v):
        d = 0
        for a, b in zip(fx.src, fx.tgt):
            if a not in alive or b not in alive:
                continue
            if not fx.directed or mode == "all":
                d += (a == v) + (b == v)
            elif mode == "out":
                d += a == v
            else:
                d += b == v
        return d

    core = [0] * fx.n
    alive = set(range(fx.n))
    k = 0
    while alive:
        k += 1
        changed = True
        while changed:
            low = {v for v in alive if degree(alive, v) < k}
            changed = bool(low)
            alive -= low
        for v in alive:
            core[v] = k
    return [float(c) for c in core]


def simple_adjacency(fx) -> list[set[int]]:
    adj = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            adj[a].add(b)
            adj[b].add(a)
    return adj


def maximal_cliques(fx) -> list[tuple[int, ...]]:
    adj = simple_adjacency(fx)
    out = []

    def grow(clique, common):
        # `common` holds every node adjacent to all members, of any index.
        if not common:
            out.append(tuple(clique))
        last = clique[-1]
        for v in sorted(common):
            if v > last:
                grow(clique + [v], common & adj[v])

    for v in range(fx.n):
        grow([v], set(adj[v]))
    return out


def modularity(fx, membership, gamma: Fraction, weighted: bool) -> float:
    n = fx.n
    a = [[Fraction(0)] * n for _ in range(n)]
    for i, (s, t) in enumerate(zip(fx.src, fx.tgt)):
        w = Fraction(fx.weights[i]) if weighted else Fraction(1)
        if fx.directed:
            a[s][t] += w
        elif s == t:
            a[s][s] += 2 * w
        else:
            a[s][t] += w
            a[t][s] += w
    k_out = [sum(row) for row in a]
    k_in = [sum(a[i][j] for i in range(n)) for j in range(n)]
    total = sum(k_out)  # m for directed, 2m for undirected
    if total == 0:
        return float("nan")
    q = Fraction(0)
    for i in range(n):
        for j in range(n):
            if membership[i] == membership[j]:
                q += a[i][j] - gamma * k_out[i] * k_in[j] / total
    return float(q / total)


def engine_louvain(fx, spec):
    # Imported lazily: only Louvain needs the engine (for its partition).
    from engine import Graph, _load_bindings
    _load_bindings()
    binding = sys.modules["bindings_subgroups"]
    g = Graph(fx.n, fx.src, fx.tgt, fx.weights, fx.directed)
    membership, _ = binding.louvain(g, spec.get("resolution", 1.0),
                                    spec.get("weighted", False), spec.get("seed", 1))
    return membership


def compute(spec: dict, fx):
    e = spec["engine"]
    fn = e["fn"]
    if fn == "coreness":
        return coreness(fx, e["mode"])
    if fn == "max_cliques":
        return [c for c in maximal_cliques(fx) if len(c) >= e.get("min_size", 0)]
    if fn == "clique_number":
        return [max((len(c) for c in maximal_cliques(fx)), default=0)]
    gamma = Fraction(e.get("resolution", 1.0))
    weighted = e.get("weighted", False)
    if fn == "modularity":
        assert e["partition"] == "blocks3"
        blocks = [3 * v // fx.n for v in range(fx.n)]
        return [modularity(fx, blocks, gamma, weighted)]
    if fn == "louvain_modularity":
        return [modularity(fx, engine_louvain(fx, e), gamma, weighted)]
    return None


def main() -> None:
    manifest = load_manifest()
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        if fx.n > MAX_NODES:
            continue
        for variant, spec in manifest.items():
            if spec["engine"]["fn"] not in MEASURES:
                continue
            if ORACLE not in spec["oracles"] or not applies(spec, fx):
                continue
            values = compute(spec, fx)
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (subgroups): wrote {count} golden files")


if __name__ == "__main__":
    main()
