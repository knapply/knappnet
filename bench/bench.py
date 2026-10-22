#!/usr/bin/env python3
"""Times the engine against NetworkX on random graphs.

    python3 bench/bench.py

igraph comparisons run from bench/bench_igraph.R. Results depend on the
machine, so report the core count alongside any numbers.
"""

from __future__ import annotations

import os
import random
import sys
import time
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
from engine import Graph  # noqa: E402


def random_graph(n: int, avg_degree: float, directed: bool, seed: int):
    rng = random.Random(seed)
    m = int(n * avg_degree / (1 if directed else 2))
    src = [rng.randrange(n) for _ in range(m)]
    tgt = [rng.randrange(n) for _ in range(m)]
    keep = [(a, b) for a, b in zip(src, tgt) if a != b]
    keep = list(dict.fromkeys(keep))
    return [a for a, _ in keep], [b for _, b in keep]


def timed(fn, repeat: int = 3) -> float:
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def main() -> None:
    print(f"cores: {os.cpu_count()}")
    print(f"{'graph':<28}{'measure':<14}{'engine (s)':>12}{'networkx (s)':>14}{'speedup':>10}")
    for n, deg, directed in [(2_000, 8, False), (2_000, 8, True), (5_000, 10, False)]:
        src, tgt = random_graph(n, deg, directed, seed=n)
        g = Graph(n, src, tgt, None, directed)
        G = nx.DiGraph() if directed else nx.Graph()
        G.add_nodes_from(range(n))
        G.add_edges_from(zip(src, tgt))
        label = f"n={n}, m={len(src)}, {'dir' if directed else 'und'}"
        cases = [("betweenness", lambda: g.betweenness(), lambda: nx.betweenness_centrality(G))]
        if directed:
            cases.append(("triad census", lambda: g.triad_census(), lambda: nx.triadic_census(G)))
        for name, ours, theirs in cases:
            a = timed(ours)
            b = timed(theirs, repeat=1)
            print(f"{label:<28}{name:<14}{a:>12.3f}{b:>14.3f}{b / a:>9.0f}x")


if __name__ == "__main__":
    main()
