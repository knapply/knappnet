#!/usr/bin/env python3
"""Times NetworkX functions in NetworkX itself and with backend="knappnet".

    python3 bench/networkx_backend.py              # full run
    python3 bench/networkx_backend.py --quick      # small graphs, short budget
    python3 bench/networkx_backend.py --budget 60  # seconds per NetworkX call

Needs the Python package installed (python -m pip install ./python). For
each function and graph it reports:

- nx: NetworkX's own implementation (backend="networkx");
- first: backend="knappnet" on a graph that has not been converted yet
  (conversion plus the engine call);
- cached: the same call again, reusing the conversion NetworkX caches in
  G.__networkx_cache__;
- convert: the conversion alone;
- speedup: nx / cached and nx / first.

Graphs are G(n, m) random graphs (largest connected component, so path
measures are defined), from small to large and sparse to dense, with
uniform random weights in [1, 10); directed functions use a directed graph
with the same n and m. A NetworkX call that exceeds the budget is stopped
(Unix only), and larger graphs are then skipped for that function. Results
of the two implementations are compared, and a mismatch is flagged.
"""

from __future__ import annotations

import argparse
import math
import os
import platform
import random
import signal
import statistics
import sys
import time
import warnings
from collections import Counter
from collections.abc import Iterator

import networkx as nx

warnings.simplefilter("ignore")
nx.config.warnings_to_ignore.add("cache")


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout


def run_limited(f, budget):
    """Seconds f() takes, or None if it ran past budget (where supported)."""
    can_stop = hasattr(signal, "setitimer")
    if can_stop:
        signal.signal(signal.SIGALRM, _alarm)
        signal.setitimer(signal.ITIMER_REAL, budget)
    try:
        t0 = time.perf_counter()
        result = f()
        return time.perf_counter() - t0, result
    except Timeout:
        return None, None
    finally:
        if can_stop:
            signal.setitimer(signal.ITIMER_REAL, 0)


def median_time(f, min_total=0.2, max_reps=7):
    """Median seconds per call, repeating short calls."""
    times = []
    total = 0.0
    while len(times) < max_reps and (total < min_total or len(times) < 3):
        t0 = time.perf_counter()
        f()
        dt = time.perf_counter() - t0
        times.append(dt)
        total += dt
        if dt > 5:
            break
    return statistics.median(times)


def materialize(x):
    if isinstance(x, Iterator):
        return list(x)
    if isinstance(x, tuple):
        return tuple(materialize(v) for v in x)
    return x


def agree(a, b):
    """Loose comparison of results (values to 1e-6, collections unordered)."""
    a, b = materialize(a), materialize(b)
    if isinstance(b, dict):
        return set(a) == set(b) and all(agree(a[k], b[k]) for k in b)
    if isinstance(b, tuple):
        return all(agree(x, y) for x, y in zip(a, b))
    if isinstance(b, list):
        if b and isinstance(b[0], tuple) and len(b[0]) == 2 and isinstance(b[0][1], dict):
            return agree(dict(a), dict(b))
        def norm(c):
            return Counter(frozenset(x) if isinstance(x, (set, frozenset, list)) else x
                           for x in c)

        return norm(a) == norm(b)
    if isinstance(b, (set, frozenset)):
        return set(a) == set(b)
    if isinstance(b, float) or isinstance(a, float):
        if math.isnan(b):
            return math.isnan(a)
        return math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-9)
    return a == b


def make_graph(n, avg_degree, directed, seed):
    m = int(n * avg_degree / (1 if directed else 2))
    G = nx.gnm_random_graph(n, m, seed=seed, directed=directed)
    comp = (max(nx.weakly_connected_components(G), key=len) if directed
            else max(nx.connected_components(G), key=len))
    G = G.subgraph(comp).copy()
    rng = random.Random(seed)
    for _, _, d in G.edges(data=True):
        d["weight"] = rng.uniform(1, 10)
    return G


# name -> (callable taking (G, backend), needs directed graph, largest n)
FUNCTIONS = {
    "betweenness_centrality": (lambda G, b: nx.betweenness_centrality(G, backend=b), False, None),
    "betweenness_centrality (weighted)": (
        lambda G, b: nx.betweenness_centrality(G, weight="weight", backend=b), False, None),
    "edge_betweenness_centrality": (
        lambda G, b: nx.edge_betweenness_centrality(G, backend=b), False, None),
    "closeness_centrality": (lambda G, b: nx.closeness_centrality(G, backend=b), False, None),
    "harmonic_centrality": (lambda G, b: nx.harmonic_centrality(G, backend=b), False, None),
    "pagerank": (lambda G, b: nx.pagerank(G, tol=1e-10, backend=b), False, None),
    "eigenvector_centrality_numpy": (
        lambda G, b: nx.eigenvector_centrality_numpy(G, backend=b), False, None),
    "hits": (lambda G, b: nx.hits(G, max_iter=10000, tol=1e-10, backend=b), True, None),
    "clustering": (lambda G, b: nx.clustering(G, backend=b), False, None),
    "clustering (weighted)": (lambda G, b: nx.clustering(G, weight="weight", backend=b),
                              False, None),
    "transitivity": (lambda G, b: nx.transitivity(G, backend=b), False, None),
    "core_number": (lambda G, b: nx.core_number(G, backend=b), False, None),
    "find_cliques": (lambda G, b: list(nx.find_cliques(G, backend=b)), False, None),
    "connected_components": (lambda G, b: list(nx.connected_components(G, backend=b)),
                             False, None),
    "strongly_connected_components": (
        lambda G, b: list(nx.strongly_connected_components(G, backend=b)), True, None),
    "triadic_census": (lambda G, b: nx.triadic_census(G, backend=b), True, None),
    "biconnected_components": (lambda G, b: list(nx.biconnected_components(G, backend=b)),
                               False, None),
    "all_pairs_shortest_path_length": (
        lambda G, b: dict(nx.all_pairs_shortest_path_length(G, backend=b)), False, 3000),
    "average_shortest_path_length": (
        lambda G, b: nx.average_shortest_path_length(G, backend=b), False, None),
    "diameter": (lambda G, b: nx.diameter(G, backend=b), False, None),
    "global_efficiency": (lambda G, b: nx.global_efficiency(G, backend=b), False, None),
    "constraint": (lambda G, b: nx.constraint(G, backend=b), False, None),
    "degree_assortativity_coefficient": (
        lambda G, b: nx.degree_assortativity_coefficient(G, backend=b), False, None),
    "louvain_communities": (
        lambda G, b: nx.community.modularity(
            G, nx.community.louvain_communities(G, seed=1, backend=b)), False, None),
}

# Louvain's partitions differ by design (different random streams); compare
# modularity loosely instead.
LOOSE = {"louvain_communities": 0.02}


def tiers(quick):
    if quick:
        return [("sparse", 500, 8), ("sparse", 2000, 8), ("medium", 1000, 32)]
    return [("sparse", 1000, 8), ("sparse", 10000, 8), ("sparse", 100000, 8),
            ("medium", 1000, 32), ("medium", 10000, 32),
            ("dense", 1000, 200), ("dense", 3000, 600)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--budget", type=float, default=None)
    parser.add_argument("--only", nargs="*", default=None)
    args = parser.parse_args()
    budget = args.budget or (5.0 if args.quick else 60.0)
    import knappnet  # noqa: F401  (fails early if the package is missing)
    from knappnet._graph import convert_from_nx

    print(f"# NetworkX {nx.__version__} vs backend='knappnet' ({platform.platform()}, "
          f"Python {platform.python_version()}, {os.cpu_count()} CPUs, "
          f"threads={nx.config.backends.knappnet.threads or 'all'}); "
          f"NetworkX budget {budget:g} s")
    header = ("| function | graph | n | m | nx s | first s | cached s | convert s | "
              "speedup cached | speedup first | agree |")
    print(header)
    print("|" + "---|" * 11)
    graphs = {}
    skipped = set()
    failures = 0
    for name, (fn, directed, max_n) in FUNCTIONS.items():
        if args.only and not any(o in name for o in args.only):
            continue
        for kind, n, deg in tiers(args.quick):
            if max_n and n > max_n:
                continue
            key = (kind, n, directed)
            if graphs.get(key) is None:
                graphs[key] = make_graph(n, deg, directed, seed=n + deg)
            G = graphs[key]
            G.__networkx_cache__.clear()
            weight = "weight" if "weighted" in name or name == "pagerank" else None
            t_conv = median_time(lambda: convert_from_nx(
                G, edge_attrs={weight: 1} if weight else None))
            G.__networkx_cache__.clear()
            t0 = time.perf_counter()
            res_kn = materialize(fn(G, "knappnet"))
            t_first = time.perf_counter() - t0
            t_cached = median_time(lambda: materialize(fn(G, "knappnet")))
            if name in skipped:
                t_nx, res_nx = None, None
            else:
                t_nx, res_nx = run_limited(lambda: materialize(fn(G, "networkx")), budget)
                if t_nx is None:
                    skipped.add(name)
            if t_nx is None:
                ok = "-"
                nx_s = f">{budget:g}"
                sp_c = f">{budget / t_cached:,.0f}x"
                sp_f = f">{budget / t_first:,.0f}x"
            else:
                if name in LOOSE:
                    good = abs(res_kn - res_nx) <= LOOSE[name]
                else:
                    good = agree(res_kn, res_nx)
                ok = "yes" if good else "NO"
                failures += not good
                nx_s = f"{t_nx:.3g}"
                sp_c = f"{t_nx / t_cached:,.1f}x"
                sp_f = f"{t_nx / t_first:,.1f}x"
            print(f"| {name} | {kind}{' directed' if directed else ''} | {G.number_of_nodes()} "
                  f"| {G.number_of_edges()} | {nx_s} | {t_first:.3g} | {t_cached:.3g} | "
                  f"{t_conv:.3g} | {sp_c} | {sp_f} | {ok} |", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
