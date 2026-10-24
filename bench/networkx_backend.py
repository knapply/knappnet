#!/usr/bin/env python3
"""Times NetworkX functions in NetworkX itself and with each NetworkX backend.

    python3 bench/networkx_backend.py              # full run
    python3 bench/networkx_backend.py --quick      # small graphs, short budget
    python3 bench/networkx_backend.py --budget 60  # seconds per NetworkX call
    python3 bench/networkx_backend.py --backends knappnet cugraph

Needs the Python package installed (python -m pip install ./python); the
other backends are optional (tools/bench_backends.sh installs them): knappnet,
parallel (nx-parallel), graphblas (graphblas-algorithms), and cugraph
(nx-cugraph, needs a GPU). One table row per function, graph, and backend:

- networkx: NetworkX's own implementation (backend="networkx");
- first: the backend on a graph that has not been converted yet (conversion
  plus the call);
- cached: the same call again, reusing the conversion NetworkX caches in
  G.__networkx_cache__;
- convert: the backend's conversion alone, where it exposes one;
- x cached, x first: NetworkX's time divided by the backend's;
- agree: the backend's result matches NetworkX's.

A backend that is not installed, does not implement the function, or does
not take the arguments is reported as n/a.

Graphs are G(n, m) random graphs (largest connected component, so path
measures are defined), from small to large and sparse to dense, with
uniform random weights in [1, 10); directed functions use a directed graph
with the same n and m. A call that exceeds the budget is stopped (Unix
only), and larger graphs are then skipped for NetworkX. Results are
compared with NetworkX's, and a mismatch is flagged.
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
    if hasattr(b, "shape"):
        return a.shape == b.shape and bool(((a - b) ** 2).max() < 1e-9 if a.size else True)
    if isinstance(b, nx.Graph):
        return (type(a) is type(b) and set(a) == set(b)
                and {frozenset(e) for e in a.edges} == {frozenset(e) for e in b.edges})
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


BETWEENNESS = {"max_n": 20000}  # O(nm): the 100,000-node graph takes minutes

# name -> (callable taking (G, backend), needs a directed graph, limits)
FUNCTIONS = {
    "betweenness_centrality": (
        lambda G, b: nx.betweenness_centrality(G, backend=b), False, BETWEENNESS),
    "betweenness_centrality (weighted)": (
        lambda G, b: nx.betweenness_centrality(G, weight="weight", backend=b), False,
        BETWEENNESS),
    "edge_betweenness_centrality": (
        lambda G, b: nx.edge_betweenness_centrality(G, backend=b), False, BETWEENNESS),
    "closeness_centrality": (lambda G, b: nx.closeness_centrality(G, backend=b), False, {}),
    "harmonic_centrality": (lambda G, b: nx.harmonic_centrality(G, backend=b), False, {}),
    # NetworkX stops when the L1 change is below n * tol; tol=1e-13 gives
    # agreement to about 1e-7 relative with the engine (which solves to 1e-12).
    "pagerank": (lambda G, b: nx.pagerank(G, tol=1e-13, max_iter=10000, backend=b), False, {}),
    "eigenvector_centrality_numpy": (
        lambda G, b: nx.eigenvector_centrality_numpy(G, backend=b), False, {}),
    "hits": (lambda G, b: nx.hits(G, max_iter=10000, tol=1e-10, backend=b), True, {}),
    "clustering": (lambda G, b: nx.clustering(G, backend=b), False, {}),
    "clustering (weighted)": (lambda G, b: nx.clustering(G, weight="weight", backend=b),
                              False, {}),
    "transitivity": (lambda G, b: nx.transitivity(G, backend=b), False, {}),
    "core_number": (lambda G, b: nx.core_number(G, backend=b), False, {}),
    # Dense graphs have too many maximal cliques to hold in memory.
    "find_cliques": (lambda G, b: list(nx.find_cliques(G, backend=b)), False,
                     {"max_degree": 32}),
    "connected_components": (lambda G, b: list(nx.connected_components(G, backend=b)),
                             False, {}),
    "strongly_connected_components": (
        lambda G, b: list(nx.strongly_connected_components(G, backend=b)), True, {}),
    "triadic_census": (lambda G, b: nx.triadic_census(G, backend=b), True, {}),
    "biconnected_components": (lambda G, b: list(nx.biconnected_components(G, backend=b)),
                               False, {}),
    # n^2 entries in NetworkX's dict of dicts.
    "all_pairs_shortest_path_length": (
        lambda G, b: dict(nx.all_pairs_shortest_path_length(G, backend=b)), False,
        {"max_n": 3000}),
    "average_shortest_path_length": (
        lambda G, b: nx.average_shortest_path_length(G, backend=b), False, {}),
    "diameter": (lambda G, b: nx.diameter(G, backend=b), False, {}),
    "global_efficiency": (lambda G, b: nx.global_efficiency(G, backend=b), False, {}),
    "constraint": (lambda G, b: nx.constraint(G, backend=b), False, {}),
    "degree_assortativity_coefficient": (
        lambda G, b: nx.degree_assortativity_coefficient(G, backend=b), False, {}),
    "louvain_communities": (
        lambda G, b: nx.community.modularity(
            G, nx.community.louvain_communities(G, seed=1, backend=b)), False, {}),
    # Shortest paths.
    "single_source_bellman_ford_path_length (weighted)": (
        lambda G, b: nx.single_source_bellman_ford_path_length(
            G, next(iter(G)), weight="weight", backend=b), False, {}),
    "floyd_warshall (weighted)": (
        lambda G, b: nx.floyd_warshall(G, weight="weight", backend=b), False,
        {"max_n": 1500}),
    # Centrality extras.
    "load_centrality": (lambda G, b: nx.load_centrality(G, backend=b), False, BETWEENNESS),
    "betweenness_centrality_subset": (
        lambda G, b: nx.betweenness_centrality_subset(
            G, list(G)[:50], list(G)[-50:], backend=b), False, {"max_n": 20000}),
    "current_flow_betweenness_centrality": (
        lambda G, b: nx.current_flow_betweenness_centrality(G, backend=b), False,
        {"max_n": 1000}),
    "second_order_centrality": (
        lambda G, b: nx.second_order_centrality(G, backend=b), False, {"max_n": 1000}),
    # Distance indices.
    "hyper_wiener_index": (lambda G, b: nx.hyper_wiener_index(G, backend=b), False,
                           {"max_n": 20000}),
    "local_efficiency": (lambda G, b: nx.local_efficiency(G, backend=b), False,
                         {"max_n": 3000}),
    # Cores, mixing, cuts, colouring.
    "k_truss": (lambda G, b: nx.k_truss(G, 4, backend=b), False, {}),
    "onion_layers": (lambda G, b: nx.onion_layers(G, backend=b), False, {}),
    "average_neighbor_degree": (lambda G, b: nx.average_neighbor_degree(G, backend=b),
                                False, {}),
    "square_clustering": (lambda G, b: nx.square_clustering(G, backend=b), False,
                          {"max_n": 3000}),
    "conductance": (lambda G, b: nx.conductance(G, list(G)[::2], backend=b), False, {}),
    "greedy_color": (lambda G, b: nx.greedy_color(G, "largest_first", backend=b), False,
                     {}),
    # Similarity and isomorphism.
    "jaccard_coefficient": (lambda G, b: list(nx.jaccard_coefficient(G, backend=b)), False,
                            {"max_n": 3000}),
    "weisfeiler_lehman_graph_hash": (
        lambda G, b: nx.weisfeiler_lehman_graph_hash(G, backend=b), False, {}),
    "simrank_similarity": (lambda G, b: nx.simrank_similarity(G, backend=b), False,
                           {"max_n": 300}),
    "vf2pp_is_isomorphic": (
        lambda G, b: nx.vf2pp_is_isomorphic(G, G.copy(), backend=b), False, {}),
    # Flow and connectivity.
    "maximum_flow_value": (
        lambda G, b: nx.maximum_flow_value(G, next(iter(G)), list(G)[-1],
                                           capacity="weight", backend=b), True, {}),
    "edge_connectivity": (lambda G, b: nx.edge_connectivity(G, backend=b), False,
                          {"max_n": 3000}),
    "gomory_hu_tree": (lambda G, b: nx.gomory_hu_tree(G, capacity="weight", backend=b),
                       False, {"max_n": 300}),
    "condensation": (lambda G, b: nx.condensation(G, backend=b), True, {}),
    # Trees, DAGs and cycles.
    "minimum_spanning_tree (weighted)": (
        lambda G, b: nx.minimum_spanning_tree(G, weight="weight", backend=b), False, {}),
    "cycle_basis": (lambda G, b: nx.cycle_basis(G, backend=b), False, {}),
    "transitive_reduction (DAG)": (
        lambda G, b: nx.transitive_reduction(nx.condensation(G, backend=b), backend=b),
        True, {"max_n": 3000}),
}

# Louvain's partitions differ by design (different random streams); compare
# modularity loosely instead.
LOOSE = {"louvain_communities": 0.02}


def tiers(quick):
    """(label, n, average degree), in increasing number of edges."""
    if quick:
        return [("sparse", 500, 8), ("medium", 1000, 32), ("sparse", 2000, 8)]
    return [("sparse", 1000, 8), ("medium", 1000, 32), ("sparse", 10000, 8),
            ("dense", 1000, 200), ("medium", 10000, 32), ("sparse", 100000, 8),
            ("dense", 3000, 600)]


def fmt(t):
    return f"{t:.3g}"


EXTRA_BACKENDS = ("knappnet", "parallel", "graphblas", "cugraph")


def dispatchable(name):
    """NetworkX's dispatchable function a benchmark entry exercises."""
    base = name.split(" (")[0]
    return getattr(nx, base, None) or getattr(nx.community, base, None)


def implements(backend, name):
    """Whether the backend registers an implementation of the function."""
    if backend == "knappnet":
        return True
    f = dispatchable(name)
    return f is not None and backend in getattr(f, "backends", {})


def backend_converter(backend):
    """The backend's convert_from_nx, where it exposes one."""
    try:
        return nx.utils.backends.backends[backend].load().convert_from_nx
    except Exception:  # not installed, or no conversion exposed
        return None


def installed(backend):
    try:
        return backend in nx.utils.backends.backends and (
            nx.utils.backends.backends[backend].load() is not None)
    except Exception:
        return False


def conversion_time(convert, G, weight):
    if convert is None:
        return None
    try:
        def once():
            G.__networkx_cache__.clear()
            convert(G, edge_attrs={weight: 1} if weight else None)

        return median_time(once)
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--budget", type=float, default=None)
    parser.add_argument("--only", nargs="*", default=None)
    parser.add_argument("--backends", nargs="*", default=list(EXTRA_BACKENDS),
                        help="backends to time next to NetworkX itself "
                             f"(default: {' '.join(EXTRA_BACKENDS)}); one that is not "
                             "installed, or does not implement a function, is reported n/a")
    args = parser.parse_args()
    budget = args.budget or (5.0 if args.quick else 60.0)
    import knappnet  # noqa: F401  (fails early if the package is missing)

    versions = []
    for b in args.backends:
        status = "installed" if installed(b) else "not installed"
        versions.append(f"{b} ({status})")
    print(f"# NetworkX {nx.__version__} ({platform.platform()}, "
          f"Python {platform.python_version()}, {os.cpu_count()} CPUs, "
          f"knappnet threads={nx.config.backends.knappnet.threads or 'all'}); "
          f"NetworkX budget {budget:g} s; backends: {', '.join(versions)}")
    print("| function | graph | n | m | backend | first s | cached s | convert s | "
          "x cached | x first | agree |")
    print("|" + "---|" * 11)
    graphs = {}
    timed_out = {}  # function -> (n, m) of the graph on which NetworkX ran out of time
    failures = 0
    for name, (fn, directed, limits) in FUNCTIONS.items():
        if args.only and not any(o in name for o in args.only):
            continue
        for kind, n, deg in tiers(args.quick):
            if n > limits.get("max_n", n) or deg > limits.get("max_degree", deg):
                continue
            key = (kind, n, directed)
            if graphs.get(key) is None:
                graphs[key] = make_graph(n, deg, directed, seed=n + deg)
            G = graphs[key]
            size = (G.number_of_nodes(), G.number_of_edges())
            weight = "weight" if "weighted" in name or name == "pagerank" else None
            head = f"| {name} | {kind}{' dir.' if directed else ''} | {size[0]} | {size[1]} |"
            G.__networkx_cache__.clear()
            fail = timed_out.get(name)
            if fail and size[0] >= fail[0] and size[1] >= fail[1]:
                t_nx, res_nx = None, None  # larger than a graph NetworkX could not finish
            else:
                t_nx, res_nx = run_limited(lambda: materialize(fn(G, "networkx")), budget)
                if t_nx is None and not fail:
                    timed_out[name] = size
            print(f"{head} networkx | {f'{t_nx:.3g}' if t_nx else f'>{budget:g}'} | - | - "
                  f"| 1x | 1x | - |", flush=True)
            for backend in args.backends:
                if not installed(backend) or not implements(backend, name):
                    print(f"{head} {backend} | n/a | n/a | n/a | n/a | n/a | n/a |", flush=True)
                    continue
                t_conv = conversion_time(backend_converter(backend), G, weight)
                G.__networkx_cache__.clear()
                try:
                    t_first, res = run_limited(lambda: materialize(fn(G, backend)), budget)
                    if t_first is None:
                        print(f"{head} {backend} | >{budget:g} | - | "
                              f"{fmt(t_conv) if t_conv else '-'} | - | - | - |", flush=True)
                        continue
                    t_cached = median_time(lambda: materialize(fn(G, backend)))
                except Exception:  # an argument or input the backend does not take
                    print(f"{head} {backend} | n/a | n/a | n/a | n/a | n/a | n/a |", flush=True)
                    continue
                if t_nx is None:
                    ok = "-"
                    sp_c = f">{budget / t_cached:,.0f}x" if t_cached < budget else "-"
                    sp_f = f">{budget / t_first:,.0f}x" if t_first < budget else "-"
                else:
                    if name in LOOSE:
                        good = abs(res - res_nx) <= LOOSE[name]
                    else:
                        good = agree(res, res_nx)
                    ok = "yes" if good else "NO"
                    failures += not good
                    sp_c = f"{t_nx / t_cached:,.1f}x"
                    sp_f = f"{t_nx / t_first:,.1f}x"
                print(f"{head} {backend} | {fmt(t_first)} | {fmt(t_cached)} | "
                      f"{fmt(t_conv) if t_conv else '-'} | {sp_c} | {sp_f} | {ok} |",
                      flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
