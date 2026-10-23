#!/usr/bin/env python3
"""Benchmarks the engine against NetworkX, python-igraph, and NetworKit,
and checks that every pair of tools being timed gives the same answer.

    python3 bench/bench.py            # full run (NetworkX takes a while)
    python3 bench/bench.py --quick    # small networks, for CI

The engine is called through its C interface (the same code the R package
runs). python-igraph and NetworKit are optional: comparisons with a tool
that is not installed are skipped. Writes bench/results/py-<platform>.csv,
prints a Markdown table, and exits with status 1 if any pair of tools
disagrees.

Timings are medians of repeated runs (operations under 0.2 s are looped so
each timed run lasts at least 0.2 s). Graphs are built before timing; the
cost of building one is timed separately. NetworkX and igraph run on one
core; NetworKit uses every core by default, and is timed that way.
"""

from __future__ import annotations

import csv
import ctypes
import math
import os
import platform
import random
import statistics
import sys
import time
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "validation"))
sys.path.insert(0, str(ROOT / "validation" / "bindings"))
from engine import ERR_LEN, Graph, load  # noqa: E402
import subgroups  # noqa: E402
import topography  # noqa: E402

try:
    import igraph as ig
except ImportError:
    ig = None
try:
    import networkit as nk
except ImportError:
    nk = None

QUICK = "--quick" in sys.argv[1:]
BUDGET = 2.0 if QUICK else 15.0
SCALE = 0.2 if QUICK else 1.0


def timing(f) -> float:
    """Median seconds per call of f()."""
    t0 = time.perf_counter()
    f()
    first = time.perf_counter() - t0
    if first >= 0.2:
        reps = max(0, min(4, int((BUDGET - first) / first)))
        times = [first]
        for _ in range(reps):
            t0 = time.perf_counter()
            f()
            times.append(time.perf_counter() - t0)
        return statistics.median(times)
    inner = max(1, math.ceil(0.2 / max(first, 1e-5)))
    times = []
    for _ in range(5):
        t0 = time.perf_counter()
        for _ in range(inner):
            f()
        times.append((time.perf_counter() - t0) / inner)
    return statistics.median(times)


def close(a, b, rel=1e-8, abs_=1e-10) -> bool:
    a, b = list(a), list(b)
    if len(a) != len(b):
        return False
    for x, y in zip(a, b):
        if math.isnan(x) and math.isnan(y):
            continue
        if not math.isclose(x, y, rel_tol=rel, abs_tol=abs_):
            return False
    return True


def same_partition(a, b) -> bool:
    return len(set(a)) == len(set(b)) == len(set(zip(a, b)))


# ---- Engine calls ---------------------------------------------------------

LIB = load()
topography._declare()
subgroups._lib()
LIB.net_structural_holes.argtypes = [ctypes.c_void_p, ctypes.c_int32, ctypes.c_int32,
                                     ctypes.c_int32, ctypes.POINTER(ctypes.c_double),
                                     ctypes.c_char_p, ctypes.c_size_t]


def _out(n):
    return (ctypes.c_double * max(n, 1))(), ctypes.create_string_buffer(ERR_LEN)


def _ok(status, err):
    if status:
        raise RuntimeError(err.value.decode())


def e_betweenness(g, threads, weighted=False):
    return g.betweenness(weighted=weighted, threads=threads)


MODES = {"out": 1, "in": 2, "all": 3}


def e_closeness(g, threads, harmonic=False, wf=False, mode="all", weighted=False):
    """Closeness (Wasserman-Faust when `wf`, as NetworkX) or harmonic
    centrality. NetworkX measures distances to a node, which is mode "in"
    for directed graphs."""
    out, err = _out(g.n)
    _ok(LIB.net_closeness(g.handle, int(harmonic), MODES[mode], int(weighted), 0, int(wf),
                          1 if wf else 0, threads, out, err, ERR_LEN), err)
    return list(out)[:g.n]


def e_eigenvector(g, threads, unit=False):
    out, err = _out(g.n)
    value, connected = ctypes.c_double(), ctypes.c_int32()
    _ok(LIB.net_eigenvector(g.handle, 1, 0, int(unit), 1 if unit else 0, 0.0, 0, threads, out,
                            ctypes.byref(value), ctypes.byref(connected), err, ERR_LEN), err)
    return list(out)[:g.n]


def e_constraint(g, threads):
    out, err = _out(g.n)
    _ok(LIB.net_structural_holes(g.handle, 0, 0, threads, out, err, ERR_LEN), err)
    return list(out)[:g.n]


def e_components(g):
    member = (ctypes.c_int32 * max(g.n, 1))()
    count = ctypes.c_int64()
    err = ctypes.create_string_buffer(ERR_LEN)
    _ok(LIB.net_components(g.handle, 0, member, ctypes.byref(count), err, ERR_LEN), err)
    return list(member)[:g.n]


# ---- Graphs ---------------------------------------------------------------

def gnm(n, m, directed, seed):
    rng = random.Random(seed)
    edges = set()
    while len(edges) < m:
        a, b = rng.randrange(n), rng.randrange(n)
        if a == b:
            continue
        if not directed and a > b:
            a, b = b, a
        edges.add((a, b))
    return sorted(edges)


def preferential(n, k, seed):
    """Barabasi-Albert preferential attachment, k edges per new node."""
    return sorted(tuple(sorted(e)) for e in nx.barabasi_albert_graph(n, k, seed=seed).edges())


def weights_for(m, integer, seed):
    """Whole numbers 1 to 5 (the engine's Dial path) or reals in [0.1, 3]
    (Dijkstra's algorithm)."""
    rng = random.Random(seed)
    return [float(rng.randint(1, 5)) if integer else rng.uniform(0.1, 3.0) for _ in range(m)]


class Net:
    def __init__(self, label, n, edges, directed, weights=None):
        self.label, self.n, self.edges, self.directed = label, n, edges, directed
        src, tgt = [a for a, _ in edges], [b for _, b in edges]
        self.g = Graph(n, src, tgt, weights, directed)
        self.nx = nx.DiGraph() if directed else nx.Graph()
        self.nx.add_nodes_from(range(n))
        if weights is None:
            self.nx.add_edges_from(edges)
        else:
            self.nx.add_weighted_edges_from((a, b, w) for (a, b), w in zip(edges, weights))
        self.ig = None
        if ig:
            self.ig = ig.Graph(n=n, edges=edges, directed=directed,
                               edge_attrs={"weight": weights} if weights is not None else None)
        self.nk = None
        if nk and weights is None:
            self.nk = nk.Graph(n, directed=directed)
            for a, b in edges:
                self.nk.addEdge(a, b)


rows = []
notes = {}


def record(net, measure, ours1, ours, tool, theirs, agree, note=None):
    """`note` (optional) is printed after the row and again if the tools
    disagree, for example both modularity values for Louvain."""
    rows.append({"graph": net.label, "measure": measure, "knappnet_1_thread": ours1,
                 "knappnet": ours, "tool": tool, "tool_seconds": theirs,
                 "speedup_1_thread": theirs / ours1, "speedup": theirs / ours,
                 "agree": agree})
    if note:
        notes[(net.label, measure, tool)] = note
    print(f"{net.label:<24} {measure:<24} {tool:<9} {ours1:9.4f} {ours:9.4f} {theirs:9.4f} "
          f"{theirs / ours1:7.1f}x {theirs / ours:7.1f}x {'agree' if agree else 'DISAGREE'}"
          f"{'  ' + note if note else ''}", flush=True)


def both(f):
    return timing(lambda: f(1)), timing(lambda: f(0))


def compare(net, measure, ours_fn, threaded, tools, note=None):
    """tools: list of (name, available, run, check) where check(ours, theirs).
    note(ours, theirs), if given, describes the two results for the output."""
    one, all_ = both(ours_fn) if threaded else (timing(lambda: ours_fn(0)),) * 2
    ours = ours_fn(0)
    for name, available, run, check in tools:
        if not available:
            continue
        theirs_t = timing(run)
        theirs = run()
        record(net, measure, one, all_, name, theirs_t, bool(check(ours, theirs)),
               note(ours, theirs) if note else None)


def main() -> int:
    cpu = platform.processor() or platform.machine()
    if sys.platform.startswith("linux"):
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    print(f"Python {platform.python_version()}; {cpu}; {os.cpu_count()} cores; "
          f"networkx {nx.__version__}; igraph {ig.__version__ if ig else 'not installed'}; "
          f"networkit {nk.__version__ if nk else 'not installed'}")
    if nk:
        nk.engineering.setNumberOfThreads(os.cpu_count() or 1)

    n_mid = int(5000 * SCALE)
    n_big = int(100000 * SCALE)
    und = Net(f"G(n={n_mid}, m={5 * n_mid})", n_mid, gnm(n_mid, 5 * n_mid, False, 1), False)
    dire = Net(f"directed G({n_mid}, {5 * n_mid})", n_mid, gnm(n_mid, 5 * n_mid, True, 2), True)
    pa_edges = preferential(n_big, 3, 3)
    pa = Net(f"scale-free n={n_big}, m={len(pa_edges)}", n_big, pa_edges, False)

    print(f"{'graph':<24} {'measure':<24} {'tool':<9} {'1 thread':>9} {'all':>9} {'tool':>9}")

    # Building a network from an edge list.
    src, tgt = [a for a, _ in pa.edges], [b for _, b in pa.edges]
    ours = timing(lambda: Graph(pa.n, src, tgt, None, False).degree())

    def build_nx():
        G = nx.Graph()
        G.add_nodes_from(range(pa.n))
        G.add_edges_from(pa.edges)
        return [d for _, d in G.degree()]
    record(pa, "build + degree", ours, ours, "networkx", timing(build_nx),
           close(Graph(pa.n, src, tgt, None, False).degree(), build_nx()))
    if ig:
        def build_ig():
            return ig.Graph(n=pa.n, edges=pa.edges).degree()
        record(pa, "build + degree", ours, ours, "igraph", timing(build_ig),
               close(Graph(pa.n, src, tgt, None, False).degree(), build_ig()))

    for net in (und, dire):
        compare(net, "betweenness", lambda k, net=net: e_betweenness(net.g, k), True, [
            ("networkx", True,
             lambda net=net: list(nx.betweenness_centrality(net.nx, normalized=False).values()),
             close),
            ("igraph", ig is not None, lambda net=net: net.ig.betweenness(), close),
            ("networkit", nk is not None,
             lambda net=net: nk.centrality.Betweenness(net.nk).run().scores(),
             # NetworKit counts each unordered pair twice in undirected graphs.
             lambda o, t, net=net: close(o, [x / (1 if net.directed else 2) for x in t])),
        ])

    compare(und, "closeness", lambda k: e_closeness(und.g, k, wf=True), True, [
        ("networkx", True, lambda: list(nx.closeness_centrality(und.nx).values()), close),
    ])
    compare(und, "harmonic centrality", lambda k: e_closeness(und.g, k, harmonic=True), True, [
        ("networkx", True, lambda: list(nx.harmonic_centrality(und.nx).values()), close),
        ("igraph", ig is not None,
         lambda: und.ig.harmonic_centrality(normalized=False), close),
        ("networkit", nk is not None,
         lambda: nk.centrality.HarmonicCloseness(und.nk, normalized=False).run().scores(), close),
    ])

    # Directed closeness and harmonic centrality: NetworkX measures
    # distances to each node (the engine's mode "in").
    compare(dire, "closeness", lambda k: e_closeness(dire.g, k, wf=True, mode="in"), True, [
        ("networkx", True, lambda: list(nx.closeness_centrality(dire.nx).values()), close),
    ])
    compare(dire, "harmonic centrality",
            lambda k: e_closeness(dire.g, k, harmonic=True, mode="in"), True, [
                ("networkx", True, lambda: list(nx.harmonic_centrality(dire.nx).values()), close),
                ("igraph", ig is not None,
                 lambda: dire.ig.harmonic_centrality(mode="in", normalized=False), close),
            ])

    # Weighted shortest paths: whole-number weights take the engine's Dial
    # bucket queue, real weights Dijkstra's algorithm.
    n_w = int(2000 * SCALE)
    for directed in (False, True):
        for integer in (True, False):
            edges = gnm(n_w, 5 * n_w, directed, 10 + 2 * directed + integer)
            w = weights_for(len(edges), integer, 20 + 2 * directed + integer)
            kind = "integer (Dial)" if integer else "real (Dijkstra)"
            net = Net(f"{'directed, ' if directed else ''}{kind} weights, {n_w} / {len(edges)}",
                      n_w, edges, directed, w)
            mode = "in" if directed else "all"
            compare(net, "betweenness", lambda k, net=net: e_betweenness(net.g, k, True), True, [
                ("networkx", True,
                 lambda net=net: list(nx.betweenness_centrality(
                     net.nx, normalized=False, weight="weight").values()), close),
                ("igraph", ig is not None,
                 lambda net=net: net.ig.betweenness(weights="weight"), close),
            ])
            compare(net, "closeness",
                    lambda k, net=net, mode=mode: e_closeness(net.g, k, wf=True, mode=mode,
                                                              weighted=True), True, [
                        ("networkx", True,
                         lambda net=net: list(nx.closeness_centrality(
                             net.nx, distance="weight").values()), close),
                    ])
            compare(net, "harmonic centrality",
                    lambda k, net=net, mode=mode: e_closeness(net.g, k, harmonic=True, mode=mode,
                                                              weighted=True), True, [
                        ("networkx", True,
                         lambda net=net: list(nx.harmonic_centrality(
                             net.nx, distance="weight").values()), close),
                        ("igraph", ig is not None,
                         lambda net=net, mode=mode: net.ig.harmonic_centrality(
                             mode=mode, weights="weight", normalized=False), close),
                    ])

    def unit(v):
        s = math.sqrt(sum(x * x for x in v))
        return [abs(x) / s for x in v]
    compare(pa, "eigenvector centrality", lambda k: e_eigenvector(pa.g, k, unit=True), True, [
        ("networkx", True,
         lambda: list(nx.eigenvector_centrality_numpy(pa.nx).values()),
         lambda o, t: close(o, unit(t), 1e-6, 1e-9)),
        ("igraph", ig is not None, lambda: pa.ig.eigenvector_centrality(scale=True),
         lambda o, t: close(o, unit(t), 1e-6, 1e-9)),
        # NetworKit's power iteration stops at a looser tolerance.
        ("networkit", nk is not None,
         lambda: nk.centrality.EigenvectorCentrality(pa.nk).run().scores(),
         lambda o, t: close(o, unit(t), 1e-3, 1e-6)),
    ])

    small_n = int(2000 * SCALE)
    small = Net(f"G(n={small_n}, m={5 * small_n})", small_n,
                gnm(small_n, 5 * small_n, False, 4), False)
    compare(small, "constraint", lambda k: e_constraint(small.g, k), True, [
        ("networkx", True, lambda: list(nx.constraint(small.nx).values()), close),
        ("igraph", ig is not None, lambda: small.ig.constraint(), close),
    ])

    compare(pa, "k-cores", lambda k: subgroups.coreness(pa.g), False, [
        ("networkx", True, lambda: list(nx.core_number(pa.nx).values()), close),
        ("igraph", ig is not None, lambda: pa.ig.coreness(), close),
        ("networkit", nk is not None,
         lambda: nk.centrality.CoreDecomposition(pa.nk).run().scores(), close),
    ])

    def nx_components():
        label = [0] * pa.n
        for i, comp in enumerate(nx.connected_components(pa.nx)):
            for v in comp:
                label[v] = i
        return label
    compare(pa, "components", lambda k: e_components(pa.g), False, [
        ("networkx", True, nx_components, same_partition),
        ("igraph", ig is not None, lambda: pa.ig.connected_components().membership,
         same_partition),
        ("networkit", nk is not None,
         lambda: nk.components.ConnectedComponents(pa.nk).run().getPartition().getVector(),
         same_partition),
    ])

    compare(pa, "maximal cliques", lambda k: subgroups.max_cliques(pa.g, 1, k), True, [
        ("networkx", True, lambda: list(nx.find_cliques(pa.nx)),
         lambda o, t: len(o) == len(t)),
        ("igraph", ig is not None, lambda: pa.ig.maximal_cliques(),
         lambda o, t: len(o) == len(t)),
        ("networkit", nk is not None,
         lambda: nk.clique.MaximalCliques(pa.nk).run().getCliques(),
         lambda o, t: len(o) == len(t)),
    ])

    # Louvain is a randomized heuristic, so tools find different partitions.
    # Agreement means comparable quality: the engine's modularity is at least
    # 98% of the other tool's (finding a higher modularity is fine).
    q_ours = subgroups.louvain(pa.g, 1.0, False, 1)[1]

    def near(q):
        return q_ours >= 0.98 * q
    compare(pa, "Louvain communities", lambda k: subgroups.louvain(pa.g, 1.0, False, 1), False, [
        ("networkx", True,
         lambda: nx.community.modularity(pa.nx, nx.community.louvain_communities(pa.nx, seed=1)),
         lambda o, t: near(t)),
        ("igraph", ig is not None, lambda: pa.ig.community_multilevel().modularity,
         lambda o, t: near(t)),
        ("networkit", nk is not None,
         lambda: nk.community.Modularity().getQuality(
             nk.community.detectCommunities(pa.nk, algo=nk.community.PLM(pa.nk, False),
                                            inspect=False), pa.nk),
         lambda o, t: near(t)),
    ], note=lambda o, t: f"modularity {q_ours:.4f} (knappnet) vs {t:.4f}")

    names = nx.algorithms.triads.TRIAD_NAMES
    compare(dire, "triad census", lambda k: dire.g.triad_census(threads=k), True, [
        ("networkx", True, lambda: nx.triadic_census(dire.nx),
         lambda o, t: close(o, [t[name] for name in names])),
        ("igraph", ig is not None, lambda: list(dire.ig.triad_census()), close),
    ])

    out_dir = ROOT / "bench" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    name = f"py-{sys.platform}{'-quick' if QUICK else ''}.csv"
    with open(out_dir / name, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    def fs(s):
        if s >= 1:
            return f"{s:.2f} s"
        if s >= 1e-3:
            return f"{s * 1e3:.1f} ms"
        return f"{s * 1e6:.0f} µs"

    def fx(x):
        return f"{x:.0f}x" if x >= 10 else f"{x:.1f}x"
    print("\n| Network | Measure | Compared with | knappnet, 1 thread | knappnet, all cores | "
          "Other tool | Speedup, 1 thread | Speedup, all cores |")
    print("|---|---|---|---:|---:|---:|---:|---:|")
    for r in rows:
        print(f"| {r['graph']} | {r['measure']} | {r['tool']} | {fs(r['knappnet_1_thread'])} | "
              f"{fs(r['knappnet'])} | {fs(r['tool_seconds'])} | {fx(r['speedup_1_thread'])} | "
              f"{fx(r['speedup'])} |")
    bad = [r for r in rows if not r["agree"]]
    if bad:
        print("\nResults differ for:")
        for r in bad:
            note = notes.get((r["graph"], r["measure"], r["tool"]))
            print(f"  {r['graph']}: {r['measure']} vs {r['tool']}" + (f" ({note})" if note else ""))
        return 1
    print("\nAll results agree.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
