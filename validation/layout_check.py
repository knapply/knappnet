#!/usr/bin/env python3
"""Check network layouts by their properties.

A layout is not unique (any rotation, reflection, or translation of a good
layout is equally good, and different local optima can be equally
acceptable), so layouts cannot be compared value by value with another
tool. This script checks what must hold instead, on every fixture in
validation/fixtures/:

1. Finite coordinates for every method.
2. Reproducibility: identical coordinates, bit for bit, from repeated calls
   with the same seed and from 1, 2, and 3 threads.
3. Separation: in stress, force-directed, and tree layouts, nodes of
   different weak components are at least the packing gap apart (the mean
   shortest edge length between distinct adjacent nodes for stress, 1
   otherwise).
4. Quality: on the largest component of each fixture, the stress of the
   engine's default layout is no worse than that of NetworkX's
   kamada_kawai_layout(), which minimizes the same energy (Kamada and Kawai
   1989), up to a 2% allowance. Both are measured here with the same
   formula, sum over pairs of (||x_i - x_j|| - d_ij)^2 / d_ij^2 with d from
   NetworkX's shortest paths, after scaling each layout by the factor that
   minimizes its stress (NetworkX rescales its output to [-1, 1]). The
   engine's own stress function must agree with this formula.
5. Sparse stress against stress on all pairs, reported as the ratio of
   their stresses: on fixtures with n/5 pivots (a deliberately coarse
   approximation; with every node a pivot the two models coincide), and on
   larger generated graphs with the default 100 pivots.
6. Simple layouts: the circle lies on the unit circle, tree edges join
   consecutive levels in trees, two-mode rows follow the node types.

    python3 validation/layout_check.py

Prints the numbers and exits non-zero on any failure.
"""

from __future__ import annotations

import ctypes
import math
import sys
import time
from pathlib import Path

import networkx as nx
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import fixture_names, load_fixture  # noqa: E402
from engine import ERR_LEN, EngineError, Graph, load  # noqa: E402

# The engine's stress must be no more than this factor times NetworkX's.
ALLOWANCE = 1.02

failures: list[str] = []


def check(ok: bool, what: str) -> None:
    if not ok:
        failures.append(what)
        print(f"FAIL  {what}")


# ---- Engine binding -------------------------------------------------------

def _bind():
    lib = load()
    h, i32, i64, u64, f64 = (ctypes.c_void_p, ctypes.c_int32, ctypes.c_int64,
                             ctypes.c_uint64, ctypes.c_double)
    f64p = ctypes.POINTER(ctypes.c_double)
    cp, sz = ctypes.c_char_p, ctypes.c_size_t
    lib.net_layout_stress.argtypes = [h, i32, i32, i32, i64, i64, f64, u64, i32, f64p, f64p,
                                      cp, sz]
    lib.net_layout_fr.argtypes = [h, i32, i32, i64, i32, f64, u64, i32, f64p, f64p, cp, sz]
    lib.net_layout_simple.argtypes = [h, i32, i64, f64p, f64p, cp, sz]
    lib.net_layout_tree.argtypes = [h, ctypes.POINTER(i64), i64, f64p, f64p, cp, sz]
    lib.net_layout_twomode.argtypes = [h, ctypes.POINTER(i32), i64, f64p, f64p, cp, sz]
    lib.net_layout_stress_value.argtypes = [h, f64p, f64p, i32, i32, i32, f64p, cp, sz]
    for fn in ("net_layout_stress", "net_layout_fr", "net_layout_simple", "net_layout_tree",
               "net_layout_twomode", "net_layout_stress_value"):
        getattr(lib, fn).restype = i32
    return lib


LIB = None


def _run(fn, g: Graph, *args) -> np.ndarray:
    x = (ctypes.c_double * max(g.n, 1))()
    y = (ctypes.c_double * max(g.n, 1))()
    err = ctypes.create_string_buffer(ERR_LEN)
    if fn(g.handle, *args, x, y, err, ERR_LEN):
        raise EngineError(err.value.decode())
    return np.column_stack([np.array(x[:g.n]), np.array(y[:g.n])]).reshape(g.n, 2)


def stress_layout(g, algorithm=0, weighted=False, pivots=0, max_iter=0, tol=0.0, seed=1,
                  threads=0):
    return _run(LIB.net_layout_stress, g, algorithm, int(weighted), 0, pivots, max_iter, tol,
                seed, threads)


def fr_layout(g, weighted=False, iterations=0, approximation=0, seed=1, threads=0):
    return _run(LIB.net_layout_fr, g, int(weighted), 0, iterations, approximation, 0.0, seed,
                threads)


def simple_layout(g, kind, param):
    return _run(LIB.net_layout_simple, g, kind, param)


def tree_layout(g, roots=()):
    r = (ctypes.c_int64 * max(len(roots), 1))(*roots)
    return _run(LIB.net_layout_tree, g, r, len(roots))


def twomode_layout(g, types, sweeps=0):
    t = (ctypes.c_int32 * max(len(types), 1))(*types)
    return _run(LIB.net_layout_twomode, g, t, sweeps)


def engine_stress(g, xy, weighted=False) -> float:
    x = (ctypes.c_double * max(g.n, 1))(*xy[:, 0])
    y = (ctypes.c_double * max(g.n, 1))(*xy[:, 1])
    out = (ctypes.c_double * 1)()
    err = ctypes.create_string_buffer(ERR_LEN)
    if LIB.net_layout_stress_value(g.handle, x, y, int(weighted), 0, 0, out, err, ERR_LEN):
        raise EngineError(err.value.decode())
    return out[0]


# ---- Independent reference computations -----------------------------------

def skeleton(n, src, tgt, weights) -> nx.Graph:
    """Undirected simple graph; parallel edges keep their smallest weight and
    self-loops are dropped, as the engine's layouts read a network."""
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for e, (a, b) in enumerate(zip(src, tgt)):
        if a == b:
            continue
        w = 1.0 if weights is None else weights[e]
        if G.has_edge(a, b):
            G[a][b]["weight"] = min(G[a][b]["weight"], w)
        else:
            G.add_edge(a, b, weight=w)
    return G


def distance_matrix(G: nx.Graph, weighted: bool) -> np.ndarray:
    nodes = list(G.nodes())
    index = {v: i for i, v in enumerate(nodes)}
    D = np.full((len(nodes), len(nodes)), np.inf)
    lengths = (nx.all_pairs_dijkstra_path_length(G, weight="weight") if weighted
               else nx.all_pairs_shortest_path_length(G))
    for a, row in lengths:
        for b, d in row.items():
            D[index[a], index[b]] = d
    return D


def stress(xy: np.ndarray, D: np.ndarray, rescale: bool) -> float:
    """sum over pairs i < j with finite d of ((e - d) / d)^2, after optimal
    scaling when `rescale`."""
    iu = np.triu_indices(len(xy), 1)
    d = D[iu]
    e = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1))[iu]
    ok = np.isfinite(d) & (d > 0)
    d, e = d[ok], e[ok]
    if rescale and (e * e / (d * d)).sum() > 0:
        e = e * (e / d).sum() / (e * e / (d * d)).sum()
    return float((((e - d) / d) ** 2).sum())


def largest_component(fx):
    G = skeleton(fx.n, fx.src, fx.tgt, fx.weights)
    keep = sorted(max(nx.connected_components(G), key=lambda c: (len(c), -min(c))))
    index = {v: i for i, v in enumerate(keep)}
    src, tgt, w = [], [], []
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a in index and b in index:
            src.append(index[a])
            tgt.append(index[b])
            if fx.weights is not None:
                w.append(fx.weights[e])
    return len(keep), src, tgt, (w if fx.weights is not None else None)


# ---- Checks -----------------------------------------------------------------

def methods(g, weighted):
    yield "stress", lambda t, s=1: stress_layout(g, 0, weighted, seed=s, threads=t)
    yield "sparse stress", lambda t, s=1: stress_layout(g, 2, weighted, pivots=8, seed=s,
                                                        threads=t)
    yield "fr exact", lambda t, s=1: fr_layout(g, weighted, 100, 1, seed=s, threads=t)
    yield "fr barnes-hut", lambda t, s=1: fr_layout(g, weighted, 100, 2, seed=s, threads=t)
    yield "tree", lambda t, s=1: tree_layout(g)


def check_fixtures() -> None:
    print("== Finite, reproducible, thread-independent, components separated")
    for name in fixture_names():
        fx = load_fixture(name)
        for weighted in ([False, True] if fx.weighted else [False]):
            if weighted and min(fx.weights, default=1.0) <= 0:
                continue
            g = Graph(fx.n, fx.src, fx.tgt, fx.weights if weighted else None, fx.directed)
            G = skeleton(fx.n, fx.src, fx.tgt, fx.weights if weighted else None)
            comp = {}
            for k, c in enumerate(nx.connected_components(G)):
                for v in c:
                    comp[v] = k
            lens = [d["weight"] for _, _, d in G.edges(data=True)]
            stress_gap = float(np.mean(lens)) if lens else 1.0
            for method, run in methods(g, weighted):
                tag = f"{name}{' weighted' if weighted else ''} {method}"
                a = run(1)
                check(bool(np.isfinite(a).all()), f"{tag}: finite coordinates")
                check(np.array_equal(a, run(1)), f"{tag}: repeated call")
                check(np.array_equal(a, run(2)), f"{tag}: 2 threads")
                check(np.array_equal(a, run(3)), f"{tag}: 3 threads")
                gap = stress_gap if "stress" in method else 1.0
                worst = math.inf
                for u in range(fx.n):
                    for v in range(u + 1, fx.n):
                        if comp[u] != comp[v]:
                            worst = min(worst, float(np.hypot(*(a[u] - a[v]))))
                check(worst >= gap * (1 - 1e-9),
                      f"{tag}: components {worst:.4g} apart, gap {gap:.4g}")
            if fx.n >= 4 and G.number_of_edges() >= 3:
                fa, fb = fr_layout(g, weighted, 50, 0, seed=1), fr_layout(g, weighted, 50, 0,
                                                                         seed=2)
                check(not np.array_equal(fa, fb), f"{name}: seeds give different FR starts")
    print(f"   {len(fixture_names())} fixtures checked")


def check_quality() -> None:
    print("== Stress on the largest component (lower is better)")
    print(f"   {'fixture':<24}{'n':>5}{'engine':>12}{'networkx kk':>14}{'ratio':>8}"
          f"{'sparse k=n/5':>15}")
    for name in fixture_names():
        fx = load_fixture(name)
        n, src, tgt, w = largest_component(fx)
        if n < 3:
            continue
        for weighted in ([False, True] if fx.weighted else [False]):
            wts = w if weighted else None
            g = Graph(n, src, tgt, wts, fx.directed)
            G = skeleton(n, src, tgt, wts)
            D = distance_matrix(G, weighted)
            ours = stress_layout(g, 0, weighted)
            mine = stress(ours, D, rescale=False)
            reported = engine_stress(g, ours, weighted)
            check(abs(mine - reported) <= 1e-9 * max(1.0, mine),
                  f"{name}: engine stress {reported} vs formula {mine}")
            kk = nx.kamada_kawai_layout(G, weight="weight" if weighted else None)
            kxy = np.array([kk[v] for v in range(n)])
            s_ours = stress(ours, D, rescale=True)
            s_kk = stress(kxy, D, rescale=True)
            sparse = stress_layout(g, 2, weighted, pivots=max(3, n // 5))
            s_sparse = stress(sparse, D, rescale=True)
            ratio = s_ours / s_kk if s_kk > 0 else (1.0 if s_ours < 1e-12 else math.inf)
            label = name + (" (w)" if weighted else "")
            print(f"   {label:<24}{n:>5}{s_ours:>12.4f}{s_kk:>14.4f}{ratio:>8.3f}"
                  f"{(s_sparse / s_ours if s_ours > 0 else 1.0):>15.3f}")
            check(s_ours <= ALLOWANCE * s_kk + 1e-9,
                  f"{label}: stress {s_ours:.6g} worse than NetworkX {s_kk:.6g}")


def check_sparse_on_larger_graphs() -> None:
    print("== Sparse stress (100 pivots) against stress on all pairs, generated graphs")
    graphs = {
        "grid 40x40": nx.convert_node_labels_to_integers(nx.grid_2d_graph(40, 40)),
        "watts-strogatz 1500": nx.connected_watts_strogatz_graph(1500, 6, 0.05, seed=3),
        "barabasi-albert 1500": nx.barabasi_albert_graph(1500, 2, seed=3),
        "random tree 1500": nx.random_labeled_tree(1500, seed=3),
    }
    for label, G in graphs.items():
        n = G.number_of_nodes()
        src = [a for a, b in G.edges()]
        tgt = [b for a, b in G.edges()]
        g = Graph(n, src, tgt, None, False)
        D = distance_matrix(G, False)
        t0 = time.perf_counter()
        full = stress_layout(g, 1)
        t1 = time.perf_counter()
        sparse = stress_layout(g, 2)
        t2 = time.perf_counter()
        s_full, s_sparse = stress(full, D, True), stress(sparse, D, True)
        print(f"   {label:<22} n={n:<5} all pairs {s_full:12.1f} ({t1 - t0:.2f} s)   "
              f"sparse {s_sparse:12.1f} ({t2 - t1:.2f} s)   ratio {s_sparse / s_full:.3f}")
        check(bool(np.isfinite(full).all() and np.isfinite(sparse).all()),
              f"{label}: finite coordinates")
        check(s_sparse <= 1.5 * s_full, f"{label}: sparse stress far worse than all pairs")


def check_simple_layouts() -> None:
    print("== Circle, star, grid, tree, and two-mode layouts")
    for name in fixture_names():
        fx = load_fixture(name)
        g = Graph(fx.n, fx.src, fx.tgt, None, fx.directed)
        c = simple_layout(g, 0, 0)
        if fx.n > 1:
            check(bool(np.allclose(np.hypot(c[:, 0], c[:, 1]), 1.0)), f"{name}: circle radius")
        s = simple_layout(g, 1, -1)
        r = np.hypot(s[:, 0], s[:, 1])
        check(int((r < 1e-12).sum()) == min(fx.n, 1) and
              bool(np.allclose(r[r > 1e-12], 1.0)), f"{name}: star")
        grid = simple_layout(g, 2, 0)
        check(bool(np.isfinite(grid).all()), f"{name}: grid")
        G = skeleton(fx.n, fx.src, fx.tgt, None)
        t = tree_layout(g)
        if fx.n > 0 and nx.is_tree(G):
            check(all(abs(abs(t[a, 1] - t[b, 1]) - 1) < 1e-9 for a, b in G.edges()),
                  f"{name}: tree edges join consecutive levels")
        types = [v % 2 for v in range(fx.n)]
        two = twomode_layout(g, types)
        top = [two[v, 1] for v in range(fx.n) if types[v] == 0]
        bottom = [two[v, 1] for v in range(fx.n) if types[v] == 1]
        if top and bottom:
            check(min(top) > max(bottom), f"{name}: two-mode rows")
    # A two-mode network with a crossing-free drawing must get one.
    B = nx.Graph([(0, 5), (1, 4), (2, 3)])
    g = Graph(6, [0, 1, 2], [5, 4, 3], None, False)
    two = twomode_layout(g, [0, 0, 0, 1, 1, 1])
    crossings = sum(1 for (a, b) in B.edges() for (c, d) in B.edges()
                    if (two[a, 0] - two[c, 0]) * (two[b, 0] - two[d, 0]) < 0) // 2
    check(crossings == 0, f"two-mode barycenter ordering leaves {crossings} crossings")


def main() -> int:
    global LIB
    LIB = _bind()
    check_fixtures()
    check_simple_layouts()
    check_quality()
    check_sparse_on_larger_graphs()
    print(f"{len(failures)} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    if {"-h", "--help"} & set(sys.argv[1:]):
        print(__doc__.strip())
        sys.exit(0)
    sys.exit(main())
