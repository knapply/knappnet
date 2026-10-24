#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for
matchings, greedy colouring, maximal independent sets, square clustering,
and generalized degree (see ../bindings/matching.py for the derived inputs
and encodings).

    python3 validation/oracles/matching_networkx_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture, write_values,  # noqa: E402
                    write_version)

ORACLE = "networkx"


def simple(fx, prefer_max=True):
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        w = fx.weights[i] if fx.weighted else 1.0
        if a == b and not fx.directed:
            g.add_edge(a, b, weight=w)
            continue
        if g.has_edge(a, b):
            old = g[a][b]["weight"]
            w = max(old, w) if prefer_max else min(old, w)
        g.add_edge(a, b, weight=w)
    return g


def natural(fx):
    pairs = list(zip(fx.src, fx.tgt))
    keys = {p if fx.directed else frozenset(p) for p in pairs}
    multi = len(keys) < len(pairs)
    cls = {(True, True): nx.MultiDiGraph, (True, False): nx.DiGraph,
           (False, True): nx.MultiGraph, (False, False): nx.Graph}[(fx.directed, multi)]
    g = cls()
    g.add_nodes_from(range(fx.n))
    g.add_edges_from(pairs)
    return g


def summary(g, m, weighted):
    total = sum(g[u][v]["weight"] if weighted else 1.0 for u, v in m)
    return [(len(m), total, int(nx.is_matching(g, m)))]


def compute(fx, e):
    fn = e["fn"]
    weighted = bool(e.get("weighted"))
    if fn == "max_weight_matching":
        g = simple(fx)
        noloops = g.copy()
        noloops.remove_edges_from(nx.selfloop_edges(noloops))
        m = nx.max_weight_matching(noloops, maxcardinality=bool(e.get("max_cardinality")),
                                   weight="weight" if weighted else None)
        out = summary(g, m, weighted)
        if weighted and not e.get("max_cardinality"):
            out = [(-1,) + out[0][1:]]
        return out
    if fn == "min_weight_matching":
        g = simple(fx, prefer_max=False)
        g.remove_edges_from(nx.selfloop_edges(g))
        m = nx.min_weight_matching(g, weight="weight" if weighted else None)
        return summary(g, m, weighted)
    if fn == "maximal_matching":
        g = simple(fx)
        return [(min(u, v), max(u, v)) for u, v in nx.maximal_matching(g)]
    if fn == "check_matching":
        g = simple(fx)
        if e["pairs"] == "candidate":
            pairs = {(fx.src[i], fx.tgt[i]) for i in range(0, len(fx.src), 2)}
        else:
            pairs = nx.maximal_matching(g)
        return [int(nx.is_matching(g, pairs)), int(nx.is_maximal_matching(g, pairs)),
                int(nx.is_perfect_matching(g, pairs))]
    if fn == "greedy_color":
        c = nx.greedy_color(natural(fx), strategy=e["strategy"])
        return [c[v] for v in range(fx.n)]
    if fn == "mis_check":
        g = simple(fx)
        looped = {a for a, b in zip(fx.src, fx.tgt) if a == b}
        start = [v for v in range(fx.n) if v not in looped][:1]
        got = nx.maximal_independent_set(g, start, seed=11)
        independent = all(not (set(g[v]) & got - {v}) for v in got)
        maximal = all(v in got or set(g[v]) & got for v in g)
        return [(int(independent), int(maximal), int(set(start) <= got))]
    if fn == "square_clustering":
        c = nx.square_clustering(natural(fx))
        return [c[v] for v in range(fx.n)]
    if fn == "generalized_degree":
        g = nx.Graph()
        g.add_nodes_from(range(fx.n))
        g.add_edges_from(zip(fx.src, fx.tgt))
        gd = nx.generalized_degree(g)
        return [tuple(x for k in sorted(gd[v]) for x in (k, gd[v][k])) for v in range(fx.n)]
    return None


def main() -> None:
    spec = load_area("matching")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(fx, v["engine"])
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    print(f"networkx {nx.__version__} (matching): wrote {count} golden files")
    write_version(ORACLE, f"networkx {nx.__version__}")


if __name__ == "__main__":
    main()
