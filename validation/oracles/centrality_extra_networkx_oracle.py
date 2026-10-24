#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for the
centrality variants: betweenness for sources and targets and from sampled
sources, percolation and load centrality, edge load, current-flow
betweenness (exact, subset, approximate) and closeness, resistance
distance, effective graph resistance, Kemeny's constant, second-order
centrality, communicability and communicability betweenness, dispersion,
VoteRank, Laplacian centrality, trophic levels, and reaching centrality.

Networks are MultiGraph or MultiDiGraph objects with the fixture weights in
"weight" (simple Graph objects for the functions NetworkX defines only for
graphs; the manifest marks multigraph fixtures there). Node sets, states,
dispersion pairs, and sampled sources repeat validation/bindings/
centrality_extra.py. A NetworkX error is written as -999, like NaN; Inf is
-1000.

    python3 validation/oracles/centrality_extra_networkx_oracle.py
"""

from __future__ import annotations

import math
import random
import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture, write_values,  # noqa: E402
                    write_version)

ORACLE = "networkx"
UNDEFINED = -999.0
INFINITE = -1000.0


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


def mark(x):
    x = float(x)
    if math.isnan(x):
        return UNDEFINED
    if math.isinf(x):
        return INFINITE
    return x


def multi(fx, weighted=True):
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if weighted and fx.weighted:
            g.add_edge(a, b, weight=fx.weights[i])
        else:
            g.add_edge(a, b)
    return g


def simple(fx, weighted=True):
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if weighted and fx.weighted:
            g.add_edge(a, b, weight=fx.weights[i])
        else:
            g.add_edge(a, b)
    return g


def edge_keys(fx):
    """The (u, v, key) NetworkX gives each input edge of a multigraph."""
    seen = {}
    out = []
    for a, b in zip(fx.src, fx.tgt):
        pair = (a, b) if fx.directed else (min(a, b), max(a, b))
        k = seen.get(pair, 0)
        seen[pair] = k + 1
        out.append((a, b, k))
    return out


def per_edge(fx, result, multigraph):
    out = []
    keys = edge_keys(fx)
    for (a, b, k) in keys:
        if multigraph:
            cands = [(a, b, k)] if fx.directed else [(a, b, k), (b, a, k)]
        else:
            cands = [(a, b)] if fx.directed else [(a, b), (b, a)]
        val = next((result[c] for c in cands if c in result), None)
        out.append(mark(val) if val is not None else 0.0)
    return out


def guarded(fn):
    try:
        return fn()
    except (nx.NetworkXError, nx.NetworkXException, ZeroDivisionError):
        return None


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    n = fx.n
    weight = "weight" if e.get("weighted") else None
    if fn == "cx_betweenness_subset":
        g = multi(fx)
        normalized = e.get("normalized", False)
        edges = e.get("output", "nodes") == "edges"
        if "sample" in e:
            k = min(int(e["sample"]), n)
            seed = int(e.get("seed", 1))
            if edges:
                r = nx.edge_betweenness_centrality(g, k=k, normalized=normalized, weight=weight,
                                                   seed=seed)
            else:
                r = nx.betweenness_centrality(g, k=k, normalized=normalized, weight=weight,
                                              endpoints=e.get("endpoints", False), seed=seed)
        elif e.get("rescale") == "networkx_subset":
            src = node_subset(n, e.get("sources", "all"))
            tgt = node_subset(n, e.get("targets", "all"))
            f = (nx.edge_betweenness_centrality_subset if edges
                 else nx.betweenness_centrality_subset)
            r = f(g, src, tgt, normalized=normalized, weight=weight)
        else:
            if edges:
                r = nx.edge_betweenness_centrality(g, normalized=normalized, weight=weight)
            else:
                r = nx.betweenness_centrality(g, normalized=normalized, weight=weight,
                                              endpoints=e.get("endpoints", False))
        if edges:
            return per_edge(fx, r, True)
        return [mark(r[v]) for v in range(n)]
    if fn == "cx_percolation":
        g = multi(fx)
        r = guarded(lambda: nx.percolation_centrality(
            g, states=dict(enumerate(states(n))), weight=weight))
        return [UNDEFINED] * n if r is None else [mark(r[v]) for v in range(n)]
    if fn == "cx_load":
        g = multi(fx)
        cutoff = e.get("cutoff")
        r = nx.load_centrality(g, normalized=e.get("normalized", False), weight=weight,
                               cutoff=cutoff)
        return [mark(r[v]) for v in range(n)]
    if fn == "cx_edge_load":
        g = multi(fx)
        r = nx.edge_load_centrality(g, cutoff=e.get("cutoff", False))
        return [(mark(r[(a, b)]), mark(r[(b, a)])) for a, b in zip(fx.src, fx.tgt)]
    if fn in ("cx_current_flow", "cx_current_flow_subset", "cx_current_flow_approx",
              "cx_current_flow_closeness"):
        g = simple(fx)
        normalized = e.get("normalized", True)
        edges = e.get("output", "nodes") == "edges"
        if fn == "cx_current_flow":
            if edges:
                r = nx.edge_current_flow_betweenness_centrality(g, normalized=normalized,
                                                                weight=weight)
                return per_edge(fx, r, False)
            r = nx.current_flow_betweenness_centrality(g, normalized=normalized, weight=weight)
        elif fn == "cx_current_flow_subset":
            src = node_subset(n, e.get("sources", "even"))
            tgt = node_subset(n, e.get("targets", "thirds"))
            if edges:
                r = nx.edge_current_flow_betweenness_centrality_subset(
                    g, src, tgt, normalized=normalized, weight=weight)
                return per_edge(fx, r, False)
            r = guarded(lambda: nx.current_flow_betweenness_centrality_subset(
                g, src, tgt, normalized=normalized, weight=weight))
            if r is None:
                return [UNDEFINED] * n
        elif fn == "cx_current_flow_approx":
            r = nx.approximate_current_flow_betweenness_centrality(
                g, normalized=normalized, weight=weight, epsilon=float(e.get("epsilon", 0.5)),
                seed=int(e.get("seed", 1)))
        else:
            r = nx.current_flow_closeness_centrality(g, weight=weight)
        return [mark(r[v]) for v in range(n)]
    if fn == "cx_resistance":
        g = multi(fx)
        r = nx.resistance_distance(g, weight=weight, invert_weight=True)
        return [tuple(mark(r[i][j]) for j in range(n)) for i in range(n)]
    if fn == "cx_effective_resistance":
        g = multi(fx)
        return [mark(nx.effective_graph_resistance(g, weight=weight, invert_weight=True))]
    if fn == "cx_kemeny":
        g = multi(fx)
        return [mark(nx.kemeny_constant(g, weight=weight))]
    if fn == "cx_second_order":
        g = simple(fx, weighted=bool(weight))
        r = nx.second_order_centrality(g, weight="weight")
        return [mark(r[v]) for v in range(n)]
    if fn == "cx_communicability":
        r = nx.communicability(simple(fx))
        return [tuple(mark(r[i][j]) for j in range(n)) for i in range(n)]
    if fn == "cx_communicability_betweenness":
        r = nx.communicability_betweenness_centrality(simple(fx))
        return [mark(r[v]) for v in range(n)]
    if fn == "cx_dispersion":
        g = multi(fx)
        r = nx.dispersion(g, normalized=e.get("normalized", True),
                          alpha=float(e.get("alpha", 1.0)), b=float(e.get("b", 0.0)),
                          c=float(e.get("c", 0.0)))
        pairs = set()
        for a, b in zip(fx.src, fx.tgt):
            pairs.add((a, b))
            if not fx.directed:
                pairs.add((b, a))
        return [mark(r[u][v]) for u, v in sorted(pairs)]
    if fn == "cx_voterank":
        g = multi(fx)
        return [float(v) for v in nx.voterank(g)]
    if fn == "cx_laplacian":
        g = multi(fx)
        r = guarded(lambda: nx.laplacian_centrality(
            g, normalized=e.get("normalized", True), weight=weight))
        return [UNDEFINED] * n if r is None else [mark(r[v]) for v in range(n)]
    if fn == "cx_trophic":
        what = e["what"]
        g = multi(fx, weighted=bool(weight)) if what == "levels" else simple(
            fx, weighted=bool(weight))
        if what == "levels":
            r = guarded(lambda: nx.trophic_levels(g, weight="weight"))
            return [UNDEFINED] * n if r is None else [mark(r[v]) for v in range(n)]
        if what == "differences":
            r = guarded(lambda: nx.trophic_differences(g, weight="weight"))
            m = len(fx.src)
            if r is None:
                return [UNDEFINED] * m
            return [mark(r[(a, b)]) for a, b in zip(fx.src, fx.tgt)]
        r = guarded(lambda: nx.trophic_incoherence_parameter(
            g, weight="weight", cannibalism=e.get("cannibalism", False)))
        return [UNDEFINED if r is None else mark(r)]
    if fn == "cx_reaching":
        g = simple(fx) if weight else multi(fx)
        normalized = e.get("normalized", True)
        if e.get("what", "local") == "global":
            r = guarded(lambda: nx.global_reaching_centrality(g, weight=weight,
                                                              normalized=normalized))
            return [UNDEFINED if r is None else mark(r)]
        out = []
        for v in range(n):
            r = guarded(lambda: nx.local_reaching_centrality(g, v, weight=weight,
                                                             normalized=normalized))
            out.append(UNDEFINED if r is None else mark(r))
        return out
    raise ValueError(f"no NetworkX oracle for {fn}")


def main() -> int:
    manifest = load_area("centrality_extra")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
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
    write_version(ORACLE, f"networkx {nx.__version__}")
    print(f"networkx (centrality_extra): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
