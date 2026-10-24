#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for
maximum flows, minimum cuts, minimum-cost flows, and Gomory-Hu trees.

NetworkX's flow functions refuse multigraphs, so the flow graph merges
parallel edges into one edge whose capacity is their sum and drops
self-loops (which carry no flow); minimum-cost flows use the multigraph,
which NetworkX's network simplex accepts. The flow problem of a fixture is
the one documented in validation/bindings/flow.py. Results that do not
exist are written as -1.

    python3 validation/oracles/flow_networkx_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture,  # noqa: E402
                    write_values, write_version)

ORACLE = "networkx"
MISSING = -1.0


def problem(fx, multi=False):
    if multi:
        if fx.n < 4:
            return None
        s, t = [0, 1], [fx.n - 2, fx.n - 1]
    else:
        if fx.n < 2:
            return None
        s, t = [0], [fx.n - 1]
    m = len(fx.src)
    cap = list(fx.weights) if fx.weighted else [1.0] * m
    cost = [float((7 * e) % 11 - 3) for e in range(m)]
    return s, t, cap, cost


def flow_graph(fx, cap):
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        if g.has_edge(a, b):
            g[a][b]["capacity"] += cap[e]
        else:
            g.add_edge(a, b, capacity=cap[e])
    return g


def with_terminals(g, s, t):
    """(graph, source, sink): super-terminals joined by arcs without a
    capacity attribute (infinite in NetworkX) when there are several."""
    if len(s) == 1:
        return g, s[0], t[0]
    h = nx.DiGraph(g)
    for v in s:
        h.add_edge("S", v)
    for v in t:
        h.add_edge(v, "T")
    return h, "S", "T"


def simple_graph(fx):
    g = nx.DiGraph() if fx.directed else nx.Graph()
    g.add_nodes_from(range(fx.n))
    g.add_edges_from((a, b) for a, b in zip(fx.src, fx.tgt) if a != b)
    return g


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    m = len(fx.src)
    if fn in ("gomory_hu", "gomory_hu_values"):
        cap = list(fx.weights) if fx.weighted else [1.0] * m
        g = flow_graph(fx, cap)
        if fn == "gomory_hu":
            t = nx.gomory_hu_tree(g, capacity="capacity")
            return [(min(u, v), max(u, v), w) for u, v, w in t.edges(data="weight")]
        return [nx.minimum_cut_value(g, u, v, capacity="capacity")
                for u in range(fx.n) for v in range(u + 1, fx.n)]
    if fn == "max_flow_min_cost":
        p = problem(fx)
        if p is None:
            return [MISSING, MISSING]
        s, t, cap, cost = p
        value = nx.maximum_flow_value(flow_graph(fx, cap), s[0], t[0], capacity="capacity")
        mg = nx.MultiDiGraph()
        mg.add_nodes_from(range(fx.n))
        for k, (a, b) in enumerate(zip(fx.src, fx.tgt)):
            mg.add_edge(a, b, capacity=cap[k], weight=cost[k])
        mg.nodes[s[0]]["demand"] = -value
        mg.nodes[t[0]]["demand"] = value
        return [value, nx.min_cost_flow_cost(mg)]
    if fn == "min_node_cut":
        if fx.n < 2:
            return [MISSING]
        s, t = 0, fx.n - 1
        if not e.get("weighted"):
            g = simple_graph(fx)
            if g.has_edge(s, t):
                return [MISSING]
            return [float(len(nx.connectivity.minimum_st_node_cut(g, s, t)))]
        g = simple_graph(fx)
        if g.has_edge(s, t):
            return [MISSING]
        h = nx.DiGraph()
        for v in range(fx.n):
            h.add_edge(f"{v}A", f"{v}B", capacity=float(v % 3 + 1))
        for a, b in g.edges():
            h.add_edge(f"{a}B", f"{b}A")
            if not fx.directed:
                h.add_edge(f"{b}B", f"{a}A")
        return [nx.minimum_cut_value(h, f"{s}B", f"{t}A", capacity="capacity")]
    p = problem(fx, e.get("multi", False))
    if p is None:
        return [MISSING, MISSING] if fn == "max_flow_check" else [MISSING]
    s, t, cap, _ = p
    g, src, snk = with_terminals(flow_graph(fx, cap), s, t)
    if fn == "max_flow":
        return [nx.maximum_flow_value(g, src, snk, capacity="capacity")]
    if fn == "max_flow_check":
        return [nx.maximum_flow_value(g, src, snk, capacity="capacity"), 1.0]
    if fn == "min_cut":
        if e["side"] == "largest":
            side = nx.minimum_cut(g, src, snk, capacity="capacity")[1][0]
        else:
            r = nx.flow.edmonds_karp(g, src, snk, capacity="capacity")
            side = {src}
            stack = [src]
            while stack:
                u = stack.pop()
                for v, attr in r[u].items():
                    if v not in side and attr["flow"] < attr["capacity"]:
                        side.add(v)
                        stack.append(v)
        return [1.0 if v in side else 0.0 for v in range(fx.n)]
    raise ValueError(fn)


def main() -> int:
    manifest = load_area("flow")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, spec in manifest.items():
            ospec = spec["oracles"].get(ORACLE)
            if ospec is None or not applies(spec, fx) or name in ospec.get("undefined_on", []):
                continue
            write_values(ORACLE, name, variant, compute(spec, fx))
            count += 1
    write_version(ORACLE, f"networkx {nx.__version__}")
    print(f"networkx (flow): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
