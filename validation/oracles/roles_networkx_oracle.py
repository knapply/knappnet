#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for reach,
proximity prestige, edge betweenness, and information centrality
(validation/manifest.d/roles.yaml).

- Reach counts the nodes within k steps found by
  `single_source_shortest_path_length(cutoff = k)` (NetworkX has no reach
  function); weighted reach uses `single_source_dijkstra_path_length`.
- Proximity prestige is `closeness_centrality()` on the directed graph,
  which NetworkX computes over paths *to* each node with the Wasserman-Faust
  scaling; on undirected graphs it is the same index.
- Edge betweenness is `edge_betweenness_centrality()` on the multigraph,
  which keys parallel edges separately and so agrees with the engine's
  one-value-per-edge convention.
- Information centrality is `information_centrality()` on the simple
  undirected graph with tie weights summed over parallel edges, multiplied
  by the number of nodes (NetworkX returns Stephenson and Zelen's I_i / n);
  NetworkX refuses disconnected graphs, which are listed as undefined in the
  manifest.

Undefined values are written as empty lines.

    python3 validation/oracles/roles_networkx_oracle.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import networkx as nx
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (GOLDEN, MANIFEST_DIR, applies, fixture_names,  # noqa: E402
                    load_fixture, write_values)

ORACLE = "networkx"


def multigraph(fx, weighted: bool, strength: bool = False):
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        w = fx.weights[i] if weighted else 1.0
        g.add_edge(a, b, weight=(1.0 / w if strength else w))
    return g


def oriented(g, mode: str):
    """The graph whose out-paths are the paths `mode` asks for."""
    if not g.is_directed():
        return g
    if mode == "out":
        return g
    if mode == "in":
        return g.reverse(copy=True)
    return g.to_undirected(as_view=False)


def reach(fx, e):
    k = math.inf if e["k"] == "inf" else float(e["k"])
    weighted = e.get("weighted", False)
    g = oriented(multigraph(fx, weighted, e.get("weights_as") == "strength"), e.get("mode", "out"))
    out = []
    for v in range(fx.n):
        if weighted:
            lengths = nx.single_source_dijkstra_path_length(g, v, weight="weight")
            # Floating-point sums: admit lengths within the engine's tolerance.
            count = sum(1 for u, d in lengths.items() if u != v and d <= k * (1 + 1e-10))
        else:
            cutoff = None if k == math.inf else int(k)
            lengths = nx.single_source_shortest_path_length(g, v, cutoff=cutoff)
            count = len(lengths) - 1
        if e.get("normalized"):
            out.append(math.nan if fx.n < 2 else count / (fx.n - 1))
        else:
            out.append(float(count))
    return out


def proximity_prestige(fx, e):
    g = multigraph(fx, e.get("weighted", False), e.get("weights_as") == "strength")
    c = nx.closeness_centrality(g, distance="weight" if e.get("weighted") else None)
    return [c[v] for v in range(fx.n)]


def edge_betweenness(fx, e):
    g = multigraph(fx, e.get("weighted", False), e.get("weights_as") == "strength")
    eb = nx.edge_betweenness_centrality(g, normalized=e.get("normalized", False),
                                        weight="weight" if e.get("weighted") else None)
    out = []
    seen = {}
    for a, b in zip(fx.src, fx.tgt):
        # Multigraph keys number parallel edges in insertion order, and an
        # undirected key may be stored under either endpoint order.
        pair = (a, b) if fx.directed else (min(a, b), max(a, b))
        key = seen.get(pair, 0)
        seen[pair] = key + 1
        if (a, b, key) in eb:
            out.append(eb[(a, b, key)])
        else:
            out.append(eb[(b, a, key)])
    return out


def information_centrality(fx, e):
    weighted = e.get("weighted", False)
    g = nx.Graph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        w = fx.weights[i] if weighted else 1.0
        if fx.directed and weighted:
            # The engine sums each direction and keeps the larger total.
            g.add_edge(a, b)
            d = g[a][b]
            d.setdefault("fwd", 0.0)
            d.setdefault("bwd", 0.0)
            if (a, b) == (min(a, b), max(a, b)):
                d["fwd"] += w
            else:
                d["bwd"] += w
        elif g.has_edge(a, b):
            if weighted:
                g[a][b]["weight"] += w
        else:
            g.add_edge(a, b, weight=w)
    if fx.directed and weighted:
        for a, b, d in g.edges(data=True):
            d["weight"] = max(d.get("fwd", 0.0), d.get("bwd", 0.0))
    if fx.n < 2 or not nx.is_connected(g):
        return None
    c = nx.information_centrality(g, weight="weight")
    return [c[v] * fx.n for v in range(fx.n)]


def compute(fx, e):
    fn = e["fn"]
    if fn == "reach":
        return reach(fx, e)
    if fn == "closeness":
        return proximity_prestige(fx, e)
    if fn == "edge_betweenness":
        return edge_betweenness(fx, e)
    if fn == "information_centrality":
        return information_centrality(fx, e)
    raise ValueError(f"no NetworkX oracle for {fn!r}")


def main() -> None:
    spec = yaml.safe_load((MANIFEST_DIR / "roles.yaml").read_text())
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(fx, v["engine"])
            if values is None:
                continue
            write_values(ORACLE, name, variant,
                         [() if isinstance(x, float) and math.isnan(x) else x for x in values])
            count += 1
    print(f"networkx {nx.__version__} (roles): wrote {count} golden files")
    (GOLDEN / ORACLE).mkdir(parents=True, exist_ok=True)
    (GOLDEN / ORACLE / "VERSION").write_text(f"networkx {nx.__version__}\n", newline="")


if __name__ == "__main__":
    main()
