#!/usr/bin/env python3
"""Golden values for cohesive subgroups from NetworkX (Hagberg, Schult,
and Swart 2008).

- Coreness: `nx.core_number` on the simple graph (Graph or DiGraph) with
  self-loops removed, since it rejects multigraphs and loops. The manifest
  records the resulting difference on fixtures with loops or parallel
  edges.
- Maximal cliques: `nx.find_cliques` on the undirected simple projection.
- Modularity: `nx.community.modularity` on the multigraph, so parallel
  edges and loops count as in igraph and the engine.
- Louvain: the modularity NetworkX assigns to the engine's partition, and
  (for partition quality) the mean modularity of NetworkX's own
  `louvain_communities` partitions over seeds 1 to 20.

    python3 validation/oracles/subgroups_networkx_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import GOLDEN, applies, fixture_names, load_fixture, write_values  # noqa: E402
from common import load_manifest  # noqa: E402

ORACLE = "networkx"
MEASURES = {"coreness", "max_cliques", "clique_number", "modularity", "louvain_modularity",
            "louvain_quality"}


def multigraph(fx):
    g = nx.MultiDiGraph() if fx.directed else nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, weight=fx.weights[i] if fx.weighted else 1.0)
    return g


def undirected_simple(fx):
    g = nx.Graph()
    g.add_nodes_from(range(fx.n))
    g.add_edges_from((a, b) for a, b in zip(fx.src, fx.tgt) if a != b)
    return g


def compute(spec: dict, fx):
    e = spec["engine"]
    fn = e["fn"]
    if fn == "coreness":
        g = nx.DiGraph() if fx.directed else nx.Graph()
        g.add_nodes_from(range(fx.n))
        g.add_edges_from((a, b) for a, b in zip(fx.src, fx.tgt) if a != b)
        core = nx.core_number(g)
        return [core[v] for v in range(fx.n)]
    if fn in ("max_cliques", "clique_number"):
        cliques = [tuple(sorted(c)) for c in nx.find_cliques(undirected_simple(fx))]
        if fn == "clique_number":
            return [max((len(c) for c in cliques), default=0)]
        return [c for c in cliques if len(c) >= e.get("min_size", 0)]
    weight = "weight" if e.get("weighted", False) else None
    resolution = e.get("resolution", 1.0)
    if fn == "modularity":
        blocks = [3 * v // fx.n for v in range(fx.n)]
        groups = [{v for v in range(fx.n) if blocks[v] == c} for c in sorted(set(blocks))]
        return [nx.community.modularity(multigraph(fx), groups, weight=weight,
                                        resolution=resolution)]
    if fn == "louvain_modularity":
        from engine import Graph, _load_bindings
        _load_bindings()
        binding = sys.modules["bindings_subgroups"]
        g = Graph(fx.n, fx.src, fx.tgt, fx.weights, fx.directed)
        membership, _ = binding.louvain(g, resolution, e.get("weighted", False),
                                        e.get("seed", 1))
        groups = [{v for v in range(fx.n) if membership[v] == c} for c in set(membership)]
        return [nx.community.modularity(multigraph(fx), groups, weight=weight,
                                        resolution=resolution)]
    if fn == "louvain_quality":
        g = multigraph(fx)
        qs = []
        for seed in range(1, e["seeds"] + 1):
            groups = nx.community.louvain_communities(g, weight=weight, resolution=resolution,
                                                      seed=seed)
            qs.append(nx.community.modularity(g, groups, weight=weight, resolution=resolution))
        return [sum(qs) / len(qs)]
    return None


def main() -> None:
    manifest = load_manifest()
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, spec in manifest.items():
            if spec["engine"]["fn"] not in MEASURES:
                continue
            if ORACLE not in spec["oracles"] or not applies(spec, fx):
                continue
            values = compute(spec, fx)
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    (GOLDEN / ORACLE).mkdir(parents=True, exist_ok=True)
    (GOLDEN / ORACLE / "VERSION").write_text(f"networkx {nx.__version__}\n", newline="")
    print(f"networkx {nx.__version__} (subgroups): wrote {count} golden files")


if __name__ == "__main__":
    main()
