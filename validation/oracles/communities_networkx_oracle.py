#!/usr/bin/env python3
"""Golden values for community detection from NetworkX (Hagberg, Schult,
and Swart 2008).

- Semi-synchronous label propagation: the partition of
  `nx.community.label_propagation_communities` (deterministic) on the
  simple undirected graph.
- Partition quality: the mean modularity, over seeds 1 to `seeds`, of the
  partitions of `asyn_lpa_communities` (label propagation; a node keeps
  a dominant label, and the method stops after a sweep without change,
  which is the engine's retention variant),
  `fast_label_propagation_communities`, `asyn_fluidc` (fluid
  communities), and `leiden_communities` (with the modularity metric),
  each evaluated with `nx.community.modularity` on the multigraph.
- Partition comparison: NetworkX has none of these measures, so none are
  written.

    python3 validation/oracles/communities_networkx_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_fixture, load_manifest,  # noqa: E402
                    write_values, write_version)

ORACLE = "networkx"
MEASURES = {"label_propagation", "fluid_communities", "leiden"}


def multigraph(fx):
    g = nx.MultiGraph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, weight=fx.weights[i] if fx.weighted else 1.0)
    return g


def simple(fx, weighted=False):
    """Simple undirected graph; with `weighted`, parallel edges add up."""
    g = nx.Graph()
    g.add_nodes_from(range(fx.n))
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a == b:
            continue
        w = fx.weights[i] if (weighted and fx.weighted) else 1.0
        if g.has_edge(a, b):
            if weighted:
                g[a][b]["weight"] += w
        else:
            g.add_edge(a, b, weight=w)
    return g


def labels_of(groups, n):
    lab = [0] * n
    for c, members in enumerate(groups):
        for v in members:
            lab[v] = c
    return lab


def mean_modularity(fx, run, seeds, weighted):
    g = multigraph(fx)
    weight = "weight" if weighted else None
    if g.number_of_edges() == 0:
        return [-1.0]
    qs = []
    for s in range(1, seeds + 1):
        groups = [set(c) for c in run(s)]
        qs.append(nx.community.modularity(g, groups, weight=weight))
    return [sum(qs) / len(qs)]


def compute(spec: dict, fx):
    e = spec["engine"]
    fn = e["fn"]
    weighted = e.get("weighted", False)
    weight = "weight" if weighted else None
    if fn == "label_propagation":
        variant = e.get("variant", "dominance")
        output = e.get("output", "membership")
        if output == "membership" and variant == "semi_synchronous":
            groups = nx.community.label_propagation_communities(simple(fx))
            return labels_of(groups, fx.n)
        if output == "quality":
            g = simple(fx, weighted)
            if variant == "retention":
                def run(s):
                    return nx.community.asyn_lpa_communities(g, weight=weight, seed=s)
            elif variant == "fast":
                def run(s):
                    return nx.community.fast_label_propagation_communities(g, weight=weight,
                                                                           seed=s)
            else:
                return None
            return mean_modularity(fx, run, e["seeds"], weighted)
        return None
    if fn == "fluid_communities" and e.get("output") == "quality":
        g = simple(fx)
        if fx.n < 2 or not nx.is_connected(g) or e["k"] > fx.n:
            return [-1.0]

        def run(s):
            return nx.community.asyn_fluidc(g, e["k"], seed=s)
        return mean_modularity(fx, run, e["seeds"], False)
    if fn == "leiden" and e.get("output") == "mean_quality" \
            and e.get("objective", "modularity") == "modularity":
        g = simple(fx, weighted)

        def run(s):
            return nx.community.leiden_communities(g, weight=weight, seed=s,
                                                   metric="modularity",
                                                   resolution=e.get("resolution", 1.0))
        return mean_modularity(fx, run, e["seeds"], weighted)
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
            if name in spec["oracles"][ORACLE].get("undefined_on", []):
                continue
            values = compute(spec, fx)
            if values is None:
                continue
            write_values(ORACLE, name, variant, values)
            count += 1
    write_version(ORACLE, f"networkx {nx.__version__}")
    print(f"networkx (communities): wrote {count} golden files")


if __name__ == "__main__":
    main()
