#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for
spanning trees, branchings and arborescences, topological orders,
transitive closure and reduction, longest paths, dominators, and lowest
common ancestors.

Graphs are built from the fixture file in order (nodes 0..n - 1, then the
edges), as a multigraph where NetworkX accepts one (spanning trees,
branchings, longest paths) and as a simple graph elsewhere. The derived
networks and encodings are those of validation/bindings/trees.py.

    python3 validation/oracles/trees_networkx_oracle.py
"""

from __future__ import annotations

import itertools
import sys
from collections import deque
from pathlib import Path

import networkx as nx

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (applies, fixture_names, load_area, load_fixture,  # noqa: E402
                    write_values, write_version)

ORACLE = "networkx"
MISSING = -1.0


def weights_of(fx, weighted):
    return list(fx.weights) if weighted and fx.weighted else [1.0] * len(fx.src)


def multigraph(fx, weighted):
    g = (nx.MultiDiGraph if fx.directed else nx.MultiGraph)()
    g.add_nodes_from(range(fx.n))
    w = weights_of(fx, weighted)
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        g.add_edge(a, b, key=e, weight=w[e])
    return g


def dag(fx, weighted=False, multi=False):
    g = nx.MultiDiGraph() if multi else nx.DiGraph()
    g.add_nodes_from(range(fx.n))
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        if a < b:
            if weighted and fx.weighted:
                g.add_edge(a, b, weight=fx.weights[e])
            else:
                g.add_edge(a, b)
    return g


def simple_digraph(fx, loops=True):
    g = nx.DiGraph()
    g.add_nodes_from(range(fx.n))
    g.add_edges_from((a, b) for a, b in zip(fx.src, fx.tgt) if loops or a != b)
    return g


def bfs_tree(fx):
    adj = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            adj[a].add(b)
            adj[b].add(a)
    parent = {0: None} if fx.n else {}
    q = deque([0] if fx.n else [])
    arcs = []
    while q:
        v = q.popleft()
        for w in sorted(adj[v]):
            if w not in parent:
                parent[w] = v
                arcs.append((v, w))
                q.append(w)
    return arcs, sorted(parent)


def edge_tuple(fx, u, v, w):
    if not fx.directed:
        u, v = min(u, v), max(u, v)
    return (u, v, w)


def all_pairs(nodes):
    return [(u, v) for i, u in enumerate(nodes) for v in nodes[i:]]


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    n = fx.n
    weighted = e.get("weighted", False)
    if fn == "spanning_tree":
        g = multigraph(fx, weighted)
        f = nx.maximum_spanning_edges if e.get("maximum") else nx.minimum_spanning_edges
        return [edge_tuple(fx, u, v, d["weight"])
                for u, v, k, d in f(g, algorithm="kruskal", keys=True, data=True)]
    if fn == "spanning_weight":
        g = multigraph(fx, weighted)
        f = nx.maximum_spanning_tree if e.get("maximum") else nx.minimum_spanning_tree
        return [f(g).size(weight="weight")]
    if fn == "branching":
        g = multigraph(fx, weighted)
        f = {"max_branching": nx.maximum_branching, "min_branching": nx.minimum_branching,
             "minimal": nx.algorithms.tree.branchings.minimal_branching,
             "max_arborescence": nx.maximum_spanning_arborescence,
             "min_arborescence": nx.minimum_spanning_arborescence}[e["kind"]]
        try:
            b = f(g, attr="weight")
        except (nx.NetworkXException, nx.NetworkXPointlessConcept):
            return [MISSING]
        return [b.size(weight="weight")]
    if fn == "topological_sort":
        return list(nx.topological_sort(dag(fx)))
    if fn == "generations":
        label = [0] * n
        for i, gen in enumerate(nx.topological_generations(dag(fx))):
            for v in gen:
                label[v] = i
        return label
    if fn == "lexicographic":
        return list(nx.lexicographical_topological_sort(dag(fx)))
    if fn == "all_topological":
        return [tuple(o) for o in itertools.islice(nx.all_topological_sorts(dag(fx)), 50)]
    if fn == "is_dag":
        g = multigraph(fx, False)
        return [1.0 if nx.is_directed_acyclic_graph(g) else 0.0]
    if fn == "closure":
        refl = {"cycles": False, "never": None, "always": True}[e["reflexive"]]
        t = nx.transitive_closure(simple_digraph(fx), reflexive=refl)
        return sorted(t.edges())
    if fn == "reduction":
        return sorted(nx.transitive_reduction(dag(fx)).edges())
    if fn == "longest_path":
        return nx.dag_longest_path(dag(fx, weighted, multi=True))
    if fn == "longest_length":
        return [nx.dag_longest_path_length(dag(fx, weighted, multi=True))]
    if fn == "idom":
        if n == 0:
            return []
        idom = nx.immediate_dominators(simple_digraph(fx), 0)
        return [0 if v == 0 else idom.get(v, MISSING) for v in range(n)]
    if fn == "frontiers":
        df = nx.dominance_frontiers(simple_digraph(fx), 0)
        return [tuple(sorted(df[v])) if v in df else (MISSING,) for v in range(n)]
    if fn == "dag_lca":
        found = dict(nx.all_pairs_lowest_common_ancestor(dag(fx), pairs=all_pairs(list(range(n)))))
        return [found.get(p, MISSING) for p in all_pairs(list(range(n)))]
    if fn == "dag_lca_check":
        found = dict(nx.all_pairs_lowest_common_ancestor(dag(fx), pairs=all_pairs(list(range(n)))))
        return [len(found), 1.0]
    if fn == "tree_lca":
        arcs, nodes = bfs_tree(fx)
        t = nx.DiGraph()
        t.add_nodes_from(range(n))
        t.add_edges_from(arcs)
        pairs = all_pairs(nodes)
        found = dict(nx.tree_all_pairs_lowest_common_ancestor(t, root=0, pairs=pairs))
        return [found[p] for p in pairs]
    raise ValueError(fn)


def main() -> int:
    manifest = load_area("trees")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        tags = set(fx.tags) | {"directed" if fx.directed else "undirected"}
        for variant, spec in manifest.items():
            ospec = spec["oracles"].get(ORACLE)
            if (ospec is None or not applies(spec, fx) or name in ospec.get("undefined_on", [])
                    or tags & set(ospec.get("undefined_when", []))):
                continue
            write_values(ORACLE, name, variant, compute(spec, fx))
            count += 1
    write_version(ORACLE, f"networkx {nx.__version__}")
    print(f"networkx (trees): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
