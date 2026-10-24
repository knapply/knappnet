#!/usr/bin/env python3
"""Golden values from NetworkX (Hagberg, Schult, and Swart 2008) for VF2++
isomorphism, induced subgraph isomorphism, and monomorphism, the quick
isomorphism checks, and tree isomorphism (see ../bindings/isomorphism.py
for the derived inputs and encodings).

    python3 validation/oracles/isomorphism_networkx_oracle.py
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import networkx as nx
from networkx.algorithms import isomorphism as iso

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bindings"))
from common import (applies, fixture_names, load_area, load_fixture, write_values,  # noqa: E402
                    write_version)

ORACLE = "networkx"
CAP = 5000
PATTERNS = {
    "path3": (3, [(0, 1), (1, 2)]),
    "triangle": (3, [(0, 1), (1, 2), (2, 0)]),
    "star3": (4, [(0, 1), (0, 2), (0, 3)]),
    "cycle4": (4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    "feedforward": (3, [(0, 1), (1, 2), (0, 2)]),
}


def build(n, pairs, directed):
    multi = len({(a, b) if directed else (min(a, b), max(a, b)) for a, b in pairs}) < len(pairs)
    cls = {(True, True): nx.MultiDiGraph, (True, False): nx.DiGraph,
           (False, True): nx.MultiGraph, (False, False): nx.Graph}[(directed, multi)]
    g = cls()
    for v in range(n):
        g.add_node(v, label=v % 2)
    g.add_edges_from(pairs)
    return g


def rev(fx):
    n = fx.n
    return n, [(n - 1 - a, n - 1 - b) for a, b in zip(fx.src, fx.tgt)][::-1]


def perturbed(fx):
    return fx.n, list(zip(fx.src, fx.tgt))[1:] + [(0, fx.n - 1)]


def valid_edges(m, a, b):
    """m maps nodes of tree a onto nodes of tree b: every edge maps onto an
    edge and the map is a bijection."""
    if sorted(m.values()) != sorted(b.nodes):
        return False
    return all(b.has_edge(m[u], m[v]) for u, v in a.edges)


def tree(fx, depth_first):
    nbrs = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            nbrs[a].add(b)
            nbrs[b].add(a)
    name = {0: 0}
    edges = []
    if depth_first:
        stack, seen = [0], {0}
        while stack:
            v = stack[-1]
            nxt = next((w for w in sorted(nbrs[v]) if w not in seen), None)
            if nxt is None:
                stack.pop()
                continue
            seen.add(nxt)
            name[nxt] = len(name)
            edges.append((name[v], name[nxt]))
            stack.append(nxt)
    else:
        queue = [0]
        for v in queue:
            for w in sorted(nbrs[v]):
                if w not in name:
                    name[w] = len(name)
                    edges.append((name[v], name[w]))
                    queue.append(w)
    g = nx.Graph()
    g.add_nodes_from(range(len(name)))
    g.add_edges_from(edges)
    return g


def compute(fx, e):
    fn = e["fn"]
    if fn == "vf2pp":
        g = build(fx.n, list(zip(fx.src, fx.tgt)), fx.directed)
        if e["problem"] == "iso":
            h = build(*rev(fx), fx.directed)
            gen = nx.vf2pp_all_isomorphisms
        else:
            h = build(*PATTERNS[e["pattern"]], fx.directed)
            gen = (nx.vf2pp_all_subgraph_isomorphisms if e["problem"] == "induced"
                   else nx.vf2pp_all_monomorphisms)
        it = gen(g, h, node_label="label" if e.get("labels") else None)
        count = sum(1 for _ in itertools.islice(it, CAP))
        return [(count, 1 if count else 0)]
    if fn == "could_be_iso":
        g = build(fx.n, list(zip(fx.src, fx.tgt)), fx.directed)
        other = rev(fx) if e["against"] == "rev" else perturbed(fx)
        h = build(*other, fx.directed)
        return [int(nx.could_be_isomorphic(g, h, properties=e["properties"]))]
    if fn == "tree_iso":
        t1, t3 = tree(fx, False), tree(fx, True)
        k = len(t1)
        t2 = nx.relabel_nodes(t1, {v: k - 1 - v for v in t1})

        def res(pairs, a, b):
            if not pairs:
                return (0, 0)
            return (1, int(valid_edges(dict(pairs), a, b)))
        free_rev = res(iso.tree_isomorphism(t1, t2), t1, t2)
        free_dfs = res(iso.tree_isomorphism(t1, t3), t1, t3) if len(t3) == k else (0, 0)
        rooted = res(iso.rooted_tree_isomorphism(t1, 0, t2, k - 1), t1, t2)
        wrong = res(iso.rooted_tree_isomorphism(t1, 0, t2, k - 2), t1, t2) if k > 2 else (0, 0)
        return [free_rev, free_dfs, rooted, wrong]
    return None


def undefined(fx, spec):
    o = spec["oracles"][ORACLE]
    tags = set(fx.tags) | {"directed" if fx.directed else "undirected"}
    return fx.name in o.get("undefined_on", []) or tags & set(o.get("undefined_when", []))


def main() -> None:
    spec = load_area("isomorphism")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx) or undefined(fx, v):
                continue
            values = compute(fx, v["engine"])
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    print(f"networkx {nx.__version__} (isomorphism): wrote {count} golden files")
    write_version(ORACLE, f"networkx {nx.__version__}")


if __name__ == "__main__":
    main()
