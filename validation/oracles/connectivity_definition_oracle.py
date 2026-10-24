#!/usr/bin/env python3
"""Definition oracle for connectivity.

Everything is evaluated from the definitions on the simple graph, by brute
force, without flows: a local node connectivity is the size of the
smallest node set (other than s and t) whose removal leaves no path from s
to t, plus one if s -> t is an arc (Menger's theorem); a local or global
edge connectivity is the least capacity of a cut over every node subset;
the global node connectivity is the size of the smallest node set whose
removal leaves a network that is not strongly connected (n - 1 for a
complete graph); k-components are the maximal node sets whose induced
subgraph has node connectivity at least k, found among all subsets;
k-edge-components group nodes with pairwise edge connectivity at least k
(both ways); local bridges and semiconnectedness come from breadth-first
search. Encodings follow validation/bindings/connectivity.py.

    python3 validation/oracles/connectivity_definition_oracle.py
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
MAX_NODES = 20
PAIRS_MAX_NODES = 12
KCOMP_MAX_NODES = 10
MISSING = -1.0


def adjacency(n, directed, edges):
    out = [set() for _ in range(n)]
    for a, b in edges:
        if a == b:
            continue
        out[a].add(b)
        if not directed:
            out[b].add(a)
    return out


def reach(adj, s, removed=frozenset(), skip_arc=None):
    seen = {s}
    stack = [s]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w in removed or w in seen or (v, w) == skip_arc:
                continue
            seen.add(w)
            stack.append(w)
    return seen


def local_kappa(adj, n, s, t, directed):
    adjacent = t in adj[s]
    skip = (s, t)
    others = [v for v in range(n) if v not in (s, t)]
    for k in range(len(others) + 1):
        for comb in itertools.combinations(others, k):
            removed = frozenset(comb)
            # A path from s that uses the edge s-t uses it as its first arc,
            # so skipping the arc s -> t removes the edge from every path.
            seen = reach(adj, s, removed, skip)
            if t not in seen:
                return k + (1 if adjacent else 0)
    raise AssertionError("unreachable")


def strongly_connected(adj, nodes, removed):
    nodes = [v for v in nodes if v not in removed]
    if len(nodes) <= 1:
        return True
    r = nodes[0]
    fwd = reach(adj, r, removed)
    if any(v not in fwd for v in nodes):
        return False
    return all(r in reach(adj, v, removed) for v in nodes)


def global_kappa(adj, nodes, directed):
    nodes = list(nodes)
    k = len(nodes)
    sub = set(nodes)
    complete = all(w in adj[v] for v in nodes for w in nodes if v != w)
    if complete:
        return k - 1
    for size in range(k - 1):
        for comb in itertools.combinations(nodes, size):
            removed = frozenset(comb) | frozenset(set(range(len(adj))) - sub)
            if not strongly_connected(adj, nodes, removed):
                return size
    return k - 1


def cut_vector(n, directed, adj):
    idx = np.arange(1 << n, dtype=np.int64)
    bits = [((idx >> v) & 1).astype(bool) for v in range(n)]
    total = np.zeros(1 << n)
    for a in range(n):
        for b in adj[a]:
            if directed:
                total += bits[a] & ~bits[b]
            elif a < b:
                total += bits[a] ^ bits[b]
    return idx, bits, total


def lam(bits, total, s, t):
    return float(total[bits[s] & ~bits[t]].min())


def k_components(adj, n):
    kappa = {}
    for size in range(2, n + 1):
        for comb in itertools.combinations(range(n), size):
            kappa[comb] = global_kappa(adj, comb, False)
    out = []
    kmax = max(kappa.values(), default=0)
    for k in range(1, kmax + 1):
        sets = [set(c) for c, kk in kappa.items() if kk >= k and len(c) > k]
        maximal = [s for s in sets if not any(s < o for o in sets)]
        for s in maximal:
            out.append(tuple(float(k * 1000 + v) for v in sorted(s)))
    return out


def compute(spec, fx):
    e = spec["engine"]
    fn = e["fn"]
    n = fx.n
    adj = adjacency(n, fx.directed, zip(fx.src, fx.tgt))
    if fn == "local_bridges":
        und = adjacency(n, False, zip(fx.src, fx.tgt))
        out = []
        for a, b in zip(fx.src, fx.tgt):
            if a == b or und[a] & und[b]:
                out.append(-1.0)
                continue
            # Distance from a to b without the edge a-b.
            dist = {a: 0}
            frontier = [a]
            found = -2.0
            while frontier and found == -2.0:
                nxt = []
                for v in frontier:
                    for w in und[v]:
                        if {v, w} == {a, b} or w in dist:
                            continue
                        dist[w] = dist[v] + 1
                        if w == b:
                            found = float(dist[w])
                        nxt.append(w)
                frontier = nxt
            out.append(found)
        return out
    if fn == "semiconnected":
        if n == 0:
            return None
        r = [reach(adj, v) for v in range(n)]
        ok = all(v in r[u] or u in r[v] for u in range(n) for v in range(n))
        return [1.0 if ok else 0.0]
    if fn == "attracting_flags":
        r = [reach(adj, v) for v in range(n)]
        return [1.0 if all(v in r[w] for w in r[v]) else 0.0 for v in range(n)]
    if fn == "strong_partition":
        r = [reach(adj, v) for v in range(n)]
        label, out = {}, []
        for v in range(n):
            key = frozenset(w for w in r[v] if v in r[w])
            out.append(label.setdefault(key, len(label)))
        return out
    if n > MAX_NODES:
        return None
    if fn in ("local_connectivity", "disjoint_paths_check"):
        if n < 2:
            return [MISSING] if fn == "local_connectivity" else [MISSING, MISSING]
        if e["kind"] == "node":
            v = local_kappa(adj, n, 0, n - 1, fx.directed)
        else:
            _, bits, total = cut_vector(n, fx.directed, adj)
            v = lam(bits, total, 0, n - 1)
        return [v] if fn == "local_connectivity" else [v, 1.0]
    if fn in ("connectivity", "node_cut_check", "edge_cut_check"):
        kind = e.get("kind", "node" if fn == "node_cut_check" else "edge")
        if kind == "node":
            v = global_kappa(adj, range(n), fx.directed)
        else:
            if n <= 1:
                v = 0.0
            else:
                _, _, total = cut_vector(n, fx.directed, adj)
                v = float(total[1:-1].min())
        return [v] if fn == "connectivity" else [v, 1.0]
    if fn == "st_node_cut_check":
        if n < 2:
            return [MISSING, MISSING]
        if n - 1 in adj[0]:
            return [MISSING, MISSING]
        return [local_kappa(adj, n, 0, n - 1, fx.directed), 1.0]
    if fn in ("all_pairs_node", "average_node"):
        if n > PAIRS_MAX_NODES:
            return None
        m = [[0.0] * n for _ in range(n)]
        for u in range(n):
            for v in range(n):
                if u != v:
                    m[u][v] = local_kappa(adj, n, u, v, fx.directed)
        if fn == "all_pairs_node":
            return [x for row in m for x in row]
        pairs = [(u, v) for u in range(n) for v in range(n)
                 if u != v and (fx.directed or u < v)]
        return [sum(m[u][v] for u, v in pairs) / len(pairs) if pairs else 0.0]
    if fn == "k_components":
        if n > KCOMP_MAX_NODES:
            return None
        return k_components(adj, n)
    if fn == "k_edge_components":
        k = e["k"]
        if n == 0:
            return []
        _, bits, total = cut_vector(n, fx.directed, adj) if n > 1 else (None, None, None)
        label = [-1] * n
        nxt = 0
        for u in range(n):
            if label[u] >= 0:
                continue
            label[u] = nxt
            for v in range(u + 1, n):
                if label[v] < 0 and lam(bits, total, u, v) >= k and lam(bits, total, v, u) >= k:
                    label[v] = nxt
            nxt += 1
        return label
    raise ValueError(fn)


def main() -> int:
    manifest = load_area("connectivity")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, spec in manifest.items():
            ospec = spec["oracles"].get(ORACLE)
            if ospec is None or not applies(spec, fx) or name in ospec.get("undefined_on", []):
                continue
            values = compute(spec, fx)
            if values is None:
                continue
            write_values(ORACLE, name, variant, values)
            count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (connectivity): wrote {count} golden files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
