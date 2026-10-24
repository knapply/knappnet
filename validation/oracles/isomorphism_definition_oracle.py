#!/usr/bin/env python3
"""Definition oracle for isomorphism, induced subgraph isomorphism,
monomorphism, the quick isomorphism checks, and tree isomorphism. Shares
no code with the engine:

- mappings: injective node maps built one node at a time, each partial map
  checked against the definition (equal edge multiplicities for every pair
  of mapped nodes, or no more in the pattern for monomorphisms; equal
  labels), counted up to 5000; isomorphisms by every permutation on
  networks of at most 8 nodes;
- quick checks: per-node degree, triangle, and maximal-clique counts by
  brute force over node triples and subsets;
- trees: canonical parenthesis strings of every rooting.

    python3 validation/oracles/isomorphism_definition_oracle.py
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
CAP = 5000
PATTERNS = {
    "path3": (3, [(0, 1), (1, 2)]),
    "triangle": (3, [(0, 1), (1, 2), (2, 0)]),
    "star3": (4, [(0, 1), (0, 2), (0, 3)]),
    "cycle4": (4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    "feedforward": (3, [(0, 1), (1, 2), (0, 2)]),
}


def counts(pairs, directed):
    c = {}
    for a, b in pairs:
        k = (a, b) if directed else frozenset((a, b))
        c[k] = c.get(k, 0) + 1
    return c


def mult(c, a, b, directed):
    return c.get((a, b) if directed else frozenset((a, b)), 0)


def count_maps(sg, fg, problem, directed, labels):
    (ns, ps), (nf, pf) = sg, fg
    cs, cf = counts(ps, directed), counts(pf, directed)
    mono = problem == "mono"
    total = 0
    m = []

    def ok_pair(u, w):
        s = mult(cs, u, w, directed)
        f = mult(cf, m[u], m[w], directed)
        return s <= f if mono else s == f

    def rec():
        nonlocal total
        if total >= CAP:
            return
        u = len(m)
        if u == ns:
            total += 1
            return
        for x in range(nf):
            if x in m or (labels and x % 2 != u % 2):
                continue
            m.append(x)
            if all(ok_pair(u, w) and (not directed or ok_pair(w, u)) for w in range(u + 1)):
                rec()
            m.pop()
    rec()
    return min(total, CAP)


def count_isos(fx, other, labels):
    n = fx.n
    pf = list(zip(fx.src, fx.tgt))
    cs, cf = counts(other[1], fx.directed), counts(pf, fx.directed)
    total = 0
    for perm in itertools.permutations(range(n)):
        if labels and any(u % 2 != perm[u] % 2 for u in range(n)):
            continue
        if all(mult(cs, u, w, fx.directed) == mult(cf, perm[u], perm[w], fx.directed)
               for u in range(n) for w in range(n)):
            total += 1
            if total >= CAP:
                break
    return total


def props(n, pairs, directed, which):
    deg = [0] * n
    for a, b in pairs:
        deg[a] += 1
        deg[b] += 1
    adj = [set() for _ in range(n)]
    for a, b in pairs:
        if a != b:
            adj[a].add(b)
            adj[b].add(a)
    rows = [[] for _ in range(n)]
    tables = []
    for p in which:
        if p == "d":
            for v in range(n):
                rows[v].append(deg[v])
        elif p == "t":
            tri = [0] * n
            for a, b, c in itertools.combinations(range(n), 3):
                if b in adj[a] and c in adj[a] and c in adj[b]:
                    tri[a] += 1
                    tri[b] += 1
                    tri[c] += 1
            for v in range(n):
                rows[v].append(tri[v])
        elif p == "c":
            cl = [0] * n
            for c in maximal_cliques(adj):
                for v in c:
                    cl[v] += 1
            for v in range(n):
                rows[v].append(cl[v])
        tables.append(sorted(tuple(r) for r in rows))
    return tables


def maximal_cliques(adj):
    """Maximal cliques by growing every clique (exhaustive; small graphs)."""
    n = len(adj)
    found = set()

    def grow(c, cand):
        ext = [v for v in cand if all(v in adj[u] for u in c)]
        if not ext:
            found.add(frozenset(c))
            return
        for v in ext:
            grow(c | {v}, [w for w in ext if w > v])
    for v in range(n):
        grow({v}, [w for w in range(n) if w != v])
    return [c for c in found if not any(c < d for d in found)]


def canon(adj, root):
    def enc(v, parent):
        return "(" + "".join(sorted(enc(w, v) for w in adj[v] if w != parent)) + ")"
    return enc(root, -1)


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
    adj = [set() for _ in range(len(name))]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def compute(fx, e):
    fn = e["fn"]
    pf = list(zip(fx.src, fx.tgt))
    if fn == "vf2pp":
        if e["problem"] == "iso":
            n = fx.n
            other = (n, [(n - 1 - a, n - 1 - b) for a, b in pf][::-1])
            c = count_isos(fx, other, e.get("labels"))
        else:
            c = count_maps(PATTERNS[e["pattern"]], (fx.n, pf), e["problem"], fx.directed,
                           e.get("labels"))
        return [(c, 1 if c else 0)]
    if fn == "could_be_iso":
        n = fx.n
        other = ([(n - 1 - a, n - 1 - b) for a, b in pf][::-1] if e["against"] == "rev"
                 else pf[1:] + [(0, n - 1)])
        a = props(n, pf, fx.directed, e["properties"])
        b = props(n, other, fx.directed, e["properties"])
        return [int(all(x == y for x, y in zip(a, b)))]
    if fn == "tree_iso":
        t1, t3 = tree(fx, False), tree(fx, True)
        k = len(t1)
        t2 = [set() for _ in range(k)]
        for v in range(k):
            for w in t1[v]:
                t2[k - 1 - v].add(k - 1 - w)

        def free(a, b):
            if len(a) != len(b):
                return 0
            ca = canon(a, 0)
            return int(any(canon(b, r) == ca for r in range(len(b))))

        def rooted(r2):
            return int(canon(t1, 0) == canon(t2, r2))
        fr = free(t1, t2)
        fd = free(t1, t3)
        ro = rooted(k - 1)
        wr = rooted(k - 2) if k > 2 else 0
        return [(fr, fr), (fd, fd), (ro, ro), (wr, wr)]
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
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (isomorphism): wrote {count} golden files")


if __name__ == "__main__":
    main()
