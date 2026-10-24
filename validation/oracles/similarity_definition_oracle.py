#!/usr/bin/env python3
"""Definition oracle for Weisfeiler-Lehman hashes, SimRank, Panther, and
graph edit distance. Shares no code with the engine:

- hashes: the relabelling written out from Shervashidze et al. (2011) and
  NetworkX's documentation (label + sorted neighbour labels, BLAKE2b from
  Python's hashlib, the graph hash over the tuple of sorted label counts);
- SimRank: Jeh and Widom's pairwise formula over in-neighbour pairs
  (math.fsum), iterated with the stopping rule of each variant; the
  fixed point by dense iteration to 1e-15;
- Panther: the exact probability that a random path visits both nodes,
  by dynamic programming over (node, seen source, seen target);
- edit distance: the minimum over every partial injective node mapping.

    python3 validation/oracles/similarity_definition_oracle.py
"""

from __future__ import annotations

import hashlib
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"
C = 0.9


def hex_ints(text):
    return tuple(int(text[i:i + 8], 16) for i in range(0, len(text), 8))


def blake(text, size):
    return hashlib.blake2b(text.encode("ascii"), digest_size=size).hexdigest()


def structure(fx):
    """Successor and predecessor sets (neighbour sets when undirected) and
    the label of every arc."""
    succ = [set() for _ in range(fx.n)]
    pred = [set() for _ in range(fx.n)]
    elabel = {}
    for e, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        succ[a].add(b)
        elabel[(a, b)] = str(e % 2)
        if fx.directed:
            pred[b].add(a)
        else:
            succ[b].add(a)
            elabel[(b, a)] = str(e % 2)
    return succ, pred, elabel


def wl(fx, e):
    succ, pred, elabel = structure(fx)
    size = int(e.get("digest_size", 16))
    iterations = int(e.get("iterations", 3))
    use_e = bool(e.get("edge_labels"))
    if e.get("node_labels"):
        labels = [str(v % 3) for v in range(fx.n)]
    elif use_e:
        labels = [""] * fx.n
    else:
        iterations -= 1
        if fx.directed:
            labels = [f"{len(pred[v])}_{len(succ[v])}" for v in range(fx.n)]
        else:
            labels = [str(len(succ[v]) + (v in succ[v])) for v in range(fx.n)]
    history = [[] for _ in range(fx.n)]
    if e.get("include_initial"):
        for v in range(fx.n):
            history[v].append(blake(labels[v], size))
    if iterations < int(e.get("iterations", 3)):
        for v in range(fx.n):
            history[v].append(blake(labels[v], size))
    counts = []
    for _ in range(iterations):
        new = []
        for v in range(fx.n):
            if fx.directed:
                s = sorted((elabel[(v, u)] if use_e else "s_") + labels[u] for u in succ[v])
                p = sorted((elabel[(u, v)] if use_e else "p_") + labels[u] for u in pred[v])
                text = labels[v] + "".join(s) + "".join(p)
            else:
                text = labels[v] + "".join(sorted(
                    (elabel[(v, u)] if use_e else "") + labels[u] for u in succ[v]))
            new.append(blake(text, size))
        labels = new
        for v in range(fx.n):
            history[v].append(labels[v])
        counts.extend(sorted(Counter(labels).items()))
    if e["fn"] == "wl_hash":
        return [hex_ints(blake(str(tuple(counts)), size))]
    return [hex_ints("".join(history[v])) for v in range(fx.n)]


def in_neighbours(fx, weighted):
    w = [dict() for _ in range(fx.n)]
    for i, (a, b) in enumerate(zip(fx.src, fx.tgt)):
        x = fx.weights[i] if weighted else 1.0
        pairs = [(a, b)] if fx.directed or a == b else [(a, b), (b, a)]
        for s, t in pairs:
            w[t][s] = w[t].get(s, 0.0) + x if weighted else 1.0
    return w


def simrank(fx, e):
    n = fx.n
    weighted = bool(e.get("weighted"))
    tol = float(e.get("tolerance", 1e-4))
    inw = in_neighbours(fx, weighted)
    if tol < 1e-11:
        w = np.zeros((n, n))
        for t in range(n):
            for s, x in inw[t].items():
                w[s, t] = x
        col = w.sum(axis=0)
        col[col == 0] = 1
        w = w / col
        s = np.eye(n)
        for _ in range(100000):
            new = C * (w.T @ s @ w)
            np.fill_diagonal(new, 1.0)
            if np.max(np.abs(new - s), initial=0) < 1e-15:
                s = new
                break
            s = new
        return [float(s[u, v]) for u in range(n) for v in range(n)]
    tot = [sum(inw[v].values()) for v in range(n)]
    s = [[1.0 if u == v else 0.0 for v in range(n)] for u in range(n)]
    relative = e.get("stop", "allclose") == "relative"
    for _ in range(1000):
        new = [[1.0] * n for _ in range(n)]
        for u in range(n):
            for v in range(n):
                if u == v:
                    continue
                if not inw[u] or not inw[v]:
                    new[u][v] = 0.0
                    continue
                total = math.fsum(xa * xb * s[a][b] for a, xa in inw[u].items()
                                  for b, xb in inw[v].items())
                new[u][v] = C * total / (tot[u] * tot[v])
        if relative:
            close = all(abs(new[u][v] - s[u][v]) <= tol * (1 + abs(s[u][v]))
                        for u in range(n) for v in range(n))
        else:
            close = all(abs(s[u][v] - new[u][v]) <= tol + 1e-5 * abs(new[u][v])
                        for u in range(n) for v in range(n))
        s = new
        if close:
            break
    return [s[u][v] for u in range(n) for v in range(n)]


def hub(fx):
    nbrs = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        if a != b:
            nbrs[a].add(b)
            nbrs[b].add(a)
    return max(range(fx.n), key=lambda v: (len(nbrs[v]), -v))


def panther(fx, steps=5):
    nb = [set() for _ in range(fx.n)]
    for a, b in zip(fx.src, fx.tgt):
        nb[a].add(b)
        nb[b].add(a)
    kept = [v for v in range(fx.n) if nb[v]]
    s = hub(fx)
    out = [-1.0] * fx.n
    for t in kept:
        if t == s:
            continue
        # P[(v, seen_s, seen_t)] after each step.
        p = {}
        for v in kept:
            key = (v, v == s, v == t)
            p[key] = p.get(key, 0.0) + 1.0 / len(kept)
        for _ in range(steps):
            q = {}
            for (v, a, b), x in p.items():
                for u in nb[v]:
                    key = (u, a or u == s, b or u == t)
                    q[key] = q.get(key, 0.0) + x / len(nb[v])
            p = q
        out[t] = math.fsum(x for (v, a, b), x in p.items() if a and b)
    return out


def ged(fx, labels, max_nodes=7):
    if fx.n > max_nodes:
        return None
    n1 = fx.n
    e1 = list(zip(fx.src, fx.tgt))
    e2 = e1[1:]
    present = set(e2) | (set() if fx.directed else {(b, a) for a, b in e2})
    if (0, n1 - 1) not in present:
        e2.append((0, n1 - 1))
    e2.append((n1, min(1, n1 - 1)))
    n2 = n1 + 1

    def edge_map(edges):
        m = {}
        for i, (a, b) in enumerate(edges):
            key = (a, b) if fx.directed else (min(a, b), max(a, b))
            m[key] = i % 2
        return m
    E1, E2 = edge_map(e1), edge_map(e2)
    node_lab = labels in ("nodes", "both")
    edge_lab = labels == "both"
    best = [math.inf]

    def key(a, b):
        return (a, b) if fx.directed else (min(a, b), max(a, b))

    def cost(f):
        c = 0
        for a in range(n1):
            if f[a] is None:
                c += 1
            elif node_lab and a % 2 != f[a] % 2:
                c += 1
        c += n2 - sum(1 for x in f if x is not None)
        covered = set()
        for (a, b), lab in E1.items():
            x, y = f[a], f[b]
            if x is not None and y is not None and key(x, y) in E2:
                covered.add(key(x, y))
                if edge_lab and E2[key(x, y)] != lab:
                    c += 1
            else:
                c += 1
        c += len(E2) - len(covered)
        return c

    def rec(a, f, used):
        if a == n1:
            best[0] = min(best[0], cost(f))
            return
        for x in [None] + [x for x in range(n2) if x not in used]:
            f.append(x)
            if x is not None:
                used.add(x)
            rec(a + 1, f, used)
            if x is not None:
                used.discard(x)
            f.pop()
    rec(0, [], set())
    return [float(best[0])]


def undefined(fx, spec):
    o = spec["oracles"][ORACLE]
    tags = set(fx.tags) | {"directed" if fx.directed else "undirected"}
    return fx.name in o.get("undefined_on", []) or tags & set(o.get("undefined_when", []))


def main() -> None:
    spec = load_area("similarity")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx) or undefined(fx, v):
                continue
            e = v["engine"]
            if e["fn"] in ("wl_hash", "wl_subgraph"):
                values = wl(fx, e)
            elif e["fn"] == "simrank":
                values = simrank(fx, e)
            elif e["fn"] == "panther":
                values = panther(fx)
            elif e["fn"] == "ged":
                values = ged(fx, e.get("labels", "none"))
            else:
                values = None
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (similarity): wrote {count} golden files")


if __name__ == "__main__":
    main()
