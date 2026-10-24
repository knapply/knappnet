#!/usr/bin/env python3
"""Definition oracle for k-shells, k-crusts, k-coronas, k-trusses, and
onion layers, written from the definitions (repeated deletion of nodes
below degree k, of edges in too few triangles, and layer-by-layer
removal). Shares no code with the engine. Degrees count each parallel
edge and a loop (twice undirected, once in and once out directed), as
the engine's core numbers do; trusses and onion layers use the simple
undirected graph without loops.

    python3 validation/oracles/cores_definition_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"


def degree_within(fx, alive):
    d = {v: 0 for v in alive}
    for a, b in zip(fx.src, fx.tgt):
        if a in alive and b in alive:
            d[a] += 1
            d[b] += 1
    return d


def core_numbers(fx):
    core = [0] * fx.n
    k = 0
    while True:
        alive = set(range(fx.n))
        while True:
            d = degree_within(fx, alive)
            drop = {v for v in alive if d[v] < k}
            if not drop:
                break
            alive -= drop
        if not alive:
            break
        for v in alive:
            core[v] = k
        k += 1
    return core


def simple_pairs(fx):
    return {frozenset((a, b)) for a, b in zip(fx.src, fx.tgt) if a != b}


def truss(fx, k):
    edges = set(simple_pairs(fx))
    while True:
        nbrs = {}
        for e in edges:
            a, b = tuple(e)
            nbrs.setdefault(a, set()).add(b)
            nbrs.setdefault(b, set()).add(a)
        drop = set()
        for e in edges:
            a, b = tuple(e)
            if len(nbrs[a] & nbrs[b]) < k - 2:
                drop.add(e)
        if not drop:
            return edges
        edges -= drop


def onion(fx):
    pairs = simple_pairs(fx)
    nbrs = {v: set() for v in range(fx.n)}
    for e in pairs:
        a, b = tuple(e)
        nbrs[a].add(b)
        nbrs[b].add(a)
    layer = {}
    current = 1
    iso = [v for v in range(fx.n) if not nbrs[v]]
    for v in iso:
        layer[v] = 1
    if iso:
        current = 2
    alive = {v for v in range(fx.n) if nbrs[v]}
    core = 1
    while alive:
        deg = {v: len(nbrs[v] & alive) for v in alive}
        core = max(core, min(deg.values()))
        this = {v for v in alive if deg[v] <= core}
        for v in this:
            layer[v] = current
        alive -= this
        current += 1
    return [layer[v] for v in range(fx.n)]


def compute(fx, e):
    fn = e["fn"]
    if fn == "core_subset":
        c = core_numbers(fx)
        which, k = e["which"], e.get("k")
        if k is None:
            if fx.n == 0:
                return []
            k = max(c) - (1 if which == "crust" else 0)
        if which == "core":
            return [v for v in range(fx.n) if c[v] >= k]
        if which == "shell":
            return [v for v in range(fx.n) if c[v] == k]
        if which == "crust":
            return [v for v in range(fx.n) if c[v] <= k]
        succ = {v: set() for v in range(fx.n)}
        for a, b in zip(fx.src, fx.tgt):
            succ[a].add(b)
            if not fx.directed:
                succ[b].add(a)
        return [v for v in range(fx.n)
                if c[v] == k and sum(1 for w in succ[v] if c[w] >= k) == k]
    if fn == "truss":
        number = {e: 2 for e in simple_pairs(fx)}
        k = 3
        while True:
            t = truss(fx, k)
            if not t:
                break
            for e in t:
                number[e] = k
            k += 1
        return [0 if a == b else number[frozenset((a, b))] for a, b in zip(fx.src, fx.tgt)]
    if fn == "k_truss":
        t = truss(fx, int(e["k"]))
        if e["part"] == "edges":
            return [i for i, (a, b) in enumerate(zip(fx.src, fx.tgt))
                    if a != b and frozenset((a, b)) in t]
        return sorted({v for x in t for v in x})
    if fn == "onion":
        return onion(fx)
    return None


def main() -> None:
    spec = load_area("cores")
    count = 0
    for name in fixture_names():
        fx = load_fixture(name)
        for variant, v in spec.items():
            if ORACLE not in v["oracles"] or not applies(v, fx):
                continue
            values = compute(fx, v["engine"])
            if values is not None:
                write_values(ORACLE, name, variant, values)
                count += 1
    write_version(ORACLE, DEFINITION_VERSION)
    print(f"definition (cores): wrote {count} golden files")


if __name__ == "__main__":
    main()
