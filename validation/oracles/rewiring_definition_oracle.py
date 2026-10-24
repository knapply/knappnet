#!/usr/bin/env python3
"""Definition oracle for edge swaps, reference graphs, and the rich-club
coefficient: the contract a rewired graph must satisfy (degrees kept,
simple, connected where required, as many edges), or [-1] where the
procedure cannot run (too few nodes or edges, a disconnected input, or,
found by exhaustive search, no swap at all); and the rich-club
coefficient by direct counting. Shares no code with the engine.

    python3 validation/oracles/rewiring_definition_oracle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (DEFINITION_VERSION, applies, fixture_names, load_area,  # noqa: E402
                    load_fixture, write_values, write_version)

ORACLE = "definition"


def pairs_of(fx):
    if fx.directed:
        return set(zip(fx.src, fx.tgt))
    return {(min(a, b), max(a, b)) for a, b in zip(fx.src, fx.tgt)}


def neighbours(fx, pairs):
    nb = {v: set() for v in range(fx.n)}
    for a, b in pairs:
        nb[a].add(b)
        if not fx.directed:
            nb[b].add(a)
    return nb


def connected(fx, pairs):
    nb = {v: set() for v in range(fx.n)}
    for a, b in pairs:
        nb[a].add(b)
        nb[b].add(a)
    seen, stack = {0}, [0]
    while stack:
        for w in nb[stack.pop()]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return fx.n > 0 and len(seen) == fx.n


def any_double_swap(fx, pairs):
    nb = neighbours(fx, pairs)
    ends = [(a, b) for a, b in pairs] + [(b, a) for a, b in pairs if a != b]
    for u, v in ends:
        for x, y in ends:
            if u != x and v != y and x not in nb[u] and y not in nb[v]:
                return True
    return False


def any_directed_swap(fx, pairs):
    nb = neighbours(fx, pairs)
    for a in range(fx.n):
        for b in nb[a]:
            for c in nb[b]:
                for d in nb[c]:
                    if a != b and b != c and c != d and c not in nb[a] and d not in nb[b] \
                            and b not in nb[c]:
                        return True
    return False


def rich_club(fx):
    if any(a == b for a, b in zip(fx.src, fx.tgt)):
        return [-1]
    pairs = pairs_of(fx)
    deg = {v: 0 for v in range(fx.n)}
    for a, b in pairs:
        deg[a] += 1
        deg[b] += 1
    if not pairs:
        return []
    out = []
    k = 0
    while True:
        rich = {v for v in range(fx.n) if deg[v] > k}
        nk = len(rich)
        if nk <= 1:
            break
        ek = sum(1 for a, b in pairs if a in rich and b in rich)
        out.append(2 * ek / (nk * (nk - 1)))
        k += 1
    return out


def compute(fx, e):
    fn = e["fn"]
    if fn == "rich_club":
        return rich_club(fx)
    if fn != "swap_contract":
        return None
    p = e["procedure"]
    pairs = pairs_of(fx)
    m = len(pairs)
    if fx.n < 4:
        return [-1]
    if p == "directed":
        if m < 3 or not any_directed_swap(fx, pairs):
            return [-1]
        return [(1, 1, 1, m)]
    if p == "double":
        if m < 2 or not any_double_swap(fx, pairs):
            return [-1]
        return [(1, 1, m)]
    if p == "connected":
        if not connected(fx, pairs):
            return [-1]
        return [(1, 1, m)]
    if m < 2:
        return [-1]
    return [(1, 1, m)]


def main() -> None:
    spec = load_area("rewiring")
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
    print(f"definition (rewiring): wrote {count} golden files")


if __name__ == "__main__":
    main()
