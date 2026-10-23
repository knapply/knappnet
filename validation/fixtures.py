#!/usr/bin/env python3
"""Generate the validation fixture corpus.

Each fixture is written to validation/fixtures/<name>/ as:
  edges.csv  - columns `from,to` (and `weight` when weighted), 0-based indices
  meta.json  - node count, directedness, tags, and provenance

Fixtures are deterministic: rerunning this script reproduces them exactly.
Classic datasets come from NetworkX's bundled copies; see docs/references.bib
for the original sources.

    python3 validation/fixtures.py
"""

from __future__ import annotations

import csv
import json
import random
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path

import networkx as nx

OUT = Path(__file__).resolve().parent / "fixtures"


@dataclass
class Fixture:
    name: str
    n: int
    edges: list[tuple[int, int]]
    directed: bool
    weights: list[float] | None = None
    source: str = "synthetic"
    tags: list[str] = field(default_factory=list)
    # Two-mode fixtures: 0 (first mode) or 1 (second mode) for every node,
    # written to meta.json as `type`.
    types: list[int] | None = None

    def finalize(self) -> "Fixture":
        tags = set(self.tags)
        if self.weights is not None:
            tags.add("weighted")
            if all(float(w).is_integer() for w in self.weights):
                tags.add("integer_weights")
        if any(a == b for a, b in self.edges):
            tags.add("loops")
        seen = set()
        for a, b in self.edges:
            key = (a, b) if self.directed else (min(a, b), max(a, b))
            if key in seen:
                tags.add("multi_edges")
            seen.add(key)
        g = nx.Graph()
        g.add_nodes_from(range(self.n))
        g.add_edges_from(self.edges)
        if self.n > 0 and not nx.is_connected(g):
            tags.add("disconnected")
        self.tags = sorted(tags)
        return self


def from_networkx(name: str, g: nx.Graph, source: str, weight: str | None = None) -> Fixture:
    index = {v: i for i, v in enumerate(g.nodes())}
    edges, weights = [], []
    for a, b, data in g.edges(data=True):
        edges.append((index[a], index[b]))
        if weight:
            weights.append(float(data[weight]))
    return Fixture(name, g.number_of_nodes(), edges, g.is_directed(),
                   weights if weight else None, source).finalize()


def gnp(name: str, n: int, p: float, directed: bool, seed: int,
        weighted: bool = False) -> Fixture:
    rng = random.Random(seed)
    edges = []
    for a in range(n):
        for b in range(n):
            if a == b or (not directed and b <= a):
                continue
            if rng.random() < p:
                edges.append((a, b))
    # Integer weights avoid floating-point ties being resolved differently by
    # implementations that compare path lengths exactly.
    weights = [float(rng.randint(1, 5)) for _ in edges] if weighted else None
    return Fixture(name, n, edges, directed, weights,
                   f"G(n={n}, p={p}), seed {seed}").finalize()


def uniform(name: str, n: int, p: float, directed: bool, seed: int, w: float) -> Fixture:
    """Every edge has the same weight, which must give unweighted results."""
    base = gnp(name, n, p, directed, seed)
    base.weights = [w] * len(base.edges)
    base.source += f", every weight {w}"
    return base.finalize()


def real_weights(name: str, n: int, p: float, directed: bool, seed: int) -> Fixture:
    """Non-integer weights, which take the general (Dijkstra) path."""
    base = gnp(name, n, p, directed, seed)
    rng = random.Random(seed + 1000)
    # Distinct-looking reals with few exact ties between path sums.
    base.weights = [round(rng.uniform(0.1, 3.0), 6) for _ in base.edges]
    base.source += ", real weights in [0.1, 3]"
    return base.finalize()


def messy(name: str, n: int, m: int, directed: bool, seed: int) -> Fixture:
    """Random multigraph with self-loops and repeated edges."""
    rng = random.Random(seed)
    edges = [(rng.randrange(n), rng.randrange(n)) for _ in range(m)]
    edges += edges[: m // 5]  # guarantee repeats
    return Fixture(name, n, edges, directed, None,
                   f"random multigraph n={n}, m={len(edges)}, seed {seed}").finalize()


def davis_twomode() -> Fixture:
    """Davis's Southern Women as a two-mode network: women first (type 0),
    then events (type 1)."""
    g = nx.davis_southern_women_graph()
    f = from_networkx("davis_twomode", g,
                      "Davis, Gardner, and Gardner (1941); davis1941deep; women type 0, events type 1")
    f.types = [int(g.nodes[v]["bipartite"]) for v in g.nodes()]
    f.tags = sorted(set(f.tags) | {"twomode"})
    return f


def random_twomode(name: str, n1: int, n2: int, p: float, seed: int,
                   weighted: bool = False) -> Fixture:
    """A connected random bipartite graph. Node order is shuffled so the two
    modes are interleaved, and each edge's orientation is random; the first
    seed from `seed` upward that gives a connected graph is used."""
    for s in range(seed, seed + 1000):
        rng = random.Random(s)
        types = [0] * n1 + [1] * n2
        rng.shuffle(types)
        first = [v for v, t in enumerate(types) if t == 0]
        second = [v for v, t in enumerate(types) if t == 1]
        edges = []
        for a in first:
            for b in second:
                if rng.random() < p:
                    edges.append((a, b) if rng.random() < 0.5 else (b, a))
        g = nx.Graph()
        g.add_nodes_from(range(n1 + n2))
        g.add_edges_from(edges)
        if nx.is_connected(g):
            weights = [float(rng.randint(1, 4)) for _ in edges] if weighted else None
            source = f"random bipartite graph, {n1} + {n2} nodes, p={p}, seed {s}"
            return Fixture(name, n1 + n2, edges, False, weights, source,
                           tags=["twomode"], types=types).finalize()
    raise RuntimeError(f"no connected graph for {name}")


def build() -> list[Fixture]:
    fx = [
        from_networkx("karate", nx.karate_club_graph(), "Zachary (1977); zachary1977information"),
        from_networkx("florentine", nx.florentine_families_graph(),
                      "Breiger and Pattison (1986); breiger1986cumulated. See also padgett1993robust"),
        from_networkx("southern_women", nx.davis_southern_women_graph(),
                      "Davis, Gardner, and Gardner (1941); davis1941deep"),
        from_networkx("kite", nx.krackhardt_kite_graph(), "Krackhardt (1990); krackhardt1990assessing"),
        from_networkx("les_miserables", nx.les_miserables_graph(),
                      "Knuth (1993); knuth1993stanford", weight="weight"),
        gnp("gnp_dir_sparse", 60, 0.05, True, 11),
        gnp("gnp_dir_dense", 40, 0.5, True, 12),
        gnp("gnp_und_sparse", 80, 0.04, False, 13),
        gnp("gnp_und_dense", 30, 0.6, False, 14),
        gnp("weighted_dir", 50, 0.1, True, 15, weighted=True),
        gnp("weighted_und", 50, 0.1, False, 16, weighted=True),
        uniform("uniform_weights", 40, 0.12, False, 21, 2.5),
        real_weights("real_weights_dir", 45, 0.1, True, 22),
        real_weights("real_weights_und", 45, 0.1, False, 23),
        messy("multi_loops_dir", 20, 60, True, 17),
        messy("multi_loops_und", 20, 60, False, 18),
        Fixture("empty", 5, [], True, tags=["edge_case"]).finalize(),
        Fixture("single_node", 1, [], False, tags=["edge_case"]).finalize(),
        Fixture("pair", 2, [(0, 1)], True, tags=["edge_case"]).finalize(),
        Fixture("two_triangles_isolate", 7, [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)],
                False, tags=["edge_case"]).finalize(),
        Fixture("path10", 10, [(i, i + 1) for i in range(9)], False).finalize(),
        Fixture("star8", 8, [(0, i) for i in range(1, 8)], False).finalize(),
        Fixture("cycle7_dir", 7, [(i, (i + 1) % 7) for i in range(7)], True).finalize(),
        Fixture("complete6", 6, [(a, b) for a in range(6) for b in range(a + 1, 6)],
                False).finalize(),
        Fixture("complete5_dir", 5, [(a, b) for a in range(5) for b in range(5) if a != b],
                True).finalize(),
        davis_twomode(),
        random_twomode("twomode_random", 12, 9, 0.3, 31),
        random_twomode("twomode_weighted", 9, 6, 0.45, 41, weighted=True),
    ]
    return fx


def write(fixtures: list[Fixture]) -> None:
    shutil.rmtree(OUT, ignore_errors=True)
    for f in fixtures:
        d = OUT / f.name
        d.mkdir(parents=True)
        with open(d / "edges.csv", "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["from", "to"] + (["weight"] if f.weights is not None else []))
            for i, (a, b) in enumerate(f.edges):
                w.writerow([a, b] + ([f.weights[i]] if f.weights is not None else []))
        meta = {"n": f.n, "directed": f.directed, "weighted": f.weights is not None,
                "tags": f.tags, "source": f.source}
        if f.types is not None:
            meta["type"] = f.types
        (d / "meta.json").write_text(json.dumps(meta, indent=2) + "\n", newline="")


if __name__ == "__main__":
    if {"-h", "--help"} & set(sys.argv[1:]):
        print(__doc__.strip())
        sys.exit(0)
    fixtures = build()
    write(fixtures)
    print(f"Wrote {len(fixtures)} fixtures to {OUT}")
