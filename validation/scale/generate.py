#!/usr/bin/env python3
"""Generate the scale tier of the validation corpus.

The main fixtures (validation/fixtures/) are deliberately small, so the
slow definition oracles can run on all of them. This tier checks the same
measures on larger, realistic networks (about 600 to 50,000 nodes)
against NetworkX, igraph, and sna. The networks are not committed: they
are generated deterministically from fixed seeds, in the same format as
the main fixtures (edges.csv and meta.json, with tags computed the same
way), so they can be regenerated anywhere.

    python3 validation/scale/generate.py [output directory]

Generators (all from NetworkX; Hagberg, Schult, and Swart 2008):
- LFR benchmark graphs with planted communities (Lancichinetti, Fortunato,
  and Radicchi 2008), nx.LFR_benchmark_graph;
- stochastic block models (Holland, Laskey, and Leinhardt 1983),
  nx.stochastic_block_model, directed;
- power-law graphs with tunable clustering (Holme and Kim 2002),
  nx.powerlaw_cluster_graph;
- random geometric graphs (spatial networks; Penrose 2003),
  nx.random_geometric_graph, weighted by Euclidean distance;
- G(n, m) random graphs (Erdos and Renyi 1959), nx.gnm_random_graph;
- a network with several components of different kinds and isolates;
- random bipartite (two-mode) graphs, nx.bipartite.random_graph, with the
  node types recorded in meta.json as `type` like the main two-mode
  fixtures.

Weights are whole numbers (the engine's Dial path) or reals (Dijkstra's
algorithm). The output depends only on the seeds and the NetworkX version
recorded in each fixture's meta.json.

LFR graphs as NetworkX generates them can contain self-loops; these are
removed, so that the measures whose loop conventions differ between tools
are compared on these networks too (loops are covered by the main
fixtures).

Besides the tags validation/fixtures.py computes, a network with real
weights is tagged `near_ties` when two different shortest-path lengths
from some node to another differ by less than the engine's relative
tolerance of 1e-10 (crates/engine/src/metrics/betweenness.rs) without
being equal. The engine and igraph count such paths as tied; NetworkX
compares lengths exactly, so its shortest-path counts can differ there.
Random geometric graphs weighted by distance have them: three nearly
collinear points make a two-step path almost exactly as long as the
direct edge. The check takes all-pairs shortest paths, so it is made
only on networks of up to NEAR_TIES_MAX_NODES nodes; no shortest-path
counting measure runs on larger ones (see check_scale.py's budgets).
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

HERE = Path(__file__).resolve().parent
DEFAULT_OUT = HERE / "work" / "fixtures"


@dataclass
class Fixture:
    name: str
    n: int
    edges: list[tuple[int, int]]
    directed: bool
    weights: list[float] | None = None
    source: str = ""
    tags: list[str] = field(default_factory=list)
    # Node types of a two-mode network (0 or 1 per node), written to
    # meta.json as `type` like the main two-mode fixtures.
    types: list[int] | None = None

    def finalize(self) -> "Fixture":
        """Tags computed exactly as validation/fixtures.py computes them, so
        the manifest's documented differences apply in the same way."""
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
        if "weighted" in tags and "integer_weights" not in tags and has_near_ties(self):
            tags.add("near_ties")
        self.tags = sorted(tags)
        return self


TIE_TOLERANCE = 1e-10
NEAR_TIES_MAX_NODES = 5000


def has_near_ties(f: "Fixture") -> bool:
    """Whether some edge (u, v, w) and source s have d(s, u) + w within
    TIE_TOLERANCE (relative) of d(s, v) without the two being equal to
    rounding error: then a path through u is tied with the shortest paths
    to v under the engine's tolerance but not under exact comparison."""
    if f.n > NEAR_TIES_MAX_NODES:
        return False
    import numpy as np
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import dijkstra

    src = np.array([a for a, _ in f.edges])
    tgt = np.array([b for _, b in f.edges])
    w = np.array(f.weights, dtype=float)
    if not f.directed:
        src, tgt, w = np.concatenate([src, tgt]), np.concatenate([tgt, src]), np.concatenate([w, w])
    keep = src != tgt
    src, tgt, w = src[keep], tgt[keep], w[keep]
    # The lightest of parallel edges (csr_matrix would add them).
    best = {}
    for a, b, x in zip(src.tolist(), tgt.tolist(), w.tolist()):
        if x < best.get((a, b), float("inf")):
            best[(a, b)] = x
    src = np.array([k[0] for k in best], dtype=int)
    tgt = np.array([k[1] for k in best], dtype=int)
    w = np.array(list(best.values()))
    adj = csr_matrix((w, (src, tgt)), shape=(f.n, f.n))
    for start in range(0, f.n, 500):
        d = dijkstra(adj, directed=True, indices=np.arange(start, min(f.n, start + 500)))
        alt, cur = d[:, src] + w, d[:, tgt]
        with np.errstate(invalid="ignore"):
            gap = np.abs(alt - cur)
            scale = np.maximum(alt, cur)
            near = (gap <= TIE_TOLERANCE * scale) & (gap > 1e-14 * scale)
        if np.any(near & np.isfinite(cur)):
            return True
    return False


def from_graph(name: str, g, source: str, weights: str | None = None,
               seed: int = 0) -> Fixture:
    """A fixture from a NetworkX graph whose nodes are relabeled 0..n-1 in
    sorted order. `weights` is None, "integer" (1 to 10), "real" (0.1 to
    3), or "distance" (Euclidean distance between node positions)."""
    nodes = sorted(g.nodes())
    index = {v: i for i, v in enumerate(nodes)}
    edges = sorted((index[a], index[b]) if g.is_directed() or index[a] <= index[b]
                   else (index[b], index[a]) for a, b in g.edges())
    w = None
    rng = random.Random(seed)
    if weights == "integer":
        w = [float(rng.randint(1, 10)) for _ in edges]
    elif weights == "real":
        w = [rng.uniform(0.1, 3.0) for _ in edges]
    elif weights == "distance":
        pos = {index[v]: p for v, p in nx.get_node_attributes(g, "pos").items()}
        w = [max(1e-6, ((pos[a][0] - pos[b][0]) ** 2 + (pos[a][1] - pos[b][1]) ** 2) ** 0.5)
             for a, b in edges]
    if weights:
        source += f", {weights} weights (seed {seed})"
    return Fixture(name, len(nodes), edges, g.is_directed(), w, source).finalize()


def lfr(n: int, avg: int, seed: int):
    g = nx.LFR_benchmark_graph(n, 3, 1.5, 0.1, average_degree=avg, max_degree=5 * avg,
                               min_community=30, max_community=150, seed=seed)
    g.remove_edges_from(list(nx.selfloop_edges(g)))
    return g


def sbm(blocks: int, size: int, p_in: float, p_out: float, seed: int):
    p = [[p_in if i == j else p_out for j in range(blocks)] for i in range(blocks)]
    return nx.stochastic_block_model([size] * blocks, p, directed=True, seed=seed)


def twomode(name: str, n1: int, n2: int, p: float, seed: int,
            weights: str | None = None) -> Fixture:
    """A random bipartite graph with the two modes interleaved in the node
    order, as in the main two-mode fixtures."""
    g = nx.bipartite.random_graph(n1, n2, p, seed=seed)
    order = list(g.nodes())
    random.Random(seed).shuffle(order)
    h = nx.relabel_nodes(g, {v: i for i, v in enumerate(order)})
    f = from_graph(name, h, f"random bipartite graph, {n1} + {n2} nodes, p = {p}, seed {seed}",
                   weights=weights, seed=seed + 100)
    f.types = [int(g.nodes[v]["bipartite"]) for v in order]
    f.tags = sorted(set(f.tags) | {"twomode"})
    return f


def components_and_isolates(seed: int) -> Fixture:
    """Disjoint pieces of different structure plus isolated nodes."""
    parts = [
        nx.powerlaw_cluster_graph(500, 3, 0.3, seed=seed),
        nx.random_geometric_graph(300, 0.1, seed=seed + 1),
        nx.cycle_graph(50),
        nx.path_graph(40),
        nx.star_graph(30),
        nx.complete_graph(12),
    ]
    g = nx.disjoint_union_all(parts)
    g.add_nodes_from(range(g.number_of_nodes(), g.number_of_nodes() + 60))
    return from_graph("components_isolates", g,
                      f"disjoint union of six graphs and 60 isolates, seed {seed}")


def build() -> list[Fixture]:
    fx = [
        # Small enough for every measure, including sna's matrix routines.
        from_graph("sbm_dir_600", sbm(4, 150, 0.04, 0.004, 31),
                   "directed stochastic block model, 4 x 150 nodes, p = 0.04 / 0.004, seed 31"),
        # About 1,000 to 1,500 nodes: shortest-path measures against
        # NetworkX (pure Python, O(nm) per measure).
        from_graph("lfr_und_1500", lfr(1500, 8, 32),
                   "LFR benchmark, n = 1500, mean degree 8, mu = 0.1, seed 32, self-loops removed"),
        from_graph("rgg_dist_1000", nx.random_geometric_graph(1000, 0.05, seed=34),
                   "random geometric graph, n = 1000, radius 0.05, seed 34",
                   weights="distance"),
        from_graph("plc_int_1000", nx.powerlaw_cluster_graph(1000, 3, 0.3, seed=35),
                   "power-law cluster graph, n = 1000, m = 3, p = 0.3, seed 35",
                   weights="integer", seed=135),
        from_graph("gnm_dir_int_1000",
                   nx.gnm_random_graph(1000, 4000, seed=36, directed=True),
                   "directed G(n = 1000, m = 4000), seed 36", weights="integer", seed=136),
        from_graph("sbm_dir_real_1200", sbm(3, 400, 0.012, 0.0012, 37),
                   "directed stochastic block model, 3 x 400 nodes, p = 0.012 / 0.0012, seed 37",
                   weights="real", seed=137),
        components_and_isolates(38),
        twomode("twomode_int_1000", 650, 350, 0.012, 45, weights="integer"),
        # About 5,000 nodes: near-linear measures, and igraph for all.
        from_graph("lfr_und_5000", lfr(5000, 10, 39),
                   "LFR benchmark, n = 5000, mean degree 10, mu = 0.1, seed 39, self-loops removed"),
        from_graph("sbm_dir_5000", sbm(5, 1000, 0.005, 0.0003, 40),
                   "directed stochastic block model, 5 x 1000 nodes, seed 40"),
        from_graph("rgg_dist_5000", nx.random_geometric_graph(5000, 0.025, seed=41),
                   "random geometric graph, n = 5000, radius 0.025, seed 41",
                   weights="distance"),
        # Large: linear-time measures only.
        from_graph("plc_und_50000", nx.powerlaw_cluster_graph(50000, 3, 0.3, seed=42),
                   "power-law cluster graph, n = 50000, m = 3, p = 0.3, seed 42"),
        from_graph("gnm_dir_real_20000",
                   nx.gnm_random_graph(20000, 80000, seed=43, directed=True),
                   "directed G(n = 20000, m = 80000), seed 43", weights="real", seed=143),
    ]
    return fx


def write(fixtures: list[Fixture], out: Path) -> None:
    shutil.rmtree(out, ignore_errors=True)
    for f in fixtures:
        d = out / f.name
        d.mkdir(parents=True)
        with open(d / "edges.csv", "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["from", "to"] + (["weight"] if f.weights is not None else []))
            for i, (a, b) in enumerate(f.edges):
                w.writerow([a, b] + ([repr(f.weights[i])] if f.weights is not None else []))
        meta = {"n": f.n, "directed": f.directed, "weighted": f.weights is not None,
                "tags": f.tags, "source": f"{f.source}; networkx {nx.__version__}"}
        if f.types is not None:
            meta["type"] = f.types
        (d / "meta.json").write_text(json.dumps(meta, indent=2) + "\n", newline="")


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUT
    fixtures = build()
    write(fixtures, out)
    for f in fixtures:
        print(f"{f.name:<22} n = {f.n:>6}  m = {len(f.edges):>7}  "
              f"{'directed' if f.directed else 'undirected'}  {', '.join(f.tags)}")
    print(f"Wrote {len(fixtures)} scale fixtures to {out}")


if __name__ == "__main__":
    main()
