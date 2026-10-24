"""Every implemented function against NetworkX's own implementation.

Each case calls a NetworkX function twice, once with backend="knappnet" and
once with backend="networkx", on graphs that exercise the conventions that
differ between tools: self-loops, direction, weights, isolated nodes,
disconnected graphs, and node labels of mixed types. A call the backend
declines (NetworkX would then run its own implementation) is skipped, but
the cases below are chosen so that most calls run on the engine, and
test_coverage_of_cases checks that every implemented function ran at least
once.
"""

from __future__ import annotations

import itertools
import math
import random
from collections.abc import Iterator

import networkx as nx
import pytest

import knappnet.backend as kb

RAN = set()


def random_graph(n, p, directed, seed, loops=False, weights=None, labels="int"):
    rng = random.Random(seed)
    G = nx.DiGraph() if directed else nx.Graph()
    names = list(range(n))
    if labels == "str":
        names = [f"v{i}" for i in range(n)]
    elif labels == "mixed":
        names = [i if i % 2 else f"n{i}" for i in range(n)]
    rng.shuffle(names)
    G.add_nodes_from(names)
    for u, v in itertools.product(names, repeat=2):
        if u == v and not loops:
            continue
        if not directed and str(u) > str(v):
            continue
        if rng.random() < p:
            if weights == "int":
                G.add_edge(u, v, weight=rng.randint(1, 5))
            elif weights == "float":
                G.add_edge(u, v, weight=rng.uniform(0.5, 3.0))
            else:
                G.add_edge(u, v)
    return G


def graphs():
    out = []
    seed = 0
    for directed in (False, True):
        for n, p in ((1, 0.0), (2, 1.0), (7, 0.3), (12, 0.25), (30, 0.12), (40, 0.3)):
            for loops in (False, True):
                for weights in (None, "int", "float"):
                    seed += 1
                    labels = ("int", "str", "mixed")[seed % 3]
                    G = random_graph(n, p, directed, seed, loops, weights, labels)
                    out.append(pytest.param(G, id=f"{'D' if directed else 'U'}{n}-p{p}"
                                            f"-{'loops' if loops else 'noloops'}-{weights}"))
    out.append(pytest.param(nx.karate_club_graph(), id="karate"))
    out.append(pytest.param(nx.les_miserables_graph(), id="lesmis"))
    out.append(pytest.param(nx.florentine_families_graph(), id="florentine"))
    out.append(pytest.param(nx.path_graph(9), id="path"))
    out.append(pytest.param(nx.cycle_graph(8, create_using=nx.DiGraph), id="dicycle"))
    out.append(pytest.param(nx.complete_graph(6), id="complete"))
    out.append(pytest.param(nx.star_graph(7), id="star"))
    out.append(pytest.param(nx.barbell_graph(5, 2), id="barbell"))
    out.append(pytest.param(nx.gnp_random_graph(60, 0.08, seed=3, directed=True), id="gnp60d"))
    M = nx.MultiGraph([(0, 1), (0, 1), (1, 2), (2, 2), (2, 3)])
    out.append(pytest.param(M, id="multigraph"))
    return out


GRAPHS = graphs()


def same(a, b, path="result"):
    """Compares a backend result with NetworkX's (unordered where NetworkX's
    order is arbitrary)."""
    if isinstance(b, Iterator) or isinstance(a, Iterator):
        a, b = list(a), list(b)
        if a and isinstance(b[0], tuple) and len(b[0]) == 2 and isinstance(b[0][1], dict):
            return same(dict(a), dict(b), path)
        return same_collection(a, b, path)
    if isinstance(b, dict) or type(b).__name__ == "dict_values":
        if not isinstance(b, dict):
            return same_collection(list(a), list(b), path)
        assert set(a) == set(b), f"{path}: keys differ: {set(a) ^ set(b)}"
        for k in b:
            same(a[k], b[k], f"{path}[{k!r}]")
        return None
    if isinstance(b, tuple) and isinstance(a, tuple):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            same(x, y, f"{path}[{i}]")
        return None
    if isinstance(b, (set, frozenset)):
        assert set(a) == set(b), f"{path}: {a} != {b}"
        return None
    if isinstance(b, list):
        return same_collection(a, b, path)
    if isinstance(b, bool) or isinstance(a, bool):
        assert a == b, f"{path}: {a} != {b}"
        return None
    if isinstance(b, (int, float)):
        if isinstance(b, float) and math.isnan(b):
            assert math.isnan(a), f"{path}: {a} != nan"
            return None
        assert math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-9), f"{path}: {a} != {b}"
        return None
    assert a == b, f"{path}: {a!r} != {b!r}"
    return None


def same_collection(a, b, path):
    def norm(x):
        if isinstance(x, (set, frozenset, list)):
            return frozenset(x)
        if isinstance(x, tuple):
            return frozenset(x)
        return x

    a, b = list(a), list(b)
    assert len(a) == len(b), f"{path}: {len(a)} items != {len(b)}"
    assert sorted(map(repr, map(norm, a))) == sorted(map(repr, map(norm, b))) or \
        {norm(x) for x in a} == {norm(x) for x in b}, f"{path}: {a} != {b}"
    return None


def run(fn, G, *args, **kwargs):
    name = getattr(fn, "name", fn.__name__)
    try:
        expected = fn(G, *args, backend="networkx", **kwargs)
        if isinstance(expected, Iterator):
            expected = list(expected)
    except nx.PowerIterationFailedConvergence:
        # NetworkX stopped before converging; the engine iterates to its own
        # tolerance and returns the converged scores (python/docs/differences.md).
        return
    except Exception:  # NetworkX raises: the backend must decline
        with pytest.raises(Exception):
            res = fn(G, *args, backend="knappnet", **kwargs)
            if isinstance(res, Iterator):
                list(res)
        return
    try:
        got = fn(G, *args, backend="knappnet", **kwargs)
    except NotImplementedError:
        return
    RAN.add(name)
    same(got, expected, name)


def weight_of(G):
    return "weight" if any("weight" in d for _, _, d in G.edges(data=True)) else None


def node_sample(G, k, seed=1):
    nodes = list(G)
    random.Random(seed).shuffle(nodes)
    return nodes[:k]


@pytest.mark.parametrize("G", GRAPHS)
def test_centrality(G):
    w = weight_of(G)
    run(nx.degree_centrality, G)
    if G.is_directed():
        run(nx.in_degree_centrality, G)
        run(nx.out_degree_centrality, G)
    for weight in {None, w}:
        for normalized in (True, False):
            run(nx.betweenness_centrality, G, normalized=normalized, weight=weight)
            run(nx.edge_betweenness_centrality, G, normalized=normalized, weight=weight)
        for wf in (True, False):
            run(nx.closeness_centrality, G, distance=weight, wf_improved=wf)
        run(nx.harmonic_centrality, G, distance=weight)
        run(nx.harmonic_centrality, G, nbunch=node_sample(G, 3), distance=weight)
        run(nx.pagerank, G, weight=weight, tol=1e-12)
        run(nx.eigenvector_centrality_numpy, G, weight=weight)
        run(nx.katz_centrality_numpy, G, alpha=0.05, weight=weight)
        if not G.is_multigraph():
            run(nx.eigenvector_centrality, G, weight=weight, max_iter=1000, tol=1e-10)
            run(nx.katz_centrality, G, alpha=0.05, weight=weight, tol=1e-12)
            run(nx.current_flow_closeness_centrality, G, weight=weight)
            run(nx.constraint, G, weight=weight)
            run(nx.effective_size, G, weight=weight)
            run(nx.effective_size, G, nodes=list(G), weight=weight)
    pers = {v: i % 3 + 1 for i, v in enumerate(G)}
    run(nx.pagerank, G, personalization=pers, tol=1e-12)
    run(nx.pagerank, G, personalization=pers, dangling={v: 1 for v in G}, tol=1e-12)
    run(nx.pagerank, G, alpha=0.6, tol=1e-12)
    run(nx.katz_centrality_numpy, G, alpha=0.05, beta={v: 1 + i % 2 for i, v in enumerate(G)})
    if not G.is_multigraph():
        run(nx.hits, G, max_iter=10000, tol=1e-14)
        run(nx.hits, G, max_iter=10000, tol=1e-14, normalized=False)


def group_betweenness_definition(G, group, normalized, weight):
    """Sum over pairs s, t outside the group of the share of shortest s-t
    paths through a group member (Everett and Borgatti 1999)."""
    group = set(group)
    outside = [v for v in G if v not in group]
    total = 0.0
    for s in outside:
        for t in outside:
            if s == t or (not G.is_directed() and str(s) > str(t)):
                continue
            if not nx.has_path(G, s, t):
                continue
            paths = list(nx.all_shortest_paths(G, s, t, weight=weight))
            through = sum(1 for p in paths if group & set(p[1:-1]))
            total += through / len(paths)
    k = len(outside)
    if normalized:
        pairs = k * (k - 1) if G.is_directed() else k * (k - 1) / 2
        return total / pairs if pairs else 0.0
    return total


@pytest.mark.parametrize("G", GRAPHS)
def test_group_centrality(G):
    if len(G) < 4 or G.is_multigraph():
        return
    w = weight_of(G)
    S = node_sample(G, 2)
    run(nx.group_degree_centrality, G, S)
    if G.is_directed():
        run(nx.group_in_degree_centrality, G, S)
        run(nx.group_out_degree_centrality, G, S)
    for weight in {None, w}:
        run(nx.group_closeness_centrality, G, S, weight=weight)
        # NetworkX 3.7's group betweenness depends on set iteration order
        # (PYTHONHASHSEED) for some groups, so the backend is checked against
        # a direct evaluation of the definition instead.
        S3 = node_sample(G, 3, 2)
        for normalized in (True, False):
            got = nx.group_betweenness_centrality(G, S3, normalized=normalized,
                                                  weight=weight, backend="knappnet")
            RAN.add("group_betweenness_centrality")
            assert got == pytest.approx(group_betweenness_definition(
                G, S3, normalized, weight), rel=1e-9, abs=1e-12)
        got = nx.group_betweenness_centrality(G, [S, S3], weight=weight, backend="knappnet")
        assert got == pytest.approx([group_betweenness_definition(G, S, True, weight),
                                     group_betweenness_definition(G, S3, True, weight)])


@pytest.mark.parametrize("G", GRAPHS)
def test_structure(G):
    w = weight_of(G)
    run(nx.transitivity, G)
    for weight in {None, w}:
        run(nx.clustering, G, weight=weight)
        run(nx.average_clustering, G, weight=weight)
        run(nx.average_clustering, G, weight=weight, count_zeros=False)
        run(nx.clustering, G, nodes=node_sample(G, 3), weight=weight)
        if len(G):
            run(nx.clustering, G, nodes=next(iter(G)), weight=weight)
    run(nx.overall_reciprocity, G)
    run(nx.reciprocity, G)
    run(nx.triadic_census, G)
    run(nx.core_number, G)
    run(nx.find_cliques, G)
    run(nx.articulation_points, G)
    run(nx.bridges, G)
    run(nx.has_bridges, G)
    run(nx.biconnected_components, G)
    run(nx.is_biconnected, G)
    run(nx.global_efficiency, G)


@pytest.mark.parametrize("G", GRAPHS)
def test_components(G):
    for fn in (nx.connected_components, nx.number_connected_components, nx.is_connected,
               nx.strongly_connected_components, nx.number_strongly_connected_components,
               nx.is_strongly_connected, nx.weakly_connected_components,
               nx.number_weakly_connected_components, nx.is_weakly_connected):
        run(fn, G)
    if G.is_directed():
        # Components come in NetworkX's order (condensation numbers by it).
        got = list(nx.strongly_connected_components(G, backend="knappnet"))
        assert got == list(nx.strongly_connected_components(G, backend="networkx"))
    elif len(G):
        got = list(nx.connected_components(G, backend="knappnet"))
        assert got == list(nx.connected_components(G, backend="networkx"))


@pytest.mark.parametrize("G", GRAPHS)
def test_paths(G):
    w = weight_of(G)
    run(nx.all_pairs_shortest_path_length, G)
    run(nx.all_pairs_shortest_path_length, G, cutoff=2)
    if w:
        run(nx.all_pairs_dijkstra_path_length, G, weight=w)
        run(nx.all_pairs_dijkstra_path_length, G, weight=w, cutoff=3.5)
    for s in node_sample(G, 3):
        run(nx.single_source_shortest_path_length, G, s)
        run(nx.single_source_shortest_path_length, G, s, cutoff=1)
        run(nx.single_target_shortest_path_length, G, s)
        run(nx.descendants, G, s)
        run(nx.ancestors, G, s)
        for weight in {None, w}:
            run(nx.shortest_path_length, G, source=s, weight=weight)
            run(nx.shortest_path_length, G, target=s, weight=weight)
            for t in node_sample(G, 2, seed=5):
                run(nx.shortest_path_length, G, source=s, target=t, weight=weight)
                run(nx.has_path, G, s, t)
        if w:
            run(nx.single_source_dijkstra_path_length, G, s, weight=w)
    for weight in {None, w}:
        run(nx.shortest_path_length, G, weight=weight)
        run(nx.average_shortest_path_length, G, weight=weight)
        for fn in (nx.eccentricity, nx.diameter, nx.radius, nx.center, nx.periphery):
            run(fn, G, weight=weight)


@pytest.mark.parametrize("G", GRAPHS)
def test_assortativity(G):
    if G.is_multigraph():
        return
    w = weight_of(G)
    H = G.copy()
    for i, v in enumerate(H):
        H.nodes[v]["group"] = ("a", "b", "c")[i % 3]
        H.nodes[v]["value"] = (i * 7) % 5
    for weight in {None, w}:
        run(nx.degree_assortativity_coefficient, H, weight=weight)
        run(nx.degree_pearson_correlation_coefficient, H, weight=weight)
    run(nx.attribute_assortativity_coefficient, H, "group")
    run(nx.numeric_assortativity_coefficient, H, "value")


@pytest.mark.parametrize("G", GRAPHS)
def test_communities(G):
    if G.is_multigraph() or len(G) == 0:
        return
    w = weight_of(G)
    nodes = list(G)
    blocks = [set(nodes[i::3]) for i in range(3) if nodes[i::3]]
    for weight in {None, w}:
        for res in (1, 0.5):
            run(nx.community.modularity, G, blocks, weight=weight, resolution=res)
    if not G.is_directed():
        # Louvain's partition depends on the random stream, so the backend's
        # result is checked for being a partition with at least the
        # modularity of a single community per component, not node by node.
        for weight in {None, w}:
            comms = nx.community.louvain_communities(G, weight=weight, seed=1,
                                                     backend="knappnet")
            RAN.add("louvain_communities")
            assert nx.community.is_partition(G, comms)
            if G.number_of_edges():
                q = nx.community.modularity(G, comms, weight=weight)
                q0 = nx.community.modularity(G, list(nx.connected_components(G)),
                                             weight=weight)
                assert q >= q0 - 1e-12


def test_bipartite():
    B = nx.complete_bipartite_graph(3, 4)
    B.remove_edge(0, 4)
    B.add_node(99)
    top = [0, 1, 2]
    run(nx.bipartite.degree_centrality, B, top)
    run(nx.bipartite.betweenness_centrality, B, top)
    run(nx.bipartite.closeness_centrality, B, top)
    run(nx.bipartite.closeness_centrality, B, top, normalized=False)


def test_coverage_of_cases():
    missing = set(kb.IMPLEMENTATIONS) - RAN
    assert not missing, f"never ran on the engine: {sorted(missing)}"
