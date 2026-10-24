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
import numpy as np
import pytest

import knappnet.backend as kb

RAN = set()
nx.config.warnings_to_ignore.add("cache")


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


def same(a, b, path="result", exact=False):
    """Compares a backend result with NetworkX's (unordered where NetworkX's
    order is arbitrary; with ``exact``, lists and generators must come in the
    same order and graphs must have the same adjacency order)."""
    if isinstance(b, np.ndarray):
        assert isinstance(a, np.ndarray) and a.shape == b.shape, f"{path}: shape"
        assert np.allclose(a, b, equal_nan=True, rtol=1e-6, atol=1e-9), f"{path}: {a} != {b}"
        return None
    if isinstance(b, Iterator) or isinstance(a, Iterator):
        a, b = list(a), list(b)
        if exact:
            return same(a, b, path, exact)
        if a and isinstance(b[0], tuple) and len(b[0]) == 2 and (
                isinstance(b[0][1], dict)
                or (isinstance(b[0][1], tuple) and any(isinstance(x, dict) for x in b[0][1]))):
            return same(dict(a), dict(b), path)
        return same_collection(a, b, path)
    if isinstance(b, nx.Graph):
        return same_graph(a, b, path, exact)
    if isinstance(b, dict) or type(b).__name__ == "dict_values":
        if not isinstance(b, dict):
            return same_collection(list(a), list(b), path)
        assert set(a) == set(b), f"{path}: keys differ: {set(a) ^ set(b)}"
        for k in b:
            same(a[k], b[k], f"{path}[{k!r}]", exact)
        return None
    if isinstance(b, tuple) and isinstance(a, tuple):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            same(x, y, f"{path}[{i}]", exact)
        return None
    if isinstance(b, (set, frozenset)):
        assert set(a) == set(b), f"{path}: {a} != {b}"
        return None
    if isinstance(b, list):
        if exact:
            assert isinstance(a, list) and len(a) == len(b), f"{path}: {a} != {b}"
            for i, (x, y) in enumerate(zip(a, b)):
                same(x, y, f"{path}[{i}]", exact)
            return None
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


def same_graph(a, b, path, exact=False):
    """Same class, graph attributes, nodes (in order) and edges with their data."""
    assert type(a) is type(b), f"{path}: {type(a).__name__} != {type(b).__name__}"
    assert a.graph == b.graph, f"{path}: graph attributes differ"
    assert list(a.nodes(data=True)) == list(b.nodes(data=True)), f"{path}: nodes differ"
    if exact:
        assert list(a.edges(data=True)) == list(b.edges(data=True)), f"{path}: edges differ"
    else:
        key = lambda e: (repr(sorted(map(repr, e[:2]))) if not a.is_directed() else repr(e[:2]), repr(e[2]))
        assert sorted(a.edges(data=True), key=key) == sorted(b.edges(data=True), key=key), \
            f"{path}: edges differ"
    return None


def same_collection(a, b, path):
    def norm(x):
        if isinstance(x, (set, frozenset, list)):
            return frozenset(x)
        if isinstance(x, tuple):
            return frozenset(repr(y) if isinstance(y, dict) else y for y in x)
        return x

    a, b = list(a), list(b)
    assert len(a) == len(b), f"{path}: {len(a)} items != {len(b)}"
    assert sorted(map(repr, map(norm, a))) == sorted(map(repr, map(norm, b))) or \
        {norm(x) for x in a} == {norm(x) for x in b}, f"{path}: {a} != {b}"
    return None


def run(fn, G, *args, exact=False, **kwargs):
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
    same(got, expected, name, exact)


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
        run(nx.node_connected_component, G, s)
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
        run(nx.wiener_index, G, weight=weight)
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


# ---------------------------------------------------------------------------
# Shortest paths, centrality extras, distance indices


def graph_in(param):
    return param.values[0] if hasattr(param, "values") else param.v


def small_graphs(limit=12):
    return [g for g in GRAPHS if g.id != "multigraph" and len(graph_in(g)) <= limit]


SMALL = small_graphs()


def in_input_order(G):
    """The same graph with adjacency in node order (the engine's order=1
    traversals reproduce NetworkX only for such graphs)."""
    H = G.__class__()
    H.add_nodes_from(G.nodes(data=True))
    pos = {v: i for i, v in enumerate(G)}
    H.add_edges_from(sorted(G.edges(data=True), key=lambda e: (pos[e[0]], pos[e[1]])))
    return H


ORDERED = [pytest.param(in_input_order(graph_in(g)), id=g.id + "-ordered") for g in SMALL]


def ends(G):
    nodes = list(G)
    return (nodes[0], nodes[-1]) if nodes else (None, None)


@pytest.mark.parametrize("G", SMALL)
def test_shortest_paths(G):
    s, t = ends(G)
    if s is None:
        return
    nodes = list(G)
    for weight in {None, weight_of(G)}:
        run(nx.multi_source_dijkstra_path_length, G, nodes[:2], weight=weight)
        run(nx.multi_source_dijkstra_path_length, G, nodes[:2], cutoff=2.5, weight=weight)
        run(nx.dijkstra_path_length, G, s, t, weight=weight)
        run(nx.single_source_bellman_ford_path_length, G, s, weight=weight)
        run(nx.bellman_ford_path_length, G, s, t, weight=weight)
        run(nx.all_pairs_bellman_ford_path_length, G, weight=weight)
        run(nx.negative_edge_cycle, G, weight=weight)
        run(nx.floyd_warshall, G, weight=weight)
        run(nx.floyd_warshall_numpy, G, weight=weight)
        run(nx.floyd_warshall_predecessor_and_distance, G, weight=weight)
        run(nx.single_source_dijkstra_path, G, s, weight=weight)
        run(nx.single_source_dijkstra, G, s, weight=weight)
        run(nx.single_source_dijkstra, G, s, target=t, weight=weight)
        run(nx.dijkstra_path, G, s, t, weight=weight)
        run(nx.dijkstra_predecessor_and_distance, G, s, weight=weight)
        run(nx.bellman_ford_predecessor_and_distance, G, s, weight=weight)
        run(nx.bellman_ford_path, G, s, t, weight=weight)
        run(nx.single_source_bellman_ford_path, G, s, weight=weight)
        run(nx.single_source_bellman_ford, G, s, weight=weight)
        run(nx.single_source_bellman_ford, G, s, target=t, weight=weight)
        run(nx.bidirectional_dijkstra, G, s, t, weight=weight)
        run(nx.all_pairs_dijkstra_path, G, weight=weight)
        run(nx.all_pairs_dijkstra, G, weight=weight)
        run(nx.all_pairs_bellman_ford_path, G, weight=weight)
        run(nx.shortest_path, G, weight=weight)
        run(nx.shortest_path, G, s, weight=weight)
        run(nx.shortest_path, G, target=t, weight=weight)
        run(nx.shortest_path, G, s, t, weight=weight)
    run(nx.single_source_shortest_path, G, s)
    run(nx.single_source_shortest_path, G, s, cutoff=2)
    run(nx.single_target_shortest_path, G, s)
    run(nx.all_pairs_shortest_path, G)
    run(nx.predecessor, G, s)
    run(nx.predecessor, G, s, cutoff=2)
    run(nx.bidirectional_shortest_path, G, s, t)


@pytest.mark.parametrize("G", SMALL)
def test_centrality_extras(G):
    nodes = list(G)
    if not nodes:
        return
    even, thirds = nodes[::2], nodes[::3]
    for weight in {None, weight_of(G)}:
        for norm in (True, False):
            run(nx.betweenness_centrality_subset, G, even, thirds, normalized=norm, weight=weight)
            run(nx.edge_betweenness_centrality_subset, G, even, thirds, normalized=norm, weight=weight)
            run(nx.load_centrality, G, normalized=norm, weight=weight)
            run(nx.current_flow_betweenness_centrality, G, normalized=norm, weight=weight)
            run(nx.edge_current_flow_betweenness_centrality, G, normalized=norm, weight=weight)
            run(nx.current_flow_betweenness_centrality_subset, G, even, thirds,
                normalized=norm, weight=weight)
            run(nx.edge_current_flow_betweenness_centrality_subset, G, even, thirds,
                normalized=norm, weight=weight)
        states = {v: ((37 * i) % 11) / 10 for i, v in enumerate(nodes)}
        run(nx.percolation_centrality, G, states=states, weight=weight)
        run(nx.percolation_centrality, G, weight=weight)
        run(nx.load_centrality, G, v=nodes[0], weight=weight)
        run(nx.load_centrality, G, cutoff=2, weight=weight)
        run(nx.resistance_distance, G, weight=weight)
        run(nx.resistance_distance, G, nodes[0], weight=weight)
        run(nx.resistance_distance, G, nodeB=nodes[-1], weight=weight)
        run(nx.resistance_distance, G, nodes[0], nodes[-1], weight=weight)
        run(nx.resistance_distance, G, weight=weight, invert_weight=False)
        run(nx.effective_graph_resistance, G, weight=weight)
        run(nx.kemeny_constant, G, weight=weight)
        run(nx.second_order_centrality, G, weight=weight or "weight")
        run(nx.laplacian_centrality, G, weight=weight)
        run(nx.laplacian_centrality, G, normalized=False, weight=weight)
        run(nx.trophic_levels, G, weight=weight or "weight")
        run(nx.trophic_differences, G, weight=weight or "weight")
        run(nx.trophic_incoherence_parameter, G, weight=weight or "weight")
        run(nx.trophic_incoherence_parameter, G, weight=weight or "weight", cannibalism=True)
        run(nx.global_reaching_centrality, G, weight=weight)
        run(nx.global_reaching_centrality, G, weight=weight, normalized=False)
        run(nx.local_reaching_centrality, G, nodes[0], weight=weight)
        run(nx.local_reaching_centrality, G, nodes[0], weight=weight, normalized=False)
    run(nx.edge_load_centrality, G)
    run(nx.edge_load_centrality, G, cutoff=2)
    run(nx.communicability_betweenness_centrality, G)
    run(nx.communicability, G)
    run(nx.communicability_exp, G)
    run(nx.subgraph_centrality, G)
    run(nx.subgraph_centrality, G, normalized=True)
    run(nx.subgraph_centrality_exp, G)
    run(nx.estrada_index, G)
    run(nx.dispersion, G)
    run(nx.dispersion, G, normalized=False)
    run(nx.dispersion, G, alpha=2, b=1, c=0.5)
    run(nx.dispersion, G, u=nodes[0])
    run(nx.dispersion, G, v=nodes[0])
    run(nx.voterank, G)
    run(nx.voterank, G, 3)


@pytest.mark.parametrize("G", SMALL)
def test_distance_indices(G):
    nodes = list(G)
    for weight in {None, weight_of(G)}:
        run(nx.schultz_index, G, weight=weight)
        run(nx.gutman_index, G, weight=weight)
        run(nx.hyper_wiener_index, G, weight=weight)
        run(nx.harmonic_diameter, G, weight=weight)
        run(nx.closeness_vitality, G, weight=weight)
        if nodes:
            run(nx.closeness_vitality, G, node=nodes[0], weight=weight)
    run(nx.local_efficiency, G)


# ---------------------------------------------------------------------------
# Cores, mixing, cuts, matching, colouring


@pytest.mark.parametrize("G", [g for g in SMALL] + ORDERED)
def test_cores_and_mixing(G):
    G = G.copy()
    rng = random.Random(len(G))
    nodes = list(G)
    for v in nodes:
        G.nodes[v]["color"] = rng.choice(["r", "g", "b", 1, 2.0])
        G.nodes[v]["val"] = rng.choice([1, 2, 3.5, 7])
    sub, sub2 = nodes[::2], nodes[1::3]
    for k in (None, 0, 1, 2, 3):
        for f in (nx.k_core, nx.k_shell, nx.k_crust):
            run(f, G, k)
    for k in (0, 1, 2, 3):
        run(nx.k_corona, G, k)
    for k in (0, 2, 3, 4):
        run(nx.k_truss, G, k)
    run(nx.onion_layers, G)
    for weight in {None, weight_of(G)}:
        for x in ("in", "out"):
            for y in ("in", "out"):
                run(nx.degree_assortativity_coefficient, G, x=x, y=y, weight=weight)
                if sub:
                    run(nx.degree_assortativity_coefficient, G, x=x, y=y, weight=weight,
                        nodes=sub)
        for src in ("in", "out", "in+out"):
            for tgt in ("in", "out", "in+out"):
                run(nx.average_neighbor_degree, G, source=src, target=tgt, weight=weight)
                run(nx.average_degree_connectivity, G, source=src, target=tgt, weight=weight)
                if sub:
                    run(nx.average_neighbor_degree, G, source=src, target=tgt, nodes=sub,
                        weight=weight)
        if sub:
            run(nx.cut_size, G, sub, weight=weight)
            run(nx.cut_size, G, sub, sub2, weight=weight)
            run(nx.volume, G, sub, weight=weight)
            for f in (nx.normalized_cut_size, nx.conductance, nx.edge_expansion,
                      nx.mixing_expansion):
                run(f, G, sub, weight=weight)
                run(f, G, sub, sub2, weight=weight)
    for x in ("in", "out"):
        for y in ("in", "out"):
            run(nx.degree_mixing_dict, G, x=x, y=y)
            run(nx.degree_mixing_dict, G, x=x, y=y, normalized=True)
            if sub:
                run(nx.degree_mixing_dict, G, x=x, y=y, nodes=sub)
    run(nx.attribute_assortativity_coefficient, G, "color")
    run(nx.numeric_assortativity_coefficient, G, "val")
    run(nx.attribute_mixing_dict, G, "color")
    run(nx.attribute_mixing_dict, G, "color", normalized=True)
    if sub:
        run(nx.attribute_assortativity_coefficient, G, "color", nodes=sub)
        run(nx.numeric_assortativity_coefficient, G, "val", nodes=sub)
        run(nx.attribute_mixing_dict, G, "color", nodes=sub)
        run(nx.node_expansion, G, sub)
        run(nx.boundary_expansion, G, sub)
    mapping = {"r": 0, "g": 1, "b": 2, 1: 3, 2.0: 4}
    run(nx.attribute_mixing_matrix, G, "color", mapping=mapping)
    run(nx.attribute_mixing_matrix, G, "color", mapping=mapping, normalized=False)
    run(nx.degree_mixing_matrix, G, mapping={i: i for i in range(8)})
    run(nx.degree_mixing_matrix, G, mapping={i: i for i in range(8)}, normalized=False)
    run(nx.square_clustering, G)
    if nodes:
        run(nx.square_clustering, G, nodes[0])
        run(nx.square_clustering, G, sub)
        run(nx.generalized_degree, G)
        run(nx.generalized_degree, G, nodes[0])
        run(nx.generalized_degree, G, sub)
    run(nx.maximal_matching, G, exact=True)
    if not G.is_directed():
        mm = nx.maximal_matching(G)
        for m in (mm, set(list(mm)[:1]), {(u, v) for u, v in G.edges()}, set()):
            run(nx.is_matching, G, m)
            run(nx.is_maximal_matching, G, m)
            run(nx.is_perfect_matching, G, m)
        run(nx.is_maximal_matching, G, {(v, u) for u, v in mm})
        run(nx.is_perfect_matching, G, {(v, u) for u, v in mm})
    for strategy in ("largest_first", "saturation_largest_first", "DSATUR"):
        run(nx.greedy_color, G, strategy)


# ---------------------------------------------------------------------------
# Link prediction, hashes, similarity, isomorphism


@pytest.mark.parametrize("G", [g for g in small_graphs(10)])
def test_similarity_measures(G):
    G = G.copy()
    nodes = list(G)
    rng = random.Random(len(nodes))
    for v in nodes:
        G.nodes[v]["community"] = rng.choice([0, 1, "a"])
        G.nodes[v]["lab"] = rng.choice(["x", "y", 3])
    for e in G.edges:
        G.edges[e]["el"] = rng.choice(["p", "q"])
    if not nodes:
        return
    pairs = [(rng.choice(nodes), rng.choice(nodes)) for _ in range(5)]
    pairs = [p for p in pairs if p[0] != p[1]]
    for f in (nx.resource_allocation_index, nx.jaccard_coefficient, nx.adamic_adar_index,
              nx.preferential_attachment, nx.cn_soundarajan_hopcroft,
              nx.ra_index_soundarajan_hopcroft, nx.within_inter_cluster,
              nx.common_neighbor_centrality):
        run(f, G)
        run(f, G, pairs)
    run(nx.common_neighbor_centrality, G, pairs, alpha=0.3)
    run(nx.within_inter_cluster, G, delta=0.5)
    run(nx.cn_soundarajan_hopcroft, G, community="lab")
    for it in (1, 3):
        for ds in (1, 16):
            run(nx.weisfeiler_lehman_graph_hash, G, iterations=it, digest_size=ds)
            run(nx.weisfeiler_lehman_graph_hash, G, edge_attr="el", node_attr="lab",
                iterations=it, digest_size=ds)
            run(nx.weisfeiler_lehman_subgraph_hashes, G, iterations=it, digest_size=ds)
            run(nx.weisfeiler_lehman_subgraph_hashes, G, edge_attr="el", node_attr="lab",
                iterations=it, digest_size=ds, include_initial_labels=True)
    run(nx.simrank_similarity, G)
    run(nx.simrank_similarity, G, importance_factor=0.6, tolerance=1e-6)
    run(nx.simrank_similarity, G, nodes[0])
    run(nx.simrank_similarity, G, nodes[0], nodes[-1])


@pytest.mark.parametrize("directed, n, seed", list(itertools.product(
    (False, True), (0, 1, 3, 5, 6), range(4))))
def test_isomorphism_pairs(directed, n, seed):
    G1 = random_graph(n, 0.4, directed, seed)
    G2 = random_graph(n, 0.4, directed, seed + 100, labels="mixed")
    perm = list(G1)
    random.Random(seed).shuffle(perm)
    G3 = nx.relabel_nodes(G1, dict(zip(G1, perm)))
    for A, B in ((G1, G2), (G1, G3)):
        for v in A:
            A.nodes[v]["lab"] = "x" if hash(str(v)) % 2 else "y"
        for v in B:
            B.nodes[v]["lab"] = "x" if len(str(v)) % 2 else "y"
        run(nx.is_isomorphic, A, B)
        run(nx.could_be_isomorphic, A, B)
        run(nx.could_be_isomorphic, A, B, properties="dc")
        run(nx.fast_could_be_isomorphic, A, B)
        run(nx.faster_could_be_isomorphic, A, B)
        run(nx.vf2pp_is_isomorphic, A, B)
        run(nx.vf2pp_is_isomorphic, A, B, node_label="lab")
        run(nx.vf2pp_subgraph_is_isomorphic, A, B)
        run(nx.vf2pp_subgraph_is_isomorphic, A, B, node_label="lab")
        run(nx.vf2pp_is_monomorphic, A, B)
        run(nx.vf2pp_is_monomorphic, A, B, node_label="lab")
        if n <= 5:
            run(nx.graph_edit_distance, A, B)
            run(nx.graph_edit_distance, A, B, upper_bound=3)


# ---------------------------------------------------------------------------
# Flow, connectivity, components, trees, DAGs and cycles


def with_capacities(G, seed):
    G = G.copy()
    rng = random.Random(seed)
    for _, _, d in G.edges(data=True):
        if "weight" in d:
            d["capacity"] = d["weight"] if rng.random() < 0.85 else 0
        elif rng.random() < 0.7:
            d["capacity"] = rng.randint(1, 4)
    for v in G:
        G.nodes[v]["demand"] = 0
    nodes = list(G)
    if G.is_directed() and len(nodes) > 1:
        a, b = rng.sample(nodes, 2)
        G.nodes[a]["demand"], G.nodes[b]["demand"] = -3, 3
        for _, _, d in G.edges(data=True):
            d["weight"] = d.get("weight", rng.randint(1, 5))
    return G


FLOW = [pytest.param(with_capacities(graph_in(g), i), id=g.id)
        for i, g in enumerate(small_graphs(10))]


@pytest.mark.parametrize("G", FLOW)
def test_flow_and_connectivity(G):
    nodes = list(G)
    if len(nodes) >= 2:
        s, t = nodes[0], nodes[-1]
        for f in (nx.maximum_flow_value, nx.minimum_cut_value, nx.minimum_cut):
            run(f, G, s, t)
            run(f, G, s, t, capacity="weight")
            run(f, G, t, s)
        run(nx.connectivity.minimum_st_edge_cut, G, s, t)
        run(nx.connectivity.minimum_st_node_cut, G, s, t)
        run(nx.minimum_node_cut, G, s, t)
        run(nx.minimum_edge_cut, G, s, t)
        run(nx.node_connectivity, G, s, t)
        run(nx.edge_connectivity, G, s, t)
        run(nx.connectivity.local_node_connectivity, G, s, t)
        run(nx.connectivity.local_edge_connectivity, G, s, t)
        # NetworkX's network simplex can cycle forever on float data.
        if all(isinstance(d.get(k, 0), int) for _, _, d in G.edges(data=True)
               for k in ("weight", "capacity")):
            run(nx.min_cost_flow_cost, G)
        if not G.is_directed():
            run(nx.gomory_hu_tree, G)
            run(nx.gomory_hu_tree, G, capacity="weight")
        run(nx.connectivity.is_locally_k_edge_connected, G, s, t, 2)
        run(nx.connectivity.is_locally_k_edge_connected, G, s, t, 1)
    run(nx.node_connectivity, G)
    run(nx.edge_connectivity, G)
    run(nx.average_node_connectivity, G)
    run(nx.connectivity.all_pairs_node_connectivity, G)
    run(nx.connectivity.all_pairs_node_connectivity, G, nodes[:3])
    run(nx.local_bridges, G)
    run(nx.local_bridges, G, with_span=False)
    run(nx.condensation, G)
    run(nx.attracting_components, G)
    run(nx.number_attracting_components, G)
    run(nx.is_attracting_component, G)
    run(nx.is_semiconnected, G)
    if not G.is_directed():
        run(nx.k_components, G)
    for k in (1, 2):
        run(nx.k_edge_components, G, k)
    for k in (3, 4):
        run(nx.k_edge_components, G, k)
    run(nx.connectivity.bridge_components, G)
    for k in (1, 2, 3):
        run(nx.connectivity.is_k_edge_connected, G, k)


def random_dag(G, seed):
    rng = random.Random(seed)
    nodes = list(G)
    rank = {v: rng.random() for v in nodes}
    H = nx.DiGraph()
    H.add_nodes_from(G.nodes(data=True))
    for u, v, d in G.edges(data=True):
        if u != v:
            a, b = (u, v) if rank[u] < rank[v] else (v, u)
            if not H.has_edge(a, b):
                H.add_edge(a, b, **d)
    return H


def random_tree(n, seed):
    rng = random.Random(seed)
    names = [f"t{i}" if seed % 2 else i for i in range(n)]
    T = nx.DiGraph()
    T.add_nodes_from(names)
    for i in range(1, n):
        if rng.random() < 0.9:
            T.add_edge(names[rng.randrange(i)], names[i])
    return T


STRUCT = (SMALL + ORDERED
          + [pytest.param(random_dag(graph_in(g), i), id=g.id + "-dag")
             for i, g in enumerate(SMALL)]
          + [pytest.param(in_input_order(random_dag(graph_in(g), i)), id=g.id + "-dag-ordered")
             for i, g in enumerate(SMALL)]
          + [pytest.param(random_tree(n, s), id=f"tree{n}-{s}")
             for n in (1, 2, 5, 12) for s in range(3)])


@pytest.mark.parametrize("G", STRUCT)
def test_trees_dags_and_cycles(G):
    G = G.copy()
    rng = random.Random(len(G))
    for _, _, d in G.edges(data=True):
        if rng.random() < 0.8 and "weight" not in d:
            d["weight"] = rng.choice([1, 2, 3, 2.5, -1])
    nodes = list(G)
    run(nx.topological_sort, G, exact=True)
    run(nx.topological_generations, G, exact=True)
    run(nx.lexicographical_topological_sort, G, exact=True)
    if len(G) <= 8:
        run(nx.all_topological_sorts, G)
    run(nx.is_directed_acyclic_graph, G)
    run(nx.algorithms.dag.has_cycle, G)
    run(nx.transitive_closure, G)
    run(nx.transitive_closure, G, reflexive=None)
    run(nx.transitive_closure, G, reflexive=True)
    run(nx.transitive_reduction, G)
    for w in ("weight", None):
        run(nx.dag_longest_path_length, G, weight=w)
        run(nx.dag_longest_path, G, weight=w)
    run(nx.dag_longest_path_length, G, default_weight=2)
    if nodes:
        for root in (nodes[0], nodes[-1]):
            run(nx.immediate_dominators, G, root)
            run(nx.dominance_frontiers, G, root)
        run(nx.lowest_common_ancestor, G, nodes[0], nodes[-1])
        run(nx.lowest_common_ancestor, G, nodes[-1], nodes[0], default="x")
        run(nx.all_pairs_lowest_common_ancestor, G)
        run(nx.all_pairs_lowest_common_ancestor, G,
            [(nodes[0], nodes[-1]), (nodes[-1], nodes[0]), (nodes[0], nodes[0])])
        run(nx.cycle_basis, G, nodes[len(nodes) // 2])
    run(nx.cycle_basis, G)
    run(nx.girth, G)
    for orientation in (None, "original", "reverse", "ignore"):
        run(nx.find_cycle, G, orientation=orientation)
        if nodes:
            run(nx.find_cycle, G, nodes[0], orientation=orientation)
    if not G.is_directed() and len(G):
        run(nx.rich_club_coefficient, G, normalized=False)
    for algorithm in ("kruskal", "prim"):
        for w in ("weight", None):
            run(nx.minimum_spanning_tree, G, weight=w, algorithm=algorithm)
            run(nx.maximum_spanning_tree, G, weight=w, algorithm=algorithm)
            run(nx.minimum_spanning_edges, G, weight=w, algorithm=algorithm, exact=True)
            run(nx.maximum_spanning_edges, G, weight=w, algorithm=algorithm, data=False, exact=True)
    run(nx.algorithms.tree.mst.kruskal_mst_edges, G, True, exact=True)
    run(nx.algorithms.tree.mst.kruskal_mst_edges, G, False, data=False, exact=True)


def test_coverage_of_cases():
    missing = set(kb.IMPLEMENTATIONS) - RAN
    assert not missing, f"never ran on the engine: {sorted(missing)}"
