# Differences from NetworkX

knappnet aims to return exactly what NetworkX returns. Every implemented
function is compared with NetworkX on graphs with self-loops, weights,
direction, isolated nodes, disconnected parts, and node labels of mixed
types (`python/tests`), and NetworkX's own test suite passes with the
backend (see [testing](testing.md)). Where knappnet cannot reproduce
NetworkX's result for some input, it declines the call and NetworkX runs
it. What remains is listed here.

## Floating-point results

Sums are accumulated in a different order, so floating-point results can
differ from NetworkX's in the last digits.

PageRank, eigenvector centrality, Katz centrality (power iteration), and
HITS are solved to the engine's tolerance (1e-12), not to NetworkX's `tol`.
They agree with NetworkX to within NetworkX's tolerance and are closer to
the exact values. Where NetworkX would stop with
`PowerIterationFailedConvergence` after `max_iter` iterations, knappnet
returns the converged scores; a `max_iter` below NetworkX's default is
left to NetworkX, so that deliberately small limits still raise.

## Ties between real-valued path lengths

With real-valued weights, the engine treats two path lengths as equal when
they agree to a relative tolerance of 1e-10; NetworkX compares them
exactly. Betweenness (node, edge, and group) can therefore differ where
paths are equal in exact arithmetic but not in floating point (0.1 + 0.2
versus 0.3): knappnet counts both as shortest paths, NetworkX only one.
Integer weights and unweighted graphs are unaffected.

## Order of generated items

Functions that return dicts return them in node order. Functions that
yield collections yield the same collections as NetworkX:

- connected and weakly connected components, and strongly connected
  components, come in NetworkX's order (`condensation` numbers components
  by it);
- articulation points, bridges, biconnected components, and maximal
  cliques come in a different order, and the nodes of a clique in a
  different order within it;
- the dicts of `single_source_shortest_path_length` and related functions
  list nodes in node order rather than in order of distance.

## Number types

Unweighted distances and eccentricities are integers, as in NetworkX.
Weighted ones are floats even when every weight is an integer (NetworkX
then returns integers); they compare equal. `harmonic_centrality` gives
`0.0` where NetworkX gives `0` for a node no other node reaches.

## Louvain communities

`louvain_communities` uses the engine's random number stream: a seed
reproduces knappnet's partition, not the one NetworkX returns for that
seed. Both are Louvain partitions of the same quality on average (the
validation suite compares their mean modularity over seeds). Five
NetworkX tests that expect NetworkX's partitions for particular seeds are
marked as known differences.

## Group betweenness centrality

In NetworkX 3.7, `group_betweenness_centrality` can return different values
for the same group in different Python processes, because its result
depends on the iteration order of a set (and so on `PYTHONHASHSEED`) for
some groups of three or more nodes. knappnet returns the value of the
definition (the share of shortest paths between pairs outside the group
that pass through it, summed over those pairs), which its tests check
against a direct enumeration of shortest paths.

## Effective size

With `nodes=None` and SciPy installed, NetworkX 3.7 computes
`effective_size` with a vectorized formula that, for weighted or directed
graphs, divides each tie by the strongest tie of the other endpoint, and so
disagrees with its own per-node computation (used when `nodes` is given)
and with Burt's definition. knappnet leaves those calls to NetworkX and
runs the others, which agree with NetworkX's per-node computation.

## Functions not run on the engine

Some calls are always left to NetworkX because the engine would give a
different answer; the [coverage list](../../docs/drafts/networkx_coverage.md)
gives every case. In particular `label_propagation_communities` and
`greedy_modularity_communities` are not wired, though the engine has both
algorithms: NetworkX's results depend on set iteration order and on its
tie-breaking, which the engine does not reproduce.
