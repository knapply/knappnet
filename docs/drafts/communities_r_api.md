# R wrappers for community detection (draft)

Drafts of the R functions and help pages for the community detection
engine (`crates/engine/src/metrics/partition/`, C ABI in
`crates/capi/src/communities.rs`). They are written against the `net_`
naming of the redesigned R API and should be adapted to its helpers when
that branch lands. Every function takes a network `x`, validates its
arguments with the existing helpers (`check_netframe()`,
`resolve_additive_weights()`, `resolve_threads()`, `check_seed()`,
`check_choice()`, `check_flag()`, `engine_call()`), and returns results
aligned with `x$nodes`. Memberships are integer vectors (1-based in R:
add 1 to the engine's 0-based labels) carrying attributes for the scalar
results, as `cluster_louvain()` does today.

C glue routines (one per function, in `r/src/communities.c`) call the
`net_*` functions declared in the "Community detection" section of
`netcapi.h`. Registration lines for `init.c` are listed at the end.

## net_walktrap

```r
net_walktrap(x, steps = 4L, weighted = NULL, threads = 0L)
```

- `steps`: length of the random walks (t), a positive integer.
- `weighted`: use edge weights (must be positive); `NULL` uses them when
  the network has them.
- `threads`: threads for the walk vectors (results do not depend on it).

Value: membership of the level of highest modularity, with attributes
`modularity`, `merges` (a two-column integer matrix in igraph's
numbering, 1-based leaves as in `igraph::cluster_walktrap()$merges`) and
`modularity_path` (modularity before any merge and after each).

Rd draft:

```
\name{net_walktrap}
\alias{net_walktrap}
\title{Walktrap communities}
\description{
Agglomerative communities from random-walk distances (Pons and Latapy
2006): groups are merged by Ward's criterion on the distance between the
distributions of walks of \code{steps} steps, and the level of highest
modularity is returned.
}
\usage{
net_walktrap(x, steps = 4L, weighted = NULL, threads = 0L)
}
\arguments{
\item{x}{An undirected network.}
\item{steps}{Length of the random walks.}
\item{weighted}{Use edge weights (positive)? \code{NULL} uses them when present.}
\item{threads}{Number of threads (0 = all cores); results do not depend on it.}
}
\details{
Every node gets a loop of weight equal to the mean weight of its other
edges (1 for an isolated node), as igraph documents. Parallel edges add
their weights. Only adjacent groups merge, so a network with c components
has n - c merges. Among merges with equal Ward distance (to a relative
1e-9) the pair with the smaller group names is merged first; igraph
documents no tie rule, so on networks with structurally equivalent nodes
the merge order can differ from igraph's. Directed networks are rejected
(igraph ignores directions). The walk vectors are dense: memory grows as
8 n^2 bytes.
}
\value{Integer membership with attributes \code{modularity}, \code{merges},
and \code{modularity_path}.}
\references{
Pons, P. and Latapy, M. (2006). Computing communities in large networks
using random walks. \emph{Journal of Graph Algorithms and Applications},
10(2), 191--218. \doi{10.7155/jgaa.00124}
}
\examples{
x <- netframe(data.frame(from = c(1, 1, 2, 4, 4, 5, 3), to = c(2, 3, 3, 5, 6, 6, 4)),
              directed = FALSE)
net_walktrap(x, threads = 1)
}
```

## net_label_prop

```r
net_label_prop(x, variant = c("dominance", "retention", "fast", "semi_synchronous"),
               mode = c("all", "out", "in"), weighted = NULL, initial = NULL,
               fixed = NULL, seed = NULL)
```

- `variant`: update and stopping rule. `"dominance"` (igraph's default):
  nodes draw among their dominant labels, stop when every label is
  dominant. `"retention"` (NetworkX's `asyn_lpa_communities()`): keep a
  dominant label, stop after a sweep without change. `"fast"` (Traag and
  Subelj 2023; igraph's `lpa_variant = "fast"`, NetworkX's
  `fast_label_propagation_communities()`). `"semi_synchronous"` (Cordasco
  and Gargano 2010; deterministic, matches NetworkX's
  `label_propagation_communities()`; undirected networks only).
- `mode`: on directed networks, `"out"` propagates labels along edges (a
  node reads its in-neighbors), `"in"` against them, `"all"` ignores
  direction (igraph's advice for directed networks).
- `initial`: starting labels aligned with `x$nodes` (`NA` = unlabeled);
  `fixed`: logical, nodes whose label never changes (must be labeled).
- `seed`: reproducible on any platform and thread count.

Value: integer membership.

Examples must show each variant on one small network, naming the tool
each reproduces:

```r
x <- netframe(data.frame(from = c(1, 1, 2, 3, 4, 4, 5, 6), to = c(2, 3, 3, 4, 5, 6, 6, 7)),
              directed = FALSE)
net_label_prop(x, "dominance", seed = 1)         # igraph cluster_label_prop()
net_label_prop(x, "retention", seed = 1)         # NetworkX asyn_lpa_communities()
net_label_prop(x, "fast", seed = 1)              # NetworkX fast_label_propagation_communities()
net_label_prop(x, "semi_synchronous")            # NetworkX label_propagation_communities()
```

References: Raghavan, Albert and Kumara (2007) \doi{10.1103/PhysRevE.76.036106};
Cordasco and Gargano (2010) \doi{10.1109/BASNA.2010.5730298}; Traag and
Subelj (2023) \doi{10.1038/s41598-023-29610-z}.

## net_infomap

```r
net_infomap(x, weighted = NULL, node_weights = NULL, trials = 10L, seed = NULL)
net_map_equation(x, membership, weighted = NULL, node_weights = NULL)
```

- `node_weights`: teleportation weights for directed networks (numeric,
  aligned with `x$nodes`; no effect on undirected networks).
- `trials`: number of independent searches; the shortest description is
  kept.

Value: membership with attribute `codelength` (bits per step);
`net_map_equation()` returns the codelength of a given partition.

Details to state: directed networks use a random walk that teleports with
probability 0.15 (always from nodes without out-links) and teleportation
steps are encoded, as in PageRank; undirected networks have no
teleportation; an undirected network without edges has codelength `NA`.
References: Rosvall and Bergstrom (2008) \doi{10.1073/pnas.0706851105};
Rosvall, Axelsson and Bergstrom (2009) \doi{10.1140/epjst/e2010-01179-1}.

## net_leiden

```r
net_leiden(x, objective = c("modularity", "cpm"), resolution = 1, beta = 0.01,
           n_iterations = -1L, weighted = NULL, node_weights = NULL,
           initial = NULL, seed = NULL)
```

- `objective`: `"modularity"` (resolution divided by 2m, as igraph does)
  or `"cpm"` (Constant Potts Model, igraph's default
  `objective_function = "CPM"`; note the different default here, chosen
  because modularity is what most users expect; state it in the Rd).
- `beta`: randomness of the refinement step.
- `n_iterations`: iterations; negative runs until stable.
- `initial`: starting membership.

Value: membership with attributes `quality` and `iterations`. Every group
is connected (guaranteed by the method). Example showing both objectives:

```r
net_leiden(x, "modularity", seed = 1)               # igraph objective_function = "modularity"
net_leiden(x, "cpm", resolution = 0.1, seed = 1)    # igraph default objective ("CPM")
```

References: Traag, Waltman and van Eck (2019) \doi{10.1038/s41598-019-41695-z};
Traag, Van Dooren and Nesterov (2011) \doi{10.1103/PhysRevE.84.016114}.

## net_spinglass

```r
net_spinglass(x, spins = 25L, parupdate = FALSE, start_temp = 1, stop_temp = 0.01,
              cool_fact = 0.99, update_rule = c("config", "simple"), gamma = 1,
              implementation = c("orig", "neg"), gamma_minus = 1, weighted = NULL,
              seed = NULL)
```

Arguments follow `igraph::cluster_spinglass()`. `implementation = "neg"`
uses the Hamiltonian of Traag and Bruggeman (2009) and accepts negative
weights (`weighted` must then be resolved without the non-negativity
check, so use `resolve_weighting()` rather than
`resolve_additive_weights()`). Connected networks only; direction ignored.

Value: membership with attributes `modularity` (`NA` with negative
weights) and `energy`. Examples: one call per `update_rule` and per
`implementation`, noting that `"config"`/`"orig"` is igraph's default.
References: Reichardt and Bornholdt (2006) \doi{10.1103/PhysRevE.74.016110};
Traag and Bruggeman (2009) \doi{10.1103/PhysRevE.80.036115}; Kirkpatrick,
Gelatt and Vecchi (1983) \doi{10.1126/science.220.4598.671}.

## net_optimal

```r
net_optimal(x, resolution = 1, weighted = NULL)
```

Exact maximum modularity by branch and bound. Errors for components of
more than 60 nodes or after 10^9 search steps (a few minutes) rather
than return a partition that might not be optimal. Directed networks use
directed modularity. Value: membership with attribute `modularity`.
References: Brandes et al. (2008) \doi{10.1109/TKDE.2007.190689};
Verfaillie, Lemaitre and Schiex (1996) (Russian doll search; AAAI-96,
181--187, no DOI).

## net_fluid

```r
net_fluid(x, k, seed = NULL)
```

Connected undirected networks; `1 <= k <= n`; weights are not used.
Value: membership with exactly `k` groups. Reference: Pares et al. (2018)
\doi{10.1007/978-3-319-72150-7_19}.

## net_voronoi

```r
net_voronoi(x, lengths = NULL, weighted = NULL, mode = c("out", "in", "all"),
            radius = NULL, threads = 0L)
```

- `lengths`: name of an edge column (or numeric vector aligned with the
  edges) giving edge lengths; `NULL` for unit lengths. Lengths are
  divided by the edge clustering coefficient.
- `weighted`: weights for strengths (density) and modularity.
- `radius`: `NULL` chooses the radius of highest modularity among every
  distinct distance (igraph uses Brent's method instead, so it can choose
  a radius of lower modularity).

Value: membership with attributes `generators` (node indices, in the
order taken), `modularity`, and `radius`. References: Deritei et al.
(2014) \doi{10.1088/1367-2630/16/6/063007}; Molnar et al. (2024)
\doi{10.1038/s41598-024-58624-4}; Radicchi et al. (2004)
\doi{10.1073/pnas.0400054101}.

## net_compare_partitions

```r
net_compare_partitions(a, b, method = c("vi", "nmi", "split_join", "rand",
                                        "adjusted_rand", "ami"))
```

Two label vectors of equal length (any atomic type; `NA` is an error).
Matches `igraph::compare()` for the methods it offers (`"split_join"` is
the sum of the two distances of `igraph::split_join_distance()`); `"ami"`
is the adjusted mutual information of Vinh, Epps and Bailey (2010) with
the arithmetic-mean normalization. Rand indices are `NA` for fewer than
two nodes (igraph errors). References: Meila (2007)
\doi{10.1016/j.jmva.2006.11.013}; Danon et al. (2005)
\doi{10.1088/1742-5468/2005/09/P09008}; van Dongen (2000) (CWI report
INS-R0012); Rand (1971) \doi{10.1080/01621459.1971.10482356}; Hubert and
Arabie (1985) \doi{10.1007/BF01908075}; Vinh, Epps and Bailey (2010)
(JMLR 11, 2837--2854).

## net_cut_dendrogram

```r
net_cut_dendrogram(merges, k = NULL, steps = NULL, x = NULL, weighted = NULL)
```

`merges` is a two-column matrix in igraph's numbering (1-based), such as
the `merges` attribute of `net_walktrap()`, or the converted dendrogram
of `cluster_fast_greedy()` / `cluster_girvan_newman()` (via
`net_merges_from_pairs`, C routine `net_merges_from_pairs`). Exactly one
of `k` (number of groups, as `igraph::cut_at(no = )`), `steps` (number of
merges to perform, as `igraph::cut_at(steps = )`), or neither with `x`
given (cut at the level of highest modularity of `x`) is used. Cutting an
incomplete dendrogram into fewer groups than it allows is an error, as in
igraph. Value: integer membership (with attribute `modularity` when cut
at the maximum). An `as.hclust()` conversion can be written in R from the
merges (negate leaves, shift clusters by n) without engine support.

## Shared-file edits needed when the wrappers are added

`r/src/init.c` (declarations and `R_CallMethodDef` entries; argument
counts follow the glue routines):

```
{"C_walktrap", (DL_FUNC) &C_walktrap, 4},
{"C_label_propagation", (DL_FUNC) &C_label_propagation, 7},
{"C_infomap", (DL_FUNC) &C_infomap, 5},
{"C_map_equation", (DL_FUNC) &C_map_equation, 4},
{"C_leiden", (DL_FUNC) &C_leiden, 9},
{"C_spinglass", (DL_FUNC) &C_spinglass, 12},
{"C_optimal_modularity", (DL_FUNC) &C_optimal_modularity, 3},
{"C_fluid_communities", (DL_FUNC) &C_fluid_communities, 3},
{"C_voronoi", (DL_FUNC) &C_voronoi, 6},
{"C_compare_partitions", (DL_FUNC) &C_compare_partitions, 3},
{"C_cut_dendrogram", (DL_FUNC) &C_cut_dendrogram, 5},
{"C_dendrogram_modularity", (DL_FUNC) &C_dendrogram_modularity, 5},
{"C_merges_from_pairs", (DL_FUNC) &C_merges_from_pairs, 4},
```

`NAMESPACE`: `export(net_walktrap, net_label_prop, net_infomap,
net_map_equation, net_leiden, net_spinglass, net_optimal, net_fluid,
net_voronoi, net_compare_partitions, net_cut_dendrogram)`.

pkgdown reference topics: all of the above (one Rd each, with
`net_map_equation` documented with `net_infomap`).
