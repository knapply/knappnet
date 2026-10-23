# knappnet

Fast, validated network analysis for R, with Python to follow.

A network is two data frames, one for nodes and one for edges, with any
attributes you like. Every measure runs in a compiled, multithreaded engine
with no dependencies, and every measure is checked against [igraph](https://r.igraph.org/),
[sna](https://cran.r-project.org/package=sna), [NetworkX](https://networkx.org/), and a direct evaluation of its definition before it
ships. Where those tools follow different conventions, the difference is
documented and selectable with one argument, never silent.

**Documentation:** [Get started](https://knapply.github.io/knappnet/articles/knappnet.html) ·
[Reference](https://knapply.github.io/knappnet/reference/index.html) ·
Coming from [igraph](https://knapply.github.io/knappnet/articles/coming-from-igraph.html),
[statnet](https://knapply.github.io/knappnet/articles/coming-from-statnet.html),
[UCINET, ORA, Pajek, and Gephi](https://knapply.github.io/knappnet/articles/coming-from-ucinet-ora-pajek-gephi.html), or
[NetworkX](https://knapply.github.io/knappnet/articles/coming-from-networkx.html) ·
[Conventions](https://knapply.github.io/knappnet/articles/conventions.html) ·
[Validation](https://knapply.github.io/knappnet/articles/validation.html) ·
[Performance](https://knapply.github.io/knappnet/articles/performance.html)

> **Status:** early development. Interfaces may change before the first
> release.

## Installation

Not yet released. Once binaries are published, installation will not require
Rust or any other build tools. To build from source today you need R (4.1 or
newer) and a Rust toolchain:

<!-- norun -->
```r
# install.packages("remotes")
remotes::install_github("knapply/knappnet", subdir = "r")
```

## A first network

```r
library(knappnet)

edges <- data.frame(
  from = c("ana", "ana", "ben", "cy", "cy"),
  to   = c("ben", "cy",  "cy",  "dee", "eli")
)
g <- netframe(edges, directed = FALSE)
g
#> <netframe> undirected network with 5 nodes and 5 edges, unweighted

g$nodes$degree <- degree(g)
g$nodes$betweenness <- betweenness(g)
g$nodes
#>   name degree betweenness
#> 1  ana      2           0
#> 2  ben      2           0
#> 3   cy      4           5
#> 4  dee      1           0
#> 5  eli      1           0
```

A `netframe` is a list holding the two tables, `g$nodes` and `g$edges`,
plus a note of which columns play which role (identifier, endpoints,
weight). Nodes that only appear in the edge list are added to the node
table automatically. Every measure returns a vector aligned with the rows
of `g$nodes`, so results can be stored as columns and analyzed with
whatever data tools you already use.

## What is here

Centrality (degree, strength, betweenness, closeness, harmonic,
eigenvector, PageRank, Katz, Bonacich power, hubs and authorities, reach,
proximity prestige, information centrality), whole-network measures
(density, reciprocity, transitivity and clustering, distances, components,
cut points and bridges, centralization, core-periphery, small worlds),
cohesion and subgroups (k-cores, cliques, n-cliques, n-clans, k-plexes,
lambda sets, factions, Louvain, Girvan-Newman, leading eigenvector, greedy
modularity, hierarchical clustering), positions and brokerage (structural
holes, Gould-Fernandez brokerage, structural, regular, and automorphic
equivalence, CONCOR, blockmodels), mixing and homophily (assortativity,
mixing matrices, the E-I index), ego networks, two-mode networks, key
players and group centrality, QAP and MRQAP, CUG tests and random graphs,
GeoPackage files,
temporal networks, geospatial networks with [sf](https://r-spatial.github.io/sf/), multidimensional
scaling, layouts and plots, conversions to and from [igraph](https://r.igraph.org/), network, and
[tidygraph](https://tidygraph.data-imaginist.com/), and nine file formats (GraphML, GEXF, GML, [Pajek](http://mrvar.fdv.uni-lj.si/pajek/), [UCINET](https://sites.google.com/site/ucinetsoftware/home) DL,
JSON node-link, DyNetML, CSV, DOT). The
[reference index](https://knapply.github.io/knappnet/reference/index.html)
lists every function by area, and each help page states the conventions
the function follows and the tools it matches.

## One switch for every tool's conventions

[igraph](https://r.igraph.org/), [NetworkX](https://networkx.org/), [sna](https://cran.r-project.org/package=sna), and [UCINET](https://sites.google.com/site/ucinetsoftware/home) often give different answers to the same
question, because they make different choices that their documentation
rarely puts side by side. Closeness on a disconnected network is the
classic case: each tool treats unreachable nodes differently. Here `a`,
`b`, and `c` form a path, and `d` and `e` are a separate pair:

```r
split <- netframe(data.frame(from = c("a", "b", "d"), to = c("b", "c", "e")),
                  directed = FALSE)
data.frame(
  node     = split$nodes$name,
  igraph   = round(closeness(split), 3),
  networkx = round(closeness(split, convention = "networkx"), 3),
  sna      = closeness(split, convention = "sna"),
  ucinet   = round(closeness(split, convention = "ucinet"), 3)
)
#>   node igraph networkx sna ucinet
#> 1    a  0.333    0.333   0  0.444
#> 2    b  0.500    0.500   0  0.500
#> 3    d  1.000    0.250   0  0.400
#> 4    c  0.333    0.333   0  0.444
#> 5    e  1.000    0.250   0  0.400
```

igraph counts only the nodes each node can reach, so `d` and `e` look
perfectly central; NetworkX scales by the share of the network reached
(Wasserman and Faust 1994); sna gives 0 to any node that cannot reach
everyone; UCINET counts each missing distance as one more than the longest
observed one. None is wrong, but a published result should say which was
used, and replicating a study done in another tool should take one
argument, not a reimplementation. `convention` is available on degree,
betweenness, closeness, harmonic and eigenvector centrality, density, and
centralization, and `?conventions` tabulates every difference. Without it,
results follow igraph.

## Validated, not trusted

Every measure is compared with every reference tool that implements it,
and with a slow, deliberately simple implementation of its textbook
definition, on a corpus of classic datasets, random networks across
densities and weightings, and edge cases (empty networks, self-loops,
parallel edges), on Linux, macOS, and Windows, on every change. Thousands
of comparisons pass; each of the few documented differences is explained
in the validation manifest, and the definition oracles show that several
of them come from the reference tools rather than from knappnet.
[validation/REPORT.md](validation/REPORT.md) has the latest results,
[validation/README.md](validation/README.md) explains the framework, and
the output shown in this README is checked against R on every change.

## Benchmarks


Measured on GitHub's standard machines, Linux (AMD EPYC 7763, 2 cores)
and macOS (Apple M1, 3 cores), with R 4.6.1, [igraph](https://r.igraph.org/) 2.3.4, [sna](https://cran.r-project.org/package=sna) 2.8,
Python 3.12, [NetworkX](https://networkx.org/) 3.7, [python-igraph](https://python.igraph.org/) 1.0.0, and [NetworKit](https://networkit.github.io/) 11.2.2.
Times are medians of repeated runs, with networks built beforehand
(building is timed separately). igraph, sna, and NetworkX use one core,
so the one-core column is the like-for-like comparison. Every pair of
timings is also a correctness check: the run fails unless both tools give
the same answer. To reproduce, run `Rscript bench/bench.R` and
`python bench/bench.py`.

From R, on Linux:

| Measure | Network | knappnet, 1 core | knappnet, 2 cores | Other tool | Speedup, 1 core | Speedup, 2 cores |
|---|---|---:|---:|---:|---:|---:|
| Build network + degree | scale-free, 100,000 nodes, 300,000 edges | 55 ms | 55 ms | 282 ms (igraph) | 5.2x | 5.2x |
| Betweenness | random, 10,000 nodes, 50,000 edges | 11.3 s | 7.7 s | 14.7 s (igraph) | 1.3x | 1.9x |
| Betweenness | directed random, 10,000 / 50,000 | 6.6 s | 4.6 s | 9.8 s (igraph) | 1.5x | 2.1x |
| Closeness | random, 10,000 / 50,000 | 2.8 s | 1.9 s | 4.8 s (igraph) | 1.7x | 2.5x |
| Harmonic centrality | random, 10,000 / 50,000 | 2.8 s | 1.9 s | 4.7 s (igraph) | 1.7x | 2.5x |
| Constraint | random, 10,000 / 50,000 | 6.1 ms | 4.5 ms | 12.4 ms (igraph) | 2.0x | 2.7x |
| Triad census | directed random, 10,000 / 50,000 | 37 ms | 27 ms | 55 ms (igraph) | 1.5x | 2.1x |
| k-cores | scale-free, 100,000 / 300,000 | 8.0 ms | 8.0 ms | 11.1 ms (igraph) | 1.4x | 1.4x |
| Maximal cliques | scale-free, 100,000 / 300,000 | 439 ms | 423 ms | 5.5 s (igraph) | 12x | 13x |
| Louvain communities | scale-free, 100,000 / 300,000 | 323 ms | 323 ms | 23.7 s (igraph) | 73x | 73x |
| Eigenvector centrality | scale-free, 100,000 / 300,000 | 537 ms | 540 ms | 172 ms (igraph) | 0.3x | 0.3x |
| Components | scale-free, 100,000 / 300,000 | 1.2 ms | 1.2 ms | 0.4 ms (igraph) | 0.3x | 0.3x |
| Read GraphML | random, 10,000 / 50,000 | 98 ms | 98 ms | 37 ms (igraph) | 0.4x | 0.4x |
| Betweenness | directed random, 2,000 / 10,000 | 236 ms | 153 ms | 1.6 s (sna) | 6.9x | 11x |
| Triad census | directed random, 2,000 / 10,000 | 6.7 ms | 5.1 ms | 239 s (sna) | 35,000x | 46,000x |
| Structural equivalence | directed random, 2,000 / 10,000 | 462 ms | 349 ms | 154 s (sna) | 332x | 440x |
| Brokerage | directed random, 2,000 / 10,000 | 0.7 ms | 0.8 ms | 100 ms (sna) | 138x | 131x |
| QAP correlation | 100 nodes, 1,000 permutations | 37 ms | 34 ms | 1.1 s (sna) | 30x | 33x |
| MRQAP (double semi-partialling) | 100 nodes, 1,000 permutations | 865 ms | 699 ms | 41 s (sna) | 48x | 59x |

On the 3-core Mac, single-core results are similar and the multicore
gains grow: betweenness, closeness, and harmonic centrality run 3.7 to
4.4 times as fast as igraph, and betweenness 155 times as fast as sna.

From Python, on Linux (speedups on 1 core / 2 cores; NetworKit always
uses every core):

| Measure | vs NetworkX | vs python-igraph | vs NetworKit |
|---|---:|---:|---:|
| Betweenness (random, 5,000 / 25,000) | 55x / 78x | 1.2x / 1.7x | 2.5x / 3.6x |
| Closeness | 33x / 49x | | |
| Harmonic centrality | 39x / 58x | 1.4x / 2.2x | 6.0x / 9.0x |
| Triad census (directed) | 96x / 132x | 1.3x / 1.8x | |
| Constraint (2,000 / 10,000) | 16x / 17x | 1.5x / 1.6x | |
| Louvain (scale-free, 100,000 / 300,000) | 723x | 39x | 1.7x |
| Maximal cliques | 6.2x / 6.6x | 5.9x / 6.2x | 0.6x / 0.6x |
| k-cores | 28x | 0.5x | 0.4x |
| Components | 4.7x | 0.3x | 0.3x |
| Eigenvector centrality | 1.6x / 1.7x | 0.8x / 0.9x | 0.3x / 0.3x |

The Python comparisons call the engine through the test harness, which
copies inputs and results between Python lists and the engine; that copy
dominates the fastest, linear-time measures and will be avoided by the
Python package.

Where knappnet is slower, and what is planned:

- **Eigenvector centrality** uses power iteration to a tight tolerance
  (10⁻¹²). igraph uses ARPACK's Lanczos-type eigensolver, and NetworKit
  stops its iteration at a looser tolerance. A Lanczos-type solver is
  planned.
- **Components** takes about a millisecond on 100,000 nodes, but igraph
  is faster still.
- **Reading GraphML** builds a full document tree before converting it; a
  streaming reader is planned.
- **NetworKit**, a parallel C++ library, is faster on k-cores and maximal
  cliques, and close on Louvain.

## How it works

```
R package (r/)            data frames in, vectors and data frames out
   │  .Call
C glue (r/src/*.c)        converts R vectors; no logic of its own
   │  C ABI
crates/capi               a small, stable C interface (netcapi.h)
   │
crates/engine             graph storage, algorithms, file formats (Rust)
```

All computation happens in one engine written in Rust with no third-party
dependencies, published on its own as a Rust crate. Languages reach it only
through the C interface, so the planned Python package will run exactly the
same code, and a result validated once is validated for every language.
Networks are stored in compressed sparse row form; the algorithm for a
path-based measure follows the weights (breadth-first search, Dial's bucket
queue, or Dijkstra's algorithm); per-node work is spread over all cores and
combined in a fixed order, and randomized methods use a counter-based
generator, so results are reproducible for a given seed and thread count.
The [performance article](https://knapply.github.io/knappnet/articles/performance.html)
explains these choices, and [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)
and [CONTRIBUTING.md](CONTRIBUTING.md) describe the layers and how to add
a measure.

## Citing

If you use knappnet in published work, please cite it using the metadata in
[CITATION.cff](CITATION.cff). The methods implemented here come from the
literature listed in [docs/references.bib](docs/references.bib) and
indexed in [docs/references.md](docs/references.md); please cite the
original authors of the methods you use as well.

## License

MIT. See [LICENSE](LICENSE).
