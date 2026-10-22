# knappnet

Fast, validated network analysis for R, with Python to follow.

knappnet is a network analysis engine built around three commitments:

1. **Rigorously accurate.** Every measure is checked against established
   reference implementations (igraph, sna, NetworkX) and against a direct
   evaluation of its mathematical definition. Where those tools follow
   different conventions, the differences are documented and selectable,
   never silent.
2. **Fast.** Algorithms run in a compiled, multithreaded core with a
   cache-friendly compressed sparse layout. Speed never comes at the cost of
   accuracy: approximate methods are always opt-in and clearly labeled.
3. **Easy.** A network is just two data frames, one for nodes and one for
   edges, with any attributes you like. Installation never requires a
   compiler for end users.

The scope covers social network analysis first, then networks of every kind:
computer, transportation, infrastructure, biological, and more. The design
already accounts for multilayer, multiplex, multimodal, temporal, and spatial
networks, and for statistical models such as QAP, exponential random graph
models, and stochastic actor-oriented models.

> **Status:** early development. The engine, R package, file formats,
> conversions, and validation framework are in place for degree, strength,
> betweenness, and the triad census. Interfaces may change before the first
> release.

## Example

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

Results are always aligned with the rows of the node table, so they can be
added as columns directly.

## Weighted and binary networks

A network without a weight column is binary: an edge exists or it doesn't,
and every algorithm uses its fastest unweighted form. With weights, you say
once what they mean. Distances (travel time, cost, latency) keep nodes
apart:

```r
roads <- data.frame(
  from    = c("a", "b", "a"),
  to      = c("b", "c", "c"),
  minutes = c(5, 5, 20)
)
trips <- netframe(roads, directed = FALSE, weight = "minutes")
trips
#> <netframe> undirected network with 3 nodes and 3 edges, weighted by `minutes` (distances)
betweenness(trips)                    # a to c is quicker through b
#> [1] 0 1 0
betweenness(trips, weighted = FALSE)  # ignoring travel time
#> [1] 0 0 0
```

Strengths (contact frequency, tie intensity) bring nodes together, so
path-based measures use 1/weight as the length, while measures that add up
weights use them as given:

```r
calls <- data.frame(
  from  = c("a", "b", "a"),
  to    = c("b", "c", "c"),
  calls = c(1, 1, 20)
)
phone <- netframe(calls, directed = FALSE, weight = "calls", weight_type = "strength")
betweenness(phone)                    # a and c talk directly, so b brokers nothing
#> [1] 0 0 0
strength(phone)
#> [1] 21  2 21
```

The engine picks the fastest exact algorithm for the weights it is given:
breadth-first search for binary networks (and for weights that are all
equal), an exact bucket queue for small whole-number weights, and Dijkstra's
algorithm otherwise.

## Bring your existing work

Switching should cost nothing. Networks convert both ways with igraph,
tidygraph, statnet's network, adjacency matrices (dense or sparse), and
plain edge lists, keeping node, edge, and network attributes:

```r
karate <- as_netframe(igraph::make_graph("Zachary"))
karate
#> <netframe> undirected network with 34 nodes and 78 edges, unweighted
round(head(betweenness(karate)), 1)
#> [1] 231.1  28.5  75.9   6.3   0.3  15.8

ig <- to_igraph(karate)               # also to_network(), to_tbl_graph(), as.matrix()
igraph::ecount(ig)
#> [1] 78
```

Files from other tools read and write directly:

```r
path <- tempfile(fileext = ".graphml")
write_network(trips, path)
read_network(path, weight = "minutes")$edges
#>   from to minutes
#> 1    a  b       5
#> 2    b  c       5
#> 3    a  c      20
```

| Format | Read | Write | Common source |
|---|---|---|---|
| GraphML | yes | yes | igraph, NetworkX, Gephi, yEd, Cytoscape |
| GEXF | yes | yes | Gephi |
| GML | yes | yes | igraph, NetworkX, Cytoscape |
| Pajek `.net` | yes | yes | Pajek, igraph, NetworkX |
| UCINET DL | yes | yes | UCINET, igraph |
| JSON node-link | yes | yes | NetworkX, D3 |
| DyNetML (`.xml`) | yes | yes | ORA (multimode meta-networks) |
| CSV/TSV edge lists | yes | yes | spreadsheets, any tool |
| DOT | no | yes | Graphviz |

Every format is checked in both directions against NetworkX's and igraph's
own readers and writers, except DyNetML, which neither reads; it is tested
against its specification and ORA's file layout.

## Installation

Not yet released. Once binaries are published, installation will not require
Rust or any other build tools. To build from source today you need R (4.1 or
newer) and a Rust toolchain:

<!-- norun -->
```r
# install.packages("remotes")
remotes::install_github("knapply/knappnet", subdir = "r")
```

## Validation

The `validation/` directory holds a corpus of fixtures (classic datasets,
random graphs across densities and weightings, and edge cases such as empty
graphs, self-loops, and parallel edges), golden results from each reference
implementation, and a manifest describing how every measure maps onto each
tool. `validation/REPORT.md` summarizes the latest comparison. The output
shown in this README is checked against R on every change by
`tools/check_readme.R`. See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for
details.

## Citing

If you use knappnet in published work, please cite it using the metadata in
[CITATION.cff](CITATION.cff). The algorithms implemented here come from the
literature listed in [docs/references.bib](docs/references.bib); please cite
the original authors of the methods you use as well.

## License

MIT. See [LICENSE](LICENSE).
