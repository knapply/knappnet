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

> **Status:** early development. The core engine, R bindings, and validation
> framework are in place for degree, strength, betweenness, and the triad
> census. Interfaces may change before the first release.

## Example

```r
library(knappnet)

edges <- data.frame(
  from = c("ana", "ana", "ben", "cy", "cy"),
  to   = c("ben", "cy",  "cy",  "dee", "eli")
)
g <- netframe(edges, directed = FALSE)

g$nodes$degree <- degree(g)
g$nodes$betweenness <- betweenness(g)
g$nodes
```

Results are always aligned with the rows of the node table, so they can be
added as columns directly.

## Bring your existing work

Switching should cost nothing. Networks convert both ways with igraph,
tidygraph, statnet's network, adjacency matrices (dense or sparse), and
plain edge lists, keeping node, edge, and network attributes:

```r
g  <- as_netframe(my_igraph)        # also network, tbl_graph, matrix, Matrix
ig <- to_igraph(g)                   # also to_network(), to_tbl_graph(), as.matrix()
```

Files from other tools read and write directly:

```r
g <- read_network("from_gephi.gexf")
write_network(g, "for_pajek.net")
```

| Format | Read | Write | Common source |
|---|---|---|---|
| GraphML | yes | yes | igraph, NetworkX, Gephi, yEd, Cytoscape |
| GEXF | yes | yes | Gephi |
| GML | yes | yes | igraph, NetworkX, Cytoscape |
| Pajek `.net` | yes | yes | Pajek, igraph, NetworkX |
| UCINET DL | yes | yes | UCINET, igraph |
| JSON node-link | yes | yes | NetworkX, D3 |
| CSV/TSV edge lists | yes | yes | spreadsheets, any tool |
| DOT | no | yes | Graphviz |

Every format is checked in both directions against NetworkX's and igraph's
own readers and writers.

## Installation

Not yet released. Once binaries are published, installation will not require
Rust or any other build tools. To build from source today you need R (4.1 or
newer) and a Rust toolchain:

```r
# install.packages("remotes")
remotes::install_github("knapply/knappnet", subdir = "r")
```

## Validation

The `validation/` directory holds a corpus of fixtures (classic datasets,
random graphs across densities, and edge cases such as empty graphs,
self-loops, and parallel edges), golden results from each reference
implementation, and a manifest describing how every measure maps onto each
tool. `validation/REPORT.md` summarizes the latest comparison. See
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for details.

## Citing

If you use knappnet in published work, please cite it using the metadata in
[CITATION.cff](CITATION.cff). The algorithms implemented here come from the
literature listed in [docs/references.bib](docs/references.bib); please cite
the original authors of the methods you use as well.

## License

MIT. See [LICENSE](LICENSE).
