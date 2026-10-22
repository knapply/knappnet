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
networks, and for statistical models such as exponential random graph models
and stochastic actor-oriented models.

> **Status:** early development. Available now: degree, strength,
> betweenness, closeness, harmonic and eigenvector centrality,
> centralization, density, components, fragmentation, the triad census,
> k-cores, cliques, modularity and Louvain
> communities, Burt's structural holes, Gould-Fernandez brokerage,
> structural equivalence, CONCOR, QAP correlation and MRQAP regression,
> conversions to and from igraph, network, tidygraph, and matrices, and nine
> file formats. Interfaces may change before the first release.

**Contents:** [A first network](#a-first-network) ·
[Weighted and binary networks](#weighted-and-binary-networks) ·
[Bring your existing work](#bring-your-existing-work) ·
[Conventions](#one-switch-for-every-tools-conventions) ·
[Cohesive subgroups](#cohesive-subgroups) ·
[Brokerage and structural holes](#brokerage-and-structural-holes) ·
[Positions and roles](#positions-and-roles) ·
[Testing hypotheses with QAP](#testing-hypotheses-with-qap) ·
[Files](#files) · [How it works](#how-it-works) ·
[Installation](#installation) · [Citing](#citing)

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

Most social data is directed. Here, `advice` records who goes to whom for
advice in a small office, with each person's department in the node table:

```r
staff <- data.frame(
  name = c("ann", "bob", "cat", "dan", "eve", "fay"),
  dept = c("sales", "sales", "eng", "eng", "eng", "legal")
)
advice <- netframe(
  data.frame(from = c("ann", "bob", "cat", "dan", "eve", "bob", "fay", "ann"),
             to   = c("bob", "cat", "dan", "eve", "cat", "fay", "dan", "cat")),
  nodes = staff
)
advice
#> <netframe> directed network with 6 nodes and 8 edges, unweighted
degree(advice, mode = "in")    # how often each person is asked
#> [1] 0 1 3 2 1 1
betweenness(advice)
#> [1] 0 1 4 5 2 1
```

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

The same data read the wrong way gives a different answer, which is why the
meaning is stated once on the network instead of being guessed by each
function. Measures defined only for binary ties (cliques, brokerage, the
triad census) ignore weights and say so in their help pages.

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

## One switch for every tool's conventions

igraph, NetworkX, sna, and UCINET often give different answers to the same
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

The network as a whole has measures too:

```r
round(graph_density(karate), 3)
#> [1] 0.139
components(split)$no
#> [1] 2
round(centralization(karate, "betweenness"), 3)
#> [1] 0.406
```

## Cohesive subgroups

Zachary's (1977) karate club split into two factions after a dispute
between the instructor (node 1) and the club's administrator (node 34).
The k-core decomposition (Seidman 1983) shows its dense center; maximal
cliques (Bron and Kerbosch 1973) list its tightest groups:

```r
table(coreness(karate))
#>
#>  1  2  3  4 
#>  1 11 12 10 
clique_number(karate)
#> [1] 5
max_cliques(karate, min_size = 5)
#> [[1]]
#> [1] 1 2 3 4 8
#>
#> [[2]]
#> [1]  1  2  3  4 14
#>
```

Louvain community detection (Blondel et al. 2008) maximizes modularity
(Newman and Girvan 2004). It visits nodes in a random order, so it takes a
seed; with the same seed it gives the same answer on every computer and
with any number of threads:

```r
groups <- cluster_louvain(karate, seed = 42)
table(groups)
#> groups
#>  1  2  3  4 
#> 12  5 13  4 
round(attr(groups, "modularity"), 3)
#> [1] 0.416
round(modularity(karate, groups), 3)  # any partition can be scored
#> [1] 0.416
```

## Brokerage and structural holes

Burt's (1992) constraint is low for people whose contacts are not connected
to each other. In the karate club, the two leaders span the most structural
holes:

```r
holes <- data.frame(node = karate$nodes$name, constraint = round(constraint(karate), 3))
head(holes[order(holes$constraint), ], 5)
#>    node constraint
#> 1     1      0.155
#> 34   34      0.156
#> 3     3      0.185
#> 32   32      0.245
#> 33   33      0.249
```

`effective_size()` and `efficiency()` complete Burt's measures. In the
first network, `cy` connects otherwise unconnected people:

```r
round(constraint(g), 3)
#> [1] 0.953 0.953 0.406 1.000 1.000
effective_size(g)
#> [1] 1.0 1.0 3.5 1.0 1.0
```

Gould and Fernandez (1989) brokerage counts, for each person, the pairs
they connect who are not tied directly, by the groups involved. A
gatekeeper controls access into their own group; a liaison connects two
groups neither of which is their own:

```r
brokerage(advice, "dept")
#>   node coordinator itinerant representative gatekeeper liaison total
#> 1  ann           0         0              0          0       0     0
#> 2  bob           0         0              1          0       0     1
#> 3  cat           1         0              0          2       0     3
#> 4  dan           1         0              0          1       0     2
#> 5  eve           1         0              0          0       0     1
#> 6  fay           0         0              0          0       1     1
```

## Positions and roles

Two people are structurally equivalent when they send ties to, and receive
ties from, the same others (Lorrain and White 1971). The Hamming distance
counts the ties on which two people differ:

```r
structural_equivalence(advice, "hamming")
#>     ann bob cat dan eve fay
#> ann   0   1   4   5   2   4
#> bob   1   0   3   6   3   3
#> cat   4   3   0   5   4   2
#> dan   5   6   5   0   3   3
#> eve   2   3   4   3   0   4
#> fay   4   3   2   3   4   0
```

CONCOR (Breiger, Boorman, and Arabie 1975) splits the network into
positions by repeatedly correlating these profiles. Each extra level of
`depth` splits every position again. CONCOR can be numerically fragile, so
knappnet warns when rounding error may have decided a split:

```r
concor(advice)
#> [1] 1 1 2 2 1 2
concor(advice, depth = 2)
#> [1] 1 1 3 4 2 3
```

## Testing hypotheses with QAP

Ties in a network are not independent observations, so ordinary standard
errors don't apply. The quadratic assignment procedure (Krackhardt 1987,
1988) builds the null distribution by permuting node labels instead. Here
we ask whether advice follows friendship or department lines:

```r
friends <- netframe(
  data.frame(from = c("ann", "ann", "bob", "cat", "dan"),
             to   = c("bob", "eve", "eve", "fay", "eve")),
  nodes = staff, directed = FALSE
)
same_dept <- outer(staff$dept, staff$dept, "==") * 1

test <- qap_cor(advice, friends, exact = TRUE)
round(test$statistic, 3)
#> [1] -0.107
round(test$p.value, 3)
#>   Pr(<=r)   Pr(>=r) Pr(>=|r|) 
#>     0.422     0.883     0.611 
```

Multiple regression QAP uses the double semi-partialling method of Dekker,
Krackhardt, and Snijders (2007), which stays accurate when predictors are
correlated:

```r
fit <- mrqap(advice, list(friends = friends, same_dept = same_dept), exact = TRUE)
round(fit$coefficients, 3)
#> (Intercept)     friends   same_dept 
#>       0.229      -0.171       0.357 
round(fit$p.value[, "Pr(>=|b|)"], 3)
#> (Intercept)     friends   same_dept 
#>       0.100       0.214       0.031 
```

With six people there are only 6! = 720 permutations, so `exact = TRUE`
evaluates all of them and the p-values are exact. Larger networks use
random permutations (`reps = 1000` by default) drawn from a seeded
generator; `print()` and `summary()` give the full tables. Matrices,
netframes, and a mix of both are accepted; inputs with node names are
aligned by name.

## Files

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

ORA's DyNetML format holds multimode meta-networks: several classes of
nodes (agents, organizations, tasks, resources) and several networks among
them. Node classes come from the `class` and `type` columns, and networks
from an edge column called `network`:

```r
org <- netframe(
  data.frame(from    = c("ann", "bob", "ann", "bob", "cat"),
             to      = c("bob", "cat", "sales", "sales", "eng"),
             network = c("advice", "advice", "member", "member", "member")),
  nodes = data.frame(name  = c("ann", "bob", "cat", "sales", "eng"),
                     class = c("people", "people", "people", "units", "units"),
                     type  = c("Agent", "Agent", "Agent", "Organization", "Organization"))
)
path <- tempfile(fileext = ".xml")
write_network(org, path, format = "dynetml")
back <- read_network(path)
back$nodes
#>      id  class         type
#> 1   ann people        Agent
#> 2   bob people        Agent
#> 3   cat people        Agent
#> 4 sales  units Organization
#> 5   eng  units Organization
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

## How it works

### Layers

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
The R package carries its own copy of the engine source and builds offline,
as CRAN requires.

### From data frames to compressed sparse rows

When a measure is first computed, the R package passes the engine just
what it needs: each edge's endpoints as integer positions in the node
table and, if relevant, the weights. The engine stores the network in
compressed sparse row (CSR) form: every node's neighbors sit next to each
other in one array, sorted, with a second array marking where each node's
run begins (directed networks also keep the reverse, CSC, for in-ties).
Walking a node's neighbors is a sequential read of memory, which is what
modern processors do fastest, and node positions are 32-bit integers, so
the structure stays compact. Attributes never leave R; they stay
in your data frames.

### Choosing the algorithm for the weights

Whether a network is weighted, and what its weights mean, decides the
algorithm, so the fastest exact method is always used:

| Network | Shortest paths | Cost per source |
|---|---|---|
| Binary, or all weights equal | breadth-first search | O(m) |
| Whole-number weights up to 1024 | bucket queue (Dial 1969) | O(m + D) |
| Any other positive weights | Dijkstra (1959), binary heap | O(m log n) |

Here m is the number of edges, n the number of nodes, and D the longest
shortest-path length. Whole-number lengths are added in exact integer
arithmetic. Betweenness then follows Brandes (2001, 2008). With real-valued weights,
two paths of the same length can sum to slightly different floating-point
values, which silently drops shortest paths in some tools; the engine
treats lengths within a relative 10⁻¹⁰ as equal, and the validation suite
contains cases built to catch this. Dense binary networks can also be held
as bit matrices, where testing a tie or intersecting two neighborhoods
handles 64 nodes per machine word.

### Parallel and reproducible

Most measures repeat independent work for every node (one shortest-path
search per source, one neighborhood per node), which the engine spreads
over all cores; `threads` limits it. Partial results are combined in a
fixed order, so the output is identical with 1 thread or 64. Randomized
methods (Louvain, QAP permutations) use a counter-based random number
generator (Philox; Salmon et al. 2011): every task gets its own stream,
computed from the seed and the task's number, so results do not depend on
which thread ran what. When you don't give a seed, one is drawn from R's
generator, so `set.seed()` works as usual.

### Conventions, stated rather than assumed

Tools disagree on details that change results: whether undirected
betweenness counts each pair once or twice, whether self-loops add to
degree, how isolates enter normalization, whether parallel edges count
separately. Each help page states the choices a measure makes and which
tools it matches, `convention` reproduces each tool's choices exactly, and
the validation manifest records every known disagreement with the reason
for it. UCINET's conventions are taken from its documentation and checked
against its published output.

### Validated, not trusted

The `validation/` directory holds a corpus of fixtures (classic datasets,
random graphs across densities and weightings, and edge cases such as empty
graphs, self-loops, and parallel edges) and golden results from three kinds
of oracle:

- **Reference tools:** igraph and sna in R, and NetworkX in Python.
- **Definitions:** slow, deliberately simple implementations written from
  the textbook formula and sharing no code with the engine, using exact
  arithmetic where possible.
- **Exact enumeration:** for QAP, all permutations of small networks, so
  p-values are checked exactly rather than within simulation error.

Every measure is compared with each oracle that implements it, on every
fixture, on Linux, macOS, and Windows, on every change. Over a thousand comparisons currently
pass, and each of the few documented differences is explained; the
definition oracles show that several of them come from reference tools
rather than from knappnet. `validation/REPORT.md` has the latest results,
and the output shown in this README is regenerated and compared with R on
every change by `tools/check_readme.R`. See
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for details.

### Room to grow

The engine keeps algorithms separate from storage, so new network types
(multilayer, temporal, spatial) are new column roles on the same two
tables rather than new object classes. Inner loops avoid
pointer-chasing data structures so that they can later gain GPU backends,
not tied to one vendor, checked against the same oracles.

## Installation

Not yet released. Once binaries are published, installation will not require
Rust or any other build tools. To build from source today you need R (4.1 or
newer) and a Rust toolchain:

<!-- norun -->
```r
# install.packages("remotes")
remotes::install_github("knapply/knappnet", subdir = "r")
```

## Citing

If you use knappnet in published work, please cite it using the metadata in
[CITATION.cff](CITATION.cff). The methods implemented here come from the
literature listed in [docs/references.bib](docs/references.bib) and
summarized in [docs/references.md](docs/references.md); please cite the
original authors of the methods you use as well.

## License

MIT. See [LICENSE](LICENSE).
