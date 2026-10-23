# Design decisions

Decisions about the data model and the R API that future work must stay
consistent with. Each entry gives the decision, the alternatives that were
considered, and the reasoning, with sources by BibTeX key
(`docs/references.bib`). The user-facing statement of the data model is
the article "What knappnet means by a network"
(`r/vignettes/articles/definitions.Rmd`); the work the decisions imply is
tracked in [ROADMAP.md](ROADMAP.md).

## October 2026: the data model and the API

These thirteen decisions were taken together in October 2026, on the basis
of a definitions memo (every term's meaning in the literature and in
[igraph](https://r.igraph.org/), network/networkDynamic, [tidygraph](https://tidygraph.data-imaginist.com/), [NetworkX](https://networkx.org/), [graph-tool](https://graph-tool.skewed.de/), [NetworKit](https://networkit.github.io/),
[UCINET](https://sites.google.com/site/ucinetsoftware/home), [Pajek](http://mrvar.fdv.uni-lj.si/pajek/), [ORA](https://www.casos.cs.cmu.edu/projects/ora/), [RSiena](https://www.stats.ox.ac.uk/~snijders/siena/), multinet, pymnet, Raphtory, teneto, and
sfnetworks) and of three surveys of the teaching sources the package is
meant to supersede: `rawlings2023network`, `shafie2026r4sna`,
`wiederstein2021network`, `almeidasilva2024network`, `epskamp2012qgraph`,
`borgatti2022analyzing`, `kadushin2012understanding`, and
`cunningham2016understanding`.

### 1. Every exported R function is prefixed `net_`

**Decision.** Every exported function is named `net_<verb or noun>`, in
the style of stringr's `str_`, with names short enough for interactive
use (standard short forms for long names). No name collides with any
known CRAN or Bioconductor export.

**Alternatives.** Keep the current unprefixed names (`degree()`,
`betweenness()`, `closeness()`, `components()`, `brokerage()`), which are
shared with [sna](https://cran.r-project.org/package=sna) and [igraph](https://r.igraph.org/) so that whichever package is attached last
wins (the [statnet](https://statnet.org/) article has to say so); family prefixes only
(`cluster_*`, `plot_*`, `random_*`, as now) with bare names for measures;
tidygraph-style `centrality_*`, `group_*`, `node_*` families.

**Reasoning.** The surveys report that readers of the teaching sources
lose time to name clashes and to inconsistent naming across packages
(`graph_from_*` versus `make_*` versus `sample_*`, dotted and underscored
names), and that a common prefix makes a package discoverable from the
prompt. A single prefix removes every masking warning at once. The
existing family prefixes survive inside the name (`net_cluster_louvain`,
`net_plot`).

### 2. `print()` is a labelled card plus table previews

**Decision.** Printing a netframe shows a card that states, in order:
size; directed or undirected; one-mode, two-mode, or k-mode with mode
labels and counts; weighted (column, declared meaning, range) or not;
multiplex (the relations) or multilayer; static or temporal, and how (edge
spells, contacts, node spells, waves); spatial (node or edge geometry,
CRS); loops and parallel edges; components, isolates, and density where
they are cheap; followed by the first rows of the node and edge tables.

**Alternatives.** The current one-line sentence ("undirected network with
5 nodes and 5 edges, unweighted"); [igraph](https://r.igraph.org/)'s coded header (`UNW-`);
[tidygraph](https://tidygraph.data-imaginist.com/)'s two-table preview without a summary.

**Reasoning.** All three surveys single out print output: one tutorial
exists only to decode igraph's header, and every survey asks that the
weight meaning, the two-mode, temporal, spatial, and signed status, and
the counts of loops and parallel edges be visible without a helper. The
table preview makes attribute discovery need no accessor function.

### 3. Spatial, temporal, and spatio-temporal networks are S3 subclasses

**Decision.** `c("spatial_netframe", "netframe")`,
`c("temporal_netframe", "netframe")`, and
`c("spatiotemporal_netframe", "spatial_netframe", "temporal_netframe",
"netframe")`.

**Alternatives.** One class with roles only (the current design, where
`is_temporal()` inspects the roles); separate, unrelated classes.

**Reasoning.** Subclasses let `print()`, `plot()`, and `summary()`
dispatch on the kind of network while every static measure keeps working
through inheritance, and the spatio-temporal class inherits both
behaviours. The roles still carry the structure; the class is a label on
top of them.

### 4. Validation in the engine, in four groups

**Decision.** The constructor validates in the Rust engine (so the Python
package applies the same rules), modelled on the checks of
`as.network.data.frame()` in the `network` package (`butts2008network`):
structure (endpoint columns atomic and of one type, no missing
endpoints, every endpoint in the node table with the missing ids listed,
unique node ids with duplicates listed, matching id types); declared
properties (a declared `loops = FALSE`, `multiple = FALSE`, or
between-modes-only rule lists the violating rows; contradictory
declarations are errors, never silent overrides); attribute sanity
(weights numeric and finite, positive when distances; spells ordered and
within node spells; geometry and CRS consistent; mode column complete);
and advisory notes (nodes added from edges, isolates, ids that look
0-based, ids coerced between types, counts of loops and parallel edges
when none were declared), silenced with `quiet = TRUE`.

**Alternatives.** Validation in R only; [igraph](https://r.igraph.org/)'s and [sna](https://cran.r-project.org/package=sna)'s approach of
accepting whatever is given and failing later inside a measure; [UCINET](https://sites.google.com/site/ucinetsoftware/home)'s
approach of transforming silently with a notice.

**Reasoning.** The surveys record that misaligned attribute tables,
endpoints missing from the node table, 0-based ids, and coerced id types
cost readers whole sections of trial and error, and that a message naming
the rows is what they wanted. Doing it once in the engine keeps the two
language bindings identical.

### 5. Weights must say what they mean

**Decision.** `netframe(weight = "w")` is an error unless `weight_type`
is `"strength"` or `"distance"`. There is no default meaning.

**Alternatives.** Default to `"distance"` ([igraph](https://r.igraph.org/)'s and [NetworkX](https://networkx.org/)'s path
functions); default to `"strength"` (what [UCINET](https://sites.google.com/site/ucinetsoftware/home), [Pajek](http://mrvar.fdv.uni-lj.si/pajek/), and every
teaching source assume); a third `"similarity"` type for correlation
networks.

**Reasoning.** `newman2001scientific`, `newman2004analysis`, and
`opsahl2010node` establish that "weight" means two different things and
that path-based measures must invert strengths. Every teaching source
reads weights as strengths while igraph's path functions read them as
distances, and one of them spends a section on closeness computed with
the wrong direction. A wrong default is an invisible error; a required
declaration costs one argument. A similarity type is deferred until a
measure needs it.

### 6. Modes are a node role; bipartite is a graph property

**Decision.** A `mode` node role holds any discrete column, so two-mode
and k-mode networks share one mechanism. [igraph](https://r.igraph.org/)'s logical `type` is still
accepted as input and written on conversion to igraph. An
`is_twomode()`-style function reports the declared modes; an
`is_bipartite()`-style function tests 2-colourability of the structure.

**Alternatives.** Keep the logical `type` role only (igraph's convention,
the current design); `network`'s count of first-mode vertices; [NetworkX](https://networkx.org/)'s
explicit node sets with no declaration; one function for both questions.

**Reasoning.** `wasserman1994social` (ch. 8), `breiger1974duality`, and
`borgatti1997network` define two-mode data by declared node sets;
`newman2018networks` defines bipartite by 2-colourability; no tool keeps
the two apart, and igraph's R `is_bipartite()` checks only for the
attribute. `fararo1984tripartite` motivates k modes. A discrete `mode`
column also removes a collision with the `type` column that the DyNetML
reader writes for [ORA](https://www.casos.cs.cmu.edu/projects/ora/) node classes, and matches ORA's node classes and
the multimode structure of the Noordin data.

### 7. Parallel edges are defined relative to roles

**Decision.** Two edges are parallel when they have the same endpoints
and the same values of the distinguishing roles (relation, time spell).
`a-b friendship` and `a-b advice` are not parallel. Simplification,
simplicity tests, and edge multiplicity follow this rule.

**Alternatives.** Endpoints only ([igraph](https://r.igraph.org/), [NetworkX](https://networkx.org/) `MultiGraph`,
`network`'s `multiple` flag, [Pajek](http://mrvar.fdv.uni-lj.si/pajek/)'s multiple lines), which would make
every multiplex edge table a multigraph; forbidding parallel edges
(NetworkX `Graph`, [NetworKit](https://networkit.github.io/)).

**Reasoning.** An edge list naturally holds repeated pairs, so networks
stay multigraphs by default as in igraph (`newman2018networks` for the
term), but a long multiplex table with one edge per relation should be
simple within each relation, as analysts expect. The rule also keeps
"multiplex" for several relations and "parallel" for repeated pairs,
avoiding `network`'s use of "multiplex edges" for the latter
(`butts2008network`).

### 8. Multiplex networks are a `relation` edge role; multilayer comes later

**Decision.** Several relations live in one long edge table under a
`relation` edge role (the Noordin data already has one). Measures take
`relation =` to select or combine relations; a split helper returns a
list of one-relation networks. A multilayer network in the sense of
`kivela2014multilayer` (node-layers with interlayer edges) is a later,
separate subclass.

**Alternatives.** Stacked matrices ([UCINET](https://sites.google.com/site/ucinetsoftware/home)) or one network per relation
([igraph](https://r.igraph.org/), `network`, two of the teaching sources); a multilayer-only model
with actor-layer vertices (multinet, `magnani2021analysis`; pymnet,
`nurmi2024pymnet`) used for multiplex data too; one row per actor with
`layer_from`/`layer_to` on edges.

**Reasoning.** `kivela2014multilayer` identify a multiplex network with
an edge-coloured multigraph, which is exactly a relation column;
`battiston2014structural` give its measures. The long table needs no new
structure and is what the DyNetML reader already produces. Node-layers
are needed only for interlayer edges between different actors and for
per-layer node attributes, so they wait for multilayer measures; the
rejected `layer_from`/`layer_to` form would give such attributes no home.
The multilevel design of `lazega2008catching` and `lazega2016multilevel`
is covered by `mode` plus `relation` with within-mode edges allowed.

### 9. Panel waves are time spells

**Decision.** Wave k of a panel is the spell `[k, k+1)` in a temporal
netframe. Helpers build the waves from a list of networks and split them
back.

**Alternatives.** A separate `netframe_panel` class holding a list of
netframes with a presence table (the user-facing shape of [RSiena](https://www.stats.ox.ac.uk/~snijders/siena/),
`ripley2026manual`, and networkDynamic's `network.list`,
`butts2026networkdynamic`); a three-dimensional array (teneto,
`thompson2017static`).

**Reasoning.** One engine path then serves panels and event data:
windows, snapshots, and rolling measures apply unchanged. Half-open
spells (`holme2012temporal` for the vocabulary; networkDynamic for the
convention) make wave boundaries unambiguous. The list-of-networks shape
remains the import and export form for RSiena and [ORA](https://www.casos.cs.cmu.edu/projects/ora/)
(`snijders2010introduction`).

### 10. Endpoint columns are named explicitly

**Decision.** `source = "source"` and `target = "target"` are explicit
arguments whose defaults are column names. Column names and positions are
never guessed.

**Alternatives.** Positional defaults (the current `from = 1L, to = 2L`;
[igraph](https://r.igraph.org/)'s `graph_from_data_frame()`, which also makes a numeric third
column the weight); recognising `from`/`to`, `sender`/`receiver`, and
`source`/`target` with a message; `source`/`target` as fixed, unnamed
conventions.

**Reasoning.** The surveys list position-based meaning as a recurring
source of silent errors (renamed Source/Target/Weight columns, a numeric
third column becoming a weight) and ask that nothing be decided by
position. Named defaults keep the common case to zero arguments and make
every other case one explicit argument.

### 11. Signed networks and missing data: adopt the recommendations provisionally

**Decision.** The owner has no preference; the research recommendations
are adopted provisionally and implemented after the core changes above. A
`sign` edge role holds -1 and +1 and weights stay non-negative
(magnitudes); `as_signed()` splits a signed value column. Missing data
get explicit roles: a node role for non-respondents and an edge role for
unobserved ties; measures refuse unresolved missingness unless an
`na_action` is given.

**Alternatives.** Sign inside the weight ([Pajek](http://mrvar.fdv.uni-lj.si/pajek/), [UCINET](https://sites.google.com/site/ucinetsoftware/home), qgraph,
`epskamp2012qgraph`); no sign concept ([igraph](https://r.igraph.org/), [NetworkX](https://networkx.org/)); missing ties as
rows in `edges` versus a separate `missing` table; treating unobserved as
absent (igraph, NetworkX, Pajek, [ORA](https://www.casos.cs.cmu.edu/projects/ora/)).

**Reasoning.** `harary1953notion` and `cartwright1956structural` make
sign categorical; keeping it apart from magnitude means path-based weight
checks are never tripped by a negative tie and a network can be signed
without being weighted. The values -1/+1 are signnet's convention
(`shafie2026r4sna`, ch. 7), so objects convert without recoding.
`kossinets2006effects`, `robins2004missing`, and `handcock2010modeling`
distinguish missing from absent; `network`'s `na` attribute and [RSiena](https://www.stats.ox.ac.uk/~snijders/siena/)'s
codes (`ripley2026manual`) show the statistical tools do too, and the
textbooks' operational notion is the non-respondent (`rawlings2023network`,
ch. 4), so the node role comes first. Whether missing ties live in
`edges` or in a separate table, and the default `na_action`, are still
open.

### 12. Hypergraphs and simplicial complexes are out of scope

**Decision.** Not first-class. The documented exchange is a two-mode
netframe with `mode = c("node", "hyperedge")`; the clique expansion is
the ordinary two-mode projection.

**Alternatives.** A third table of hyperedges; `network`'s `hyper = TRUE`
edges with several heads and tails.

**Reasoning.** `battiston2020networks` and `bianconi2021higher` define
the objects; `torres2021why` explain why a hypergraph is not a two-mode
network even though the incidence structure is the same. An edge with two
endpoints is the whole data model; XGI (`landry2023xgi`) exists for
higher-order data.

### 13. ERGM, TERGM, and SAOM stay on the `ergm-saom` branch

**Decision.** The statistical models live on the `ergm-saom` branch until
their API is redesigned; they are not on `main`. The documented hand-off
is `to_network()` for [ergm](https://cran.r-project.org/package=ergm) and tergm and the planned [RSiena](https://www.stats.ox.ac.uk/~snijders/siena/) export.

**Alternatives.** Merge the current implementations; drop them and bridge
only.

**Reasoning.** The surveys show that the modelling chapters' difficulty is
data preparation, which the bridges solve, and that the models deserve an
API designed with the rest of the package rather than one carried over.
