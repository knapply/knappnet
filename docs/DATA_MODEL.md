# The data model

The single specification of every data structure the package supports.
The R and Python bindings, the engine's `frame` module
(`crates/engine/src/frame/`), the C ABI section "Data model and
validation" of `crates/capi/include/netcapi.h`, and the documentation
follow it. Where this page and an older page disagree, this page wins;
the decisions it implements are in [DECISIONS.md](DECISIONS.md) (items
1-19), the vocabulary in the article "What knappnet means by a network"
(`r/vignettes/articles/definitions.Rmd`). Sources are cited by BibTeX key
(`docs/references.bib`).

Contents: [the two tables](#1-the-two-tables), then one section per
configuration ([directed](#2-direction), [two-mode and
k-mode](#3-two-mode-and-k-mode-networks), [weighted](#4-weighted-networks),
[signed](#5-signed-networks), [multiplex](#6-multiplex-networks),
[multilayer](#7-multilayer-networks), [temporal](#8-temporal-networks),
[node spells and node changes](#9-node-spells-and-node-changes),
[spatial](#10-spatial-networks), [spatio-temporal](#11-spatio-temporal-networks),
[ego](#12-ego-networks), [missing data](#13-missing-data),
[loops and parallel edges](#14-loops-and-parallel-edges),
[network-level attributes](#15-network-level-attributes),
[hypergraphs](#16-hypergraphs-and-simplicial-complexes-out-of-scope)),
then [how configurations combine](#17-combining-configurations), the
[print card](#18-the-print-card), and the catalogues of
[roles](#19-role-catalogue), [classes](#20-class-catalogue), and
[diagnostics](#21-diagnostic-catalogue), and the [open
questions](#22-open-questions).

## 1. The two tables

A network is a pair of tables (`wasserman1994social`, chapters 2 and 3,
separate the structural variables, the ties, from the composition
variables, the attributes; the edge list with a vertex table is the input
form of `butts2008network` and of igraph's `graph_from_data_frame()`,
`csardi2006igraph`):

| Table | One row per | Required roles | Optional roles |
|---|---|---|---|
| `nodes` | node | `id` | `mode`, `node_start`, `node_end`, `x`, `y`, `node_missing`, `ego` |
| `edges` | edge | `source`, `target` | `weight`, `relation`, `sign`, `start`, `end`, `wave`, `edge_missing`, `source_layer`, `target_layer`, `ego` |

Three optional long tables, attached to the network rather than to a row
of it, hold what one row per node or per edge cannot (decisions 15 and
19):

| Table | One row per | Required roles | Optional |
|---|---|---|---|
| `layers` | layer | `layer` (the layer identifier) | any layer attribute |
| `node_spells` | presence spell of a node | `id`, `start`, `end` | nothing |
| `node_changes` | change of a node attribute | `id`, `time` | the attribute columns whose values take effect at `time` |

A **role** names the column that carries a structural meaning; every
other column is an attribute and stays with the binding. Roles are set
by name in the constructor (`net_frame(edges, nodes, source = "source",
target = "target", id = "id", ...)`); nothing is read from column
positions (decision 10). The engine receives only the role columns, as
typed columns of one of four kinds (integer, double, text, logical, each
with missing values), validates them, and builds the graph; attributes
never cross into the engine.

**Identifiers.** `id`, `source`, and `target` are integer, double with
whole values, or text. The two endpoint columns share one type, and the
identifier column matches it; an integer column and a double column are
matched as numbers with a note (`N_IDS_COERCED`). Missing identifiers and
empty strings are errors. Node identifiers are unique. An endpoint absent
from the node table is an error (`E_UNKNOWN_NODE`) unless the frame
declares `add_nodes`, or has no node table, in which case the node is
appended in order of first appearance and reported (`N_NODES_ADDED`); an
added node has no attributes, so a frame with a `mode` or `ego` role
cannot add nodes.

**Declarations.** Alongside the roles, the constructor takes declared
properties: `directed` (required), and optionally `loops`, `multiple`,
`between_modes_only`, `egos_are_nodes`, `couplings_only`,
`overlapping_spells`, each `TRUE`, `FALSE`, or undeclared, and
`add_nodes`. A declaration is enforced: the rows that violate it are
errors (as in `as.network.data.frame()`, `butts2008network`), and two
declarations that contradict each other are an error
(`E_DECLARATION_CONFLICT`), never a silent override. An undeclared
property is observed and reported as a note (decision 4).

**Validation.** Every check runs in the engine in O(n + m) time (plus
O(s log s) to sort s node spells) and returns every error and note it can
establish, each with a stable code (section 21), the table it concerns,
a message, the first 50 offending 0-based rows, and the total count. The
result is the same for every thread count. The checks fall into four
groups: structure, declared properties, attribute sanity, and advisory
notes. A frame with an error cannot be built; notes are shown unless
`quiet = TRUE`.

**Construction.** Node `i` of the built graph is row `i` of the node
table, followed by the added nodes in order of first appearance. The
graph's edges are the **observed** rows of the edge table in row order
(section 13); when no row is flagged missing, edge `e` is row `e`.
The engine returns both mappings so that every result aligns with the
tables.

**Class.** `net_frame` (decision 16). Every configuration below says
which class it adds; section 20 gives the fixed order.

## 2. Direction

**Definition.** Direction is a property of the network: every edge is an
ordered pair (directed) or an unordered pair (undirected)
(`wasserman1994social`, chapter 4; `newman2018networks`, chapter 6).
There are no mixed networks.

**Representation.** The `directed` declaration. In an undirected network
the endpoint columns are kept as given and never swapped; `(a, b)` and
`(b, a)` are the same edge for every check and measure.

**Validation.** None beyond the structure checks. Print: `directed: yes`
or `no`.

## 3. Two-mode and k-mode networks

**Definition.** A mode is a declared, disjoint set of nodes of one kind.
A two-mode (affiliation) network has two modes and every edge joins
them (`wasserman1994social`, chapter 8; `breiger1974duality`;
`borgatti1997network`); a k-mode (k-partite) network has k modes with the
same rule (`fararo1984tripartite`); a multimode network declares modes
but allows ties within a mode (the multilevel design of
`lazega2008catching` and `lazega2016multilevel`, and an ORA meta-network,
`tsvetovat2004dynetml`). Being bipartite as a graph (2-colourable,
`newman2018networks`) is a structural property that needs no declaration
and is tested by a separate function, not by validation (decision 6).

**Representation.**

| Table | Role | Type | Constraints |
|---|---|---|---|
| nodes | `mode` | text, integer, or logical | no missing values; every node has one |

Declaration `between_modes_only`: `TRUE` makes every edge join two
different modes; `FALSE` allows within-mode edges (multimode);
undeclared observes. igraph's logical `type` column is accepted as a
mode column as is and written back on conversion to igraph. Mode labels
are numbered by first appearance in the node table; the first mode is
the row mode of the incidence matrix.

**Validation.** `E_MODE_TYPE` (double column), `E_MODE_MISSING` (missing
values, or nodes added from the edges), `E_MODE_UNDECLARED`
(`between_modes_only` without a mode role), `E_MODE_WITHIN` (declared
between-modes-only, but the listed edges join one mode),
`E_DECLARATION_CONFLICT` (`between_modes_only = TRUE` with
`loops = TRUE`: a loop stays within a mode). With layers (section 7),
coupling edges are exempt from the between-modes rule, since they join
one actor to itself.

**Class.** None; `is_twomode()`-style functions read the declaration.

**Print.** `mode: two-mode (18 women, 14 events; between modes only)`,
`mode: 3 modes (agent 20, organisation 5, location 4; within-mode ties
allowed)`, or `mode: one-mode (34 actors)`.

**Example.** A two-mode network of people and events.

nodes:

| id | mode |
|---|---|
| ann | person |
| bob | person |
| e1 | event |
| e2 | event |

edges (`between_modes_only = TRUE`):

| source | target |
|---|---|
| ann | e1 |
| ann | e2 |
| bob | e2 |

## 4. Weighted networks

**Definition.** A network is weighted when one numeric edge column is
declared the weight; otherwise it is binary. The weight's meaning is
declared once and has no default (decision 5): `strength` (frequency,
intensity, capacity; larger brings nodes closer; `barrat2004architecture`)
or `distance` (cost, length, latency; larger keeps them apart;
`newman2001scientific`, `opsahl2010node`). Path-based measures use the
weight as a length under `distance` and `1/w` under `strength`; additive
measures use weights as given.

**Representation.**

| Table | Role | Type | Constraints |
|---|---|---|---|
| edges | `weight` + `weight_type` | integer or double | finite, not missing; positive when distances; not negative when the network is signed |

**Validation.** `E_WEIGHT_TYPE`, `E_WEIGHT_TYPE_MISSING` (no declared
meaning), `E_WEIGHT_NOT_FINITE`, `E_DISTANCE_NOT_POSITIVE`,
`E_WEIGHT_NEGATIVE` (signed network).

**Class.** None. **Print.** `weights: weight (strength, 1-7)`.

**Example.** edges, `weight = "years", weight_type = "strength"`:

| source | target | years |
|---|---|---|
| ann | bob | 4 |
| ann | cy | 2 |
| bob | cy | 7 |

## 5. Signed networks

**Definition.** Every edge is positive or negative, optionally with a
magnitude (`harary1953notion`, `cartwright1956structural`). Sign is
categorical and separate from weight, so a negative tie is never read as
a negative length (decision 11).

**Representation.**

| Table | Role | Type | Constraints |
|---|---|---|---|
| edges | `sign` | integer `-1`/`1`, or logical (`TRUE` positive) | not missing |

The values `-1`/`1` are signnet's convention (`shafie2026r4sna`, chapter
7). A `weight`, if present, is the magnitude and must not be negative. A
helper (`as_signed()`) splits a signed value column into `sign` and
`weight`, which is how Pajek and UCINET data arrive
(`denooy2018exploratory`, `borgatti2002ucinet`).

**Validation.** `E_SIGN_TYPE`, `E_SIGN_VALUE`, `E_WEIGHT_NEGATIVE`.

**Class.** None. **Print.** `signs: sign (12 positive, 3 negative)`.

**Example.** edges, `sign = "sign"`:

| source | target | sign |
|---|---|---|
| ann | bob | 1 |
| ann | cy | -1 |
| bob | cy | 1 |

## 6. Multiplex networks

**Definition.** Several relations on one node set, each edge in exactly
one relation (`wasserman1994social`, chapter 3); formally an
edge-coloured multigraph and the special case of a multilayer network
whose only interlayer edges are couplings (`kivela2014multilayer`;
`battiston2014structural` give the measures). Decision 8.

**Representation.**

| Table | Role | Type | Constraints |
|---|---|---|---|
| edges | `relation` | text, integer, or logical | not missing |

One long edge table; relation labels are numbered by first appearance.
Measures take `relation =` to select or combine; a split helper returns
one network per relation. Two edges of different relations are never
parallel (section 14).

**Validation.** `E_RELATION_TYPE`, `E_RELATION_MISSING`.

**Class.** None. **Print.** `relations: friendship (42), advice (36)
(multiplex)`.

**Example.** edges, `relation = "relation"`:

| source | target | relation |
|---|---|---|
| ann | bob | friendship |
| ann | bob | advice |
| bob | cy | advice |

## 7. Multilayer networks

**Definition.** Following `kivela2014multilayer` and
`dedomenico2013mathematical`: a set of actors, a set of layers,
node-layers (an actor's presence in a layer), intralayer edges within one
layer, and interlayer edges between layers. An interlayer edge joining an
actor to itself in another layer is a **coupling**; one joining two
different actors across layers is a plain interlayer edge. The network is
node-aligned when every actor is in every layer; this model is always
node-aligned (section 22). Mucha et al. call layers "slices"
(`mucha2010community`); multinet (`magnani2021analysis`) and pymnet
(`nurmi2024pymnet`) use the same actor-by-layer vocabulary. Decision 15.

**Representation.** Actors stay in the node table, one row each; every
edge names the layer of each endpoint.

| Table | Role | Type | Constraints |
|---|---|---|---|
| edges | `source_layer` | text, integer, or logical | set together with `target_layer`; not missing |
| edges | `target_layer` | same type as `source_layer` | not missing |
| layers (optional) | `layer` | text, integer, or logical | unique, not missing; every edge layer is one of them |

Classification of an edge with endpoints `(a, b)` and layers `(la, lb)`:

| `la == lb` | `a == b` | Class |
|---|---|---|
| yes | no | intralayer |
| yes | yes | intralayer self-loop |
| no | yes | coupling |
| no | no | interlayer |

Layer labels are the layer table's identifiers in row order when the
table is given (unused layers keep a count of 0), otherwise the layer
values of the edge table in order of first appearance (source before
target within a row). Layer attributes live in the layer table and stay
with the binding; the engine receives only its identifier column.

Declaration `couplings_only`: `TRUE` makes every interlayer edge a
coupling (the multiplex special case); `FALSE` allows any; undeclared
observes. A multiplex network stored with a `relation` role needs no
layer columns; it is lifted to layer form by a helper that copies the
relation into both layer columns and adds the couplings asked for
(categorical: all pairs of layers; ordinal: adjacent layers). The
`relation` and layer roles may coexist (layers of time windows, relations
of tie type); the parallel-edge key then holds both.

**Validation.** `E_LAYER_ROLE_PAIR` (one layer column without the
other), `E_LAYER_TYPE` (double), `E_LAYER_MISSING` (missing layer
values), `E_LAYER_UNKNOWN` (a layer absent from the layer table),
`E_TABLE_ROLE_UNDECLARED` (a layer table without its `layer` role),
`E_LAYER_ID_MISSING`, `E_LAYER_ID_DUPLICATE`, `E_INTERLAYER_DECLARED`
(declared `couplings_only`, but the listed interlayer edges join
different actors). Loops and parallel edges are layer-aware (section 14):
a coupling is not a loop, and the parallel-edge key includes both
layers.

**Class.** Adds `multilayer_net_frame`.

**Print.** `layers: work, family, online (3); 120 intralayer, 60
coupling, 2 interlayer` and, when a layer table is given, `(layer
table: 3 x 2)`.

**Example.** Two actors in layers `work` and `online`, with one
coupling.

nodes:

| id |
|---|
| ann |
| bob |

edges (`source_layer = "from_layer", target_layer = "to_layer"`):

| source | target | from_layer | to_layer |
|---|---|---|---|
| ann | bob | work | work |
| ann | bob | online | online |
| ann | ann | work | online |

layers (`layer = "layer"`):

| layer | medium |
|---|---|
| work | office |
| online | chat |

Row 3 is a coupling, not a self-loop; rows 1 and 2 are not parallel
because their layers differ. Subgraphs per layer (all actors, the
intralayer edges of one layer) are built by the engine on request and
returned with their node and edge row mappings.

## 8. Temporal networks

**Definition.** Edges carry times, as instants or intervals, so the
structure is defined at every time (`holme2012temporal`,
`masuda2016guide`). In the vocabulary of `latapy2018stream`, edges with a
start only form a contact sequence, edges with a start and an end a link
stream, and a link stream whose nodes also have presence spells a stream
graph (the "interval graph" of `holme2012temporal` is avoided because
graph theory uses the term for the intersection graph of
`lekkerkerker1962representation`). A panel (longitudinal) network is the
same node set observed at discrete waves (`snijders2010introduction`).

**Representation.** Three forms, one engine path.

| Form | Table | Roles | Meaning |
|---|---|---|---|
| contacts | edges | `start` | an instantaneous contact at `start` |
| edge spells | edges | `start`, `end` | the half-open spell `[start, end)`; `end == start` is a contact |
| waves | edges | `wave` | integer wave `k` is the spell `[k, k + 1)` (decision 9) |

| Table | Role | Type | Constraints |
|---|---|---|---|
| edges | `start` | integer or double | not missing, finite |
| edges | `end` | integer or double | needs `start`; not missing, finite; `end >= start` |
| edges | `wave` | integer, or double with whole values | not missing; excludes `start` and `end` |

Times are doubles in the binding's unit (seconds since an epoch, days,
wave numbers); the binding records the class and time zone to print and
restore. Half-open spells are the convention of networkDynamic
(`butts2026networkdynamic`), so consecutive spells and window boundaries
are unambiguous; `-Inf` and `Inf` are not allowed as edge times (a
censored spell is an open question, section 22). An edge is active at
time `t` when `start <= t < end`, or `t == start` for a contact. An edge
may be active only while both endpoints are present (section 9; Pajek's
consistency rule, `denooy2018exploratory`).

**Validation.** `E_TIME_TYPE`, `E_TIME_MISSING`, `E_TIME_NOT_FINITE`,
`E_SPELL_ORDER`, `E_WAVE_TYPE`, `E_WAVE_MISSING`, `E_DECLARATION_CONFLICT`
(`end` without `start`; `wave` with `start` or `end`),
`E_SPELL_OUTSIDE_NODE`.

**Class.** Adds `temporal_net_frame`.

**Print.** `time: edge spells, 2010-01-04 to 2012-06-30`, `time:
contacts (1 203 events, 2010-01-04 to 2010-01-11)`, or `time: 3 waves`.

**Example.** edges, `start = "from", end = "to"` (days):

| source | target | from | to |
|---|---|---|---|
| ann | bob | 0 | 10 |
| ann | bob | 12 | 15 |
| bob | cy | 3 | 3 |

Rows 1 and 2 are two spells of one tie, not parallel edges; row 3 is a
contact.

## 9. Node spells and node changes

**Definition.** A node is present during one or more spells and does not
exist for measures outside them (networkDynamic's vertex activity,
`butts2026networkdynamic`; RSiena's composition change,
`ripley2026manual`; the stream graphs of `latapy2018stream`). A node
attribute may change value at known instants (networkDynamic's
temporally extended attributes; RSiena's changing covariates). Decision
19 (provisional).

**Representation.** The common case, one presence spell per node, is two
node roles; the general case is a long table. The two cannot both be
set. A node absent from the `node_spells` table, or with a missing
`node_start` or `node_end`, is present throughout (unbounded) on that
side.

| Table | Role | Type | Constraints |
|---|---|---|---|
| nodes | `node_start` | integer or double | missing means `-Inf`; `node_start <= node_end` |
| nodes | `node_end` | integer or double | missing means `Inf` |
| node_spells | `id` | the node identifier type | every value is a node |
| node_spells | `start` | integer or double | missing means `-Inf` |
| node_spells | `end` | integer or double | missing means `Inf`; `end >= start`; spells of one node do not overlap unless `overlapping_spells = TRUE` |
| node_changes | `id` | the node identifier type | every value is a node |
| node_changes | `time` | integer or double | not missing, finite |

Spells of one node are sorted by start; two spells overlap when the
later one starts before the earlier one ends (`[0, 5)` and `[5, 9)` do
not overlap). A node change at `time` sets the attribute columns of that
row from `time` onwards; the node table holds the value before the first
change. The attribute columns of `node_changes` stay with the binding;
the engine validates only `id` and `time`. An edge spell must lie inside
one presence spell of each endpoint; a contact at `t` is inside `[a, b)`
when `a <= t < b`, or `t == a == b`.

**Validation.** On the roles: `E_TIME_TYPE`, `E_NODE_SPELL_ORDER`. On
`node_spells`: `E_TABLE_ROLE_UNDECLARED` (no `id`, `start`, or `end`
role), `E_TIME_TYPE`, `E_NODE_SPELL_ORDER`, `E_UNKNOWN_NODE`,
`E_MISSING_ENDPOINT` (a missing identifier), `E_NODE_SPELL_OVERLAP`,
`E_DECLARATION_CONFLICT` (both roles and the table). On `node_changes`:
`E_TABLE_ROLE_UNDECLARED`, `E_TIME_TYPE`, `E_TIME_MISSING`,
`E_TIME_NOT_FINITE`, `E_UNKNOWN_NODE`, `E_MISSING_ENDPOINT`. Across
tables: `E_SPELL_OUTSIDE_NODE` lists the edges active while an endpoint
is absent.

**Class.** Adds `temporal_net_frame` (node spells alone make a network
temporal).

**Print.** `time: edge spells, 0 to 15; node spells (table, 5 spells for
3 nodes); 4 attribute changes`.

**Example.** node_spells (`id = "id", start = "from", end = "to"`):

| id | from | to |
|---|---|---|
| ann | 0 | 5 |
| ann | 8 | 20 |
| bob | 0 | 20 |

node_changes (`id = "id", time = "time"`; `status` stays in the binding):

| id | time | status |
|---|---|---|
| ann | 8 | member |
| bob | 12 | staff |

With the edges of section 8, row 1 (`ann-bob`, `[0, 10)`) is an error:
`ann` is absent during `[5, 8)`.

## 10. Spatial networks

**Definition.** Nodes, and possibly edges, have geometry in one
coordinate reference system; adjacency still comes from the edge table,
never from geometry (`barthelemy2011spatial`). Planarity is a graph
property tested separately.

**Representation.** Geometry stays with the binding (sf in R,
`pebesma2018sf`; shapely in Python). The engine receives representative
coordinates and the geometry flags.

| Table | Role | Type | Constraints |
|---|---|---|---|
| nodes | `x` | integer or double | set together with `y`; finite |
| nodes | `y` | integer or double | finite; in `[-90, 90]` when the coordinates are longitude and latitude |

Geometry flags: whether each table carries geometry, each table's
coordinate reference system (an EPSG code, WKT, or proj string), and
what the representative coordinates are (`unknown`, `lonlat`, `planar`).
The two tables' systems must agree. Edge geometry need not touch the node
points (sfnetworks requires it, `vandermeer2026sfnetworks`; here an
opt-in check does).

**Validation.** `E_COORD_ROLE_PAIR`, `E_COORD_TYPE`, `E_COORD_NOT_FINITE`,
`E_LATITUDE_RANGE`, `E_CRS_MISMATCH`.

**Class.** Adds `spatial_net_frame`.

**Print.** `space: node points (EPSG:4326), edge lines`.

**Example.** nodes, `x = "lon", y = "lat"`, coordinates `lonlat`:

| id | lon | lat |
|---|---|---|
| ann | -76.6 | 39.3 |
| bob | -77.0 | 38.9 |

## 11. Spatio-temporal networks

**Definition.** Spatial and temporal roles together, with fixed
geometries and timed edges or nodes (`williams2016spatiotemporal`).
Trajectories are not represented: a node whose position changes is
modelled as node changes (section 9) with a geometry attribute.

**Representation, validation.** The union of sections 8 to 10; no new
role.

**Class.** `spatiotemporal_net_frame`, which inherits from both
`spatial_net_frame` and `temporal_net_frame` (decision 17).

**Print.** The `time` and `space` lines of the two.

## 12. Ego networks

**Definition.** The ego network of a node at order k is the node (ego),
every node within k steps (the alters), and the edges among them
(`marsden2002egocentric`; `everett2005ego`); the order-1.5 network drops
ego and its ties. An egocentric dataset is a collection of such networks,
one per sampled ego, with no ties across egos (`burt1984network`;
`perry2018egocentric`). Decision 18.

**Representation.** One `net_frame` holds every ego's network; nodes and
edges carry the ego they belong to, so every measure runs per ego.

| Table | Role | Type | Constraints |
|---|---|---|---|
| nodes | `ego` | text or integer | not missing; one ego per node row |
| edges | `ego` (optional) | same as the node role | not missing; equals the ego of both endpoints |

An ego's own row, when egos are nodes, carries itself as `ego`. An alter
named by two egos is two node rows with distinct identifiers (as egor's
per-ego alter identifiers). Without an edge `ego` role the edge's ego is
its source's; the target must agree. Ego labels are numbered by first
appearance in the node table.

Declaration `egos_are_nodes`: `TRUE` requires every ego identifier to be
a node identifier (ego included, as `make_ego_graph()` of
`csardi2006igraph` and `ego_graph()` of `hagberg2008exploring` do);
`FALSE` says egos are not rows (order-1.5 data, egor's default); undeclared
observes and the summary reports how many egos are nodes. An ego
identifier is matched to a node identifier by value, so the two columns
should have one type.

The bridge to egor's three tables (egos, alters, alter-alter ties, keyed
by the ego identifier; `rawlings2023network`, chapter 6; `shafie2026r4sna`)
is: egos become node rows with `ego = id` (when `egos_are_nodes`), alters
node rows with their ego, ego-alter and alter-alter ties edge rows with
the ego; and back.

**Validation.** `E_EGO_TYPE` (double or logical column), `E_EGO_MISSING`
(missing values, or nodes added from the edges), `E_EGO_UNKNOWN`
(declared `egos_are_nodes`, but the listed node rows name an ego that is
not a node), `E_EGO_SELF` (a node that is an ego belongs to another
ego's network), `E_EGO_MISMATCH` (the listed edges join nodes of
different egos, or carry an ego that differs from their endpoints'),
`E_DECLARATION_CONFLICT` (an edge `ego` role without the node role).

**Class.** None (section 22). Ego subgraphs (the nodes and edges of one
ego) are built by the engine on request with their row mappings.

**Print.** `egos: 25 (egos are nodes; 4-12 alters each)` or `egos: 25
(alters only)`.

**Example.** `egos_are_nodes = TRUE`, `ego = "ego"` on both tables.

nodes:

| id | ego |
|---|---|
| ann | ann |
| a1 | ann |
| a2 | ann |
| bob | bob |
| b1 | bob |

edges:

| source | target | ego |
|---|---|---|
| ann | a1 | ann |
| ann | a2 | ann |
| a1 | a2 | ann |
| bob | b1 | bob |

## 13. Missing data

**Definition.** An absent tie was observed not to exist (no row). A
missing tie is a pair whose status was not observed. A non-respondent is
a node whose own reports were not collected; ties to it reported by
others are known. (`kossinets2006effects`; `robins2004missing`;
`handcock2010modeling`; network's `na` attribute, `butts2008network`;
RSiena's codes, `ripley2026manual`.) Decision 14.

**Representation.** Two logical flags; unobserved pairs are stored only
when named, never materialized.

| Table | Role | Type | Constraints |
|---|---|---|---|
| edges | `edge_missing` | logical | not missing; `TRUE` marks a pair whose status is unobserved |
| nodes | `node_missing` | logical | not missing; `TRUE` marks a non-respondent |

A row with `edge_missing = TRUE` is an unobserved pair, neither a tie nor
a non-tie: it is kept in the edge table, counted, carries its relation,
layers, and spell like any row (so it says which relation, layer, or
period is unobserved), and is **excluded from the built graph**. The
graph used by every measure holds only observed ties, and the engine
returns the graph-edge-to-row mapping and, separately, the missing pairs
(node indices and rows). A pair cannot be both observed and unobserved:
a flagged row with the same key (section 14) as an unflagged row is an
error. Two flagged rows with one key are parallel rows and counted as
such.

A non-respondent's unreported ties are implied by rule, not written out:
in a directed network, the out-ties of a non-respondent are unobserved;
in an undirected network, the pairs between two non-respondents are
unobserved and a pair with one respondent is known from that respondent's
report (`rawlings2023network`, chapter 4). Measures refuse a network with
unresolved missingness (either flag set to `TRUE` anywhere) unless an
`na_action` says how to resolve it: treat as absent, treat as present,
or exclude the pairs from denominators; this is the measures' contract,
not validation's, which only counts.

**Validation.** `E_MISSING_ROLE_TYPE` (not logical, or missing values),
`E_MISSING_CONFLICT` (a pair observed and unobserved at once),
`N_MISSING_NODES`, `N_MISSING_EDGES` (counts, so a binding can warn).

**Class.** None.

**Print.** `missing: 3 non-respondents, 5 unobserved pairs` (the line is
omitted when both are zero).

**Example.** edges, `edge_missing = "na"`:

| source | target | na |
|---|---|---|
| ann | bob | FALSE |
| ann | cy | TRUE |
| bob | cy | FALSE |

The built graph has two edges (rows 1 and 3); the pair `ann-cy` is
reported as unobserved.

## 14. Loops and parallel edges

**Definition.** A self-loop joins a node to itself. Two edges are
parallel when they have the same endpoints (ordered when directed) and
the same values of every role that distinguishes edges: the relation,
the time spell, and the layers (decision 7; `newman2018networks`,
chapter 6, for multigraphs; igraph's endpoint-only rule,
`csardi2006igraph`, and network's `multiple` flag, `butts2008network`,
are the alternatives). A network with neither is simple. Networks are
multigraphs by default: loops and parallel edges are kept, counted, and
each measure documents how it treats them.

**The key of an edge** is `(a, b, la, lb, relation, start, end)` where
`(a, la)` and `(b, lb)` are the endpoints with their layers; in an
undirected network the two node-layer pairs are ordered so that `(a, la)
<= (b, lb)`. Missing roles contribute nothing. A loop is an edge with
`a == b` and `la == lb`; a coupling (section 7) is not a loop.

**Declarations.** `loops = FALSE` makes every loop an error; `multiple =
FALSE` makes every repeated key an error; `TRUE` allows them silently;
undeclared reports a note. `simplify()` groups parallel edges, combines
their weights (sum, mean, min, max, first, count; igraph's
`edge.attr.comb` and Pajek's "Remove Multiple Lines" offer the same), and
returns an index mapping so the binding aggregates other attributes.

**Validation.** `E_LOOP_DECLARED`, `E_PARALLEL_DECLARED`, `N_LOOPS`,
`N_PARALLEL`.

**Print.** `simple: no loops, no parallel edges` or `simple: 2 loops, 5
parallel edges`.

**Example.** Undirected, `relation = "rel"`:

| source | target | rel |
|---|---|---|
| ann | bob | friend |
| bob | ann | friend |
| ann | bob | advice |
| cy | cy | friend |

Rows 1 and 2 are parallel; row 3 is not; row 4 is a loop.

## 15. Network-level attributes

**Definition.** Named metadata about the network as a whole with no
structural meaning (title, source, collection date, design).

**Representation.** A named list in the binding (`graph`); the engine
never sees it. Structural declarations (`directed`, the roles, the
declared properties, the geometry flags) are fields of the object, not
entries of the list (network stores its flags as attributes,
`butts2008network`; igraph and NetworkX keep graph attributes apart).
Reserved optional keys with a documented meaning: `title`, `source`,
`design` (boundary, sampling, instrument, fixed-choice limit;
`marsden1990network`, `handcock2010modeling`), `relations` (a data frame
of relation-level attributes), `time_class` and `tz` (how times print),
`ego_order` (for extracted ego networks).

**Validation.** None. **Print.** Not shown.

## 16. Hypergraphs and simplicial complexes (out of scope)

A hypergraph is a node set with hyperedges that are node subsets of any
size; a simplicial complex is a hypergraph closed under subsets
(`battiston2020networks`, `bianconi2021higher`). Neither is first-class
(decision 12): an edge has two endpoints. The documented exchange is a
two-mode network with `mode = c("node", "hyperedge")` and
`between_modes_only = TRUE`, one hyperedge node per hyperedge and one
edge per membership; the clique expansion is the ordinary two-mode
projection. `torres2021why` explain why the two objects differ although
the incidence structure is the same; XGI (`landry2023xgi`) is the
reference library for higher-order data.

**Example.** The hyperedges `{ann, bob, cy}` and `{bob, cy}`:

nodes:

| id | mode |
|---|---|
| ann | node |
| bob | node |
| cy | node |
| h1 | hyperedge |
| h2 | hyperedge |

edges:

| source | target |
|---|---|
| ann | h1 |
| bob | h1 |
| cy | h1 |
| bob | h2 |
| cy | h2 |

## 17. Combining configurations

Every configuration is a set of roles and declarations, so any
combination exists. The matrix says which pairs need a rule beyond the
union of their checks (a dot means none).

| | mode | weight | sign | relation | layers | time | node spells | space | ego | missing |
|---|---|---|---|---|---|---|---|---|---|---|
| **mode** | | . | . | . | A | . | . | . | . | . |
| **weight** | . | | B | . | . | . | . | . | . | . |
| **sign** | . | B | | . | . | . | . | . | . | . |
| **relation** | . | . | . | | C | . | . | . | . | D |
| **layers** | A | . | . | C | | E | F | . | . | D |
| **time** | . | . | . | . | E | | G | H | . | D |
| **node spells** | . | . | . | . | F | G | | . | . | . |
| **space** | . | . | . | . | . | H | . | | . | . |
| **ego** | . | . | . | . | . | . | . | . | | I |
| **missing** | . | . | . | D | D | D | . | . | I | |

- **A** (mode + layers): modes belong to actors, so an actor has one
  mode in every layer; `between_modes_only` applies to intralayer and
  interlayer edges alike, and coupling edges are exempt.
- **B** (weight + sign): the weight is the magnitude and must not be
  negative (`E_WEIGHT_NEGATIVE`).
- **C** (relation + layers): both may be set; the parallel-edge key holds
  both; the multiplex special case (relations as layers with couplings)
  is produced by a helper, not required.
- **D** (missing + relation, layers, or time): an unobserved pair carries
  its relation, layers, and spell, and is unobserved only there; the
  conflict check (`E_MISSING_CONFLICT`) compares full keys.
- **E** (layers + time): every edge row has its own spell, so a tie has
  one spell per layer; node spells are per actor, not per node-layer
  (section 22).
- **F** (layers + node spells): a coupling edge, like any edge, must lie
  inside the presence of its (single) actor.
- **G** (time + node spells): every edge spell or contact lies inside one
  presence spell of each endpoint (`E_SPELL_OUTSIDE_NODE`); waves are
  compared in wave units (`[k, k + 1)` inside the node's spell).
- **H** (time + space): `spatiotemporal_net_frame`; geometry is fixed.
- **I** (ego + missing): flags are per row, so a non-respondent or an
  unobserved pair belongs to one ego's network; measures per ego resolve
  missingness per ego.

Classes compose in the fixed order of section 20; `multilayer` +
`temporal`, for instance, is `c("multilayer_net_frame",
"temporal_net_frame", "net_frame")`.

## 18. The print card

The card states, in this order, one line per item, omitting lines whose
configuration is absent (decision 2):

```
# net_frame: 34 nodes, 78 edges
#   directed:  no        mode: one-mode (34 actors)
#   weights:   weight (strength, 1-7)
#   signs:     sign (70 positive, 8 negative)
#   relations: friendship (42), advice (36) (multiplex)
#   layers:    work, online (2); 70 intralayer, 8 coupling, 0 interlayer
#   time:      edge spells, 0 to 15; node spells (table, 40 spells); 4 changes
#   space:     node points (EPSG:4326)
#   egos:      25 (egos are nodes)
#   missing:   3 non-respondents, 5 unobserved pairs
#   simple:    no loops, no parallel edges
#   components: 1    isolates: 0    density: 0.139
#
# Nodes (34 x 3)              Edges (78 x 4)
#   ...                         ...
```

Every number on the card comes from one engine call (the summary), which
reports: node and edge counts, the number of observed (graph) edges,
direction, mode labels and sizes with the between-modes declaration, the
weight column, meaning, and range, the sign counts, relation labels and
sizes, layer labels and intralayer sizes with the intralayer, coupling,
and interlayer totals and whether a layer table was given, the time kind
(static, contacts, edge spells, waves) and range, whether node spells
come from roles or a table and how many spells and changes there are,
the geometry flags and coordinate reference system, the number of egos
and how many are nodes with nodes and edges per ego, the counts of
non-respondents and unobserved pairs, loops, parallel edges, isolates,
weak components, and density (loops excluded from the denominator; NaN
when undefined). Components and density are computed on the observed
graph.

## 19. Role catalogue

Node roles (C codes in parentheses, positions in the engine's list):

| Role | Type | Meaning |
|---|---|---|
| `id` (0) | integer, double, text | the node identifier; unique, not missing |
| `mode` (1) | text, integer, logical | the node's mode (section 3) |
| `node_start` (2) | numeric | start of the node's presence spell; missing is `-Inf` |
| `node_end` (3) | numeric | end of the presence spell; missing is `Inf` |
| `x` (4) | numeric | representative x coordinate or longitude |
| `y` (5) | numeric | representative y coordinate or latitude |
| `node_missing` (6) | logical | `TRUE` marks a non-respondent |
| `ego` (7) | text, integer | the ego whose network the node belongs to |

Edge roles:

| Role | Type | Meaning |
|---|---|---|
| `source` (0) | as `id` | the source endpoint (tail of a directed edge) |
| `target` (1) | as `id` | the target endpoint |
| `weight` (2) | numeric | the weight, with a declared `weight_type` |
| `relation` (3) | text, integer, logical | the relation of a multiplex network |
| `sign` (4) | integer `-1`/`1` or logical | the sign |
| `start` (5) | numeric | start of the edge's spell, or the time of a contact |
| `end` (6) | numeric | end of the spell |
| `wave` (7) | integer | the observation wave; `[k, k + 1)` |
| `edge_missing` (8) | logical | `TRUE` marks an unobserved pair |
| `source_layer` (9) | text, integer, logical | the layer of the source endpoint |
| `target_layer` (10) | text, integer, logical | the layer of the target endpoint |
| `ego` (11) | as the node `ego` | the ego whose network the edge belongs to |

Layer table roles: `layer` (0): the layer identifier (text, integer,
logical; unique, not missing).

Node spells table roles: `id` (0), `start` (1), `end` (2).

Node changes table roles: `id` (0), `time` (1).

Declarations (C codes): `directed` (0), `loops` (1), `multiple` (2),
`between_modes_only` (3), `add_nodes` (4), `egos_are_nodes` (5),
`couplings_only` (6), `overlapping_spells` (7). Values: 0 false, 1 true,
-1 undeclared (1 to 3 and 5 to 7 only).

Table codes: `nodes` (0), `edges` (1), `layers` (2), `node_spells` (3),
`node_changes` (4).

## 20. Class catalogue

| Class | Added when | Inherits from |
|---|---|---|
| `net_frame` | always | |
| `temporal_net_frame` | an edge time role (`start`, `wave`) or node spells (roles or table) is set | `net_frame` |
| `spatial_net_frame` | coordinates or geometry are set | `net_frame` |
| `spatiotemporal_net_frame` | both of the above | `spatial_net_frame`, `temporal_net_frame`, `net_frame` |
| `multilayer_net_frame` | the layer roles are set | `net_frame` |

The class vector lists, in this fixed order, each class whose condition
holds: `multilayer_net_frame`, `spatiotemporal_net_frame`,
`spatial_net_frame`, `temporal_net_frame`, `net_frame`. So a multilayer
temporal network is `c("multilayer_net_frame", "temporal_net_frame",
"net_frame")`, a spatio-temporal one `c("spatiotemporal_net_frame",
"spatial_net_frame", "temporal_net_frame", "net_frame")`, and a
multilayer spatio-temporal one has all five. Modes, weights, signs,
relations, egos, and missing flags add no class: they are roles that
every method reads through the object.

## 21. Diagnostic catalogue

Codes beginning `E_` are errors (the frame cannot be built), `N_` notes.
Each diagnostic names its table (`nodes`, `edges`, `layers`,
`node_spells`, `node_changes`, or the network as a whole) and lists the
first 50 offending rows with the total count.

Structure:

| Code | Table | Condition |
|---|---|---|
| `E_ENDPOINTS_UNDECLARED` | network | no source and target roles |
| `E_ID_UNDECLARED` | network | a node table without an identifier role |
| `E_TABLE_ROLE_UNDECLARED` | network | a layers, node_spells, or node_changes table without its required roles (new) |
| `E_ENDPOINT_TYPE` | network | endpoints not integer, double, or text, or of two types |
| `E_ID_TYPE` | network | identifier type differs from the endpoints' |
| `E_MISSING_ENDPOINT` | edges, node_spells, node_changes | a missing endpoint or identifier |
| `E_EMPTY_ID` | nodes, edges | an empty text identifier |
| `E_MISSING_ID` | nodes | a missing node identifier |
| `E_DUPLICATE_ID` | nodes | a repeated node identifier |
| `E_UNKNOWN_NODE` | edges, node_spells, node_changes | an identifier absent from the node table |

Declarations:

| Code | Table | Condition |
|---|---|---|
| `E_LOOP_DECLARED` | edges | loops present but declared absent |
| `E_PARALLEL_DECLARED` | edges | parallel edges present but declared absent |
| `E_MODE_WITHIN` | edges | within-mode edges under `between_modes_only` |
| `E_MODE_UNDECLARED` | network | `between_modes_only` without a mode role |
| `E_INTERLAYER_DECLARED` | edges | non-coupling interlayer edges under `couplings_only` (new) |
| `E_EGO_UNKNOWN` | nodes | an ego that is not a node under `egos_are_nodes` (new) |
| `E_DECLARATION_CONFLICT` | network | contradictory declarations or roles |

Weights, modes, relations, signs:

| Code | Table | Condition |
|---|---|---|
| `E_WEIGHT_TYPE` | network | weight column not numeric |
| `E_WEIGHT_TYPE_MISSING` | network | no declared weight meaning |
| `E_WEIGHT_NOT_FINITE` | edges | missing, NaN, or infinite weight |
| `E_DISTANCE_NOT_POSITIVE` | edges | a distance that is not positive |
| `E_WEIGHT_NEGATIVE` | edges | a negative magnitude in a signed network |
| `E_MODE_TYPE` | network | mode column is double |
| `E_MODE_MISSING` | nodes, edges | a node without a mode (edges: added nodes) |
| `E_RELATION_TYPE` | network | relation column is double |
| `E_RELATION_MISSING` | edges | an edge without a relation |
| `E_SIGN_TYPE` | network | sign column not integer or logical |
| `E_SIGN_VALUE` | edges | a sign other than -1 or 1, or missing |

Layers (all new):

| Code | Table | Condition |
|---|---|---|
| `E_LAYER_ROLE_PAIR` | network | one layer role without the other |
| `E_LAYER_TYPE` | network | a layer column is double, or the two differ in type |
| `E_LAYER_MISSING` | edges | a missing layer value |
| `E_LAYER_UNKNOWN` | edges | a layer absent from the layer table |
| `E_LAYER_ID_MISSING` | layers | a missing layer identifier |
| `E_LAYER_ID_DUPLICATE` | layers | a repeated layer identifier |

Egos (all new):

| Code | Table | Condition |
|---|---|---|
| `E_EGO_TYPE` | network | ego column is double or logical, or the edge and node columns differ in type |
| `E_EGO_MISSING` | nodes, edges | a missing ego (edges: added nodes, or missing edge egos) |
| `E_EGO_SELF` | nodes | a node that is an ego belongs to another ego |
| `E_EGO_MISMATCH` | edges | an edge whose endpoints, or whose own ego, disagree |

Time:

| Code | Table | Condition |
|---|---|---|
| `E_TIME_TYPE` | network | a time column not numeric |
| `E_TIME_MISSING` | edges, node_changes | a missing time |
| `E_TIME_NOT_FINITE` | edges, node_changes | an infinite time |
| `E_SPELL_ORDER` | edges | an edge ending before it starts |
| `E_NODE_SPELL_ORDER` | nodes, node_spells | a node leaving before it arrives |
| `E_NODE_SPELL_OVERLAP` | node_spells | two spells of one node overlap (new) |
| `E_SPELL_OUTSIDE_NODE` | edges | an edge active while an endpoint is absent |
| `E_WAVE_TYPE` | network | wave column not whole numbers |
| `E_WAVE_MISSING` | edges | an edge without a wave |

Space:

| Code | Table | Condition |
|---|---|---|
| `E_COORD_TYPE` | network | a coordinate column not numeric |
| `E_COORD_ROLE_PAIR` | network | `x` without `y` or the reverse |
| `E_COORD_NOT_FINITE` | nodes | a missing or infinite coordinate |
| `E_LATITUDE_RANGE` | nodes | a latitude outside [-90, 90] |
| `E_CRS_MISMATCH` | network | the two tables' systems differ |

Missing data:

| Code | Table | Condition |
|---|---|---|
| `E_MISSING_ROLE_TYPE` | network | a flag column not logical, or with missing values |
| `E_MISSING_CONFLICT` | edges | a pair observed and unobserved at once (new) |

Notes:

| Code | Table | Condition |
|---|---|---|
| `N_NODES_ADDED` | edges | nodes added from the endpoints |
| `N_ISOLATES` | nodes | nodes without edges |
| `N_ZERO_BASED_IDS` | network | integer identifiers 0 to n - 1 |
| `N_IDS_COERCED` | network | integer and double identifiers matched as numbers |
| `N_LOOPS` | edges | loops present, none declared |
| `N_PARALLEL` | edges | parallel edges present, none declared |
| `N_MISSING_NODES` | nodes | non-respondents present |
| `N_MISSING_EDGES` | edges | unobserved pairs present |

## 22. Open questions

1. **Node-layer presence and attributes.** The model keeps one row per
   actor, so an actor is present in every layer and has one value of
   each attribute across layers. Per-layer presence (an actor absent from
   a layer) and per-layer attributes would need a node-layer table
   (`id`, `layer`, attributes), the exact form of `kivela2014multilayer`;
   not added until a multilayer measure needs it.
2. **Censored spells.** Edge times must be finite. networkDynamic allows
   `-Inf` and `Inf` as "onset unknown" and "terminus unknown"; node
   spells allow them here, edge spells do not. Allowing them for edges
   would make the time range of the print card infinite.
3. **An `ego_net_frame` class.** Decision 18 makes ego data a role, not
   a class, so no class is added. If per-ego printing and plotting want
   dispatch, the class would go first in the order of section 20.
4. **The default `na_action`.** Measures refuse unresolved missingness;
   whether `"error"` stays the only default, and how `exclude` affects
   each denominator, is decided measure by measure.
5. **Non-respondents who report ties.** A directed edge whose source is a
   non-respondent is accepted without comment; a note might help.
6. **Node changes outside presence.** A change at a time when the node is
   absent is accepted; validation could note it.
7. **Structural zeros and ones.** Pairs that cannot or must be tied by
   design have no role; a rule or a pair table (decision 11) is still
   open.
8. **Several aspects.** Layers are one categorical column; a layer that
   is a combination of aspects (relation by time) is one label with
   attributes in the layer table.
9. **Structural zeros and ones.** Pairs that cannot be tied (a person
   cannot advise someone outside their unit) or must be tied (a
   reporting line) are distinct from unobserved pairs. RSiena codes them
   10 and 11 in adjacency matrices; ERGM fitting treats them as
   constraints. Proposed: an edge-level `fixed` role (`"zero"`, `"one"`)
   plus a pair-rule form for whole blocks (e.g. no ties across modes or
   outside a group), carried to models as constraints and excluded from
   descriptive denominators where a measure documents that.
