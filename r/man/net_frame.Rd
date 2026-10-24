\name{net_frame}
\alias{net_frame}
\alias{net_diagnostics}
\alias{print.net_frame}
\title{Create a network from node and edge data frames}
\description{
A network is a pair of data frames: \code{nodes}, with one row per node, and
\code{edges}, with one row per edge. Any other columns are attributes and are
carried along untouched. Column roles say which columns hold node
identifiers, edge endpoints, and, optionally, weights, modes, relations,
and signs. The tables are validated in the engine before anything is
computed: every problem is reported at once, with the rows concerned, and
advisory notes say what was noticed. Printing a net_frame shows a card of
its properties and the first rows of both tables.
}
\usage{
net_frame(edges, nodes = NULL, directed = TRUE, source = "source",
          target = "target", id = "id", weight = NULL, weight_type = NULL,
          mode = NULL, between_modes_only = NULL, relation = NULL, sign = NULL,
          start = NULL, end = NULL, wave = NULL, node_start = NULL,
          node_end = NULL, coords = NULL, crs = NULL, source_layer = NULL,
          target_layer = NULL, layers = NULL, layer = "layer",
          couplings_only = NULL, ego = NULL, edge_ego = NULL,
          egos_are_nodes = NULL, node_missing = NULL, edge_missing = NULL,
          node_spells = NULL, node_changes = NULL, spell_start = "start",
          spell_end = "end", change_time = "time", overlapping_spells = NULL,
          loops = NULL, multiple = NULL, add_nodes = FALSE, quiet = FALSE,
          graph = list())

net_diagnostics(x, threads = 0L)

\method{print}{net_frame}(x, ..., n = 6L, width = getOption("width"))
}
\arguments{
\item{edges}{A data frame of edges.}
\item{nodes}{A data frame of nodes, or \code{NULL} to derive the node set from
the edge endpoints in order of first appearance. The derived table has one
column, called \code{id} (the name given by \code{id}).}
\item{directed}{Whether edges are directed: \code{TRUE} or \code{FALSE}.}
\item{source, target}{The names of the columns of \code{edges} holding the
endpoints. The defaults are \code{"source"} and \code{"target"}; column
names and positions are never guessed, so other names are given here.}
\item{id}{The name of the column of \code{nodes} holding node identifiers
(default \code{"id"}), or the name of the identifier column to create when
\code{nodes} is \code{NULL}.}
\item{weight}{Optional name of a column of \code{edges} holding edge weights.
Without one, the network is binary: an edge either exists or it does not,
and algorithms use their faster unweighted forms.}
\item{weight_type}{What the weights mean, required whenever \code{weight} is
given: \code{"strength"} means larger weights bring nodes closer, as with
contact frequency, tie intensity, or capacity, and path-based measures
then use 1/weight as the length (Newman 2001; Opsahl, Agneessens, and
Skvoretz 2010); \code{"distance"} (the convention of igraph's and
NetworkX's path functions) means larger weights keep nodes farther apart,
as with travel time, cost, or latency, and distances must be positive.
There is no default: a weight column without a declared meaning is an
error. Measures that add up weights, such as \code{\link{net_strength}},
use them as given either way.}
\item{mode}{Optional name of a column of \code{nodes} holding each node's
mode (a discrete column: text, integer, logical, or a factor). Modes
declared here may have ties within them; \code{\link{net_twomode}} declares
modes whose ties must cross modes.}
\item{relation}{Optional name of a column of \code{edges} holding the
relation each edge belongs to, for a multiplex network kept in one long
edge table (see \code{\link{net_combine_relations}}).}
\item{sign}{Optional name of a column of \code{edges} holding the sign of
each tie: \code{-1} and \code{1}, or a logical column where \code{TRUE} is
positive. Weights of a signed network are magnitudes and must not be
negative. Provisional: the role is validated and stored; no measure uses it
yet.}
\item{between_modes_only}{Declaration for a \code{mode} column:
\code{TRUE} states that every edge joins nodes of different modes (a
two-mode or k-partite network) and lists the edges that do not;
\code{FALSE} allows ties within modes; \code{NULL} leaves it undeclared.
\code{\link{net_twomode}} sets it to \code{TRUE}.}
\item{start, end}{Optional names of columns of \code{edges} holding when
each edge is active: \code{start} alone makes every edge an instantaneous
contact; \code{start} and \code{end} make it the half-open spell
\eqn{[start, end)}. Numbers, \code{Date}, or \code{POSIXct}; every time
column of the network has one class. See \code{\link{net_temporal}}.}
\item{wave}{Optional name of a column of \code{edges} holding whole wave
numbers of a panel: wave \eqn{k} is the spell \eqn{[k, k + 1)}. Excludes
\code{start} and \code{end}. See \code{\link{net_panel}}.}
\item{node_start, node_end}{Optional names of columns of \code{nodes}
holding one presence spell per node; a missing value leaves the spell
open on that side. Edges must fall within the presence of both endpoints.
For several spells per node use \code{node_spells} instead.}
\item{coords}{Optional names of two numeric columns of \code{nodes}:
longitude and latitude, or x and y. See \code{\link{net_spatial}}.}
\item{crs}{What \code{coords} hold: \code{"lonlat"} (the default when
\code{coords} is given; latitudes must lie in [-90, 90]) or
\code{"planar"}.}
\item{source_layer, target_layer}{Optional names of two columns of
\code{edges} holding the layer of each endpoint, for a multilayer network:
equal layers make an intralayer edge, an actor joined to itself across
layers is a coupling, and other edges across layers are interlayer edges
(Kivela et al. 2014). Given together. See \code{\link{net_split_layers}}.}
\item{layers}{Optional data frame with one row per layer and any layer
attributes; every edge layer must be one of its identifiers.}
\item{layer}{The name of the identifier column of \code{layers} (default
\code{"layer"}).}
\item{couplings_only}{Declaration for a multilayer network: \code{TRUE}
states that every edge across layers is a coupling (the multiplex case)
and lists the interlayer edges that join different actors; \code{FALSE}
allows any; \code{NULL} leaves it undeclared.}
\item{ego, edge_ego}{Optional names of columns of \code{nodes} and
\code{edges} holding the ego whose network each row belongs to, for
egocentric data held in one net_frame. Without \code{edge_ego} an edge
belongs to its source's ego, and both endpoints must belong to the same
ego. See \code{\link{net_ego_frame}}.}
\item{egos_are_nodes}{Declaration for \code{ego}: \code{TRUE} requires
every ego to be a node row (whose own \code{ego} value is itself);
\code{FALSE} states that egos are not rows (alters only); \code{NULL}
leaves it undeclared.}
\item{node_missing}{Optional name of a logical column of \code{nodes}
marking non-respondents (\code{TRUE}): nodes whose own reports were not
collected. Their unreported ties are unobserved by rule and are not
written out. See \code{\link{net_resolve_missing}}.}
\item{edge_missing}{Optional name of a logical column of \code{edges}
marking unobserved pairs (\code{TRUE}): rows whose tie status was not
observed. They are kept in the table but are neither ties nor non-ties,
and are left out of the graph that measures use.}
\item{node_spells}{Optional data frame of presence spells, with columns
\code{id} (the name given by \code{id}), \code{spell_start}, and
\code{spell_end}; a node may have several. A node not in the table is
present throughout. Excludes \code{node_start} and \code{node_end}. See
\code{\link{net_node_spells}}.}
\item{node_changes}{Optional data frame of node attribute changes, with
columns \code{id}, \code{change_time}, and the attribute columns: a row
sets those attributes of its node from that time on, and \code{nodes}
holds the values before the first change. See \code{\link{net_nodes_at}}.}
\item{spell_start, spell_end}{The names of the start and end columns of
\code{node_spells}.}
\item{change_time}{The name of the time column of \code{node_changes}.}
\item{overlapping_spells}{Declaration for \code{node_spells}:
\code{TRUE} allows the spells of one node to overlap; otherwise
overlapping spells are an error.}
\item{loops, multiple}{Declarations: \code{FALSE} states that the network
has no self-loops (no parallel edges), and any violating row is an error;
\code{TRUE} states that they are allowed; \code{NULL} (the default) leaves
the property undeclared, in which case it is observed and reported as a
note.}
\item{add_nodes}{Whether endpoints absent from \code{nodes} are appended to
the node table, in order of first appearance, with \code{NA} attributes
(reported as a note). By default they are an error that lists them.}
\item{quiet}{Whether to suppress the advisory notes.}
\item{graph}{A named list of network-level attributes (for example a title or
data source), carried through conversions and file formats that support
them.}
\item{x}{A net_frame.}
\item{threads}{Number of threads for the validation; \code{0} uses all
available cores.}
\item{...}{Ignored.}
\item{n}{Rows of each table to preview.}
\item{width}{The console width the preview is laid out for.}
}
\details{
\subsection{Validation}{
The checks run in the engine, so the R and Python packages apply the same
rules, and every problem is reported in one error rather than the first
one found. They fall into four groups, modelled on the checks of
\code{as.network.data.frame()} in the network package (Butts 2008):
\itemize{
\item \emph{Structure}: the endpoint columns are integer, double, or text
and of one type, the identifier column matches it, no endpoint or
identifier is missing or an empty string, identifiers are unique, and
every endpoint is in the node table (the missing ones are listed).
\item \emph{Declared properties}: a declared \code{loops = FALSE},
\code{multiple = FALSE}, \code{between_modes_only = TRUE},
\code{couplings_only = TRUE}, or \code{egos_are_nodes = TRUE} lists the
rows that violate it; contradictory
declarations are errors, never silent overrides.
\item \emph{Attribute sanity}: weights are numeric and finite (positive
when distances) and their meaning is declared; the mode and relation
columns are discrete and complete; signs are -1 and 1; times are finite,
spells end no earlier than they start and fall within the endpoints'
presence spells; coordinates are finite and latitudes lie in [-90, 90];
sf node and edge tables share one coordinate reference system.
\item \emph{Advisory notes}, shown with \code{message()} unless
\code{quiet = TRUE}: nodes added from the edges, isolates, integer
identifiers that run from 0 to n - 1 (which may be positions in a 1-based
table), identifiers of integer and double type matched as numbers, and
the counts of self-loops and parallel edges when neither was declared.
}
\code{net_diagnostics()} runs the checks again on the current tables and
returns them as a data frame with columns \code{code} (a stable code such
as \code{E_UNKNOWN_NODE} or \code{N_ISOLATES}), \code{severity}
(\code{"error"} or \code{"note"}), \code{table} (\code{"nodes"},
\code{"edges"}, or \code{NA} for the network as a whole), \code{message},
\code{count} (rows affected in all), and \code{rows} (a list column with
the first 50 rows, 1-based).
}
\subsection{Self-loops and parallel edges}{
Loops and parallel edges are kept and counted. Two edges are parallel when
they have the same endpoints (ordered in a directed network) and the same
relation and time spell, so a multiplex or temporal network with one edge
per relation and spell is simple (see \code{\link{net_simplify}}). Each
measure documents how it counts them; \code{convention} selects the
treatment of other tools (see \code{\link{conventions}}).
}
\subsection{Identifiers}{
Identifiers may be character, factor (matched by their labels), integer,
double with whole values, or logical (read as 0 and 1), or a classed vector
such as \code{Date}; \code{integer64} identifiers (bit64) are matched as
text. Without \code{nodes}, the \code{source} and \code{target} columns
must have the same class.
}
\subsection{Editing the tables}{
Add columns to \code{x$nodes} and \code{x$edges} freely. The edge index is
computed here; to add, remove, or reorder rows build a new net_frame, or use
the transformation functions (\code{\link{net_subnetwork}},
\code{\link{net_filter_nodes}}, \code{\link{net_dichotomize}}, ...), which
rebuild it. Measures stop when the number of nodes or edges has changed;
other changes to rows cannot be detected and would misalign results.
Changed weights are picked up (the engine graph is rebuilt and validated
again), so editing \code{x$edges$w} in place works.
}
\subsection{The printed card}{
\code{print()} shows, in order: the class and size; whether edges are
directed; the modes (one-mode, two-mode, or k-mode with the mode labels
and counts); the weights (column, declared meaning, range) or
\code{unweighted}; the relations of a multiplex network; the sign column;
time (static, or contacts, edge spells, or waves with their range and
whether nodes have presence spells) and space (coordinate columns or sf
geometry with its CRS); self-loops and parallel edges, with
\code{(declared)} when their absence was declared; nodes and edges marked
unobserved; and the number of weak components, isolates, and the density.
Then come the first \code{n} rows of the node and edge tables with each
column's type abbreviated (\code{<int>}, \code{<dbl>}, \code{<chr>},
\code{<lgl>}, \code{<fct>}, \code{<date>}, \code{<dttm>}, \code{<POINT>},
...), side by side when the console is wide enough and one after the other
otherwise, followed by how many more rows and columns there are. Every
line starts with \code{#}, so printed output can be pasted into a script.
The output uses only ASCII characters. Temporal and spatial networks add
a line each (see \code{\link{net_temporal}} and \code{\link{net_spatial}}).
\code{\link[=summary.net_frame]{summary()}} gives more.
}
\subsection{Threads}{
Measures that run in parallel take \code{threads}. \code{0} (the default)
uses all available cores, except that only 2 are used when the environment
variable \env{_R_CHECK_LIMIT_CORES_} is set (as it is by
\verb{R CMD check --as-cran}), and at most \env{OMP_THREAD_LIMIT} when that is
set. Results are reproducible for a given number of threads: integer results
never depend on it, and floating-point results can differ in the last bits
between thread counts. Computations run in compiled code and cannot be
interrupted until they return.
}
\subsection{Spatial tables}{
\code{nodes} and \code{edges} may be sf data frames (simple features). Their
geometry columns stay in place and are recorded in \code{x$roles}
(\code{node_geometry}, \code{edge_geometry}, \code{crs}); both tables must
then have the same coordinate reference system, and the network gets the
class \code{spatial_net_frame}. Connectivity still comes from \code{source}
and \code{target}, never from the geometry. See \code{\link{net_to_sf}}.
}
}
\value{An object of class \code{net_frame}: a list with the two tables
(\code{nodes}, \code{edges}), the optional long tables (\code{layers},
\code{node_spells}, \code{node_changes}, or \code{NULL}), \code{directed},
\code{roles} (the column names of each role, with \code{weight_type}),
\code{declared} (the declarations that were made), \code{graph}, and the
engine's edge index. Temporal, spatial, and multilayer networks are
subclasses, composed in a fixed order:
\code{c("temporal_net_frame", "net_frame")},
\code{c("spatial_net_frame", "net_frame")},
\code{c("spatiotemporal_net_frame", "spatial_net_frame",
"temporal_net_frame", "net_frame")}, and
\code{c("multilayer_net_frame", "net_frame")} (which may precede the
others, as in \code{c("multilayer_net_frame", "temporal_net_frame",
"net_frame")}). Modes, weights, signs, relations, egos, and missing-data
flags add no class.
}
\references{
Butts, C. T. (2008). network: A package for managing relational data in R.
\emph{Journal of Statistical Software}, 24(2), 1-36.
\doi{10.18637/jss.v024.i02}

Kivela, M., Arenas, A., Barthelemy, M., Gleeson, J. P., Moreno, Y., &
Porter, M. A. (2014). Multilayer networks. \emph{Journal of Complex
Networks}, 2(3), 203-271. \doi{10.1093/comnet/cnu016}

Newman, M. E. J. (2001). Scientific collaboration networks. II. Shortest
paths, weighted networks, and centrality. \emph{Physical Review E}, 64,
016132. \doi{10.1103/PhysRevE.64.016132}

Opsahl, T., Agneessens, F., & Skvoretz, J. (2010). Node centrality in
weighted networks: Generalizing degree and shortest paths. \emph{Social
Networks}, 32(3), 245-251. \doi{10.1016/j.socnet.2010.03.006}

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods
and Applications}. Cambridge University Press.
\doi{10.1017/CBO9780511815478}
}
\seealso{
\code{\link{as_net_frame}} to convert from other packages, \code{\link{net_read}}
for files, \code{\link{net_twomode}}, \code{\link{net_temporal}}, and
\code{\link{net_spatial}} for the other kinds of network, and
\code{\link{net_simplify}} for loops and parallel edges,
\code{\link{net_split_layers}} for multilayer networks,
\code{\link{net_ego_frame}} for ego networks,
\code{\link{net_resolve_missing}} for missing data,
\code{\link{net_node_spells}} for node presence and changes, and
\code{\link{net_as_signed}} for signed values. The articles
\dQuote{Network types} and \dQuote{Definitions} show every configuration
with its tables, card, and plot.
}
\examples{
edges <- data.frame(source = c("a", "b", "c"), target = c("b", "c", "a"))
g <- net_frame(edges)
g

# The same edges read as directed and as undirected.
op <- par(mfrow = c(1, 2))
plot(g, main = "directed", seed = 1)
plot(net_frame(edges, directed = FALSE), main = "undirected", seed = 1)
par(op)

# Weights must say what they mean.
w <- data.frame(source = c("a", "a", "b"), target = c("b", "c", "c"), n = c(3, 1, 2))
net_frame(w, directed = FALSE, weight = "n", weight_type = "strength")

# Every problem is reported at once, with the rows.
nodes <- data.frame(id = c("a", "b", "b"))
try(net_frame(data.frame(source = c("a", NA), target = c("q", "b")), nodes))

# Self-loops and parallel edges are kept, reported, and drawn.
h <- net_frame(data.frame(source = c("a", "a", "a", "b"), target = c("a", "b", "b", "c")))
net_diagnostics(h)
net_is_simple(h)
plot(h, seed = 1)

# A signed network: negative ties are drawn dashed.
s <- data.frame(source = c("a", "a", "b", "c"), target = c("b", "c", "c", "d"),
                sign = c(1, -1, 1, -1))
plot(net_frame(s, directed = FALSE, sign = "sign"), seed = 1)

# A multiplex network: edges coloured by relation.
krack <- net_example("krackhardt")
krack
plot(krack, seed = 1)

# Unobserved pairs are kept in the table and drawn dotted.
u <- data.frame(source = c("a", "a", "b"), target = c("b", "c", "c"),
                unknown = c(FALSE, TRUE, FALSE))
net_frame(u, directed = FALSE, edge_missing = "unknown")
}
