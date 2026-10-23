\name{net_frame}
\alias{net_frame}
\alias{netframe}
\alias{net_diagnostics}
\alias{print.netframe}
\title{Create a network from node and edge data frames}
\description{
A network is a pair of data frames: \code{nodes}, with one row per node, and
\code{edges}, with one row per edge. Any other columns are attributes and are
carried along untouched. Column roles say which columns hold node
identifiers, edge endpoints, and, optionally, weights, modes, relations,
and signs. The tables are validated in the engine before anything is
computed: every problem is reported at once, with the rows concerned, and
advisory notes say what was noticed. Printing a netframe shows a card of
its properties and the first rows of both tables.
}
\usage{
net_frame(edges, nodes = NULL, directed = TRUE, source = "source",
          target = "target", id = "id", weight = NULL, weight_type = NULL,
          mode = NULL, relation = NULL, sign = NULL, loops = NULL,
          multiple = NULL, add_nodes = FALSE, quiet = FALSE, graph = list())

net_diagnostics(x, threads = 0L)

\method{print}{netframe}(x, ..., n = 6L, width = getOption("width"))
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
\item{x}{A netframe.}
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
\code{multiple = FALSE}, or (for \code{\link{net_twomode}}) the rule that
edges join different modes lists the rows that violate it; contradictory
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
computed here; to add, remove, or reorder rows build a new netframe, or use
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
\code{\link[=summary.netframe]{summary()}} gives more.
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
class \code{spatial_netframe}. Connectivity still comes from \code{source}
and \code{target}, never from the geometry. See \code{\link{net_to_sf}}.
}
}
\value{An object of class \code{netframe}: a list with the two tables
(\code{nodes}, \code{edges}), \code{directed}, \code{roles} (the column
names of each role, with \code{weight_type}), \code{declared}
(\code{loops}, \code{multiple}, and \code{between_modes_only} when
declared), \code{graph}, and the engine's edge index. Temporal and spatial
networks are subclasses: \code{c("temporal_netframe", "netframe")},
\code{c("spatial_netframe", "netframe")}, and
\code{c("spatiotemporal_netframe", "spatial_netframe", "temporal_netframe",
"netframe")}.
}
\references{
Butts, C. T. (2008). network: A package for managing relational data in R.
\emph{Journal of Statistical Software}, 24(2), 1-36.
\doi{10.18637/jss.v024.i02}

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
\code{\link{net_from}} to convert from other packages, \code{\link{net_read}}
for files, \code{\link{net_twomode}}, \code{\link{net_temporal}}, and
\code{\link{net_spatial}} for the other kinds of network, and
\code{\link{net_simplify}} for loops and parallel edges.
}
\examples{
edges <- data.frame(source = c("a", "b", "c"), target = c("b", "c", "a"))
g <- net_frame(edges)
g

# Weights must say what they mean.
w <- data.frame(source = c("a", "a", "b"), target = c("b", "c", "c"), n = c(3, 1, 2))
net_frame(w, directed = FALSE, weight = "n", weight_type = "strength")

# Every problem is reported at once, with the rows.
nodes <- data.frame(id = c("a", "b", "b"))
try(net_frame(data.frame(source = c("a", NA), target = c("q", "b")), nodes))

# Self-loops and parallel edges are kept and reported.
h <- net_frame(data.frame(source = c("a", "a", "a"), target = c("a", "b", "b")))
net_diagnostics(h)
net_is_simple(h)
}
