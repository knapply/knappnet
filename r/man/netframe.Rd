\name{netframe}
\alias{netframe}
\alias{print.netframe}
\title{Create a network from node and edge data frames}
\usage{
netframe(edges, nodes = NULL, directed = TRUE, from = 1L, to = 2L,
         id = 1L, weight = NULL, weight_type = c("distance", "strength"),
         graph = list())

\method{print}{netframe}(x, ...)
}
\arguments{
\item{edges}{A data frame of edges.}
\item{nodes}{A data frame of nodes, or \code{NULL} to derive the node set from
the edge endpoints in order of first appearance.}
\item{directed}{Whether edges are directed: \code{TRUE} or \code{FALSE}.}
\item{from, to}{Columns of \code{edges} holding the endpoints, by name or
position.}
\item{id}{Column of \code{nodes} holding node identifiers, by name or
position.}
\item{weight}{Optional column of \code{edges} holding edge weights. Without
one, the network is binary: an edge either exists or it does not, and
algorithms use their faster unweighted forms.}
\item{weight_type}{What the weights mean. \code{"distance"} (the default, and
the convention of igraph and NetworkX) means larger weights keep nodes
farther apart, as with travel time, cost, or latency. \code{"strength"} means
larger weights bring nodes closer, as with contact frequency or tie
intensity; path-based measures then use 1/weight as the length (Newman 2001;
Opsahl, Agneessens, and Skvoretz 2010). Measures that add up weights, such as
\code{\link{strength}}, use them as given either way.}
\item{graph}{A named list of network-level attributes (for example a title or
data source), carried through conversions and file formats that support
them.}
\item{x}{A netframe.}
\item{...}{Ignored.}
}
\value{An object of class \code{netframe}. Printing it gives the size of
the network and whether it is weighted, and adds a line for self-loops and
parallel edges and for missing or non-finite weights when there are any.}
\description{
A network is a pair of data frames: \code{nodes}, with one row per node, and
\code{edges}, with one row per edge. Any other columns are attributes and are
carried along untouched. Column roles say which columns hold node
identifiers, edge endpoints, and (optionally) edge weights.
}
\section{Missing values and invalid input}{
\itemize{
\item Edge endpoints and node identifiers must not be missing (\code{NA}) or
empty strings; the error names the rows. Identifiers may be character,
factor (matched by their labels), integer, double, or logical, or a classed
vector such as \code{Date}; \code{integer64} identifiers (bit64) are matched
as text. Without \code{nodes}, classed \code{from} and \code{to} columns must
have the same class.
\item Node identifiers must be unique, and every endpoint must appear in
\code{nodes}.
\item Networks without edges, and isolated nodes, are allowed.
\item Self-loops and parallel edges are allowed and kept, and \code{print()}
reports how many there are. Each measure documents how it counts them;
\code{convention} selects the treatment of other tools (see
\code{\link{conventions}}).
\item Weights must be numeric. Missing, \code{NaN}, and infinite weights are
accepted here, so the data can still be stored, converted, and written to
files, and \code{print()} reports them; but every measure stops with an error
naming the rows, even measures that ignore weights, until they are removed
or replaced. Zero and negative weights are accepted and checked by each
measure that uses weights: path-based measures (betweenness, closeness,
harmonic centrality) and structural holes need them strictly positive;
modularity, Louvain, and eigenvector centrality need them non-negative;
strength, structural equivalence, and CONCOR take any finite value. Errors
name the row of \code{x$edges} that holds the offending weight.
\item Arguments are checked before any computation: \code{x} must be a
netframe, logical options \code{TRUE} or \code{FALSE}, \code{threads} a
non-negative whole number, and choices one of the listed values. Errors name
the argument.
\item The edge index is computed once, here. Add columns to \code{x$nodes}
and \code{x$edges} freely, but to add, remove, or reorder rows build a new
netframe. Measures stop when the number of edges has changed or an edge
points past the last node; other changes to rows, such as reordering,
cannot be detected and would misalign results.
}
}

\section{Spatial tables}{
\code{nodes} and \code{edges} may be sf data frames (simple features). Their
geometry columns stay in place and are recorded in \code{x$roles}
(\code{node_geometry}, \code{edge_geometry}, \code{crs}); both tables must
then have the same coordinate reference system. Connectivity still comes
from \code{from} and \code{to}, never from the geometry. See
\code{\link{to_sf}} for what the geometries are used for.
}

\section{Threads}{
Measures that run in parallel take \code{threads}. \code{0} (the default)
uses all available cores, except that only 2 are used when the environment
variable \env{_R_CHECK_LIMIT_CORES_} is set (as it is by
\verb{R CMD check --as-cran}), and at most \env{OMP_THREAD_LIMIT} when that is
set. Results are reproducible for a given number of threads: integer results
never depend on it, and floating-point results can differ in the last bits
between thread counts. Computations run in compiled code and cannot be
interrupted until they return.
}

\references{
Newman, M. E. J. (2001). Scientific collaboration networks. II. Shortest
paths, weighted networks, and centrality. \emph{Physical Review E}, 64,
016132. \doi{10.1103/PhysRevE.64.016132}

Opsahl, T., Agneessens, F., & Skvoretz, J. (2010). Node centrality in
weighted networks: Generalizing degree and shortest paths. \emph{Social
Networks}, 32(3), 245-251. \doi{10.1016/j.socnet.2010.03.006}
}
\examples{
edges <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"))
g <- netframe(edges)
g

# Self-loops, parallel edges, and missing weights are reported.
h <- netframe(data.frame(from = c("a", "a", "a", "b"), to = c("a", "b", "b", "c"),
                         w = c(1, 2, NA, 1)), weight = "w")
h
}
