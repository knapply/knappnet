\name{twomode}
\alias{twomode}
\alias{as_twomode}
\alias{is_twomode}
\alias{incidence_matrix}
\title{Two-mode networks}
\usage{
twomode(edges, nodes = NULL, type = "type", first = NULL, directed = FALSE,
  from = 1L, to = 2L, id = 1L, weight = NULL,
  weight_type = c("distance", "strength"), graph = list())

as_twomode(x, type = NULL, first = NULL, weighted = NULL)

is_twomode(x)

incidence_matrix(x, attr = NULL, sparse = FALSE)
}
\arguments{
\item{edges}{A data frame of edges.}
\item{nodes}{A data frame of nodes, or \code{NULL} to take the first mode
from the \code{from} column and the second from the \code{to} column.}
\item{type}{For \code{twomode()} and \code{as_twomode()}, the node column
holding the type, by name. When \code{twomode()} builds the node table
itself, this is the name of the logical column it creates.}
\item{first}{The value of \code{type} that marks the first mode; see
Details.}
\item{directed}{Whether edges are directed. Two-mode measures ignore
direction.}
\item{from, to}{Columns of \code{edges} holding the endpoints, by name or
position.}
\item{id}{Column of \code{nodes} holding node identifiers, by name or
position.}
\item{weight}{Optional column of \code{edges} holding edge weights.}
\item{weight_type}{What the weights mean: \code{"distance"} or
\code{"strength"}; see \code{\link{netframe}}.}
\item{graph}{A named list of network-level attributes.}
\item{x}{For \code{as_twomode()}, the object to convert; for the other
functions, a two-mode \code{\link{netframe}}.}
\item{weighted}{For matrices, whether values are edge weights. \code{NULL}
treats any value other than 0 and 1 as a weight.}
\item{attr}{For \code{incidence_matrix()}, an edge column whose values fill
the matrix; \code{NULL} counts edges.}
\item{sparse}{For \code{incidence_matrix()}, whether to return a sparse
\code{Matrix}.}
}
\value{
\code{twomode()} and \code{as_twomode()} return a \code{\link{netframe}}
with a node type role. \code{is_twomode()} returns \code{TRUE} or
\code{FALSE}. \code{incidence_matrix()} returns a matrix with first-mode
nodes as rows and second-mode nodes as columns; parallel edges add up.
}
\description{
A two-mode (bipartite, affiliation) network has two sets of nodes, such as
people and the events they attend, or people and the organizations they
belong to, with ties only between the sets. \code{twomode()} creates one
from data frames, \code{as_twomode()} converts a netframe, an igraph or
network object, an incidence matrix, or a network read from a file, and
\code{is_twomode()} tests for one. \code{incidence_matrix()} returns the
rectangular matrix with one row per node of the first mode and one column
per node of the second.
}
\details{
A two-mode network is a \code{\link{netframe}} whose node table has a
logical column, recorded as the node type role (\code{x$roles$type}):
\code{FALSE} for the first mode and \code{TRUE} for the second, igraph's
convention. Every edge must join the two modes; edges within a mode are an
error.

The type may be given as a logical column, or as a column with exactly two
distinct values (for example \code{"person"} and \code{"event"}, the two
node classes of a DyNetML file, or the 1 and 2 of a bipartite network from
the network package). \code{first} names the value of the first mode; by
default it is the first factor level, or the smallest value (text is
compared byte by byte, so the choice does not depend on the locale). A
logical column is then added: it is called \code{type} (so that
\code{\link{to_igraph}} gives igraph's bipartite functions the \code{type}
attribute they look for) unless the node table already has a column with
that name, in which case it is called \code{is_} followed by the second
mode's value. The original column is kept, and the two values are recorded
in \code{x$roles$type_labels}.

Without a node table, \code{twomode()} takes the first mode from the
\code{from} column of the edges and the second from the \code{to} column,
so a list of (person, event) pairs is all it needs. A value appearing in
both columns is an error.

\code{as_twomode()} finds the type column itself when \code{type} is
\code{NULL}: a logical \code{type} column, else a \code{class} column with
exactly two values (as \code{\link{read_network}} gives for a DyNetML file
with two node classes), else a two-valued \code{type} or \code{mode}
column. A matrix (or sparse \code{Matrix}) is read as an incidence matrix:
rows are the first mode and columns the second, with identifiers from the
row and column names (\code{r1}, \code{r2}, \ldots and \code{c1},
\code{c2}, \ldots if there are none); values other than 0 and 1 become edge
weights unless \code{weighted = FALSE}. A data frame is passed to
\code{twomode()}.

Two-mode measures (\code{\link{project}}, \code{\link{twomode_degree}} and
its relatives) read the network as an affiliation matrix and ignore the
direction of edges, so \code{twomode()} creates undirected networks by
default.
}
\references{
Breiger, R. L. (1974). The duality of persons and groups. \emph{Social
Forces}, 53(2), 181-190. \doi{10.2307/2576011}

Borgatti, S. P., & Everett, M. G. (1997). Network analysis of 2-mode data.
\emph{Social Networks}, 19(3), 243-269.
\doi{10.1016/S0378-8733(96)00301-2}

Borgatti, S. P., & Halgin, D. S. (2011). Analyzing affiliation networks. In
J. Scott & P. J. Carrington (Eds.), \emph{The SAGE Handbook of Social
Network Analysis} (pp. 417-433). SAGE. \doi{10.4135/9781446294413.n28}
}
\examples{
attended <- data.frame(
  person = c("ann", "ann", "bob", "bob", "cat", "dan"),
  event  = c("gala", "fair", "gala", "show", "show", "fair")
)
g <- twomode(attended)
g
g$nodes
incidence_matrix(g)
is_twomode(g)
}
