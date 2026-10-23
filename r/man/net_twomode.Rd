\name{net_twomode}
\alias{net_twomode}
\alias{net_as_twomode}
\alias{net_is_twomode}
\alias{net_is_bipartite}
\alias{net_incidence}
\title{Two-mode networks}
\description{
A two-mode (affiliation) network has two declared sets of nodes, such as
people and the events they attend, with ties only between the sets.
\code{net_twomode()} creates one from data frames, \code{net_as_twomode()}
converts a net_frame, an igraph or network object, an incidence matrix, or
a network read from a file, \code{net_is_twomode()} tests whether a
network declares two modes, \code{net_is_bipartite()} tests the structure
(whether the nodes can be 2-coloured), and \code{net_incidence()} returns
the rectangular matrix with one row per node of the first mode and one
column per node of the second.
}
\usage{
net_twomode(edges, nodes = NULL, mode = "mode", first = NULL, directed = FALSE,
            source = "source", target = "target", id = "id", weight = NULL,
            weight_type = NULL, quiet = FALSE, graph = list())

net_as_twomode(x, mode = NULL, first = NULL, weighted = NULL, quiet = FALSE)

net_is_twomode(x)

net_is_bipartite(x)

net_incidence(x, attr = NULL, sparse = FALSE)
}
\arguments{
\item{edges}{A data frame of edges.}
\item{nodes}{A data frame of nodes, or \code{NULL} to take the first mode
from the \code{source} column and the second from the \code{target}
column.}
\item{mode}{The name of the node column holding each node's mode. A
logical column called \code{type} (igraph's convention, \code{FALSE} for
the first mode) is accepted when there is no \code{mode} column. When
\code{net_twomode()} builds the node table itself, this is the name of the
column it creates, whose values are the names of the \code{source} and
\code{target} columns.}
\item{first}{The mode that counts as the first (the rows of the incidence
matrix, the \code{"first"} projection); by default the mode of the first
node in the table. For \code{net_as_twomode()} on a two-mode network,
\code{first} alone changes which mode is first.}
\item{directed}{Whether edges are directed. Two-mode measures ignore
direction.}
\item{source, target}{The names of the columns of \code{edges} holding the
endpoints.}
\item{id}{The name of the column of \code{nodes} holding node identifiers,
or of the column to create without \code{nodes}.}
\item{weight}{Optional name of a column of \code{edges} holding edge weights.}
\item{weight_type}{What the weights mean, \code{"strength"} or
\code{"distance"}, required with \code{weight}; see \code{\link{net_frame}}.}
\item{quiet}{Whether to suppress the advisory notes of validation.}
\item{graph}{A named list of network-level attributes.}
\item{x}{For \code{net_as_twomode()}, the object to convert; for the other
functions, a \code{\link[=net_frame]{net_frame}} (two-mode for
\code{net_incidence()}).}
\item{weighted}{For matrices, whether values are edge weights (read as
strengths). \code{NULL} treats any value other than 0 and 1 as a weight.}
\item{attr}{For \code{net_incidence()}, an edge column whose values fill
the matrix; \code{NULL} counts edges.}
\item{sparse}{For \code{net_incidence()}, whether to return a sparse
\code{Matrix}.}
}
\details{
A mode is a declared, disjoint set of nodes of one kind (Wasserman and
Faust 1994, ch. 8; Borgatti and Everett 1997). A two-mode network is a
\code{\link[=net_frame]{net_frame}} whose node table has a mode column
(role \code{x$roles$mode}, any discrete column: text, integer, logical, or
factor) and whose edges are declared to join different modes only
(\code{x$declared$between_modes_only}); an edge within a mode is an error
that names the rows. A column with more than two modes gives a k-partite
network, which the two-mode measures refuse. Modes declared with
\code{net_frame(mode = )} instead allow ties within a mode (a multimode or
multilevel network).

Being bipartite is a property of the structure, not a declaration: a
network is bipartite when its nodes can be split into two sets with every
edge between them (Newman 2018, ch. 6), which \code{net_is_bipartite()}
tests by 2-colouring each component (a self-loop makes a network
non-bipartite). Every two-mode network is bipartite, but a bipartite
network declares nothing: declare its modes with \code{net_as_twomode()}.

The modes are labelled by the values of the mode column, in order of
first appearance; \code{first} chooses which is the first mode.
\code{\link{net_to_igraph}} writes igraph's logical \code{type} attribute
(\code{FALSE} for the first mode) for every two-mode network, and
\code{\link{as_net_frame}} reads it back as the mode role.

Without a node table, \code{net_twomode()} takes the first mode from the
\code{source} column and the second from the \code{target} column, so a
list of (person, event) pairs is all it needs; the mode column it creates
holds the two column names. A value appearing in both columns is an error.

\code{net_as_twomode()} reads a matrix (or sparse \code{Matrix}) as an
incidence matrix: rows are the first mode (labelled \code{row}) and
columns the second (\code{column}), with identifiers from the row and
column names (\code{r1}, \code{r2}, \ldots and \code{c1}, \code{c2},
\ldots if there are none); values other than 0 and 1 become strengths
unless \code{weighted = FALSE}. A data frame is passed to
\code{net_twomode()}. A net_frame keeps its tables and declares the modes.

Two-mode measures (\code{\link{net_project}}, \code{\link{net_twomode_degree}}
and its relatives) read the network as an affiliation matrix and ignore
the direction of edges, so \code{net_twomode()} creates undirected networks
by default.
}
\value{
\code{net_twomode()} and \code{net_as_twomode()} return a net_frame with a
mode role and edges declared to join different modes only.
\code{net_is_twomode()} and \code{net_is_bipartite()} return \code{TRUE} or
\code{FALSE}. \code{net_incidence()} returns a matrix with first-mode
nodes as rows and second-mode nodes as columns; parallel edges add up.
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

Newman, M. E. J. (2018). \emph{Networks} (2nd ed.). Oxford University
Press. \doi{10.1093/oso/9780198805090.001.0001}

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods
and Applications}. Cambridge University Press.
\doi{10.1017/CBO9780511815478}
}
\examples{
attended <- data.frame(
  person = c("ann", "ann", "bob", "bob", "cat", "dan"),
  event  = c("gala", "fair", "gala", "show", "show", "fair")
)
g <- net_twomode(attended, source = "person", target = "event")
g
g$nodes
net_incidence(g)
net_is_twomode(g)

# A cycle of even length is bipartite but declares no modes.
ring <- net_frame(data.frame(source = 1:4, target = c(2:4, 1)), directed = FALSE)
net_is_bipartite(ring)
net_is_twomode(ring)
}
