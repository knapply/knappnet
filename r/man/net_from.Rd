\name{net_from}
\alias{net_from}
\alias{net_from.default}
\alias{net_from.netframe}
\alias{net_from.data.frame}
\alias{net_from.igraph}
\alias{net_from.network}
\alias{net_from.matrix}
\alias{net_from.Matrix}
\alias{net_from.sf}
\alias{net_from.sfnetwork}
\alias{net_to_igraph}
\alias{net_to_network}
\alias{net_to_tbl_graph}
\alias{as.matrix.netframe}
\alias{as.igraph.netframe}
\alias{as.network.netframe}
\alias{as_tbl_graph.netframe}
\title{Convert networks to and from other packages}
\description{
\code{net_from()} turns network objects from other packages, adjacency
matrices, and edge lists into a \code{\link[=net_frame]{netframe}}. \code{net_to_igraph()},
\code{net_to_network()}, and \code{net_to_tbl_graph()} convert back, and
\code{as.matrix()} gives an adjacency matrix. Node, edge, and network-level
attributes are carried across in both directions, so existing code and data
keep working.

The usual generics of other packages also work on a netframe:
\code{igraph::as.igraph()}, \code{network::as.network()}, and
\code{tidygraph::as_tbl_graph()}.
}
\usage{
net_from(x, ...)

\method{net_from}{default}(x, ...)

\method{net_from}{netframe}(x, ...)

\method{net_from}{data.frame}(x, ...)

\method{net_from}{igraph}(x, weight_type = "distance", quiet = FALSE, ...)

\method{net_from}{network}(x, weight_type = "strength", quiet = FALSE, ...)

\method{net_from}{matrix}(x, directed = NULL, weighted = NULL, weight_type = "strength",
         quiet = FALSE, ...)

\method{net_from}{Matrix}(x, directed = NULL, weighted = NULL, weight_type = "strength",
         quiet = FALSE, ...)

\method{net_from}{sf}(x, nodes = NULL, ..., edges = NULL)

\method{net_from}{sfnetwork}(x, weight_type = "distance", quiet = FALSE, ...)

net_to_igraph(x)

net_to_tbl_graph(x)

net_to_network(x)

\method{as.matrix}{netframe}(x, attr = NULL, sparse = FALSE, ...)
}
\arguments{
\item{x}{An object to convert.}
\item{nodes}{For sf edge tables, the node table (a data frame or sf data
frame), as in \code{\link{net_frame}}.}
\item{edges}{For sf node tables, the edge table: when given, \code{x} is
taken as the node table instead.}
\item{...}{Passed on to methods. For data frames, arguments of
\code{\link{net_frame}}.}
\item{directed}{For matrices, whether edges are directed. \code{NULL} infers
it from symmetry.}
\item{weighted}{For matrices, whether values are edge weights. \code{NULL}
treats any value other than 0 and 1 as a weight.}
\item{weight_type}{What the weights of the converted object mean (see
\code{\link{net_frame}}), when it has any. The defaults follow the source:
\code{"distance"} for igraph and sfnetworks objects, whose path functions
read weights as lengths, and \code{"strength"} for network objects and
valued matrices, which UCINET and sna read as tie strengths. Give the
other value when your data mean the opposite.}
\item{quiet}{Whether to suppress the advisory notes of validation.}
\item{attr}{For \code{as.matrix()}, an edge column whose values fill the
matrix; \code{NULL} counts edges.}
\item{sparse}{For \code{as.matrix()}, whether to return a sparse
\code{Matrix}.}
}
\details{
\itemize{
\item \strong{igraph} (including \strong{tidygraph} \code{tbl_graph}
objects): vertex names become node identifiers in a column called
\code{name} (vertex positions are used if the graph has no names). An
edge attribute called \code{weight} becomes the weight column, and a
logical vertex attribute \code{type} declares the modes of a two-mode
network (see \code{\link{net_twomode}}); converting back writes both.
\item \strong{network} (statnet): vertex names become node identifiers in
a column called \code{name}. For bipartite networks a \code{mode} node
column (\code{"first"} or \code{"second"}) declares the two modes, and
converting a two-mode netframe back restores the bipartite structure
when the first mode comes first in the node table.
\item \strong{matrix} or sparse \strong{Matrix}: row names become
identifiers in a column called \code{id}. Values other than 0 and 1
become edge weights. Direction is inferred from symmetry unless
\code{directed} is given; an undirected matrix contributes only its upper
triangle. Missing values (\code{NA}) are read as 0, no tie.
\item \strong{data.frame}: treated as an edge list, as in
\code{\link{net_frame}}.
\item \strong{sf} data frames: an edge list whose geometry is kept (see
\code{\link{net_to_sf}}), or the node table when \code{edges} is given.
\item \strong{sfnetwork} (sfnetworks): see \code{\link{net_to_sfnetwork}}.
}
\code{as.matrix()} sums parallel edges. In undirected networks the matrix is
symmetric and a self-loop contributes once to the diagonal.
}
\value{
\code{net_from()} returns a \code{\link[=net_frame]{netframe}}. The
\code{net_to_*()} functions return objects of the corresponding package,
and \code{as.matrix()} returns an adjacency matrix.
}
\examples{
m <- matrix(c(0, 2, 0,
              2, 0, 1,
              0, 1, 0), 3, byrow = TRUE,
            dimnames = list(letters[1:3], letters[1:3]))
g <- net_from(m)
g$edges
as.matrix(g, attr = "weight")
}
