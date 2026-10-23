\name{as_netframe}
\alias{as_netframe}
\alias{as_netframe.default}
\alias{as_netframe.netframe}
\alias{as_netframe.data.frame}
\alias{as_netframe.igraph}
\alias{as_netframe.network}
\alias{as_netframe.matrix}
\alias{as_netframe.Matrix}
\alias{as_netframe.sf}
\alias{as_netframe.sfnetwork}
\alias{to_igraph}
\alias{to_network}
\alias{to_tbl_graph}
\alias{as.matrix.netframe}
\alias{as.igraph.netframe}
\alias{as.network.netframe}
\alias{as_tbl_graph.netframe}
\title{Convert networks to and from other packages}
\description{
\code{as_netframe()} turns network objects from other packages, adjacency
matrices, and edge lists into a \code{\link{netframe}}. \code{to_igraph()},
\code{to_network()}, and \code{to_tbl_graph()} convert back, and
\code{as.matrix()} gives an adjacency matrix. Node, edge, and network-level
attributes are carried across in both directions, so existing code and data
keep working.

The usual generics of other packages also work on a netframe:
\code{igraph::as.igraph()}, \code{network::as.network()}, and
\code{tidygraph::as_tbl_graph()}.
}
\usage{
as_netframe(x, ...)

\method{as_netframe}{default}(x, ...)

\method{as_netframe}{netframe}(x, ...)

\method{as_netframe}{data.frame}(x, ...)

\method{as_netframe}{igraph}(x, ...)

\method{as_netframe}{network}(x, ...)

\method{as_netframe}{matrix}(x, directed = NULL, weighted = NULL, ...)

\method{as_netframe}{Matrix}(x, directed = NULL, weighted = NULL, ...)

\method{as_netframe}{sf}(x, nodes = NULL, ..., edges = NULL)

\method{as_netframe}{sfnetwork}(x, ...)

to_igraph(x)

to_tbl_graph(x)

to_network(x)

\method{as.matrix}{netframe}(x, attr = NULL, sparse = FALSE, ...)
}
\arguments{
\item{x}{An object to convert.}
\item{nodes}{For sf edge tables, the node table (a data frame or sf data
frame), as in \code{\link{netframe}}.}
\item{edges}{For sf node tables, the edge table: when given, \code{x} is
taken as the node table instead.}
\item{...}{Passed on to methods. For data frames, arguments of
\code{\link{netframe}}.}
\item{directed}{For matrices, whether edges are directed. \code{NULL} infers
it from symmetry.}
\item{weighted}{For matrices, whether values are edge weights. \code{NULL}
treats any value other than 0 and 1 as a weight.}
\item{attr}{For \code{as.matrix()}, an edge column whose values fill the
matrix; \code{NULL} counts edges.}
\item{sparse}{For \code{as.matrix()}, whether to return a sparse
\code{Matrix}.}
}
\details{
\itemize{
\item \strong{igraph} (including \strong{tidygraph} \code{tbl_graph}
objects): vertex names become node identifiers (vertex positions are used
if the graph has no names). An edge attribute called \code{weight} becomes
the weight column.
\item \strong{network} (statnet): vertex names become node identifiers. For
bipartite networks a \code{mode} node column (1 or 2) records the two sets,
and converting back restores the bipartite structure.
\item \strong{matrix} or sparse \strong{Matrix}: row names become
identifiers. Values other than 0 and 1 become edge weights. Direction is
inferred from symmetry unless \code{directed} is given; an undirected matrix
contributes only its upper triangle. Missing values (\code{NA}) are read as 0,
no tie.
\item \strong{data.frame}: treated as an edge list, as in
\code{\link{netframe}}.
\item \strong{sf} data frames: an edge list whose geometry is kept (see
\code{\link{to_sf}}), or the node table when \code{edges} is given.
\item \strong{sfnetwork} (sfnetworks): see \code{\link{to_sfnetwork}}.
}
\code{as.matrix()} sums parallel edges. In undirected networks the matrix is
symmetric and a self-loop contributes once to the diagonal.
}
\value{
\code{as_netframe()} returns a \code{\link{netframe}}. The \code{to_*()}
functions return objects of the corresponding package, and
\code{as.matrix()} returns an adjacency matrix.
}
\examples{
m <- matrix(c(0, 2, 0,
              2, 0, 1,
              0, 1, 0), 3, byrow = TRUE,
            dimnames = list(letters[1:3], letters[1:3]))
g <- as_netframe(m)
g$edges
as.matrix(g, attr = "weight")
}
