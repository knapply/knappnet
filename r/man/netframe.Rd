\name{netframe}
\alias{netframe}
\alias{print.netframe}
\title{Create a network from node and edge data frames}
\usage{
netframe(edges, nodes = NULL, directed = TRUE, from = 1L, to = 2L,
         id = 1L, weight = NULL, graph = list())

\method{print}{netframe}(x, ...)
}
\arguments{
\item{edges}{A data frame of edges.}
\item{nodes}{A data frame of nodes, or \code{NULL} to derive the node set from
the edge endpoints in order of first appearance.}
\item{directed}{Whether edges are directed.}
\item{from, to}{Columns of \code{edges} holding the endpoints, by name or
position.}
\item{id}{Column of \code{nodes} holding node identifiers, by name or
position.}
\item{weight}{Optional column of \code{edges} holding edge weights.}
\item{graph}{A named list of network-level attributes (for example a title or
data source), carried through conversions and file formats that support
them.}
\item{x}{A netframe.}
\item{...}{Ignored.}
}
\value{An object of class \code{netframe}.}
\description{
A network is a pair of data frames: \code{nodes}, with one row per node, and
\code{edges}, with one row per edge. Any other columns are attributes and are
carried along untouched. Column roles say which columns hold node
identifiers, edge endpoints, and (optionally) edge weights.
}
\examples{
edges <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"))
g <- netframe(edges)
g
}
