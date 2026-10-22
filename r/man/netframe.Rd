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
\item{directed}{Whether edges are directed.}
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
\value{An object of class \code{netframe}.}
\description{
A network is a pair of data frames: \code{nodes}, with one row per node, and
\code{edges}, with one row per edge. Any other columns are attributes and are
carried along untouched. Column roles say which columns hold node
identifiers, edge endpoints, and (optionally) edge weights.
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
}
