\name{graph_density}
\alias{graph_density}
\title{Edge density}
\description{
The number of edges divided by the number of possible edges:
\eqn{n(n-1)} in directed and \eqn{n(n-1)/2} in undirected networks, or
\eqn{n^2} and \eqn{n(n+1)/2} when self-loops count as possible edges
(Wasserman and Faust 1994).
}
\usage{
graph_density(x, loops = FALSE, convention = c("igraph", "networkx", "sna", "ucinet"))
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{loops}{Whether self-loops count as possible edges.}
\item{convention}{How edges are counted; see \code{\link{conventions}}.
\code{"igraph"} and \code{"networkx"} count every edge, so parallel edges and
loops can push the density above 1; \code{"sna"} and \code{"ucinet"} count
distinct ties between distinct nodes (the adjacency matrix without its
diagonal).}
}
\details{
The name avoids masking \code{stats::density()}, the kernel density
estimator, which many scripts use alongside network functions.

With the default convention every edge counts, including parallel edges and
self-loops; this is igraph's \code{edge_density()} and, for networks with
more than one node, NetworkX's \code{density()}. With \code{convention =
"sna"} or \code{"ucinet"}, parallel edges count once and self-loops are
ignored, as in sna's \code{gden(diag = FALSE)} and UCINET's overall density
with the diagonal excluded. Weights are not used.

A network with one node has no possible edges between distinct nodes, so
its density is \code{NaN} (as in igraph and sna), except under the NetworkX
convention, which returns 0.
}
\value{A number.}
\references{
Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods and
Applications}. Cambridge University Press. \doi{10.1017/CBO9780511815478}
}
\examples{
g <- netframe(data.frame(from = c(1, 2, 3), to = c(2, 3, 4)), directed = FALSE)
graph_density(g)

# A repeated edge and a self-loop: igraph and NetworkX count every edge,
# while sna and UCINET count distinct ties between different nodes.
h <- netframe(data.frame(from = c(1, 1, 2, 3), to = c(2, 2, 3, 3)), directed = FALSE)
graph_density(h)                          # igraph
graph_density(h, convention = "networkx") # NetworkX
graph_density(h, convention = "sna")      # sna
graph_density(h, convention = "ucinet")   # UCINET

# A single node: NetworkX reports 0, the others NaN.
one <- netframe(data.frame(from = integer(0), to = integer(0)), nodes = data.frame(id = 1))
graph_density(one)
graph_density(one, convention = "networkx")
}
