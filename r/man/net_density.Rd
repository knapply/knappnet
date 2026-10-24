\name{net_density}
\alias{net_density}
\title{Edge density}
\description{
The number of edges divided by the number of possible edges:
\eqn{n(n-1)} in directed and \eqn{n(n-1)/2} in undirected networks, or
\eqn{n^2} and \eqn{n(n+1)/2} when self-loops count as possible edges
(Wasserman and Faust 1994).
}
\usage{
net_density(x, loops = FALSE, convention = c("igraph", "networkx", "sna", "ucinet"),
            relation = NULL)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{loops}{Whether self-loops count as possible edges.}
\item{relation}{For a multiplex network (see \code{\link{net_combine_relations}}),
the relation or relations whose edges count; \code{NULL} counts every edge.}
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
g <- net_frame(data.frame(source = c(1, 2, 3), target = c(2, 3, 4)), directed = FALSE)
net_density(g)

# A repeated edge and a self-loop: igraph and NetworkX count every edge,
# while sna and UCINET count distinct ties between different nodes.
h <- net_frame(data.frame(source = c(1, 1, 2, 3), target = c(2, 2, 3, 3)), directed = FALSE)
net_density(h)                          # igraph
net_density(h, convention = "networkx") # NetworkX
net_density(h, convention = "sna")      # sna
net_density(h, convention = "ucinet")   # UCINET

# A single node: NetworkX reports 0, the others NaN.
one <- net_frame(data.frame(source = integer(0), target = integer(0)), nodes = data.frame(id = 1))
net_density(one)
net_density(one, convention = "networkx")
}
