\name{betweenness}
\alias{betweenness}
\title{Betweenness centrality}
\usage{
betweenness(x, weighted = NULL, normalized = FALSE, threads = 0L,
            weights_as = NULL)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{weighted}{Whether to use edge weights. \code{NULL} (the default) uses
them when the network has a weight column.}
\item{normalized}{Whether to normalize by the number of pairs.}
\item{threads}{Number of threads; \code{0} uses all available cores.}
\item{weights_as}{Override the network's \code{weight_type}:
\code{"distance"} or \code{"strength"}.}
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\description{
Shortest-path betweenness (Freeman 1977), computed with the algorithm of
Brandes (2001) in parallel.
}
\details{
Endpoints are excluded and, in undirected networks, each pair of nodes is
counted once. Parallel edges count as distinct shortest paths.

By default, weighted networks use their weights and binary networks do not.
Weights are read as the network declares them (see \code{weight_type} in
\code{\link{netframe}}): distances as given, or strengths converted to
lengths 1/weight. They must be positive. The fastest exact algorithm is
chosen automatically: breadth-first search for binary networks and for
weights that are all equal, a bucket queue with exact integer arithmetic for
small whole-number lengths (Dial 1969), and Dijkstra's algorithm (Dijkstra
1959) otherwise.

Normalized values divide by the number of pairs not involving the
node: \eqn{(n-1)(n-2)} for directed and \eqn{(n-1)(n-2)/2} for undirected
networks. These conventions match igraph and NetworkX.
}
\references{
Freeman, L. C. (1977). A set of measures of centrality based on betweenness.
\emph{Sociometry}, 40(1), 35-41. \doi{10.2307/3033543}

Brandes, U. (2001). A faster algorithm for betweenness centrality.
\emph{Journal of Mathematical Sociology}, 25(2), 163-177.
\doi{10.1080/0022250X.2001.9990249}

Dial, R. B. (1969). Algorithm 360: Shortest-path forest with topological
ordering. \emph{Communications of the ACM}, 12(11), 632-633.
\doi{10.1145/363269.363610}

Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
\emph{Numerische Mathematik}, 1, 269-271. \doi{10.1007/BF01386390}
}
\examples{
g <- netframe(data.frame(from = 1:3, to = 2:4), directed = FALSE)
betweenness(g)
}
