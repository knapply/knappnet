\name{betweenness}
\alias{betweenness}
\title{Betweenness centrality}
\usage{
betweenness(x, weighted = FALSE, normalized = FALSE, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{weighted}{Whether to use edge weights as distances.}
\item{normalized}{Whether to normalize by the number of pairs.}
\item{threads}{Number of threads; \code{0} uses all available cores.}
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\description{
Shortest-path betweenness (Freeman 1977), computed with the algorithm of
Brandes (2001) in parallel.
}
\details{
Endpoints are excluded and, in undirected networks, each pair of nodes is
counted once. Parallel edges count as distinct shortest paths. With
\code{weighted = TRUE}, weights are treated as distances and must be
positive. Normalized values divide by the number of pairs not involving the
node: \eqn{(n-1)(n-2)} for directed and \eqn{(n-1)(n-2)/2} for undirected
networks. These conventions match igraph and NetworkX.
}
\references{
Freeman, L. C. (1977). A set of measures of centrality based on betweenness.
\emph{Sociometry}, 40(1), 35-41. \doi{10.2307/3033543}

Brandes, U. (2001). A faster algorithm for betweenness centrality.
\emph{Journal of Mathematical Sociology}, 25(2), 163-177.
\doi{10.1080/0022250X.2001.9990249}
}
\examples{
g <- netframe(data.frame(from = 1:3, to = 2:4), directed = FALSE)
betweenness(g)
}
