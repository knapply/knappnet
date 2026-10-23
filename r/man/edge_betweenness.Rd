\name{edge_betweenness}
\alias{edge_betweenness}
\title{Edge betweenness}
\description{
The number of shortest paths that pass through each edge, summed over
pairs of nodes with each path weighted by the share of the pair's shortest
paths it represents (Girvan and Newman 2002, after Freeman 1977). Edges
with high betweenness join parts of the network that have few other
connections; removing them in order of betweenness is the Girvan-Newman
community method.
}
\usage{
edge_betweenness(x, weighted = NULL, normalized = FALSE, threads = 0L,
                 weights_as = NULL)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{weighted}{Whether to use edge weights. \code{NULL} (the default) uses
them when the network has a weight column.}
\item{normalized}{Whether to divide by the number of pairs of nodes:
\eqn{n(n-1)} in directed and \eqn{n(n-1)/2} in undirected networks.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
\item{weights_as}{Override the network's \code{weight_type}:
\code{"distance"} or \code{"strength"}.}
}
\details{
Values are computed with Brandes' (2008) accumulation for edges, one
search per node, run in parallel. In undirected networks each unordered
pair of nodes is counted once. Parallel edges are distinct edges and
distinct shortest paths, so two parallel edges each carry half the paths
one edge would; self-loops never lie on a shortest path and score 0.
Weights are lengths, as in \code{\link{betweenness}}: distances, or
strengths converted to 1/weight, and must be positive.

Unnormalized values match igraph's \code{edge_betweenness()} and NetworkX's
\code{edge_betweenness_centrality(normalized = False)}; normalized values
match NetworkX's default. sna has no edge betweenness function.
}
\value{A numeric vector aligned with the rows of \code{x$edges}.}
\references{
Freeman, L. C. (1977). A set of measures of centrality based on
betweenness. \emph{Sociometry}, 40(1), 35-41. \doi{10.2307/3033543}

Girvan, M., & Newman, M. E. J. (2002). Community structure in social and
biological networks. \emph{Proceedings of the National Academy of
Sciences}, 99(12), 7821-7826. \doi{10.1073/pnas.122653799}

Brandes, U. (2008). On variants of shortest-path betweenness centrality
and their generic computation. \emph{Social Networks}, 30(2), 136-145.
\doi{10.1016/j.socnet.2007.11.001}
}
\seealso{\code{\link{betweenness}}.}
\examples{
g <- example_network("florentine")
eb <- edge_betweenness(g)
g$edges[order(-eb), ][1:5, ]  # the five edges most shortest paths cross
}
