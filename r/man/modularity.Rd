\name{modularity}
\alias{modularity}
\title{Modularity of a partition}
\usage{
modularity(x, membership, resolution = 1, weighted = NULL)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{membership}{Group of each node, aligned with the rows of
\code{x$nodes}. Any atomic vector (integer, character, factor); only
equality of values matters. Missing values are not allowed.}
\item{resolution}{Resolution parameter \eqn{\gamma} (non-negative).}
\item{weighted}{Whether to use edge weights. \code{NULL} (the default) uses
them when the network has a weight column.}
}
\value{A single number; \code{NaN} for a network without edges (or with
total weight zero).}
\description{
Modularity (Newman and Girvan 2004; Newman 2006) compares the weight of
edges inside groups with its expectation under a random network with the
same degrees.
}
\details{
For an undirected network with total edge weight \eqn{m} and node
strengths \eqn{k_i},
\deqn{Q = \frac{1}{2m} \sum_{ij} \left[A_{ij} - \gamma \frac{k_i k_j}{2m}\right] \delta(c_i, c_j),}{Q = 1/(2m) sum_ij [A_ij - gamma k_i k_j / (2m)] delta(c_i, c_j),}
and for a directed network (Leicht and Newman 2008)
\deqn{Q = \frac{1}{m} \sum_{ij} \left[A_{ij} - \gamma \frac{k_i^{out} k_j^{in}}{m}\right] \delta(c_i, c_j).}{Q = 1/m sum_ij [A_ij - gamma k_i^out k_j^in / m] delta(c_i, c_j).}
The resolution parameter \eqn{\gamma} follows Reichardt and Bornholdt
(2006); \eqn{\gamma = 1} gives classical modularity.

Weights (Newman 2004) are used as given, whatever the network's
\code{weight_type}, and must be non-negative. Parallel edges each count
and their weights add. In an undirected network a self-loop of weight
\eqn{w} adds \eqn{2w} to \eqn{A_{ii}} and to its node's strength; in a
directed network it adds \eqn{w} to \eqn{A_{ii}} and to out- and
in-strength.

These conventions match \code{igraph::modularity()} and NetworkX's
\code{community.modularity()} on multigraphs, except that NetworkX returns
0 for a network without edges, where modularity is undefined.
}
\references{
Newman, M. E. J., & Girvan, M. (2004). Finding and evaluating community
structure in networks. \emph{Physical Review E}, 69(2), 026113.
\doi{10.1103/PhysRevE.69.026113}

Newman, M. E. J. (2004). Analysis of weighted networks. \emph{Physical
Review E}, 70(5), 056131. \doi{10.1103/PhysRevE.70.056131}

Newman, M. E. J. (2006). Modularity and community structure in networks.
\emph{Proceedings of the National Academy of Sciences}, 103(23),
8577-8582. \doi{10.1073/pnas.0601602103}

Reichardt, J., & Bornholdt, S. (2006). Statistical mechanics of community
detection. \emph{Physical Review E}, 74(1), 016110.
\doi{10.1103/PhysRevE.74.016110}

Leicht, E. A., & Newman, M. E. J. (2008). Community structure in directed
networks. \emph{Physical Review Letters}, 100(11), 118703.
\doi{10.1103/PhysRevLett.100.118703}
}
\seealso{\code{\link{cluster_louvain}}}
\examples{
g <- netframe(data.frame(from = c(1, 1, 2, 3, 4, 4, 5), to = c(2, 3, 3, 4, 5, 6, 6)),
              directed = FALSE)
modularity(g, c(1, 1, 1, 2, 2, 2))
modularity(g, c("a", "a", "a", "b", "b", "b"), resolution = 0.5)
}
