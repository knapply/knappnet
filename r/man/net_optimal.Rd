\name{net_optimal}
\alias{net_optimal}
\title{Exact maximum modularity}
\description{
Finds a partition of maximum modularity exactly, by branch and bound, for
networks small enough to solve; it fails rather than return a partition
that might not be optimal. Finding the optimum is NP-hard (Brandes et al.
2008), so this is a benchmark for the heuristics of
\code{\link{net_louvain}} and \code{\link{net_leiden}} rather than a method
for large networks.
}
\usage{
net_optimal(x, resolution = 1, weighted = NULL)
}
\arguments{
\item{x}{A \code{\link{net_frame}}, directed or undirected, in which every
component has at most 60 nodes.}
\item{resolution}{Resolution parameter \eqn{\gamma} (non-negative); see
\code{\link{net_modularity}}.}
\item{weighted}{Whether to use edge weights, as given. \code{NULL} (the
default) uses them when the network has a weight column. They must not be
negative.}
}
\details{
Modularity is as \code{\link{net_modularity}} defines it, directed for
directed networks (igraph's \code{cluster_optimal()} accepts both).
Parallel edges add their weights and self-loops only add a constant. The
groups of an optimal partition are connected, so each component is solved
on its own, nodes are assigned in order of strongest attachment to those
already placed, and the suffixes of that order are solved in turn
(Russian doll search; Verfaillie, Lemaitre, and Schiex 1996) so that each
optimum bounds the next. The result is exact: an error is raised for a
component of more than 60 nodes, or after \eqn{10^9} search steps (a few
minutes), instead of a partition that may not be optimal. Running time
depends on structure more than size: the karate club and the Southern women
network take well under a second, a dense random network of 30 nodes about
40 seconds, and sparse random networks of 40 to 50 nodes with little
community structure may not finish.

When several partitions are optimal, which one is returned is unspecified.
Nodes without edges are groups of their own. A network without edges has
modularity \code{NaN} and one group per node. Unlike igraph's, which
solves the integer program of Brandes et al. with GLPK, the search needs no
solver; the validation suite compares the optimum with igraph's on every
network where both finish.
}
\value{An integer vector of group labels aligned with the rows of
\code{x$nodes}, numbered 1, 2, ... in order of first appearance, with the
maximum modularity in attribute \code{"modularity"}.}
\references{
Brandes, U., Delling, D., Gaertler, M., Goerke, R., Hoefer, M., Nikoloski,
Z., & Wagner, D. (2008). On modularity clustering. \emph{IEEE Transactions
on Knowledge and Data Engineering}, 20(2), 172--188.
\doi{10.1109/TKDE.2007.190689}

Verfaillie, G., Lemaitre, M., & Schiex, T. (1996). Russian doll search for
solving constraint optimization problems. In \emph{Proceedings of the
Thirteenth National Conference on Artificial Intelligence (AAAI-96)},
181--187. AAAI Press.
}
\seealso{\code{\link{net_louvain}}, \code{\link{net_leiden}},
\code{\link{net_modularity}}}
\examples{
g <- net_example("florentine")
best <- net_optimal(g)
best
net_plot_partition(g, best, seed = 1)

# How close is a heuristic?
attr(net_louvain(g, seed = 1), "modularity") / attr(best, "modularity")
}
