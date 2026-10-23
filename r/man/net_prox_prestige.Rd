\name{net_prox_prestige}
\alias{net_prox_prestige}
\title{Proximity prestige}
\description{
The prestige of a node in a directed network from the nodes that can reach
it, its influence domain (Lin 1976; Wasserman and Faust 1994, section
5.3.3). With \eqn{I_i} the number of nodes that reach node \eqn{i} (other
than itself), \eqn{n} the number of nodes, and \eqn{d(j, i)} the
shortest-path length from \eqn{j} to \eqn{i},
\deqn{P_P(i) = \frac{I_i / (n - 1)}{\sum_{j \in I_i} d(j, i) / I_i},}{P_P(i) = (I_i / (n - 1)) / (sum_j d(j, i) / I_i),}
the share of nodes that reach \eqn{i} divided by their mean distance from
it. It is 0 for a node that no other node reaches and 1 for a node every
other node points to directly.
}
\usage{
net_prox_prestige(x, weighted = NULL, weights_as = NULL, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{weighted}{Whether to use edge weights. \code{NULL} (the default) uses
them when the network has a weight column.}
\item{weights_as}{Override the network's \code{weight_type}: \code{"distance"}
or \code{"strength"}.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{net_frame}}).}
}
\details{
Lengths are handled as in \code{\link{net_closeness}}: without weights every
edge has length 1; with weights they are distances, or strengths converted
to lengths 1/weight, and must be positive. Self-loops never lie on a
shortest path and parallel edges do not change path lengths. In an
undirected network every node that reaches \eqn{i} is also reached by it,
so the measure reduces to the Wasserman and Faust closeness index of
\code{net_closeness(x, convention = "networkx")}; for a directed network it is
that index computed over paths \emph{to} each node, which is NetworkX's
\code{closeness_centrality()} on a directed graph. Neither igraph nor sna
has a function for it.
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\references{
Lin, N. (1976). \emph{Foundations of Social Research}. McGraw-Hill.

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods
and Applications}. Cambridge University Press.
\doi{10.1017/CBO9780511815478}
}
\seealso{\code{\link{net_closeness}}, \code{\link{net_reach}}.}
\examples{
# 1 -> 2 -> 3, and 4 -> 3: node 3 is reached by everyone, 1 by no one.
g <- net_frame(data.frame(source = c(1, 2, 4), target = c(2, 3, 3)))
net_prox_prestige(g)
}
