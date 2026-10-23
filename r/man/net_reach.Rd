\name{net_reach}
\alias{net_reach}
\alias{net_dw_reach}
\title{Reach centrality}
\description{
\code{net_reach()} counts the other nodes within \code{k} steps of each node
(UCINET's reach centrality, and the reach of Borgatti's (2006) key player
problem). \code{net_dw_reach()} is UCINET's average reciprocal
distance (ARD): the mean of \eqn{1/d} over the other nodes, with \eqn{1/d}
taken as 0 for nodes that cannot be reached, which is normalized harmonic
centrality (see \code{\link{net_harmonic}}).
}
\usage{
net_reach(x, k = 2, mode = c("out", "in", "all"), weighted = NULL, weights_as = NULL,
      normalized = FALSE, threads = 0L)

net_dw_reach(x, mode = c("out", "in", "all"), weighted = NULL,
                        weights_as = NULL, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{k}{Largest distance counted: a number of steps, or a path length when
weights are used. \code{Inf} counts every reachable node (see
\code{\link{net_reachability}}).}
\item{mode}{For directed networks: \code{"out"} counts the nodes a node
reaches, \code{"in"} the nodes that reach it, and \code{"all"} ignores
direction. Ignored for undirected networks.}
\item{weighted}{Whether to use edge weights. \code{NULL} (the default) uses
them when the network has a weight column.}
\item{weights_as}{Override the network's \code{weight_type}: \code{"distance"}
or \code{"strength"}.}
\item{normalized}{Whether to divide by \eqn{n - 1}, the largest possible
count (UCINET reports reach as this proportion).}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{net_frame}}).}
}
\details{
Without weights every edge has length 1, so \code{k} counts steps. With
weights (distances, or strengths converted to lengths 1/weight, as in
\code{\link{net_betweenness}}) \code{k} is a path length and a node counts when
its shortest-path length is at most \code{k}. Self-loops never lie on a
shortest path and parallel edges do not change path lengths. The node
itself is never counted. Searches from each node run in parallel.

\code{net_reach(x, k = 1)} is the number of distinct neighbours and
\code{net_reach(x, k = Inf)} is \code{\link{net_reachability}}; UCINET's 2-step
reach is the default \code{k = 2}. igraph's \code{net_ego_size(order = k)}
counts the node itself as well, and NetworkX has no reach function.

\code{net_dw_reach()} equals
\code{net_harmonic(x, mode, normalized = TRUE)}, the normalized
harmonic centrality of igraph, sna's \code{net_closeness(cmode = "suminvdir")},
and UCINET's reciprocal closeness.
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\references{
Borgatti, S. P. (2006). Identifying sets of key players in a social
network. \emph{Computational and Mathematical Organization Theory}, 12(1),
21-34. \doi{10.1007/s10588-006-7084-x}

Everton, S. F. (2012). \emph{Disrupting Dark Networks}. Cambridge
University Press. \doi{10.1017/CBO9781139136877}
}
\seealso{\code{\link{net_harmonic}}, \code{\link{net_closeness}},
\code{\link{net_reachability}}.}
\examples{
g <- net_example("florentine")
net_reach(g, k = 2)
net_reach(g, k = 2, normalized = TRUE)
net_dw_reach(g)

# In a directed network, out-reach counts the nodes a node reaches and
# in-reach the nodes that reach it.
d <- net_frame(data.frame(source = c(1, 2, 3), target = c(2, 3, 4)))
net_reach(d, k = 2, mode = "out")
net_reach(d, k = 2, mode = "in")
}
