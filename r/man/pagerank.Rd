\name{pagerank}
\alias{pagerank}
\title{PageRank}
\description{
The stationary distribution of a random walk that, at each step, follows one
of the current node's out-ties with probability \code{damping} and otherwise
jumps to a node chosen at random (Brin and Page 1998).
}
\usage{
pagerank(x, damping = 0.85, personalized = NULL, directed = TRUE, weighted = NULL,
         dangling = c("personalized", "uniform"),
         convention = c("igraph", "networkx"), tol = 1e-12, max_iter = 10000L,
         threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{damping}{Probability of following a tie rather than jumping, between 0
and 1.}
\item{personalized}{Jump probabilities: a numeric vector with one nonnegative
value per node (aligned with \code{x$nodes}), or the name of such a column of
\code{x$nodes}. It is rescaled to sum to 1. \code{NULL} jumps to every node
with equal probability.}
\item{directed}{Whether to follow the direction of ties in directed networks.
\code{FALSE} treats every tie as reciprocal.}
\item{weighted}{Whether to choose ties in proportion to their weights, used as
given whatever the network's \code{weight_type}. \code{NULL} (the default)
uses weights when the network has a weight column. They must not be negative.}
\item{dangling}{Where a walker at a node without out-ties (or whose out-ties
all have weight 0) jumps: \code{"personalized"} by the personalization
vector, as igraph and NetworkX do, or \code{"uniform"} to any node with equal
probability.}
\item{convention}{\code{"igraph"} or \code{"networkx"}; they differ only in
how an undirected self-loop counts (see Details).}
\item{tol}{Convergence tolerance: iteration stops when the scores, which sum
to 1, change by at most \code{tol} in total.}
\item{max_iter}{Iteration cap. Not converging within it is an error.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\details{
Scores are computed by power iteration from equal scores. The error after
stopping is at most \code{damping / (1 - damping)} times \code{tol}. Each step
runs in parallel on large networks.

Parallel edges add up, as their number or their weights. In undirected
networks a self-loop counts as two tie ends under the igraph convention (so a
walker at the node returns to it with probability 2/degree) and as one under
NetworkX's; in directed networks it counts once. The treatment of dangling
nodes follows Langville and Meyer (2004). sna and UCINET have no PageRank, so
those conventions are not available.

Results match igraph's \code{page_rank()} and NetworkX's \code{pagerank()} to
within their convergence tolerances.
}
\value{A numeric vector aligned with the rows of \code{x$nodes}, summing to
1.}
\references{
Brin, S., & Page, L. (1998). The anatomy of a large-scale hypertextual Web
search engine. \emph{Computer Networks and ISDN Systems}, 30(1-7), 107-117.
\doi{10.1016/S0169-7552(98)00110-X}

Langville, A. N., & Meyer, C. D. (2004). Deeper inside PageRank.
\emph{Internet Mathematics}, 1(3), 335-380.
\doi{10.1080/15427951.2004.10129091}
}
\examples{
# An undirected network with a self-loop at node 4: igraph counts the loop
# as two tie ends, NetworkX as one.
g <- netframe(data.frame(from = c(1, 2, 3, 1, 4), to = c(2, 3, 4, 3, 4)),
              directed = FALSE)
pagerank(g)                          # igraph
pagerank(g, convention = "networkx") # NetworkX

# Personalized PageRank: every random jump returns to node 1.
d <- netframe(data.frame(from = c(1, 1, 2, 3, 4), to = c(2, 3, 3, 1, 3)))
pagerank(d, personalized = c(1, 0, 0, 0))
}
