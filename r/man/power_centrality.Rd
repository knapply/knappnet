\name{power_centrality}
\alias{power_centrality}
\title{Bonacich power centrality}
\description{
Bonacich's (1987) family of measures \eqn{c(\alpha, \beta) = \alpha (I -
\beta A)^{-1} A 1}{c(alpha, beta) = alpha (I - beta A)^-1 A 1}: a weighted
count of the walks starting at each node, with a walk of length k weighted by
\eqn{\beta^{k-1}}.
}
\usage{
power_centrality(x, exponent = 1, mode = c("out", "in", "all"), weighted = NULL,
                 loops = FALSE, rescale = FALSE,
                 convention = c("igraph", "sna"), tol = 1e-12,
                 max_iter = 10000L, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{exponent}{The attenuation \eqn{\beta}. Positive values make a node
powerful when its contacts are well connected; negative values when its
contacts are weak (Bonacich 1987).}
\item{mode}{For directed networks: \code{"out"} counts walks along out-ties
(Bonacich, igraph, sna), \code{"in"} along in-ties, \code{"all"} ignores
direction.}
\item{weighted}{Whether to use edge weights, as given whatever the network's
\code{weight_type}. \code{NULL} uses them when the network has a weight
column. They must not be negative.}
\item{loops}{Whether self-loops enter the diagonal of the adjacency matrix
(the number of loops, or their weights). igraph and sna ignore the diagonal
by default.}
\item{rescale}{Whether to divide the scores by their sum instead of scaling
them so that their squares sum to the number of nodes.}
\item{convention}{\code{"igraph"} counts every edge, so parallel edges add
up; \code{"sna"} reads the network as an adjacency matrix, so parallel edges
collapse into one tie (keeping the largest weight).}
\item{tol}{Convergence tolerance of the iterative solver, relative to the
largest score.}
\item{max_iter}{Iteration cap. Not converging within it is an error.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\details{
\eqn{\alpha} is chosen so that the squared scores sum to the number of nodes,
as in igraph's \code{power_centrality()} and sna's \code{bonpow()}. When
\eqn{|\beta| \lambda_{max} < 1} (checked with Collatz's bounds) the walk
series converges and the scores are found by iterating \eqn{y \leftarrow A 1 +
\beta A y}{y <- A 1 + beta A y}. Otherwise, as with the common default
\eqn{\beta = 1}, the series diverges, but the linear system usually still has
a solution; it is then solved by Gaussian elimination with partial pivoting,
which is limited to 5,000 nodes. If \eqn{1/\beta} is an eigenvalue of
\eqn{A} the system is singular and an error is reported. Without ties every
score is \code{NaN} (\eqn{\alpha} would be 0/0).

In directed networks the default counts walks along out-ties, so a node
without out-ties scores 0 (Bonacich 1987). Weights are used as given.
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\references{
Bonacich, P. (1987). Power and centrality: A family of measures.
\emph{American Journal of Sociology}, 92(5), 1170-1182. \doi{10.1086/228631}

Golub, G. H., & Van Loan, C. F. (2013). \emph{Matrix Computations} (4th
ed.). Johns Hopkins University Press.
}
\examples{
# A star around node 1 with a pendant at node 5, and the tie 1-2 recorded
# twice: igraph counts both edges, sna reads one tie.
g <- netframe(data.frame(from = c(1, 1, 1, 1, 2), to = c(2, 2, 3, 4, 5)),
              directed = FALSE)
power_centrality(g, exponent = 0.3)                     # igraph
power_centrality(g, exponent = 0.3, convention = "sna") # sna

# A negative exponent rewards ties to weakly connected nodes.
power_centrality(g, exponent = -0.3)
}
