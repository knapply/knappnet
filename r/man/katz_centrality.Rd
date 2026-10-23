\name{katz_centrality}
\alias{katz_centrality}
\title{Katz centrality}
\usage{
katz_centrality(x, alpha = 0.1, beta = 1, mode = c("in", "out", "all"),
                weighted = NULL, normalized = NULL,
                convention = c("igraph", "networkx"), tol = 1e-12,
                max_iter = 10000L, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{alpha}{The attenuation factor. It must be below \eqn{1/\lambda_{max}},
the reciprocal of the largest eigenvalue of the adjacency matrix, for the
walk counts to converge; larger values are an error.}
\item{beta}{The exogenous term: one number for every node, a vector with one
value per node, or the name of a column of \code{x$nodes}.}
\item{mode}{For directed networks: \code{"in"} counts walks ending at the node
(NetworkX, igraph), \code{"out"} walks starting there, and \code{"all"}
ignores direction.}
\item{weighted}{Whether to use edge weights, as given whatever the network's
\code{weight_type}. \code{NULL} uses them when the network has a weight
column. They must not be negative.}
\item{normalized}{Whether to scale the scores to Euclidean length 1.
\code{NULL} uses the convention's default: \code{TRUE} for NetworkX,
\code{FALSE} for igraph.}
\item{convention}{\code{"igraph"} (\code{alpha_centrality()}) or
\code{"networkx"} (\code{katz_centrality()}); see Details.}
\item{tol}{Convergence tolerance of the iteration, relative to the largest
score.}
\item{max_iter}{Iteration cap. Not converging within it is an error.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\description{
Counts the walks that end at each node, with a walk of length k weighted by
\eqn{\alpha^k} (Katz 1953): the solution of
\eqn{x = \alpha A^T x + \beta}{x = alpha t(A) x + beta}.
}
\details{
Expanding the solution gives \eqn{x = \beta (1 + \alpha A^T 1 + \alpha^2
(A^T)^2 1 + \dots)}{x = beta (1 + alpha t(A) 1 + alpha^2 t(A)^2 1 + ...)}, a
weighted count of walks; this is NetworkX's \code{katz_centrality()} and
igraph's \code{alpha_centrality()} (Bonacich and Lloyd 2001). Katz's original
index leaves out the walks of length 0, which gives \eqn{x/\beta - 1} for
constant \eqn{\beta}. The series converges only for
\eqn{\alpha < 1/\lambda_{max}}; knappnet checks this with rigorous bounds on
\eqn{\lambda_{max}} (Collatz 1942) and stops with an error otherwise, where
NetworkX's direct solver returns numbers that are no longer walk counts.
Scores are found by iterating \eqn{x \leftarrow \alpha A^T x + \beta}{x <-
alpha t(A) x + beta}, in parallel on large networks.

Parallel edges add up. Under the igraph convention self-loops are ignored and
scores are not normalized; under NetworkX's a self-loop adds its weight once
to the diagonal and scores have unit length. Weights are used as given and
must not be negative.
}
\references{
Katz, L. (1953). A new status index derived from sociometric analysis.
\emph{Psychometrika}, 18(1), 39-43. \doi{10.1007/BF02289026}

Bonacich, P., & Lloyd, P. (2001). Eigenvector-like measures of centrality for
asymmetric relations. \emph{Social Networks}, 23(3), 191-201.
\doi{10.1016/S0378-8733(01)00038-7}

Collatz, L. (1942). Einschließungssatz für die charakteristischen Zahlen von
Matrizen. \emph{Mathematische Zeitschrift}, 48(1), 221-226.
\doi{10.1007/BF01180013}
}
\examples{
d <- netframe(data.frame(from = c(1, 1, 2, 3, 4), to = c(2, 3, 3, 1, 3)))
katz_centrality(d, alpha = 0.2)                          # igraph: walk counts
katz_centrality(d, alpha = 0.2, convention = "networkx") # NetworkX: unit length
}
