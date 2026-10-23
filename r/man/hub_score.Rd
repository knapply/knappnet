\name{hub_score}
\alias{hub_score}
\alias{authority_score}
\title{Hub and authority scores}
\description{
Kleinberg's (1999) hubs and authorities: a good hub points to good
authorities, and a good authority is pointed to by good hubs. Authority
scores are the leading eigenvector of \eqn{A^T A}{t(A) A} and hub scores that
of \eqn{A A^T}{A t(A)}.
}
\usage{
hub_score(x, weighted = NULL, scale = NULL, convention = c("igraph", "networkx"),
          tol = 1e-12, max_iter = 10000L, threads = 0L)

authority_score(x, weighted = NULL, scale = NULL,
                convention = c("igraph", "networkx"), tol = 1e-12,
                max_iter = 10000L, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{weighted}{Whether to use edge weights, as given whatever the network's
\code{weight_type}. \code{NULL} (the default) uses them when the network has a
weight column. They must not be negative.}
\item{scale}{\code{"max"} (largest score 1, igraph), \code{"sum"} (scores sum
to 1, NetworkX), or \code{"unit"} (Euclidean length 1). \code{NULL} uses the
convention's choice.}
\item{convention}{\code{"igraph"} or \code{"networkx"}: the scaling, and how
an undirected self-loop counts (twice or once on the diagonal).}
\item{tol}{Convergence tolerance (relative residual of the eigenvector).}
\item{max_iter}{Cap on matrix-vector products. Not converging within it is an
error.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\details{
The authority vector is computed with the Lanczos method on \eqn{A^T A}{t(A)
A}, which is never formed, and the hub vector is \eqn{A a}{A a}. Parallel
edges add up; weights are used as given and must not be negative.

The scores are unique only when the largest singular value of \eqn{A} is
simple. When it is not (for example in bipartite or disconnected networks, or
a directed cycle), a warning is given and the result is the limit of
Kleinberg's iteration \eqn{a \leftarrow A^T h}{a <- t(A) h}, \eqn{h \leftarrow
A a}{h <- A a} from equal hub scores, which is also what NetworkX's
\code{hits()} returns; igraph's answer then depends on its starting vector.
Without ties every score is \code{NaN}.
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\references{
Kleinberg, J. M. (1999). Authoritative sources in a hyperlinked environment.
\emph{Journal of the ACM}, 46(5), 604-632. \doi{10.1145/324133.324140}
}
\examples{
d <- netframe(data.frame(from = c(1, 1, 2, 3, 4), to = c(3, 4, 3, 4, 3)))
hub_score(d)                          # igraph: largest score 1
hub_score(d, convention = "networkx") # NetworkX: scores sum to 1
hub_score(d, scale = "unit")          # Euclidean length 1

authority_score(d)                          # igraph
authority_score(d, convention = "networkx") # NetworkX
}
