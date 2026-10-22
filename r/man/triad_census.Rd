\name{triad_census}
\alias{triad_census}
\title{Triad census}
\usage{
triad_census(x, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{threads}{Number of threads; \code{0} uses all available cores.}
}
\value{A named numeric vector of 16 counts in the standard order.}
\description{
Counts every set of three nodes by its type among the 16 directed triad types
of Holland and Leinhardt (1976), using the algorithm of Batagelj and Mrvar
(2001).
}
\details{
Self-loops are ignored and parallel edges are counted once. In an undirected
network every edge is treated as mutual, as in igraph. Weights are ignored:
the census is defined on binary ties.
}
\references{
Holland, P. W., & Leinhardt, S. (1976). Local structure in social networks.
\emph{Sociological Methodology}, 7, 1-45. \doi{10.2307/270703}

Batagelj, V., & Mrvar, A. (2001). A subquadratic triad census algorithm for
large sparse networks with small maximum degree. \emph{Social Networks},
23(3), 237-243. \doi{10.1016/S0378-8733(01)00035-1}
}
\examples{
g <- netframe(data.frame(from = c(1, 2, 3), to = c(2, 3, 1)))
triad_census(g)
}
