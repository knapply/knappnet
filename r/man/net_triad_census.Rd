\name{net_triad_census}
\alias{net_triad_census}
\title{Triad census}
\description{
Counts every set of three nodes by its type among the 16 directed triad types
of Holland and Leinhardt (1976), using the algorithm of Batagelj and Mrvar
(2001).
}
\usage{
net_triad_census(x, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{net_frame}}).}
}
\details{
Self-loops are ignored and parallel edges are counted once. In an undirected
network every edge is treated as mutual, as in igraph. Weights are ignored:
the census is defined on binary ties.
}
\value{A named numeric vector of 16 counts in the standard order.}
\references{
Holland, P. W., & Leinhardt, S. (1976). Local structure in social networks.
\emph{Sociological Methodology}, 7, 1-45. \doi{10.2307/270703}

Batagelj, V., & Mrvar, A. (2001). A subquadratic triad census algorithm for
large sparse networks with small maximum degree. \emph{Social Networks},
23(3), 237-243. \doi{10.1016/S0378-8733(01)00035-1}
}
\examples{
g <- net_frame(data.frame(source = c(1, 2, 3), target = c(2, 3, 1)))
net_triad_census(g)
}
