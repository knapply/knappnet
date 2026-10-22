\name{concor}
\alias{concor}
\title{CONCOR blockmodel}
\usage{
concor(x, depth = 1L, weighted = NULL, tolerance = 1e-7, max_iter = 100L,
       threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{depth}{Number of rounds of splitting (at least 1).}
\item{weighted}{Whether ties are total edge weights rather than 0 or 1.
\code{NULL} (the default) uses weights when the network has a weight
column.}
\item{tolerance}{Convergence tolerance, between 0 and 1.}
\item{max_iter}{Maximum number of iterations per split.}
\item{threads}{Number of threads; \code{0} uses all available cores.}
}
\value{An integer vector of block numbers (starting at 1) aligned with the
rows of \code{x$nodes}.}
\description{
Partitions nodes into positions with CONCOR (Breiger, Boorman, and Arabie
1975): the matrix of profile correlations is correlated with itself
repeatedly until its entries reach +1 or -1, which splits the nodes in two.
Each block is split again, \code{depth} times in all.
}
\details{
The starting correlations are those of
\code{structural_equivalence(x, "correlation")}: each node's row and column
of the adjacency matrix, compared over the other nodes, leaving out the
diagonal and the pair's own ties. Each block after the first split starts
from its submatrix of these correlations, so \code{depth} rounds give at
most \eqn{2^{depth}} blocks of nodes with non-constant profiles.

Nodes whose profile is constant (isolates, for example) cannot be
correlated with other nodes; they are set aside and given blocks of their
own, one per distinct profile.

Iteration stops when every entry has absolute value at least
\code{1 - tolerance}, when an iteration changes no entry by more than
\eqn{10^{-12}} (a fixed point other than +1 or -1, which some highly
symmetric networks have), or after \code{max_iter} iterations. Nodes are
then split by the sign of their correlation with the block's first node:
positive with it, zero or negative apart. A block whose correlations are all
positive stays whole. Blocks are numbered in depth-first order of the
splits, the side holding the lower-numbered node first, followed by the
blocks of constant profiles.

CONCOR can be numerically fragile. A warning is given when a split did not
converge, or when some node's correlation with the block's first node came
within \eqn{10^{-6}} of zero during the iteration: rounding error may then
decide that node's block, so another implementation or platform may give a
different partition.
}
\references{
Breiger, R. L., Boorman, S. A., & Arabie, P. (1975). An algorithm for
clustering relational data with applications to social network analysis and
comparison with multidimensional scaling. \emph{Journal of Mathematical
Psychology}, 12(3), 328-383. \doi{10.1016/0022-2496(75)90028-0}

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods
and Applications}. Cambridge University Press.
\doi{10.1017/CBO9780511815478}
}
\examples{
e <- expand.grid(from = 1:3, to = 4:6)
g <- netframe(e, directed = FALSE)
concor(g)
}
