\name{blockmodel}
\alias{blockmodel}
\alias{block_summary}
\alias{print.blockmodel}
\title{Blockmodel of a partition}
\description{
Reduces a network to its positions: given a partition of the nodes into
blocks, \code{blockmodel()} computes the density of ties within and
between every pair of blocks, the binary image matrix that codes each
block as present or absent at a cutoff, and how well the image fits the
network (White, Boorman, and Breiger 1976; Wasserman and Faust 1994, ch.
10). \code{block_summary()} lists the same quantities as one row per pair
of blocks.
}
\usage{
blockmodel(x, membership, weighted = NULL, loops = FALSE, cutoff = NULL)

block_summary(bm)

\method{print}{blockmodel}(x, ...)
}
\arguments{
\item{x}{A \code{\link{netframe}}; for \code{print}, a blockmodel.}
\item{membership}{The block of each node: a vector aligned with the rows of
\code{x$nodes}, or the name of a column of \code{x$nodes}, for example the
result of \code{\link{concor}} or \code{\link{regular_equivalence}}. Only
equality of values matters; blocks are listed in order of first
appearance.}
\item{weighted}{Whether to use total edge weights as tie values, so that
block densities are mean tie values. \code{NULL} uses weights when the
network has a weight column; they are used as given, whatever the
network's \code{weight_type}.}
\item{loops}{Whether self-loops and the diagonal of the adjacency matrix
count: with \code{FALSE} (the default, and the usual practice) the
diagonal blocks have \eqn{n_r (n_r - 1)} pairs, with \code{TRUE}
\eqn{n_r^2}.}
\item{cutoff}{The density at or above which a block is a 1-block in the
image. \code{NULL} uses the density of the whole network, the
alpha-density criterion of Arabie, Boorman, and Levitt (1978) that
UCINET's Blockmodel routine and Pajek apply by default.}
\item{bm}{A blockmodel from \code{blockmodel()}.}
\item{...}{Passed to \code{print()} for the density matrix.}
}
\details{
The density of block \eqn{(r, s)} is the sum of the ties from nodes in
\eqn{r} to nodes in \eqn{s} divided by the number of such pairs. Without
weights a tie is 1 when there is an edge from \eqn{i} to \eqn{j}, so
parallel edges count once; with weights it is the total weight of those
edges, as in sna's \code{blockmodel(block.content = "density")} on a valued
matrix. An undirected edge is a tie in both directions, so the density
matrix of an undirected network is symmetric. A one-node block without
loops has no pairs on the diagonal: its density and image are \code{NaN}.

The fit is the squared Pearson correlation between the ties and the
density of the block each pair of nodes falls in, over the pairs that
count (Wasserman and Faust 1994, ch. 10): the share of the variance in the
ties that the block densities account for, 1 for a perfect blockmodel and
0 when every block has the same density. It is \code{NaN} when the ties or
the densities have no variance.

The validation suite compares block densities with
\code{sna::blockmodel(block.content = "density")}, which also leaves the
diagonal out; UCINET's density blockmodel documents the same definition.
}
\value{\code{blockmodel()} returns a list of class \code{"blockmodel"} with
\code{membership} (the block labels, aligned with \code{x$nodes}),
\code{blocks} (the distinct labels, in order), \code{sizes}, the
block-by-block matrices \code{density}, \code{ties} (sums of tie values),
\code{pairs}, and \code{image} (1, 0, or \code{NaN}), the \code{cutoff}
used, \code{overall_density}, and \code{r_squared}. \code{block_summary()}
returns a data frame with one row per ordered pair of blocks: \code{from},
\code{to}, \code{size_from}, \code{size_to}, \code{ties}, \code{pairs},
\code{density}, and \code{image}.}
\references{
White, H. C., Boorman, S. A., & Breiger, R. L. (1976). Social structure
from multiple networks. I. Blockmodels of roles and positions.
\emph{American Journal of Sociology}, 81(4), 730-780. \doi{10.1086/226141}

Arabie, P., Boorman, S. A., & Levitt, P. R. (1978). Constructing
blockmodels: How and why. \emph{Journal of Mathematical Psychology},
17(1), 21-63. \doi{10.1016/0022-2496(78)90034-2}

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods
and Applications}. Cambridge University Press.
\doi{10.1017/CBO9780511815478}
}
\seealso{\code{\link{concor}}, \code{\link{regular_equivalence}},
\code{\link{structural_equivalence}}.}
\examples{
g <- example_network("florentine")
bm <- blockmodel(g, concor(g, depth = 2))
bm
bm$image
bm$r_squared
block_summary(bm)
}
