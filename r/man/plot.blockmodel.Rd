\name{plot.blockmodel}
\alias{plot.blockmodel}
\title{Plot a blockmodel}
\description{
Draws the block densities of a \code{\link{net_blockmodel}} as a heat map,
with the 1-blocks of the image matrix outlined, and, when the network is
given, its adjacency matrix permuted so that nodes of the same block are
adjacent, with lines at the block boundaries: the two displays used to
present blockmodels since White, Boorman, and Breiger (1976).
}
\usage{
\method{plot}{blockmodel}(x, network = NULL, which = c("both", "density", "matrix"),
     values = TRUE, ...)
}
\arguments{
\item{x}{A \code{\link{net_blockmodel}}.}
\item{network}{The \code{\link{net_frame}} the blockmodel was computed on,
for the permuted adjacency matrix. \code{NULL} draws only the densities.}
\item{which}{\code{"both"} (when \code{network} is given),
\code{"density"}, or \code{"matrix"}.}
\item{values}{Whether to write each density in its cell.}
\item{...}{Further arguments passed to \code{\link[graphics]{image}} for
both panels, such as \code{main} (which replaces the panel titles) or
\code{cex_labels}.}
}
\details{
Densities are shaded from light (0) to dark (1, or the largest mean tie
value for a weighted blockmodel); blocks whose density reaches the cutoff
(the 1-blocks of \code{x$image}) are outlined. Blocks with no pairs (a
one-node block without loops) are blank. With both panels, the permuted
adjacency matrix is drawn on the left and the densities on the right.
}
\value{\code{x}, invisibly.}
\references{
White, H. C., Boorman, S. A., & Breiger, R. L. (1976). Social structure
from multiple networks. I. Blockmodels of roles and positions.
\emph{American Journal of Sociology}, 81(4), 730-780. \doi{10.1086/226141}
}
\seealso{\code{\link{net_blockmodel}}, \code{\link{net_plot_matrix}}.}
\examples{
g <- net_example("florentine")
bm <- net_blockmodel(g, net_concor(g, depth = 2))
plot(bm)
plot(bm, network = g)
}
