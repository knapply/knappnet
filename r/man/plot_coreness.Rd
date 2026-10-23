\name{plot_coreness}
\alias{plot_coreness}
\title{Draw the k-core shells of a network}
\description{
Places nodes on concentric rings by \code{\link{coreness}}, the innermost
ring holding the nodes of the highest core, and colors them by coreness,
the layered view of core structure that Seidman (1983) describes.
}
\usage{
plot_coreness(x, mode = c("all", "out", "in"), seed = NULL, threads = 0L, ...)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{mode}{For directed networks: \code{"all"}, \code{"out"}, or
\code{"in"}, as in \code{\link{coreness}}.}
\item{seed}{Seed for the layout whose angles order the nodes within each
ring (see Details), as in \code{\link{layout_network}}.}
\item{threads}{Number of threads for that layout; \code{0} uses all
available cores.}
\item{...}{Further arguments passed to \code{\link{plot.netframe}}, such as
\code{node_label} or \code{main}.}
}
\details{
Each distinct coreness value gets a ring, the highest in the center (a
single node of the highest core sits at the center itself). Within a ring
nodes are evenly spaced, ordered by their angle in a stress layout of the
whole network (\code{\link{layout_network}}), so that nodes close in the
network tend to be close on the ring. Rings are drawn as dotted circles
and the coreness values are listed in a legend inside the plot.
}
\value{The layout used, invisibly: a matrix with columns \code{x} and
\code{y}.}
\references{
Seidman, S. B. (1983). Network structure and minimum degree. \emph{Social
Networks}, 5(3), 269-287. \doi{10.1016/0378-8733(83)90028-X}
}
\seealso{\code{\link{coreness}}.}
\examples{
g <- example_network("karate")
plot_coreness(g, seed = 1, threads = 2, main = "Karate club: k-cores")
}
