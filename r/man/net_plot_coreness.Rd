\name{net_plot_coreness}
\alias{net_plot_coreness}
\title{Plot the k-core shells of a network}
\description{
Places nodes on concentric rings by \code{\link{net_coreness}}, the innermost
ring holding the nodes of the highest core, and colors them by coreness,
the layered view of core structure that Seidman (1983) describes.
}
\usage{
net_plot_coreness(x, mode = c("all", "out", "in"), seed = NULL, threads = 0L, ...)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{mode}{For directed networks: \code{"all"}, \code{"out"}, or
\code{"in"}, as in \code{\link{net_coreness}}.}
\item{seed}{Seed for the layout whose angles order the nodes within each
ring (see Details), as in \code{\link{net_layout}}.}
\item{threads}{Number of threads for that layout, a whole number; \code{0}
uses all available cores, or 2 during \verb{R CMD check} (see \sQuote{Threads}
in \code{\link{net_frame}}).}
\item{...}{Further arguments passed to \code{\link{plot.net_frame}}, such as
\code{node_label} or \code{main}.}
}
\details{
Each distinct coreness value gets a ring, the highest in the center (a
single node of the highest core sits at the center itself). Within a ring
nodes are evenly spaced, ordered by their angle in a stress layout of the
whole network (\code{\link{net_layout}}), so that nodes close in the
network tend to be close on the ring. Rings are drawn as dotted circles
and the coreness values are listed in a legend inside the plot.
}
\value{The layout used, invisibly: a matrix with columns \code{x} and
\code{y}.}
\references{
Seidman, S. B. (1983). Network structure and minimum degree. \emph{Social
Networks}, 5(3), 269-287. \doi{10.1016/0378-8733(83)90028-X}
}
\seealso{\code{\link{net_coreness}}.}
\examples{
g <- net_example("karate")
net_plot_coreness(g, seed = 1, threads = 2, main = "Karate club: k-cores")
}
