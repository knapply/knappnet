\name{net_plot_ego}
\alias{net_plot_ego}
\title{Plot an ego network}
\description{
Extracts the ego network of one node with \code{\link{net_ego}} and
draws it with ego enlarged and colored apart from its alters.
}
\usage{
net_plot_ego(x, node, order = 1L, mode = c("all", "out", "in"), ...)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{node}{The identifier of the ego, as in the identifier column of
\code{x$nodes}.}
\item{order}{Alters are the nodes within this many steps of ego.}
\item{mode}{For directed networks, \code{"all"}, \code{"out"}, or
\code{"in"}; see \code{\link{net_ego}}.}
\item{...}{Further arguments passed to \code{\link{plot.netframe}}, such as
\code{layout}, \code{seed}, or \code{main}.}
}
\value{The ego network (a \code{\link{net_frame}}), invisibly.}
\seealso{\code{\link{net_ego}}, \code{\link{net_ego_summary}}.}
\examples{
g <- net_example("florentine")
net_plot_ego(g, "Medici", seed = 1)
net_plot_ego(g, "Medici", order = 2, seed = 1, main = "Within two steps")
}
