\name{plot_ego}
\alias{plot_ego}
\title{Draw an ego network}
\usage{
plot_ego(x, node, order = 1L, mode = c("all", "out", "in"), ...)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{node}{The identifier of the ego, as in the identifier column of
\code{x$nodes}.}
\item{order}{Alters are the nodes within this many steps of ego.}
\item{mode}{For directed networks, \code{"all"}, \code{"out"}, or
\code{"in"}; see \code{\link{ego_network}}.}
\item{...}{Further arguments passed to \code{\link{plot.netframe}}, such as
\code{layout}, \code{seed}, or \code{main}.}
}
\value{The ego network (a \code{\link{netframe}}), invisibly.}
\description{
Extracts the ego network of one node with \code{\link{ego_network}} and
draws it with ego enlarged and colored apart from its alters.
}
\seealso{\code{\link{ego_network}}, \code{\link{ego_summary}}.}
\examples{
g <- example_network("florentine")
plot_ego(g, "Medici", seed = 1)
plot_ego(g, "Medici", order = 2, seed = 1, main = "Within two steps")
}
