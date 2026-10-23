\name{net_plot_degree}
\alias{net_plot_degree}
\title{Plot the degree distribution}
\description{
Draws how many nodes have each degree, as a bar chart of counts or as the
complementary cumulative distribution \eqn{P(D \ge d)}{P(D >= d)}, which
on logarithmic axes is the usual display for judging heavy tails (Newman
2003; Clauset, Shalizi, and Newman 2009).
}
\usage{
net_plot_degree(x, mode = c("all", "out", "in"), loops = TRUE,
                         type = c("histogram", "ccdf"), log = FALSE,
                         convention = c("igraph", "networkx", "sna", "ucinet"),
                         ...)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{mode, loops, convention}{As in \code{\link{net_degree}}.}
\item{type}{\code{"histogram"} for counts per degree, or \code{"ccdf"} for
the share of nodes with at least each degree.}
\item{log}{Whether to use logarithmic axes. Nodes of degree 0 cannot be
shown on a logarithmic degree axis and are left out, with a note in the
subtitle.}
\item{...}{Further arguments passed to \code{\link[graphics]{barplot}}
(histogram on linear axes) or \code{\link[graphics]{plot}} (otherwise),
such as \code{main}.}
}
\details{
Degrees are counted by \code{\link{net_degree}}, so parallel edges each count
and a self-loop adds 2 in an undirected network, unless
\code{loops = FALSE}. The histogram has one bar per degree value from 0 to
the largest (not binned intervals); with more than 60 distinct values, or
with \code{log = TRUE}, it is drawn as points instead. The CCDF is a step
function over the observed degrees; a straight line on logarithmic axes
suggests a power law, but a fitted model and a goodness-of-fit test
(Clauset et al. 2009) are needed to claim one.
}
\value{A data frame with one row per observed degree, invisibly:
\code{degree}, \code{count}, and \code{ccdf} (the share of nodes with at
least that degree).}
\references{
Newman, M. E. J. (2003). The structure and function of complex networks.
\emph{SIAM Review}, 45(2), 167-256. \doi{10.1137/S003614450342480}

Clauset, A., Shalizi, C. R., & Newman, M. E. J. (2009). Power-law
distributions in empirical data. \emph{SIAM Review}, 51(4), 661-703.
\doi{10.1137/070710111}
}
\seealso{\code{\link{net_degree}}, \code{\link{net_plot_scores}}.}
\examples{
g <- net_example("karate")
net_plot_degree(g)
net_plot_degree(g, type = "ccdf", log = TRUE)
}
