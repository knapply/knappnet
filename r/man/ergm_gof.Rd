\name{ergm_gof}
\alias{ergm_gof}
\alias{print.ergm_gof}
\alias{plot.ergm_gof}
\title{Goodness of fit of an exponential random graph model}
\description{
Simulates networks from a fitted model and compares distributions of
structural features that need not be in the model with the observed
network (Hunter, Goodreau, and Handcock 2008): the degree distribution,
edgewise and dyadwise shared partners, geodesic distances, the triad
census, and the model's own statistics.
}
\usage{
ergm_gof(object, nsim = 100, seed = NULL,
         stats = c("degree", "esp", "dsp", "distance", "triad", "model"),
         burnin = object$control$burnin, interval = object$control$interval,
         threads = object$control$threads)

\method{print}{ergm_gof}(x, digits = max(3L, getOption("digits") - 3L), ...)

\method{plot}{ergm_gof}(x, which = NULL, ask = NULL, ...)
}
\arguments{
\item{object}{A fitted model from \code{\link{ergm_fit}}.}
\item{nsim}{Number of networks to simulate.}
\item{seed}{Seed for the simulation, a whole number between 0 and
\eqn{2^{53}}; \code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}()} makes results reproducible.}
\item{stats}{Which statistics to compare: any of \code{"degree"} (in- and
out-degree for directed networks), \code{"esp"}, \code{"dsp"},
\code{"distance"}, \code{"triad"} (directed networks only), and
\code{"model"}.}
\item{burnin, interval, threads}{Sampler settings, as in
\code{\link{simulate.ergm_fit}}.}
\item{x}{An \code{ergm_gof} object.}
\item{digits}{Digits to print.}
\item{which}{Names of the statistics to plot; default all.}
\item{ask}{Whether to prompt before each page when several statistics are
plotted on an interactive device.}
\item{...}{For \code{plot()}, further arguments to
\code{\link[graphics]{boxplot}}.}
}
\details{
For each statistic the table reports the observed value, the minimum, mean,
and maximum over the simulated networks, and a Monte Carlo p-value: twice
the smaller of the proportions of simulated values at or below and at or
above the observed value, capped at 1, as \code{gof()} in \pkg{ergm}
reports it.

The statistics are those of \code{gof()}'s defaults: the number of nodes
with each degree \eqn{0, \ldots, n - 1}; the number of ties (edgewise) or
dyads (dyadwise) with each number of shared partners \eqn{0, \ldots, n -
2}, where shared partners of a directed tie \eqn{i \to j}{i -> j} are the
outgoing two-path partners \eqn{i \to k \to j}{i -> k -> j} and dyadwise
counts are over ordered pairs in directed networks; the number of pairs
(unordered for undirected, ordered for directed networks) at each geodesic
distance \eqn{1, \ldots, n - 1}, with unreachable pairs in the \code{Inf}
cell; and, for directed networks, the triad census in the order of
\code{\link{triad_census}}.
}
\value{
An object of class \code{"ergm_gof"}: a list with one element per
statistic (\code{degree} or \code{idegree} and \code{odegree}, \code{esp},
\code{dsp}, \code{distance}, \code{triadcensus}, \code{model}), each a
list of \code{observed}, \code{simulated} (an nsim-by-k matrix), and
\code{summary} (a data frame with columns \code{obs}, \code{min},
\code{mean}, \code{max}, and \code{p}); plus \code{nsim} and \code{seed}.
\code{plot()} draws one panel per statistic: box plots of the simulated
values with the observed values as a line, like \code{plot.gof()} in
\pkg{ergm}.
}
\references{
Hunter, D. R., Goodreau, S. M., & Handcock, M. S. (2008). Goodness of fit
of social network models. \emph{Journal of the American Statistical
Association}, 103(481), 248-258. \doi{10.1198/016214507000000446}
}
\seealso{\code{\link{ergm_fit}}, \code{\link{mcmc_diagnostics}}}
\examples{
flo <- example_network("florentine")
fit <- ergm_fit(flo ~ edges, estimate = "MPLE")
g <- ergm_gof(fit, nsim = 50, seed = 1, burnin = 1000, interval = 20, threads = 2)
g
plot(g, which = "degree")
}
