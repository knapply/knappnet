\name{simulate.ergm_fit}
\alias{simulate.ergm_fit}
\alias{ergm_simulate}
\title{Simulate networks from an exponential random graph model}
\description{
Draws networks from a fitted model, or from a formula and coefficients, by
the Metropolis-Hastings sampler of \code{\link{ergm_fit}} (tie/no-tie
proposals by default; Morris, Handcock, and Hunter 2008). One chain starts
from the observed network, discards \code{burnin} proposals, and then
records a network every \code{interval} proposals, as \code{simulate()} in
the \pkg{ergm} package does.
}
\usage{
\method{simulate}{ergm_fit}(object, nsim = 1, seed = NULL, ..., coef = NULL,
         burnin = object$control$burnin, interval = object$control$interval,
         proposal = object$control$proposal, threads = object$control$threads,
         output = c("netframe", "stats"))

ergm_simulate(formula, coef, nsim = 1, seed = NULL, burnin = 10000,
              interval = 100, proposal = c("tnt", "random"),
              constraints = NULL, threads = 0L,
              output = c("netframe", "stats"))
}
\arguments{
\item{object}{A fitted model from \code{\link{ergm_fit}}.}
\item{formula}{A model formula, as in \code{\link{ergm_fit}}.}
\item{coef}{Coefficients, one per statistic in formula order. For a fitted
model, \code{NULL} uses the estimates.}
\item{nsim}{Number of networks to draw.}
\item{seed}{Seed, a whole number between 0 and \eqn{2^{53}}; \code{NULL} draws
one from R's random number generator, so \code{\link{set.seed}()} makes
results reproducible.}
\item{burnin, interval, proposal, threads}{Sampler settings; see
\code{\link{ergm_control}}. For a fitted model the defaults are its
control settings.}
\item{constraints}{Fixed dyads, as in \code{\link{ergm_fit}}.}
\item{output}{\code{"netframe"} returns the networks; \code{"stats"} only
the statistics of each draw.}
\item{...}{Not used.}
}
\details{
Under an \code{edges}-only model every dyad is an independent Bernoulli
tie, so the sampled tie counts have a known distribution; this, and the
agreement of sampled means with exact expectations on tiny networks, is
how the sampler is validated. Fixed dyads (\code{constraints}) keep their
observed state in every draw.
}
\value{
A list of \code{\link{netframe}}s with the same nodes as the model's
network, with the sampled statistics (one row per network) in attribute
\code{"stats"} and the seed in attribute \code{"seed"}; or, with
\code{output = "stats"}, the matrix of statistics alone.
}
\references{
Morris, M., Handcock, M. S., & Hunter, D. R. (2008). Specification of
exponential-family random graph models: Terms and computational aspects.
\emph{Journal of Statistical Software}, 24(4), 1-24.
\doi{10.18637/jss.v024.i04}
}
\seealso{\code{\link{ergm_fit}}, \code{\link{ergm_gof}}}
\examples{
flo <- example_network("florentine")
nets <- ergm_simulate(flo ~ edges + triangle, coef = c(-1.7, 0.2), nsim = 3,
                      seed = 1, burnin = 2000, interval = 50, threads = 2)
attr(nets, "stats")
nets[[1]]

fit <- ergm_fit(flo ~ edges, estimate = "MPLE")
simulate(fit, nsim = 2, seed = 2, output = "stats", threads = 2)
}
