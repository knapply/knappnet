\name{saom_simulate}
\alias{saom_simulate}
\title{Simulate network panels from a stochastic actor-oriented model}
\description{
Simulates the change from each wave of \code{data} to the next under given
parameters, as \code{\link{saom_fit}} does internally: for one time unit
per period with the given rates, or conditional on the observed number of
changes (\code{rates = NULL}). Useful for simulation studies and for
generating panels from a known model.
}
\usage{
saom_simulate(data, effects, theta, rates = NULL, n = 1, seed = NULL)
}
\arguments{
\item{data}{An \code{\link{saom_data}} object; each period is simulated
from the observed network at its start.}
\item{effects}{The effects, as in \code{\link{saom_fit}}.}
\item{theta}{Effect parameters, one per effect, in the order of
\code{saom_effects(effects)}.}
\item{rates}{Rate parameters, one per period (or one value recycled), or
\code{NULL} to simulate until the observed number of changes of each
period is reached, as conditional estimation does.}
\item{n}{Number of panels to simulate.}
\item{seed}{Seed, a whole number between 0 and \eqn{2^{53}}; \code{NULL} draws
one from R's random number generator, so \code{\link{set.seed}()} makes
results reproducible. Simulation \eqn{s} uses stream \eqn{s - 1} of the seed.}
}
\value{
A list of \code{n} simulated panels. Each panel is a list with one
\code{\link{netframe}} per period, the network at the end of the period, on
the nodes of \code{data}; its attribute \code{"statistics"} holds the
effect statistics of each period's network (one row per period) and
\code{"times"} the simulated period lengths (the stopping times under rate
1 when \code{rates} is \code{NULL}, otherwise 1).
}
\seealso{
\code{\link{saom_fit}}, \code{\link{saom_gof}}
}
\examples{
w1 <- netframe(data.frame(from = c("a", "b", "c"), to = c("b", "c", "a")),
               nodes = data.frame(name = c("a", "b", "c", "d")))
sims <- saom_simulate(saom_data(list(w1, w1)), ~ density + recip, theta = c(-1, 1),
                      rates = 2, n = 2, seed = 7)
sims[[1]][[1]]
attr(sims[[1]], "statistics")
}
