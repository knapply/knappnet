\name{plot.saom_fit}
\alias{plot.saom_fit}
\title{Plot the estimation diagnostics of a SAOM}
\description{
\code{which = "trajectory"} draws the parameter values after every phase 2
iteration of the Robbins-Monro algorithm, with the subphases separated by
dashed lines and the final estimate as a horizontal line: a well-behaved
estimation settles down within each subphase. \code{which = "deviations"}
draws boxplots of the phase 3 deviations of the simulated statistics from
their targets, divided by their standard deviations; they should be
centered on 0.
}
\usage{
\method{plot}{saom_fit}(x, which = c("trajectory", "deviations"), ...)
}
\arguments{
\item{x}{An \code{\link{saom_fit}}.}
\item{which}{\code{"trajectory"} or \code{"deviations"}.}
\item{...}{Passed to the plotting functions.}
}
\value{
\code{x}, invisibly.
}
\seealso{
\code{\link{saom_fit}}, \code{\link{saom_gof}}
}
\examples{
set.seed(2)
start <- netframe(data.frame(from = sample(15, 30, TRUE), to = sample(15, 30, TRUE)),
                  nodes = data.frame(id = 1:15))
panel <- saom_simulate(saom_data(list(start, start)), ~ density + recip,
                       theta = c(-1.5, 1), rates = 3, seed = 1)[[1]]
fit <- saom_fit(saom_data(c(list(start), panel)), ~ density + recip,
                control = saom_control(n3 = 100, nsub = 2, seed = 2, threads = 2))
plot(fit)
plot(fit, which = "deviations")
}
