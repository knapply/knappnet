\name{mcmc_diagnostics}
\alias{mcmc_diagnostics}
\alias{print.ergm_mcmc}
\alias{plot.ergm_mcmc}
\title{Markov chain diagnostics for a Monte Carlo fit}
\description{
Summarizes the final Monte Carlo sample of an \code{\link{ergm_fit}}: for
each model statistic the sampled values minus the observed value (centered
on zero when the estimate solves the estimating equation), Geweke's (1992)
convergence diagnostic, the effective sample size, and autocorrelations,
as \code{mcmc.diagnostics()} in \pkg{ergm} reports them.
}
\usage{
mcmc_diagnostics(object, lags = c(1, 2, 5, 10, 20))

\method{print}{ergm_mcmc}(x, digits = max(3L, getOption("digits") - 3L), ...)

\method{plot}{ergm_mcmc}(x, which = NULL, ask = NULL, ...)
}
\arguments{
\item{object}{A fitted model from \code{\link{ergm_fit}} with
\code{estimate = "MCMLE"}.}
\item{lags}{Lags, in samples, at which to report autocorrelation.}
\item{x}{An \code{ergm_mcmc} object.}
\item{digits}{Digits to print.}
\item{which}{Statistics to plot, by name or position; default all.}
\item{ask}{Whether to prompt before each page on an interactive device.}
\item{...}{Not used.}
}
\details{
The sample holds \code{chains} chains of \code{nsim} draws each (see
\code{\link{ergm_control}}). Geweke's statistic compares the mean of the
first 10 percent of a chain with the mean of its last 50 percent, divided by
the square root of the sum of the two variances; each variance is estimated
by batch means (Geyer 1992) rather than the spectral density estimate of
the original paper. The statistic is approximately standard normal for a
chain in equilibrium, and the p-value is two-sided. The statistics of
independent chains are combined as \eqn{\sum_c z_c / \sqrt{C}}. The
effective sample size is the number of draws divided by the integrated
autocorrelation time estimated from the batch means, summed over chains.
Autocorrelations are averaged over chains.
}
\value{
An object of class \code{"ergm_mcmc"}: a list with \code{summary} (a data
frame with one row per statistic: \code{mean} and \code{sd} of the centered
sample, \code{geweke_z}, \code{geweke_p}, and \code{ess}),
\code{autocorrelation} (a lags-by-statistics matrix), \code{sample} (the
centered draws, chains stacked), \code{chains}, and \code{nsim}.
\code{plot()} draws, for each statistic, a trace plot with chain boundaries
marked, a density plot, and the autocorrelation function of the first
chain, with base graphics.
}
\references{
Geweke, J. (1992). Evaluating the accuracy of sampling-based approaches to
the calculation of posterior moments. In J. M. Bernardo, J. O. Berger, A.
P. Dawid, & A. F. M. Smith (Eds.), \emph{Bayesian Statistics 4} (pp.
169-193). Oxford University Press.

Geyer, C. J. (1992). Practical Markov chain Monte Carlo. \emph{Statistical
Science}, 7(4), 473-483. \doi{10.1214/ss/1177011137}
}
\seealso{\code{\link{ergm_fit}}, \code{\link{ergm_gof}}}
\examples{
flo <- example_network("florentine")
ctl <- ergm_control(seed = 1, threads = 2, burnin = 2000, interval = 20,
                    nsim = 500, loglik = FALSE)
fit <- ergm_fit(flo ~ edges + gwesp(0.5, fixed = TRUE), control = ctl)
d <- mcmc_diagnostics(fit)
d
plot(d, which = "edges")
}
