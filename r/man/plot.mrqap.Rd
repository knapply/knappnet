\name{plot.mrqap}
\alias{plot.mrqap}
\title{Plot a multiple regression QAP}
\description{
Draws either a coefficient plot, each observed statistic against the
interval of its permutation distribution, or the permutation distribution
of every coefficient as a histogram with the observed statistic marked, as
in the diagnostics Dekker, Krackhardt, and Snijders (2007) use to compare
QAP variants.
}
\usage{
\method{plot}{mrqap}(x, which = c("coefficients", "distributions"), level = 0.95,
     intercept = FALSE, breaks = "Sturges", ...)
}
\arguments{
\item{x}{An \code{\link{net_mrqap}} result.}
\item{which}{\code{"coefficients"} for the coefficient plot or
\code{"distributions"} for the histograms.}
\item{level}{Coverage of the permutation intervals in the coefficient
plot, between 0 and 1.}
\item{intercept}{Whether to include the intercept, whose scale is often
very different from the slopes.}
\item{breaks}{For \code{"distributions"}, passed to
\code{\link[graphics]{hist}}.}
\item{...}{Further arguments passed to \code{\link[graphics]{plot}}
(coefficient plot) or \code{\link[graphics]{hist}} (distributions), such
as \code{main}.}
}
\details{
The statistics shown are the ones the permutation test uses
(\code{x$statistic}): t-statistics (the default of \code{\link{net_mrqap}},
recommended by Dekker et al. 2007) or raw coefficients. In the coefficient
plot each row has a horizontal segment from the \code{(1 - level) / 2} to
the \code{1 - (1 - level) / 2} quantile of the permutation distribution
(under the null hypothesis) and a point at the observed statistic; a point
outside its segment is significant at level \code{1 - level} in the
two-sided test, whose p-value \code{Pr(>=|b|)} is written beside it. Under
double semi-partialling the permutation distributions are centred near 0
for every coefficient; under Y-permutation they need not be.
}
\value{\code{x}, invisibly.}
\references{
Dekker, D., Krackhardt, D., & Snijders, T. A. B. (2007). Sensitivity of
MRQAP tests to collinearity and autocorrelation conditions.
\emph{Psychometrika}, 72(4), 563-581. \doi{10.1007/s11336-007-9016-1}
}
\seealso{\code{\link{net_mrqap}}, \code{\link{plot.qap_cor}}.}
\examples{
set.seed(2)
x1 <- matrix(rnorm(144), 12)
x2 <- matrix(rnorm(144), 12)
y <- 0.8 * x1 + matrix(rnorm(144), 12)
fit <- net_mrqap(y, list(same = x1, other = x2), reps = 300, seed = 1, threads = 2)
plot(fit)
plot(fit, which = "distributions")
}
