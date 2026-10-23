\name{saom_gof}
\alias{saom_gof}
\alias{print.saom_gof}
\alias{plot.saom_gof}
\title{Goodness of fit of a stochastic actor-oriented model}
\description{
Compares a statistic of the observed networks (summed over the waves after
the first) with its distribution over the networks simulated in phase 3 of
the estimation, as proposed by Lospinoso and Snijders (2019) and
implemented in RSiena's \code{sienaGOF()}. The discrepancy is the
Mahalanobis distance of the observed statistic from the mean of the
simulated ones, using a generalized (Moore-Penrose) inverse of their
covariance matrix, and the p-value is the proportion of simulations whose
own distance is at least as large. A small p-value means the model does
not reproduce this feature of the data.
}
\usage{
saom_gof(fit, statistic = c("indegree", "outdegree", "triad_census", "geodesic"),
         levels = NULL, cumulative = TRUE)

\method{plot}{saom_gof}(x, violin = TRUE, ...)
}
\arguments{
\item{fit}{An \code{\link{saom_fit}} estimated with
\code{keep_simulations = TRUE} (the default).}
\item{statistic}{Which auxiliary statistic: the in-degree or out-degree
distribution, the triad census (Holland and Leinhardt 1976), or the
distribution of geodesic distances between ordered pairs of actors.}
\item{levels}{For the degree distributions, the degrees to report (default
0 to 8); for geodesics, the distances (default 1 to 5 and \code{Inf} for
unreachable pairs). Ignored for the triad census.}
\item{cumulative}{Whether degree and geodesic counts are cumulative (actors
with degree at most \eqn{k}; pairs at distance at most \eqn{d}), as in
RSiena.}
\item{x}{An \code{"saom_gof"} object.}
\item{violin}{Whether to draw the simulated distributions as violins
(kernel densities) rather than boxplots.}
\item{...}{Passed to \code{\link[graphics]{plot}}.}
}
\value{
An object of class \code{"saom_gof"} with \code{observed}, \code{simulated}
(one row per phase 3 simulation), \code{p} (the Monte Carlo p-value),
\code{distance} (the observed Mahalanobis distance), and \code{distances}
(one per simulation). \code{plot()} draws the simulated distributions with
the observed values as a red line and returns \code{x} invisibly.
}
\references{
Lospinoso, J., and Snijders, T. A. B. (2019). Goodness of fit for
stochastic actor-oriented models. \emph{Methodological Innovations}, 12(3).
\doi{10.1177/2059799119884282}

Holland, P. W., and Leinhardt, S. (1976). Local structure in social
networks. \emph{Sociological Methodology}, 7, 1-45. \doi{10.2307/270703}
}
\seealso{
\code{\link{saom_fit}}, \code{\link{triad_census}}
}
\examples{
set.seed(3)
start <- netframe(data.frame(from = sample(15, 30, TRUE), to = sample(15, 30, TRUE)),
                  nodes = data.frame(id = 1:15))
panel <- saom_simulate(saom_data(list(start, start)), ~ density + recip,
                       theta = c(-1.5, 1), rates = 3, seed = 1)[[1]]
fit <- saom_fit(saom_data(c(list(start), panel)), ~ density + recip,
                control = saom_control(n3 = 100, nsub = 2, seed = 2, threads = 2))
gof <- saom_gof(fit, "indegree")
gof
plot(gof)
plot(saom_gof(fit, "triad_census"))
}
