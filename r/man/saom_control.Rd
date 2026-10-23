\name{saom_control}
\alias{saom_control}
\title{Settings of the SAOM estimation algorithm}
\description{
Settings of the three-phase Robbins-Monro algorithm used by
\code{\link{saom_fit}}, with RSiena's defaults.
}
\usage{
saom_control(conditional = TRUE, nsub = 4, n3 = 1000, n1 = NULL, firstg = 0.2,
             reduceg = 0.5, diagonalize = 0.2, truncation = 5, theta_bound = 50,
             seed = NULL, threads = 0L, keep_simulations = TRUE)
}
\arguments{
\item{conditional}{Whether to estimate conditionally on the observed number
of changes per period (the default, as in RSiena): each period's simulation
runs until as many tie variables differ from the start of the period as were
observed to differ, and the rate parameters are estimated afterwards from
the simulated period lengths. \code{FALSE} simulates one time unit per
period and estimates the rates as parameters of the Robbins-Monro
iteration.}
\item{nsub}{Number of subphases in phase 2.}
\item{n3}{Number of simulations in phase 3, which give the standard errors,
convergence diagnostics, and goodness-of-fit statistics. 1000 is adequate
for exploration; use more for publication.}
\item{n1}{Number of phase 1 simulations; \code{NULL} for \eqn{7 + 3p}, with
\eqn{p} the number of parameters.}
\item{firstg}{Initial gain of the Robbins-Monro steps.}
\item{reduceg}{Factor by which the gain shrinks after each subphase.}
\item{diagonalize}{Weight in \eqn{[0, 1]} of the diagonal of the derivative
matrix in the step matrix (RSiena's default 0.2; 1 uses only the
diagonal).}
\item{truncation}{Phase 2 deviations larger than this many standard
deviations (estimated in phase 1) are truncated.}
\item{theta_bound}{Estimation stops with an error when a parameter exceeds
this in absolute value, a sign of divergence.}
\item{seed}{Seed for the simulations, a whole number between 0 and
\eqn{2^{53}}; \code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}()} makes results reproducible.}
\item{threads}{Number of threads for phases 1 and 3, a whole number; \code{0}
uses all available cores, or 2 during \verb{R CMD check} (see \sQuote{Threads}
in \code{\link{netframe}}). Results do not depend on the thread count except
in the last bits of floating-point sums.}
\item{keep_simulations}{Whether to keep the auxiliary statistics of every
phase 3 simulation for \code{\link{saom_gof}} (\eqn{n_3 (3n + 16)} numbers
for \eqn{n} actors).}
}
\value{
A list of class \code{"saom_control"}.
}
\seealso{
\code{\link{saom_fit}}
}
\examples{
saom_control(n3 = 500, seed = 1, threads = 2)
}
