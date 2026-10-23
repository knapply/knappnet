\name{saom_fit}
\alias{saom_fit}
\alias{print.saom_fit}
\alias{summary.saom_fit}
\alias{print.summary.saom_fit}
\alias{coef.saom_fit}
\alias{vcov.saom_fit}
\title{Fit a stochastic actor-oriented model}
\usage{
saom_fit(data, effects, control = saom_control(), initial = NULL)
}
\arguments{
\item{data}{An \code{\link{saom_data}} object.}
\item{effects}{The effects of the evaluation function: a one-sided formula
such as \code{~ density + recip + transTrip + simX("smoke")}, a character
vector of structural effects, or the data frame returned by
\code{\link{saom_effects}}, which documents every effect and its
definition.}
\item{control}{Algorithm settings from \code{\link{saom_control}}.}
\item{initial}{Initial parameter values: \code{NULL} for the standard
initial values, or a numeric vector with one value per effect (rates first,
one per period, when \code{control$conditional} is \code{FALSE});
\code{NA} keeps the standard value for that parameter.}
}
\value{
An object of class \code{"saom_fit"}: a list with \code{coefficients},
\code{se}, \code{vcov}, \code{tconv} (convergence t-ratios),
\code{tconv_max} (the overall maximum convergence ratio), \code{rates} (a
data frame with the rate parameter and standard error of each period),
\code{effects}, \code{targets} (the observed statistics), \code{means}
(their phase 3 means), \code{derivative}, \code{statistic_vcov},
\code{deviations} (phase 3 deviations of the simulated statistics from the
targets, one row per simulation), \code{phase2_deviations},
\code{trajectory} (phase 2 parameter values by iteration), \code{times}
(phase 3 simulated period lengths), \code{aux} (auxiliary statistics for
\code{\link{saom_gof}}), \code{iterations}, \code{seed}, \code{control}, and
the \code{call}. \code{summary()} adds \code{converged}. \code{coef()} and
\code{vcov()} return the estimates and their covariance matrix.
}
\description{
Estimates a stochastic actor-oriented model (SAOM; Snijders 2001; Snijders,
van de Bunt, and Steglich 2010; Snijders 2017) of the changes between the
waves of a network panel by the method of moments, using the three-phase
Robbins-Monro stochastic approximation algorithm of RSiena.
}
\details{
The model is a continuous-time Markov chain of \emph{ministeps} (Snijders
2001): actor \eqn{i} gets opportunities to change at rate
\eqn{\rho_m}{rho_m} during period \eqn{m}, and at each opportunity adds or
drops one outgoing tie, or changes nothing, choosing the resulting network
\eqn{x'} with probability proportional to \eqn{\exp f_i(x')}{exp f_i(x')},
where the evaluation function
\eqn{f_i(x) = \sum_k \beta_k s_{ik}(x)}{f_i(x) = sum_k beta_k s_ik(x)} is a
weighted sum of the effects in \code{\link{saom_effects}}. The method of
moments chooses \eqn{\beta}{beta} so that the statistics
\eqn{\sum_m \sum_i s_{ik}}{sum_m sum_i s_ik} of networks simulated from
each wave to the next equal, on average, the observed statistics of the
next wave (with changing covariates at their values at the start of each
period).

Estimation follows RSiena's algorithm (Ripley et al. 2026, chapter 7, and
the document \emph{Siena algorithms}): phase 1 estimates the derivative of
the expected statistics with respect to the parameters by the
score-function method (Schweinberger and Snijders 2007) from \eqn{7 + 3p}
simulations; phase 2 runs \code{nsub} subphases of Robbins-Monro steps
\eqn{\theta \leftarrow \theta - a D^{-1} (S - s)}{theta <- theta - a D^-1 (S - s)},
each subphase ending once the deviations have crossed zero and setting
\eqn{\theta}{theta} to its average over the subphase, with the gain
\eqn{a} halved between subphases (Robbins and Monro 1951; Snijders 2001,
appendix); phase 3 simulates \code{n3} times at the final
\eqn{\theta}{theta} to obtain the convergence t-ratios
\eqn{(\bar S_k - s_k) / sd(S_k)}{(mean S_k - s_k) / sd(S_k)}, which should
be below 0.1 in absolute value, the overall maximum convergence ratio
\eqn{\sqrt{(\bar S - s)' \Sigma^{-1} (\bar S - s)}}{sqrt((mean S - s)' Sigma^-1 (mean S - s))},
which should be below 0.25, and the standard errors by the delta method,
\eqn{cov(\hat\theta) = D^{-1} \Sigma D^{-T}}{cov(theta) = D^-1 Sigma D^-T}.
Conditional estimation (the default) simulates each period until the
observed number of tie changes is reached and estimates the rate of each
period as the mean simulated period length in phase 3, reporting the
standard deviation of those lengths as its standard error (Snijders 2001,
section 4.2; RSiena may compute this standard error differently).

Standard initial values follow Snijders (2005, equation 11.41) in the form
RSiena uses: the out-degree parameter starts at
\eqn{\log(c^+ d / (c^- (1 - d))) / 2}{log(c+ d / (c- (1 - d))) / 2}, with
\eqn{c^+}{c+} ties created, \eqn{c^-}{c-} ties dissolved, and \eqn{d} the
mean density; other effects start at 0 and, in unconditional estimation,
each rate at twice the observed number of changes per actor. The
simulations use the engine's seeded, counter-based random number generator,
so a seed reproduces the estimate exactly for any number of threads.

The deterministic parts of the model (target statistics, covariate
centering, change statistics) are validated against RSiena and against a
direct evaluation of the definitions; the estimates are compared with
RSiena's on the same data, within the standard errors both report, in the
package tests.
}
\section{Not implemented}{
Behavior dynamics and co-evolution, several dependent networks,
non-directed networks, creation and endowment effects, rate effects beyond
a constant rate per period, missing data, structural zeros and ones,
composition change, maximum likelihood and Bayesian estimation, score-type
tests, and time heterogeneity tests. For these, see RSiena.
}
\references{
Snijders, T. A. B. (2001). The statistical evaluation of social network
dynamics. \emph{Sociological Methodology}, 31(1), 361-395.
\doi{10.1111/0081-1750.00099}

Snijders, T. A. B. (2005). Models for longitudinal network data. In P. J.
Carrington, J. Scott, and S. Wasserman (Eds.), \emph{Models and Methods in
Social Network Analysis} (pp. 215-247). Cambridge University Press.
\doi{10.1017/CBO9780511811395.011}

Snijders, T. A. B., van de Bunt, G. G., and Steglich, C. E. G. (2010).
Introduction to stochastic actor-based models for network dynamics.
\emph{Social Networks}, 32(1), 44-60. \doi{10.1016/j.socnet.2009.02.004}

Snijders, T. A. B. (2017). Stochastic actor-oriented models for network
dynamics. \emph{Annual Review of Statistics and Its Application}, 4,
343-363. \doi{10.1146/annurev-statistics-060116-054035}

Schweinberger, M., and Snijders, T. A. B. (2007). Markov models for digraph
panel data: Monte Carlo-based derivative estimation. \emph{Computational
Statistics & Data Analysis}, 51(9), 4465-4483.
\doi{10.1016/j.csda.2006.07.014}

Robbins, H., and Monro, S. (1951). A stochastic approximation method.
\emph{The Annals of Mathematical Statistics}, 22(3), 400-407.
\doi{10.1214/aoms/1177729586}

Ripley, R. M., Snijders, T. A. B., Boda, Z., and Voros, A. (2026).
\emph{Manual for RSiena}. University of Oxford, Department of Statistics and
Nuffield College; University of Groningen, Department of Sociology.
\url{https://www.stats.ox.ac.uk/~snijders/siena/RSiena_Manual.pdf}
}
\seealso{
\code{\link{saom_data}}, \code{\link{saom_effects}},
\code{\link{saom_control}}, \code{\link{saom_gof}},
\code{\link{saom_simulate}}, \code{\link{plot.saom_fit}}.
}
\examples{
# A small panel simulated from known parameters, then re-estimated.
set.seed(1)
start <- netframe(data.frame(from = sample(20, 40, TRUE), to = sample(20, 40, TRUE)),
                  nodes = data.frame(id = 1:20, smoke = rep(0:1, 10)))
panel <- saom_simulate(saom_data(list(start, start), covariates = "smoke"),
                       ~ density + recip + simX("smoke"), theta = c(-2, 1.5, 0.5),
                       rates = 4, seed = 1)[[1]]
data <- saom_data(c(list(start), panel), covariates = "smoke")
fit <- saom_fit(data, ~ density + recip + simX("smoke"),
                control = saom_control(n3 = 200, nsub = 2, seed = 2, threads = 2))
fit
summary(fit)
coef(fit)

# Unconditional estimation treats the rates as ordinary parameters.
saom_fit(data, ~ density + recip,
         control = saom_control(conditional = FALSE, n3 = 200, nsub = 2, seed = 3,
                                threads = 2))
}
