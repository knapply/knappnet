\name{tergm_fit}
\alias{tergm_fit}
\alias{print.tergm_fit}
\alias{summary.tergm_fit}
\alias{print.summary.tergm_fit}
\alias{coef.tergm_fit}
\alias{confint.tergm_fit}
\title{Temporal exponential random graph models by pooled pseudolikelihood}
\usage{
tergm_fit(formula, boot = 200, seed = NULL, level = 0.95, threads = 0L,
          control = ergm_control())

\method{summary}{tergm_fit}(object, ...)
\method{coef}{tergm_fit}(object, ...)
\method{confint}{tergm_fit}(object, parm, level = NULL, ...)
}
\arguments{
\item{formula}{A formula whose left-hand side is a list of
\code{\link{netframe}}s in time order (for example from
\code{\link{snapshots}}), all with the same nodes in the same order, and
whose right-hand side sums the model terms of \code{\link{ergm_fit}} and
the temporal terms below.}
\item{boot}{Number of bootstrap replicates.}
\item{seed}{Seed for the bootstrap, a whole number in [0, 2^53);
\code{NULL} draws one from R's random number generator.}
\item{level}{Confidence level of the percentile intervals; for
\code{confint()}, \code{NULL} keeps the fit's level.}
\item{threads}{Number of threads; \code{0} uses all available cores, or 2
during \verb{R CMD check}.}
\item{control}{Pseudolikelihood settings (\code{mple_max_iter},
\code{mple_tol}) from \code{\link{ergm_control}}.}
\item{object}{A \code{tergm_fit} object.}
\item{parm}{Coefficients to report, by name or position; default all.}
\item{...}{Not used.}
}
\value{
An object of class \code{"tergm_fit"}: a list with \code{coefficients};
\code{boot}, a boot-by-p matrix of replicate estimates (\code{NA} where a
replicate did not converge); \code{ci}, the percentile intervals;
\code{level}; \code{transitions}; \code{dyads} per transition;
\code{deviance} and \code{null_deviance} of the pooled fit; \code{seed};
\code{names}; and the \code{call}. \code{summary()} gives the coefficients
with their intervals; \code{coef()} and \code{confint()} extract them.
}
\description{
Fits a discrete-time temporal ERGM (TERGM; Hanneke, Fu, and Xing 2010) to a
sequence of networks on the same nodes by maximum pseudolikelihood pooled
over the transitions, with bootstrap confidence intervals from resampling
the transitions (Desmarais and Cranmer 2012), the method of the
\pkg{btergm} package (Leifeld, Cranmer, and Desmarais 2018).
}
\details{
Each network \eqn{y_t} (\eqn{t = 2, \ldots, T}) is modeled conditional on
the previous one,
\deqn{P(Y_t = y_t \mid y_{t-1}) = \exp\{\theta' g(y_t, y_{t-1})\} / c(\theta, y_{t-1}),}{P(Y_t = y_t | y_{t-1}) = exp(theta' g(y_t, y_{t-1})) / c(theta, y_{t-1}),}
with the ERGM terms of \code{\link{ergm_fit}} evaluated on \eqn{y_t} and
temporal terms that read \eqn{y_{t-1}} as a dyadic covariate:
\describe{
\item{\code{memory(type = "stability")}}{The number of dyads in the same
state as at the previous time (ties that persist and non-ties that
persist). Its change statistic is \eqn{2 y^{t-1}_{ij} - 1}, so the term is
\code{edgecov(2 * y_prev - 1)}. \code{type = "autoregression"} counts
persisting ties (\code{edgecov(y_prev)}), \code{"innovation"} newly formed
ties (\code{edgecov(1 - y_prev)}), and \code{"loss"} dissolved ties
(\code{edgecov(-y_prev)}).}
\item{\code{delrecip}}{Delayed reciprocity: ties \eqn{i \to j}{i -> j}
whose reverse \eqn{j \to i}{j -> i} existed at the previous time,
\code{edgecov(t(y_prev))}; directed networks only.}
}
Memory terms use the previous network's ties as 0/1 regardless of weights,
with self-loops ignored.

The pooled pseudolikelihood stacks the change statistics of every free
dyad of every transition and fits one logistic regression (Strauss and
Ikeda 1990; see \code{\link{ergm_fit}}); the point estimates equal those of
\code{btergm()}. Because dyads within a network are dependent, the
logistic standard errors are not used. Instead the transitions are
resampled with replacement \code{boot} times, the pooled fit is repeated on
each resample, and confidence intervals are the percentiles of the
replicate estimates; replicates whose fit does not converge are dropped
and counted in \code{summary()}. With few transitions the intervals are
coarse, which is a property of the method. Replicate \eqn{r} draws from
stream \eqn{r} of \code{seed}, so results do not depend on \code{threads}.
}
\references{
Hanneke, S., Fu, W., & Xing, E. P. (2010). Discrete temporal models of
social networks. \emph{Electronic Journal of Statistics}, 4, 585-605.
\doi{10.1214/09-EJS548}

Desmarais, B. A., & Cranmer, S. J. (2012). Statistical mechanics of
networks: Estimation and uncertainty. \emph{Physica A}, 391(4), 1865-1876.
\doi{10.1016/j.physa.2011.10.018}

Leifeld, P., Cranmer, S. J., & Desmarais, B. A. (2018). Temporal exponential
random graph models with btergm: Estimation and bootstrap confidence
intervals. \emph{Journal of Statistical Software}, 83(6), 1-36.
\doi{10.18637/jss.v083.i06}

Strauss, D., & Ikeda, M. (1990). Pseudolikelihood estimation for social
networks. \emph{Journal of the American Statistical Association}, 85(409),
204-212. \doi{10.1080/01621459.1990.10475327}
}
\seealso{\code{\link{ergm_fit}}, \code{\link{snapshots}}}
\examples{
# Three waves of a friendship network: the second and third are the first
# with some ties rewired, so ties mostly persist.
set.seed(1)
n <- 12
wave <- function(prev, k) {
  e <- prev
  drop <- sample(nrow(e), k)
  e <- e[-drop, ]
  add <- data.frame(from = sample(n, k, TRUE), to = sample(n, k, TRUE))
  add <- add[add$from != add$to, ]
  unique(rbind(e, add))
}
e1 <- data.frame(from = sample(n, 30, TRUE), to = sample(n, 30, TRUE))
e1 <- unique(e1[e1$from != e1$to, ])
e2 <- wave(e1, 6)
e3 <- wave(e2, 6)
nodes <- data.frame(id = 1:n, group = rep(c("a", "b"), 6))
nets <- lapply(list(e1, e2, e3), netframe, nodes = nodes, directed = TRUE)

fit <- tergm_fit(nets ~ edges + mutual + memory("stability") + delrecip +
                   nodematch("group"), boot = 100, seed = 1, threads = 2)
summary(fit)
confint(fit)
}
