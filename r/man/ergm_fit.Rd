\name{ergm_fit}
\alias{ergm_fit}
\alias{ergm_stats}
\alias{ergm_control}
\alias{print.ergm_fit}
\alias{summary.ergm_fit}
\alias{print.summary.ergm_fit}
\alias{coef.ergm_fit}
\alias{vcov.ergm_fit}
\alias{logLik.ergm_fit}
\title{Exponential random graph models}
\description{
\code{ergm_fit()} fits an exponential random graph model (ERGM) to a
network from a formula of model terms written as in the statnet
\pkg{ergm} package. \code{ergm_stats()} evaluates the terms on the observed
network, and \code{ergm_control()} collects the estimation settings. The
functions carry an \code{ergm_} prefix rather than statnet's names
(\code{ergm()}, \code{gof()}, \code{mcmc.diagnostics()}) so that both
packages can be attached without masking each other.
}
\usage{
ergm_fit(formula, estimate = c("MCMLE", "MPLE", "exact"), constraints = NULL,
         control = ergm_control(), start = NULL)

ergm_stats(formula)

ergm_control(seed = NULL, threads = 0L, burnin = 10000, interval = 100,
             nsim = 1000, chains = 2L, proposal = c("tnt", "random"),
             max_iter = 60L, margin = 0.05, confidence = 0.95, loglik = TRUE,
             loglik_steps = 16L, mple_max_iter = 50L, mple_tol = 1e-12)

\method{summary}{ergm_fit}(object, ...)
\method{coef}{ergm_fit}(object, ...)
\method{vcov}{ergm_fit}(object, ...)
\method{logLik}{ergm_fit}(object, ...)
}
\arguments{
\item{formula}{A formula whose left-hand side is a \code{\link{netframe}}
and whose right-hand side sums model terms (see Details), for example
\code{g ~ edges + mutual + gwesp(0.5, fixed = TRUE) + nodematch("group")}.}
\item{estimate}{\code{"MCMLE"} (Monte Carlo maximum likelihood, the
default), \code{"MPLE"} (maximum pseudolikelihood), or \code{"exact"}
(maximum likelihood by enumerating every network; at most 22 free dyads).
A model whose terms are all dyad independent is fitted exactly by
\code{"MPLE"}, and \code{"MCMLE"} returns that fit.}
\item{constraints}{\code{NULL}, or an n-by-n logical or 0/1 matrix (rows
and columns in the order of \code{x$nodes}) marking dyads that keep their
observed state: structural zeros and ones. Fixed dyads are never toggled,
enter no estimate, and are left out of the count of dyads.}
\item{control}{Estimation settings from \code{ergm_control()}.}
\item{start}{Starting coefficients for \code{"MCMLE"} and \code{"exact"};
\code{NULL} uses the pseudolikelihood estimate when it is finite and zero
otherwise.}
\item{seed}{Seed for the Markov chains, a whole number between 0 and
\eqn{2^{53}}; \code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}()} makes results reproducible.}
\item{threads}{Number of threads for parallel chains, a whole number; \code{0}
uses all available cores, or 2 during \verb{R CMD check} (see \sQuote{Threads}
in \code{\link{netframe}}).}
\item{burnin, interval, nsim, chains}{Markov chain settings: proposals
discarded before sampling, proposals between samples, samples per chain,
and the number of chains run in parallel. Each chain has its own random
stream, so results do not depend on \code{threads}.}
\item{proposal}{\code{"tnt"} (tie/no-tie, the default and statnet's) or
\code{"random"} (a uniformly random dyad).}
\item{max_iter}{Maximum number of Monte Carlo MLE iterations.}
\item{margin}{How far inside the convex hull of the sampled statistics the
pseudo-observation of each step must lie, as a fraction (statnet's
\code{MCMLE.steplength.margin}).}
\item{confidence}{Level of the chi-squared test that declares the
estimating equation solved.}
\item{loglik}{Whether to estimate the log-likelihood of a Monte Carlo fit by
path sampling, so that \code{logLik()}, AIC, and BIC are available.}
\item{loglik_steps}{Number of steps of the path.}
\item{mple_max_iter, mple_tol}{Iteration limit and relative convergence
tolerance of the logistic regression behind the pseudolikelihood (R's
\code{glm()} uses 1e-8; the tighter default makes the estimate at least as
close to the optimum).}
\item{object}{A fitted model.}
\item{...}{Not used.}
}
\details{
An ERGM gives a network \eqn{y} on \eqn{n} nodes the probability
\deqn{P(Y = y) = \exp(\theta' g(y)) / c(\theta),}{P(Y = y) = exp(theta' g(y)) / c(theta),}
where \eqn{g(y)} is a vector of network statistics, \eqn{\theta} their
coefficients, and \eqn{c(\theta)} sums \eqn{\exp(\theta' g(y))}{exp(theta'
g(y))} over every network on the same nodes (Holland and Leinhardt 1981;
Frank and Strauss 1986; Wasserman and Pattison 1996; Snijders et al. 2006).

\strong{Terms.} Names and parameterizations follow statnet's \pkg{ergm}
package (Hunter et al. 2008; Morris, Handcock, and Hunter 2008) exactly, so
\code{ergm_stats()} agrees with \code{summary()} there and estimates are
directly comparable. Attributes are columns of \code{x$nodes}, named as
strings.
\describe{
\item{\code{edges}}{Number of ties.}
\item{\code{mutual}}{Pairs of nodes tied in both directions (directed).}
\item{\code{nodecov(attr)}, \code{nodeicov(attr)}, \code{nodeocov(attr)}}{
Sum over ties of a numeric attribute at both endpoints, at the head, or at
the tail.}
\item{\code{nodefactor(attr, levels = -1)}, \code{nodeifactor()},
\code{nodeofactor()}}{For each kept level of a categorical attribute, the
number of tie endpoints (heads, tails) at nodes with that level.
\code{levels} are positions among the sorted levels (negative to drop; the
default drops the first), level names, or \code{NULL} for all.}
\item{\code{nodematch(attr, diff = FALSE, levels = NULL)}}{Ties whose
endpoints share a level (uniform homophily); with \code{diff = TRUE} one
statistic per level (differential homophily).}
\item{\code{absdiff(attr, pow = 1)}}{Sum over ties of
\eqn{|x_i - x_j|^{pow}}.}
\item{\code{edgecov(x, attrname = NULL)}}{Sum over ties of a dyadic
covariate: an n-by-n matrix, or a netframe on the same nodes (whose edge
column \code{attrname}, or its tie count, fills the matrix). For undirected
ties the cell (min, max) is read, so an asymmetric matrix should be
symmetrized first.}
\item{\code{kstar(k)}, \code{istar(k)}, \code{ostar(k)}}{Number of
\eqn{k}-stars (in- and out-stars in directed networks); a vector \code{k}
gives one statistic per value.}
\item{\code{triangle}, \code{ttriple}, \code{ctriple}}{Triangles. In a
directed network \code{triangle} counts transitive plus cyclic triples,
which \code{ttriple} and \code{ctriple} count separately.}
\item{\code{gwesp(decay, fixed = TRUE)}, \code{gwdsp()}, \code{gwnsp()}}{
Geometrically weighted edgewise, dyadwise, and non-edgewise shared partners
(Snijders et al. 2006; Hunter 2007) with fixed decay \eqn{\alpha}:
\deqn{e^{\alpha} \sum_{k \ge 1} \{1 - (1 - e^{-\alpha})^k\} EP_k,}{exp(a) sum_k (1 - (1 - exp(-a))^k) EP_k,}
where \eqn{EP_k} is the number of ties (dyads, non-ties) with exactly
\eqn{k} shared partners. Shared partners of a directed tie
\eqn{i \to j}{i -> j} are the outgoing two-path partners \eqn{k} with
\eqn{i \to k \to j}{i -> k -> j}, statnet's default \code{type = "OTP"};
dyadwise statistics count ordered pairs in directed networks. Only
\code{fixed = TRUE} is supported (curved terms are not estimated).}
\item{\code{gwdegree(decay, fixed = TRUE)}, \code{gwidegree()},
\code{gwodegree()}}{Geometrically weighted degree, the same weighting
applied to the number of nodes of each degree.}
\item{\code{degree(d)}, \code{idegree(d)}, \code{odegree(d)}}{Nodes of
degree exactly \code{d} (a vector gives several statistics).}
\item{\code{isolates}}{Nodes with no ties.}
}
Networks are binary and simple: parallel edges collapse to one tie,
weights are ignored, and self-loops are an error, as in statnet's
\code{network} objects built from a sociomatrix. Categorical levels are
sorted (or a factor's levels). Statistics are validated to agree with
\code{summary()} in \pkg{ergm} on every network of the validation corpus.

\strong{Pseudolikelihood.} The maximum pseudolikelihood estimate (MPLE;
Strauss and Ikeda 1990) is the logistic regression of each free dyad's
state on its change statistics \eqn{\delta_{ij} = g(y^+_{ij}) -
g(y^-_{ij})}, fitted by iteratively reweighted least squares. When every
term is dyad independent (\code{edges}, the \code{nodecov}, \code{nodefactor},
\code{nodematch}, \code{absdiff}, and \code{edgecov} families) the dyads
are independent and the MPLE is the exact maximum likelihood estimate with
valid standard errors. For dependent terms the MPLE is a rough estimate
whose logistic standard errors are not valid (van Duijn, Gile, and Handcock
2009); \code{summary()} says so. A statistic that separates ties from
non-ties (as any count at its minimum or maximum does) has no finite
estimate, and a warning reports that the regression did not converge.

\strong{Monte Carlo maximum likelihood.} The Monte Carlo MLE (Geyer and
Thompson 1992) maximizes the estimated log-likelihood ratio
\deqn{l(\theta) - l(\theta_0) = (\theta - \theta_0)' g(y) - \log \frac{1}{m} \sum_k \exp\{(\theta - \theta_0)' g(Y_k)\}}{l(theta) - l(theta_0) = (theta - theta_0)' g(y) - log mean_k exp((theta - theta_0)' g(Y_k))}
over networks \eqn{Y_k} simulated at \eqn{\theta_0} by Metropolis-Hastings
with tie/no-tie proposals (Morris, Handcock, and Hunter 2008). Because the
approximation is good only near \eqn{\theta_0}, the fit iterates the
stepping algorithm of Hummel, Hunter, and Handcock (2012): starting from
the MPLE, each iteration samples at the current coefficients, finds the
largest step toward the observed statistics whose target lies inside the
convex hull of the sample (shrunk by \code{margin}; a linear program), and
takes a Newton step toward that target. Once the observed statistics lie in
the hull, Newton's method on the importance-sampled ratio gives the
estimate, a new sample is drawn at it, and the fit stops when the sampled
means agree with the observed statistics within Monte Carlo error (a
Hotelling-type chi-squared test at \code{confidence}). If \code{max_iter}
iterations pass without convergence the fit is returned with a warning.

Standard errors follow Hunter and Handcock (2006): the covariance of the
estimate is \eqn{I^{-1} + I^{-1} V I^{-1}}, where \eqn{I} is the Fisher
information estimated by the sample covariance of the statistics at the
estimate and \eqn{V} the Monte Carlo covariance of their mean, estimated by
batch means within each chain. \code{summary()} reports the share of each
variance due to the second term as \code{MCMC \%}, as \pkg{ergm} does.

The log-likelihood is estimated by path sampling (Gelman and Meng 1998;
Hunter and Handcock 2006, Section 4) from the null model, whose
log-likelihood is \eqn{-D \log 2} for \eqn{D} free dyads:
\eqn{l(\hat\theta) - l(0) = \hat\theta' g(y) - \int_0^1 \hat\theta' E_{t
\hat\theta}[g(Y)] dt}, with the integral evaluated by the midpoint rule
over \code{loglik_steps} samples. Its Monte Carlo standard error is
reported, and AIC and BIC (with the number of dyads as sample size, as in
\pkg{ergm}) follow from it.

A model that is degenerate at the current coefficients (the sampler sits at
the empty or the complete network, so the statistics never change;
Handcock 2003; Schweinberger 2011) stops with an error that says so.

\strong{Exact likelihood.} With \code{estimate = "exact"} every network on
the free dyads is enumerated in Gray-code order, one change statistic per
step, giving the exact normalizing constant, log-likelihood, and MLE by
Newton's method. This is possible for at most 22 free dyads (about 7
undirected or 5 directed nodes without constraints) and is the reference
against which the Monte Carlo MLE is validated: on every such network the
Monte Carlo estimate recovers the exact one within Monte Carlo error. The
exact MLE does not exist when the observed statistics lie on the boundary
of their convex hull (Handcock 2003), and an error says so.

\strong{Conventions and reproducibility.} Undirected ties are unordered
pairs and directed ties ordered pairs; the count of dyads \eqn{D} is
\eqn{n(n-1)/2} or \eqn{n(n-1)} minus the fixed dyads. Chain \eqn{c} draws
from stream \eqn{c} of \code{seed}, so a fit is reproducible for a given
seed on any number of threads and any platform.
}
\value{
\code{ergm_fit()} returns an object of class \code{"ergm_fit"}: a list with
\code{coefficients}; \code{vcov}, the full covariance matrix, and
\code{mcmc_vcov}, its Monte Carlo part (\code{NULL} unless
\code{estimate = "MCMLE"}); \code{estimate}; \code{stats}, the observed
statistics; \code{loglik} and \code{loglik_se} (\code{NULL} when the
log-likelihood is not available); \code{pseudo_loglik}; \code{null_loglik};
\code{dyads}; \code{mple}, the pseudolikelihood fit; \code{converged} and
\code{iterations}; \code{sample}, the statistics of the final Monte Carlo
sample (for \code{\link{mcmc_diagnostics}}); \code{model}; \code{control};
and the \code{call}.

\code{summary()} gives a coefficient table with standard errors, the
percentage of each variance due to Monte Carlo error (\code{MCMC \%}), z
values and p-values, and the deviance, AIC, and BIC when the
log-likelihood is available. \code{logLik()} returns the (exact,
path-sampled, or, for dyad-independent models, logistic) log-likelihood
with attributes \code{df}, \code{nobs} (the dyads), and \code{se}, and
stops when it is not available. \code{ergm_stats()} returns the observed
statistics as a named numeric vector, like statnet's
\code{summary(formula)}. \code{ergm_control()} returns a list.
}
\references{
Holland, P. W., & Leinhardt, S. (1981). An exponential family of
probability distributions for directed graphs. \emph{Journal of the American
Statistical Association}, 76(373), 33-50.
\doi{10.1080/01621459.1981.10477598}

Frank, O., & Strauss, D. (1986). Markov graphs. \emph{Journal of the
American Statistical Association}, 81(395), 832-842.
\doi{10.1080/01621459.1986.10478342}

Wasserman, S., & Pattison, P. (1996). Logit models and logistic regressions
for social networks: I. An introduction to Markov graphs and p*.
\emph{Psychometrika}, 61(3), 401-425. \doi{10.1007/BF02294547}

Snijders, T. A. B., Pattison, P. E., Robins, G. L., & Handcock, M. S.
(2006). New specifications for exponential random graph models.
\emph{Sociological Methodology}, 36(1), 99-153.
\doi{10.1111/j.1467-9531.2006.00176.x}

Hunter, D. R. (2007). Curved exponential family models for social networks.
\emph{Social Networks}, 29(2), 216-230. \doi{10.1016/j.socnet.2006.08.005}

Strauss, D., & Ikeda, M. (1990). Pseudolikelihood estimation for social
networks. \emph{Journal of the American Statistical Association}, 85(409),
204-212. \doi{10.1080/01621459.1990.10475327}

van Duijn, M. A. J., Gile, K. J., & Handcock, M. S. (2009). A framework for
the comparison of maximum pseudo-likelihood and maximum likelihood
estimation of exponential family random graph models. \emph{Social
Networks}, 31(1), 52-62. \doi{10.1016/j.socnet.2008.10.003}

Geyer, C. J., & Thompson, E. A. (1992). Constrained Monte Carlo maximum
likelihood for dependent data. \emph{Journal of the Royal Statistical
Society: Series B}, 54(3), 657-683. \doi{10.1111/j.2517-6161.1992.tb01443.x}

Hunter, D. R., & Handcock, M. S. (2006). Inference in curved exponential
family models for networks. \emph{Journal of Computational and Graphical
Statistics}, 15(3), 565-583. \doi{10.1198/106186006X133069}

Hummel, R. M., Hunter, D. R., & Handcock, M. S. (2012). Improving
simulation-based algorithms for fitting ERGMs. \emph{Journal of
Computational and Graphical Statistics}, 21(4), 920-939.
\doi{10.1080/10618600.2012.679224}

Hunter, D. R., Handcock, M. S., Butts, C. T., Goodreau, S. M., & Morris, M.
(2008). ergm: A package to fit, simulate and diagnose exponential-family
models for networks. \emph{Journal of Statistical Software}, 24(3), 1-29.
\doi{10.18637/jss.v024.i03}

Morris, M., Handcock, M. S., & Hunter, D. R. (2008). Specification of
exponential-family random graph models: Terms and computational aspects.
\emph{Journal of Statistical Software}, 24(4), 1-24.
\doi{10.18637/jss.v024.i04}

Gelman, A., & Meng, X.-L. (1998). Simulating normalizing constants: From
importance sampling to bridge sampling to path sampling. \emph{Statistical
Science}, 13(2), 163-185. \doi{10.1214/ss/1028905934}

Handcock, M. S. (2003). Assessing degeneracy in statistical models of
social networks. Working Paper 39, Center for Statistics and the Social
Sciences, University of Washington.

Schweinberger, M. (2011). Instability, sensitivity, and degeneracy of
discrete exponential families. \emph{Journal of the American Statistical
Association}, 106(496), 1361-1370. \doi{10.1198/jasa.2011.tm10747}
}
\seealso{\code{\link{simulate.ergm_fit}}, \code{\link{predict.ergm_fit}},
\code{\link{ergm_gof}}, \code{\link{mcmc_diagnostics}},
\code{\link{tergm_fit}}}
\examples{
# Padgett's Florentine marriage network with the families' wealth
# (thousands of lira), as in the ergm package's flomarriage data.
flo <- example_network("florentine")
wealth <- c(Acciaiuoli = 10, Albizzi = 36, Barbadori = 55, Bischeri = 44,
            Castellani = 20, Ginori = 32, Guadagni = 8, Lamberteschi = 42,
            Medici = 103, Pazzi = 48, Peruzzi = 49, Ridolfi = 27,
            Salviati = 10, Strozzi = 146, Tornabuoni = 48)
flo$nodes$wealth <- wealth[flo$nodes$name]

# Observed statistics, as summary(flomarriage ~ ...) in ergm.
ergm_stats(flo ~ edges + triangle + nodecov("wealth") + gwesp(0.5, fixed = TRUE))

# A dyad-independent model: the pseudolikelihood is the exact likelihood.
fit <- ergm_fit(flo ~ edges + nodecov("wealth"), estimate = "MPLE")
summary(fit)

# A dependent model by Monte Carlo maximum likelihood (small settings so
# the example is quick; the defaults are larger).
ctl <- ergm_control(seed = 1, threads = 2, burnin = 2000, interval = 20,
                    nsim = 500, loglik_steps = 8)
fit2 <- ergm_fit(flo ~ edges + gwesp(0.5, fixed = TRUE), control = ctl)
summary(fit2)
logLik(fit2)

# The same model by exact enumeration is possible only on tiny networks:
tiny <- netframe(data.frame(from = c(1, 1, 2, 2, 3, 4),
                            to   = c(2, 3, 3, 4, 5, 6)), directed = FALSE)
exact <- ergm_fit(tiny ~ edges + gwesp(0.5, fixed = TRUE), estimate = "exact")
coef(exact)
}
