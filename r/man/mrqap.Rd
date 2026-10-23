\name{mrqap}
\alias{mrqap}
\alias{print.mrqap}
\alias{summary.mrqap}
\alias{print.summary.mrqap}
\title{Multiple regression QAP}
\description{
Regresses one network, or any \eqn{n \times n}{n x n} matrix of dyadic data,
on others by ordinary least squares over the dyads, with p-values from the
quadratic assignment procedure: double semi-partialling (Dekker,
Krackhardt, and Snijders 2007) or permutation of the response (Krackhardt
1988).
}
\usage{
mrqap(y, x, reps = 1000, method = c("dsp", "y"), seed = NULL,
      intercept = TRUE, mode = c("digraph", "graph"), diag = FALSE,
      statistic = c("t", "beta"), exact = FALSE, attr = NULL,
      tol = 1e-7, threads = 0L)
}
\arguments{
\item{y}{The response: a square numeric matrix or a
\code{\link{netframe}}.}
\item{x}{The predictors: a matrix or netframe, or a list of them. List
names become coefficient names (\code{x1}, \code{x2}, \ldots otherwise).}
\item{reps}{Number of random permutations.}
\item{method}{\code{"dsp"} (double semi-partialling, the default) or
\code{"y"} (permutation of the response).}
\item{seed}{Seed for the permutations, a whole number between 0 and
\eqn{2^{53}}; \code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}()} makes results reproducible.}
\item{intercept}{Whether to include an intercept.}
\item{mode}{\code{"digraph"} uses every ordered pair of distinct nodes;
\code{"graph"} uses each unordered pair once, read from the lower triangle,
for undirected or symmetric data.}
\item{diag}{Whether to use the diagonal (self-dyads).}
\item{statistic}{The statistic the permutation test uses: \code{"t"} (the
default, as Dekker et al. 2007 recommend) or \code{"beta"}, the raw
coefficient.}
\item{exact}{Whether to evaluate all \eqn{n!} permutations instead of
\code{reps} random ones (at most 10 nodes).}
\item{attr}{For netframe inputs, the edge column that fills each matrix.
\code{NULL} uses each network's weight column if it has one, and counts
edges otherwise.}
\item{tol}{Tolerance for detecting a rank-deficient design: a column that
keeps at most this fraction of its norm after the preceding columns are
projected out.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\details{
\strong{Estimates.} Each matrix is vectorized over the dyads (see
\code{mode} and \code{diag}, whose conventions are those of
\code{\link{qap_cor}} and of sna) and the response is regressed on the
predictors, with an intercept unless \code{intercept = FALSE}, by
Householder QR. Coefficients, standard errors, t-statistics, and R-squared
equal those of \code{lm()} on the vectorized dyads and of
\code{sna::netlm()}. R-squared, adjusted R-squared, and \eqn{F} follow
\code{summary.lm()}. The classical p-values of \code{lm()} are not reported,
since dyads that share a node are not independent.

\strong{Y-permutation} (\code{method = "y"}; Krackhardt 1988) permutes the
rows and columns of \code{y} together and refits the model; this is
\code{sna::netlm(nullhyp = "qapy")}. It is exact for the null hypothesis
that \code{y} is unrelated to all predictors jointly, but can be misleading
for one coefficient when predictors are correlated.

\strong{Double semi-partialling} (\code{method = "dsp"}; Dekker et al.
2007) tests each coefficient \eqn{k} by replacing predictor \eqn{k} with its
residuals from a regression on the other columns of the design (the
intercept included), permuting that residual matrix, and refitting
\code{y} on the other columns and the permuted residuals. The statistic is
the t-statistic of the permuted column; Dekker et al. show that this pivotal
statistic keeps the test valid under collinearity and skewed data. This is
\code{sna::netlm(nullhyp = "qap")} (the default there), which is the same as
\code{nullhyp = "qapspp"}, including its treatment of the intercept, which
is tested like any other column. As in sna, a model with a single column
uses Y-permutation, since there is nothing to partial out; the
\code{method} element reports the method used.

\strong{Missing values} (\code{NA}) are handled by listwise deletion of
dyads: each fit, observed or permuted, uses the dyads observed in every
matrix that enters it, and missing cells move with a matrix when it is
permuted. For DSP, the residuals of predictor \eqn{k} come from a
regression on the other predictors only, so they are defined wherever every
predictor is observed. This matches \code{sna::netlm()} for Y-permutation;
sna's DSP also matches when only \code{y} has missing values; with a
missing value in a predictor, sna's DSP fills the residual matrix from a
shorter vector by recycling, so results there are not comparable. A
permuted fit that is rank deficient, or has no more dyads than coefficients,
is left out of every proportion, and the \code{counts} report how many draws
were valid. Infinite values are an error.

\strong{Permutations and p-values} work as in \code{\link{qap_cor}}: draw
\eqn{i} uses stream \eqn{i} of \code{seed}, so results do not depend on the
number of threads; \code{exact = TRUE} evaluates every permutation once.
The p-values are the proportions of draws whose statistic \eqn{s^*}{s*}
satisfies \eqn{s^* \le s}{s* <= s} (\code{Pr(<=b)}), \eqn{s^* \ge s}{s* >= s}
(\code{Pr(>=b)}), and \eqn{|s^*| \ge |s|}{|s*| >= |s|} (\code{Pr(>=|b|)}),
where \eqn{s} is the observed statistic, as reported by
\code{sna::netlm()}. The observed statistic is not counted among the
draws, which is sna's convention; \code{(counts + 1) / (valid + 1)} gives
the form that counts it (Phipson and Smyth 2010). Values within
\eqn{10^{-10} \max(1, |s|)}{1e-10 * max(1, |s|)} of \eqn{s} count as ties.

The Monte Carlo p-values estimate the same quantities as
\code{sna::netlm()} but use a different random number generator, so they
agree with sna within Monte Carlo error, not exactly. Exact enumeration is
validated against an independent implementation for small networks.
}
\value{
An object of class \code{"mrqap"}: a list with \code{coefficients},
\code{std.error}, and \code{t}; \code{p.value}, a matrix with columns
\code{Pr(<=b)}, \code{Pr(>=b)}, and \code{Pr(>=|b|)}; \code{counts}, the
numbers of draws behind them and the number of valid draws; \code{dyads},
\code{df.residual}, \code{rss}, \code{mss}, \code{r.squared},
\code{adj.r.squared}, \code{sigma}, and \code{fstatistic}, as in
\code{summary.lm()}; the \code{method} actually used; \code{draws};
\code{dist}, a draws-by-coefficients matrix of the permuted statistics
(\code{NA} where undefined); the settings used, including \code{seed}; and
the \code{call}. \code{coef()} returns the coefficients. \code{summary()}
gives a regression table.
}
\references{
Krackhardt, D. (1988). Predicting with networks: Nonparametric multiple
regression analysis of dyadic data. \emph{Social Networks}, 10(4),
359-381. \doi{10.1016/0378-8733(88)90004-4}

Dekker, D., Krackhardt, D., & Snijders, T. A. B. (2007). Sensitivity of
MRQAP tests to collinearity and autocorrelation conditions.
\emph{Psychometrika}, 72(4), 563-581. \doi{10.1007/s11336-007-9016-1}

Krackhardt, D. (1987). QAP partialling as a test of spuriousness.
\emph{Social Networks}, 9(2), 171-186. \doi{10.1016/0378-8733(87)90012-8}

Hubert, L. J. (1987). \emph{Assignment Methods in Combinatorial Data
Analysis}. Marcel Dekker.

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero: Calculating exact p-values when permutations are randomly drawn.
\emph{Statistical Applications in Genetics and Molecular Biology}, 9(1),
Article 39. \doi{10.2202/1544-6115.1585}
}
\seealso{\code{\link{qap_cor}}}
\examples{
set.seed(1)
n <- 12
x1 <- matrix(rnorm(n * n), n, n)
x2 <- matrix(rbinom(n * n, 1, 0.4), n, n)
y <- 0.5 * x1 - 0.3 * x2 + matrix(rnorm(n * n), n, n)
fit <- mrqap(y, list(trust = x1, same_unit = x2), reps = 500)
fit
summary(fit)
}
