\name{qap_cor}
\alias{qap_cor}
\alias{print.qap_cor}
\alias{summary.qap_cor}
\alias{print.summary.qap_cor}
\title{QAP correlation test}
\usage{
qap_cor(y, x, reps = 1000, seed = NULL, diag = FALSE,
        mode = c("digraph", "graph"), exact = FALSE, attr = NULL,
        threads = 0L)
}
\arguments{
\item{y, x}{Square numeric matrices of the same size, or
\code{\link{netframe}}s (see \code{attr}).}
\item{reps}{Number of random permutations.}
\item{seed}{Seed for the permutations, a whole number in [0, 2^53).
\code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}} makes results reproducible.}
\item{diag}{Whether to use the diagonal (self-dyads).}
\item{mode}{\code{"digraph"} uses every ordered pair of distinct nodes;
\code{"graph"} uses each unordered pair once, read from the lower triangle,
for undirected or symmetric data.}
\item{exact}{Whether to evaluate all \eqn{n!} permutations instead of
\code{reps} random ones (at most 10 nodes).}
\item{attr}{For netframe inputs, the edge column that fills the matrix.
\code{NULL} uses the network's weight column if it has one, and counts
edges otherwise.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\value{
An object of class \code{"qap_cor"}: a list with the observed correlation
\code{statistic}; \code{p.value}, the three proportions described below;
\code{counts}, the numbers of draws behind them; \code{valid}, the number of
draws with a defined correlation; \code{draws}; \code{dyads}, the number of
dyads used; \code{dist}, the correlation of every draw in order; the
settings used, including \code{seed}; and the \code{call}. \code{summary()}
adds quantiles of the reference distribution.
}
\description{
Tests the Pearson correlation between two networks, or any two
\eqn{n \times n}{n x n} matrices of dyadic data, with the quadratic
assignment procedure (QAP; Hubert and Schultz 1976; Hubert 1987;
Krackhardt 1987). The reference distribution is formed by relabeling the
nodes of \code{x}: its rows and columns are permuted together, which keeps
the dependence among dyads that share a node.
}
\details{
\strong{Dyads.} With \code{mode = "digraph"}, every ordered pair
\eqn{(i, j)}, \eqn{i \ne j}{i != j}, is a dyad. With \code{mode = "graph"},
each unordered pair is used once, read from the lower triangle; the upper
triangle is ignored. The diagonal is used only with \code{diag = TRUE}.
These are the conventions of \code{sna::gcor()}. A netframe becomes its
adjacency matrix (see \code{\link{as.matrix.netframe}}), so absent ties are
0; an undirected netframe gives a symmetric matrix, for which
\code{mode = "graph"} is appropriate. When every matrix has row and column
names, they are aligned by name.

\strong{Missing values.} \code{NA} cells are handled by listwise deletion:
each correlation, observed or permuted, uses the dyads observed in both
matrices. Missing cells move with \code{x} when it is permuted, so the dyads
used can differ between draws. This matches \code{sna::gcor()}, which uses
\code{use = "complete.obs"}. Infinite values are an error.

\strong{Permutations.} Draw \eqn{i} is the permutation produced by stream
\eqn{i} of a counter-based random number generator (Philox) seeded with
\code{seed}, so results are identical for any number of threads. Draws are
uniform over all \eqn{n!} permutations and may repeat. With
\code{exact = TRUE}, every permutation is evaluated exactly once (including
the identity, the observed labeling), giving the exact permutation
distribution.

\strong{p-values} are proportions of the draws, in the forms sna reports:
\code{Pr(<=r)} is the proportion with \eqn{r^* \le r}{r* <= r},
\code{Pr(>=r)} the proportion with \eqn{r^* \ge r}{r* >= r}, and
\code{Pr(>=|r|)} the proportion with \eqn{|r^*| \ge |r|}{|r*| >= |r|}, where
\eqn{r} is observed and \eqn{r^*}{r*} a draw. The observed statistic is not
counted among the draws, which is sna's convention; with \code{reps}
random draws, \code{(counts + 1) / (valid + 1)} gives the form that counts
it (Phipson and Smyth 2010), which is never zero. Values within
\eqn{10^{-10} \max(1, |r|)}{1e-10 * max(1, |r|)} of \eqn{r} count as ties,
so that relabelings giving the same correlation in exact arithmetic are not
split by rounding. Draws whose correlation is undefined (a constant
permuted vector after deletion) are left out of every proportion.

The correlation matches \code{sna::gcor()} exactly. The p-values estimate
the same quantities as \code{sna::qaptest(dat, gcor, g1 = 1, g2 = 2)},
which permutes both graphs independently; since only the relative labeling
matters, the distributions are the same, but individual Monte Carlo draws
differ. Exact enumeration is validated against an independent
implementation for small networks.
}
\references{
Hubert, L., & Schultz, J. (1976). Quadratic assignment as a general data
analysis strategy. \emph{British Journal of Mathematical and Statistical
Psychology}, 29(2), 190-241. \doi{10.1111/j.2044-8317.1976.tb00714.x}

Hubert, L. J. (1987). \emph{Assignment Methods in Combinatorial Data
Analysis}. Marcel Dekker.

Krackhardt, D. (1987). QAP partialling as a test of spuriousness.
\emph{Social Networks}, 9(2), 171-186. \doi{10.1016/0378-8733(87)90012-8}

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero: Calculating exact p-values when permutations are randomly drawn.
\emph{Statistical Applications in Genetics and Molecular Biology}, 9(1),
Article 39. \doi{10.2202/1544-6115.1585}
}
\seealso{\code{\link{mrqap}}}
\examples{
set.seed(1)
y <- matrix(rbinom(100, 1, 0.3), 10, 10)
x <- y + matrix(rnorm(100), 10, 10)
qap_cor(y, x, reps = 500)
summary(qap_cor(y[1:6, 1:6], x[1:6, 1:6], exact = TRUE))
}
