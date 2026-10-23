\name{cug_test}
\alias{cug_test}
\alias{print.cug_test}
\alias{summary.cug_test}
\alias{print.summary.cug_test}
\alias{plot.cug_test}
\title{Conditional uniform graph tests}
\description{
Compares a network statistic with its distribution over random networks
that share some feature of the observed network (Anderson, Butts, and
Carley 1999), as sna's \code{cug.test()} does. The null model is uniform
over networks with the same number of nodes (\code{conditioning = "size"},
that is \eqn{G(n, 1/2)}), the same number of edges (\code{"edges"},
\eqn{G(n, m)}), or the same dyad census (\code{"dyad.census"}, the U|MAN
distribution of Holland and Leinhardt 1976: the same numbers of mutual,
asymmetric, and null dyads).
}
\usage{
cug_test(x, statistic, conditioning = c("size", "edges", "dyad.census"),
         reps = 1000L, ..., seed = NULL, threads = 0L)

\method{print}{cug_test}(x, ...)

\method{summary}{cug_test}(object, ...)

\method{plot}{cug_test}(x, ..., breaks = "Sturges")
}
\arguments{
\item{x}{A \code{\link{netframe}} (for \code{print} and \code{plot}, a
\code{cug_test} object).}
\item{statistic}{A function of a netframe returning one number.}
\item{conditioning}{What the random networks share with the observed one:
\code{"size"}, \code{"edges"}, or \code{"dyad.census"}.}
\item{reps}{Number of random networks.}
\item{...}{Further arguments passed to \code{statistic}. For \code{plot},
arguments passed to \code{\link[graphics]{plot.histogram}} and on to
\code{\link[graphics]{title}}, such as \code{main}; ignored by
\code{print} and \code{summary}.}
\item{seed}{Seed for the random networks, a whole number between 0 and
\eqn{2^{53}}; \code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}()} makes results reproducible.}
\item{threads}{Number of threads for drawing the random networks, a whole
number; \code{0} uses all available cores, or 2 during \verb{R CMD check} (see
\sQuote{Threads} in \code{\link{netframe}}).}
\item{object}{A \code{cug_test}.}
\item{breaks}{Breakpoints for the histogram, as in
\code{\link[graphics]{hist}}.}
}
\details{
\code{statistic} is any function that takes a netframe (and \code{...}) and
returns a single number, such as \code{\link{transitivity}},
\code{\link{centralization}}, or \code{function(g) mean(degree(g))}. It is
applied to the observed network as given and to \code{reps} random networks
on \code{1, 2, ..., n} nodes, drawn in parallel with the generators of
\code{\link{random_gnp}}, \code{\link{random_gnm}}, and
\code{\link{random_uman}} from a seeded counter-based generator, so the
null distribution for a given \code{seed} does not depend on
\code{threads}; the statistic itself runs in R, one network at a time. All
random networks are drawn before the statistic runs, so memory grows with
\code{reps} times the number of edges.

The random networks are simple and directed if and only if the observed
network is. The observed network's number of edges and dyad census are
counted on its simple reading (parallel edges once, self-loops ignored), as
sna does on an adjacency matrix with \code{diag = FALSE}; a network with
self-loops or parallel edges gets a warning, since the statistic sees them
while the null networks never have them. In an undirected network every tie
is a mutual dyad, so conditioning on the dyad census is the same as
conditioning on the number of edges and is treated so. Weights and
attributes are not carried to the random networks.

The p-values are the proportions of random networks whose statistic is at
least (\code{p_greater}) and at most (\code{p_less}) the observed one, as
sna reports them (\code{Pr(X >= Obs)} and \code{Pr(X <= Obs)}). A
statistic that is \code{NA} or \code{NaN} on some random networks is
dropped from the distribution for those networks, and the number dropped is
reported. The null generators are checked exactly by the validation suite
(\code{validation/cug_check.py}: \eqn{G(n, m)} edge counts, U|MAN dyad
censuses, configuration-model degrees), and the p-values against sna's
\code{cug.test()} on the same statistics within Monte Carlo error.
}
\value{
An object of class \code{cug_test}: a list with \code{observed} (the
statistic on \code{x}), \code{null} (its values on the random networks),
\code{p_greater} and \code{p_less}, \code{conditioning}, \code{reps},
\code{dropped} (random networks on which the statistic was undefined), and
\code{statistic} (the expression given). \code{summary()} adds the mean,
standard deviation, and quantiles of the null distribution. \code{plot()}
draws a histogram of the null distribution with the observed value marked
and the p-values in a legend, and returns \code{x} invisibly.
}
\references{
Anderson, B. S., Butts, C., & Carley, K. (1999). The interaction of size
and density with graph-level indices. \emph{Social Networks}, 21(3),
239-267. \doi{10.1016/S0378-8733(99)00011-8}

Holland, P. W., & Leinhardt, S. (1976). Local structure in social networks.
\emph{Sociological Methodology}, 7, 1-45. \doi{10.2307/270703}
}
\seealso{\code{\link{random_gnp}}, \code{\link{qap_cor}} for permutation
tests of association between networks}
\examples{
g <- example_network("karate")
# Is the karate club more clustered than random networks of its size and
# density?
ct <- cug_test(g, transitivity, conditioning = "edges", reps = 200, seed = 1, threads = 2)
ct
summary(ct)
plot(ct)

# The three null models: conditioning only on size (density 1/2) makes
# almost any sparse network look unusual, conditioning on edges or on the
# dyad census (the same thing for undirected networks) is the usual choice.
cug_test(g, transitivity, "size", reps = 200, seed = 1, threads = 2)$p_greater
cug_test(g, transitivity, "edges", reps = 200, seed = 1, threads = 2)$p_greater
cug_test(g, transitivity, "dyad.census", reps = 200, seed = 1, threads = 2)$p_greater

# A directed network: the dyad census fixes reciprocity, so a statistic
# driven by mutual ties is no longer unusual under it.
d <- random_uman(15, mutual = 20, asym = 10, seed = 2, threads = 2)
cug_test(d, reciprocity, "edges", reps = 200, seed = 1, threads = 2)$p_greater
cug_test(d, reciprocity, "dyad.census", reps = 200, seed = 1, threads = 2)$p_greater

# Extra arguments go to the statistic.
cug_test(g, centralization, "edges", measure = "degree", reps = 100, seed = 1, threads = 2)
}
