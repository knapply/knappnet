\name{net_small_world}
\alias{net_small_world}
\alias{print.small_world}
\alias{plot.small_world}
\title{Small-world quotient}
\description{
Watts and Strogatz (1998) call a network a small world when its clustering
is much higher than that of a random network of the same size and density
while its mean distance is about as short. The small-world quotient (Uzzi
and Spiro 2005; the "small-world-ness" of Humphries and Gurney 2008)
combines the two ratios:
\deqn{Q = \frac{C / C_{random}}{L / L_{random}}}{Q = (C / C_random) / (L / L_random)}
Here \eqn{C_{random}} and \eqn{L_{random}} are the means over \code{reps}
simulated \eqn{G(n, m)} networks with the same number of nodes and edges
(see \code{\link{net_random_gnm}}), rather than the approximations
\eqn{C_{random} \approx k/n} and \eqn{L_{random} \approx \ln n / \ln k}
used in some textbooks, which fail for small or sparse networks.
}
\usage{
net_small_world(x, local = FALSE, reps = 100L, seed = NULL, threads = 0L)

\method{print}{small_world}(x, ...)

\method{plot}{small_world}(x, ...)
}
\arguments{
\item{x}{A \code{\link{net_frame}} (for \code{print} and \code{plot}, a
\code{small_world} object).}
\item{local}{Whether to use the mean local clustering coefficient
(\code{TRUE}) or global transitivity (\code{FALSE}).}
\item{reps}{Number of random networks.}
\item{seed}{Seed for the random networks, a whole number between 0 and
\eqn{2^{53}}; \code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}()} makes results reproducible.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{net_frame}}).}
\item{...}{For \code{plot}, \code{main}; otherwise ignored.}
}
\details{
Everything is computed on the undirected simple graph underlying the
network (directions ignored, parallel edges once, self-loops ignored),
without weights, and the random networks are undirected and simple with as
many edges as that graph. \eqn{C} is the global transitivity (three times
the triangles over the connected triples; \code{local = FALSE}, igraph's
\code{net_transitivity(type = "global")}) or the mean local clustering
coefficient of Watts and Strogatz with nodes of fewer than two neighbors
counted as 0 (\code{local = TRUE}, NetworkX's \code{net_avg_clustering()}
and igraph's \code{net_transitivity(type = "average", isolates = "zero")}).
\eqn{L} is the mean distance over pairs of nodes joined by a path (igraph's
\code{net_mean_distance(unconnected = TRUE)}). Parts that are undefined (no
connected triples, no connected pairs) give \code{NaN}.

Replicate \eqn{r} uses random stream \code{(seed, r)} and replicates run in
parallel, so the result for a given \code{seed} does not depend on
\code{threads}. The observed values are checked by the validation suite
against igraph and NetworkX; the random means depend on the seed and are
checked for reproducibility and against the known expectation of
\eqn{G(n, m)} (clustering near the density).
}
\value{
A list of class \code{small_world} with \code{clustering} and
\code{mean_distance} of the network, \code{random_clustering} and
\code{random_mean_distance} (means over the random networks),
\code{quotient}, \code{reps}, and \code{local}. \code{plot} draws the
observed and random clustering and mean distance side by side and returns
\code{x} invisibly.
}
\references{
Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
'small-world' networks. \emph{Nature}, 393(6684), 440-442.
\doi{10.1038/30918}

Uzzi, B., & Spiro, J. (2005). Collaboration and creativity: The small world
problem. \emph{American Journal of Sociology}, 111(2), 447-504.
\doi{10.1086/432782}

Humphries, M. D., & Gurney, K. (2008). Network 'small-world-ness': A
quantitative method for determining canonical network equivalence.
\emph{PLoS ONE}, 3(4), e0002051. \doi{10.1371/journal.pone.0002051}
}
\seealso{\code{\link{net_transitivity}}, \code{\link{net_mean_distance}},
\code{\link{net_random_gnm}}}
\examples{
g <- net_example("karate")
sw <- net_small_world(g, reps = 50, seed = 1, threads = 2)
sw
plot(sw)

# Global transitivity (igraph's default) and mean local clustering
# (NetworkX's average_clustering) give different quotients.
net_small_world(g, local = FALSE, reps = 50, seed = 1, threads = 2)$quotient  # transitivity
net_small_world(g, local = TRUE, reps = 50, seed = 1, threads = 2)$quotient   # average clustering
}
