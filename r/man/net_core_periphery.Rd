\name{net_core_periphery}
\alias{net_core_periphery}
\alias{print.core_periphery}
\alias{plot.core_periphery}
\title{Core-periphery structure}
\description{
Fits the core-periphery models of Borgatti and Everett (1999), as in
UCINET's Network > Core/Periphery > Categorical and Continuous routines. In
the discrete model every node is in the core or the periphery; the ideal
pattern has every core-core pair tied and every periphery-periphery pair
untied, and the core is chosen to maximize the Pearson correlation between
the adjacency matrix and the pattern. In the continuous model every node has
a coreness \eqn{c_i \ge 0} and the pattern is \eqn{c_i c_j}, so that ties
are expected between core nodes, less so between core and peripheral nodes,
and least between peripheral nodes; coreness is fitted by least squares
(MINRES, UCINET's method) or by maximizing the correlation.
}
\usage{
net_core_periphery(x, method = c("discrete", "continuous"),
               links = c("missing", "ones", "zeros"),
               fit = c("minres", "correlation"), weighted = NULL,
               starts = 20L, max_stall = 200L, tenure = NULL, max_iter = 1000L,
               tol = 1e-10, seed = NULL, threads = 0L)

\method{print}{core_periphery}(x, ...)

\method{plot}{core_periphery}(x, what = NULL, layout = NULL, ..., seed = NULL)
}
\arguments{
\item{x}{A \code{\link{net_frame}} (for \code{print} and \code{plot}, a
\code{core_periphery} object).}
\item{method}{\code{"discrete"} (categorical core and periphery) or
\code{"continuous"} (a coreness score for every node).}
\item{links}{For the discrete model, the pattern value of core-periphery
pairs: \code{"missing"}, \code{"ones"}, or \code{"zeros"}; see Details.}
\item{fit}{For the continuous model: \code{"minres"} (least squares,
UCINET's method) or \code{"correlation"}; see Details.}
\item{weighted}{Whether to use edge weights (total weight between two nodes)
as the pair values instead of binary ties; \code{NULL} (the default) uses them
when the network has a weight column.}
\item{starts}{Number of starting points (random cores for the discrete
model, points for the continuous correlation fit).}
\item{max_stall}{Steps without improvement after which a discrete start
ends.}
\item{tenure}{Tabu tenure for the discrete search; \code{NULL} chooses
\code{max(7, n / 4)}.}
\item{max_iter}{Iteration limit per start (sweeps, for MINRES) for the
continuous model.}
\item{tol}{Convergence tolerance of MINRES: the largest change in coreness
over a sweep, relative to the largest coreness, below which sweeps stop.}
\item{seed}{Seed for the starts, a whole number between 0 and
\eqn{2^{53}}; \code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}()} makes results reproducible. For \code{plot}, the
seed of the layout.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{net_frame}}).}
\item{what}{For \code{plot}: \code{"matrix"} draws the adjacency matrix with
nodes ordered core first (by decreasing coreness for the continuous model)
and the core boundary marked (for the continuous model, between the nodes
above and below the mean coreness); \code{"network"} draws the network with
core nodes (or coreness) coloured. \code{NULL} chooses \code{"matrix"} for
the discrete and \code{"network"} for the continuous model.}
\item{layout}{For \code{plot(what = "network")}: node positions, as in
\code{\link{plot.netframe}}.}
\item{...}{For \code{plot}: further arguments passed to
\code{\link[graphics]{title}} (and, with \code{what = "network"}, to
\code{\link{plot.netframe}}), such as \code{main}. Ignored by \code{print}.}
}
\details{
Pairs of distinct nodes enter the fit: each unordered pair once in
undirected networks and both ordered pairs in directed ones. A pair's value
is 1 for a tie (parallel edges count once) or, with \code{weighted = TRUE},
the total weight of the edges between the two nodes in that direction;
weights must be non-negative. Self-loops are ignored. A network whose pairs
all have the same value (no ties, or every tie) has no core-periphery
structure and is an error.

\strong{Discrete model.} \code{links} says what the pattern expects of
core-periphery pairs: \code{"missing"} leaves them out of the correlation
(Borgatti and Everett's equation 4, the UCINET default), \code{"ones"}
expects them tied (equation 2), and \code{"zeros"} expects them untied
(equation 3). The search is tabu search (Glover 1989) over single-node moves
from \code{starts} random cores, as in \code{\link{net_factions}}; UCINET's
routine also uses a combinatorial search whose details are not published.
The best core over all starts is returned, and the validation suite checks
it against exhaustive search on small networks. The correlation is undefined
for cores whose pattern has no variance (for example a core of one node with
\code{links = "missing"} in a two-node network); such cores are never
returned, and if no core has a defined correlation the result is
\code{NaN}.

\strong{Continuous model, \code{fit = "minres"}.} Coreness minimizes the
residual sum of squares \eqn{\sum_{i \ne j} (a_{ij} - c_i c_j)^2} over
\eqn{c \ge 0}: the minimum residual method of factor analysis (Comrey 1962)
with one factor and the diagonal left out, which Borgatti and Everett
suggest and Boyd, Fitzgerald, Mahutga, and Smith (2010) describe for this
model; it is what UCINET's continuous routine computes. Each coreness in
turn is set to the value that minimizes the residuals given the others
(clamped at 0), and sweeps continue until no coreness changes by more than
\code{tol} relative to the largest, or \code{max_iter} sweeps have been
made. The start is coreness proportional to degree. On data such as a star,
whose leaves share no ties, the residual sum has no minimizer (it keeps
falling as the hub's coreness grows); the sweeps then stop at
\code{max_iter}, and the unit-norm coreness, which still converges, is
reported.

\strong{Continuous model, \code{fit = "correlation"}.} Coreness maximizes
the Pearson correlation between the data and the pattern, the criterion
Borgatti and Everett state, by the quasi-Newton method of Broyden,
Fletcher, Goldfarb, and Shanno with an Armijo line search (Nocedal and
Wright 2006, chapters 3 and 6) on the exact gradient, over \eqn{t} with
\eqn{c = t^2} so that coreness stays non-negative (without the constraint
the product pattern also fits two separate dense groups with coreness of
opposite signs). The correlation has several local maxima in general, so
the search runs from \code{starts} points (degree plus one, then random
points from the seeded generator) and keeps the best. On many networks the
supremum is only approached as all coreness values become equal, where the
pattern reduces to the additive \eqn{c_i + c_j} up to scale; the returned
coreness is then nearly uniform and the correlation is that of the best
additive fit. MINRES does not have this degeneracy, which is why it is the
default.

Both fits are scaled to unit Euclidean norm (the correlation does not change
under scaling, and UCINET reports MINRES coreness so). The validation suite
checks MINRES against scipy's bounded quasi-Newton minimizer and the
correlation fit against scipy's BFGS from many starts. Starts run in
parallel, and the result for a given \code{seed} does not depend on
\code{threads}.
}
\value{
A list of class \code{core_periphery} with \code{method}; for the discrete
model \code{links} and \code{core}, a logical vector aligned with the rows of
\code{x$nodes} (\code{TRUE} for core nodes); for the continuous model
\code{fit}, \code{coreness}, a numeric vector aligned with \code{x$nodes}
(non-negative, unit Euclidean norm), \code{residual} (the MINRES residual
sum of squares, \code{NaN} for the correlation fit), and \code{iterations};
\code{correlation}, the correlation between the data and the fitted
pattern; and \code{network}, the netframe fitted (used by \code{plot}).
\code{plot} returns \code{x} invisibly.
}
\references{
Borgatti, S. P., & Everett, M. G. (1999). Models of core/periphery
structures. \emph{Social Networks}, 21(4), 375-395.
\doi{10.1016/S0378-8733(99)00019-2}

Boyd, J. P., Fitzgerald, W. J., Mahutga, M. C., & Smith, D. A. (2010).
Computing continuous core/periphery structures for social relations data
with MINRES/SVD. \emph{Social Networks}, 32(2), 125-137.
\doi{10.1016/j.socnet.2009.09.003}

Comrey, A. L. (1962). The minimum residual method of factor analysis.
\emph{Psychological Reports}, 11(1), 15-18. \doi{10.2466/pr0.1962.11.1.15}

Glover, F. (1989). Tabu search, part I. \emph{ORSA Journal on Computing},
1(3), 190-206. \doi{10.1287/ijoc.1.3.190}

Nocedal, J., & Wright, S. J. (2006). \emph{Numerical Optimization} (2nd
ed.). Springer. \doi{10.1007/978-0-387-40065-5}
}
\seealso{\code{\link{net_factions}}, \code{\link{net_coreness}} (k-cores, a
different notion)}
\examples{
# Node 2 ties to most others, 5 and 7 to 2 and each other; the rest hang on.
g <- net_frame(data.frame(source = c(1, 1, 2, 2, 2, 2, 3, 4, 5),
                         target = c(2, 7, 4, 5, 6, 7, 5, 5, 7)), directed = FALSE)
cp <- net_core_periphery(g, seed = 1, threads = 2)
cp
g$nodes$name[cp$core]
plot(cp)

# The three treatments of core-periphery pairs give different cores.
net_core_periphery(g, links = "missing", seed = 1, threads = 2)$core  # UCINET default (eq. 4)
net_core_periphery(g, links = "ones", seed = 1, threads = 2)$core     # expected tied (eq. 2)
net_core_periphery(g, links = "zeros", seed = 1, threads = 2)$core    # expected untied (eq. 3)

# Continuous coreness: UCINET's least-squares fit, and the correlation fit,
# which here collapses toward equal coreness.
karate <- net_example("karate")
cc <- net_core_periphery(karate, method = "continuous", threads = 2)        # MINRES (UCINET)
round(sort(cc$coreness, decreasing = TRUE)[1:5], 3)
cr <- net_core_periphery(karate, method = "continuous", fit = "correlation",
                     seed = 1, threads = 2)                            # maximum correlation
round(range(cr$coreness), 3)
plot(cc, seed = 1)
}
