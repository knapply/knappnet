\name{twomode_degree}
\alias{twomode_degree}
\alias{twomode_closeness}
\alias{twomode_betweenness}
\alias{twomode_eigenvector}
\alias{twomode_density}
\title{Centrality and density of two-mode networks}
\description{
Degree, closeness, betweenness, and eigenvector centrality of the nodes of a
two-mode network, normalized as Borgatti and Everett (1997) and Borgatti and
Halgin (2011) describe, and two-mode density.
}
\usage{
twomode_degree(x, normalized = TRUE)

twomode_closeness(x, normalized = TRUE,
  unreachable = c("wasserman_faust", "ignore"), threads = 0L)

twomode_betweenness(x, normalized = TRUE, threads = 0L)

twomode_eigenvector(x, tol = 1e-12, max_iter = 10000L, threads = 0L)

twomode_density(x)
}
\arguments{
\item{x}{A two-mode \code{\link{netframe}} (see \code{\link{twomode}}).}
\item{normalized}{Whether to apply the two-mode normalization.}
\item{unreachable}{How closeness treats nodes that cannot be reached:
\code{"wasserman_faust"} or \code{"ignore"}; see Details.}
\item{tol}{Convergence tolerance for eigenvector scores: iteration stops when
no score (scaled so the largest is 1) changes by more than \code{tol}.}
\item{max_iter}{Iteration cap. Not converging within it is an error.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\details{
The scores are those of the bipartite graph; only the normalization differs
from the one-mode case, because a node can be tied only to the other mode.
For a node whose own mode has \eqn{n_o} nodes, with \eqn{n_i} nodes in the
other mode:
\describe{
  \item{degree}{\eqn{d / n_i}.}
  \item{closeness}{\eqn{(n_i + 2(n_o - 1)) / S}, where \eqn{S} is the sum of
    the node's distances to the others; the numerator is the smallest
    possible sum (distance 1 to every node of the other mode, 2 to every
    other node of its own).}
  \item{betweenness}{\eqn{b / b_{max}}{b / b_max} with
    \eqn{b_{max} = [n_i^2 (s+1)^2 + n_i (s+1)(2t - s - 1) - t(2s - t + 3)] / 2}{b_max = [n_i^2 (s+1)^2 + n_i (s+1)(2t - s - 1) - t(2s - t + 3)] / 2},
    \eqn{s = (n_o - 1) \;\mathrm{div}\; n_i}{s = (n_o - 1) div n_i} and
    \eqn{t = (n_o - 1) \bmod n_i}{t = (n_o - 1) mod n_i}, the largest
    betweenness a node can have given the sizes of the two modes (Borgatti
    and Everett 1997). The raw betweenness counts each unordered pair
    once.}
  \item{eigenvector}{the principal eigenvector of the bipartite graph, which
    is the pair of leading singular vectors of the incidence matrix
    (Borgatti and Everett 1997; Borgatti and Halgin 2011). The two modes'
    halves of the unit eigenvector have the same length, so Borgatti and
    Everett's normalization (dividing by \eqn{\sqrt{1/2}}{sqrt(1/2)}) gives
    each mode Euclidean length 1, which is what is returned. Their further
    two-mode normalization of eigenvector scores has no closed form and is
    not provided.}
  \item{density}{the number of ties divided by \eqn{n_1 n_2}, the number of
    possible ties between the modes.}
}
With \code{normalized = FALSE} degree, betweenness, and closeness are the
raw scores \eqn{d}, \eqn{b}, and \eqn{1/S}.

Ties are binary and undirected: weights are ignored, parallel edges count
once, and direction does not matter. When the network is not connected,
\eqn{S} sums the distances to the nodes a node reaches. With
\code{unreachable = "wasserman_faust"} (the default) closeness is then
multiplied by \eqn{r/(n-1)}, the share of the other nodes it reaches, and a
node that reaches no other scores 0; with \code{"ignore"} it is not
rescaled and a node that reaches no other gets \code{NaN}. Eigenvector
scores are then unique only within the component with the largest
eigenvalue (a warning is given); every score is \code{NaN} in a network
without ties. Normalizations that divide by zero (an empty mode, or
\eqn{b_{max} = 0}{b_max = 0} when the other mode has one node) give
\code{NaN}.

Degree, closeness, and betweenness match NetworkX's
\code{bipartite.degree_centrality()}, \code{bipartite.closeness_centrality()}
(\code{normalized = TRUE} for the Wasserman-Faust rule, \code{FALSE} for
\code{"ignore"}), and \code{bipartite.betweenness_centrality()}, and density
matches \code{bipartite.density()}. Raw scores match igraph's
\code{degree()}, \code{closeness()}, and \code{betweenness()} on the
simplified network, and eigenvector scores match igraph's
\code{eigen_centrality()} rescaled within each mode. UCINET reports the same
normalized scores as percentages.
}
\value{
\code{twomode_density()} returns a number; the others return a numeric
vector aligned with the rows of \code{x$nodes}.
}
\references{
Borgatti, S. P., & Everett, M. G. (1997). Network analysis of 2-mode data.
\emph{Social Networks}, 19(3), 243-269.
\doi{10.1016/S0378-8733(96)00301-2}

Borgatti, S. P., & Halgin, D. S. (2011). Analyzing affiliation networks. In
J. Scott & P. J. Carrington (Eds.), \emph{The SAGE Handbook of Social
Network Analysis} (pp. 417-433). SAGE. \doi{10.4135/9781446294413.n28}

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods
and Applications}. Cambridge University Press.
\doi{10.1017/CBO9780511815478}
}
\examples{
attended <- data.frame(
  person = c("ann", "ann", "bob", "bob", "cat", "dan"),
  event  = c("gala", "fair", "gala", "show", "show", "fair")
)
g <- twomode(attended)
data.frame(g$nodes,
           degree = twomode_degree(g),
           closeness = twomode_closeness(g),
           betweenness = twomode_betweenness(g))
twomode_density(g)
}
