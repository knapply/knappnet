\name{flow_betweenness}
\alias{flow_betweenness}
\title{Flow betweenness}
\usage{
flow_betweenness(x, scale = c("raw", "normalized", "fractional"),
                 capacity = NULL, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{scale}{How scores are scaled; see Details. \code{"raw"} and
\code{"normalized"} are Freeman, Borgatti, and White's measure as UCINET
reports it (and sna's \code{cmode = "rawflow"} and \code{"normflow"});
\code{"fractional"} is sna's \code{cmode = "fracflow"}.}
\item{capacity}{Whether edge weights are capacities. \code{NULL} (the
default) uses them when the network's weights are strengths
(\code{weight_type = "strength"} in \code{\link{netframe}}), since tie
strength is how much can flow through a tie; distances are not capacities,
so networks weighted by distance use unit capacities unless
\code{capacity = TRUE}.}
\item{threads}{Number of threads; \code{0} uses all available cores (see
"Threads" in \code{\link{netframe}}).}
}
\value{
A numeric vector aligned with the rows of \code{x$nodes}.
}
\description{
Freeman, Borgatti, and White's (1991) flow betweenness: how much of the
maximum flow between every pair of other nodes must pass through each
node. Unlike shortest-path \code{\link{betweenness}()}, it credits nodes on
every route between a pair, not only the shortest ones, which suits
information or resources that travel by any available path, as in covert
networks that route around exposed members.
}
\details{
Let \eqn{f(i, j)} be the maximum flow from \eqn{i} to \eqn{j} (the
largest number of edge-disjoint paths, without capacities) and
\eqn{f_v(i, j)} the maximum flow once \eqn{v} is removed. For each node
\eqn{v}, over pairs of other nodes:
\describe{
  \item{\code{"raw"}}{\eqn{\sum_{i, j} (f(i, j) - f_v(i, j))}: the part of
    each pair's maximum flow that cannot avoid \eqn{v}.}
  \item{\code{"normalized"}}{The raw score divided by
    \eqn{\sum_{i, j} f(i, j)}, the total maximum flow between pairs not
    involving \eqn{v} (UCINET reports this as a percentage);
    \code{NaN} when that total is 0.}
  \item{\code{"fractional"}}{\eqn{\sum_{i, j} (f(i, j) - f_v(i, j)) / f(i, j)}
    over pairs with \eqn{f(i, j) > 0} (Koschutzki et al. 2005): the share
    of each pair's flow that depends on \eqn{v}, summed, in the manner of
    shortest-path betweenness.}
}
Pairs are ordered in directed networks and unordered in undirected ones,
as in sna (\code{gmode = "graph"}).

Conventions: without capacities every edge carries one unit; with them,
edge weights are capacities and must be non-negative. Parallel edges add
their capacities. An undirected edge carries flow in either direction up to
its capacity; a directed one only along its arc. Nodes have unlimited
capacity, and self-loops carry nothing.

Algorithm: maximum flows use Dinic's (1970) algorithm. A node that
carries no flow in a maximum flow for a pair loses nothing when removed,
since the same flow avoids it, so only nodes carrying flow are rerun with
the node removed. The cost is still up to \eqn{O(n^3)} maximum flows in
the worst case, run in parallel over source nodes; networks of a few
hundred nodes take seconds to minutes.

Validation: matches sna's \code{flowbet()} for all three scales, and a
brute-force evaluation of the definition (one maximum flow per pair and
removed node, by the Edmonds-Karp algorithm in exact arithmetic). sna
reads an adjacency matrix, so on networks with parallel edges it uses one
tie per pair where this function adds the capacities.
}
\references{
Freeman, L. C., Borgatti, S. P., & White, D. R. (1991). Centrality in
valued graphs: A measure of betweenness based on network flow.
\emph{Social Networks}, 13(2), 141-154. \doi{10.1016/0378-8733(91)90017-N}

Koschutzki, D., Lehmann, K. A., Peeters, L., Richter, S.,
Tenfelde-Podehl, D., & Zlotowski, O. (2005). Centrality indices. In U.
Brandes & T. Erlebach (Eds.), \emph{Network Analysis: Methodological
Foundations} (Lecture Notes in Computer Science 3418, pp. 16-61).
Springer. \doi{10.1007/978-3-540-31955-9_3}

Dinic, E. A. (1970). Algorithm for solution of a problem of maximum flow
in a network with power estimation. \emph{Soviet Mathematics Doklady}, 11,
1277-1280. Translated from \emph{Doklady Akademii Nauk SSSR}, 194(4),
754-757.

Ford, L. R., & Fulkerson, D. R. (1956). Maximal flow through a network.
\emph{Canadian Journal of Mathematics}, 8, 399-404.
\doi{10.4153/CJM-1956-045-5}

Everton, S. F. (2012). \emph{Disrupting Dark Networks}. Cambridge
University Press. \doi{10.1017/CBO9781139136877}
}
\seealso{
\code{\link{betweenness}}, \code{\link{key_players}},
\code{\link{group_centrality}}
}
\examples{
# A hypothetical covert network: two cells whose leaders (a1, b1) are
# linked through a courier, and a financier tied to both cells.
cells <- data.frame(
  from = c("a1", "a1", "a2", "b1", "b1", "b2", "a1", "courier", "fin", "fin"),
  to   = c("a2", "a3", "a3", "b2", "b3", "b3", "courier", "b1", "a2", "b2")
)
g <- netframe(cells, directed = FALSE)

# Freeman, Borgatti, and White's totals (UCINET; sna cmode = "rawflow").
flow_betweenness(g, threads = 1)
# Divided by the flow between other pairs (UCINET's normalized scores, as
# proportions; sna cmode = "normflow").
flow_betweenness(g, scale = "normalized", threads = 1)
# Shares of each pair's flow, summed (sna cmode = "fracflow").
flow_betweenness(g, scale = "fractional", threads = 1)

# Shortest-path betweenness credits only the shortest routes.
betweenness(g)
}
