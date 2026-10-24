\name{net_ego_summary}
\alias{net_ego_summary}
\alias{net_ego_size}
\alias{net_ego_density}
\alias{net_ego_betweenness}
\title{Egocentric measures}
\description{
Measures of every node's ego network, computed together in parallel and
returned aligned with the rows of \code{x$nodes}. \code{net_ego_summary()} gives
them all in a data frame; \code{net_ego_size()}, \code{net_ego_density()}, and
\code{net_ego_betweenness()} return one measure as a vector.
}
\usage{
net_ego_summary(x, order = 1L, mode = c("all", "out", "in"), betweenness = TRUE,
  threads = 0L)

net_ego_size(x, order = 1L, mode = c("all", "out", "in"), include_ego = TRUE,
  threads = 0L)

net_ego_density(x, order = 1L, mode = c("all", "out", "in"), threads = 0L)

net_ego_betweenness(x, order = 1L, mode = c("all", "out", "in"), threads = 0L)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{order}{Alters are the nodes within this many steps of ego. \code{1}
gives the first-order neighborhood.}
\item{mode}{For directed networks, \code{"all"}, \code{"out"}, or
\code{"in"}; see \code{\link{net_ego}}.}
\item{betweenness}{Whether to compute ego betweenness, the costliest of the
measures.}
\item{include_ego}{Whether ego itself is counted (igraph's
\code{mindist = 0}) or only its alters (\code{mindist = 1}).}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{net_frame}}).}
}
\details{
For a node with \eqn{k} alters (the nodes within \code{order} steps; see
\code{\link{net_ego}}):
\describe{
  \item{\code{size}}{\eqn{k}, the number of alters.}
  \item{\code{ties}}{The number of ties among the alters, leaving out ego's
    own ties.}
  \item{\code{pairs}}{The number of pairs of alters: \eqn{k(k-1)/2}, or
    \eqn{k(k-1)} ordered pairs in a directed network.}
  \item{\code{density}}{\code{ties / pairs}, the density of the ego network
    without ego; \code{NaN} for fewer than two alters.}
  \item{\code{components}}{The number of weak components among the alters
    once ego is removed; isolated alters count as components.}
  \item{\code{ego_betweenness}}{The betweenness of ego within its own ego
    network (Everett and Borgatti 2005; Marsden 2002): for every pair of
    alters, the share of the shortest paths between them, inside the ego
    network, that pass through ego. Each unordered pair counts once in an
    undirected network and each ordered pair in a directed one. Not
    normalized; the largest possible value is \code{pairs}.}
}

\code{net_ego_size()} is igraph's \code{net_ego_size()}: with
\code{include_ego = TRUE} (the default, igraph's \code{mindist = 0}) it
counts ego as well, so it is \code{size + 1}.

Ties are binary: weights are ignored, parallel edges count once, and
self-loops are ignored, as in UCINET's ego network measures. In directed
networks \code{mode} chooses the alters (see \code{\link{net_ego}}),
every tie among them counts with its direction, and ego betweenness follows
directed paths. Burt's structural-hole measures, which are also egocentric,
are \code{\link{net_constraint}}, \code{\link{net_effective_size}}, and
\code{\link{net_efficiency}}.

Results match igraph (\code{net_ego_size()}, and \code{ecount()},
\code{edge_density()}, \code{count_components()}, and \code{net_betweenness()}
of \code{induced_subgraph(simplify(g), ego(g, ...))}) and NetworkX
(\code{ego_graph()} on the simple graph).
}
\value{
\code{net_ego_summary()} returns a data frame with one row per node, in the
order of \code{x$nodes}, and columns \code{node} (the identifier),
\code{size}, \code{ties}, \code{pairs}, \code{density}, \code{components},
and (if requested) \code{ego_betweenness}. The other functions return a
numeric vector aligned with the rows of \code{x$nodes}.
}
\references{
Everett, M., & Borgatti, S. P. (2005). Ego network betweenness.
\emph{Social Networks}, 27(1), 31-38. \doi{10.1016/j.socnet.2004.11.007}

Marsden, P. V. (2002). Egocentric and sociocentric measures of network
centrality. \emph{Social Networks}, 24(4), 407-422.
\doi{10.1016/S0378-8733(02)00016-3}

Perry, B. L., Pescosolido, B. A., & Borgatti, S. P. (2018).
\emph{Egocentric Network Analysis: Foundations, Methods, and Models}.
Cambridge University Press. \doi{10.1017/9781316443255}
}
\examples{
g <- net_frame(data.frame(source = c("a", "a", "a", "b", "d"),
                         target   = c("b", "c", "d", "c", "e")),
              directed = FALSE)
net_ego_summary(g)
net_ego_size(g, order = 2)
}
