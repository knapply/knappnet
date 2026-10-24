\name{net_infomap}
\alias{net_infomap}
\alias{net_map_equation}
\title{Infomap communities and the map equation}
\description{
\code{net_infomap()} finds the partition that gives the shortest
description of a random walker's path on the network (Rosvall and
Bergstrom 2008), the one that minimizes the map equation;
\code{net_map_equation()} evaluates the map equation for a partition you
supply.
}
\usage{
net_infomap(x, weighted = NULL, node_weights = NULL, trials = 10L, seed = NULL)

net_map_equation(x, membership, weighted = NULL, node_weights = NULL)
}
\arguments{
\item{x}{A \code{\link{net_frame}}, directed or undirected.}
\item{weighted}{Whether to use edge weights, as given. \code{NULL} (the
default) uses them when the network has a weight column. They must not be
negative.}
\item{node_weights}{Teleportation weights of directed networks: a numeric
vector with one value per node (or the name of such a column of
\code{x$nodes}), non-negative and not all zero. \code{NULL} gives every
node the same weight. Undirected networks have no teleportation and ignore
them.}
\item{trials}{Number of independent searches, a positive whole number; the
partition of shortest description is kept.}
\item{seed}{Seed for the random node orders, a whole number between 0 and
\eqn{2^{53}}; \code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}()} makes results reproducible.}
\item{membership}{For \code{net_map_equation()}, the partition to evaluate:
a vector with one label per node (any atomic type, without missing values)
or the name of a column of \code{x$nodes}.}
}
\details{
With \eqn{p_a} the visit rate of node \eqn{a} and \eqn{q_i} the rate at
which the walker leaves group \eqn{i}, the map equation is
\deqn{L = q \log_2 q - 2 \sum_i q_i \log_2 q_i - \sum_a p_a \log_2 p_a +
\sum_i (q_i + p_i) \log_2 (q_i + p_i),}{L = plogp(q) - 2 sum_i plogp(q_i) -
sum_a plogp(p_a) + sum_i plogp(q_i + p_i),}
with \eqn{q = \sum_i q_i}, \eqn{p_i} the visit rate of group \eqn{i}, and
\eqn{\mathrm{plogp}(x) = x \log_2 x}{plogp(x) = x log2 x}; it is
measured in bits per step.

Undirected networks: the walker moves along an edge in proportion to its
weight, so \eqn{p_a} is the strength of \eqn{a} over twice the total
weight (a self-loop counts twice) and there is no teleportation. Directed
networks (edge directions are used, as in igraph's
\code{cluster_infomap()}): the walker follows an out-edge in proportion to
its weight with probability 0.85 and teleports with probability 0.15
(always from a node without out-edges), landing on node \eqn{a} in
proportion to its node weight; teleportation steps are part of the code,
as in PageRank, instead of the recorded-teleportation variant of equation 9
of Rosvall and Bergstrom (2008). Parallel edges add their weights.

The search is the core algorithm of Rosvall, Axelsson, and Bergstrom
(2009): nodes are moved to the neighboring group that most shortens the
code, groups are merged into nodes and the process repeats, and
single-node and submodule movements then refine the result until it stops
improving. A move must shorten the code by more than \eqn{10^{-10}} bits.
The search is a heuristic and its trials are random, so results differ from
other implementations in detail (the validation suite compares code
lengths); the same \code{seed} gives the same result on every platform.

An undirected network without edges has no random walk: every node is its
own group and the code length is \code{NA}.
}
\value{
\code{net_infomap()} returns an integer vector of group labels aligned with
the rows of \code{x$nodes}, numbered 1, 2, ... in order of first appearance,
with the code length of the partition in bits per step in attribute
\code{"codelength"}.

\code{net_map_equation()} returns the code length of \code{membership} in
bits per step.
}
\references{
Rosvall, M., & Bergstrom, C. T. (2008). Maps of random walks on complex
networks reveal community structure. \emph{Proceedings of the National
Academy of Sciences}, 105(4), 1118--1123. \doi{10.1073/pnas.0706851105}

Rosvall, M., Axelsson, D., & Bergstrom, C. T. (2009). The map equation.
\emph{The European Physical Journal Special Topics}, 178(1), 13--23.
\doi{10.1140/epjst/e2010-01179-1}
}
\seealso{\code{\link{net_louvain}}, \code{\link{net_plot_partition}}}
\examples{
g <- net_example("karate")
groups <- net_infomap(g, seed = 1)
groups
net_plot_partition(g, groups, seed = 1)

# The code length of any partition, here the two factions.
net_map_equation(g, "faction")
net_map_equation(g, groups)

# Directed networks use the direction of the edges.
d <- net_example("krackhardt", relation = "advice")
net_infomap(d, seed = 1)
}
