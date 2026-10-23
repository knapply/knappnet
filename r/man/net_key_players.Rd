\name{net_key_players}
\alias{net_key_players}
\alias{net_key_player_score}
\title{Key players}
\description{
Borgatti's (2006) key player problems. KPP-Neg asks which \eqn{k} nodes,
if removed, would most fragment the network: the question behind
disrupting a covert network by arresting or isolating its members. KPP-Pos
asks which \eqn{k} nodes are best placed to reach the rest of it: whom to
recruit, monitor, or inform through. \code{net_key_players()} searches for the
best set; \code{net_key_player_score()} evaluates any set.
}
\usage{
net_key_players(x, k, measure = c("fragmentation", "distance_fragmentation",
                              "m_reach", "distance_reach"),
            m = 1, mode = c("out", "in", "all"), weighted = NULL,
            weights_as = NULL, method = c("greedy", "exact"), restarts = 10L,
            seed = NULL, threads = 0L)

net_key_player_score(x, set, measure = c("fragmentation", "distance_fragmentation",
                                     "m_reach", "distance_reach"),
                 m = 1, mode = c("out", "in", "all"), weighted = NULL,
                 weights_as = NULL)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{k}{Number of key players to find.}
\item{set}{Identifiers of the nodes in the set to score (values of the
identifier column of \code{x$nodes}); may be empty.}
\item{measure}{The criterion; see Details. \code{"fragmentation"} and
\code{"distance_fragmentation"} solve KPP-Neg, \code{"m_reach"} and
\code{"distance_reach"} KPP-Pos.}
\item{m}{Distance limit for \code{"m_reach"}: the number of steps in a
binary network, or the largest path length when weights are used. Ignored
by the other criteria.}
\item{mode}{For directed networks, the direction of paths: \code{"out"}
follows the arcs (the set reaches others; \eqn{i} reaches \eqn{j}),
\code{"in"} reverses them, and \code{"all"} ignores direction. Ignored for
undirected networks.}
\item{weighted}{Whether to use edge weights as path lengths. \code{NULL} (the
default) uses them when the network has a weight column.}
\item{weights_as}{Override the network's \code{weight_type}: \code{"distance"}
or \code{"strength"}.}
\item{method}{\code{"greedy"} (Borgatti's swap heuristic with random
restarts) or \code{"exact"} (every set of size \code{k}; refused when there
are more than 50 million sets).}
\item{restarts}{Number of random starting sets for the greedy search.}
\item{seed}{Seed for the random starting sets, a whole number between 0 and
\eqn{2^{53}}; \code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}()} makes results reproducible.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{net_frame}}). Results do not depend on it.}
}
\details{
With \eqn{S} the set and \eqn{n} the number of nodes of the network being
measured, the criteria are (Borgatti 2006, equations 3, 9, 12, and 14):
\describe{
  \item{\code{"fragmentation"}}{\eqn{F = 1 - \sum_{i \ne j} r_{ij} / (n(n-1))}
    in the network without \eqn{S}, where \eqn{r_{ij} = 1} if \eqn{i}
    reaches \eqn{j}: the share of ordered pairs of remaining nodes that can
    no longer reach each other. For undirected networks this is
    \code{\link{net_fragmentation}()} of the remaining network.}
  \item{\code{"distance_fragmentation"}}{\eqn{DF = 1 - \sum_{i \ne j} (1/d_{ij}) / (n(n-1))}
    in the network without \eqn{S}, with \eqn{1/\infty = 0}: as \eqn{F},
    but pairs that stay connected count less the farther apart they are.}
  \item{\code{"m_reach"}}{The number of nodes outside \eqn{S} within
    distance \code{m} of at least one member. Borgatti defines it as a
    count, and it is reported as one.}
  \item{\code{"distance_reach"}}{\eqn{DR = \sum_j (1/d_{Sj}) / n}, where
    \eqn{d_{Sj}} is the distance from the nearest member of \eqn{S} to
    \eqn{j}, and members count as distance 1. It is 1 when every node is a
    member or adjacent to one.}
}
Larger values are better for every criterion. Borgatti writes \eqn{F} and
\eqn{DF} over the unordered pairs of an undirected network with a factor
of 2; the ordered-pair forms above are the same there and extend them to
directed paths.

\code{method = "greedy"} is Borgatti's (2006, figure 6) heuristic: start
from a random set, evaluate every swap of a member for a non-member, make
the swap that improves the criterion most, and repeat until no swap
improves it. Each of the \code{restarts} starting sets comes from its own
random stream of \code{seed}, the restarts run in parallel, and the best
result is returned, so the answer is the same for any number of threads.
The search can stop at a local optimum; more restarts make that less
likely, and \code{method = "exact"} settles it on small problems. Ties
between sets go to the set whose sorted node rows come first, and scores
within a relative \eqn{10^{-9}} of each other count as ties.

Conventions: distances are shortest-path lengths. Without weights they
are numbers of steps, found by breadth-first search. With weights, they
are distances, or strengths converted to lengths \eqn{1/w} (see
\code{weight_type} in \code{\link{net_frame}}); they must be positive.
Self-loops and parallel edges do not change distances, so they do not
matter. KPP-Neg criteria are \code{NaN} when fewer than two nodes remain.

Cost: each KPP-Neg evaluation measures every pair of the remaining
network, so one greedy step costs \eqn{O(k (n - k) n m)}; KPP-Pos
evaluations need a single search from the set.

Validation: every criterion is checked for every single-node set, and
the best value over all sets of size 2 or 3 is checked, against a
brute-force evaluation of the definitions. For single nodes on undirected
networks, \code{"distance_fragmentation"} equals the fragmentation
centrality of the keyplayer package (\code{keyplayer::fragment()} with
\code{binary = TRUE}; An and Liu 2016). keyplayer's group scores use a
different construction (contracting the group into one node), so they are
not comparable to these criteria.
}
\value{
\code{net_key_players()} returns a list with \code{set} (the identifiers of
the key players, in node order), \code{index} (their rows in
\code{x$nodes}), \code{score} (the criterion's value for the set),
\code{measure}, and \code{evaluations} (the number of sets evaluated).
\code{net_key_player_score()} returns a number.
}
\references{
Borgatti, S. P. (2006). Identifying sets of key players in a social network.
\emph{Computational and Mathematical Organization Theory}, 12(1), 21-34.
\doi{10.1007/s10588-006-7084-x}

Everton, S. F. (2012). \emph{Disrupting Dark Networks}. Cambridge
University Press. \doi{10.1017/CBO9781139136877}

An, W., & Liu, Y.-H. (2016). keyplayer: An R package for locating key
players in social networks. \emph{The R Journal}, 8(1), 257-268.
\doi{10.32614/RJ-2016-018}
}
\seealso{
\code{\link{net_fragmentation}}, \code{\link{net_group_centrality}},
\code{\link{net_flow_betweenness}}; \code{\link{net_example}} for the
networks in the examples.
}
\examples{
# The Jemaah Islamiyah cell behind the 2002 Bali bombings (Koschade 2006).
ji <- net_example("jemaah_islamiyah")

# Whom to remove to break the cell apart (KPP-Neg): the greedy search and
# exact enumeration both pick Muklas, the overall leader, and Samudra, the
# field commander.
net_key_players(ji, k = 2, measure = "fragmentation", seed = 1, threads = 2)
net_key_players(ji, k = 2, measure = "distance_fragmentation", method = "exact",
            threads = 2)

# Whom to watch to reach the whole cell quickly (KPP-Pos).
net_key_players(ji, k = 1, measure = "distance_reach", method = "exact", threads = 2)
net_key_players(ji, k = 2, measure = "m_reach", m = 1, method = "exact", threads = 2)

# Scoring a set chosen by hand: the command team.
command <- ji$nodes$name[ji$nodes$role == "command team"]
net_key_player_score(ji, command, measure = "fragmentation")
}
