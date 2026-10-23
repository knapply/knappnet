\name{n_cliques}
\alias{n_cliques}
\alias{n_clans}
\alias{k_plexes}
\title{n-cliques, n-clans, and k-plexes}
\description{
Relaxations of the clique (see \code{\link{max_cliques}}) that admit groups
whose members are not all directly tied:
\describe{
  \item{n-clique}{A maximal set of nodes every two of which are joined by a
    path of length at most \code{n} in the whole network (Luce 1950; Alba
    1973). The paths may run through nodes outside the set. 1-cliques are
    cliques.}
  \item{n-clan}{An n-clique whose own subgraph has diameter at most
    \code{n}, so that the short paths stay inside the group (Mokken 1979).}
  \item{k-plex}{A maximal set of nodes each of which is adjacent to all but
    at most \code{k - 1} of the other members (Seidman and Foster 1978).
    1-plexes are cliques.}
}
}
\usage{
n_cliques(x, n = 2L, min_size = 1L, threads = 0L)

n_clans(x, n = 2L, min_size = 1L, threads = 0L)

k_plexes(x, k = 2L, min_size = NULL, max_steps = 1e7, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{n}{Largest path length between members: 1 for cliques, 2 for
2-cliques, and so on.}
\item{k}{Each member may have at most \code{k - 1} non-neighbors among the
other members (\code{k >= 1}).}
\item{min_size}{Smallest group to report. For k-plexes \code{NULL} means
\code{2k - 1}.}
\item{max_steps}{Limit on the number of search steps for k-plexes.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{netframe}}).}
}
\details{
All three use the undirected simple graph underlying the network:
directions, parallel edges, self-loops, and weights are ignored, as in
UCINET's Network > Subgroups > N-Cliques, N-Clans, and K-Plex routines for
symmetric data. Members are listed in node order and groups in
lexicographic order of member positions.

n-cliques are the maximal cliques of the graph joining two nodes when their
distance is at most \code{n}, which is built with one breadth-first search
per node (in parallel) and searched with the Bron-Kerbosch algorithm (Bron
and Kerbosch 1973); n-clans are the n-cliques that pass the diameter test.
With \code{min_size} at most 1 an isolated node is an n-clique and n-clan of
its own. \code{min_size} only drops smaller groups; it never reports
non-maximal subsets.

k-plexes are listed with a Bron-Kerbosch set enumeration (every subset of a
k-plex is a k-plex, which makes the scheme exact), after dropping the nodes
with fewer than \code{min_size - k} neighbors, which cannot belong to a
k-plex of \code{min_size} nodes (Seidman and Foster 1978). Any set of at
most \code{k} nodes is a k-plex, so small values of \code{min_size} report
many trivial sets; \code{min_size} below \code{2k - 1} also admits
disconnected k-plexes. The default, \code{2k - 1}, is the smallest size at
which every k-plex is connected with diameter at most 2. The number of
maximal k-plexes can grow exponentially with the network, so the search
stops with an error after \code{max_steps} steps rather than running
indefinitely; raise the limit, \code{k}, or \code{min_size} as needed.
Whether the limit is reached does not depend on \code{threads}.

The validation suite checks every group against its definition by
exhaustive enumeration, and n-cliques against NetworkX's
\code{find_cliques()} on the distance-\code{n} power graph.
}
\value{
A list of vectors of node identifiers (values of the identifier column of
\code{x$nodes}), one per group, each in node order; the groups are in
lexicographic order of node positions.
}
\references{
Luce, R. D. (1950). Connectivity and generalized cliques in sociometric
group structure. \emph{Psychometrika}, 15(2), 169-190.
\doi{10.1007/BF02289199}

Alba, R. D. (1973). A graph-theoretic definition of a sociometric clique.
\emph{Journal of Mathematical Sociology}, 3(1), 113-126.
\doi{10.1080/0022250X.1973.9989826}

Mokken, R. J. (1979). Cliques, clubs and clans. \emph{Quality and
Quantity}, 13(2), 161-173. \doi{10.1007/BF00139635}

Seidman, S. B., & Foster, B. L. (1978). A graph-theoretic generalization of
the clique concept. \emph{Journal of Mathematical Sociology}, 6(1),
139-154. \doi{10.1080/0022250X.1978.9989883}

Bron, C., & Kerbosch, J. (1973). Algorithm 457: Finding all cliques of an
undirected graph. \emph{Communications of the ACM}, 16(9), 575-577.
\doi{10.1145/362342.362367}
}
\seealso{\code{\link{max_cliques}}, \code{\link{lambda_sets}},
\code{\link{factions}}}
\examples{
# A hexagon 1-2-3-4-5-6-1: no triangle, so no clique beyond pairs.
g <- netframe(data.frame(from = 1:6, to = c(2:6, 1)), directed = FALSE)
max_cliques(g)
# Nodes two steps apart are within distance 2 through an outside node, so
# the alternating triples {1, 3, 5} and {2, 4, 6} are 2-cliques; they have
# no ties among themselves, so they are not 2-clans.
n_cliques(g, n = 2, threads = 2)
n_clans(g, n = 2, threads = 2)
# In a 2-plex every member may miss one other member: three consecutive
# nodes qualify.
k_plexes(g, k = 2, threads = 2)

jemaah <- example_network("jemaah_islamiyah")
length(n_cliques(jemaah, n = 2, min_size = 4, threads = 2))
length(k_plexes(jemaah, k = 2, min_size = 5, threads = 2))
}
