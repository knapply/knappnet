\name{net_girvan_newman}
\alias{net_girvan_newman}
\alias{net_leading_eigen}
\alias{net_fast_greedy}
\title{Divisive, spectral, and greedy community detection}
\description{
Three classic partitioning methods for undirected networks, each returning
the partition with its modularity (see \code{\link{net_modularity}}):
\describe{
  \item{\code{net_girvan_newman()}}{Repeatedly removes the edge with the
    highest edge betweenness, recomputing betweenness after every removal,
    so that the network falls apart into components (Girvan and Newman
    2002); of the partitions met along the way, the one with the highest
    modularity is returned (Newman and Girvan 2004). Edge betweenness is
    accumulated with Brandes' (2001) algorithm, recomputed only in the
    component that lost an edge.}
  \item{\code{net_leading_eigen()}}{Divides the network in two by the
    signs of the leading eigenvector of the modularity matrix, and each
    group again by its generalized modularity matrix, until no division
    raises modularity (Newman 2006).}
  \item{\code{net_fast_greedy()}}{Starts with every node in its own
    group and repeatedly merges the two connected groups whose merger raises
    modularity the most, returning the best partition along the way
    (Clauset, Newman, and Moore 2004; Newman 2004).}
}
}
\usage{
net_girvan_newman(x, threads = 0L)

net_leading_eigen(x, weighted = NULL)

net_fast_greedy(x, resolution = 1, weighted = NULL)
}
\arguments{
\item{x}{An undirected \code{\link{net_frame}}.}
\item{threads}{Number of threads for the betweenness computations, a whole
number; \code{0} uses all available cores, or 2 during \verb{R CMD check} (see
\sQuote{Threads} in \code{\link{net_frame}}).}
\item{weighted}{Whether to use edge weights, as given (larger means a stronger
tie). \code{NULL} (the default) uses them when the network has a weight
column. They must not be negative.}
\item{resolution}{Resolution parameter \eqn{\gamma} (non-negative); see
\code{\link{net_modularity}}.}
}
\details{
All three accept only undirected networks, as the igraph functions they
correspond to (\code{cluster_edge_betweenness()},
\code{net_leading_eigen()}, and \code{net_fast_greedy()}) do; a
directed network is an error rather than being silently symmetrized (see
\code{\link{net_symmetrize}}). Groups are numbered in order of their first
node. A network without edges puts every node in its own group with
modularity \code{NaN}.

\strong{Girvan-Newman} reads the network as a simple graph (parallel edges
count once, self-loops are ignored) and counts hops, without weights; the
reported modularity is that of the simple graph, as NetworkX's
\code{girvan_newman()} with \code{net_modularity()} on a \code{Graph} gives.
The edge removed at each step is the one with the largest betweenness;
values within a relative \eqn{10^{-9}} of the largest count as tied, and
the tie goes to the edge whose endpoints come first in node order. igraph
removes the tied edge with the lowest edge id and NetworkX the first in its
iteration order, so on networks with tied betweenness the three can split
differently. The chosen partition is the first of highest modularity in
removal order. Cost is \eqn{O(m^2 n)} in the worst case; searches from each
source run in parallel.

\strong{Leading eigenvector} and \strong{greedy merging} treat weights (when
used) as in \code{\link{net_modularity}}: they are added as given and must be
non-negative, parallel edges add up, and a self-loop of weight \eqn{w} adds
\eqn{2w} to its node's diagonal entry and degree. The leading eigenvector
of each group's matrix comes from a dense symmetric eigendecomposition
(Householder tridiagonalization and implicit QL; Martin, Reinsch, and
Wilkinson 1968; Bowdler et al. 1968), which costs \eqn{O(n^3)} time and
\eqn{O(n^2)} memory for the first division, so networks above 3000 nodes
are refused. A group is left undivided when the leading eigenvalue is not
positive (larger than \eqn{10^{-10}} times a bound on the eigenvalues) or
when the division does not raise modularity. The eigenvector is scaled so
that its largest component is positive; nodes with positive components form
one side and the others the other. Results match igraph's
\code{net_leading_eigen()} except where an eigenvalue is degenerate (the
network has a symmetry) or a component is numerically zero, where the
division is not unique.

Greedy merging compares the change in modularity exactly when weights are
whole numbers and \code{resolution = 1}, and breaks ties toward the groups
with the lowest-numbered members. igraph and NetworkX break ties by their
own internal order, so on networks with tied mergers the partitions can
differ while having similar modularity; the validation suite verifies that
any such difference is a pure tie-break. Mergers continue past the first
that lowers modularity until no two groups are connected, as in igraph, and
the first partition of highest modularity is returned. \code{resolution} is
the parameter of Reichardt and Bornholdt (2006), as in NetworkX's
\code{greedy_modularity_communities()}.
}
\value{
An integer vector of group labels aligned with the rows of \code{x$nodes},
numbered 1, 2, ... in order of first appearance, with the partition's
modularity in attribute \code{"modularity"}.
\code{net_girvan_newman()} also records the edges removed, in order, as
a two-column data frame of node identifiers in attribute \code{"removed"},
the modularity before any removal and after each removal in attribute
\code{"path"}, and the divisive hierarchy as an \code{hclust} object in
attribute \code{"tree"}, whose heights count removals from the end (the
last edge removed joins two nodes at height 1, components never connected
join at height \eqn{m + 1}), so that \code{plot()} of it draws the
dendrogram. \code{net_leading_eigen()} records the number of divisions
in attribute \code{"splits"}.
}
\references{
Girvan, M., & Newman, M. E. J. (2002). Community structure in social and
biological networks. \emph{Proceedings of the National Academy of
Sciences}, 99(12), 7821-7826. \doi{10.1073/pnas.122653799}

Newman, M. E. J., & Girvan, M. (2004). Finding and evaluating community
structure in networks. \emph{Physical Review E}, 69(2), 026113.
\doi{10.1103/PhysRevE.69.026113}

Brandes, U. (2001). A faster algorithm for betweenness centrality.
\emph{Journal of Mathematical Sociology}, 25(2), 163-177.
\doi{10.1080/0022250X.2001.9990249}

Newman, M. E. J. (2006). Modularity and community structure in networks.
\emph{Proceedings of the National Academy of Sciences}, 103(23), 8577-8582.
\doi{10.1073/pnas.0601602103}

Clauset, A., Newman, M. E. J., & Moore, C. (2004). Finding community
structure in very large networks. \emph{Physical Review E}, 70(6), 066111.
\doi{10.1103/PhysRevE.70.066111}

Newman, M. E. J. (2004). Fast algorithm for detecting community structure
in networks. \emph{Physical Review E}, 69(6), 066133.
\doi{10.1103/PhysRevE.69.066133}

Reichardt, J., & Bornholdt, S. (2006). Statistical mechanics of community
detection. \emph{Physical Review E}, 74(1), 016110.
\doi{10.1103/PhysRevE.74.016110}

Martin, R. S., Reinsch, C., & Wilkinson, J. H. (1968). Householder's
tridiagonalization of a symmetric matrix. \emph{Numerische Mathematik},
11(3), 181-195. \doi{10.1007/BF02161841}

Bowdler, H., Martin, R. S., Reinsch, C., & Wilkinson, J. H. (1968). The QR
and QL algorithms for symmetric matrices. \emph{Numerische Mathematik},
11(4), 293-306. \doi{10.1007/BF02166681}
}
\seealso{\code{\link{net_modularity}}, \code{\link{net_louvain}},
\code{\link{net_factions}}}
\examples{
g <- net_example("karate")
gn <- net_girvan_newman(g, threads = 2)
table(gn)
attr(gn, "modularity")
# Modularity after each removal; the maximum is the returned partition.
plot(attr(gn, "path"), type = "l", xlab = "edges removed", ylab = "modularity")
# The divisive hierarchy as a dendrogram, and the partition on the network.
plot(attr(gn, "tree"), main = "Girvan-Newman")
net_plot_partition(g, gn, seed = 1, main = "Girvan-Newman partition")

le <- net_leading_eigen(g)
attr(le, "modularity")
fg <- net_fast_greedy(g)
attr(fg, "modularity")
# A finer partition at higher resolution.
table(net_fast_greedy(g, resolution = 2))
}
