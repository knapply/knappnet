\name{cluster_louvain}
\alias{cluster_louvain}
\title{Louvain community detection}
\description{
Finds groups of nodes by greedy optimization of modularity with the
method of Blondel, Guillaume, Lambiotte, and Lefebvre (2008).
}
\usage{
cluster_louvain(x, resolution = 1, weighted = NULL, seed = NULL)
}
\arguments{
\item{x}{An undirected \code{\link{netframe}}.}
\item{resolution}{Resolution parameter \eqn{\gamma} (non-negative); see
\code{\link{modularity}}.}
\item{weighted}{Whether to use edge weights. \code{NULL} (the default) uses
them when the network has a weight column.}
\item{seed}{Seed for the random node order: a whole number between 0 and
\eqn{2^{53}}. \code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}()} makes results reproducible.}
}
\details{
Each level repeatedly moves single nodes to the neighboring group that
most increases modularity (staying put on ties) until no move improves it,
then merges each group into one node and repeats on the merged network,
until a level makes no move. Within a level, nodes are processed from a
queue: after a node moves, only its neighbors outside its new group are
revisited, as in the fast local move procedure of Traag, Waltman, and van
Eck (2019). This finds the same kind of local optimum as repeated sweeps
over every node, with far less work on large networks. Nodes are first
queued in a random order drawn
from a counter-based generator (Philox; Salmon et al. 2011), so the same
\code{seed} gives the same result on every platform.

Only undirected networks are supported, as in igraph; a directed network
is an error rather than being silently symmetrized. Weights, self-loops,
parallel edges, and the resolution parameter are treated as in
\code{\link{modularity}}, and the reported modularity is
\code{modularity(x, labels, resolution)} of the returned partition. A
network without edges puts every node in its own group, with modularity
\code{NaN}.

Like every implementation of this heuristic, results depend on the node
order and are local optima. Partitions generally differ from those of
\code{igraph::cluster_louvain()} and NetworkX's
\code{louvain_communities()}; the validation suite checks that their
modularity is comparable.
}
\value{An integer vector of group labels aligned with the rows of
\code{x$nodes}, numbered 1, 2, ... in order of first appearance, with the
modularity of the partition in attribute \code{"modularity"}.}
\references{
Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
Fast unfolding of communities in large networks. \emph{Journal of
Statistical Mechanics: Theory and Experiment}, 2008(10), P10008.
\doi{10.1088/1742-5468/2008/10/P10008}

Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to
Leiden: Guaranteeing well-connected communities. \emph{Scientific
Reports}, 9, 5233. \doi{10.1038/s41598-019-41695-z}

Newman, M. E. J., & Girvan, M. (2004). Finding and evaluating community
structure in networks. \emph{Physical Review E}, 69(2), 026113.
\doi{10.1103/PhysRevE.69.026113}

Reichardt, J., & Bornholdt, S. (2006). Statistical mechanics of community
detection. \emph{Physical Review E}, 74(1), 016110.
\doi{10.1103/PhysRevE.74.016110}

Salmon, J. K., Moraes, M. A., Dror, R. O., & Shaw, D. E. (2011). Parallel
random numbers: As easy as 1, 2, 3. In \emph{Proceedings of the 2011
International Conference for High Performance Computing, Networking,
Storage and Analysis (SC '11)}. \doi{10.1145/2063384.2063405}
}
\seealso{\code{\link{modularity}}}
\examples{
g <- netframe(data.frame(from = c(1, 1, 2, 3, 4, 4, 5), to = c(2, 3, 3, 4, 5, 6, 6)),
              directed = FALSE)
groups <- cluster_louvain(g, seed = 1)
groups
attr(groups, "modularity")
}
