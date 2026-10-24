\name{net_leiden}
\alias{net_leiden}
\title{Leiden communities}
\description{
Finds groups of nodes with the Leiden algorithm (Traag, Waltman, and van
Eck 2019), which improves on the Louvain method by refining each group
before aggregating, so that every group it returns is connected, and
converges to partitions that are locally optimal at several scales. The
quality function is either modularity or the Constant Potts Model.
}
\usage{
net_leiden(x, objective = c("modularity", "cpm"), resolution = 1, beta = 0.01,
           n_iterations = -1L, weighted = NULL, node_weights = NULL,
           initial = NULL, seed = NULL)
}
\arguments{
\item{x}{An undirected \code{\link{net_frame}}.}
\item{objective}{Quality function: \code{"modularity"} or \code{"cpm"} (the
Constant Potts Model); see \sQuote{Objectives}.}
\item{resolution}{Resolution parameter \eqn{\gamma} (non-negative); larger
values give smaller groups.}
\item{beta}{Randomness of the refinement step, a positive number; the
default follows the paper.}
\item{n_iterations}{Number of iterations, a whole number; a negative number
(the default) iterates until an iteration changes nothing (at most 10,000).}
\item{weighted}{Whether to use edge weights, as given. \code{NULL} (the
default) uses them when the network has a weight column. They must not be
negative.}
\item{node_weights}{Node weights \eqn{n_i}: a numeric vector with one
value per node (or the name of such a column of \code{x$nodes}), finite
and non-negative. \code{NULL} uses the weighted degree for modularity and
1 for the CPM.}
\item{initial}{A starting partition: a vector with one label per node (any
atomic type, without missing values) or the name of a column of
\code{x$nodes}. \code{NULL} starts with every node alone.}
\item{seed}{Seed for the random node orders and the refinement, a whole
number between 0 and \eqn{2^{53}}; \code{NULL} draws one from R's random
number generator, so \code{\link{set.seed}()} makes results reproducible.}
}
\details{
\strong{Objectives.} Following igraph's \code{cluster_leiden()}, the
method maximizes
\deqn{Q = \frac{1}{2m} \sum_{ij} (A_{ij} - \gamma n_i n_j) \delta(c_i, c_j),}{
Q = 1/(2m) sum_ij (A_ij - gamma n_i n_j) delta(c_i, c_j),}
with \eqn{m} the total edge weight and \eqn{c_i} the group of node
\eqn{i}. With \code{objective = "modularity"} the node weights are the
weighted degrees and \eqn{\gamma} is divided by \eqn{2m}, as igraph does,
so that \eqn{Q} is the modularity of \code{\link{net_modularity}} with that
resolution. With \code{objective = "cpm"} the node weights are 1 (or those
given) and \eqn{\gamma} is the density threshold of the Constant Potts
Model (Traag, Van Dooren, and Nesterov 2011): groups are kept together only
if their internal density exceeds \eqn{\gamma}. \strong{The default differs
from igraph}, whose \code{cluster_leiden()} uses
\code{objective_function = "CPM"}; modularity is what most users expect.

Self-loops add twice their weight to the node's degree and never affect a
move; parallel edges add their weights. A network without edges puts every
node in its own group with quality \code{NaN}. Moves must raise the quality
by more than zero, and ties between equally good groups go to the first
neighbor group in node order. Every group is connected whenever
\eqn{\gamma > 0}; the paper proves this holds after each iteration.

The same \code{seed} gives the same result on every platform, but the
random numbers differ from igraph's, so partitions differ from
\code{igraph::cluster_leiden()}; the validation suite compares the mean
quality over many seeds and recomputes the quality of each result from its
definition.
}
\value{An integer vector of group labels aligned with the rows of
\code{x$nodes}, numbered 1, 2, ... in order of first appearance, with the
attributes
\describe{
  \item{\code{"quality"}}{The value of \eqn{Q} above for the partition.}
  \item{\code{"iterations"}}{Number of iterations run.}
}
}
\references{
Traag, V. A., Waltman, L., & van Eck, N. J. (2019). From Louvain to
Leiden: Guaranteeing well-connected communities. \emph{Scientific
Reports}, 9, 5233. \doi{10.1038/s41598-019-41695-z}

Traag, V. A., Van Dooren, P., & Nesterov, Y. (2011). Narrow scope for
resolution-limit-free community detection. \emph{Physical Review E},
84(1), 016114. \doi{10.1103/PhysRevE.84.016114}
}
\seealso{\code{\link{net_louvain}}, \code{\link{net_modularity}},
\code{\link{net_plot_partition}}}
\examples{
g <- net_example("karate")
groups <- net_leiden(g, seed = 1)               # igraph objective_function = "modularity"
groups
net_plot_partition(g, groups, seed = 1)

# Constant Potts Model: igraph's default objective ("CPM"); its resolution
# is a density, so values far below 1 are needed.
net_leiden(g, "cpm", resolution = 0.1, seed = 1)

# A finer resolution for modularity.
max(net_leiden(g, resolution = 2, seed = 1))
}
