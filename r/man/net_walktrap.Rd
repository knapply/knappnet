\name{net_walktrap}
\alias{net_walktrap}
\title{Walktrap communities}
\description{
Agglomerative communities from random-walk distances (Pons and Latapy 2006):
starting from single nodes, the two adjacent groups whose merger least
increases the mean squared distance between the walks that start in a node
and the walks that start in its group (Ward's criterion) are merged, and
the level of highest modularity is returned together with the whole
dendrogram.
}
\usage{
net_walktrap(x, steps = 4L, weighted = NULL, threads = 0L)
}
\arguments{
\item{x}{An undirected \code{\link{net_frame}}.}
\item{steps}{Length \eqn{t} of the random walks, a positive whole number.
Short walks make the distances local; the default follows the paper and
igraph.}
\item{weighted}{Whether to use edge weights, as given (larger means a
stronger tie). \code{NULL} (the default) uses them when the network has a
weight column. They must be positive.}
\item{threads}{Number of threads for the walk vectors, a whole number;
\code{0} uses all available cores, or 2 during \verb{R CMD check} (see
\sQuote{Threads} in \code{\link{net_frame}}). Results do not depend on it.}
}
\details{
Self-loops of the network add their weight to the diagonal of the
adjacency matrix, and parallel edges add their weights. As in the paper,
every node also gets a loop of its own, of weight equal to the mean weight
of its other incident edges (1 without weights); a node without edges gets
a loop of weight 1, as igraph documents. Directed networks are an error
(igraph silently ignores directions); see \code{\link{net_symmetrize}}.

Every distance is evaluated from the group vectors, so each merge is the
exact minimum up to rounding. Only adjacent groups merge, so a network with
\eqn{c} components has \eqn{n - c} merges and its dendrogram is
incomplete (igraph completes it by joining the components, the paper does
not). Among merges with equal Ward distance (to a relative
\eqn{10^{-9}}) the pair with the smaller group names is made first; igraph
documents no tie rule. Exact ties are common, for example between
structurally equivalent nodes such as the leaves of a star, and on such
networks (Zachary's karate club among them) the merge order, and with it
the partition, can differ from igraph's \code{cluster_walktrap()}. The
validation suite compares the engine with an exact-arithmetic
implementation of the paper everywhere, and with igraph on the networks
without ties, where merges and modularity agree. The walk vectors are
dense: memory grows as \eqn{8 n^2} bytes, so very large networks are out
of reach.

The \code{"merges"} attribute uses igraph's numbering: the \eqn{n} nodes are
clusters 1 to \eqn{n} and row \eqn{i} joins two clusters into cluster
\eqn{n + i}. It carries the number of nodes in its attribute
\code{"leaves"}. Use \code{\link{net_cut_dendrogram}} to cut the tree at
another level and \code{\link{net_to_hclust}} to plot it.
}
\value{An integer vector of group labels aligned with the rows of
\code{x$nodes}, numbered 1, 2, ... in order of first appearance, with the
attributes
\describe{
  \item{\code{"modularity"}}{Modularity of the partition (see
    \code{\link{net_modularity}}).}
  \item{\code{"merges"}}{A two-column integer matrix with the merges, in
    igraph's numbering.}
  \item{\code{"modularity_path"}}{Modularity before any merge and after
    each merge (one more value than there are merges).}
}
}
\references{
Pons, P., & Latapy, M. (2006). Computing communities in large networks
using random walks. \emph{Journal of Graph Algorithms and Applications},
10(2), 191--218. \doi{10.7155/jgaa.00124}
}
\seealso{\code{\link{net_cut_dendrogram}}, \code{\link{net_to_hclust}},
\code{\link{net_louvain}}, \code{\link{net_plot_partition}}}
\examples{
g <- net_example("karate")
groups <- net_walktrap(g, threads = 1)
groups
net_plot_partition(g, groups, seed = 1)

# Longer walks smooth the distances over more of the network.
attr(net_walktrap(g, steps = 8, threads = 1), "modularity")

# The dendrogram, cut into three groups.
net_cut_dendrogram(attr(groups, "merges"), k = 3)
}
