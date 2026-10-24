\name{net_cut_dendrogram}
\alias{net_cut_dendrogram}
\alias{net_to_hclust}
\title{Cut and convert a dendrogram of communities}
\description{
\code{net_cut_dendrogram()} turns the merges of a hierarchical community
method, such as the \code{"merges"} attribute of \code{\link{net_walktrap}},
into a partition, cutting it at a number of groups, after a number of
merges, or at the level of highest modularity. \code{net_to_hclust()}
converts the merges into an \code{\link[stats]{hclust}} object, to plot or
to use with the functions of \pkg{stats}.
}
\usage{
net_cut_dendrogram(merges, k = NULL, steps = NULL, x = NULL, weighted = NULL)

net_to_hclust(merges, x = NULL)
}
\arguments{
\item{merges}{A two-column matrix of whole numbers in igraph's numbering:
the \eqn{n} nodes are clusters 1 to \eqn{n} and row \eqn{i} joins two
existing clusters into cluster \eqn{n + i}. The number of nodes is read
from the attribute \code{"leaves"} that \code{net_walktrap()} sets, from
\code{x}, or else taken to be one more than the number of rows (a complete
dendrogram).}
\item{k}{Number of groups to cut into, a positive whole number.}
\item{steps}{Number of merges to perform, a non-negative whole number.}
\item{x}{A \code{\link{net_frame}} (undirected, for cutting) whose nodes
are the leaves. \code{net_cut_dendrogram()} needs it only to cut at the
level of highest modularity, when neither \code{k} nor \code{steps} is
given; \code{net_to_hclust()} uses its node identifiers as labels.}
\item{weighted}{Whether to use edge weights, as given, in the modularity
used to cut at its maximum. \code{NULL} (the default) uses them when the
network has a weight column.}
}
\details{
Give at most one of \code{k} and \code{steps}; with neither, \code{x} must
be given and the tree is cut at the first level of highest modularity (the
modularity of \code{\link{net_modularity}} with resolution 1). A
dendrogram of a disconnected network from a method that merges only
connected groups is incomplete, with \eqn{n - c} merges for \eqn{c}
components; cutting it into fewer than \eqn{c} groups is an error, as in
igraph's \code{cut_at()}. \code{steps} corresponds to igraph's
\code{community_to_membership()} and \code{k} to \code{cut_at(no = k)}.

\code{net_to_hclust()} builds the \code{hclust} tree by negating leaves and
shifting cluster numbers, as R's convention requires. A merge list has no
distances, so the heights are the merge order, 1, 2, ..., which draws the
tree with equally spaced levels and cuts it with \code{\link[stats]{cutree}}
exactly as \code{net_cut_dendrogram()} does; an incomplete dendrogram is
completed above its last merge by joining the remaining clusters one after
another.
}
\value{
\code{net_cut_dendrogram()} returns an integer vector of group labels
aligned with the leaves, numbered 1, 2, ... in order of first appearance.
When cut at the maximum modularity it has the attributes
\code{"modularity"} and \code{"modularity_path"} (the modularity of the
network after 0, 1, ..., all merges).

\code{net_to_hclust()} returns an object of class \code{"hclust"}.
}
\seealso{\code{\link{net_walktrap}}, \code{\link{net_hclust}},
\code{\link{net_cut_clusters}}}
\examples{
g <- net_example("karate")
tree <- attr(net_walktrap(g, threads = 1), "merges")

net_cut_dendrogram(tree, k = 3)
net_cut_dendrogram(tree, steps = 30)
best <- net_cut_dendrogram(tree, x = g)
attr(best, "modularity")

# As an hclust object, drawn and cut with the functions of stats.
hc <- net_to_hclust(tree, g)
plot(hc, cex = 0.6, main = "Walktrap dendrogram")
stats::cutree(hc, k = 3)
}
