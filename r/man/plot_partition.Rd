\name{plot_partition}
\alias{plot_partition}
\title{Plot a partition of the network}
\description{
Draws the network with \code{\link{plot.netframe}}, coloring each node by
its group in a partition: communities from \code{\link{cluster_louvain}},
positions from \code{\link{concor}} or \code{\link{regular_equivalence}},
components, or any node attribute.
}
\usage{
plot_partition(x, membership, title = "group", ...)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{membership}{The group of each node: a vector aligned with the rows of
\code{x$nodes} (only equality of values matters), or the name of a column
of \code{x$nodes}.}
\item{title}{The name of the partition, used for the legend.}
\item{...}{Further arguments passed to \code{\link{plot.netframe}}, such as
\code{layout}, \code{seed}, \code{node_size}, or \code{main}.}
}
\details{
Groups are colored with the qualitative palette of
\code{\link{plot.netframe}} in order of their first appearance in
\code{x$nodes}, so group numbers from \code{\link{cluster_louvain}}
(numbered in order of first appearance) keep their order. More than 15
groups are listed in the legend as a count. When \code{membership} carries
a modularity (as \code{cluster_louvain()}'s result does) and no
\code{main} is given, the title reports the number of groups and the
modularity.
}
\value{The layout used, invisibly, as \code{\link{plot.netframe}} returns.}
\seealso{\code{\link{cluster_louvain}}, \code{\link{plot_matrix}} for the
same partition as a sorted adjacency matrix, \code{\link{plot.blockmodel}}.}
\examples{
g <- example_network("karate")
groups <- cluster_louvain(g, seed = 1)
plot_partition(g, groups, title = "community", seed = 1)
plot_partition(g, "faction", seed = 1, main = "Factions")
}
