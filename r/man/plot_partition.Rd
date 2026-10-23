\name{plot_partition}
\alias{plot_partition}
\title{Draw a partition on a network}
\usage{
plot_partition(x, membership, layout = NULL, ..., legend = TRUE, seed = NULL)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{membership}{Group of each node, aligned with the rows of
\code{x$nodes}: any atomic vector, as \code{\link{modularity}} accepts.}
\item{layout}{Node positions, as in \code{\link{plot.netframe}}.}
\item{...}{Further arguments passed to \code{\link{plot.netframe}} and on
to \code{\link[graphics]{title}}, such as \code{main}, \code{node_size}, or
\code{node_label}.}
\item{legend}{Whether to draw a legend of the groups.}
\item{seed}{Seed for the layout, as in \code{\link{layout_network}}.}
}
\value{\code{x}, invisibly.}
\description{
Draws the network with nodes coloured by group, for the results of
\code{\link{factions}}, \code{\link{cluster_girvan_newman}},
\code{\link{cluster_leading_eigen}}, \code{\link{cluster_fast_greedy}},
\code{\link{cluster_louvain}}, \code{\link{cut_clusters}}, or any other
partition aligned with the rows of \code{x$nodes}. Colours come from the
Okabe-Ito palette (\code{\link[grDevices]{palette.colors}}), which is safe
for colour-blind readers, for up to eight groups, and from the HCL "Dark 3"
palette beyond that. Groups are listed in sorted order for numeric and
logical memberships and in order of first appearance otherwise.
}
\seealso{\code{\link{plot.netframe}}}
\examples{
g <- example_network("karate")
plot_partition(g, factions(g, k = 2, seed = 1, threads = 2), seed = 1,
               main = "Two factions")
plot_partition(g, g$nodes$faction, seed = 1, main = "The actual split")
}
