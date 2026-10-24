\name{net_distances}
\alias{net_distances}
\alias{net_mean_distance}
\alias{net_diameter}
\alias{net_eccentricity}
\alias{net_radius}
\title{Shortest-path distances and their summaries}
\description{
\code{net_distances()} gives the length of a shortest path between every pair of
nodes. \code{net_mean_distance()} averages those lengths and \code{net_diameter()} is
the longest of them. A node's eccentricity is its largest distance to the
nodes it reaches, and the radius is the smallest eccentricity.
}
\usage{
net_distances(x, from = NULL, to = NULL, mode = c("out", "in", "all"),
          weighted = NULL, weights_as = NULL, threads = 0L)

net_mean_distance(x, mode = c("out", "in", "all"), weighted = NULL,
              weights_as = NULL, unconnected = TRUE, threads = 0L)

net_diameter(x, mode = c("out", "in", "all"), weighted = NULL,
         weights_as = NULL, unconnected = TRUE, threads = 0L)

net_eccentricity(x, mode = c("all", "out", "in"), weighted = NULL,
             weights_as = NULL, threads = 0L)

net_radius(x, mode = c("all", "out", "in"), weighted = NULL,
       weights_as = NULL, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{from, to}{Nodes to measure paths from and to, as values of the node
identifier column; \code{NULL} uses every node.}
\item{mode}{For directed networks: \code{"out"} follows edges in their
direction (row \eqn{i}, column \eqn{j} of \code{net_distances()} holds the
distance from \eqn{i} to \eqn{j}), \code{"in"} against it (the distance from
\eqn{j} to \eqn{i}), and \code{"all"} ignores direction. Ignored for
undirected networks. The defaults follow igraph: \code{"out"} for
\code{net_distances()}, \code{net_mean_distance()}, and \code{net_diameter()}, and
\code{"all"} for \code{net_eccentricity()} and \code{net_radius()}.}
\item{weighted}{Whether to use edge weights. \code{NULL} (the default) uses
them when the network has a weight column.}
\item{weights_as}{Override the network's \code{weight_type}: \code{"distance"}
or \code{"strength"}.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{net_frame}}).}
\item{unconnected}{What to do when some node does not reach every other:
\code{TRUE} averages (or maximizes) over the pairs joined by a path, and
\code{FALSE} returns \code{Inf}.}
}
\details{
A node is at distance 0 from itself and at distance \code{Inf} from nodes it
cannot reach. Self-loops never shorten a path and parallel edges count
through the shortest of them. Weights are handled as in
\code{\link{net_betweenness}}: distances as given, or strengths converted to
lengths 1/weight, and they must be positive; when every weight is the same,
lengths are that weight times the number of steps. The fastest exact
algorithm is chosen automatically (breadth-first search, Dial's bucket queue
with exact integer arithmetic, or Dijkstra's algorithm), with one search per
source node, run in parallel. \code{is.finite(net_distances(x))} is the
reachability matrix.

\code{net_mean_distance()} is the characteristic path length of Watts and
Strogatz (1998), averaged over ordered pairs of distinct nodes. As in
igraph's \code{net_mean_distance()} and \code{net_diameter()},
\code{unconnected = TRUE} uses only the pairs joined by a path, and
\code{unconnected = FALSE} gives \code{Inf} when some node does not reach
every other (for directed networks with \code{mode = "out"}, when the
network is not strongly connected). The mean is \code{NaN} when no pair is
joined; the diameter is 0 for a single node and \code{NaN} for a network
without nodes. NetworkX's \code{average_shortest_path_length()} and
\code{net_diameter()} agree on connected networks and refuse others.

Eccentricity ignores the nodes a node cannot reach, so an isolated node has
eccentricity 0, as in igraph's \code{net_eccentricity()} and \code{net_radius()}
(see Hage and Harary 1995 for these measures). NetworkX's
\code{net_eccentricity()} agrees on connected networks and refuses others.

Every function here matches igraph, and \code{net_distances()} also NetworkX and
sna's \code{geodist()}, on every network in the validation suite.
}
\value{
\code{net_distances()} returns a numeric matrix with one row per node in
\code{from} and one column per node in \code{to}, with node identifiers as
dimnames. \code{net_mean_distance()}, \code{net_diameter()}, and \code{net_radius()}
return a number; \code{net_eccentricity()} returns a numeric vector aligned with
the rows of \code{x$nodes}.
}
\references{
Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
\emph{Numerische Mathematik}, 1(1), 269-271. \doi{10.1007/BF01386390}

Dial, R. B. (1969). Algorithm 360: Shortest-path forest with topological
ordering. \emph{Communications of the ACM}, 12(11), 632-633.
\doi{10.1145/363269.363610}

Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of 'small-world'
networks. \emph{Nature}, 393(6684), 440-442. \doi{10.1038/30918}

Hage, P., & Harary, F. (1995). Eccentricity and centrality in networks.
\emph{Social Networks}, 17(1), 57-63. \doi{10.1016/0378-8733(94)00248-9}
}
\examples{
g <- net_frame(data.frame(source = c("a", "b", "c"), target = c("b", "c", "d"),
                         minutes = c(5, 5, 20)), weight = "minutes", weight_type = "distance")
net_distances(g)
net_distances(g, weighted = FALSE, mode = "all")
net_distances(g, from = "a", to = c("c", "d"))
net_mean_distance(g)
net_diameter(g, weighted = FALSE)
net_eccentricity(g, weighted = FALSE)
net_radius(g, weighted = FALSE)
}
