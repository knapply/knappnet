\name{distances}
\alias{distances}
\alias{mean_distance}
\alias{diameter}
\alias{eccentricity}
\alias{radius}
\title{Shortest-path distances, mean distance, diameter, eccentricity, and radius}
\description{
\code{distances()} gives the length of a shortest path between every pair of
nodes. \code{mean_distance()} averages those lengths and \code{diameter()} is
the longest of them. A node's eccentricity is its largest distance to the
nodes it reaches, and the radius is the smallest eccentricity.
}
\usage{
distances(x, from = NULL, to = NULL, mode = c("out", "in", "all"),
          weighted = NULL, weights_as = NULL, threads = 0L)

mean_distance(x, mode = c("out", "in", "all"), weighted = NULL,
              weights_as = NULL, unconnected = TRUE, threads = 0L)

diameter(x, mode = c("out", "in", "all"), weighted = NULL,
         weights_as = NULL, unconnected = TRUE, threads = 0L)

eccentricity(x, mode = c("all", "out", "in"), weighted = NULL,
             weights_as = NULL, threads = 0L)

radius(x, mode = c("all", "out", "in"), weighted = NULL,
       weights_as = NULL, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{from, to}{Nodes to measure paths from and to, as values of the node
identifier column; \code{NULL} uses every node.}
\item{mode}{For directed networks: \code{"out"} follows edges in their
direction (row \eqn{i}, column \eqn{j} of \code{distances()} holds the
distance from \eqn{i} to \eqn{j}), \code{"in"} against it (the distance from
\eqn{j} to \eqn{i}), and \code{"all"} ignores direction. Ignored for
undirected networks. The defaults follow igraph: \code{"out"} for
\code{distances()}, \code{mean_distance()}, and \code{diameter()}, and
\code{"all"} for \code{eccentricity()} and \code{radius()}.}
\item{weighted}{Whether to use edge weights. \code{NULL} (the default) uses
them when the network has a weight column.}
\item{weights_as}{Override the network's \code{weight_type}:
\code{"distance"} or \code{"strength"}.}
\item{threads}{Number of threads; \code{0} uses all available cores.}
\item{unconnected}{What to do when some node does not reach every other:
\code{TRUE} averages (or maximizes) over the pairs joined by a path, and
\code{FALSE} returns \code{Inf}.}
}
\details{
A node is at distance 0 from itself and at distance \code{Inf} from nodes it
cannot reach. Self-loops never shorten a path and parallel edges count
through the shortest of them. Weights are handled as in
\code{\link{betweenness}}: distances as given, or strengths converted to
lengths 1/weight, and they must be positive; when every weight is the same,
lengths are that weight times the number of steps. The fastest exact
algorithm is chosen automatically (breadth-first search, Dial's bucket queue
with exact integer arithmetic, or Dijkstra's algorithm), with one search per
source node, run in parallel. \code{is.finite(distances(x))} is the
reachability matrix.

\code{mean_distance()} is the characteristic path length of Watts and
Strogatz (1998), averaged over ordered pairs of distinct nodes. As in
igraph's \code{mean_distance()} and \code{diameter()},
\code{unconnected = TRUE} uses only the pairs joined by a path, and
\code{unconnected = FALSE} gives \code{Inf} when some node does not reach
every other (for directed networks with \code{mode = "out"}, when the
network is not strongly connected). The mean is \code{NaN} when no pair is
joined; the diameter is 0 for a single node and \code{NaN} for a network
without nodes. NetworkX's \code{average_shortest_path_length()} and
\code{diameter()} agree on connected networks and refuse others.

Eccentricity ignores the nodes a node cannot reach, so an isolated node has
eccentricity 0, as in igraph's \code{eccentricity()} and \code{radius()}
(see Hage and Harary 1995 for these measures). NetworkX's
\code{eccentricity()} agrees on connected networks and refuses others.

Every function here matches igraph, and \code{distances()} also NetworkX and
sna's \code{geodist()}, on every network in the validation suite.
}
\value{
\code{distances()} returns a numeric matrix with one row per node in
\code{from} and one column per node in \code{to}, with node identifiers as
dimnames. \code{mean_distance()}, \code{diameter()}, and \code{radius()}
return a number; \code{eccentricity()} returns a numeric vector aligned with
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
g <- netframe(data.frame(from = c("a", "b", "c"), to = c("b", "c", "d"),
                         minutes = c(5, 5, 20)), weight = "minutes")
distances(g)
distances(g, weighted = FALSE, mode = "all")
distances(g, from = "a", to = c("c", "d"))
mean_distance(g)
diameter(g, weighted = FALSE)
eccentricity(g, weighted = FALSE)
radius(g, weighted = FALSE)
}
