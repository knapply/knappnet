\name{net_cut_points}
\alias{net_cut_points}
\alias{net_bridges}
\title{Cut points and bridges}
\description{
A cut point (articulation point) is a node whose removal increases the
number of components; a bridge is an edge whose removal does.
}
\usage{
net_cut_points(x)

net_bridges(x)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
}
\details{
Both use the undirected reading of the network in which every edge stays an
edge: directions are ignored, so two parallel edges, or two arcs in opposite
directions, between the same nodes are two edges and neither is a bridge;
self-loops are ignored. Weights do not matter. Both come from one
depth-first search that computes low points (Hopcroft and Tarjan 1973),
written without recursion so that very deep networks cannot exhaust the
stack, in time linear in the size of the network.

Results match NetworkX's \code{articulation_points()} and \code{net_bridges()} on
an undirected multigraph, igraph's \code{articulation_points()}, and sna's
\code{cutpoints()} on the symmetrized network, on every network in the
validation suite. igraph's \code{net_bridges()} agrees on every validation
network whose undirected reading has no parallel edges; its documentation
does not say how it treats parallel edges.
}
\value{
\code{net_cut_points()} returns the identifiers (values of
\code{x$nodes[[id]]}) of the cut points, in the order of \code{x$nodes}.
\code{net_bridges()} returns the row numbers of the bridges in \code{x$edges},
in increasing order.
}
\references{
Hopcroft, J., & Tarjan, R. (1973). Algorithm 447: Efficient algorithms for
graph manipulation. \emph{Communications of the ACM}, 16(6), 372-378.
\doi{10.1145/362248.362272}
}
\examples{
# Two triangles joined through node "c".
g <- net_frame(data.frame(source = c("a", "b", "c", "c", "d", "e"),
                         target   = c("b", "c", "a", "d", "e", "c")),
              directed = FALSE)
net_cut_points(g)
net_bridges(g)
# A path: every edge is a bridge.
p <- net_frame(data.frame(source = c(1, 2, 3), target = c(2, 3, 4)), directed = FALSE)
net_cut_points(p)
p$edges[net_bridges(p), ]
}
