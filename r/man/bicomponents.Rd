\name{bicomponents}
\alias{bicomponents}
\title{Bicomponents (blocks)}
\description{
A bicomponent is a maximal set of nodes that stays connected when any one
of its nodes is removed; two bicomponents share at most one node, a cut
point (Wasserman and Faust 1994, section 4.2). Found with the depth-first
search of Hopcroft and Tarjan (1973) in time linear in the size of the
network.
}
\usage{
bicomponents(x)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
}
\details{
The undirected simple graph underlying the network is used: an edge in
either direction joins two nodes, parallel edges count once, self-loops
and weights are ignored. A bridge and its two endpoints form a
bicomponent of two nodes, so every edge lies in exactly one bicomponent;
isolated nodes belong to none. This matches igraph's
\code{biconnected_components()} and NetworkX's
\code{biconnected_components()}, and UCINET's Network > Regions >
Bi-Component routine. See \code{\link{cut_points}} for the cut points
themselves.
}
\value{
A list of vectors of node identifiers (values of the identifier column of
\code{x$nodes}), one per bicomponent, each in node order; the bicomponents
are in lexicographic order of node positions.
}
\references{
Hopcroft, J., & Tarjan, R. (1973). Algorithm 447: Efficient algorithms for
graph manipulation. \emph{Communications of the ACM}, 16(6), 372-378.
\doi{10.1145/362248.362272}

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods
and Applications}. Cambridge University Press.
\doi{10.1017/CBO9780511815478}
}
\seealso{\code{\link{cut_points}}, \code{\link{components}}}
\examples{
# Two triangles joined by a bridge 3 - 4, and a pendant 6 - 7.
g <- netframe(data.frame(from = c(1, 1, 2, 3, 4, 4, 5, 6),
                         to = c(2, 3, 3, 4, 5, 6, 6, 7)), directed = FALSE)
bicomponents(g)
cut_points(g)
}
