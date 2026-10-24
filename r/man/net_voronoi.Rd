\name{net_voronoi}
\alias{net_voronoi}
\title{Voronoi communities}
\description{
Partitions a network into the graph Voronoi cells of generator nodes of
high local density (Deritei et al. 2014), in the extension to weighted and
directed networks of Molnar et al. (2024) that picks the radius of highest
modularity. This is the method of igraph's \code{cluster_voronoi()}.
}
\usage{
net_voronoi(x, lengths = NULL, weighted = NULL, mode = c("out", "in", "all"),
            radius = NULL, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{net_frame}}, directed or undirected.}
\item{lengths}{Edge lengths: a non-negative numeric vector with one value
per row of \code{x$edges}, or the name of such a column of \code{x$edges}.
\code{NULL} gives every edge length 1. Lengths are independent of
\code{weighted}.}
\item{weighted}{Whether to use edge weights, as given, for the strengths
(in the local density) and for modularity. \code{NULL} (the default) uses
them when the network has a weight column. They must not be negative.}
\item{mode}{For directed networks, the direction of distances:
\code{"out"} measures them from a generator along the edges (the generator
reaches the node), \code{"in"} against them (the node reaches the
generator), and \code{"all"} ignores directions. Undirected networks ignore
it.}
\item{radius}{The radius \eqn{r} that decides which nodes become
generators, a non-negative number (\code{Inf} allowed). \code{NULL} (the
default) chooses the radius of highest modularity.}
\item{threads}{Number of threads for the shortest-path searches, a whole
number; \code{0} uses all available cores, or 2 during \verb{R CMD check}
(see \sQuote{Threads} in \code{\link{net_frame}}). Results do not depend
on it.}
}
\details{
The length of an edge is its given length divided by its edge clustering
coefficient (Radicchi et al. 2004), \eqn{C_{ij} = (z_{ij} + 1) / \min(k_i
- 1, k_j - 1)} with \eqn{z_{ij}} the number of common neighbors and
\eqn{k} the degree in the underlying simple undirected network, so edges in
triangles are shorter. When \eqn{\min(k_i - 1, k_j - 1) = 0} (an edge to a
node of degree one) the length is 0, the limit of the formula, and the node
always joins its neighbor's cell. The local relative density of node
\eqn{i} is \eqn{\rho_i = s_i m_i / (m_i + k_i^{out})}{rho_i = s_i m_i /
(m_i + kout_i)}, with \eqn{s_i} its strength, \eqn{m_i} the number of edges
among \eqn{i} and its neighbors and \eqn{k_i^{out}}{kout_i} the number of
edges leaving that set. For a given radius, the densest node not yet
excluded becomes a generator, every node within distance less than
\eqn{r} of it is excluded, and this repeats until every node is a generator
or excluded; each node joins the cell of its nearest generator.

Degrees, neighborhoods, and clustering coefficients ignore direction,
parallel edges, and self-loops; among parallel arcs the shortest is used,
and self-loops are left out of strengths. Ties: generators are taken in
decreasing density, ties going to the lower node index, and a node at equal
distance from several generators joins the one taken first (the paper
chooses at random). Distances (and modularities) within a relative
\eqn{10^{-9}} are taken as equal. A node without edges has density 0.

The generator set changes only where \eqn{r} crosses a shortest-path
distance, so when \code{radius} is \code{NULL} the candidates are every
distinct finite distance between two nodes (and infinity, which leaves one
generator per set of mutually reachable nodes); with more than 1,000
distinct distances, 1,000 evenly spaced by rank are tried. The candidate of
highest modularity is used, ties going to the smallest radius. igraph (and
the paper) instead search with Brent's method, a local optimizer, so where
they differ this choice has the higher modularity, and the validation suite
lists those networks rather than expecting a match.
}
\value{An integer vector of group labels aligned with the rows of
\code{x$nodes}, numbered 1, 2, ... in order of first appearance, with the
attributes
\describe{
  \item{\code{"generators"}}{Positions in \code{x$nodes} of the generator
    nodes, in the order they were taken.}
  \item{\code{"modularity"}}{Modularity of the partition (directed for
    directed networks; \code{NaN} without edges).}
  \item{\code{"radius"}}{The radius used.}
}
}
\references{
Deritei, D., Lazar, Z. I., Papp, I., Jarai-Szabo, F., Sumi, R., Varga, L.,
Ravasz Regan, E., & Ercsey-Ravasz, M. (2014). Community detection by graph
Voronoi diagrams. \emph{New Journal of Physics}, 16(6), 063007.
\doi{10.1088/1367-2630/16/6/063007}

Molnar, B., Marton, I.-B., Horvat, S., & Ercsey-Ravasz, M. (2024).
Community detection in directed weighted networks using Voronoi
partitioning. \emph{Scientific Reports}, 14, 8124.
\doi{10.1038/s41598-024-58624-4}

Radicchi, F., Castellano, C., Cecconi, F., Loreto, V., & Parisi, D. (2004).
Defining and identifying communities in networks. \emph{Proceedings of the
National Academy of Sciences}, 101(9), 2658--2663.
\doi{10.1073/pnas.0400054101}
}
\seealso{\code{\link{net_louvain}}, \code{\link{net_modularity}},
\code{\link{net_plot_partition}}}
\examples{
g <- net_example("karate")
groups <- net_voronoi(g, threads = 1)
groups
attr(groups, "generators")
net_plot_partition(g, groups, seed = 1)

# Directed networks: the direction of the distances.
d <- net_example("krackhardt", relation = "advice")
net_voronoi(d, mode = "out", threads = 1)   # igraph cluster_voronoi(mode = "out")
net_voronoi(d, mode = "in", threads = 1)    # igraph mode = "in"
net_voronoi(d, mode = "all", threads = 1)   # igraph mode = "all"

# A fixed radius instead of the best one.
net_voronoi(g, radius = 2, threads = 1)
}
