\name{net_fluid}
\alias{net_fluid}
\title{Fluid communities}
\description{
Finds exactly \code{k} communities by letting \code{k} fluids expand and
compete for the nodes of a network (Pares et al. 2018).
}
\usage{
net_fluid(x, k, seed = NULL)
}
\arguments{
\item{x}{A connected, undirected \code{\link{net_frame}}.}
\item{k}{Number of communities, a whole number from 1 to the number of
nodes.}
\item{seed}{Seed for the starting nodes and the sweeps, a whole number
between 0 and \eqn{2^{53}}; \code{NULL} draws one from R's random number
generator, so \code{\link{set.seed}()} makes results reproducible.}
}
\details{
The \code{k} communities start on \code{k} distinct random nodes, each with
density 1; the density of a community is 1 over its number of nodes. In
every sweep the nodes are visited in a new random order, and a node joins
the community with the largest total density over the node and its
neighbors, keeping its own community when that is among the largest and
otherwise choosing uniformly among them. The method stops after a sweep in
which no node changes (at most 100,000 sweeps). No community ever
disappears, so the result has exactly \code{k} groups.

Weights are not used. Parallel edges count once and self-loops are ignored
(igraph's \code{cluster_fluid_communities()} requires a simple network;
this is the same neighborhood). Directed networks and disconnected networks
are errors rather than being silently symmetrized or split. Densities are
compared exactly as fractions; NetworkX's \code{asyn_fluidc()} compares
them to a tolerance of \eqn{10^{-4}}, which makes no difference unless
communities exceed about 100 nodes. Unlike igraph, which requires
\eqn{k < n}, \eqn{k = n} is accepted. The same \code{seed} gives the same
result on every platform, but not the same partition as the other tools.
}
\value{An integer vector of group labels aligned with the rows of
\code{x$nodes}, numbered 1, 2, ... in order of first appearance, with
exactly \code{k} distinct values.}
\references{
Pares, F., Garcia-Gasulla, D., Vilalta, A., Moreno, J., Ayguade, E.,
Labarta, J., Cortes, U., & Suzumura, T. (2018). Fluid communities: A
competitive, scalable and diverse community detection algorithm. In C.
Cherifi, H. Cherifi, M. Karsai, & M. Musolesi (Eds.), \emph{Complex
Networks & Their Applications VI}, Studies in Computational Intelligence
689, 229--240. Springer. \doi{10.1007/978-3-319-72150-7_19}
}
\seealso{\code{\link{net_louvain}}, \code{\link{net_compare_partitions}},
\code{\link{net_plot_partition}}}
\examples{
g <- net_example("karate")
groups <- net_fluid(g, k = 2, seed = 1)
table(groups, g$nodes$faction)
net_plot_partition(g, groups, seed = 1)

# Four communities.
table(net_fluid(g, k = 4, seed = 1))
}
