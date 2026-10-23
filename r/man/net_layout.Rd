\name{net_layout}
\alias{net_layout}
\alias{net_layout_stress}
\title{Node coordinates for drawing a network}
\description{
\code{net_layout()} places the nodes of a network in the plane.
\code{net_layout_stress()} measures how well a layout represents the network's
shortest-path distances.
}
\usage{
net_layout(x, method = c("auto", "stress", "sparse_stress", "fr", "circle",
                             "star", "grid", "tree", "twomode"),
               seed = NULL, weighted = NULL, weights_as = NULL, types = NULL,
               root = NULL, center = NULL, max_iter = NULL, tol = NULL,
               pivots = NULL, approximation = c("auto", "exact", "barnes_hut"),
               threads = 0L)

net_layout_stress(x, layout, weighted = NULL, weights_as = NULL, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{method}{The layout algorithm; see Details.}
\item{seed}{Seed for the random parts of a layout, a whole number between 0
and \eqn{2^{53}}; \code{NULL} draws one from R's random number generator, so
\code{\link{set.seed}()} makes results reproducible.}
\item{weighted}{Whether to use edge weights (stress and force-directed layouts
only). \code{NULL} (the default) uses them when the network has a weight
column.}
\item{weights_as}{Override the network's \code{weight_type}: \code{"distance"}
or \code{"strength"}.}
\item{types}{For \code{method = "twomode"}: which mode each node belongs
to, as the name of a column of \code{x$nodes} or a vector with one value
per node, with at most two distinct values. \code{NULL} uses the node type
of a network made by \code{\link{net_twomode}} or \code{\link{net_as_twomode}}, or
else a column named \code{type}. Nodes of the first value (\code{FALSE}, or
the first factor level or smallest value) form the top row.}
\item{root}{For \code{method = "tree"}: identifiers of root nodes, at most
one per component. \code{NULL} uses each component's center.}
\item{center}{For \code{method = "star"}: the identifier of the center
node. \code{NULL} uses the first node with the most neighbors.}
\item{max_iter}{Iteration cap for \code{"stress"}, \code{"sparse_stress"},
and \code{"fr"} (the number of iterations for \code{"fr"}). \code{NULL}
uses the defaults: 500 for stress on all pairs, 100 for sparse stress, and
500 for \code{"fr"}.}
\item{tol}{Convergence tolerance for stress layouts. \code{NULL} uses the
defaults (see Details).}
\item{pivots}{Number of pivots for sparse stress (default 100) and for the
starting layout of stress on all pairs (default 50). \code{NULL} uses the
defaults.}
\item{approximation}{For \code{"fr"}: \code{"auto"} computes repulsion
exactly for components of up to 1,000 nodes and with a Barnes-Hut quadtree
beyond, \code{"exact"} always computes it exactly, and \code{"barnes_hut"}
always approximates it.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{net_frame}}). Results do not depend on the number of threads.}
\item{layout}{A layout of \code{x}: a numeric matrix or data frame with two
columns and one row per node.}
}
\details{
Methods:
\describe{
  \item{\code{"auto"}}{Stress majorization (\code{"stress"}) for components
    of up to 2,000 nodes and sparse stress (\code{"sparse_stress"}) for
    larger ones, so the default stays fast on very large networks.}
  \item{\code{"stress"}}{Stress majorization (Gansner, Koren, and North
    2005), which places nodes so that their distances in the plane match
    their shortest-path distances as closely as possible. It minimizes
    stress, \eqn{\sum_{i<j} (\|x_i - x_j\| - d_{ij})^2 / d_{ij}^2}{sum over
    pairs of (||x_i - x_j|| - d_ij)^2 / d_ij^2}, the energy of Kamada and
    Kawai (1989), starting from pivot multidimensional scaling (Brandes and
    Pich 2007). It needs all shortest-path distances, so time per iteration
    and memory grow with the square of the component size; components of
    more than 10,000 nodes are refused. Iteration stops when stress falls by
    less than a relative 1e-4 in an iteration. Components of at most 100
    nodes are also laid out from four random starts, keeping the layout
    with the least stress, with a tolerance of 1e-5.}
  \item{\code{"sparse_stress"}}{The sparse stress model of Ortmann,
    Klimenta, and Brandes (2017): stress over edges and over distances to a
    random sample of pivot nodes, each pivot weighted by the part of the
    network it represents. Time per iteration grows linearly with the
    network. It approximates the stress objective; it never changes any
    measure.}
  \item{\code{"fr"}}{The force-directed algorithm of Fruchterman and
    Reingold (1991), from random starting positions. Repulsion between all
    pairs is computed exactly or, for large components, with the Barnes and
    Hut (1986) quadtree approximation (see \code{approximation}), which only
    changes the simulated forces. With weights, an edge's attraction is
    multiplied by its strength (the weight, or 1/weight for distances).}
  \item{\code{"circle"}}{Nodes on a circle in the order of
    \code{x$nodes}.}
  \item{\code{"star"}}{One node in the center, the others on a circle.}
  \item{\code{"grid"}}{Nodes on a square grid in the order of
    \code{x$nodes}, row by row from the top.}
  \item{\code{"tree"}}{A layered drawing of a breadth-first spanning tree
    of each component: depth increases downwards, leaves are one unit
    apart, and each parent is centered over its children. The default root
    is a center of the component (the middle of the path between the ends
    of a double breadth-first sweep).}
  \item{\code{"twomode"}}{Two rows for a two-mode (bipartite) network,
    given by \code{types}, ordered by the barycenter heuristic of Sugiyama,
    Tagawa, and Toda (1981) to reduce edge crossings.}
}

All methods ignore edge direction and self-loops and count parallel edges
once. In stress layouts, weights are lengths: distances as given, or
strengths converted to 1/weight; with several parallel edges the shortest
counts. Weights must be positive.

Stress, force-directed, and tree layouts treat each weak component
separately and pack the components in rows, largest first, separated by
one typical edge length; isolated nodes end up in rows of their own. Each
stress or force-directed component is centered and rotated so its longest
extent is horizontal. Coordinates are in units of path length for stress
layouts (hops, without weights) and of the ideal edge length for
\code{"fr"}.

Layouts are identical for a given \code{seed} on any number of threads.
Stress layouts use the seed only for a tiny perturbation of the starting
layout, the random restarts of small components, and the pivots of sparse
stress.

\code{net_layout_stress()} returns the stress of a layout,
\eqn{\sum_{i<j} (\|x_i - x_j\| - d_{ij})^2 / d_{ij}^2}{sum over pairs of
(||x_i - x_j|| - d_ij)^2 / d_ij^2} over pairs of nodes in the same
component, with shortest-path lengths \eqn{d_{ij}}{d_ij} computed as for
stress layouts. It depends on the scale of the layout: compare layouts at
the scale that minimizes it, or use layouts in path-length units.

Stress layouts of the validation networks have stress no higher than
NetworkX's \code{kamada_kawai_layout()} (within 2\%, and usually lower);
see \file{validation/layout_check.py} in the source repository.
}
\value{\code{net_layout()} returns a numeric matrix with columns
\code{x} and \code{y} and one row per node, aligned with the rows of
\code{x$nodes}. \code{net_layout_stress()} returns a number.}
\references{
Gansner, E. R., Koren, Y., & North, S. (2005). Graph drawing by stress
majorization. In \emph{Graph Drawing (GD 2004)}, Lecture Notes in Computer
Science 3383, 239-250. \doi{10.1007/978-3-540-31843-9_25}

Kamada, T., & Kawai, S. (1989). An algorithm for drawing general undirected
graphs. \emph{Information Processing Letters}, 31(1), 7-15.
\doi{10.1016/0020-0190(89)90102-6}

Brandes, U., & Pich, C. (2007). Eigensolver methods for progressive
multidimensional scaling of large data. In \emph{Graph Drawing (GD 2006)},
Lecture Notes in Computer Science 4372, 42-53.
\doi{10.1007/978-3-540-70904-6_6}

Ortmann, M., Klimenta, M., & Brandes, U. (2017). A sparse stress model.
\emph{Journal of Graph Algorithms and Applications}, 21(5), 791-821.
\doi{10.7155/jgaa.00440}

Fruchterman, T. M. J., & Reingold, E. M. (1991). Graph drawing by
force-directed placement. \emph{Software: Practice and Experience}, 21(11),
1129-1164. \doi{10.1002/spe.4380211102}

Barnes, J., & Hut, P. (1986). A hierarchical O(N log N) force-calculation
algorithm. \emph{Nature}, 324, 446-449. \doi{10.1038/324446a0}

Sugiyama, K., Tagawa, S., & Toda, M. (1981). Methods for visual
understanding of hierarchical system structures. \emph{IEEE Transactions on
Systems, Man, and Cybernetics}, 11(2), 109-125.
\doi{10.1109/TSMC.1981.4308636}
}
\seealso{\code{\link{plot.netframe}}}
\examples{
# A triangle with a pendant node, a separate pair, and an isolate.
g <- net_frame(data.frame(source = c(1, 2, 3, 3, 5), target = c(2, 3, 1, 4, 6)),
              nodes = data.frame(id = 1:7), directed = FALSE)
xy <- net_layout(g, seed = 1, threads = 2)
xy
net_layout_stress(g, xy, threads = 2)
net_layout(g, "fr", seed = 1, threads = 2)
net_layout(g, "circle")
net_layout(g, "tree", root = c(3, 5))

# Two rows for a two-mode network.
attended <- data.frame(person = c("ann", "ann", "bob", "cat"),
                       event  = c("gala", "fair", "gala", "fair"))
net_layout(net_twomode(attended), "twomode")
}
