\name{project}
\alias{project}
\title{One-mode projections of a two-mode network}
\usage{
project(x, which = c("first", "second", "both"), weighted = TRUE,
  method = c("count", "newman", "product"), loops = FALSE, threads = 0L)
}
\arguments{
\item{x}{A two-mode \code{\link{netframe}} (see \code{\link{twomode}}).}
\item{which}{\code{"first"}, \code{"second"}, or \code{"both"} modes.}
\item{weighted}{Whether the result has a \code{weight} column. Weights are
tie strengths (\code{weight_type = "strength"}), since sharing more partners
brings two nodes closer.}
\item{method}{How ties are weighted; see Details.}
\item{loops}{Whether to keep the diagonal as self-loops.}
\item{threads}{Number of threads; \code{0} uses all available cores.}
}
\value{
An undirected \code{\link{netframe}} whose nodes are the rows of
\code{x$nodes} in the chosen mode (all columns kept) and whose edges join
each pair of nodes with a shared partner once, ordered by node position.
For \code{which = "both"}, a list with elements \code{first} and
\code{second}.
}
\description{
Joins two nodes of the same mode when they share a partner in the other
mode: two people who attended the same event, or two events with a common
participant (Breiger 1974). \code{which = "both"} returns both projections,
the dual-projection view of Everett and Borgatti (2013).
}
\details{
Tie weights, chosen by \code{method}, are
\describe{
  \item{\code{"count"}}{the number of shared partners, the off-diagonal of
    \eqn{BB^T}{B B'} for the binary incidence matrix \eqn{B}. This is
    igraph's \code{bipartite_projection(multiplicity = TRUE)} and NetworkX's
    \code{weighted_projected_graph()}.}
  \item{\code{"newman"}}{Newman's (2001) collaboration weight
    \eqn{\sum_e 1/(d_e - 1)}{sum_e 1/(d_e - 1)} over shared partners
    \eqn{e}, where \eqn{d_e} is the number of partners of \eqn{e}, so that
    a large event ties each pair of its participants more weakly. This is
    NetworkX's \code{collaboration_weighted_projected_graph()}.}
  \item{\code{"product"}}{\eqn{\sum_e w_{ie} w_{je}}{sum_e w_ie w_je} over
    shared partners, the off-diagonal of \eqn{WW^T}{W W'} for the valued
    incidence matrix \eqn{W} (see \code{\link{incidence_matrix}}). It needs
    edge weights, which are used as given whatever the network's
    \code{weight_type}.}
}
For \code{"count"} and \code{"newman"} parallel edges count once (shared
partners are distinct partners); for \code{"product"} their weights add up,
as in \code{\link{incidence_matrix}}. Direction is ignored.

With \code{loops = TRUE} the diagonal is kept as self-loops: each node's
number of partners for \code{"count"} and its sum of squared weights for
\code{"product"}, so the projection holds all of \eqn{BB^T}{B B'}, as the
dual projection approach uses. Newman's weights have no diagonal.

The projection is computed with a sparse accumulator for each node, in
parallel, in time proportional to the number of two-step paths.
}
\references{
Breiger, R. L. (1974). The duality of persons and groups. \emph{Social
Forces}, 53(2), 181-190. \doi{10.2307/2576011}

Newman, M. E. J. (2001). Scientific collaboration networks. II. Shortest
paths, weighted networks, and centrality. \emph{Physical Review E}, 64(1),
016132. \doi{10.1103/PhysRevE.64.016132}

Everett, M. G., & Borgatti, S. P. (2013). The dual-projection approach for
two-mode networks. \emph{Social Networks}, 35(2), 204-210.
\doi{10.1016/j.socnet.2012.05.004}
}
\examples{
attended <- data.frame(
  person = c("ann", "ann", "bob", "bob", "cat", "dan"),
  event  = c("gala", "fair", "gala", "show", "show", "fair")
)
g <- twomode(attended)
project(g)$edges
project(g, "second", method = "newman")$edges
}
