\name{assortativity}
\alias{assortativity}
\alias{assortativity_nominal}
\alias{assortativity_degree}
\title{Assortativity coefficients}
\usage{
assortativity(x, values, values_in = NULL, directed = TRUE, normalized = TRUE,
              convention = c("igraph", "networkx", "homophily"))

assortativity_nominal(x, types, directed = TRUE, normalized = TRUE,
                      convention = c("igraph", "networkx", "homophily"))

assortativity_degree(x, directed = TRUE, weighted = FALSE, normalized = TRUE,
                     convention = c("igraph", "networkx", "homophily"))
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{values}{Numeric node values: a vector with one value per node or the
name of a column of \code{x$nodes}.}
\item{values_in}{Values at the receiving end of ties in directed networks,
given the same way; \code{NULL} uses \code{values}.}
\item{directed}{Whether to respect the direction of ties in directed
networks. \code{FALSE} treats every tie as reciprocal.}
\item{normalized}{Whether to return the coefficient (\code{TRUE}) or its
numerator: the covariance of the values at the two ends of ties for
\code{assortativity()}, and the modularity of the categories for
\code{assortativity_nominal()}.}
\item{convention}{How an undirected self-loop counts: \code{"igraph"} and
\code{"homophily"} treat it as two arcs (each undirected edge is a
reciprocal pair), \code{"networkx"} as one. \code{"homophily"} also ignores
direction in \code{assortativity_nominal()}, as homophily's
\code{assort_discrete()} does.}
\item{types}{Node categories: a vector with one value per node or the name of
a column of \code{x$nodes}. Missing values are not allowed.}
\item{weighted}{For \code{assortativity_degree()}: use strengths (sums of
weights) instead of degrees, as NetworkX does with \code{weight}.}
}
\value{A number.}
\description{
Newman's (2002, 2003) assortativity: the tendency of ties to join nodes with
similar values. \code{assortativity()} uses numeric node values,
\code{assortativity_nominal()} categories, and \code{assortativity_degree()}
degrees.
}
\details{
Every tie is read as arcs: an arc of a directed network once (or, with
\code{directed = FALSE}, in both directions), an undirected edge in both
directions. Parallel edges each count.

For categories, with \eqn{e_{ij}} the share of arcs from category \eqn{i} to
category \eqn{j}, \eqn{a_i = \sum_j e_{ij}} and \eqn{b_j = \sum_i e_{ij}},
\deqn{r = \frac{\sum_i e_{ii} - \sum_i a_i b_i}{1 - \sum_i a_i b_i}}{r =
(sum_i e_ii - sum_i a_i b_i) / (1 - sum_i a_i b_i)}
(Newman 2003, eq. 2). For numeric values, \eqn{r} is the Pearson correlation
of the values at the source and target of every arc; directed networks may
use different values at each end (\code{values_in}). Degree assortativity
correlates the out-degree of the source with the in-degree of the target, or
total degree when direction is ignored; using degree minus one (Newman's
"excess degree") gives the same correlation. Degrees count every edge end, so
an undirected self-loop adds 2.

Results are \code{NaN} when undefined: without ties, or when the values at one
end do not vary. They match igraph's \code{assortativity()},
\code{assortativity_nominal()}, and \code{assortativity_degree()}, NetworkX's
\code{numeric_assortativity_coefficient()},
\code{attribute_assortativity_coefficient()}, and
\code{degree_assortativity_coefficient()} (which return 0 instead of
\code{NaN} without ties), and homophily's \code{assort_scalar()},
\code{assort_discrete()} (which ignores direction, so it equals
\code{directed = FALSE} here), and \code{assort_degree()}.
}
\references{
Newman, M. E. J. (2002). Assortative mixing in networks. \emph{Physical
Review Letters}, 89(20), 208701. \doi{10.1103/PhysRevLett.89.208701}

Newman, M. E. J. (2003). Mixing patterns in networks. \emph{Physical Review
E}, 67(2), 026126. \doi{10.1103/PhysRevE.67.026126}
}
\seealso{
\code{\link{mixing_matrix}}, \code{\link{ei_index}}
}
\examples{
# An undirected network with a self-loop at node 5: igraph and homophily
# count the loop as two arcs, NetworkX as one.
g <- netframe(data.frame(from = c(1, 1, 2, 3, 4, 5), to = c(2, 3, 3, 4, 5, 5)),
              nodes = data.frame(id = 1:5, team = c("a", "a", "a", "b", "b"),
                                 age = c(30, 32, 35, 50, 52)),
              directed = FALSE)
assortativity_nominal(g, "team")                           # igraph
assortativity_nominal(g, "team", convention = "networkx")  # NetworkX
assortativity_nominal(g, "team", convention = "homophily") # homophily

assortativity(g, "age")                           # igraph
assortativity(g, "age", convention = "networkx")  # NetworkX
assortativity(g, "age", convention = "homophily") # homophily (as igraph)

assortativity_degree(g)                           # igraph
assortativity_degree(g, convention = "networkx")  # NetworkX
assortativity_degree(g, convention = "homophily") # homophily (as igraph)

# In a directed network homophily's assort_discrete() ignores direction.
d <- netframe(data.frame(from = c(1, 1, 1, 2, 3, 3), to = c(2, 3, 5, 4, 2, 4)),
              nodes = data.frame(id = 1:5, team = c("a", "a", "a", "b", "b")))
assortativity_nominal(d, "team")                           # igraph, NetworkX
assortativity_nominal(d, "team", convention = "homophily") # homophily
}
