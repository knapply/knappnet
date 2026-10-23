\name{regular_equivalence}
\alias{regular_equivalence}
\title{Regular equivalence (CATREGE)}
\description{
Nodes are regularly equivalent when they have the same kinds of ties to
the same kinds of nodes: for every node tied to one of them, some node of
the same class is tied in the same way to the other (White and Reitz
1983). Unlike structural equivalence, the alters need not be the same
nodes. The largest such equivalence, the maximal regular equivalence, is
found by CATREGE (Borgatti and Everett 1993): start with every node in one
class (or with \code{seed}), and repeatedly split each class by the set of
(type of tie, class of alter) pairs its members have, until no class
splits.
}
\usage{
regular_equivalence(x, seed = NULL, weighted = NULL)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{seed}{A partition to refine: a vector aligned with the rows of
\code{x$nodes}, or the name of a column of \code{x$nodes}. Nodes in
different seed classes are never equivalent; only equality of values
matters. \code{NULL} starts from a single class.}
\item{weighted}{Whether to treat total edge weights as categories of tie:
nodes with ties of different weights are then not equivalent. \code{NULL}
uses weights when the network has a weight column. Weights are categories,
not magnitudes, so \code{weight_type} does not matter.}
}
\details{
The type of tie between \eqn{i} and \eqn{k} is the pair (tie from \eqn{i}
to \eqn{k}, tie from \eqn{k} to \eqn{i}), so the network and its transpose
are treated as one multiplex relation and a mutual tie differs from a
one-way tie, as in sna's \code{redist(method = "catrege")} and UCINET's
CATREGE. In an undirected network every tie has the same type. Self-loops
are ignored and parallel edges count once (with weights, their weights add
up).

In an undirected network without a seed partition the maximal regular
equivalence is degenerate: every node with at least one tie ends in one
class and isolated nodes in another, as Borgatti and Everett note. Classes
are numbered from 1 in order of their first node.

The validation suite compares the result with the classes at distance 0
in sna's \code{redist(method = "catrege")}, with and without
\code{seed.partition}, on every network it holds; UCINET's CATREGE
computes the same partition. REGE, the iterative algorithm for valued data that scores
the degree of regular equivalence between pairs, is not provided.
}
\value{An integer vector of class numbers (starting at 1) aligned with the
rows of \code{x$nodes}, with attribute \code{rounds}, the number of
refinement rounds that split a class.}
\references{
White, D. R., & Reitz, K. P. (1983). Graph and semigroup homomorphisms on
networks of relations. \emph{Social Networks}, 5(2), 193-234.
\doi{10.1016/0378-8733(83)90025-4}

Borgatti, S. P., & Everett, M. G. (1993). Two algorithms for computing
regular equivalence. \emph{Social Networks}, 15(4), 361-376.
\doi{10.1016/0378-8733(93)90012-A}

Everett, M. G., & Borgatti, S. P. (1994). Regular equivalence: General
theory. \emph{Journal of Mathematical Sociology}, 19(1), 29-52.
\doi{10.1080/0022250X.1994.9990134}
}
\seealso{\code{\link{structural_equivalence}},
\code{\link{automorphic_equivalence}}, \code{\link{blockmodel}}.}
\examples{
# A two-level hierarchy: 1 commands 2 and 3, who each command two others.
g <- netframe(data.frame(from = c(1, 1, 2, 2, 3, 3), to = c(2, 3, 4, 5, 6, 7)))
regular_equivalence(g)  # three classes: top, middle, bottom

# Undirected networks need a seed partition to give a useful answer.
u <- example_network("florentine")
regular_equivalence(u)
regular_equivalence(u, seed = concor(u))
}
