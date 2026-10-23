\name{net_simplify}
\alias{net_simplify}
\alias{net_is_simple}
\alias{net_multiplicity}
\title{Self-loops, parallel edges, and simplification}
\description{
A network is simple when it has neither self-loops nor parallel edges.
\code{net_multiplicity()} counts, for every edge, the edges that share its
key; \code{net_is_simple()} tests for simplicity; and
\code{net_simplify()} collapses each group of parallel edges into one,
combining the weights, and optionally drops self-loops.
}
\usage{
net_simplify(x, combine = c("sum", "mean", "min", "max", "first", "count"),
             remove_loops = FALSE)

net_is_simple(x)

net_multiplicity(x)
}
\arguments{
\item{x}{A \code{\link[=net_frame]{net_frame}}.}
\item{combine}{How the weights of parallel edges are combined:
\code{"sum"} (the default; igraph's \code{simplify()} does the same),
\code{"mean"}, \code{"min"}, \code{"max"}, \code{"first"} (the first
edge's weight), or \code{"count"} (the number of edges collapsed, which
becomes the weight of an unweighted network as a strength, in a column
called \code{count}).}
\item{remove_loops}{Whether to drop self-loops as well.}
}
\details{
Two edges are parallel when they have the same endpoints (as an ordered
pair in a directed network, unordered otherwise) \emph{and} the same values
of every role that distinguishes edges: the relation of a multiplex
network and the time spell of a temporal one. A friendship tie and an
advice tie between the same two people are not parallel, and neither are
two spells of one tie. This differs from igraph, NetworkX's
\code{MultiGraph}, and the network package, which compare endpoints only,
and keeps "multiplex" for several relations and "parallel" for repeated
pairs (Newman 2018, ch. 6, for the terms).

\code{net_simplify()} keeps the first row of each group of parallel edges
with its attributes, except for the weight column, which receives the
combined value; the kept edges stay in their original order. The result
declares \code{multiple = FALSE} (and \code{loops = FALSE} when loops were
removed), so the property is checked when the tables are edited.
}
\value{
\code{net_multiplicity()} returns an integer vector aligned with the rows
of \code{x$edges}: 1 for an edge with no parallel, and the size of its
group otherwise. \code{net_is_simple()} returns \code{TRUE} or
\code{FALSE}. \code{net_simplify()} returns a net_frame.
}
\references{
Newman, M. E. J. (2018). \emph{Networks} (2nd ed.). Oxford University
Press. \doi{10.1093/oso/9780198805090.001.0001}
}
\seealso{
\code{\link{net_frame}} for the \code{loops} and \code{multiple}
declarations and \code{\link{net_combine_relations}} for collapsing the
relations of a multiplex network.
}
\examples{
e <- data.frame(source = c("a", "a", "a", "b", "b"), target = c("b", "b", "a", "c", "c"),
                w = c(1, 2, 5, 1, 1), relation = c("x", "x", "x", "x", "y"))
g <- net_frame(e, directed = FALSE, weight = "w", weight_type = "strength",
               relation = "relation")
net_multiplicity(g)          # b-c x and b-c y are not parallel
net_is_simple(g)
s <- net_simplify(g, remove_loops = TRUE)
s$edges
net_is_simple(s)
}
