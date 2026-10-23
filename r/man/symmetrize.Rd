\name{symmetrize}
\alias{symmetrize}
\alias{transpose_network}
\title{Symmetrize or transpose a directed network}
\description{
\code{symmetrize()} turns a directed network into an undirected one by
combining the two directions of every pair, as UCINET's Transform >
Symmetrize does with the method of the same name (Hanneman and Riddle 2005,
chapter 6):
\describe{
  \item{\code{"max"}}{A pair is tied if either direction is, with the larger
    weight (Boolean OR for binary ties).}
  \item{\code{"min"}}{A pair is tied only if both directions are, with the
    smaller weight (Boolean AND).}
  \item{\code{"sum"}}{The weights of the two directions add up (1 + 1 = 2
    for a reciprocated binary tie); the result is weighted.}
  \item{\code{"average"}}{Half the sum (an unreciprocated binary tie becomes
    1/2); the result is weighted.}
  \item{\code{"upper"}, \code{"lower"}}{The tie from the node that comes
    first (\code{"upper"}) or last (\code{"lower"}) in \code{x$nodes}
    decides, as copying the upper or lower triangle of the adjacency matrix
    onto the other.}
}
\code{transpose_network()} reverses every edge (UCINET: Transform >
Transpose).
}
\usage{
symmetrize(x, method = c("max", "min", "sum", "average", "upper", "lower"))

transpose_network(x)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{method}{How the two directions combine; see Description.}
}
\details{
An absent tie counts as weight 0, so for weighted networks \code{"max"},
\code{"sum"}, and \code{"average"} tie a pair when either direction is tied
and \code{"min"} when both are. Parallel edges in one direction collapse
first, keeping the largest weight (\code{"min"}, \code{"max"},
\code{"upper"}, \code{"lower"}) or adding up (\code{"sum"},
\code{"average"}). Self-loops are kept once per node, with the same rule.
The result's edge table has \code{from}, \code{to}, and, for weighted input
or the \code{"sum"} and \code{"average"} methods, \code{weight} (with the
network's \code{weight_type}, or \code{"strength"} for counts of binary
ties); other edge columns are dropped, since the two directions can carry
different values. Each pair appears once with \code{from} before \code{to}
in node order. \code{transpose_network()} keeps every column. An undirected
network is returned unchanged by \code{transpose_network()} and is an error
for \code{symmetrize()}.
}
\value{A netframe: undirected for \code{symmetrize()}.}
\references{
Hanneman, R. A., & Riddle, M. (2005). \emph{Introduction to Social Network
Methods}. University of California, Riverside.
}
\seealso{\code{\link{dichotomize}}, \code{\link{reciprocity}}}
\examples{
# 1 -> 2, 2 -> 1 (reciprocated), 1 -> 3, 3 -> 2.
d <- netframe(data.frame(from = c(1, 2, 1, 3), to = c(2, 1, 3, 2)))
symmetrize(d, "max")$edges      # any tie: 1-2, 1-3, 2-3
symmetrize(d, "min")$edges      # both directions: 1-2 only
symmetrize(d, "sum")$edges      # 1-2 weighs 2, the others 1
symmetrize(d, "average")$edges  # 1-2 weighs 1, the others 1/2
symmetrize(d, "upper")$edges    # the earlier node decides: 1-2, 1-3, but not 2-3 (3 -> 2)
symmetrize(d, "lower")$edges    # the later node decides: 1-2, 2-3, but not 1-3 (1 -> 3)

transpose_network(d)$edges
}
