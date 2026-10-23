\name{mixing_matrix}
\alias{mixing_matrix}
\title{Mixing matrices}
\description{
Counts ties by the categories of the nodes at their two ends. The categories
can come from one node attribute (the classical mixing matrix) or from two
different ones, such as node identifiers and a group attribute, which gives
each node's ties to each group.
}
\usage{
mixing_matrix(x, dim1, dim2 = NULL, direction = c("all", "out", "in"))
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{dim1}{Row categories: the name of a column of \code{x$nodes}, or a
vector with one value per node.}
\item{dim2}{Column categories, given the same way. \code{NULL} uses
\code{dim1}.}
\item{direction}{For directed networks: \code{"out"} counts ties from
\code{dim1} (rows) to \code{dim2} (columns), \code{"in"} ties received (see
Details), and \code{"all"} both. Ignored for undirected networks.}
}
\details{
The layout follows the R package homophily (Knapp 2019), so that results
carry over:
\itemize{
  \item Undirected networks: every edge counts in both directions. When
    \code{dim1} and \code{dim2} are the same attribute, the matrix is folded
    into its upper triangle: the diagonal holds twice the number of edges
    within each category, entry \eqn{(i, j)} for \eqn{i < j} the number of
    edges between categories \eqn{i} and \eqn{j} counted from both sides
    (twice the number of edges), and the lower triangle is zero. The matrix
    sums to twice the number of edges.
  \item Directed networks with one attribute: rows are the categories of
    senders and columns those of receivers, so each arc counts once;
    \code{direction = "in"} transposes the matrix.
  \item Directed networks with two attributes: \code{"out"} counts each arc
    by \code{dim1} of its sender and \code{dim2} of its receiver;
    \code{"in"} gives a matrix with \code{dim2} categories as rows and
    \code{dim1} categories as columns, counting each arc by \code{dim2} of
    its sender and \code{dim1} of its receiver; \code{"all"} counts every
    arc in both directions.
}
Categories appear in order of first appearance in \code{x$nodes}. Ties
involving a node whose category is missing are not counted. Parallel edges
each count, and a self-loop counts like any other tie (twice in an
undirected network). Weights are not used. homophily also folds a square
matrix built from two different attributes; knappnet folds only when both
dimensions are the same attribute.
}
\value{A numeric matrix with dimnames.}
\references{
Newman, M. E. J. (2003). Mixing patterns in networks. \emph{Physical Review
E}, 67(2), 026126. \doi{10.1103/PhysRevE.67.026126}

Knapp, B. G. (2019). homophily: Measuring network homophily. R package,
\url{https://github.com/knapply/homophily}.
}
\seealso{
\code{\link{ei_index}}, \code{\link{assortativity_nominal}}
}
\examples{
staff <- data.frame(name = c("ann", "bob", "cat", "dan"),
                    dept = c("sales", "sales", "eng", "eng"))
g <- netframe(data.frame(from = c("ann", "ann", "bob", "cat"),
                         to = c("bob", "cat", "dan", "dan")),
              nodes = staff, directed = FALSE)
mixing_matrix(g, "dept")
mixing_matrix(g, "name", "dept")
}
