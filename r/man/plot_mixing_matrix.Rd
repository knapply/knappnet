\name{plot_mixing_matrix}
\alias{plot_mixing_matrix}
\title{Plot a mixing matrix}
\description{
Draws the ties between categories of nodes as a heat map with the counts
written in the cells: the mixing matrix of \code{\link{mixing_matrix}}, or
any square or rectangular table of counts.
}
\usage{
plot_mixing_matrix(x, dim1 = NULL, dim2 = NULL, direction = c("all", "out", "in"),
                   values = TRUE, ...)
}
\arguments{
\item{x}{A \code{\link{netframe}}, or a numeric matrix as returned by
\code{\link{mixing_matrix}}.}
\item{dim1, dim2, direction}{For a netframe, passed to
\code{\link{mixing_matrix}}.}
\item{values}{Whether to write the counts in the cells.}
\item{...}{Further arguments passed to \code{\link[graphics]{image}}, such
as \code{main}.}
}
\details{
Darker cells hold more ties. For an undirected network with one attribute
the matrix is folded into its upper triangle (see
\code{\link{mixing_matrix}}), so the lower triangle is 0; the diagonal
holds the ties within each category. Compare the diagonal with the rest to
judge homophily, or use \code{\link{assortativity_nominal}} and
\code{\link{ei_index}} for one number.
}
\value{The matrix drawn, invisibly.}
\seealso{\code{\link{mixing_matrix}}, \code{\link{plot_matrix}}.}
\examples{
g <- example_network("karate")
plot_mixing_matrix(g, "faction", main = "Ties by faction")
m <- mixing_matrix(g, "faction")
plot_mixing_matrix(m)
}
