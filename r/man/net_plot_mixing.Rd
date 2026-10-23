\name{net_plot_mixing}
\alias{net_plot_mixing}
\title{Plot a mixing matrix}
\description{
Draws the ties between categories of nodes as a heat map with the counts
written in the cells: the mixing matrix of \code{\link{net_mixing}}, or
any square or rectangular table of counts.
}
\usage{
net_plot_mixing(x, dim1 = NULL, dim2 = NULL, direction = c("all", "out", "in"),
                   values = TRUE, ...)
}
\arguments{
\item{x}{A \code{\link{net_frame}}, or a numeric matrix as returned by
\code{\link{net_mixing}}.}
\item{dim1, dim2, direction}{For a netframe, passed to
\code{\link{net_mixing}}.}
\item{values}{Whether to write the counts in the cells.}
\item{...}{Further arguments passed to \code{\link[graphics]{image}}, such
as \code{main}.}
}
\details{
Darker cells hold more ties. For an undirected network with one attribute
the matrix is folded into its upper triangle (see
\code{\link{net_mixing}}), so the lower triangle is 0; the diagonal
holds the ties within each category. Compare the diagonal with the rest to
judge homophily, or use \code{\link{net_assort_nominal}} and
\code{\link{net_ei_index}} for one number.
}
\value{The matrix drawn, invisibly.}
\seealso{\code{\link{net_mixing}}, \code{\link{net_plot_matrix}}.}
\examples{
g <- net_example("karate")
net_plot_mixing(g, "faction", main = "Ties by faction")
m <- net_mixing(g, "faction")
net_plot_mixing(m)
}
