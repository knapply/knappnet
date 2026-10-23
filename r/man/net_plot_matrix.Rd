\name{net_plot_matrix}
\alias{net_plot_matrix}
\title{Plot an adjacency or incidence matrix}
\description{
Draws a network as its adjacency matrix (or, for a two-mode network, its
incidence matrix) with one shaded cell per tie, the way many network texts
present blockmodels and permuted matrices. Rows and columns can be sorted
by a partition of the nodes, with lines between the groups.
}
\usage{
net_plot_matrix(x, order = NULL, attr = NULL, values = FALSE, key = NULL,
            labels = NULL, ...)
}
\arguments{
\item{x}{A \code{\link{net_frame}}, or a numeric matrix.}
\item{order}{How to arrange rows and columns: \code{NULL} keeps the order
of \code{x$nodes}; a partition (the name of a column of \code{x$nodes}, or
a vector with one value per node, such as the result of
\code{\link{net_louvain}} or \code{\link{net_concor}}) sorts the nodes by
group, in order of the groups' first appearance, and draws a line between
groups. For a matrix, a vector with one value per row (and column, if
square).}
\item{attr}{For a net_frame, an edge column whose values fill the matrix;
\code{NULL} counts edges (see \code{\link{as.matrix.net_frame}}).}
\item{values}{Whether to write each cell's value in it.}
\item{key}{Whether to draw a color key. \code{NULL} draws one when the
matrix has more than two distinct values.}
\item{labels}{Whether to label rows and columns with node identifiers.
\code{NULL} labels them when there are at most 80.}
\item{...}{Further arguments passed to \code{\link[graphics]{image}}, such
as \code{main}.}
}
\details{
Cells are shaded from light (0) to dark (the largest value) with a
single-hue sequential palette, so a binary network shows ties as dark
cells on a light background. In a directed network row \eqn{i}, column
\eqn{j} is the tie from \eqn{i} to \eqn{j}; undirected networks give a
symmetric matrix. Parallel edges add up and a self-loop fills the
diagonal. A two-mode network (\code{\link{net_twomode}}) is drawn as its
\code{\link{net_incidence}}, first mode in rows, with the partition
applied to each mode's nodes.
}
\value{The matrix as drawn (after any reordering), invisibly.}
\seealso{\code{\link{plot.blockmodel}} for block densities,
\code{\link{net_plot_mixing}} for ties between categories.}
\examples{
g <- net_example("karate")
net_plot_matrix(g, main = "Karate club")
net_plot_matrix(g, order = "faction", main = "Sorted by faction")
net_plot_matrix(g, order = net_louvain(g, seed = 1))

women <- net_example("southern_women")
net_plot_matrix(women, main = "Southern women: incidence matrix")
}
