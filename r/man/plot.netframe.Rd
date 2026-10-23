\name{plot.netframe}
\alias{plot.netframe}
\title{Plot a network}
\usage{
\method{plot}{netframe}(x, layout = NULL, node_size = NULL, node_color = NULL,
     node_label = NULL, edge_width = NULL, edge_color = NULL, ...,
     legend = TRUE, seed = NULL)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{layout}{Node positions: \code{NULL} for the default of
\code{\link{layout_network}}, the name of a \code{layout_network()} method
(such as \code{"fr"} or \code{"circle"}), the names of two numeric columns
of \code{x$nodes} (for example \code{c("lon", "lat")}), a matrix or data
frame with two columns and one row per node (as returned by
\code{layout_network()}), or a function that takes the network and returns
such a matrix.}
\item{node_size, node_color, node_label, edge_width, edge_color}{Appearance
of nodes and edges. Each can be the name of a column of \code{x$nodes} (for
node arguments) or \code{x$edges} (for edge arguments), a vector with one
value per node or edge, or a single value for all. Columns and vectors are
mapped to sizes, widths, or colors (see Details) and get a legend; wrap a
vector in \code{\link{I}()} to use its values as they are (sizes as
\code{cex} multipliers, widths as \code{lwd}, colors as R colors).
\code{NULL} uses the defaults (see Details). \code{node_label = TRUE}
labels nodes with their identifiers and \code{FALSE} draws no labels.}
\item{...}{Further arguments passed to \code{\link[graphics]{title}}, such
as \code{main}.}
\item{legend}{Whether to draw legends for mapped columns and vectors.}
\item{seed}{Seed for the layout, as in \code{\link{layout_network}}, when
\code{layout} is \code{NULL} or a method name.}
}
\value{The layout used, invisibly: a matrix with columns \code{x} and
\code{y} and one row per node, which can be passed back as \code{layout}.}
\description{
Draws a network with base R graphics. With no arguments beyond the
network, nodes are placed by \code{\link{layout_network}}, sized by degree,
and labelled when the network is small; edges are semi-transparent, with
arrows when the network is directed. Node and edge appearance can be set
from columns of \code{x$nodes} and \code{x$edges}.
}
\details{
Defaults:
\itemize{
  \item Node size grows with degree (as area), and nodes are smaller in
    larger networks.
  \item Nodes are a single blue with a white border.
  \item Nodes are labelled with their identifiers in networks of at most
    40 nodes.
  \item Edges are gray and more transparent the more edges there are. In
    directed networks, arrows end at the border of the target node.
  \item Edge width is constant, except in networks whose weights are
    strengths (\code{weight_type = "strength"} in \code{\link{netframe}}),
    where it grows with the weight. Distance weights do not change the
    width, since a thick line would suggest a strong tie.
}

Mapping: numeric values give sizes and widths that grow with the value
(node area grows linearly) and colors from the sequential palette
\code{grDevices::hcl.colors(, "viridis")}. Factors, character vectors, and
logical vectors give colors from the qualitative palette
\code{grDevices::hcl.colors(, "Dark 3")}, in the order of the factor levels
or of the sorted values. Missing values are drawn in light gray (colors) or
at the smallest size or width. Sizes and widths must be numeric.

Legends are drawn in the right margin. Coordinates keep their aspect
ratio, so distances in the plot are proportional to distances in the
layout. Self-loops are not drawn, and parallel edges are drawn on top of
each other.

The plot is reproducible for a given \code{seed} (or after
\code{\link{set.seed}}).
}
\seealso{\code{\link{layout_network}}}
\examples{
g <- netframe(data.frame(from = c("a", "a", "b", "c", "c"),
                         to   = c("b", "c", "c", "d", "e")),
              directed = FALSE)
g$nodes$group <- c("x", "x", "y", "y", "y")
plot(g, node_color = "group", seed = 1)
}
