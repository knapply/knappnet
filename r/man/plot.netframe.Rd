\name{plot.netframe}
\alias{plot.netframe}
\title{Plot a network}
\description{
Draws a network with base R graphics. With no arguments beyond the
network, nodes are placed by \code{\link{layout_network}}, sized by degree,
and labelled when the network is small; edges are semi-transparent, with
arrows when the network is directed. Node and edge appearance can be set
from columns of \code{x$nodes} and \code{x$edges}.
}
\usage{
\method{plot}{netframe}(x, layout = NULL, node_size = NULL, node_color = NULL,
     node_label = NULL, edge_width = NULL, edge_color = NULL, ...,
     legend = TRUE, seed = NULL)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{layout}{Node positions: \code{NULL} for the default of
\code{\link{layout_network}} (or the node locations of a spatial network;
see Details), \code{"geo"} for the node locations, the name of a
\code{layout_network()} method (such as \code{"fr"} or \code{"circle"}),
the names of two numeric columns of \code{x$nodes} (for example
\code{c("lon", "lat")}), a matrix or data frame with two columns and one
row per node (as returned by \code{layout_network()}), or a function that
takes the network and returns such a matrix.}
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
\details{
Defaults:
\itemize{
  \item Node size grows with degree (as area), and nodes are smaller in
    larger networks.
  \item Nodes are a single blue with a white border, except in two-mode
    networks (\code{\link{twomode}}), where each mode has its own color.
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
(node area grows linearly) and colors from the sequential viridis palette
(\code{grDevices::hcl.colors(, "viridis")}; light for small and dark for
large values). Factors, character vectors, and logical vectors give colors
from the qualitative palette of Okabe and Ito (2008), as
\code{grDevices::palette.colors(, "Okabe-Ito")} ships it with black left
out, in the order of the factor levels or of the sorted values; beyond
eight categories the HCL palette \code{"Dark 3"} of Zeileis et al. (2020)
is used instead, which stays distinguishable in sequence but is not
guaranteed to be safe for every kind of color vision deficiency. Every
plot in the package uses these two palettes, so a category keeps its color
from one plot to the next. Missing values are drawn in light gray (colors)
or at the smallest size or width. Sizes and widths must be numeric.

Legends are drawn in the right margin. Coordinates keep their aspect
ratio, so distances in the plot are proportional to distances in the
layout. Self-loops are not drawn, and parallel edges are drawn on top of
each other.

The plot is reproducible for a given \code{seed} (or after
\code{\link{set.seed}}).

Spatial networks (those with sf node geometry, or coordinate columns
recorded by \code{\link{spatial_network}} or
\code{\link{proximity_network}}) are drawn at their node locations by
\code{\link{plot_spatial}} when \code{layout} is \code{NULL} or
\code{"geo"}, with the same appearance arguments; any other \code{layout}
draws them like any network. \code{"geo"} is an error for a network without
locations.
}
\value{The layout used, invisibly: a matrix with columns \code{x} and
\code{y} and one row per node, which can be passed back as \code{layout}.
For a spatial network drawn at its locations, \code{x} invisibly, as
\code{\link{plot_spatial}} returns.}
\references{
Okabe, M., & Ito, K. (2008). Color Universal Design (CUD): How to make
figures and presentations that are friendly to colorblind people.
\url{https://jfly.uni-koeln.de/color/}

Zeileis, A., Fisher, J. C., Hornik, K., Ihaka, R., McWhite, C. D., Murrell,
P., Stauffer, R., & Wilke, C. O. (2020). colorspace: A toolbox for
manipulating and assessing colors and palettes. \emph{Journal of
Statistical Software}, 96(1), 1-49. \doi{10.18637/jss.v096.i01}
}
\seealso{\code{\link{layout_network}}, \code{\link{plot_spatial}}, and the
other plotting functions: \code{\link{plot_matrix}},
\code{\link{plot_centrality}}, \code{\link{plot_partition}},
\code{\link{plot_degree_distribution}}.}
\examples{
g <- netframe(data.frame(from  = c("a", "a", "b", "c", "c", "d"),
                         to    = c("b", "c", "c", "d", "e", "e"),
                         calls = c(1, 4, 2, 8, 3, 5)),
              directed = FALSE)
g$nodes$group <- c("x", "x", "y", "y", "y")
plot(g, seed = 1)
plot(g, node_color = "group", edge_width = "calls", seed = 1, main = "Calls")
xy <- plot(g, layout = "circle", node_label = FALSE)
plot(g, layout = xy, node_size = I(c(3, 1, 1, 1, 1)))

attended <- data.frame(person = c("ann", "ann", "bob", "cat"),
                       event  = c("gala", "fair", "gala", "fair"))
plot(twomode(attended), layout = "twomode")
}
