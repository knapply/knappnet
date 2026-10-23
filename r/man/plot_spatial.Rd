\name{plot_spatial}
\alias{plot_spatial}
\alias{plot_edge_lengths}
\title{Plot a spatial network on its map}
\description{
Draws a spatial network with base graphics at its true locations: nodes
and edges at their geometries when the tables are sf (see
\code{\link{to_sf}}), and at the representative points otherwise. Polygons
are filled with the node color, lines drawn in it, points drawn as points,
and any other geometry as sf draws it; edges follow their own geometry when
they have one and are straight segments between representative points
otherwise, with arrowheads at the end of directed edges. Planar coordinates
keep an aspect ratio of 1; longitude and latitude get 1 / cos(mean
latitude), so shapes are not stretched. Networks without sf geometry
(coordinates recorded by \code{\link{spatial_network}}) need no sf.

\code{plot_edge_lengths()} draws the distribution of
\code{\link{edge_lengths}} as a histogram, in the unit of the lengths, and,
when \code{scale} is given, the \code{\link{distance_decay}} curve with that
scale on a second axis, so the choice of scale can be judged against the
lengths it will weight.

Colors come from palettes that are safe for color-blind readers: the
Okabe-Ito palette for categories (eight colors; more categories recycle
them, with a warning) and viridis for numbers.
}
\usage{
plot_spatial(x, node_color = NULL, node_size = NULL, node_label = NULL,
             edge_width = NULL, edge_color = NULL, ..., legend = TRUE)

plot_edge_lengths(x, scale = NULL, form = c("exponential", "gaussian", "power"),
                  exponent = 2, ...)
}
\arguments{
\item{x}{A spatial \code{\link{netframe}}: one with sf node geometry, or
one returned by \code{\link{spatial_network}} or
\code{\link{proximity_network}}.}
\item{node_color, edge_color}{Colors: the name of a column of \code{x$nodes}
or \code{x$edges} (categories or numbers, mapped to a palette with a
legend), a vector with one value per node or edge, or a single color.}
\item{node_size}{Node sizes for point nodes: a numeric column of
\code{x$nodes}, a vector with one value per node, or a single \code{cex};
the default grows with degree.}
\item{node_label}{Labels: \code{NULL} labels networks of at most 40 nodes
with their identifiers, \code{TRUE} always, \code{FALSE} never, or a column
of \code{x$nodes} or a vector with one label per node.}
\item{edge_width}{Line widths for edges: a numeric column of
\code{x$edges}, a vector with one value per edge, or a single \code{lwd}.}
\item{legend}{Whether to draw a legend for a mapped \code{node_color}.}
\item{...}{Further arguments to \code{\link[graphics]{title}}, such as
\code{main}; for \code{plot_edge_lengths()}, to \code{\link[graphics]{hist}}.}
\item{scale, form, exponent}{For \code{plot_edge_lengths()}, the
\code{\link{distance_decay}} curve to overlay; \code{NULL} draws none.}
}
\value{\code{x}, invisibly.}
\seealso{\code{\link{plot.netframe}} for abstract layouts,
\code{\link{to_sf}}, \code{\link{spatial_network}}.}
\examples{
nodes <- data.frame(name = c("London", "Paris", "Berlin", "Madrid"),
                    lon = c(-0.1276, 2.3522, 13.4050, -3.7038),
                    lat = c(51.5072, 48.8566, 52.5200, 40.4168),
                    region = c("north", "north", "north", "south"))
edges <- data.frame(from = c("London", "Paris", "Paris"),
                    to = c("Paris", "Berlin", "Madrid"))
g <- spatial_network(netframe(edges, nodes, directed = FALSE), radius = 6371.0088)
plot_spatial(g, node_color = "region", main = "Four cities")
plot_edge_lengths(g, scale = 500)
if (requireNamespace("sf", quietly = TRUE)) {
  sq <- function(x, y) sf::st_polygon(list(rbind(c(x, y), c(x + 1, y), c(x + 1, y + 1),
                                                 c(x, y + 1), c(x, y))))
  regions <- sf::st_sf(region = c("sw", "se", "nw", "ne"), pop = c(10, 40, 20, 30),
                       geometry = sf::st_sfc(sq(0, 0), sq(1, 0), sq(0, 1), sq(1, 1)))
  plot_spatial(from_polygons(regions, id = "region"), node_color = "pop")
}
}
