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

Appearance arguments are read exactly as by \code{\link{plot.netframe}}
(columns, vectors, single values, or vectors wrapped in \code{I()}), with
the same palettes: Okabe-Ito colors for categories and viridis for
numbers.
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
\item{node_size, node_color, node_label, edge_width, edge_color}{Appearance
of nodes and edges, as in \code{\link{plot.netframe}}: the name of a column
of \code{x$nodes} or \code{x$edges}, a vector with one value per node or
edge, a single value for all, or a vector wrapped in \code{\link{I}()} to
use its values as they are. Node sizes apply to point nodes; polygons and
lines are filled or drawn in the node color.}
\item{legend}{Whether to draw legends for mapped columns and vectors.}
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
