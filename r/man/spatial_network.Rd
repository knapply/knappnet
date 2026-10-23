\name{spatial_network}
\alias{spatial_network}
\alias{edge_lengths}
\alias{distance_decay}
\alias{group_geometry}
\title{Geospatial networks}
\usage{
spatial_network(x, coords = c("lon", "lat"), crs = c("lonlat", "planar"),
                radius = 6371008.7714)

edge_lengths(x)

distance_decay(d, scale, form = c("exponential", "gaussian", "power"),
               exponent = 2)

group_geometry(x, groups)
}
\arguments{
\item{x}{A \code{\link{netframe}}; for every function but
\code{spatial_network()}, one returned by \code{spatial_network()} or
\code{\link{proximity_network}}.}
\item{coords}{The two columns of \code{x$nodes} holding the coordinates, by
name or position: longitude then latitude for \code{crs = "lonlat"}, x
then y for \code{"planar"}.}
\item{crs}{\code{"lonlat"} for longitude and latitude in degrees on a
sphere (distances are great-circle distances in the units of
\code{radius}), or \code{"planar"} for Euclidean distances in the units of
the coordinates.}
\item{radius}{Radius of the sphere for \code{crs = "lonlat"}. The default
is the mean Earth radius of the Geodetic Reference System 1980,
6371008.7714 m (Moritz 2000), so distances are in metres; use 6371.0088
for kilometres.}
\item{d}{Distances, such as the result of \code{edge_lengths()}.}
\item{scale}{The distance at which the weight has fallen to \eqn{e^{-1}}
(exponential and gaussian) or \eqn{2^{-e}} (power).}
\item{form}{The decay function; see Details.}
\item{exponent}{The exponent \eqn{e} of the \code{"power"} form.}
\item{groups}{The group of each node: a vector aligned with the rows of
\code{x$nodes}, or the name of a column of \code{x$nodes}.}
}
\value{\code{spatial_network()} returns \code{x} with the coordinate
columns and coordinate system recorded in \code{x$roles} (\code{x},
\code{y}, \code{crs}, \code{radius}). \code{edge_lengths()} returns a
numeric vector aligned with the rows of \code{x$edges};
\code{distance_decay()} a numeric vector the length of \code{d}.
\code{group_geometry()} returns a data frame with one row per group:
\code{group}, \code{size}, the centroid coordinates (named after the
coordinate columns), \code{mean_distance}, and \code{standard_distance}.}
\description{
A spatial network is a netframe whose nodes have coordinates: longitude
and latitude in degrees, or planar x and y. \code{spatial_network()}
records which node columns hold them; the other functions here use that
record. Spatial networks in general are reviewed by Barthelemy (2011).
}
\details{
Great-circle distances use the haversine formula (Sinnott 1984), which is
accurate to well under a millimetre on the Earth except between nearly
antipodal points, where it loses about half its digits (an error of up to
about 0.2 m). A sphere is an approximation: distances on it differ from
ellipsoidal geodesics by up to about 0.5\%. Coordinates must be finite and
latitudes within \eqn{[-90, 90]}.

\code{edge_lengths()} gives the distance between the endpoints of every
edge (0 for a self-loop); \code{summary()} or \code{hist()} of it describe
the edge-length distribution, and the lengths can serve as distance
weights in a new netframe. \code{distance_decay()} turns lengths into tie
weights that fall with distance, as in gravity models of spatial
interaction: \code{"exponential"} gives \eqn{\exp(-d / s)}{exp(-d / s)},
\code{"gaussian"} \eqn{\exp(-(d / s)^2)}{exp(-(d / s)^2)}, and
\code{"power"} \eqn{(1 + d / s)^{-e}}, with \code{scale} \eqn{s} and
\code{exponent} \eqn{e}; each is 1 at distance 0.

\code{group_geometry()} gives, for each group of nodes, its size, its
centroid, the mean distance of its members from the centroid, and the
standard distance (the root mean squared distance from the centroid; Bachi
1963). The centroid of planar points is the mean of their coordinates; for
longitude and latitude it is the geographic midpoint, the normalized mean
of the members' unit vectors converted back to longitude and latitude,
which is undefined (\code{NaN}) when that mean vanishes, as for two
antipodal points.
}
\references{
Barthelemy, M. (2011). Spatial networks. \emph{Physics Reports}, 499(1-3),
1-101. \doi{10.1016/j.physrep.2010.11.002}

Sinnott, R. W. (1984). Virtues of the haversine. \emph{Sky and Telescope},
68(2), 158-159.

Moritz, H. (2000). Geodetic Reference System 1980. \emph{Journal of
Geodesy}, 74(1), 128-133. \doi{10.1007/s001900050278}

Bachi, R. (1963). Standard distance measures and related methods for
spatial analysis. \emph{Papers of the Regional Science Association},
10(1), 83-132. \doi{10.1007/BF01934680}
}
\seealso{\code{\link{proximity_network}} to build a network from locations.}
\examples{
nodes <- data.frame(name = c("London", "Paris", "Berlin", "Madrid"),
                    lon = c(-0.1276, 2.3522, 13.4050, -3.7038),
                    lat = c(51.5072, 48.8566, 52.5200, 40.4168),
                    region = c("north", "north", "north", "south"))
edges <- data.frame(from = c("London", "Paris", "Paris"),
                    to = c("Paris", "Berlin", "Madrid"))
g <- spatial_network(netframe(edges, nodes, directed = FALSE),
                     coords = c("lon", "lat"), radius = 6371.0088)  # kilometres
edge_lengths(g)
distance_decay(edge_lengths(g), scale = 500)
group_geometry(g, "region")
}
