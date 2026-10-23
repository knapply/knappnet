\name{net_spatial}
\alias{net_spatial}
\alias{net_is_spatial}
\alias{spatial_netframe}
\alias{print.spatial_netframe}
\alias{net_edge_lengths}
\alias{net_distance_decay}
\alias{net_group_geometry}
\title{Geospatial networks}
\description{
A spatial network is a netframe whose nodes have coordinates: longitude
and latitude in degrees, or planar x and y. \code{net_spatial()}
records which node columns hold them; the other functions here use that
record. Networks whose node table is an sf data frame (see
\code{\link{net_to_sf}}) need no such record: their geometries, of any type,
are reduced to representative points, and
\code{net_spatial(x, representative = )} chooses how. Either way the
network gets the class \code{spatial_netframe}, which
\code{net_is_spatial()} tests. Spatial networks in general are reviewed by
Barthelemy (2011).
}
\usage{
net_spatial(x, coords = c("lon", "lat"), crs = c("lonlat", "planar"),
            radius = 6371008.7714,
            representative = c("point_on_surface", "centroid"), quiet = FALSE)

net_is_spatial(x)

\method{print}{spatial_netframe}(x, ..., n = 6L, width = getOption("width"))

net_edge_lengths(x)

net_distance_decay(d, scale, form = c("exponential", "gaussian", "power"),
               exponent = 2)

net_group_geometry(x, groups)
}
\arguments{
\item{x}{A \code{\link[=net_frame]{netframe}}; for every function but
\code{net_spatial()} and \code{net_is_spatial()}, a spatial one returned by
\code{net_spatial()} or \code{\link{net_proximity}}, or one with sf node
geometry (see \code{\link{net_to_sf}}).}
\item{coords}{The names of the two columns of \code{x$nodes} holding the
coordinates: longitude then latitude for \code{crs = "lonlat"}, x then y
for \code{"planar"}. Not used, and an error if given, when \code{x$nodes}
carries sf geometry.}
\item{representative}{For sf node geometries other than points, how each
is reduced to one point: \code{\link[sf]{st_point_on_surface}} (the
default; always on the geometry) or \code{\link[sf]{st_centroid}}. See
\code{\link{net_to_sf}}.}
\item{crs}{\code{"lonlat"} for longitude and latitude in degrees on a
sphere (distances are great-circle distances in the units of
\code{radius}), or \code{"planar"} for Euclidean distances in the units of
the coordinates.}
\item{radius}{Radius of the sphere for \code{crs = "lonlat"}. The default
is the mean Earth radius of the Geodetic Reference System 1980,
6371008.7714 m (Moritz 2000), so distances are in metres; use 6371.0088
for kilometres.}
\item{quiet}{Whether to suppress the advisory notes of validation.}
\item{..., n, width}{As in \code{\link[=net_frame]{print.netframe}}.}
\item{d}{Distances, such as the result of \code{net_edge_lengths()}.}
\item{scale}{The distance at which the weight has fallen to \eqn{e^{-1}}
(exponential and gaussian) or \eqn{2^{-e}} (power).}
\item{form}{The decay function; see Details.}
\item{exponent}{The exponent \eqn{e} of the \code{"power"} form.}
\item{groups}{The group of each node: a vector aligned with the rows of
\code{x$nodes}, or the name of a column of \code{x$nodes}.}
}
\details{
Great-circle distances use the haversine formula (Sinnott 1984), which is
accurate to well under a millimetre on the Earth except between nearly
antipodal points, where it loses about half its digits (an error of up to
about 0.2 m). A sphere is an approximation: distances on it differ from
ellipsoidal geodesics by up to about 0.5\%. Coordinates must be finite and
latitudes within \eqn{[-90, 90]}.

\code{net_edge_lengths()} gives the distance between the endpoints of every
edge (0 for a self-loop); \code{summary()} or \code{hist()} of it describe
the edge-length distribution, and the lengths can serve as distance
weights in a new netframe. \code{net_distance_decay()} turns lengths into tie
weights that fall with distance, as in gravity models of spatial
interaction: \code{"exponential"} gives \eqn{\exp(-d / s)}{exp(-d / s)},
\code{"gaussian"} \eqn{\exp(-(d / s)^2)}{exp(-(d / s)^2)}, and
\code{"power"} \eqn{(1 + d / s)^{-e}}, with \code{scale} \eqn{s} and
\code{exponent} \eqn{e}; each is 1 at distance 0.

\code{net_group_geometry()} gives, for each group of nodes, its size, its
centroid, the mean distance of its members from the centroid, and the
standard distance (the root mean squared distance from the centroid; Bachi
1963). The centroid of planar points is the mean of their coordinates; for
longitude and latitude it is the geographic midpoint, the normalized mean
of the members' unit vectors converted back to longitude and latitude,
which is undefined (\code{NaN}) when that mean vanishes, as for two
antipodal points.

\strong{Classes.} A netframe whose node table has coordinate roles, or
whose node or edge table is an sf data frame, has class
\code{c("spatial_netframe", "netframe")}, or
\code{c("spatiotemporal_netframe", "spatial_netframe",
"temporal_netframe", "netframe")} when its edges also carry times. The
coordinate columns are validated by the engine (finite, latitudes within
[-90, 90]) when the network is built. Transformations that keep the
coordinate columns keep the class; dropping them drops it. The printed
card names the coordinate columns or the geometry types and the CRS, and
adds a line on how distances are measured.
}
\value{\code{net_spatial()} returns \code{x} with the coordinate
columns and coordinate system recorded in \code{x$roles} (\code{x},
\code{y}, \code{crs}, \code{radius}; for sf nodes \code{representative}
and \code{radius} only, since \code{crs} then holds the sf CRS) and the
class \code{spatial_netframe}. \code{net_is_spatial()} returns
\code{TRUE} or \code{FALSE}.
\code{net_edge_lengths()} returns a numeric vector aligned with the rows of
\code{x$edges} (from the edge geometries when the edge table has them, with
the unit in \code{attr(, "unit")}; see \code{\link{net_to_sf}});
\code{net_distance_decay()} a numeric vector the length of \code{d}.
\code{net_group_geometry()} returns a data frame with one row per group:
\code{group}, \code{size}, the centroid coordinates (named after the
coordinate columns), \code{mean_distance}, and \code{standard_distance}.}
\references{
Barthelemy, M. (2011). Spatial networks. \emph{Physics Reports}, 499(1-3),
1-101. \doi{10.1016/j.physrep.2010.11.002}

Sinnott, R. W. (1984). Virtues of the haversine. \emph{Sky and Telescope},
68(2), 159.

Moritz, H. (2000). Geodetic Reference System 1980. \emph{Journal of
Geodesy}, 74(1), 128-133. \doi{10.1007/s001900050278}

Bachi, R. (1963). Standard distance measures and related methods for
spatial analysis. \emph{Papers of the Regional Science Association},
10(1), 83-132. \doi{10.1007/BF01934680}
}
\seealso{\code{\link{net_proximity}} to build a network from locations.}
\examples{
nodes <- data.frame(id = c("London", "Paris", "Berlin", "Madrid"),
                    lon = c(-0.1276, 2.3522, 13.4050, -3.7038),
                    lat = c(51.5072, 48.8566, 52.5200, 40.4168),
                    region = c("north", "north", "north", "south"))
edges <- data.frame(source = c("London", "Paris", "Paris"),
                    target = c("Paris", "Berlin", "Madrid"))
g <- net_spatial(net_frame(edges, nodes, directed = FALSE),
                 coords = c("lon", "lat"), radius = 6371.0088)  # kilometres
g
net_edge_lengths(g)
net_distance_decay(net_edge_lengths(g), scale = 500)
net_group_geometry(g, "region")
}
