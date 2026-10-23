\name{to_sf}
\alias{to_sf}
\alias{check_geometry}
\alias{transform_network}
\alias{spatial_filter}
\alias{spatial_join}
\title{Spatial networks with sf}
\usage{
to_sf(x, what = c("nodes", "edges"), implicit = TRUE)

check_geometry(x, tolerance = 0)

transform_network(x, crs)

spatial_filter(x, y, predicate = sf::st_intersects, ...)

spatial_join(x, y, join = sf::st_intersects, ...)
}
\arguments{
\item{x}{A \code{\link{netframe}} with sf nodes or edges; for
\code{to_sf(what = "nodes")} also one with coordinate columns from
\code{\link{spatial_network}}.}
\item{what}{Which table to return: \code{"nodes"} or \code{"edges"}.}
\item{implicit}{For \code{what = "edges"}, whether edges without geometry
get straight-line geometries between their endpoints.}
\item{tolerance}{Largest distance, in the units of the CRS (metres for
longitude and latitude), at which an edge geometry's end still counts as
touching its node.}
\item{crs}{Target coordinate reference system, as accepted by
\code{\link[sf]{st_transform}} (an EPSG code, a string, or an sf
\code{crs}).}
\item{y}{Features to filter or join by: an sf or sfc object, an sf
\code{bbox}, or a spatial netframe (its node geometries are used).}
\item{predicate, join}{A geometric binary predicate from sf, such as
\code{\link[sf]{st_intersects}} (the default), \code{\link[sf]{st_within}},
or \code{\link[sf]{st_is_within_distance}}, called as
\code{predicate(nodes, y, ...)}.}
\item{...}{Further arguments passed to \code{predicate} or \code{join}, such
as \code{dist} for \code{\link[sf]{st_is_within_distance}}.}
}
\value{\code{to_sf()} returns an sf data frame. \code{check_geometry()}
returns a logical vector aligned with the rows of \code{x$edges}. The other
functions return a netframe: \code{transform_network()} and
\code{spatial_join()} with the same nodes and edges, \code{spatial_filter()}
the induced subnetwork.}
\description{
A spatial network is a \code{\link{netframe}} whose node table, edge table,
or both are simple-feature data frames from the sf package (Pebesma 2018).
Nodes and edges may carry any geometry type (POINT, LINESTRING, POLYGON,
their MULTI forms, or GEOMETRYCOLLECTION) and a coordinate reference
system; edges may also have no geometry at all, or geometry for only some
rows (empty geometries for the rest). Connectivity always comes from the
\code{from} and \code{to} columns, never from the geometry, so a line that
happens to cross a node does not join it. Geometry is handled entirely by
sf (GEOS and s2); the engine only ever sees numbers derived from it, such
as edge lengths and representative points.

\code{\link{netframe}} and \code{\link{as_netframe}} notice sf tables on
their own and record the geometry columns and CRS in \code{x$roles}
(\code{node_geometry}, \code{edge_geometry}, and \code{crs}, which then
holds an sf \code{crs} object). Nodes and edges that both have geometry
must share a CRS; use \code{\link[sf]{st_transform}} or
\code{\link[sf]{st_set_crs}} on one table first if they do not. The
functions here read those roles:
\itemize{
\item \code{to_sf()} returns the node or edge table as an sf data frame
with all its attributes. Nodes recorded with coordinate columns by
\code{\link{spatial_network}} become POINT features (WGS 84 for longitude
and latitude, no CRS for planar coordinates). Edges without geometry
become straight LINESTRINGs between the representative points of their
endpoints when \code{implicit = TRUE} (sfnetworks calls these spatially
implicit edges); with \code{implicit = FALSE} the edge table is returned as
it is, and an error is raised if it has no geometry.
\code{\link[sf]{st_as_sf}} on a netframe does the same.
\item \code{check_geometry()} is the optional consistency check between
edge geometries and node geometries. It is never required: a network is
valid however its geometries lie. For each edge with geometry it reports
whether the geometry starts within \code{tolerance} of the \code{from}
node's geometry and ends within \code{tolerance} of the \code{to} node's
geometry (in either orientation for undirected networks), \code{NA} for
edges without geometry. Distances are \code{\link[sf]{st_distance}}, so a
line ending anywhere inside a polygon node is at distance 0 from it.
\item \code{transform_network()} applies \code{\link[sf]{st_transform}} to
both tables at once and updates the recorded CRS. Columns computed earlier
from the geometry, such as lengths, are not recomputed.
\item \code{spatial_filter()} keeps the nodes whose geometry satisfies
\code{predicate} with any feature of \code{y} (any sf or sfc object, a
bbox, or another spatial netframe) and returns the induced subnetwork:
those nodes and every edge between two of them, with all columns and
roles, as \code{\link{ego_network}} does.
\item \code{spatial_join()} attaches the attributes of the feature of
\code{y} that each node's geometry matches under \code{join} (an
\code{\link[sf]{st_join}} with \code{left = TRUE}). The node table keeps
one row per node: a node that matches several features keeps the first
and a warning says how many nodes had more; a node that matches none gets
\code{NA}. Columns of \code{y} whose names clash with node columns get the
suffix \code{_y}.
}
}
\section{Representative points}{
Every point-based function of \code{\link{spatial_network}}
(\code{edge_lengths()} for edges without geometry,
\code{\link{group_geometry}}, \code{\link{proximity_network}}, the
map drawn by \code{\link{plot_spatial}}) reduces each node
geometry to one point. POINT geometries are used as they are. Other
geometries are reduced by \code{\link[sf]{st_point_on_surface}} by default,
which always lies on the geometry (a crescent's centroid does not), or by
\code{\link[sf]{st_centroid}} when
\code{spatial_network(x, representative = "centroid")} has been called. For
longitude and latitude, point on surface is computed on the planar
coordinates (sf warns of this; the warning is silenced here because the
point still lies on the geometry), while centroids are computed on the
sphere by s2. A node whose geometry is empty has no representative point
and those functions stop with an error naming the row.
}
\section{Edge lengths}{
When the edge table has geometry, \code{\link{edge_lengths}} returns
\code{\link[sf]{st_length}} of each edge geometry: a geodesic length for
longitude and latitude (computed on the sphere by s2 when
\code{\link[sf]{sf_use_s2}} is on, otherwise on the ellipsoid), and a
planar length in the units of the CRS otherwise. Edges with empty geometry
get the distance between the representative points of their endpoints
(\code{\link[sf]{st_distance}}, in the same units). The result is a plain
numeric vector with the unit recorded in \code{attr(, "unit")} (for example
\code{"m"}), empty when the CRS is unknown. Without edge geometry,
\code{\link{edge_lengths}} uses the representative points and the
haversine formula on a sphere of radius \code{x$roles$radius}, as
documented in \code{\link{spatial_network}}; the two differ by the small
difference between s2's Earth radius and the GRS 80 mean radius, and by
the ellipsoid when s2 is off.
}
\references{
Pebesma, E. (2018). Simple features for R: Standardized support for
spatial vector data. \emph{The R Journal}, 10(1), 439-446.
\doi{10.32614/RJ-2018-009}

Barthelemy, M. (2011). Spatial networks. \emph{Physics Reports}, 499(1-3),
1-101. \doi{10.1016/j.physrep.2010.11.002}
}
\seealso{\code{\link{from_lines}} and \code{\link{from_polygons}} to build
spatial networks from sf data, \code{\link{to_sfnetwork}} for sfnetworks,
\code{\link{read_network}} for GeoPackage files,
\code{\link{spatial_network}} and \code{\link{proximity_network}} for the
point-based functions.}
\examples{
if (requireNamespace("sf", quietly = TRUE)) {
  pts <- sf::st_as_sf(data.frame(name = c("a", "b", "c"), x = c(0, 1, 1), y = c(0, 0, 1)),
                      coords = c("x", "y"), crs = 32633)  # UTM zone 33N, metres
  g <- netframe(data.frame(from = c("a", "b", "c"), to = c("b", "c", "a")), pts,
                directed = FALSE)
  g
  edge_lengths(g)
  to_sf(g, "edges")
  h <- netframe(to_sf(g, "edges"), pts, directed = FALSE)  # explicit edge geometry
  check_geometry(h)
  transform_network(g, 4326)$roles$crs$epsg
  box <- sf::st_bbox(c(xmin = -0.5, ymin = -0.5, xmax = 0.5, ymax = 0.5), crs = 32633)
  spatial_filter(g, box)
}
}
