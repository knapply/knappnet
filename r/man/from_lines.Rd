\name{from_lines}
\alias{from_lines}
\alias{from_polygons}
\title{Spatial networks from lines and polygons}
\description{
\code{from_lines()} builds a network from LINESTRING features such as
roads, rivers, or pipelines: each line becomes an edge, and its two
endpoints become nodes, shared with every other line that ends at the same
place. \code{from_polygons()} builds the contiguity network of regions:
each polygon becomes a node, and two polygons are joined when they share a
boundary. Both are standard constructions in spatial analysis (Barthelemy
2011; Bivand, Pebesma, and Gomez-Rubio 2013).
}
\usage{
from_lines(lines, tolerance = 0, directed = FALSE)

from_polygons(polygons, contiguity = c("queen", "rook"), id = NULL)
}
\arguments{
\item{lines}{An sf data frame or sfc of LINESTRING geometries (cast
MULTILINESTRINGs with \code{\link[sf]{st_cast}} first).}
\item{tolerance}{Distance within which line endpoints are the same node;
see the snapping rule.}
\item{directed}{Whether edges run from each line's first vertex to its last
(\code{TRUE}) or are undirected (\code{FALSE}, the default).}
\item{polygons}{An sf data frame or sfc of POLYGON or MULTIPOLYGON
geometries.}
\item{contiguity}{\code{"queen"} (shared boundary point) or \code{"rook"}
(shared boundary line).}
\item{id}{Column of \code{polygons} holding region identifiers, by name or
position, or \code{NULL} to number the regions 1, 2, ... in a new first
column \code{id} (an existing \code{id} column is renamed \code{id_attr}).}
}
\value{A netframe whose node table is an sf data frame; see
\code{\link{to_sf}}.}
\section{Snapping rule of from_lines()}{
The first and last vertex of every line are collected, in the order of the
lines (start, then end). With \code{tolerance = 0}, two endpoints are the
same node exactly when their coordinates are equal, as in
\code{sfnetworks::as_sfnetwork()}. With a positive \code{tolerance},
endpoints at most \code{tolerance} apart
(\code{\link[sf]{st_is_within_distance}}, in the units of the CRS: metres
for longitude and latitude) are the same node, and so, transitively, are
endpoints linked through others (single linkage, the connected components
of the within-tolerance relation). A tolerance that is small next to the
shortest line keeps chains from merging distinct junctions. Each node is
placed at the first endpoint that formed it and numbered by first
appearance; line geometries are kept as they are, so with a positive
tolerance an edge can end up to \code{tolerance} from its node, which
\code{\link{check_geometry}} with the same tolerance accepts.

Nodes get an integer \code{id} and POINT geometry. Edges keep every column
of \code{lines} (columns called \code{from}, \code{to}, or \code{length}
are renamed \code{from_attr}, \code{to_attr}, and \code{length_attr}), gain
\code{from} and \code{to}, and a \code{length} column from
\code{\link[sf]{st_length}} (see the Edge lengths section of
\code{\link{to_sf}}) that is the network's weight column, as distances, so
path-based measures use travel distance along the lines. A line that
starts and ends at the same node (a ring road) is a self-loop and is kept.
The network is undirected unless \code{directed = TRUE}, in which case each
edge runs from a line's first vertex to its last.
}
\section{Contiguity rule of from_polygons()}{
Two polygons are neighbours under \code{"queen"} contiguity when their
boundaries share at least one point, and under \code{"rook"} contiguity
when they share a boundary line, not just corners, as a rook and a queen
move on a chessboard (Anselin 1988). The relations are evaluated exactly by
\code{\link[sf]{st_relate}} with the DE-9IM patterns \code{F***T****}
(queen) and \code{F***1****} (rook) of Egenhofer and Franzosa (1991) and
Clementini, Di Felice, and van Oosterom (1993): interiors disjoint,
boundaries intersecting in anything (queen) or in a line (rook).
\code{spdep::poly2nb()} instead compares vertices: polygons are queen
neighbours when a vertex of one lies within \code{snap} (by default
\code{sqrt(.Machine$double.eps)}) of a vertex of the other, and rook
neighbours when more than one does. The two rules agree when shared
boundaries have shared vertices, as in most maps; they differ when a
shared boundary line ends at a vertex of only one polygon (a T-junction,
rook here but not in spdep) and when boundaries are drawn with slightly
different coordinates. Such polygons do not touch here, and should be snapped first (\code{\link[sf]{st_snap}}, or
rounding with \code{\link[sf]{st_set_precision}}). Overlapping polygons,
whose interiors intersect, are not neighbours under either rule.
Coordinates are treated as planar, also for longitude and latitude (as sf
and spdep do). The network is undirected, each pair listed once with the
polygon that comes first in \code{polygons} first, and a polygon with no
neighbours is an isolate. The edges have no geometry;
\code{\link{edge_lengths}} then measures between representative points
(see \code{\link{to_sf}}).
}
\references{
Anselin, L. (1988). \emph{Spatial Econometrics: Methods and Models}.
Kluwer. \doi{10.1007/978-94-015-7799-1}

Egenhofer, M. J., & Franzosa, R. D. (1991). Point-set topological spatial
relations. \emph{International Journal of Geographical Information
Systems}, 5(2), 161-174. \doi{10.1080/02693799108927841}

Clementini, E., Di Felice, P., & van Oosterom, P. (1993). A small set of
formal topological relationships suitable for end-user interaction. In
\emph{Advances in Spatial Databases} (pp. 277-295). Springer.
\doi{10.1007/3-540-56869-7_16}

Bivand, R. S., Pebesma, E., & Gomez-Rubio, V. (2013). \emph{Applied
Spatial Data Analysis with R} (2nd ed.). Springer.
\doi{10.1007/978-1-4614-7618-4}

Barthelemy, M. (2011). Spatial networks. \emph{Physics Reports}, 499(1-3),
1-101. \doi{10.1016/j.physrep.2010.11.002}
}
\seealso{\code{\link{to_sf}}, \code{\link{proximity_network}},
\code{\link{to_sfnetwork}}.}
\examples{
if (requireNamespace("sf", quietly = TRUE)) {
  # A small road grid: two streets crossing two avenues.
  seg <- function(x0, y0, x1, y1) sf::st_linestring(rbind(c(x0, y0), c(x1, y1)))
  roads <- sf::st_sf(name = c("1st", "2nd", "A", "B"),
                     geometry = sf::st_sfc(seg(0, 0, 2, 0), seg(0, 1, 2, 1),
                                           seg(0, 0, 0, 1), seg(2, 0, 2, 1)))
  g <- from_lines(roads)
  g$nodes
  g$edges
  edge_lengths(g)

  # Four square regions in a row, each touching the next.
  sq <- function(x) sf::st_polygon(list(rbind(c(x, 0), c(x + 1, 0), c(x + 1, 1),
                                              c(x, 1), c(x, 0))))
  regions <- sf::st_sf(region = c("w", "x", "y", "z"),
                       geometry = sf::st_sfc(sq(0), sq(1), sq(2), sq(3)))
  from_polygons(regions, id = "region")$edges
}
}
