\name{net_proximity}
\alias{net_proximity}
\title{Proximity networks from locations}
\description{
Builds a network that joins located points to the points near them: every
pair within a distance (a fixed-radius or disk network), or each point to
its \code{k} nearest neighbours. Neighbours are found with a k-d tree
(Bentley 1975; Friedman, Bentley, and Finkel 1977), so networks of 100,000
points take seconds, and every candidate's distance is then computed
exactly, so the tree changes the speed but never the result.
}
\usage{
net_proximity(nodes, coords = c("lon", "lat"), crs = c("lonlat", "planar"),
                  radius = 6371008.7714, within = NULL, k = NULL,
                  knn = c("directed", "either", "mutual"), id = 1L,
                  threads = 0L, exact = FALSE)
}
\arguments{
\item{nodes}{A data frame of points, a spatial net_frame (its nodes and
coordinate record are used), or an sf data frame of any geometry type (see
Details).}
\item{coords}{The two columns of \code{nodes} holding the coordinates, by
name or position: longitude then latitude for \code{crs = "lonlat"}, x
then y for \code{"planar"}.}
\item{crs}{\code{"lonlat"} for longitude and latitude in degrees on a
sphere, or \code{"planar"} for Euclidean distances; see
\code{\link{net_spatial}}.}
\item{radius}{Radius of the sphere for \code{crs = "lonlat"}; the default
gives metres.}
\item{within}{Join every pair of distinct points at most this far apart
(in the units of the distances). Give \code{within} or \code{k}, not
both.}
\item{k}{Join each point to its \code{k} nearest other points.}
\item{knn}{For \code{k} nearest neighbours: \code{"directed"} keeps an arc
from each point to each of its neighbours; \code{"either"} makes an
undirected edge wherever at least one of the two points lists the other;
\code{"mutual"} only where both do.}
\item{id}{Column of \code{nodes} holding node identifiers, as in
\code{\link{net_frame}}. For sf nodes, the default is the first column other
than the geometry; a table with no other column gets a new column \code{id}
numbering the features.}
\item{threads}{Number of threads, a whole number; \code{0} uses all available
cores, or 2 during \verb{R CMD check} (see \sQuote{Threads} in
\code{\link{net_frame}}).}
\item{exact}{For sf nodes only: whether to measure distances between the
geometries themselves with \code{\link[sf]{st_distance}} (so two polygons
that touch are at distance 0) rather than between their representative
points. Exact distances with \code{k} compute every pair, so they are
limited to 5000 features.}
}
\details{
Within-distance networks are undirected, contain one edge per pair at
distance at most \code{within} (so coincident points are joined), and
list each edge with the point that comes first in \code{nodes} first.
Nearest neighbours are listed in order of increasing distance, ties broken
by position in \code{nodes}, so the result is reproducible; a point with
fewer than \code{k} other points gets them all. The edges carry a
\code{distance} column, which is the network's weight column (as
distances), so path-based measures use spatial lengths unless asked not
to.

For an sf node table (or a net_frame with sf node geometry), when
\code{coords} is not given, the geometries are reduced to representative
points (see \code{\link{net_to_sf}}), longitude and latitude are recognized
from the CRS (with the default Earth radius, so distances are in metres),
and planar coordinates give distances in the units of the CRS. The result
keeps the sf node table and its geometry roles. With \code{exact = TRUE}
the distances are those between the geometries themselves.

For longitude and latitude the tree searches the points on the unit
sphere, where the chord between two points increases with their
great-circle distance, so it finds the same neighbours as the exact
metric.
}
\value{A spatial net_frame (see \code{\link{net_spatial}}) with edge
columns \code{from}, \code{to}, and \code{distance}.}
\references{
Bentley, J. L. (1975). Multidimensional binary search trees used for
associative searching. \emph{Communications of the ACM}, 18(9), 509-517.
\doi{10.1145/361002.361007}

Friedman, J. H., Bentley, J. L., & Finkel, R. A. (1977). An algorithm for
finding best matches in logarithmic expected time. \emph{ACM Transactions
on Mathematical Software}, 3(3), 209-226. \doi{10.1145/355744.355745}
}
\seealso{\code{\link{net_spatial}}, \code{\link{net_edge_lengths}}.}
\examples{
set.seed(1)
pts <- data.frame(id = 1:50, x = runif(50), y = runif(50))
g <- net_proximity(pts, coords = c("x", "y"), crs = "planar", within = 0.2)
g
h <- net_proximity(pts, coords = c("x", "y"), crs = "planar", k = 3, knn = "mutual")
summary(net_edge_lengths(h))
}
