\name{to_sfnetwork}
\alias{to_sfnetwork}
\title{Convert to and from sfnetworks}
\usage{
to_sfnetwork(x, edges_as_lines = NULL, ...)
}
\arguments{
\item{x}{A \code{\link{netframe}} with node geometry or coordinate columns
(see \code{\link{to_sf}}).}
\item{edges_as_lines}{Whether edges without geometry should get
straight-line geometries; \code{NULL} keeps the edge table as it is.}
\item{...}{Further arguments to \code{sfnetworks::sfnetwork()}, such as
\code{force = TRUE}.}
}
\value{\code{to_sfnetwork()} returns an \code{sfnetwork};
\code{\link{as_netframe}} on an \code{sfnetwork} a netframe whose node table
(and edge table, when edges are spatially explicit) is an sf data frame.}
\description{
\code{to_sfnetwork()} turns a spatial netframe into an \code{sfnetwork}
object of the sfnetworks package (van der Meer et al. 2026);
\code{\link{as_netframe}} on an \code{sfnetwork} converts back, and so does
\code{sfnetworks::as_sfnetwork()} on a netframe. Node, edge, and
network-level attributes are carried across.
}
\details{
sfnetworks represents nodes as POINT features, so node geometries of other
types are replaced by their representative points (see \code{\link{to_sf}})
with a message. Node identifiers become the \code{name} node column that
tidygraph and igraph use (an existing \code{name} attribute is renamed
\code{name_attr}), and come back from it. Edge geometries are passed on
when the edge table has them; edges without geometry get straight lines
when \code{edges_as_lines = TRUE}, and remain spatially implicit by
default. \code{sfnetworks::sfnetwork()} checks that every edge line starts
and ends at its nodes' points and stops otherwise; pass \code{force = TRUE}
through \code{...} to skip that check for networks whose geometries do not
meet it (knappnet does not require it; see \code{\link{check_geometry}}).
Columns of class \code{units}, such as sfnetworks' \code{weight}, are
converted to plain numbers on the way back.
}
\references{
van der Meer, L., Abad, L., Gilardi, A., & Lovelace, R. (2026). sfnetworks:
Tidy geospatial networks. R package.
\doi{10.32614/CRAN.package.sfnetworks}
}
\seealso{\code{\link{to_sf}}, \code{\link{as_netframe}}.}
\examples{
if (requireNamespace("sfnetworks", quietly = TRUE)) {
  pts <- sf::st_as_sf(data.frame(name = c("a", "b", "c"), x = c(0, 1, 1), y = c(0, 0, 1)),
                      coords = c("x", "y"), crs = 32633)
  g <- netframe(data.frame(from = c("a", "b"), to = c("b", "c")), pts, directed = FALSE)
  net <- to_sfnetwork(g, edges_as_lines = TRUE)
  net
  as_netframe(net)
}
}
