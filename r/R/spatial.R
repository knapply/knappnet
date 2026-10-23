#' Geospatial networks
#'
#' A spatial network is a netframe whose nodes have coordinates:
#' longitude and latitude in degrees, or planar x and y. `spatial_network()`
#' records which node columns hold them; the other functions here use that
#' record. Spatial networks in general are reviewed by Barthelemy (2011).
#'
#' @param x A [netframe()]; for every function but `spatial_network()` and
#'   `proximity_network()`, one returned by `spatial_network()`.
#' @param coords The two columns of `x$nodes` (or of `nodes`) holding the
#'   coordinates, by name or position: longitude then latitude for
#'   `crs = "lonlat"`, x then y for `"planar"`.
#' @param crs `"lonlat"` for longitude and latitude in degrees on a sphere
#'   (distances are great-circle distances in the units of `radius`), or
#'   `"planar"` for Euclidean distances in the units of the coordinates.
#' @param radius Radius of the sphere for `crs = "lonlat"`. The default is
#'   the mean Earth radius of the Geodetic Reference System 1980,
#'   6371008.7714 m (Moritz 2000), so distances are in metres; use
#'   6371.0088 for kilometres.
#'
#' @details
#' Great-circle distances use the haversine formula (Sinnott 1984), which
#' is accurate to well under a millimetre on the Earth except between
#' nearly antipodal points, where it loses about half its digits (an error
#' of up to about 0.2 m). A sphere is an approximation: distances on it
#' differ from ellipsoidal geodesics by up to about 0.5%. Coordinates must
#' be finite and latitudes within \[-90, 90\].
#'
#' `edge_lengths()` gives the distance between the endpoints of every
#' edge (0 for a self-loop); `summary()` or `hist()` of it describe the
#' edge-length distribution, and the lengths can serve as distance weights
#' in a new netframe. `distance_decay()` turns lengths into tie weights
#' that fall with distance, as in gravity models of spatial interaction:
#' `"exponential"` gives \eqn{\exp(-d / s)}, `"gaussian"`
#' \eqn{\exp(-(d / s)^2)}, and `"power"` \eqn{(1 + d / s)^{-e}}, with
#' `scale` \eqn{s} and `exponent` \eqn{e}; each is 1 at distance 0.
#'
#' `group_geometry()` gives, for each group of nodes, its size, its
#' centroid, the mean distance of its members from the centroid, and the
#' standard distance (the root mean squared distance from the centroid;
#' Bachi 1963). The centroid of planar points is the mean of their
#' coordinates; for longitude and latitude it is the geographic midpoint,
#' the normalized mean of the members' unit vectors converted back to
#' longitude and latitude, which is undefined (`NaN`) when that mean
#' vanishes, as for two antipodal points.
#'
#' @references
#' Barthelemy, M. (2011). Spatial networks. *Physics Reports*, 499(1-3),
#' 1-101. \doi{10.1016/j.physrep.2010.11.002}
#'
#' Sinnott, R. W. (1984). Virtues of the haversine. *Sky and Telescope*,
#' 68(2), 159.
#'
#' Moritz, H. (2000). Geodetic Reference System 1980. *Journal of Geodesy*,
#' 74(1), 128-133. \doi{10.1007/s001900050278}
#'
#' Bachi, R. (1963). Standard distance measures and related methods for
#' spatial analysis. *Papers of the Regional Science Association*, 10(1),
#' 83-132. \doi{10.1007/BF01934680}
#'
#' @return `spatial_network()` returns `x` with the coordinate columns and
#'   coordinate system recorded in `x$roles` (`x`, `y`, `crs`, `radius`).
#'   `edge_lengths()` returns a numeric vector aligned with the rows of
#'   `x$edges`; `distance_decay()` a numeric vector the length of `d`.
#'   `group_geometry()` returns a data frame with one row per group:
#'   `group`, `size`, the centroid coordinates (named after the coordinate
#'   columns), `mean_distance`, and `standard_distance`.
#' @seealso [proximity_network()] to build a network from locations.
#' @examples
#' nodes <- data.frame(name = c("London", "Paris", "Berlin", "Madrid"),
#'                     lon = c(-0.1276, 2.3522, 13.4050, -3.7038),
#'                     lat = c(51.5072, 48.8566, 52.5200, 40.4168),
#'                     region = c("north", "north", "north", "south"))
#' edges <- data.frame(from = c("London", "Paris", "Paris"),
#'                     to = c("Paris", "Berlin", "Madrid"))
#' g <- spatial_network(netframe(edges, nodes, directed = FALSE),
#'                      coords = c("lon", "lat"), radius = 6371.0088)  # kilometres
#' edge_lengths(g)
#' distance_decay(edge_lengths(g), scale = 500)
#' group_geometry(g, "region")
#' @export
spatial_network <- function(x, coords = c("lon", "lat"), crs = c("lonlat", "planar"),
                            radius = 6371008.7714) {
  check_netframe(x)
  crs <- match.arg(crs)
  if (length(coords) != 2L) {
    stop("`coords` must name two columns of `x$nodes`: longitude and latitude, or x and y.",
         call. = FALSE)
  }
  cx <- resolve_column(x$nodes, coords[[1L]], "coords[1]", "nodes")
  cy <- resolve_column(x$nodes, coords[[2L]], "coords[2]", "nodes")
  if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) || !is.finite(radius) ||
      radius <= 0) {
    stop("`radius` must be a single positive number.", call. = FALSE)
  }
  check_coordinates(x$nodes[[cx]], x$nodes[[cy]], cx, cy, crs, "x$nodes")
  x$roles$x <- cx
  x$roles$y <- cy
  x$roles$crs <- crs
  x$roles$radius <- as.double(radius)
  x
}

check_coordinates <- function(cx_values, cy_values, cx, cy, crs, table) {
  for (col in list(list(cx, cx_values), list(cy, cy_values))) {
    v <- col[[2L]]
    if (!is.numeric(v)) {
      stop(sprintf("Coordinate column `%s` of `%s` must be numeric, not %s.", col[[1L]],
                   table, describe_class(v)), call. = FALSE)
    }
    bad <- which(!is.finite(v))
    if (length(bad)) {
      stop(sprintf("Coordinate column `%s` of `%s` has missing or non-finite values in %s.",
                   col[[1L]], table, rows_text(bad)), call. = FALSE)
    }
  }
  if (crs == "lonlat") {
    bad <- which(cy_values < -90 | cy_values > 90)
    if (length(bad)) {
      stop(sprintf("Latitude column `%s` of `%s` has values outside [-90, 90] in %s.", cy,
                   table, rows_text(bad)), call. = FALSE)
    }
  }
  invisible(NULL)
}

# Node coordinates of a spatial network: an n-by-2 matrix (longitude and
# latitude, or x and y), whether they are longitude and latitude, the
# sphere radius, and the column names. Every point-based function reads
# coordinates through this helper, so that node geometries of other kinds
# can later be reduced to representative points in one place.
node_coordinates <- function(x) {
  check_netframe(x)
  if (is.null(x$roles$x) || is.null(x$roles$y) || is.null(x$roles$crs)) {
    stop("`x` has no coordinates; record them with spatial_network(x, coords = ...).",
         call. = FALSE)
  }
  cx <- x$roles$x
  cy <- x$roles$y
  if (!all(c(cx, cy) %in% names(x$nodes))) {
    stop(sprintf("Coordinate columns `%s` and `%s` are no longer in `x$nodes`.", cx, cy),
         call. = FALSE)
  }
  check_coordinates(x$nodes[[cx]], x$nodes[[cy]], cx, cy, x$roles$crs, "x$nodes")
  coords <- cbind(as.double(x$nodes[[cx]]), as.double(x$nodes[[cy]]))
  colnames(coords) <- c(cx, cy)
  list(coords = coords, lonlat = identical(x$roles$crs, "lonlat"),
       radius = if (is.null(x$roles$radius)) 6371008.7714 else x$roles$radius)
}

# The engine's view of node_coordinates(): separate vectors and codes.
spatial_coords <- function(x) {
  p <- node_coordinates(x)
  list(x = p$coords[, 1L], y = p$coords[, 2L], crs = if (p$lonlat) 0L else 1L,
       radius = p$radius, names = colnames(p$coords))
}

#' @rdname spatial_network
#' @export
edge_lengths <- function(x) {
  p <- spatial_coords(x)
  ptr <- engine_graph(x)
  engine_call("edge_lengths", .Call(C_edge_lengths, ptr, p$x, p$y, p$crs, p$radius))
}

#' @rdname spatial_network
#' @param d Distances, such as the result of `edge_lengths()`.
#' @param scale The distance at which the weight has fallen to \eqn{e^{-1}}
#'   (exponential and gaussian) or \eqn{2^{-e}} (power).
#' @param form The decay function; see Details.
#' @param exponent The exponent \eqn{e} of the `"power"` form.
#' @export
distance_decay <- function(d, scale, form = c("exponential", "gaussian", "power"),
                           exponent = 2) {
  form <- match.arg(form)
  if (!is.numeric(d)) stop("`d` must be numeric.", call. = FALSE)
  if (!is.numeric(scale) || length(scale) != 1L || is.na(scale) || scale <= 0) {
    stop("`scale` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(exponent) || length(exponent) != 1L || is.na(exponent) || exponent <= 0) {
    stop("`exponent` must be a single positive number.", call. = FALSE)
  }
  if (any(d < 0, na.rm = TRUE)) stop("`d` must not be negative.", call. = FALSE)
  r <- d / scale
  switch(form,
    exponential = exp(-r),
    gaussian = exp(-r * r),
    power = (1 + r)^(-exponent)
  )
}

#' @rdname spatial_network
#' @param groups The group of each node: a vector aligned with the rows of
#'   `x$nodes`, or the name of a column of `x$nodes`.
#' @export
group_geometry <- function(x, groups) {
  p <- spatial_coords(x)
  grp <- node_groups(x, groups, "groups")
  k <- length(grp$labels)
  m <- engine_call("group_geometry",
                   .Call(C_group_geometry, p$x, p$y, grp$codes, p$crs, p$radius, k),
                   edge_rows = FALSE)
  out <- data.frame(group = grp$labels, size = as.integer(m[, 1L]), stringsAsFactors = FALSE)
  out[[p$names[[1L]]]] <- m[, 2L]
  out[[p$names[[2L]]]] <- m[, 3L]
  out$mean_distance <- m[, 4L]
  out$standard_distance <- m[, 5L]
  out
}

#' Proximity networks from locations
#'
#' Builds a network that joins located points to the points near them:
#' every pair within a distance (a fixed-radius or disk network), or each
#' point to its `k` nearest neighbours. Neighbours are found with a k-d
#' tree (Bentley 1975; Friedman, Bentley, and Finkel 1977), so networks of
#' 100,000 points take seconds, and every candidate's distance is then
#' computed exactly, so the tree changes the speed but never the result.
#'
#' @inheritParams spatial_network
#' @param nodes A data frame of points, or a spatial netframe (its nodes and
#'   coordinate record are used).
#' @param within Join every pair of distinct points at most this far apart
#'   (in the units of the distances). Give `within` or `k`, not both.
#' @param k Join each point to its `k` nearest other points.
#' @param knn For `k` nearest neighbours: `"directed"` keeps an arc from
#'   each point to each of its neighbours; `"either"` makes an undirected
#'   edge wherever at least one of the two points lists the other;
#'   `"mutual"` only where both do.
#' @param id Column of `nodes` holding node identifiers, as in
#'   [netframe()].
#' @param threads Number of threads; see [netframe()].
#'
#' @details
#' Within-distance networks are undirected, contain one edge per pair at
#' distance at most `within` (so coincident points are joined), and list
#' each edge with the point that comes first in `nodes` first. Nearest
#' neighbours are listed in order of increasing distance, ties broken by
#' position in `nodes`, so the result is reproducible; a point with fewer
#' than `k` other points gets them all. The edges carry a `distance`
#' column, which is the network's weight column (as distances), so
#' path-based measures use spatial lengths unless asked not to.
#'
#' @references
#' Bentley, J. L. (1975). Multidimensional binary search trees used for
#' associative searching. *Communications of the ACM*, 18(9), 509-517.
#' \doi{10.1145/361002.361007}
#'
#' Friedman, J. H., Bentley, J. L., & Finkel, R. A. (1977). An algorithm
#' for finding best matches in logarithmic expected time. *ACM Transactions
#' on Mathematical Software*, 3(3), 209-226. \doi{10.1145/355744.355745}
#'
#' @return A spatial netframe (see [spatial_network()]) with edge columns
#'   `from`, `to`, and `distance`.
#' @examples
#' set.seed(1)
#' pts <- data.frame(id = 1:50, x = runif(50), y = runif(50))
#' g <- proximity_network(pts, coords = c("x", "y"), crs = "planar", within = 0.2)
#' g
#' h <- proximity_network(pts, coords = c("x", "y"), crs = "planar", k = 3, knn = "mutual")
#' summary(edge_lengths(h))
#' @export
proximity_network <- function(nodes, coords = c("lon", "lat"), crs = c("lonlat", "planar"),
                              radius = 6371008.7714, within = NULL, k = NULL,
                              knn = c("directed", "either", "mutual"), id = 1L,
                              threads = 0L) {
  knn <- match.arg(knn)
  threads <- resolve_threads(threads)
  if (inherits(nodes, "netframe")) {
    if (missing(coords) && !is.null(nodes$roles$x) && !is.null(nodes$roles$y)) {
      coords <- c(nodes$roles$x, nodes$roles$y)
      if (missing(crs)) crs <- nodes$roles$crs
      if (missing(radius)) radius <- nodes$roles$radius
    }
    if (missing(id)) id <- nodes$roles$id
    nodes <- nodes$nodes
  }
  crs <- match.arg(crs)
  if (!is.data.frame(nodes)) {
    stop(sprintf("`nodes` must be a data frame or a netframe, not %s.", describe_class(nodes)),
         call. = FALSE)
  }
  if (length(coords) != 2L) {
    stop("`coords` must name two columns of `nodes`.", call. = FALSE)
  }
  cx <- resolve_column(nodes, coords[[1L]], "coords[1]", "nodes")
  cy <- resolve_column(nodes, coords[[2L]], "coords[2]", "nodes")
  id <- resolve_column(nodes, id, "id", "nodes")
  if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) || !is.finite(radius) ||
      radius <= 0) {
    stop("`radius` must be a single positive number.", call. = FALSE)
  }
  check_coordinates(nodes[[cx]], nodes[[cy]], cx, cy, crs, "nodes")
  if (is.null(within) == is.null(k)) {
    stop("Give exactly one of `within` (a distance) or `k` (a number of neighbours).",
         call. = FALSE)
  }
  if (!is.null(within) && (!is.numeric(within) || length(within) != 1L || is.na(within) ||
                           within < 0)) {
    stop("`within` must be a single non-negative distance.", call. = FALSE)
  }
  if (!is.null(k) && (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1 ||
                      k != floor(k))) {
    stop("`k` must be a positive whole number.", call. = FALSE)
  }
  r <- engine_call("proximity_network",
                   .Call(C_proximity, as.double(nodes[[cx]]), as.double(nodes[[cy]]),
                         if (crs == "lonlat") 0L else 1L, as.double(radius),
                         if (is.null(within)) 0 else as.double(within),
                         if (is.null(k)) 0 else as.double(k), threads),
                   edge_rows = FALSE)
  from <- r[[1L]] + 1
  to <- r[[2L]] + 1
  distance <- r[[3L]]
  directed <- FALSE
  if (!is.null(k)) {
    if (knn == "directed") {
      directed <- TRUE
    } else {
      a <- pmin(from, to)
      b <- pmax(from, to)
      key <- paste(a, b)
      if (knn == "either") {
        keep <- !duplicated(key)
      } else {
        keep <- duplicated(key)
      }
      from <- a[keep]
      to <- b[keep]
      distance <- distance[keep]
    }
  }
  ids <- nodes[[id]]
  edges <- data.frame(from = ids[from], to = ids[to], distance = distance,
                      stringsAsFactors = FALSE)
  g <- netframe(edges, nodes, directed = directed, id = id, weight = "distance",
                weight_type = "distance")
  spatial_network(g, coords = c(cx, cy), crs = crs, radius = radius)
}
