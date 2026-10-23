# Geospatial networks with coordinate columns: recording coordinates,
# edge lengths, distance decay, group geometry, and proximity networks.
# Networks whose tables are sf data frames are handled in sf.R, and the
# map plot in plot-spatial.R.

# Geospatial networks (man/net_spatial.Rd).
net_spatial <- function(x, coords = c("lon", "lat"), crs = c("lonlat", "planar"),
                        radius = 6371008.7714,
                        representative = c("point_on_surface", "centroid"), quiet = FALSE) {
  check_netframe(x)
  crs <- match.arg(crs)
  representative <- match.arg(representative)
  radius <- check_number(radius, "radius", lower = 0, open = TRUE)
  quiet <- check_flag(quiet, "quiet")
  if (has_node_geometry(x)) {
    if (!missing(coords)) {
      stop(sprintf(paste0("`x$nodes` already carries geometry in column `%s`, from which ",
                          "coordinates are taken; drop it with sf::st_drop_geometry() ",
                          "to use coordinate columns instead."), x$roles$node_geometry),
           call. = FALSE)
    }
    x$roles$representative <- representative
    x$roles$radius <- as.double(radius)
    return(x)
  }
  if (length(coords) != 2L) {
    stop("`coords` must name two columns of `x$nodes`: longitude and latitude, or x and y.",
         call. = FALSE)
  }
  cx <- resolve_column(x$nodes, coords[[1L]], "coords[1]", "nodes")
  cy <- resolve_column(x$nodes, coords[[2L]], "coords[2]", "nodes")
  roles <- x$roles
  roles$x <- cx
  roles$y <- cy
  roles$crs <- crs
  roles$radius <- as.double(radius)
  out <- new_netframe(x$edges, x$nodes, x$directed, roles, x$declared, x$graph)
  build_netframe(out, quiet, "net_spatial")
}

# See man/net_spatial.Rd.
net_is_spatial <- function(x) inherits(x, "spatial_netframe")

check_spatial <- function(x) {
  check_netframe(x)
  if (!inherits(x, "spatial_netframe")) {
    stop(paste0("`x` must be a spatial network (class `spatial_netframe`), with node ",
                "coordinates or sf geometry; create one with net_spatial(x, coords = ...) ",
                "or from sf tables."), call. = FALSE)
  }
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
  if (has_node_geometry(x)) {
    require_pkg("sf", "A netframe with sf node geometry")
    return(geometry_coordinates(x))
  }
  if (is.null(x$roles$x) || is.null(x$roles$y) || is.null(x$roles$crs)) {
    stop("`x` has no coordinates; record them with net_spatial(x, coords = ...).",
         call. = FALSE)
  }
  cx <- x$roles$x
  cy <- x$roles$y
  if (!all(c(cx, cy) %in% names(x$nodes))) {
    stop(sprintf("Coordinate columns `%s` and `%s` are no longer in `x$nodes`.", cx, cy),
         call. = FALSE)
  }
  kind <- crs_kind(x)
  check_coordinates(x$nodes[[cx]], x$nodes[[cy]], cx, cy, kind, "x$nodes")
  coords <- cbind(as.double(x$nodes[[cx]]), as.double(x$nodes[[cy]]))
  colnames(coords) <- c(cx, cy)
  list(coords = coords, lonlat = identical(kind, "lonlat"),
       radius = if (is.null(x$roles$radius)) 6371008.7714 else x$roles$radius)
}

# The engine's view of node_coordinates(): separate vectors and codes.
spatial_coords <- function(x) {
  p <- node_coordinates(x)
  list(x = p$coords[, 1L], y = p$coords[, 2L], crs = if (p$lonlat) 0L else 1L,
       radius = p$radius, names = colnames(p$coords))
}

# See man/net_spatial.Rd.
net_edge_lengths <- function(x) {
  check_netframe(x)
  if (has_edge_geometry(x)) {
    require_pkg("sf", "A netframe with sf edge geometry")
    return(geometry_edge_lengths(x))
  }
  p <- spatial_coords(x)
  ptr <- engine_graph(x)
  engine_call("net_edge_lengths", .Call(C_edge_lengths, ptr, p$x, p$y, p$crs, p$radius))
}

# See man/net_spatial.Rd.
net_distance_decay <- function(d, scale, form = c("exponential", "gaussian", "power"),
                           exponent = 2) {
  form <- match.arg(form)
  if (!is.numeric(d)) stop("`d` must be numeric.", call. = FALSE)
  scale <- check_number(scale, "scale", lower = 0, open = TRUE)
  exponent <- check_number(exponent, "exponent", lower = 0, open = TRUE)
  if (any(d < 0, na.rm = TRUE)) stop("`d` must not be negative.", call. = FALSE)
  r <- d / scale
  switch(form,
    exponential = exp(-r),
    gaussian = exp(-r * r),
    power = (1 + r)^(-exponent)
  )
}

# See man/net_spatial.Rd.
net_group_geometry <- function(x, groups) {
  p <- spatial_coords(x)
  grp <- node_groups(x, groups, "groups")
  k <- length(grp$labels)
  m <- engine_call("net_group_geometry",
                   .Call(C_group_geometry, p$x, p$y, grp$codes, p$crs, p$radius, k),
                   edge_rows = FALSE)
  out <- data.frame(group = grp$labels, size = as.integer(m[, 1L]), stringsAsFactors = FALSE)
  out[[p$names[[1L]]]] <- m[, 2L]
  out[[p$names[[2L]]]] <- m[, 3L]
  out$mean_distance <- m[, 4L]
  out$standard_distance <- m[, 5L]
  out
}

# Proximity networks from locations (man/net_proximity.Rd).
net_proximity <- function(nodes, coords = c("lon", "lat"), crs = c("lonlat", "planar"),
                              radius = 6371008.7714, within = NULL, k = NULL,
                              knn = c("directed", "either", "mutual"), id = 1L,
                              threads = 0L, exact = FALSE) {
  knn <- match.arg(knn)
  threads <- resolve_threads(threads)
  check_flag(exact, "exact")
  check_proximity_args(within, k)
  sf_nodes <- NULL
  if (inherits(nodes, "netframe")) {
    if (missing(coords) && has_node_geometry(nodes)) {
      sf_nodes <- nodes$nodes
    } else if (missing(coords) && !is.null(nodes$roles$x) && !is.null(nodes$roles$y)) {
      coords <- c(nodes$roles$x, nodes$roles$y)
      if (missing(crs)) crs <- nodes$roles$crs
      if (missing(radius)) radius <- nodes$roles$radius
    }
    if (missing(id)) id <- nodes$roles$id
    nodes <- nodes$nodes
  } else if (inherits(nodes, "sf") && missing(coords)) {
    sf_nodes <- nodes
    if (missing(id)) id <- NULL
  }
  if (!is.null(sf_nodes)) {
    require_pkg("sf", "net_proximity() with sf nodes")
    return(sf_proximity_network(sf_nodes, within, k, knn, id, exact, threads))
  }
  if (exact) {
    stop("`exact = TRUE` needs sf nodes; see ?proximity_network.", call. = FALSE)
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
  radius <- check_number(radius, "radius", lower = 0, open = TRUE)
  check_coordinates(nodes[[cx]], nodes[[cy]], cx, cy, crs, "nodes")
  r <- proximity_pairs(as.double(nodes[[cx]]), as.double(nodes[[cy]]), crs == "lonlat",
                       radius, within, k, threads)
  r <- symmetrize_knn(r, k, knn)
  ids <- nodes[[id]]
  edges <- data.frame(source = ids[r$from], target = ids[r$to], distance = r$distance,
                      stringsAsFactors = FALSE)
  g <- net_frame(edges, nodes, directed = r$directed, id = id, weight = "distance",
                 weight_type = "distance", quiet = TRUE)
  net_spatial(g, coords = c(cx, cy), crs = crs, radius = radius, quiet = TRUE)
}

check_proximity_args <- function(within, k) {
  if (is.null(within) == is.null(k)) {
    stop("Give exactly one of `within` (a distance) or `k` (a number of neighbours).",
         call. = FALSE)
  }
  if (!is.null(within)) check_number(within, "within", lower = 0)
  if (!is.null(k)) check_count(k, "k", min = 1L)
  invisible(NULL)
}

# The engine's neighbour search: pairs within a distance or each point's
# k nearest neighbours, as 1-based positions with their distances.
proximity_pairs <- function(px, py, lonlat, radius, within, k, threads) {
  r <- engine_call("net_proximity",
                   .Call(C_proximity, px, py, if (lonlat) 0L else 1L, as.double(radius),
                         if (is.null(within)) 0 else as.double(within),
                         if (is.null(k)) 0 else as.double(k), threads),
                   edge_rows = FALSE)
  list(from = r[[1L]] + 1L, to = r[[2L]] + 1L, distance = r[[3L]])
}

# Turns k-nearest-neighbour arcs into the edges `knn` asks for; pairs from
# `within` are already undirected.
symmetrize_knn <- function(r, k, knn) {
  r$directed <- FALSE
  if (is.null(k)) return(r)
  if (knn == "directed") {
    r$directed <- TRUE
    return(r)
  }
  a <- pmin(r$from, r$to)
  b <- pmax(r$from, r$to)
  key <- paste(a, b)
  keep <- if (knn == "either") !duplicated(key) else duplicated(key)
  list(from = a[keep], to = b[keep], distance = r$distance[keep], directed = FALSE)
}
