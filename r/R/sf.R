# Spatial networks whose node or edge tables are sf data frames
# (pebesma2018sf). All geometry work is delegated to sf; the engine receives
# derived numbers only. Polygon contiguity uses the DE-9IM relations of
# egenhofer1991pointset and clementini1993smallset with the queen and rook
# definitions of anselin1988spatial, matching spdep::poly2nb
# (bivand2013applied, bivand2018comparing). Line snapping at tolerance 0
# matches sfnetworks (vandermeer2026sfnetworks); the two-table layout of
# nodes and edges follows OSMnx's graph_to_gdfs (boeing2017osmnx) and
# sfnetworks. sf, sfnetworks, and spdep are suggested packages: every use is
# guarded by requireNamespace() and the package works without them.

#' Spatial networks with sf
#'
#' A spatial network is a [netframe()] whose node table, edge table, or
#' both are simple-feature data frames from the sf package (Pebesma 2018).
#' Nodes and edges may carry any geometry type (POINT, LINESTRING,
#' POLYGON, their MULTI forms, or GEOMETRYCOLLECTION) and a coordinate
#' reference system; edges may also have no geometry at all, or geometry
#' for only some rows (empty geometries for the rest). Connectivity always
#' comes from the `from` and `to` columns, never from the geometry, so a
#' line that happens to cross a node does not join it. Geometry is handled
#' entirely by sf (GEOS and s2); the engine only ever sees numbers derived
#' from it, such as edge lengths and representative points.
#'
#' [netframe()] and [as_netframe()] notice sf tables on their own and
#' record the geometry columns and CRS in `x$roles` (`node_geometry`,
#' `edge_geometry`, and `crs`, which then holds an sf `crs` object). Nodes
#' and edges that both have geometry must share a CRS; use
#' [sf::st_transform()] or [sf::st_set_crs()] on one table first if they
#' do not. The functions here read those roles:
#'
#' * `to_sf()` returns the node or edge table as an sf data frame with all
#'   its attributes. Nodes recorded with coordinate columns by
#'   [spatial_network()] become POINT features (WGS 84 for longitude and
#'   latitude, no CRS for planar coordinates). Edges without geometry
#'   become straight LINESTRINGs between the representative points of
#'   their endpoints when `implicit = TRUE` (sfnetworks calls these
#'   spatially implicit edges); with `implicit = FALSE` the edge table is
#'   returned as it is, and an error is raised if it has no geometry.
#'   [sf::st_as_sf()] on a netframe does the same.
#' * `check_geometry()` is the optional consistency check between edge
#'   geometries and node geometries. It is never required: a network is
#'   valid however its geometries lie. For each edge with geometry it
#'   reports whether the geometry starts within `tolerance` of the `from`
#'   node's geometry and ends within `tolerance` of the `to` node's
#'   geometry (in either orientation for undirected networks), `NA` for
#'   edges without geometry. Distances are [sf::st_distance()], so a line
#'   ending anywhere inside a polygon node is at distance 0 from it.
#' * `transform_network()` applies [sf::st_transform()] to both tables at
#'   once and updates the recorded CRS. Columns computed earlier from the
#'   geometry, such as lengths, are not recomputed.
#' * `spatial_filter()` keeps the nodes whose geometry satisfies
#'   `predicate` with any feature of `y` (any sf or sfc object, a bbox, or
#'   another spatial netframe) and returns the induced subnetwork: those
#'   nodes and every edge between two of them, with all columns and roles,
#'   as [ego_network()] does.
#' * `spatial_join()` attaches the attributes of the feature of `y` that
#'   each node's geometry matches under `join` (an [sf::st_join()] with
#'   `left = TRUE`). The node table keeps one row per node: a node that
#'   matches several features keeps the first and a warning says how many
#'   nodes had more; a node that matches none gets `NA`. Columns of `y`
#'   whose names clash with node columns get the suffix `_y`.
#'
#' @section Representative points:
#' Every point-based function of [spatial_network()] (`edge_lengths()` for
#' edges without geometry, [group_geometry()], [proximity_network()], the
#' map drawn by [plot_spatial()]) reduces each node geometry to one
#' point. POINT geometries are used as they are. Other geometries are
#' reduced by [sf::st_point_on_surface()] by default, which always lies on
#' the geometry (a crescent's centroid does not), or by
#' [sf::st_centroid()] when `spatial_network(x, representative =
#' "centroid")` has been called. For longitude and latitude, point on
#' surface is computed on the planar coordinates (sf warns of this; the
#' warning is silenced here because the point still lies on the geometry),
#' while centroids are computed on the sphere by s2. A node whose geometry
#' is empty has no representative point and those functions stop with an
#' error naming the row.
#'
#' @section Edge lengths:
#' When the edge table has geometry, [edge_lengths()] returns
#' [sf::st_length()] of each edge geometry: a geodesic length for
#' longitude and latitude (computed on the sphere by s2 when
#' [sf::sf_use_s2()] is on, otherwise on the ellipsoid), and a planar
#' length in the units of the CRS otherwise. Edges with empty geometry get
#' the distance between the representative points of their endpoints
#' ([sf::st_distance()], in the same units). The result is a plain numeric
#' vector with the unit recorded in `attr(, "unit")` (for example `"m"`),
#' empty when the CRS is unknown. Without edge geometry, [edge_lengths()]
#' uses the representative points and the haversine formula on a sphere
#' of radius `x$roles$radius`, as documented in [spatial_network()]; the
#' two differ by the small difference between s2's Earth radius and the
#' GRS 80 mean radius, and by the ellipsoid when s2 is off.
#'
#' @param x A [netframe()] with sf nodes or edges; for `to_sf(what =
#'   "nodes")` also one with coordinate columns from [spatial_network()].
#' @param what Which table to return: `"nodes"` or `"edges"`.
#' @param implicit For `what = "edges"`, whether edges without geometry get
#'   straight-line geometries between their endpoints.
#' @param tolerance Largest distance, in the units of the CRS (metres for
#'   longitude and latitude), at which an edge geometry's end still counts
#'   as touching its node.
#' @param crs Target coordinate reference system, as accepted by
#'   [sf::st_transform()] (an EPSG code, a string, or an sf `crs`).
#' @param y Features to filter or join by: an sf or sfc object, an sf
#'   `bbox`, or a spatial netframe (its node geometries are used).
#' @param predicate,join A geometric binary predicate from sf, such as
#'   [sf::st_intersects()] (the default), [sf::st_within()], or
#'   [sf::st_is_within_distance()], called as `predicate(nodes, y, ...)`.
#' @param ... Further arguments passed to `predicate` or `join`, such as
#'   `dist` for [sf::st_is_within_distance()].
#'
#' @return `to_sf()` returns an sf data frame. `check_geometry()` returns a
#'   logical vector aligned with the rows of `x$edges`. The other functions
#'   return a netframe: `transform_network()` and `spatial_join()` with the
#'   same nodes and edges, `spatial_filter()` the induced subnetwork.
#'
#' @references
#' Pebesma, E. (2018). Simple features for R: Standardized support for
#' spatial vector data. *The R Journal*, 10(1), 439-446.
#' \doi{10.32614/RJ-2018-009}
#'
#' Barthelemy, M. (2011). Spatial networks. *Physics Reports*, 499(1-3),
#' 1-101. \doi{10.1016/j.physrep.2010.11.002}
#'
#' @seealso [from_lines()] and [from_polygons()] to build spatial networks
#'   from sf data, [to_sfnetwork()] for sfnetworks, [read_network()] for
#'   GeoPackage files, [spatial_network()] and [proximity_network()] for
#'   the point-based functions.
#' @examples
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   pts <- sf::st_as_sf(data.frame(name = c("a", "b", "c"), x = c(0, 1, 1), y = c(0, 0, 1)),
#'                       coords = c("x", "y"), crs = 32633)  # UTM zone 33N, metres
#'   g <- netframe(data.frame(from = c("a", "b", "c"), to = c("b", "c", "a")), pts,
#'                 directed = FALSE)
#'   g
#'   edge_lengths(g)
#'   to_sf(g, "edges")
#'   h <- netframe(to_sf(g, "edges"), pts, directed = FALSE)  # explicit edge geometry
#'   check_geometry(h)
#'   transform_network(g, 4326)$roles$crs$epsg
#'   box <- sf::st_bbox(c(xmin = -0.5, ymin = -0.5, xmax = 0.5, ymax = 0.5), crs = 32633)
#'   spatial_filter(g, box)
#' }
#' @export
to_sf <- function(x, what = c("nodes", "edges"), implicit = TRUE) {
  check_netframe(x)
  what <- match.arg(what)
  check_flag(implicit, "implicit")
  require_sf("to_sf()")
  if (what == "nodes") return(node_table_sf(x))
  edge_table_sf(x, implicit)
}

# The node table as sf: as it is when it has geometry, or POINT features
# built from the coordinate columns recorded by spatial_network().
node_table_sf <- function(x) {
  if (has_node_geometry(x)) {
    nodes <- x$nodes
    if (!inherits(nodes, "sf")) nodes <- sf::st_sf(nodes, sf_column_name = x$roles$node_geometry)
    return(nodes)
  }
  if (is.null(x$roles$x) || is.null(x$roles$y)) {
    stop(paste0("`x$nodes` has no geometry. Give sf nodes to netframe(), or record ",
                "coordinate columns with spatial_network(x, coords = ...)."), call. = FALSE)
  }
  p <- node_coordinates(x)
  crs <- if (p$lonlat) sf::st_crs(4326) else sf::st_crs(NA)
  sf::st_as_sf(x$nodes, coords = colnames(p$coords), crs = crs, remove = FALSE)
}

# The edge table as sf. Without geometry, or for rows whose geometry is
# empty, straight lines between the endpoints' representative points.
edge_table_sf <- function(x, implicit) {
  edges <- x$edges
  if (has_edge_geometry(x)) {
    if (!inherits(edges, "sf")) edges <- sf::st_sf(edges, sf_column_name = x$roles$edge_geometry)
    if (!implicit) return(edges)
    geom <- sf::st_geometry(edges)
    empty <- which(sf::st_is_empty(geom))
    if (length(empty)) {
      geom[empty] <- straight_lines(x, empty)
      sf::st_geometry(edges) <- geom
    }
    return(edges)
  }
  if (!implicit) {
    stop(paste0("`x$edges` has no geometry; use `implicit = TRUE` for straight lines ",
                "between the endpoints."), call. = FALSE)
  }
  geom <- straight_lines(x, seq_len(nrow(edges)))
  col <- if ("geometry" %in% names(edges)) "edge_geometry" else "geometry"
  edges[[col]] <- geom
  sf::st_sf(edges, sf_column_name = col)
}

# Straight LINESTRINGs between the representative points of the endpoints
# of the edges at 1-based positions `rows`, in the network's CRS.
straight_lines <- function(x, rows) {
  pts <- node_points_any(x)
  xy <- point_matrix(pts, "`x$nodes`")
  src <- x$index$src[rows] + 1L
  tgt <- x$index$tgt[rows] + 1L
  lines <- lapply(seq_along(rows), function(k) {
    sf::st_linestring(rbind(xy[src[k], ], xy[tgt[k], ]))
  })
  sf::st_sfc(lines, crs = sf::st_crs(pts))
}

#' @rdname to_sf
#' @export
check_geometry <- function(x, tolerance = 0) {
  check_netframe(x)
  require_sf("check_geometry()")
  check_tolerance(tolerance)
  if (!has_edge_geometry(x)) {
    return(rep(NA, nrow(x$edges)))
  }
  if (!has_node_geometry(x)) {
    stop("`x$nodes` has no geometry to check the edges against.", call. = FALSE)
  }
  geom <- sf::st_geometry(x$edges)
  nodes <- sf::st_geometry(x$nodes)
  ends <- line_ends(geom)
  out <- rep(NA, length(geom))
  has <- which(!sf::st_is_empty(geom))
  if (!length(has)) return(out)
  src <- x$index$src[has] + 1L
  tgt <- x$index$tgt[has] + 1L
  near <- function(a, b) {
    strip_units(sf::st_distance(a, b, by_element = TRUE)) <= tolerance
  }
  forward <- near(ends$start[has], nodes[src]) & near(ends$end[has], nodes[tgt])
  if (!x$directed) {
    forward <- forward | (near(ends$start[has], nodes[tgt]) & near(ends$end[has], nodes[src]))
  }
  out[has] <- forward
  out
}

check_tolerance <- function(tolerance) {
  if (!is.numeric(tolerance) || length(tolerance) != 1L || is.na(tolerance) || tolerance < 0) {
    stop("`tolerance` must be a single non-negative distance.", call. = FALSE)
  }
  invisible(tolerance)
}

#' @rdname to_sf
#' @export
transform_network <- function(x, crs) {
  check_netframe(x)
  require_sf("transform_network()")
  if (!has_node_geometry(x) && !has_edge_geometry(x)) {
    stop("`x` has no geometry to transform; see ?to_sf for spatial networks.", call. = FALSE)
  }
  if (has_node_geometry(x)) x$nodes <- sf::st_transform(node_table_sf(x), crs)
  if (has_edge_geometry(x)) x$edges <- sf::st_transform(edge_table_sf(x, implicit = FALSE), crs)
  x$roles$crs <- sf::st_crs(if (has_node_geometry(x)) x$nodes else x$edges)
  x
}

#' @rdname to_sf
#' @export
spatial_filter <- function(x, y, predicate = sf::st_intersects, ...) {
  check_netframe(x)
  require_sf("spatial_filter()")
  nodes <- need_node_geometry(x, "spatial_filter()")
  target <- filter_geometry(y)
  if (!is.function(predicate)) {
    stop("`predicate` must be a function such as sf::st_intersects.", call. = FALSE)
  }
  hits <- predicate(nodes, target, ...)
  if (length(hits) != length(nodes)) {
    stop("`predicate` must return one element per node, as sf's binary predicates do.",
         call. = FALSE)
  }
  keep <- which(lengths(hits) > 0L)
  induced_netframe(x, keep)
}

filter_geometry <- function(y) {
  if (inherits(y, "netframe")) {
    return(need_node_geometry(y, "spatial_filter()"))
  }
  if (inherits(y, "bbox")) return(sf::st_as_sfc(y))
  if (inherits(y, "sf") || inherits(y, "sfc")) return(sf::st_geometry(y))
  stop(sprintf("`y` must be an sf or sfc object, a bbox, or a spatial netframe, not %s.",
               describe_class(y)), call. = FALSE)
}

#' @rdname to_sf
#' @export
spatial_join <- function(x, y, join = sf::st_intersects, ...) {
  check_netframe(x)
  require_sf("spatial_join()")
  need_node_geometry(x, "spatial_join()")
  if (!inherits(y, "sf")) {
    stop(sprintf("`y` must be an sf data frame, not %s.", describe_class(y)), call. = FALSE)
  }
  if (!is.function(join)) {
    stop("`join` must be a function such as sf::st_intersects.", call. = FALSE)
  }
  nodes <- node_table_sf(x)
  key <- ".knappnet_row"
  nodes[[key]] <- seq_len(nrow(nodes))
  joined <- sf::st_join(nodes, y, join = join, left = TRUE, suffix = c("", "_y"), ...)
  rows <- joined[[key]]
  dup <- duplicated(rows)
  if (any(dup)) {
    warning(sprintf(paste0("%s matched more than one feature of `y`; the first match ",
                           "was kept for each"), plural(length(unique(rows[dup])), "node")),
            call. = FALSE)
    joined <- joined[!dup, , drop = FALSE]
    rows <- rows[!dup]
  }
  joined <- joined[order(rows), , drop = FALSE]
  joined[[key]] <- NULL
  rownames(joined) <- NULL
  if (nrow(joined) != nrow(x$nodes)) {
    stop("spatial_join(): the join did not return one row per node.", call. = FALSE)
  }
  x$nodes <- joined
  x
}

# --- Geometry roles -----------------------------------------------------

# Records the geometry columns and CRS of sf node and edge tables in the
# roles of a new netframe. Called by netframe(); a netframe built from
# plain data frames passes through untouched, so nothing here needs sf.
record_geometry <- function(x) {
  ng <- if (inherits(x$nodes, "sf")) attr(x$nodes, "sf_column")
  eg <- if (inherits(x$edges, "sf")) attr(x$edges, "sf_column")
  if (is.null(ng) && is.null(eg)) return(x)
  require_sf("A netframe with sf nodes or edges")
  if (!is.null(ng) && !inherits(x$nodes[[ng]], "sfc")) {
    stop(sprintf("The geometry column `%s` of `nodes` is not an sfc list column.", ng),
         call. = FALSE)
  }
  if (!is.null(eg) && !inherits(x$edges[[eg]], "sfc")) {
    stop(sprintf("The geometry column `%s` of `edges` is not an sfc list column.", eg),
         call. = FALSE)
  }
  crs_n <- if (!is.null(ng)) sf::st_crs(x$nodes)
  crs_e <- if (!is.null(eg)) sf::st_crs(x$edges)
  if (!is.null(crs_n) && !is.null(crs_e) && !same_crs(crs_n, crs_e)) {
    stop(paste0("`nodes` and `edges` have different coordinate reference systems (",
                crs_label(crs_n), " and ", crs_label(crs_e), "). Transform one table ",
                "with sf::st_transform(), or set a missing CRS with sf::st_set_crs(), ",
                "before building the network."), call. = FALSE)
  }
  x$roles$node_geometry <- ng
  x$roles$edge_geometry <- eg
  x$roles$crs <- if (!is.null(crs_n)) crs_n else crs_e
  # Coordinate columns recorded earlier would now be ignored.
  x$roles$x <- NULL
  x$roles$y <- NULL
  x
}

same_crs <- function(a, b) {
  if (is.na(a) && is.na(b)) return(TRUE)
  if (is.na(a) || is.na(b)) return(FALSE)
  isTRUE(a == b)
}

crs_label <- function(crs) {
  if (is.na(crs)) return("no CRS")
  name <- tryCatch(crs$Name, error = function(e) NULL)
  epsg <- tryCatch(crs$epsg, error = function(e) NA)
  if (!is.null(name) && !is.na(name) && nzchar(name)) {
    if (!is.na(epsg)) sprintf("%s (EPSG:%d)", name, as.integer(epsg)) else name
  } else if (!is.na(epsg)) {
    sprintf("EPSG:%d", as.integer(epsg))
  } else {
    "a user-defined CRS"
  }
}

has_node_geometry <- function(x) !is.null(x$roles$node_geometry)
has_edge_geometry <- function(x) !is.null(x$roles$edge_geometry)

# Whether node positions are available for layouts: sf node geometry or
# coordinate columns recorded by spatial_network().
has_positions <- function(x) {
  has_node_geometry(x) || (!is.null(x$roles$x) && !is.null(x$roles$y))
}

require_sf <- function(what) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop(sprintf("%s needs the sf package; install it with install.packages(\"sf\").", what),
         call. = FALSE)
  }
}

need_node_geometry <- function(x, fn) {
  if (!has_node_geometry(x)) {
    stop(sprintf(paste0("%s needs node geometries: build the network from an sf node ",
                        "table, or convert coordinate columns with to_sf(x, \"nodes\")."),
                 fn), call. = FALSE)
  }
  sf::st_geometry(x$nodes)
}

# "lonlat" or "planar" for any recorded CRS: the strings spatial_network()
# stores, or an sf crs (longitude/latitude when sf says so; an unknown CRS
# counts as planar).
crs_kind <- function(x) {
  crs <- x$roles$crs
  if (is.character(crs)) return(crs)
  if (inherits(crs, "crs")) {
    require_sf("A netframe with sf geometry")
    return(if (isTRUE(sf::st_is_longlat(crs))) "lonlat" else "planar")
  }
  "planar"
}

# Representative points for any spatial netframe: from sf node geometry,
# or POINT features built from the coordinate columns.
node_points_any <- function(x) {
  if (has_node_geometry(x)) node_points(x) else sf::st_geometry(node_table_sf(x))
}

# One representative POINT per node: POINT geometries as they are, others
# reduced by point on surface or centroid (see ?to_sf).
node_points <- function(x) {
  geom <- need_node_geometry(x, "This function")
  method <- x$roles$representative
  if (is.null(method)) method <- "point_on_surface"
  is_point <- as.character(sf::st_geometry_type(geom, by_geometry = TRUE)) == "POINT"
  if (all(is_point)) return(geom)
  empty <- sf::st_is_empty(geom)
  out <- geom
  idx <- which(!is_point & !empty)
  if (length(idx)) {
    out[idx] <- if (method == "centroid") {
      sf::st_centroid(geom[idx])
    } else {
      quiet_longlat(sf::st_point_on_surface(geom[idx]))
    }
  }
  # Empty geometries of other types become empty points, which
  # point_matrix() reports by row.
  none <- which(!is_point & empty)
  if (length(none)) out[none] <- sf::st_point()
  sf::st_cast(out, "POINT")
}

# An n-by-2 matrix of coordinates from POINT geometries, stopping on empty
# points (which have no coordinates) so rows cannot misalign.
point_matrix <- function(pts, table) {
  empty <- which(sf::st_is_empty(pts))
  if (length(empty)) {
    stop(sprintf("%s has empty geometries in %s; they have no location.", table,
                 rows_text(empty)), call. = FALSE)
  }
  if (!length(pts)) return(matrix(numeric(), 0L, 2L))
  co <- sf::st_coordinates(pts)
  co[, 1:2, drop = FALSE]
}

# node_coordinates() for networks with node geometry: representative
# points as a matrix, plus the longitude/latitude flag and sphere radius.
geometry_coordinates <- function(x) {
  pts <- node_points(x)
  lonlat <- crs_kind(x) == "lonlat"
  coords <- point_matrix(pts, "`x$nodes`")
  colnames(coords) <- if (lonlat) c("lon", "lat") else c("x", "y")
  check_coordinates(coords[, 1L], coords[, 2L], colnames(coords)[1L], colnames(coords)[2L],
                    if (lonlat) "lonlat" else "planar", "x$nodes")
  list(coords = coords, lonlat = lonlat,
       radius = if (is.null(x$roles$radius)) 6371008.7714 else x$roles$radius)
}

# edge_lengths() for networks with edge geometry; see ?to_sf.
geometry_edge_lengths <- function(x) {
  geom <- sf::st_geometry(x$edges)
  len <- sf::st_length(geom)
  unit <- unit_string(len)
  len <- strip_units(len)
  empty <- which(sf::st_is_empty(geom))
  if (length(empty)) {
    if (!has_positions(x)) {
      stop(sprintf(paste0("edge_lengths(): %s of `x$edges` %s empty geometry and `x$nodes` ",
                          "has no locations to measure between."),
                   rows_text(empty), if (length(empty) == 1L) "has" else "have"), call. = FALSE)
    }
    pts <- node_points_any(x)
    point_matrix(pts, "`x$nodes`")
    d <- sf::st_distance(pts[x$index$src[empty] + 1L], pts[x$index$tgt[empty] + 1L],
                         by_element = TRUE)
    len[empty] <- strip_units(d)
  }
  structure(as.double(len), unit = unit)
}

# The unit of a `units` vector as text ("m", "US_survey_foot", "m/s"), or
# "" for a plain number. Reads the symbolic units directly, so the units
# package is not needed beyond what sf already loads.
unit_string <- function(v) {
  u <- attr(v, "units")
  if (is.null(u)) return("")
  num <- paste(u$numerator, collapse = "*")
  den <- paste(u$denominator, collapse = "*")
  if (!nzchar(den)) return(num)
  if (!nzchar(num)) num <- "1"
  paste0(num, "/", den)
}

strip_units <- function(v) {
  if (inherits(v, "units")) {
    v <- unclass(v)
    attr(v, "units") <- NULL
  }
  as.double(v)
}

# First and last vertex of each (MULTI)LINESTRING as POINT geometries,
# aligned with `geom`; empty geometries give empty points.
line_ends <- function(geom) {
  m <- length(geom)
  crs <- sf::st_crs(geom)
  start <- matrix(NA_real_, m, 2L)
  end <- start
  if (m) {
    co <- sf::st_coordinates(geom)
    if (nrow(co)) {
      fid <- co[, ncol(co)]
      first <- !duplicated(fid)
      last <- !duplicated(fid, fromLast = TRUE)
      start[fid[first], ] <- co[first, 1:2]
      end[fid[last], ] <- co[last, 1:2]
    }
  }
  list(start = points_from_matrix(start, crs), end = points_from_matrix(end, crs))
}

points_from_matrix <- function(xy, crs) {
  df <- data.frame(x = xy[, 1L], y = xy[, 2L])
  sf::st_geometry(sf::st_as_sf(df, coords = c("x", "y"), crs = crs, na.fail = FALSE))
}

# One line for print.netframe().
spatial_summary <- function(x) {
  if (!has_node_geometry(x) && !has_edge_geometry(x)) {
    if (!is.null(x$roles$x) && !is.null(x$roles$y)) {
      return(sprintf("spatial: coordinates in `%s`, `%s` (%s)\n", x$roles$x, x$roles$y,
                     if (identical(x$roles$crs, "lonlat")) "longitude and latitude" else "planar"))
    }
    return("")
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    return("spatial: sf geometry (install sf to use it)\n")
  }
  describe <- function(geom) {
    types <- as.character(sf::st_geometry_type(geom, by_geometry = TRUE))
    empty <- sf::st_is_empty(geom)
    kinds <- unique(types[!empty])
    label <- if (!length(kinds)) "no geometry" else if (length(kinds) <= 2L) {
      paste(kinds, collapse = " and ")
    } else {
      "mixed geometry"
    }
    if (any(empty) && length(kinds)) {
      label <- sprintf("%s, %d without geometry", label, sum(empty))
    }
    label
  }
  parts <- c(
    if (has_node_geometry(x)) sprintf("nodes %s", describe(sf::st_geometry(x$nodes))),
    if (has_edge_geometry(x)) sprintf("edges %s", describe(sf::st_geometry(x$edges)))
  )
  sprintf("spatial: %s; CRS %s\n", paste(parts, collapse = ", "), crs_label(x$roles$crs))
}

# --- Builders -------------------------------------------------------------

#' Spatial networks from lines and polygons
#'
#' `from_lines()` builds a network from LINESTRING features such as roads,
#' rivers, or pipelines: each line becomes an edge, and its two endpoints
#' become nodes, shared with every other line that ends at the same place.
#' `from_polygons()` builds the contiguity network of regions: each
#' polygon becomes a node, and two polygons are joined when they share a
#' boundary. Both are standard constructions in spatial analysis (Barthelemy
#' 2011; Bivand, Pebesma, and Gomez-Rubio 2013).
#'
#' @section Snapping rule of `from_lines()`:
#' The first and last vertex of every line are collected, in the order of
#' the lines (start, then end). With `tolerance = 0`, two endpoints are the
#' same node exactly when their coordinates are equal, as in
#' `sfnetworks::as_sfnetwork()`. With a positive `tolerance`, endpoints at
#' most `tolerance` apart ([sf::st_is_within_distance()], in the units of
#' the CRS: metres for longitude and latitude) are the same node, and so,
#' transitively, are endpoints linked through others (single linkage, the
#' connected components of the within-tolerance relation). A tolerance
#' that is small next to the shortest line keeps chains from merging
#' distinct junctions. Each node is placed at the first endpoint that
#' formed it and numbered by first appearance; line geometries are kept as
#' they are, so with a positive tolerance an edge can end up to
#' `tolerance` from its node, which `check_geometry()` with the same
#' tolerance accepts.
#'
#' Nodes get an integer `id` and POINT geometry. Edges keep every column
#' of `lines` (columns called `from`, `to`, or `length` are renamed
#' `from_attr`, `to_attr`, and `length_attr`), gain `from` and `to`, and a
#' `length` column from [sf::st_length()] (see the Edge lengths section
#' of [to_sf()]) that is the network's weight column, as distances, so
#' path-based measures use travel distance along the lines. A line that
#' starts and ends at the same node (a ring road) is a self-loop and is
#' kept. The network is undirected unless `directed = TRUE`, in which case
#' each edge runs from a line's first vertex to its last.
#'
#' @section Contiguity rule of `from_polygons()`:
#' Two polygons are neighbours under `"queen"` contiguity when their
#' boundaries share at least one point, and under `"rook"` contiguity
#' when they share a boundary line, not just corners, as a rook and a
#' queen move on a chessboard (Anselin 1988). The relations are evaluated
#' exactly by [sf::st_relate()] with the DE-9IM patterns `F***T****`
#' (queen) and `F***1****` (rook) of Egenhofer and Franzosa (1991) and
#' Clementini, Di Felice, and van Oosterom (1993): interiors disjoint,
#' boundaries intersecting in anything (queen) or in a line (rook). That
#' is the definition of `spdep::poly2nb()` with `snap = 0`; polygons whose
#' shared boundaries are drawn with slightly different coordinates do not
#' touch, and should be snapped first ([sf::st_snap()], or rounding with
#' [sf::st_set_precision()]). Overlapping polygons, whose interiors
#' intersect, are not neighbours under either rule. Coordinates are
#' treated as planar, also for longitude and latitude (as sf and spdep
#' do). The network is undirected, each pair listed once with the polygon
#' that comes first in `polygons` first, and a polygon with no neighbours
#' is an isolate. The edges have no geometry; `edge_lengths()` then
#' measures between representative points (see [to_sf()]).
#'
#' @param lines An sf data frame or sfc of LINESTRING geometries (cast
#'   MULTILINESTRINGs with [sf::st_cast()] first).
#' @param tolerance Distance within which line endpoints are the same
#'   node; see the snapping rule.
#' @param directed Whether edges run from each line's first vertex to its
#'   last (`TRUE`) or are undirected (`FALSE`, the default).
#' @param polygons An sf data frame or sfc of POLYGON or MULTIPOLYGON
#'   geometries.
#' @param contiguity `"queen"` (shared boundary point) or `"rook"` (shared
#'   boundary line).
#' @param id Column of `polygons` holding region identifiers, by name or
#'   position, or `NULL` to number the regions 1, 2, ... in a new first
#'   column `id` (an existing `id` column is renamed `id_attr`).
#'
#' @return A netframe whose node table is an sf data frame; see
#'   [to_sf()].
#'
#' @references
#' Anselin, L. (1988). *Spatial Econometrics: Methods and Models*. Kluwer.
#' \doi{10.1007/978-94-015-7799-1}
#'
#' Egenhofer, M. J., & Franzosa, R. D. (1991). Point-set topological
#' spatial relations. *International Journal of Geographical Information
#' Systems*, 5(2), 161-174. \doi{10.1080/02693799108927841}
#'
#' Clementini, E., Di Felice, P., & van Oosterom, P. (1993). A small set
#' of formal topological relationships suitable for end-user interaction.
#' In *Advances in Spatial Databases* (pp. 277-295). Springer.
#' \doi{10.1007/3-540-56869-7_16}
#'
#' Bivand, R. S., Pebesma, E., & Gomez-Rubio, V. (2013). *Applied Spatial
#' Data Analysis with R* (2nd ed.). Springer.
#' \doi{10.1007/978-1-4614-7618-4}
#'
#' Barthelemy, M. (2011). Spatial networks. *Physics Reports*, 499(1-3),
#' 1-101. \doi{10.1016/j.physrep.2010.11.002}
#'
#' @seealso [to_sf()], [proximity_network()], [to_sfnetwork()].
#' @examples
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   # A small road grid: two streets crossing two avenues.
#'   seg <- function(x0, y0, x1, y1) sf::st_linestring(rbind(c(x0, y0), c(x1, y1)))
#'   roads <- sf::st_sf(name = c("1st", "2nd", "A", "B"),
#'                      geometry = sf::st_sfc(seg(0, 0, 2, 0), seg(0, 1, 2, 1),
#'                                            seg(0, 0, 0, 1), seg(2, 0, 2, 1)))
#'   g <- from_lines(roads)
#'   g$nodes
#'   g$edges
#'   edge_lengths(g)
#'
#'   # Four square regions in a row, each touching the next.
#'   sq <- function(x) sf::st_polygon(list(rbind(c(x, 0), c(x + 1, 0), c(x + 1, 1),
#'                                               c(x, 1), c(x, 0))))
#'   regions <- sf::st_sf(region = c("w", "x", "y", "z"),
#'                        geometry = sf::st_sfc(sq(0), sq(1), sq(2), sq(3)))
#'   from_polygons(regions, id = "region")$edges
#' }
#' @export
from_lines <- function(lines, tolerance = 0, directed = FALSE) {
  require_sf("from_lines()")
  check_tolerance(tolerance)
  check_flag(directed, "directed")
  lines <- as_sf_table(lines, "lines")
  geom <- sf::st_geometry(lines)
  types <- as.character(sf::st_geometry_type(geom, by_geometry = TRUE))
  bad <- which(types != "LINESTRING")
  if (length(bad)) {
    stop(sprintf(paste0("`lines` must hold LINESTRING geometries; %s %s %s. Cast with ",
                        "sf::st_cast(lines, \"LINESTRING\")."), rows_text(bad),
                 if (length(bad) == 1L) "is" else "are",
                 paste(unique(types[bad]), collapse = ", ")), call. = FALSE)
  }
  empty <- which(sf::st_is_empty(geom))
  if (length(empty)) {
    stop(sprintf("`lines` has empty geometries in %s.", rows_text(empty)), call. = FALSE)
  }
  m <- length(geom)
  ends <- line_ends(geom)
  # Endpoints in line order: start of line 1, end of line 1, start of line 2, ...
  pts <- c(ends$start, ends$end)[as.vector(rbind(seq_len(m), m + seq_len(m)))]
  node_of <- snap_points(pts, tolerance)
  node_pts <- pts[!duplicated(node_of)]
  n <- length(node_pts)
  nodes <- sf::st_sf(id = seq_len(n), geometry = node_pts)
  from <- node_of[seq(1L, by = 2L, length.out = m)]
  to <- node_of[seq(2L, by = 2L, length.out = m)]
  attrs <- sf::st_drop_geometry(lines)
  names(attrs) <- rename_reserved(names(attrs), c("from", "to", "length"))
  len <- strip_units(sf::st_length(geom))
  edges <- cbind(data.frame(from = from, to = to), attrs, data.frame(length = len))
  edges <- sf::st_sf(edges, geometry = geom)
  netframe(edges, nodes, directed = directed, from = "from", to = "to", id = "id",
           weight = "length", weight_type = "distance")
}

# Groups POINT geometries into nodes: equal coordinates at tolerance 0,
# single linkage within `tolerance` otherwise. Returns group numbers in
# order of first appearance.
snap_points <- function(pts, tolerance) {
  k <- length(pts)
  if (k == 0L) return(integer())
  if (tolerance == 0) {
    co <- sf::st_coordinates(pts)
    key <- paste(sprintf("%.17g", co[, 1L]), sprintf("%.17g", co[, 2L]))
    return(match(key, unique(key)))
  }
  near <- sf::st_is_within_distance(pts, pts, dist = tolerance)
  i <- rep(seq_len(k), lengths(near))
  j <- unlist(near, use.names = FALSE)
  keep <- i < j
  pairs <- data.frame(from = i[keep], to = j[keep])
  g <- netframe(pairs, nodes = data.frame(id = seq_len(k)), directed = FALSE)
  member <- components(g, "weak")$membership
  match(member, unique(member))
}

# Runs sf code whose notice that longitude/latitude coordinates are
# treated as planar is documented behaviour here.
quiet_longlat <- function(expr) {
  withCallingHandlers(expr,
    message = function(m) {
      if (grepl("longitude/latitude", conditionMessage(m))) invokeRestart("muffleMessage")
    },
    warning = function(w) {
      if (grepl("longitude/latitude", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}

rename_reserved <- function(nm, reserved) {
  nm[nm %in% reserved] <- paste0(nm[nm %in% reserved], "_attr")
  nm
}

as_sf_table <- function(obj, name) {
  if (inherits(obj, "sfc")) return(sf::st_sf(geometry = obj))
  if (!inherits(obj, "sf")) {
    stop(sprintf("`%s` must be an sf data frame or an sfc geometry list, not %s.", name,
                 describe_class(obj)), call. = FALSE)
  }
  obj
}

#' @rdname from_lines
#' @export
from_polygons <- function(polygons, contiguity = c("queen", "rook"), id = NULL) {
  require_sf("from_polygons()")
  contiguity <- match.arg(contiguity)
  polygons <- as_sf_table(polygons, "polygons")
  geom <- sf::st_geometry(polygons)
  types <- as.character(sf::st_geometry_type(geom, by_geometry = TRUE))
  bad <- which(!types %in% c("POLYGON", "MULTIPOLYGON"))
  if (length(bad)) {
    stop(sprintf("`polygons` must hold POLYGON or MULTIPOLYGON geometries; %s %s %s.",
                 rows_text(bad), if (length(bad) == 1L) "is" else "are",
                 paste(unique(types[bad]), collapse = ", ")), call. = FALSE)
  }
  if (is.null(id)) {
    attrs <- sf::st_drop_geometry(polygons)
    names(attrs) <- rename_reserved(names(attrs), "id")
    nodes <- sf::st_sf(cbind(data.frame(id = seq_len(nrow(polygons))), attrs),
                       geometry = geom)
    id <- "id"
  } else {
    nodes <- polygons
    id <- resolve_column(sf::st_drop_geometry(nodes), id, "id", "polygons")
  }
  pattern <- if (contiguity == "queen") "F***T****" else "F***1****"
  # st_relate is planar; sf says so for longitude and latitude, as documented.
  rel <- quiet_longlat(sf::st_relate(geom, geom, pattern = pattern))
  i <- rep(seq_along(rel), lengths(rel))
  j <- unlist(rel, use.names = FALSE)
  keep <- i < j
  i <- i[keep]
  j <- j[keep]
  o <- order(i, j)
  ids <- nodes[[id]]
  edges <- data.frame(from = ids[i[o]], to = ids[j[o]], stringsAsFactors = FALSE)
  netframe(edges, nodes, directed = FALSE, from = "from", to = "to", id = id)
}

# proximity_network() for sf node tables: representative points through
# the engine's k-d tree, or exact geometry distances when `exact = TRUE`.
sf_proximity_network <- function(nodes, within, k, knn, id, exact, threads) {
  if (is.null(id)) {
    attrs <- sf::st_drop_geometry(nodes)
    if (ncol(attrs)) {
      id <- names(attrs)[[1L]]
    } else {
      nodes <- sf::st_sf(data.frame(id = seq_len(nrow(nodes))), geometry = sf::st_geometry(nodes))
      id <- "id"
    }
  } else {
    id <- resolve_column(sf::st_drop_geometry(nodes), id, "id", "nodes")
  }
  holder <- netframe(data.frame(from = nodes[[id]][0], to = nodes[[id]][0]), nodes,
                     directed = FALSE, id = id)
  geom <- sf::st_geometry(nodes)
  lonlat <- crs_kind(holder) == "lonlat"
  if (exact) {
    r <- exact_proximity(geom, within, k)
  } else {
    xy <- point_matrix(node_points(holder), "`nodes`")
    r <- proximity_pairs(xy[, 1L], xy[, 2L], lonlat, 6371008.7714, within, k, threads)
  }
  r <- symmetrize_knn(r, k, knn)
  ids <- nodes[[id]]
  edges <- data.frame(from = ids[r$from], to = ids[r$to], distance = r$distance,
                      stringsAsFactors = FALSE)
  netframe(edges, nodes, directed = r$directed, id = id, weight = "distance",
           weight_type = "distance")
}

# Pairs within a distance, or k nearest neighbours, by sf::st_distance on
# the geometries themselves (polygons at distance 0 when they touch).
exact_proximity <- function(geom, within, k) {
  n <- length(geom)
  if (!is.null(within)) {
    near <- sf::st_is_within_distance(geom, geom, dist = within)
    i <- rep(seq_len(n), lengths(near))
    j <- unlist(near, use.names = FALSE)
    keep <- i < j
    i <- i[keep]
    j <- j[keep]
    d <- if (length(i)) strip_units(sf::st_distance(geom[i], geom[j], by_element = TRUE)) else numeric()
    return(list(from = i, to = j, distance = d))
  }
  if (n > 5000L) {
    stop(paste0("`exact = TRUE` with `k` computes every pairwise distance, which is ",
                "impractical beyond 5000 features; use `exact = FALSE`."), call. = FALSE)
  }
  d <- strip_units(sf::st_distance(geom))
  dim(d) <- c(n, n)
  from <- integer()
  to <- integer()
  dist <- numeric()
  for (i in seq_len(n)) {
    others <- setdiff(seq_len(n), i)
    o <- others[order(d[i, others])]
    pick <- utils::head(o, k)
    from <- c(from, rep(i, length(pick)))
    to <- c(to, pick)
    dist <- c(dist, d[i, pick])
  }
  list(from = from, to = to, distance = dist)
}

# --- sfnetworks -------------------------------------------------------------

#' Convert to and from sfnetworks
#'
#' `to_sfnetwork()` turns a spatial netframe into an `sfnetwork` object of
#' the sfnetworks package (van der Meer et al. 2026); `as_netframe()` on an
#' `sfnetwork` converts back, and so does `sfnetworks::as_sfnetwork()` on
#' a netframe. Node, edge, and network-level attributes are carried across.
#'
#' sfnetworks represents nodes as POINT features, so node geometries of
#' other types are replaced by their representative points (see
#' [to_sf()]) with a message. Node identifiers become the `name` node
#' column that tidygraph and igraph use (an existing `name` attribute is
#' renamed `name_attr`), and come back from it. Edge geometries are passed
#' on when the edge table has them; edges without geometry get straight
#' lines when `edges_as_lines = TRUE`, and remain spatially implicit by
#' default. `sfnetworks::sfnetwork()` checks that every edge line starts
#' and ends at its nodes' points and stops otherwise; pass `force = TRUE`
#' through `...` to skip that check for networks whose geometries do not
#' meet it (knappnet does not require it; see [check_geometry()]). Columns
#' of class `units`, such as sfnetworks' `weight`, are converted to plain
#' numbers on the way back.
#'
#' @param x A [netframe()] with node geometry or coordinate columns, or an
#'   `sfnetwork`.
#' @param edges_as_lines Whether edges without geometry should get
#'   straight-line geometries; `NULL` keeps the edge table as it is.
#' @param ... Further arguments to `sfnetworks::sfnetwork()`, such as
#'   `force = TRUE`.
#'
#' @return `to_sfnetwork()` returns an `sfnetwork`;
#'   `as_netframe()` a netframe whose node table (and edge table, when
#'   edges are spatially explicit) is an sf data frame.
#' @references
#' van der Meer, L., Abad, L., Gilardi, A., & Lovelace, R. (2026).
#' sfnetworks: Tidy geospatial networks. R package.
#' \doi{10.32614/CRAN.package.sfnetworks}
#' @seealso [to_sf()], [as_netframe()].
#' @examples
#' if (requireNamespace("sfnetworks", quietly = TRUE)) {
#'   pts <- sf::st_as_sf(data.frame(name = c("a", "b", "c"), x = c(0, 1, 1), y = c(0, 0, 1)),
#'                       coords = c("x", "y"), crs = 32633)
#'   g <- netframe(data.frame(from = c("a", "b"), to = c("b", "c")), pts, directed = FALSE)
#'   net <- to_sfnetwork(g, edges_as_lines = TRUE)
#'   net
#'   as_netframe(net)
#' }
#' @export
to_sfnetwork <- function(x, edges_as_lines = NULL, ...) {
  check_netframe(x)
  require_pkg("sfnetworks")
  require_sf("to_sfnetwork()")
  check_flag_or_null(edges_as_lines, "edges_as_lines")
  nodes <- node_table_sf(x)
  geom <- sf::st_geometry(nodes)
  if (!all(as.character(sf::st_geometry_type(geom, by_geometry = TRUE)) == "POINT")) {
    message("Node geometries were replaced by representative points, which sfnetworks needs.")
    sf::st_geometry(nodes) <- node_points(x)
  }
  id <- x$roles$id
  attrs <- sf::st_drop_geometry(nodes)
  ids <- attrs[[id]]
  attrs <- attrs[setdiff(names(attrs), id)]
  names(attrs) <- rename_reserved(names(attrs), "name")
  node_df <- cbind(data.frame(name = ids, stringsAsFactors = FALSE), attrs)
  nodes <- sf::st_sf(node_df, geometry = sf::st_geometry(nodes))

  explicit <- has_edge_geometry(x)
  edges <- if (explicit || isTRUE(edges_as_lines)) {
    edge_table_sf(x, implicit = TRUE)
  } else {
    x$edges
  }
  eattrs <- if (inherits(edges, "sf")) sf::st_drop_geometry(edges) else as.data.frame(edges)
  eattrs <- eattrs[setdiff(names(eattrs), c(x$roles$from, x$roles$to))]
  names(eattrs) <- rename_reserved(names(eattrs), c("from", "to"))
  edge_df <- cbind(data.frame(from = x$index$src + 1L, to = x$index$tgt + 1L), eattrs)
  if (inherits(edges, "sf")) {
    edge_df <- sf::st_sf(edge_df, geometry = sf::st_geometry(edges))
  }
  net <- sfnetworks::sfnetwork(nodes, edge_df, directed = x$directed, node_key = "name",
                               edges_as_lines = inherits(edge_df, "sf"), ...)
  if (length(x$graph)) {
    cls <- class(net)
    for (a in names(x$graph)) net <- igraph::set_graph_attr(net, a, x$graph[[a]])
    class(net) <- cls
  }
  net
}

#' @rdname as_netframe
#' @export
as_netframe.sfnetwork <- function(x, ...) {
  require_pkg("sfnetworks")
  require_sf("as_netframe() for sfnetwork objects")
  nodes <- sf::st_as_sf(x, active = "nodes")
  n <- nrow(nodes)
  if ("name" %in% names(nodes) && !anyDuplicated(nodes$name)) {
    ids <- nodes$name
  } else {
    ids <- seq_len(n)
    attrs <- sf::st_drop_geometry(nodes)
    names(attrs) <- rename_reserved(names(attrs), "name")
    nodes <- sf::st_sf(cbind(data.frame(name = ids), attrs), geometry = sf::st_geometry(nodes))
  }
  nodes <- strip_unit_columns(nodes)
  edges <- tryCatch(suppressMessages(sf::st_as_sf(x, active = "edges")),
                    error = function(e) NULL)
  if (is.null(edges)) {
    edges <- as.data.frame(suppressMessages(tidygraph::as_tibble(x, active = "edges")))
  }
  edges <- strip_unit_columns(edges)
  edges$from <- ids[edges$from]
  edges$to <- ids[edges$to]
  weight <- if (is.numeric(edges$weight)) "weight"
  # sfnetworks keeps its own bookkeeping in graph attributes; carry over scalars only.
  graph <- Filter(function(v) is.atomic(v) && length(v) == 1L, igraph::graph_attr(x))
  netframe(edges, nodes, directed = igraph::is_directed(x), from = "from", to = "to",
           id = "name", weight = weight, graph = graph)
}

strip_unit_columns <- function(df) {
  for (col in names(df)) {
    if (inherits(df[[col]], "units")) df[[col]] <- strip_units(df[[col]])
  }
  df
}

#' @rdname as_netframe
#' @param nodes For sf edge tables, the node table (a data frame or sf
#'   data frame), as in [netframe()].
#' @param edges For sf node tables, the edge table: when given, `x` is
#'   taken as the node table instead.
#' @export
as_netframe.sf <- function(x, nodes = NULL, ..., edges = NULL) {
  if (!is.null(edges)) return(netframe(edges, nodes = x, ...))
  netframe(x, nodes = nodes, ...)
}

# Methods for generics in suggested packages, registered when they load.
st_as_sf.netframe <- function(x, what = c("nodes", "edges"), implicit = TRUE, ...) {
  to_sf(x, what = what, implicit = implicit)
}
as_sfnetwork.netframe <- function(x, ...) to_sfnetwork(x, ...)

# --- GeoPackage ------------------------------------------------------------

# A spatial network as a GeoPackage: layers `nodes` and `edges` (either
# may be an aspatial table) and a small aspatial layer `network` holding
# the direction, the column roles, and network-level attributes as text.
write_gpkg <- function(x, file, weight) {
  require_sf("Writing a GeoPackage")
  nodes <- if (has_node_geometry(x)) node_table_sf(x) else as.data.frame(x$nodes)
  edges <- if (has_edge_geometry(x)) edge_table_sf(x, implicit = FALSE) else as.data.frame(x$edges)
  check_gpkg_columns(nodes, "nodes")
  check_gpkg_columns(edges, "edges")
  graph <- Filter(function(v) is.atomic(v) && length(v) == 1L, x$graph)
  meta <- data.frame(
    key = c("directed", "id", "from", "to", "weight", "weight_type", paste0("graph.", names(graph))),
    value = c(if (x$directed) "true" else "false", x$roles$id, x$roles$from, x$roles$to,
              if (is.null(weight)) "" else weight,
              if (is.null(weight)) "" else if (is.null(x$roles$weight_type)) "distance" else x$roles$weight_type,
              vapply(graph, function(v) as.character(v), "")),
    stringsAsFactors = FALSE
  )
  path <- path.expand(file)
  if (file.exists(path)) unlink(path)
  sf::st_write(nodes, path, layer = "nodes", driver = "GPKG", quiet = TRUE)
  sf::st_write(edges, path, layer = "edges", driver = "GPKG", quiet = TRUE)
  sf::st_write(meta, path, layer = "network", driver = "GPKG", quiet = TRUE)
  invisible(file)
}

check_gpkg_columns <- function(df, what) {
  plain <- if (inherits(df, "sf")) sf::st_drop_geometry(df) else df
  for (col in names(plain)) {
    v <- plain[[col]]
    if (is.list(v)) {
      stop(sprintf("Column `%s` of `%s` is a list; convert it to text or numbers before writing.",
                   col, what), call. = FALSE)
    }
  }
  invisible(NULL)
}

read_gpkg <- function(file, directed, weight, weight_type) {
  require_sf("Reading a GeoPackage")
  path <- path.expand(file)
  layers <- sf::st_layers(path)$name
  missing <- setdiff(c("nodes", "edges"), layers)
  if (length(missing)) {
    stop(sprintf(paste0("GeoPackage `%s` has no layer %s; a network needs layers `nodes` ",
                        "and `edges` (it has: %s)."), file,
                 paste0("`", missing, "`", collapse = " or "),
                 paste0("`", layers, "`", collapse = ", ")), call. = FALSE)
  }
  nodes <- read_gpkg_layer(path, "nodes")
  edges <- read_gpkg_layer(path, "edges")
  meta <- list()
  if ("network" %in% layers) {
    m <- read_gpkg_layer(path, "network")
    if (all(c("key", "value") %in% names(m))) {
      meta <- stats::setNames(as.list(as.character(m$value)), as.character(m$key))
    }
  }
  pick <- function(key, default) {
    v <- meta[[key]]
    if (is.null(v) || is.na(v) || !nzchar(v)) default else v
  }
  if (is.null(directed)) directed <- identical(pick("directed", "true"), "true")
  id <- pick("id", names(nodes)[[1L]])
  from <- pick("from", "from")
  to <- pick("to", "to")
  if (is.null(weight)) {
    w <- pick("weight", NULL)
    weight <- if (!is.null(w) && w %in% names(edges) && is.numeric(edges[[w]])) w else pick_weight(edges, NULL)
    wt <- pick("weight_type", NULL)
    if (!is.null(wt) && wt %in% c("distance", "strength")) weight_type <- wt
  } else if (!weight %in% names(edges)) {
    stop(sprintf("Edge column `%s` not found in the file.", weight), call. = FALSE)
  }
  graph <- meta[startsWith(names(meta), "graph.")]
  names(graph) <- substring(names(graph), 7L)
  netframe(edges, nodes, directed = directed, from = from, to = to, id = id, weight = weight,
           weight_type = weight_type, graph = graph)
}

# st_read() returns a data frame for an aspatial layer and says so; the
# notice is expected here.
read_gpkg_layer <- function(path, layer) {
  out <- withCallingHandlers(
    sf::st_read(path, layer = layer, quiet = TRUE, stringsAsFactors = FALSE),
    message = function(m) {
      if (grepl("no simple feature geometries", conditionMessage(m))) invokeRestart("muffleMessage")
    },
    warning = function(w) {
      if (grepl("no simple feature geometries", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
  # A layer written from a table without geometry comes back aspatial.
  if (inherits(out, "sf")) out else as.data.frame(out, stringsAsFactors = FALSE)
}

# --- Drawing with sf ----------------------------------------------------------

# What plot_spatial() draws from sf tables: each node's geometry kind
# ("point", "line", "polygon", or "other"), the node and edge geometries,
# which edges have their own geometry, and the bounding box of it all.
spatial_plot_geometry <- function(x) {
  require_sf("Plotting a netframe with sf geometry")
  out <- list(nodes = NULL, kind = rep("point", nrow(x$nodes)), edges = NULL, explicit = NULL,
              bbox = NULL)
  boxes <- list()
  if (has_node_geometry(x)) {
    geom <- sf::st_geometry(x$nodes)
    type <- as.character(sf::st_geometry_type(geom, by_geometry = TRUE))
    kind <- rep("other", length(geom))
    kind[type %in% c("POINT", "MULTIPOINT")] <- "point"
    kind[type %in% c("LINESTRING", "MULTILINESTRING")] <- "line"
    kind[type %in% c("POLYGON", "MULTIPOLYGON")] <- "polygon"
    kind[sf::st_is_empty(geom)] <- "point"
    out$nodes <- geom
    out$kind <- kind
    boxes <- c(boxes, list(sf::st_bbox(geom)))
  }
  if (has_edge_geometry(x)) {
    geom <- sf::st_geometry(x$edges)
    out$edges <- geom
    out$explicit <- !sf::st_is_empty(geom)
    if (any(out$explicit)) boxes <- c(boxes, list(sf::st_bbox(geom)))
  }
  if (length(boxes)) {
    b <- do.call(rbind, lapply(boxes, as.numeric))
    bbox <- c(min(b[, 1L]), min(b[, 2L]), max(b[, 3L]), max(b[, 4L]))
    if (all(is.finite(bbox))) out$bbox <- bbox
  }
  out
}

# Polygons filled with the node color, lines in it, and any other
# non-point geometry as sf draws it; points are drawn afterwards, on top.
draw_node_geometry <- function(geo, fill) {
  if (is.null(geo$nodes)) return(invisible(NULL))
  geom <- geo$nodes
  poly <- which(geo$kind == "polygon")
  if (length(poly)) {
    graphics::plot(sf::st_cast(geom[poly], "MULTIPOLYGON"), add = TRUE, col = fill[poly],
                   border = "white", lwd = 0.6)
  }
  line <- which(geo$kind == "line")
  if (length(line)) {
    graphics::plot(sf::st_cast(geom[line], "MULTILINESTRING"), add = TRUE, col = fill[line],
                   lwd = 2)
  }
  other <- which(geo$kind == "other")
  if (length(other)) {
    graphics::plot(geom[other], add = TRUE, col = fill[other], border = "white")
  }
  invisible(NULL)
}

# Edges with their own geometry, plus an arrowhead on the last segment of
# each directed one.
draw_edge_geometry <- function(geo, col, lwd, directed, n) {
  idx <- which(geo$explicit)
  if (!length(idx)) return(invisible(NULL))
  geom <- geo$edges[idx]
  graphics::plot(geom, add = TRUE, col = col[idx], lwd = lwd[idx])
  if (!directed) return(invisible(NULL))
  co <- sf::st_coordinates(geom)
  if (nrow(co) < 2L) return(invisible(NULL))
  fid <- co[, ncol(co)]
  last <- which(!duplicated(fid, fromLast = TRUE))
  prev <- last - 1L
  ok <- prev >= 1L
  ok[ok] <- fid[prev[ok]] == fid[last[ok]] &
    (co[last[ok], 1L] != co[prev[ok], 1L] | co[last[ok], 2L] != co[prev[ok], 2L])
  if (!any(ok)) return(invisible(NULL))
  which_edge <- idx[fid[last[ok]]]
  graphics::arrows(co[prev[ok], 1L], co[prev[ok], 2L], co[last[ok], 1L], co[last[ok], 2L],
                   length = if (n <= 100) 0.08 else 0.05, angle = 22, col = col[which_edge],
                   lwd = lwd[which_edge])
  invisible(NULL)
}
