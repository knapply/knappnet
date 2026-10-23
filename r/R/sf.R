# Spatial networks whose node or edge tables are sf data frames
# (pebesma2018sf). All geometry work is delegated to sf; the engine receives
# derived numbers only. Polygon contiguity uses the DE-9IM relations of
# egenhofer1991pointset and clementini1993smallset with the queen and rook
# definitions of anselin1988spatial, compared with spdep::poly2nb
# (bivand2013applied, bivand2018comparing). Line snapping at tolerance 0
# matches sfnetworks (vandermeer2026sfnetworks); the two-table layout of
# nodes and edges follows OSMnx's graph_to_gdfs (boeing2017osmnx) and
# sfnetworks. sf, sfnetworks, and spdep are suggested packages: every use is
# guarded by requireNamespace() and the package works without them.

# Spatial networks with sf (man/to_sf.Rd).
to_sf <- function(x, what = c("nodes", "edges"), implicit = TRUE) {
  check_netframe(x)
  what <- match.arg(what)
  check_flag(implicit, "implicit")
  require_pkg("sf", "to_sf()")
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

# See man/to_sf.Rd.
check_geometry <- function(x, tolerance = 0) {
  check_netframe(x)
  require_pkg("sf", "check_geometry()")
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

# See man/to_sf.Rd.
transform_network <- function(x, crs) {
  check_netframe(x)
  require_pkg("sf", "transform_network()")
  if (!has_node_geometry(x) && !has_edge_geometry(x)) {
    stop("`x` has no geometry to transform; see ?to_sf for spatial networks.", call. = FALSE)
  }
  if (has_node_geometry(x)) x$nodes <- sf::st_transform(node_table_sf(x), crs)
  if (has_edge_geometry(x)) x$edges <- sf::st_transform(edge_table_sf(x, implicit = FALSE), crs)
  x$roles$crs <- sf::st_crs(if (has_node_geometry(x)) x$nodes else x$edges)
  x
}

# See man/to_sf.Rd.
spatial_filter <- function(x, y, predicate = sf::st_intersects, ...) {
  check_netframe(x)
  require_pkg("sf", "spatial_filter()")
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

# See man/to_sf.Rd.
spatial_join <- function(x, y, join = sf::st_intersects, ...) {
  check_netframe(x)
  require_pkg("sf", "spatial_join()")
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
  require_pkg("sf", "A netframe with sf nodes or edges")
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
    require_pkg("sf", "A netframe with sf geometry")
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

# Spatial networks from lines and polygons (man/from_lines.Rd).
from_lines <- function(lines, tolerance = 0, directed = FALSE) {
  require_pkg("sf", "from_lines()")
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

# See man/from_lines.Rd.
from_polygons <- function(polygons, contiguity = c("queen", "rook"), id = NULL) {
  require_pkg("sf", "from_polygons()")
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

# Convert to and from sfnetworks (man/to_sfnetwork.Rd).
to_sfnetwork <- function(x, edges_as_lines = NULL, ...) {
  check_netframe(x)
  require_pkg("sfnetworks")
  require_pkg("sf", "to_sfnetwork()")
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
  attrs <- x$graph
  if (!is.null(x$roles$weight)) {
    # Recorded so that as_netframe() restores the weight role.
    attrs$knappnet_weight <- x$roles$weight
    attrs$knappnet_weight_type <- x$roles$weight_type
  }
  if (length(attrs)) {
    cls <- class(net)
    for (a in names(attrs)) net <- igraph::set_graph_attr(net, a, attrs[[a]])
    class(net) <- cls
  }
  net
}

# See man/as_netframe.Rd.
as_netframe.sfnetwork <- function(x, ...) {
  require_pkg("sfnetworks")
  require_pkg("sf", "as_netframe() for sfnetwork objects")
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
  # sfnetworks keeps its own bookkeeping in graph attributes; carry over scalars only.
  graph <- Filter(function(v) is.atomic(v) && length(v) == 1L, igraph::graph_attr(x))
  # to_sfnetwork() records the weight column and its type; otherwise a
  # numeric column called `weight` is the weight, as in igraph.
  weight <- graph[["knappnet_weight"]]
  weight_type <- graph[["knappnet_weight_type"]]
  graph[c("knappnet_weight", "knappnet_weight_type")] <- NULL
  if (is.null(weight) || !weight %in% names(edges) || !is.numeric(edges[[weight]])) {
    weight <- if ("weight" %in% names(edges) && is.numeric(edges[["weight"]])) "weight"
    weight_type <- NULL
  }
  if (is.null(weight_type) || !weight_type %in% c("distance", "strength")) weight_type <- "distance"
  netframe(edges, nodes, directed = igraph::is_directed(x), from = "from", to = "to",
           id = "name", weight = weight, weight_type = weight_type, graph = graph)
}

strip_unit_columns <- function(df) {
  for (col in names(df)) {
    if (inherits(df[[col]], "units")) df[[col]] <- strip_units(df[[col]])
  }
  df
}

# See man/as_netframe.Rd.
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
# A role recorded in the GeoPackage metadata layer: "" when unset.
role_text <- function(v) if (is.null(v)) "" else as.character(v)

write_gpkg <- function(x, file, weight) {
  require_pkg("sf", "Writing a GeoPackage")
  nodes <- if (has_node_geometry(x)) node_table_sf(x) else as.data.frame(x$nodes)
  edges <- if (has_edge_geometry(x)) edge_table_sf(x, implicit = FALSE) else as.data.frame(x$edges)
  check_gpkg_columns(nodes, "nodes")
  check_gpkg_columns(edges, "edges")
  graph <- Filter(function(v) is.atomic(v) && length(v) == 1L, x$graph)
  meta <- data.frame(
    key = c("directed", "id", "from", "to", "weight", "weight_type", if (length(graph)) paste0("graph.", names(graph))),
    value = c(if (x$directed) "true" else "false", role_text(x$roles$id),
              role_text(x$roles$from), role_text(x$roles$to),
              if (is.null(weight)) "" else weight,
              if (is.null(weight)) "" else if (is.null(x$roles$weight_type)) "distance" else x$roles$weight_type,
              vapply(graph, function(v) as.character(v), "", USE.NAMES = FALSE)),
    stringsAsFactors = FALSE
  )
  path <- path.expand(file)
  if (file.exists(path)) unlink(path)
  # delete_dsn replaces a file that could not be removed (on Windows a file
  # read earlier in the session can stay locked).
  sf::st_write(nodes, path, layer = "nodes", driver = "GPKG", quiet = TRUE,
               delete_dsn = file.exists(path))
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
  require_pkg("sf", "Reading a GeoPackage")
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
  meta <- stats::setNames(list(), character())
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
  require_pkg("sf", "Plotting a netframe with sf geometry")
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
