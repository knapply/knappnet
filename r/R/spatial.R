# Geospatial networks (man/spatial_network.Rd).
spatial_network <- function(x, coords = c("lon", "lat"), crs = c("lonlat", "planar"),
                            radius = 6371008.7714,
                            representative = c("point_on_surface", "centroid")) {
  check_netframe(x)
  crs <- match.arg(crs)
  representative <- match.arg(representative)
  radius <- check_number(radius, "radius", lower = 0, open = TRUE)
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
  if (has_node_geometry(x)) {
    require_pkg("sf", "A netframe with sf node geometry")
    return(geometry_coordinates(x))
  }
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

# See man/spatial_network.Rd.
edge_lengths <- function(x) {
  check_netframe(x)
  if (has_edge_geometry(x)) {
    require_pkg("sf", "A netframe with sf edge geometry")
    return(geometry_edge_lengths(x))
  }
  p <- spatial_coords(x)
  ptr <- engine_graph(x)
  engine_call("edge_lengths", .Call(C_edge_lengths, ptr, p$x, p$y, p$crs, p$radius))
}

# See man/spatial_network.Rd.
distance_decay <- function(d, scale, form = c("exponential", "gaussian", "power"),
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

# See man/spatial_network.Rd.
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

# Proximity networks from locations (man/proximity_network.Rd).
proximity_network <- function(nodes, coords = c("lon", "lat"), crs = c("lonlat", "planar"),
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
    require_pkg("sf", "proximity_network() with sf nodes")
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
  edges <- data.frame(from = ids[r$from], to = ids[r$to], distance = r$distance,
                      stringsAsFactors = FALSE)
  g <- netframe(edges, nodes, directed = r$directed, id = id, weight = "distance",
                weight_type = "distance")
  spatial_network(g, coords = c(cx, cy), crs = crs, radius = radius)
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
  r <- engine_call("proximity_network",
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

# Plot a spatial network on its map (man/plot_spatial.Rd).
plot_spatial <- function(x, node_color = NULL, node_size = NULL, node_label = NULL,
                         edge_width = NULL, edge_color = NULL, ..., legend = TRUE) {
  check_netframe(x)
  check_flag(legend, "legend")
  if (!has_positions(x)) {
    stop(paste0("`x` has no node locations: give sf nodes to netframe(), or record ",
                "coordinate columns with spatial_network()."), call. = FALSE)
  }
  n <- nrow(x$nodes)
  m <- nrow(x$edges)
  p <- node_coordinates(x)
  xy <- p$coords
  geo <- if (has_node_geometry(x) || has_edge_geometry(x)) spatial_plot_geometry(x) else NULL

  fill <- spatial_colors(node_color, x$nodes, n, "node_color", "#0072B2")
  alpha <- min(0.7, max(0.1, 10 / sqrt(max(m, 1))))
  ecol <- spatial_colors(edge_color, x$edges, m, "edge_color",
                         grDevices::adjustcolor("gray35", alpha.f = alpha))
  cex <- spatial_sizes(node_size, x, n)
  lwd <- spatial_widths(edge_width, x$edges, m)
  labels <- spatial_labels(node_label, x, n)

  dots <- list(...)
  key <- if (legend && !is.null(fill$legend)) fill$legend else NULL
  right <- 0.4
  legend_cex <- 0.75
  if (!is.null(key)) {
    text <- c(key$title, key$labels)
    widths <- tryCatch(graphics::strwidth(text, units = "inches", cex = legend_cex),
                       error = function(e) 0.09 * legend_cex * nchar(text))
    right <- right + (max(widths) + 0.45) / (graphics::par("csi") * graphics::par("mex"))
  }
  top <- if (!is.null(dots$main)) 2.2 else 0.4
  old <- graphics::par(mar = c(0.4, 0.4, top, right), xpd = NA)
  on.exit(graphics::par(old), add = TRUE)
  graphics::plot.new()
  if (n == 0) {
    if (length(dots)) do.call(graphics::title, dots)
    return(invisible(x))
  }

  xr <- range(xy[, 1])
  yr <- range(xy[, 2])
  if (!is.null(geo$bbox)) {
    xr <- range(xr, geo$bbox[c(1L, 3L)])
    yr <- range(yr, geo$bbox[c(2L, 4L)])
  }
  span <- max(diff(xr), diff(yr))
  pad <- if (span > 0) 0.05 * span else 1
  low <- if (any(nzchar(labels))) 0.04 * max(span, 1) else 0
  asp <- if (p$lonlat) 1 / cos(mean(yr) * pi / 180) else 1
  graphics::plot.window(xlim = xr + c(-pad, pad), ylim = yr + c(-pad - low, pad), asp = asp)
  usr <- graphics::par("usr")
  per_inch <- diff(usr[1:2]) / graphics::par("pin")[1]
  # Radius of a pch 21 circle: 0.375 of the symbol size.
  radius <- 0.375 * cex * graphics::par("cex") * graphics::par("ps") / 72 * per_inch

  if (!is.null(geo)) draw_node_geometry(geo, fill$col)
  if (m > 0) {
    from <- x$index$src + 1L
    to <- x$index$tgt + 1L
    keep <- from != to
    if (!is.null(geo$edges)) {
      keep <- keep & !geo$explicit
      draw_edge_geometry(geo, ecol$col, lwd, x$directed, n)
    }
    draw_segments(xy, from[keep], to[keep], ecol$col[keep], lwd[keep], x$directed, radius,
                  per_inch, n)
  }
  points_only <- if (is.null(geo)) rep(TRUE, n) else geo$kind == "point"
  if (any(points_only)) {
    graphics::points(xy[points_only, 1], xy[points_only, 2], pch = 21,
                     bg = fill$col[points_only], col = if (n <= 2000) "white" else NA,
                     cex = cex[points_only], lwd = 0.6)
  }
  if (any(nzchar(labels))) {
    label_cex <- if (n <= 15) 0.9 else if (n <= 40) 0.75 else 0.6
    below <- radius + 0.3 * graphics::strheight("M", cex = label_cex)
    graphics::text(xy[, 1], xy[, 2] - below, labels, adj = c(0.5, 1), cex = label_cex,
                   col = "gray15")
  }
  if (length(dots)) do.call(graphics::title, dots)
  if (!is.null(key)) {
    graphics::legend(x = usr[2] + 0.15 * per_inch, y = usr[4], legend = key$labels,
                     title = key$title, bty = "n", cex = legend_cex, xjust = 0, yjust = 1,
                     title.adj = 0, pch = 22, pt.bg = key$col, col = "gray40", pt.cex = 1.6)
  }
  invisible(x)
}

# Straight edges between node positions, shortened to the node borders and
# tipped with arrowheads when directed.
draw_segments <- function(xy, from, to, col, lwd, directed, radius, per_inch, n) {
  if (!length(from)) return(invisible(NULL))
  x0 <- xy[from, 1]
  y0 <- xy[from, 2]
  x1 <- xy[to, 1]
  y1 <- xy[to, 2]
  if (!directed) {
    graphics::segments(x0, y0, x1, y1, col = col, lwd = lwd)
    return(invisible(NULL))
  }
  dx <- x1 - x0
  dy <- y1 - y0
  len <- sqrt(dx^2 + dy^2)
  start <- radius[from]
  end <- radius[to] + 0.01 * per_inch
  ok <- len > start + end
  if (any(ok)) {
    ux <- dx[ok] / len[ok]
    uy <- dy[ok] / len[ok]
    graphics::arrows(x0[ok] + ux * start[ok], y0[ok] + uy * start[ok],
                     x1[ok] - ux * end[ok], y1[ok] - uy * end[ok],
                     length = if (n <= 100) 0.08 else 0.05, angle = 22, col = col[ok],
                     lwd = lwd[ok])
  }
  invisible(NULL)
}

# How an appearance argument was given: nothing, a column name, one value
# per element, or a single value.
spatial_aes <- function(value, table, k, arg) {
  if (is.null(value) || k == 0L) return(list(kind = "default"))
  if (is.character(value) && length(value) == 1L && value %in% names(table)) {
    return(list(kind = "mapped", values = table[[value]], title = value))
  }
  if (length(value) == 1L) return(list(kind = "constant", value = value))
  if (length(value) == k) return(list(kind = "mapped", values = value, title = NULL))
  stop(sprintf("`%s` must be a column name, one value, or one value per element.", arg),
       call. = FALSE)
}

# Okabe-Ito for categories, viridis for numbers; missing values light gray.
spatial_colors <- function(value, table, k, arg, default) {
  a <- spatial_aes(value, table, k, arg)
  if (a$kind == "default") return(list(col = rep(default, k)))
  if (a$kind == "constant") {
    ok <- tryCatch({
      grDevices::col2rgb(a$value)
      TRUE
    }, error = function(e) FALSE)
    if (!ok) stop(sprintf("`%s` is neither a column nor a color.", arg), call. = FALSE)
    return(list(col = rep(a$value, k)))
  }
  v <- a$values
  missing <- "gray80"
  if (is.numeric(v)) {
    ok <- is.finite(v)
    col <- rep(missing, k)
    if (!any(ok)) {
      return(list(col = col, legend = list(title = a$title, labels = "NA", col = missing)))
    }
    rng <- range(v[ok])
    pal <- grDevices::hcl.colors(100, "viridis")
    pick <- function(z) {
      if (diff(rng) > 0) pal[1L + round(99 * (z - rng[1]) / diff(rng))] else pal[50L]
    }
    col[ok] <- pick(v[ok])
    breaks <- pretty(rng, 4)
    breaks <- breaks[breaks >= rng[1] & breaks <= rng[2]]
    if (length(breaks) < 2L) breaks <- unique(rng)
    return(list(col = col, legend = list(title = a$title, labels = format(breaks),
                                         col = pick(breaks))))
  }
  f <- if (is.factor(v)) droplevels(v) else factor(v, levels = sort(unique(v[!is.na(v)])))
  # Okabe-Ito without black first, so the default node color stays distinct.
  pal <- unname(grDevices::palette.colors(8L, palette = "Okabe-Ito")[c(2:8, 1L)])
  if (nlevels(f) > 8L) {
    warning(sprintf("`%s` has %d categories; the eight Okabe-Ito colors are recycled.", arg,
                    nlevels(f)), call. = FALSE)
  }
  key <- rep_len(pal, nlevels(f))
  col <- key[as.integer(f)]
  col[is.na(f)] <- missing
  labels <- levels(f)
  if (anyNA(f)) {
    labels <- c(labels, "NA")
    key <- c(key, missing)
  }
  list(col = col, legend = list(title = a$title, labels = labels, col = key))
}

spatial_sizes <- function(value, x, n) {
  a <- spatial_aes(value, x$nodes, n, "node_size")
  base <- if (n > 0) min(2.5, max(0.25, 4 / n^0.25)) else 1
  if (a$kind == "default") {
    d <- if (n > 0) degree(x, loops = FALSE) else numeric()
    return(spatial_scale(d, 0.6 * base, 1.6 * base))
  }
  if (a$kind == "constant") {
    if (!is.numeric(a$value) || !is.finite(a$value) || a$value < 0) {
      stop("`node_size` must be a non-negative number.", call. = FALSE)
    }
    return(rep(a$value, n))
  }
  if (!is.numeric(a$values)) stop("`node_size` must be numeric.", call. = FALSE)
  spatial_scale(a$values, 0.6 * base, 1.6 * base)
}

spatial_widths <- function(value, edges, m) {
  a <- spatial_aes(value, edges, m, "edge_width")
  if (a$kind == "default") return(rep(if (m <= 200) 1 else if (m <= 5000) 0.7 else 0.5, m))
  if (a$kind == "constant") {
    if (!is.numeric(a$value) || !is.finite(a$value) || a$value < 0) {
      stop("`edge_width` must be a non-negative number.", call. = FALSE)
    }
    return(rep(a$value, m))
  }
  if (!is.numeric(a$values)) stop("`edge_width` must be numeric.", call. = FALSE)
  spatial_scale(a$values, 0.4, 3, area = FALSE)
}

# Maps values onto [lo, hi], by square root so that areas grow linearly.
spatial_scale <- function(v, lo, hi, area = TRUE) {
  ok <- is.finite(v)
  out <- rep(lo, length(v))
  if (!any(ok)) return(out)
  rng <- range(v[ok])
  t <- if (diff(rng) > 0) (v[ok] - rng[1]) / diff(rng) else rep(0.5, sum(ok))
  out[ok] <- lo + (hi - lo) * (if (area) sqrt(t) else t)
  out
}

spatial_labels <- function(value, x, n) {
  ids <- function() as.character(x$nodes[[x$roles$id]])
  if (is.null(value)) return(if (n > 0 && n <= 40) ids() else rep("", n))
  if (isFALSE(value)) return(rep("", n))
  if (isTRUE(value)) return(ids())
  a <- spatial_aes(value, x$nodes, n, "node_label")
  labels <- as.character(rep_len(if (a$kind == "constant") a$value else a$values, n))
  labels[is.na(labels)] <- ""
  labels
}

# See man/plot_spatial.Rd.
plot_edge_lengths <- function(x, scale = NULL, form = c("exponential", "gaussian", "power"),
                              exponent = 2, ...) {
  form <- match.arg(form)
  len <- edge_lengths(x)
  unit <- attr(len, "unit")
  len <- as.vector(len)
  dots <- list(...)
  if (is.null(dots$main)) dots$main <- "Edge lengths"
  if (is.null(dots$xlab)) {
    dots$xlab <- if (!is.null(unit) && nzchar(unit)) sprintf("Edge length (%s)", unit) else "Edge length"
  }
  if (is.null(dots$col)) dots$col <- "#56B4E9"
  if (is.null(dots$border)) dots$border <- "white"
  if (!length(len)) {
    graphics::plot.new()
    do.call(graphics::title, dots[intersect(names(dots), c("main", "xlab", "ylab", "sub"))])
    return(invisible(x))
  }
  h <- do.call(graphics::hist, c(list(x = len), dots))
  if (!is.null(scale)) {
    d <- seq(0, max(h$breaks), length.out = 200L)
    old <- graphics::par(new = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    graphics::plot(d, distance_decay(d, scale, form, exponent), type = "l", lwd = 2,
                   col = "#D55E00", axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1),
                   xlim = range(h$breaks))
    graphics::axis(4, col = "#D55E00", col.axis = "#D55E00")
    graphics::mtext(sprintf("Decay weight (%s, scale %s)", form, format(scale)), side = 4,
                    line = 2.5, col = "#D55E00", cex = 0.9)
  }
  invisible(x)
}
