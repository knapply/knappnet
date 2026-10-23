# Plot a spatial network on its map (man/net_plot_spatial.Rd). Appearance
# arguments are resolved by the helpers in utils-plot.R that plot.net_frame()
# uses; sf geometries are drawn by the helpers in sf.R.
net_plot_spatial <- function(x, node_color = NULL, node_size = NULL, node_label = NULL,
                         edge_width = NULL, edge_color = NULL, ..., legend = TRUE) {
  check_net_frame(x)
  legend <- check_flag(legend, "legend")
  if (!has_positions(x)) {
    stop(paste0("`x` has no node locations: give sf nodes to net_frame(), or record ",
                "coordinate columns with net_spatial()."), call. = FALSE)
  }
  n <- nrow(x$nodes)
  m <- nrow(x$edges)
  p <- node_coordinates(x)
  xy <- p$coords
  geo <- if (has_node_geometry(x) || has_edge_geometry(x)) spatial_plot_geometry(x) else NULL

  size <- plot_node_size(x, node_size, n)
  fill <- plot_color(node_color, x$nodes, n, "node_color", "nodes", plot_blue)
  labels <- plot_labels(x, node_label, n)
  width <- plot_edge_width(x, edge_width, m)
  alpha <- min(0.7, max(0.1, 10 / sqrt(max(m, 1))))
  ecol <- plot_color(edge_color, x$edges, m, "edge_color", "edges",
                     grDevices::adjustcolor("gray35", alpha.f = alpha))
  if (ecol$mapped) ecol$col <- grDevices::adjustcolor(ecol$col, alpha.f = max(alpha, 0.7))

  legends <- if (legend) {
    Filter(Negate(is.null), list(fill$legend, size$legend, ecol$legend, width$legend))
  } else {
    list()
  }
  dots <- list(...)
  top <- if (!is.null(dots$main)) 2.2 else 0.4
  legend_cex <- 0.75
  right <- 0.4 + legend_margin(legends, legend_cex)
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
  # Longitude and latitude are drawn with the local aspect ratio of a
  # plate carree map, so that distances look right near the data.
  asp <- if (p$lonlat) 1 / cos(mean(yr) * pi / 180) else 1
  graphics::plot.window(xlim = xr + c(-pad, pad), ylim = yr + c(-pad - low, pad), asp = asp)
  usr <- graphics::par("usr")
  per_inch <- diff(usr[1:2]) / graphics::par("pin")[1]
  radius <- symbol_radius(size$cex, per_inch)

  if (!is.null(geo)) draw_node_geometry(geo, fill$col)
  if (m > 0) {
    from <- x$index$src + 1L
    to <- x$index$tgt + 1L
    keep <- from != to
    if (!is.null(geo$edges)) {
      keep <- keep & !geo$explicit
      draw_edge_geometry(geo, ecol$col, width$lwd, x$directed, n)
    }
    draw_edges(xy, from[keep], to[keep], ecol$col[keep], width$lwd[keep], x$directed,
               radius, per_inch, n)
  }
  points_only <- if (is.null(geo)) rep(TRUE, n) else geo$kind == "point"
  if (any(points_only)) {
    graphics::points(xy[points_only, 1], xy[points_only, 2], pch = 21,
                     bg = fill$col[points_only], col = if (n <= 2000) "white" else NA,
                     cex = size$cex[points_only], lwd = 0.6)
  }
  draw_node_labels(xy, labels, radius, n)
  if (length(dots)) do.call(graphics::title, dots)
  draw_legends(legends, usr, per_inch, legend_cex)
  invisible(x)
}

# Histogram of edge lengths with an optional decay curve (man/net_plot_spatial.Rd).
net_plot_edge_lengths <- function(x, scale = NULL, form = c("exponential", "gaussian", "power"),
                              exponent = 2, ...) {
  form <- match.arg(form)
  len <- net_edge_lengths(x)
  unit <- attr(len, "unit")
  len <- as.vector(len)
  dots <- list(...)
  if (is.null(dots$main)) dots$main <- "Edge lengths"
  if (is.null(dots$xlab)) {
    dots$xlab <- if (!is.null(unit) && nzchar(unit)) sprintf("Edge length (%s)", unit) else "Edge length"
  }
  if (is.null(dots$col)) dots$col <- plot_blue
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
    graphics::plot(d, net_distance_decay(d, scale, form, exponent), type = "l", lwd = 2,
                   col = plot_accent, axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1),
                   xlim = range(h$breaks))
    graphics::axis(4, col = plot_accent, col.axis = plot_accent)
    graphics::mtext(sprintf("Decay weight (%s, scale %s)", form, format(scale)), side = 4,
                    line = 2.5, col = plot_accent, cex = 0.9)
  }
  invisible(x)
}
