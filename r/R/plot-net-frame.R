# Plot a network (man/plot.net_frame.Rd). The appearance arguments are
# resolved by the helpers in utils-plot.R, which net_plot_spatial() shares.
plot.net_frame <- function(x, layout = NULL, node_size = NULL, node_color = NULL,
                          node_label = NULL, edge_width = NULL, edge_color = NULL, ...,
                          legend = TRUE, seed = NULL) {
  check_net_frame(x)
  legend <- check_flag(legend, "legend")
  geo <- identical(layout, "geo")
  if ((is.null(layout) || geo) && has_positions(x)) {
    return(net_plot_spatial(x, node_color = node_color, node_size = node_size,
                        node_label = node_label, edge_width = edge_width,
                        edge_color = edge_color, ..., legend = legend))
  }
  if (geo) {
    stop(paste0("`layout = \"geo\"` needs node locations: give sf nodes to net_frame(), or ",
                "record coordinate columns with net_spatial()."), call. = FALSE)
  }
  n <- nrow(x$nodes)
  m <- nrow(x$edges)
  xy <- plot_layout(x, layout, seed)

  size <- plot_node_size(x, node_size, n)
  fill <- if (is.null(node_color) && n > 0 && net_is_twomode(x)) {
    twomode_fill(x)
  } else {
    plot_color(node_color, x$nodes, n, "node_color", "nodes", plot_blue)
  }
  labels <- plot_labels(x, node_label, n)
  width <- plot_edge_width(x, edge_width, m)
  alpha <- min(0.6, max(0.05, 10 / sqrt(max(m, 1))))
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
    return(invisible(xy))
  }

  span <- max(diff(range(xy[, 1])), diff(range(xy[, 2])))
  pad <- if (span > 0) 0.05 * span else 1
  low <- if (any(nzchar(labels))) 0.04 * max(span, 1) else 0
  graphics::plot.window(xlim = range(xy[, 1]) + c(-pad, pad),
                        ylim = range(xy[, 2]) + c(-pad - low, pad), asp = 1)
  usr <- graphics::par("usr")
  per_inch <- diff(usr[1:2]) / graphics::par("pin")[1]
  radius <- symbol_radius(size$cex, per_inch)

  if (m > 0) {
    from <- x$index$src + 1L
    to <- x$index$tgt + 1L
    keep <- from != to
    draw_edges(xy, from[keep], to[keep], ecol$col[keep], width$lwd[keep], x$directed,
               radius, per_inch, n)
  }
  graphics::points(xy[, 1], xy[, 2], pch = 21, bg = fill$col,
                   col = if (n <= 2000) "white" else NA, cex = size$cex, lwd = 0.6)
  draw_node_labels(xy, labels, radius, n)
  if (length(dots)) do.call(graphics::title, dots)
  draw_legends(legends, usr, per_inch, legend_cex)
  invisible(xy)
}
