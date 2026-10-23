# Plot a network (man/plot.net_frame.Rd). The appearance arguments are
# resolved by the helpers in utils-plot.R, which net_plot_spatial() shares.
# The defaults make each configuration of the data model visible: modes
# in colours (two-mode networks in two rows), relations and layers as edge
# colours, negative ties dashed, unobserved pairs dotted gray, self-loops
# as small circles, and parallel edges as arcs.
plot.net_frame <- function(x, layout = NULL, node_size = NULL, node_color = NULL,
                           node_label = NULL, edge_width = NULL, edge_color = NULL,
                           edge_linetype = NULL, curved = NULL, ..., legend = TRUE,
                           seed = NULL) {
  check_net_frame(x)
  legend <- check_flag(legend, "legend")
  check_flag_or_null(curved, "curved")
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
  if (is.null(layout) && n > 0 && net_is_twomode(x)) layout <- "twomode"
  xy <- plot_layout(x, layout, seed)

  size <- plot_node_size(x, node_size, n)
  fill <- if (is.null(node_color) && n > 0 && net_is_twomode(x)) {
    twomode_fill(x)
  } else if (is.null(node_color) && n > 0 && !is.null(x$roles$mode)) {
    plot_color(x$roles$mode, x$nodes, n, "node_color", "nodes", plot_blue)
  } else {
    plot_color(node_color, x$nodes, n, "node_color", "nodes", plot_blue)
  }
  labels <- plot_labels(x, node_label, n)
  width <- plot_edge_width(x, edge_width, m)
  alpha <- min(0.6, max(0.05, 10 / sqrt(max(m, 1))))
  ecol <- if (is.null(edge_color) && m > 0) {
    default_edge_color(x, m, alpha)
  } else {
    plot_color(edge_color, x$edges, m, "edge_color", "edges",
               grDevices::adjustcolor("gray35", alpha.f = alpha))
  }
  if (ecol$mapped) ecol$col <- grDevices::adjustcolor(ecol$col, alpha.f = max(alpha, 0.7))
  lty <- plot_edge_linetype(x, edge_linetype, m)
  # Unobserved pairs are drawn dotted and light gray, whatever else is mapped.
  unobserved <- edge_unobserved(x)
  if (any(unobserved)) {
    ecol$col[unobserved] <- plot_missing
    lty$lty[unobserved] <- 3L
    lty$legend <- c(lty$legend, list(list(type = "line", title = "unobserved",
                                          labels = "pair not observed", col = plot_missing,
                                          lwd = 1.5, lty = 3L)))
  }

  legends <- if (legend) {
    Filter(Negate(is.null), c(list(fill$legend, size$legend, ecol$legend, width$legend),
                              lty$legend))
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
    curve <- edge_curvature(from, to, curved)
    draw_edges(xy, from, to, ecol$col, width$lwd, x$directed, radius, per_inch, n,
               lty = lty$lty, curve = curve)
  }
  graphics::points(xy[, 1], xy[, 2], pch = 21, bg = fill$col,
                   col = if (n <= 2000) "white" else NA, cex = size$cex, lwd = 0.6)
  draw_node_labels(xy, labels, radius, n)
  if (length(dots)) do.call(graphics::title, dots)
  draw_legends(legends, usr, per_inch, legend_cex)
  invisible(xy)
}

# Default edge colours: by relation in a multiplex network, by layer in a
# multilayer network (couplings and interlayer edges gray), gray otherwise.
default_edge_color <- function(x, m, alpha) {
  gray <- grDevices::adjustcolor("gray35", alpha.f = alpha)
  if (inherits(x, "multilayer_net_frame")) {
    el <- net_edge_layers(x)
    labels <- frame_summary_cached(x)$layers
    layer <- factor(el$source_layer, levels = labels)
    pal <- plot_qualitative(length(labels))
    col <- pal[as.integer(layer)]
    col[el$class != "intralayer"] <- "gray55"
    lab <- c(labels, if (any(el$class != "intralayer")) "between layers")
    key <- c(pal, if (any(el$class != "intralayer")) "gray55")
    return(list(col = col, mapped = TRUE,
                legend = list(type = "line", title = "layer", labels = lab, col = key,
                              lwd = 2)))
  }
  if (!is.null(x$roles$relation)) {
    return(plot_color(x$roles$relation, x$edges, m, "edge_color", "edges", gray))
  }
  list(col = rep(gray, m), mapped = FALSE)
}

# Line types for edges: a column, a vector, one value, or I(); by default
# negative ties of a signed network are dashed.
plot_edge_linetype <- function(x, value, m) {
  a <- plot_aes(value, x$edges, m, "edge_linetype", "edges")
  if (a$kind == "default") {
    lty <- rep(1L, m)
    out <- list(lty = lty, legend = list())
    if (!is.null(x$roles$sign) && m > 0) {
      s <- x$edges[[x$roles$sign]]
      neg <- if (is.logical(s)) s %in% FALSE else s %in% -1
      out$lty[neg] <- 2L
      if (any(neg)) {
        out$legend <- list(list(type = "line", title = "sign",
                                labels = c("positive", "negative"), col = "gray35",
                                lwd = 1.5, lty = c(1L, 2L)))
      }
    }
    return(out)
  }
  if (a$kind == "constant" || a$kind == "asis") {
    return(list(lty = rep_len(a$values, m), legend = list()))
  }
  v <- a$values
  f <- if (is.factor(v)) droplevels(v) else factor(v, levels = sort(unique(v[!is.na(v)])))
  types <- c(1L, 2L, 3L, 4L, 5L, 6L)
  k <- nlevels(f)
  if (k > length(types)) {
    stop(sprintf("`edge_linetype` maps at most %d categories; this column has %d.",
                 length(types), k), call. = FALSE)
  }
  lty <- types[as.integer(f)]
  lty[is.na(lty)] <- 1L
  list(lty = lty, legend = list(list(type = "line", title = a$title, labels = levels(f),
                                     col = "gray35", lwd = 1.5, lty = types[seq_len(k)])))
}

# Rows of the edge table flagged as unobserved pairs.
edge_unobserved <- function(x) {
  col <- x$roles$edge_missing
  if (is.null(col)) return(rep(FALSE, nrow(x$edges)))
  x$edges[[col]] %in% TRUE
}

# Curvature of each edge: parallel edges (and, in a directed network,
# reciprocated ones) are spread into arcs so that each stays visible.
# `curved` NULL curves only where edges share a pair of nodes.
edge_curvature <- function(from, to, curved) {
  m <- length(from)
  out <- rep(0, m)
  if (isFALSE(curved) || !m) return(out)
  a <- pmin(from, to)
  b <- pmax(from, to)
  key <- paste(a, b)
  size <- as.vector(table(key)[key])
  rank <- stats::ave(seq_len(m), key, FUN = seq_along)
  multi <- size > 1L & from != to
  step <- 0.25
  out[multi] <- (rank[multi] - (size[multi] + 1) / 2) * step
  # Even groups have no straight edge; odd groups keep one in the middle.
  flip <- from > to
  out[flip] <- -out[flip]
  if (isTRUE(curved)) out[!multi & from != to] <- 0.15
  out
}
