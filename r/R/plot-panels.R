# Small multiples: one panel per layer of a multilayer network, per time
# window or wave of a temporal network, and per ego of an egocentric
# dataset. Panels of one network share its layout, so that a node keeps
# its place from panel to panel and only the ties change.

# Layers as panels (man/net_plot_panels.Rd).
net_plot_layers <- function(x, layout = NULL, ..., ncol = NULL, seed = NULL) {
  check_multilayer(x)
  xy <- plot_layout(x, layout, seed)
  layers <- net_split_layers(x)
  plot_panels(layers, names(layers), xy, ncol, ...)
  invisible(xy)
}

# Snapshots over time as panels (man/net_plot_panels.Rd).
net_plot_snapshots <- function(x, width = NULL, breaks = NULL, layout = NULL, ...,
                               ncol = NULL, seed = NULL) {
  check_temporal(x)
  static <- static_net_frame(x)
  xy <- plot_layout(static, layout, seed)
  if (!is.null(x$roles$wave) && is.null(width) && is.null(breaks)) {
    nets <- net_waves(x, nodes = "all")
    titles <- paste("wave", names(nets))
    windows <- lapply(as.numeric(names(nets)), function(k) c(k, k + 1))
  } else {
    if (is.null(width) && is.null(breaks)) {
      tt <- temporal_times(x)
      last <- if (is.null(tt$end)) tt$start else tt$end
      if (!length(tt$start)) stop("`x` has no edges to draw.", call. = FALSE)
      lo <- min(tt$start)
      hi <- max(last)
      b <- if (hi > lo) seq(lo, hi, length.out = 5L) else c(lo, lo + 1)
      b[length(b)] <- b[length(b)] + max(1e-9, 1e-9 * abs(b[length(b)]))
      breaks <- from_seconds(b, x)
    }
    nets <- net_snapshots(x, width = width, breaks = breaks, nodes = "all")
    w <- attr(nets, "windows")
    titles <- sprintf("%s to %s", format(w$start), format(w$end))
    windows <- lapply(seq_len(nrow(w)), function(k) {
      c(time_arg(w$start[k], x, "start"), time_arg(w$end[k], x, "end"))
    })
  }
  timed_nodes <- !is.null(x$roles$node_start) || !is.null(x$roles$node_end) ||
    !is.null(x$node_spells)
  dots <- list(...)
  plot_panels(nets, titles, xy, ncol, ..., node_color_fun = if (timed_nodes &&
                                                                   is.null(dots$node_color)) {
    function(k) {
      present <- node_present(x, windows[[k]][1L], windows[[k]][2L])
      I(ifelse(present, plot_blue, "gray88"))
    }
  })
  invisible(xy)
}

# Ego networks as panels (man/net_plot_panels.Rd).
net_plot_egos <- function(x, egos = NULL, order = 1L, ..., ncol = NULL, max_panels = 12L,
                          seed = NULL) {
  check_net_frame(x)
  max_panels <- check_count(max_panels, "max_panels", min = 1L)
  if (is.null(x$roles$ego)) {
    if (is.null(egos)) egos <- utils::head(node_ids(x), max_panels)
    x <- net_ego_frame(x, egos = egos, order = order, threads = 1L)
  }
  nets <- net_split_egos(x)
  if (!is.null(egos) && !is.null(x$roles$ego)) {
    keep <- names(nets) %in% as.character(egos)
    if (any(keep)) nets <- nets[keep]
  }
  if (length(nets) > max_panels) nets <- nets[seq_len(max_panels)]
  dots <- list(...)
  plot_panels(nets, names(nets), NULL, ncol, ..., seed = seed,
              node_color_fun = if (is.null(dots$node_color)) {
                function(k) {
                  g <- nets[[k]]
                  is_ego <- node_ids(g) == names(nets)[[k]]
                  I(ifelse(is_ego, plot_accent, plot_blue))
                }
              })
  invisible(nets)
}

# Draws one network per panel in a grid of `ncol` columns: with the shared
# layout `xy` (aligned with every network's nodes), or each network's own
# layout when `xy` is NULL. Legends are left out, since the panels share
# their encodings; `node_color_fun(k)`, when given, colours panel k.
plot_panels <- function(nets, titles, xy, ncol, ..., node_color_fun = NULL, seed = NULL) {
  k <- length(nets)
  if (!k) stop("There is nothing to draw.", call. = FALSE)
  if (is.null(ncol)) ncol <- min(k, if (k <= 4L) k else ceiling(sqrt(k)))
  ncol <- check_count(ncol, "ncol", min = 1L)
  nrow <- ceiling(k / ncol)
  old <- graphics::par(mfrow = c(nrow, ncol))
  on.exit(graphics::par(old), add = TRUE)
  for (i in seq_len(k)) {
    args <- c(list(nets[[i]], main = titles[[i]], legend = FALSE), list(...))
    if (!is.null(xy)) args$layout <- xy
    if (is.null(xy) && is.null(args$seed)) args$seed <- seed
    if (!is.null(node_color_fun)) args$node_color <- node_color_fun(i)
    if (is.null(args$cex.main)) args$cex.main <- 0.9
    do.call(plot.net_frame, args)
  }
  invisible(NULL)
}

# The static network of a temporal one: every edge, no time roles.
static_net_frame <- function(x) {
  rebuild_net_frame(x, x$edges, x$nodes, start = NULL, end = NULL, duration = NULL,
                    wave = NULL, node_start = NULL, node_end = NULL,
                    tables = list(node_spells = NULL, node_changes = NULL))
}
