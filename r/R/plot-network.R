# Draw a network with nodes sized and colored by a score (man/plot_centrality.Rd).
plot_centrality <- function(x, scores, title = "score", color = TRUE, ...) {
  check_netframe(x)
  color <- check_flag(color, "color")
  if (is.character(scores) && length(scores) == 1L && scores %in% names(x$nodes) &&
      identical(title, "score")) {
    title <- scores
  }
  scores <- node_column(x, scores, "scores")
  if (!is.numeric(scores)) {
    stop(paste0("`scores` must be a numeric vector with one value per node, or the name ",
                "of a numeric column of `x$nodes`."), call. = FALSE)
  }
  check_string(title, "title")
  w <- with_node_column(x, as.numeric(scores), title)
  invisible(plot(w$x, node_size = w$name, node_color = if (color) w$name else NULL, ...))
}

# Draw a network colored by a partition (man/plot_partition.Rd).
plot_partition <- function(x, membership, title = "group", ...) {
  check_netframe(x)
  if (is.character(membership) && length(membership) == 1L && membership %in% names(x$nodes) &&
      identical(title, "group")) {
    title <- membership
  }
  check_string(title, "title")
  q <- attr(membership, "modularity")
  g <- plot_groups(x, membership, "membership")
  w <- with_node_column(x, g, title)
  dots <- list(...)
  if (is.null(dots$main) && !is.null(q) && is.numeric(q) && length(q) == 1L) {
    dots$main <- sprintf("%d groups, modularity %s", nlevels(g), format(q, digits = 3))
  }
  invisible(do.call(plot, c(list(x = w$x, node_color = w$name), dots)))
}

# Draw the k-core shells of a network (man/plot_coreness.Rd).
plot_coreness <- function(x, mode = c("all", "out", "in"), seed = NULL, threads = 0L, ...) {
  check_netframe(x)
  mode <- match.arg(mode)
  n <- nrow(x$nodes)
  k <- coreness(x, mode = mode)
  dots <- list(...)
  xy <- if (n > 0L) {
    layout_network(x, seed = seed, threads = threads)
  } else {
    cbind(x = numeric(), y = numeric())
  }
  levels_ <- sort(unique(k), decreasing = TRUE)
  out <- matrix(0, n, 2L, dimnames = list(NULL, c("x", "y")))
  radii <- numeric(length(levels_))
  for (i in seq_along(levels_)) {
    members <- which(k == levels_[i])
    r <- if (i == 1L && length(members) == 1L) 0 else i
    radii[i] <- r
    if (r == 0) next
    ang <- atan2(xy[members, 2] - mean(xy[, 2]), xy[members, 1] - mean(xy[, 1]))
    members <- members[order(ang)]
    theta <- seq(0, 2 * pi, length.out = length(members) + 1L)[-1L] + pi / 2
    out[members, 1] <- r * cos(theta)
    out[members, 2] <- r * sin(theta)
  }
  # The legend is drawn here rather than by plot.netframe() so that the
  # margins are known and the rings can be added in the same coordinates.
  mapped <- map_colors(k)
  args <- c(list(x = x, layout = out, node_color = I(mapped$col), legend = FALSE), dots)
  if (is.null(args$node_size)) args$node_size <- 1
  drawn <- do.call(plot, args)
  if (n > 0L) {
    top <- if (!is.null(dots$main)) 2.2 else 0.4
    old <- graphics::par(mar = c(0.4, 0.4, top, 0.4), xpd = NA)
    on.exit(graphics::par(old), add = TRUE)
    theta <- seq(0, 2 * pi, length.out = 181L)
    for (r in radii[radii > 0]) {
      graphics::lines(r * cos(theta), r * sin(theta), col = "gray80", lty = 3)
    }
    graphics::legend("topleft", legend = mapped$labels, title = "coreness", pch = 21,
                     pt.bg = mapped$key, col = "gray40", pt.cex = 1.6, bty = "n", cex = 0.75,
                     title.adj = 0)
  }
  invisible(drawn)
}

# Draw an ego network (man/plot_ego.Rd).
plot_ego <- function(x, node, order = 1L, mode = c("all", "out", "in"), ...) {
  check_netframe(x)
  mode <- match.arg(mode)
  e <- ego_network(x, node, order = order, mode = mode, include_ego = TRUE)
  is_ego <- node_ids(e) == as.character(node)
  role <- factor(ifelse(is_ego, "ego", "alter"), levels = c("ego", "alter"))
  w <- with_node_column(e, role, "role")
  dots <- list(...)
  args <- c(list(x = w$x, node_color = w$name), dots)
  if (is.null(args$node_size)) args$node_size <- I(ifelse(is_ego, 2.4, 1.3))
  do.call(plot, args)
  invisible(e)
}

# Robustness of a network to node removal (man/plot_robustness.Rd).
plot_robustness <- function(x, nodes = NULL, strategy = c("degree", "betweenness", "random"),
                            measure = c("largest_component", "fragmentation"), steps = 25,
                            mode = c("weak", "strong"), seed = NULL, threads = 0L,
                            label = "selected nodes", ...) {
  check_netframe(x)
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  if (!is.null(strategy)) {
    if (!is.character(strategy) || anyNA(strategy) ||
        !all(strategy %in% c("degree", "betweenness", "random"))) {
      stop("`strategy` must be some of \"degree\", \"betweenness\", and \"random\", or NULL.",
           call. = FALSE)
    }
    strategy <- unique(strategy)
  }
  steps <- check_number(steps, "steps", lower = 2)
  if (!is.character(label) || length(label) != 1L || is.na(label) || label %in% strategy) {
    stop("`label` must be a single string that is not the name of a strategy.", call. = FALSE)
  }
  n <- nrow(x$nodes)
  if (n < 2L) stop("Robustness needs a network with at least two nodes.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  evaluate <- function(removed) {
    keep <- setdiff(seq_len(n), removed)
    sub <- induced_netframe(x, keep)
    if (measure == "largest_component") {
      comp <- components(sub, mode)
      if (comp$no == 0L) 0 else max(comp$csize) / n
    } else {
      if (length(keep) < 2L) NaN else fragmentation(sub, mode)
    }
  }
  curve_for <- function(removal, counts) {
    vapply(counts, function(k) evaluate(removal[seq_len(k)]), numeric(1))
  }
  counts <- unique(round(seq(0, n - 1, length.out = min(steps, n))))
  curves <- list()
  for (s in strategy) {
    removal <- switch(s,
      degree = order(degree(x), decreasing = TRUE),
      betweenness = order(betweenness(x, threads = threads), decreasing = TRUE),
      random = sample.int(n)
    )
    curves[[s]] <- data.frame(curve = s, removed = counts, fraction = counts / n,
                              value = curve_for(removal, counts), stringsAsFactors = FALSE)
  }
  if (!is.null(nodes)) {
    pos <- node_positions(x, nodes, "nodes", empty_ok = FALSE, distinct = TRUE)
    k <- 0:length(pos)
    curves[[label]] <- data.frame(curve = label, removed = k, fraction = k / n,
                                  value = curve_for(pos, k), stringsAsFactors = FALSE)
  }
  if (!length(curves)) stop("Give `nodes`, a `strategy`, or both.", call. = FALSE)
  out <- do.call(rbind, unname(curves))
  rownames(out) <- NULL
  dots <- list(...)
  ylab <- if (measure == "largest_component") "largest component (share of nodes)" else "fragmentation"
  if (is.null(dots$main)) dots$main <- "Robustness to node removal"
  if (is.null(dots$xlab)) dots$xlab <- "nodes removed (share)"
  if (is.null(dots$ylab)) dots$ylab <- ylab
  if (is.null(dots$xlim)) dots$xlim <- c(0, max(out$fraction))
  if (is.null(dots$ylim)) dots$ylim <- c(0, 1)
  names_ <- names(curves)
  cols <- plot_qualitative(length(names_))
  if (!is.null(nodes)) cols[length(cols)] <- plot_accent
  old <- graphics::par(mar = c(4.6, 4.6, 4.1, 1.1))
  on.exit(graphics::par(old), add = TRUE)
  do.call(graphics::plot, c(list(x = 0, y = 0, type = "n"), dots))
  graphics::grid(col = "gray90", lty = 1)
  for (i in seq_along(names_)) {
    cv <- curves[[i]]
    graphics::lines(cv$fraction, cv$value, col = cols[i], lwd = 2)
    graphics::points(cv$fraction, cv$value, col = cols[i], pch = 16, cex = 0.7)
  }
  graphics::legend("topright", legend = names_, col = cols, lwd = 2, pch = 16, bty = "n",
                   cex = 0.85)
  invisible(out)
}
