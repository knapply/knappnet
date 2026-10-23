# Plot methods for the results of cohesive subgroup, topography, and
# conditional uniform graph routines (partitions are drawn by
# plot_partition() in plot-network.R). All use base graphics, colour-blind
# safe palettes (Okabe and Ito's qualitative palette through
# grDevices::palette.colors(), and viridis for sequential values), return
# their input invisibly, and pass `...` on to graphics::title().

# See man/cug_test.Rd.
plot.cug_test <- function(x, ..., breaks = "Sturges") {
  pal <- plot_qualitative(8L)
  null <- x$null
  if (!length(null)) stop("The null distribution is empty; nothing to plot.", call. = FALSE)
  dots <- list(...)
  if (is.null(dots$main)) dots$main <- sprintf("CUG test of %s (conditioning on %s)",
                                               x$statistic, x$conditioning)
  if (is.null(dots$xlab)) dots$xlab <- x$statistic
  xlim <- range(c(null, x$observed), finite = TRUE)
  xlim <- xlim + c(-1, 1) * 0.05 * max(diff(xlim), 1e-8)
  h <- graphics::hist(null, breaks = breaks, plot = FALSE)
  do.call(graphics::plot, c(list(h, col = pal[[2L]], border = "white", xlim = xlim), dots))
  graphics::abline(v = x$observed, col = pal[[5L]], lwd = 2.5)
  graphics::legend(if (x$observed > stats::median(null)) "topleft" else "topright",
                   bty = "n", cex = 0.85, lty = c(1, NA, NA), lwd = c(2.5, NA, NA),
                   col = c(pal[[5L]], NA, NA),
                   legend = c(sprintf("observed = %s", format(x$observed, digits = 4)),
                              sprintf("Pr(X >= obs) = %s", format(x$p_greater, digits = 3)),
                              sprintf("Pr(X <= obs) = %s", format(x$p_less, digits = 3))))
  invisible(x)
}

# See man/core_periphery.Rd.
plot.core_periphery <- function(x, what = NULL, layout = NULL, ..., seed = NULL) {
  if (is.null(what)) what <- if (x$method == "discrete") "matrix" else "network"
  what <- check_choice(what, c("matrix", "network"), "what")
  g <- x$network
  pal <- plot_qualitative(8L)
  if (what == "network") {
    if (x$method == "discrete") {
      plot(g, layout = layout, node_color = I(ifelse(x$core, pal[[5L]], pal[[1L]])), ...,
           legend = FALSE, seed = seed)
      graphics::legend("topleft", legend = c("core", "periphery"), pch = 21,
                       pt.bg = pal[c(5L, 1L)], col = "white", pt.cex = 1.6, bty = "n",
                       cex = 0.8)
    } else {
      plot(g, layout = layout, node_color = x$coreness, node_size = x$coreness, ..., seed = seed)
    }
    return(invisible(x))
  }
  n <- nrow(g$nodes)
  if (n == 0L) stop("The network has no nodes; nothing to plot.", call. = FALSE)
  order <- if (x$method == "discrete") {
    order(!x$core, seq_len(n))
  } else {
    order(-x$coreness, seq_len(n))
  }
  weighted <- !is.null(g$roles$weight)
  a <- as.matrix(g, attr = g$roles$weight)[order, order, drop = FALSE]
  diag(a) <- 0
  if (!weighted) a <- (a > 0) * 1
  dots <- list(...)
  if (is.null(dots$main)) {
    dots$main <- sprintf("%s core-periphery fit (correlation %s)",
                         if (x$method == "discrete") "Discrete" else "Continuous",
                         format(x$correlation, digits = 3))
  }
  labelled <- n <= 40
  old <- graphics::par(mar = c(0.5, if (labelled) 5 else 1, if (labelled) 6.5 else 3, 1),
                       xpd = NA)
  on.exit(graphics::par(old), add = TRUE)
  # Row i at the top, as a matrix is printed.
  colors <- if (weighted) {
    plot_matrix_palette(64L)
  } else {
    c("white", pal[[2L]])
  }
  graphics::image(seq_len(n), seq_len(n), t(a)[, n:1, drop = FALSE], col = colors,
                  axes = FALSE, xlab = "", ylab = "", useRaster = FALSE)
  graphics::box(col = "gray50")
  if (labelled) {
    labels <- as.character(g$nodes[[g$roles$id]])[order]
    graphics::axis(2, at = n:1, labels = labels, las = 2, tick = FALSE, cex.axis = 0.7,
                   line = -0.5)
    graphics::axis(3, at = seq_len(n), labels = labels, las = 2, tick = FALSE,
                   cex.axis = 0.7, line = -0.5)
  }
  boundary <- if (x$method == "discrete") sum(x$core) else sum(x$coreness > mean(x$coreness))
  if (boundary > 0 && boundary < n) {
    graphics::abline(v = boundary + 0.5, h = n - boundary + 0.5, col = pal[[5L]], lwd = 2)
  }
  do.call(graphics::title, c(dots, list(line = if (labelled) 5 else 1.5)))
  invisible(x)
}

# See man/small_world.Rd.
plot.small_world <- function(x, ...) {
  pal <- plot_qualitative(8L)
  dots <- list(...)
  old <- graphics::par(mfrow = c(1L, 2L), mar = c(3, 4, 4, 1), oma = c(0, 0, 2, 0))
  on.exit(graphics::par(old), add = TRUE)
  panel <- function(obs, rnd, ylab) {
    v <- c(observed = obs, random = rnd)
    graphics::barplot(v, col = pal[c(2L, 1L)], border = NA, ylab = ylab,
                      ylim = c(0, max(v[is.finite(v)], 0) * 1.15))
    graphics::text(c(0.7, 1.9), v, format(v, digits = 3), pos = 3, cex = 0.8, xpd = NA)
  }
  panel(x$clustering, x$random_clustering,
        if (x$local) "mean local clustering" else "transitivity")
  panel(x$mean_distance, x$random_mean_distance, "mean distance")
  main <- if (is.null(dots$main)) {
    sprintf("Small-world quotient %s (%d random networks)", format(x$quotient, digits = 3),
            x$reps)
  } else {
    dots$main
  }
  graphics::mtext(main, outer = TRUE, cex = 1.1, font = 2)
  invisible(x)
}

# See man/lambda_sets.Rd.
plot.lambda_sets <- function(x, ...) {
  m <- x$connectivity
  n <- nrow(m)
  if (n < 2L) stop("Fewer than two nodes; nothing to plot.", call. = FALSE)
  top <- max(m)
  # Single linkage on top - lambda recovers the nested lambda sets exactly,
  # since lambda(a, c) >= min(lambda(a, b), lambda(b, c)).
  d <- stats::as.dist(top - m)
  tree <- stats::hclust(d, method = "single")
  dots <- list(...)
  if (is.null(dots$main)) dots$main <- "Lambda sets"
  if (is.null(dots$ylab)) dots$ylab <- "edge connectivity"
  args <- c(list(tree, yaxt = "n", sub = "", xlab = ""), dots)
  do.call(graphics::plot, args)
  at <- sort(unique(c(0, tree$height)))
  graphics::axis(2, at = at, labels = top - at, las = 1)
  invisible(x)
}
