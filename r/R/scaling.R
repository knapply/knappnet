# Multidimensional scaling: classical scaling with its print and plot
# methods, stress, and distance matrices for scaling a network.

# Classical multidimensional scaling (man/net_mds.Rd).
net_mds <- function(d, k = 2) {
  d <- as_dissimilarity(d)
  n <- nrow(d)
  k <- check_count(k, "k", min = 1L)
  if (n < 2L) {
    stop("Classical scaling needs at least two points.", call. = FALSE)
  }
  if (k > n - 1L) {
    stop(sprintf("`k` must be at most %d, one less than the number of points.", n - 1L),
         call. = FALSE)
  }
  r <- engine_call("net_mds", .Call(C_classical_mds, d, n, k), edge_rows = FALSE)
  dims <- r[[2L]]
  points <- r[[1L]][, seq_len(dims), drop = FALSE]
  dimnames(points) <- list(rownames(d), if (dims) paste0("Dim", seq_len(dims)))
  stress <- if (dims) {
    engine_call("net_mds", .Call(C_mds_stress, d, n, points, dims), edge_rows = FALSE)
  } else {
    NaN
  }
  structure(list(points = points, eigenvalues = r[[3L]], gof = r[[4L]], stress = stress),
            class = "classical_mds")
}

# See man/net_mds.Rd.
print.classical_mds <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sprintf("<classical_mds> %d points in %d dimension%s\n", nrow(x$points), ncol(x$points),
              if (ncol(x$points) == 1L) "" else "s"))
  cat(sprintf("Stress: %s; goodness of fit: %s, %s\n", format(x$stress, digits = digits),
              format(x$gof[[1L]], digits = digits), format(x$gof[[2L]], digits = digits)))
  cat("Leading eigenvalues:\n")
  print(utils::head(x$eigenvalues, 5L), digits = digits, ...)
  invisible(x)
}

# Plot a multidimensional scaling (man/plot.classical_mds.Rd).
plot.classical_mds <- function(x, which = c("configuration", "shepard", "scree"), d = NULL,
                               dims = c(1, 2), labels = NULL, ...) {
  which <- match.arg(which)
  check_flag_or_null(labels, "labels")
  dots <- list(...)
  pts <- x$points
  n <- nrow(pts)
  if (which == "configuration") {
    if (ncol(pts) == 0L) {
      stop("The scaling has no dimensions with a positive eigenvalue to draw.", call. = FALSE)
    }
    if (!is.numeric(dims) || length(dims) != 2L || anyNA(dims) || any(dims < 1) ||
        any(dims != floor(dims))) {
      stop("`dims` must be two dimension numbers.", call. = FALSE)
    }
    if (ncol(pts) == 1L) {
      pts <- cbind(pts, 0)
      colnames(pts)[2L] <- ""
      dims <- c(1, 2)
    }
    if (any(dims > ncol(pts))) {
      stop(sprintf("`dims` must be at most %d, the number of dimensions kept.", ncol(pts)),
           call. = FALSE)
    }
    px <- pts[, dims[1L]]
    py <- pts[, dims[2L]]
    show <- if (is.null(labels)) n <= 60 else labels
    nm <- rownames(pts)
    if (show && is.null(nm)) nm <- as.character(seq_len(n))
    if (is.null(dots$xlab)) dots$xlab <- colnames(pts)[dims[1L]]
    if (is.null(dots$ylab)) dots$ylab <- colnames(pts)[dims[2L]]
    if (is.null(dots$main)) dots$main <- "Classical scaling"
    if (is.null(dots$asp)) dots$asp <- 1
    if (is.null(dots$pch)) dots$pch <- 21
    if (is.null(dots$bg)) dots$bg <- plot_blue
    if (is.null(dots$col)) dots$col <- "white"
    if (is.null(dots$cex)) dots$cex <- 1.4
    old <- graphics::par(mar = c(4.6, 4.6, 4.1, 1.1))
    on.exit(graphics::par(old), add = TRUE)
    do.call(graphics::plot, c(list(x = px, y = py), dots))
    graphics::abline(h = 0, v = 0, col = "gray85", lty = 3)
    if (show) {
      graphics::text(px, py, nm, pos = 3, offset = 0.6, cex = if (n <= 20) 0.8 else 0.65,
                     col = "gray15")
    }
    graphics::mtext(sprintf("stress %s", format(x$stress, digits = 3)), side = 3,
                    line = 0.3, cex = 0.8, col = "gray30")
  } else if (which == "shepard") {
    if (is.null(d)) {
      stop("The Shepard diagram needs the dissimilarities: pass `d`.", call. = FALSE)
    }
    d <- as_dissimilarity(d)
    if (nrow(d) != n) {
      stop(sprintf("`d` has %d points but the scaling has %d.", nrow(d), n), call. = FALSE)
    }
    if (ncol(pts) == 0L) {
      stop("The scaling has no dimensions with a positive eigenvalue to draw.", call. = FALSE)
    }
    fitted <- as.matrix(stats::dist(pts))
    lower <- lower.tri(d)
    delta <- d[lower]
    fit_d <- fitted[lower]
    if (is.null(dots$xlab)) dots$xlab <- "Dissimilarity"
    if (is.null(dots$ylab)) dots$ylab <- "Distance in the configuration"
    if (is.null(dots$main)) dots$main <- "Shepard diagram"
    if (is.null(dots$pch)) dots$pch <- 16
    if (is.null(dots$col)) {
      dots$col <- grDevices::adjustcolor(plot_blue, alpha.f = if (length(delta) > 500) 0.3 else 0.7)
    }
    lim <- range(c(delta, fit_d, 0), na.rm = TRUE)
    if (is.null(dots$xlim)) dots$xlim <- lim
    if (is.null(dots$ylim)) dots$ylim <- lim
    old <- graphics::par(mar = c(4.6, 4.6, 4.1, 1.1))
    on.exit(graphics::par(old), add = TRUE)
    do.call(graphics::plot, c(list(x = delta, y = fit_d), dots))
    graphics::abline(0, 1, col = plot_accent, lwd = 1.5)
    graphics::mtext(sprintf("stress %s over %d pairs", format(x$stress, digits = 3),
                            length(delta)), side = 3, line = 0.3, cex = 0.8, col = "gray30")
  } else {
    ev <- x$eigenvalues
    k <- ncol(pts)
    cols <- ifelse(ev < 0, plot_accent, plot_gray)
    cols[seq_len(k)] <- plot_blue
    if (is.null(dots$main)) dots$main <- "Eigenvalues of the scaling"
    if (is.null(dots$ylab)) dots$ylab <- "Eigenvalue"
    if (is.null(dots$xlab)) dots$xlab <- "Dimension"
    if (is.null(dots$col)) dots$col <- cols
    if (is.null(dots$border)) dots$border <- NA
    if (is.null(dots$names.arg)) {
      dots$names.arg <- if (length(ev) <= 30) as.character(seq_along(ev)) else rep("", length(ev))
    }
    old <- graphics::par(mar = c(4.6, 4.6, 4.1, 1.1))
    on.exit(graphics::par(old), add = TRUE)
    do.call(graphics::barplot, c(list(height = ev), dots))
    graphics::abline(h = 0, col = "gray40")
    graphics::legend("topright", bty = "n", cex = 0.8, fill = c(plot_blue, plot_gray, plot_accent),
                     border = NA, legend = c(sprintf("kept (%d)", k), "positive, not kept",
                                             "negative (non-Euclidean)"))
  }
  invisible(x)
}

# See man/net_mds.Rd.
net_mds_stress <- function(d, points) {
  d <- as_dissimilarity(d)
  n <- nrow(d)
  if (!is.numeric(points) || !is.matrix(points)) {
    stop("`points` must be a numeric matrix with one row per point.", call. = FALSE)
  }
  if (nrow(points) != n) {
    stop(sprintf("`points` has %d rows but `d` has %d points.", nrow(points), n),
         call. = FALSE)
  }
  if (anyNA(points) || any(!is.finite(points))) {
    stop("`points` must be finite.", call. = FALSE)
  }
  storage.mode(points) <- "double"
  engine_call("net_mds_stress", .Call(C_mds_stress, d, n, points, ncol(points)), edge_rows = FALSE)
}

# A dissimilarity matrix as a square double matrix. Structural checks
# (symmetry, zero diagonal, finiteness) are made by the engine.
as_dissimilarity <- function(d) {
  if (inherits(d, "dist")) d <- as.matrix(d)
  if (!is.numeric(d) || !is.matrix(d) || nrow(d) != ncol(d)) {
    stop("`d` must be a square numeric matrix of dissimilarities, or a \"dist\" object.",
         call. = FALSE)
  }
  if (anyNA(d)) {
    stop("`d` must not contain missing values.", call. = FALSE)
  }
  storage.mode(d) <- "double"
  d
}

# Distances between nodes for scaling (man/net_mds_distances.Rd).
net_mds_distances <- function(x, method = c("geodesic", "euclidean", "hamming", "correlation"),
                          unreachable = "max_plus_one", weighted = NULL, weights_as = NULL,
                          threads = 0L) {
  check_netframe(x)
  method <- match.arg(method)
  if (method == "geodesic") {
    d <- net_distances(x, mode = "all", weighted = weighted, weights_as = weights_as,
                   threads = threads)
    if (any(is.infinite(d))) {
      finite <- d[is.finite(d)]
      longest <- if (length(finite)) max(finite) else 0
      fill <- if (is.numeric(unreachable)) {
        if (length(unreachable) != 1L || is.na(unreachable) || unreachable < 0) {
          stop("`unreachable` must be \"max_plus_one\", \"diameter\", or a non-negative number.",
               call. = FALSE)
        }
        unreachable
      } else {
        switch(check_choice(unreachable, c("max_plus_one", "diameter"), "unreachable"),
               max_plus_one = longest + 1, diameter = longest)
      }
      d[is.infinite(d)] <- fill
    }
    return(d)
  }
  m <- net_struct_equiv(x, method = method, weighted = weighted, threads = threads)
  if (method == "correlation") {
    m <- 1 - m
    diag(m) <- 0
  }
  m
}
