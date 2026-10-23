# Heat maps of adjacency, incidence, and mixing matrices, drawn by
# draw_heatmap() in utils-plot.R.

# Plot an adjacency or incidence matrix (man/net_plot_matrix.Rd).
net_plot_matrix <- function(x, order = NULL, attr = NULL, values = FALSE, key = NULL,
                        labels = NULL, ...) {
  values <- check_flag(values, "values")
  check_flag_or_null(key, "key")
  check_flag_or_null(labels, "labels")
  if (inherits(x, "netframe")) {
    check_netframe(x)
    if (net_is_twomode(x)) {
      m <- net_incidence(x, attr = attr)
      row_lines <- col_lines <- integer()
      if (!is.null(order)) {
        g <- plot_groups(x, order, "order")
        second <- twomode_types(x) == 1L
        r <- group_order(g[!second])
        c_ <- group_order(g[second])
        m <- m[r$order, c_$order, drop = FALSE]
        row_lines <- r$lines
        col_lines <- c_$lines
      }
    } else {
      m <- as.matrix(x, attr = attr)
      row_lines <- col_lines <- integer()
      if (!is.null(order)) {
        g <- plot_groups(x, order, "order")
        o <- group_order(g)
        m <- m[o$order, o$order, drop = FALSE]
        row_lines <- col_lines <- o$lines
      }
    }
  } else {
    if (!is.matrix(x) || !(is.numeric(x) || is.logical(x))) {
      stop("`x` must be a netframe or a numeric matrix.", call. = FALSE)
    }
    m <- x
    storage.mode(m) <- "double"
    row_lines <- col_lines <- integer()
    if (!is.null(order)) {
      if (!is.atomic(order) || length(order) != nrow(m) || anyNA(order)) {
        stop("`order` must have one value per row of the matrix, without missing values.",
             call. = FALSE)
      }
      g <- factor(as.character(order), levels = unique(as.character(order)))
      o <- group_order(g)
      if (nrow(m) == ncol(m)) {
        m <- m[o$order, o$order, drop = FALSE]
        col_lines <- o$lines
      } else {
        m <- m[o$order, , drop = FALSE]
      }
      row_lines <- o$lines
    }
  }
  if (is.null(key)) key <- length(unique(m[is.finite(m)])) > 2L
  rl <- if (is.null(rownames(m))) as.character(seq_len(nrow(m))) else rownames(m)
  cl <- if (is.null(colnames(m))) as.character(seq_len(ncol(m))) else colnames(m)
  show <- if (isTRUE(labels)) Inf else if (isFALSE(labels)) -1 else 80
  if (nrow(m) > show) rl <- rep("", nrow(m))
  if (ncol(m) > show) cl <- rep("", ncol(m))
  finite <- m[is.finite(m)]
  top <- if (length(finite)) max(1, finite) else 1
  draw_heatmap(m, row_labels = rl, col_labels = cl, values = values, key = key,
               range = c(0, top), row_lines = row_lines, col_lines = col_lines, ...)
  invisible(m)
}

# The order that sorts nodes by group (stable), and the 0-based positions
# after which a group ends.
group_order <- function(g) {
  o <- order(as.integer(g))
  sizes <- tabulate(as.integer(g), nbins = nlevels(g))
  sizes <- sizes[sizes > 0L]
  ends <- cumsum(sizes)
  list(order = o, lines = ends[-length(ends)])
}

# Plot a mixing matrix (man/net_plot_mixing.Rd).
net_plot_mixing <- function(x, dim1 = NULL, dim2 = NULL, direction = c("all", "out", "in"),
                               values = TRUE, ...) {
  values <- check_flag(values, "values")
  if (inherits(x, "netframe")) {
    if (is.null(dim1)) stop("`dim1` is needed to compute the mixing matrix.", call. = FALSE)
    m <- net_mixing(x, dim1, dim2, direction = match.arg(direction))
  } else {
    if (!is.matrix(x) || !is.numeric(x)) {
      stop("`x` must be a netframe or a numeric matrix.", call. = FALSE)
    }
    m <- x
  }
  finite <- m[is.finite(m)]
  top <- if (length(finite)) max(1, finite) else 1
  draw_heatmap(m, values = values, digits = 0L, range = c(0, top), ...)
  invisible(m)
}
