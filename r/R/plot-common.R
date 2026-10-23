# Helpers shared by the plot methods and plotting functions. Every plot in
# the package draws with base graphics, uses the same palettes as
# plot.netframe() (hcl.colors "Dark 3" for categories, "viridis" for scores,
# the single-hue "Blues 3" for matrices), restores graphical parameters it
# changes, and returns its input or the plotted data invisibly.

# Qualitative colors for k categories, as in plot.netframe().
plot_qualitative <- function(k) {
  if (k <= 0L) return(character())
  grDevices::hcl.colors(k, "Dark 3")
}

# Sequential colors from light to dark for matrices and heat maps.
plot_sequential <- function(n) grDevices::hcl.colors(n, "Blues 3", rev = TRUE)

# Accent colors: observed statistics and highlights (vermilion), the main
# series (the blue of plot.netframe()), and reference lines.
plot_accent <- "#D55E00"
plot_blue <- "#2C6FAC"
plot_gray <- "gray55"

# Width of the widest label in lines of text, for margins; strwidth() needs
# an open device, so fall back to a character-count estimate without one.
label_lines <- function(labels, cex = 1) {
  if (!length(labels)) return(0)
  labels <- as.character(labels)
  inches <- tryCatch(max(graphics::strwidth(labels, units = "inches", cex = cex)),
                     error = function(e) 0.09 * cex * max(nchar(labels)))
  csi <- tryCatch(graphics::par("csi"), error = function(e) 0.2)
  inches / csi
}

# The `main` title and other arguments meant for title(), taken out of
# `...` so the rest can go to the drawing call.
split_dots <- function(dots, names) {
  take <- intersect(names(dots), names)
  list(taken = dots[take], rest = dots[setdiff(names(dots), take)])
}

# Draws a numeric matrix as a heat map with row 1 at the top, as a matrix
# is printed: rows are labelled on the left and columns at the bottom.
# `col` maps the range of the finite values from light to dark; cells that
# are NA or NaN are left blank. `lines` are 0-based positions after which a
# boundary is drawn (between rows/columns k and k + 1), and `values` writes
# each cell's value in it.
draw_heatmap <- function(m, row_labels = rownames(m), col_labels = colnames(m),
                         col = plot_sequential(64), range = NULL, values = FALSE,
                         digits = 2L, row_lines = integer(), col_lines = integer(),
                         key = TRUE, main = NULL, xlab = "", ylab = "", cex_labels = NULL,
                         ...) {
  nr <- nrow(m)
  nc <- ncol(m)
  if (nr == 0L || nc == 0L) {
    stop("The matrix to plot has no rows or no columns.", call. = FALSE)
  }
  if (is.null(row_labels)) row_labels <- as.character(seq_len(nr))
  if (is.null(col_labels)) col_labels <- as.character(seq_len(nc))
  show_row <- nr <= 80
  show_col <- nc <= 80
  if (is.null(cex_labels)) {
    cex_labels <- if (max(nr, nc) <= 20) 0.8 else if (max(nr, nc) <= 50) 0.6 else 0.45
  }
  z <- t(m[rev(seq_len(nr)), , drop = FALSE])
  storage.mode(z) <- "double"
  finite <- z[is.finite(z)]
  if (is.null(range)) {
    range <- if (length(finite)) base::range(finite) else c(0, 1)
  }
  if (!(diff(range) > 0)) range <- range + c(-0.5, 0.5)
  breaks <- seq(range[1], range[2], length.out = length(col) + 1L)
  z[is.finite(z) & z < range[1]] <- range[1]
  z[is.finite(z) & z > range[2]] <- range[2]
  left <- if (show_row) label_lines(row_labels, cex_labels) + 1 else 1
  bottom <- if (show_col) label_lines(col_labels, cex_labels) + 1 else 1
  if (nzchar(ylab)) left <- left + 1.5
  if (nzchar(xlab)) bottom <- bottom + 1.5
  right <- if (isTRUE(key)) 4.5 else 1
  top <- if (!is.null(main)) 3 else 1
  old <- graphics::par(mar = c(max(bottom, 2), max(left, 2), top, right))
  on.exit(graphics::par(old), add = TRUE)
  # Cell boundaries rather than centers, so that one-row and one-column
  # matrices draw too.
  graphics::image(seq(0.5, nc + 0.5, by = 1), seq(0.5, nr + 0.5, by = 1), z, col = col,
                  breaks = breaks, axes = FALSE, xlab = "", ylab = "", main = main, ...)
  if (show_col) {
    graphics::axis(1, at = seq_len(nc), labels = col_labels, las = 2, tick = FALSE,
                   cex.axis = cex_labels, line = -0.5)
  }
  if (show_row) {
    graphics::axis(2, at = rev(seq_len(nr)), labels = row_labels, las = 1, tick = FALSE,
                   cex.axis = cex_labels, line = -0.5)
  }
  if (nzchar(xlab)) graphics::mtext(xlab, side = 1, line = bottom - 1.2)
  if (nzchar(ylab)) graphics::mtext(ylab, side = 2, line = left - 1.2)
  if (isTRUE(values)) {
    cx <- rep(seq_len(nc), times = nr)
    cy <- rep(rev(seq_len(nr)), each = nc)
    v <- as.vector(t(m))
    lab <- ifelse(is.finite(v), format(round(v, digits), trim = TRUE), "")
    shade <- rep("gray15", length(v))
    dark <- is.finite(v) & (v - range[1]) / diff(range) > 0.6
    shade[dark] <- "white"
    graphics::text(cx, cy, lab, cex = cex_labels, col = shade)
  }
  if (length(row_lines)) {
    graphics::abline(h = nr - row_lines + 0.5, col = plot_accent, lwd = 1.5)
  }
  if (length(col_lines)) {
    graphics::abline(v = col_lines + 0.5, col = plot_accent, lwd = 1.5)
  }
  graphics::box(col = "gray40")
  if (isTRUE(key)) draw_color_key(col, range)
  invisible(m)
}

# A vertical color key in the right margin of the current plot.
draw_color_key <- function(col, range) {
  usr <- graphics::par("usr")
  old <- graphics::par(xpd = NA)
  on.exit(graphics::par(old), add = TRUE)
  w <- diff(usr[1:2])
  h <- diff(usr[3:4])
  x0 <- usr[2] + 0.04 * w
  x1 <- usr[2] + 0.08 * w
  ys <- seq(usr[3] + 0.1 * h, usr[4] - 0.1 * h, length.out = length(col) + 1L)
  graphics::rect(x0, ys[-length(ys)], x1, ys[-1L], col = col, border = NA)
  graphics::rect(x0, ys[1L], x1, ys[length(ys)], border = "gray40")
  ticks <- pretty(range, 4)
  ticks <- ticks[ticks >= range[1] & ticks <= range[2]]
  if (length(ticks) < 2L) ticks <- range
  at <- ys[1L] + (ticks - range[1]) / diff(range) * (ys[length(ys)] - ys[1L])
  graphics::segments(x1, at, x1 + 0.01 * w, at, col = "gray40")
  graphics::text(x1 + 0.015 * w, at, format(ticks, trim = TRUE), adj = c(0, 0.5), cex = 0.7)
}

# A node grouping (a column name or a vector aligned with x$nodes) as a
# factor with levels in order of first appearance, for coloring.
plot_groups <- function(x, groups, name) {
  g <- node_groups(x, groups, name)
  factor(as.character(g$values), levels = as.character(g$labels))
}

# Adds a column to the node table under a unique name so plot.netframe()
# can map it and label the legend with `title`.
with_node_column <- function(x, values, title) {
  name <- make.unique(c(names(x$nodes), title))[ncol(x$nodes) + 1L]
  x$nodes[[name]] <- values
  list(x = x, name = name)
}
