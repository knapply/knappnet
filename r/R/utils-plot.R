# Helpers shared by every plot in the package. All plots draw with base
# graphics, use the palettes below, restore the graphical parameters they
# change, and return their input or the plotted data invisibly.
#
# Colours (one scheme for the whole package):
# * Categories: the eight hues of Okabe and Ito's colour-blind-safe palette
#   (okabe2008color), as grDevices::palette.colors() ships them, with black
#   left out so that labels stay readable. Beyond eight categories the HCL
#   "Dark 3" palette of Zeileis et al. (zeileis2020colorspace) is used; it
#   stays distinguishable in sequence but is not guaranteed safe for every
#   kind of colour vision deficiency, which is why plots with many groups
#   should prefer a layout or facets over colour.
# * Numbers: viridis (vanderwalt2015viridis), the perceptually uniform,
#   colour-blind-safe sequential palette, light for small and dark for
#   large values.
# * Matrices and heat maps: the single-hue sequential "Blues 3" palette
#   from the same HCL family (zeileis2020colorspace), light for small and
#   dark for large cells, so that cells read as densities rather than as
#   categories.
# * Accents: Okabe-Ito blue for the main series and vermilion for observed
#   statistics and highlights; mid gray for reference series and light gray
#   for missing values.

# Colours for k categories.
plot_qualitative <- function(k) {
  if (k <= 0L) return(character())
  if (k <= 8L) {
    return(unname(grDevices::palette.colors(9L, "Okabe-Ito")[-1L])[seq_len(k)])
  }
  grDevices::hcl.colors(k, "Dark 3")
}

# n colours for numeric values, light to dark.
plot_sequential <- function(n) grDevices::hcl.colors(n, "viridis")

# n colours for matrix cells, light to dark.
plot_matrix_palette <- function(n) grDevices::hcl.colors(n, "Blues 3", rev = TRUE)

plot_blue <- "#0072B2"
plot_accent <- "#D55E00"
plot_gray <- "gray55"
plot_missing <- "gray80"

# ---- Layout of the page ----------------------------------------------------------

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

# Right margin (in lines) that fits the widest text of `legends`, each a
# list with `title` and `labels`.
legend_margin <- function(legends, cex) {
  if (!length(legends)) return(0)
  text <- unlist(lapply(legends, function(l) c(l$title, l$labels)))
  widths <- tryCatch(graphics::strwidth(text, units = "inches", cex = cex),
                     error = function(e) 0.09 * cex * nchar(text))
  (max(widths) + 0.45) / (graphics::par("csi") * graphics::par("mex"))
}

# Draws `legends` down the right margin of a network plot, starting at the
# top of the plot region. Each legend is a list with `type` ("fill",
# "size", or "line"), `title`, `labels`, and the colours, sizes, or widths
# of its keys.
draw_legends <- function(legends, usr, per_inch, cex) {
  left <- usr[2] + 0.15 * per_inch
  ytop <- usr[4]
  for (l in legends) {
    args <- list(x = left, y = ytop, legend = l$labels, title = l$title, bty = "n",
                 cex = cex, xjust = 0, yjust = 1, title.adj = 0)
    args <- switch(l$type,
      fill = c(args, list(pch = 21, pt.bg = l$col, col = "gray40", pt.cex = 1.6)),
      size = c(args, list(pch = 21, pt.bg = "gray70", col = "gray40", pt.cex = l$cex)),
      line = c(args, list(lty = 1, lwd = l$lwd, col = l$col))
    )
    box <- do.call(graphics::legend, args)
    ytop <- ytop - box$rect$h - 0.1 * per_inch
  }
  invisible(NULL)
}

# ---- Heat maps -------------------------------------------------------------------

# Draws a numeric matrix as a heat map with row 1 at the top, as a matrix
# is printed: rows are labelled on the left and columns at the bottom.
# `col` maps the range of the finite values from light to dark; cells that
# are NA or NaN are left blank. `lines` are 0-based positions after which a
# boundary is drawn (between rows/columns k and k + 1), and `values` writes
# each cell's value in it.
draw_heatmap <- function(m, row_labels = rownames(m), col_labels = colnames(m),
                         col = plot_matrix_palette(64), range = NULL, values = FALSE,
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

# ---- Node groupings for overlays --------------------------------------------------

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

# ---- Appearance arguments of network plots ---------------------------------------
# Shared by plot.netframe() and net_plot_spatial(). Each appearance argument
# can be the name of a column of the node or edge table, a vector with one
# value per node or edge, a single value for all, or a vector wrapped in
# I() to be used as it is.

# How an appearance argument was given: "default", "constant", "asis", or
# "mapped" (from a column or a vector), with its values.
plot_aes <- function(value, table, k, arg, what) {
  if (is.null(value) || k == 0) return(list(kind = "default"))
  if (inherits(value, "AsIs")) {
    value <- unclass(value)
    if (length(value) == 1L) value <- rep(value, k)
    if (length(value) != k) {
      stop(sprintf("`%s` wrapped in I() must have one value per %s.", arg,
                   sub("s$", "", what)), call. = FALSE)
    }
    return(list(kind = "asis", values = value))
  }
  if (is.character(value) && length(value) == 1L && value %in% names(table)) {
    return(list(kind = "mapped", values = table[[value]], title = value))
  }
  if (length(value) == k && k != 1L) {
    return(list(kind = "mapped", values = value, title = NULL))
  }
  if (length(value) == 1L) return(list(kind = "constant", values = value))
  stop(sprintf("`%s` must be a column of `x$%s`, one value, or one value per %s.", arg, what,
               sub("s$", "", what)), call. = FALSE)
}

# Colours for nodes or edges: `col` per element, whether they were mapped,
# and the legend for a mapping.
plot_color <- function(value, table, k, arg, what, default) {
  a <- plot_aes(value, table, k, arg, what)
  if (a$kind == "default") return(list(col = rep(default, k), mapped = FALSE))
  if (a$kind == "constant") {
    ok <- tryCatch({
      grDevices::col2rgb(a$values)
      TRUE
    }, error = function(e) FALSE)
    if (!ok) {
      stop(sprintf("`%s` is neither a column of `x$%s` nor a color.", arg, what),
           call. = FALSE)
    }
    return(list(col = rep(a$values, k), mapped = FALSE))
  }
  if (a$kind == "asis") return(list(col = a$values, mapped = FALSE))
  mapped <- map_colors(a$values)
  type <- if (what == "edges") "line" else "fill"
  list(col = mapped$col, mapped = TRUE,
       legend = list(type = type, title = a$title, labels = mapped$labels,
                     col = mapped$key, lwd = 2))
}

# Colours for a vector of values: numbers through the sequential palette,
# anything else through the qualitative palette in the order of the factor
# levels or of the sorted values; missing values light gray. Returns the
# colours with the legend's labels and keys.
map_colors <- function(v) {
  missing <- plot_missing
  if (is.numeric(v)) {
    ok <- is.finite(v)
    pal <- plot_sequential(100)
    if (!any(ok)) {
      return(list(col = rep(missing, length(v)), labels = "NA", key = missing))
    }
    rng <- range(v[ok])
    pick <- function(z) {
      if (diff(rng) > 0) pal[1L + round(99 * (z - rng[1]) / diff(rng))] else pal[50L]
    }
    col <- rep(missing, length(v))
    col[ok] <- pick(v[ok])
    breaks <- pretty(rng, 4)
    breaks <- breaks[breaks >= rng[1] & breaks <= rng[2]]
    if (length(breaks) < 2L) breaks <- unique(rng)
    labels <- format(breaks)
    key <- pick(breaks)
    if (!all(ok)) {
      labels <- c(labels, "NA")
      key <- c(key, missing)
    }
    return(list(col = col, labels = labels, key = key))
  }
  f <- if (is.factor(v)) droplevels(v) else factor(v, levels = sort(unique(v[!is.na(v)])))
  k <- nlevels(f)
  pal <- plot_qualitative(k)
  col <- pal[as.integer(f)]
  col[is.na(f)] <- missing
  labels <- levels(f)
  key <- pal
  if (k > 15) {
    labels <- c(labels[1:14], sprintf("(%d more)", k - 14))
    key <- c(key[1:14], NA)
  }
  if (anyNA(f)) {
    labels <- c(labels, "NA")
    key <- c(key, missing)
  }
  list(col = col, labels = labels, key = key)
}

# Linear map of numeric values onto [lo, hi]; with `area`, onto the square
# root so that areas grow linearly. Missing values get `lo`.
map_numeric <- function(v, lo, hi, area) {
  ok <- is.finite(v)
  out <- rep(lo, length(v))
  if (!any(ok)) return(list(values = out, rng = c(NA_real_, NA_real_)))
  rng <- range(v[ok])
  f <- function(z) {
    t <- if (diff(rng) > 0) (z - rng[1]) / diff(rng) else rep(0.5, length(z))
    lo + (hi - lo) * (if (area) sqrt(t) else t)
  }
  out[ok] <- f(v[ok])
  list(values = out, rng = rng, f = f)
}

numeric_key <- function(map) {
  breaks <- pretty(map$rng, 3)
  breaks <- breaks[breaks >= map$rng[1] & breaks <= map$rng[2]]
  if (length(breaks) < 2L) breaks <- unique(map$rng)
  breaks
}

# Base symbol size for n nodes: smaller in larger networks.
plot_base_size <- function(n) if (n > 0) min(2.5, max(0.25, 4 / n^0.25)) else 1

# Base line width for m edges.
plot_base_width <- function(m) if (m <= 200) 1 else if (m <= 5000) 0.7 else 0.5

# Node sizes as `cex`, by degree when nothing is given.
plot_node_size <- function(x, value, n, base = plot_base_size(n)) {
  a <- plot_aes(value, x$nodes, n, "node_size", "nodes")
  lo <- 0.6 * base
  hi <- 1.6 * base
  if (a$kind == "default") {
    d <- if (n > 0) net_degree(x, loops = FALSE) else numeric()
    return(list(cex = map_numeric(d, lo, hi, area = TRUE)$values))
  }
  if (!is.numeric(a$values)) {
    stop("`node_size` must be numeric.", call. = FALSE)
  }
  if (a$kind == "constant" || a$kind == "asis") {
    cex <- rep_len(a$values, n)
    if (any(!is.finite(cex) | cex < 0)) {
      stop("`node_size` values used as given must be finite and not negative.",
           call. = FALSE)
    }
    return(list(cex = cex))
  }
  map <- map_numeric(a$values, lo, hi, area = TRUE)
  out <- list(cex = map$values)
  if (!is.na(map$rng[1])) {
    breaks <- numeric_key(map)
    out$legend <- list(type = "size", title = a$title, labels = format(breaks),
                       cex = map$f(breaks))
  }
  out
}

# Edge widths as `lwd`: constant, except that strength weights widen edges
# when nothing is given.
plot_edge_width <- function(x, value, m, base = plot_base_width(m)) {
  a <- plot_aes(value, x$edges, m, "edge_width", "edges")
  if (a$kind == "default") {
    weight <- x$roles$weight
    if (!is.null(weight) && identical(x$roles$weight_type, "strength") && m > 0) {
      return(list(lwd = map_numeric(x$edges[[weight]], 0.5 * base, 3 * base,
                                    area = FALSE)$values))
    }
    return(list(lwd = rep(base, m)))
  }
  if (!is.numeric(a$values)) {
    stop("`edge_width` must be numeric.", call. = FALSE)
  }
  if (a$kind == "constant" || a$kind == "asis") {
    lwd <- rep_len(a$values, m)
    if (any(!is.finite(lwd) | lwd < 0)) {
      stop("`edge_width` values used as given must be finite and not negative.",
           call. = FALSE)
    }
    return(list(lwd = lwd))
  }
  map <- map_numeric(a$values, 0.4 * base, 3 * base, area = FALSE)
  out <- list(lwd = map$values)
  if (!is.na(map$rng[1])) {
    breaks <- numeric_key(map)
    out$legend <- list(type = "line", title = a$title, labels = format(breaks),
                       lwd = map$f(breaks), col = "gray35")
  }
  out
}

# Node labels: identifiers in small networks, none in large ones, or as
# given.
plot_labels <- function(x, value, n) {
  if (is.null(value)) return(if (n > 0 && n <= 40) node_ids(x) else rep("", n))
  if (isFALSE(value)) return(rep("", n))
  if (isTRUE(value)) return(node_ids(x))
  a <- plot_aes(value, x$nodes, n, "node_label", "nodes")
  labels <- as.character(rep_len(a$values, n))
  labels[is.na(labels)] <- ""
  labels
}

# Default node colors of a two-mode network: one per mode, labelled with
# the mode labels ("first" and "second" for a logical column).
twomode_fill <- function(x) {
  labels <- x$roles$mode_levels
  if (length(labels) != 2L || identical(sort(labels), c("FALSE", "TRUE"))) {
    labels <- c("first", "second")
  }
  second <- twomode_types(x)
  pal <- plot_qualitative(2L)
  list(col = pal[second + 1L], mapped = TRUE,
       legend = list(type = "fill", title = "mode", labels = labels, col = pal, lwd = 2))
}

# Node positions for plot.netframe(): a layout method, two node columns, a
# matrix or data frame, or a function of the network.
plot_layout <- function(x, layout, seed) {
  n <- nrow(x$nodes)
  methods <- c("auto", "stress", "sparse_stress", "fr", "circle", "star", "grid", "tree",
               "twomode")
  if (is.null(layout)) {
    return(net_layout(x, seed = seed))
  }
  if (is.character(layout) && length(layout) == 1L) {
    if (!layout %in% methods) {
      stop(sprintf("Unknown layout method `%s`; use one of %s.", layout,
                   paste(sprintf("\"%s\"", methods), collapse = ", ")), call. = FALSE)
    }
    return(net_layout(x, method = layout, seed = seed))
  }
  if (is.character(layout) && length(layout) == 2L) {
    missing <- setdiff(layout, names(x$nodes))
    if (length(missing)) {
      stop(sprintf("Layout column `%s` not found in the node table.", missing[1]),
           call. = FALSE)
    }
    cols <- lapply(layout, function(col) x$nodes[[col]])
    if (!all(vapply(cols, is.numeric, logical(1)))) {
      stop("Layout columns must be numeric.", call. = FALSE)
    }
    out <- layout_matrix(cbind(cols[[1]], cols[[2]]), n)
    dimnames(out) <- list(NULL, c("x", "y"))
    return(out)
  }
  if (is.function(layout)) layout <- layout(x)
  out <- layout_matrix(layout, n)
  dimnames(out) <- list(NULL, c("x", "y"))
  out
}

# Straight edges between node positions, shortened to the node borders and
# tipped with arrowheads when directed. `radius` is each node's symbol
# radius in user units and `per_inch` the user units per inch.
draw_edges <- function(xy, from, to, col, lwd, directed, radius, per_inch, n) {
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

# Radius of a pch 21 circle in user units: 0.375 of the symbol size (cex
# times the point size, in points).
symbol_radius <- function(cex, per_inch) {
  0.375 * cex * graphics::par("cex") * graphics::par("ps") / 72 * per_inch
}

# Draws node labels below the symbols.
draw_node_labels <- function(xy, labels, radius, n) {
  if (!any(nzchar(labels))) return(invisible(NULL))
  label_cex <- if (n <= 15) 0.9 else if (n <= 40) 0.75 else 0.6
  below <- radius + 0.3 * graphics::strheight("M", cex = label_cex)
  graphics::text(xy[, 1], xy[, 2] - below, labels, adj = c(0.5, 1), cex = label_cex,
                 col = "gray15")
  invisible(NULL)
}
