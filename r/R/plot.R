#' Plot a network
#'
#' Draws a network with base R graphics. With no arguments beyond the
#' network, nodes are placed by [layout_network()], sized by degree, and
#' labelled when the network is small; edges are semi-transparent, with
#' arrows when the network is directed. Node and edge appearance can be set
#' from columns of `x$nodes` and `x$edges`.
#'
#' @param x A [netframe()].
#' @param layout Node positions: `NULL` for [layout_network()]'s default,
#'   the name of a [layout_network()] method (such as `"fr"` or `"circle"`),
#'   the names of two numeric columns of `x$nodes` (for example
#'   `c("lon", "lat")`), a matrix or data frame with two columns and one row
#'   per node (as returned by [layout_network()]), or a function that takes
#'   the network and returns such a matrix.
#' @param node_size,node_color,node_label,edge_width,edge_color Appearance
#'   of nodes and edges. Each can be the name of a column of `x$nodes` (for
#'   node arguments) or `x$edges` (for edge arguments), a vector with one
#'   value per node or edge, or a single value for all. Columns and vectors
#'   are mapped to sizes, widths, or colors (see Details) and get a legend;
#'   wrap a vector in [I()] to use its values as they are (sizes as `cex`
#'   multipliers, widths as `lwd`, colors as R colors). `NULL` uses the
#'   defaults (see Details). `node_label = TRUE` labels nodes with their
#'   identifiers and `FALSE` draws no labels.
#' @param ... Further arguments passed to [graphics::title()], such as
#'   `main`.
#' @param legend Whether to draw legends for mapped columns and vectors.
#' @param seed Seed for the layout, as in [layout_network()], when `layout`
#'   is `NULL` or a method name.
#'
#' @details
#' Defaults:
#' \itemize{
#'   \item Node size grows with degree (as area), and nodes are smaller in
#'     larger networks.
#'   \item Nodes are a single blue with a white border, except in two-mode
#'     networks ([twomode()]), where each mode has its own color.
#'   \item Nodes are labelled with their identifiers in networks of at most
#'     40 nodes.
#'   \item Edges are gray and more transparent the more edges there are.
#'     In directed networks, arrows end at the border of the target node.
#'   \item Edge width is constant, except in networks whose weights are
#'     strengths (`weight_type = "strength"` in [netframe()]), where it
#'     grows with the weight. Distance weights do not change the width,
#'     since a thick line would suggest a strong tie.
#' }
#'
#' Mapping: numeric values give sizes and widths that grow with the value
#' (node area grows linearly) and colors from the sequential palette
#' `grDevices::hcl.colors(, "viridis")`. Factors, character vectors, and
#' logical vectors give colors from the qualitative palette
#' `grDevices::hcl.colors(, "Dark 3")`, in the order of the factor levels
#' or of the sorted values. Missing values are drawn in light gray (colors)
#' or at the smallest size or width. Sizes and widths must be numeric.
#'
#' Legends are drawn in the right margin. Coordinates keep their aspect
#' ratio, so distances in the plot are proportional to distances in the
#' layout. Self-loops are not drawn, and parallel edges are drawn on top of
#' each other.
#'
#' The plot is reproducible for a given `seed` (or after [set.seed()]).
#'
#' @return The layout used, invisibly: a matrix with columns `x` and `y`
#'   and one row per node, which can be passed back as `layout`.
#' @seealso [layout_network()]
#' @examples
#' g <- netframe(data.frame(from  = c("a", "a", "b", "c", "c", "d"),
#'                          to    = c("b", "c", "c", "d", "e", "e"),
#'                          calls = c(1, 4, 2, 8, 3, 5)),
#'               directed = FALSE)
#' g$nodes$group <- c("x", "x", "y", "y", "y")
#' plot(g, seed = 1)
#' plot(g, node_color = "group", edge_width = "calls", seed = 1, main = "Calls")
#' xy <- plot(g, layout = "circle", node_label = FALSE)
#' plot(g, layout = xy, node_size = I(c(3, 1, 1, 1, 1)))
#'
#' attended <- data.frame(person = c("ann", "ann", "bob", "cat"),
#'                        event  = c("gala", "fair", "gala", "fair"))
#' plot(twomode(attended), layout = "twomode")
#' @export
plot.netframe <- function(x, layout = NULL, node_size = NULL, node_color = NULL,
                          node_label = NULL, edge_width = NULL, edge_color = NULL, ...,
                          legend = TRUE, seed = NULL) {
  check_netframe(x)
  n <- nrow(x$nodes)
  m <- nrow(x$edges)
  xy <- plot_layout(x, layout, seed)

  base <- if (n > 0) min(2.5, max(0.25, 4 / n^0.25)) else 1
  size <- plot_node_size(x, node_size, n, base)
  fill <- if (is.null(node_color) && n > 0 && is_twomode(x)) {
    twomode_fill(x)
  } else {
    plot_color(node_color, x$nodes, n, "node_color", "nodes", "#2C6FAC")
  }
  labels <- plot_labels(x, node_label, n)
  lwd_base <- if (m <= 200) 1 else if (m <= 5000) 0.7 else 0.5
  width <- plot_edge_width(x, edge_width, m, lwd_base)
  alpha <- min(0.6, max(0.05, 10 / sqrt(max(m, 1))))
  ecol <- plot_color(edge_color, x$edges, m, "edge_color", "edges",
                     grDevices::adjustcolor("gray35", alpha.f = alpha))
  if (ecol$mapped) ecol$col <- grDevices::adjustcolor(ecol$col, alpha.f = max(alpha, 0.7))

  legends <- if (isTRUE(legend)) {
    Filter(Negate(is.null), list(fill$legend, size$legend, ecol$legend, width$legend))
  } else {
    list()
  }
  dots <- list(...)
  top <- if (!is.null(dots$main)) 2.2 else 0.4
  right <- 0.4
  legend_cex <- 0.75
  if (length(legends)) {
    text <- unlist(lapply(legends, function(l) c(l$title, l$labels)))
    widths <- tryCatch(graphics::strwidth(text, units = "inches", cex = legend_cex),
                       error = function(e) 0.09 * legend_cex * nchar(text))
    inches <- max(widths) + 0.45
    right <- right + inches / (graphics::par("csi") * graphics::par("mex"))
  }
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
  # Radius of a pch 21 circle: 0.375 of the symbol size (cex times the
  # point size, in points).
  radius <- 0.375 * size$cex * graphics::par("cex") * graphics::par("ps") / 72 * per_inch

  if (m > 0) {
    from <- x$index$src + 1L
    to <- x$index$tgt + 1L
    keep <- from != to
    from <- from[keep]
    to <- to[keep]
    col <- ecol$col[keep]
    lwd <- width$lwd[keep]
    x0 <- xy[from, 1]
    y0 <- xy[from, 2]
    x1 <- xy[to, 1]
    y1 <- xy[to, 2]
    if (x$directed) {
      dx <- x1 - x0
      dy <- y1 - y0
      len <- sqrt(dx^2 + dy^2)
      gap <- 0.01 * per_inch
      start <- radius[from]
      end <- radius[to] + gap
      ok <- len > start + end
      ux <- dx[ok] / len[ok]
      uy <- dy[ok] / len[ok]
      head <- if (n <= 100) 0.08 else 0.05
      if (any(ok)) {
        graphics::arrows(x0[ok] + ux * start[ok], y0[ok] + uy * start[ok],
                         x1[ok] - ux * end[ok], y1[ok] - uy * end[ok],
                         length = head, angle = 22, col = col[ok], lwd = lwd[ok])
      }
    } else if (length(from)) {
      graphics::segments(x0, y0, x1, y1, col = col, lwd = lwd)
    }
  }

  border <- if (n <= 2000) "white" else NA
  graphics::points(xy[, 1], xy[, 2], pch = 21, bg = fill$col, col = border,
                   cex = size$cex, lwd = 0.6)
  if (any(nzchar(labels))) {
    label_cex <- if (n <= 15) 0.9 else if (n <= 40) 0.75 else 0.6
    below <- radius + 0.3 * graphics::strheight("M", cex = label_cex)
    graphics::text(xy[, 1], xy[, 2] - below, labels, adj = c(0.5, 1), cex = label_cex,
                   col = "gray15")
  }
  if (length(dots)) do.call(graphics::title, dots)

  if (length(legends)) {
    left <- usr[2] + 0.15 * per_inch
    ytop <- usr[4]
    for (l in legends) {
      args <- list(x = left, y = ytop, legend = l$labels, title = l$title, bty = "n",
                   cex = legend_cex, xjust = 0, yjust = 1, title.adj = 0)
      args <- switch(l$type,
        fill = c(args, list(pch = 21, pt.bg = l$col, col = "gray40", pt.cex = 1.6)),
        size = c(args, list(pch = 21, pt.bg = "gray70", col = "gray40", pt.cex = l$cex)),
        line = c(args, list(lty = 1, lwd = l$lwd, col = l$col))
      )
      box <- do.call(graphics::legend, args)
      ytop <- ytop - box$rect$h - 0.1 * per_inch
    }
  }
  invisible(xy)
}

# Default node colors of a two-mode network: one per mode, labelled with
# the two values recorded by twomode() or as_twomode() ("first" and
# "second" for a logical type column).
twomode_fill <- function(x) {
  labels <- x$roles$type_labels
  if (length(labels) != 2L || identical(labels, c("FALSE", "TRUE"))) {
    labels <- c("first", "second")
  }
  second <- x$nodes[[x$roles$type]]
  pal <- grDevices::hcl.colors(2L, "Dark 3")
  list(col = pal[second + 1L], mapped = TRUE,
       legend = list(type = "fill", title = "mode", labels = labels, col = pal, lwd = 2))
}

plot_layout <- function(x, layout, seed) {
  n <- nrow(x$nodes)
  methods <- c("auto", "stress", "sparse_stress", "fr", "circle", "star", "grid", "tree",
               "twomode")
  if (is.null(layout)) {
    return(layout_network(x, seed = seed))
  }
  if (is.character(layout) && length(layout) == 1L) {
    if (!layout %in% methods) {
      stop(sprintf("Unknown layout method `%s`; use one of %s.", layout,
                   paste(sprintf("\"%s\"", methods), collapse = ", ")), call. = FALSE)
    }
    return(layout_network(x, method = layout, seed = seed))
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

map_colors <- function(v) {
  missing <- "gray80"
  if (is.numeric(v)) {
    ok <- is.finite(v)
    pal <- grDevices::hcl.colors(100, "viridis")
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
  pal <- if (k > 0) grDevices::hcl.colors(k, "Dark 3") else character()
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

plot_node_size <- function(x, value, n, base) {
  a <- plot_aes(value, x$nodes, n, "node_size", "nodes")
  lo <- 0.6 * base
  hi <- 1.6 * base
  if (a$kind == "default") {
    d <- if (n > 0) degree(x, loops = FALSE) else numeric()
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

plot_edge_width <- function(x, value, m, base) {
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

plot_labels <- function(x, value, n) {
  ids <- function() as.character(x$nodes[[x$roles$id]])
  if (is.null(value)) return(if (n > 0 && n <= 40) ids() else rep("", n))
  if (isFALSE(value)) return(rep("", n))
  if (isTRUE(value)) return(ids())
  a <- plot_aes(value, x$nodes, n, "node_label", "nodes")
  labels <- as.character(rep_len(a$values, n))
  labels[is.na(labels)] <- ""
  labels
}
