#' Plot an adjacency or incidence matrix
#'
#' Draws a network as its adjacency matrix (or, for a two-mode network, its
#' incidence matrix) with one shaded cell per tie, the way many network
#' texts present blockmodels and permuted matrices. Rows and columns can be
#' sorted by a partition of the nodes, with lines between the groups.
#'
#' @param x A [netframe()], or a numeric matrix.
#' @param order How to arrange rows and columns: `NULL` keeps the order of
#'   `x$nodes`; a partition (the name of a column of `x$nodes`, or a vector
#'   with one value per node, such as the result of [cluster_louvain()] or
#'   [concor()]) sorts the nodes by group, in order of the groups' first
#'   appearance, and draws a line between groups. For a matrix, a vector
#'   with one value per row (and column, if square).
#' @param attr For a netframe, an edge column whose values fill the matrix;
#'   `NULL` counts edges (see [as.matrix.netframe()]).
#' @param values Whether to write each cell's value in it.
#' @param key Whether to draw a color key.
#' @param labels Whether to label rows and columns with node identifiers.
#'   `NULL` labels them when there are at most 80.
#' @param ... Further arguments passed to [graphics::image()], such as
#'   `main`.
#'
#' @details
#' Cells are shaded from light (0) to dark (the largest value) with a
#' single-hue sequential palette, so a binary network shows ties as dark
#' cells on a light background. In a directed network row \eqn{i}, column
#' \eqn{j} is the tie from \eqn{i} to \eqn{j}; undirected networks give a
#' symmetric matrix. Parallel edges add up and a self-loop fills the
#' diagonal. A two-mode network ([twomode()]) is drawn as its
#' [incidence_matrix()], first mode in rows, with the partition applied to
#' each mode's nodes.
#'
#' @return The matrix as drawn (after any reordering), invisibly.
#' @seealso [plot.blockmodel()] for block densities, [plot_mixing_matrix()]
#'   for ties between categories.
#' @examples
#' g <- example_network("karate")
#' plot_matrix(g, main = "Karate club")
#' plot_matrix(g, order = "faction", main = "Sorted by faction")
#' plot_matrix(g, order = cluster_louvain(g, seed = 1))
#'
#' women <- example_network("southern_women")
#' plot_matrix(women, main = "Southern women: incidence matrix")
#' @export
plot_matrix <- function(x, order = NULL, attr = NULL, values = FALSE, key = NULL,
                        labels = NULL, ...) {
  values <- check_flag(values, "values")
  check_flag_or_null(key, "key")
  check_flag_or_null(labels, "labels")
  if (inherits(x, "netframe")) {
    check_netframe(x)
    if (is_twomode(x)) {
      m <- incidence_matrix(x, attr = attr)
      row_lines <- col_lines <- integer()
      if (!is.null(order)) {
        g <- plot_groups(x, order, "order")
        second <- x$nodes[[x$roles$type]]
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

#' Plot a mixing matrix
#'
#' Draws the ties between categories of nodes as a heat map with the counts
#' written in the cells: the mixing matrix of [mixing_matrix()], or any
#' square or rectangular table of counts.
#'
#' @param x A [netframe()], or a numeric matrix as returned by
#'   [mixing_matrix()].
#' @param dim1,dim2,direction For a netframe, passed to [mixing_matrix()].
#' @param values Whether to write the counts in the cells.
#' @param ... Further arguments passed to [graphics::image()], such as
#'   `main`.
#'
#' @details
#' Darker cells hold more ties. For an undirected network with one
#' attribute the matrix is folded into its upper triangle (see
#' [mixing_matrix()]), so the lower triangle is 0; the diagonal holds the
#' ties within each category. Compare the diagonal with the rest to judge
#' homophily, or use [assortativity_nominal()] and [ei_index()] for one
#' number.
#'
#' @return The matrix drawn, invisibly.
#' @seealso [mixing_matrix()], [plot_matrix()].
#' @examples
#' g <- example_network("karate")
#' plot_mixing_matrix(g, "faction", main = "Ties by faction")
#' m <- mixing_matrix(g, "faction")
#' plot_mixing_matrix(m)
#' @export
plot_mixing_matrix <- function(x, dim1 = NULL, dim2 = NULL, direction = c("all", "out", "in"),
                               values = TRUE, ...) {
  values <- check_flag(values, "values")
  if (inherits(x, "netframe")) {
    if (is.null(dim1)) stop("`dim1` is needed to compute the mixing matrix.", call. = FALSE)
    m <- mixing_matrix(x, dim1, dim2, direction = match.arg(direction))
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

#' Plot a blockmodel
#'
#' Draws the block densities of a [blockmodel()] as a heat map, with the
#' 1-blocks of the image matrix outlined, and, when the network is given,
#' its adjacency matrix permuted so that nodes of the same block are
#' adjacent, with lines at the block boundaries: the two displays used to
#' present blockmodels since White, Boorman, and Breiger (1976).
#'
#' @param x A [blockmodel()].
#' @param network The [netframe()] the blockmodel was computed on, for the
#'   permuted adjacency matrix. `NULL` draws only the densities.
#' @param which `"both"` (when `network` is given), `"density"`, or
#'   `"matrix"`.
#' @param values Whether to write each density in its cell.
#' @param ... Further arguments passed to [graphics::image()] for both
#'   panels, such as `main` (which replaces the panel titles) or
#'   `cex_labels`.
#'
#' @details
#' Densities are shaded from light (0) to dark (1, or the largest mean tie
#' value for a weighted blockmodel); blocks whose density reaches the
#' cutoff (the 1-blocks of `x$image`) are outlined. Blocks with no pairs
#' (a one-node block without loops) are blank. With both panels, the
#' permuted adjacency matrix is drawn on the left and the densities on the
#' right.
#'
#' @references
#' White, H. C., Boorman, S. A., & Breiger, R. L. (1976). Social structure
#' from multiple networks. I. Blockmodels of roles and positions.
#' *American Journal of Sociology*, 81(4), 730-780. \doi{10.1086/226141}
#'
#' @return `x`, invisibly.
#' @seealso [blockmodel()], [plot_matrix()].
#' @examples
#' g <- example_network("florentine")
#' bm <- blockmodel(g, concor(g, depth = 2))
#' plot(bm)
#' plot(bm, network = g)
#' @export
plot.blockmodel <- function(x, network = NULL, which = c("both", "density", "matrix"),
                            values = TRUE, ...) {
  which <- match.arg(which)
  values <- check_flag(values, "values")
  if (!is.null(network)) {
    check_netframe(network)
    if (nrow(network$nodes) != length(x$membership)) {
      stop("`network` must have one node per element of `x$membership`.", call. = FALSE)
    }
  } else if (which != "density") {
    if (which == "matrix") {
      stop("The permuted adjacency matrix needs `network`.", call. = FALSE)
    }
    which <- "density"
  }
  if (which == "both") {
    old <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(old), add = TRUE)
  }
  dots <- list(...)
  if (which != "density") {
    args <- c(list(x = network, order = x$membership, key = FALSE), dots)
    if (is.null(args$main)) args$main <- "Permuted adjacency matrix"
    do.call(plot_matrix, args)
  }
  if (which != "matrix") {
    d <- x$density
    finite <- d[is.finite(d)]
    top <- if (length(finite)) max(1, finite) else 1
    args <- c(list(m = d, values = values, range = c(0, top), digits = 2L), dots)
    if (is.null(args$main)) args$main <- "Block densities"
    do.call(draw_heatmap, args)
    k <- nrow(d)
    ones <- which(is.finite(x$image) & x$image == 1, arr.ind = TRUE)
    if (nrow(ones)) {
      graphics::rect(ones[, 2] - 0.5, k - ones[, 1] + 0.5, ones[, 2] + 0.5,
                     k - ones[, 1] + 1.5, border = plot_accent, lwd = 2.5)
    }
  }
  invisible(x)
}
