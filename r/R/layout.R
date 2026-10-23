#' Node coordinates for drawing a network
#'
#' `layout_network()` places the nodes of a network in the plane.
#' `layout_stress()` measures how well a layout represents the network's
#' shortest-path distances.
#'
#' @param x A [netframe()].
#' @param method The layout algorithm; see Details.
#' @param seed Seed for the random parts of a layout (a whole number in
#'   \[0, 2^53)). `NULL` draws one from R's random number generator, so
#'   [set.seed()] makes layouts reproducible.
#' @param weighted Whether to use edge weights (stress and force-directed
#'   layouts only). `NULL` uses them when the network has a weight column.
#' @param weights_as Override the network's `weight_type`: `"distance"` or
#'   `"strength"`.
#' @param types For `method = "twomode"`: which mode each node belongs to,
#'   as the name of a column of `x$nodes` or a vector with one value per
#'   node, with exactly two distinct values. `NULL` uses a column named
#'   `type`. Nodes of the first value (`FALSE`, or the first level or
#'   sorted value) form the top row.
#' @param root For `method = "tree"`: identifiers of root nodes, at most
#'   one per component. `NULL` uses each component's center.
#' @param center For `method = "star"`: the identifier of the center node.
#'   `NULL` uses the first node with the most neighbors.
#' @param max_iter Iteration cap for `"stress"`, `"sparse_stress"`, and
#'   `"fr"` (number of iterations for `"fr"`). `NULL` uses the defaults: 500
#'   for stress on all pairs, 100 for sparse stress, 500 for `"fr"`.
#' @param tol Convergence tolerance for stress layouts. `NULL` uses the
#'   defaults (see Details).
#' @param pivots Number of pivots for sparse stress (default 100) and for
#'   the starting layout of stress on all pairs (default 50). `NULL` uses
#'   the defaults.
#' @param approximation For `"fr"`: `"auto"` computes repulsion exactly for
#'   components of up to 1,000 nodes and with a Barnes-Hut quadtree beyond,
#'   `"exact"` always computes it exactly, and `"barnes_hut"` always
#'   approximates it.
#' @param threads Number of threads; `0` uses all available cores. Results
#'   do not depend on the number of threads.
#'
#' @details
#' Methods:
#' \describe{
#'   \item{`"auto"`}{Stress majorization (`"stress"`) for components of up
#'     to 2,000 nodes and sparse stress (`"sparse_stress"`) for larger
#'     ones, so the default stays fast on very large networks.}
#'   \item{`"stress"`}{Stress majorization (Gansner, Koren, and North
#'     2005), which places nodes so that their distances in the plane match
#'     their shortest-path distances as closely as possible. It minimizes
#'     stress, \eqn{\sum_{i<j} (\|x_i - x_j\| - d_{ij})^2 / d_{ij}^2}, the
#'     energy of Kamada and Kawai (1989), starting from pivot
#'     multidimensional scaling (Brandes and Pich 2007). It needs all
#'     shortest-path distances, so time per iteration and memory grow with
#'     the square of the component size; components of more than 10,000
#'     nodes are refused. Iteration stops when stress falls by less than a
#'     relative 1e-4 per iteration. Components of at most 100 nodes are
#'     also laid out from four random starts, keeping the layout with the
#'     least stress, with a tolerance of 1e-5.}
#'   \item{`"sparse_stress"`}{The sparse stress model of Ortmann, Klimenta,
#'     and Brandes (2017): stress over edges and over distances to a random
#'     sample of pivot nodes, each pivot weighted by the part of the network
#'     it represents. Time per iteration grows linearly with the network.
#'     It approximates the stress objective; it never changes any measure.}
#'   \item{`"fr"`}{The force-directed algorithm of Fruchterman and Reingold
#'     (1991), from random starting positions. Repulsion between all pairs
#'     is computed exactly or, for large components, with the Barnes and Hut
#'     (1986) quadtree approximation (see `approximation`), which only
#'     changes the simulated forces. With weights, an edge's attraction is
#'     multiplied by its strength (the weight, or 1/weight for distances).}
#'   \item{`"circle"`}{Nodes on a circle in the order of `x$nodes`.}
#'   \item{`"star"`}{One node in the center, the others on a circle.}
#'   \item{`"grid"`}{Nodes on a square grid in the order of `x$nodes`, row
#'     by row from the top.}
#'   \item{`"tree"`}{A layered drawing of a breadth-first spanning tree of
#'     each component: depth increases downwards, leaves are one unit
#'     apart, and each parent is centered over its children. The default
#'     root is a center of the component (the middle of the path between
#'     the ends of a double breadth-first sweep).}
#'   \item{`"twomode"`}{Two rows for a two-mode (bipartite) network, given
#'     by `types`, ordered by the barycenter heuristic of Sugiyama, Tagawa,
#'     and Toda (1981) to reduce edge crossings.}
#' }
#'
#' All methods ignore edge direction and self-loops and count parallel
#' edges once. In stress layouts, weights are lengths: distances as given,
#' or strengths converted to 1/weight; with several parallel edges the
#' shortest counts. Weights must be positive.
#'
#' Stress, force-directed, and tree layouts treat each weak component
#' separately and pack the components in rows, largest first, separated by
#' one typical edge length; isolated nodes end up in rows of their own.
#' Each stress or force-directed component is centered and rotated so its
#' longest extent is horizontal. Coordinates are in units of path length
#' for stress layouts (hops, without weights) and of the ideal edge length
#' for `"fr"`.
#'
#' Layouts are identical for a given `seed` on any number of threads. Stress
#' layouts use the seed only for a tiny perturbation of the starting
#' layout, the random restarts of small components, and the pivots of
#' sparse stress.
#'
#' `layout_stress()` returns the stress of a layout, \eqn{\sum_{i<j} (\|x_i
#' - x_j\| - d_{ij})^2 / d_{ij}^2} over pairs of nodes in the same
#' component, with shortest-path lengths \eqn{d_{ij}} computed as for
#' stress layouts. It depends on the scale of the layout: compare layouts
#' at the scale that minimizes it, or use layouts in path-length units.
#'
#' Stress layouts of the validation networks have stress no higher than
#' NetworkX's `kamada_kawai_layout()` (within 2\%, and usually lower); see
#' `validation/layout_check.py` in the source repository.
#'
#' @references
#' Gansner, E. R., Koren, Y., & North, S. (2005). Graph drawing by stress
#' majorization. In *Graph Drawing (GD 2004)*, Lecture Notes in Computer
#' Science 3383, 239-250. \doi{10.1007/978-3-540-31843-9_25}
#'
#' Kamada, T., & Kawai, S. (1989). An algorithm for drawing general
#' undirected graphs. *Information Processing Letters*, 31(1), 7-15.
#' \doi{10.1016/0020-0190(89)90102-6}
#'
#' Brandes, U., & Pich, C. (2007). Eigensolver methods for progressive
#' multidimensional scaling of large data. In *Graph Drawing (GD 2006)*,
#' Lecture Notes in Computer Science 4372, 42-53.
#' \doi{10.1007/978-3-540-70904-6_6}
#'
#' Ortmann, M., Klimenta, M., & Brandes, U. (2017). A sparse stress model.
#' *Journal of Graph Algorithms and Applications*, 21(5), 791-821.
#' \doi{10.7155/jgaa.00440}
#'
#' Fruchterman, T. M. J., & Reingold, E. M. (1991). Graph drawing by
#' force-directed placement. *Software: Practice and Experience*, 21(11),
#' 1129-1164. \doi{10.1002/spe.4380211102}
#'
#' Barnes, J., & Hut, P. (1986). A hierarchical O(N log N) force-calculation
#' algorithm. *Nature*, 324, 446-449. \doi{10.1038/324446a0}
#'
#' Sugiyama, K., Tagawa, S., & Toda, M. (1981). Methods for visual
#' understanding of hierarchical system structures. *IEEE Transactions on
#' Systems, Man, and Cybernetics*, 11(2), 109-125.
#' \doi{10.1109/TSMC.1981.4308636}
#'
#' @return `layout_network()` returns a numeric matrix with columns `x` and
#'   `y` and one row per node, aligned with the rows of `x$nodes`.
#'   `layout_stress()` returns a number.
#' @seealso [plot.netframe()]
#' @export
layout_network <- function(x, method = c("auto", "stress", "sparse_stress", "fr", "circle",
                                         "star", "grid", "tree", "twomode"),
                           seed = NULL, weighted = NULL, weights_as = NULL, types = NULL,
                           root = NULL, center = NULL, max_iter = NULL, tol = NULL,
                           pivots = NULL, approximation = c("auto", "exact", "barnes_hut"),
                           threads = 0L) {
  check_netframe(x)
  method <- match.arg(method)
  approximation <- match.arg(approximation)
  ptr <- engine_graph(x)
  xy <- switch(method,
    auto = , stress = , sparse_stress = {
      w <- resolve_weighting(x, weighted, weights_as)
      algorithm <- switch(method, auto = 0L, stress = 1L, sparse_stress = 2L)
      .Call(C_layout_stress, ptr, algorithm, w$weighted, w$code,
            layout_count(pivots, "pivots"), layout_count(max_iter, "max_iter"),
            layout_tol(tol), layout_seed(seed), as.integer(threads))
    },
    fr = {
      w <- resolve_weighting(x, weighted, weights_as)
      code <- switch(approximation, auto = 0L, exact = 1L, barnes_hut = 2L)
      .Call(C_layout_fr, ptr, w$weighted, w$code, layout_count(max_iter, "max_iter"), code,
            0, layout_seed(seed), as.integer(threads))
    },
    circle = .Call(C_layout_simple, ptr, 0L, 0),
    star = {
      param <- if (is.null(center)) -1 else node_index(x, center, "center", one = TRUE) - 1
      .Call(C_layout_simple, ptr, 1L, as.double(param))
    },
    grid = .Call(C_layout_simple, ptr, 2L, 0),
    tree = {
      roots <- if (is.null(root)) double() else node_index(x, root, "root") - 1
      .Call(C_layout_tree, ptr, as.double(roots))
    },
    twomode = .Call(C_layout_twomode, ptr, twomode_types(x, types), 0)
  )
  n <- nrow(x$nodes)
  matrix(xy, nrow = n, ncol = 2L, dimnames = list(NULL, c("x", "y")))
}

#' @rdname layout_network
#' @param layout A layout of `x`: a numeric matrix or data frame with two
#'   columns and one row per node.
#' @export
layout_stress <- function(x, layout, weighted = NULL, weights_as = NULL, threads = 0L) {
  check_netframe(x)
  xy <- layout_matrix(layout, nrow(x$nodes))
  w <- resolve_weighting(x, weighted, weights_as)
  .Call(C_layout_stress_value, engine_graph(x), as.double(xy), w$weighted, w$code,
        as.integer(threads))
}

layout_seed <- function(seed) {
  if (is.null(seed)) {
    # Two draws give 52 random bits, so set.seed() controls the result.
    return(floor(stats::runif(1L) * 2^26) * 2^26 + floor(stats::runif(1L) * 2^26))
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed < 0 ||
      seed >= 2^53 || seed != floor(seed)) {
    stop("`seed` must be a whole number in [0, 2^53).", call. = FALSE)
  }
  as.double(seed)
}

layout_count <- function(value, what) {
  if (is.null(value)) return(0)
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value < 1 ||
      value != floor(value)) {
    stop(sprintf("`%s` must be a positive whole number.", what), call. = FALSE)
  }
  as.double(value)
}

layout_tol <- function(tol) {
  if (is.null(tol)) return(0)
  if (!is.numeric(tol) || length(tol) != 1L || !isTRUE(tol > 0)) {
    stop("`tol` must be a positive number.", call. = FALSE)
  }
  as.double(tol)
}

# 1-based row indices of node identifiers.
node_index <- function(x, ids, what, one = FALSE) {
  if (one && length(ids) != 1L) {
    stop(sprintf("`%s` must be a single node identifier.", what), call. = FALSE)
  }
  idx <- match(ids, x$nodes[[x$roles$id]])
  if (anyNA(idx)) {
    stop(sprintf("`%s` contains identifiers that are not nodes of `x`: %s.", what,
                 paste(format(utils::head(ids[is.na(idx)], 5L)), collapse = ", ")),
         call. = FALSE)
  }
  idx
}

# Node types for a two-mode layout as 0 (top row) and 1 (bottom row).
twomode_types <- function(x, types) {
  n <- nrow(x$nodes)
  if (is.null(types)) {
    if (!"type" %in% names(x$nodes)) {
      stop("A two-mode layout needs node types: supply `types =` or a `type` column ",
           "in the node table.", call. = FALSE)
    }
    types <- x$nodes$type
  } else if (is.character(types) && length(types) == 1L && types %in% names(x$nodes)) {
    types <- x$nodes[[types]]
  } else if (is.character(types) && length(types) == 1L && n != 1L) {
    stop(sprintf("Column `%s` not found in the node table.", types), call. = FALSE)
  }
  if (length(types) != n) {
    stop("`types` must have one value per node.", call. = FALSE)
  }
  if (anyNA(types)) stop("`types` must not contain missing values.", call. = FALSE)
  values <- if (is.factor(types)) levels(droplevels(types)) else sort(unique(types))
  if (length(values) > 2L) {
    stop("`types` must have exactly two distinct values for a two-mode layout.",
         call. = FALSE)
  }
  as.integer(match(as.character(types), as.character(values)) - 1L)
}

# A layout given as a matrix or data frame, checked against the node count.
layout_matrix <- function(layout, n) {
  if (is.data.frame(layout)) {
    cols <- if (all(c("x", "y") %in% names(layout))) c("x", "y") else seq_len(min(2L, ncol(layout)))
    layout <- as.matrix(layout[, cols, drop = FALSE])
  }
  if (!is.matrix(layout) || !is.numeric(layout) || ncol(layout) != 2L ||
      nrow(layout) != n) {
    stop(sprintf("A layout must be a numeric matrix or data frame with 2 columns and %d rows.",
                 n), call. = FALSE)
  }
  if (!all(is.finite(layout))) {
    stop("Layout coordinates must be finite.", call. = FALSE)
  }
  storage.mode(layout) <- "double"
  layout
}
