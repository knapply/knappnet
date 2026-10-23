#' Regular equivalence (CATREGE)
#'
#' Nodes are regularly equivalent when they have the same kinds of ties to
#' the same kinds of nodes: for every node tied to one of them, some node
#' of the same class is tied in the same way to the other (White and Reitz
#' 1983). Unlike structural equivalence, the alters need not be the same
#' nodes. The largest such equivalence, the maximal regular equivalence,
#' is found by CATREGE (Borgatti and Everett 1993): start with every node
#' in one class (or with `seed`), and repeatedly split each class by the
#' set of (type of tie, class of alter) pairs its members have, until no
#' class splits.
#'
#' @inheritParams betweenness
#' @param seed A partition to refine: a vector aligned with the rows of
#'   `x$nodes`, or the name of a column of `x$nodes`. Nodes in different
#'   seed classes are never equivalent; only equality of values matters.
#'   `NULL` starts from a single class.
#' @param weighted Whether to treat total edge weights as categories of
#'   tie: nodes with ties of different weights are then not equivalent.
#'   `NULL` uses weights when the network has a weight column. Weights are
#'   categories, not magnitudes, so `weight_type` does not matter.
#'
#' @details
#' The type of tie between i and k is the pair (tie from i to k, tie from
#' k to i), so the network and its transpose are treated as one multiplex
#' relation and a mutual tie differs from a one-way tie, as in sna's
#' `redist(method = "catrege")` and UCINET's CATREGE. In an undirected
#' network every tie has the same type. Self-loops are ignored and
#' parallel edges count once (with weights, their weights add up).
#'
#' In an undirected network without a seed partition the maximal regular
#' equivalence is degenerate: every node with at least one tie ends in one
#' class and isolated nodes in another, as Borgatti and Everett note.
#' Classes are numbered from 1 in order of their first node.
#'
#' The result matches the maximal regular equivalence computed by sna's
#' `redist(method = "catrege")` with `seed.partition` and by UCINET's
#' CATREGE; sna's `equiv.clust()` then cuts the resulting distance matrix.
#'
#' @references
#' White, D. R., & Reitz, K. P. (1983). Graph and semigroup homomorphisms
#' on networks of relations. *Social Networks*, 5(2), 193-234.
#' \doi{10.1016/0378-8733(83)90025-4}
#'
#' Borgatti, S. P., & Everett, M. G. (1993). Two algorithms for computing
#' regular equivalence. *Social Networks*, 15(4), 361-376.
#' \doi{10.1016/0378-8733(93)90012-A}
#'
#' Everett, M. G., & Borgatti, S. P. (1994). Regular equivalence: General
#' theory. *Journal of Mathematical Sociology*, 19(1), 29-52.
#' \doi{10.1080/0022250X.1994.9990134}
#'
#' @return An integer vector of class numbers (starting at 1) aligned with
#'   the rows of `x$nodes`, with attribute `rounds`, the number of
#'   refinement rounds that split a class.
#' @seealso [structural_equivalence()], [automorphic_equivalence()],
#'   [blockmodel()].
#' @examples
#' # A two-level hierarchy: 1 commands 2 and 3, who each command two others.
#' g <- netframe(data.frame(from = c(1, 1, 2, 2, 3, 3), to = c(2, 3, 4, 5, 6, 7)))
#' regular_equivalence(g)  # three classes: top, middle, bottom
#' @export
regular_equivalence <- function(x, seed = NULL, weighted = NULL) {
  check_netframe(x)
  weighted <- positions_weighted(x, weighted)
  codes <- if (is.null(seed)) NULL else node_groups(x, seed, "seed")$codes
  ptr <- engine_graph(x)
  r <- engine_call("regular_equivalence", .Call(C_regular_equivalence, ptr, codes, weighted))
  structure(r[[1L]] + 1L, rounds = r[[2L]])
}

#' Automorphic equivalence
#'
#' Nodes are automorphically equivalent when some automorphism of the
#' network, a relabelling of the nodes that maps the network onto itself,
#' maps one onto the other (Everett 1985; Borgatti and Everett 1992): they
#' are interchangeable, and have the same value on every measure that
#' depends on the structure alone. The classes are the orbits of the
#' automorphism group, computed exactly by an individualization-refinement
#' search (McKay and Piperno 2014).
#'
#' @inheritParams betweenness
#' @param max_steps A cap on the number of refined colourings the search
#'   examines. The search is exponential in the worst case (highly
#'   symmetric networks such as strongly regular graphs), so it stops with
#'   an error, rather than run indefinitely, once it has examined this many;
#'   results are never approximate.
#'
#' @details
#' The equivalence is defined on binary ties: weights are ignored and
#' parallel edges count once. Direction is respected: an automorphism must
#' preserve every arc. A self-loop is a property of its node, so a node
#' with a loop can only map to a node with a loop. Classes are numbered
#' from 1 in order of their first node.
#'
#' Automorphic equivalence implies regular equivalence and is implied by
#' structural equivalence. Networks of a few hundred nodes are usually
#' quick; igraph and NetworkX have no function for it, and UCINET's
#' routine (MAXSIM) computes an approximation based on distance profiles
#' rather than the exact orbits.
#'
#' @references
#' Everett, M. G. (1985). Role similarity and complexity in social
#' networks. *Social Networks*, 7(4), 353-359.
#' \doi{10.1016/0378-8733(85)90013-9}
#'
#' Borgatti, S. P., & Everett, M. G. (1992). Notions of position in social
#' network analysis. *Sociological Methodology*, 22, 1-35.
#' \doi{10.2307/270991}
#'
#' McKay, B. D., & Piperno, A. (2014). Practical graph isomorphism, II.
#' *Journal of Symbolic Computation*, 60, 94-112.
#' \doi{10.1016/j.jsc.2013.09.003}
#'
#' @return An integer vector of class numbers (starting at 1) aligned with
#'   the rows of `x$nodes`.
#' @examples
#' g <- netframe(data.frame(from = c(1, 1, 2, 2, 3, 3), to = c(2, 3, 4, 5, 6, 7)))
#' automorphic_equivalence(g)  # {1}, {2, 3}, {4, 5, 6, 7}
#' @export
automorphic_equivalence <- function(x, max_steps = 1e6) {
  check_netframe(x)
  if (!is.numeric(max_steps) || length(max_steps) != 1L || is.na(max_steps) ||
      max_steps < 1 || max_steps != floor(max_steps) || max_steps > 2^53) {
    stop("`max_steps` must be a positive whole number.", call. = FALSE)
  }
  ptr <- engine_graph(x)
  engine_call("automorphic_equivalence",
              .Call(C_automorphic_equivalence, ptr, as.double(max_steps))) + 1L
}

#' Blockmodel of a partition
#'
#' Reduces a network to its positions: given a partition of the nodes into
#' blocks, `blockmodel()` computes the density of ties within and between
#' every pair of blocks, the binary image matrix that codes each block as
#' present or absent at a cutoff, and how well the image fits the network
#' (White, Boorman, and Breiger 1976; Wasserman and Faust 1994, ch. 10).
#' `block_summary()` lists the same quantities as one row per pair of
#' blocks.
#'
#' @inheritParams betweenness
#' @param membership The block of each node: a vector aligned with the rows
#'   of `x$nodes`, or the name of a column of `x$nodes`, for example the
#'   result of [concor()] or [regular_equivalence()]. Only equality of
#'   values matters; blocks are listed in order of first appearance.
#' @param weighted Whether to use total edge weights as tie values, so that
#'   block densities are mean tie values. `NULL` uses weights when the
#'   network has a weight column; they are used as given, whatever the
#'   network's `weight_type`.
#' @param loops Whether self-loops and the diagonal of the adjacency matrix
#'   count: with `FALSE` (the default, and the usual practice) the diagonal
#'   blocks have \eqn{n_r (n_r - 1)} pairs, with `TRUE` \eqn{n_r^2}.
#' @param cutoff The density at or above which a block is a 1-block in the
#'   image. `NULL` uses the density of the whole network, the alpha-density
#'   criterion of Arabie, Boorman, and Levitt (1978) that UCINET's
#'   Blockmodel routine and Pajek apply by default.
#' @param bm A blockmodel from `blockmodel()`.
#'
#' @details
#' The density of block (r, s) is the sum of the ties from nodes in r to
#' nodes in s divided by the number of such pairs. Without weights a tie
#' is 1 when there is an edge from i to j, so parallel edges count once;
#' with weights it is the total weight of those edges, as in sna's
#' `blockmodel(block.content = "density")` on a valued matrix. An
#' undirected edge is a tie in both directions, so the density matrix of
#' an undirected network is symmetric. A one-node block without loops has
#' no pairs on the diagonal: its density and image are `NaN`.
#'
#' The fit is the squared Pearson correlation between the ties and the
#' density of the block each pair of nodes falls in, over the pairs that
#' count (Wasserman and Faust 1994, ch. 10): the share of the variance in
#' the ties that the block densities account for, 1 for a perfect
#' blockmodel and 0 when every block has the same density. It is `NaN`
#' when the ties or the densities have no variance. UCINET reports the
#' same correlation for its density blockmodels.
#'
#' @references
#' White, H. C., Boorman, S. A., & Breiger, R. L. (1976). Social structure
#' from multiple networks. I. Blockmodels of roles and positions.
#' *American Journal of Sociology*, 81(4), 730-780. \doi{10.1086/226141}
#'
#' Arabie, P., Boorman, S. A., & Levitt, P. R. (1978). Constructing
#' blockmodels: How and why. *Journal of Mathematical Psychology*, 17(1),
#' 21-63. \doi{10.1016/0022-2496(78)90034-2}
#'
#' Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and
#' Applications*. Cambridge University Press.
#' \doi{10.1017/CBO9780511815478}
#'
#' @return `blockmodel()` returns a list of class `"blockmodel"` with
#'   `membership` (the block labels, aligned with `x$nodes`), `blocks` (the
#'   distinct labels, in order), `sizes`, the block-by-block matrices
#'   `density`, `ties` (sums of tie values), `pairs`, and `image` (1, 0, or
#'   `NaN`), the `cutoff` used, `overall_density`, and `r_squared`.
#'   `block_summary()` returns a data frame with one row per ordered pair
#'   of blocks: `from`, `to`, `size_from`, `size_to`, `ties`, `pairs`,
#'   `density`, and `image`.
#' @seealso [concor()], [regular_equivalence()], [structural_equivalence()].
#' @examples
#' g <- example_network("florentine")
#' bm <- blockmodel(g, concor(g, depth = 2))
#' bm$density
#' bm$image
#' bm$r_squared
#' block_summary(bm)
#' @export
blockmodel <- function(x, membership, weighted = NULL, loops = FALSE, cutoff = NULL) {
  check_netframe(x)
  weighted <- positions_weighted(x, weighted)
  loops <- check_flag(loops, "loops")
  if (!is.null(cutoff) && (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff))) {
    stop("`cutoff` must be a single number, or NULL for the overall density.", call. = FALSE)
  }
  grp <- node_groups(x, membership, "membership")
  k <- length(grp$labels)
  ptr <- engine_graph(x)
  r <- engine_call("blockmodel", .Call(C_block_model, ptr, grp$codes, k, weighted, loops))
  names_ <- as.character(grp$labels)
  density <- r[[1L]]
  ties <- r[[2L]]
  pairs <- r[[3L]]
  dimnames(density) <- dimnames(ties) <- dimnames(pairs) <- list(names_, names_)
  overall <- r[[4L]][[1L]]
  if (is.null(cutoff)) cutoff <- overall
  image <- ifelse(is.nan(density), NaN, as.numeric(density >= cutoff))
  dim(image) <- dim(density)
  dimnames(image) <- dimnames(density)
  structure(
    list(
      membership = grp$values,
      blocks = grp$labels,
      sizes = tabulate(grp$codes + 1L, nbins = k),
      density = density,
      ties = ties,
      pairs = pairs,
      image = image,
      cutoff = cutoff,
      overall_density = overall,
      r_squared = r[[4L]][[2L]]
    ),
    class = "blockmodel"
  )
}

#' @rdname blockmodel
#' @export
block_summary <- function(bm) {
  if (!inherits(bm, "blockmodel")) {
    stop("`bm` must be the result of blockmodel().", call. = FALSE)
  }
  k <- length(bm$blocks)
  r <- rep(seq_len(k), times = k)
  s <- rep(seq_len(k), each = k)
  data.frame(
    from = bm$blocks[r],
    to = bm$blocks[s],
    size_from = bm$sizes[r],
    size_to = bm$sizes[s],
    ties = as.vector(bm$ties),
    pairs = as.vector(bm$pairs),
    density = as.vector(bm$density),
    image = as.vector(bm$image),
    stringsAsFactors = FALSE
  )
}

#' @export
print.blockmodel <- function(x, ...) {
  k <- length(x$blocks)
  cat(sprintf("<blockmodel> %d blocks of %d nodes; cutoff %s, R-squared %s\n", k,
              length(x$membership), format(x$cutoff, digits = 4),
              format(x$r_squared, digits = 4)))
  cat("Block densities:\n")
  print(round(x$density, 4), ...)
  invisible(x)
}

# A grouping of the nodes given as a vector aligned with x$nodes or the
# name of a column of x$nodes: the values, their distinct labels in order
# of first appearance, and 0-based codes into those labels.
node_groups <- function(x, groups, name) {
  n <- nrow(x$nodes)
  if (is.character(groups) && length(groups) == 1L && groups %in% names(x$nodes)) {
    groups <- x$nodes[[groups]]
  }
  if (is.factor(groups)) groups <- as.character(groups)
  if (!is.atomic(groups)) {
    stop(sprintf("`%s` must be an atomic vector (integer, character, factor, ...).", name),
         call. = FALSE)
  }
  if (length(groups) != n) {
    stop(sprintf("`%s` must have one value per node (%d), or name a column of `x$nodes`.",
                 name, n), call. = FALSE)
  }
  if (anyNA(groups)) {
    stop(sprintf("`%s` must not contain missing values.", name), call. = FALSE)
  }
  labels <- unique(groups)
  list(values = groups, labels = labels, codes = match(groups, labels) - 1L)
}
