#' Coreness (k-cores)
#'
#' The coreness of a node is the largest k such that the node belongs to
#' the k-core, the largest subgraph in which every node has degree at least
#' k (Seidman 1983). Computed with the O(m) algorithm of Batagelj and
#' Zaversnik (2003).
#'
#' Parallel edges each count. A self-loop adds 2 to its node's degree in an
#' undirected network, and 1 to out-degree and 1 to in-degree in a directed
#' one. Weights are ignored. These match igraph's `coreness()`.
#'
#' @inheritParams degree
#' @param mode For directed networks, `"all"` (in-degree plus out-degree),
#'   `"out"`, or `"in"`. Ignored for undirected networks.
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
coreness <- function(x, mode = c("all", "out", "in")) {
  check_netframe(x)
  mode <- match.arg(mode)
  code <- switch(mode, out = 1L, `in` = 2L, all = 3L)
  ptr <- engine_graph(x)
  engine_call("coreness", .Call(C_coreness, ptr, code))
}

#' Maximal cliques and the clique number
#'
#' Lists the maximal cliques of the undirected simple graph underlying the
#' network, using the Bron-Kerbosch algorithm with pivoting over a
#' degeneracy ordering.
#'
#' @inheritParams betweenness
#' @param min_size Smallest clique size to return.
#'
#' @return `max_cliques()`: a list of vectors of node identifiers (values
#'   of `x$nodes[[id]]`). `clique_number()`: a single integer.
#' @export
max_cliques <- function(x, min_size = 1L, threads = 0L) {
  check_netframe(x)
  min_size <- check_count(min_size, "min_size")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  idx <- engine_call("max_cliques", .Call(C_max_cliques, ptr, min_size, threads))
  ids <- x$nodes[[x$roles$id]]
  lapply(idx, function(i) ids[i])
}

#' @rdname max_cliques
#' @export
clique_number <- function(x, threads = 0L) {
  check_netframe(x)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("clique_number", .Call(C_clique_number, ptr, threads))
}

#' Modularity of a partition
#'
#' Newman-Girvan modularity, with weights, a resolution parameter, and the
#' directed form of Leicht and Newman (2008).
#'
#' @inheritParams betweenness
#' @param membership Group of each node, aligned with the rows of
#'   `x$nodes`. Any atomic vector (integer, character, factor); only
#'   equality of values matters. Missing values are not allowed.
#' @param resolution Resolution parameter (non-negative).
#'
#' @return A single number; `NaN` for a network without edges.
#' @export
modularity <- function(x, membership, resolution = 1, weighted = NULL) {
  check_netframe(x)
  if (length(membership) != nrow(x$nodes)) {
    stop(sprintf("`membership` has length %d, but the network has %d nodes.",
                 length(membership), nrow(x$nodes)), call. = FALSE)
  }
  if (anyNA(membership)) {
    stop("`membership` must not contain missing values.", call. = FALSE)
  }
  labels <- match(membership, unique(membership)) - 1L
  resolution <- check_resolution(resolution)
  weighted <- resolve_additive_weights(x, weighted)
  ptr <- engine_graph(x)
  engine_call("modularity", .Call(C_modularity, ptr, labels, resolution, weighted))
}

#' Louvain community detection
#'
#' Finds groups by greedy modularity optimization with the method of
#' Blondel et al. (2008).
#'
#' @inheritParams modularity
#' @param seed Seed for the random node order: a whole number between 0
#'   and 2^53. `NULL` draws one from R's random number generator, so
#'   [set.seed()] makes results reproducible.
#'
#' @return An integer vector of group labels aligned with the rows of
#'   `x$nodes`, numbered 1, 2, ... in order of first appearance, with the
#'   partition's modularity in attribute `"modularity"`.
#' @export
cluster_louvain <- function(x, resolution = 1, weighted = NULL, seed = NULL) {
  check_netframe(x)
  if (x$directed) {
    stop("Louvain community detection needs an undirected network; ",
         "rebuild it with netframe(..., directed = FALSE).", call. = FALSE)
  }
  resolution <- check_resolution(resolution)
  weighted <- resolve_additive_weights(x, weighted)
  if (is.null(seed)) seed <- floor(stats::runif(1L) * 2^31)
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
      seed < 0 || seed > 2^53 || seed != floor(seed)) {
    stop("`seed` must be a whole number between 0 and 2^53.", call. = FALSE)
  }
  ptr <- engine_graph(x)
  res <- engine_call("cluster_louvain",
                     .Call(C_louvain, ptr, resolution, weighted, as.double(seed)))
  structure(res[[1L]], modularity = res[[2L]])
}

# Weights that add up (modularity, strength) are used as given, whatever
# the network's weight_type.
resolve_additive_weights <- function(x, weighted) {
  check_flag_or_null(weighted, "weighted")
  has <- !is.null(x$roles$weight)
  weighted <- if (is.null(weighted)) has else isTRUE(weighted)
  if (weighted) require_weights(x)
  weighted
}

check_resolution <- function(resolution) {
  if (!is.numeric(resolution) || length(resolution) != 1L || is.na(resolution) ||
      !is.finite(resolution) || resolution < 0) {
    stop("`resolution` must be a single finite, non-negative number.", call. = FALSE)
  }
  as.double(resolution)
}

check_count <- function(value, what) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value < 0 ||
      value != floor(value) || value > .Machine$integer.max) {
    stop(sprintf("`%s` must be a single non-negative whole number.", what),
         call. = FALSE)
  }
  as.integer(value)
}
