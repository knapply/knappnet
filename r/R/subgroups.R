# Coreness (k-cores) (man/coreness.Rd).
coreness <- function(x, mode = c("all", "out", "in")) {
  check_netframe(x)
  mode <- match.arg(mode)
  code <- switch(mode, out = 1L, `in` = 2L, all = 3L)
  ptr <- engine_graph(x)
  engine_call("coreness", .Call(C_coreness, ptr, code))
}

# Maximal cliques and the clique number (man/max_cliques.Rd).
max_cliques <- function(x, min_size = 1L, threads = 0L) {
  check_netframe(x)
  min_size <- check_count(min_size, "min_size")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  idx <- engine_call("max_cliques", .Call(C_max_cliques, ptr, min_size, threads))
  ids <- x$nodes[[x$roles$id]]
  lapply(idx, function(i) ids[i])
}

# See man/max_cliques.Rd.
clique_number <- function(x, threads = 0L) {
  check_netframe(x)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("clique_number", .Call(C_clique_number, ptr, threads))
}

# Modularity of a partition (man/modularity.Rd).
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

# Louvain community detection (man/cluster_louvain.Rd).
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
