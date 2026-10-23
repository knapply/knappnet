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
  labels <- node_groups(x, membership, "membership")$codes
  resolution <- check_resolution(resolution)
  weighted <- resolve_additive_weights(x, weighted)
  ptr <- engine_graph(x)
  engine_call("modularity", .Call(C_modularity, ptr, labels, resolution, weighted))
}

# Louvain community detection (man/cluster_louvain.Rd).
cluster_louvain <- function(x, resolution = 1, weighted = NULL, seed = NULL) {
  check_netframe(x)
  require_undirected(x, "Louvain community detection")
  resolution <- check_resolution(resolution)
  weighted <- resolve_additive_weights(x, weighted)
  seed <- check_seed(seed)
  ptr <- engine_graph(x)
  res <- engine_call("cluster_louvain", .Call(C_louvain, ptr, resolution, weighted, seed))
  structure(res[[1L]], modularity = res[[2L]])
}
