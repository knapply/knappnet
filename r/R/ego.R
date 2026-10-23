# Ego networks (man/ego_network.Rd).
ego_network <- function(x, node, order = 1L, mode = c("all", "out", "in"),
                        include_ego = TRUE) {
  check_netframe(x)
  if (length(node) != 1L) {
    stop("`node` must be a single node identifier; use ego_networks() for several.",
         call. = FALSE)
  }
  ego_networks(x, node, order = order, mode = match.arg(mode),
               include_ego = include_ego, threads = 1L)[[1L]]
}

# See man/ego_network.Rd.
ego_networks <- function(x, nodes = NULL, order = 1L, mode = c("all", "out", "in"),
                         include_ego = TRUE, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  order <- check_count(order, "order")
  include_ego <- check_flag(include_ego, "include_ego")
  egos <- if (is.null(nodes)) seq_len(nrow(x$nodes)) else node_positions(x, nodes, "nodes")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  members <- engine_call("ego_networks",
                         .Call(C_ego_members, ptr, as.integer(egos - 1L), order,
                               mode_code(mode), if (include_ego) 0L else 1L, threads))
  out <- lapply(members, function(keep) induced_netframe(x, keep))
  names(out) <- node_ids(x)[egos]
  out
}

# Egocentric measures (man/ego_summary.Rd).
ego_summary <- function(x, order = 1L, mode = c("all", "out", "in"), betweenness = TRUE,
                        threads = 0L) {
  check_netframe(x)
  betweenness <- check_flag(betweenness, "betweenness")
  res <- ego_measures(x, order, match.arg(mode), betweenness, threads, "ego_summary")
  out <- data.frame(node = x$nodes[[x$roles$id]], size = res[[1]], ties = res[[2]],
                    pairs = res[[3]], density = res[[4]], components = res[[5]],
                    stringsAsFactors = FALSE)
  if (betweenness) out$ego_betweenness <- res[[6]]
  out
}

# See man/ego_summary.Rd.
ego_size <- function(x, order = 1L, mode = c("all", "out", "in"), include_ego = TRUE,
                     threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  order <- check_count(order, "order")
  include_ego <- check_flag(include_ego, "include_ego")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("ego_size", .Call(C_ego_sizes, ptr, NULL, order, mode_code(mode),
                                if (include_ego) 0L else 1L, threads))
}

# See man/ego_summary.Rd.
ego_density <- function(x, order = 1L, mode = c("all", "out", "in"), threads = 0L) {
  check_netframe(x)
  ego_measures(x, order, match.arg(mode), FALSE, threads, "ego_density")[[4]]
}

# See man/ego_summary.Rd.
ego_betweenness <- function(x, order = 1L, mode = c("all", "out", "in"), threads = 0L) {
  check_netframe(x)
  ego_measures(x, order, match.arg(mode), TRUE, threads, "ego_betweenness")[[6]]
}

ego_measures <- function(x, order, mode, betweenness, threads, fn) {
  order <- check_count(order, "order", min = 1L)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call(fn, .Call(C_ego_summary, ptr, order, mode_code(mode), betweenness, threads))
}
