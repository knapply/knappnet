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
  egos <- if (is.null(nodes)) seq_len(nrow(x$nodes)) else node_positions(x, nodes)
  members <- .Call(C_ego_members, engine_graph(x), as.integer(egos - 1L), order,
                   mode_code(mode), if (isTRUE(include_ego)) 0L else 1L, resolve_threads(threads))
  out <- lapply(members, function(keep) induced_netframe(x, keep))
  names(out) <- as.character(x$nodes[[x$roles$id]][egos])
  out
}

# Egocentric measures (man/ego_summary.Rd).
ego_summary <- function(x, order = 1L, mode = c("all", "out", "in"), betweenness = TRUE,
                        threads = 0L) {
  check_netframe(x)
  res <- ego_measures(x, order, match.arg(mode), isTRUE(betweenness), threads)
  out <- data.frame(node = x$nodes[[x$roles$id]], size = res[[1]], ties = res[[2]],
                    pairs = res[[3]], density = res[[4]], components = res[[5]],
                    stringsAsFactors = FALSE)
  if (isTRUE(betweenness)) out$ego_betweenness <- res[[6]]
  out
}

# See man/ego_summary.Rd.
ego_size <- function(x, order = 1L, mode = c("all", "out", "in"), include_ego = TRUE,
                     threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  order <- check_count(order, "order")
  .Call(C_ego_sizes, engine_graph(x), NULL, order, mode_code(mode),
        if (isTRUE(include_ego)) 0L else 1L, resolve_threads(threads))
}

# See man/ego_summary.Rd.
ego_density <- function(x, order = 1L, mode = c("all", "out", "in"), threads = 0L) {
  check_netframe(x)
  ego_measures(x, order, match.arg(mode), FALSE, threads)[[4]]
}

# See man/ego_summary.Rd.
ego_betweenness <- function(x, order = 1L, mode = c("all", "out", "in"), threads = 0L) {
  check_netframe(x)
  ego_measures(x, order, match.arg(mode), TRUE, threads)[[6]]
}

ego_measures <- function(x, order, mode, betweenness, threads) {
  order <- check_count(order, "order")
  if (order < 1L) stop("`order` must be at least 1.", call. = FALSE)
  .Call(C_ego_summary, engine_graph(x), order, mode_code(mode), betweenness,
        resolve_threads(threads))
}


# The subgraph induced by the nodes at 1-based positions `keep` (in
# ascending order): those nodes and every edge between two of them, with
# all columns and roles (including a two-mode network's node type).
induced_netframe <- function(x, keep) {
  inside <- logical(nrow(x$nodes))
  inside[keep] <- TRUE
  e <- inside[x$index$src + 1L] & inside[x$index$tgt + 1L]
  nodes <- x$nodes[keep, , drop = FALSE]
  edges <- x$edges[e, , drop = FALSE]
  rownames(nodes) <- NULL
  rownames(edges) <- NULL
  weight_type <- if (is.null(x$roles$weight_type)) "distance" else x$roles$weight_type
  out <- netframe(edges, nodes, directed = x$directed, from = x$roles$from,
                  to = x$roles$to, id = x$roles$id, weight = x$roles$weight,
                  weight_type = weight_type, graph = x$graph)
  extra <- setdiff(names(x$roles), names(out$roles))
  if (length(extra)) out$roles[extra] <- x$roles[extra]
  out
}
