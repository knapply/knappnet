# Key players (man/key_players.Rd).
key_players <- function(x, k, measure = c("fragmentation", "distance_fragmentation",
                                          "m_reach", "distance_reach"),
                        m = 1, mode = c("out", "in", "all"), weighted = NULL,
                        weights_as = NULL, method = c("greedy", "exact"), restarts = 10L,
                        seed = NULL, threads = 0L) {
  check_netframe(x)
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  method <- match.arg(method)
  n <- nrow(x$nodes)
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1 || k > n || k != floor(k)) {
    stop(sprintf("`k` must be a whole number between 1 and the number of nodes (%d).", n),
         call. = FALSE)
  }
  restarts <- check_count(restarts, "restarts", min = 1L)
  m <- check_number(m, "m", lower = 0, finite = FALSE)
  w <- resolve_weighting(x, weighted, weights_as)
  seed <- check_seed(seed)
  threads <- resolve_threads(threads)
  res <- engine_call("key_players",
    .Call(C_key_players, engine_graph(x), as.double(k), kpp_code(measure), m, mode_code(mode),
          w$weighted, w$code, if (method == "exact") 1L else 0L, as.double(restarts), seed,
          threads))
  index <- res[[1]] + 1L
  list(set = x$nodes[[x$roles$id]][index], index = index, score = res[[2]],
       measure = measure, evaluations = res[[3]])
}

# See man/key_players.Rd.
key_player_score <- function(x, set, measure = c("fragmentation", "distance_fragmentation",
                                                 "m_reach", "distance_reach"),
                             m = 1, mode = c("out", "in", "all"), weighted = NULL,
                             weights_as = NULL) {
  check_netframe(x)
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  m <- check_number(m, "m", lower = 0, finite = FALSE)
  index <- node_positions(x, set, "set", distinct = TRUE)
  w <- resolve_weighting(x, weighted, weights_as)
  engine_call("key_player_score",
    .Call(C_kpp_score, engine_graph(x), as.integer(index - 1L), kpp_code(measure), m,
          mode_code(mode), w$weighted, w$code))
}

# Group centrality (man/group_centrality.Rd).
group_centrality <- function(x, group, measure = c("degree", "closeness", "betweenness"),
                             mode = c("out", "in", "all"), normalized = TRUE,
                             weighted = NULL, weights_as = NULL, threads = 0L) {
  check_netframe(x)
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  normalized <- check_flag(normalized, "normalized")
  check_flag_or_null(weighted, "weighted")
  if (measure == "degree") {
    if (isTRUE(weighted)) {
      stop("Group degree is defined on binary ties; use `weighted = FALSE`.", call. = FALSE)
    }
    weighted <- FALSE
  }
  index <- node_positions(x, group, "group", empty_ok = FALSE, distinct = TRUE)
  w <- resolve_weighting(x, weighted, weights_as)
  threads <- resolve_threads(threads)
  code <- switch(measure, degree = 0L, closeness = 1L, betweenness = 2L)
  engine_call("group_centrality",
    .Call(C_group_centrality, engine_graph(x), as.integer(index - 1L), code, mode_code(mode),
          w$weighted, w$code, normalized, threads))
}

# Flow betweenness (man/flow_betweenness.Rd).
flow_betweenness <- function(x, scale = c("raw", "normalized", "fractional"),
                             capacity = NULL, threads = 0L) {
  check_netframe(x)
  scale <- match.arg(scale)
  check_flag_or_null(capacity, "capacity")
  if (is.null(capacity)) {
    capacity <- !is.null(x$roles$weight) && identical(x$roles$weight_type, "strength")
  }
  if (capacity) require_weights(x)
  threads <- resolve_threads(threads)
  code <- switch(scale, raw = 0L, normalized = 1L, fractional = 2L)
  engine_call("flow_betweenness",
    .Call(C_flow_betweenness, engine_graph(x), capacity, code, threads))
}

kpp_code <- function(measure) {
  switch(measure, fragmentation = 0L, distance_fragmentation = 1L, m_reach = 2L,
         distance_reach = 3L)
}
