#' Key players
#'
#' Finds the set of `k` nodes that, removed, most fragment the network
#' (KPP-Neg) or that reach most of it (KPP-Pos), following Borgatti (2006).
#' `key_player_score()` evaluates the same criteria for a given set. See
#' the help page for definitions and conventions.
#'
#' @param x A [netframe()].
#' @param k Number of key players.
#' @param set Node identifiers of the set to score.
#' @param measure `"fragmentation"`, `"distance_fragmentation"`, `"m_reach"`,
#'   or `"distance_reach"`.
#' @param m Distance limit for m-reach.
#' @param mode Direction of paths in directed networks.
#' @param weighted,weights_as Whether and how edge weights are used as path
#'   lengths.
#' @param method `"greedy"` or `"exact"`.
#' @param restarts Number of greedy restarts.
#' @param seed Seed of the random starting sets.
#' @param threads Number of threads.
#'
#' @return `key_players()`: a list with `set`, `index`, `score`, `measure`,
#'   and `evaluations`. `key_player_score()`: a number.
#' @export
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
  if (!is.numeric(restarts) || length(restarts) != 1L || is.na(restarts) || restarts < 1 ||
      restarts != floor(restarts)) {
    stop("`restarts` must be a positive whole number.", call. = FALSE)
  }
  m <- check_reach_limit(m)
  w <- resolve_weighting(x, weighted, weights_as)
  seed <- qap_seed(seed)
  threads <- resolve_threads(threads)
  res <- engine_call("key_players",
    .Call(C_key_players, engine_graph(x), as.double(k), kpp_code(measure), m, mode_code(mode),
          w$weighted, w$code, if (method == "exact") 1L else 0L, as.double(restarts), seed,
          threads))
  index <- res[[1]] + 1L
  list(set = x$nodes[[x$roles$id]][index], index = index, score = res[[2]],
       measure = measure, evaluations = res[[3]])
}

#' @rdname key_players
#' @export
key_player_score <- function(x, set, measure = c("fragmentation", "distance_fragmentation",
                                                 "m_reach", "distance_reach"),
                             m = 1, mode = c("out", "in", "all"), weighted = NULL,
                             weights_as = NULL) {
  check_netframe(x)
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  m <- check_reach_limit(m)
  index <- set_positions(x, set, "set", allow_empty = TRUE)
  w <- resolve_weighting(x, weighted, weights_as)
  engine_call("key_player_score",
    .Call(C_kpp_score, engine_graph(x), as.integer(index - 1L), kpp_code(measure), m,
          mode_code(mode), w$weighted, w$code))
}

#' Group centrality
#'
#' Degree, closeness, and betweenness of a group of nodes taken together
#' (Everett and Borgatti 1999). See the help page for definitions and
#' conventions.
#'
#' @inheritParams key_players
#' @param group Node identifiers of the group's members.
#' @param measure `"degree"`, `"closeness"`, or `"betweenness"`.
#' @param mode For directed networks: which ties (degree) or the direction
#'   of paths (closeness). Ignored by betweenness.
#' @param normalized Whether to normalize.
#'
#' @return A number.
#' @export
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
  index <- set_positions(x, group, "group", allow_empty = FALSE)
  w <- resolve_weighting(x, weighted, weights_as)
  threads <- resolve_threads(threads)
  code <- switch(measure, degree = 0L, closeness = 1L, betweenness = 2L)
  engine_call("group_centrality",
    .Call(C_group_centrality, engine_graph(x), as.integer(index - 1L), code, mode_code(mode),
          w$weighted, w$code, normalized, threads))
}

#' Flow betweenness
#'
#' The part of the maximum flow between every pair of other nodes that must
#' pass through each node (Freeman, Borgatti, and White 1991). See the help
#' page for definitions and conventions.
#'
#' @inheritParams key_players
#' @param scale `"raw"`, `"normalized"`, or `"fractional"`.
#' @param capacity Whether edge weights are capacities. `NULL` uses them
#'   when the network's weights are strengths.
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
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

check_reach_limit <- function(m) {
  if (!is.numeric(m) || length(m) != 1L || is.na(m) || m < 0) {
    stop("`m` must be a single non-negative number.", call. = FALSE)
  }
  as.double(m)
}

# Positions (1-based) of a set of distinct node identifiers in x$nodes.
set_positions <- function(x, nodes, what, allow_empty) {
  if (!is.atomic(nodes) || (!allow_empty && !length(nodes))) {
    stop(sprintf("`%s` must be a vector of node identifiers.", what), call. = FALSE)
  }
  ids <- x$nodes[[x$roles$id]]
  if (is.factor(ids)) ids <- as.character(ids)
  if (is.factor(nodes)) nodes <- as.character(nodes)
  if (anyNA(nodes)) {
    stop(sprintf("`%s` must not contain missing values (NA).", what), call. = FALSE)
  }
  pos <- match(nodes, ids)
  if (anyNA(pos)) {
    shown <- utils::head(unique(nodes[is.na(pos)]), 5L)
    stop(sprintf("%d node(s) in `%s` are not in the network, for example: %s.",
                 length(unique(nodes[is.na(pos)])), what,
                 paste(format(shown), collapse = ", ")), call. = FALSE)
  }
  if (anyDuplicated(pos)) {
    stop(sprintf("`%s` lists a node more than once.", what), call. = FALSE)
  }
  pos
}
