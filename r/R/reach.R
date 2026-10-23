# Reach centrality (man/reach.Rd).
reach <- function(x, k = 2, mode = c("out", "in", "all"), weighted = NULL, weights_as = NULL,
                  normalized = FALSE, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  k <- check_number(k, "k", lower = 0, finite = FALSE)
  normalized <- check_flag(normalized, "normalized")
  threads <- resolve_threads(threads)
  w <- resolve_weighting(x, weighted, weights_as)
  ptr <- engine_graph(x)
  engine_call("reach", .Call(C_reach, ptr, as.double(k), mode_code(mode), w$weighted, w$code,
                             normalized, threads))
}

# See man/reach.Rd.
distance_weighted_reach <- function(x, mode = c("out", "in", "all"), weighted = NULL,
                                    weights_as = NULL, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  threads <- resolve_threads(threads)
  w <- resolve_weighting(x, weighted, weights_as)
  ptr <- engine_graph(x)
  engine_call("distance_weighted_reach",
              .Call(C_closeness, ptr, TRUE, mode_code(mode), w$weighted, w$code, TRUE, 0L,
                    threads))
}

# Proximity prestige (man/proximity_prestige.Rd).
proximity_prestige <- function(x, weighted = NULL, weights_as = NULL, threads = 0L) {
  check_netframe(x)
  threads <- resolve_threads(threads)
  w <- resolve_weighting(x, weighted, weights_as)
  ptr <- engine_graph(x)
  engine_call("proximity_prestige",
              .Call(C_closeness, ptr, FALSE, mode_code("in"), w$weighted, w$code, TRUE, 1L,
                    threads))
}

# Edge betweenness (man/edge_betweenness.Rd).
edge_betweenness <- function(x, weighted = NULL, normalized = FALSE, threads = 0L,
                             weights_as = NULL) {
  check_netframe(x)
  normalized <- check_flag(normalized, "normalized")
  threads <- resolve_threads(threads)
  w <- resolve_weighting(x, weighted, weights_as)
  ptr <- engine_graph(x)
  engine_call("edge_betweenness",
              .Call(C_edge_betweenness, ptr, w$weighted, w$code, normalized, threads))
}

# Information centrality (man/information_centrality.Rd).
information_centrality <- function(x, weighted = NULL, convention = c("sna", "networkx")) {
  check_netframe(x)
  convention <- match.arg(convention)
  weighted <- resolve_additive_weights(x, weighted)
  ptr <- engine_graph(x)
  out <- engine_call("information_centrality", .Call(C_information_centrality, ptr, weighted))
  if (convention == "networkx") {
    n <- nrow(x$nodes)
    if (n > 1L && (components(x, "weak")$no > 1L || anyNA(out))) {
      stop("The network is not connected; NetworkX does not compute information ",
           "centrality for such networks. Use `convention = \"sna\"`, which gives 0 to ",
           "isolates, or analyse components separately.", call. = FALSE)
    }
    out <- out / n
  }
  out
}

# Bridging centrality (man/bridging_centrality.Rd).
bridging_centrality <- function(x, threads = 0L) {
  check_netframe(x)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  b <- engine_call("bridging_centrality",
                   .Call(C_betweenness, ptr, FALSE, 0L, FALSE, threads))
  b * engine_call("bridging_centrality", .Call(C_bridging_coefficient, ptr))
}

# See man/bridging_centrality.Rd.
bridging_coefficient <- function(x) {
  check_netframe(x)
  ptr <- engine_graph(x)
  engine_call("bridging_coefficient", .Call(C_bridging_coefficient, ptr))
}

# Centralization of any scores (man/centralize.Rd).
centralize <- function(scores, theoretical_max = NULL) {
  if (!is.numeric(scores)) {
    stop("`scores` must be a numeric vector.", call. = FALSE)
  }
  if (!is.null(theoretical_max) &&
      (!is.numeric(theoretical_max) || length(theoretical_max) != 1L ||
       is.na(theoretical_max) || theoretical_max <= 0)) {
    stop("`theoretical_max` must be a single positive number, or NULL.", call. = FALSE)
  }
  if (!length(scores) || anyNA(scores)) return(NaN)
  raw <- sum(max(scores) - scores)
  if (is.null(theoretical_max)) raw else raw / theoretical_max
}
