# Node-level centrality: degree and strength, betweenness, closeness and
# harmonic centrality, eigenvector centrality, reach, proximity prestige,
# edge betweenness, and information and bridging centrality. Spectral
# measures (PageRank, Katz, HITS, Bonacich power) live in spectral.R, group
# centrality and key players in resilience.R, and two-mode centrality in
# twomode.R. Each function follows the same shape: check_netframe(), its own
# argument checks, resolve_threads(), engine_graph() or
# engine_graph_simple() for the convention, then one engine_call().

# Degree and strength (man/degree.Rd).
degree <- function(x, mode = c("all", "out", "in"), loops = TRUE, weighted = FALSE,
                   convention = c("igraph", "networkx", "sna", "ucinet")) {
  check_netframe(x)
  mode <- match.arg(mode)
  convention <- match.arg(convention)
  loops <- check_flag(loops, "loops")
  weighted <- check_flag(weighted, "weighted")
  if (weighted) require_weights(x)
  code <- switch(mode, out = 1L, `in` = 2L, all = 3L)
  simple <- convention %in% c("sna", "ucinet")
  ptr <- if (simple) engine_graph_simple(x) else engine_graph(x)
  if (simple) loops <- FALSE
  engine_call(if (weighted) "strength" else "degree",
              .Call(C_degree, ptr, code, loops, weighted), edge_rows = !simple)
}

# See man/degree.Rd.
strength <- function(x, mode = c("all", "out", "in"), loops = TRUE) {
  degree(x, mode = match.arg(mode), loops = loops, weighted = TRUE)
}

# Betweenness centrality (man/betweenness.Rd).
betweenness <- function(x, weighted = NULL, normalized = FALSE, threads = 0L,
                        weights_as = NULL,
                        convention = c("igraph", "networkx", "sna", "ucinet")) {
  check_netframe(x)
  convention <- match.arg(convention)
  check_flag_or_null(weighted, "weighted")
  normalized <- check_flag(normalized, "normalized")
  threads <- resolve_threads(threads)
  binary <- convention %in% c("sna", "ucinet")
  if (binary && isTRUE(weighted)) {
    stop(sprintf("The %s convention uses binary ties; use `weighted = FALSE`.", convention),
         call. = FALSE)
  }
  w <- resolve_weighting(x, if (binary) FALSE else weighted, weights_as)
  ptr <- switch(convention,
    igraph = engine_graph(x),
    # NetworkX keeps the shortest of several parallel edges: the smallest
    # distance, or the largest strength.
    networkx = engine_graph_simple(x, keep = if (w$code == 1L) "max" else "min"),
    engine_graph_simple(x)
  )
  engine_call("betweenness",
              .Call(C_betweenness, ptr, w$weighted, w$code, normalized, threads),
              edge_rows = convention == "igraph")
}

# Closeness and harmonic centrality (man/closeness.Rd).
closeness <- function(x, mode = NULL, weighted = NULL, normalized = NULL, threads = 0L,
                      weights_as = NULL,
                      convention = c("igraph", "networkx", "sna", "ucinet")) {
  closeness_common(x, FALSE, mode, weighted, normalized, threads, weights_as,
                   match.arg(convention), "closeness")
}

# See man/closeness.Rd.
harmonic_centrality <- function(x, mode = NULL, weighted = NULL, normalized = NULL,
                                threads = 0L, weights_as = NULL,
                                convention = c("igraph", "networkx", "sna", "ucinet")) {
  closeness_common(x, TRUE, mode, weighted, normalized, threads, weights_as,
                   match.arg(convention), "harmonic_centrality")
}

closeness_common <- function(x, harmonic, mode, weighted, normalized, threads, weights_as,
                             convention, fn) {
  check_netframe(x)
  if (is.null(mode)) mode <- if (convention == "networkx") "in" else "out"
  mode <- check_choice(mode, c("out", "in", "all"), "mode")
  check_flag_or_null(weighted, "weighted")
  check_flag_or_null(normalized, "normalized")
  threads <- resolve_threads(threads)
  if (is.null(normalized)) {
    normalized <- !(convention == "igraph" || (harmonic && convention == "networkx"))
  }
  binary <- convention %in% c("sna", "ucinet")
  if (binary && isTRUE(weighted)) {
    stop(sprintf("The %s convention computes closeness on binary ties; use `weighted = FALSE`.",
                 convention), call. = FALSE)
  }
  w <- resolve_weighting(x, if (binary) FALSE else weighted, weights_as)
  unreachable <- switch(convention, igraph = 0L, networkx = 1L, sna = 2L, ucinet = 3L)
  ptr <- if (binary) engine_graph_simple(x) else engine_graph(x)
  engine_call(fn, .Call(C_closeness, ptr, harmonic, mode_code(mode), w$weighted, w$code,
                        normalized, unreachable, threads), edge_rows = !binary)
}

# Eigenvector centrality (man/eigenvector_centrality.Rd).
eigenvector_centrality <- function(x, mode = NULL, weighted = NULL, scale = NULL,
                                   convention = c("igraph", "networkx", "sna", "ucinet"),
                                   tol = 1e-12, max_iter = 10000L, threads = 0L) {
  check_netframe(x)
  convention <- match.arg(convention)
  if (is.null(mode)) {
    mode <- switch(convention, igraph = "in", networkx = "in", sna = "out", ucinet = "all")
  }
  mode <- check_choice(mode, c("in", "out", "all"), "mode")
  if (is.null(scale)) scale <- if (convention == "igraph") "max" else "unit"
  scale <- check_choice(scale, c("max", "unit"), "scale")
  threads <- resolve_threads(threads)
  weighted <- resolve_additive_weights(x, weighted)
  check_iteration(tol, max_iter)

  loops <- switch(convention, igraph = 0L, networkx = 1L, sna = 2L, ucinet = 2L)
  ptr <- switch(convention,
    sna = engine_graph_simple(x),
    ucinet = engine_graph_simple(x, undirected = TRUE),
    engine_graph(x)
  )
  if (convention == "ucinet") {
    csize <- components(x, "weak")$csize
    if (sum(csize > 1L) > 1L) {
      stop("More than one component has ties. UCINET then reports lambda-squared ",
           "scores computed within components, which are not implemented; use ",
           "another convention or analyse components separately.", call. = FALSE)
    }
  }
  res <- engine_call("eigenvector_centrality",
    .Call(C_eigenvector, ptr, mode_code(mode), weighted, if (scale == "max") 0L else 1L,
          loops, as.double(tol), as.double(max_iter), threads),
    edge_rows = convention %in% c("igraph", "networkx"))
  if (!res[[3]] && length(res[[1]]) > 1L) {
    if (convention == "networkx") {
      stop("The network is not ", if (x$directed && mode != "all") "strongly ",
           "connected; NetworkX does not compute eigenvector centrality for such networks.",
           call. = FALSE)
    }
    warning("The network is not ", if (x$directed && mode != "all") "strongly ",
            "connected, so eigenvector centrality is zero outside the dominant ",
            "component(s) and may not be unique.", call. = FALSE)
  }
  res[[1]]
}

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
