# Centrality from matrix iterations: PageRank, Katz, Bonacich power, and
# hub and authority scores. Eigenvector centrality itself is in
# centrality.R.

# PageRank (man/pagerank.Rd).
pagerank <- function(x, damping = 0.85, personalized = NULL, directed = TRUE, weighted = NULL,
                     dangling = c("personalized", "uniform"),
                     convention = c("igraph", "networkx"), tol = 1e-12, max_iter = 10000L,
                     threads = 0L) {
  check_netframe(x)
  dangling <- match.arg(dangling)
  convention <- match.arg(convention)
  damping <- check_number(damping, "damping", lower = 0, upper = 1)
  if (!is.null(personalized)) {
    personalized <- node_values(x, personalized, "personalized")
  }
  directed <- check_flag(directed, "directed")
  weighted <- resolve_additive_weights(x, weighted)
  check_iteration(tol, max_iter)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("pagerank",
    .Call(C_pagerank, ptr, as.double(damping), directed, weighted, personalized,
          if (dangling == "personalized") 0L else 1L,
          if (convention == "igraph") 0L else 1L, as.double(tol), as.double(max_iter),
          threads))
}

# Katz centrality (man/katz_centrality.Rd).
katz_centrality <- function(x, alpha = 0.1, beta = 1, mode = c("in", "out", "all"),
                            weighted = NULL, normalized = NULL,
                            convention = c("igraph", "networkx"), tol = 1e-12,
                            max_iter = 10000L, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  convention <- match.arg(convention)
  alpha <- check_number(alpha, "alpha", lower = 0)
  if (is.character(beta) || length(beta) != 1L) beta <- node_values(x, beta, "beta")
  if (!is.numeric(beta) || anyNA(beta)) {
    stop("`beta` must be numeric without missing values.", call. = FALSE)
  }
  check_flag_or_null(normalized, "normalized")
  if (is.null(normalized)) normalized <- convention == "networkx"
  weighted <- resolve_additive_weights(x, weighted)
  check_iteration(tol, max_iter)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("katz_centrality",
    .Call(C_katz, ptr, as.double(alpha), as.double(beta), mode_code(mode), weighted,
          if (convention == "igraph") 2L else 1L, normalized, as.double(tol),
          as.double(max_iter), threads))
}

# Bonacich power centrality (man/power_centrality.Rd).
power_centrality <- function(x, exponent = 1, mode = c("out", "in", "all"), weighted = NULL,
                             loops = FALSE, rescale = FALSE,
                             convention = c("igraph", "sna"), tol = 1e-12,
                             max_iter = 10000L, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  convention <- match.arg(convention)
  exponent <- check_number(exponent, "exponent")
  loops <- check_flag(loops, "loops")
  rescale <- check_flag(rescale, "rescale")
  if (convention == "sna" && loops) {
    stop("The sna convention reads the network without self-loops; use `loops = FALSE`.",
         call. = FALSE)
  }
  weighted <- resolve_additive_weights(x, weighted)
  check_iteration(tol, max_iter)
  threads <- resolve_threads(threads)
  ptr <- if (convention == "sna") engine_graph_simple(x) else engine_graph(x)
  engine_call("power_centrality",
    .Call(C_power_centrality, ptr, as.double(exponent), mode_code(mode), weighted,
          if (loops) 1L else 2L, rescale, as.double(tol), as.double(max_iter), threads),
    edge_rows = convention == "igraph")
}

# Hub and authority scores (man/hub_score.Rd).
hub_score <- function(x, weighted = NULL, scale = NULL, convention = c("igraph", "networkx"),
                      tol = 1e-12, max_iter = 10000L, threads = 0L) {
  hits_common(x, "hub_score", weighted, scale, match.arg(convention), tol, max_iter,
              threads)
}

# See man/hub_score.Rd.
authority_score <- function(x, weighted = NULL, scale = NULL,
                            convention = c("igraph", "networkx"), tol = 1e-12,
                            max_iter = 10000L, threads = 0L) {
  hits_common(x, "authority_score", weighted, scale, match.arg(convention), tol, max_iter,
              threads)
}

hits_common <- function(x, fn, weighted, scale, convention, tol, max_iter, threads) {
  check_netframe(x)
  if (is.null(scale)) scale <- if (convention == "igraph") "max" else "sum"
  scale <- check_choice(scale, c("max", "sum", "unit"), "scale")
  weighted <- resolve_additive_weights(x, weighted)
  check_iteration(tol, max_iter)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  res <- engine_call(fn,
    .Call(C_hits, ptr, weighted, if (convention == "igraph") 0L else 1L,
          switch(scale, max = 0L, sum = 1L, unit = 2L), as.double(tol),
          as.double(max_iter), threads))
  if (!res[[4]] && nrow(x$nodes) > 1L && !all(is.nan(res[[1]]))) {
    warning("The largest singular value of the adjacency matrix is not simple, so hub ",
            "and authority scores are not unique; the result is the limit of ",
            "Kleinberg's iteration from equal scores.", call. = FALSE)
  }
  if (fn == "hub_score") res[[1]] else res[[2]]
}
