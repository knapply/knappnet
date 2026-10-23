# PageRank (man/pagerank.Rd).
pagerank <- function(x, damping = 0.85, personalized = NULL, directed = TRUE, weighted = NULL,
                     dangling = c("personalized", "uniform"),
                     convention = c("igraph", "networkx"), tol = 1e-12, max_iter = 10000L,
                     threads = 0L) {
  check_netframe(x)
  dangling <- match.arg(dangling)
  convention <- match.arg(convention)
  if (!is.numeric(damping) || length(damping) != 1L || is.na(damping) ||
      damping < 0 || damping > 1) {
    stop("`damping` must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.null(personalized)) {
    personalized <- node_values(x, personalized, "personalized")
  }
  directed <- check_flag(directed, "directed")
  weighted <- positions_weighted(x, weighted)
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
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha < 0) {
    stop("`alpha` must be a single nonnegative number.", call. = FALSE)
  }
  if (is.character(beta) || length(beta) != 1L) beta <- node_values(x, beta, "beta")
  if (!is.numeric(beta) || anyNA(beta)) {
    stop("`beta` must be numeric without missing values.", call. = FALSE)
  }
  check_flag_or_null(normalized, "normalized")
  if (is.null(normalized)) normalized <- convention == "networkx"
  weighted <- positions_weighted(x, weighted)
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
  if (!is.numeric(exponent) || length(exponent) != 1L || !is.finite(exponent)) {
    stop("`exponent` must be a single finite number.", call. = FALSE)
  }
  loops <- check_flag(loops, "loops")
  rescale <- check_flag(rescale, "rescale")
  if (convention == "sna" && loops) {
    stop("The sna convention reads the network without self-loops; use `loops = FALSE`.",
         call. = FALSE)
  }
  weighted <- positions_weighted(x, weighted)
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
  weighted <- positions_weighted(x, weighted)
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

check_iteration <- function(tol, max_iter) {
  if (!is.numeric(tol) || length(tol) != 1L || !isTRUE(tol > 0)) {
    stop("`tol` must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !isTRUE(max_iter >= 1) ||
      max_iter != floor(max_iter)) {
    stop("`max_iter` must be a positive whole number.", call. = FALSE)
  }
}

# Numeric node values: a vector aligned with x$nodes, or the name of a
# column of x$nodes.
node_values <- function(x, values, what) {
  if (is.character(values) && length(values) == 1L) {
    if (!values %in% names(x$nodes)) {
      stop(sprintf("`%s` names no column of `x$nodes`.", what), call. = FALSE)
    }
    values <- x$nodes[[values]]
  }
  if (!is.numeric(values) || length(values) != nrow(x$nodes) || anyNA(values)) {
    stop(sprintf(paste0("`%s` must be a numeric vector with one value per node (%d) and ",
                        "no missing values, or the name of such a column of `x$nodes`."),
                 what, nrow(x$nodes)), call. = FALSE)
  }
  as.double(values)
}
