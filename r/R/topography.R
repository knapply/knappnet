# Edge density (man/graph_density.Rd).
graph_density <- function(x, loops = FALSE, convention = c("igraph", "networkx", "sna", "ucinet")) {
  check_netframe(x)
  convention <- match.arg(convention)
  loops <- check_flag(loops, "loops")
  if (convention %in% c("sna", "ucinet")) {
    if (loops) {
      stop("The sna and UCINET conventions ignore self-loops; use `loops = FALSE`.",
           call. = FALSE)
    }
    ptr <- engine_graph_simple(x)
    return(engine_call("graph_density", .Call(C_density, ptr, FALSE), edge_rows = FALSE))
  }
  if (convention == "networkx") {
    if (loops) {
      stop("NetworkX's density() has no option for self-loops; use `loops = FALSE`.",
           call. = FALSE)
    }
    if (nrow(x$nodes) == 1L) return(0)
  }
  ptr <- engine_graph(x)
  engine_call("graph_density", .Call(C_density, ptr, loops))
}

# Connected components (man/components.Rd).
components <- function(x, mode = c("weak", "strong")) {
  check_netframe(x)
  mode <- match.arg(mode)
  ptr <- engine_graph(x)
  m <- engine_call("components",
                   .Call(C_components, ptr, if (mode == "strong") 1L else 0L)) + 1L
  csize <- tabulate(m, nbins = if (length(m)) max(m) else 0L)
  list(membership = m, csize = csize, no = length(csize))
}

# See man/components.Rd.
fragmentation <- function(x, mode = c("weak", "strong")) {
  check_netframe(x)
  mode <- match.arg(mode)
  ptr <- engine_graph(x)
  engine_call("fragmentation", .Call(C_fragmentation, ptr, if (mode == "strong") 1L else 0L))
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

# Freeman centralization (man/centralization.Rd).
centralization <- function(x, measure = c("degree", "betweenness", "closeness", "eigenvector"),
                           mode = NULL,
                           loops = TRUE, normalized = TRUE,
                           convention = c("igraph", "networkx", "sna", "ucinet"),
                           threads = 0L) {
  check_netframe(x)
  measure <- match.arg(measure)
  convention <- match.arg(convention)
  if (convention == "networkx") {
    stop("NetworkX has no centralization index; use the igraph, sna, or UCINET convention.",
         call. = FALSE)
  }
  if (measure == "eigenvector" && convention != "igraph") {
    stop("Eigenvector centralization is provided for the igraph convention only.",
         call. = FALSE)
  }
  if (is.null(mode)) {
    mode <- if (measure != "eigenvector" && (measure == "closeness" || convention == "ucinet")) {
      "out"
    } else {
      "all"
    }
  }
  mode <- check_choice(mode, c("all", "out", "in"), "mode")
  loops <- check_flag(loops, "loops")
  normalized <- check_flag(normalized, "normalized")
  threads <- resolve_threads(threads)
  code <- switch(measure, degree = 0L, betweenness = 1L, closeness = 2L, eigenvector = 3L)
  conv <- switch(convention, igraph = 0L, sna = 2L, ucinet = 3L)
  ptr <- engine_graph(x)
  if (measure == "eigenvector" && nrow(x$nodes) > 1L) {
    csize <- components(x, "weak")$csize
    if (sum(csize > 1L) > 1L) {
      warning("More than one component has ties, so eigenvector centrality may not be ",
              "unique and neither is its centralization.", call. = FALSE)
    }
  }
  engine_call("centralization", .Call(C_centralization, ptr, code, mode_code(mode), loops,
                                      normalized, conv, threads))
}
