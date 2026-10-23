# Mixing matrices (man/mixing_matrix.Rd).
mixing_matrix <- function(x, dim1, dim2 = NULL, direction = c("all", "out", "in")) {
  check_netframe(x)
  direction <- match.arg(direction)
  d1 <- node_categories(x, dim1, "dim1", missing_ok = TRUE)
  same <- is.null(dim2)
  if (!same) {
    d2 <- node_categories(x, dim2, "dim2", missing_ok = TRUE)
    same <- identical(d1$values, d2$values)
  }
  if (same) d2 <- d1
  ptr <- engine_graph(x)
  res <- engine_call("mixing_matrix",
    .Call(C_mixing_matrix, ptr, d1$codes, length(d1$levels), d2$codes, length(d2$levels),
          same, mode_code(direction)))
  k1 <- length(d1$levels)
  k2 <- length(d2$levels)
  if (res[[2]]) {
    out <- matrix(res[[1]], nrow = k2, ncol = k1, byrow = TRUE,
                  dimnames = list(d2$levels, d1$levels))
  } else {
    out <- matrix(res[[1]], nrow = k1, ncol = k2, byrow = TRUE,
                  dimnames = list(d1$levels, d2$levels))
  }
  out
}

# Assortativity coefficients (man/assortativity.Rd).
assortativity <- function(x, values, values_in = NULL, directed = TRUE, normalized = TRUE,
                          convention = c("igraph", "networkx", "homophily")) {
  check_netframe(x)
  convention <- match.arg(convention)
  v <- node_values(x, values, "values")
  vin <- if (is.null(values_in)) NULL else node_values(x, values_in, "values_in")
  directed <- check_flag(directed, "directed")
  normalized <- check_flag(normalized, "normalized")
  if (!is.null(vin) && !x$directed) {
    warning("`values_in` is ignored for undirected networks.", call. = FALSE)
    vin <- NULL
  }
  ptr <- engine_graph(x)
  engine_call("assortativity",
    .Call(C_assortativity_values, ptr, v, vin, directed, assort_loops(convention),
          normalized))
}

# See man/assortativity.Rd.
assortativity_nominal <- function(x, types, directed = TRUE, normalized = TRUE,
                                  convention = c("igraph", "networkx", "homophily")) {
  check_netframe(x)
  convention <- match.arg(convention)
  t <- node_categories(x, types, "types", missing_ok = FALSE)
  directed <- check_flag(directed, "directed")
  normalized <- check_flag(normalized, "normalized")
  if (convention == "homophily") directed <- FALSE
  ptr <- engine_graph(x)
  engine_call("assortativity_nominal",
    .Call(C_assortativity_nominal, ptr, t$codes, length(t$levels), directed,
          assort_loops(convention), normalized))
}

# See man/assortativity.Rd.
assortativity_degree <- function(x, directed = TRUE, weighted = FALSE, normalized = TRUE,
                                 convention = c("igraph", "networkx", "homophily")) {
  check_netframe(x)
  convention <- match.arg(convention)
  directed <- check_flag(directed, "directed")
  weighted <- check_flag(weighted, "weighted")
  normalized <- check_flag(normalized, "normalized")
  if (weighted) require_weights(x)
  ptr <- engine_graph(x)
  engine_call("assortativity_degree",
    .Call(C_assortativity_degree, ptr, weighted, directed, assort_loops(convention),
          normalized))
}

assort_loops <- function(convention) if (convention == "networkx") 1L else 0L

# E-I index (man/ei_index.Rd).
ei_index <- function(x, groups, scope = c("global", "group", "node"),
                     direction = c("all", "out", "in"),
                     convention = c("krackhardt", "homophily")) {
  check_netframe(x)
  scope <- match.arg(scope)
  direction <- match.arg(direction)
  convention <- match.arg(convention)
  g <- node_categories(x, groups, "groups", missing_ok = FALSE)
  ptr <- engine_graph(x)
  warn_parallel(x)
  res <- engine_call("ei_index",
    .Call(C_ei_index, ptr, g$codes, length(g$levels),
          switch(scope, global = 0L, group = 1L, node = 2L), mode_code(direction),
          if (convention == "krackhardt") 0L else 1L))
  out <- res[[1]]
  if (scope == "group") names(out) <- g$levels
  attr(out, "external") <- res[[2]]
  attr(out, "internal") <- res[[3]]
  out
}

# Within-group tie shares (man/homophily.Rd).
homophily <- function(x, groups, direction = c("all", "out", "in")) {
  check_netframe(x)
  direction <- match.arg(direction)
  g <- node_categories(x, groups, "groups", missing_ok = FALSE)
  k <- length(g$levels)
  ptr <- engine_graph(x)
  v <- engine_call("homophily", .Call(C_group_shares, ptr, g$codes, k, mode_code(direction)))
  m <- matrix(v, nrow = k, ncol = 5L, byrow = TRUE)
  data.frame(group = g$levels, size = as.integer(m[, 1]), internal = m[, 2],
             external = m[, 3], observed = m[, 4], expected = m[, 5],
             stringsAsFactors = FALSE)
}

warn_parallel <- function(x) {
  src <- x$index$src
  tgt <- x$index$tgt
  key <- if (x$directed) paste(src, tgt) else paste(pmin(src, tgt), pmax(src, tgt))
  if (anyDuplicated(key)) {
    warning("The network has parallel edges; each counts as a separate tie.",
            call. = FALSE)
  }
}
