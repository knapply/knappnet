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

# Triad census (man/triad_census.Rd).
triad_census <- function(x, threads = 0L) {
  check_netframe(x)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  out <- engine_call("triad_census", .Call(C_triad_census, ptr, threads))
  names(out) <- c("003", "012", "102", "021D", "021U", "021C", "111D", "111U",
                  "030T", "030C", "201", "120D", "120U", "120C", "210", "300")
  out
}
