# Structural holes: constraint, effective size, efficiency, and hierarchy (man/constraint.Rd).
constraint <- function(x, weighted = NULL, threads = 0L) {
  structural_holes(x, 0L, weighted, threads, "constraint")
}

# See man/constraint.Rd.
effective_size <- function(x, weighted = NULL, threads = 0L) {
  structural_holes(x, 1L, weighted, threads, "effective_size")
}

# See man/constraint.Rd.
efficiency <- function(x, weighted = NULL, threads = 0L) {
  structural_holes(x, 2L, weighted, threads, "efficiency")
}

# See man/constraint.Rd.
hierarchy <- function(x, weighted = NULL, threads = 0L) {
  structural_holes(x, 3L, weighted, threads, "hierarchy")
}

structural_holes <- function(x, measure, weighted, threads, fn) {
  check_netframe(x)
  weighted <- resolve_additive_weights(x, weighted)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call(fn, .Call(C_structural_holes, ptr, measure, weighted, threads))
}

# Gould-Fernandez brokerage roles (man/brokerage.Rd).
brokerage <- function(x, groups, threads = 0L) {
  check_netframe(x)
  threads <- resolve_threads(threads)
  codes <- node_groups(x, groups, "groups")$codes
  ptr <- engine_graph(x)
  m <- engine_call("brokerage", .Call(C_brokerage, ptr, codes, threads))
  out <- data.frame(node = x$nodes[[x$roles$id]], stringsAsFactors = FALSE)
  roles <- c("coordinator", "itinerant", "representative", "gatekeeper",
             "liaison", "total")
  for (r in seq_along(roles)) out[[roles[[r]]]] <- m[, r]
  out
}

# Structural equivalence (man/structural_equivalence.Rd).
structural_equivalence <- function(x, method = c("euclidean", "correlation", "hamming"),
                                   weighted = NULL, threads = 0L) {
  check_netframe(x)
  method <- match.arg(method)
  code <- switch(method, euclidean = 0L, correlation = 1L, hamming = 2L)
  weighted <- resolve_additive_weights(x, weighted)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  m <- engine_call("structural_equivalence",
                   .Call(C_structural_equivalence, ptr, code, weighted, threads))
  ids <- as.character(x$nodes[[x$roles$id]])
  dimnames(m) <- list(ids, ids)
  m
}

# CONCOR blockmodel (man/concor.Rd).
concor <- function(x, depth = 1L, weighted = NULL, tolerance = 1e-7, max_iter = 100L,
                   threads = 0L) {
  check_netframe(x)
  depth <- check_count(depth, "depth", min = 1L)
  tolerance <- check_number(tolerance, "tolerance", lower = 0, upper = 1, open = TRUE)
  max_iter <- check_count(max_iter, "max_iter")
  weighted <- resolve_additive_weights(x, weighted)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  r <- engine_call("concor", .Call(C_concor, ptr, depth, weighted, tolerance,
                                   max_iter, threads))
  if (r[[2]][[1]] > 0L) {
    warning(sprintf("%d CONCOR split(s) did not converge to +1/-1 within `max_iter` iterations.",
                    r[[2]][[1]]), call. = FALSE)
  }
  if (r[[2]][[2]] > 0L) {
    warning(sprintf(paste0("%d CONCOR split(s) had a correlation within 1e-6 of zero; ",
                           "rounding error may have decided some nodes' blocks."),
                    r[[2]][[2]]), call. = FALSE)
  }
  r[[1]]
}
