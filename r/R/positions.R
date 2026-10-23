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
  weighted <- positions_weighted(x, weighted)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call(fn, .Call(C_structural_holes, ptr, measure, weighted, threads))
}

# Gould-Fernandez brokerage roles (man/brokerage.Rd).
brokerage <- function(x, groups, threads = 0L) {
  check_netframe(x)
  threads <- resolve_threads(threads)
  n <- nrow(x$nodes)
  if (is.character(groups) && length(groups) == 1L && groups %in% names(x$nodes)) {
    groups <- x$nodes[[groups]]
  }
  if (length(groups) != n) {
    stop(sprintf("`groups` must have one value per node (%d), or name a column of `x$nodes`.",
                 n), call. = FALSE)
  }
  if (anyNA(groups)) {
    stop("`groups` must not contain missing values.", call. = FALSE)
  }
  if (!is.atomic(groups)) {
    stop("`groups` must be an atomic vector (integer, character, factor, ...).",
         call. = FALSE)
  }
  codes <- match(groups, unique(groups)) - 1L
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
  weighted <- positions_weighted(x, weighted)
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
  if (!is.numeric(depth) || length(depth) != 1L || is.na(depth) || depth < 1 ||
      depth != floor(depth) || depth > .Machine$integer.max) {
    stop("`depth` must be a whole number of at least 1.", call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L || !isTRUE(tolerance > 0) ||
      !isTRUE(tolerance < 1)) {
    stop("`tolerance` must be a number between 0 and 1.", call. = FALSE)
  }
  max_iter <- check_count(max_iter, "max_iter")
  weighted <- positions_weighted(x, weighted)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  r <- engine_call("concor", .Call(C_concor, ptr, as.integer(depth), weighted,
                                   as.double(tolerance), max_iter, threads))
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

# Whether a measure that reads weights as given should use them.
positions_weighted <- function(x, weighted) {
  check_flag_or_null(weighted, "weighted")
  weighted <- if (is.null(weighted)) !is.null(x$roles$weight) else isTRUE(weighted)
  if (weighted) require_weights(x)
  weighted
}
