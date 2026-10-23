# Positions and roles: structural holes, brokerage roles, structural,
# regular, and automorphic equivalence, CONCOR, and blockmodels with their
# print and plot methods.

# Structural holes: constraint, effective size, efficiency, and hierarchy (man/net_constraint.Rd).
net_constraint <- function(x, weighted = NULL, threads = 0L) {
  structural_holes(x, 0L, weighted, threads, "net_constraint")
}

# See man/net_constraint.Rd.
net_effective_size <- function(x, weighted = NULL, threads = 0L) {
  structural_holes(x, 1L, weighted, threads, "net_effective_size")
}

# See man/net_constraint.Rd.
net_efficiency <- function(x, weighted = NULL, threads = 0L) {
  structural_holes(x, 2L, weighted, threads, "net_efficiency")
}

# See man/net_constraint.Rd.
net_hierarchy <- function(x, weighted = NULL, threads = 0L) {
  structural_holes(x, 3L, weighted, threads, "net_hierarchy")
}

structural_holes <- function(x, measure, weighted, threads, fn) {
  check_net_frame(x)
  weighted <- resolve_additive_weights(x, weighted)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call(fn, .Call(C_structural_holes, ptr, measure, weighted, threads))
}

# Gould-Fernandez brokerage roles (man/net_brokerage.Rd).
net_brokerage <- function(x, groups, threads = 0L) {
  check_net_frame(x)
  threads <- resolve_threads(threads)
  codes <- node_groups(x, groups, "groups")$codes
  ptr <- engine_graph(x)
  m <- engine_call("net_brokerage", .Call(C_brokerage, ptr, codes, threads))
  out <- data.frame(node = x$nodes[[x$roles$id]], stringsAsFactors = FALSE)
  roles <- c("coordinator", "itinerant", "representative", "gatekeeper",
             "liaison", "total")
  for (r in seq_along(roles)) out[[roles[[r]]]] <- m[, r]
  out
}

# Structural equivalence (man/net_struct_equiv.Rd).
net_struct_equiv <- function(x, method = c("euclidean", "correlation", "hamming"),
                                   weighted = NULL, threads = 0L) {
  check_net_frame(x)
  method <- match.arg(method)
  code <- switch(method, euclidean = 0L, correlation = 1L, hamming = 2L)
  weighted <- resolve_additive_weights(x, weighted)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  m <- engine_call("net_struct_equiv",
                   .Call(C_structural_equivalence, ptr, code, weighted, threads))
  ids <- as.character(x$nodes[[x$roles$id]])
  dimnames(m) <- list(ids, ids)
  m
}

# CONCOR blockmodel (man/net_concor.Rd).
net_concor <- function(x, depth = 1L, weighted = NULL, tolerance = 1e-7, max_iter = 100L,
                   threads = 0L) {
  check_net_frame(x)
  depth <- check_count(depth, "depth", min = 1L)
  tolerance <- check_number(tolerance, "tolerance", lower = 0, upper = 1, open = TRUE)
  max_iter <- check_count(max_iter, "max_iter")
  weighted <- resolve_additive_weights(x, weighted)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  r <- engine_call("net_concor", .Call(C_concor, ptr, depth, weighted, tolerance,
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

# Regular equivalence (CATREGE) (man/net_regular_equiv.Rd).
net_regular_equiv <- function(x, seed = NULL, weighted = NULL) {
  check_net_frame(x)
  weighted <- resolve_additive_weights(x, weighted)
  codes <- if (is.null(seed)) NULL else node_groups(x, seed, "seed")$codes
  ptr <- engine_graph(x)
  r <- engine_call("net_regular_equiv", .Call(C_regular_equivalence, ptr, codes, weighted))
  structure(r[[1L]] + 1L, rounds = r[[2L]])
}

# Automorphic equivalence (man/net_auto_equiv.Rd).
net_auto_equiv <- function(x, max_steps = 1e6) {
  check_net_frame(x)
  if (!is.numeric(max_steps) || length(max_steps) != 1L || is.na(max_steps) ||
      max_steps < 1 || max_steps != floor(max_steps) || max_steps > 2^53) {
    stop("`max_steps` must be a positive whole number.", call. = FALSE)
  }
  ptr <- engine_graph(x)
  engine_call("net_auto_equiv",
              .Call(C_automorphic_equivalence, ptr, as.double(max_steps))) + 1L
}

# Blockmodel of a partition (man/net_blockmodel.Rd).
net_blockmodel <- function(x, membership, weighted = NULL, loops = FALSE, cutoff = NULL) {
  check_net_frame(x)
  weighted <- resolve_additive_weights(x, weighted)
  loops <- check_flag(loops, "loops")
  if (!is.null(cutoff) && (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff))) {
    stop("`cutoff` must be a single number, or NULL for the overall density.", call. = FALSE)
  }
  grp <- node_groups(x, membership, "membership")
  k <- length(grp$labels)
  ptr <- engine_graph(x)
  r <- engine_call("net_blockmodel", .Call(C_block_model, ptr, grp$codes, k, weighted, loops))
  names_ <- as.character(grp$labels)
  density <- r[[1L]]
  ties <- r[[2L]]
  pairs <- r[[3L]]
  dimnames(density) <- dimnames(ties) <- dimnames(pairs) <- list(names_, names_)
  overall <- r[[4L]][[1L]]
  if (is.null(cutoff)) cutoff <- overall
  image <- ifelse(is.nan(density), NaN, as.numeric(density >= cutoff))
  dim(image) <- dim(density)
  dimnames(image) <- dimnames(density)
  structure(
    list(
      membership = grp$values,
      blocks = grp$labels,
      sizes = tabulate(grp$codes + 1L, nbins = k),
      density = density,
      ties = ties,
      pairs = pairs,
      image = image,
      cutoff = cutoff,
      overall_density = overall,
      r_squared = r[[4L]][[2L]]
    ),
    class = "blockmodel"
  )
}

# See man/net_blockmodel.Rd.
net_block_summary <- function(bm) {
  if (!inherits(bm, "blockmodel")) {
    stop("`bm` must be the result of net_blockmodel().", call. = FALSE)
  }
  k <- length(bm$blocks)
  r <- rep(seq_len(k), times = k)
  s <- rep(seq_len(k), each = k)
  data.frame(
    from = bm$blocks[r],
    to = bm$blocks[s],
    size_from = bm$sizes[r],
    size_to = bm$sizes[s],
    ties = as.vector(bm$ties),
    pairs = as.vector(bm$pairs),
    density = as.vector(bm$density),
    image = as.vector(bm$image),
    stringsAsFactors = FALSE
  )
}

# See man/net_blockmodel.Rd.
print.blockmodel <- function(x, ...) {
  k <- length(x$blocks)
  cat(sprintf("<blockmodel> %d blocks of %d nodes; cutoff %s, R-squared %s\n", k,
              length(x$membership), format(x$cutoff, digits = 4),
              format(x$r_squared, digits = 4)))
  cat("Block densities:\n")
  print(round(x$density, 4), ...)
  invisible(x)
}

# Plot a blockmodel (man/plot.blockmodel.Rd).
plot.blockmodel <- function(x, network = NULL, which = c("both", "density", "matrix"),
                            values = TRUE, ...) {
  which <- match.arg(which)
  values <- check_flag(values, "values")
  if (!is.null(network)) {
    check_net_frame(network)
    if (nrow(network$nodes) != length(x$membership)) {
      stop("`network` must have one node per element of `x$membership`.", call. = FALSE)
    }
  } else if (which != "density") {
    if (which == "matrix") {
      stop("The permuted adjacency matrix needs `network`.", call. = FALSE)
    }
    which <- "density"
  }
  if (which == "both") {
    old <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(old), add = TRUE)
  }
  dots <- list(...)
  if (which != "density") {
    args <- c(list(x = network, order = x$membership, key = FALSE), dots)
    if (is.null(args$main)) args$main <- "Permuted adjacency matrix"
    do.call(net_plot_matrix, args)
  }
  if (which != "matrix") {
    d <- x$density
    finite <- d[is.finite(d)]
    top <- if (length(finite)) max(1, finite) else 1
    args <- c(list(m = d, values = values, range = c(0, top), digits = 2L), dots)
    if (is.null(args$main)) args$main <- "Block densities"
    do.call(draw_heatmap, args)
    k <- nrow(d)
    ones <- which(is.finite(x$image) & x$image == 1, arr.ind = TRUE)
    if (nrow(ones)) {
      graphics::rect(ones[, 2] - 0.5, k - ones[, 1] + 0.5, ones[, 2] + 0.5,
                     k - ones[, 1] + 1.5, border = plot_accent, lwd = 2.5)
    }
  }
  invisible(x)
}
