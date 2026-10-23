# Regular equivalence (CATREGE) (man/regular_equivalence.Rd).
regular_equivalence <- function(x, seed = NULL, weighted = NULL) {
  check_netframe(x)
  weighted <- resolve_additive_weights(x, weighted)
  codes <- if (is.null(seed)) NULL else node_groups(x, seed, "seed")$codes
  ptr <- engine_graph(x)
  r <- engine_call("regular_equivalence", .Call(C_regular_equivalence, ptr, codes, weighted))
  structure(r[[1L]] + 1L, rounds = r[[2L]])
}

# Automorphic equivalence (man/automorphic_equivalence.Rd).
automorphic_equivalence <- function(x, max_steps = 1e6) {
  check_netframe(x)
  if (!is.numeric(max_steps) || length(max_steps) != 1L || is.na(max_steps) ||
      max_steps < 1 || max_steps != floor(max_steps) || max_steps > 2^53) {
    stop("`max_steps` must be a positive whole number.", call. = FALSE)
  }
  ptr <- engine_graph(x)
  engine_call("automorphic_equivalence",
              .Call(C_automorphic_equivalence, ptr, as.double(max_steps))) + 1L
}

# Blockmodel of a partition (man/blockmodel.Rd).
blockmodel <- function(x, membership, weighted = NULL, loops = FALSE, cutoff = NULL) {
  check_netframe(x)
  weighted <- resolve_additive_weights(x, weighted)
  loops <- check_flag(loops, "loops")
  if (!is.null(cutoff) && (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff))) {
    stop("`cutoff` must be a single number, or NULL for the overall density.", call. = FALSE)
  }
  grp <- node_groups(x, membership, "membership")
  k <- length(grp$labels)
  ptr <- engine_graph(x)
  r <- engine_call("blockmodel", .Call(C_block_model, ptr, grp$codes, k, weighted, loops))
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

# See man/blockmodel.Rd.
block_summary <- function(bm) {
  if (!inherits(bm, "blockmodel")) {
    stop("`bm` must be the result of blockmodel().", call. = FALSE)
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

# See man/blockmodel.Rd.
print.blockmodel <- function(x, ...) {
  k <- length(x$blocks)
  cat(sprintf("<blockmodel> %d blocks of %d nodes; cutoff %s, R-squared %s\n", k,
              length(x$membership), format(x$cutoff, digits = 4),
              format(x$r_squared, digits = 4)))
  cat("Block densities:\n")
  print(round(x$density, 4), ...)
  invisible(x)
}
