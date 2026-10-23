# Classical multidimensional scaling (man/classical_mds.Rd).
classical_mds <- function(d, k = 2) {
  d <- as_dissimilarity(d)
  n <- nrow(d)
  k <- check_count(k, "k", min = 1L)
  if (n < 2L) {
    stop("Classical scaling needs at least two points.", call. = FALSE)
  }
  if (k > n - 1L) {
    stop(sprintf("`k` must be at most %d, one less than the number of points.", n - 1L),
         call. = FALSE)
  }
  r <- engine_call("classical_mds", .Call(C_classical_mds, d, n, k), edge_rows = FALSE)
  dims <- r[[2L]]
  points <- r[[1L]][, seq_len(dims), drop = FALSE]
  dimnames(points) <- list(rownames(d), if (dims) paste0("Dim", seq_len(dims)))
  stress <- if (dims) {
    engine_call("classical_mds", .Call(C_mds_stress, d, n, points, dims), edge_rows = FALSE)
  } else {
    NaN
  }
  structure(list(points = points, eigenvalues = r[[3L]], gof = r[[4L]], stress = stress),
            class = "classical_mds")
}

# See man/classical_mds.Rd.
print.classical_mds <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sprintf("<classical_mds> %d points in %d dimension%s\n", nrow(x$points), ncol(x$points),
              if (ncol(x$points) == 1L) "" else "s"))
  cat(sprintf("Stress: %s; goodness of fit: %s, %s\n", format(x$stress, digits = digits),
              format(x$gof[[1L]], digits = digits), format(x$gof[[2L]], digits = digits)))
  cat("Leading eigenvalues:\n")
  print(utils::head(x$eigenvalues, 5L), digits = digits, ...)
  invisible(x)
}

# See man/classical_mds.Rd.
mds_stress <- function(d, points) {
  d <- as_dissimilarity(d)
  n <- nrow(d)
  if (!is.numeric(points) || !is.matrix(points)) {
    stop("`points` must be a numeric matrix with one row per point.", call. = FALSE)
  }
  if (nrow(points) != n) {
    stop(sprintf("`points` has %d rows but `d` has %d points.", nrow(points), n),
         call. = FALSE)
  }
  if (anyNA(points) || any(!is.finite(points))) {
    stop("`points` must be finite.", call. = FALSE)
  }
  storage.mode(points) <- "double"
  engine_call("mds_stress", .Call(C_mds_stress, d, n, points, ncol(points)), edge_rows = FALSE)
}

# A dissimilarity matrix as a square double matrix. Structural checks
# (symmetry, zero diagonal, finiteness) are made by the engine.
as_dissimilarity <- function(d) {
  if (inherits(d, "dist")) d <- as.matrix(d)
  if (!is.numeric(d) || !is.matrix(d) || nrow(d) != ncol(d)) {
    stop("`d` must be a square numeric matrix of dissimilarities, or a \"dist\" object.",
         call. = FALSE)
  }
  if (anyNA(d)) {
    stop("`d` must not contain missing values.", call. = FALSE)
  }
  storage.mode(d) <- "double"
  d
}

# Distances between nodes for scaling (man/mds_distances.Rd).
mds_distances <- function(x, method = c("geodesic", "euclidean", "hamming", "correlation"),
                          unreachable = "max_plus_one", weighted = NULL, weights_as = NULL,
                          threads = 0L) {
  check_netframe(x)
  method <- match.arg(method)
  if (method == "geodesic") {
    d <- distances(x, mode = "all", weighted = weighted, weights_as = weights_as,
                   threads = threads)
    if (any(is.infinite(d))) {
      finite <- d[is.finite(d)]
      longest <- if (length(finite)) max(finite) else 0
      fill <- if (is.numeric(unreachable)) {
        if (length(unreachable) != 1L || is.na(unreachable) || unreachable < 0) {
          stop("`unreachable` must be \"max_plus_one\", \"diameter\", or a non-negative number.",
               call. = FALSE)
        }
        unreachable
      } else {
        switch(check_choice(unreachable, c("max_plus_one", "diameter"), "unreachable"),
               max_plus_one = longest + 1, diameter = longest)
      }
      d[is.infinite(d)] <- fill
    }
    return(d)
  }
  m <- structural_equivalence(x, method = method, weighted = weighted, threads = threads)
  if (method == "correlation") {
    m <- 1 - m
    diag(m) <- 0
  }
  m
}
