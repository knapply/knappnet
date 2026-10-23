# Core-periphery structure (man/core_periphery.Rd).
core_periphery <- function(x, method = c("discrete", "continuous"),
                           links = c("missing", "ones", "zeros"),
                           fit = c("minres", "correlation"), weighted = NULL,
                           starts = 20L, max_stall = 200L, tenure = NULL, max_iter = 1000L,
                           tol = 1e-10, seed = NULL, threads = 0L) {
  check_netframe(x)
  method <- match.arg(method)
  links <- match.arg(links)
  fit <- match.arg(fit)
  weighted <- resolve_additive_weights(x, weighted)
  starts <- check_count(starts, "starts")
  max_stall <- check_count(max_stall, "max_stall")
  tenure <- if (is.null(tenure)) 0L else check_count(tenure, "tenure")
  max_iter <- check_count(max_iter, "max_iter")
  tol <- check_number(tol, "tol", lower = 0)
  seed <- check_seed(seed)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  if (method == "discrete") {
    code <- switch(links, missing = 0L, ones = 1L, zeros = 2L)
    res <- engine_call("core_periphery",
                       .Call(C_core_periphery, ptr, code, weighted, starts, max_stall, tenure,
                             seed, threads))
    out <- list(method = "discrete", links = links, core = res[[1L]], correlation = res[[2L]],
                network = x)
  } else {
    res <- engine_call("core_periphery",
                       .Call(C_core_periphery_continuous, ptr,
                             if (fit == "minres") 0L else 1L, weighted, max_iter,
                             tol, starts, seed, threads))
    out <- list(method = "continuous", fit = fit, coreness = res[[1L]],
                correlation = res[[2L]], residual = res[[3L]], iterations = res[[4L]],
                network = x)
  }
  structure(out, class = "core_periphery")
}

# See man/core_periphery.Rd.
print.core_periphery <- function(x, ...) {
  if (x$method == "discrete") {
    cat(sprintf("Discrete core-periphery fit (core-periphery pairs %s): %d core, %d periphery\n",
                switch(x$links, missing = "left out", ones = "expected tied",
                       zeros = "expected untied"),
                sum(x$core), sum(!x$core)))
  } else {
    cat(sprintf("Continuous core-periphery fit (%s): coreness of %d nodes (unit norm)\n",
                if (x$fit == "minres") "MINRES" else "maximum correlation",
                length(x$coreness)))
    if (x$fit == "minres") {
      cat(sprintf("Residual sum of squares: %s\n", format(x$residual, digits = 4)))
    }
  }
  cat(sprintf("Correlation with the ideal pattern: %s\n", format(x$correlation, digits = 4)))
  invisible(x)
}

# Small-world quotient (man/small_world.Rd).
small_world <- function(x, local = FALSE, reps = 100L, seed = NULL, threads = 0L) {
  check_netframe(x)
  local <- check_flag(local, "local")
  reps <- check_count(reps, "reps", min = 1L)
  seed <- check_seed(seed)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  v <- engine_call("small_world", .Call(C_small_world, ptr, local, reps, seed, threads))
  structure(list(clustering = v[[1L]], mean_distance = v[[2L]], random_clustering = v[[3L]],
                 random_mean_distance = v[[4L]], quotient = v[[5L]], reps = reps,
                 local = local, n = nrow(x$nodes)),
            class = "small_world")
}

# See man/small_world.Rd.
print.small_world <- function(x, ...) {
  cat(sprintf("Small-world quotient: %s\n", format(x$quotient, digits = 4)))
  cat(sprintf("  %s: %s observed, %s in %d random networks (ratio %s)\n",
              if (x$local) "mean local clustering" else "transitivity",
              format(x$clustering, digits = 4), format(x$random_clustering, digits = 4),
              x$reps, format(x$clustering / x$random_clustering, digits = 4)))
  cat(sprintf("  mean distance: %s observed, %s in random networks (ratio %s)\n",
              format(x$mean_distance, digits = 4), format(x$random_mean_distance, digits = 4),
              format(x$mean_distance / x$random_mean_distance, digits = 4)))
  invisible(x)
}

# Compactness and distance-weighted fragmentation (man/compactness.Rd).
compactness <- function(x, mode = c("out", "in", "all"), weighted = NULL, weights_as = NULL,
                        threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  w <- resolve_weighting(x, weighted, weights_as)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("compactness", .Call(C_compactness, ptr, mode_code(mode), w$weighted, w$code,
                                   threads))
}

# See man/compactness.Rd.
distance_weighted_fragmentation <- function(x, mode = c("out", "in", "all"), weighted = NULL,
                                            weights_as = NULL, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  w <- resolve_weighting(x, weighted, weights_as)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  1 - engine_call("distance_weighted_fragmentation",
                  .Call(C_compactness, ptr, mode_code(mode), w$weighted, w$code, threads))
}

# Average degree and the degree distribution (man/average_degree.Rd).
average_degree <- function(x, mode = c("all", "out", "in"), loops = TRUE) {
  check_netframe(x)
  mode <- match.arg(mode)
  loops <- check_flag(loops, "loops")
  d <- degree(x, mode = mode, loops = loops)
  if (!length(d)) return(NaN)
  mean(d)
}

# See man/average_degree.Rd.
degree_distribution <- function(x, mode = c("all", "out", "in"), loops = TRUE,
                                cumulative = FALSE) {
  check_netframe(x)
  mode <- match.arg(mode)
  loops <- check_flag(loops, "loops")
  cumulative <- check_flag(cumulative, "cumulative")
  d <- degree(x, mode = mode, loops = loops)
  values <- sort(unique(d))
  count <- vapply(values, function(v) sum(d == v), integer(1))
  out <- data.frame(degree = values, count = count, proportion = count / length(d))
  if (cumulative) out$at_least <- rev(cumsum(rev(out$proportion)))
  rownames(out) <- NULL
  out
}
