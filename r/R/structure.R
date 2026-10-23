# Measures of the network as a whole: density, average degree and the
# degree distribution, components and fragmentation, the triad census,
# centralization, compactness, the small-world quotient, and core-periphery
# structure, with the print and plot methods of the result classes.

# Edge density (man/net_density.Rd).
net_density <- function(x, loops = FALSE, convention = c("igraph", "networkx", "sna", "ucinet"),
                          relation = NULL) {
  check_net_frame(x)
  convention <- match.arg(convention)
  loops <- check_flag(loops, "loops")
  x <- relation_subset(x, relation, "net_density")
  if (convention %in% c("sna", "ucinet")) {
    if (loops) {
      stop("The sna and UCINET conventions ignore self-loops; use `loops = FALSE`.",
           call. = FALSE)
    }
    ptr <- engine_graph_simple(x)
    return(engine_call("net_density", .Call(C_density, ptr, FALSE), edge_rows = FALSE))
  }
  if (convention == "networkx") {
    if (loops) {
      stop("NetworkX's density() has no option for self-loops; use `loops = FALSE`.",
           call. = FALSE)
    }
    if (nrow(x$nodes) == 1L) return(0)
  }
  ptr <- engine_graph(x)
  engine_call("net_density", .Call(C_density, ptr, loops))
}

# Average degree and the degree distribution (man/net_avg_degree.Rd).
net_avg_degree <- function(x, mode = c("all", "out", "in"), loops = TRUE) {
  check_net_frame(x)
  mode <- match.arg(mode)
  loops <- check_flag(loops, "loops")
  d <- net_degree(x, mode = mode, loops = loops)
  if (!length(d)) return(NaN)
  mean(d)
}

# See man/net_avg_degree.Rd.
net_degree_dist <- function(x, mode = c("all", "out", "in"), loops = TRUE,
                                cumulative = FALSE) {
  check_net_frame(x)
  mode <- match.arg(mode)
  loops <- check_flag(loops, "loops")
  cumulative <- check_flag(cumulative, "cumulative")
  d <- net_degree(x, mode = mode, loops = loops)
  values <- sort(unique(d))
  count <- vapply(values, function(v) sum(d == v), integer(1))
  out <- data.frame(degree = values, count = count, proportion = count / length(d))
  if (cumulative) out$at_least <- rev(cumsum(rev(out$proportion)))
  rownames(out) <- NULL
  out
}

# Connected components (man/net_components.Rd).
net_components <- function(x, mode = c("weak", "strong")) {
  check_net_frame(x)
  mode <- match.arg(mode)
  ptr <- engine_graph(x)
  m <- engine_call("net_components",
                   .Call(C_components, ptr, if (mode == "strong") 1L else 0L)) + 1L
  csize <- tabulate(m, nbins = if (length(m)) max(m) else 0L)
  list(membership = m, csize = csize, no = length(csize))
}

# See man/net_components.Rd.
net_fragmentation <- function(x, mode = c("weak", "strong")) {
  check_net_frame(x)
  mode <- match.arg(mode)
  ptr <- engine_graph(x)
  engine_call("net_fragmentation", .Call(C_fragmentation, ptr, if (mode == "strong") 1L else 0L))
}

# Triad census (man/net_triad_census.Rd).
net_triad_census <- function(x, threads = 0L) {
  check_net_frame(x)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  out <- engine_call("net_triad_census", .Call(C_triad_census, ptr, threads))
  names(out) <- c("003", "012", "102", "021D", "021U", "021C", "111D", "111U",
                  "030T", "030C", "201", "120D", "120U", "120C", "210", "300")
  out
}

# Freeman centralization (man/net_centralization.Rd).
net_centralization <- function(x, measure = c("degree", "betweenness", "closeness", "eigenvector"),
                           mode = NULL,
                           loops = TRUE, normalized = TRUE,
                           convention = c("igraph", "networkx", "sna", "ucinet"),
                           threads = 0L) {
  check_net_frame(x)
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
    csize <- net_components(x, "weak")$csize
    if (sum(csize > 1L) > 1L) {
      warning("More than one component has ties, so eigenvector centrality may not be ",
              "unique and neither is its centralization.", call. = FALSE)
    }
  }
  engine_call("net_centralization", .Call(C_centralization, ptr, code, mode_code(mode), loops,
                                      normalized, conv, threads))
}

# Centralization of any scores (man/net_centralize.Rd).
net_centralize <- function(scores, theoretical_max = NULL) {
  if (!is.numeric(scores)) {
    stop("`scores` must be a numeric vector.", call. = FALSE)
  }
  if (!is.null(theoretical_max) &&
      (!is.numeric(theoretical_max) || length(theoretical_max) != 1L ||
       is.na(theoretical_max) || theoretical_max <= 0)) {
    stop("`theoretical_max` must be a single positive number, or NULL.", call. = FALSE)
  }
  if (!length(scores) || anyNA(scores)) return(NaN)
  raw <- sum(max(scores) - scores)
  if (is.null(theoretical_max)) raw else raw / theoretical_max
}

# Compactness and distance-weighted fragmentation (man/net_compactness.Rd).
net_compactness <- function(x, mode = c("out", "in", "all"), weighted = NULL, weights_as = NULL,
                        threads = 0L) {
  check_net_frame(x)
  mode <- match.arg(mode)
  w <- resolve_weighting(x, weighted, weights_as)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("net_compactness", .Call(C_compactness, ptr, mode_code(mode), w$weighted, w$code,
                                   threads))
}

# See man/net_compactness.Rd.
net_dw_fragmentation <- function(x, mode = c("out", "in", "all"), weighted = NULL,
                                            weights_as = NULL, threads = 0L) {
  check_net_frame(x)
  mode <- match.arg(mode)
  w <- resolve_weighting(x, weighted, weights_as)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  1 - engine_call("net_dw_fragmentation",
                  .Call(C_compactness, ptr, mode_code(mode), w$weighted, w$code, threads))
}

# Small-world quotient (man/net_small_world.Rd).
net_small_world <- function(x, local = FALSE, reps = 100L, seed = NULL, threads = 0L) {
  check_net_frame(x)
  local <- check_flag(local, "local")
  reps <- check_count(reps, "reps", min = 1L)
  seed <- check_seed(seed)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  v <- engine_call("net_small_world", .Call(C_small_world, ptr, local, reps, seed, threads))
  structure(list(clustering = v[[1L]], mean_distance = v[[2L]], random_clustering = v[[3L]],
                 random_mean_distance = v[[4L]], quotient = v[[5L]], reps = reps,
                 local = local, n = nrow(x$nodes)),
            class = "small_world")
}

# See man/net_small_world.Rd.
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

# See man/net_small_world.Rd.
plot.small_world <- function(x, ...) {
  pal <- plot_qualitative(8L)
  dots <- list(...)
  old <- graphics::par(mfrow = c(1L, 2L), mar = c(3, 4, 4, 1), oma = c(0, 0, 2, 0))
  on.exit(graphics::par(old), add = TRUE)
  panel <- function(obs, rnd, ylab) {
    v <- c(observed = obs, random = rnd)
    graphics::barplot(v, col = pal[c(2L, 1L)], border = NA, ylab = ylab,
                      ylim = c(0, max(v[is.finite(v)], 0) * 1.15))
    graphics::text(c(0.7, 1.9), v, format(v, digits = 3), pos = 3, cex = 0.8, xpd = NA)
  }
  panel(x$clustering, x$random_clustering,
        if (x$local) "mean local clustering" else "transitivity")
  panel(x$mean_distance, x$random_mean_distance, "mean distance")
  main <- if (is.null(dots$main)) {
    sprintf("Small-world quotient %s (%d random networks)", format(x$quotient, digits = 3),
            x$reps)
  } else {
    dots$main
  }
  graphics::mtext(main, outer = TRUE, cex = 1.1, font = 2)
  invisible(x)
}

# Core-periphery structure (man/net_core_periphery.Rd).
net_core_periphery <- function(x, method = c("discrete", "continuous"),
                           links = c("missing", "ones", "zeros"),
                           fit = c("minres", "correlation"), weighted = NULL,
                           starts = 20L, max_stall = 200L, tenure = NULL, max_iter = 1000L,
                           tol = 1e-10, seed = NULL, threads = 0L) {
  check_net_frame(x)
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
    res <- engine_call("net_core_periphery",
                       .Call(C_core_periphery, ptr, code, weighted, starts, max_stall, tenure,
                             seed, threads))
    out <- list(method = "discrete", links = links, core = res[[1L]], correlation = res[[2L]],
                network = x)
  } else {
    res <- engine_call("net_core_periphery",
                       .Call(C_core_periphery_continuous, ptr,
                             if (fit == "minres") 0L else 1L, weighted, max_iter,
                             tol, starts, seed, threads))
    out <- list(method = "continuous", fit = fit, coreness = res[[1L]],
                correlation = res[[2L]], residual = res[[3L]], iterations = res[[4L]],
                network = x)
  }
  structure(out, class = "core_periphery")
}

# See man/net_core_periphery.Rd.
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

# See man/net_core_periphery.Rd.
plot.core_periphery <- function(x, what = NULL, layout = NULL, ..., seed = NULL) {
  if (is.null(what)) what <- if (x$method == "discrete") "matrix" else "network"
  what <- check_choice(what, c("matrix", "network"), "what")
  g <- x$network
  pal <- plot_qualitative(8L)
  if (what == "network") {
    if (x$method == "discrete") {
      plot(g, layout = layout, node_color = I(ifelse(x$core, pal[[5L]], pal[[1L]])), ...,
           legend = FALSE, seed = seed)
      graphics::legend("topleft", legend = c("core", "periphery"), pch = 21,
                       pt.bg = pal[c(5L, 1L)], col = "white", pt.cex = 1.6, bty = "n",
                       cex = 0.8)
    } else {
      plot(g, layout = layout, node_color = x$coreness, node_size = x$coreness, ..., seed = seed)
    }
    return(invisible(x))
  }
  n <- nrow(g$nodes)
  if (n == 0L) stop("The network has no nodes; nothing to plot.", call. = FALSE)
  order <- if (x$method == "discrete") {
    order(!x$core, seq_len(n))
  } else {
    order(-x$coreness, seq_len(n))
  }
  weighted <- !is.null(g$roles$weight)
  a <- as.matrix(g, attr = g$roles$weight)[order, order, drop = FALSE]
  diag(a) <- 0
  if (!weighted) a <- (a > 0) * 1
  dots <- list(...)
  if (is.null(dots$main)) {
    dots$main <- sprintf("%s core-periphery fit (correlation %s)",
                         if (x$method == "discrete") "Discrete" else "Continuous",
                         format(x$correlation, digits = 3))
  }
  labelled <- n <= 40
  old <- graphics::par(mar = c(0.5, if (labelled) 5 else 1, if (labelled) 6.5 else 3, 1),
                       xpd = NA)
  on.exit(graphics::par(old), add = TRUE)
  # Row i at the top, as a matrix is printed.
  colors <- if (weighted) {
    plot_matrix_palette(64L)
  } else {
    c("white", pal[[2L]])
  }
  graphics::image(seq_len(n), seq_len(n), t(a)[, n:1, drop = FALSE], col = colors,
                  axes = FALSE, xlab = "", ylab = "", useRaster = FALSE)
  graphics::box(col = "gray50")
  if (labelled) {
    labels <- as.character(g$nodes[[g$roles$id]])[order]
    graphics::axis(2, at = n:1, labels = labels, las = 2, tick = FALSE, cex.axis = 0.7,
                   line = -0.5)
    graphics::axis(3, at = seq_len(n), labels = labels, las = 2, tick = FALSE,
                   cex.axis = 0.7, line = -0.5)
  }
  boundary <- if (x$method == "discrete") sum(x$core) else sum(x$coreness > mean(x$coreness))
  if (boundary > 0 && boundary < n) {
    graphics::abline(v = boundary + 0.5, h = n - boundary + 0.5, col = pal[[5L]], lwd = 2)
  }
  do.call(graphics::title, c(dots, list(line = if (labelled) 5 else 1.5)))
  invisible(x)
}
