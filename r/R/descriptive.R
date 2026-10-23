# Reciprocity (man/reciprocity.Rd).
reciprocity <- function(x, mode = c("default", "ratio", "dyadic", "correlation")) {
  check_netframe(x)
  mode <- match.arg(mode)
  code <- switch(mode, default = 0L, ratio = 1L, dyadic = 2L, correlation = 3L)
  ptr <- engine_graph(x)
  engine_call("reciprocity", .Call(C_reciprocity, ptr, code))
}

# Transitivity and clustering coefficients (man/transitivity.Rd).
transitivity <- function(x, convention = c("igraph", "networkx"), threads = 0L) {
  check_netframe(x)
  convention <- match.arg(convention)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("transitivity",
              .Call(C_transitivity, ptr, undefined_code(convention), threads))
}

# See man/transitivity.Rd.
clustering <- function(x, weighted = FALSE, method = NULL,
                       convention = c("igraph", "networkx"), threads = 0L) {
  clustering_common(x, FALSE, weighted, method, match.arg(convention), threads)
}

# See man/transitivity.Rd.
average_clustering <- function(x, weighted = FALSE, method = NULL,
                               convention = c("igraph", "networkx"), threads = 0L) {
  clustering_common(x, TRUE, weighted, method, match.arg(convention), threads)
}

clustering_common <- function(x, average, weighted, method, convention, threads) {
  check_netframe(x)
  weighted <- check_flag(weighted, "weighted")
  weighting <- 0L
  if (weighted) {
    require_weights(x)
    if (is.null(method)) method <- if (convention == "networkx") "onnela" else "barrat"
    method <- check_choice(method, c("barrat", "onnela"), "method")
    weighting <- if (method == "barrat") 1L else 2L
  }
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call(if (average) "average_clustering" else "clustering",
              .Call(C_clustering, ptr, average, weighting, undefined_code(convention),
                    threads))
}

# The engine's code for an undefined ratio: NaN (igraph) or 0 (NetworkX).
undefined_code <- function(convention) if (convention == "networkx") 1L else 0L

# Shortest-path distances, mean distance, diameter, eccentricity, and radius (man/distances.Rd).
distances <- function(x, from = NULL, to = NULL, mode = c("out", "in", "all"),
                      weighted = NULL, weights_as = NULL, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  w <- resolve_weighting(x, weighted, weights_as)
  rows <- node_positions(x, from, "from")
  cols <- node_positions(x, to, "to")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  d <- engine_call("distances",
                   .Call(C_distances, ptr, if (is.null(rows)) NULL else rows - 1L,
                         mode_code(mode), w$weighted, w$code, threads))
  labels <- node_ids(x)
  dimnames(d) <- list(if (is.null(rows)) labels else labels[rows], labels)
  if (!is.null(cols)) d <- d[, cols, drop = FALSE]
  d
}

# See man/distances.Rd.
mean_distance <- function(x, mode = c("out", "in", "all"), weighted = NULL,
                          weights_as = NULL, unconnected = TRUE, threads = 0L) {
  path_summary(x, NULL, match.arg(mode), weighted, weights_as, unconnected, threads,
               "mean_distance")[[6L]]
}

# See man/distances.Rd.
diameter <- function(x, mode = c("out", "in", "all"), weighted = NULL,
                     weights_as = NULL, unconnected = TRUE, threads = 0L) {
  path_summary(x, NULL, match.arg(mode), weighted, weights_as, unconnected, threads,
               "diameter")[[7L]]
}

# See man/distances.Rd.
eccentricity <- function(x, mode = c("all", "out", "in"), weighted = NULL,
                         weights_as = NULL, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  w <- resolve_weighting(x, weighted, weights_as)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("eccentricity",
              .Call(C_eccentricity, ptr, mode_code(mode), w$weighted, w$code, threads))
}

# See man/distances.Rd.
radius <- function(x, mode = c("all", "out", "in"), weighted = NULL,
                   weights_as = NULL, threads = 0L) {
  e <- eccentricity(x, match.arg(mode), weighted, weights_as, threads)
  if (length(e) == 0L) NaN else min(e)
}

# Totals over the searches from `sources` (0-based positions, or NULL for
# every node): sources, pairs, total, longest, incomplete, mean, diameter.
path_summary <- function(x, sources, mode, weighted, weights_as, unconnected, threads, fn) {
  check_netframe(x)
  w <- resolve_weighting(x, weighted, weights_as)
  unconnected <- check_flag(unconnected, "unconnected")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call(fn, .Call(C_path_summary, ptr, sources, mode_code(mode), w$weighted, w$code,
                        unconnected, threads))
}

# Reachability (man/reachability.Rd).
reachability <- function(x, mode = c("out", "in", "all"), threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("reachability", .Call(C_reachability, ptr, mode_code(mode), threads))
}

# Cut points and bridges (man/cut_points.Rd).
cut_points <- function(x) {
  check_netframe(x)
  ptr <- engine_graph(x)
  r <- engine_call("cut_points", .Call(C_cut_structure, ptr))
  x$nodes[[x$roles$id]][r[[1L]]]
}

# See man/cut_points.Rd.
bridges <- function(x) {
  check_netframe(x)
  ptr <- engine_graph(x)
  which(engine_call("bridges", .Call(C_cut_structure, ptr))[[2L]])
}

# Summarize a network (man/summary.netframe.Rd).
summary.netframe <- function(object, paths = NULL, weighted = NULL, threads = 0L, ...) {
  x <- object
  check_netframe(x)
  check_flag_or_null(paths, "paths")
  w <- resolve_weighting(x, weighted, NULL)
  threads <- resolve_threads(threads)
  n <- nrow(x$nodes)
  ptr <- engine_graph(x)
  counts <- engine_call("summary", .Call(C_edge_counts, ptr))
  comp <- components(x, "weak")
  largest <- if (comp$no > 0L) which.max(comp$csize) else NA_integer_
  largest_size <- if (comp$no > 0L) comp$csize[[largest]] else 0L
  if (is.null(paths)) paths <- largest_size <= summary_path_limit
  mean_dist <- diam <- NA_real_
  if (paths && largest_size > 0L) {
    src <- which(comp$membership == largest) - 1L
    s <- engine_call("summary", .Call(C_path_summary, ptr, src, mode_code("out"), w$weighted,
                                      w$code, TRUE, threads))
    mean_dist <- s[[6L]]
    diam <- s[[7L]]
  }
  deg <- degree(x)
  deg <- if (n > 0L) {
    c(mean = mean(deg), min = min(deg), max = max(deg))
  } else {
    c(mean = NaN, min = NaN, max = NaN)
  }
  structure(
    list(
      nodes = n,
      edges = nrow(x$edges),
      directed = x$directed,
      weight = x$roles$weight,
      weight_type = x$roles$weight_type,
      loops = counts[[1L]],
      multiple = counts[[2L]],
      isolates = counts[[3L]],
      density = graph_density(x),
      components = comp$no,
      largest_component = largest_size,
      reciprocity = if (x$directed) reciprocity(x) else NA_real_,
      transitivity = transitivity(x, threads = threads),
      paths = paths,
      path_weighted = w$weighted,
      mean_distance = mean_dist,
      diameter = diam,
      degree = deg,
      centralization = centralization(x, "degree", threads = threads)
    ),
    class = "summary.netframe"
  )
}

# Largest component for which summary() computes distances by default.
summary_path_limit <- 5000L

# See man/summary.netframe.Rd.
print.summary.netframe <- function(x, digits = 3L, ...) {
  num <- function(v) format(v, digits = digits)
  int <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
  weighting <- if (is.null(x$weight)) {
    "no"
  } else {
    sprintf("yes, by `%s` (%ss)", x$weight, x$weight_type)
  }
  components <- if (x$components > 1L) {
    sprintf("%s (largest has %s nodes)", int(x$components), int(x$largest_component))
  } else {
    int(x$components)
  }
  rows <- c(
    "Nodes" = int(x$nodes),
    "Edges" = int(x$edges),
    "Directed" = if (x$directed) "yes" else "no",
    "Weighted" = weighting,
    "Self-loops" = int(x$loops),
    "Multiple edges" = int(x$multiple),
    "Isolates" = int(x$isolates),
    "Density" = num(x$density),
    "Components" = components
  )
  if (x$directed) rows <- c(rows, "Reciprocity" = num(x$reciprocity))
  rows <- c(rows, "Transitivity" = num(x$transitivity))
  if (x$paths) {
    rows <- c(rows, "Mean distance" = num(x$mean_distance), "Diameter" = num(x$diameter))
  } else {
    rows <- c(rows, "Mean distance" = "not computed", "Diameter" = "not computed")
  }
  rows <- c(rows,
    "Degree (mean, min, max)" = paste(num(x$degree[["mean"]]), num(x$degree[["min"]]),
                                      num(x$degree[["max"]]), sep = ", "),
    "Degree centralization" = num(x$centralization)
  )
  cat("<netframe summary>\n")
  cat(paste0(format(names(rows)), "  ", rows, "\n"), sep = "")
  notes <- character(0)
  if (x$paths && x$components > 1L) {
    notes <- c(notes, sprintf(
      "Mean distance and diameter are for the largest component (%s of %s nodes).",
      int(x$largest_component), int(x$nodes)))
  }
  if (x$paths && x$directed) {
    notes <- c(notes, "Distances follow edge directions; unreachable pairs are left out.")
  }
  if (x$paths && isTRUE(x$path_weighted)) {
    notes <- c(notes, sprintf("Distances use the weights as %ss.", x$weight_type))
  }
  if (!x$paths) {
    notes <- c(notes, sprintf(paste0(
      "Mean distance and diameter were skipped because the largest component has %s ",
      "nodes; use summary(x, paths = TRUE) to compute them."), int(x$largest_component)))
  }
  if (length(notes)) cat(paste0("Note: ", notes, "\n"), sep = "")
  invisible(x)
}
