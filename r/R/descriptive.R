# Reciprocity (man/reciprocity.Rd).
reciprocity <- function(x, mode = c("default", "ratio", "dyadic", "correlation")) {
  check_netframe(x)
  mode <- match.arg(mode)
  code <- switch(mode, default = 0L, ratio = 1L, dyadic = 2L, correlation = 3L)
  .Call(C_reciprocity, engine_graph(x), code)
}

# Transitivity and clustering coefficients (man/transitivity.Rd).
transitivity <- function(x, convention = c("igraph", "networkx"), threads = 0L) {
  check_netframe(x)
  convention <- match.arg(convention)
  .Call(C_transitivity, engine_graph(x), undefined_code(convention), resolve_threads(threads))
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
  weighting <- 0L
  if (isTRUE(weighted)) {
    require_weights(x)
    if (is.null(method)) method <- if (convention == "networkx") "onnela" else "barrat"
    method <- match.arg(method, c("barrat", "onnela"))
    weighting <- if (method == "barrat") 1L else 2L
  }
  .Call(C_clustering, engine_graph(x), average, weighting, undefined_code(convention),
        resolve_threads(threads))
}

undefined_code <- function(convention) if (convention == "networkx") 1L else 0L

# Shortest-path distances, mean distance, diameter, eccentricity, and radius (man/distances.Rd).
distances <- function(x, from = NULL, to = NULL, mode = c("out", "in", "all"),
                      weighted = NULL, weights_as = NULL, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  w <- resolve_weighting(x, weighted, weights_as)
  ids <- x$nodes[[x$roles$id]]
  rows <- node_positions(x, from, "from")
  cols <- node_positions(x, to, "to")
  d <- .Call(C_distances, engine_graph(x), if (is.null(rows)) NULL else rows - 1L,
             mode_code(mode), w$weighted, w$code, resolve_threads(threads))
  labels <- as.character(ids)
  dimnames(d) <- list(if (is.null(rows)) labels else labels[rows], labels)
  if (!is.null(cols)) d <- d[, cols, drop = FALSE]
  d
}

# See man/distances.Rd.
mean_distance <- function(x, mode = c("out", "in", "all"), weighted = NULL,
                          weights_as = NULL, unconnected = TRUE, threads = 0L) {
  path_summary(x, NULL, match.arg(mode), weighted, weights_as, unconnected, threads)[[6L]]
}

# See man/distances.Rd.
diameter <- function(x, mode = c("out", "in", "all"), weighted = NULL,
                     weights_as = NULL, unconnected = TRUE, threads = 0L) {
  path_summary(x, NULL, match.arg(mode), weighted, weights_as, unconnected, threads)[[7L]]
}

# See man/distances.Rd.
eccentricity <- function(x, mode = c("all", "out", "in"), weighted = NULL,
                         weights_as = NULL, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  w <- resolve_weighting(x, weighted, weights_as)
  .Call(C_eccentricity, engine_graph(x), mode_code(mode), w$weighted, w$code,
        resolve_threads(threads))
}

# See man/distances.Rd.
radius <- function(x, mode = c("all", "out", "in"), weighted = NULL,
                   weights_as = NULL, threads = 0L) {
  e <- eccentricity(x, match.arg(mode), weighted, weights_as, threads)
  if (length(e) == 0L) NaN else min(e)
}

# Totals over the searches from `sources` (0-based positions, or NULL for
# every node): sources, pairs, total, longest, incomplete, mean, diameter.
path_summary <- function(x, sources, mode, weighted, weights_as, unconnected, threads) {
  check_netframe(x)
  w <- resolve_weighting(x, weighted, weights_as)
  if (!is.logical(unconnected) || length(unconnected) != 1L || is.na(unconnected)) {
    stop("`unconnected` must be TRUE or FALSE.", call. = FALSE)
  }
  .Call(C_path_summary, engine_graph(x), sources, mode_code(mode), w$weighted, w$code,
        unconnected, resolve_threads(threads))
}

# Positions (1-based) of node identifiers in x$nodes, or NULL for NULL.
# Shared by every function that takes node identifiers.
node_positions <- function(x, nodes, what = "nodes") {
  if (is.null(nodes)) return(NULL)
  ids <- x$nodes[[x$roles$id]]
  if (is.factor(ids)) ids <- as.character(ids)
  if (is.factor(nodes)) nodes <- as.character(nodes)
  pos <- match(nodes, ids)
  if (anyNA(pos)) {
    shown <- utils::head(unique(nodes[is.na(pos)]), 5L)
    stop(sprintf("%d node(s) in `%s` are not in the network, for example: %s.",
                 length(unique(nodes[is.na(pos)])), what,
                 paste(format(shown), collapse = ", ")), call. = FALSE)
  }
  pos
}

# Reachability (man/reachability.Rd).
reachability <- function(x, mode = c("out", "in", "all"), threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  .Call(C_reachability, engine_graph(x), mode_code(mode), resolve_threads(threads))
}

# Cut points and bridges (man/cut_points.Rd).
cut_points <- function(x) {
  check_netframe(x)
  r <- .Call(C_cut_structure, engine_graph(x))
  x$nodes[[x$roles$id]][r[[1L]]]
}

# See man/cut_points.Rd.
bridges <- function(x) {
  check_netframe(x)
  which(.Call(C_cut_structure, engine_graph(x))[[2L]])
}

# Summarize a network (man/summary.netframe.Rd).
summary.netframe <- function(object, paths = NULL, weighted = NULL, threads = 0L, ...) {
  x <- object
  check_netframe(x)
  n <- nrow(x$nodes)
  ptr <- engine_graph(x)
  counts <- .Call(C_edge_counts, ptr)
  comp <- components(x, "weak")
  largest <- if (comp$no > 0L) which.max(comp$csize) else NA_integer_
  largest_size <- if (comp$no > 0L) comp$csize[[largest]] else 0L
  w <- resolve_weighting(x, weighted, NULL)
  if (is.null(paths)) paths <- largest_size <= summary_path_limit
  if (!is.logical(paths) || length(paths) != 1L || is.na(paths)) {
    stop("`paths` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }
  mean_dist <- diam <- NA_real_
  if (paths && largest_size > 0L) {
    src <- which(comp$membership == largest) - 1L
    s <- .Call(C_path_summary, ptr, src, mode_code("out"), w$weighted, w$code, TRUE,
               resolve_threads(threads))
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
