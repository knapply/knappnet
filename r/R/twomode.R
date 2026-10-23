# Two-mode (affiliation) networks: construction through the mode role with
# edges declared to join different modes only, the bipartite test,
# incidence matrices, projections, and two-mode centrality and density.
# A mode is a declared, disjoint set of nodes (wasserman1994social, ch. 8;
# borgatti1997network); being bipartite (2-colourable) is a property of
# the structure (newman2018networks, ch. 6), tested separately.

# Two-mode networks (man/net_twomode.Rd).
net_twomode <- function(edges, nodes = NULL, mode = "mode", first = NULL, directed = FALSE,
                        source = "source", target = "target", id = "id", weight = NULL,
                        weight_type = NULL, quiet = FALSE, graph = list()) {
  check_tables(edges, nodes)
  mode <- check_string(mode, "mode")
  id <- check_string(id, "id")
  src <- resolve_column(edges, source, "source", "edges")
  tgt <- resolve_column(edges, target, "target", "edges")
  if (is.null(nodes)) {
    a <- edges[[src]]
    b <- edges[[tgt]]
    if (is.factor(a) || is.factor(b)) {
      a <- as.character(a)
      b <- as.character(b)
    }
    a <- unique(a)
    b <- unique(b)
    both <- intersect(a, b)
    if (length(both)) {
      stop(sprintf(paste0(
        "%d identifier(s) appear in both `%s` and `%s`, for example %s; in a two-mode ",
        "network every node belongs to one mode. Supply `nodes` with a mode column."),
        length(both), src, tgt, paste(format(utils::head(both, 3L)), collapse = ", ")),
        call. = FALSE)
    }
    if (mode == id) {
      stop("`mode` and `id` must name different columns of the node table.", call. = FALSE)
    }
    nodes <- data.frame(c(a, b), stringsAsFactors = FALSE)
    names(nodes) <- id
    # The modes are labelled by the endpoint columns they came from.
    nodes[[mode]] <- rep(c(src, tgt), c(length(a), length(b)))
    if (!is.null(first)) {
      stop("`first` applies only when `nodes` is given; without it the `source` column ",
           "is the first mode.", call. = FALSE)
    }
  } else {
    id <- resolve_column(nodes, id, "id", "nodes")
    mode <- resolve_mode_column(nodes, mode)
  }
  roles <- list(id = id, source = src, target = tgt,
                weight = resolve_optional(edges, weight, "weight", "edges"),
                weight_type = check_weight_type(weight, weight_type), mode = mode)
  x <- new_net_frame(edges, nodes, check_flag(directed, "directed"), roles,
                    list(between_modes_only = TRUE), check_graph(graph))
  x <- build_net_frame(x, check_flag(quiet, "quiet"), "net_twomode")
  set_first_mode(x, first)
}

# The mode column: `mode` when present, else igraph's logical `type`
# column, which is accepted as an alias.
resolve_mode_column <- function(nodes, mode) {
  if (mode %in% names(nodes)) return(mode)
  if (mode == "mode" && "type" %in% names(nodes) && is.logical(nodes$type)) return("type")
  resolve_column(nodes, mode, "mode", "nodes")
}

# See man/net_twomode.Rd.
net_as_twomode <- function(x, mode = NULL, first = NULL, weighted = NULL, quiet = FALSE) {
  if (is.matrix(x) || (isS4(x) && methods::is(x, "Matrix"))) {
    return(incidence_to_twomode(x, weighted, quiet))
  }
  if (is.data.frame(x) && !inherits(x, "net_frame")) {
    return(net_twomode(x, mode = if (is.null(mode)) "mode" else mode, first = first,
                       quiet = quiet))
  }
  x <- as_net_frame(x)
  if (is.null(mode)) {
    if (net_is_twomode(x)) return(set_first_mode(x, first))
    mode <- if (!is.null(x$roles$mode)) x$roles$mode else "mode"
  }
  mode <- resolve_mode_column(x$nodes, mode)
  declared <- x$declared
  declared$between_modes_only <- TRUE
  out <- rebuild_net_frame(x, x$edges, x$nodes, mode = mode, declared = declared,
                          quiet = quiet, fn = "net_as_twomode")
  set_first_mode(out, first)
}

# See man/net_twomode.Rd.
net_is_twomode <- function(x) {
  if (!inherits(x, "net_frame") || is.null(x$roles$mode)) return(FALSE)
  isTRUE(x$declared$between_modes_only) && length(mode_labels(x)) == 2L
}

# See man/net_twomode.Rd.
net_is_bipartite <- function(x) {
  check_net_frame(x)
  require_resolved(x)
  n <- nrow(x$nodes)
  src <- x$index$src + 1L
  tgt <- x$index$tgt + 1L
  if (any(src == tgt)) return(FALSE)
  # Breadth-first 2-colouring over the undirected adjacency lists.
  adj <- split(c(tgt, src), factor(c(src, tgt), levels = seq_len(n)))
  color <- integer(n)
  for (root in seq_len(n)) {
    if (color[[root]] != 0L) next
    color[[root]] <- 1L
    queue <- root
    while (length(queue)) {
      v <- queue[[1L]]
      queue <- queue[-1L]
      for (u in adj[[v]]) {
        if (color[[u]] == 0L) {
          color[[u]] <- 3L - color[[v]]
          queue <- c(queue, u)
        } else if (color[[u]] == color[[v]]) {
          return(FALSE)
        }
      }
    }
  }
  TRUE
}

# See man/net_twomode.Rd.
net_incidence <- function(x, attr = NULL, sparse = FALSE) {
  check_twomode(x)
  second <- twomode_types(x) == 1L
  m <- as.matrix(x, attr = attr, sparse = sparse)
  rows <- which(!second)
  cols <- which(second)
  # A directed network may store a tie either way; an undirected one is
  # already symmetric.
  if (x$directed) {
    m[rows, cols, drop = FALSE] + t(m[cols, rows, drop = FALSE])
  } else {
    m[rows, cols, drop = FALSE]
  }
}

# One-mode projections of a two-mode network (man/net_project.Rd).
net_project <- function(x, which = c("first", "second", "both"), weighted = TRUE,
                        method = c("count", "newman", "product"), loops = FALSE,
                        threads = 0L) {
  check_twomode(x)
  which <- match.arg(which)
  method <- match.arg(method)
  if (which == "both") {
    return(list(
      first = net_project(x, "first", weighted, method, loops, threads),
      second = net_project(x, "second", weighted, method, loops, threads)
    ))
  }
  weighted <- check_flag(weighted, "weighted")
  loops <- check_flag(loops, "loops")
  if (method == "product") require_weights(x)
  if (method == "newman" && loops) {
    stop("Newman's collaboration weights have no diagonal; use `loops = FALSE`.",
         call. = FALSE)
  }
  second <- which == "second"
  code <- switch(method, count = 0L, newman = 1L, product = 2L)
  threads <- resolve_threads(threads)
  types <- twomode_types(x)
  ptr <- engine_graph(x)
  res <- engine_call("net_project", .Call(C_project, ptr, types, second, code, loops, threads))
  in_mode <- types == as.integer(second)
  nodes <- x$nodes[in_mode, , drop = FALSE]
  rownames(nodes) <- NULL
  ids <- x$nodes[[x$roles$id]]
  edges <- data.frame(source = ids[res[[1]]], target = ids[res[[2]]], stringsAsFactors = FALSE)
  if (weighted) edges$weight <- res[[3]]
  net_frame(edges, nodes, directed = FALSE, id = x$roles$id,
            weight = if (weighted) "weight", weight_type = "strength",
            loops = if (!loops) FALSE, quiet = TRUE, graph = x$graph)
}

# Centrality and density of two-mode networks (man/net_twomode_degree.Rd).
net_twomode_degree <- function(x, normalized = TRUE) {
  twomode_centrality(x, 0L, normalized, 0L, 1L, "net_twomode_degree")
}

# See man/net_twomode_degree.Rd.
net_twomode_closeness <- function(x, normalized = TRUE,
                                  unreachable = c("wasserman_faust", "ignore"), threads = 0L) {
  unreachable <- match.arg(unreachable)
  if (unreachable == "wasserman_faust" && !isTRUE(normalized)) {
    stop("The Wasserman-Faust rule applies to normalized closeness; use ",
         "`unreachable = \"ignore\"` with `normalized = FALSE`.", call. = FALSE)
  }
  twomode_centrality(x, 1L, normalized, if (unreachable == "ignore") 0L else 1L, threads,
                     "net_twomode_closeness")
}

# See man/net_twomode_degree.Rd.
net_twomode_betweenness <- function(x, normalized = TRUE, threads = 0L) {
  twomode_centrality(x, 2L, normalized, 0L, threads, "net_twomode_betweenness")
}

# See man/net_twomode_degree.Rd.
net_twomode_eigen <- function(x, tol = 1e-12, max_iter = 10000L, threads = 0L) {
  check_twomode(x)
  check_iteration(tol, max_iter)
  threads <- resolve_threads(threads)
  types <- twomode_types(x)
  ptr <- engine_graph(x)
  res <- engine_call("net_twomode_eigen",
                     .Call(C_twomode_eigenvector, ptr, types, as.double(tol),
                           as.double(max_iter), threads))
  if (!res[[3]] && nrow(x$edges) > 0L) {
    warning("The network is not connected, so eigenvector scores are zero outside the ",
            "dominant component(s) and may not be unique.", call. = FALSE)
  }
  res[[1]]
}

# See man/net_twomode_degree.Rd.
net_twomode_density <- function(x) {
  check_twomode(x)
  types <- twomode_types(x)
  ptr <- engine_graph(x)
  engine_call("net_twomode_density", .Call(C_twomode_density, ptr, types))
}

twomode_centrality <- function(x, measure, normalized, unreachable, threads, fn) {
  check_twomode(x)
  normalized <- check_flag(normalized, "normalized")
  threads <- resolve_threads(threads)
  types <- twomode_types(x)
  ptr <- engine_graph(x)
  engine_call(fn, .Call(C_twomode_centrality, ptr, types, measure, normalized, unreachable,
                        threads))
}

check_twomode <- function(x) {
  check_net_frame(x)
  if (is.null(x$roles$mode) || !isTRUE(x$declared$between_modes_only)) {
    stop("`x` is not a two-mode network; create one with net_twomode() or net_as_twomode().",
         call. = FALSE)
  }
  k <- length(mode_labels(x))
  if (k != 2L) {
    stop(sprintf(paste0("`x` has %d modes (column `%s` of `x$nodes`); two-mode measures ",
                        "need exactly two."), k, x$roles$mode), call. = FALSE)
  }
  invisible(TRUE)
}

# The mode labels in the engine's order (first appearance in the node
# table).
mode_labels <- function(x) frame_summary_cached(x)$modes

# The mode of each node as the engine's 0-based code by first appearance.
mode_codes <- function(x) {
  engine_graph(x, missing_ok = TRUE)
  engine_call("net_frame", .Call(C_frame_codes, x$cache$frame, 0L, as.double(nrow(x$nodes))))
}

# Node types as the two-mode engine functions expect them: 0 for the first
# mode, 1 for the second.
twomode_types <- function(x) {
  labels <- mode_labels(x)
  first <- x$roles$mode_levels
  first <- if (is.null(first)) labels[[1L]] else first[[1L]]
  as.integer(labels[mode_codes(x) + 1L] != first)
}

# Records which mode is the first (`mode_levels`: the first mode's label,
# then the others in order of appearance).
set_first_mode <- function(x, first) {
  labels <- mode_labels(x)
  if (is.null(first)) {
    if (is.null(x$roles$mode_levels) || !all(x$roles$mode_levels %in% labels)) {
      x$roles$mode_levels <- labels
    }
    return(x)
  }
  if (length(first) != 1L || is.na(first) || !as.character(first) %in% labels) {
    stop(sprintf("`first` must be one of the modes in `%s`: %s.", x$roles$mode,
                 paste(labels, collapse = ", ")), call. = FALSE)
  }
  first <- as.character(first)
  x$roles$mode_levels <- c(first, setdiff(labels, first))
  x
}

# An incidence matrix as a two-mode network: rows are the first mode.
incidence_to_twomode <- function(m, weighted, quiet = FALSE) {
  rows <- rownames(m)
  cols <- colnames(m)
  if (is.null(rows)) rows <- paste0("r", seq_len(nrow(m)))
  if (is.null(cols)) cols <- paste0("c", seq_len(ncol(m)))
  if (anyDuplicated(c(rows, cols))) {
    stop("Row and column names of an incidence matrix must be distinct identifiers.",
         call. = FALSE)
  }
  if (is.matrix(m)) {
    if (!is.numeric(m) && !is.logical(m)) {
      stop("An incidence matrix must be numeric or logical.", call. = FALSE)
    }
    m[is.na(m)] <- 0
    idx <- which(m != 0, arr.ind = TRUE)
    i <- idx[, 1]
    j <- idx[, 2]
    v <- as.numeric(m[idx])
  } else {
    require_pkg("Matrix")
    s <- methods::as(m, "TsparseMatrix")
    i <- s@i + 1L
    j <- s@j + 1L
    v <- if (methods::.hasSlot(s, "x")) as.numeric(s@x) else rep(1, length(i))
    keep <- v != 0
    i <- i[keep]
    j <- j[keep]
    v <- v[keep]
  }
  o <- order(i, j)
  i <- i[o]
  j <- j[o]
  v <- v[o]
  weighted <- if (is.null(weighted)) any(v != 1) else isTRUE(weighted)
  edges <- data.frame(source = rows[i], target = cols[j], stringsAsFactors = FALSE)
  if (weighted) edges$weight <- v
  nodes <- data.frame(id = c(rows, cols),
                      mode = rep(c("row", "column"), c(length(rows), length(cols))),
                      stringsAsFactors = FALSE)
  net_twomode(edges, nodes, mode = "mode", weight = if (weighted) "weight",
              weight_type = if (weighted) "strength", quiet = quiet)
}
