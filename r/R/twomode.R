# Two-mode (bipartite) networks: construction, incidence matrices,
# projections, and two-mode centrality and density.

# Two-mode networks (man/twomode.Rd).
twomode <- function(edges, nodes = NULL, type = "type", first = NULL, directed = FALSE,
                    from = 1L, to = 2L, id = 1L, weight = NULL,
                    weight_type = c("distance", "strength"), graph = list()) {
  weight_type <- match.arg(weight_type)
  if (!is.data.frame(edges)) {
    stop("`edges` must be a data frame.", call. = FALSE)
  }
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("`type` must be a single column name.", call. = FALSE)
  }
  if (is.null(nodes)) {
    if (type == "name") {
      stop("Without `nodes`, identifiers go in a column called `name`; choose another `type`.",
           call. = FALSE)
    }
    from_col <- resolve_column(edges, from, "from", "edges")
    to_col <- resolve_column(edges, to, "to", "edges")
    a <- unique(as.vector(edges[[from_col]]))
    b <- unique(as.vector(edges[[to_col]]))
    both <- intersect(a, b)
    if (length(both)) {
      stop(sprintf(paste0(
        "%d identifier(s) appear in both `%s` and `%s`, for example %s; in a two-mode ",
        "network every node belongs to one mode. Supply `nodes` with a type column."),
        length(both), from_col, to_col, paste(format(utils::head(both, 3L)), collapse = ", ")),
        call. = FALSE)
    }
    nodes <- data.frame(name = c(a, b), stringsAsFactors = FALSE)
    nodes[[type]] <- rep(c(FALSE, TRUE), c(length(a), length(b)))
    id <- "name"
    if (!is.null(first)) {
      stop("`first` applies only when `nodes` is given; without it the `from` column ",
           "is the first mode.", call. = FALSE)
    }
  }
  x <- netframe(edges, nodes, directed = directed, from = from, to = to, id = id,
                weight = weight, weight_type = weight_type, graph = graph)
  set_twomode(x, type, first)
}

# See man/twomode.Rd.
as_twomode <- function(x, type = NULL, first = NULL, weighted = NULL) {
  if (is.matrix(x) || (isS4(x) && methods::is(x, "Matrix"))) {
    return(incidence_to_twomode(x, weighted))
  }
  if (is.data.frame(x) && !inherits(x, "netframe")) {
    return(twomode(x, type = if (is.null(type)) "type" else type, first = first))
  }
  x <- as_netframe(x)
  if (is.null(type)) {
    if (is_twomode(x)) return(x)
    type <- find_type_column(x$nodes)
  }
  set_twomode(x, type, first)
}

# See man/twomode.Rd.
is_twomode <- function(x) {
  if (!inherits(x, "netframe") || is.null(x$roles$type)) return(FALSE)
  v <- x$nodes[[x$roles$type]]
  is.logical(v) && !anyNA(v)
}

# See man/twomode.Rd.
incidence_matrix <- function(x, attr = NULL, sparse = FALSE) {
  check_twomode(x)
  second <- x$nodes[[x$roles$type]]
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

# One-mode projections of a two-mode network (man/project.Rd).
project <- function(x, which = c("first", "second", "both"), weighted = TRUE,
                    method = c("count", "newman", "product"), loops = FALSE, threads = 0L) {
  check_twomode(x)
  which <- match.arg(which)
  method <- match.arg(method)
  if (which == "both") {
    return(list(
      first = project(x, "first", weighted, method, loops, threads),
      second = project(x, "second", weighted, method, loops, threads)
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
  res <- engine_call("project", .Call(C_project, ptr, types, second, code, loops, threads))
  in_mode <- x$nodes[[x$roles$type]] == second
  nodes <- x$nodes[in_mode, , drop = FALSE]
  rownames(nodes) <- NULL
  ids <- x$nodes[[x$roles$id]]
  edges <- data.frame(from = ids[res[[1]]], to = ids[res[[2]]], stringsAsFactors = FALSE)
  if (weighted) edges$weight <- res[[3]]
  netframe(edges, nodes, directed = FALSE, from = "from", to = "to", id = x$roles$id,
           weight = if (weighted) "weight", weight_type = "strength",
           graph = x$graph)
}

# Centrality and density of two-mode networks (man/twomode_degree.Rd).
twomode_degree <- function(x, normalized = TRUE) {
  twomode_centrality(x, 0L, normalized, 0L, 1L, "twomode_degree")
}

# See man/twomode_degree.Rd.
twomode_closeness <- function(x, normalized = TRUE,
                              unreachable = c("wasserman_faust", "ignore"), threads = 0L) {
  unreachable <- match.arg(unreachable)
  if (unreachable == "wasserman_faust" && !isTRUE(normalized)) {
    stop("The Wasserman-Faust rule applies to normalized closeness; use ",
         "`unreachable = \"ignore\"` with `normalized = FALSE`.", call. = FALSE)
  }
  twomode_centrality(x, 1L, normalized, if (unreachable == "ignore") 0L else 1L, threads,
                     "twomode_closeness")
}

# See man/twomode_degree.Rd.
twomode_betweenness <- function(x, normalized = TRUE, threads = 0L) {
  twomode_centrality(x, 2L, normalized, 0L, threads, "twomode_betweenness")
}

# See man/twomode_degree.Rd.
twomode_eigenvector <- function(x, tol = 1e-12, max_iter = 10000L, threads = 0L) {
  check_twomode(x)
  check_iteration(tol, max_iter)
  threads <- resolve_threads(threads)
  types <- twomode_types(x)
  ptr <- engine_graph(x)
  res <- engine_call("twomode_eigenvector",
                     .Call(C_twomode_eigenvector, ptr, types, as.double(tol),
                           as.double(max_iter), threads))
  if (!res[[3]] && nrow(x$edges) > 0L) {
    warning("The network is not connected, so eigenvector scores are zero outside the ",
            "dominant component(s) and may not be unique.", call. = FALSE)
  }
  res[[1]]
}

# See man/twomode_degree.Rd.
twomode_density <- function(x) {
  check_twomode(x)
  types <- twomode_types(x)
  ptr <- engine_graph(x)
  engine_call("twomode_density", .Call(C_twomode_density, ptr, types))
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
  check_netframe(x)
  if (is.null(x$roles$type)) {
    stop("`x` is not a two-mode network; create one with twomode() or as_twomode().",
         call. = FALSE)
  }
  v <- x$nodes[[x$roles$type]]
  if (!is.logical(v) || anyNA(v)) {
    stop(sprintf(paste0("The node type column `%s` must be logical without missing values; ",
                        "recreate the network with as_twomode()."), x$roles$type),
         call. = FALSE)
  }
}

# Node types as the engine expects them: 0 for the first mode, 1 for the
# second.
twomode_types <- function(x) as.integer(x$nodes[[x$roles$type]])

# Records `type` (a node column) as the two-mode node type, converting a
# two-valued column to a new logical column, and checks every edge.
set_twomode <- function(x, type, first) {
  if (!is.character(type) || length(type) != 1L || !type %in% names(x$nodes)) {
    stop(sprintf("Node column `%s` not found.", format(type)), call. = FALSE)
  }
  v <- x$nodes[[type]]
  if (anyNA(v)) {
    stop(sprintf("Node type column `%s` has missing values.", type), call. = FALSE)
  }
  if (is.logical(v)) {
    if (!is.null(first)) {
      stop("`first` applies only to a type column that is not logical; a logical ",
           "column always has FALSE for the first mode.", call. = FALSE)
    }
    col <- type
    labels <- c("FALSE", "TRUE")
    is_second <- v
  } else {
    # Radix sorting orders text by bytes, the same in every locale.
    values <- if (is.factor(v)) levels(droplevels(v)) else sort(unique(v), method = "radix")
    if (length(values) != 2L) {
      stop(sprintf(paste0("Node type column `%s` must be logical or have exactly two ",
                          "distinct values; it has %d."), type, length(values)),
           call. = FALSE)
    }
    key <- if (is.factor(v)) as.character(v) else v
    if (!is.null(first)) {
      if (length(first) != 1L || !as.character(first) %in% as.character(values)) {
        stop(sprintf("`first` must be one of the values of `%s`: %s.", type,
                     paste(format(values), collapse = ", ")), call. = FALSE)
      }
      first_value <- values[as.character(values) == as.character(first)]
    } else {
      first_value <- values[[1L]]
    }
    second_value <- values[as.character(values) != as.character(first_value)]
    is_second <- as.character(key) != as.character(first_value)
    labels <- c(as.character(first_value), as.character(second_value))
    col <- if ("type" %in% names(x$nodes)) {
      make.unique(c(names(x$nodes), make.names(paste0("is_", labels[[2L]]))))[ncol(x$nodes) + 1L]
    } else {
      "type"
    }
    x$nodes[[col]] <- is_second
  }
  bad <- which(is_second[x$index$src + 1L] == is_second[x$index$tgt + 1L])
  if (length(bad)) {
    k <- bad[[1L]]
    ids <- x$nodes[[x$roles$id]]
    a <- x$index$src[[k]] + 1L
    stop(sprintf(paste0(
      "%d edge(s) join two nodes of the same mode, for example edge %d between `%s` and ",
      "`%s` (both %s); in a two-mode network every edge joins the two modes."),
      length(bad), k, format(ids[[a]]), format(ids[[x$index$tgt[[k]] + 1L]]),
      if (is_second[[a]]) labels[[2L]] else labels[[1L]]), call. = FALSE)
  }
  x$roles$type <- col
  x$roles$type_labels <- labels
  x
}

# A logical `type` column, else a two-valued `class` (DyNetML node
# classes), `type`, or `mode` (network's bipartite networks) column.
find_type_column <- function(nodes) {
  two <- function(col) {
    col %in% names(nodes) && !anyNA(nodes[[col]]) &&
      length(unique(as.character(nodes[[col]]))) == 2L
  }
  if ("type" %in% names(nodes) && is.logical(nodes$type)) return("type")
  for (col in c("class", "type", "mode")) if (two(col)) return(col)
  stop(paste0("Cannot tell the two modes apart. Name the node column that holds them with ",
              "`type =` (a logical column, or one with exactly two distinct values)."),
       call. = FALSE)
}

# An incidence matrix as a two-mode network: rows are the first mode.
incidence_to_twomode <- function(m, weighted) {
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
  edges <- data.frame(from = rows[i], to = cols[j], stringsAsFactors = FALSE)
  if (weighted) edges$weight <- v
  nodes <- data.frame(name = c(rows, cols), type = rep(c(FALSE, TRUE), c(length(rows),
                                                                         length(cols))),
                      stringsAsFactors = FALSE)
  twomode(edges, nodes, type = "type", from = "from", to = "to", id = "name",
          weight = if (weighted) "weight")
}
