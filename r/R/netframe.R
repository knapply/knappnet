#' Create a network from node and edge data frames
#'
#' A network is a pair of data frames: `nodes`, with one row per node, and
#' `edges`, with one row per edge. Any other columns are attributes and are
#' carried along untouched. Column roles say which columns hold node
#' identifiers, edge endpoints, and (optionally) edge weights.
#'
#' @param edges A data frame of edges.
#' @param nodes A data frame of nodes, or `NULL` to derive the node set from
#'   the edge endpoints in order of first appearance.
#' @param directed Whether edges are directed.
#' @param from,to Columns of `edges` holding the endpoints, by name or
#'   position.
#' @param id Column of `nodes` holding node identifiers, by name or position.
#' @param weight Optional column of `edges` holding edge weights.
#' @param graph A named list of network-level attributes (for example a
#'   title or data source), carried through conversions and file formats
#'   that support them.
#'
#' @return An object of class `netframe`.
#' @export
netframe <- function(edges, nodes = NULL, directed = TRUE,
                     from = 1L, to = 2L, id = 1L, weight = NULL, graph = list()) {
  if (!is.data.frame(edges)) {
    stop("`edges` must be a data frame.", call. = FALSE)
  }
  from <- resolve_column(edges, from, "from", "edges")
  to <- resolve_column(edges, to, "to", "edges")
  if (!is.null(weight)) {
    weight <- resolve_column(edges, weight, "weight", "edges")
    if (!is.numeric(edges[[weight]])) {
      stop(sprintf("Weight column `%s` must be numeric.", weight), call. = FALSE)
    }
  }
  if (anyNA(edges[[from]]) || anyNA(edges[[to]])) {
    stop("Edge endpoints must not contain missing values.", call. = FALSE)
  }

  if (is.null(nodes)) {
    ids <- unique(c(as.vector(edges[[from]]), as.vector(edges[[to]])))
    nodes <- data.frame(name = ids, stringsAsFactors = FALSE)
    id <- "name"
  } else {
    if (!is.data.frame(nodes)) {
      stop("`nodes` must be a data frame or NULL.", call. = FALSE)
    }
    id <- resolve_column(nodes, id, "id", "nodes")
    if (anyNA(nodes[[id]])) {
      stop("Node identifiers must not contain missing values.", call. = FALSE)
    }
    if (anyDuplicated(nodes[[id]])) {
      dup <- nodes[[id]][anyDuplicated(nodes[[id]])]
      stop(sprintf("Node identifiers must be unique; `%s` appears more than once.",
                   format(dup)), call. = FALSE)
    }
  }

  src <- match(edges[[from]], nodes[[id]])
  tgt <- match(edges[[to]], nodes[[id]])
  missing <- c(edges[[from]][is.na(src)], edges[[to]][is.na(tgt)])
  if (length(missing)) {
    shown <- utils::head(unique(missing), 5L)
    stop(sprintf(
      "%d edge endpoint(s) are not in `nodes`, for example: %s.",
      length(unique(missing)), paste(format(shown), collapse = ", ")
    ), call. = FALSE)
  }

  structure(
    list(
      nodes = nodes,
      edges = edges,
      directed = isTRUE(directed),
      roles = list(id = id, from = from, to = to, weight = weight),
      graph = as.list(graph),
      index = list(src = src - 1L, tgt = tgt - 1L),
      cache = new.env(parent = emptyenv())
    ),
    class = "netframe"
  )
}

#' @export
print.netframe <- function(x, ...) {
  cat(sprintf(
    "<netframe> %s network with %d nodes and %d edges%s\n",
    if (x$directed) "directed" else "undirected",
    nrow(x$nodes), nrow(x$edges),
    if (is.null(x$roles$weight)) "" else sprintf(" (weights: `%s`)", x$roles$weight)
  ))
  invisible(x)
}

resolve_column <- function(df, col, role, what) {
  if (is.numeric(col) && length(col) == 1L && col >= 1 && col <= ncol(df)) {
    return(names(df)[[col]])
  }
  if (is.character(col) && length(col) == 1L && col %in% names(df)) {
    return(col)
  }
  stop(sprintf("Column for `%s` not found in `%s`.", role, what), call. = FALSE)
}

# Returns the engine graph for `x`, building it on first use. The handle is
# an external pointer, which does not survive saveRDS(); a dead handle is
# rebuilt transparently.
engine_graph <- function(x) {
  ptr <- x$cache$ptr
  if (is.null(ptr) || !.Call(C_graph_alive, ptr)) {
    w <- if (is.null(x$roles$weight)) NULL else as.double(x$edges[[x$roles$weight]])
    ptr <- .Call(C_graph_new, nrow(x$nodes), x$index$src, x$index$tgt, w, x$directed)
    x$cache$ptr <- ptr
  }
  ptr
}

check_netframe <- function(x) {
  if (!inherits(x, "netframe")) {
    stop("`x` must be a netframe; create one with netframe().", call. = FALSE)
  }
}
