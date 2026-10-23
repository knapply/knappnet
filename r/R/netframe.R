# The netframe class: the constructor, which checks the tables and builds
# the edge index, and its print method. The argument checks it uses are in
# utils-checks.R, the engine handle in utils-engine.R, and summary() in
# descriptive.R.

# Create a network from node and edge data frames (man/netframe.Rd).
netframe <- function(edges, nodes = NULL, directed = TRUE,
                     from = 1L, to = 2L, id = 1L, weight = NULL,
                     weight_type = c("distance", "strength"), graph = list()) {
  weight_type <- match.arg(weight_type)
  if (!is.data.frame(edges)) {
    stop(sprintf("`edges` must be a data frame, not %s.", describe_class(edges)),
         call. = FALSE)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("`directed` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(graph)) graph <- list()
  if (!is.list(graph) && !is.atomic(graph)) {
    stop("`graph` must be a named list of network-level attributes.", call. = FALSE)
  }
  graph <- as.list(graph)
  if (length(graph) && (is.null(names(graph)) || anyNA(names(graph)) ||
                        !all(nzchar(names(graph))))) {
    stop("Every element of `graph` must be named, as in `graph = list(title = \"...\")`.",
         call. = FALSE)
  }
  from <- resolve_column(edges, from, "from", "edges")
  to <- resolve_column(edges, to, "to", "edges")
  if (!is.null(weight)) {
    weight <- resolve_column(edges, weight, "weight", "edges")
    if (!is.numeric(edges[[weight]])) {
      stop(sprintf("Weight column `%s` must be numeric, not %s.", weight,
                   describe_class(edges[[weight]])), call. = FALSE)
    }
  }
  src_ids <- check_ids(edges[[from]], from, "edges", "Edge endpoints")
  tgt_ids <- check_ids(edges[[to]], to, "edges", "Edge endpoints")

  if (is.null(nodes)) {
    if ((is_classed(src_ids) || is_classed(tgt_ids)) &&
        !identical(class(src_ids), class(tgt_ids))) {
      stop(sprintf(paste0(
        "Columns `%s` and `%s` of `edges` hold different types (%s and %s); ",
        "convert them to the same type first."),
        from, to, describe_class(src_ids), describe_class(tgt_ids)), call. = FALSE)
    }
    ids <- unique(c(src_ids, tgt_ids))
    nodes <- data.frame(name = ids, stringsAsFactors = FALSE)
    id <- "name"
    node_ids <- ids
  } else {
    if (!is.data.frame(nodes)) {
      stop(sprintf("`nodes` must be a data frame or NULL, not %s.", describe_class(nodes)),
           call. = FALSE)
    }
    id <- resolve_column(nodes, id, "id", "nodes")
    node_ids <- check_ids(nodes[[id]], id, "nodes", "Node identifiers")
    if (anyDuplicated(node_ids)) {
      dup <- nodes[[id]][anyDuplicated(node_ids)]
      stop(sprintf(paste0("Node identifiers must be unique; `%s` appears more than once ",
                          "in column `%s` of `nodes`."), format(dup), id), call. = FALSE)
    }
  }

  src <- match(src_ids, node_ids)
  tgt <- match(tgt_ids, node_ids)
  missing <- c(edges[[from]][is.na(src)], edges[[to]][is.na(tgt)])
  if (length(missing)) {
    shown <- utils::head(unique(missing), 5L)
    stop(sprintf(paste0(
      "%d edge endpoint(s) are not in `nodes`, for example: %s. Add them to `nodes`, ",
      "or use `nodes = NULL` to take the node set from the edges."),
      length(unique(missing)), paste(format(shown), collapse = ", ")
    ), call. = FALSE)
  }

  out <- structure(
    list(
      nodes = nodes,
      edges = edges,
      directed = isTRUE(directed),
      roles = list(id = id, from = from, to = to, weight = weight,
                   weight_type = if (!is.null(weight)) weight_type),
      graph = graph,
      index = list(src = src - 1L, tgt = tgt - 1L),
      cache = new.env(parent = emptyenv())
    ),
    class = "netframe"
  )
  # sf node or edge tables also record their geometry columns and CRS.
  record_geometry(out)
}

# See man/netframe.Rd.
print.netframe <- function(x, ...) {
  cat(sprintf(
    "<netframe> %s network with %d nodes and %d edges%s\n",
    if (x$directed) "directed" else "undirected",
    nrow(x$nodes), nrow(x$edges),
    if (is.null(x$roles$weight)) ", unweighted"
    else sprintf(", weighted by `%s` (%ss)", x$roles$weight, x$roles$weight_type)
  ))
  # Only unusual features get a line, so ordinary networks print one line.
  s <- edge_summary(x)
  notes <- c(
    if (s$loops) plural(s$loops, "self-loop"),
    if (s$parallel) plural(s$parallel, "parallel edge")
  )
  if (length(notes)) {
    cat(sprintf("including %s; each measure's help page says how they count\n",
                paste(notes, collapse = " and ")))
  }
  if (!is.null(x$roles$start)) cat(temporal_summary(x))
  cat(spatial_summary(x))
  if (s$bad_weights) {
    cat(sprintf(paste0("%s in weight column `%s` %s missing or not finite; ",
                       "measures will refuse the network until %s fixed\n"),
                plural(s$bad_weights, "value"), x$roles$weight,
                if (s$bad_weights == 1L) "is" else "are",
                if (s$bad_weights == 1L) "it is" else "they are"))
  }
  invisible(x)
}

# Counts of self-loops, parallel edges (edges repeating an earlier edge
# between the same ordered pair of nodes, or unordered pair in undirected
# networks), and missing or non-finite weights.
edge_summary <- function(x) {
  src <- x$index$src
  tgt <- x$index$tgt
  n <- nrow(x$nodes)
  a <- if (x$directed) src else pmin(src, tgt)
  b <- if (x$directed) tgt else pmax(src, tgt)
  # a * n + b is exact in double precision while n^2 < 2^53.
  key <- if (n < 2^26) as.double(a) * n + b else paste(a, b)
  w <- if (is.null(x$roles$weight)) NULL else x$edges[[x$roles$weight]]
  list(loops = sum(src == tgt), parallel = sum(duplicated(key)),
       bad_weights = if (is.null(w)) 0L else sum(!is.finite(w)))
}
