#' Read and write network files
#'
#' `read_network()` reads a file into a [netframe()]; `write_network()`
#' writes one out. The format is taken from the file extension unless
#' `format` is given.
#'
#' | Format | Extensions | Read | Write | Common source |
#' |---|---|---|---|---|
#' | GraphML | `.graphml`, `.xml` | yes | yes | igraph, NetworkX, Gephi, yEd, Cytoscape |
#' | GEXF | `.gexf` | yes | yes | Gephi |
#' | GML | `.gml` | yes | yes | igraph, NetworkX, Cytoscape |
#' | Pajek | `.net`, `.paj` | yes | yes | Pajek, igraph, NetworkX |
#' | UCINET DL | `.dl` | yes | yes | UCINET, igraph |
#' | JSON node-link | `.json` | yes | yes | NetworkX, D3 |
#' | Edge list | `.csv`, `.tsv` | yes | yes | spreadsheets, any tool |
#' | DOT | `.dot`, `.gv` | no | yes | Graphviz |
#'
#' Attributes keep their types where the format records them (GraphML,
#' GEXF, JSON) and are inferred otherwise. Node identifiers that are all
#' whole numbers are returned as integers. Formats differ in what they can
#' hold: Pajek and UCINET DL store only structure, node labels, and a
#' single edge value, and DOT is for drawing.
#'
#' @param file Path to the file.
#' @param format Format name (`"graphml"`, `"gexf"`, `"gml"`, `"pajek"`,
#'   `"dl"`, `"json"`, `"csv"`, `"tsv"`, or `"dot"`), or `NULL` to use the
#'   file extension.
#' @param directed For formats that do not record direction (UCINET DL and
#'   edge lists), whether edges are directed. Defaults to `TRUE`.
#' @param nodes For edge lists only: a data frame of nodes, or the path to
#'   a CSV file of nodes whose first column is the node identifier.
#' @param x A [netframe()].
#' @param weight When reading, the edge column holding weights; by default
#'   a numeric column called `weight` is used if there is one. When writing,
#'   the edge column written as the link value in formats with a single
#'   positional value (Pajek, UCINET DL) and as the GEXF edge weight;
#'   defaults to the network's weight column.
#' @param weight_type When reading, what the weights mean (`"distance"` or
#'   `"strength"`); see [netframe()].
#' @param nodes_file For edge lists only: optional path for a CSV file of
#'   nodes, so isolates and node attributes are kept.
#'
#' @return `read_network()` returns a [netframe()]. `write_network()`
#'   returns `file`, invisibly.
#'
#' @examples
#' g <- netframe(data.frame(from = c("a", "b"), to = c("b", "c"), w = c(1, 2)),
#'               weight = "w", directed = FALSE)
#' path <- tempfile(fileext = ".graphml")
#' write_network(g, path)
#' read_network(path)
#' @export
read_network <- function(file, format = NULL, directed = NULL, nodes = NULL,
                         weight = NULL, weight_type = c("distance", "strength")) {
  weight_type <- match.arg(weight_type)
  if (!is.character(file) || length(file) != 1L) {
    stop("`file` must be a single path.", call. = FALSE)
  }
  if (!file.exists(file)) {
    stop(sprintf("File `%s` does not exist.", file), call. = FALSE)
  }
  fmt <- resolve_format(format, file, write = FALSE)

  if (fmt %in% c("csv", "tsv")) {
    sep <- if (fmt == "tsv") "\t" else ","
    edges <- utils::read.csv(file, sep = sep, stringsAsFactors = FALSE,
                             check.names = FALSE, encoding = "UTF-8")
    if (is.character(nodes)) {
      nodes <- utils::read.csv(nodes, stringsAsFactors = FALSE, check.names = FALSE,
                               encoding = "UTF-8")
    }
    return(netframe(edges, nodes, directed = if (is.null(directed)) TRUE else directed,
                    weight = pick_weight(edges, weight), weight_type = weight_type))
  }

  res <- .Call(C_read_network, path.expand(file), fmt,
               if (is.null(directed)) NULL else isTRUE(directed))
  node_df <- as_frame(res$nodes)
  edge_df <- as_frame(res$edges)
  ids <- integer_ids(node_df$id, edge_df$from, edge_df$to)
  node_df$id <- ids[[1]]
  edge_df$from <- ids[[2]]
  edge_df$to <- ids[[3]]
  netframe(edge_df, node_df, directed = res$directed, from = "from", to = "to",
           id = "id", weight = pick_weight(edge_df, weight), weight_type = weight_type,
           graph = res$graph)
}

pick_weight <- function(edges, weight) {
  if (!is.null(weight)) {
    if (!weight %in% names(edges)) {
      stop(sprintf("Edge column `%s` not found in the file.", weight), call. = FALSE)
    }
    return(weight)
  }
  if ("weight" %in% names(edges) && is.numeric(edges$weight)) "weight"
}

#' @rdname read_network
#' @export
write_network <- function(x, file, format = NULL, weight = NULL, nodes_file = NULL) {
  check_netframe(x)
  if (!is.character(file) || length(file) != 1L) {
    stop("`file` must be a single path.", call. = FALSE)
  }
  fmt <- resolve_format(format, file, write = TRUE)
  weight <- if (is.null(weight)) x$roles$weight else weight
  if (!is.null(weight) && !weight %in% names(x$edges)) {
    stop(sprintf("Edge column `%s` not found.", weight), call. = FALSE)
  }
  # GEXF stores the edge weight in an attribute that must be called
  # `weight`; every other format keeps the column's own name.
  tables <- export_tables(x, if (fmt == "gexf") weight)

  if (fmt %in% c("csv", "tsv")) {
    sep <- if (fmt == "tsv") "\t" else ","
    utils::write.table(tables$edges, file, sep = sep, row.names = FALSE,
                       qmethod = "double", fileEncoding = "UTF-8")
    if (!is.null(nodes_file)) {
      utils::write.table(tables$nodes, nodes_file, sep = sep, row.names = FALSE,
                         qmethod = "double", fileEncoding = "UTF-8")
    }
    return(invisible(file))
  }

  weight_col <- if (is.null(weight)) NULL else if (fmt == "gexf") "weight" else weight
  .Call(C_write_network, x$directed, tables$graph, as.list(tables$nodes),
        as.list(tables$edges), path.expand(file), fmt, weight_col)
  invisible(file)
}

format_names <- c(
  graphml = "graphml", xml = "graphml", gexf = "gexf", gml = "gml",
  pajek = "pajek", net = "pajek", paj = "pajek", dl = "dl", ucinet = "dl",
  json = "json", csv = "csv", tsv = "tsv", txt = "tsv", dot = "dot", gv = "dot"
)

resolve_format <- function(format, file, write) {
  key <- if (is.null(format)) tolower(tools::file_ext(file)) else tolower(format)
  fmt <- unname(format_names[key])
  if (is.na(fmt) || !nzchar(key)) {
    stop(sprintf(
      "Cannot tell the format of `%s`. Use `format =` with one of: %s.",
      file, paste(unique(format_names), collapse = ", ")
    ), call. = FALSE)
  }
  if (!write && fmt == "dot") {
    stop("DOT files can be written but not read.", call. = FALSE)
  }
  fmt
}

as_frame <- function(cols) {
  df <- as.data.frame(cols, stringsAsFactors = FALSE, optional = TRUE)
  names(df) <- names(cols)
  if (!length(cols)) df <- data.frame()
  df
}

# Identifiers that are all whole numbers become integers, consistently for
# nodes and edge endpoints so they still match.
integer_ids <- function(...) {
  ids <- list(...)
  pattern <- "^-?(0|[1-9][0-9]{0,8})$"
  ok <- all(vapply(ids, function(v) is.character(v) && all(grepl(pattern, v)), TRUE))
  if (ok) lapply(ids, as.integer) else ids
}

# Prepares node and edge tables for the engine: structural columns first
# under the names the engine expects, the weight column renamed to `weight`
# when `weight` is given, and R-specific types converted to plain vectors.
export_tables <- function(x, weight) {
  nodes <- x$nodes
  edges <- x$edges
  id <- x$roles$id
  node_out <- c(list(id = as.character(nodes[[id]])),
                export_columns(nodes[setdiff(names(nodes), id)], reserved = "id"))
  attrs <- edges[setdiff(names(edges), c(x$roles$from, x$roles$to))]
  if (!is.null(weight) && weight != "weight") {
    if ("weight" %in% names(attrs)) {
      names(attrs)[names(attrs) == "weight"] <- "weight_attr"
    }
    names(attrs)[names(attrs) == weight] <- "weight"
  }
  edge_out <- c(list(from = as.character(edges[[x$roles$from]]),
                     to = as.character(edges[[x$roles$to]])),
                export_columns(attrs, reserved = c("from", "to")))
  graph <- Filter(function(v) is.atomic(v) && length(v) == 1L, x$graph)
  graph <- lapply(graph, function(v) export_columns(list(v = v), character())[[1]])
  list(graph = graph,
       nodes = as.data.frame(node_out, stringsAsFactors = FALSE, optional = TRUE),
       edges = as.data.frame(edge_out, stringsAsFactors = FALSE, optional = TRUE))
}

export_columns <- function(cols, reserved) {
  out <- lapply(names(cols), function(name) {
    v <- cols[[name]]
    if (is.factor(v)) v <- as.character(v)
    if (inherits(v, "Date")) v <- format(v, "%Y-%m-%d")
    if (inherits(v, "POSIXt")) v <- format(v, "%Y-%m-%dT%H:%M:%S%z")
    if (is.list(v)) {
      stop(sprintf("Column `%s` is a list; convert it to text or numbers before writing.",
                   name), call. = FALSE)
    }
    if (!is.logical(v) && !is.numeric(v) && !is.character(v)) v <- as.character(v)
    if (is.numeric(v) && !is.integer(v)) v <- as.double(v)
    v
  })
  nm <- names(cols)
  nm[nm %in% reserved] <- paste0(nm[nm %in% reserved], "_attr")
  stats::setNames(out, nm)
}
