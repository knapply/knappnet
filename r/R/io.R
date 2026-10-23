# Reading and writing network files. Text formats and CSV are handled
# here; GeoPackage files need sf and live in sf.R.

# Read and write network files (man/read_network.Rd).
read_network <- function(file, format = NULL, directed = NULL, nodes = NULL,
                         weight = NULL, weight_type = c("distance", "strength")) {
  weight_type <- match.arg(weight_type)
  check_path(file, "file")
  if (!file.exists(file)) {
    stop(sprintf("File `%s` does not exist.", file), call. = FALSE)
  }
  check_flag_or_null(directed, "directed")
  check_string(weight, "weight", null_ok = TRUE)
  check_string(format, "format", null_ok = TRUE)
  fmt <- resolve_format(format, file, write = FALSE)
  if (fmt == "gpkg") return(read_gpkg(file, directed, weight, weight_type))

  if (fmt %in% c("csv", "tsv")) {
    sep <- if (fmt == "tsv") "\t" else ","
    edges <- utils::read.csv(file, sep = sep, stringsAsFactors = FALSE,
                             check.names = FALSE, encoding = "UTF-8")
    if (is.character(nodes)) {
      check_path(nodes, "nodes")
      if (!file.exists(nodes)) {
        stop(sprintf("Nodes file `%s` does not exist.", nodes), call. = FALSE)
      }
      nodes <- utils::read.csv(nodes, stringsAsFactors = FALSE, check.names = FALSE,
                               encoding = "UTF-8")
    }
    return(netframe(edges, nodes, directed = if (is.null(directed)) TRUE else directed,
                    weight = pick_weight(edges, weight), weight_type = weight_type))
  }

  # A .xml file may be GraphML or DyNetML; the engine reads its root
  # element to decide.
  if (is.null(format) && identical(tolower(tools::file_ext(file)), "xml")) fmt <- NULL
  path <- path.expand(file)
  dir <- if (is.null(directed)) NULL else directed
  res <- engine_call("read_network", .Call(C_read_network, path, fmt, dir))
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

# See man/read_network.Rd.
write_network <- function(x, file, format = NULL, weight = NULL, nodes_file = NULL) {
  check_netframe(x)
  check_path(file, "file")
  check_string(format, "format", null_ok = TRUE)
  check_string(weight, "weight", null_ok = TRUE)
  if (!is.null(nodes_file)) check_path(nodes_file, "nodes_file")
  fmt <- resolve_format(format, file, write = TRUE)
  weight <- if (is.null(weight)) x$roles$weight else weight
  if (!is.null(weight) && !weight %in% names(x$edges)) {
    stop(sprintf("Edge column `%s` not found.", weight), call. = FALSE)
  }
  if (fmt == "gpkg") return(write_gpkg(x, file, weight))
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
  node_cols <- as.list(tables$nodes)
  edge_cols <- as.list(tables$edges)
  path <- path.expand(file)
  engine_call("write_network", .Call(C_write_network, x$directed, tables$graph, node_cols,
                                     edge_cols, path, fmt, weight_col))
  invisible(file)
}

format_names <- c(
  graphml = "graphml", xml = "graphml", gexf = "gexf", gml = "gml",
  pajek = "pajek", net = "pajek", paj = "pajek", dl = "dl", ucinet = "dl",
  json = "json", dynetml = "dynetml", ora = "dynetml", csv = "csv", tsv = "tsv",
  txt = "tsv", dot = "dot", gv = "dot", gpkg = "gpkg", geopackage = "gpkg"
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
