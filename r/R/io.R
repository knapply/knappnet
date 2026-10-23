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
#' | DyNetML | `.xml` | yes | yes | ORA (CASOS, Carnegie Mellon) |
#' | Edge list | `.csv`, `.tsv` | yes | yes | spreadsheets, any tool |
#' | DOT | `.dot`, `.gv` | no | yes | Graphviz |
#' | GeoPackage | `.gpkg` | yes | yes | QGIS, GDAL, sf (spatial networks; needs sf) |
#'
#' Attributes keep their types where the format records them (GraphML,
#' GEXF, JSON) and are inferred otherwise. Node identifiers that are all
#' whole numbers are returned as integers. Formats differ in what they can
#' hold: Pajek and UCINET DL store only structure, node labels, and a
#' single edge value, and DOT is for drawing.
#'
#' Files are read and written as UTF-8, and file paths may contain any
#' character, so identifiers and attributes in any script round-trip on
#' every platform. The exception is writing edge lists, which uses
#' [utils::write.table()]: in R before 4.2 on Windows, and in other locales
#' that are not UTF-8, it can only write characters of the native encoding.
#'
#' GraphML and DyNetML files both use `.xml`. When reading a `.xml` file
#' without `format`, the root element decides (`<graphml>` or a DyNetML
#' `<DynamicNetwork>`, `<DynamicMetaNetwork>`, `<MetaNetwork>`, or
#' `<MetaMatrix>`). Writing a `.xml` file without `format` produces GraphML;
#' use `format = "dynetml"` (or `"ora"`) to write DyNetML.
#'
#' @section GeoPackage:
#' A GeoPackage (`.gpkg`) holds a spatial network (see [to_sf()]) as two
#' layers, `nodes` and `edges`, each written with [sf::st_write()] when it
#' has geometry and as a plain table otherwise, plus a small table layer
#' `network` with the direction, the identifier, endpoint, and weight
#' columns, and network-level attributes (as text). Reading needs only
#' the `nodes` and `edges` layers: without `network`, the first node
#' column is the identifier, the edge columns are `from` and `to`, the
#' network is directed unless `directed` says otherwise, and the weight
#' is found as for other formats. Geometries and their CRS round-trip;
#' attribute types follow GDAL's GeoPackage driver (logical columns come
#' back as integers). Requires the sf package.
#'
#' @section DyNetML:
#' DyNetML (Tsvetovat, Reminga, and Carley 2004) is the XML format of ORA
#' for multimode meta-networks: node classes such as agents, organizations,
#' locations, tasks, and resources, with several networks among them, and
#' optionally one meta-network per time period. It is read into one
#' network:
#'
#' * Nodes get `class` (the node class id) and `type` (its ORA type, such
#'   as `Agent`), `title`, and one column per node property, typed from the
#'   file's property identities (`number` as double; `text`,
#'   `categoryText`, `URI`, and `date` as text).
#' * The same node id may appear in two classes. Only when it does are
#'   identifiers written `class::id`, with the original in `node_id`.
#' * Edges get `network` (the network id) and `weight` (the link value;
#'   links without a value get 1).
#' * Each network says whether it is directed. A file mixing directed and
#'   undirected networks is read as directed, each undirected link stored
#'   as two opposite edges, with a logical `directed` column marking them.
#' * A file with several meta-networks (time periods) gets a `period`
#'   column on nodes and edges, and identifiers `period::id`; each period's
#'   nodes are separate rows.
#'
#' Ports, survey (AnthroPac) data, and properties on networks are not
#' supported and give an error naming the line. Writing produces one
#' meta-network per `period` value, node classes from the `class` and
#' `type` columns (a single `Agent` class without them), and one network
#' per `network` value (or per pair of node classes); other node columns
#' become node properties. Logical node columns are written as ORA
#' category text and read back as text, and integer columns read back as
#' double.
#'
#' @references
#' Tsvetovat, M., Reminga, J., and Carley, K. M. (2004). DyNetML:
#' Interchange format for rich social network data. Technical report
#' CMU-ISRI-04-105, Carnegie Mellon University. \doi{10.1184/R1/6622022.v1}
#'
#' Altman, N., Carley, K. M., and Reminga, J. (2020). ORA User's Guide
#' 2020. Technical report CMU-ISR-20-110, Carnegie Mellon University.
#'
#' @param file Path to the file.
#' @param format Format name (`"graphml"`, `"gexf"`, `"gml"`, `"pajek"`,
#'   `"dl"`, `"json"`, `"dynetml"`, `"csv"`, `"tsv"`, `"gpkg"`, or `"dot"`),
#'   or `NULL` to use the file extension.
#' @param directed For formats that do not record direction (UCINET DL and
#'   edge lists), whether edges are directed. Defaults to `TRUE`.
#' @param nodes For edge lists only: a data frame of nodes, or the path to
#'   a CSV file of nodes whose first column is the node identifier.
#' @param x A [netframe()].
#' @param weight When reading, the edge column holding weights; by default
#'   a numeric column called `weight` is used if there is one. When writing,
#'   the edge column written as the link value in formats with a single
#'   positional value (Pajek, UCINET DL, DyNetML) and as the GEXF edge weight;
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
  check_path(file, "file")
  if (!file.exists(file)) {
    stop(sprintf("File `%s` does not exist.", file), call. = FALSE)
  }
  check_flag_or_null(directed, "directed")
  check_string_or_null(weight, "weight")
  check_string_or_null(format, "format")
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

check_path <- function(value, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop(sprintf("`%s` must be a single file path.", name), call. = FALSE)
  }
}

check_string_or_null <- function(value, name) {
  if (!is.null(value) && (!is.character(value) || length(value) != 1L || is.na(value))) {
    stop(sprintf("`%s` must be a single string or NULL.", name), call. = FALSE)
  }
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
  check_path(file, "file")
  check_string_or_null(format, "format")
  check_string_or_null(weight, "weight")
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
