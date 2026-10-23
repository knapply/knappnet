# Conversion to and from other packages' objects (igraph, tidygraph,
# network, dense and sparse adjacency matrices). Methods for sf and
# sfnetworks are in sf.R. Methods for other packages' generics are
# registered in NAMESPACE with S3method(pkg::generic, netframe) and take
# effect when that package is loaded.

# Convert networks to and from other packages (man/net_from.Rd).
net_from <- function(x, ...) {
  UseMethod("net_from")
}

# See man/net_from.Rd.
net_from.default <- function(x, ...) {
  if (isS4(x) && methods::is(x, "Matrix")) {
    return(net_from.Matrix(x, ...))
  }
  stop(sprintf(paste0(
    "Cannot convert an object of class `%s` to a netframe. Supported: igraph, ",
    "tbl_graph, network, matrix, Matrix, and data.frame (an edge list)."),
    class(x)[[1]]), call. = FALSE)
}

# See man/net_from.Rd.
net_from.netframe <- function(x, ...) x

# See man/net_from.Rd.
net_from.data.frame <- function(x, ...) net_frame(x, ...)

# See man/net_from.Rd.
net_from.igraph <- function(x, weight_type = "distance", quiet = FALSE, ...) {
  require_pkg("igraph")
  n <- igraph::vcount(x)
  vattrs <- igraph::vertex_attr(x)
  named <- "name" %in% names(vattrs)
  ids <- if (named) vattrs$name else seq_len(n)
  vattrs$name <- NULL
  nodes <- data.frame(name = ids, stringsAsFactors = FALSE)
  for (a in names(vattrs)) nodes[[a]] <- plain_column(vattrs[[a]])

  ends <- igraph::ends(x, igraph::E(x), names = FALSE)
  edges <- data.frame(source = ids[ends[, 1]], target = ids[ends[, 2]],
                      stringsAsFactors = FALSE)
  eattrs <- igraph::edge_attr(x)
  for (a in names(eattrs)) {
    col <- if (a %in% c("source", "target")) paste0(a, "_attr") else a
    edges[[col]] <- plain_column(eattrs[[a]])
  }
  weight <- if (is.numeric(eattrs$weight)) "weight"
  if (!is.null(weight)) check_weight_type(weight, weight_type)
  graph <- igraph::graph_attr(x)
  # igraph's logical `type` attribute declares the modes of a two-mode
  # network.
  if (is.logical(nodes$type) && !anyNA(nodes$type) && nrow(nodes)) {
    return(net_twomode(edges, nodes, mode = "type", directed = igraph::is_directed(x),
                       id = "name", weight = weight, weight_type = weight_type, quiet = quiet,
                       graph = graph))
  }
  net_frame(edges, nodes, directed = igraph::is_directed(x), id = "name", weight = weight,
            weight_type = weight_type, quiet = quiet, graph = graph)
}

# See man/net_from.Rd.
net_from.network <- function(x, weight_type = "strength", quiet = FALSE, ...) {
  require_pkg("network")
  n <- network::network.size(x)
  vnames <- network::network.vertex.names(x)
  ids <- if (length(vnames) == n && !anyDuplicated(vnames)) vnames else seq_len(n)
  nodes <- data.frame(name = ids, stringsAsFactors = FALSE)
  for (a in setdiff(network::list.vertex.attributes(x), c("vertex.names", "na"))) {
    nodes[[a]] <- plain_column(network::get.vertex.attribute(x, a, unlist = FALSE))
  }
  b <- network::get.network.attribute(x, "bipartite")
  bipartite <- is.numeric(b) && length(b) == 1L && b > 0
  if (bipartite) nodes$mode <- ifelse(seq_len(n) <= b, "first", "second")

  # Read edges in storage order so edge attributes line up.
  eids <- network::valid.eids(x)
  tails <- vapply(x$mel[eids], function(e) as.numeric(e$outl[[1]]), 0)
  heads <- vapply(x$mel[eids], function(e) as.numeric(e$inl[[1]]), 0)
  edges <- data.frame(source = ids[tails], target = ids[heads], stringsAsFactors = FALSE)
  for (a in setdiff(network::list.edge.attributes(x), "na")) {
    vals <- network::get.edge.attribute(x, a, unlist = FALSE)[eids]
    col <- if (a %in% c("source", "target")) paste0(a, "_attr") else a
    edges[[col]] <- plain_column(vals)
  }
  graph <- list()
  for (a in setdiff(network::list.network.attributes(x),
                    c("bipartite", "directed", "hyper", "loops", "multiple", "n", "mnext"))) {
    v <- network::get.network.attribute(x, a)
    if (is.atomic(v) && length(v) == 1L) graph[[a]] <- v
  }
  weight <- if (is.numeric(edges$weight)) "weight"
  if (!is.null(weight)) check_weight_type(weight, weight_type)
  if (bipartite) {
    return(net_twomode(edges, nodes, mode = "mode", directed = network::is.directed(x),
                       id = "name", weight = weight, weight_type = weight_type, quiet = quiet,
                       graph = graph))
  }
  net_frame(edges, nodes, directed = network::is.directed(x), id = "name", weight = weight,
            weight_type = weight_type, quiet = quiet, graph = graph)
}

# See man/net_from.Rd.
net_from.matrix <- function(x, directed = NULL, weighted = NULL, weight_type = "strength",
                               quiet = FALSE, ...) {
  check_flag_or_null(directed, "directed")
  check_flag_or_null(weighted, "weighted")
  if (nrow(x) != ncol(x)) {
    stop("An adjacency matrix must be square.", call. = FALSE)
  }
  if (!is.numeric(x) && !is.logical(x)) {
    stop("An adjacency matrix must be numeric or logical.", call. = FALSE)
  }
  x[is.na(x)] <- 0
  idx <- which(x != 0, arr.ind = TRUE)
  matrix_to_netframe(idx[, 1], idx[, 2], as.numeric(x[idx]), nrow(x), rownames(x),
                     symmetric = isSymmetric(unname(x)), directed, weighted, weight_type,
                     quiet)
}

# See man/net_from.Rd.
net_from.Matrix <- function(x, directed = NULL, weighted = NULL, weight_type = "strength",
                               quiet = FALSE, ...) {
  require_pkg("Matrix")
  check_flag_or_null(directed, "directed")
  check_flag_or_null(weighted, "weighted")
  if (nrow(x) != ncol(x)) {
    stop("An adjacency matrix must be square.", call. = FALSE)
  }
  t <- methods::as(x, "TsparseMatrix")
  i <- t@i + 1L
  j <- t@j + 1L
  v <- if (methods::.hasSlot(t, "x")) as.numeric(t@x) else rep(1, length(i))
  # Symmetric storage keeps one triangle; expand it.
  if (methods::is(t, "symmetricMatrix")) {
    off <- i != j
    i2 <- c(i, j[off])
    j <- c(j, i[off])
    i <- i2
    v <- c(v, v[off])
  }
  # As for dense matrices, missing values are read as 0 (no tie).
  keep <- !is.na(v) & v != 0
  matrix_to_netframe(i[keep], j[keep], v[keep], nrow(x), rownames(x),
                     symmetric = Matrix::isSymmetric(x), directed, weighted, weight_type, quiet)
}

matrix_to_netframe <- function(i, j, v, n, ids, symmetric, directed, weighted, weight_type,
                               quiet) {
  if (is.null(ids)) ids <- seq_len(n)
  if (anyDuplicated(ids)) {
    stop("Row names of the adjacency matrix must be unique.", call. = FALSE)
  }
  directed <- if (is.null(directed)) !symmetric else isTRUE(directed)
  if (!directed) {
    keep <- i <= j
    i <- i[keep]
    j <- j[keep]
    v <- v[keep]
  }
  o <- order(i, j)
  i <- i[o]
  j <- j[o]
  v <- v[o]
  weighted <- if (is.null(weighted)) any(v != 1) else isTRUE(weighted)
  if (weighted) check_weight_type("weight", weight_type)
  edges <- data.frame(source = ids[i], target = ids[j], stringsAsFactors = FALSE)
  if (weighted) edges$weight <- v
  nodes <- data.frame(id = ids, stringsAsFactors = FALSE)
  net_frame(edges, nodes, directed = directed, weight = if (weighted) "weight",
            weight_type = if (weighted) weight_type, quiet = quiet)
}

# See man/net_from.Rd.
net_to_igraph <- function(x) {
  check_netframe(x)
  require_pkg("igraph")
  g <- igraph::make_empty_graph(n = nrow(x$nodes), directed = x$directed)
  if (nrow(x$edges)) {
    g <- igraph::add_edges(g, rbind(x$index$src, x$index$tgt) + 1L)
  }
  igraph::V(g)$name <- as.character(x$nodes[[x$roles$id]])
  for (a in setdiff(names(x$nodes), x$roles$id)) {
    col <- if (a == "name") "name_attr" else a
    g <- igraph::set_vertex_attr(g, col, value = x$nodes[[a]])
  }
  for (a in setdiff(names(x$edges), c(x$roles$source, x$roles$target))) {
    g <- igraph::set_edge_attr(g, a, value = x$edges[[a]])
  }
  # igraph declares two-mode networks with a logical `type` attribute.
  if (net_is_twomode(x) && !identical(x$roles$mode, "type")) {
    g <- igraph::set_vertex_attr(g, "type", value = twomode_types(x) == 1L)
  }
  if (!is.null(x$roles$weight) && x$roles$weight != "weight" &&
      !"weight" %in% names(x$edges)) {
    g <- igraph::set_edge_attr(g, "weight", value = x$edges[[x$roles$weight]])
  }
  for (a in names(x$graph)) g <- igraph::set_graph_attr(g, a, x$graph[[a]])
  g
}

# See man/net_from.Rd.
net_to_tbl_graph <- function(x) {
  require_pkg("tidygraph")
  tidygraph::as_tbl_graph(net_to_igraph(x))
}

# See man/net_from.Rd.
net_to_network <- function(x) {
  check_netframe(x)
  require_pkg("network")
  n <- nrow(x$nodes)
  src <- x$index$src + 1L
  tgt <- x$index$tgt + 1L
  pairs <- if (x$directed) paste(src, tgt) else paste(pmin(src, tgt), pmax(src, tgt))
  bipartite <- if (net_is_twomode(x)) bipartite_count(twomode_types(x))
  net <- network::network.initialize(n, directed = x$directed, loops = any(src == tgt),
                                     multiple = anyDuplicated(pairs) > 0,
                                     bipartite = if (is.null(bipartite)) FALSE else bipartite)
  network::set.vertex.attribute(net, "vertex.names", as.character(x$nodes[[x$roles$id]]))
  for (a in setdiff(names(x$nodes), x$roles$id)) {
    network::set.vertex.attribute(net, a, x$nodes[[a]])
  }
  if (length(src)) {
    attrs <- setdiff(names(x$edges), c(x$roles$source, x$roles$target))
    names_list <- rep(list(attrs), length(src))
    vals_list <- lapply(seq_along(src), function(k) lapply(attrs, function(a) x$edges[[a]][[k]]))
    network::add.edges(net, tail = src, head = tgt,
                       names.eval = if (length(attrs)) names_list,
                       vals.eval = if (length(attrs)) vals_list)
  }
  for (a in names(x$graph)) network::set.network.attribute(net, a, x$graph[[a]])
  net
}

# network's bipartite flag requires the first `b` vertices to be one mode;
# only then is it set (the edges of a two-mode netframe already cross
# modes).
bipartite_count <- function(types) {
  b <- sum(types == 0L)
  if (b > 0 && all(types[seq_len(b)] == 0L)) b else NULL
}

# See man/net_from.Rd.
as.matrix.netframe <- function(x, attr = NULL, sparse = FALSE, ...) {
  check_netframe(x)
  check_flag(sparse, "sparse")
  n <- nrow(x$nodes)
  i <- x$index$src + 1L
  j <- x$index$tgt + 1L
  if (!is.null(attr) && (!is.atomic(attr) || length(attr) != 1L || is.na(attr))) {
    stop("`attr` must name one edge column, or be NULL.", call. = FALSE)
  }
  if (!is.null(attr) && is.null(x$edges[[attr]])) {
    stop(sprintf("Edge column `%s` not found.", attr), call. = FALSE)
  }
  v <- if (is.null(attr)) rep(1, length(i)) else as.numeric(x$edges[[attr]])
  if (!x$directed) {
    off <- i != j
    i2 <- c(i, j[off])
    j <- c(j, i[off])
    i <- i2
    v <- c(v, v[off])
  }
  ids <- as.character(x$nodes[[x$roles$id]])
  if (sparse) {
    require_pkg("Matrix")
    return(Matrix::sparseMatrix(i = i, j = j, x = v, dims = c(n, n),
                                dimnames = list(ids, ids)))
  }
  m <- matrix(0, n, n, dimnames = list(ids, ids))
  if (length(i)) {
    cell <- (j - 1) * n + i
    sums <- rowsum(v, cell)
    m[as.numeric(rownames(sums))] <- sums[, 1]
  }
  m
}

plain_column <- function(v) {
  if (is.list(v)) {
    flat <- vapply(v, function(e) length(e) == 1L && is.atomic(e), TRUE)
    if (all(flat)) return(unlist(v, use.names = FALSE))
  }
  v
}

# Methods for other packages' generics, registered when those packages load.
as.igraph.netframe <- function(x, ...) net_to_igraph(x)
as.network.netframe <- function(x, ...) net_to_network(x)
as_tbl_graph.netframe <- function(x, ...) net_to_tbl_graph(x)
