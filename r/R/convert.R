#' Convert networks to and from other packages
#'
#' `as_netframe()` turns network objects from other packages, adjacency
#' matrices, and edge lists into a [netframe()]. `to_igraph()`,
#' `to_network()`, and `to_tbl_graph()` convert back, and `as.matrix()`
#' gives an adjacency matrix. Node, edge, and network-level attributes are
#' carried across in both directions, so existing code and data keep
#' working.
#'
#' The usual generics of other packages also work on a netframe:
#' `igraph::as.igraph()`, `network::as.network()`, and
#' `tidygraph::as_tbl_graph()`.
#'
#' Notes on each source:
#' * **igraph** (including **tidygraph** `tbl_graph` objects): vertex names
#'   become node identifiers (vertex positions are used if the graph has no
#'   names). An edge attribute called `weight` becomes the weight column.
#' * **network** (statnet): vertex names become node identifiers. For
#'   bipartite networks a `mode` node column (1 or 2) records the two sets.
#' * **matrix** or sparse **Matrix**: row names become identifiers. Values
#'   other than 0 and 1 become edge weights. Direction is inferred from
#'   symmetry unless `directed` is given; an undirected matrix contributes
#'   only its upper triangle.
#' * **data.frame**: treated as an edge list, as in [netframe()].
#'
#' @param x An object to convert.
#' @param ... Passed on to methods. For data frames, arguments of
#'   [netframe()].
#' @param directed For matrices, whether edges are directed. `NULL` infers
#'   it from symmetry.
#' @param weighted For matrices, whether values are edge weights. `NULL`
#'   treats any value other than 0 and 1 as a weight.
#'
#' @return `as_netframe()` returns a [netframe()]. The `to_*()` functions
#'   return objects of the corresponding package.
#'
#' @examples
#' m <- matrix(c(0, 2, 0,
#'               2, 0, 1,
#'               0, 1, 0), 3, byrow = TRUE, dimnames = list(letters[1:3], letters[1:3]))
#' g <- as_netframe(m)
#' g$edges
#' as.matrix(g, attr = "weight")
#' @export
as_netframe <- function(x, ...) {
  UseMethod("as_netframe")
}

#' @rdname as_netframe
#' @export
as_netframe.default <- function(x, ...) {
  if (isS4(x) && methods::is(x, "Matrix")) {
    return(as_netframe.Matrix(x, ...))
  }
  stop(sprintf(paste0(
    "Cannot convert an object of class `%s` to a netframe. Supported: igraph, ",
    "tbl_graph, network, matrix, Matrix, and data.frame (an edge list)."),
    class(x)[[1]]), call. = FALSE)
}

#' @rdname as_netframe
#' @export
as_netframe.netframe <- function(x, ...) x

#' @rdname as_netframe
#' @export
as_netframe.data.frame <- function(x, ...) netframe(x, ...)

#' @rdname as_netframe
#' @export
as_netframe.igraph <- function(x, ...) {
  require_pkg("igraph")
  n <- igraph::vcount(x)
  vattrs <- igraph::vertex_attr(x)
  named <- "name" %in% names(vattrs)
  ids <- if (named) vattrs$name else seq_len(n)
  vattrs$name <- NULL
  nodes <- data.frame(name = ids, stringsAsFactors = FALSE)
  for (a in names(vattrs)) nodes[[a]] <- plain_column(vattrs[[a]])

  ends <- igraph::ends(x, igraph::E(x), names = FALSE)
  edges <- data.frame(from = ids[ends[, 1]], to = ids[ends[, 2]], stringsAsFactors = FALSE)
  eattrs <- igraph::edge_attr(x)
  for (a in names(eattrs)) {
    col <- if (a %in% c("from", "to")) paste0(a, "_attr") else a
    edges[[col]] <- plain_column(eattrs[[a]])
  }
  weight <- if (is.numeric(eattrs$weight)) "weight"
  netframe(edges, nodes, directed = igraph::is_directed(x), from = "from", to = "to",
           id = "name", weight = weight, graph = igraph::graph_attr(x))
}

#' @rdname as_netframe
#' @export
as_netframe.network <- function(x, ...) {
  require_pkg("network")
  n <- network::network.size(x)
  vnames <- network::network.vertex.names(x)
  ids <- if (length(vnames) == n && !anyDuplicated(vnames)) vnames else seq_len(n)
  nodes <- data.frame(name = ids, stringsAsFactors = FALSE)
  for (a in setdiff(network::list.vertex.attributes(x), c("vertex.names", "na"))) {
    nodes[[a]] <- plain_column(network::get.vertex.attribute(x, a, unlist = FALSE))
  }
  b <- network::get.network.attribute(x, "bipartite")
  if (is.numeric(b) && length(b) == 1L && b > 0) {
    nodes$mode <- ifelse(seq_len(n) <= b, 1L, 2L)
  }

  # Read edges in storage order so edge attributes line up.
  eids <- network::valid.eids(x)
  tails <- vapply(x$mel[eids], function(e) as.numeric(e$outl[[1]]), 0)
  heads <- vapply(x$mel[eids], function(e) as.numeric(e$inl[[1]]), 0)
  edges <- data.frame(from = ids[tails], to = ids[heads], stringsAsFactors = FALSE)
  for (a in setdiff(network::list.edge.attributes(x), "na")) {
    vals <- network::get.edge.attribute(x, a, unlist = FALSE)[eids]
    col <- if (a %in% c("from", "to")) paste0(a, "_attr") else a
    edges[[col]] <- plain_column(vals)
  }
  graph <- list()
  for (a in setdiff(network::list.network.attributes(x),
                    c("bipartite", "directed", "hyper", "loops", "multiple", "n", "mnext"))) {
    v <- network::get.network.attribute(x, a)
    if (is.atomic(v) && length(v) == 1L) graph[[a]] <- v
  }
  weight <- if (is.numeric(edges$weight)) "weight"
  netframe(edges, nodes, directed = network::is.directed(x), from = "from", to = "to",
           id = "name", weight = weight, graph = graph)
}

#' @rdname as_netframe
#' @export
as_netframe.matrix <- function(x, directed = NULL, weighted = NULL, ...) {
  if (nrow(x) != ncol(x)) {
    stop("An adjacency matrix must be square.", call. = FALSE)
  }
  if (!is.numeric(x) && !is.logical(x)) {
    stop("An adjacency matrix must be numeric or logical.", call. = FALSE)
  }
  x[is.na(x)] <- 0
  idx <- which(x != 0, arr.ind = TRUE)
  matrix_to_netframe(idx[, 1], idx[, 2], as.numeric(x[idx]), nrow(x), rownames(x),
                     symmetric = isSymmetric(unname(x)), directed, weighted)
}

#' @rdname as_netframe
#' @export
as_netframe.Matrix <- function(x, directed = NULL, weighted = NULL, ...) {
  require_pkg("Matrix")
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
  keep <- v != 0
  matrix_to_netframe(i[keep], j[keep], v[keep], nrow(x), rownames(x),
                     symmetric = Matrix::isSymmetric(x), directed, weighted)
}

matrix_to_netframe <- function(i, j, v, n, ids, symmetric, directed, weighted) {
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
  edges <- data.frame(from = ids[i], to = ids[j], stringsAsFactors = FALSE)
  if (weighted) edges$weight <- v
  nodes <- data.frame(name = ids, stringsAsFactors = FALSE)
  netframe(edges, nodes, directed = directed, from = "from", to = "to", id = "name",
           weight = if (weighted) "weight")
}

#' @rdname as_netframe
#' @export
to_igraph <- function(x) {
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
  for (a in setdiff(names(x$edges), c(x$roles$from, x$roles$to))) {
    g <- igraph::set_edge_attr(g, a, value = x$edges[[a]])
  }
  if (!is.null(x$roles$weight) && x$roles$weight != "weight" &&
      !"weight" %in% names(x$edges)) {
    g <- igraph::set_edge_attr(g, "weight", value = x$edges[[x$roles$weight]])
  }
  for (a in names(x$graph)) g <- igraph::set_graph_attr(g, a, x$graph[[a]])
  g
}

#' @rdname as_netframe
#' @export
to_tbl_graph <- function(x) {
  require_pkg("tidygraph")
  tidygraph::as_tbl_graph(to_igraph(x))
}

#' @rdname as_netframe
#' @export
to_network <- function(x) {
  check_netframe(x)
  require_pkg("network")
  n <- nrow(x$nodes)
  src <- x$index$src + 1L
  tgt <- x$index$tgt + 1L
  pairs <- if (x$directed) paste(src, tgt) else paste(pmin(src, tgt), pmax(src, tgt))
  bipartite <- if ("mode" %in% names(x$nodes)) bipartite_count(x$nodes$mode, src, tgt)
  net <- network::network.initialize(n, directed = x$directed, loops = any(src == tgt),
                                     multiple = anyDuplicated(pairs) > 0,
                                     bipartite = if (is.null(bipartite)) FALSE else bipartite)
  network::set.vertex.attribute(net, "vertex.names", as.character(x$nodes[[x$roles$id]]))
  for (a in setdiff(names(x$nodes), c(x$roles$id, if (!is.null(bipartite)) "mode"))) {
    network::set.vertex.attribute(net, a, x$nodes[[a]])
  }
  if (length(src)) {
    attrs <- setdiff(names(x$edges), c(x$roles$from, x$roles$to))
    names_list <- rep(list(attrs), length(src))
    vals_list <- lapply(seq_along(src), function(k) lapply(attrs, function(a) x$edges[[a]][[k]]))
    network::add.edges(net, tail = src, head = tgt,
                       names.eval = if (length(attrs)) names_list,
                       vals.eval = if (length(attrs)) vals_list)
  }
  for (a in names(x$graph)) network::set.network.attribute(net, a, x$graph[[a]])
  net
}

# network's bipartite flag requires the first `b` vertices to be one mode
# and every edge to cross modes; only then is it set.
bipartite_count <- function(mode, src, tgt) {
  b <- sum(mode == 1, na.rm = TRUE)
  ok <- !anyNA(mode) && all(mode[seq_len(b)] == 1) && all(mode[src] != mode[tgt])
  if (ok && b > 0) b else NULL
}

#' @rdname as_netframe
#' @param attr For `as.matrix()`, an edge column whose values fill the
#'   matrix; `NULL` counts edges.
#' @param sparse For `as.matrix()`, whether to return a sparse `Matrix`.
#' @details `as.matrix()` sums parallel edges. In undirected networks the
#'   matrix is symmetric and a self-loop contributes once to the diagonal.
#' @export
as.matrix.netframe <- function(x, attr = NULL, sparse = FALSE, ...) {
  n <- nrow(x$nodes)
  i <- x$index$src + 1L
  j <- x$index$tgt + 1L
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

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("The %s package is required; install it with install.packages(\"%s\").",
                 pkg, pkg), call. = FALSE)
  }
}

# Methods for other packages' generics, registered when those packages load.
as.igraph.netframe <- function(x, ...) to_igraph(x)
as.network.netframe <- function(x, ...) to_network(x)
as_tbl_graph.netframe <- function(x, ...) to_tbl_graph(x)
