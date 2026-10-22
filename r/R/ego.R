#' Ego networks
#'
#' `ego_network()` extracts the ego network of one node: the node (ego),
#' every node within `order` steps of it (its alters), and all ties among
#' them. `ego_networks()` does the same for several nodes, or for every
#' node, at once.
#'
#' The result is the subgraph induced by ego and its alters, as igraph's
#' `make_ego_graph()` and NetworkX's `ego_graph()` return: every edge whose
#' two endpoints are both in the ego network is kept, with all its columns,
#' and so is every column of the node table. Nodes and edges keep their
#' order in `x`. The weight column and its meaning carry over, and so does
#' the node type of a two-mode network (see [twomode()]).
#'
#' In directed networks `mode` chooses the alters: `"out"` takes the nodes
#' ego reaches along outgoing ties (whom ego names), `"in"` the nodes that
#' reach ego (who names ego), and `"all"` ignores direction. Ties among the
#' chosen nodes are kept whatever their direction. Undirected networks
#' ignore `mode`. Steps are counted without weights.
#'
#' @param x A [netframe()].
#' @param node The identifier of the ego, as in the identifier column of
#'   `x$nodes`.
#' @param nodes Identifiers of the egos; `NULL` (the default) uses every
#'   node.
#' @param order Alters are the nodes within this many steps of ego. `1`
#'   gives the first-order neighborhood.
#' @param mode For directed networks, `"all"`, `"out"`, or `"in"`; see
#'   Details.
#' @param include_ego Whether ego itself is part of the result (igraph's
#'   `mindist = 0`) or only its alters (`mindist = 1`).
#' @param threads Number of threads; `0` uses all available cores.
#'
#' @references
#' Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and
#' Applications*. Cambridge University Press.
#' \doi{10.1017/CBO9780511815478}
#'
#' Perry, B. L., Pescosolido, B. A., & Borgatti, S. P. (2018). *Egocentric
#' Network Analysis: Foundations, Methods, and Models*. Cambridge
#' University Press. \doi{10.1017/9781316443255}
#'
#' @seealso [ego_summary()] for measures of every ego network at once.
#'
#' @return `ego_network()` returns a [netframe()]. `ego_networks()` returns
#'   a list of netframes named by the egos' identifiers.
#'
#' @examples
#' g <- netframe(data.frame(from = c("a", "a", "b", "c"), to = c("b", "c", "c", "d")),
#'               directed = FALSE)
#' ego_network(g, "c")
#' ego_network(g, "c", include_ego = FALSE)$nodes
#' lengths(lapply(ego_networks(g), function(e) e$nodes$name))
#' @export
ego_network <- function(x, node, order = 1L, mode = c("all", "out", "in"),
                        include_ego = TRUE) {
  check_netframe(x)
  if (length(node) != 1L) {
    stop("`node` must be a single node identifier; use ego_networks() for several.",
         call. = FALSE)
  }
  ego_networks(x, node, order = order, mode = match.arg(mode),
               include_ego = include_ego, threads = 1L)[[1L]]
}

#' @rdname ego_network
#' @export
ego_networks <- function(x, nodes = NULL, order = 1L, mode = c("all", "out", "in"),
                         include_ego = TRUE, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  order <- check_count(order, "order")
  egos <- if (is.null(nodes)) seq_len(nrow(x$nodes)) else node_positions(x, nodes)
  members <- .Call(C_ego_members, engine_graph(x), as.integer(egos - 1L), order,
                   mode_code(mode), if (isTRUE(include_ego)) 0L else 1L, resolve_threads(threads))
  out <- lapply(members, function(keep) induced_netframe(x, keep))
  names(out) <- as.character(x$nodes[[x$roles$id]][egos])
  out
}

#' Egocentric measures
#'
#' Measures of every node's ego network, computed together in parallel and
#' returned aligned with the rows of `x$nodes`. `ego_summary()` gives them
#' all in a data frame; `ego_size()`, `ego_density()`, and
#' `ego_betweenness()` return one measure as a vector.
#'
#' For a node with \eqn{k} alters (the nodes within `order` steps; see
#' [ego_network()]):
#' \describe{
#'   \item{`size`}{\eqn{k}, the number of alters.}
#'   \item{`ties`}{The number of ties among the alters, leaving out ego's
#'     own ties.}
#'   \item{`pairs`}{The number of pairs of alters: \eqn{k(k-1)/2}, or
#'     \eqn{k(k-1)} ordered pairs in a directed network.}
#'   \item{`density`}{`ties / pairs`, the density of the ego network
#'     without ego; `NaN` for fewer than two alters.}
#'   \item{`components`}{The number of weak components among the alters
#'     once ego is removed; isolated alters count as components.}
#'   \item{`ego_betweenness`}{The betweenness of ego within its own ego
#'     network (Everett and Borgatti 2005; Marsden 2002): for every pair of
#'     alters, the share of the shortest paths between them, inside the ego
#'     network, that pass through ego. Each unordered pair counts once in an
#'     undirected network and each ordered pair in a directed one. Not
#'     normalized; the largest possible value is `pairs`.}
#' }
#'
#' `ego_size()` is igraph's `ego_size()`: with `include_ego = TRUE` (the
#' default, igraph's `mindist = 0`) it counts ego as well, so it is
#' `size + 1`.
#'
#' Ties are binary: weights are ignored, parallel edges count once, and
#' self-loops are ignored, as in UCINET's ego network measures. In directed
#' networks `mode` chooses the alters (see [ego_network()]), every tie
#' among them counts with its direction, and ego betweenness follows
#' directed paths. Burt's structural-hole measures, which are also
#' egocentric, are [constraint()], [effective_size()], and [efficiency()].
#'
#' Results match igraph (`ego_size()`, and `ecount()`, `edge_density()`,
#' `count_components()`, and `betweenness()` of
#' `induced_subgraph(simplify(g), ego(g, ...))`) and NetworkX (`ego_graph()`
#' on the simple graph).
#'
#' @inheritParams ego_network
#' @param betweenness Whether to compute ego betweenness, the costliest of
#'   the measures.
#'
#' @references
#' Everett, M., & Borgatti, S. P. (2005). Ego network betweenness. *Social
#' Networks*, 27(1), 31-38. \doi{10.1016/j.socnet.2004.11.007}
#'
#' Marsden, P. V. (2002). Egocentric and sociocentric measures of network
#' centrality. *Social Networks*, 24(4), 407-422.
#' \doi{10.1016/S0378-8733(02)00016-3}
#'
#' Perry, B. L., Pescosolido, B. A., & Borgatti, S. P. (2018). *Egocentric
#' Network Analysis: Foundations, Methods, and Models*. Cambridge
#' University Press. \doi{10.1017/9781316443255}
#'
#' @return `ego_summary()` returns a data frame with one row per node, in
#'   the order of `x$nodes`, and columns `node` (the identifier), `size`,
#'   `ties`, `pairs`, `density`, `components`, and (if requested)
#'   `ego_betweenness`. The other functions return a numeric vector aligned
#'   with the rows of `x$nodes`.
#'
#' @examples
#' g <- netframe(data.frame(from = c("a", "a", "a", "b", "d"),
#'                          to   = c("b", "c", "d", "c", "e")),
#'               directed = FALSE)
#' ego_summary(g)
#' ego_size(g, order = 2)
#' @export
ego_summary <- function(x, order = 1L, mode = c("all", "out", "in"), betweenness = TRUE,
                        threads = 0L) {
  check_netframe(x)
  res <- ego_measures(x, order, match.arg(mode), isTRUE(betweenness), threads)
  out <- data.frame(node = x$nodes[[x$roles$id]], size = res[[1]], ties = res[[2]],
                    pairs = res[[3]], density = res[[4]], components = res[[5]],
                    stringsAsFactors = FALSE)
  if (isTRUE(betweenness)) out$ego_betweenness <- res[[6]]
  out
}

#' @rdname ego_summary
#' @export
ego_size <- function(x, order = 1L, mode = c("all", "out", "in"), include_ego = TRUE,
                     threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  order <- check_count(order, "order")
  .Call(C_ego_sizes, engine_graph(x), NULL, order, mode_code(mode),
        if (isTRUE(include_ego)) 0L else 1L, resolve_threads(threads))
}

#' @rdname ego_summary
#' @export
ego_density <- function(x, order = 1L, mode = c("all", "out", "in"), threads = 0L) {
  check_netframe(x)
  ego_measures(x, order, match.arg(mode), FALSE, threads)[[4]]
}

#' @rdname ego_summary
#' @export
ego_betweenness <- function(x, order = 1L, mode = c("all", "out", "in"), threads = 0L) {
  check_netframe(x)
  ego_measures(x, order, match.arg(mode), TRUE, threads)[[6]]
}

ego_measures <- function(x, order, mode, betweenness, threads) {
  order <- check_count(order, "order")
  if (order < 1L) stop("`order` must be at least 1.", call. = FALSE)
  .Call(C_ego_summary, engine_graph(x), order, mode_code(mode), betweenness,
        resolve_threads(threads))
}

# Positions (1-based) of node identifiers in x$nodes.
node_positions <- function(x, nodes) {
  ids <- x$nodes[[x$roles$id]]
  if (is.factor(ids)) ids <- as.character(ids)
  if (is.factor(nodes)) nodes <- as.character(nodes)
  pos <- match(nodes, ids)
  if (anyNA(pos)) {
    shown <- utils::head(unique(nodes[is.na(pos)]), 5L)
    stop(sprintf("%d node(s) not found in `x$nodes`, for example: %s.",
                 length(unique(nodes[is.na(pos)])), paste(format(shown), collapse = ", ")),
         call. = FALSE)
  }
  pos
}

# The subgraph induced by the nodes at 1-based positions `keep` (in
# ascending order): those nodes and every edge between two of them, with
# all columns and roles (including a two-mode network's node type).
induced_netframe <- function(x, keep) {
  inside <- logical(nrow(x$nodes))
  inside[keep] <- TRUE
  e <- inside[x$index$src + 1L] & inside[x$index$tgt + 1L]
  nodes <- x$nodes[keep, , drop = FALSE]
  edges <- x$edges[e, , drop = FALSE]
  rownames(nodes) <- NULL
  rownames(edges) <- NULL
  weight_type <- if (is.null(x$roles$weight_type)) "distance" else x$roles$weight_type
  out <- netframe(edges, nodes, directed = x$directed, from = x$roles$from,
                  to = x$roles$to, id = x$roles$id, weight = x$roles$weight,
                  weight_type = weight_type, graph = x$graph)
  extra <- setdiff(names(x$roles), names(out$roles))
  if (length(extra)) out$roles[extra] <- x$roles[extra]
  out
}
