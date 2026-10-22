#' Degree and strength
#'
#' Counts edge endpoints at each node. With `weighted = TRUE`, sums edge
#' weights instead (often called strength).
#'
#' Parallel edges each count. In undirected networks a self-loop adds 2 to
#' its node; in directed networks it adds 1 to out-degree and 1 to
#' in-degree. These match igraph and NetworkX. Set `loops = FALSE` to ignore
#' self-loops, which matches the sna default (`diag = FALSE`).
#'
#' @param x A [netframe()].
#' @param mode For directed networks, `"all"`, `"out"`, or `"in"`. Ignored
#'   for undirected networks.
#' @param loops Whether self-loops count.
#' @param weighted Whether to sum edge weights.
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
degree <- function(x, mode = c("all", "out", "in"), loops = TRUE, weighted = FALSE) {
  check_netframe(x)
  mode <- match.arg(mode)
  if (weighted) require_weights(x)
  code <- switch(mode, out = 1L, `in` = 2L, all = 3L)
  .Call(C_degree, engine_graph(x), code, isTRUE(loops), isTRUE(weighted))
}

#' Betweenness centrality
#'
#' Shortest-path betweenness (Freeman 1977), computed with the algorithm of
#' Brandes (2001) in parallel.
#'
#' Endpoints are excluded and, in undirected networks, each pair of nodes is
#' counted once. Parallel edges count as distinct shortest paths. With
#' `weighted = TRUE`, weights are treated as distances and must be positive.
#' Normalized values divide by the number of pairs not involving the node:
#' \eqn{(n-1)(n-2)} for directed and \eqn{(n-1)(n-2)/2} for undirected
#' networks. These conventions match igraph and NetworkX.
#'
#' @inheritParams degree
#' @param weighted Whether to use edge weights as distances.
#' @param normalized Whether to normalize by the number of pairs.
#' @param threads Number of threads; `0` uses all available cores.
#'
#' @references
#' Freeman, L. C. (1977). A set of measures of centrality based on
#' betweenness. *Sociometry*, 40(1), 35-41. \doi{10.2307/3033543}
#'
#' Brandes, U. (2001). A faster algorithm for betweenness centrality.
#' *Journal of Mathematical Sociology*, 25(2), 163-177.
#' \doi{10.1080/0022250X.2001.9990249}
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
betweenness <- function(x, weighted = FALSE, normalized = FALSE, threads = 0L) {
  check_netframe(x)
  if (weighted) require_weights(x)
  .Call(C_betweenness, engine_graph(x), isTRUE(weighted), isTRUE(normalized),
        as.integer(threads))
}

#' Triad census
#'
#' Counts every set of three nodes by its type among the 16 directed triad
#' types of Holland and Leinhardt (1976), using the algorithm of Batagelj
#' and Mrvar (2001).
#'
#' Self-loops are ignored and parallel edges are counted once. In an
#' undirected network every edge is treated as mutual, as in igraph.
#'
#' @inheritParams betweenness
#'
#' @references
#' Holland, P. W., & Leinhardt, S. (1976). Local structure in social
#' networks. *Sociological Methodology*, 7, 1-45. \doi{10.2307/270703}
#'
#' Batagelj, V., & Mrvar, A. (2001). A subquadratic triad census algorithm
#' for large sparse networks with small maximum degree. *Social Networks*,
#' 23(3), 237-243. \doi{10.1016/S0378-8733(01)00035-1}
#'
#' @return A named numeric vector of 16 counts in the standard order.
#' @export
triad_census <- function(x, threads = 0L) {
  check_netframe(x)
  out <- .Call(C_triad_census, engine_graph(x), as.integer(threads))
  names(out) <- c("003", "012", "102", "021D", "021U", "021C", "111D", "111U",
                  "030T", "030C", "201", "120D", "120U", "120C", "210", "300")
  out
}

require_weights <- function(x) {
  if (is.null(x$roles$weight)) {
    stop("This network has no weights; supply `weight =` to netframe().",
         call. = FALSE)
  }
}
