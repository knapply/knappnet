#' Reach centrality
#'
#' `reach()` counts the other nodes within `k` steps of each node (UCINET's
#' reach centrality, and the reach of Borgatti's (2006) key player
#' problem). `distance_weighted_reach()` is UCINET's average reciprocal
#' distance (ARD): the mean of 1/d over the other nodes, with 1/d taken as
#' 0 for nodes that cannot be reached, which is normalized harmonic
#' centrality (see [harmonic_centrality()]).
#'
#' @inheritParams betweenness
#' @param k Largest distance counted: a number of steps, or a path length
#'   when weights are used. `Inf` counts every reachable node (see
#'   [reachability()]).
#' @param mode For directed networks: `"out"` counts the nodes a node
#'   reaches, `"in"` the nodes that reach it, and `"all"` ignores
#'   direction. Ignored for undirected networks.
#' @param normalized Whether to divide by \eqn{n - 1}, the largest possible
#'   count (UCINET reports reach as this proportion).
#'
#' @details
#' Without weights every edge has length 1, so `k` counts steps. With
#' weights (distances, or strengths converted to lengths 1/weight, as in
#' [betweenness()]) `k` is a path length and a node counts when its
#' shortest-path length is at most `k`. Self-loops never lie on a shortest
#' path and parallel edges do not change path lengths. The node itself is
#' never counted. Searches from each node run in parallel.
#'
#' `reach(x, k = 1)` is the number of distinct neighbours and
#' `reach(x, k = Inf)` is [reachability()]; UCINET's 2-step reach is the
#' default `k = 2`.
#'
#' @references
#' Borgatti, S. P. (2006). Identifying sets of key players in a social
#' network. *Computational and Mathematical Organization Theory*, 12(1),
#' 21-34. \doi{10.1007/s10588-006-7084-x}
#'
#' Everton, S. F. (2012). *Disrupting Dark Networks*. Cambridge University
#' Press. \doi{10.1017/CBO9781139136877}
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @seealso [harmonic_centrality()], [closeness()], [reachability()].
#' @export
reach <- function(x, k = 2, mode = c("out", "in", "all"), weighted = NULL, weights_as = NULL,
                  normalized = FALSE, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 0) {
    stop("`k` must be a single non-negative number (Inf counts every reachable node).",
         call. = FALSE)
  }
  normalized <- check_flag(normalized, "normalized")
  threads <- resolve_threads(threads)
  w <- resolve_weighting(x, weighted, weights_as)
  ptr <- engine_graph(x)
  engine_call("reach", .Call(C_reach, ptr, as.double(k), mode_code(mode), w$weighted, w$code,
                             normalized, threads))
}

#' @rdname reach
#' @export
distance_weighted_reach <- function(x, mode = c("out", "in", "all"), weighted = NULL,
                                    weights_as = NULL, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  threads <- resolve_threads(threads)
  w <- resolve_weighting(x, weighted, weights_as)
  ptr <- engine_graph(x)
  engine_call("distance_weighted_reach",
              .Call(C_closeness, ptr, TRUE, mode_code(mode), w$weighted, w$code, TRUE, 0L,
                    threads))
}

#' Proximity prestige
#'
#' The prestige of a node in a directed network from the nodes that can
#' reach it, its influence domain (Lin 1976; Wasserman and Faust 1994,
#' section 5.3.3). With \eqn{I_i} the number of nodes that reach node i
#' (other than itself), \eqn{n} the number of nodes, and \eqn{d(j, i)} the
#' shortest-path length from j to i,
#' \deqn{P_P(i) = \frac{I_i / (n - 1)}{\sum_{j \in I_i} d(j, i) / I_i},}{P_P(i) = (I_i / (n - 1)) / (sum_j d(j, i) / I_i),}
#' the share of nodes that reach i divided by their mean distance from it.
#' It is 0 for a node that no other node reaches and 1 for a node every
#' other node points to directly.
#'
#' @inheritParams betweenness
#'
#' @details
#' Lengths are handled as in [closeness()]: without weights every edge has
#' length 1; with weights they are distances, or strengths converted to
#' lengths 1/weight, and must be positive. Self-loops never lie on a
#' shortest path and parallel edges do not change path lengths. In an
#' undirected network every node that reaches i is also reached by it, so
#' the measure reduces to the Wasserman and Faust closeness index of
#' `closeness(x, convention = "networkx")`; for a directed network it is
#' that index computed over paths *to* each node, which is NetworkX's
#' `closeness_centrality()` on a directed graph. Neither igraph nor sna
#' has a function for it.
#'
#' @references
#' Lin, N. (1976). *Foundations of Social Research*. McGraw-Hill.
#'
#' Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and
#' Applications*. Cambridge University Press.
#' \doi{10.1017/CBO9780511815478}
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
proximity_prestige <- function(x, weighted = NULL, weights_as = NULL, threads = 0L) {
  check_netframe(x)
  threads <- resolve_threads(threads)
  w <- resolve_weighting(x, weighted, weights_as)
  ptr <- engine_graph(x)
  engine_call("proximity_prestige",
              .Call(C_closeness, ptr, FALSE, mode_code("in"), w$weighted, w$code, TRUE, 1L,
                    threads))
}

#' Edge betweenness
#'
#' The number of shortest paths that pass through each edge, summed over
#' pairs of nodes with each path weighted by the share of the pair's
#' shortest paths it represents (Girvan and Newman 2002, after Freeman
#' 1977). Edges with high betweenness join parts of the network that have
#' few other connections; removing them in order of betweenness is the
#' Girvan-Newman community method.
#'
#' @inheritParams betweenness
#' @param normalized Whether to divide by the number of pairs of nodes:
#'   \eqn{n(n-1)} in directed and \eqn{n(n-1)/2} in undirected networks.
#'
#' @details
#' Values are computed with Brandes' (2008) accumulation for edges, one
#' search per node, run in parallel. In undirected networks each unordered
#' pair of nodes is counted once. Parallel edges are distinct edges and
#' distinct shortest paths, so two parallel edges each carry half the
#' paths one edge would; self-loops never lie on a shortest path and score
#' 0. Weights are lengths, as in [betweenness()]: distances, or strengths
#' converted to 1/weight, and must be positive.
#'
#' Unnormalized values match igraph's `edge_betweenness()` and NetworkX's
#' `edge_betweenness_centrality(normalized = FALSE)`; normalized values
#' match NetworkX's default. sna has no edge betweenness function.
#'
#' @references
#' Freeman, L. C. (1977). A set of measures of centrality based on
#' betweenness. *Sociometry*, 40(1), 35-41. \doi{10.2307/3033543}
#'
#' Girvan, M., & Newman, M. E. J. (2002). Community structure in social and
#' biological networks. *Proceedings of the National Academy of Sciences*,
#' 99(12), 7821-7826. \doi{10.1073/pnas.122653799}
#'
#' Brandes, U. (2008). On variants of shortest-path betweenness centrality
#' and their generic computation. *Social Networks*, 30(2), 136-145.
#' \doi{10.1016/j.socnet.2007.11.001}
#'
#' @return A numeric vector aligned with the rows of `x$edges`.
#' @export
edge_betweenness <- function(x, weighted = NULL, normalized = FALSE, threads = 0L,
                             weights_as = NULL) {
  check_netframe(x)
  normalized <- check_flag(normalized, "normalized")
  threads <- resolve_threads(threads)
  w <- resolve_weighting(x, weighted, weights_as)
  ptr <- engine_graph(x)
  engine_call("edge_betweenness",
              .Call(C_edge_betweenness, ptr, w$weighted, w$code, normalized, threads))
}

#' Information centrality
#'
#' Stephenson and Zelen's (1989) centrality from the information carried
#' by every path between a node and the others, not only the shortest.
#' With \eqn{A} the symmetric adjacency matrix, \eqn{B = L + J} (the
#' Laplacian plus a matrix of ones) and \eqn{C = B^{-1}},
#' \deqn{I_i = \frac{1}{c_{ii} + (T - 2R)/n},}{I_i = 1 / (c_ii + (T - 2R)/n),}
#' where \eqn{T} is the trace of \eqn{C} and \eqn{R} any of its row sums.
#'
#' @inheritParams betweenness
#' @param weighted Whether to use edge weights as tie strengths
#'   (conductances: a stronger tie carries more information), whatever the
#'   network's `weight_type`. `NULL` uses them when the network has a
#'   weight column. They must be positive.
#' @param convention `"sna"` returns \eqn{I_i} as above, the value of
#'   `sna::infocent()` and UCINET; `"networkx"` returns \eqn{I_i / n}, the
#'   value of NetworkX's `information_centrality()`, which computes
#'   current-flow closeness (Brandes and Fleischer 2005 show the two agree
#'   up to that factor).
#'
#' @details
#' The measure is defined for undirected ties: a directed network is
#' symmetrized with \eqn{a_{ij} = \max(a_{ij}, a_{ji})}, sna's default for
#' `infocent()`. Self-loops are ignored. Without weights parallel edges
#' count once (the adjacency matrix has a 1); with weights
#' \eqn{a_{ij}} is the total weight of the edges between i and j.
#'
#' Isolated nodes score 0 and are left out of the matrix, as in sna, so
#' \eqn{n} is the number of nodes with ties. The formula needs those nodes
#' to be connected (otherwise \eqn{B} is singular); when they are not,
#' every node with ties gets `NaN`. NetworkX refuses such networks
#' altogether, so under the `"networkx"` convention they are an error, and
#' \eqn{n} is the number of nodes in the network, with no isolates allowed.
#'
#' @references
#' Stephenson, K., & Zelen, M. (1989). Rethinking centrality: Methods and
#' examples. *Social Networks*, 11(1), 1-37.
#' \doi{10.1016/0378-8733(89)90016-6}
#'
#' Brandes, U., & Fleischer, D. (2005). Centrality measures based on current
#' flow. In *STACS 2005* (Lecture Notes in Computer Science 3404, pp.
#' 533-544). Springer. \doi{10.1007/978-3-540-31856-9_44}
#'
#' Butts, C. T. (2008). Social network analysis with sna. *Journal of
#' Statistical Software*, 24(6), 1-51. \doi{10.18637/jss.v024.i06}
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @examples
#' g <- netframe(data.frame(from = c(1, 1, 1, 2), to = c(2, 3, 4, 3)),
#'               directed = FALSE)
#' information_centrality(g)                           # sna::infocent(), UCINET
#' information_centrality(g, convention = "networkx")  # nx.information_centrality()
#' @export
information_centrality <- function(x, weighted = NULL, convention = c("sna", "networkx")) {
  check_netframe(x)
  convention <- match.arg(convention)
  weighted <- positions_weighted(x, weighted)
  ptr <- engine_graph(x)
  out <- engine_call("information_centrality", .Call(C_information_centrality, ptr, weighted))
  if (convention == "networkx") {
    n <- nrow(x$nodes)
    if (n > 1L && (components(x, "weak")$no > 1L || anyNA(out))) {
      stop("The network is not connected; NetworkX does not compute information ",
           "centrality for such networks. Use `convention = \"sna\"`, which gives 0 to ",
           "isolates, or analyse components separately.", call. = FALSE)
    }
    out <- out / n
  }
  out
}

#' Bridging centrality
#'
#' Hwang et al.'s (2008) bridging centrality finds nodes that lie between
#' densely connected regions: it is betweenness times the bridging
#' coefficient
#' \deqn{BC(v) = \frac{1 / d(v)}{\sum_{u \in N(v)} 1 / d(u)},}{BC(v) = (1 / d(v)) / sum_u 1 / d(u),}
#' where \eqn{d} is the number of distinct neighbours and the sum runs
#' over the neighbours of v. The coefficient is large for nodes of low
#' degree whose neighbours have high degree.
#'
#' @inheritParams betweenness
#'
#' @details
#' The bridging coefficient ignores direction, weights, self-loops, and
#' parallel edges, and is `NaN` for isolated nodes. `bridging_centrality()`
#' multiplies it by unnormalized shortest-path betweenness computed on
#' binary ties ([betweenness()] with `weighted = FALSE`), which respects
#' direction in directed networks, as Hwang et al. define it on the
#' undirected networks they study. Neither igraph, NetworkX, nor sna
#' implements the measure.
#'
#' @references
#' Hwang, W., Kim, T., Ramanathan, M., & Zhang, A. (2008). Bridging
#' centrality: Graph mining from element level to group level. In
#' *Proceedings of the 14th ACM SIGKDD International Conference on
#' Knowledge Discovery and Data Mining* (pp. 336-344). ACM.
#' \doi{10.1145/1401890.1401934}
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
bridging_centrality <- function(x, threads = 0L) {
  check_netframe(x)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  b <- engine_call("bridging_centrality",
                   .Call(C_betweenness, ptr, FALSE, 0L, FALSE, threads))
  b * engine_call("bridging_centrality", .Call(C_bridging_coefficient, ptr))
}

#' @rdname bridging_centrality
#' @export
bridging_coefficient <- function(x) {
  check_netframe(x)
  ptr <- engine_graph(x)
  engine_call("bridging_coefficient", .Call(C_bridging_coefficient, ptr))
}

#' Centralization of any scores
#'
#' Freeman's (1978) centralization index for a vector of node scores:
#' \eqn{\sum_i (c_{\max} - c_i)}, the sum of differences between the
#' largest score and every score, divided by `theoretical_max`, the largest
#' value that sum can take for a network with the same number of nodes
#' under the measure's definition. [centralization()] computes the index
#' with the theoretical maxima of degree, betweenness, closeness, and
#' eigenvector centrality; this function serves any other measure whose
#' maximum is known, or reports the unnormalized sum when it is not.
#'
#' @param scores A numeric vector of node scores.
#' @param theoretical_max The largest possible value of the sum of
#'   differences for `length(scores)` nodes, or `NULL` to return the sum
#'   itself. It depends on the measure, its normalization, and the
#'   network's direction; see the table in [centralization()].
#'
#' @references
#' Freeman, L. C. (1978). Centrality in social networks: Conceptual
#' clarification. *Social Networks*, 1(3), 215-239.
#' \doi{10.1016/0378-8733(78)90021-7}
#'
#' @return A number: `NaN` if `scores` is empty or contains `NaN`.
#' @examples
#' g <- netframe(data.frame(from = c(1, 1, 1, 1), to = c(2, 3, 4, 5)), directed = FALSE)
#' d <- degree(g)
#' centralize(d, theoretical_max = (5 - 1) * (5 - 2))  # a star: 1
#' centralization(g, "degree")                        # the same
#' @export
centralize <- function(scores, theoretical_max = NULL) {
  if (!is.numeric(scores)) {
    stop("`scores` must be a numeric vector.", call. = FALSE)
  }
  if (!is.null(theoretical_max) &&
      (!is.numeric(theoretical_max) || length(theoretical_max) != 1L ||
       is.na(theoretical_max) || theoretical_max <= 0)) {
    stop("`theoretical_max` must be a single positive number, or NULL.", call. = FALSE)
  }
  if (!length(scores) || anyNA(scores)) return(NaN)
  raw <- sum(max(scores) - scores)
  if (is.null(theoretical_max)) raw else raw / theoretical_max
}
