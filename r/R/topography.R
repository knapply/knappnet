#' Edge density
#'
#' The number of edges divided by the number of possible edges: \eqn{n(n-1)}
#' in directed and \eqn{n(n-1)/2} in undirected networks, or \eqn{n^2} and
#' \eqn{n(n+1)/2} when self-loops count as possible edges.
#'
#' @param x A [netframe()].
#' @param loops Whether self-loops count as possible edges.
#' @param convention How edges are counted; see [conventions]. `"igraph"`
#'   and `"networkx"` count every edge, so parallel edges and loops can push
#'   the density above 1; `"sna"` and `"ucinet"` count distinct ties between
#'   distinct nodes (the adjacency matrix without its diagonal).
#'
#' @return A number.
#' @export
graph_density <- function(x, loops = FALSE, convention = c("igraph", "networkx", "sna", "ucinet")) {
  check_netframe(x)
  convention <- match.arg(convention)
  loops <- check_flag(loops, "loops")
  if (convention %in% c("sna", "ucinet")) {
    if (loops) {
      stop("The sna and UCINET conventions ignore self-loops; use `loops = FALSE`.",
           call. = FALSE)
    }
    ptr <- engine_graph_simple(x)
    return(engine_call("graph_density", .Call(C_density, ptr, FALSE), edge_rows = FALSE))
  }
  if (convention == "networkx") {
    if (loops) {
      stop("NetworkX's density() has no option for self-loops; use `loops = FALSE`.",
           call. = FALSE)
    }
    if (nrow(x$nodes) == 1L) return(0)
  }
  ptr <- engine_graph(x)
  engine_call("graph_density", .Call(C_density, ptr, loops))
}

#' Connected components
#'
#' `components()` finds the weak or strong components of a network.
#' `fragmentation()` is the proportion of ordered pairs of distinct nodes
#' that lie in different components (Borgatti 2006):
#' \deqn{F = 1 - \frac{\sum_k s_k (s_k - 1)}{n (n - 1)}}{F = 1 - sum_k s_k (s_k - 1) / (n (n - 1))}
#' where \eqn{s_k} are the component sizes.
#'
#' In a weak component every node can reach every other when the direction
#' of edges is ignored; in a strong component every node reaches every
#' other along directed paths. Undirected networks have one kind of
#' component, whichever `mode` is chosen. Isolated nodes are components of
#' size 1. Weights, parallel edges, and self-loops do not matter.
#'
#' Strong components are found with Tarjan's (1972) algorithm, written
#' without recursion, in time linear in the size of the network; weak
#' components with union-find over the edge list (Tarjan 1975), in
#' near-linear time. Components are numbered in order
#' of their first node in `x$nodes`. Membership matches igraph's
#' `components()` and NetworkX's component functions up to the numbering.
#'
#' @inheritParams graph_density
#' @param mode `"weak"` or `"strong"`.
#'
#' @references
#' Tarjan, R. (1972). Depth-first search and linear graph algorithms.
#' *SIAM Journal on Computing*, 1(2), 146-160. \doi{10.1137/0201010}
#'
#' Tarjan, R. E. (1975). Efficiency of a good but not linear set union
#' algorithm. *Journal of the ACM*, 22(2), 215-225.
#'
#' Borgatti, S. P. (2006). Identifying sets of key players in a social
#' network. *Computational and Mathematical Organization Theory*, 12(1),
#' 21-34. \doi{10.1007/s10588-006-7084-x}
#'
#' @return `components()` returns a list with `membership` (the component of
#'   each node, aligned with the rows of `x$nodes`), `csize` (the size of
#'   each component), and `no` (the number of components), as igraph does.
#'   `fragmentation()` returns a number, `NaN` for networks with fewer than
#'   two nodes.
#' @export
components <- function(x, mode = c("weak", "strong")) {
  check_netframe(x)
  mode <- match.arg(mode)
  ptr <- engine_graph(x)
  m <- engine_call("components",
                   .Call(C_components, ptr, if (mode == "strong") 1L else 0L)) + 1L
  csize <- tabulate(m, nbins = if (length(m)) max(m) else 0L)
  list(membership = m, csize = csize, no = length(csize))
}

#' @rdname components
#' @export
fragmentation <- function(x, mode = c("weak", "strong")) {
  check_netframe(x)
  mode <- match.arg(mode)
  ptr <- engine_graph(x)
  engine_call("fragmentation", .Call(C_fragmentation, ptr, if (mode == "strong") 1L else 0L))
}

#' Closeness and harmonic centrality
#'
#' Closeness is the reciprocal of the sum of a node's shortest-path lengths
#' to the other nodes (Sabidussi 1966; Freeman 1978). Harmonic centrality
#' sums the reciprocals of those lengths instead, so unreachable nodes
#' simply contribute 0 (Marchiori and Latora 2000; Rochat 2009).
#'
#' @inheritParams betweenness
#' @param mode For directed networks: `"out"` uses paths from the node,
#'   `"in"` paths to it, and `"all"` ignores direction. `NULL` uses the
#'   convention's default: `"in"` for NetworkX and `"out"` otherwise.
#' @param normalized Whether to normalize (see Details). `NULL` uses the
#'   convention's default: unnormalized under igraph's convention and for
#'   NetworkX's harmonic centrality, normalized otherwise. NetworkX's
#'   closeness exists only in normalized form.
#' @param convention Whose definition to follow; see Details and
#'   [conventions].
#'
#' @details
#' Let \eqn{n} be the number of nodes, \eqn{r} the number of other nodes a
#' node reaches, and \eqn{S} the sum of the lengths to them. Tools differ in
#' how they handle nodes that cannot be reached:
#' \describe{
#'   \item{`"igraph"`}{Only reachable nodes count: \eqn{1/S}, or \eqn{r/S}
#'     normalized. A node that reaches no other gets `NaN`.}
#'   \item{`"networkx"`}{The Wasserman and Faust (1994, p. 201) index
#'     \eqn{(r/S)(r/(n-1))}, which is always normalized; 0 if \eqn{r = 0}.}
#'   \item{`"sna"`}{\eqn{(n-1)/S} (or \eqn{1/S} unnormalized), and 0 for any
#'     node that does not reach every other node.}
#'   \item{`"ucinet"`}{\eqn{(n-1)/S'}, where \eqn{S'} counts every missing
#'     length as one more than the longest shortest path in the network,
#'     UCINET's default ("max observed distance plus 1").}
#' }
#' All four agree on strongly connected networks. With the sna and UCINET
#' conventions the network is read as an adjacency matrix (parallel edges
#' collapse, self-loops are ignored) and weights are not used, as in those
#' tools.
#'
#' Harmonic centrality is \eqn{\sum 1/d}, divided by \eqn{n-1} when
#' normalized; this is igraph's `harmonic_centrality()`, NetworkX's
#' `harmonic_centrality()` (never normalized, paths to the node), sna's
#' `closeness(cmode = "suminvdir")`, and UCINET's reciprocal closeness (both
#' normalized).
#'
#' Self-loops never lie on a shortest path, and parallel edges do not change
#' path lengths (with weights the shortest one counts). Weights are handled
#' as in [betweenness()]: distances, or strengths converted to lengths
#' 1/weight, and they must be positive. Searches from each node run in
#' parallel.
#'
#' @references
#' Sabidussi, G. (1966). The centrality index of a graph. *Psychometrika*,
#' 31(4), 581-603. \doi{10.1007/BF02289527}
#'
#' Freeman, L. C. (1978). Centrality in social networks: Conceptual
#' clarification. *Social Networks*, 1(3), 215-239.
#' \doi{10.1016/0378-8733(78)90021-7}
#'
#' Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and
#' Applications*. Cambridge University Press.
#' \doi{10.1017/CBO9780511815478}
#'
#' Marchiori, M., & Latora, V. (2000). Harmony in the small-world.
#' *Physica A*, 285(3-4), 539-546. \doi{10.1016/S0378-4371(00)00311-3}
#'
#' Rochat, Y. (2009). Closeness centrality extended to unconnected graphs:
#' The harmonic centrality index. Applications of Social Network Analysis
#' (ASNA), Zurich.
#'
#' Boldi, P., & Vigna, S. (2014). Axioms for centrality. *Internet
#' Mathematics*, 10(3-4), 222-262. \doi{10.1080/15427951.2013.865686}
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
closeness <- function(x, mode = NULL, weighted = NULL, normalized = NULL, threads = 0L,
                      weights_as = NULL,
                      convention = c("igraph", "networkx", "sna", "ucinet")) {
  closeness_common(x, FALSE, mode, weighted, normalized, threads, weights_as,
                   match.arg(convention), "closeness")
}

#' @rdname closeness
#' @export
harmonic_centrality <- function(x, mode = NULL, weighted = NULL, normalized = NULL,
                                threads = 0L, weights_as = NULL,
                                convention = c("igraph", "networkx", "sna", "ucinet")) {
  closeness_common(x, TRUE, mode, weighted, normalized, threads, weights_as,
                   match.arg(convention), "harmonic_centrality")
}

closeness_common <- function(x, harmonic, mode, weighted, normalized, threads, weights_as,
                             convention, fn) {
  check_netframe(x)
  if (is.null(mode)) mode <- if (convention == "networkx") "in" else "out"
  mode <- check_choice(mode, c("out", "in", "all"), "mode")
  check_flag_or_null(weighted, "weighted")
  check_flag_or_null(normalized, "normalized")
  threads <- resolve_threads(threads)
  if (is.null(normalized)) {
    normalized <- !(convention == "igraph" || (harmonic && convention == "networkx"))
  }
  binary <- convention %in% c("sna", "ucinet")
  if (binary && isTRUE(weighted)) {
    stop(sprintf("The %s convention computes closeness on binary ties; use `weighted = FALSE`.",
                 convention), call. = FALSE)
  }
  w <- resolve_weighting(x, if (binary) FALSE else weighted, weights_as)
  unreachable <- switch(convention, igraph = 0L, networkx = 1L, sna = 2L, ucinet = 3L)
  ptr <- if (binary) engine_graph_simple(x) else engine_graph(x)
  engine_call(fn, .Call(C_closeness, ptr, harmonic, mode_code(mode), w$weighted, w$code,
                        normalized, unreachable, threads), edge_rows = !binary)
}

#' Eigenvector centrality
#'
#' Each node's score is proportional to the sum of its neighbors' scores:
#' the principal eigenvector of the adjacency matrix (Bonacich 1972).
#'
#' @inheritParams closeness
#' @param mode For directed networks: `"in"` sums the scores of nodes with
#'   ties *to* the node (igraph with `directed = TRUE`, NetworkX), `"out"`
#'   the nodes it has ties to (sna), and `"all"` ignores direction (igraph's
#'   R default, `directed = FALSE`, and UCINET, which symmetrizes). `NULL`
#'   uses the convention's choice.
#' @param weighted Whether to use edge weights, as given (larger means a
#'   stronger tie), whatever the network's `weight_type`. `NULL` uses them
#'   when the network has a weight column. They must not be negative.
#' @param scale `"max"` (largest score 1, igraph) or `"unit"` (Euclidean
#'   norm 1, NetworkX, sna, and UCINET). `NULL` uses the convention's
#'   choice.
#' @param tol Convergence tolerance: iteration stops when no score (scaled
#'   so the largest is 1) changes by more than `tol`.
#' @param max_iter Iteration cap. Not converging within it is an error.
#'
#' @details
#' Scores are computed by power iteration on \eqn{A + I} (von Mises and
#' Pollaczek-Geiringer 1929), starting from equal scores. The shift by the
#' identity changes no eigenvector but guarantees convergence on bipartite
#' and other periodic networks, where plain iteration on \eqn{A} cycles.
#' Each iteration costs time linear in the size of the network and runs in
#' parallel on large networks.
#'
#' Parallel edges add up. In undirected networks a self-loop adds 2 to its
#' diagonal entry (two edge endpoints) under the igraph convention and 1
#' under NetworkX's; directed loops add 1. The sna and UCINET conventions
#' read the network as an adjacency matrix: parallel edges collapse
#' (keeping the largest weight) and self-loops are ignored.
#'
#' The eigenvector is unique (up to scale) and positive only when the
#' network is connected (strongly connected, for directed networks with
#' `mode` `"in"` or `"out"`); otherwise a warning is given. For undirected
#' networks (and `mode = "all"`) each component is then computed
#' separately: components whose largest eigenvalue is below the network's
#' score 0, and components that share it each contribute their unit
#' eigenvector \eqn{u} times \eqn{\sum u}, the limit of power iteration
#' from equal scores. For directed networks the result is the limit of the
#' iteration; if components sharing the largest eigenvalue are joined by a
#' directed path, the iteration does not converge and an error is reported.
#' Convergence is also slow when the two largest eigenvalues are close.
#' NetworkX refuses disconnected networks, and so does the `"networkx"`
#' convention. When more than one component has ties, UCINET switches to
#' a different measure (lambda-squared scores computed within components),
#' which is not implemented here, so the `"ucinet"` convention stops.
#'
#' Networks without cycles have no positive eigenvalue. Following igraph,
#' nodes whose score feeds no other node's score (nodes without outgoing
#' ties for `mode = "in"`, without incoming ties for `"out"`, and every node
#' in a network without edges) then score 1 before scaling and all others 0.
#'
#' Results match igraph's `eigen_centrality()`, NetworkX's
#' `eigenvector_centrality_numpy()`, and sna's `evcent()` to within their
#' convergence tolerances on connected networks.
#'
#' @references
#' Bonacich, P. (1972). Factoring and weighting approaches to status scores
#' and clique identification. *Journal of Mathematical Sociology*, 2(1),
#' 113-120. \doi{10.1080/0022250X.1972.9989806}
#'
#' Bonacich, P. (1987). Power and centrality: A family of measures.
#' *American Journal of Sociology*, 92(5), 1170-1182. \doi{10.1086/228631}
#'
#' Mises, R. v., & Pollaczek-Geiringer, H. (1929). Praktische Verfahren der
#' Gleichungsauflösung. *Zeitschrift für Angewandte Mathematik und
#' Mechanik*, 9(2), 152-164. \doi{10.1002/zamm.19290090206}
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
eigenvector_centrality <- function(x, mode = NULL, weighted = NULL, scale = NULL,
                                   convention = c("igraph", "networkx", "sna", "ucinet"),
                                   tol = 1e-12, max_iter = 10000L, threads = 0L) {
  check_netframe(x)
  convention <- match.arg(convention)
  if (is.null(mode)) {
    mode <- switch(convention, igraph = "in", networkx = "in", sna = "out", ucinet = "all")
  }
  mode <- check_choice(mode, c("in", "out", "all"), "mode")
  if (is.null(scale)) scale <- if (convention == "igraph") "max" else "unit"
  scale <- check_choice(scale, c("max", "unit"), "scale")
  check_flag_or_null(weighted, "weighted")
  threads <- resolve_threads(threads)
  weighted <- if (is.null(weighted)) !is.null(x$roles$weight) else isTRUE(weighted)
  if (weighted) require_weights(x)
  if (!is.numeric(tol) || length(tol) != 1L || !isTRUE(tol > 0)) {
    stop("`tol` must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !isTRUE(max_iter >= 1) ||
      max_iter != floor(max_iter)) {
    stop("`max_iter` must be a positive whole number.", call. = FALSE)
  }

  loops <- switch(convention, igraph = 0L, networkx = 1L, sna = 2L, ucinet = 2L)
  ptr <- switch(convention,
    sna = engine_graph_simple(x),
    ucinet = engine_graph_simple(x, undirected = TRUE),
    engine_graph(x)
  )
  if (convention == "ucinet") {
    csize <- components(x, "weak")$csize
    if (sum(csize > 1L) > 1L) {
      stop("More than one component has ties. UCINET then reports lambda-squared ",
           "scores computed within components, which are not implemented; use ",
           "another convention or analyse components separately.", call. = FALSE)
    }
  }
  res <- engine_call("eigenvector_centrality",
    .Call(C_eigenvector, ptr, mode_code(mode), weighted, if (scale == "max") 0L else 1L,
          loops, as.double(tol), as.double(max_iter), threads),
    edge_rows = convention %in% c("igraph", "networkx"))
  if (!res[[3]] && length(res[[1]]) > 1L) {
    if (convention == "networkx") {
      stop("The network is not ", if (x$directed && mode != "all") "strongly ",
           "connected; NetworkX does not compute eigenvector centrality for such networks.",
           call. = FALSE)
    }
    warning("The network is not ", if (x$directed && mode != "all") "strongly ",
            "connected, so eigenvector centrality is zero outside the dominant ",
            "component(s) and may not be unique.", call. = FALSE)
  }
  res[[1]]
}

#' Freeman centralization
#'
#' How far the distribution of a centrality is from the most unequal one
#' possible: the sum of differences between the largest score and every
#' score, divided by the largest value that sum can take in a network of
#' the same size, which a star attains (Freeman 1978).
#'
#' @inheritParams closeness
#' @param measure `"degree"`, `"betweenness"`, `"closeness"`, or
#'   `"eigenvector"`.
#' @param mode For directed networks: which degree, or which direction of
#'   paths for closeness. Ignored for betweenness and undirected networks,
#'   and always `"all"` for eigenvector centrality. `NULL` uses igraph's
#'   defaults (`"all"` for degree and eigenvector centrality, `"out"` for
#'   closeness), or out-degree under the UCINET convention, which reports
#'   only in- and out-degree centralization.
#' @param loops Whether self-loops count toward degree (igraph convention
#'   only).
#' @param normalized Whether to divide by the theoretical maximum.
#' @param convention Whose theoretical maxima and closeness definition to
#'   use: `"igraph"`, `"sna"`, or `"ucinet"`. NetworkX has no
#'   centralization index, so `"networkx"` is an error.
#'
#' @details
#' Degree uses edge counts and betweenness is unnormalized shortest-path
#' betweenness; weights are not used. Closeness centralization uses
#' normalized closeness as each tool defines it (see [closeness()]); under
#' the igraph convention it is `NaN` if some node reaches no other, and
#' under the UCINET convention unless every node reaches every other, as in
#' those tools.
#'
#' Eigenvector centralization follows igraph's `centr_eigen()`: scores
#' scaled so the largest is 1 (see [eigenvector_centrality()]), direction
#' ignored, and the theoretical maximum \eqn{n - 2}, attained by a network
#' of one edge and \eqn{n - 2} isolated nodes (scores 1, 1, and 0), which
#' igraph's documentation gives as the most centralized network for this
#' measure; a star then scores less than 1 (0.75 for ten nodes). It is
#' `NaN` for fewer than three nodes. When more than one component has
#' ties the eigenvector depends on more than the network (see
#' [eigenvector_centrality()]; a warning is given), and so does its
#' centralization. sna and UCINET compute eigenvector
#' scores with unit norm and use maxima of their own that could not be
#' verified from their documentation, so only the igraph convention is
#' provided for this measure.
#'
#' The theoretical maxima are:
#' \tabular{llll}{
#'   Measure \tab igraph \tab sna \tab UCINET \cr
#'   degree, undirected \tab \eqn{(n-1)(n-2)}, or \eqn{(n-1)n} with loops \tab \eqn{(n-1)(n-2)} \tab \eqn{(n-1)(n-2)} \cr
#'   in- or out-degree \tab \eqn{(n-1)^2}, or \eqn{(n-1)n} with loops \tab \eqn{(n-1)^2} \tab \eqn{(n-1)^2} \cr
#'   total degree \tab \eqn{2(n-1)(n-2)}, or \eqn{2(n-1)^2} with loops \tab \eqn{2(n-1)(n-2)} \tab (not reported) \cr
#'   betweenness \tab \eqn{(n-1)^2(n-2)}, halved if undirected \tab same \tab same \cr
#'   closeness, undirected or `"all"` \tab \eqn{(n-1)(n-2)/(2n-3)} \tab same \tab same \cr
#'   closeness, in or out \tab \eqn{(n-1)^2/n} \tab same \tab \eqn{(n-1)(n-2)/(2n-3)} \cr
#'   eigenvector (largest score 1) \tab \eqn{n-2} \tab (not provided) \tab (not provided)
#' }
#' With the sna and UCINET conventions the network is read as an adjacency
#' matrix: parallel edges collapse and self-loops are ignored. The igraph
#' convention matches `centr_degree()`, `centr_betw()`, and `centr_clo()`.
#' UCINET reports centralization as a percentage in some routines; this
#' function always returns a proportion. UCINET versions before 6.849 used
#' \eqn{(n-1)(n-2)} for directed degree centralization.
#'
#' @references
#' Freeman, L. C. (1978). Centrality in social networks: Conceptual
#' clarification. *Social Networks*, 1(3), 215-239.
#' \doi{10.1016/0378-8733(78)90021-7}
#'
#' @return A number.
#' @export
centralization <- function(x, measure = c("degree", "betweenness", "closeness", "eigenvector"),
                           mode = NULL,
                           loops = TRUE, normalized = TRUE,
                           convention = c("igraph", "networkx", "sna", "ucinet"),
                           threads = 0L) {
  check_netframe(x)
  measure <- match.arg(measure)
  convention <- match.arg(convention)
  if (convention == "networkx") {
    stop("NetworkX has no centralization index; use the igraph, sna, or UCINET convention.",
         call. = FALSE)
  }
  if (measure == "eigenvector" && convention != "igraph") {
    stop("Eigenvector centralization is provided for the igraph convention only.",
         call. = FALSE)
  }
  if (is.null(mode)) {
    mode <- if (measure != "eigenvector" && (measure == "closeness" || convention == "ucinet")) {
      "out"
    } else {
      "all"
    }
  }
  mode <- check_choice(mode, c("all", "out", "in"), "mode")
  loops <- check_flag(loops, "loops")
  normalized <- check_flag(normalized, "normalized")
  threads <- resolve_threads(threads)
  code <- switch(measure, degree = 0L, betweenness = 1L, closeness = 2L, eigenvector = 3L)
  conv <- switch(convention, igraph = 0L, sna = 2L, ucinet = 3L)
  ptr <- engine_graph(x)
  if (measure == "eigenvector" && nrow(x$nodes) > 1L) {
    csize <- components(x, "weak")$csize
    if (sum(csize > 1L) > 1L) {
      warning("More than one component has ties, so eigenvector centrality may not be ",
              "unique and neither is its centralization.", call. = FALSE)
    }
  }
  engine_call("centralization", .Call(C_centralization, ptr, code, mode_code(mode), loops,
                                      normalized, conv, threads))
}

mode_code <- function(mode) switch(mode, out = 1L, `in` = 2L, all = 3L)

# The network as an adjacency matrix, the way sna and UCINET read it:
# parallel edges collapse into one tie (keeping the largest weight, or the
# smallest with keep = "min") and self-loops are dropped. Cached alongside
# the engine graph.
engine_graph_simple <- function(x, undirected = FALSE, keep = c("max", "min")) {
  keep <- match.arg(keep)
  key <- paste("simple", isTRUE(undirected), keep, sep = "_")
  # engine_graph() first: it empties the cache when the handle is stale.
  base <- engine_graph(x)
  ptr <- x$cache[[key]]
  if (is.null(ptr) || !.Call(C_graph_alive, ptr)) {
    ptr <- engine_call("netframe", .Call(C_graph_simplify, base, isTRUE(undirected),
                                         if (keep == "max") 0L else 1L))
    assign(key, ptr, envir = x$cache)
  }
  ptr
}
