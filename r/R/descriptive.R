#' Reciprocity
#'
#' How far ties in a directed network are returned. Reciprocity is a
#' property of dyads, so the network is read as an adjacency matrix:
#' parallel arcs count once and self-loops are ignored.
#'
#' @param x A [netframe()].
#' @param mode Which index; see Details.
#'
#' @details
#' With \eqn{M} mutual, \eqn{A} asymmetric, and \eqn{N} null dyads among
#' the \eqn{n(n-1)/2} pairs of nodes:
#' \describe{
#'   \item{`"default"`}{The share of arcs that are reciprocated,
#'     \eqn{2M / (2M + A)} (Newman, Forrest, and Balthrop 2002). igraph's
#'     `reciprocity()`, sna's `grecip(measure = "edgewise")`, and
#'     NetworkX's `overall_reciprocity()`.}
#'   \item{`"ratio"`}{Mutual dyads over connected dyads, \eqn{M / (M + A)}.
#'     igraph's `mode = "ratio"` and sna's `"dyadic.nonnull"`.}
#'   \item{`"dyadic"`}{Symmetric dyads (mutual or null) over all dyads,
#'     \eqn{(M + N) / (M + A + N)}. sna's `"dyadic"`.}
#'   \item{`"correlation"`}{Garlaschelli and Loffredo's (2004)
#'     \eqn{\rho = (r - p) / (1 - p)}{rho = (r - p) / (1 - p)}, where
#'     \eqn{r} is the `"default"` index and \eqn{p} the density: the
#'     correlation between \eqn{x_{ij}}{x_ij} and \eqn{x_{ji}}{x_ji} over
#'     ordered pairs of distinct nodes, which is sna's `"correlation"`. It
#'     is 0 when ties are returned as often as chance would have it,
#'     positive when more often, and negative when less.}
#' }
#' A ratio whose denominator is 0 (a network without arcs, or a complete one
#' for `"correlation"`) is `NaN`; sna reports 1 for the correlation in those
#' cases. An undirected network is read as a
#' symmetric directed one, so every tie is reciprocated. Weights are
#' ignored.
#'
#' @references
#' Newman, M. E. J., Forrest, S., & Balthrop, J. (2002). Email networks and
#' the spread of computer viruses. *Physical Review E*, 66(3), 035101.
#' \doi{10.1103/PhysRevE.66.035101}
#'
#' Garlaschelli, D., & Loffredo, M. I. (2004). Patterns of link reciprocity
#' in directed networks. *Physical Review Letters*, 93(26), 268701.
#' \doi{10.1103/PhysRevLett.93.268701}
#'
#' Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and
#' Applications*. Cambridge University Press.
#' \doi{10.1017/CBO9780511815478}
#'
#' @return A number.
#' @export
reciprocity <- function(x, mode = c("default", "ratio", "dyadic", "correlation")) {
  check_netframe(x)
  mode <- match.arg(mode)
  code <- switch(mode, default = 0L, ratio = 1L, dyadic = 2L, correlation = 3L)
  .Call(C_reciprocity, engine_graph(x), code)
}

#' Transitivity and clustering coefficients
#'
#' `transitivity()` is the share of connected triples that are closed: three
#' times the number of triangles divided by the number of paths of length
#' two (Newman 2003). `clustering()` gives each node's local clustering
#' coefficient, the share of pairs of its neighbors that are themselves
#' adjacent (Watts and Strogatz 1998), and `average_clustering()` their
#' mean.
#'
#' @inheritParams reciprocity
#' @param weighted Whether to compute a weighted local coefficient (see
#'   Details). The binary coefficient is the default, whatever the network's
#'   weights, as in igraph and NetworkX.
#' @param method The weighted coefficient: `"barrat"` (Barrat et al. 2004,
#'   igraph) or `"onnela"` (Onnela et al. 2005, NetworkX). `NULL` uses the
#'   convention's choice.
#' @param convention How a coefficient that is 0 / 0 is reported, and which
#'   weighted coefficient is the default: `"igraph"` (the default) or
#'   `"networkx"`; see Details and [conventions].
#' @param threads Number of threads; `0` uses all available cores.
#'
#' @details
#' All three use the simple undirected reading of the network: directions
#' are ignored, parallel edges count once, and self-loops are ignored, as in
#' igraph's `transitivity()`. A node's degree \eqn{k} is its number of
#' distinct neighbors other than itself. Triangles are listed once each
#' after ordering nodes by degree (Chiba and Nishizeki 1985), in parallel.
#'
#' The coefficient of a node with fewer than two neighbors, and the
#' transitivity of a network without connected triples, are 0 / 0. Under
#' the igraph convention they are `NaN`, and `average_clustering()` leaves
#' those nodes out (igraph's `transitivity(type = "average", isolates =
#' "NaN")`). Under the NetworkX convention they are 0, and the average
#' includes them as 0 (NetworkX's `clustering()`, `transitivity()`, and
#' `average_clustering()`; igraph's `isolates = "zero"`). The two averages
#' can differ greatly in networks with many leaves (Kaiser 2008).
#'
#' Weighted coefficients use, for each pair of adjacent nodes, the total
#' weight of the edges joining them in either direction, read as tie
#' strengths whatever the network's `weight_type`; weights must be
#' positive. With \eqn{s_i} the sum of node \eqn{i}'s weights, Barrat's
#' coefficient is \deqn{\frac{1}{s_i (k_i - 1)} \sum_{j,h}
#' \frac{w_{ij} + w_{ih}}{2} a_{ij} a_{ih} a_{jh}}{1 / (s_i (k_i - 1))
#' sum_(j,h) (w_ij + w_ih) / 2 a_ij a_ih a_jh} (igraph's `type =
#' "barrat"`), and Onnela's is \deqn{\frac{1}{k_i (k_i - 1)} \sum_{j,h}
#' (\hat w_{ij} \hat w_{ih} \hat w_{jh})^{1/3}}{1 / (k_i (k_i - 1)) sum_(j,h)
#' (w_ij w_ih w_jh)^(1/3)} with every weight divided by the largest one
#' (NetworkX's weighted `clustering()`), both summing over ordered pairs of
#' neighbors. On networks without parallel edges these are the published
#' definitions; with equal weights both reduce to the binary coefficient.
#' Results match igraph and NetworkX on every network in the validation
#' suite.
#'
#' @references
#' Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
#' 'small-world' networks. *Nature*, 393(6684), 440-442.
#' \doi{10.1038/30918}
#'
#' Newman, M. E. J. (2003). The structure and function of complex networks.
#' *SIAM Review*, 45(2), 167-256. \doi{10.1137/S003614450342480}
#'
#' Barrat, A., Barthélemy, M., Pastor-Satorras, R., & Vespignani, A.
#' (2004). The architecture of complex weighted networks. *Proceedings of
#' the National Academy of Sciences*, 101(11), 3747-3752.
#' \doi{10.1073/pnas.0400087101}
#'
#' Onnela, J.-P., Saramäki, J., Kertész, J., & Kaski, K. (2005). Intensity
#' and coherence of motifs in weighted complex networks. *Physical Review
#' E*, 71(6), 065103. \doi{10.1103/PhysRevE.71.065103}
#'
#' Kaiser, M. (2008). Mean clustering coefficients: The role of isolated
#' nodes and leafs on clustering measures for small-world networks. *New
#' Journal of Physics*, 10(8), 083042. \doi{10.1088/1367-2630/10/8/083042}
#'
#' Chiba, N., & Nishizeki, T. (1985). Arboricity and subgraph listing
#' algorithms. *SIAM Journal on Computing*, 14(1), 210-223.
#' \doi{10.1137/0214017}
#'
#' @return `transitivity()` and `average_clustering()` return a number;
#'   `clustering()` a numeric vector aligned with the rows of `x$nodes`.
#' @export
transitivity <- function(x, convention = c("igraph", "networkx"), threads = 0L) {
  check_netframe(x)
  convention <- match.arg(convention)
  .Call(C_transitivity, engine_graph(x), undefined_code(convention), resolve_threads(threads))
}

#' @rdname transitivity
#' @export
clustering <- function(x, weighted = FALSE, method = NULL,
                       convention = c("igraph", "networkx"), threads = 0L) {
  clustering_common(x, FALSE, weighted, method, match.arg(convention), threads)
}

#' @rdname transitivity
#' @export
average_clustering <- function(x, weighted = FALSE, method = NULL,
                               convention = c("igraph", "networkx"), threads = 0L) {
  clustering_common(x, TRUE, weighted, method, match.arg(convention), threads)
}

clustering_common <- function(x, average, weighted, method, convention, threads) {
  check_netframe(x)
  weighting <- 0L
  if (isTRUE(weighted)) {
    require_weights(x)
    if (is.null(method)) method <- if (convention == "networkx") "onnela" else "barrat"
    method <- match.arg(method, c("barrat", "onnela"))
    weighting <- if (method == "barrat") 1L else 2L
  }
  .Call(C_clustering, engine_graph(x), average, weighting, undefined_code(convention),
        resolve_threads(threads))
}

undefined_code <- function(convention) if (convention == "networkx") 1L else 0L

#' Shortest-path distances, mean distance, diameter, eccentricity, and radius
#'
#' `distances()` gives the length of a shortest path between every pair of
#' nodes. `mean_distance()` averages those lengths and `diameter()` is the
#' longest of them. A node's eccentricity is its largest distance to the
#' nodes it reaches, and the radius is the smallest eccentricity.
#'
#' @inheritParams betweenness
#' @param from,to Nodes to measure paths from and to, as values of the node
#'   identifier column; `NULL` uses every node.
#' @param mode For directed networks: `"out"` follows edges in their
#'   direction (row \eqn{i}, column \eqn{j} holds the distance from
#'   \eqn{i} to \eqn{j}), `"in"` against it (the distance from \eqn{j} to
#'   \eqn{i}), and `"all"` ignores direction. Ignored for undirected
#'   networks. The defaults follow igraph: `"out"` for `distances()`,
#'   `mean_distance()`, and `diameter()`, and `"all"` for `eccentricity()`
#'   and `radius()`.
#' @param unconnected What to do when some node does not reach every other:
#'   `TRUE` averages (or maximizes) over the pairs joined by a path, and
#'   `FALSE` returns `Inf`.
#'
#' @details
#' A node is at distance 0 from itself and at distance `Inf` from nodes it
#' cannot reach. Self-loops never shorten a path and parallel edges count
#' through the shortest of them. Weights are handled as in
#' [betweenness()]: distances as given, or strengths converted to lengths
#' 1/weight, and they must be positive; when every weight is the same,
#' lengths are that weight times the number of steps. The fastest exact
#' algorithm is chosen automatically (breadth-first search, Dial's bucket
#' queue, or Dijkstra's algorithm), with one search per source node run in
#' parallel. `is.finite(distances(x))` is the reachability matrix.
#'
#' `mean_distance()` is the characteristic path length of Watts and
#' Strogatz (1998), averaged over ordered pairs of distinct nodes. As in
#' igraph's `mean_distance()` and `diameter()`, `unconnected = TRUE` uses
#' only pairs joined by a path, and `unconnected = FALSE` gives `Inf` when
#' some node does not reach every other (in directed networks with
#' `mode = "out"`, when the network is not strongly connected). The mean
#' is `NaN` when no pair is joined; the diameter is 0 for a single node and
#' `NaN` for a network without nodes. NetworkX's
#' `average_shortest_path_length()` and `diameter()` agree on connected
#' networks and refuse others.
#'
#' Eccentricity ignores nodes that cannot be reached, so an isolated node
#' has eccentricity 0, as in igraph's `eccentricity()` and `radius()`
#' (Hage and Harary 1995). NetworkX's `eccentricity()` agrees on connected
#' networks and refuses others.
#'
#' Every function here matches igraph, and `distances()` also NetworkX and
#' sna's `geodist()`, on every network in the validation suite.
#'
#' @references
#' Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
#' *Numerische Mathematik*, 1, 269-271. \doi{10.1007/BF01386390}
#'
#' Dial, R. B. (1969). Algorithm 360: Shortest-path forest with topological
#' ordering. *Communications of the ACM*, 12(11), 632-633.
#' \doi{10.1145/363269.363610}
#'
#' Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
#' 'small-world' networks. *Nature*, 393(6684), 440-442.
#' \doi{10.1038/30918}
#'
#' Hage, P., & Harary, F. (1995). Eccentricity and centrality in networks.
#' *Social Networks*, 17(1), 57-63. \doi{10.1016/0378-8733(94)00248-9}
#'
#' @return `distances()` returns a numeric matrix with one row per node in
#'   `from` and one column per node in `to`, with node identifiers as
#'   dimnames. `mean_distance()`, `diameter()`, and `radius()` return a
#'   number; `eccentricity()` a numeric vector aligned with the rows of
#'   `x$nodes`.
#' @export
distances <- function(x, from = NULL, to = NULL, mode = c("out", "in", "all"),
                      weighted = NULL, weights_as = NULL, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  w <- resolve_weighting(x, weighted, weights_as)
  ids <- x$nodes[[x$roles$id]]
  rows <- node_positions(x, from, "from")
  cols <- node_positions(x, to, "to")
  d <- .Call(C_distances, engine_graph(x), if (is.null(rows)) NULL else rows - 1L,
             mode_code(mode), w$weighted, w$code, resolve_threads(threads))
  labels <- as.character(ids)
  dimnames(d) <- list(if (is.null(rows)) labels else labels[rows], labels)
  if (!is.null(cols)) d <- d[, cols, drop = FALSE]
  d
}

#' @rdname distances
#' @export
mean_distance <- function(x, mode = c("out", "in", "all"), weighted = NULL,
                          weights_as = NULL, unconnected = TRUE, threads = 0L) {
  path_summary(x, NULL, match.arg(mode), weighted, weights_as, unconnected, threads)[[6L]]
}

#' @rdname distances
#' @export
diameter <- function(x, mode = c("out", "in", "all"), weighted = NULL,
                     weights_as = NULL, unconnected = TRUE, threads = 0L) {
  path_summary(x, NULL, match.arg(mode), weighted, weights_as, unconnected, threads)[[7L]]
}

#' @rdname distances
#' @export
eccentricity <- function(x, mode = c("all", "out", "in"), weighted = NULL,
                         weights_as = NULL, threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  w <- resolve_weighting(x, weighted, weights_as)
  .Call(C_eccentricity, engine_graph(x), mode_code(mode), w$weighted, w$code,
        resolve_threads(threads))
}

#' @rdname distances
#' @export
radius <- function(x, mode = c("all", "out", "in"), weighted = NULL,
                   weights_as = NULL, threads = 0L) {
  e <- eccentricity(x, match.arg(mode), weighted, weights_as, threads)
  if (length(e) == 0L) NaN else min(e)
}

# Totals over the searches from `sources` (0-based positions, or NULL for
# every node): sources, pairs, total, longest, incomplete, mean, diameter.
path_summary <- function(x, sources, mode, weighted, weights_as, unconnected, threads) {
  check_netframe(x)
  w <- resolve_weighting(x, weighted, weights_as)
  if (!is.logical(unconnected) || length(unconnected) != 1L || is.na(unconnected)) {
    stop("`unconnected` must be TRUE or FALSE.", call. = FALSE)
  }
  .Call(C_path_summary, engine_graph(x), sources, mode_code(mode), w$weighted, w$code,
        unconnected, resolve_threads(threads))
}

# Positions (1-based) of node identifiers in x$nodes, or NULL for NULL.
# Shared by every function that takes node identifiers.
node_positions <- function(x, nodes, what = "nodes") {
  if (is.null(nodes)) return(NULL)
  ids <- x$nodes[[x$roles$id]]
  if (is.factor(ids)) ids <- as.character(ids)
  if (is.factor(nodes)) nodes <- as.character(nodes)
  pos <- match(nodes, ids)
  if (anyNA(pos)) {
    shown <- utils::head(unique(nodes[is.na(pos)]), 5L)
    stop(sprintf("%d node(s) in `%s` are not in the network, for example: %s.",
                 length(unique(nodes[is.na(pos)])), what,
                 paste(format(shown), collapse = ", ")), call. = FALSE)
  }
  pos
}

#' Reachability
#'
#' The number of other nodes each node can reach along a path.
#'
#' @inheritParams betweenness
#' @param mode For directed networks: `"out"` counts the nodes a node
#'   reaches, `"in"` the nodes that reach it, and `"all"` the nodes joined
#'   to it when directions are ignored (the size of its weak component,
#'   less one). Ignored for undirected networks.
#'
#' @details
#' Weights, parallel edges, and self-loops do not matter. Undirected
#' networks and `mode = "all"` take time linear in the size of the network;
#' otherwise one breadth-first search runs from every node, in parallel.
#' The full reachability matrix is
#' `is.finite(distances(x, mode = mode, weighted = FALSE))`. Results match
#' igraph's `ego_size()` with an unlimited order (less one) and NetworkX's
#' `descendants()` and `ancestors()`.
#'
#' @return A numeric vector aligned with the rows of `x$nodes`.
#' @export
reachability <- function(x, mode = c("out", "in", "all"), threads = 0L) {
  check_netframe(x)
  mode <- match.arg(mode)
  .Call(C_reachability, engine_graph(x), mode_code(mode), resolve_threads(threads))
}

#' Cut points and bridges
#'
#' A cut point (articulation point) is a node whose removal increases the
#' number of components; a bridge is an edge whose removal does.
#'
#' @inheritParams betweenness
#'
#' @details
#' Both use the undirected reading of the network in which every edge stays
#' an edge: directions are ignored, so two parallel edges, or two arcs in
#' opposite directions, between the same nodes are two edges and neither is
#' a bridge; self-loops are ignored. Weights do not matter. Both come from
#' one depth-first search that computes low points (Hopcroft and Tarjan
#' 1973), written without recursion so that very deep networks cannot
#' exhaust the stack, in time linear in the size of the network.
#'
#' Results match NetworkX's `articulation_points()` and `bridges()` on an
#' undirected multigraph, and igraph's `articulation_points()` and sna's
#' `cutpoints()` (on the symmetrized network). igraph's `bridges()` agrees
#' on every validation network whose undirected reading has no parallel
#' edges; how it treats parallel edges is not documented.
#'
#' @references
#' Hopcroft, J., & Tarjan, R. (1973). Algorithm 447: Efficient algorithms
#' for graph manipulation. *Communications of the ACM*, 16(6), 372-378.
#' \doi{10.1145/362248.362272}
#'
#' @return `cut_points()` returns the identifiers (values of
#'   `x$nodes[[id]]`) of the cut points, in the order of `x$nodes`.
#'   `bridges()` returns the row numbers of the bridges in `x$edges`, in
#'   increasing order.
#' @export
cut_points <- function(x) {
  check_netframe(x)
  r <- .Call(C_cut_structure, engine_graph(x))
  x$nodes[[x$roles$id]][r[[1L]]]
}

#' @rdname cut_points
#' @export
bridges <- function(x) {
  check_netframe(x)
  which(.Call(C_cut_structure, engine_graph(x))[[2L]])
}

#' Summarize a network
#'
#' The standard descriptive table of a network: its size, ties, components,
#' cohesion, distances, and degree distribution.
#'
#' @param object A [netframe()].
#' @param paths Whether to compute the mean distance and diameter, which
#'   need a shortest-path search from every node of the largest component.
#'   `NULL` (the default) computes them when that component has at most
#'   5,000 nodes and otherwise skips them, saying so.
#' @param weighted Whether the mean distance and diameter use edge weights
#'   (which must then be positive). `NULL` uses them when the network has a
#'   weight column; see [distances()].
#' @param threads Number of threads; `0` uses all available cores.
#' @param ... Not used.
#'
#' @details
#' The summary reports:
#' \describe{
#'   \item{Nodes, edges}{Rows of `x$nodes` and `x$edges`.}
#'   \item{Self-loops, multiple edges}{Edges from a node to itself, and
#'     edges repeating an earlier edge between the same (ordered, if
#'     directed) pair of nodes, as igraph's `sum(which_multiple())`.}
#'   \item{Isolates}{Nodes without an edge to another node.}
#'   \item{Density}{[graph_density()].}
#'   \item{Components}{The number of weak components ([components()]) and
#'     the size of the largest.}
#'   \item{Reciprocity}{[reciprocity()], for directed networks.}
#'   \item{Transitivity}{[transitivity()].}
#'   \item{Mean distance, diameter}{[mean_distance()] and [diameter()]
#'     with `unconnected = TRUE`, following edge directions. When the
#'     network has more than one component they are computed on the
#'     largest weak component (the first, if several are equally large),
#'     and the printout says so.}
#'   \item{Degree}{The mean, smallest, and largest [degree()] (all edges,
#'     self-loops counted twice in undirected networks).}
#'   \item{Degree centralization}{[centralization()] of degree with its
#'     defaults.}
#' }
#' Everything except the mean distance and diameter takes time roughly
#' linear in the size of the network, so summaries of networks with
#' hundreds of thousands of nodes are quick.
#'
#' @return An object of class `summary.netframe`: a list holding the values
#'   above, with a `print()` method.
#' @export
summary.netframe <- function(object, paths = NULL, weighted = NULL, threads = 0L, ...) {
  x <- object
  check_netframe(x)
  n <- nrow(x$nodes)
  ptr <- engine_graph(x)
  counts <- .Call(C_edge_counts, ptr)
  comp <- components(x, "weak")
  largest <- if (comp$no > 0L) which.max(comp$csize) else NA_integer_
  largest_size <- if (comp$no > 0L) comp$csize[[largest]] else 0L
  w <- resolve_weighting(x, weighted, NULL)
  if (is.null(paths)) paths <- largest_size <= summary_path_limit
  if (!is.logical(paths) || length(paths) != 1L || is.na(paths)) {
    stop("`paths` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }
  mean_dist <- diam <- NA_real_
  if (paths && largest_size > 0L) {
    src <- which(comp$membership == largest) - 1L
    s <- .Call(C_path_summary, ptr, src, mode_code("out"), w$weighted, w$code, TRUE,
               resolve_threads(threads))
    mean_dist <- s[[6L]]
    diam <- s[[7L]]
  }
  deg <- degree(x)
  deg <- if (n > 0L) {
    c(mean = mean(deg), min = min(deg), max = max(deg))
  } else {
    c(mean = NaN, min = NaN, max = NaN)
  }
  structure(
    list(
      nodes = n,
      edges = nrow(x$edges),
      directed = x$directed,
      weight = x$roles$weight,
      weight_type = x$roles$weight_type,
      loops = counts[[1L]],
      multiple = counts[[2L]],
      isolates = counts[[3L]],
      density = graph_density(x),
      components = comp$no,
      largest_component = largest_size,
      reciprocity = if (x$directed) reciprocity(x) else NA_real_,
      transitivity = transitivity(x, threads = threads),
      paths = paths,
      path_weighted = w$weighted,
      mean_distance = mean_dist,
      diameter = diam,
      degree = deg,
      centralization = centralization(x, "degree", threads = threads)
    ),
    class = "summary.netframe"
  )
}

# Largest component for which summary() computes distances by default.
summary_path_limit <- 5000L

#' @rdname summary.netframe
#' @param x A `summary.netframe` object.
#' @param digits Significant digits for proportions and means.
#' @export
print.summary.netframe <- function(x, digits = 3L, ...) {
  num <- function(v) format(v, digits = digits)
  int <- function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)
  weighting <- if (is.null(x$weight)) {
    "no"
  } else {
    sprintf("yes, by `%s` (%ss)", x$weight, x$weight_type)
  }
  components <- if (x$components > 1L) {
    sprintf("%s (largest has %s nodes)", int(x$components), int(x$largest_component))
  } else {
    int(x$components)
  }
  rows <- c(
    "Nodes" = int(x$nodes),
    "Edges" = int(x$edges),
    "Directed" = if (x$directed) "yes" else "no",
    "Weighted" = weighting,
    "Self-loops" = int(x$loops),
    "Multiple edges" = int(x$multiple),
    "Isolates" = int(x$isolates),
    "Density" = num(x$density),
    "Components" = components
  )
  if (x$directed) rows <- c(rows, "Reciprocity" = num(x$reciprocity))
  rows <- c(rows, "Transitivity" = num(x$transitivity))
  if (x$paths) {
    rows <- c(rows, "Mean distance" = num(x$mean_distance), "Diameter" = num(x$diameter))
  } else {
    rows <- c(rows, "Mean distance" = "not computed", "Diameter" = "not computed")
  }
  rows <- c(rows,
    "Degree (mean, min, max)" = paste(num(x$degree[["mean"]]), num(x$degree[["min"]]),
                                      num(x$degree[["max"]]), sep = ", "),
    "Degree centralization" = num(x$centralization)
  )
  cat("<netframe summary>\n")
  cat(paste0(format(names(rows)), "  ", rows, "\n"), sep = "")
  notes <- character(0)
  if (x$paths && x$components > 1L) {
    notes <- c(notes, sprintf(
      "Mean distance and diameter are for the largest component (%s of %s nodes).",
      int(x$largest_component), int(x$nodes)))
  }
  if (x$paths && x$directed) {
    notes <- c(notes, "Distances follow edge directions; unreachable pairs are left out.")
  }
  if (x$paths && isTRUE(x$path_weighted)) {
    notes <- c(notes, sprintf("Distances use the weights as %ss.", x$weight_type))
  }
  if (!x$paths) {
    notes <- c(notes, sprintf(paste0(
      "Mean distance and diameter were skipped because the largest component has %s ",
      "nodes; use summary(x, paths = TRUE) to compute them."), int(x$largest_component)))
  }
  if (length(notes)) cat(paste0("Note: ", notes, "\n"), sep = "")
  invisible(x)
}
