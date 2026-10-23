#' Mixing matrices
#'
#' Counts ties by the categories of the nodes at their two ends. The
#' categories can come from one node attribute (the classical mixing
#' matrix) or from two different ones, such as node identifiers and a group
#' attribute, which gives each node's ties to each group.
#'
#' @param x A [netframe()].
#' @param dim1 Row categories: the name of a column of `x$nodes`, or a
#'   vector with one value per node.
#' @param dim2 Column categories, given the same way. `NULL` uses `dim1`.
#' @param direction For directed networks: `"out"` counts ties from `dim1`
#'   (rows) to `dim2` (columns), `"in"` ties received (see Details), and
#'   `"all"` both. Ignored for undirected networks.
#'
#' @details
#' The layout follows the R package homophily (Knapp 2019), so that results
#' carry over:
#' \itemize{
#'   \item Undirected networks: every edge counts in both directions. When
#'     `dim1` and `dim2` are the same attribute, the matrix is folded into
#'     its upper triangle: the diagonal holds twice the number of edges
#'     within each category, entry \eqn{(i, j)} for \eqn{i < j} the number
#'     of edges between categories \eqn{i} and \eqn{j} counted from both
#'     sides (twice the number of edges), and the lower triangle is zero.
#'     The matrix sums to twice the number of edges.
#'   \item Directed networks with one attribute: rows are the categories of
#'     senders and columns those of receivers, so each arc counts once;
#'     `direction = "in"` transposes the matrix.
#'   \item Directed networks with two attributes: `"out"` counts each arc
#'     by `dim1` of its sender and `dim2` of its receiver; `"in"` gives a
#'     matrix with `dim2` categories as rows and `dim1` categories as
#'     columns, counting each arc by `dim2` of its sender and `dim1` of its
#'     receiver; `"all"` counts every arc in both directions.
#' }
#' Categories appear in order of first appearance in `x$nodes`. Ties
#' involving a node whose category is missing are not counted. Parallel
#' edges each count, and a self-loop counts like any other tie (twice in an
#' undirected network). Weights are not used. homophily also folds a
#' square matrix built from two different attributes; knappnet folds only
#' when both dimensions are the same attribute.
#'
#' @references
#' Newman, M. E. J. (2003). Mixing patterns in networks. *Physical Review
#' E*, 67(2), 026126. \doi{10.1103/PhysRevE.67.026126}
#'
#' Knapp, B. G. (2019). homophily: Measuring network homophily. R package,
#' \url{https://github.com/knapply/homophily}.
#'
#' @return A numeric matrix with dimnames.
#' @export
mixing_matrix <- function(x, dim1, dim2 = NULL, direction = c("all", "out", "in")) {
  check_netframe(x)
  direction <- match.arg(direction)
  d1 <- node_categories(x, dim1, "dim1", missing_ok = TRUE)
  same <- is.null(dim2)
  if (!same) {
    d2 <- node_categories(x, dim2, "dim2", missing_ok = TRUE)
    same <- identical(d1$values, d2$values)
  }
  if (same) d2 <- d1
  ptr <- engine_graph(x)
  res <- engine_call("mixing_matrix",
    .Call(C_mixing_matrix, ptr, d1$codes, length(d1$levels), d2$codes, length(d2$levels),
          same, mode_code(direction)))
  k1 <- length(d1$levels)
  k2 <- length(d2$levels)
  if (res[[2]]) {
    out <- matrix(res[[1]], nrow = k2, ncol = k1, byrow = TRUE,
                  dimnames = list(d2$levels, d1$levels))
  } else {
    out <- matrix(res[[1]], nrow = k1, ncol = k2, byrow = TRUE,
                  dimnames = list(d1$levels, d2$levels))
  }
  out
}

#' Assortativity coefficients
#'
#' Newman's (2002, 2003) assortativity: the tendency of ties to join nodes
#' with similar values. `assortativity()` uses numeric node values,
#' `assortativity_nominal()` categories, and `assortativity_degree()`
#' degrees.
#'
#' @inheritParams mixing_matrix
#' @param values Numeric node values: a vector with one value per node or
#'   the name of a column of `x$nodes`.
#' @param values_in Values at the receiving end of ties in directed
#'   networks, given the same way; `NULL` uses `values`.
#' @param types Node categories: a vector with one value per node or the
#'   name of a column of `x$nodes`. Missing values are not allowed.
#' @param directed Whether to respect the direction of ties in directed
#'   networks. `FALSE` treats every tie as reciprocal.
#' @param normalized Whether to return the coefficient (`TRUE`) or its
#'   numerator: the covariance of the values at the two ends of ties for
#'   `assortativity()`, and the modularity of the categories for
#'   `assortativity_nominal()`.
#' @param weighted For `assortativity_degree()`: use strengths (sums of
#'   weights) instead of degrees, as NetworkX does with `weight`.
#' @param convention How an undirected self-loop counts: `"igraph"` and
#'   `"homophily"` treat it as two arcs (each undirected edge is a
#'   reciprocal pair), `"networkx"` as one. `"homophily"` also ignores
#'   direction in `assortativity_nominal()`, as homophily's
#'   `assort_discrete()` does.
#'
#' @details
#' Every tie is read as arcs: an arc of a directed network once (or, with
#' `directed = FALSE`, in both directions), an undirected edge in both
#' directions. Parallel edges each count.
#'
#' For categories, with \eqn{e_{ij}} the share of arcs from category
#' \eqn{i} to category \eqn{j}, \eqn{a_i = \sum_j e_{ij}} and
#' \eqn{b_j = \sum_i e_{ij}},
#' \deqn{r = \frac{\sum_i e_{ii} - \sum_i a_i b_i}{1 - \sum_i a_i b_i}}{r =
#' (sum_i e_ii - sum_i a_i b_i) / (1 - sum_i a_i b_i)}
#' (Newman 2003, eq. 2). For numeric values, \eqn{r} is the Pearson
#' correlation of the values at the source and target of every arc;
#' directed networks may use different values at each end (`values_in`).
#' Degree assortativity correlates the out-degree of the source with the
#' in-degree of the target, or total degree when direction is ignored;
#' using degree minus one (Newman's "excess degree") gives the same
#' correlation. Degrees count every edge end, so an undirected self-loop
#' adds 2.
#'
#' Results are `NaN` when undefined: without ties, or when the values at
#' one end do not vary. They match igraph's `assortativity()`,
#' `assortativity_nominal()`, and `assortativity_degree()`, NetworkX's
#' `numeric_assortativity_coefficient()`,
#' `attribute_assortativity_coefficient()`, and
#' `degree_assortativity_coefficient()` (which return 0 instead of `NaN`
#' without ties), and homophily's `assort_scalar()`, `assort_discrete()`
#' (which ignores direction, so it equals `directed = FALSE` here), and
#' `assort_degree()`.
#'
#' @references
#' Newman, M. E. J. (2002). Assortative mixing in networks. *Physical
#' Review Letters*, 89(20), 208701. \doi{10.1103/PhysRevLett.89.208701}
#'
#' Newman, M. E. J. (2003). Mixing patterns in networks. *Physical Review
#' E*, 67(2), 026126. \doi{10.1103/PhysRevE.67.026126}
#'
#' @return A number.
#' @export
assortativity <- function(x, values, values_in = NULL, directed = TRUE, normalized = TRUE,
                          convention = c("igraph", "networkx", "homophily")) {
  check_netframe(x)
  convention <- match.arg(convention)
  v <- node_values(x, values, "values")
  vin <- if (is.null(values_in)) NULL else node_values(x, values_in, "values_in")
  directed <- check_flag(directed, "directed")
  normalized <- check_flag(normalized, "normalized")
  if (!is.null(vin) && !x$directed) {
    warning("`values_in` is ignored for undirected networks.", call. = FALSE)
    vin <- NULL
  }
  ptr <- engine_graph(x)
  engine_call("assortativity",
    .Call(C_assortativity_values, ptr, v, vin, directed, assort_loops(convention),
          normalized))
}

#' @rdname assortativity
#' @export
assortativity_nominal <- function(x, types, directed = TRUE, normalized = TRUE,
                                  convention = c("igraph", "networkx", "homophily")) {
  check_netframe(x)
  convention <- match.arg(convention)
  t <- node_categories(x, types, "types", missing_ok = FALSE)
  directed <- check_flag(directed, "directed")
  normalized <- check_flag(normalized, "normalized")
  if (convention == "homophily") directed <- FALSE
  ptr <- engine_graph(x)
  engine_call("assortativity_nominal",
    .Call(C_assortativity_nominal, ptr, t$codes, length(t$levels), directed,
          assort_loops(convention), normalized))
}

#' @rdname assortativity
#' @export
assortativity_degree <- function(x, directed = TRUE, weighted = FALSE, normalized = TRUE,
                                 convention = c("igraph", "networkx", "homophily")) {
  check_netframe(x)
  convention <- match.arg(convention)
  directed <- check_flag(directed, "directed")
  weighted <- check_flag(weighted, "weighted")
  normalized <- check_flag(normalized, "normalized")
  if (weighted) require_weights(x)
  ptr <- engine_graph(x)
  engine_call("assortativity_degree",
    .Call(C_assortativity_degree, ptr, weighted, directed, assort_loops(convention),
          normalized))
}

assort_loops <- function(convention) if (convention == "networkx") 1L else 0L

#' E-I index
#'
#' Krackhardt and Stern's (1988) E-I index compares the ties that cross
#' group boundaries (external, E) with those inside groups (internal, I):
#' \deqn{\frac{E - I}{E + I}}{(E - I) / (E + I)}
#' It ranges from -1 (every tie inside a group) to 1 (every tie between
#' groups).
#'
#' @inheritParams mixing_matrix
#' @param groups Group of each node: the name of a column of `x$nodes`, or
#'   a vector with one value per node. Missing values are not allowed.
#' @param scope `"global"` for the whole network, `"group"` for each group,
#'   or `"node"` for each node.
#' @param direction For groups and nodes of directed networks: count ties
#'   sent (`"out"`), received (`"in"`), or both (`"all"`).
#' @param convention How group-level ties are counted: `"krackhardt"` (see
#'   Details) or `"homophily"`, which reproduces the R package homophily.
#'
#' @details
#' Network: E is the number of ties between nodes of different groups and I
#' the number within groups. Group (Krackhardt and Stern 1988): I is the
#' number of ties among the group's members and E the number of ties
#' between a member and a non-member (with `direction = "out"` or `"in"`,
#' only those sent or received by members). Node: the node's own ties to
#' members of its group (I) and to others (E); in directed networks ties
#' sent and received both count with `direction = "all"`, so a mutual dyad
#' counts twice.
#'
#' Ties are binary: weights are ignored and each parallel edge counts as a
#' separate tie (a warning is given when there are parallel edges, as
#' homophily does). A self-loop is an internal tie. An index is `NaN` where
#' there are no ties.
#'
#' `convention = "homophily"` matches homophily's `ei_index()`. Network
#' and node indices agree with the default except that homophily counts an
#' undirected self-loop twice at its node. For groups, homophily takes I
#' from the diagonal of its folded mixing matrix (see [mixing_matrix()])
#' and E from the group's row of the upper triangle, so a tie between two
#' groups counts as external only for the group that appears first in
#' `x$nodes`. UCINET also reports the index for groups and nodes; its
#' documentation could not be consulted, so agreement with UCINET is not
#' claimed for those levels.
#'
#' @references
#' Krackhardt, D., & Stern, R. N. (1988). Informal networks and
#' organizational crises: An experimental simulation. *Social Psychology
#' Quarterly*, 51(2), 123-140. \doi{10.2307/2786835}
#'
#' Knapp, B. G. (2019). homophily: Measuring network homophily. R package,
#' \url{https://github.com/knapply/homophily}.
#'
#' @return For `"global"`, a number. For `"group"`, a vector named by group
#'   in order of first appearance. For `"node"`, a vector aligned with the
#'   rows of `x$nodes`. External and internal tie counts are attached as
#'   attributes `"external"` and `"internal"`.
#' @export
ei_index <- function(x, groups, scope = c("global", "group", "node"),
                     direction = c("all", "out", "in"),
                     convention = c("krackhardt", "homophily")) {
  check_netframe(x)
  scope <- match.arg(scope)
  direction <- match.arg(direction)
  convention <- match.arg(convention)
  g <- node_categories(x, groups, "groups", missing_ok = FALSE)
  ptr <- engine_graph(x)
  warn_parallel(x)
  res <- engine_call("ei_index",
    .Call(C_ei_index, ptr, g$codes, length(g$levels),
          switch(scope, global = 0L, group = 1L, node = 2L), mode_code(direction),
          if (convention == "krackhardt") 0L else 1L))
  out <- res[[1]]
  if (scope == "group") names(out) <- g$levels
  attr(out, "external") <- res[[2]]
  attr(out, "internal") <- res[[3]]
  out
}

#' Within-group tie shares
#'
#' For each group, the share of its members' ties that stay inside the
#' group, against the share expected if ties ignored groups: the baseline
#' homophily of McPherson, Smith-Lovin, and Cook (2001).
#'
#' @inheritParams ei_index
#' @param direction For directed networks: ties sent (`"out"`), received
#'   (`"in"`), or both (`"all"`) by members.
#'
#' @details
#' Ties are counted from the members' side: a tie between two members
#' counts once for each, so `internal` is twice the number of internal ties
#' in an undirected network (and with `direction = "all"` in a directed
#' one). `observed` is internal / (internal + external), `NaN` for a group
#' without ties. `expected` is \eqn{(n_g - 1)/(n - 1)}, the share of a
#' member's possible contacts who belong to its group, which is the
#' expected share when every tie is equally likely; `observed` above
#' `expected` indicates ties concentrated within the group. Self-loops are
#' ignored, parallel edges each count, and weights are not used.
#'
#' @references
#' McPherson, M., Smith-Lovin, L., & Cook, J. M. (2001). Birds of a
#' feather: Homophily in social networks. *Annual Review of Sociology*, 27,
#' 415-444. \doi{10.1146/annurev.soc.27.1.415}
#'
#' @return A data frame with one row per group (in order of first
#'   appearance): `group`, `size`, `internal`, `external`, `observed`, and
#'   `expected`.
#' @export
homophily <- function(x, groups, direction = c("all", "out", "in")) {
  check_netframe(x)
  direction <- match.arg(direction)
  g <- node_categories(x, groups, "groups", missing_ok = FALSE)
  k <- length(g$levels)
  ptr <- engine_graph(x)
  v <- engine_call("homophily", .Call(C_group_shares, ptr, g$codes, k, mode_code(direction)))
  m <- matrix(v, nrow = k, ncol = 5L, byrow = TRUE)
  data.frame(group = g$levels, size = as.integer(m[, 1]), internal = m[, 2],
             external = m[, 3], observed = m[, 4], expected = m[, 5],
             stringsAsFactors = FALSE)
}

# Node categories: a column of x$nodes or a vector aligned with the nodes.
# Returns the values, 0-based integer codes in order of first appearance
# (NA where missing), and the category labels.
node_categories <- function(x, values, what, missing_ok) {
  if (is.character(values) && length(values) == 1L && nrow(x$nodes) != 1L) {
    if (!values %in% names(x$nodes)) {
      stop(sprintf("`%s` names no column of `x$nodes`.", what), call. = FALSE)
    }
    values <- x$nodes[[values]]
  } else if (is.character(values) && length(values) == 1L && values %in% names(x$nodes)) {
    values <- x$nodes[[values]]
  }
  if (!is.atomic(values) || length(values) != nrow(x$nodes)) {
    stop(sprintf("`%s` must have one value per node (%d), or name a column of `x$nodes`.",
                 what, nrow(x$nodes)), call. = FALSE)
  }
  if (!missing_ok && anyNA(values)) {
    stop(sprintf("`%s` must not contain missing values.", what), call. = FALSE)
  }
  plain <- if (is.factor(values)) as.character(values) else values
  levels <- unique(plain[!is.na(plain)])
  list(values = values, codes = match(plain, levels) - 1L, levels = as.character(levels))
}

warn_parallel <- function(x) {
  src <- x$index$src
  tgt <- x$index$tgt
  key <- if (x$directed) paste(src, tgt) else paste(pmin(src, tgt), pmax(src, tgt))
  if (anyDuplicated(key)) {
    warning("The network has parallel edges; each counts as a separate tie.",
            call. = FALSE)
  }
}
