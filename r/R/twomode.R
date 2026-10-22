#' Two-mode networks
#'
#' A two-mode (bipartite, affiliation) network has two sets of nodes, such
#' as people and the events they attend, or people and the organizations
#' they belong to, with ties only between the sets. `twomode()` creates
#' one from data frames, `as_twomode()` converts a netframe, an igraph or
#' network object, an incidence matrix, or a network read from a file, and
#' `is_twomode()` tests for one. `incidence_matrix()` returns the
#' rectangular matrix with one row per node of the first mode and one
#' column per node of the second.
#'
#' A two-mode network is a [netframe()] whose node table has a logical
#' column, recorded as the node type role (`x$roles$type`): `FALSE` for the
#' first mode and `TRUE` for the second, igraph's convention. Every edge
#' must join the two modes; edges within a mode are an error.
#'
#' The type may be given as a logical column, or as a column with exactly
#' two distinct values (for example `"person"` and `"event"`, the two node
#' classes of a DyNetML file, or the 1 and 2 of a bipartite network from
#' the network package). `first` names the value of the first mode; by
#' default it is the first factor level, or the smallest value (text is
#' compared byte by byte, so the choice does not depend on the locale). A logical
#' column is then added: it is called `type` (so that [to_igraph()] gives
#' igraph's bipartite functions the `type` attribute they look for) unless
#' the node table already has a column with that name, in which case it is
#' called `is_` followed by the second mode's value. The original column
#' is kept, and the two values are recorded in `x$roles$type_labels`.
#'
#' Without a node table, `twomode()` takes the first mode from the `from`
#' column of the edges and the second from the `to` column, so a list of
#' (person, event) pairs is all it needs. A value appearing in both
#' columns is an error.
#'
#' `as_twomode()` finds the type column itself when `type` is `NULL`: a
#' logical `type` column, else a `class` column with exactly two values (as
#' [read_network()] gives for a DyNetML file with two node classes), else a
#' two-valued `type` or `mode` column. A matrix (or sparse `Matrix`) is read
#' as an incidence matrix: rows are the first mode and columns the second,
#' with identifiers from the row and column names (`r1`, `r2`, ... and
#' `c1`, `c2`, ... if there are none); values other than 0 and 1 become
#' edge weights unless `weighted = FALSE`. A data frame is passed to
#' `twomode()`.
#'
#' Two-mode measures ([project()], [twomode_degree()] and its relatives)
#' read the network as an affiliation matrix and ignore the direction of
#' edges, so `twomode()` creates undirected networks by default.
#'
#' @inheritParams netframe
#' @param type For `twomode()` and `as_twomode()`, the node column holding
#'   the type, by name. When `twomode()` builds the node table itself, this
#'   is the name of the logical column it creates.
#' @param first The value of `type` that marks the first mode; see Details.
#' @param x For `as_twomode()`, the object to convert; for the other
#'   functions, a two-mode [netframe()].
#' @param weighted For matrices, whether values are edge weights. `NULL`
#'   treats any value other than 0 and 1 as a weight.
#' @param attr For `incidence_matrix()`, an edge column whose values fill
#'   the matrix; `NULL` counts edges.
#' @param sparse For `incidence_matrix()`, whether to return a sparse
#'   `Matrix`.
#'
#' @references
#' Breiger, R. L. (1974). The duality of persons and groups. *Social
#' Forces*, 53(2), 181-190. \doi{10.2307/2576011}
#'
#' Borgatti, S. P., & Everett, M. G. (1997). Network analysis of 2-mode
#' data. *Social Networks*, 19(3), 243-269.
#' \doi{10.1016/S0378-8733(96)00301-2}
#'
#' Borgatti, S. P., & Halgin, D. S. (2011). Analyzing affiliation networks.
#' In J. Scott & P. J. Carrington (Eds.), *The SAGE Handbook of Social
#' Network Analysis* (pp. 417-433). SAGE. \doi{10.4135/9781446294413.n28}
#'
#' @return `twomode()` and `as_twomode()` return a [netframe()] with a node
#'   type role. `is_twomode()` returns `TRUE` or `FALSE`.
#'   `incidence_matrix()` returns a matrix with first-mode nodes as rows and
#'   second-mode nodes as columns; parallel edges add up.
#'
#' @examples
#' attended <- data.frame(
#'   person = c("ann", "ann", "bob", "bob", "cat", "dan"),
#'   event  = c("gala", "fair", "gala", "show", "show", "fair")
#' )
#' g <- twomode(attended)
#' g
#' g$nodes
#' incidence_matrix(g)
#' is_twomode(g)
#' @export
twomode <- function(edges, nodes = NULL, type = "type", first = NULL, directed = FALSE,
                    from = 1L, to = 2L, id = 1L, weight = NULL,
                    weight_type = c("distance", "strength"), graph = list()) {
  weight_type <- match.arg(weight_type)
  if (!is.data.frame(edges)) {
    stop("`edges` must be a data frame.", call. = FALSE)
  }
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("`type` must be a single column name.", call. = FALSE)
  }
  if (is.null(nodes)) {
    if (type == "name") {
      stop("Without `nodes`, identifiers go in a column called `name`; choose another `type`.",
           call. = FALSE)
    }
    from_col <- resolve_column(edges, from, "from", "edges")
    to_col <- resolve_column(edges, to, "to", "edges")
    a <- unique(as.vector(edges[[from_col]]))
    b <- unique(as.vector(edges[[to_col]]))
    both <- intersect(a, b)
    if (length(both)) {
      stop(sprintf(paste0(
        "%d identifier(s) appear in both `%s` and `%s`, for example %s; in a two-mode ",
        "network every node belongs to one mode. Supply `nodes` with a type column."),
        length(both), from_col, to_col, paste(format(utils::head(both, 3L)), collapse = ", ")),
        call. = FALSE)
    }
    nodes <- data.frame(name = c(a, b), stringsAsFactors = FALSE)
    nodes[[type]] <- rep(c(FALSE, TRUE), c(length(a), length(b)))
    id <- "name"
    if (!is.null(first)) {
      stop("`first` applies only when `nodes` is given; without it the `from` column ",
           "is the first mode.", call. = FALSE)
    }
  }
  x <- netframe(edges, nodes, directed = directed, from = from, to = to, id = id,
                weight = weight, weight_type = weight_type, graph = graph)
  set_twomode(x, type, first)
}

#' @rdname twomode
#' @export
as_twomode <- function(x, type = NULL, first = NULL, weighted = NULL) {
  if (is.matrix(x) || (isS4(x) && methods::is(x, "Matrix"))) {
    return(incidence_to_twomode(x, weighted))
  }
  if (is.data.frame(x) && !inherits(x, "netframe")) {
    return(twomode(x, type = if (is.null(type)) "type" else type, first = first))
  }
  x <- as_netframe(x)
  if (is.null(type)) {
    if (is_twomode(x)) return(x)
    type <- find_type_column(x$nodes)
  }
  set_twomode(x, type, first)
}

#' @rdname twomode
#' @export
is_twomode <- function(x) {
  if (!inherits(x, "netframe") || is.null(x$roles$type)) return(FALSE)
  v <- x$nodes[[x$roles$type]]
  is.logical(v) && !anyNA(v)
}

#' @rdname twomode
#' @export
incidence_matrix <- function(x, attr = NULL, sparse = FALSE) {
  check_twomode(x)
  second <- x$nodes[[x$roles$type]]
  m <- as.matrix(x, attr = attr, sparse = sparse)
  rows <- which(!second)
  cols <- which(second)
  # A directed network may store a tie either way; an undirected one is
  # already symmetric.
  if (x$directed) {
    m[rows, cols, drop = FALSE] + t(m[cols, rows, drop = FALSE])
  } else {
    m[rows, cols, drop = FALSE]
  }
}

#' One-mode projections of a two-mode network
#'
#' Joins two nodes of the same mode when they share a partner in the other
#' mode: two people who attended the same event, or two events with a
#' common participant (Breiger 1974). `which = "both"` returns both
#' projections, the dual-projection view of Everett and Borgatti (2013).
#'
#' Tie weights, chosen by `method`, are
#' \describe{
#'   \item{`"count"`}{the number of shared partners, the off-diagonal of
#'     \eqn{BB^T} for the binary incidence matrix \eqn{B}. This is igraph's
#'     `bipartite_projection(multiplicity = TRUE)` and NetworkX's
#'     `weighted_projected_graph()`.}
#'   \item{`"newman"`}{Newman's (2001) collaboration weight
#'     \eqn{\sum_e 1/(d_e - 1)} over shared partners \eqn{e}, where
#'     \eqn{d_e} is the number of partners of \eqn{e}, so that a large event
#'     ties each pair of its participants more weakly. This is NetworkX's
#'     `collaboration_weighted_projected_graph()`.}
#'   \item{`"product"`}{\eqn{\sum_e w_{ie} w_{je}} over shared partners,
#'     the off-diagonal of \eqn{WW^T} for the valued incidence matrix
#'     \eqn{W} (see [incidence_matrix()]). It needs edge weights, which are
#'     used as given whatever the network's `weight_type`.}
#' }
#' For `"count"` and `"newman"` parallel edges count once (shared partners
#' are distinct partners); for `"product"` their weights add up, as in
#' [incidence_matrix()]. Direction is ignored.
#'
#' With `loops = TRUE` the diagonal is kept as self-loops: each node's
#' number of partners for `"count"` and its sum of squared weights for
#' `"product"`, so the projection holds all of \eqn{BB^T}, as the dual
#' projection approach uses. Newman's weights have no diagonal.
#'
#' The projection is computed with a sparse accumulator for each node, in
#' parallel, in time proportional to the number of two-step paths.
#'
#' @param x A two-mode [netframe()] (see [twomode()]).
#' @param which `"first"`, `"second"`, or `"both"` modes.
#' @param weighted Whether the result has a `weight` column. Weights are
#'   tie strengths (`weight_type = "strength"`), since sharing more
#'   partners brings two nodes closer.
#' @param method How ties are weighted; see Details.
#' @param loops Whether to keep the diagonal as self-loops.
#' @param threads Number of threads; `0` uses all available cores.
#'
#' @references
#' Breiger, R. L. (1974). The duality of persons and groups. *Social
#' Forces*, 53(2), 181-190. \doi{10.2307/2576011}
#'
#' Newman, M. E. J. (2001). Scientific collaboration networks. II.
#' Shortest paths, weighted networks, and centrality. *Physical Review E*,
#' 64(1), 016132. \doi{10.1103/PhysRevE.64.016132}
#'
#' Everett, M. G., & Borgatti, S. P. (2013). The dual-projection approach
#' for two-mode networks. *Social Networks*, 35(2), 204-210.
#' \doi{10.1016/j.socnet.2012.05.004}
#'
#' @return An undirected [netframe()] whose nodes are the rows of `x$nodes`
#'   in the chosen mode (all columns kept) and whose edges join each pair of
#'   nodes with a shared partner once, ordered by node position. For
#'   `which = "both"`, a list with elements `first` and `second`.
#'
#' @examples
#' attended <- data.frame(
#'   person = c("ann", "ann", "bob", "bob", "cat", "dan"),
#'   event  = c("gala", "fair", "gala", "show", "show", "fair")
#' )
#' g <- twomode(attended)
#' project(g)$edges
#' project(g, "second", method = "newman")$edges
#' @export
project <- function(x, which = c("first", "second", "both"), weighted = TRUE,
                    method = c("count", "newman", "product"), loops = FALSE, threads = 0L) {
  check_twomode(x)
  which <- match.arg(which)
  method <- match.arg(method)
  if (which == "both") {
    return(list(
      first = project(x, "first", weighted, method, loops, threads),
      second = project(x, "second", weighted, method, loops, threads)
    ))
  }
  if (method == "product") require_weights(x)
  if (method == "newman" && isTRUE(loops)) {
    stop("Newman's collaboration weights have no diagonal; use `loops = FALSE`.",
         call. = FALSE)
  }
  second <- which == "second"
  code <- switch(method, count = 0L, newman = 1L, product = 2L)
  res <- .Call(C_project, engine_graph(x), twomode_types(x), second, code, isTRUE(loops),
               as.integer(threads))
  in_mode <- x$nodes[[x$roles$type]] == second
  nodes <- x$nodes[in_mode, , drop = FALSE]
  rownames(nodes) <- NULL
  ids <- x$nodes[[x$roles$id]]
  edges <- data.frame(from = ids[res[[1]]], to = ids[res[[2]]], stringsAsFactors = FALSE)
  if (isTRUE(weighted)) edges$weight <- res[[3]]
  netframe(edges, nodes, directed = FALSE, from = "from", to = "to", id = x$roles$id,
           weight = if (isTRUE(weighted)) "weight", weight_type = "strength",
           graph = x$graph)
}

#' Centrality and density of two-mode networks
#'
#' Degree, closeness, betweenness, and eigenvector centrality of the nodes
#' of a two-mode network, normalized as Borgatti and Everett (1997) and
#' Borgatti and Halgin (2011) describe, and two-mode density.
#'
#' The scores are those of the bipartite graph; only the normalization
#' differs from the one-mode case, because a node can be tied only to the
#' other mode. For a node whose own mode has \eqn{n_o} nodes, with
#' \eqn{n_i} nodes in the other mode:
#' \describe{
#'   \item{degree}{\eqn{d / n_i}.}
#'   \item{closeness}{\eqn{(n_i + 2(n_o - 1)) / S}, where \eqn{S} is the sum
#'     of the node's distances to the others; the numerator is the smallest
#'     possible sum (distance 1 to every node of the other mode, 2 to every
#'     other node of its own).}
#'   \item{betweenness}{\eqn{b / b_{max}} with
#'     \eqn{b_{max} = \frac{1}{2}[n_i^2 (s+1)^2 + n_i (s+1)(2t - s - 1) - t(2s - t + 3)]},
#'     \eqn{s = (n_o - 1) \mathrm{div}\ n_i} and
#'     \eqn{t = (n_o - 1) \bmod n_i}, the largest betweenness a node can
#'     have given the sizes of the two modes (Borgatti and Everett 1997).
#'     The raw betweenness counts each unordered pair once.}
#'   \item{eigenvector}{the principal eigenvector of the bipartite graph,
#'     which is the pair of leading singular vectors of the incidence
#'     matrix (Borgatti and Everett 1997; Borgatti and Halgin 2011). The
#'     two modes' halves of the unit eigenvector have the same length, so
#'     Borgatti and Everett's normalization (dividing by
#'     \eqn{\sqrt{1/2}}) gives each mode Euclidean length 1, which is what
#'     is returned. Their further two-mode normalization of eigenvector
#'     scores has no closed form and is not provided.}
#'   \item{density}{the number of ties divided by \eqn{n_1 n_2}, the number
#'     of possible ties between the modes.}
#' }
#' With `normalized = FALSE` degree, betweenness, and closeness are the raw
#' scores \eqn{d}, \eqn{b}, and \eqn{1/S}.
#'
#' Ties are binary and undirected: weights are ignored, parallel edges
#' count once, and direction does not matter. When the network is not
#' connected, \eqn{S} sums the distances to the nodes a node reaches. With
#' `unreachable = "wasserman_faust"` (the default) closeness is then
#' multiplied by \eqn{r/(n-1)}, the share of the other nodes it reaches,
#' and a node that reaches no other scores 0; with `"ignore"` it is not
#' rescaled and a node that reaches no other gets `NaN`. Eigenvector scores
#' are then unique only within the component with the largest eigenvalue
#' (a warning is given); every score is `NaN` in a network without ties.
#' Normalizations that divide by zero (an empty mode, or \eqn{b_{max} = 0}
#' when the other mode has one node) give `NaN`.
#'
#' Degree, closeness, and betweenness match NetworkX's
#' `bipartite.degree_centrality()`, `bipartite.closeness_centrality()`
#' (`normalized = TRUE` for the Wasserman-Faust rule, `FALSE` for
#' `"ignore"`), and `bipartite.betweenness_centrality()`, and density
#' matches `bipartite.density()`. Raw scores match igraph's `degree()`,
#' `closeness()`, and `betweenness()` on the simplified network, and
#' eigenvector scores match igraph's `eigen_centrality()` rescaled within
#' each mode. UCINET reports the same normalized scores as percentages.
#'
#' @param x A two-mode [netframe()] (see [twomode()]).
#' @param normalized Whether to apply the two-mode normalization.
#' @param unreachable How closeness treats nodes that cannot be reached:
#'   `"wasserman_faust"` or `"ignore"`; see Details.
#' @param tol Convergence tolerance for eigenvector scores: iteration stops
#'   when no score (scaled so the largest is 1) changes by more than `tol`.
#' @param max_iter Iteration cap. Not converging within it is an error.
#' @param threads Number of threads; `0` uses all available cores.
#'
#' @references
#' Borgatti, S. P., & Everett, M. G. (1997). Network analysis of 2-mode
#' data. *Social Networks*, 19(3), 243-269.
#' \doi{10.1016/S0378-8733(96)00301-2}
#'
#' Borgatti, S. P., & Halgin, D. S. (2011). Analyzing affiliation networks.
#' In J. Scott & P. J. Carrington (Eds.), *The SAGE Handbook of Social
#' Network Analysis* (pp. 417-433). SAGE. \doi{10.4135/9781446294413.n28}
#'
#' Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and
#' Applications*. Cambridge University Press.
#' \doi{10.1017/CBO9780511815478}
#'
#' @return `twomode_density()` returns a number; the others return a
#'   numeric vector aligned with the rows of `x$nodes`.
#'
#' @examples
#' attended <- data.frame(
#'   person = c("ann", "ann", "bob", "bob", "cat", "dan"),
#'   event  = c("gala", "fair", "gala", "show", "show", "fair")
#' )
#' g <- twomode(attended)
#' data.frame(g$nodes,
#'            degree = twomode_degree(g),
#'            closeness = twomode_closeness(g),
#'            betweenness = twomode_betweenness(g))
#' twomode_density(g)
#' @export
twomode_degree <- function(x, normalized = TRUE) {
  twomode_centrality(x, 0L, normalized, 0L, 1L)
}

#' @rdname twomode_degree
#' @export
twomode_closeness <- function(x, normalized = TRUE,
                              unreachable = c("wasserman_faust", "ignore"), threads = 0L) {
  unreachable <- match.arg(unreachable)
  if (unreachable == "wasserman_faust" && !isTRUE(normalized)) {
    stop("The Wasserman-Faust rule applies to normalized closeness; use ",
         "`unreachable = \"ignore\"` with `normalized = FALSE`.", call. = FALSE)
  }
  twomode_centrality(x, 1L, normalized, if (unreachable == "ignore") 0L else 1L, threads)
}

#' @rdname twomode_degree
#' @export
twomode_betweenness <- function(x, normalized = TRUE, threads = 0L) {
  twomode_centrality(x, 2L, normalized, 0L, threads)
}

#' @rdname twomode_degree
#' @export
twomode_eigenvector <- function(x, tol = 1e-12, max_iter = 10000L, threads = 0L) {
  check_twomode(x)
  if (!is.numeric(tol) || length(tol) != 1L || !isTRUE(tol > 0)) {
    stop("`tol` must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || !isTRUE(max_iter >= 1)) {
    stop("`max_iter` must be a positive whole number.", call. = FALSE)
  }
  res <- .Call(C_twomode_eigenvector, engine_graph(x), twomode_types(x), as.double(tol),
               as.double(max_iter), as.integer(threads))
  if (!res[[3]] && nrow(x$edges) > 0L) {
    warning("The network is not connected, so eigenvector scores are zero outside the ",
            "dominant component(s) and may not be unique.", call. = FALSE)
  }
  res[[1]]
}

#' @rdname twomode_degree
#' @export
twomode_density <- function(x) {
  check_twomode(x)
  .Call(C_twomode_density, engine_graph(x), twomode_types(x))
}

twomode_centrality <- function(x, measure, normalized, unreachable, threads) {
  check_twomode(x)
  .Call(C_twomode_centrality, engine_graph(x), twomode_types(x), measure,
        isTRUE(normalized), unreachable, as.integer(threads))
}

check_twomode <- function(x) {
  check_netframe(x)
  if (is.null(x$roles$type)) {
    stop("`x` is not a two-mode network; create one with twomode() or as_twomode().",
         call. = FALSE)
  }
  v <- x$nodes[[x$roles$type]]
  if (!is.logical(v) || anyNA(v)) {
    stop(sprintf(paste0("The node type column `%s` must be logical without missing values; ",
                        "recreate the network with as_twomode()."), x$roles$type),
         call. = FALSE)
  }
}

# Node types as the engine expects them: 0 for the first mode, 1 for the
# second.
twomode_types <- function(x) as.integer(x$nodes[[x$roles$type]])

# Records `type` (a node column) as the two-mode node type, converting a
# two-valued column to a new logical column, and checks every edge.
set_twomode <- function(x, type, first) {
  if (!is.character(type) || length(type) != 1L || !type %in% names(x$nodes)) {
    stop(sprintf("Node column `%s` not found.", format(type)), call. = FALSE)
  }
  v <- x$nodes[[type]]
  if (anyNA(v)) {
    stop(sprintf("Node type column `%s` has missing values.", type), call. = FALSE)
  }
  if (is.logical(v)) {
    if (!is.null(first)) {
      stop("`first` applies only to a type column that is not logical; a logical ",
           "column always has FALSE for the first mode.", call. = FALSE)
    }
    col <- type
    labels <- c("FALSE", "TRUE")
    is_second <- v
  } else {
    # Radix sorting orders text by bytes, the same in every locale.
    values <- if (is.factor(v)) levels(droplevels(v)) else sort(unique(v), method = "radix")
    if (length(values) != 2L) {
      stop(sprintf(paste0("Node type column `%s` must be logical or have exactly two ",
                          "distinct values; it has %d."), type, length(values)),
           call. = FALSE)
    }
    key <- if (is.factor(v)) as.character(v) else v
    if (!is.null(first)) {
      if (length(first) != 1L || !as.character(first) %in% as.character(values)) {
        stop(sprintf("`first` must be one of the values of `%s`: %s.", type,
                     paste(format(values), collapse = ", ")), call. = FALSE)
      }
      first_value <- values[as.character(values) == as.character(first)]
    } else {
      first_value <- values[[1L]]
    }
    second_value <- values[as.character(values) != as.character(first_value)]
    is_second <- as.character(key) != as.character(first_value)
    labels <- c(as.character(first_value), as.character(second_value))
    col <- if ("type" %in% names(x$nodes)) {
      make.unique(c(names(x$nodes), make.names(paste0("is_", labels[[2L]]))))[ncol(x$nodes) + 1L]
    } else {
      "type"
    }
    x$nodes[[col]] <- is_second
  }
  bad <- which(is_second[x$index$src + 1L] == is_second[x$index$tgt + 1L])
  if (length(bad)) {
    k <- bad[[1L]]
    ids <- x$nodes[[x$roles$id]]
    a <- x$index$src[[k]] + 1L
    stop(sprintf(paste0(
      "%d edge(s) join two nodes of the same mode, for example edge %d between `%s` and ",
      "`%s` (both %s); in a two-mode network every edge joins the two modes."),
      length(bad), k, format(ids[[a]]), format(ids[[x$index$tgt[[k]] + 1L]]),
      if (is_second[[a]]) labels[[2L]] else labels[[1L]]), call. = FALSE)
  }
  x$roles$type <- col
  x$roles$type_labels <- labels
  x
}

# A logical `type` column, else a two-valued `class` (DyNetML node
# classes), `type`, or `mode` (network's bipartite networks) column.
find_type_column <- function(nodes) {
  two <- function(col) {
    col %in% names(nodes) && !anyNA(nodes[[col]]) &&
      length(unique(as.character(nodes[[col]]))) == 2L
  }
  if ("type" %in% names(nodes) && is.logical(nodes$type)) return("type")
  for (col in c("class", "type", "mode")) if (two(col)) return(col)
  stop(paste0("Cannot tell the two modes apart. Name the node column that holds them with ",
              "`type =` (a logical column, or one with exactly two distinct values)."),
       call. = FALSE)
}

# An incidence matrix as a two-mode network: rows are the first mode.
incidence_to_twomode <- function(m, weighted) {
  rows <- rownames(m)
  cols <- colnames(m)
  if (is.null(rows)) rows <- paste0("r", seq_len(nrow(m)))
  if (is.null(cols)) cols <- paste0("c", seq_len(ncol(m)))
  if (anyDuplicated(c(rows, cols))) {
    stop("Row and column names of an incidence matrix must be distinct identifiers.",
         call. = FALSE)
  }
  if (is.matrix(m)) {
    if (!is.numeric(m) && !is.logical(m)) {
      stop("An incidence matrix must be numeric or logical.", call. = FALSE)
    }
    m[is.na(m)] <- 0
    idx <- which(m != 0, arr.ind = TRUE)
    i <- idx[, 1]
    j <- idx[, 2]
    v <- as.numeric(m[idx])
  } else {
    require_pkg("Matrix")
    s <- methods::as(m, "TsparseMatrix")
    i <- s@i + 1L
    j <- s@j + 1L
    v <- if (methods::.hasSlot(s, "x")) as.numeric(s@x) else rep(1, length(i))
    keep <- v != 0
    i <- i[keep]
    j <- j[keep]
    v <- v[keep]
  }
  o <- order(i, j)
  i <- i[o]
  j <- j[o]
  v <- v[o]
  weighted <- if (is.null(weighted)) any(v != 1) else isTRUE(weighted)
  edges <- data.frame(from = rows[i], to = cols[j], stringsAsFactors = FALSE)
  if (weighted) edges$weight <- v
  nodes <- data.frame(name = c(rows, cols), type = rep(c(FALSE, TRUE), c(length(rows),
                                                                         length(cols))),
                      stringsAsFactors = FALSE)
  twomode(edges, nodes, type = "type", from = "from", to = "to", id = "name",
          weight = if (weighted) "weight")
}
