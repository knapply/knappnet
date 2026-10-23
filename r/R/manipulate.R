#' Combine the relations of a multiplex network
#'
#' A multiplex network stores several relations in one edge table, told
#' apart by a column such as `relation` (as [example_network("noordin")]
#' and files read from DyNetML have). `combine_relations()` builds a
#' one-relation network from them:
#' \describe{
#'   \item{`"union"`}{A pair is tied if it is tied in any of the relations
#'     (UCINET: Transform > Matrix Operations > Between Datasets >
#'     Statistical Summaries, Maximum; Boolean OR).}
#'   \item{`"intersection"`}{A pair is tied if it is tied in every one of the
#'     relations (Statistical Summaries, Minimum; Boolean AND).}
#'   \item{`"sum"`}{A pair is tied with weight equal to the number of
#'     relations in which it is tied, or the sum of its weights over the
#'     relations when the network is weighted (Statistical Summaries, Sum).}
#'   \item{`"difference"`}{A pair is tied if it is tied in the first of two
#'     relations and not in the second (Transform > Matrix Operations >
#'     Between Datasets > Boolean Combinations).}
#'   \item{`"select"`}{The edges of one relation, with all their columns
#'     (Data > Unpack).}
#' }
#'
#' @details
#' Within one relation parallel edges count once (a pair is tied or not),
#' except that `"sum"` adds the weights of parallel edges when the network
#' is weighted. Self-loops are kept as ties of a node with itself. Pairs are
#' ordered in directed networks and unordered in undirected ones. The result
#' keeps the node table and the network's direction; its edge table has
#' `from`, `to`, and `relations` (the number of selected relations in which
#' the pair is tied), plus `weight` for `"sum"`, which becomes the weight
#' column (with the network's `weight_type`, or `"strength"` for counts).
#' Other edge columns are dropped, since they cannot be merged across
#' relations; `"select"` keeps every column.
#'
#' The UCINET menu paths follow Hanneman and Riddle (2005, chapter 16),
#' which also describes Transform > Multiplex, the inverse operation that
#' stacks relations into one nominal-valued matrix.
#'
#' @inheritParams graph_density
#' @param relation The edge column that names each edge's relation, by name
#'   or position; `NULL` looks for a column named `relation` or `network`.
#' @param method How to combine; see Description.
#' @param which Which relations (values of the relation column) to combine;
#'   `NULL` means all of them. `"difference"` needs exactly two, in order,
#'   and `"select"` exactly one.
#'
#' @references
#' Hanneman, R. A., & Riddle, M. (2005). *Introduction to Social Network
#' Methods*. University of California, Riverside.
#'
#' Borgatti, S. P., Everett, M. G., & Freeman, L. C. (2002). *UCINET for
#' Windows: Software for Social Network Analysis*. Analytic Technologies.
#'
#' @return A netframe.
#' @export
combine_relations <- function(x, relation = NULL,
                              method = c("union", "intersection", "sum", "difference", "select"),
                              which = NULL) {
  check_netframe(x)
  method <- match.arg(method)
  relation <- relation_column(x, relation)
  rel <- x$edges[[relation]]
  if (anyNA(rel)) {
    stop(sprintf("Column `%s` of `x$edges` has missing values; every edge needs a relation.",
                 relation), call. = FALSE)
  }
  rel <- as.character(rel)
  all_rel <- unique(rel)
  which_given <- which
  if (is.null(which)) {
    which <- all_rel
  } else {
    if (!is.atomic(which) || anyNA(which)) {
      stop("`which` must be a vector of relation names.", call. = FALSE)
    }
    which <- as.character(which)
    bad <- setdiff(which, all_rel)
    if (length(bad)) {
      stop(sprintf("Unknown relation(s) %s; the network has %s.",
                   paste0("\"", bad, "\"", collapse = ", "),
                   paste0("\"", all_rel, "\"", collapse = ", ")), call. = FALSE)
    }
  }
  # The order of a difference matters, so the two relations must be named.
  if (method == "difference" && (is.null(which_given) || length(which) != 2L)) {
    stop("`method = \"difference\"` needs exactly two relations in `which` (first minus second).",
         call. = FALSE)
  }
  if (method == "select" && length(which) != 1L) {
    stop("`method = \"select\"` needs exactly one relation in `which`.", call. = FALSE)
  }
  if (method == "select") {
    edges <- x$edges[rel == which, , drop = FALSE]
    rownames(edges) <- NULL
    return(with_tables(x, edges, x$nodes))
  }

  n <- nrow(x$nodes)
  src <- x$index$src
  tgt <- x$index$tgt
  a <- if (x$directed) src else pmin(src, tgt)
  b <- if (x$directed) tgt else pmax(src, tgt)
  keep <- rel %in% which
  key <- a[keep] * n + b[keep]
  rel <- rel[keep]
  weighted <- !is.null(x$roles$weight)
  w <- if (weighted) x$edges[[x$roles$weight]][keep] else rep(1, sum(keep))
  # One row per (pair, relation): tied or not; weights add up within it.
  pair_rel <- paste(key, rel)
  first <- !duplicated(pair_rel)
  tied_key <- key[first]
  count <- tabulate(match(tied_key, unique(tied_key)), nbins = length(unique(tied_key)))
  pairs <- unique(tied_key)
  total <- if (weighted) as.numeric(rowsum(w, match(key, pairs))) else count
  selected <- switch(method,
    union = rep(TRUE, length(pairs)),
    intersection = count == length(which),
    sum = rep(TRUE, length(pairs)),
    difference = {
      in_first <- pairs %in% key[rel == which[[1L]]]
      in_second <- pairs %in% key[rel == which[[2L]]]
      in_first & !in_second
    }
  )
  pairs <- pairs[selected]
  ids <- x$nodes[[x$roles$id]]
  edges <- data.frame(from = ids[pairs %/% n + 1L], to = ids[pairs %% n + 1L],
                      relations = count[selected], stringsAsFactors = FALSE)
  weight <- NULL
  weight_type <- "distance"
  if (method == "sum") {
    edges$weight <- as.numeric(total[selected])
    weight <- "weight"
    weight_type <- if (weighted) x$roles$weight_type else "strength"
  }
  rownames(edges) <- NULL
  netframe(edges, x$nodes, directed = x$directed, from = "from", to = "to", id = x$roles$id,
           weight = weight, weight_type = weight_type, graph = x$graph)
}

relation_column <- function(x, relation) {
  if (is.null(relation)) {
    hit <- intersect(c("relation", "network"), names(x$edges))
    if (!length(hit)) {
      stop("No `relation` or `network` column in `x$edges`; name the column with `relation =`.",
           call. = FALSE)
    }
    return(hit[[1L]])
  }
  resolve_column(x$edges, relation, "relation", "edges")
}

#' Dichotomize a weighted network
#'
#' Keeps the edges whose weight passes a test against a threshold and
#' drops the weight role, giving a binary network (UCINET: Transform >
#' Dichotomize, with the same cut-off operators). Edges that fail the test
#' are removed; the weight column stays on the surviving edges as an
#' ordinary attribute, and every other column is kept.
#'
#' @inheritParams graph_density
#' @param threshold The cut-off value.
#' @param op The comparison applied as `weight op threshold`: `">"` (the
#'   default, UCINET's "greater than"), `">="`, `"<"`, `"<="`, `"=="`, or
#'   `"!="`.
#' @param column The numeric edge column to test, by name or position;
#'   `NULL` uses the network's weight column.
#'
#' @references
#' Hanneman, R. A., & Riddle, M. (2005). *Introduction to Social Network
#' Methods*. University of California, Riverside.
#'
#' @return An unweighted netframe.
#' @export
dichotomize <- function(x, threshold = 0, op = c(">", ">=", "<", "<=", "==", "!="),
                        column = NULL) {
  check_netframe(x)
  op <- match.arg(op)
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
    stop("`threshold` must be a single number.", call. = FALSE)
  }
  if (is.null(column)) {
    if (is.null(x$roles$weight)) {
      stop("The network has no weight column; name the column to test with `column =`.",
           call. = FALSE)
    }
    column <- x$roles$weight
  } else {
    column <- resolve_column(x$edges, column, "column", "edges")
  }
  v <- x$edges[[column]]
  if (!is.numeric(v)) {
    stop(sprintf("Column `%s` of `x$edges` must be numeric, not %s.", column,
                 describe_class(v)), call. = FALSE)
  }
  if (anyNA(v)) {
    stop(sprintf("Column `%s` of `x$edges` has missing values; remove or fill them first.",
                 column), call. = FALSE)
  }
  keep <- switch(op, ">" = v > threshold, ">=" = v >= threshold, "<" = v < threshold,
                 "<=" = v <= threshold, "==" = v == threshold, "!=" = v != threshold)
  edges <- x$edges[keep, , drop = FALSE]
  rownames(edges) <- NULL
  with_tables(x, edges, x$nodes, weight = NULL)
}

#' Symmetrize or transpose a directed network
#'
#' `symmetrize()` turns a directed network into an undirected one by
#' combining the two directions of every pair, as UCINET's Transform >
#' Symmetrize does with the method of the same name:
#' \describe{
#'   \item{`"max"`}{A pair is tied if either direction is, with the larger
#'     weight (Boolean OR for binary ties).}
#'   \item{`"min"`}{A pair is tied only if both directions are, with the
#'     smaller weight (Boolean AND).}
#'   \item{`"sum"`}{The weights of the two directions add up (1 + 1 = 2 for a
#'     reciprocated binary tie); the result is weighted.}
#'   \item{`"average"`}{Half the sum (an unreciprocated binary tie becomes
#'     1/2); the result is weighted.}
#'   \item{`"upper"`, `"lower"`}{The tie from the node that comes first
#'     (`"upper"`) or last (`"lower"`) in `x$nodes` decides, as copying the
#'     upper or lower triangle of the adjacency matrix onto the other.}
#' }
#' `transpose_network()` reverses every edge (UCINET: Transform >
#' Transpose).
#'
#' @details
#' An absent tie counts as weight 0, so for weighted networks `"max"`,
#' `"sum"`, and `"average"` tie a pair when either direction is tied and
#' `"min"` when both are. Parallel edges in one direction collapse first,
#' keeping the largest weight (`"min"`, `"max"`, `"upper"`, `"lower"`) or
#' adding up (`"sum"`, `"average"`). Self-loops are kept once per node, with
#' the same rule. The result's edge table has `from`, `to`, and, for
#' weighted input or the `"sum"` and `"average"` methods, `weight` (with the
#' network's `weight_type`, or `"strength"` for counts of binary ties);
#' other edge columns are dropped, since the two directions can carry
#' different values. Each pair appears once with `from` before `to` in node
#' order. `transpose_network()` keeps every column. An undirected network
#' is returned unchanged by `transpose_network()` and is an error for
#' `symmetrize()`.
#'
#' @inheritParams graph_density
#' @param method How the two directions combine; see Description.
#'
#' @references
#' Hanneman, R. A., & Riddle, M. (2005). *Introduction to Social Network
#' Methods*. University of California, Riverside.
#'
#' @return A netframe: undirected for `symmetrize()`.
#' @export
symmetrize <- function(x, method = c("max", "min", "sum", "average", "upper", "lower")) {
  check_netframe(x)
  method <- match.arg(method)
  if (!x$directed) {
    stop("The network is already undirected; symmetrize() needs a directed one.",
         call. = FALSE)
  }
  n <- nrow(x$nodes)
  src <- x$index$src
  tgt <- x$index$tgt
  weighted <- !is.null(x$roles$weight)
  w <- if (weighted) as.numeric(x$edges[[x$roles$weight]]) else rep(1, length(src))
  # Collapse parallel arcs (groups come back in the order of `arcs`).
  arc <- src * n + tgt
  arcs <- unique(arc)
  idx <- match(arc, arcs)
  arc_w <- if (method %in% c("sum", "average")) {
    as.numeric(rowsum(w, idx))
  } else {
    vapply(split(w, idx), max, numeric(1), USE.NAMES = FALSE)
  }
  a_src <- arcs %/% n
  a_tgt <- arcs %% n
  lo <- pmin(a_src, a_tgt)
  hi <- pmax(a_src, a_tgt)
  pair <- lo * n + hi
  pairs <- unique(pair)
  value <- function(sel) {
    # Weights of the forward (lo -> hi) and backward (hi -> lo) arcs, 0 if
    # absent; loops have one arc, read as both directions.
    fwd <- numeric(length(pairs))
    bwd <- numeric(length(pairs))
    idx <- match(pair, pairs)
    is_fwd <- a_src <= a_tgt
    fwd[idx[is_fwd]] <- arc_w[is_fwd]
    bwd[idx[!is_fwd]] <- arc_w[!is_fwd]
    loop <- pairs %/% n == pairs %% n
    bwd[loop] <- fwd[loop]
    sel(fwd, bwd)
  }
  weight <- switch(method,
    max = value(pmax),
    min = value(pmin),
    sum = value(`+`),
    average = value(function(f, b) (f + b) / 2),
    upper = value(function(f, b) f),
    lower = value(function(f, b) b)
  )
  keep <- weight > 0
  ids <- x$nodes[[x$roles$id]]
  edges <- data.frame(from = ids[pairs[keep] %/% n + 1L], to = ids[pairs[keep] %% n + 1L],
                      stringsAsFactors = FALSE)
  weight_col <- NULL
  weight_type <- "distance"
  if (weighted || method %in% c("sum", "average")) {
    edges$weight <- weight[keep]
    weight_col <- "weight"
    weight_type <- if (weighted) x$roles$weight_type else "strength"
  }
  rownames(edges) <- NULL
  netframe(edges, x$nodes, directed = FALSE, from = "from", to = "to", id = x$roles$id,
           weight = weight_col, weight_type = weight_type, graph = x$graph)
}

#' @rdname symmetrize
#' @export
transpose_network <- function(x) {
  check_netframe(x)
  if (!x$directed) return(x)
  edges <- x$edges
  from <- edges[[x$roles$from]]
  edges[[x$roles$from]] <- edges[[x$roles$to]]
  edges[[x$roles$to]] <- from
  with_tables(x, edges, x$nodes)
}

#' Subnetworks: node selection, filters, isolates, pendants, components
#'
#' Induced subnetworks keep a set of nodes and every edge between them:
#' \describe{
#'   \item{`subnetwork()`}{Keeps the nodes named in `nodes` (UCINET: Data >
#'     Extract).}
#'   \item{`filter_nodes()`}{Keeps the nodes for which a condition on the
#'     columns of `x$nodes` holds, written as in [subset()]: for example
#'     `filter_nodes(x, role == "leader")`.}
#'   \item{`remove_isolates()`}{Drops nodes with no ties to other nodes
#'     (UCINET: Data > Remove isolates).}
#'   \item{`remove_pendants()`}{Drops nodes with exactly one neighbor, once
#'     or repeatedly (UCINET: Data > Remove pendants).}
#'   \item{`largest_component()`}{Keeps the largest weak or strong component
#'     (UCINET: Data > Extract main component).}
#' }
#'
#' @details
#' Node and edge attributes, weights, direction, and network attributes are
#' kept. Isolates and pendants are judged on the undirected simple graph
#' underlying the network: a node is an isolate when it has no neighbor
#' other than itself (a self-loop does not save it) and a pendant when it
#' has exactly one such neighbor, whatever the direction and number of
#' edges between them. With `recursive = TRUE`, `remove_pendants()` repeats
#' until no pendant is left, which prunes trees to their core; otherwise
#' only the nodes that are pendants in the original network are dropped.
#' When several components are largest, the one containing the earliest
#' node in `x$nodes` is kept.
#'
#' @inheritParams graph_density
#' @param nodes Node identifiers (values of `x$nodes[[id]]`) to keep, or a
#'   logical vector with one entry per row of `x$nodes`.
#' @param condition A logical expression in the columns of `x$nodes`;
#'   missing values count as `FALSE`.
#' @param recursive Whether to keep removing pendants until none is left.
#' @param mode `"weak"` or `"strong"` components; see [components()].
#'
#' @references
#' Hanneman, R. A., & Riddle, M. (2005). *Introduction to Social Network
#' Methods*. University of California, Riverside.
#'
#' @return A netframe.
#' @export
subnetwork <- function(x, nodes) {
  check_netframe(x)
  keep <- node_selection(x, nodes)
  induced_netframe(x, keep)
}

#' @rdname subnetwork
#' @export
filter_nodes <- function(x, condition) {
  check_netframe(x)
  keep <- eval(substitute(condition), x$nodes, parent.frame())
  if (!is.logical(keep) || !(length(keep) == nrow(x$nodes) || length(keep) == 1L)) {
    stop("`condition` must give a logical value for every row of `x$nodes`.", call. = FALSE)
  }
  keep <- rep_len(keep, nrow(x$nodes)) & !is.na(keep)
  induced_netframe(x, which(keep))
}

#' @rdname subnetwork
#' @export
remove_isolates <- function(x) {
  check_netframe(x)
  induced_netframe(x, which(neighbor_counts(x) > 0L))
}

#' @rdname subnetwork
#' @export
remove_pendants <- function(x, recursive = FALSE) {
  check_netframe(x)
  recursive <- check_flag(recursive, "recursive")
  repeat {
    counts <- neighbor_counts(x)
    if (!any(counts == 1L)) return(x)
    x <- induced_netframe(x, which(counts != 1L))
    if (!recursive) return(x)
  }
}

#' @rdname subnetwork
#' @export
largest_component <- function(x, mode = c("weak", "strong")) {
  check_netframe(x)
  mode <- match.arg(mode)
  if (!nrow(x$nodes)) return(x)
  cc <- components(x, mode)
  induced_netframe(x, which(cc$membership == which.max(cc$csize)))
}

# Distinct neighbors other than the node itself, ignoring direction.
neighbor_counts <- function(x) {
  n <- nrow(x$nodes)
  src <- x$index$src
  tgt <- x$index$tgt
  keep <- src != tgt
  a <- c(src[keep], tgt[keep])
  b <- c(tgt[keep], src[keep])
  key <- unique(a * n + b)
  tabulate(key %/% n + 1L, nbins = n)
}

# Positions of the nodes a user selects by identifier or logical mask.
node_selection <- function(x, nodes) {
  n <- nrow(x$nodes)
  if (is.logical(nodes)) {
    if (length(nodes) != n || anyNA(nodes)) {
      stop("A logical `nodes` must have one TRUE or FALSE per row of `x$nodes`.",
           call. = FALSE)
    }
    return(which(nodes))
  }
  if (!is.atomic(nodes)) {
    stop("`nodes` must be a vector of node identifiers or a logical vector.", call. = FALSE)
  }
  ids <- x$nodes[[x$roles$id]]
  pos <- match(check_ids(nodes, x$roles$id, "nodes", "Node identifiers"),
               check_ids(ids, x$roles$id, "nodes", "Node identifiers"))
  if (anyNA(pos)) {
    missing <- unique(nodes[is.na(pos)])
    stop(sprintf("%d of the requested nodes are not in `x$nodes`, for example: %s.",
                 length(missing),
                 paste(format(utils::head(missing, 5L)), collapse = ", ")), call. = FALSE)
  }
  sort(unique(pos))
}

#' Recode and join node attributes
#'
#' `recode_attribute()` replaces values of a column of `x$nodes` (or
#' `x$edges`) according to a mapping, as UCINET's Transform > Recode does
#' for data values. `join_attributes()` adds the columns of a data frame to
#' `x$nodes`, matched by node identifier (UCINET: Data > Join > Columns, or
#' Data > Attribute for loading an attribute file).
#'
#' @details
#' `map` is a named vector, old value to new (`c(male = "M", female =
#' "F")`), or a function applied to the whole column. Values not named in
#' the vector are left unchanged unless `default` is given, in which case
#' they become `default`; a factor column is recoded through its labels and
#' returned as character. Recoding the identifier column rebuilds the edge
#' index, so new identifiers must stay unique.
#'
#' `join_attributes()` matches `data[[by]]` to `x$nodes[[id]]`; nodes
#' without a matching row get `NA`, rows of `data` without a matching node
#' are dropped, and duplicate keys in `data` are an error. Columns of `data`
#' that already exist in `x$nodes` replace them, except the identifier
#' column.
#'
#' @inheritParams graph_density
#' @param column The column to recode, by name or position.
#' @param map A named vector (old value = new value) or a function.
#' @param default Value for entries not named in `map` (ignored when `map`
#'   is a function); `NULL` leaves them unchanged.
#' @param table `"nodes"` or `"edges"`.
#' @param data A data frame with one row per node to join.
#' @param by The column of `data` holding node identifiers, by name or
#'   position; `NULL` uses a column with the same name as the network's
#'   identifier column, or else the first column.
#'
#' @references
#' Hanneman, R. A., & Riddle, M. (2005). *Introduction to Social Network
#' Methods*. University of California, Riverside.
#'
#' @return A netframe.
#' @export
recode_attribute <- function(x, column, map, default = NULL, table = c("nodes", "edges")) {
  check_netframe(x)
  table <- match.arg(table)
  df <- x[[table]]
  column <- resolve_column(df, column, "column", table)
  v <- df[[column]]
  if (is.factor(v)) v <- as.character(v)
  if (is.function(map)) {
    new <- map(v)
    if (length(new) != length(v)) {
      stop("`map` must return one value per row when it is a function.", call. = FALSE)
    }
  } else {
    if (!is.atomic(map) || is.null(names(map)) || anyNA(names(map)) || !all(nzchar(names(map)))) {
      stop("`map` must be a named vector (old value = new value) or a function.",
           call. = FALSE)
    }
    hit <- match(as.character(v), names(map))
    new <- if (is.null(default)) v else rep_len(default, length(v))
    if (is.character(map) && !is.character(new)) new <- as.character(new)
    new[!is.na(hit)] <- map[hit[!is.na(hit)]]
    names(new) <- NULL
  }
  df[[column]] <- new
  if (table == "nodes") {
    if (column == x$roles$id) {
      # Edge endpoints must follow the identifiers.
      edges <- x$edges
      old <- x$nodes[[x$roles$id]]
      edges[[x$roles$from]] <- new[match(edges[[x$roles$from]], old)]
      edges[[x$roles$to]] <- new[match(edges[[x$roles$to]], old)]
      return(with_tables(x, edges, df))
    }
    with_tables(x, x$edges, df)
  } else {
    with_tables(x, df, x$nodes)
  }
}

#' @rdname recode_attribute
#' @export
join_attributes <- function(x, data, by = NULL) {
  check_netframe(x)
  if (!is.data.frame(data)) {
    stop(sprintf("`data` must be a data frame, not %s.", describe_class(data)), call. = FALSE)
  }
  if (is.null(by)) {
    by <- if (x$roles$id %in% names(data)) x$roles$id else {
      if (!ncol(data)) stop("`data` has no columns.", call. = FALSE)
      names(data)[[1L]]
    }
  } else {
    by <- resolve_column(data, by, "by", "data")
  }
  keys <- check_ids(data[[by]], by, "data", "Join keys")
  if (anyDuplicated(keys)) {
    stop(sprintf("Column `%s` of `data` has duplicated keys; each node may appear once.", by),
         call. = FALSE)
  }
  ids <- check_ids(x$nodes[[x$roles$id]], x$roles$id, "nodes", "Node identifiers")
  pos <- match(ids, keys)
  nodes <- x$nodes
  for (col in setdiff(names(data), by)) {
    nodes[[col]] <- data[[col]][pos]
  }
  with_tables(x, x$edges, nodes)
}

# A netframe with the same roles and attributes as `x` but new tables (and
# possibly no weight role).
with_tables <- function(x, edges, nodes, directed = x$directed, weight = x$roles$weight) {
  weight_type <- if (is.null(x$roles$weight_type)) "distance" else x$roles$weight_type
  out <- netframe(edges, nodes, directed = directed, from = x$roles$from, to = x$roles$to,
                  id = x$roles$id, weight = weight, weight_type = weight_type,
                  graph = x$graph)
  extra <- setdiff(names(x$roles), names(out$roles))
  if (length(extra)) out$roles[extra] <- x$roles[extra]
  out
}
