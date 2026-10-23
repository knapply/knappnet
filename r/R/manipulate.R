# Transformations of a network's tables: combining relations,
# dichotomizing and symmetrizing, subnetworks and node filters, and recoding
# or joining attributes. Each returns a new netframe built through
# with_tables() or induced_netframe() (utils-nodes.R), which keep the
# roles and attributes of the input.

# Combine the relations of a multiplex network (man/combine_relations.Rd).
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
  edges <- data.frame(source = ids[pairs %/% n + 1L], target = ids[pairs %% n + 1L],
                      relations = count[selected], stringsAsFactors = FALSE)
  weight <- NULL
  weight_type <- NULL
  if (method == "sum") {
    edges$weight <- as.numeric(total[selected])
    weight <- "weight"
    weight_type <- if (weighted) x$roles$weight_type else "strength"
  }
  rownames(edges) <- NULL
  net_frame(edges, x$nodes, directed = x$directed, id = x$roles$id, weight = weight,
            weight_type = weight_type, quiet = TRUE, graph = x$graph)
}

# The relation column: the network's relation role, or a column named
# here.
relation_column <- function(x, relation) {
  if (is.null(relation)) {
    if (is.null(x$roles$relation)) {
      stop(paste0("`x` has no relation role; build it with net_frame(relation = ...) or ",
                  "name the column with `relation =`."), call. = FALSE)
    }
    return(x$roles$relation)
  }
  resolve_column(x$edges, relation, "relation", "edges")
}

# The edges of the chosen relations of a multiplex network, as a netframe
# with the same roles; `relation` NULL keeps every edge. Measures that take
# `relation =` go through this.
relation_subset <- function(x, relation, fn) {
  if (is.null(relation)) return(x)
  if (is.null(x$roles$relation)) {
    stop(sprintf("%s(): `relation` needs a network with a relation role; build it with ",
                 "net_frame(relation = ...).", fn), call. = FALSE)
  }
  rel <- as.character(x$edges[[x$roles$relation]])
  if (!is.atomic(relation) || anyNA(relation)) {
    stop(sprintf("%s(): `relation` must be a vector of relation names.", fn), call. = FALSE)
  }
  relation <- as.character(relation)
  bad <- setdiff(relation, unique(rel))
  if (length(bad)) {
    stop(sprintf("%s(): unknown relation(s) %s; the network has %s.", fn,
                 paste0("\"", bad, "\"", collapse = ", "),
                 paste0("\"", unique(rel), "\"", collapse = ", ")), call. = FALSE)
  }
  edges <- x$edges[rel %in% relation, , drop = FALSE]
  rownames(edges) <- NULL
  rebuild_netframe(x, edges, x$nodes)
}

# Split a multiplex network into one network per relation (man/net_combine_relations.Rd).
net_split_relations <- function(x, relation = NULL) {
  check_netframe(x)
  relation <- relation_column(x, relation)
  rel <- x$edges[[relation]]
  if (anyNA(rel)) {
    stop(sprintf("Column `%s` of `x$edges` has missing values; every edge needs a relation.",
                 relation), call. = FALSE)
  }
  rel <- as.character(rel)
  labels <- unique(rel)
  out <- lapply(labels, function(label) {
    edges <- x$edges[rel == label, , drop = FALSE]
    rownames(edges) <- NULL
    rebuild_netframe(x, edges, x$nodes, relation = NULL)
  })
  names(out) <- labels
  out
}

# Dichotomize a weighted network (man/dichotomize.Rd).
dichotomize <- function(x, threshold = 0, op = c(">", ">=", "<", "<=", "==", "!="),
                        column = NULL) {
  check_netframe(x)
  op <- match.arg(op)
  threshold <- check_number(threshold, "threshold", finite = FALSE)
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

# Symmetrize or transpose a directed network (man/symmetrize.Rd).
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
  edges <- data.frame(source = ids[pairs[keep] %/% n + 1L],
                      target = ids[pairs[keep] %% n + 1L], stringsAsFactors = FALSE)
  weight_col <- NULL
  weight_type <- NULL
  if (weighted || method %in% c("sum", "average")) {
    edges$weight <- weight[keep]
    weight_col <- "weight"
    weight_type <- if (weighted) x$roles$weight_type else "strength"
  }
  rownames(edges) <- NULL
  net_frame(edges, x$nodes, directed = FALSE, id = x$roles$id, weight = weight_col,
            weight_type = weight_type, quiet = TRUE, graph = x$graph)
}

# See man/symmetrize.Rd.
transpose_network <- function(x) {
  check_netframe(x)
  if (!x$directed) return(x)
  edges <- x$edges
  from <- edges[[x$roles$source]]
  edges[[x$roles$source]] <- edges[[x$roles$target]]
  edges[[x$roles$target]] <- from
  with_tables(x, edges, x$nodes)
}

# Subnetworks: node selection, filters, isolates, pendants, components (man/subnetwork.Rd).
subnetwork <- function(x, nodes) {
  check_netframe(x)
  keep <- node_selection(x, nodes)
  induced_netframe(x, keep)
}

# See man/subnetwork.Rd.
filter_nodes <- function(x, condition) {
  check_netframe(x)
  keep <- eval(substitute(condition), x$nodes, parent.frame())
  if (!is.logical(keep) || !(length(keep) == nrow(x$nodes) || length(keep) == 1L)) {
    stop("`condition` must give a logical value for every row of `x$nodes`.", call. = FALSE)
  }
  keep <- rep_len(keep, nrow(x$nodes)) & !is.na(keep)
  induced_netframe(x, which(keep))
}

# See man/subnetwork.Rd.
remove_isolates <- function(x) {
  check_netframe(x)
  induced_netframe(x, which(neighbor_counts(x) > 0L))
}

# See man/subnetwork.Rd.
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

# See man/subnetwork.Rd.
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

# Recode and join node attributes (man/recode_attribute.Rd).
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
      edges[[x$roles$source]] <- new[match(edges[[x$roles$source]], old)]
      edges[[x$roles$target]] <- new[match(edges[[x$roles$target]], old)]
      return(with_tables(x, edges, df))
    }
    with_tables(x, x$edges, df)
  } else {
    with_tables(x, df, x$nodes)
  }
}

# See man/recode_attribute.Rd.
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
