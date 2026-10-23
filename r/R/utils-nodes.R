# Helpers for arguments that refer to nodes: identifiers, positions,
# attributes given as a column name or a vector, and subsets of a network
# that keep its roles and attributes.

# Node identifiers of a net_frame as character.
node_ids <- function(x) as.character(x$nodes[[x$roles$id]])

# Positions (1-based) in `x$nodes` of node identifiers given as values of
# the identifier column; NULL stays NULL (unless nodes are required). `one`
# requires exactly one identifier, `empty_ok = FALSE` at least one, and
# `distinct` rejects a node named twice. Factors are matched by their
# labels.
node_positions <- function(x, nodes, name, one = FALSE, empty_ok = TRUE, distinct = FALSE) {
  if (is.null(nodes) && empty_ok && !one) return(NULL)
  if (!is.atomic(nodes) || is.null(nodes)) {
    stop(sprintf("`%s` must be a vector of node identifiers.", name), call. = FALSE)
  }
  if (one && length(nodes) != 1L) {
    stop(sprintf("`%s` must be a single node identifier.", name), call. = FALSE)
  }
  if (!empty_ok && !length(nodes)) {
    stop(sprintf("`%s` must name at least one node.", name), call. = FALSE)
  }
  if (anyNA(nodes)) {
    stop(sprintf("`%s` must not contain missing values (NA).", name), call. = FALSE)
  }
  ids <- x$nodes[[x$roles$id]]
  if (is.factor(ids)) ids <- as.character(ids)
  if (is.factor(nodes)) nodes <- as.character(nodes)
  pos <- match(nodes, ids)
  if (anyNA(pos)) {
    missing <- unique(nodes[is.na(pos)])
    stop(sprintf("%d node(s) in `%s` are not in the network, for example: %s.",
                 length(missing), name,
                 paste(format(utils::head(missing, 5L)), collapse = ", ")), call. = FALSE)
  }
  if (distinct && anyDuplicated(pos)) {
    stop(sprintf("`%s` lists a node more than once.", name), call. = FALSE)
  }
  pos
}

# Lists of node identifiers from the engine's lists of 1-based positions
# (cliques, bicomponents, ...).
node_sets <- function(x, idx) {
  ids <- x$nodes[[x$roles$id]]
  lapply(idx, function(i) ids[i])
}

# ---- Node attributes --------------------------------------------------------------

# A node attribute given as the name of a column of `x$nodes` or as an
# atomic vector with one value per node. A single string is read as a
# column name, except in a one-node network where it is not a column.
node_column <- function(x, values, name) {
  n <- nrow(x$nodes)
  if (is.character(values) && length(values) == 1L &&
      (values %in% names(x$nodes) || n != 1L)) {
    if (!values %in% names(x$nodes)) {
      stop(sprintf("`%s` names no column of `x$nodes`.", name), call. = FALSE)
    }
    values <- x$nodes[[values]]
  }
  if (!is.atomic(values) || length(values) != n) {
    stop(sprintf("`%s` must have one value per node (%d), or name a column of `x$nodes`.",
                 name, n), call. = FALSE)
  }
  values
}

# A grouping of the nodes (a column name or a vector aligned with
# `x$nodes`, without missing values): the values with factors as character,
# their distinct labels in order of first appearance, and 0-based codes
# into those labels.
node_groups <- function(x, groups, name) {
  groups <- node_column(x, groups, name)
  if (is.factor(groups)) groups <- as.character(groups)
  if (anyNA(groups)) {
    stop(sprintf("`%s` must not contain missing values.", name), call. = FALSE)
  }
  labels <- unique(groups)
  list(values = groups, labels = labels, codes = match(groups, labels) - 1L)
}

# node_groups() for measures that can leave nodes out: missing values get
# code NA and `levels` are the labels as character.
node_categories <- function(x, values, name, missing_ok) {
  values <- node_column(x, values, name)
  if (!missing_ok && anyNA(values)) {
    stop(sprintf("`%s` must not contain missing values.", name), call. = FALSE)
  }
  plain <- if (is.factor(values)) as.character(values) else values
  levels <- unique(plain[!is.na(plain)])
  list(values = values, codes = match(plain, levels) - 1L, levels = as.character(levels))
}

# A numeric node attribute without missing values, as a double vector.
node_values <- function(x, values, name) {
  values <- node_column(x, values, name)
  if (!is.numeric(values) || anyNA(values)) {
    stop(sprintf(paste0("`%s` must be a numeric vector with one value per node (%d) and ",
                        "no missing values, or the name of such a column of `x$nodes`."),
                 name, nrow(x$nodes)), call. = FALSE)
  }
  as.double(values)
}

# ---- Derived networks -------------------------------------------------------------

# A net_frame with the same roles, declarations, and network attributes as
# `x` but new tables (and possibly no weight role or another direction);
# see rebuild_net_frame() in net_frame.R.
with_tables <- function(x, edges, nodes, directed = x$directed, weight = x$roles$weight) {
  rebuild_net_frame(x, edges, nodes, directed = directed, weight = weight)
}

# The subgraph induced by the nodes at 1-based positions `keep` (in
# ascending order): those nodes and every edge between two of them, with
# all columns and roles.
induced_net_frame <- function(x, keep) {
  inside <- logical(nrow(x$nodes))
  inside[keep] <- TRUE
  e <- inside[x$index$src + 1L] & inside[x$index$tgt + 1L]
  nodes <- x$nodes[keep, , drop = FALSE]
  edges <- x$edges[e, , drop = FALSE]
  rownames(nodes) <- NULL
  rownames(edges) <- NULL
  with_tables(x, edges, nodes)
}
