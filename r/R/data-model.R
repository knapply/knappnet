# Helpers for the configurations of the data model (docs/DATA_MODEL.md):
# resolving missing data, signed values, multilayer networks (edge
# classes, layers as networks, lifting a multiplex network), ego networks
# held in one net_frame (building, splitting, and the bridge to egor's
# three tables), node presence spells, and node attributes at a time.

# ---- Missing data (section 13) ------------------------------------------------

# Resolve non-respondents and unobserved pairs (man/net_resolve_missing.Rd).
net_resolve_missing <- function(x, na_action) {
  check_net_frame(x)
  if (missing(na_action)) {
    stop(paste0("`na_action` must say how to resolve the missing data: \"absent\" (an ",
                "unobserved pair is no tie) or \"present\" (it is a tie). There is no ",
                "default."), call. = FALSE)
  }
  check_choice(na_action, c("absent", "present"), "na_action")
  r <- x$roles
  edges <- x$edges
  if (!is.null(r$node_missing) && na_action == "present" &&
      any(x$nodes[[r$node_missing]] %in% TRUE)) {
    stop(paste0("`na_action = \"present\"` would tie every non-respondent to every node ",
                "whose tie it did not report; use \"absent\", or impute those ties as rows ",
                "of `x$edges` and rebuild the network."), call. = FALSE)
  }
  if (!is.null(r$edge_missing) && na_action == "absent") {
    keep <- !(edges[[r$edge_missing]] %in% TRUE)
    edges <- edges[keep, , drop = FALSE]
    rownames(edges) <- NULL
  }
  rebuild_net_frame(x, edges, x$nodes, edge_missing = NULL, node_missing = NULL,
                    fn = "net_resolve_missing")
}

# ---- Signed networks (section 5) ------------------------------------------------

# A signed value column as sign and magnitude (man/net_as_signed.Rd).
net_as_signed <- function(edges, value, sign = "sign", weight = "weight") {
  if (!is.data.frame(edges)) {
    stop(sprintf("`edges` must be a data frame, not %s.", describe_class(edges)),
         call. = FALSE)
  }
  value <- role_column(edges, value, "value", "edges")
  sign <- check_string(sign, "sign")
  weight <- check_string(weight, "weight")
  v <- edges[[value]]
  if (!is.numeric(v) || is.object(v)) {
    stop(sprintf("Column `%s` of `edges` must be numeric, not %s.", value, describe_class(v)),
         call. = FALSE)
  }
  bad <- which(is.na(v) | v == 0)
  if (length(bad)) {
    stop(sprintf(paste0("Column `%s` of `edges` has zero or missing values in %s; a zero ",
                        "is neither positive nor negative, so drop those rows first."),
                 value, rows_text(bad)), call. = FALSE)
  }
  edges[[sign]] <- ifelse(v < 0, -1L, 1L)
  edges[[weight]] <- abs(v)
  edges
}

# ---- Multilayer networks (section 7) ---------------------------------------------

# See man/net_split_layers.Rd.
net_is_multilayer <- function(x) inherits(x, "multilayer_net_frame")

check_multilayer <- function(x) {
  check_net_frame(x)
  if (!inherits(x, "multilayer_net_frame")) {
    stop(paste0("`x` must be a multilayer network (class `multilayer_net_frame`), with ",
                "`source_layer` and `target_layer` edge columns; create one with ",
                "net_frame(source_layer = ..., target_layer = ...) or net_as_multilayer()."),
         call. = FALSE)
  }
}

# See man/net_split_layers.Rd.
net_edge_layers <- function(x) {
  check_multilayer(x)
  engine_graph(x, missing_ok = TRUE)
  r <- engine_call("net_edge_layers",
                   .Call(C_frame_edge_layers, x$cache$frame, as.double(nrow(x$edges))))
  labels <- frame_summary_cached(x)$layers
  data.frame(source_layer = labels[r$source_layer + 1L],
             target_layer = labels[r$target_layer + 1L],
             class = c("intralayer", "coupling", "interlayer")[r$class + 1L],
             stringsAsFactors = FALSE)
}

# Layers as networks (man/net_split_layers.Rd).
net_split_layers <- function(x) {
  check_multilayer(x)
  el <- net_edge_layers(x)
  labels <- frame_summary_cached(x)$layers
  out <- lapply(labels, function(label) {
    keep <- el$class == "intralayer" & el$source_layer == label
    edges <- x$edges[keep, , drop = FALSE]
    rownames(edges) <- NULL
    rebuild_net_frame(x, edges, x$nodes, source_layer = NULL, target_layer = NULL,
                      fn = "net_split_layers")
  })
  names(out) <- labels
  out
}

# A multiplex network in layer form (man/net_split_layers.Rd).
net_as_multilayer <- function(x, couplings = c("categorical", "ordinal", "none"),
                              relation = NULL, coupling_weight = 1) {
  check_net_frame(x)
  couplings <- match.arg(couplings)
  if (inherits(x, "multilayer_net_frame")) {
    stop("`x` is already a multilayer network.", call. = FALSE)
  }
  if (!is.null(x$roles$start) || !is.null(x$roles$wave)) {
    stop(paste0("Coupling edges would need times of their own; lift a temporal network ",
                "one period at a time, or build the layer columns yourself."), call. = FALSE)
  }
  relation <- relation_column(x, relation)
  rel <- x$edges[[relation]]
  if (anyNA(rel)) {
    stop(sprintf("Column `%s` of `x$edges` has missing values; every edge needs a relation.",
                 relation), call. = FALSE)
  }
  rel <- as.character(rel)
  labels <- unique(rel)
  if (!is.null(x$roles$weight)) {
    coupling_weight <- check_number(coupling_weight, "coupling_weight", lower = 0,
                                    open = identical(x$roles$weight_type, "distance"))
  }
  edges <- x$edges
  cols <- unique_names(c("source_layer", "target_layer"), names(edges))
  edges[[cols[[1L]]]] <- rel
  edges[[cols[[2L]]]] <- rel
  pairs <- switch(couplings,
    categorical = if (length(labels) > 1L) t(utils::combn(labels, 2L)) else NULL,
    ordinal = if (length(labels) > 1L) cbind(labels[-length(labels)], labels[-1L]) else NULL,
    none = NULL
  )
  if (!is.null(pairs) && nrow(pairs)) {
    ids <- x$nodes[[x$roles$id]]
    if (is.factor(ids)) ids <- as.character(ids)
    n <- length(ids)
    a <- rep(pairs[, 1L], each = n)
    b <- rep(pairs[, 2L], each = n)
    if (x$directed) {
      a0 <- a
      a <- c(a, b)
      b <- c(b, a0)
    }
    k <- length(a)
    who <- rep(ids, length.out = k)
    extra <- edges[rep(NA_integer_, k), , drop = FALSE]
    extra[[x$roles$source]] <- who
    extra[[x$roles$target]] <- who
    extra[[cols[[1L]]]] <- a
    extra[[cols[[2L]]]] <- b
    if (!is.null(x$roles$weight)) extra[[x$roles$weight]] <- coupling_weight
    if (!is.null(x$roles$sign)) extra[[x$roles$sign]] <- 1L
    if (!is.null(x$roles$edge_missing)) extra[[x$roles$edge_missing]] <- FALSE
    if (!is.null(x$roles$edge_ego)) extra[[x$roles$edge_ego]] <- x$nodes[[x$roles$ego]][match(who, ids)]
    edges <- rbind(edges, extra)
    rownames(edges) <- NULL
  }
  layer_col <- "layer"
  layers <- data.frame(labels, stringsAsFactors = FALSE)
  names(layers) <- layer_col
  declared <- x$declared
  declared$couplings_only <- TRUE
  out <- new_net_frame(edges, x$nodes, x$directed,
                       c(x$roles[setdiff(names(x$roles), c("relation", "node_geometry",
                                                           "edge_geometry"))],
                         list(source_layer = cols[[1L]], target_layer = cols[[2L]],
                              layer = layer_col)),
                       declared, x$graph,
                       tables = c(frame_tables(x)[c("node_spells", "node_changes")],
                                  list(layers = layers)))
  build_net_frame(out, quiet = TRUE, fn = "net_as_multilayer")
}

# Names not already in `taken`, by appending underscores.
unique_names <- function(wanted, taken) {
  out <- wanted
  for (k in seq_along(out)) {
    while (out[[k]] %in% c(taken, out[-k])) out[[k]] <- paste0(out[[k]], "_")
  }
  out
}

# ---- Ego networks (section 12) ----------------------------------------------------

# Ego networks in one net_frame (man/net_ego_frame.Rd).
net_ego_frame <- function(x, egos = NULL, order = 1L, mode = c("all", "out", "in"),
                          include_ego = TRUE, sep = ":", threads = 0L) {
  check_net_frame(x)
  mode <- match.arg(mode)
  order <- check_count(order, "order", min = 1L)
  include_ego <- check_flag(include_ego, "include_ego")
  sep <- check_string(sep, "sep")
  if (!is.null(x$roles$ego)) {
    stop("`x` already holds ego networks (it has an ego role).", call. = FALSE)
  }
  pos <- if (is.null(egos)) seq_len(nrow(x$nodes)) else node_positions(x, egos, "egos")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  members <- engine_call("net_ego_frame",
                         .Call(C_ego_members, ptr, as.integer(pos - 1L), order,
                               mode_code(mode), if (include_ego) 0L else 1L, threads))
  ids <- node_ids(x)
  src <- x$index$src + 1L
  tgt <- x$index$tgt + 1L
  ego_col <- unique_names("ego", c(names(x$nodes), names(x$edges)))
  node_col <- unique_names("node", names(x$nodes))
  node_parts <- vector("list", length(pos))
  edge_parts <- vector("list", length(pos))
  for (k in seq_along(pos)) {
    keep <- sort(as.integer(members[[k]]))
    e <- pos[[k]]
    new_id <- ifelse(keep == e, ids[keep], paste0(ids[e], sep, ids[keep]))
    nd <- x$nodes[keep, , drop = FALSE]
    nd[[node_col]] <- ids[keep]
    nd[[x$roles$id]] <- new_id
    nd[[ego_col]] <- rep(ids[e], length(keep))
    inside <- match(seq_len(nrow(x$nodes)), keep)
    rows <- which(!is.na(inside[src]) & !is.na(inside[tgt]))
    ed <- x$edges[rows, , drop = FALSE]
    ed[[x$roles$source]] <- new_id[inside[src[rows]]]
    ed[[x$roles$target]] <- new_id[inside[tgt[rows]]]
    ed[[ego_col]] <- rep(ids[e], length(rows))
    node_parts[[k]] <- nd
    edge_parts[[k]] <- ed
  }
  nodes <- do.call(rbind, node_parts)
  edges <- do.call(rbind, edge_parts)
  rownames(nodes) <- NULL
  rownames(edges) <- NULL
  if (anyDuplicated(nodes[[x$roles$id]])) {
    stop(sprintf(paste0("Combining ego and alter identifiers with `sep = \"%s\"` gives ",
                        "repeated identifiers; choose another `sep`."), sep), call. = FALSE)
  }
  declared <- x$declared
  declared$egos_are_nodes <- include_ego
  out <- rebuild_net_frame(x, edges, nodes, ego = ego_col, edge_ego = ego_col,
                           declared = declared, fn = "net_ego_frame")
  out$graph$ego_order <- order
  out
}

# See man/net_ego_frame.Rd.
net_split_egos <- function(x) {
  check_net_frame(x)
  if (is.null(x$roles$ego)) {
    stop("`x` has no ego role; build it with net_frame(ego = ...) or net_ego_frame().",
         call. = FALSE)
  }
  ego <- as.character(x$nodes[[x$roles$ego]])
  labels <- unique(ego)
  edge_ego <- ego[x$index$src + 1L]
  out <- lapply(labels, function(label) {
    nodes <- x$nodes[ego == label, , drop = FALSE]
    edges <- x$edges[edge_ego == label, , drop = FALSE]
    rownames(nodes) <- NULL
    rownames(edges) <- NULL
    rebuild_net_frame(x, edges, nodes, ego = NULL, fn = "net_split_egos")
  })
  names(out) <- labels
  out
}

# egor's three tables as one net_frame (man/net_ego_frame.Rd).
as_net_frame.egor <- function(x, directed = FALSE, sep = ":", quiet = FALSE, ...) {
  sep <- check_string(sep, "sep")
  plain <- function(tb) {
    if (!is.data.frame(tb) && is.list(tb) && !is.null(tb$variables)) tb <- tb$variables
    as.data.frame(tb, stringsAsFactors = FALSE)
  }
  ego <- plain(x$ego)
  alter <- plain(x$alter)
  aatie <- if (is.null(x$aatie)) NULL else plain(x$aatie)
  need <- list(list(ego, ".egoID", "ego"), list(alter, c(".egoID", ".altID"), "alter"))
  if (!is.null(aatie)) need <- c(need, list(list(aatie, c(".egoID", ".srcID", ".tgtID"), "aatie")))
  for (n in need) {
    miss <- setdiff(n[[2L]], names(n[[1L]]))
    if (length(miss)) {
      stop(sprintf("The %s table of the egor object has no column %s.", n[[3L]],
                   paste0("`", miss, "`", collapse = ", ")), call. = FALSE)
    }
  }
  ego_id <- as.character(ego$.egoID)
  alter_id <- paste0(alter$.egoID, sep, alter$.altID)
  nodes <- bind_rows_fill(
    cbind(data.frame(id = ego_id, ego = ego_id, stringsAsFactors = FALSE),
          ego[setdiff(names(ego), c(".egoID", "id", "ego"))]),
    cbind(data.frame(id = alter_id, ego = as.character(alter$.egoID),
                     stringsAsFactors = FALSE),
          alter[setdiff(names(alter), c(".egoID", ".altID", "id", "ego"))])
  )
  edges <- data.frame(source = as.character(alter$.egoID), target = alter_id,
                      ego = as.character(alter$.egoID), tie = "ego-alter",
                      stringsAsFactors = FALSE)
  if (!is.null(aatie) && nrow(aatie)) {
    aa <- cbind(data.frame(source = paste0(aatie$.egoID, sep, aatie$.srcID),
                           target = paste0(aatie$.egoID, sep, aatie$.tgtID),
                           ego = as.character(aatie$.egoID), tie = "alter-alter",
                           stringsAsFactors = FALSE),
                aatie[setdiff(names(aatie), c(".egoID", ".srcID", ".tgtID", "source",
                                              "target", "ego", "tie"))])
    edges <- bind_rows_fill(edges, aa)
  }
  net_frame(edges, nodes, directed = directed, ego = "ego", edge_ego = "ego",
            egos_are_nodes = TRUE, quiet = quiet)
}

# Ego networks as egor's three tables (man/net_ego_frame.Rd).
net_to_egor <- function(x) {
  check_net_frame(x)
  require_pkg("egor", "net_to_egor()")
  if (is.null(x$roles$ego)) {
    stop("`x` has no ego role; build it with net_frame(ego = ...) or net_ego_frame().",
         call. = FALSE)
  }
  ids <- node_ids(x)
  ego <- as.character(x$nodes[[x$roles$ego]])
  is_ego <- ids == ego
  if (!all(unique(ego) %in% ids[is_ego])) {
    stop("net_to_egor() needs every ego as a node row (`egos_are_nodes = TRUE`).",
         call. = FALSE)
  }
  attrs <- setdiff(names(x$nodes), c(x$roles$id, x$roles$ego))
  egos <- cbind(data.frame(.egoID = ids[is_ego], stringsAsFactors = FALSE),
                x$nodes[is_ego, attrs, drop = FALSE])
  alters <- cbind(data.frame(.egoID = ego[!is_ego], .altID = ids[!is_ego],
                             stringsAsFactors = FALSE),
                  x$nodes[!is_ego, attrs, drop = FALSE])
  src <- x$index$src + 1L
  tgt <- x$index$tgt + 1L
  aa <- !is_ego[src] & !is_ego[tgt]
  eattrs <- setdiff(names(x$edges), unlist(x$roles[intersect(names(x$roles),
                                                             edge_role_names)]))
  aaties <- cbind(data.frame(.egoID = ego[src[aa]], .srcID = ids[src[aa]],
                             .tgtID = ids[tgt[aa]], stringsAsFactors = FALSE),
                  x$edges[aa, eattrs, drop = FALSE])
  rownames(egos) <- rownames(alters) <- rownames(aaties) <- NULL
  egor::egor(alters = alters, egos = egos, aaties = aaties,
             ID.vars = list(ego = ".egoID", alter = ".altID", source = ".srcID",
                            target = ".tgtID"))
}

# Row-binds data frames with different columns, filling with NA.
bind_rows_fill <- function(...) {
  dfs <- list(...)
  cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(d) {
    for (col in setdiff(cols, names(d))) d[[col]] <- rep(NA, nrow(d))
    d[cols]
  })
  out <- do.call(rbind, dfs)
  rownames(out) <- NULL
  out
}

# ---- Node spells and changes (section 9) ------------------------------------------

# Presence spells (man/net_node_spells.Rd).
net_node_spells <- function(x) {
  check_net_frame(x)
  engine_graph(x, missing_ok = TRUE)
  n <- nrow(x$nodes)
  r <- engine_call("net_node_spells", .Call(C_frame_node_spells, x$cache$frame, as.double(n)))
  counts <- diff(r$offsets)
  ids <- x$nodes[[x$roles$id]]
  if (is.factor(ids)) ids <- as.character(ids)
  node <- rep(seq_len(n), counts)
  out <- data.frame(ids[node], stringsAsFactors = FALSE)
  names(out) <- x$roles$id
  timed <- !is.null(x$roles$time_class)
  out$start <- if (timed) from_seconds(r$start, x) else r$start
  out$end <- if (timed) from_seconds(r$end, x) else r$end
  out
}

# Node attributes at a time (man/net_node_spells.Rd).
net_nodes_at <- function(x, time) {
  check_net_frame(x)
  nodes <- x$nodes
  ch <- x$node_changes
  if (is.null(ch) || !nrow(ch)) return(nodes)
  t <- time_arg(time, x, "time")
  id <- x$roles$id
  tc <- time_values(ch[[x$roles$change_time]], x$roles$change_time, "node_changes")
  cols <- setdiff(names(ch), c(id, x$roles$change_time))
  ord <- order(tc)
  ch <- ch[ord, , drop = FALSE]
  tc <- tc[ord]
  done <- ch[tc <= t, , drop = FALSE]
  if (!nrow(done)) return(nodes)
  pos <- match(as.character(done[[id]]), as.character(nodes[[id]]))
  for (col in cols) {
    if (!col %in% names(nodes)) nodes[[col]] <- rep(done[[col]][NA_integer_], nrow(nodes))
    for (k in seq_along(pos)) nodes[[col]][pos[[k]]] <- done[[col]][[k]]
  }
  nodes
}
