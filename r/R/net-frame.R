# The net_frame class: the constructor, which hands the role columns of the
# node and edge tables (and of the optional layer, node spell, and node
# change tables) to the engine for validation and construction, the
# diagnostics it reports, the class vector, the summary, and parallel
# edges. docs/DATA_MODEL.md is the specification. The argument checks are
# in utils-checks.R, the engine handle in utils-engine.R, the print card
# in print-net-frame.R, and summary() in descriptive.R.
#
# A net_frame is a list: `nodes` and `edges` (the two tables, untouched
# except for nodes appended from the edges), `layers`, `node_spells`, and
# `node_changes` (the optional long tables, or NULL), `directed`, `roles`
# (which columns carry structural meaning, plus the geometry columns, the
# coordinate reference system, and how times are stored), `declared` (the
# declared properties: TRUE or FALSE when declared, absent otherwise),
# `graph` (network-level attributes), `index` (0-based node positions of
# every edge row's endpoints, from the engine), and `cache` (an
# environment holding the engine handles and the summary).

# Create a network from node and edge data frames (man/net_frame.Rd).
net_frame <- function(edges, nodes = NULL, directed = TRUE, source = "source",
                      target = "target", id = "id", weight = NULL, weight_type = NULL,
                      mode = NULL, between_modes_only = NULL, relation = NULL, sign = NULL,
                      start = NULL, end = NULL, wave = NULL, node_start = NULL,
                      node_end = NULL, coords = NULL, crs = NULL, source_layer = NULL,
                      target_layer = NULL, layers = NULL, layer = "layer",
                      couplings_only = NULL, ego = NULL, edge_ego = NULL,
                      egos_are_nodes = NULL, node_missing = NULL, edge_missing = NULL,
                      node_spells = NULL, node_changes = NULL, spell_start = "start",
                      spell_end = "end", change_time = "time", overlapping_spells = NULL,
                      loops = NULL, multiple = NULL, add_nodes = FALSE, quiet = FALSE,
                      graph = list()) {
  a <- as.list(environment())
  quiet <- check_flag(quiet, "quiet")
  build_net_frame(prepare_frame(a), quiet, "net_frame")
}

# An unbuilt net_frame from the arguments of net_frame() as a list `a`
# (shared by the constructors that add to them, such as net_temporal()).
prepare_frame <- function(a) {
  formals_ <- formals(net_frame)
  for (nm in names(formals_)) {
    if (!nm %in% names(a)) a[nm] <- list(eval(formals_[[nm]]))
  }
  unknown <- setdiff(names(a), names(formals_))
  if (length(unknown)) {
    stop(sprintf("Unknown argument%s: %s.", if (length(unknown) > 1L) "s" else "",
                 paste0("`", unknown, "`", collapse = ", ")), call. = FALSE)
  }
  edges <- a$edges
  nodes <- a$nodes
  directed <- a$directed
  source <- a$source
  target <- a$target
  id <- a$id
  weight <- a$weight
  weight_type <- a$weight_type
  mode <- a$mode
  between_modes_only <- a$between_modes_only
  relation <- a$relation
  sign <- a$sign
  start <- a$start
  end <- a$end
  wave <- a$wave
  node_start <- a$node_start
  node_end <- a$node_end
  coords <- a$coords
  crs <- a$crs
  source_layer <- a$source_layer
  target_layer <- a$target_layer
  layers <- a$layers
  layer <- a$layer
  couplings_only <- a$couplings_only
  ego <- a$ego
  edge_ego <- a$edge_ego
  egos_are_nodes <- a$egos_are_nodes
  node_missing <- a$node_missing
  edge_missing <- a$edge_missing
  node_spells <- a$node_spells
  node_changes <- a$node_changes
  spell_start <- a$spell_start
  spell_end <- a$spell_end
  change_time <- a$change_time
  overlapping_spells <- a$overlapping_spells
  loops <- a$loops
  multiple <- a$multiple
  add_nodes <- a$add_nodes
  graph <- a$graph
  check_tables(edges, nodes)
  directed <- check_flag(directed, "directed")
  add_nodes <- check_flag(add_nodes, "add_nodes")
  graph <- check_graph(graph)
  for (nm in c("layers", "node_spells", "node_changes")) {
    tb <- a[[nm]]
    if (!is.null(tb) && !is.data.frame(tb)) {
      stop(sprintf("`%s` must be a data frame or NULL, not %s.", nm, describe_class(tb)),
           call. = FALSE)
    }
  }
  if (!is.null(node_spells) && is.null(nodes)) {
    stop("`node_spells` needs a node table: give `nodes`.", call. = FALSE)
  }
  if (!is.null(node_changes) && is.null(nodes)) {
    stop("`node_changes` needs a node table: give `nodes`.", call. = FALSE)
  }
  roles <- list(
    id = check_string(id, "id"),
    source = role_column(edges, source, "source", "edges"),
    target = role_column(edges, target, "target", "edges"),
    weight = role_optional(edges, weight, "weight", "edges"),
    weight_type = check_weight_type(weight, weight_type),
    mode = role_optional(nodes, mode, "mode", "nodes"),
    relation = role_optional(edges, relation, "relation", "edges"),
    sign = role_optional(edges, sign, "sign", "edges"),
    start = role_optional(edges, start, "start", "edges"),
    end = role_optional(edges, end, "end", "edges"),
    wave = role_optional(edges, wave, "wave", "edges"),
    node_start = role_optional(nodes, node_start, "node_start", "nodes"),
    node_end = role_optional(nodes, node_end, "node_end", "nodes"),
    source_layer = role_optional(edges, source_layer, "source_layer", "edges"),
    target_layer = role_optional(edges, target_layer, "target_layer", "edges"),
    ego = role_optional(nodes, ego, "ego", "nodes"),
    edge_ego = role_optional(edges, edge_ego, "edge_ego", "edges"),
    node_missing = role_optional(nodes, node_missing, "node_missing", "nodes"),
    edge_missing = role_optional(edges, edge_missing, "edge_missing", "edges")
  )
  if (!is.null(nodes)) roles$id <- role_column(nodes, id, "id", "nodes")
  if (!is.null(layers)) roles$layer <- role_column(layers, layer, "layer", "layers")
  if (!is.null(node_spells)) {
    role_column(node_spells, roles$id, "id", "node_spells")
    roles$spell_start <- role_column(node_spells, spell_start, "spell_start", "node_spells")
    roles$spell_end <- role_column(node_spells, spell_end, "spell_end", "node_spells")
  }
  if (!is.null(node_changes)) {
    role_column(node_changes, roles$id, "id", "node_changes")
    roles$change_time <- role_column(node_changes, change_time, "change_time",
                                     "node_changes")
  }
  roles <- c(roles, time_roles(edges, nodes, node_spells, node_changes, roles))
  if (!is.null(coords)) {
    roles <- c(roles, coordinate_roles(nodes, coords, crs))
  } else if (!is.null(crs)) {
    stop("`crs` describes `coords`; give the coordinate columns too.", call. = FALSE)
  }
  declared <- list(
    loops = check_flag_or_null(loops, "loops"),
    multiple = check_flag_or_null(multiple, "multiple"),
    between_modes_only = check_flag_or_null(between_modes_only, "between_modes_only"),
    egos_are_nodes = check_flag_or_null(egos_are_nodes, "egos_are_nodes"),
    couplings_only = check_flag_or_null(couplings_only, "couplings_only"),
    overlapping_spells = check_flag_or_null(overlapping_spells, "overlapping_spells")
  )
  new_net_frame(edges, nodes, directed, roles, declared, graph, add_nodes,
                tables = list(layers = layers, node_spells = node_spells,
                              node_changes = node_changes))
}

# See man/net_frame.Rd.
net_diagnostics <- function(x, threads = 0L) {
  check_net_frame(x)
  threads <- resolve_threads(threads)
  frame_diagnostics(engine_frame(x), threads)
}

# ---- Construction ----------------------------------------------------------

# Checks the two tables.
check_tables <- function(edges, nodes) {
  if (!is.data.frame(edges)) {
    stop(sprintf("`edges` must be a data frame, not %s.", describe_class(edges)),
         call. = FALSE)
  }
  if (!is.null(nodes) && !is.data.frame(nodes)) {
    stop(sprintf("`nodes` must be a data frame or NULL, not %s.", describe_class(nodes)),
         call. = FALSE)
  }
  invisible(NULL)
}

# The network-level attributes as a named list.
check_graph <- function(graph) {
  if (is.null(graph)) graph <- list()
  if (!is.list(graph) && !is.atomic(graph)) {
    stop("`graph` must be a named list of network-level attributes.", call. = FALSE)
  }
  graph <- as.list(graph)
  if (length(graph) && (is.null(names(graph)) || anyNA(names(graph)) ||
                        !all(nzchar(names(graph))))) {
    stop("Every element of `graph` must be named, as in `graph = list(title = \"...\")`.",
         call. = FALSE)
  }
  graph
}

# The declared meaning of the weights: required when there is a weight
# column, NULL otherwise (the engine reports E_WEIGHT_TYPE_MISSING when
# the two disagree; the check here gives the message before any work).
check_weight_type <- function(weight, weight_type) {
  if (is.null(weight)) {
    if (!is.null(weight_type)) {
      check_choice(weight_type, c("strength", "distance"), "weight_type")
    }
    return(NULL)
  }
  if (is.null(weight_type)) {
    stop(paste0("`weight_type` must say what the weights mean: \"strength\" (larger ",
                "brings nodes closer: frequency, intensity, capacity) or \"distance\" ",
                "(larger keeps them apart: cost, length, latency). There is no default."),
         call. = FALSE)
  }
  check_choice(weight_type, c("strength", "distance"), "weight_type")
}

# A role column given by name (never by position: decision 10).
role_column <- function(df, col, role, what) {
  if (!is.character(col) || length(col) != 1L || is.na(col) || !nzchar(col)) {
    stop(sprintf("`%s` must be the name of a column of `%s`.", role, what), call. = FALSE)
  }
  if (!col %in% names(df)) {
    shown <- if (ncol(df)) {
      sprintf(" It has %s: %s.", plural(ncol(df), "column"),
              paste0("`", utils::head(names(df), 10L), "`", collapse = ", "))
    } else {
      " It has no columns."
    }
    stop(sprintf("`%s` names column `%s`, which is not in `%s`.%s", role, col, what, shown),
         call. = FALSE)
  }
  col
}

# A role column given by name, or NULL when the role is not used.
role_optional <- function(df, col, role, what) {
  if (is.null(col)) return(NULL)
  if (is.null(df)) {
    stop(sprintf("`%s` names a column of `%s`, but no `%s` were given.", role, what, what),
         call. = FALSE)
  }
  role_column(df, col, role, what)
}

# Kept for the other constructors: a column by name or position.
resolve_optional <- function(df, col, role, what) {
  if (is.null(col)) return(NULL)
  if (is.null(df)) {
    stop(sprintf("`%s` names a column of `%s`, but no `%s` were given.", role, what, what),
         call. = FALSE)
  }
  resolve_column(df, col, role, what)
}

# How times are stored: the class every time column shares ("numeric",
# "Date", or "POSIXct") and the time zone of date-times. Waves are whole
# numbers; node spells and changes then count waves too.
time_roles <- function(edges, nodes, node_spells, node_changes, roles) {
  cols <- list(
    list(edges, roles$start, "Start times", "edges"),
    list(edges, roles$end, "End times", "edges"),
    list(nodes, roles$node_start, "Node start times", "nodes"),
    list(nodes, roles$node_end, "Node end times", "nodes"),
    list(node_spells, roles$spell_start, "Spell start times", "node_spells"),
    list(node_spells, roles$spell_end, "Spell end times", "node_spells"),
    list(node_changes, roles$change_time, "Change times", "node_changes")
  )
  cols <- Filter(function(c) !is.null(c[[2L]]), cols)
  if (!is.null(roles$wave)) {
    w <- edges[[roles$wave]]
    if (!is.numeric(w) || is.object(w)) {
      stop(sprintf("Waves (column `%s` of `edges`) must be whole numbers, not %s.",
                   roles$wave, describe_class(w)), call. = FALSE)
    }
    if (!is.null(roles$start) || !is.null(roles$end)) {
      stop("Give `wave`, or `start` (with `end`), not both: wave k is the spell [k, k + 1).",
           call. = FALSE)
    }
  }
  if (!is.null(roles$end) && is.null(roles$start)) {
    stop("`end` needs `start`: an edge spell runs from `start` to `end`.", call. = FALSE)
  }
  if (!length(cols)) {
    return(if (is.null(roles$wave)) list() else list(time_class = "numeric"))
  }
  first <- cols[[1L]]
  kind <- time_kind(first[[1L]][[first[[2L]]]])
  if (is.na(kind)) {
    stop(sprintf("%s (column `%s` of `%s`) must be numbers, dates, or date-times, not %s.",
                 first[[3L]], first[[2L]], first[[4L]],
                 describe_class(first[[1L]][[first[[2L]]]])), call. = FALSE)
  }
  if (!is.null(roles$wave) && kind != "numeric") {
    stop("Node spells of a panel count waves, so they must be numbers.", call. = FALSE)
  }
  for (c in cols[-1L]) {
    check_same_kind(c[[1L]][[c[[2L]]]], kind,
                    sprintf("%s (column `%s` of `%s`)", c[[3L]], c[[2L]], c[[4L]]))
  }
  out <- list(time_class = kind)
  if (kind == "POSIXct") {
    z <- attr(as.POSIXct(first[[1L]][[first[[2L]]]]), "tzone")
    out$tz <- if (is.null(z)) "" else z[[1L]]
  }
  out
}

# The coordinate roles: two node columns and what they are.
coordinate_roles <- function(nodes, coords, crs) {
  if (is.null(nodes)) {
    stop("`coords` names columns of `nodes`, but no `nodes` were given.", call. = FALSE)
  }
  if (!is.character(coords) || length(coords) != 2L) {
    stop("`coords` must name two columns of `nodes`: longitude and latitude, or x and y.",
         call. = FALSE)
  }
  if (is.null(crs)) crs <- "lonlat"
  check_choice(crs, c("lonlat", "planar"), "crs")
  list(x = role_column(nodes, coords[[1L]], "coords[1]", "nodes"),
       y = role_column(nodes, coords[[2L]], "coords[2]", "nodes"), crs = crs)
}

# An unbuilt net_frame: the tables, roles, and declarations, with no index
# or cache yet. build_net_frame() validates it and fills the rest in.
new_net_frame <- function(edges, nodes, directed, roles, declared, graph, add_nodes = FALSE,
                          tables = list()) {
  roles <- roles[!vapply(roles, is.null, TRUE)]
  structure(
    list(
      nodes = nodes,
      edges = edges,
      layers = tables$layers,
      node_spells = tables$node_spells,
      node_changes = tables$node_changes,
      directed = isTRUE(directed),
      roles = roles,
      declared = declared[!vapply(declared, is.null, TRUE)],
      graph = graph,
      index = NULL,
      cache = new.env(parent = emptyenv())
    ),
    class = "net_frame",
    add_nodes = add_nodes
  )
}

# The long tables of a net_frame, for new_net_frame().
frame_tables <- function(x) {
  list(layers = x$layers, node_spells = x$node_spells, node_changes = x$node_changes)
}

# Validates and builds a net_frame through the engine: stops with every
# error, reports the notes (unless `quiet`), appends nodes taken from the
# edges, records the edge index, caches the graph handle and the summary,
# and sets the class. `fn` names the constructor in messages.
build_net_frame <- function(x, quiet = FALSE, fn = "net_frame", threads = 1L) {
  # sf node or edge tables record their geometry columns and CRS.
  x <- record_geometry(x)
  f <- engine_frame(x)
  d <- frame_diagnostics(f, threads)
  errors <- d[d$severity == "error", , drop = FALSE]
  if (nrow(errors)) stop(format_diagnostics(errors, fn, x), call. = FALSE)
  if (!quiet) {
    notes <- d[d$severity == "note", , drop = FALSE]
    for (k in seq_len(nrow(notes))) message("Note: ", diagnostic_text(notes[k, ], x))
  }
  ptr <- engine_call(fn, .Call(C_frame_build, f, threads))
  idx <- .Call(C_frame_index, f, as.double(nrow(x$edges)))
  if (length(idx$added_row)) x <- append_nodes(x, idx)
  if (is.null(x$nodes)) {
    nodes <- data.frame(x$edges[[x$roles$source]][0], stringsAsFactors = FALSE)
    names(nodes) <- x$roles$id
    x$nodes <- nodes
  }
  x$index <- list(src = idx$src, tgt = idx$tgt)
  attr(x, "add_nodes") <- NULL
  cache <- x$cache
  rm(list = ls(cache, all.names = TRUE), envir = cache)
  cache$frame <- f
  cache$ptr <- ptr
  cache$n <- nrow(x$nodes)
  cache$m <- nrow(x$edges)
  cache$weights <- if (is.null(x$roles$weight)) NULL else x$edges[[x$roles$weight]]
  cache$graph_rows <- if (is.null(idx$graph_rows)) NULL else idx$graph_rows + 1L
  cache$summary <- frame_summary(f, threads)
  cache$diagnostics <- d
  class(x) <- net_frame_class(x)
  x
}

# Appends the nodes the engine took from the edges (in its order) to the
# node table, creating the table when there was none.
append_nodes <- function(x, idx) {
  src <- x$edges[[x$roles$source]]
  tgt <- x$edges[[x$roles$target]]
  if (is.factor(src) || is.factor(tgt)) {
    src <- as.character(src)
    tgt <- as.character(tgt)
  }
  if (!identical(class(src), class(tgt))) {
    stop(sprintf(paste0(
      "Columns `%s` and `%s` of `edges` hold different types (%s and %s); ",
      "convert them to the same type first."),
      x$roles$source, x$roles$target, describe_class(src), describe_class(tgt)),
      call. = FALSE)
  }
  both <- c(src, tgt)
  m <- length(src)
  ids <- both[idx$added_row + 1L + ifelse(idx$added_target, m, 0L)]
  if (is.null(x$nodes)) {
    nodes <- data.frame(ids, stringsAsFactors = FALSE)
    names(nodes) <- x$roles$id
    x$nodes <- nodes
    return(x)
  }
  if (inherits(x$nodes, "sf")) {
    stop(sprintf(paste0(
      "%s of `edges` name nodes that are not in the sf node table, which cannot be ",
      "extended without their geometry; add those nodes to `nodes` first."),
      rows_text(idx$added_row + 1L)), call. = FALSE)
  }
  extra <- x$nodes[rep(NA_integer_, length(ids)), , drop = FALSE]
  extra[[x$roles$id]] <- ids
  nodes <- rbind(x$nodes, extra)
  rownames(nodes) <- NULL
  x$nodes <- nodes
  x
}

# The class vector a net_frame gets from its roles, in the fixed order of
# DATA_MODEL.md section 20: multilayer when edges carry layers; temporal
# when edges carry times or waves or nodes carry presence spells; spatial
# when nodes carry coordinates or either table carries geometry;
# spatio-temporal when both.
net_frame_class <- function(x) {
  r <- x$roles
  multilayer <- !is.null(r$source_layer) && !is.null(r$target_layer)
  temporal <- !is.null(r$start) || !is.null(r$wave) || !is.null(r$node_start) ||
    !is.null(r$node_end) || !is.null(x$node_spells)
  spatial <- (!is.null(r$x) && !is.null(r$y)) || !is.null(r$node_geometry) ||
    !is.null(r$edge_geometry)
  c(if (multilayer) "multilayer_net_frame",
    if (temporal && spatial) "spatiotemporal_net_frame",
    if (spatial) "spatial_net_frame",
    if (temporal) "temporal_net_frame",
    "net_frame")
}

# Role names of each table, used when tables change.
node_role_names <- c("id", "mode", "node_start", "node_end", "x", "y", "node_missing", "ego")
edge_role_names <- c("source", "target", "weight", "relation", "sign", "start", "end",
                     "duration", "wave", "edge_missing", "source_layer", "target_layer",
                     "edge_ego")

# A net_frame rebuilt from new tables with the same roles, declarations,
# long tables, and network attributes (roles whose columns are gone are
# dropped, so the class follows the columns). `...` overrides roles:
# `weight = NULL` drops the weight role, `directed = FALSE` changes the
# direction, `declared = list(...)` replaces the declarations, and
# `tables = list(...)` replaces long tables. Node spells and changes of
# nodes that are no longer in the node table are dropped.
rebuild_net_frame <- function(x, edges, nodes, ..., quiet = TRUE, fn = "net_frame") {
  changes <- list(...)
  roles <- x$roles
  directed <- x$directed
  declared <- x$declared
  tables <- frame_tables(x)
  for (nm in names(changes)) {
    if (nm == "directed") {
      directed <- changes[[nm]]
    } else if (nm == "declared") {
      declared <- changes[[nm]]
    } else if (nm == "tables") {
      for (t in names(changes$tables)) tables[t] <- list(changes$tables[[t]])
    } else {
      roles[nm] <- list(changes[[nm]])
    }
  }
  if (is.null(roles$weight)) roles$weight_type <- NULL
  for (nm in node_role_names) {
    if (!is.null(roles[[nm]]) && !roles[[nm]] %in% names(nodes)) roles[nm] <- list(NULL)
  }
  for (nm in edge_role_names) {
    if (!is.null(roles[[nm]]) && !roles[[nm]] %in% names(edges)) roles[nm] <- list(NULL)
  }
  if (is.null(roles$x) || is.null(roles$y)) {
    roles$x <- NULL
    roles$y <- NULL
    if (is.character(roles$crs)) roles$crs <- NULL
  }
  if (is.null(roles$source_layer) || is.null(roles$target_layer)) {
    roles$source_layer <- NULL
    roles$target_layer <- NULL
    tables$layers <- NULL
    declared$couplings_only <- NULL
  }
  if (is.null(tables$layers)) roles$layer <- NULL
  if (is.null(roles$start)) {
    roles$end <- NULL
    roles$duration <- NULL
  }
  ids <- as.character(nodes[[roles$id]])
  for (t in c("node_spells", "node_changes")) {
    tb <- tables[[t]]
    if (!is.null(tb) && roles$id %in% names(tb)) {
      keep <- as.character(tb[[roles$id]]) %in% ids
      if (!all(keep)) {
        tb <- tb[keep, , drop = FALSE]
        rownames(tb) <- NULL
      }
      tables[t] <- list(tb)
    }
  }
  if (is.null(tables$node_spells)) {
    roles$spell_start <- NULL
    roles$spell_end <- NULL
    declared$overlapping_spells <- NULL
  }
  if (is.null(tables$node_changes)) roles$change_time <- NULL
  timed <- !is.null(roles$start) || !is.null(roles$wave) || !is.null(roles$node_start) ||
    !is.null(roles$node_end) || !is.null(tables$node_spells) || !is.null(tables$node_changes)
  if (!timed) {
    roles$time_class <- NULL
    roles$tz <- NULL
  }
  if (is.null(roles$mode)) {
    roles$mode_levels <- NULL
    declared$between_modes_only <- NULL
  }
  if (is.null(roles$ego)) {
    roles$edge_ego <- NULL
    declared$egos_are_nodes <- NULL
  }
  roles$node_geometry <- NULL
  roles$edge_geometry <- NULL
  if (inherits(roles$crs, "crs")) roles$crs <- NULL
  out <- new_net_frame(edges, nodes, directed, roles, declared, x$graph, tables = tables)
  build_net_frame(out, quiet = quiet, fn = fn)
}

# ---- The engine frame -----------------------------------------------------------

# Codes of the roles in the C interface (DATA_MODEL.md section 19).
node_role_codes <- c(id = 0L, mode = 1L, node_start = 2L, node_end = 3L, x = 4L, y = 5L,
                     node_missing = 6L, ego = 7L)
edge_role_codes <- c(source = 0L, target = 1L, weight = 2L, relation = 3L, sign = 4L,
                     start = 5L, end = 6L, wave = 7L, edge_missing = 8L, source_layer = 9L,
                     target_layer = 10L, edge_ego = 11L)
time_role_names <- c("node_start", "node_end", "start", "end")

# The engine frame for `x`: the role columns of every table as typed
# buffers, the roles, the declarations, and the geometry flags. Built on
# demand; the handle is cached by build_net_frame() and engine_graph().
engine_frame <- function(x) {
  roles <- x$roles
  nodes <- x$nodes
  edges <- x$edges
  cols <- list(list(), list(), list(), list(), list())
  tables <- integer()
  codes <- integer()
  names_ <- character()
  add <- function(table, role, col, values) {
    tables <<- c(tables, table)
    codes <<- c(codes, role)
    names_ <<- c(names_, col)
    cols[[table + 1L]][[col]] <<- values
  }
  column_of <- function(df, col, role, what) {
    if (!col %in% names(df)) {
      stop(sprintf("Column `%s` (the %s role) is no longer in `x$%s`.", col, role, what),
           call. = FALSE)
    }
    df[[col]]
  }
  for (nm in if (is.null(nodes)) character() else names(node_role_codes)) {
    col <- roles[[nm]]
    if (is.null(col)) next
    v <- column_of(nodes, col, nm, "nodes")
    values <- if (nm %in% time_role_names) {
      time_values(v, col, "nodes")
    } else {
      engine_column(v, col, "nodes", identifier = nm %in% c("id", "ego"),
                    discrete = nm == "mode")
    }
    add(0L, node_role_codes[[nm]], col, values)
  }
  for (nm in names(edge_role_codes)) {
    col <- roles[[nm]]
    if (is.null(col)) next
    v <- column_of(edges, col, nm, "edges")
    values <- if (nm %in% time_role_names) {
      time_values(v, col, "edges")
    } else {
      engine_column(v, col, "edges",
                    identifier = nm %in% c("source", "target", "edge_ego"),
                    discrete = nm %in% c("relation", "sign", "wave", "source_layer",
                                         "target_layer"))
    }
    add(1L, edge_role_codes[[nm]], col, values)
  }
  # A duration column stands in for the end times: start plus duration.
  if (!is.null(roles$duration) && is.null(roles$end)) {
    col <- roles$duration
    d <- column_of(edges, col, "duration", "edges")
    d <- if (inherits(d, "difftime")) as.double(d, units = "secs") else as.double(d) * time_scale(x)
    s <- cols[[2L]][[roles$start]]
    add(1L, edge_role_codes[["end"]], col, s + d)
  }
  if (!is.null(x$layers)) {
    add(2L, 0L, roles$layer, engine_column(column_of(x$layers, roles$layer, "layer", "layers"),
                                           roles$layer, "layers", discrete = TRUE))
  }
  if (!is.null(x$node_spells)) {
    sp <- x$node_spells
    add(3L, 0L, roles$id, engine_column(column_of(sp, roles$id, "id", "node_spells"),
                                        roles$id, "node_spells", identifier = TRUE))
    add(3L, 1L, roles$spell_start,
        time_values(column_of(sp, roles$spell_start, "spell_start", "node_spells"),
                    roles$spell_start, "node_spells"))
    add(3L, 2L, roles$spell_end,
        time_values(column_of(sp, roles$spell_end, "spell_end", "node_spells"),
                    roles$spell_end, "node_spells"))
  }
  if (!is.null(x$node_changes)) {
    ch <- x$node_changes
    add(4L, 0L, roles$id, engine_column(column_of(ch, roles$id, "id", "node_changes"),
                                        roles$id, "node_changes", identifier = TRUE))
    add(4L, 1L, roles$change_time,
        time_values(column_of(ch, roles$change_time, "change_time", "node_changes"),
                    roles$change_time, "node_changes"))
  }
  weight_type <- if (is.null(roles$weight_type)) -1L else {
    if (roles$weight_type == "strength") 1L else 0L
  }
  code_of <- function(v) if (is.null(v)) -1L else as.integer(isTRUE(v))
  dec <- x$declared
  declarations <- c(code_of(dec$loops), code_of(dec$multiple),
                    code_of(dec$between_modes_only), as.integer(isTRUE(attr(x, "add_nodes"))),
                    code_of(dec$egos_are_nodes), code_of(dec$couplings_only),
                    code_of(dec$overlapping_spells))
  geometry <- NULL
  if (!is.null(roles$node_geometry) || !is.null(roles$edge_geometry)) {
    geometry <- list(!is.null(roles$node_geometry), !is.null(roles$edge_geometry),
                     crs_text(roles$crs))
  }
  coordinates <- if (!is.null(roles$x) && !is.null(roles$y)) {
    if (identical(crs_kind(x), "lonlat")) 1L else 2L
  }
  extra <- list(if (is.null(x$layers)) NULL else cols[[3L]],
                if (is.null(x$node_spells)) NULL else cols[[4L]],
                if (is.null(x$node_changes)) NULL else cols[[5L]])
  engine_call("net_frame",
    .Call(C_frame_new, x$directed, if (is.null(nodes)) 0 else as.double(nrow(nodes)),
          if (is.null(nodes)) NULL else cols[[1L]], cols[[2L]], extra,
          list(tables, codes, names_), weight_type, declarations, geometry, coordinates))
}

# Times as seconds for the engine, in the network's time class; a column
# of another type goes through as it is so the engine reports it.
time_values <- function(v, col, table) {
  k <- time_kind(v)
  if (is.na(k)) return(engine_column(v, col, table))
  to_seconds(v, k)
}

# A column as the engine takes it: integer, double, logical, or character.
# Factors and 64-bit integers (bit64) become text, dates and date-times
# seconds, and logical identifiers integers.
engine_column <- function(v, col, table, identifier = FALSE, discrete = FALSE) {
  if (inherits(v, "POSIXt")) return(as.double(unclass(as.POSIXct(v))))
  if (!is.atomic(v) || is.matrix(v) || is.null(v)) {
    stop(sprintf("Column `%s` of `%s` must be an atomic vector, not %s.", col, table,
                 describe_class(v)), call. = FALSE)
  }
  if (is.factor(v)) return(as.character(v))
  if (inherits(v, "integer64")) return(as.character(v))
  if (inherits(v, "Date")) return(as.double(unclass(v)) * 86400)
  if (inherits(v, "difftime")) return(as.double(v, units = "secs"))
  if (is.logical(v)) return(if (identifier) as.integer(v) else as.vector(v))
  if (is.integer(v)) return(as.vector(unclass(v)))
  if (is.double(v)) {
    v <- as.vector(unclass(v))
    # A whole-numbered double column is discrete enough for a mode,
    # relation, sign, wave, or layer.
    whole <- discrete && all(is.na(v) | (is.finite(v) & v == round(v) & abs(v) < 2^31))
    return(if (whole) as.integer(v) else v)
  }
  if (is.character(v)) return(as.vector(v))
  stop(sprintf(paste0("Column `%s` of `%s` must hold numbers, text, logical values, ",
                      "factors, or dates, not %s."), col, table, describe_class(v)),
       call. = FALSE)
}

# The coordinate reference system as text for the engine: the strings
# spatial networks record, or an sf crs by its input string.
crs_text <- function(crs) {
  if (is.null(crs)) return(NULL)
  if (is.character(crs)) return(crs)
  if (inherits(crs, "crs")) {
    if (is.na(crs)) return(NULL)
    s <- tryCatch(crs$input, error = function(e) NULL)
    if (is.null(s) || is.na(s)) s <- tryCatch(crs$wkt, error = function(e) NULL)
    if (is.null(s) || is.na(s)) return(NULL)
    return(as.character(s))
  }
  as.character(crs)
}

# ---- Diagnostics ------------------------------------------------------------------

table_names <- c("nodes", "edges", "layers", "node_spells", "node_changes")

# The diagnostics of frame `f` as a data frame: code, severity ("error" or
# "note"), table ("nodes", "edges", "layers", "node_spells",
# "node_changes", or NA for the network), message, count (rows affected in
# all), and rows (a list column of the first 1-based rows listed).
frame_diagnostics <- function(f, threads) {
  v <- .Call(C_frame_validate, f, as.integer(threads))
  tab <- rep(NA_character_, length(v$table))
  ok <- v$table >= 0L & v$table < length(table_names)
  tab[ok] <- table_names[v$table[ok] + 1L]
  out <- data.frame(
    code = v$code,
    severity = ifelse(v$severity == 0L, "error", "note"),
    table = tab,
    message = v$message,
    count = as.integer(v$count),
    stringsAsFactors = FALSE
  )
  out$rows <- I(lapply(v$rows, function(r) unique(as.integer(r) + 1L)))
  attr(out, "errors") <- v$errors
  attr(out, "notes") <- v$notes
  out
}

# One diagnostic as a sentence: the engine's message, the rows of the
# table it concerns, and the offending identifiers where that helps.
diagnostic_text <- function(d, x = NULL) {
  msg <- d$message
  rows <- d$rows[[1L]]
  if (length(rows) && !is.na(d$table)) {
    msg <- sprintf("%s; %s of `%s`", msg, rows_text(rows, d$count), d$table)
  }
  detail <- if (is.null(x)) NULL else diagnostic_detail(d, rows, x)
  if (!is.null(detail)) msg <- sprintf("%s (%s)", msg, detail)
  paste0(msg, ".")
}

# The identifiers behind an identifier diagnostic, up to five, and the
# argument that changes the rule where there is one.
diagnostic_detail <- function(d, rows, x) {
  shown <- function(v) {
    v <- unique(as.character(v))
    if (!length(v)) return(NULL)
    more <- if (length(v) > 5L) sprintf(", and %d more", length(v) - 5L) else ""
    paste0(paste(utils::head(v, 5L), collapse = ", "), more)
  }
  if (d$code == "E_UNKNOWN_NODE" && identical(d$table, "edges") && !is.null(x$nodes)) {
    vals <- c(x$edges[[x$roles$source]][rows], x$edges[[x$roles$target]][rows])
    ids <- x$nodes[[x$roles$id]]
    unknown <- shown(vals[!as.character(vals) %in% as.character(ids)])
    if (!is.null(unknown)) {
      return(sprintf("not in `nodes`: %s; add them to `nodes`, or set `add_nodes = TRUE`",
                     unknown))
    }
  }
  if (d$code == "E_DUPLICATE_ID" && identical(d$table, "nodes")) {
    dup <- shown(x$nodes[[x$roles$id]][rows])
    if (!is.null(dup)) return(paste("repeated:", dup))
  }
  NULL
}

# The error message of a constructor: one sentence for one problem, a list
# for several. Never partial: every error is listed.
format_diagnostics <- function(errors, fn, x = NULL) {
  lines <- vapply(seq_len(nrow(errors)), function(k) diagnostic_text(errors[k, ], x), "")
  if (length(lines) == 1L) return(sprintf("%s(): %s", fn, lines))
  sprintf("%s(): the network is invalid (%d problems):\n%s", fn, length(lines),
          paste0("  - ", lines, collapse = "\n"))
}

# ---- The summary ------------------------------------------------------------------

# The engine's derived properties of a built frame as a list (see
# net_frame_summary in netcapi.h for the fields).
frame_summary <- function(f, threads) {
  s <- .Call(C_frame_summary, f, as.integer(threads))
  v <- s$values
  time_kinds <- c("static", "contacts", "spells", "waves")
  flag <- function(k) if (is.na(v[[k]]) || v[[k]] < 0) NULL else v[[k]] != 0
  count <- function(k) if (is.na(v[[k]])) 0 else v[[k]]
  list(
    nodes = v[[1L]], edges = v[[2L]], directed = v[[3L]] != 0,
    modes = s$modes, mode_counts = s$mode_counts,
    between_modes_only = flag(5L),
    weighted = v[[6L]] != 0,
    weight = s$weight,
    weight_type = if (v[[7L]] < 0) NULL else if (v[[7L]] == 1) "strength" else "distance",
    weight_min = v[[8L]], weight_max = v[[9L]],
    relations = s$relations, relation_counts = s$relation_counts,
    signed = v[[11L]] != 0, negative = count(12L),
    time = time_kinds[[v[[13L]] + 1L]], time_start = v[[14L]], time_end = v[[15L]],
    node_spells = c("none", "roles", "table")[[count(16L) + 1L]],
    node_geometry = v[[17L]] != 0, edge_geometry = v[[18L]] != 0,
    coordinates = if (is.na(v[[19L]]) || v[[19L]] < 0) NULL else {
      c("unknown", "lonlat", "planar")[[v[[19L]] + 1L]]
    },
    crs = s$crs,
    loops = count(20L), parallel = count(21L), isolates = count(22L),
    components = count(23L), density = v[[24L]],
    missing_nodes = count(25L), missing_edges = count(26L),
    layers = s$layers, layer_counts = s$layer_counts, layer_table = isTRUE(v[[28L]] != 0),
    intralayer = count(29L), coupling = count(30L), interlayer = count(31L),
    couplings_only = flag(32L),
    egos = s$egos, ego_nodes = s$ego_nodes, ego_edges = s$ego_edges,
    egos_as_nodes = count(34L), egos_are_nodes = flag(35L),
    spell_rows = count(36L), change_rows = count(37L)
  )
}

# The cached summary of a net_frame, recomputed when the handle is stale.
frame_summary_cached <- function(x) {
  engine_graph(x, missing_ok = TRUE)
  s <- x$cache$summary
  if (is.null(s)) {
    s <- frame_summary(x$cache$frame, 1L)
    x$cache$summary <- s
  }
  s
}

# ---- Parallel edges ---------------------------------------------------------------

# Edge multiplicity, simplicity, and simplification (man/net_simplify.Rd).
net_multiplicity <- function(x) {
  check_net_frame(x)
  engine_graph(x, missing_ok = TRUE)
  engine_call("net_multiplicity",
              .Call(C_frame_multiplicity, x$cache$frame, as.double(nrow(x$edges))))
}

# See man/net_simplify.Rd.
net_is_simple <- function(x) {
  check_net_frame(x)
  engine_graph(x, missing_ok = TRUE)
  engine_call("net_is_simple", .Call(C_frame_is_simple, x$cache$frame))
}

# See man/net_simplify.Rd.
net_simplify <- function(x, combine = c("sum", "mean", "min", "max", "first", "count"),
                         remove_loops = FALSE, attributes = "first") {
  check_net_frame(x)
  combine <- match.arg(combine)
  remove_loops <- check_flag(remove_loops, "remove_loops")
  attributes <- check_attribute_rule(attributes)
  engine_graph(x, missing_ok = TRUE)
  weighted <- !is.null(x$roles$weight)
  code <- match(combine, c("sum", "mean", "min", "max", "first", "count")) - 1L
  r <- engine_call("net_simplify",
                   .Call(C_frame_simplify, x$cache$frame, as.double(nrow(x$edges)), code,
                         remove_loops, weighted || combine == "count"))
  group <- r$group
  keep <- group >= 0L & !duplicated(group)
  edges <- x$edges[keep, , drop = FALSE]
  edges <- aggregate_edge_attributes(x, edges, group, keep, attributes)
  weight <- x$roles$weight
  if (weighted) {
    edges[[weight]] <- r$weight[keep]
  } else if (combine == "count") {
    weight <- if ("count" %in% names(edges)) "count_" else "count"
    edges[[weight]] <- r$weight[keep]
  }
  rownames(edges) <- NULL
  declared <- x$declared
  declared$multiple <- FALSE
  if (remove_loops) declared$loops <- FALSE
  rebuild_net_frame(x, edges, x$nodes, weight = weight,
                    weight_type = if (weighted) x$roles$weight_type else "strength",
                    declared = declared, fn = "net_simplify")
}

# The rule for edge attributes other than the roles when parallel edges
# are merged: "first" keeps the first row's value, "concat" pastes the
# distinct values, or a function of one vector that returns one value.
check_attribute_rule <- function(attributes) {
  if (is.function(attributes)) return(attributes)
  check_choice(attributes, c("first", "concat"), "attributes")
}

# Aggregates the non-role edge columns of each group of parallel edges with
# the rule (the kept row is each group's first row, so "first" is a no-op).
aggregate_edge_attributes <- function(x, edges, group, keep, rule) {
  if (identical(rule, "first")) return(edges)
  roles <- unlist(x$roles[intersect(names(x$roles), edge_role_names)], use.names = FALSE)
  cols <- setdiff(names(x$edges), c(roles, x$roles$edge_geometry))
  if (!length(cols)) return(edges)
  g <- group[group >= 0L]
  rows <- which(group >= 0L)
  kept_groups <- group[keep]
  fun <- if (is.function(rule)) rule else function(v) {
    v <- unique(as.character(v[!is.na(v)]))
    if (!length(v)) NA_character_ else paste(v, collapse = "; ")
  }
  for (col in cols) {
    v <- x$edges[[col]]
    if (is.list(v)) next
    vals <- lapply(split(rows, factor(g, levels = kept_groups)), function(r) fun(v[r]))
    edges[[col]] <- unlist(vals, use.names = FALSE)
  }
  edges
}
