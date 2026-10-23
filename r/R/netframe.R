# The netframe class: the constructor, which hands the role columns of the
# node and edge tables to the engine for validation and construction, the
# diagnostics it reports, the subclasses, and the print card. The argument
# checks are in utils-checks.R, the engine handle in utils-engine.R, the
# print card in print-netframe.R, and summary() in descriptive.R.
#
# A netframe is a list: `nodes` and `edges` (the two tables, untouched
# except for nodes appended from the edges), `directed`, `roles` (which
# columns carry structural meaning: id, source, target, weight and
# weight_type, mode, relation, sign, start, end, duration, wave,
# node_start, node_end, x, y, node_missing, edge_missing, plus the
# geometry columns, the coordinate reference system, and how times are
# stored), `declared` (loops, multiple, between_modes_only: TRUE, FALSE,
# or NULL for not declared), `graph` (network-level attributes), `index`
# (0-based node positions of every edge's endpoints, from the engine), and
# `cache` (an environment holding the engine handles and the summary).

# Create a network from node and edge data frames (man/net_frame.Rd).
net_frame <- function(edges, nodes = NULL, directed = TRUE, source = "source",
                      target = "target", id = "id", weight = NULL, weight_type = NULL,
                      mode = NULL, relation = NULL, sign = NULL, loops = NULL,
                      multiple = NULL, add_nodes = FALSE, quiet = FALSE, graph = list()) {
  check_tables(edges, nodes)
  directed <- check_flag(directed, "directed")
  add_nodes <- check_flag(add_nodes, "add_nodes")
  quiet <- check_flag(quiet, "quiet")
  graph <- check_graph(graph)
  roles <- list(
    id = check_string(id, "id"),
    source = resolve_column(edges, source, "source", "edges"),
    target = resolve_column(edges, target, "target", "edges"),
    weight = resolve_optional(edges, weight, "weight", "edges"),
    weight_type = check_weight_type(weight, weight_type),
    mode = resolve_optional(nodes, mode, "mode", "nodes"),
    relation = resolve_optional(edges, relation, "relation", "edges"),
    sign = resolve_optional(edges, sign, "sign", "edges")
  )
  if (!is.null(nodes)) roles$id <- resolve_column(nodes, id, "id", "nodes")
  declared <- list(loops = check_flag_or_null(loops, "loops"),
                   multiple = check_flag_or_null(multiple, "multiple"),
                   between_modes_only = NULL)
  x <- new_netframe(edges, nodes, directed, roles, declared, graph, add_nodes)
  build_netframe(x, quiet, "net_frame")
}

# See man/net_frame.Rd.
net_diagnostics <- function(x, threads = 0L) {
  check_netframe(x)
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

# A role column given by name, or NULL when the role is not used.
resolve_optional <- function(df, col, role, what) {
  if (is.null(col)) return(NULL)
  if (is.null(df)) {
    stop(sprintf("`%s` names a column of `%s`, but no `%s` were given.", role, what, what),
         call. = FALSE)
  }
  resolve_column(df, col, role, what)
}

# An unbuilt netframe: the tables, roles, and declarations, with no index
# or cache yet. build_netframe() validates it and fills the rest in.
new_netframe <- function(edges, nodes, directed, roles, declared, graph, add_nodes = FALSE) {
  roles <- roles[!vapply(roles, is.null, TRUE)]
  structure(
    list(
      nodes = nodes,
      edges = edges,
      directed = isTRUE(directed),
      roles = roles,
      declared = declared[!vapply(declared, is.null, TRUE)],
      graph = graph,
      index = NULL,
      cache = new.env(parent = emptyenv())
    ),
    class = "netframe",
    add_nodes = add_nodes
  )
}

# Validates and builds a netframe through the engine: stops with every
# error, reports the notes (unless `quiet`), appends nodes taken from the
# edges, records the edge index, caches the graph handle and the summary,
# and sets the class. `fn` names the constructor in messages.
build_netframe <- function(x, quiet = FALSE, fn = "net_frame", threads = 1L) {
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
  cache$summary <- frame_summary(f, threads)
  cache$diagnostics <- d
  class(x) <- netframe_class(x)
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

# The class vector a netframe gets from its roles: temporal when edges
# carry times or waves, spatial when nodes carry coordinates or either
# table carries geometry, both when both.
netframe_class <- function(x) {
  r <- x$roles
  temporal <- !is.null(r$start) || !is.null(r$wave)
  spatial <- (!is.null(r$x) && !is.null(r$y)) || !is.null(r$node_geometry) ||
    !is.null(r$edge_geometry)
  if (temporal && spatial) {
    c("spatiotemporal_netframe", "spatial_netframe", "temporal_netframe", "netframe")
  } else if (temporal) {
    c("temporal_netframe", "netframe")
  } else if (spatial) {
    c("spatial_netframe", "netframe")
  } else {
    "netframe"
  }
}

# A netframe rebuilt from new tables with the same roles, declarations,
# and network attributes (roles whose columns are gone are dropped, so the
# class follows the columns). `...` overrides roles: `weight = NULL` drops
# the weight role, `directed = FALSE` changes the direction.
rebuild_netframe <- function(x, edges, nodes, ..., quiet = TRUE, fn = "net_frame") {
  changes <- list(...)
  roles <- x$roles
  directed <- x$directed
  declared <- x$declared
  for (nm in names(changes)) {
    if (nm == "directed") {
      directed <- changes[[nm]]
    } else if (nm == "declared") {
      declared <- changes[[nm]]
    } else {
      roles[nm] <- list(changes[[nm]])
    }
  }
  if (is.null(roles$weight)) roles$weight_type <- NULL
  node_roles <- c("id", "mode", "node_start", "node_end", "x", "y", "node_missing")
  edge_roles <- c("source", "target", "weight", "relation", "sign", "start", "end",
                  "duration", "wave", "edge_missing")
  for (nm in node_roles) {
    if (!is.null(roles[[nm]]) && !roles[[nm]] %in% names(nodes)) roles[nm] <- list(NULL)
  }
  for (nm in edge_roles) {
    if (!is.null(roles[[nm]]) && !roles[[nm]] %in% names(edges)) roles[nm] <- list(NULL)
  }
  if (is.null(roles$x) || is.null(roles$y)) {
    roles$x <- NULL
    roles$y <- NULL
    if (is.character(roles$crs)) roles$crs <- NULL
  }
  if (is.null(roles$start) && is.null(roles$wave)) {
    roles$end <- NULL
    roles$duration <- NULL
    roles$node_start <- NULL
    roles$node_end <- NULL
    roles$time_class <- NULL
    roles$tz <- NULL
  }
  if (is.null(roles$mode)) {
    roles$mode_levels <- NULL
    declared$between_modes_only <- NULL
  }
  roles$node_geometry <- NULL
  roles$edge_geometry <- NULL
  if (inherits(roles$crs, "crs")) roles$crs <- NULL
  out <- new_netframe(edges, nodes, directed, roles, declared, x$graph)
  build_netframe(out, quiet = quiet, fn = fn)
}

# ---- The engine frame -----------------------------------------------------------

# Codes of the roles in the C interface.
node_role_codes <- c(id = 0L, mode = 1L, node_start = 2L, node_end = 3L, x = 4L, y = 5L,
                     node_missing = 6L)
edge_role_codes <- c(source = 0L, target = 1L, weight = 2L, relation = 3L, sign = 4L,
                     start = 5L, end = 6L, wave = 7L, edge_missing = 8L)

# The engine frame for `x`: the role columns of both tables as typed
# buffers, the roles, the declarations, and the geometry flags. Built on
# demand; the handle is cached by build_netframe() and engine_graph().
engine_frame <- function(x) {
  roles <- x$roles
  nodes <- x$nodes
  edges <- x$edges
  node_cols <- list()
  edge_cols <- list()
  tables <- integer()
  codes <- integer()
  names_ <- character()
  add <- function(table, role, col, values) {
    tables <<- c(tables, table)
    codes <<- c(codes, role)
    names_ <<- c(names_, col)
    if (table == 0L) node_cols[[col]] <<- values else edge_cols[[col]] <<- values
  }
  for (nm in if (is.null(nodes)) character() else names(node_role_codes)) {
    col <- roles[[nm]]
    if (is.null(col)) next
    if (!col %in% names(nodes)) {
      stop(sprintf("Column `%s` (the %s role) is no longer in `x$nodes`.", col, nm),
           call. = FALSE)
    }
    v <- nodes[[col]]
    values <- if (nm %in% c("node_start", "node_end")) {
      time_values(v, col, "nodes")
    } else {
      engine_column(v, col, "nodes", identifier = nm == "id", discrete = nm == "mode")
    }
    add(0L, node_role_codes[[nm]], col, values)
  }
  for (nm in names(edge_role_codes)) {
    col <- roles[[nm]]
    if (is.null(col)) next
    if (!col %in% names(edges)) {
      stop(sprintf("Column `%s` (the %s role) is no longer in `x$edges`.", col, nm),
           call. = FALSE)
    }
    v <- edges[[col]]
    values <- if (nm %in% c("start", "end")) {
      time_values(v, col, "edges")
    } else {
      engine_column(v, col, "edges", identifier = nm %in% c("source", "target"),
                    discrete = nm %in% c("relation", "sign"))
    }
    add(1L, edge_role_codes[[nm]], col, values)
  }
  # A duration column stands in for the end times: start plus duration.
  if (!is.null(roles$duration) && is.null(roles$end)) {
    col <- roles$duration
    if (!col %in% names(edges)) {
      stop(sprintf("Column `%s` (the duration role) is no longer in `x$edges`.", col),
           call. = FALSE)
    }
    d <- edges[[col]]
    d <- if (inherits(d, "difftime")) as.double(d, units = "secs") else as.double(d) * time_scale(x)
    s <- edge_cols[[roles$start]]
    add(1L, edge_role_codes[["end"]], col, s + d)
  }
  weight_type <- if (is.null(roles$weight_type)) -1L else {
    if (roles$weight_type == "strength") 1L else 0L
  }
  code_of <- function(v) if (is.null(v)) -1L else as.integer(isTRUE(v))
  declarations <- c(code_of(x$declared$loops), code_of(x$declared$multiple),
                    code_of(x$declared$between_modes_only),
                    as.integer(isTRUE(attr(x, "add_nodes"))))
  geometry <- NULL
  if (!is.null(roles$node_geometry) || !is.null(roles$edge_geometry)) {
    geometry <- list(!is.null(roles$node_geometry), !is.null(roles$edge_geometry),
                     crs_text(roles$crs))
  }
  coordinates <- if (!is.null(roles$x) && !is.null(roles$y)) {
    if (identical(crs_kind(x), "lonlat")) 1L else 2L
  }
  engine_call("net_frame",
    .Call(C_frame_new, x$directed, if (is.null(nodes)) 0 else as.double(nrow(nodes)),
          if (is.null(nodes)) NULL else node_cols, edge_cols,
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
    # relation, or sign.
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

# The diagnostics of frame `f` as a data frame: code, severity ("error" or
# "note"), table ("nodes", "edges", or NA for the network), message, count
# (rows affected in all), and rows (a list column of the first 1-based
# rows listed).
frame_diagnostics <- function(f, threads) {
  v <- .Call(C_frame_validate, f, as.integer(threads))
  out <- data.frame(
    code = v$code,
    severity = ifelse(v$severity == 0L, "error", "note"),
    table = ifelse(v$table == 0L, "nodes", ifelse(v$table == 1L, "edges", NA_character_)),
    message = v$message,
    count = as.integer(v$count),
    stringsAsFactors = FALSE
  )
  out$rows <- I(lapply(v$rows, function(r) as.integer(r) + 1L))
  attr(out, "errors") <- v$errors
  attr(out, "notes") <- v$notes
  out
}

# One diagnostic as a sentence: the engine's message, the rows of the
# table it concerns, and the offending identifiers where that helps.
diagnostic_text <- function(d, x = NULL) {
  msg <- d$message
  rows <- d$rows[[1L]]
  if (length(rows)) {
    msg <- sprintf("%s; %s of `%s`", msg, rows_text(rows, d$count), d$table)
  }
  detail <- if (is.null(x)) NULL else diagnostic_detail(d, rows, x)
  if (!is.null(detail)) msg <- sprintf("%s (%s)", msg, detail)
  paste0(msg, ".")
}

# The identifiers behind an identifier diagnostic, up to five.
diagnostic_detail <- function(d, rows, x) {
  shown <- function(v) {
    v <- unique(as.character(v))
    if (!length(v)) return(NULL)
    more <- if (length(v) > 5L) sprintf(", and %d more", length(v) - 5L) else ""
    paste0(paste(utils::head(v, 5L), collapse = ", "), more)
  }
  if (d$code == "E_UNKNOWN_NODE" && !is.null(x$nodes)) {
    vals <- c(x$edges[[x$roles$source]][rows], x$edges[[x$roles$target]][rows])
    ids <- x$nodes[[x$roles$id]]
    unknown <- shown(vals[!as.character(vals) %in% as.character(ids)])
    if (!is.null(unknown)) return(paste("not in `nodes`:", unknown))
  }
  if (d$code == "E_DUPLICATE_ID") {
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
  list(
    nodes = v[[1L]], edges = v[[2L]], directed = v[[3L]] != 0,
    modes = s$modes, mode_counts = s$mode_counts,
    between_modes_only = if (v[[5L]] < 0) NULL else v[[5L]] != 0,
    weighted = v[[6L]] != 0,
    weight = s$weight,
    weight_type = if (v[[7L]] < 0) NULL else if (v[[7L]] == 1) "strength" else "distance",
    weight_min = v[[8L]], weight_max = v[[9L]],
    relations = s$relations, relation_counts = s$relation_counts,
    signed = v[[11L]] != 0, negative = v[[12L]],
    time = time_kinds[[v[[13L]] + 1L]], time_start = v[[14L]], time_end = v[[15L]],
    node_spells = v[[16L]] != 0,
    node_geometry = v[[17L]] != 0, edge_geometry = v[[18L]] != 0,
    coordinates = if (v[[19L]] < 0) NULL else c("unknown", "lonlat", "planar")[[v[[19L]] + 1L]],
    crs = s$crs,
    loops = v[[20L]], parallel = v[[21L]], isolates = v[[22L]], components = v[[23L]],
    density = v[[24L]], missing_nodes = v[[25L]], missing_edges = v[[26L]]
  )
}

# The cached summary of a netframe, recomputed when the handle is stale.
netframe_summary <- function(x) {
  engine_graph(x)
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
  check_netframe(x)
  engine_graph(x)
  engine_call("net_multiplicity",
              .Call(C_frame_multiplicity, x$cache$frame, as.double(nrow(x$edges))))
}

# See man/net_simplify.Rd.
net_is_simple <- function(x) {
  check_netframe(x)
  engine_graph(x)
  engine_call("net_is_simple", .Call(C_frame_is_simple, x$cache$frame))
}

# See man/net_simplify.Rd.
net_simplify <- function(x, combine = c("sum", "mean", "min", "max", "first", "count"),
                         remove_loops = FALSE) {
  check_netframe(x)
  combine <- match.arg(combine)
  remove_loops <- check_flag(remove_loops, "remove_loops")
  engine_graph(x)
  weighted <- !is.null(x$roles$weight)
  code <- match(combine, c("sum", "mean", "min", "max", "first", "count")) - 1L
  r <- engine_call("net_simplify",
                   .Call(C_frame_simplify, x$cache$frame, as.double(nrow(x$edges)), code,
                         remove_loops, weighted || combine == "count"))
  group <- r$group
  keep <- group >= 0L & !duplicated(group)
  edges <- x$edges[keep, , drop = FALSE]
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
  rebuild_netframe(x, edges, x$nodes, weight = weight,
                   weight_type = if (weighted) x$roles$weight_type else "strength",
                   declared = declared, fn = "net_simplify")
}
