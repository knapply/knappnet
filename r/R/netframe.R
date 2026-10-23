#' Create a network from node and edge data frames
#'
#' A network is a pair of data frames: `nodes`, with one row per node, and
#' `edges`, with one row per edge. Any other columns are attributes and are
#' carried along untouched. Column roles say which columns hold node
#' identifiers, edge endpoints, and (optionally) edge weights.
#'
#' @param edges A data frame of edges.
#' @param nodes A data frame of nodes, or `NULL` to derive the node set from
#'   the edge endpoints in order of first appearance.
#' @param directed Whether edges are directed: `TRUE` or `FALSE`.
#' @param from,to Columns of `edges` holding the endpoints, by name or
#'   position.
#' @param id Column of `nodes` holding node identifiers, by name or position.
#' @param weight Optional column of `edges` holding edge weights. Without
#'   one, the network is binary: an edge either exists or it does not, and
#'   algorithms use their faster unweighted forms.
#' @param weight_type What the weights mean. `"distance"` (the default, and
#'   the convention of igraph and NetworkX) means larger weights keep nodes
#'   farther apart, as with travel time, cost, or latency. `"strength"` means
#'   larger weights bring nodes closer, as with contact frequency or tie
#'   intensity; path-based measures then use 1/weight as the length (Newman
#'   2001; Opsahl, Agneessens, and Skvoretz 2010). Measures that add up
#'   weights, such as strength, use them as given either way.
#' @param graph A named list of network-level attributes (for example a
#'   title or data source), carried through conversions and file formats
#'   that support them.
#'
#' @section Missing values and invalid input:
#' * Edge endpoints and node identifiers must not be missing (`NA`) or empty
#'   strings; the error names the rows. Identifiers may be character,
#'   factor (matched by their labels), integer, double, or logical, or a
#'   classed vector such as `Date`; `integer64` identifiers (bit64) are
#'   matched as text. Without `nodes`, classed `from` and `to` columns must
#'   have the same class.
#' * Node identifiers must be unique, and every endpoint must appear in
#'   `nodes`.
#' * Networks without edges, and isolated nodes, are allowed.
#' * Self-loops and parallel edges are allowed and kept, and `print()`
#'   reports how many there are. Each measure documents how it counts them;
#'   `convention` selects the treatment of other tools.
#' * Weights must be numeric. Missing, `NaN`, and infinite weights are
#'   accepted here, so the data can still be stored, converted, and
#'   written to files, and `print()` reports them; but every measure stops
#'   with an error naming the rows, even measures that ignore weights,
#'   until they are removed or replaced. Zero and negative weights are
#'   accepted and checked by each measure that uses weights: path-based
#'   measures (betweenness, closeness, harmonic centrality) and structural
#'   holes need them strictly positive; modularity, Louvain, and
#'   eigenvector centrality need them non-negative; strength, structural
#'   equivalence, and CONCOR take any finite value. Errors name the row of
#'   `x$edges` that holds the offending weight.
#' * Arguments are checked before any computation: `x` must be a netframe,
#'   logical options `TRUE` or `FALSE`, `threads` a non-negative whole
#'   number, and choices one of the listed values. Errors name the argument.
#' * The edge index is computed once, here. Add columns to `x$nodes` and
#'   `x$edges` freely, but to add, remove, or reorder rows build a new
#'   netframe. Measures stop when the number of edges has changed or an
#'   edge points past the last node; other changes to rows, such as
#'   reordering, cannot be detected and would misalign results.
#'
#' @section Spatial tables:
#' `nodes` and `edges` may be sf data frames (simple features). Their
#' geometry columns stay in place and are recorded in `x$roles`
#' (`node_geometry`, `edge_geometry`, `crs`); both tables must then have
#' the same coordinate reference system. Connectivity still comes from
#' `from` and `to`, never from the geometry. See [to_sf()] for what the
#' geometries are used for.
#'
#' @section Threads:
#' Measures that run in parallel take `threads`. `0` (the default) uses all
#' available cores, except that only 2 are used when the environment
#' variable `_R_CHECK_LIMIT_CORES_` is set (as it is by
#' `R CMD check --as-cran`), and at most `OMP_THREAD_LIMIT` when that is
#' set. Results are reproducible for a given number of threads: integer
#' results never depend on it, and floating-point results can differ in
#' the last bits between thread counts. Computations run in compiled code
#' and cannot be interrupted until they return.
#'
#' @return An object of class `netframe`.
#' @export
netframe <- function(edges, nodes = NULL, directed = TRUE,
                     from = 1L, to = 2L, id = 1L, weight = NULL,
                     weight_type = c("distance", "strength"), graph = list()) {
  weight_type <- match.arg(weight_type)
  if (!is.data.frame(edges)) {
    stop(sprintf("`edges` must be a data frame, not %s.", describe_class(edges)),
         call. = FALSE)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("`directed` must be TRUE or FALSE.", call. = FALSE)
  }
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
  from <- resolve_column(edges, from, "from", "edges")
  to <- resolve_column(edges, to, "to", "edges")
  if (!is.null(weight)) {
    weight <- resolve_column(edges, weight, "weight", "edges")
    if (!is.numeric(edges[[weight]])) {
      stop(sprintf("Weight column `%s` must be numeric, not %s.", weight,
                   describe_class(edges[[weight]])), call. = FALSE)
    }
  }
  src_ids <- check_ids(edges[[from]], from, "edges", "Edge endpoints")
  tgt_ids <- check_ids(edges[[to]], to, "edges", "Edge endpoints")

  if (is.null(nodes)) {
    if ((is_classed(src_ids) || is_classed(tgt_ids)) &&
        !identical(class(src_ids), class(tgt_ids))) {
      stop(sprintf(paste0(
        "Columns `%s` and `%s` of `edges` hold different types (%s and %s); ",
        "convert them to the same type first."),
        from, to, describe_class(src_ids), describe_class(tgt_ids)), call. = FALSE)
    }
    ids <- unique(c(src_ids, tgt_ids))
    nodes <- data.frame(name = ids, stringsAsFactors = FALSE)
    id <- "name"
    node_ids <- ids
  } else {
    if (!is.data.frame(nodes)) {
      stop(sprintf("`nodes` must be a data frame or NULL, not %s.", describe_class(nodes)),
           call. = FALSE)
    }
    id <- resolve_column(nodes, id, "id", "nodes")
    node_ids <- check_ids(nodes[[id]], id, "nodes", "Node identifiers")
    if (anyDuplicated(node_ids)) {
      dup <- nodes[[id]][anyDuplicated(node_ids)]
      stop(sprintf(paste0("Node identifiers must be unique; `%s` appears more than once ",
                          "in column `%s` of `nodes`."), format(dup), id), call. = FALSE)
    }
  }

  src <- match(src_ids, node_ids)
  tgt <- match(tgt_ids, node_ids)
  missing <- c(edges[[from]][is.na(src)], edges[[to]][is.na(tgt)])
  if (length(missing)) {
    shown <- utils::head(unique(missing), 5L)
    stop(sprintf(paste0(
      "%d edge endpoint(s) are not in `nodes`, for example: %s. Add them to `nodes`, ",
      "or use `nodes = NULL` to take the node set from the edges."),
      length(unique(missing)), paste(format(shown), collapse = ", ")
    ), call. = FALSE)
  }

  out <- structure(
    list(
      nodes = nodes,
      edges = edges,
      directed = isTRUE(directed),
      roles = list(id = id, from = from, to = to, weight = weight,
                   weight_type = if (!is.null(weight)) weight_type),
      graph = graph,
      index = list(src = src - 1L, tgt = tgt - 1L),
      cache = new.env(parent = emptyenv())
    ),
    class = "netframe"
  )
  # sf node or edge tables also record their geometry columns and CRS.
  record_geometry(out)
}

#' @export
print.netframe <- function(x, ...) {
  cat(sprintf(
    "<netframe> %s network with %d nodes and %d edges%s\n",
    if (x$directed) "directed" else "undirected",
    nrow(x$nodes), nrow(x$edges),
    if (is.null(x$roles$weight)) ", unweighted"
    else sprintf(", weighted by `%s` (%ss)", x$roles$weight, x$roles$weight_type)
  ))
  # Only unusual features get a line, so ordinary networks print one line.
  s <- edge_summary(x)
  notes <- c(
    if (s$loops) plural(s$loops, "self-loop"),
    if (s$parallel) plural(s$parallel, "parallel edge")
  )
  if (length(notes)) {
    cat(sprintf("including %s; each measure's help page says how they count\n",
                paste(notes, collapse = " and ")))
  }
  if (!is.null(x$roles$start)) cat(temporal_summary(x))
  cat(spatial_summary(x))
  if (s$bad_weights) {
    cat(sprintf(paste0("%s in weight column `%s` %s missing or not finite; ",
                       "measures will refuse the network until %s fixed\n"),
                plural(s$bad_weights, "value"), x$roles$weight,
                if (s$bad_weights == 1L) "is" else "are",
                if (s$bad_weights == 1L) "it is" else "they are"))
  }
  invisible(x)
}

# Counts of self-loops, parallel edges (edges repeating an earlier edge
# between the same ordered pair of nodes, or unordered pair in undirected
# networks), and missing or non-finite weights.
edge_summary <- function(x) {
  src <- x$index$src
  tgt <- x$index$tgt
  n <- nrow(x$nodes)
  a <- if (x$directed) src else pmin(src, tgt)
  b <- if (x$directed) tgt else pmax(src, tgt)
  # a * n + b is exact in double precision while n^2 < 2^53.
  key <- if (n < 2^26) as.double(a) * n + b else paste(a, b)
  w <- if (is.null(x$roles$weight)) NULL else x$edges[[x$roles$weight]]
  list(loops = sum(src == tgt), parallel = sum(duplicated(key)),
       bad_weights = if (is.null(w)) 0L else sum(!is.finite(w)))
}

plural <- function(k, what) sprintf("%d %s%s", k, what, if (k == 1L) "" else "s")

resolve_column <- function(df, col, role, what) {
  if (is.numeric(col) && length(col) == 1L && !is.na(col) && col >= 1 &&
      col <= ncol(df) && col == floor(col)) {
    return(names(df)[[col]])
  }
  if (is.character(col) && length(col) == 1L && !is.na(col) && col %in% names(df)) {
    return(col)
  }
  shown <- if (ncol(df)) {
    sprintf(" It has %s: %s.", plural(ncol(df), "column"),
            paste0("`", utils::head(names(df), 10L), "`", collapse = ", "))
  } else {
    " It has no columns."
  }
  stop(sprintf(paste0("Column for `%s` not found in `%s`; give a column name or ",
                      "position.%s"), role, what, shown), call. = FALSE)
}

# Checks an identifier column and returns it in the form used for
# matching: factors as character (as before), 64-bit integers from bit64 as
# character (their storage is not meaningful as doubles), other classed
# vectors such as dates unchanged, and plain vectors without attributes.
check_ids <- function(v, col, table, what) {
  if (!is.atomic(v) || is.matrix(v)) {
    stop(sprintf("%s (column `%s` of `%s`) must be an atomic vector, not %s.",
                 what, col, table, describe_class(v)), call. = FALSE)
  }
  v <- if (inherits(v, "integer64")) {
    as.character(v)
  } else if (is.factor(v)) {
    as.vector(v)
  } else if (is.object(v)) {
    v
  } else {
    as.vector(v)
  }
  bad <- which(is.na(v))
  if (length(bad)) {
    stop(sprintf(paste0("%s must not contain missing values (NA); %s of `%s` %s one ",
                        "in column `%s`. Remove those rows or fill in the identifiers."),
                 what, rows_text(bad), table, if (length(bad) == 1L) "has" else "have",
                 col), call. = FALSE)
  }
  if (is.character(v)) {
    bad <- which(!nzchar(v))
    if (length(bad)) {
      stop(sprintf(paste0("%s must not be empty strings; %s of `%s` %s one in column ",
                          "`%s`. Remove those rows or fill in the identifiers."),
                   what, rows_text(bad), table, if (length(bad) == 1L) "has" else "have",
                   col), call. = FALSE)
    }
  }
  v
}

is_classed <- function(v) is.object(v) && !is.factor(v)

# "row 3", "rows 3 and 7", or "rows 3, 7, 9, 12, 15, ... (40 in all)".
rows_text <- function(idx) {
  k <- length(idx)
  shown <- format(utils::head(idx, 5L), scientific = FALSE, trim = TRUE)
  if (k == 1L) return(paste("row", shown))
  if (k == 2L) return(paste("rows", shown[1], "and", shown[2]))
  if (k <= 5L) {
    return(paste0("rows ", paste(shown[-k], collapse = ", "), ", and ", shown[k]))
  }
  sprintf("rows %s, ... (%d in all)", paste(shown, collapse = ", "), k)
}

describe_class <- function(v) {
  if (is.null(v)) return("NULL")
  sprintf("an object of class `%s`", class(v)[[1]])
}

# Returns the engine graph for `x`, building it on first use. The handle is
# an external pointer, which does not survive saveRDS(); a dead handle is
# rebuilt transparently. The cache is shared by copies of `x` (it is an
# environment), so the handle is also rebuilt when the weights or the
# number of nodes no longer match what it was built from, for example after
# `x$edges$w <- 2 * x$edges$w`.
engine_graph <- function(x) {
  check_index(x)
  cache <- x$cache
  w <- if (is.null(x$roles$weight)) NULL else x$edges[[x$roles$weight]]
  ptr <- cache$ptr
  fresh <- !is.null(ptr) && .Call(C_graph_alive, ptr) &&
    identical(cache$n, nrow(x$nodes)) && identical(cache$weights, w)
  if (!fresh) {
    m <- length(x$index$src)
    if (m && max(x$index$src, x$index$tgt) >= nrow(x$nodes)) index_error()
    if (!is.null(w)) {
      if (!is.numeric(w)) {
        stop(sprintf("Weight column `%s` must be numeric, not %s.", x$roles$weight,
                     describe_class(w)), call. = FALSE)
      }
      bad <- which(!is.finite(w))
      if (length(bad)) {
        stop(sprintf(paste0(
          "Weight column `%s` has missing or non-finite values (NA, NaN, or Inf) in %s ",
          "of `x$edges`. Every measure needs finite weights, even when it ignores them: ",
          "remove or replace those values, or rebuild the network without `weight =`."),
          x$roles$weight, rows_text(bad)), call. = FALSE)
      }
    }
    # Anything else cached was derived from the old handle.
    rm(list = ls(cache, all.names = TRUE), envir = cache)
    ptr <- engine_call("netframe",
      .Call(C_graph_new, nrow(x$nodes), x$index$src, x$index$tgt,
            if (is.null(w)) NULL else as.double(w), x$directed))
    cache$ptr <- ptr
    cache$n <- nrow(x$nodes)
    cache$weights <- w
  }
  ptr
}

check_netframe <- function(x) {
  if (!inherits(x, "netframe")) {
    stop(sprintf(paste0("`x` must be a netframe, not %s; create one with netframe() ",
                        "or convert with as_netframe()."), describe_class(x)),
         call. = FALSE)
  }
  check_index(x)
}

# The edge index is computed once by netframe(); rows added to or removed
# from `x$nodes` or `x$edges` afterwards would silently misalign results.
# This cheap check runs on every call; engine_graph() also checks that every
# endpoint is still a row of `x$nodes` whenever it builds a handle.
check_index <- function(x) {
  if (!is.data.frame(x$nodes) || !is.data.frame(x$edges) ||
      nrow(x$edges) != length(x$index$src)) {
    index_error()
  }
}

index_error <- function() {
  stop(paste0("Rows of `x$nodes` or `x$edges` were added or removed after the ",
              "netframe was created, so its edge index no longer matches. Rebuild ",
              "it with netframe(x$edges, x$nodes, ...)."), call. = FALSE)
}

# Runs `expr`, a .Call into the engine, and rewrites its errors for R users:
# the message names the R function, and edge numbers (0-based positions in
# the engine's edge list) become rows of `x$edges` when the handle was built
# from the edge table in its own order (`edge_rows = TRUE`); handles built by
# collapsing parallel edges number edges differently, so the number is
# dropped. Every argument of `expr` other than the .Call itself must already
# be evaluated, so that only errors from compiled code are rewritten.
engine_call <- function(fn, expr, edge_rows = TRUE) {
  tryCatch(expr, error = function(e) {
    stop(engine_message(fn, conditionMessage(e), edge_rows), call. = FALSE)
  })
}

engine_message <- function(fn, msg, edge_rows = TRUE) {
  msg <- sub("[.[:space:]]+$", "", msg)
  hit <- regmatches(msg, regexec("^edge ([0-9]+) ", msg))[[1]]
  if (length(hit)) {
    where <- if (edge_rows) {
      sprintf("the edge in row %s of `x$edges` ",
              format(as.numeric(hit[2]) + 1, scientific = FALSE))
    } else {
      "an edge "
    }
    msg <- paste0(where, substring(msg, nchar(hit[1]) + 1L))
  }
  hint <- if (grepl("weight", msg) &&
              grepl("strictly positive|non-?negative|negative weight", msg)) {
    paste0(" Use `weighted = FALSE` to ignore weights, or remove or transform ",
           "the offending weights.")
  } else {
    ""
  }
  sprintf("%s(): %s.%s", fn, msg, hint)
}

# Thread counts. `0` means all available cores, except where the
# environment limits them: R CMD check sets _R_CHECK_LIMIT_CORES_ (CRAN
# allows at most 2 cores during checks), and OMP_THREAD_LIMIT is the usual
# way machines cap threads.
resolve_threads <- function(threads) {
  if (!is.numeric(threads) || length(threads) != 1L || is.na(threads) || threads < 0 ||
      threads != floor(threads) || threads > .Machine$integer.max) {
    stop(paste0("`threads` must be a single non-negative whole number; ",
                "0 uses all available cores."), call. = FALSE)
  }
  threads <- as.integer(threads)
  if (threads == 0L) threads <- default_threads()
  threads
}

default_threads <- function() {
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && !identical(tolower(chk), "false")) return(2L)
  limit <- suppressWarnings(as.integer(Sys.getenv("OMP_THREAD_LIMIT", "")))
  if (!is.na(limit) && limit > 0L) return(limit)
  0L
}

check_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}

# For arguments whose `NULL` default means "decided elsewhere".
check_flag_or_null <- function(value, name) {
  if (!is.null(value) && (!is.logical(value) || length(value) != 1L || is.na(value))) {
    stop(sprintf("`%s` must be TRUE, FALSE, or NULL.", name), call. = FALSE)
  }
  value
}

# match.arg() for arguments that are not formals with a vector of choices.
check_choice <- function(value, choices, name) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !value %in% choices) {
    stop(sprintf("`%s` must be one of %s.", name,
                 paste0("\"", choices, "\"", collapse = ", ")), call. = FALSE)
  }
  value
}
