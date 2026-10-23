# The bridge between a net_frame and the compiled engine: the graph handle
# and its cache, the wrapper that rewrites engine errors for R users,
# thread counts, and how a measure decides whether to use weights. Every
# measure follows the same sequence: check_net_frame(x), check its own
# arguments, resolve_threads(), engine_graph(x), engine_call().

check_net_frame <- function(x) {
  if (!inherits(x, "net_frame")) {
    stop(sprintf(paste0("`x` must be a net_frame, not %s; create one with net_frame() ",
                        "or convert with as_net_frame()."), describe_class(x)),
         call. = FALSE)
  }
  check_index(x)
}

# The edge index is computed once by net_frame(); rows added to or removed
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
  stop(paste0("Rows of `x$nodes` or `x$edges` were added, removed, or reordered after ",
              "the net_frame was created, so its edge index no longer matches. Rebuild ",
              "it with net_frame(x$edges, x$nodes, ...)."), call. = FALSE)
}

# Stops unless the network is undirected; `what` names the measure.
require_undirected <- function(x, what) {
  if (x$directed) {
    stop(what, " needs an undirected network; rebuild it with ",
         "net_frame(..., directed = FALSE) or net_symmetrize() it.", call. = FALSE)
  }
  invisible(TRUE)
}

# ---- Engine handles ------------------------------------------------------------

# Returns the engine graph for `x`, building it on first use. The handle is
# an external pointer, which does not survive saveRDS(); a dead handle is
# rebuilt transparently, through the frame, so the tables are validated
# again. The cache is shared by copies of `x` (it is an environment), so
# the handle is also rebuilt when the weights or the number of nodes or
# edges no longer match what it was built from, for example after
# `x$edges$w <- 2 * x$edges$w`.
engine_graph <- function(x) {
  check_index(x)
  cache <- x$cache
  w <- if (is.null(x$roles$weight)) NULL else x$edges[[x$roles$weight]]
  ptr <- cache$ptr
  fresh <- !is.null(ptr) && .Call(C_graph_alive, ptr) &&
    !is.null(cache$frame) && .Call(C_graph_alive, cache$frame) &&
    identical(cache$n, nrow(x$nodes)) && identical(cache$m, nrow(x$edges)) &&
    identical(cache$weights, w)
  if (!fresh) {
    m <- length(x$index$src)
    if (m && max(x$index$src, x$index$tgt) >= nrow(x$nodes)) index_error()
    f <- engine_frame(x)
    d <- frame_diagnostics(f, 1L)
    errors <- d[d$severity == "error", , drop = FALSE]
    if (nrow(errors)) {
      stop(paste0(format_diagnostics(errors, "net_frame", x), "\nThe tables of `x` were ",
                  "edited after the network was built; fix them or rebuild the network."),
           call. = FALSE)
    }
    ptr <- engine_call("net_frame", .Call(C_frame_build, f, 1L))
    idx <- .Call(C_frame_index, f, as.double(nrow(x$edges)))
    if (length(idx$added_row) || !identical(idx$src, x$index$src) ||
        !identical(idx$tgt, x$index$tgt)) {
      index_error()
    }
    # Anything else cached was derived from the old handle.
    rm(list = ls(cache, all.names = TRUE), envir = cache)
    cache$frame <- f
    cache$ptr <- ptr
    cache$n <- nrow(x$nodes)
    cache$m <- nrow(x$edges)
    cache$weights <- w
  }
  ptr
}

# The network as an adjacency matrix, the way sna and UCINET read it:
# parallel edges collapse into one tie (keeping the largest weight, or the
# smallest with keep = "min") and self-loops are dropped. Cached alongside
# the engine graph.
engine_graph_simple <- function(x, undirected = FALSE, keep = c("max", "min")) {
  keep <- match.arg(keep)
  key <- paste("simple", isTRUE(undirected), keep, sep = "_")
  # engine_graph() first: it empties the cache when the handle is stale.
  base <- engine_graph(x)
  ptr <- x$cache[[key]]
  if (is.null(ptr) || !.Call(C_graph_alive, ptr)) {
    ptr <- engine_call("net_frame", .Call(C_graph_simplify, base, isTRUE(undirected),
                                         if (keep == "max") 0L else 1L))
    assign(key, ptr, envir = x$cache)
  }
  ptr
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

# The engine's code for a `mode` argument.
mode_code <- function(mode) switch(mode, out = 1L, `in` = 2L, all = 3L)

# ---- Threads ---------------------------------------------------------------------

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

# ---- Weights ---------------------------------------------------------------------

require_weights <- function(x) {
  if (is.null(x$roles$weight)) {
    stop(paste0("This network has no weights; use `weighted = FALSE`, or supply ",
                "`weight =` to net_frame()."), call. = FALSE)
  }
  invisible(TRUE)
}

# Decides whether and how a path-based measure uses weights: `weighted`
# NULL means "when the network has them", and `weights_as` overrides the
# network's weight_type. Returns the flag and the engine's code for the
# kind (0 distances, 1 strengths).
resolve_weighting <- function(x, weighted, weights_as) {
  check_flag_or_null(weighted, "weighted")
  if (!is.null(weights_as)) check_choice(weights_as, c("distance", "strength"), "weights_as")
  has <- !is.null(x$roles$weight)
  weighted <- if (is.null(weighted)) has else isTRUE(weighted)
  if (weighted) require_weights(x)
  kind <- if (is.null(weights_as)) x$roles$weight_type else weights_as
  if (is.null(kind)) kind <- "distance"
  list(weighted = weighted, code = if (kind == "strength") 1L else 0L)
}

# Whether a measure that reads weights as given (strength, modularity,
# structural equivalence) should use them: `weighted` NULL means "when the
# network has them".
resolve_additive_weights <- function(x, weighted) {
  check_flag_or_null(weighted, "weighted")
  weighted <- if (is.null(weighted)) !is.null(x$roles$weight) else isTRUE(weighted)
  if (weighted) require_weights(x)
  weighted
}
