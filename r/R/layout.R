# Node coordinates for drawing a network (man/layout_network.Rd).
layout_network <- function(x, method = c("auto", "stress", "sparse_stress", "fr", "circle",
                                         "star", "grid", "tree", "twomode"),
                           seed = NULL, weighted = NULL, weights_as = NULL, types = NULL,
                           root = NULL, center = NULL, max_iter = NULL, tol = NULL,
                           pivots = NULL, approximation = c("auto", "exact", "barnes_hut"),
                           threads = 0L) {
  check_netframe(x)
  method <- match.arg(method)
  approximation <- match.arg(approximation)
  threads <- resolve_threads(threads)
  # Arguments are checked before the engine runs, and the seed is drawn
  # only for the methods that use it.
  args <- switch(method,
    auto = , stress = , sparse_stress = list(
      w = resolve_weighting(x, weighted, weights_as),
      pivots = layout_count(pivots, "pivots"), max_iter = layout_count(max_iter, "max_iter"),
      tol = layout_tol(tol), seed = check_seed(seed)),
    fr = list(
      w = resolve_weighting(x, weighted, weights_as),
      max_iter = layout_count(max_iter, "max_iter"), seed = check_seed(seed)),
    star = list(center = if (is.null(center)) -1 else
      node_positions(x, center, "center", one = TRUE) - 1),
    tree = list(roots = if (is.null(root)) double() else node_positions(x, root, "root") - 1),
    twomode = list(types = layout_types(x, types)),
    list()
  )
  ptr <- engine_graph(x)
  xy <- engine_call("layout_network", switch(method,
    auto = , stress = , sparse_stress = {
      algorithm <- switch(method, auto = 0L, stress = 1L, sparse_stress = 2L)
      .Call(C_layout_stress, ptr, algorithm, args$w$weighted, args$w$code, args$pivots,
            args$max_iter, args$tol, args$seed, threads)
    },
    fr = {
      code <- switch(approximation, auto = 0L, exact = 1L, barnes_hut = 2L)
      .Call(C_layout_fr, ptr, args$w$weighted, args$w$code, args$max_iter, code, 0,
            args$seed, threads)
    },
    circle = .Call(C_layout_simple, ptr, 0L, 0),
    star = .Call(C_layout_simple, ptr, 1L, as.double(args$center)),
    grid = .Call(C_layout_simple, ptr, 2L, 0),
    tree = .Call(C_layout_tree, ptr, as.double(args$roots)),
    twomode = .Call(C_layout_twomode, ptr, args$types, 0)
  ))
  n <- nrow(x$nodes)
  matrix(xy, nrow = n, ncol = 2L, dimnames = list(NULL, c("x", "y")))
}

# See man/layout_network.Rd.
layout_stress <- function(x, layout, weighted = NULL, weights_as = NULL, threads = 0L) {
  check_netframe(x)
  xy <- layout_matrix(layout, nrow(x$nodes))
  w <- resolve_weighting(x, weighted, weights_as)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  engine_call("layout_stress",
              .Call(C_layout_stress_value, ptr, as.double(xy), w$weighted, w$code, threads))
}

# Layout controls whose NULL default (0 for the engine) means "choose".
layout_count <- function(value, name) {
  if (is.null(value)) 0 else as.double(check_count(value, name, min = 1L))
}

layout_tol <- function(tol) {
  if (is.null(tol)) 0 else check_number(tol, "tol", lower = 0, open = TRUE)
}

# Node types for a two-mode layout as 0 (top row) and 1 (bottom row).
layout_types <- function(x, types) {
  n <- nrow(x$nodes)
  if (is.null(types)) {
    col <- if (!is.null(x$roles$type)) x$roles$type else "type"
    if (!col %in% names(x$nodes)) {
      stop(paste0("A two-mode layout needs node types: supply `types =`, use a network ",
                  "made by twomode() or as_twomode(), or add a `type` column to the node ",
                  "table."), call. = FALSE)
    }
    types <- x$nodes[[col]]
  } else if (is.character(types) && length(types) == 1L && n != 1L) {
    if (!types %in% names(x$nodes)) {
      stop(sprintf("Column `%s` not found in the node table.", types), call. = FALSE)
    }
    types <- x$nodes[[types]]
  }
  if (length(types) != n) {
    stop("`types` must have one value per node.", call. = FALSE)
  }
  if (anyNA(types)) stop("`types` must not contain missing values.", call. = FALSE)
  # Radix sorting orders text by bytes, the same in every locale.
  values <- if (is.factor(types)) levels(droplevels(types)) else
    sort(unique(types), method = "radix")
  if (length(values) > 2L) {
    stop("`types` must have at most two distinct values for a two-mode layout.",
         call. = FALSE)
  }
  key <- if (is.factor(types)) as.character(types) else types
  as.integer(match(as.character(key), as.character(values)) - 1L)
}

# A layout given as a matrix or data frame, checked against the node count.
layout_matrix <- function(layout, n) {
  if (is.data.frame(layout)) {
    cols <- if (all(c("x", "y") %in% names(layout))) c("x", "y") else seq_len(min(2L, ncol(layout)))
    layout <- as.matrix(layout[, cols, drop = FALSE])
  }
  if (!is.matrix(layout) || !is.numeric(layout) || ncol(layout) != 2L ||
      nrow(layout) != n) {
    stop(sprintf("A layout must be a numeric matrix or data frame with 2 columns and %d rows.",
                 n), call. = FALSE)
  }
  if (!all(is.finite(layout))) {
    stop("Layout coordinates must be finite.", call. = FALSE)
  }
  storage.mode(layout) <- "double"
  layout
}
