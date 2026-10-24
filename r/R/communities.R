# Community detection beyond Louvain (Walktrap, label propagation,
# Infomap, Leiden, spinglass, exact maximum modularity, fluid and Voronoi
# communities), comparison of partitions, and dendrogram cutting and
# conversion. Memberships are integer vectors aligned with `x$nodes`,
# numbered 1, 2, ... in order of first appearance, with the scalar results
# in attributes, as net_louvain() returns them.

# Walktrap communities (man/net_walktrap.Rd).
net_walktrap <- function(x, steps = 4L, weighted = NULL, threads = 0L) {
  check_net_frame(x)
  require_undirected(x, "Walktrap community detection")
  steps <- check_count(steps, "steps", min = 1L)
  weighted <- resolve_additive_weights(x, weighted)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  res <- engine_call("net_walktrap", .Call(C_walktrap, ptr, steps, weighted, threads))
  merges <- res[[3L]]
  attr(merges, "leaves") <- nrow(x$nodes)
  structure(res[[1L]], modularity = res[[2L]], merges = merges, modularity_path = res[[4L]])
}

# Label propagation communities (man/net_label_prop.Rd).
net_label_prop <- function(x, variant = c("dominance", "retention", "fast", "semi_synchronous"),
                           mode = c("all", "out", "in"), weighted = NULL, initial = NULL,
                           fixed = NULL, seed = NULL) {
  check_net_frame(x)
  variant <- match.arg(variant)
  mode <- match.arg(mode)
  if (variant == "semi_synchronous") {
    require_undirected(x, "The semi-synchronous variant of label propagation")
  }
  weighted <- resolve_additive_weights(x, weighted)
  if (!is.null(initial)) {
    initial <- node_categories(x, initial, "initial", missing_ok = TRUE)$codes
  }
  if (!is.null(fixed)) {
    fixed <- node_column(x, fixed, "fixed")
    if (!is.logical(fixed) || anyNA(fixed)) {
      stop("`fixed` must be a logical vector without missing values, one per node.",
           call. = FALSE)
    }
  }
  seed <- check_seed(seed)
  ptr <- engine_graph(x)
  code <- match(variant, c("dominance", "retention", "fast", "semi_synchronous")) - 1L
  engine_call("net_label_prop",
              .Call(C_label_propagation, ptr, code, mode_code(mode), weighted, initial, fixed,
                    seed))
}

# Infomap communities (man/net_infomap.Rd).
net_infomap <- function(x, weighted = NULL, node_weights = NULL, trials = 10L, seed = NULL) {
  check_net_frame(x)
  weighted <- resolve_additive_weights(x, weighted)
  if (!is.null(node_weights)) node_weights <- node_values(x, node_weights, "node_weights")
  trials <- check_count(trials, "trials", min = 1L)
  seed <- check_seed(seed)
  ptr <- engine_graph(x)
  res <- engine_call("net_infomap", .Call(C_infomap, ptr, weighted, node_weights, trials, seed))
  structure(res[[1L]], codelength = res[[2L]])
}

# The map equation of a partition (man/net_infomap.Rd).
net_map_equation <- function(x, membership, weighted = NULL, node_weights = NULL) {
  check_net_frame(x)
  labels <- node_groups(x, membership, "membership")$codes
  weighted <- resolve_additive_weights(x, weighted)
  if (!is.null(node_weights)) node_weights <- node_values(x, node_weights, "node_weights")
  ptr <- engine_graph(x)
  engine_call("net_map_equation", .Call(C_map_equation, ptr, labels, weighted, node_weights))
}

# Leiden communities (man/net_leiden.Rd).
net_leiden <- function(x, objective = c("modularity", "cpm"), resolution = 1, beta = 0.01,
                       n_iterations = -1L, weighted = NULL, node_weights = NULL,
                       initial = NULL, seed = NULL) {
  check_net_frame(x)
  require_undirected(x, "Leiden community detection")
  objective <- match.arg(objective)
  resolution <- check_resolution(resolution)
  beta <- check_number(beta, "beta", lower = 0, open = TRUE)
  if (!is.numeric(n_iterations) || length(n_iterations) != 1L || is.na(n_iterations) ||
      n_iterations != floor(n_iterations) || abs(n_iterations) > .Machine$integer.max) {
    stop("`n_iterations` must be a single whole number; a negative number runs until stable.",
         call. = FALSE)
  }
  weighted <- resolve_additive_weights(x, weighted)
  if (!is.null(node_weights)) node_weights <- node_values(x, node_weights, "node_weights")
  if (!is.null(initial)) initial <- node_groups(x, initial, "initial")$codes
  seed <- check_seed(seed)
  ptr <- engine_graph(x)
  res <- engine_call("net_leiden",
                     .Call(C_leiden, ptr, if (objective == "cpm") 0L else 1L, resolution, beta,
                           as.double(n_iterations), weighted, node_weights, initial, seed))
  structure(res[[1L]], quality = res[[2L]], iterations = res[[3L]])
}

# Spinglass communities (man/net_spinglass.Rd).
net_spinglass <- function(x, spins = 25L, parupdate = FALSE, start_temp = 1, stop_temp = 0.01,
                          cool_fact = 0.99, update_rule = c("config", "simple"), gamma = 1,
                          implementation = c("orig", "neg"), gamma_minus = 1, weighted = NULL,
                          seed = NULL) {
  check_net_frame(x)
  spins <- check_count(spins, "spins", min = 1L)
  parupdate <- check_flag(parupdate, "parupdate")
  start_temp <- check_number(start_temp, "start_temp", lower = 0, open = TRUE)
  stop_temp <- check_number(stop_temp, "stop_temp", lower = 0, open = TRUE)
  if (stop_temp > start_temp) stop("`stop_temp` must not exceed `start_temp`.", call. = FALSE)
  cool_fact <- check_number(cool_fact, "cool_fact", lower = 0, upper = 1, open = TRUE)
  update_rule <- match.arg(update_rule)
  gamma <- check_number(gamma, "gamma", lower = 0)
  implementation <- match.arg(implementation)
  gamma_minus <- check_number(gamma_minus, "gamma_minus", lower = 0)
  weighted <- resolve_additive_weights(x, weighted)
  seed <- check_seed(seed)
  ptr <- engine_graph(x)
  res <- engine_call("net_spinglass",
                     .Call(C_spinglass, ptr, spins, parupdate, start_temp, stop_temp, cool_fact,
                           if (update_rule == "config") 1L else 0L, gamma,
                           implementation == "neg", gamma_minus, weighted, seed))
  structure(res[[1L]], modularity = res[[2L]], energy = res[[3L]])
}

# Exact maximum modularity (man/net_optimal.Rd).
net_optimal <- function(x, resolution = 1, weighted = NULL) {
  check_net_frame(x)
  resolution <- check_resolution(resolution)
  weighted <- resolve_additive_weights(x, weighted)
  ptr <- engine_graph(x)
  res <- engine_call("net_optimal", .Call(C_optimal_modularity, ptr, resolution, weighted))
  structure(res[[1L]], modularity = res[[2L]])
}

# Fluid communities (man/net_fluid.Rd).
net_fluid <- function(x, k, seed = NULL) {
  check_net_frame(x)
  require_undirected(x, "Fluid community detection")
  k <- check_count(k, "k", min = 1L)
  seed <- check_seed(seed)
  ptr <- engine_graph(x)
  engine_call("net_fluid", .Call(C_fluid_communities, ptr, k, seed))
}

# Voronoi communities (man/net_voronoi.Rd).
net_voronoi <- function(x, lengths = NULL, weighted = NULL, mode = c("out", "in", "all"),
                        radius = NULL, threads = 0L) {
  check_net_frame(x)
  mode <- match.arg(mode)
  weighted <- resolve_additive_weights(x, weighted)
  if (!is.null(lengths)) {
    if (is.character(lengths) && length(lengths) == 1L && !is.na(lengths)) {
      if (!lengths %in% names(x$edges)) {
        stop("`lengths` names no column of `x$edges`.", call. = FALSE)
      }
      lengths <- x$edges[[lengths]]
    }
    if (!is.numeric(lengths) || length(lengths) != nrow(x$edges) || anyNA(lengths) ||
        any(!is.finite(lengths)) || any(lengths < 0)) {
      stop(sprintf(paste0("`lengths` must be non-negative finite numbers, one per edge (%d), ",
                          "or the name of such a column of `x$edges`."), nrow(x$edges)),
           call. = FALSE)
    }
    lengths <- as.double(lengths)
  }
  radius <- if (is.null(radius)) -1 else check_number(radius, "radius", lower = 0, finite = FALSE)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  rows <- x$cache$graph_rows
  if (!is.null(lengths) && !is.null(rows)) lengths <- lengths[rows]
  res <- engine_call("net_voronoi",
                     .Call(C_voronoi, ptr, lengths, weighted, mode_code(mode), radius, threads))
  structure(res[[1L]], generators = res[[2L]], modularity = res[[3L]], radius = res[[4L]])
}

# Distances and similarities between partitions (man/net_compare_partitions.Rd).
net_compare_partitions <- function(a, b, method = c("vi", "nmi", "split_join", "rand",
                                                    "adjusted_rand", "ami")) {
  method <- match.arg(method)
  codes <- function(v, name) {
    if (!is.atomic(v) || is.null(v) || is.matrix(v)) {
      stop(sprintf("`%s` must be a vector of group labels.", name), call. = FALSE)
    }
    if (anyNA(v)) stop(sprintf("`%s` must not contain missing values.", name), call. = FALSE)
    if (is.factor(v)) v <- as.character(v)
    match(v, unique(v)) - 1L
  }
  ca <- codes(a, "a")
  cb <- codes(b, "b")
  if (length(ca) != length(cb)) {
    stop("`a` and `b` must have the same length.", call. = FALSE)
  }
  code <- match(method, c("vi", "nmi", "split_join", "rand", "adjusted_rand", "ami")) - 1L
  engine_call("net_compare_partitions", .Call(C_compare_partitions, ca, cb, code))
}

# A two-column dendrogram in igraph's numbering with integer storage;
# `leaves` (the number of leaves) is kept as an attribute when given.
check_merges <- function(merges) {
  if (!is.matrix(merges) || ncol(merges) != 2L || !is.numeric(merges) || anyNA(merges) ||
      any(merges != floor(merges)) || any(merges < 1)) {
    stop(paste0("`merges` must be a two-column matrix of positive whole numbers, such as the ",
                "`merges` attribute of net_walktrap()."), call. = FALSE)
  }
  leaves <- attr(merges, "leaves")
  out <- matrix(as.integer(merges), ncol = 2L)
  if (!is.null(leaves)) attr(out, "leaves") <- leaves
  out
}

# The number of leaves of a dendrogram: the nodes of `x` when given, else
# the recorded count, else one more than the merges (a complete tree).
dendrogram_leaves <- function(merges, x) {
  recorded <- attr(merges, "leaves")
  if (!is.null(x)) {
    n <- nrow(x$nodes)
    if (!is.null(recorded) && recorded != n) {
      stop(sprintf("`merges` describes %d leaves but `x` has %d nodes.", recorded, n),
           call. = FALSE)
    }
    return(n)
  }
  if (is.null(recorded)) nrow(merges) + 1L else as.integer(recorded)
}

# Cut a dendrogram (man/net_cut_dendrogram.Rd).
net_cut_dendrogram <- function(merges, k = NULL, steps = NULL, x = NULL, weighted = NULL) {
  merges <- check_merges(merges)
  if (!is.null(k) && !is.null(steps)) stop("Give at most one of `k` and `steps`.", call. = FALSE)
  if (is.null(k) && is.null(steps) && is.null(x)) {
    stop("Give `k`, `steps`, or `x` (to cut at the level of highest modularity).", call. = FALSE)
  }
  if (!is.null(x)) check_net_frame(x)
  n <- dendrogram_leaves(merges, x)
  a <- merges[, 1L]
  b <- merges[, 2L]
  if (is.null(k) && is.null(steps)) {
    require_undirected(x, "Cutting a dendrogram at maximum modularity")
    weighted <- resolve_additive_weights(x, weighted)
    ptr <- engine_graph(x)
    res <- engine_call("net_cut_dendrogram",
                       .Call(C_dendrogram_modularity, ptr, a, b, 1, weighted))
    return(structure(res[[2L]], modularity = res[[3L]], modularity_path = res[[1L]]))
  }
  groups <- if (is.null(k)) 0L else check_count(k, "k", min = 1L)
  stops <- if (is.null(steps)) 0L else check_count(steps, "steps")
  engine_call("net_cut_dendrogram", .Call(C_cut_dendrogram, n, a, b, stops, groups))
}

# Convert a dendrogram to an hclust object (man/net_cut_dendrogram.Rd).
net_to_hclust <- function(merges, x = NULL) {
  merges <- check_merges(merges)
  if (!is.null(x)) check_net_frame(x)
  n <- dendrogram_leaves(merges, x)
  k <- nrow(merges)
  limit <- n + seq_len(k) - 1L
  if (any(merges > limit) || anyDuplicated(as.vector(merges)) || any(merges[, 1L] == merges[, 2L])) {
    stop(paste0("`merges` is not a dendrogram: each row must join two different existing ",
                "clusters, and each cluster can be joined only once."), call. = FALSE)
  }
  labels <- if (is.null(x)) as.character(seq_len(n)) else node_ids(x)
  # A network with several components gives an incomplete dendrogram; the
  # remaining roots are joined one after another above the last merge.
  roots <- setdiff(seq_len(n + k), as.vector(merges))
  if (length(roots) > 1L) {
    cur <- roots[[1L]]
    for (r in roots[-1L]) {
      merges <- rbind(merges, c(cur, r))
      cur <- n + nrow(merges)
    }
  }
  total <- nrow(merges)
  # hclust numbers leaves -j and clusters by the row that made them, and
  # lists singletons first (in increasing order) before clusters.
  one <- function(v) if (v <= n) -v else v - n
  merge <- matrix(0L, total, 2L)
  for (i in seq_len(total)) {
    pair <- c(one(merges[i, 1L]), one(merges[i, 2L]))
    neg <- pair[pair < 0L]
    pos <- pair[pair > 0L]
    merge[i, ] <- as.integer(c(sort(neg, decreasing = TRUE), sort(pos)))
  }
  structure(list(merge = merge, height = as.double(seq_len(total)),
                 order = if (total) leaf_order(merge) else seq_len(n), labels = labels,
                 method = "merges", call = NULL, dist.method = NULL), class = "hclust")
}
