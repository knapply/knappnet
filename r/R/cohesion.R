# Bicomponents (blocks) (man/bicomponents.Rd).
bicomponents <- function(x) {
  check_netframe(x)
  ptr <- engine_graph(x)
  idx <- engine_call("bicomponents", .Call(C_bicomponents, ptr))
  ids <- x$nodes[[x$roles$id]]
  lapply(idx, function(i) ids[i])
}

# n-cliques, n-clans, and k-plexes (man/n_cliques.Rd).
n_cliques <- function(x, n = 2L, min_size = 1L, threads = 0L) {
  check_netframe(x)
  n <- check_count(n, "n", min = 1L)
  min_size <- check_count(min_size, "min_size")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  idx <- engine_call("n_cliques", .Call(C_n_cliques, ptr, n, min_size, FALSE, threads))
  ids <- x$nodes[[x$roles$id]]
  lapply(idx, function(i) ids[i])
}

# See man/n_cliques.Rd.
n_clans <- function(x, n = 2L, min_size = 1L, threads = 0L) {
  check_netframe(x)
  n <- check_count(n, "n", min = 1L)
  min_size <- check_count(min_size, "min_size")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  idx <- engine_call("n_clans", .Call(C_n_cliques, ptr, n, min_size, TRUE, threads))
  ids <- x$nodes[[x$roles$id]]
  lapply(idx, function(i) ids[i])
}

# See man/n_cliques.Rd.
k_plexes <- function(x, k = 2L, min_size = NULL, max_steps = 1e7, threads = 0L) {
  check_netframe(x)
  k <- check_count(k, "k", min = 1L)
  min_size <- if (is.null(min_size)) 2L * k - 1L else check_count(min_size, "min_size")
  if (!is.numeric(max_steps) || length(max_steps) != 1L || is.na(max_steps) || max_steps < 1) {
    stop("`max_steps` must be a number of at least 1.", call. = FALSE)
  }
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  idx <- engine_call("k_plexes",
                     .Call(C_k_plexes, ptr, k, min_size, as.double(max_steps), threads))
  ids <- x$nodes[[x$roles$id]]
  lapply(idx, function(i) ids[i])
}

# Lambda sets and edge connectivity (man/lambda_sets.Rd).
lambda_sets <- function(x) {
  check_netframe(x)
  ptr <- engine_graph(x)
  m <- engine_call("lambda_sets", .Call(C_lambda_connectivity, ptr))
  ids <- x$nodes[[x$roles$id]]
  dimnames(m) <- list(ids, ids)
  n <- nrow(m)
  sets <- list()
  seen <- character()
  top <- if (n) max(m) else 0
  for (k in seq_len(top)) {
    done <- logical(n)
    for (v in seq_len(n)) {
      if (done[v]) next
      members <- which(m[v, ] >= k | seq_len(n) == v)
      done[members] <- TRUE
      if (length(members) < 2L) next
      key <- paste(members, collapse = " ")
      if (key %in% seen) {
        # The same set at a higher level: record its true connectivity.
        sets[[match(key, seen)]]$lambda <- k
        next
      }
      seen <- c(seen, key)
      sets[[length(sets) + 1L]] <- list(lambda = k, members = ids[members])
    }
  }
  structure(list(connectivity = m, sets = sets), class = "lambda_sets")
}

# See man/lambda_sets.Rd.
print.lambda_sets <- function(x, ...) {
  cat(sprintf("Lambda sets of %d nodes: %d sets\n", nrow(x$connectivity), length(x$sets)))
  for (s in x$sets) {
    members <- format(s$members)
    shown <- if (length(members) > 12L) {
      paste0(paste(members[1:12], collapse = ", "), ", ... (", length(members), " nodes)")
    } else {
      paste(members, collapse = ", ")
    }
    cat(sprintf("  lambda = %d: %s\n", s$lambda, shown))
  }
  invisible(x)
}

# Factions (man/factions.Rd).
factions <- function(x, k = 2L, fit = c("errors", "correlation"), starts = 20L,
                     max_stall = 200L, tenure = NULL, seed = NULL, threads = 0L) {
  check_netframe(x)
  k <- check_count(k, "k")
  fit <- match.arg(fit)
  starts <- check_count(starts, "starts")
  max_stall <- check_count(max_stall, "max_stall")
  tenure <- if (is.null(tenure)) 0L else check_count(tenure, "tenure")
  seed <- check_seed(seed)
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  res <- engine_call("factions",
                     .Call(C_factions, ptr, k, if (fit == "errors") 0L else 1L, starts,
                           max_stall, tenure, seed, threads))
  structure(res[[1L]], errors = res[[2L]], correlation = res[[3L]])
}

# Girvan-Newman, leading eigenvector, and greedy modularity communities (man/cluster_girvan_newman.Rd).
cluster_girvan_newman <- function(x, threads = 0L) {
  check_netframe(x)
  require_undirected(x, "Girvan-Newman clustering")
  threads <- resolve_threads(threads)
  ptr <- engine_graph(x)
  res <- engine_call("cluster_girvan_newman", .Call(C_girvan_newman, ptr, threads))
  ids <- x$nodes[[x$roles$id]]
  removed <- data.frame(from = ids[res[[3L]]], to = ids[res[[4L]]], stringsAsFactors = FALSE)
  structure(res[[1L]], modularity = res[[2L]], removed = removed, path = res[[5L]],
            tree = removal_tree(length(ids), res[[3L]], res[[4L]], as.character(ids)))
}

# The divisive hierarchy of Girvan-Newman as an hclust object: read in
# reverse, the removals are merges. The edge removed last joins two
# singletons at height 1, the edge removed first makes the last merge at
# height m, and components that were never connected are joined above
# that, at height m + 1.
removal_tree <- function(n, from, to, labels) {
  m <- length(from)
  parent <- seq_len(n)
  find <- function(v) {
    while (parent[v] != v) v <- parent[v]
    v
  }
  # hclust's convention: singletons are negative, earlier merges positive.
  cluster <- -seq_len(n)
  merge <- matrix(0L, 0L, 2L)
  height <- numeric()
  join <- function(a, b, h) {
    ra <- find(a)
    rb <- find(b)
    if (ra == rb) return(FALSE)
    merge <<- rbind(merge, as.integer(c(cluster[ra], cluster[rb])))
    height <<- c(height, h)
    parent[rb] <<- ra
    cluster[ra] <<- nrow(merge)
    TRUE
  }
  for (i in rev(seq_len(m))) join(from[i], to[i], m - i + 1)
  roots <- unique(vapply(seq_len(n), find, integer(1)))
  if (length(roots) > 1L) {
    for (r in roots[-1L]) join(roots[[1L]], r, m + 1)
  }
  if (n < 2L) {
    return(structure(list(merge = merge, height = height, order = seq_len(n), labels = labels,
                          method = "girvan-newman"), class = "hclust"))
  }
  tree <- structure(list(merge = merge, height = height, order = integer(), labels = labels,
                         method = "girvan-newman"), class = "hclust")
  # Leaf order from the merge structure, as stats computes it.
  tree$order <- leaf_order(merge)
  tree
}

leaf_order <- function(merge) {
  # Depth-first over the merge tree with an explicit stack (no recursion,
  # so large trees cannot exhaust R's expression nesting).
  out <- integer()
  stack <- nrow(merge)
  while (length(stack)) {
    node <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    if (node < 0) {
      out <- c(out, -node)
    } else {
      stack <- c(stack, merge[node, 2L], merge[node, 1L])
    }
  }
  out
}

# See man/cluster_girvan_newman.Rd.
cluster_leading_eigen <- function(x, weighted = NULL) {
  check_netframe(x)
  require_undirected(x, "The leading eigenvector method")
  weighted <- resolve_additive_weights(x, weighted)
  ptr <- engine_graph(x)
  res <- engine_call("cluster_leading_eigen", .Call(C_leading_eigen, ptr, weighted))
  structure(res[[1L]], modularity = res[[2L]], splits = res[[3L]])
}

# See man/cluster_girvan_newman.Rd.
cluster_fast_greedy <- function(x, resolution = 1, weighted = NULL) {
  check_netframe(x)
  require_undirected(x, "Greedy modularity optimization")
  resolution <- check_resolution(resolution)
  weighted <- resolve_additive_weights(x, weighted)
  ptr <- engine_graph(x)
  res <- engine_call("cluster_fast_greedy", .Call(C_fast_greedy, ptr, resolution, weighted))
  structure(res[[1L]], modularity = res[[2L]])
}

# Hierarchical clustering of nodes (man/cluster_hierarchical.Rd).
cluster_hierarchical <- function(d, method = c("average", "single", "complete", "ward"),
                                 similarity = FALSE) {
  method <- match.arg(method)
  similarity <- check_flag(similarity, "similarity")
  if (!is.matrix(d) || !is.numeric(d) || nrow(d) != ncol(d)) {
    stop("`d` must be a square numeric matrix.", call. = FALSE)
  }
  if (anyNA(d)) stop("`d` must not contain missing values.", call. = FALSE)
  if (!isTRUE(all.equal(d, t(d), check.attributes = FALSE))) {
    stop("`d` must be symmetric.", call. = FALSE)
  }
  finite <- d[is.finite(d)]
  if (similarity) {
    if (any(!is.finite(d))) stop("A similarity matrix must be finite.", call. = FALSE)
    d <- max(finite) - d
  } else if (any(is.infinite(d))) {
    d[is.infinite(d)] <- if (length(finite)) max(finite) + 1 else 1
  }
  labels <- rownames(d)
  if (is.null(labels)) labels <- as.character(seq_len(nrow(d)))
  dist <- stats::as.dist(d)
  attr(dist, "Labels") <- labels
  stats::hclust(dist, method = if (method == "ward") "ward.D2" else method)
}

# See man/cluster_hierarchical.Rd.
cut_clusters <- function(tree, k = NULL, height = NULL) {
  if (!inherits(tree, "hclust")) {
    stop("`tree` must be the result of cluster_hierarchical() (an `hclust` object).",
         call. = FALSE)
  }
  if (is.null(k) == is.null(height)) {
    stop("Give exactly one of `k` and `height`.", call. = FALSE)
  }
  if (!is.null(k)) {
    k <- check_count(k, "k")
    if (k < 1L || k > length(tree$order)) {
      stop("`k` must be between 1 and the number of nodes.", call. = FALSE)
    }
    groups <- stats::cutree(tree, k = k)
  } else {
    if (!is.numeric(height) || length(height) != 1L || is.na(height)) {
      stop("`height` must be a single number.", call. = FALSE)
    }
    groups <- stats::cutree(tree, h = height)
  }
  out <- match(groups, unique(groups))
  names(out) <- NULL
  out
}
