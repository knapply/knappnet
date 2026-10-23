#!/usr/bin/env Rscript
# Golden values from igraph (Csardi and Nepusz 2006), sna (Butts 2008), and
# base R's stats::cmdscale() for reach, edge betweenness, eigenvector
# centralization, information centrality, regular equivalence, blockmodel
# densities, and classical scaling (validation/manifest.d/roles.yaml).
#
#   Rscript validation/oracles/r_roles.R
#
# Requires the igraph, sna, and yaml packages. Run from the repository root.
#
# Golden-file conventions (see validation/bindings/roles.py): an undefined
# value is an empty line, block densities are written one block row per
# line, classical scaling writes one point per line, and partitions are
# written as class numbers (compared up to relabelling). The fixtures carry
# no partitions, so block and seed partitions put node v (0-based) in class
# v %% 3 (`mod3`) or v %% 2 (`mod2`), as every oracle does. Classical
# scaling runs on hop distances of the undirected simple reading of the
# network with unreachable pairs set to the longest finite distance plus 1,
# as in the bindings; its coordinates are written with each dimension's
# sign fixed so that the first coordinate larger than 1e-8 times the
# dimension's largest magnitude is positive (the engine's rule), and only
# dimensions whose eigenvalue exceeds 1e-12 times the largest are kept.

suppressPackageStartupMessages({
  library(igraph)
  library(sna)
})

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "roles.yaml"))

fmt <- function(x) formatC(as.numeric(x), digits = 17, format = "g")

# `values` is a numeric vector (one value per line, NaN or NA as an empty
# line) or a matrix (one row per line).
write_values <- function(oracle, fixture, variant, values) {
  d <- file.path(golden_dir, oracle, fixture)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  lines <- if (is.matrix(values)) {
    if (nrow(values)) apply(values, 1, function(r) paste(fmt(r), collapse = " ")) else character(0)
  } else {
    values <- as.numeric(values)
    ifelse(is.na(values), "", fmt(values))
  }
  writeLines(c("value", lines), file.path(d, paste0(variant, ".csv")))
}

applies <- function(spec, meta) {
  kind <- if (is.null(spec$applies_to)) "any" else spec$applies_to
  if (kind == "directed" && !meta$directed) return(FALSE)
  if (kind == "undirected" && meta$directed) return(FALSE)
  if (isTRUE(spec$requires_weights) && !meta$weighted) return(FALSE)
  TRUE
}

# meta.json is flat and simple; parse it without a JSON dependency.
read_meta <- function(path) {
  txt <- paste(readLines(path), collapse = "")
  get <- function(key) sub(sprintf('.*"%s": *([^,}]+).*', key), "\\1", txt)
  list(n = as.integer(get("n")),
       directed = get("directed") == "true",
       weighted = get("weighted") == "true")
}

partition_for <- function(rule, n) {
  v <- seq_len(n) - 1L
  switch(rule, mod3 = v %% 3L, mod2 = v %% 2L, stop("unknown partition rule ", rule))
}

# Classes (numbered in order of first node) of the nodes at distance 0 in
# a CATREGE distance matrix.
zero_classes <- function(d) {
  d <- as.matrix(d)
  keys <- apply(d == 0, 1, function(r) paste(as.integer(r), collapse = ""))
  match(keys, unique(keys)) - 1L
}

# Hop distances of the undirected simple reading, unreachable pairs filled
# with the longest finite distance plus 1 (1 when no pair is connected).
mds_input <- function(g) {
  u <- igraph::simplify(igraph::as_undirected(g, mode = "collapse"))
  d <- igraph::distances(u, weights = NA)
  finite <- d[is.finite(d)]
  fill <- if (length(finite)) max(finite) + 1 else 1
  d[!is.finite(d)] <- fill
  unname(d)
}

fix_signs <- function(points) {
  for (j in seq_len(ncol(points))) {
    v <- points[, j]
    big <- max(abs(v))
    first <- v[which(abs(v) > 1e-8 * big)[1]]
    if (!is.na(first) && first < 0) points[, j] <- -v
  }
  points
}

igraph_values <- function(variant, spec, g, meta) {
  e <- spec$engine
  w <- if (meta$weighted) igraph::E(g)$weight else NULL
  quiet <- function(expr) suppressWarnings(expr)
  switch(e$fn,
    reach = {
      order <- if (identical(e$k, "inf")) igraph::vcount(g) else as.integer(e$k)
      mode <- if (is.null(e$mode)) "out" else e$mode
      igraph::ego_size(g, order = order, mode = mode) - 1
    },
    edge_betweenness = {
      weights <- if (isTRUE(e$weighted)) {
        if (identical(e$weights_as, "strength")) 1 / w else w
      } else {
        NA
      }
      igraph::edge_betweenness(g, directed = TRUE, weights = weights)
    },
    centralization_eigenvector = quiet(
      igraph::centr_eigen(g, directed = FALSE, scale = TRUE, normalized = TRUE)$centralization),
    NULL
  )
}

# `binary` is the 0/1 adjacency matrix (parallel edges collapsed, loops on
# the diagonal); `valued` holds total edge weights.
sna_values <- function(variant, spec, binary, valued, meta) {
  e <- spec$engine
  gmode <- if (meta$directed) "digraph" else "graph"
  switch(e$fn,
    information_centrality = {
      m <- if (isTRUE(e$weighted)) valued else binary
      sna::infocent(m, gmode = gmode, cmode = "weak", diag = FALSE)
    },
    regular_equivalence = {
      if (isTRUE(e$weighted)) return(NULL)
      seed <- if (is.null(e$seed)) NULL else partition_for(e$seed, meta$n) + 1L
      d <- sna::redist(binary, mode = gmode, method = "catrege", diag = FALSE,
                       seed.partition = seed)
      zero_classes(d)
    },
    block_density = {
      if (isTRUE(e$loops)) return(NULL)
      m <- if (isTRUE(e$weighted)) valued else binary
      classes <- partition_for(e$partition, meta$n) + 1L
      bm <- sna::blockmodel(m, classes, block.content = "density", mode = gmode, diag = FALSE)
      out <- unname(as.matrix(bm$block.model))
      out[is.na(out)] <- -1
      out
    },
    NULL
  )
}

stats_values <- function(variant, spec, g, meta) {
  e <- spec$engine
  if (e$fn != "classical_mds") return(NULL)
  k <- as.integer(e$k)
  n <- meta$n
  if (n <= k) return(matrix(numeric(0), 0, 0))
  d <- mds_input(g)
  fit <- suppressWarnings(stats::cmdscale(d, k = k, eig = TRUE))
  switch(e$output,
    points = {
      pts <- as.matrix(fit$points)
      floor <- 1e-12 * abs(fit$eig[1])
      keep <- which(fit$eig[seq_len(ncol(pts))] > floor)
      fix_signs(pts[, keep, drop = FALSE])
    },
    eigenvalues = fit$eig,
    gof = fit$GOF,
    NULL
  )
}

counts <- c(igraph = 0L, sna = 0L, stats = 0L)
for (name in sort(list.dirs(fixtures_dir, full.names = FALSE, recursive = FALSE))) {
  d <- file.path(fixtures_dir, name)
  meta <- read_meta(file.path(d, "meta.json"))
  e <- utils::read.csv(file.path(d, "edges.csv"))
  n <- meta$n
  g <- igraph::make_empty_graph(n = n, directed = meta$directed)
  if (nrow(e)) {
    g <- igraph::add_edges(g, as.vector(t(as.matrix(e[, c("from", "to")]))) + 1L)
    if (meta$weighted) igraph::E(g)$weight <- e$weight
  }
  binary <- matrix(0, n, n)
  valued <- matrix(0, n, n)
  if (nrow(e)) {
    i <- e$from + 1L
    j <- e$to + 1L
    w <- if (meta$weighted) e$weight else rep(1, nrow(e))
    binary[cbind(i, j)] <- 1
    if (!meta$directed) binary[cbind(j, i)] <- 1
    for (k in seq_along(i)) {
      valued[i[k], j[k]] <- valued[i[k], j[k]] + w[k]
      if (!meta$directed && i[k] != j[k]) valued[j[k], i[k]] <- valued[j[k], i[k]] + w[k]
    }
  }

  for (variant in names(manifest)) {
    spec <- manifest[[variant]]
    if (!applies(spec, meta)) next
    for (oracle in intersect(c("igraph", "sna", "stats"), names(spec$oracles))) {
      if (oracle == "sna" && n == 0) next
      v <- tryCatch(
        switch(oracle,
          igraph = igraph_values(variant, spec, g, meta),
          sna = sna_values(variant, spec, binary, valued, meta),
          stats = stats_values(variant, spec, g, meta)
        ),
        error = function(err) {
          message(sprintf("%s %s on %s failed: %s", oracle, variant, name, conditionMessage(err)))
          NULL
        }
      )
      if (!is.null(v)) {
        write_values(oracle, name, variant, v)
        counts[[oracle]] <- counts[[oracle]] + 1L
      }
    }
  }
}

writeLines(sprintf("igraph %s", utils::packageVersion("igraph")),
           file.path(golden_dir, "igraph", "VERSION"))
writeLines(sprintf("sna %s", utils::packageVersion("sna")),
           file.path(golden_dir, "sna", "VERSION"))
dir.create(file.path(golden_dir, "stats"), recursive = TRUE, showWarnings = FALSE)
writeLines(sprintf("R %s stats::cmdscale", getRversion()),
           file.path(golden_dir, "stats", "VERSION"))
cat(sprintf("Wrote %d igraph, %d sna, and %d stats golden files (roles).\n",
            counts[["igraph"]], counts[["sna"]], counts[["stats"]]))
