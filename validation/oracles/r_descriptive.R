#!/usr/bin/env Rscript
# Golden values from igraph (Csardi and Nepusz 2006) and sna (Butts 2008)
# for reciprocity, transitivity and clustering, distances, mean distance,
# diameter, eccentricity, radius, reachability, cut points, and bridges
# (validation/manifest.d/descriptive.yaml).
#
#   Rscript validation/oracles/r_descriptive.R
#
# Requires the igraph, sna, and yaml packages. Run from the repository root.
# NaN and NA are written as -999 and infinity as -1000, as the engine runner
# does. Where a measure is 0 / 0 on a fixture, nothing is written, and the
# manifest declares the fixture with undefined_on.

suppressPackageStartupMessages({
  library(igraph)
  library(sna)
})

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "descriptive.yaml"))

encode <- function(v) {
  v <- as.numeric(v)
  v[is.na(v)] <- -999
  v[is.infinite(v)] <- -1000
  formatC(v, digits = 17, format = "g")
}

# A vector is written one value per line; a matrix one row per line.
write_values <- function(oracle, fixture, variant, values) {
  d <- file.path(golden_dir, oracle, fixture)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  lines <- if (is.matrix(values)) {
    vapply(seq_len(nrow(values)), function(i) paste(encode(values[i, ]), collapse = " "),
           character(1))
  } else {
    encode(values)
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

# as_undirected() replaced as.undirected() in igraph 2.1.
to_undirected <- function(g, mode) {
  f <- if ("as_undirected" %in% getNamespaceExports("igraph")) {
    getExportedValue("igraph", "as_undirected")
  } else {
    getExportedValue("igraph", "as.undirected")
  }
  f(g, mode = mode)
}

# Endpoints of edges as rows (smaller first), 0-based.
edge_pairs <- function(g, edges) {
  if (length(edges) == 0L) return(matrix(numeric(0), 0L, 2L))
  e <- igraph::ends(g, edges, names = FALSE) - 1
  cbind(pmin(e[, 1], e[, 2]), pmax(e[, 1], e[, 2]))
}

no_triples <- function(su) sum(choose(igraph::degree(su), 2)) == 0

igraph_values <- function(variant, g, gu, meta) {
  w <- if (meta$weighted) igraph::E(g)$weight else NULL
  # The simple undirected reading, without weights.
  su <- igraph::simplify(to_undirected(gu, "collapse"))
  simple <- igraph::simplify(gu)
  everything <- igraph::vcount(gu)
  switch(variant,
    reciprocity = if (igraph::ecount(simple) > 0)
      igraph::reciprocity(simple, mode = "default"),
    reciprocity_ratio = if (igraph::ecount(simple) > 0)
      igraph::reciprocity(simple, mode = "ratio"),
    transitivity = if (!no_triples(su)) igraph::transitivity(su, type = "global"),
    clustering = igraph::transitivity(su, type = "local", isolates = "NaN"),
    clustering_zero = igraph::transitivity(su, type = "local", isolates = "zero"),
    average_clustering = igraph::transitivity(su, type = "average", isolates = "NaN"),
    average_clustering_zero = igraph::transitivity(su, type = "average", isolates = "zero"),
    clustering_barrat = igraph::transitivity(g, type = "barrat", weights = w,
                                             isolates = "NaN"),
    distances = igraph::distances(gu, mode = "out"),
    distances_in = igraph::distances(gu, mode = "in"),
    distances_all = igraph::distances(gu, mode = "all"),
    distances_weighted = igraph::distances(g, mode = "out", weights = w),
    distances_weighted_strength = igraph::distances(g, mode = "out", weights = 1 / w),
    mean_distance = igraph::mean_distance(gu, directed = TRUE, unconnected = TRUE),
    mean_distance_all = igraph::mean_distance(gu, directed = FALSE, unconnected = TRUE),
    mean_distance_connected = igraph::mean_distance(gu, directed = TRUE, unconnected = FALSE),
    mean_distance_weighted = igraph::mean_distance(g, weights = w, directed = TRUE,
                                                   unconnected = TRUE),
    diameter = igraph::diameter(gu, directed = TRUE, unconnected = TRUE),
    diameter_connected = igraph::diameter(gu, directed = TRUE, unconnected = FALSE),
    diameter_all = igraph::diameter(gu, directed = FALSE, unconnected = TRUE),
    diameter_weighted = igraph::diameter(g, directed = TRUE, unconnected = TRUE, weights = w),
    eccentricity = igraph::eccentricity(gu, mode = "all"),
    eccentricity_out = igraph::eccentricity(gu, mode = "out"),
    radius = igraph::radius(gu, mode = "all"),
    reachability = igraph::ego_size(gu, order = everything, mode = "out") - 1,
    reachability_in = igraph::ego_size(gu, order = everything, mode = "in") - 1,
    reachability_all = igraph::ego_size(gu, order = everything, mode = "all") - 1,
    cut_points = sort(as.integer(igraph::articulation_points(su))) - 1,
    bridges = {
      u <- to_undirected(gu, "each")
      # Compared only where the undirected reading has no parallel edges
      # (see the manifest).
      if (!any(igraph::which_multiple(u) & !igraph::which_loop(u))) {
        edge_pairs(u, igraph::bridges(u))
      }
    },
    NULL
  )
}

sna_values <- function(variant, m0, msym) {
  switch(variant,
    reciprocity = if (sum(m0) > 0) sna::grecip(m0, measure = "edgewise"),
    reciprocity_ratio = if (sum(m0) > 0) sna::grecip(m0, measure = "dyadic.nonnull"),
    reciprocity_dyadic = sna::grecip(m0, measure = "dyadic"),
    reciprocity_correlation = {
      off <- m0[row(m0) != col(m0)]
      # 0 / 0 when every value is the same (sna then reports 1); skipped.
      if (length(unique(off)) > 1L) sna::grecip(m0, measure = "correlation")
    },
    distances = sna::geodist(m0)$gdist,
    cut_points = sort(as.integer(sna::cutpoints(msym, mode = "graph"))) - 1,
    NULL
  )
}

# An oracle that fails on one fixture is reported and skipped, so the
# comparison shows it as not run instead of losing every other golden file.
attempt <- function(oracle, variant, name, expr) {
  tryCatch(expr, error = function(e) {
    message(sprintf("%s %s on %s failed: %s", oracle, variant, name, conditionMessage(e)))
    NULL
  })
}

counts <- c(igraph = 0L, sna = 0L)
for (name in sort(list.dirs(fixtures_dir, full.names = FALSE, recursive = FALSE))) {
  d <- file.path(fixtures_dir, name)
  meta <- read_meta(file.path(d, "meta.json"))
  e <- utils::read.csv(file.path(d, "edges.csv"))
  g <- igraph::make_empty_graph(n = meta$n, directed = meta$directed)
  if (nrow(e)) {
    g <- igraph::add_edges(g, as.vector(t(as.matrix(e[, c("from", "to")]))) + 1L)
    if (meta$weighted) igraph::E(g)$weight <- e$weight
  }
  # The same network without weights, for unweighted measures, so that no
  # igraph function picks up the weight attribute on its own.
  gu <- if (meta$weighted) igraph::delete_edge_attr(g, "weight") else g
  # The adjacency matrix without its diagonal, and its symmetrization.
  m0 <- matrix(0, meta$n, meta$n)
  if (nrow(e)) {
    m0[cbind(e$from + 1L, e$to + 1L)] <- 1
    if (!meta$directed) m0[cbind(e$to + 1L, e$from + 1L)] <- 1
  }
  diag(m0) <- 0
  msym <- (m0 + t(m0) > 0) * 1

  for (variant in names(manifest)) {
    spec <- manifest[[variant]]
    if (!applies(spec, meta)) next
    if ("igraph" %in% names(spec$oracles)) {
      v <- attempt("igraph", variant, name, igraph_values(variant, g, gu, meta))
      if (!is.null(v)) {
        write_values("igraph", name, variant, v)
        counts[["igraph"]] <- counts[["igraph"]] + 1L
      }
    }
    if ("sna" %in% names(spec$oracles) && meta$n > 0) {
      v <- attempt("sna", variant, name, sna_values(variant, m0, msym))
      if (!is.null(v)) {
        write_values("sna", name, variant, v)
        counts[["sna"]] <- counts[["sna"]] + 1L
      }
    }
  }
}

writeLines(sprintf("igraph %s", utils::packageVersion("igraph")),
           file.path(golden_dir, "igraph", "VERSION"))
writeLines(sprintf("sna %s", utils::packageVersion("sna")),
           file.path(golden_dir, "sna", "VERSION"))
cat(sprintf("Wrote %d igraph and %d sna golden files (descriptive).\n",
            counts[["igraph"]], counts[["sna"]]))
