#!/usr/bin/env Rscript
# Golden values from igraph (Csardi and Nepusz 2006) and sna (Butts 2008)
# for density, components, closeness, harmonic centrality, eigenvector
# centrality, and centralization (validation/manifest.d/topography.yaml).
#
#   Rscript validation/oracles/r_topography.R
#
# Requires the igraph, sna, and yaml packages. Run from the repository root.
# Undefined values (NaN or NA) are written as -1, as the engine runner does.

suppressPackageStartupMessages({
  library(igraph)
  library(sna)
})

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "topography.yaml"))

write_values <- function(oracle, fixture, variant, values) {
  values <- as.numeric(values)
  values[is.na(values)] <- -1
  d <- file.path(golden_dir, oracle, fixture)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(c("value", formatC(values, digits = 17, format = "g")),
             file.path(d, paste0(variant, ".csv")))
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

igraph_values <- function(variant, g, meta) {
  w <- if (meta$weighted) igraph::E(g)$weight else NULL
  quiet <- function(expr) suppressWarnings(expr)
  switch(variant,
    density = igraph::edge_density(g, loops = FALSE),
    density_loops = igraph::edge_density(g, loops = TRUE),
    components_weak = igraph::components(g, mode = "weak")$membership,
    components_strong = igraph::components(g, mode = "strong")$membership,
    closeness = igraph::closeness(g, mode = "out", weights = NA, normalized = FALSE),
    closeness_normalized = igraph::closeness(g, mode = "out", weights = NA, normalized = TRUE),
    closeness_in = igraph::closeness(g, mode = "in", weights = NA, normalized = TRUE),
    closeness_all = igraph::closeness(g, mode = "all", weights = NA, normalized = TRUE),
    closeness_weighted = igraph::closeness(g, mode = "out", weights = w, normalized = TRUE),
    closeness_weighted_strength = igraph::closeness(g, mode = "out", weights = 1 / w,
                                                    normalized = TRUE),
    harmonic = igraph::harmonic_centrality(g, mode = "out", weights = NA, normalized = FALSE),
    harmonic_normalized = igraph::harmonic_centrality(g, mode = "out", weights = NA,
                                                      normalized = TRUE),
    harmonic_in = igraph::harmonic_centrality(g, mode = "in", weights = NA, normalized = FALSE),
    harmonic_weighted = igraph::harmonic_centrality(g, mode = "out", weights = w,
                                                    normalized = FALSE),
    eigenvector = quiet(igraph::eigen_centrality(g, directed = TRUE, weights = NA)$vector),
    eigenvector_weighted = quiet(igraph::eigen_centrality(g, directed = TRUE,
                                                          weights = w)$vector),
    eigenvector_undirected = quiet(igraph::eigen_centrality(g, directed = FALSE,
                                                            weights = NA)$vector),
    centralization_degree = igraph::centr_degree(g, mode = "all", loops = TRUE,
                                                 normalized = TRUE)$centralization,
    centralization_degree_in = igraph::centr_degree(g, mode = "in", loops = TRUE,
                                                    normalized = TRUE)$centralization,
    centralization_betweenness = igraph::centr_betw(g, directed = TRUE,
                                                    normalized = TRUE)$centralization,
    centralization_closeness = quiet(igraph::centr_clo(g, mode = "out",
                                                       normalized = TRUE)$centralization),
    NULL
  )
}

sna_values <- function(variant, m, meta) {
  gmode <- if (meta$directed) "digraph" else "graph"
  path_mode <- if (meta$directed) "directed" else "undirected"
  switch(variant,
    density = sna::gden(m, mode = gmode, diag = FALSE),
    density_simple = sna::gden(m, mode = gmode, diag = FALSE),
    degree_sna = sna::degree(m, gmode = gmode,
                             cmode = if (meta$directed) "freeman" else "indegree"),
    betweenness_networkx = sna::betweenness(m, gmode = gmode, cmode = path_mode),
    closeness_sna = sna::closeness(m, gmode = gmode, cmode = path_mode),
    harmonic_normalized = sna::closeness(m, gmode = gmode,
                                         cmode = if (meta$directed) "suminvdir" else "suminvundir"),
    eigenvector_sna = {
      if (meta$directed) {
        v <- suppressWarnings(sna::evcent(m, gmode = "digraph", diag = FALSE))
      } else {
        v <- sna::evcent(m, gmode = "graph", diag = FALSE, use.eigen = TRUE)
      }
      v <- Re(v)
      # sna returns NA where the leading eigenvector is undefined.
      if (anyNA(v) || all(v == 0)) return(NULL)
      if (sum(v) < 0) v <- -v
      v / sqrt(sum(v^2))
    },
    centralization_degree_sna = sna::centralization(
      m, sna::degree, mode = gmode, cmode = if (meta$directed) "freeman" else "indegree"),
    centralization_betweenness = sna::centralization(m, sna::betweenness, mode = gmode,
                                                     cmode = path_mode),
    centralization_closeness_sna = sna::centralization(m, sna::closeness, mode = gmode,
                                                       cmode = path_mode),
    NULL
  )
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
  m <- matrix(0, meta$n, meta$n)
  if (nrow(e)) {
    m[cbind(e$from + 1L, e$to + 1L)] <- 1
    if (!meta$directed) m[cbind(e$to + 1L, e$from + 1L)] <- 1
  }

  for (variant in names(manifest)) {
    spec <- manifest[[variant]]
    if (!applies(spec, meta)) next
    if ("igraph" %in% names(spec$oracles)) {
      v <- igraph_values(variant, g, meta)
      if (!is.null(v)) {
        write_values("igraph", name, variant, v)
        counts[["igraph"]] <- counts[["igraph"]] + 1L
      }
    }
    if ("sna" %in% names(spec$oracles) && meta$n > 0) {
      v <- sna_values(variant, m, meta)
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
cat(sprintf("Wrote %d igraph and %d sna golden files (topography).\n",
            counts[["igraph"]], counts[["sna"]]))
