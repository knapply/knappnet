#!/usr/bin/env Rscript
# Golden values from igraph (Csardi and Nepusz 2006) and sna (Butts 2008).
#
#   Rscript validation/oracles/r_oracles.R
#
# Requires the igraph, sna, and yaml packages. Run from the repository root.

suppressPackageStartupMessages({
  library(igraph)
  library(sna)
})

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.yaml"))

write_values <- function(oracle, fixture, variant, values) {
  d <- file.path(golden_dir, oracle, fixture)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(c("value", formatC(as.numeric(values), digits = 17, format = "g")),
             file.path(d, paste0(variant, ".csv")))
}

applies <- function(spec, meta) {
  kind <- if (is.null(spec$applies_to)) "any" else spec$applies_to
  if (kind == "directed" && !meta$directed) return(FALSE)
  if (kind == "undirected" && meta$directed) return(FALSE)
  if (isTRUE(spec$requires_weights) && !meta$weighted) return(FALSE)
  TRUE
}

load_fixture <- function(name) {
  d <- file.path(fixtures_dir, name)
  meta <- jsonlite_free_read(file.path(d, "meta.json"))
  edges <- utils::read.csv(file.path(d, "edges.csv"))
  list(meta = meta, edges = edges)
}

# meta.json is flat and simple; parse it without a JSON dependency.
jsonlite_free_read <- function(path) {
  txt <- paste(readLines(path), collapse = "")
  get <- function(key) sub(sprintf('.*"%s": *([^,}]+).*', key), "\\1", txt)
  list(n = as.integer(get("n")),
       directed = get("directed") == "true",
       weighted = get("weighted") == "true")
}

igraph_values <- function(variant, g, meta) {
  directed <- meta$directed
  switch(variant,
    degree_all = igraph::degree(g, mode = "all", loops = TRUE),
    degree_all_noloops = igraph::degree(g, mode = "all", loops = FALSE),
    degree_out = igraph::degree(g, mode = "out", loops = TRUE),
    degree_in = igraph::degree(g, mode = "in", loops = TRUE),
    strength_all = igraph::strength(g, mode = "all", loops = TRUE,
                                    weights = igraph::E(g)$weight),
    betweenness = igraph::betweenness(g, directed = directed, weights = NA,
                                      normalized = FALSE),
    betweenness_normalized = igraph::betweenness(g, directed = directed, weights = NA,
                                                 normalized = TRUE),
    betweenness_weighted = igraph::betweenness(g, directed = directed,
                                               weights = igraph::E(g)$weight,
                                               normalized = FALSE),
    betweenness_weighted_strength = igraph::betweenness(g, directed = directed,
                                                        weights = 1 / igraph::E(g)$weight,
                                                        normalized = FALSE),
    triad_census = suppressWarnings(igraph::triad_census(g)),
    NULL
  )
}

sna_values <- function(variant, m, meta) {
  gmode <- if (meta$directed) "digraph" else "graph"
  switch(variant,
    degree_all = sna::degree(m, gmode = gmode,
                             cmode = if (meta$directed) "freeman" else "indegree"),
    degree_all_noloops = sna::degree(m, gmode = gmode, diag = FALSE,
                                     cmode = if (meta$directed) "freeman" else "indegree"),
    degree_out = sna::degree(m, gmode = gmode, cmode = "outdegree"),
    degree_in = sna::degree(m, gmode = gmode, cmode = "indegree"),
    betweenness = sna::betweenness(m, gmode = gmode,
                                   cmode = if (meta$directed) "directed" else "undirected"),
    triad_census = if (meta$directed) drop(sna::triad.census(m, mode = "digraph")) else NULL,
    NULL
  )
}

counts <- c(igraph = 0L, sna = 0L)
for (name in sort(list.dirs(fixtures_dir, full.names = FALSE, recursive = FALSE))) {
  fx <- load_fixture(name)
  meta <- fx$meta
  e <- fx$edges
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
cat(sprintf("Wrote %d igraph and %d sna golden files.\n", counts[["igraph"]], counts[["sna"]]))
