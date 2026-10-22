#!/usr/bin/env Rscript
# Golden values for cohesive subgroups from igraph (Csardi and Nepusz 2006).
#
#   Rscript validation/oracles/r_subgroups.R
#
# Requires the igraph and yaml packages. Run from the repository root.
# Cliques are written one per line as 0-based node indices.

suppressPackageStartupMessages(library(igraph))

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "subgroups.yaml"))

fmt <- function(x) formatC(as.numeric(x), digits = 17, format = "g")

write_values <- function(fixture, variant, values) {
  d <- file.path(golden_dir, "igraph", fixture)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  lines <- if (is.list(values)) {
    vapply(values, function(v) paste(fmt(v), collapse = " "), character(1))
  } else {
    fmt(values)
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

# The fixed partition used to validate modularity: node v (0-based) of n
# goes to group floor(3 v / n). igraph wants labels starting at 1.
blocks <- function(n) (3L * (seq_len(n) - 1L)) %/% n + 1L

edge_weights <- function(g, weighted) {
  if (weighted) igraph::E(g)$weight else rep(1, igraph::ecount(g))
}

cliques0 <- function(g, min_size) {
  cl <- suppressWarnings(igraph::max_cliques(g, min = min_size))
  lapply(cl, function(c) sort(as.integer(c)) - 1L)
}

igraph_values <- function(spec, g) {
  e <- spec$engine
  res <- if (is.null(e$resolution)) 1 else e$resolution
  weighted <- isTRUE(e$weighted)
  switch(e$fn,
    coreness = igraph::coreness(g, mode = e$mode),
    max_cliques = cliques0(g, e$min_size),
    clique_number = suppressWarnings(igraph::clique_num(g)),
    modularity = igraph::modularity(g, blocks(igraph::vcount(g)),
                                    weights = edge_weights(g, weighted),
                                    resolution = res, directed = TRUE),
    louvain_quality = {
      w <- edge_weights(g, weighted)
      q <- vapply(seq_len(e$seeds), function(s) {
        set.seed(s)
        cl <- igraph::cluster_louvain(g, weights = w, resolution = res)
        igraph::modularity(g, igraph::membership(cl), weights = w, resolution = res)
      }, numeric(1))
      mean(q)
    },
    NULL
  )
}

count <- 0L
for (name in sort(list.dirs(fixtures_dir, full.names = FALSE, recursive = FALSE))) {
  d <- file.path(fixtures_dir, name)
  meta <- read_meta(file.path(d, "meta.json"))
  e <- utils::read.csv(file.path(d, "edges.csv"))
  g <- igraph::make_empty_graph(n = meta$n, directed = meta$directed)
  if (nrow(e)) {
    g <- igraph::add_edges(g, as.vector(t(as.matrix(e[, c("from", "to")]))) + 1L)
    if (meta$weighted) igraph::E(g)$weight <- e$weight
  }
  for (variant in names(manifest)) {
    spec <- manifest[[variant]]
    if (!("igraph" %in% names(spec$oracles)) || !applies(spec, meta)) next
    v <- igraph_values(spec, g)
    if (!is.null(v)) {
      write_values(name, variant, v)
      count <- count + 1L
    }
  }
}

dir.create(file.path(golden_dir, "igraph"), recursive = TRUE, showWarnings = FALSE)
writeLines(sprintf("igraph %s", utils::packageVersion("igraph")),
           file.path(golden_dir, "igraph", "VERSION"))
cat(sprintf("Wrote %d igraph golden files for subgroups.\n", count))
