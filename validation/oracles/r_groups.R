#!/usr/bin/env Rscript
# Golden values from igraph (Csardi and Nepusz 2006) for bicomponents, edge
# connectivity, Girvan-Newman, leading eigenvector, and greedy modularity
# partitions, compactness, and the small-world inputs
# (validation/manifest.d/groups.yaml).
#
#   Rscript validation/oracles/r_groups.R
#
# Requires the igraph and yaml packages. Run from the repository root.
# Undefined values (NaN or NA) are written as -1, as the engine runner does.
# Node sets and partitions use 0-based node positions.

suppressPackageStartupMessages(library(igraph))

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "groups.yaml"))

fmt <- function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- -1
  formatC(x, digits = 17, format = "g")
}

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

# meta.json is flat and simple; parse it without a JSON dependency. The
# tags (a JSON array of strings) are read for the manifest's undefined_when.
read_meta <- function(path) {
  txt <- paste(readLines(path), collapse = "")
  get <- function(key) sub(sprintf('.*"%s": *([^,}]+).*', key), "\\1", txt)
  tags <- regmatches(txt, regexec('"tags": *\\[([^]]*)\\]', txt))[[1]]
  tags <- if (length(tags) == 2L) {
    regmatches(tags[[2]], gregexpr('"[^"]*"', tags[[2]]))[[1]]
  } else {
    character()
  }
  list(n = as.integer(get("n")),
       directed = get("directed") == "true",
       weighted = get("weighted") == "true",
       tags = gsub('"', "", tags))
}

# The undirected simple graph underlying g (directions, loops, and edge
# multiplicity dropped).
simple_undirected <- function(g) igraph::simplify(igraph::as_undirected(g, mode = "collapse"))

# The simple graph with parallel edges counted in a weight attribute, for
# routines that refuse multigraphs. Loops are dropped too; they would count
# toward degree and modularity in the engine, so fixtures with loops are
# listed as undefined for this oracle in the manifest.
collapse_weights <- function(g, weighted) {
  igraph::E(g)$weight <- if (weighted) igraph::E(g)$weight else 1
  igraph::simplify(g, remove.loops = TRUE, edge.attr.comb = list(weight = "sum"))
}

connectivity <- function(g) {
  n <- igraph::vcount(g)
  m <- matrix(0, n, n)
  if (n < 2) return(m)
  for (a in seq_len(n - 1)) {
    for (b in (a + 1):n) {
      m[a, b] <- m[b, a] <- igraph::edge_connectivity(g, source = a, target = b)
    }
  }
  m
}

mean_reciprocal <- function(g, weights) {
  n <- igraph::vcount(g)
  if (n < 2) return(NaN)
  d <- igraph::distances(g, mode = "out", weights = weights)
  mean((1 / d)[row(d) != col(d)])
}

# igraph uses an edge attribute called `weight` unless told otherwise, so
# unweighted variants pass weights = NA explicitly or drop the attribute.
unweighted <- function(g) {
  if ("weight" %in% igraph::edge_attr_names(g)) g <- igraph::delete_edge_attr(g, "weight")
  g
}

igraph_values <- function(spec, g, meta) {
  e <- spec$engine
  su <- simple_undirected(g)
  switch(e$fn,
    bicomponents = lapply(igraph::biconnected_components(su)$components,
                          function(c) sort(as.integer(c)) - 1L),
    lambda_connectivity = {
      m <- connectivity(su)
      m[upper.tri(m)][order(row(m)[upper.tri(m)], col(m)[upper.tri(m)])]
    },
    girvan_newman = {
      if (igraph::ecount(su) == 0) return(NaN)
      cl <- igraph::cluster_edge_betweenness(su, weights = NA, directed = FALSE)
      max(cl$modularity)
    },
    leading_eigen = {
      # Without weights the weight attribute is removed, so that no igraph
      # routine can pick it up.
      gl <- if (isTRUE(e$weighted)) g else unweighted(g)
      w <- if (isTRUE(e$weighted)) igraph::E(g)$weight else NULL
      cl <- tryCatch(igraph::cluster_leading_eigen(gl, weights = w), error = function(err) NULL)
      if (is.null(cl)) return(NULL)
      if (identical(e$output, "membership")) {
        as.integer(igraph::membership(cl)) - 1L
      } else {
        igraph::modularity(gl, igraph::membership(cl), weights = w)
      }
    },
    fast_greedy = {
      sg <- collapse_weights(g, isTRUE(e$weighted))
      if (igraph::ecount(sg) == 0) return(NaN)
      cl <- igraph::cluster_fast_greedy(sg, weights = igraph::E(sg)$weight)
      igraph::modularity(sg, igraph::membership(cl), weights = igraph::E(sg)$weight)
    },
    compactness = mean_reciprocal(g, if (isTRUE(e$weighted)) igraph::E(g)$weight else NA),
    small_world_observed = {
      c <- if (isTRUE(e$local)) {
        igraph::transitivity(su, type = "average", isolates = "zero")
      } else {
        igraph::transitivity(su, type = "global")
      }
      if (igraph::vcount(su) == 0) c <- NaN
      c(c, igraph::mean_distance(su, weights = NA, directed = FALSE, unconnected = TRUE))
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
    ora <- spec$oracles$igraph
    if (name %in% ora$undefined_on || any(meta$tags %in% ora$undefined_when)) next
    v <- igraph_values(spec, g, meta)
    if (!is.null(v)) {
      write_values(name, variant, v)
      count <- count + 1L
    }
  }
}

writeLines(sprintf("igraph %s", utils::packageVersion("igraph")),
           file.path(golden_dir, "igraph", "VERSION"))
cat(sprintf("Wrote %d igraph golden files (groups).\n", count))
