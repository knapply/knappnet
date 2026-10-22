#!/usr/bin/env Rscript
# Golden values for brokerage and positions from igraph (Csardi and Nepusz
# 2006) and sna (Butts 2008).
#
#   Rscript validation/oracles/r_positions.R
#
# Requires the igraph, sna, and yaml packages. Run from the repository root.
#
# Golden-file conventions (see validation/bindings/positions.py): an
# undefined value (NaN) is an empty line, brokerage writes six counts per
# node, and structural equivalence writes one matrix row per line. The
# fixtures have no group attribute, so brokerage uses node index modulo 3
# (`groups: mod3` in the manifest), 0-based as in the engine.

suppressPackageStartupMessages({
  library(igraph)
  library(sna)
})

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "positions.yaml"))

fmt <- function(x) formatC(as.numeric(x), digits = 17, format = "g")

# `values` is a numeric vector (one value per line, NaN as an empty line)
# or a matrix (one row per line).
write_values <- function(oracle, fixture, variant, values) {
  d <- file.path(golden_dir, oracle, fixture)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  lines <- if (is.matrix(values)) {
    apply(values, 1, function(r) paste(fmt(r), collapse = " "))
  } else {
    ifelse(is.nan(values) | is.na(values), "", fmt(values))
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

groups_for <- function(rule, n) {
  if (identical(rule, "mod3")) return((seq_len(n) - 1L) %% 3L)
  stop("unknown grouping rule ", rule)
}

igraph_values <- function(spec, g) {
  e <- spec$engine
  if (e$fn != "constraint") return(NULL)
  if (isTRUE(e$weighted)) {
    return(igraph::constraint(g, weights = igraph::E(g)$weight))
  }
  if ("weight" %in% igraph::edge_attr_names(g)) g <- igraph::delete_edge_attr(g, "weight")
  igraph::constraint(g)
}

# `binary` is the 0/1 adjacency matrix (parallel edges collapsed, loops on
# the diagonal); `valued` holds total edge weights.
sna_values <- function(spec, binary, valued, meta) {
  e <- spec$engine
  switch(e$fn,
    brokerage = {
      br <- sna::brokerage(binary, cl = groups_for(e$groups, meta$n))$raw.nli
      unname(br[, c("w_I", "w_O", "b_IO", "b_OI", "b_O", "t"), drop = FALSE])
    },
    structural_equivalence = {
      m <- if (isTRUE(e$weighted)) valued else binary
      unname(as.matrix(sna::sedist(m, method = e$method)))
    },
    NULL
  )
}

counts <- c(igraph = 0L, sna = 0L)
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
    for (oracle in intersect(c("igraph", "sna"), names(spec$oracles))) {
      v <- tryCatch(
        if (oracle == "igraph") igraph_values(spec, g) else sna_values(spec, binary, valued, meta),
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
cat(sprintf("Wrote %d igraph and %d sna golden files (positions).\n",
            counts[["igraph"]], counts[["sna"]]))
