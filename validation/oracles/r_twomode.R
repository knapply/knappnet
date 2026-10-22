#!/usr/bin/env Rscript
# Golden values from igraph (Csardi and Nepusz 2006) for ego networks and
# two-mode networks (validation/manifest.d/twomode.yaml).
#
#   Rscript validation/oracles/r_twomode.R
#
# Requires the igraph and yaml packages. Run from the repository root.
# Undefined values (NaN or NA) are written as -1, as the engine runner does.
# Ego measures other than sizes and neighborhoods use the simple graph
# (parallel edges collapsed, self-loops dropped), as the engine does.

suppressPackageStartupMessages(library(igraph))

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "twomode.yaml"))

fmt <- function(x) formatC(x, digits = 17, format = "g")

write_values <- function(oracle, fixture, variant, values) {
  values <- as.numeric(values)
  values[is.na(values)] <- -1
  d <- file.path(golden_dir, oracle, fixture)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(c("value", if (length(values)) fmt(values)),
             file.path(d, paste0(variant, ".csv")))
}

# One group of numbers per line, space-separated (an empty line for an
# empty group), as common.write_values does for tuples.
write_groups <- function(oracle, fixture, variant, groups) {
  lines <- vapply(groups, function(x) paste(fmt(as.numeric(x)), collapse = " "), "")
  d <- file.path(golden_dir, oracle, fixture)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(c("value", lines), file.path(d, paste0(variant, ".csv")))
}

applies <- function(spec, meta) {
  kind <- if (is.null(spec$applies_to)) "any" else spec$applies_to
  if (kind == "directed" && !meta$directed) return(FALSE)
  if (kind == "undirected" && meta$directed) return(FALSE)
  if (isTRUE(spec$requires_weights) && !meta$weighted) return(FALSE)
  TRUE
}

# meta.json is JSON, which the YAML parser reads; `type` (0 or 1 per node)
# is present only for two-mode fixtures.
read_meta <- function(path) {
  m <- yaml::read_yaml(path)
  list(n = as.integer(m$n), directed = isTRUE(m$directed), weighted = isTRUE(m$weighted),
       type = if (!is.null(m$type)) as.integer(unlist(m$type)))
}

mode_of <- function(e) if (is.null(e$mode)) "all" else e$mode
order_of <- function(e) if (is.null(e$order)) 1L else as.integer(e$order)

ego_values <- function(e, g, gs, meta) {
  n <- meta$n
  ids <- seq_len(n)
  mode <- mode_of(e)
  order <- order_of(e)
  if (e$fn == "ego_size") {
    mindist <- if (is.null(e$mindist)) 0L else as.integer(e$mindist)
    return(igraph::ego_size(g, order = order, nodes = ids, mode = mode, mindist = mindist))
  }
  if (e$fn == "ego_members") {
    mindist <- if (is.null(e$mindist)) 0L else as.integer(e$mindist)
    hood <- igraph::ego(g, order = order, nodes = ids, mode = mode, mindist = mindist)
    return(lapply(hood, function(h) sort(as.integer(h)) - 1L))
  }
  vapply(ids, function(v) {
    if (e$measure == "betweenness") {
      members <- as.integer(igraph::ego(gs, order = order, nodes = v, mode = mode)[[1]])
      sub <- igraph::induced_subgraph(gs, members)
      return(as.numeric(igraph::betweenness(sub, v = as.character(v - 1L), directed = TRUE,
                                            weights = NA)))
    }
    alters <- as.integer(igraph::ego(gs, order = order, nodes = v, mode = mode,
                                     mindist = 1L)[[1]])
    sub <- igraph::induced_subgraph(gs, alters)
    as.numeric(switch(e$measure,
      ties = igraph::ecount(sub),
      density = if (length(alters) > 1L) igraph::edge_density(sub, loops = FALSE) else NA,
      components = if (length(alters)) igraph::count_components(sub, mode = "weak") else 0,
      stop("unknown ego measure ", e$measure)
    ))
  }, 0)
}

projection_groups <- function(p) {
  if (igraph::ecount(p) == 0L) return(list())
  ends <- igraph::as_edgelist(p, names = TRUE)
  a <- as.numeric(ends[, 1])
  b <- as.numeric(ends[, 2])
  w <- igraph::E(p)$weight
  lo <- pmin(a, b)
  hi <- pmax(a, b)
  o <- order(lo, hi)
  lapply(o, function(k) c(lo[k], hi[k], w[k]))
}

twomode_values <- function(variant, e, g, gs, meta) {
  types <- meta$type
  if (e$fn == "project") {
    if (!identical(e$weight, "count") || isTRUE(e$loops)) return(NULL)
    plain <- g
    if ("weight" %in% igraph::edge_attr_names(plain)) {
      plain <- igraph::delete_edge_attr(plain, "weight")
    }
    proj <- igraph::bipartite_projection(plain, types = types == 1L, multiplicity = TRUE)
    return(projection_groups(if (isTRUE(e$second)) proj$proj2 else proj$proj1))
  }
  switch(variant,
    twomode_degree_raw = igraph::degree(gs),
    twomode_closeness_raw = igraph::closeness(gs, weights = NA, normalized = FALSE),
    twomode_betweenness_raw = igraph::betweenness(gs, weights = NA),
    twomode_eigenvector = {
      v <- igraph::eigen_centrality(gs, weights = NA)$vector
      for (t in 0:1) {
        idx <- which(types == t)
        v[idx] <- v[idx] / sqrt(sum(v[idx]^2))
      }
      v
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
  igraph::V(g)$name <- as.character(seq_len(meta$n) - 1L)
  gs <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                         edge.attr.comb = "ignore")

  for (variant in names(manifest)) {
    spec <- manifest[[variant]]
    if (!("igraph" %in% names(spec$oracles)) || !applies(spec, meta)) next
    eng <- spec$engine
    if (eng$fn %in% c("ego_size", "ego_members", "ego_summary")) {
      v <- ego_values(eng, g, gs, meta)
    } else if (is.null(meta$type)) {
      next
    } else {
      v <- twomode_values(variant, eng, g, gs, meta)
    }
    if (is.null(v)) next
    if (is.list(v)) write_groups("igraph", name, variant, v)
    else write_values("igraph", name, variant, v)
    count <- count + 1L
  }
}

dir.create(file.path(golden_dir, "igraph"), recursive = TRUE, showWarnings = FALSE)
writeLines(sprintf("igraph %s", utils::packageVersion("igraph")),
           file.path(golden_dir, "igraph", "VERSION"))
cat(sprintf("Wrote %d igraph golden files (ego and two-mode).\n", count))
