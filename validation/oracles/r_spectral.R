#!/usr/bin/env Rscript
# Golden values from igraph (Csardi and Nepusz 2006), sna (Butts 2008), and,
# when it is installed, the R package homophily (Knapp 2019) for PageRank,
# Katz and Bonacich power centrality, hubs and authorities, assortativity,
# mixing matrices, and the E-I index (validation/manifest.d/spectral.yaml).
#
#   Rscript validation/oracles/r_spectral.R
#
# Requires igraph, sna, and yaml; homophily is optional
# (remotes::install_github("knapply/homophily")). Run from the repository
# root. Undefined values (NaN or NA) are written as empty lines; a call that
# fails writes no golden file. Node attributes follow the rules in
# validation/bindings/spectral.py: group = v %% 3, value = 7 v %% 11,
# value_in = 3 v %% 5 for 0-based v; nodes are named "n1", "n2", ...

suppressPackageStartupMessages({
  library(igraph)
  library(sna)
})
have_homophily <- requireNamespace("homophily", quietly = TRUE)

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "spectral.yaml"))

fmt <- function(v) ifelse(is.na(v), "", formatC(v, digits = 17, format = "g"))

# A vector is written one value per line; a matrix one row per line.
write_values <- function(oracle, fixture, variant, values) {
  if (is.matrix(values)) {
    lines <- apply(values, 1L, function(r) paste(fmt(as.numeric(r)), collapse = " "))
    if (!nrow(values)) lines <- character(0)
  } else {
    lines <- fmt(as.numeric(values))
  }
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

read_meta <- function(path) {
  txt <- paste(readLines(path), collapse = "")
  get <- function(key) sub(sprintf('.*"%s": *([^,}]+).*', key), "\\1", txt)
  list(n = as.integer(get("n")),
       directed = get("directed") == "true",
       weighted = get("weighted") == "true")
}

# Largest eigenvalue modulus of a square matrix.
radius <- function(a) if (nrow(a)) max(Mod(eigen(a, only.values = TRUE)$values)) else 0

# Adjacency matrix with parallel edges added up and self-loops dropped.
counts_no_loops <- function(e, meta) {
  a <- matrix(0, meta$n, meta$n)
  for (k in seq_len(nrow(e))) {
    i <- e$from[k] + 1L
    j <- e$to[k] + 1L
    if (i == j) next
    a[i, j] <- a[i, j] + 1
    if (!meta$directed) a[j, i] <- a[j, i] + 1
  }
  a
}

hits_vector <- function(f, g) {
  v <- f(g, weights = NA)
  if (is.list(v)) v <- v$vector
  v
}

igraph_values <- function(variant, spec, g, e, meta) {
  n <- meta$n
  v0 <- 0:(n - 1)
  w <- if (meta$weighted) igraph::E(g)$weight else NULL
  quiet <- function(expr) suppressWarnings(expr)
  plain <- counts_no_loops(e, meta)
  fraction <- function(key) {
    lam <- radius(plain)
    f <- spec$engine[[paste0(key, "_fraction")]]
    if (lam > 0) f / lam else f
  }
  noloops <- igraph::simplify(g, remove.multiple = FALSE, remove.loops = TRUE)
  switch(variant,
    pagerank = igraph::page_rank(g, damping = 0.85, weights = NA)$vector,
    pagerank_weighted = igraph::page_rank(g, damping = 0.85, weights = w)$vector,
    pagerank_personalized = igraph::page_rank(g, damping = 0.85, weights = NA,
                                              personalized = v0 %% 3 + 1)$vector,
    pagerank_undirected = igraph::page_rank(g, damping = 0.85, directed = FALSE,
                                            weights = NA)$vector,
    katz_igraph = igraph::alpha_centrality(noloops, alpha = fraction("alpha"), exo = 1,
                                           weights = NA),
    katz_exogenous = igraph::alpha_centrality(noloops, alpha = fraction("alpha"),
                                              exo = v0 %% 3 + 1, weights = NA),
    bonacich = quiet(igraph::power_centrality(g, exponent = fraction("exponent"))),
    bonacich_negative = quiet(igraph::power_centrality(g, exponent = fraction("exponent"))),
    bonacich_unit = quiet(igraph::power_centrality(g, exponent = 1)),
    bonacich_rescale = quiet(igraph::power_centrality(g, exponent = fraction("exponent"),
                                                      rescale = TRUE)),
    hits_authorities = quiet(hits_vector(igraph::authority_score, g)),
    hits_hubs = quiet(hits_vector(igraph::hub_score, g)),
    assortativity_degree = igraph::assortativity_degree(g, directed = TRUE),
    assortativity_degree_undirected = igraph::assortativity_degree(g, directed = FALSE),
    assortativity_nominal = igraph::assortativity_nominal(g, types = v0 %% 3 + 1,
                                                          directed = TRUE),
    assortativity_nominal_undirected = igraph::assortativity_nominal(
      g, types = v0 %% 3 + 1, directed = FALSE),
    assortativity_modularity = igraph::assortativity_nominal(
      g, types = v0 %% 3 + 1, directed = TRUE, normalized = FALSE),
    assortativity_values = igraph::assortativity(g, values = (7 * v0) %% 11, directed = TRUE),
    assortativity_values_in = igraph::assortativity(g, values = (7 * v0) %% 11,
                                                    values.in = (3 * v0) %% 5,
                                                    directed = TRUE),
    assortativity_covariance = igraph::assortativity(g, values = (7 * v0) %% 11,
                                                     directed = TRUE, normalized = FALSE),
    NULL
  )
}

sna_values <- function(variant, m, meta) {
  gmode <- if (meta$directed) "digraph" else "graph"
  switch(variant,
    bonacich_sna = sna::bonpow(m, gmode = gmode, diag = FALSE, exponent = 1),
    NULL
  )
}

homophily_values <- function(variant, g, meta) {
  if (!have_homophily) return(NULL)
  n <- meta$n
  v0 <- 0:(n - 1)
  mat <- function(m) as.matrix(m)
  quiet <- function(expr) suppressWarnings(expr)
  switch(variant,
    assortativity_degree = homophily::assort_degree(g),
    assortativity_nominal_undirected = homophily::assort_discrete(g, "group"),
    assortativity_values = homophily::assort_scalar(g, (7 * v0) %% 11),
    assortativity_values_in = homophily::assort_scalar(g, (7 * v0) %% 11, (3 * v0) %% 5),
    mixing_matrix = mat(homophily::as_mixing_matrix(g, dim1 = "group")),
    mixing_matrix_in = mat(homophily::as_mixing_matrix(g, dim1 = "group", direction = "in")),
    mixing_matrix_nodes = mat(homophily::as_mixing_matrix(g, dim1 = "name", dim2 = "group")),
    mixing_matrix_nodes_out = mat(homophily::as_mixing_matrix(g, dim1 = "name", dim2 = "group",
                                                              direction = "out")),
    mixing_matrix_nodes_in = mat(homophily::as_mixing_matrix(g, dim1 = "name", dim2 = "group",
                                                             direction = "in")),
    ei_index = quiet(homophily::ei_index(g, "group", scope = "global")),
    ei_index_group_homophily = quiet(homophily::ei_index(g, "group", scope = "group")),
    ei_index_node_homophily = quiet(homophily::ei_index(g, "group", scope = "node")),
    NULL
  )
}

attempt <- function(expr) tryCatch(expr, error = function(err) NULL)

counts <- c(igraph = 0L, sna = 0L, homophily = 0L)
for (name in sort(list.dirs(fixtures_dir, full.names = FALSE, recursive = FALSE))) {
  d <- file.path(fixtures_dir, name)
  meta <- read_meta(file.path(d, "meta.json"))
  e <- utils::read.csv(file.path(d, "edges.csv"))
  g <- igraph::make_empty_graph(n = meta$n, directed = meta$directed)
  if (nrow(e)) {
    g <- igraph::add_edges(g, as.vector(t(as.matrix(e[, c("from", "to")]))) + 1L)
    if (meta$weighted) igraph::E(g)$weight <- e$weight
  }
  if (meta$n > 0) {
    v0 <- 0:(meta$n - 1)
    igraph::V(g)$name <- paste0("n", seq_len(meta$n))
    igraph::V(g)$group <- v0 %% 3
  }
  m <- matrix(0, meta$n, meta$n)
  if (nrow(e)) {
    m[cbind(e$from + 1L, e$to + 1L)] <- 1
    if (!meta$directed) m[cbind(e$to + 1L, e$from + 1L)] <- 1
  }
  diag(m) <- 0

  for (variant in names(manifest)) {
    spec <- manifest[[variant]]
    if (!applies(spec, meta) || meta$n == 0) next
    for (oracle in intersect(names(spec$oracles), c("igraph", "sna", "homophily"))) {
      v <- attempt(switch(oracle,
        igraph = igraph_values(variant, spec, g, e, meta),
        sna = sna_values(variant, m, meta),
        homophily = homophily_values(variant, g, meta)
      ))
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
if (have_homophily) {
  dir.create(file.path(golden_dir, "homophily"), showWarnings = FALSE)
  writeLines(sprintf("homophily %s", utils::packageVersion("homophily")),
             file.path(golden_dir, "homophily", "VERSION"))
}
cat(sprintf("Wrote %d igraph, %d sna, and %d homophily golden files (spectral).\n",
            counts[["igraph"]], counts[["sna"]], counts[["homophily"]]))
