#!/usr/bin/env Rscript
# Golden values from the statnet packages ergm (Hunter et al. 2008) and
# network (Butts 2008) for exponential random graph models
# (validation/manifest.d/ergm.yaml): summary statistics, maximum
# pseudolikelihood estimates and their standard errors, Monte Carlo maximum
# likelihood estimates, and bridge-sampled log-likelihoods.
#
#   Rscript validation/oracles/r_ergm.R
#
# Requires the ergm, network, and yaml packages. Run from the repository
# root. Networks are built from the binary sociomatrix (self-loops dropped,
# parallel edges collapsed). Node attributes follow the scheme of the
# manifest: x = i mod 5 + 1, a = i mod 3, w[i, j] = ((7 i + 3 j) mod 5) / 2
# for 0-based i, j, with (min, max) ordering when undirected.

suppressPackageStartupMessages({
  library(network)
  library(ergm)
})

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "ergm.yaml"))

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

# meta.json is flat and simple; parse it without a JSON dependency.
read_meta <- function(path) {
  txt <- paste(readLines(path), collapse = "")
  get <- function(key) sub(sprintf('.*"%s": *([^,}]+).*', key), "\\1", txt)
  list(n = as.integer(get("n")),
       directed = get("directed") == "true",
       weighted = get("weighted") == "true")
}

attributes_for <- function(n, directed) {
  i <- seq_len(n) - 1L
  x <- i %% 5 + 1
  a <- i %% 3
  w <- matrix(0, n, n)
  for (r in seq_len(n)) {
    for (c in seq_len(n)) {
      lo <- if (directed) r - 1L else min(r, c) - 1L
      hi <- if (directed) c - 1L else max(r, c) - 1L
      w[r, c] <- ((7 * lo + 3 * hi) %% 5) / 2
    }
  }
  list(x = x, a = a, w = w)
}

# Right-hand sides of the formulas for each term set named in the manifest.
# `w` is found in the calling environment when the formula is evaluated.
term_sets <- list(
  dyadic = "edges + nodecov('x') + nodefactor('a') + nodematch('a') + nodematch('a', diff = TRUE) + absdiff('x') + absdiff('x', pow = 2) + edgecov(w)",
  undirected = "kstar(2) + kstar(3) + triangle + gwesp(0.5, fixed = TRUE) + gwdsp(0.3, fixed = TRUE) + gwnsp(0.3, fixed = TRUE) + gwdegree(0.5, fixed = TRUE) + degree(1) + degree(2) + isolates",
  directed = "mutual + istar(2) + ostar(2) + ttriple + ctriple + triangle + gwesp(0.5, fixed = TRUE) + gwdsp(0.3, fixed = TRUE) + gwnsp(0.3, fixed = TRUE) + gwidegree(0.5, fixed = TRUE) + gwodegree(0.5, fixed = TRUE) + idegree(1) + odegree(2) + isolates + nodeicov('x') + nodeocov('x') + nodeifactor('a') + nodeofactor('a')",
  mple_dyadic = "edges + nodecov('x') + nodematch('a') + absdiff('x')",
  mple_gwesp = "edges + gwesp(0.5, fixed = TRUE)",
  mple_directed = "edges + mutual + gwesp(0.5, fixed = TRUE)",
  edges = "edges",
  edges_triangle = "edges + triangle",
  edges_gwesp = "edges + gwesp(0.5, fixed = TRUE)",
  edges_kstar2 = "edges + kstar(2)"
)

count <- 0L
for (name in sort(list.dirs(fixtures_dir, full.names = FALSE, recursive = FALSE))) {
  d <- file.path(fixtures_dir, name)
  meta <- read_meta(file.path(d, "meta.json"))
  e <- utils::read.csv(file.path(d, "edges.csv"))
  n <- meta$n
  m <- matrix(0, n, n)
  if (nrow(e)) {
    keep <- e$from != e$to
    m[cbind(e$from[keep] + 1L, e$to[keep] + 1L)] <- 1
    if (!meta$directed) m[cbind(e$to[keep] + 1L, e$from[keep] + 1L)] <- 1
  }
  attrs <- attributes_for(n, meta$directed)
  w <- attrs$w
  net <- network::network(m, directed = meta$directed, loops = FALSE)
  network::set.vertex.attribute(net, "x", attrs$x)
  network::set.vertex.attribute(net, "a", attrs$a)

  for (variant in names(manifest)) {
    spec <- manifest[[variant]]
    if (!applies(spec, meta) || !("ergm" %in% names(spec$oracles))) next
    if (name %in% unlist(spec$oracles$ergm$undefined_on)) next
    eng <- spec$engine
    if (!is.null(eng$fixtures) && !(name %in% unlist(eng$fixtures))) next
    rhs <- term_sets[[eng$terms]]
    f <- stats::as.formula(paste("net ~", rhs))
    values <- switch(eng$fn,
      ergm_stats = unname(summary(f)),
      ergm_mple = {
        fit <- suppressMessages(ergm::ergm(f, estimate = "MPLE"))
        if (identical(eng$output, "se")) sqrt(diag(stats::vcov(fit))) else stats::coef(fit)
      },
      ergm_mcmle = {
        fit <- suppressMessages(ergm::ergm(
          f, control = ergm::control.ergm(seed = 1, MCMC.samplesize = 10000,
                                          MCMC.burnin = 20000, MCMC.interval = 200)))
        if (identical(eng$output, "loglik")) as.numeric(stats::logLik(fit)) else stats::coef(fit)
      },
      NULL
    )
    if (is.null(values)) next
    write_values("ergm", name, variant, unname(values))
    count <- count + 1L
  }
}

writeLines(sprintf("ergm %s", utils::packageVersion("ergm")),
           file.path(golden_dir, "ergm", "VERSION"))
cat(sprintf("Wrote %d ergm golden files.\n", count))
