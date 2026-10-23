#!/usr/bin/env Rscript
# Golden values from the statnet packages ergm (Hunter et al. 2008) and
# network (Butts 2008) for exponential random graph models
# (validation/manifest.d/ergm.yaml): summary statistics, maximum
# pseudolikelihood estimates and their standard errors, Monte Carlo maximum
# likelihood estimates, and bridge-sampled log-likelihoods; and from btergm
# (Leifeld, Cranmer, and Desmarais 2018) for pooled temporal
# pseudolikelihood estimates, when it is installed.
#
#   Rscript validation/oracles/r_ergm.R
#
# Requires the ergm, network, and yaml packages (btergm optional). Run from
# the repository root. Networks are built from the binary sociomatrix
# (self-loops dropped, parallel edges collapsed). Node attributes follow the
# scheme of the manifest: x = i mod 5 + 1, a = i mod 3,
# w[i, j] = ((7 i + 3 j) mod 5) / 2 for 0-based i, j, with (min, max)
# ordering when undirected. The second wave of a temporal model toggles
# dyad (i, j) when (7 i + 3 j) mod 11 == 0.

suppressPackageStartupMessages({
  library(network)
  library(ergm)
})
has_btergm <- requireNamespace("btergm", quietly = TRUE)
if (has_btergm) suppressPackageStartupMessages(library(btergm))

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
  edges_kstar2 = "edges + kstar(2)",
  tergm_undirected = "edges + nodecov('x') + memory(type = 'stability')",
  tergm_directed = "edges + mutual + nodecov('x') + memory(type = 'stability') + delrecip"
)

# The sociomatrix of the second wave: dyad (i, j) (0-based, canonical order)
# toggled when (7 i + 3 j) mod 11 == 0.
second_wave <- function(m, directed) {
  n <- nrow(m)
  b <- m
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j || (!directed && j < i)) next
      if ((7 * (i - 1) + 3 * (j - 1)) %% 11 == 0) {
        b[i, j] <- 1 - b[i, j]
        if (!directed) b[j, i] <- b[i, j]
      }
    }
  }
  b
}

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
    if (!applies(spec, meta)) next
    oracle <- if ("ergm" %in% names(spec$oracles)) "ergm" else if ("btergm" %in% names(spec$oracles)) "btergm" else next
    if (name %in% unlist(spec$oracles[[oracle]]$undefined_on)) next
    eng <- spec$engine
    if (!is.null(eng$fixtures) && !(name %in% unlist(eng$fixtures))) next
    rhs <- term_sets[[eng$terms]]
    f <- stats::as.formula(paste("net ~", rhs))
    if (eng$fn == "tergm_mple") {
      if (!has_btergm) next
      net2 <- network::network(second_wave(m, meta$directed), directed = meta$directed,
                               loops = FALSE)
      network::set.vertex.attribute(net2, "x", attrs$x)
      network::set.vertex.attribute(net2, "a", attrs$a)
      waves <- list(net, net2)
      tf <- stats::as.formula(paste("waves ~", rhs))
      fit <- suppressMessages(suppressWarnings(btergm::btergm(tf, R = 10, verbose = FALSE)))
      write_values("btergm", name, variant, unname(fit@coef))
      count <- count + 1L
      next
    }
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

dir.create(file.path(golden_dir, "ergm"), recursive = TRUE, showWarnings = FALSE)
writeLines(sprintf("ergm %s", utils::packageVersion("ergm")),
           file.path(golden_dir, "ergm", "VERSION"))
if (has_btergm) {
  dir.create(file.path(golden_dir, "btergm"), recursive = TRUE, showWarnings = FALSE)
  writeLines(sprintf("btergm %s", utils::packageVersion("btergm")),
             file.path(golden_dir, "btergm", "VERSION"))
}
cat(sprintf("Wrote %d ergm and btergm golden files.\n", count))
