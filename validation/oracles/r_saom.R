#!/usr/bin/env Rscript
# Golden values from RSiena (Ripley, Snijders, Boda, Voros, and Preciado)
# for the deterministic parts of the stochastic actor-oriented model: the
# target statistics of every effect and the rate statistics (number of
# changed tie variables) per period (validation/manifest.d/saom.yaml,
# variants with an `rsiena` oracle).
#
#   Rscript validation/oracles/r_saom.R
#
# Requires the RSiena and yaml packages; exits quietly when RSiena is not
# installed. Run from the repository root.
#
# Every directed fixture is read as the three-wave panel of
# validation/saom_fixtures.py, rebuilt here with the same rules: wave 1 is
# the fixture's adjacency matrix without self-loops; wave 2 flips tie
# (i, j) (0-based) when (37 (i n + j) + 11) mod 7 == 0; wave 3 flips tie
# (i, j) of wave 2 when (53 (i n + j) + 5) mod 9 == 0. Covariates of actor
# i: `constant` = (7 i) mod 5, `changing` = (3 i + 2 w) mod 4 in wave w
# (w = 0, 1, 2), `binary` = i mod 2, and the dyadic `dyadic` =
# ((5 i + 3 j) mod 7) / 2.
#
# RSiena is used as a black box: for each effect, siena07() runs with
# nsub = 0 and a few phase 3 iterations, unconditional (cond = FALSE), and
# its `targets2` (observed statistics by period) gives the golden values,
# so centering of covariates and similarity means are RSiena's own. Rate
# rows hold the number of changed tie variables of their period.

if (!requireNamespace("RSiena", quietly = TRUE)) {
  cat("RSiena is not installed; no RSiena golden files written (saom).\n")
  quit(status = 0)
}
suppressPackageStartupMessages(library(RSiena))

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "saom.yaml"))

ORACLE <- "rsiena"

write_values <- function(oracle, fixture, variant, values) {
  d <- file.path(golden_dir, oracle, fixture)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(c("value", formatC(as.numeric(values), digits = 17, format = "g")),
             file.path(d, paste0(variant, ".csv")))
}

# meta.json is flat and simple; parse it without a JSON dependency.
read_meta <- function(path) {
  txt <- paste(readLines(path), collapse = "")
  get <- function(key) sub(sprintf('.*"%s": *([^,}]+).*', key), "\\1", txt)
  list(n = as.integer(get("n")), directed = get("directed") == "true")
}

# ---- The panel, as validation/saom_fixtures.py builds it -------------------

toggle <- function(a, mult, add, mod) {
  n <- nrow(a)
  for (i in 0:(n - 1)) {
    for (j in 0:(n - 1)) {
      if (i != j && (mult * (i * n + j) + add) %% mod == 0) {
        a[i + 1, j + 1] <- 1 - a[i + 1, j + 1]
      }
    }
  }
  a
}

panel_waves <- function(n, edges) {
  w1 <- matrix(0, n, n)
  keep <- edges$from != edges$to
  w1[cbind(edges$from[keep] + 1, edges$to[keep] + 1)] <- 1
  w2 <- toggle(w1, 37, 11, 7)
  w3 <- toggle(w2, 53, 5, 9)
  list(w1, w2, w3)
}

covariates <- function(n) {
  i <- 0:(n - 1)
  changing <- sapply(0:2, function(w) (i * 3 + 2 * w) %% 4)
  dyadic <- outer(i, i, function(a, b) ((a * 5 + b * 3) %% 7) / 2)
  list(constant = (i * 7) %% 5, changing = changing, binary = i %% 2, dyadic = dyadic)
}

# ---- RSiena ------------------------------------------------------------------

make_data <- function(waves, cov) {
  n <- nrow(waves[[1]])
  net <- sienaDependent(array(unlist(waves), dim = c(n, n, length(waves))))
  constant <- coCovar(cov$constant)
  changing <- varCovar(cov$changing)
  binary <- coCovar(cov$binary)
  dyadic <- coDyadCovar(cov$dyadic)
  sienaDataCreate(net, constant, changing, binary, dyadic)
}

# Runs RSiena with one effect added to the default specification and
# returns the row of targets2 for that effect (or the rate rows), by name.
rsiena_targets <- function(dat, engine) {
  eff <- getEffects(dat)
  fn <- engine$fn
  name <- engine$effect
  cov <- if (is.null(engine$covariate)) "" else engine$covariate
  if (fn == "saom_targets") {
    # do.call() passes the short name as a string constant, which is what
    # `character = TRUE` expects.
    if (is.null(engine$parameter)) {
      eff <- do.call(includeEffects, list(eff, name, interaction1 = cov,
                                          character = TRUE, verbose = FALSE))
    } else {
      eff <- do.call(setEffect, list(eff, name, parameter = as.numeric(engine$parameter),
                                     interaction1 = cov, character = TRUE,
                                     verbose = FALSE))
    }
  }
  alg <- sienaAlgorithmCreate(nsub = 0, n3 = 5, cond = FALSE, simOnly = TRUE,
                              seed = 1, projname = NULL)
  ans <- siena07(alg, data = dat, effects = eff, batch = TRUE, silent = TRUE)
  re <- ans$requestedEffects
  t2 <- ans$targets2
  if (fn == "saom_distances") {
    rows <- which(re$type == "rate")
    periods <- ncol(t2)
    if (length(rows) != periods) stop("unexpected number of rate rows")
    return(vapply(seq_len(periods), function(m) t2[rows[m], m], numeric(1)))
  }
  rows <- which(re$shortName == name & re$type == "eval" & re$interaction1 == cov)
  if (length(rows) != 1L) {
    stop(sprintf("effect %s (%s) matched %d rows", name, cov, length(rows)))
  }
  as.numeric(t2[rows, ])
}

meta_of <- function(name) read_meta(file.path(fixtures_dir, name, "meta.json"))

fixtures <- sort(list.files(fixtures_dir))
fixtures <- fixtures[file.exists(file.path(fixtures_dir, fixtures, "meta.json"))]
failed <- character()
for (fixture in fixtures) {
  meta <- meta_of(fixture)
  if (!meta$directed) next
  edges <- utils::read.csv(file.path(fixtures_dir, fixture, "edges.csv"))
  waves <- panel_waves(meta$n, edges)
  dat <- tryCatch(make_data(waves, covariates(meta$n)), error = function(e) {
    cat(sprintf("%s: RSiena cannot build the panel: %s\n", fixture, conditionMessage(e)))
    NULL
  })
  if (is.null(dat)) next
  for (variant in names(manifest)) {
    spec <- manifest[[variant]]
    if (is.null(spec$oracles$rsiena)) next
    values <- tryCatch(rsiena_targets(dat, spec$engine), error = function(e) {
      cat(sprintf("%s / %s: %s\n", fixture, variant, conditionMessage(e)))
      NULL
    })
    if (is.null(values)) {
      failed <- c(failed, paste(fixture, variant))
      next
    }
    write_values(ORACLE, fixture, variant, values)
  }
  cat(sprintf("%s: done\n", fixture))
}
dir.create(file.path(golden_dir, ORACLE), recursive = TRUE, showWarnings = FALSE)
writeLines(sprintf("RSiena %s", as.character(utils::packageVersion("RSiena"))),
           file.path(golden_dir, ORACLE, "VERSION"))
if (length(failed)) {
  cat(sprintf("RSiena gave no answer for %d fixture/variant pairs (not run).\n",
              length(failed)))
}
