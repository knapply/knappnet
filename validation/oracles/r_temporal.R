#!/usr/bin/env Rscript
# Golden values from tsna (Bender-deMoll and Morris) for earliest-arrival
# paths and temporal reachability over instantaneous contacts
# (validation/manifest.d/temporal.yaml, variants with a `tsna` oracle).
#
#   Rscript validation/oracles/r_temporal.R
#
# Requires the tsna, networkDynamic, network, and yaml packages; exits
# quietly when tsna is not installed. Run from the repository root.
#
# Every fixture is read as a contact sequence with the times of
# validation/derived/temporal.py, with every contact instantaneous: contact
# e (0-based) starts at (37 e + 11) mod p, p = max(3, floor(m / 3) + 2). An
# instantaneous contact is a networkDynamic point spell (onset = terminus),
# and tsna::tPath() with graph.step.time = 0 follows paths whose contact
# times are non-decreasing, the engine's default rule, from `start` (the
# earliest contact) until `end` (one past the latest contact, since tPath
# ignores contacts starting at or after `end`). Model times are
# `start + tdist`; unreachable nodes have tdist = Inf. Interval contacts are
# left out because tsna lets a path join a contact part-way through its
# spell, which the engine does not model (see ?temporal_distances).
#
# Infinities are written as -1000 and NaN as -999, as the engine runner
# does.

if (!requireNamespace("tsna", quietly = TRUE) ||
    !requireNamespace("networkDynamic", quietly = TRUE)) {
  cat("tsna is not installed; no tsna golden files written (temporal).\n")
  quit(status = 0)
}
suppressPackageStartupMessages({
  library(network)
  library(networkDynamic)
  library(tsna)
})

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "temporal.yaml"))

UNREACHABLE <- -1000
UNDEFINED <- -999

write_values <- function(oracle, fixture, variant, values) {
  values <- as.numeric(values)
  values[is.nan(values) | is.na(values)] <- UNDEFINED
  values[is.infinite(values)] <- UNREACHABLE
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

# Start times of validation/derived/temporal.py for m contacts.
contact_starts <- function(m) {
  if (m == 0L) return(numeric())
  p <- max(3, m %/% 3 + 2)
  e <- seq_len(m) - 1
  (e * 37 + 11) %% p
}

# Earliest arrival at every node from every node (row = origin), and the
# number of other nodes each node reaches, by tsna::tPath().
tsna_paths <- function(e, meta) {
  n <- meta$n
  m <- nrow(e)
  arrival <- matrix(Inf, n, n)
  if (m == 0L) {
    diag(arrival) <- 0
    return(list(arrival = arrival, reach = rep(0, n)))
  }
  t0 <- contact_starts(m)
  from <- e$from + 1L
  to <- e$to + 1L
  if (!meta$directed) {
    # One undirected edge per unordered pair, whatever the contact's order.
    tail <- pmin(from, to)
    head <- pmax(from, to)
    from <- tail
    to <- head
  }
  base <- network::network.initialize(n, directed = meta$directed, loops = TRUE,
                                      multiple = FALSE)
  spells <- data.frame(onset = t0, terminus = t0, tail = from, head = to)
  nd <- suppressMessages(networkDynamic::networkDynamic(base.net = base, edge.spells = spells,
                                                        verbose = FALSE))
  start <- min(t0)
  end <- max(t0) + 1
  reach <- numeric(n)
  for (v in seq_len(n)) {
    p <- suppressWarnings(tsna::tPath(nd, v = v, direction = "fwd", type = "earliest.arrive",
                                      start = start, end = end, graph.step.time = 0))
    d <- as.numeric(p$tdist)
    arrival[v, ] <- ifelse(is.finite(d), start + d, Inf)
    reach[v] <- sum(is.finite(d)) - 1
  }
  list(arrival = arrival, reach = reach)
}

count <- 0L
for (name in sort(list.dirs(fixtures_dir, full.names = FALSE, recursive = FALSE))) {
  d <- file.path(fixtures_dir, name)
  meta <- read_meta(file.path(d, "meta.json"))
  e <- utils::read.csv(file.path(d, "edges.csv"))
  wanted <- Filter(function(v) "tsna" %in% names(manifest[[v]]$oracles) &&
                     applies(manifest[[v]], meta), names(manifest))
  if (!length(wanted) || meta$n == 0L) next
  res <- tsna_paths(e, meta)
  for (variant in wanted) {
    fn <- manifest[[variant]]$engine$fn
    values <- switch(fn,
      temporal_paths = as.vector(t(res$arrival)),
      temporal_reach = res$reach,
      NULL
    )
    if (is.null(values)) next
    write_values("tsna", name, variant, values)
    count <- count + 1L
  }
}

dir.create(file.path(golden_dir, "tsna"), recursive = TRUE, showWarnings = FALSE)
writeLines(sprintf("tsna %s; networkDynamic %s", utils::packageVersion("tsna"),
                   utils::packageVersion("networkDynamic")),
           file.path(golden_dir, "tsna", "VERSION"))
cat(sprintf("Wrote %d tsna golden files (temporal).\n", count))
