#!/usr/bin/env Rscript
# Golden values from sna (Butts 2008) for flow betweenness (Freeman,
# Borgatti, and White 1991), and from keyplayer (An and Liu 2016) for the
# fragmentation centrality of single nodes (Borgatti 2006)
# (validation/manifest.d/resilience.yaml).
#
#   Rscript validation/oracles/r_resilience.R
#
# Requires the sna and yaml packages; keyplayer is used when installed.
# Run from the repository root. Undefined values (NaN or NA) are written as
# -1, as the engine runner does.

suppressPackageStartupMessages(library(sna))
has_keyplayer <- requireNamespace("keyplayer", quietly = TRUE)

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "resilience.yaml"))

write_values <- function(oracle, fixture, variant, values) {
  values <- as.numeric(values)
  values[is.na(values)] <- -1
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

# Whether the manifest declares the oracle undefined on this fixture.
undefined <- function(ospec, name, meta) {
  tags <- c(meta$tags, if (meta$directed) "directed" else "undirected")
  name %in% unlist(ospec$undefined_on) || any(tags %in% unlist(ospec$undefined_when))
}

read_meta <- function(path) {
  txt <- paste(readLines(path), collapse = "")
  get <- function(key) sub(sprintf('.*"%s": *([^,}]+).*', key), "\\1", txt)
  tags <- sub('.*"tags": *\\[([^]]*)\\].*', "\\1", txt)
  tags <- gsub('[" ]', "", strsplit(tags, ",")[[1]])
  list(n = as.integer(get("n")),
       directed = get("directed") == "true",
       weighted = get("weighted") == "true",
       tags = tags[nzchar(tags)])
}

sna_values <- function(variant, m, mw, meta) {
  gmode <- if (meta$directed) "digraph" else "graph"
  switch(variant,
    flow_betweenness = sna::flowbet(m, gmode = gmode, cmode = "rawflow", ignore.eval = TRUE),
    flow_betweenness_normalized = sna::flowbet(m, gmode = gmode, cmode = "normflow",
                                               ignore.eval = TRUE),
    flow_betweenness_fractional = sna::flowbet(m, gmode = gmode, cmode = "fracflow",
                                               ignore.eval = TRUE),
    flow_betweenness_capacity = sna::flowbet(mw, gmode = gmode, cmode = "rawflow",
                                             ignore.eval = FALSE),
    NULL
  )
}

keyplayer_values <- function(variant, m) {
  switch(variant,
    kpp_distance_fragmentation_single = keyplayer::fragment(m, binary = TRUE),
    NULL
  )
}

counts <- c(sna = 0L, keyplayer = 0L)
for (name in sort(list.dirs(fixtures_dir, full.names = FALSE, recursive = FALSE))) {
  d <- file.path(fixtures_dir, name)
  meta <- read_meta(file.path(d, "meta.json"))
  e <- utils::read.csv(file.path(d, "edges.csv"))
  e <- e[e$from != e$to, , drop = FALSE]
  m <- matrix(0, meta$n, meta$n)
  mw <- m
  if (nrow(e)) {
    m[cbind(e$from + 1L, e$to + 1L)] <- 1
    if (!meta$directed) m[cbind(e$to + 1L, e$from + 1L)] <- 1
    if (meta$weighted) {
      mw[cbind(e$from + 1L, e$to + 1L)] <- e$weight
      if (!meta$directed) mw[cbind(e$to + 1L, e$from + 1L)] <- e$weight
    }
  }

  for (variant in names(manifest)) {
    spec <- manifest[[variant]]
    if (!applies(spec, meta)) next
    ospec <- spec$oracles$sna
    if (!is.null(ospec) && !undefined(ospec, name, meta)) {
      v <- sna_values(variant, m, mw, meta)
      if (!is.null(v)) {
        write_values("sna", name, variant, v)
        counts[["sna"]] <- counts[["sna"]] + 1L
      }
    }
    ospec <- spec$oracles$keyplayer
    if (has_keyplayer && !is.null(ospec) && !undefined(ospec, name, meta)) {
      v <- keyplayer_values(variant, m)
      if (!is.null(v)) {
        write_values("keyplayer", name, variant, v)
        counts[["keyplayer"]] <- counts[["keyplayer"]] + 1L
      }
    }
  }
}

writeLines(sprintf("sna %s", utils::packageVersion("sna")),
           file.path(golden_dir, "sna", "VERSION"))
if (has_keyplayer) {
  writeLines(sprintf("keyplayer %s", utils::packageVersion("keyplayer")),
             file.path(golden_dir, "keyplayer", "VERSION"))
} else {
  cat("keyplayer is not installed; its golden files were not written.\n")
}
cat(sprintf("Wrote %d sna and %d keyplayer golden files (resilience).\n",
            counts[["sna"]], counts[["keyplayer"]]))
