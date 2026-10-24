#!/usr/bin/env Rscript
# Golden values for community detection and partition comparison from
# igraph (Csardi and Nepusz 2006).
#
#   Rscript validation/oracles/r_communities.R
#
# Requires the igraph and yaml packages (cluster_optimal() needs igraph
# built with GLPK). Run from the repository root. Undefined modularity
# (no edges) is written as -1, as the engine runners do. A call that
# igraph rejects, or an argument this igraph version lacks, writes nothing
# for that fixture.

suppressPackageStartupMessages(library(igraph))

root <- normalizePath(".")
fixtures_dir <- file.path(root, "validation", "fixtures")
golden_dir <- file.path(root, "validation", "golden")
manifest <- yaml::read_yaml(file.path(root, "validation", "manifest.d", "communities.yaml"))

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

# The two fixed partitions compared by the comparison measures (1-based
# labels): floor(3 v / n) and (7 v) mod 4 for v = 0, ..., n - 1.
blocks <- function(n) (3L * (seq_len(n) - 1L)) %/% n + 1L
stripes <- function(n) (7L * (seq_len(n) - 1L)) %% 4L + 1L

# NaN (modularity without edges) is written as -1, as the engine does.
q0 <- function(x) ifelse(is.nan(x) | is.na(x), -1, x)

# Modularity of a membership on graph h with weights w (NULL = unweighted
# graph h without a weight attribute).
modq <- function(h, memb, w = NULL, resolution = 1) {
  q0(igraph::modularity(h, memb, weights = w, resolution = resolution))
}

# Mean over set.seed(1..seeds) of f().
seeded_mean <- function(seeds, f) {
  mean(vapply(seq_len(seeds), function(s) {
    set.seed(s)
    f()
  }, numeric(1)))
}

# CPM quality (1/2m) sum_ij (A_ij - gamma) delta(c_i, c_j) on the
# unweighted graph h (a loop adds 2 to A_ii).
cpm_quality <- function(h, memb, gamma) {
  m <- igraph::ecount(h)
  if (m == 0) return(-1)
  ends <- igraph::ends(h, igraph::E(h), names = FALSE)
  inside <- sum(memb[ends[, 1]] == memb[ends[, 2]])
  sizes <- tabulate(memb)
  (2 * inside - gamma * sum(sizes^2)) / (2 * m)
}

# cluster_leiden()'s resolution argument was renamed in igraph 2.0.
leiden <- function(h, objective, resolution) {
  args <- names(formals(igraph::cluster_leiden))
  call <- list(h, objective_function = objective, n_iterations = -1, weights = NA)
  if ("resolution" %in% args) {
    call$resolution <- resolution
  } else {
    call$resolution_parameter <- resolution
  }
  do.call(igraph::cluster_leiden, call)
}

label_prop <- function(h, variant) {
  if (variant == "dominance") {
    return(igraph::cluster_label_prop(h, weights = NA))
  }
  if (!("lpa_variant" %in% names(formals(igraph::cluster_label_prop)))) {
    stop("this igraph version has no lpa_variant argument")
  }
  igraph::cluster_label_prop(h, weights = NA, lpa_variant = variant)
}

igraph_values <- function(spec, g, gu, meta) {
  e <- spec$engine
  weighted <- isTRUE(e$weighted)
  h <- if (weighted) g else gu
  w <- if (weighted) igraph::E(g)$weight else NULL
  n <- meta$n
  out <- if (is.null(e$output)) "" else e$output
  switch(e$fn,
    walktrap = {
      steps <- if (is.null(e$steps)) 4 else e$steps
      cl <- igraph::cluster_walktrap(h, weights = w, steps = steps)
      if (out == "membership") {
        igraph::membership(cl)
      } else if (out == "merges") {
        mg <- cl$merges - 1
        lapply(seq_len(nrow(mg)), function(i) sort(mg[i, ]))
      } else if (out == "cut") {
        igraph::cut_at(cl, no = e$groups)
      } else {
        modq(h, igraph::membership(cl), w)
      }
    },
    girvan_newman_cut = {
      cl <- igraph::cluster_edge_betweenness(igraph::simplify(gu), weights = NULL)
      igraph::cut_at(cl, no = e$groups)
    },
    label_propagation = {
      if (out != "quality") return(NULL)
      variant <- if (is.null(e$variant)) "dominance" else e$variant
      seeded_mean(e$seeds, function() modq(gu, igraph::membership(label_prop(gu, variant))))
    },
    infomap = {
      if (out != "quality") return(NULL)
      seeded_mean(e$seeds, function() {
        q0(igraph::cluster_infomap(gu, nb.trials = e$trials)$codelength)
      })
    },
    leiden = {
      if (out != "mean_quality") return(NULL)
      res <- if (is.null(e$resolution)) 1 else e$resolution
      if (identical(e$objective, "cpm")) {
        seeded_mean(e$seeds, function() {
          cpm_quality(gu, igraph::membership(leiden(gu, "CPM", res)), res)
        })
      } else {
        seeded_mean(e$seeds, function() {
          modq(gu, igraph::membership(leiden(gu, "modularity", res)), NULL, res)
        })
      }
    },
    spinglass = {
      if (out != "quality") return(NULL)
      seeded_mean(e$seeds, function() {
        modq(gu, igraph::membership(igraph::cluster_spinglass(gu, weights = NA)))
      })
    },
    optimal_modularity = {
      cl <- igraph::cluster_optimal(h, weights = w)
      modq(h, igraph::membership(cl), w)
    },
    fluid_communities = {
      if (out != "quality") return(NULL)
      s <- igraph::simplify(gu)
      seeded_mean(e$seeds, function() {
        modq(gu, igraph::membership(igraph::cluster_fluid_communities(s, e$k)))
      })
    },
    voronoi = {
      if (out != "modularity") return(NULL)
      ns <- asNamespace("igraph")
      if (exists("cluster_voronoi", envir = ns, inherits = FALSE)) {
        memb <- igraph::membership(get("cluster_voronoi", envir = ns)(gu))
      } else if (exists("community_voronoi_impl", envir = ns, inherits = FALSE)) {
        # igraph 2.3 has the C library's Voronoi communities only as an
        # internal wrapper; a radius of -1 selects it automatically.
        res <- get("community_voronoi_impl", envir = ns)(gu, mode = "out", radius = -1)
        memb <- as.integer(res$membership)
        if (length(memb) && min(memb) == 0L) memb <- memb + 1L
      } else {
        stop("this igraph version has no Voronoi communities")
      }
      modq(gu, memb)
    },
    compare_partitions = {
      a <- blocks(n)
      b <- stripes(n)
      if (e$method == "split_join") {
        sum(igraph::split_join_distance(a, b))
      } else {
        method <- c(vi = "vi", nmi = "nmi", rand = "rand",
                    adjusted_rand = "adjusted.rand")[[e$method]]
        igraph::compare(a, b, method = method)
      }
    },
    NULL
  )
}

count <- 0L
for (name in sort(list.dirs(fixtures_dir, full.names = FALSE, recursive = FALSE))) {
  d <- file.path(fixtures_dir, name)
  if (!file.exists(file.path(d, "meta.json"))) next
  meta <- read_meta(file.path(d, "meta.json"))
  e <- utils::read.csv(file.path(d, "edges.csv"))
  gu <- igraph::make_empty_graph(n = meta$n, directed = meta$directed)
  if (nrow(e)) {
    gu <- igraph::add_edges(gu, as.vector(t(as.matrix(e[, c("from", "to")]))) + 1L)
  }
  g <- gu
  if (meta$weighted && nrow(e)) igraph::E(g)$weight <- e$weight
  for (variant in names(manifest)) {
    spec <- manifest[[variant]]
    if (!("igraph" %in% names(spec$oracles)) || !applies(spec, meta)) next
    if (name %in% unlist(spec$oracles$igraph$undefined_on)) next
    v <- tryCatch(igraph_values(spec, g, gu, meta), error = function(err) NULL)
    if (!is.null(v)) {
      write_values(name, variant, v)
      count <- count + 1L
    }
  }
}

dir.create(file.path(golden_dir, "igraph"), recursive = TRUE, showWarnings = FALSE)
writeLines(sprintf("igraph %s", utils::packageVersion("igraph")),
           file.path(golden_dir, "igraph", "VERSION"))
cat(sprintf("Wrote %d igraph golden files for communities.\n", count))
