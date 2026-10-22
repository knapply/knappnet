#!/usr/bin/env Rscript
# Times this package against igraph on random graphs.
#
#   Rscript bench/bench_igraph.R
#
# Requires the package to be installed, plus igraph. Results depend on the
# machine, so report the core count alongside any numbers.

suppressPackageStartupMessages({
  library(knappnet)
  library(igraph)
})

random_edges <- function(n, avg_degree, directed, seed) {
  set.seed(seed)
  m <- as.integer(n * avg_degree / if (directed) 1 else 2)
  e <- data.frame(from = sample.int(n, m, TRUE), to = sample.int(n, m, TRUE))
  e <- e[e$from != e$to, ]
  e[!duplicated(e), ]
}

best_time <- function(expr, times = 3L) {
  f <- substitute(expr)
  env <- parent.frame()
  min(vapply(seq_len(times), function(i) system.time(eval(f, env))[["elapsed"]], 0))
}

cat(sprintf("cores: %d\n", parallel::detectCores()))
cat(sprintf("%-28s%-14s%12s%12s%10s\n", "graph", "measure", "ours (s)", "igraph (s)", "speedup"))
for (cfg in list(list(n = 5000, d = 10, dir = FALSE),
                 list(n = 5000, d = 10, dir = TRUE),
                 list(n = 20000, d = 10, dir = FALSE))) {
  e <- random_edges(cfg$n, cfg$d, cfg$dir, cfg$n)
  nodes <- data.frame(id = seq_len(cfg$n))
  g <- netframe(e, nodes, directed = cfg$dir)
  ig <- graph_from_data_frame(e, directed = cfg$dir, vertices = nodes)
  label <- sprintf("n=%d, m=%d, %s", cfg$n, nrow(e), if (cfg$dir) "dir" else "und")

  ours <- best_time(betweenness(g))
  theirs <- best_time(igraph::betweenness(ig), times = 1L)
  cat(sprintf("%-28s%-14s%12.3f%12.3f%9.1fx\n", label, "betweenness", ours, theirs, theirs / ours))

  if (cfg$dir) {
    ours <- best_time(triad_census(g))
    theirs <- best_time(igraph::triad_census(ig))
    cat(sprintf("%-28s%-14s%12.3f%12.3f%9.1fx\n", label, "triad census", ours, theirs,
                theirs / ours))
  }
}
