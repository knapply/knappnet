#!/usr/bin/env Rscript
# Benchmarks this package against igraph and sna on random networks, and
# checks that every pair of tools being timed gives the same answer: a
# fast wrong answer is not a result.
#
#   Rscript bench/bench.R            # full run (several minutes)
#   Rscript bench/bench.R --quick    # small networks, for CI
#
# Writes bench/results/r-<platform>.csv and prints a Markdown table. Exits
# with status 1 if any pair of tools disagrees.
#
# Timings are medians of repeated runs (operations under 0.2 s are looped
# so each timed run lasts at least 0.2 s). Graphs are built before timing:
# this package caches its compressed graph inside the net_frame, and igraph
# and sna objects already hold theirs, so the times compare the
# computations alone. Building a network from a data frame is timed
# separately. igraph and sna run on one core; this package is timed both
# on one thread and on every core.

# igraph and sna are used through `::` only: attaching them would mask
# this package's functions of the same names.
suppressPackageStartupMessages(library(knappnet))
for (pkg in c("igraph", "sna")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " is required.", call. = FALSE)
}

quick <- "--quick" %in% commandArgs(trailingOnly = TRUE)
budget <- if (quick) 2 else 15        # seconds per timed expression
scale <- if (quick) 0.1 else 1

now <- function() proc.time()[["elapsed"]]

# Median seconds per call of f().
timing <- function(f) {
  t0 <- now()
  f()
  first <- now() - t0
  if (first >= 0.2) {
    reps <- max(0L, min(4L, floor((budget - first) / first)))
    times <- c(first, vapply(seq_len(reps), function(i) {
      t0 <- now()
      f()
      now() - t0
    }, 0))
    return(stats::median(times))
  }
  inner <- max(1L, ceiling(0.2 / max(first, 1e-4)))
  stats::median(vapply(1:5, function(i) {
    t0 <- now()
    for (k in seq_len(inner)) f()
    (now() - t0) / inner
  }, 0))
}

same <- function(a, b, tol = 1e-8) {
  isTRUE(all.equal(unname(as.numeric(a)), unname(as.numeric(b)), tolerance = tol))
}

# Two partitions are the same if they group the nodes identically.
same_partition <- function(a, b) {
  length(unique(a)) == length(unique(b)) &&
    nrow(unique(data.frame(a, b))) == length(unique(a))
}

rows <- list()
notes <- character()
# `note` (optional) is printed after the row and, if the tools disagree,
# with the failure: for example both modularity values for Louvain.
record <- function(graph, measure, ours1, ours, tool, theirs, agree, note = NULL) {
  rows[[length(rows) + 1L]] <<- data.frame(
    graph = graph, measure = measure, knappnet_1_thread = ours1,
    knappnet = ours, tool = tool, tool_seconds = theirs,
    speedup_1_thread = theirs / ours1, speedup = theirs / ours,
    agree = agree, stringsAsFactors = FALSE)
  cat(sprintf("%-22s %-24s %-7s %9.4f %9.4f %9.4f %7.1fx %7.1fx %s%s\n", graph, measure,
              tool, ours1, ours, theirs, theirs / ours1, theirs / ours,
              if (agree) "agree" else "DISAGREE", if (is.null(note)) "" else paste0("  ", note)))
  if (!is.null(note)) notes[[paste(graph, measure, tool)]] <<- note
}

# Times a function that takes a `threads` argument on one thread and on
# all cores; functions without one are timed once.
both <- function(f, threaded = TRUE) {
  one <- timing(function() f(1L))
  all <- if (threaded) timing(function() f(0L)) else one
  c(one, all)
}

cpu <- switch(Sys.info()[["sysname"]],
  Linux = sub(".*: ", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]),
  Darwin = system("sysctl -n machdep.cpu.brand_string", intern = TRUE),
  Windows = Sys.getenv("PROCESSOR_IDENTIFIER"),
  "unknown")
cat(sprintf("%s; %s; %d cores; igraph %s; sna %s; knappnet %s\n", R.version.string, cpu,
            parallel::detectCores(), utils::packageVersion("igraph"), utils::packageVersion("sna"),
            utils::packageVersion("knappnet")))

n_big <- round(100000 * scale)
n_mid <- round(10000 * scale)
n_small <- round(2000 * scale)

set.seed(1)
und <- igraph::sample_gnm(n_mid, 5L * n_mid)
dir <- igraph::sample_gnm(n_mid, 5L * n_mid, directed = TRUE)
pa <- igraph::sample_pa(n_big, m = 3L, directed = FALSE)
small <- igraph::sample_gnm(n_small, 5L * n_small, directed = TRUE)

und_nf <- as_net_frame(und)
dir_nf <- as_net_frame(dir)
pa_nf <- as_net_frame(pa)
small_nf <- as_net_frame(small)
small_m <- as.matrix(small_nf)
invisible(lapply(list(und_nf, dir_nf, pa_nf, small_nf), degree))  # build caches

lab_und <- sprintf("G(n=%d, m=%d)", igraph::vcount(und), igraph::ecount(und))
lab_dir <- sprintf("directed G(%d, %d)", igraph::vcount(dir), igraph::ecount(dir))
lab_pa <- sprintf("scale-free n=%d, m=%d", igraph::vcount(pa), igraph::ecount(pa))
lab_small <- sprintf("directed G(%d, %d)", igraph::vcount(small), igraph::ecount(small))

cat(sprintf("%-22s %-24s %-7s %9s %9s %9s %8s %8s\n", "graph", "measure", "tool",
            "1 thread", "all", "tool", "1 thread", "all"))

# ---- Building a network ------------------------------------------------
e <- igraph::as_edgelist(pa, names = FALSE)
edges <- data.frame(source = as.integer(e[, 1]), target = as.integer(e[, 2]))
nodes <- data.frame(name = seq_len(igraph::vcount(pa)))
ours <- timing(function() net_degree(net_frame(edges, nodes, directed = FALSE)))
theirs <- timing(function() igraph::degree(igraph::graph_from_data_frame(edges, FALSE, nodes)))
record(lab_pa, "build + degree", ours, ours, "igraph", theirs,
       same(net_degree(net_frame(edges, nodes, directed = FALSE)),
            igraph::degree(igraph::graph_from_data_frame(edges, FALSE, nodes))))

# ---- Centrality --------------------------------------------------------
for (case in list(list(lab_und, und_nf, und), list(lab_dir, dir_nf, dir))) {
  t <- both(function(k) net_betweenness(case[[2]], threads = k))
  theirs <- timing(function() igraph::betweenness(case[[3]]))
  record(case[[1]], "betweenness", t[1], t[2], "igraph", theirs,
         same(net_betweenness(case[[2]]), igraph::betweenness(case[[3]])))
}

t <- both(function(k) net_closeness(und_nf, threads = k))
theirs <- timing(function() igraph::closeness(und))
record(lab_und, "closeness", t[1], t[2], "igraph", theirs,
       same(net_closeness(und_nf), igraph::closeness(und)))

t <- both(function(k) net_harmonic(und_nf, threads = k))
theirs <- timing(function() igraph::harmonic_centrality(und))
record(lab_und, "harmonic centrality", t[1], t[2], "igraph", theirs,
       same(net_harmonic(und_nf), igraph::harmonic_centrality(und)))

# ---- Directed and weighted shortest paths -------------------------------
# Directed closeness and harmonic centrality (out-distances, as igraph's
# default), and the two weighted algorithms: whole-number weights up to
# 1024 take Dial's bucket queue, other weights Dijkstra's algorithm.
t <- both(function(k) net_closeness(dir_nf, mode = "out", normalized = FALSE, threads = k))
theirs <- timing(function() igraph::closeness(dir, mode = "out", normalized = FALSE))
record(lab_dir, "closeness", t[1], t[2], "igraph", theirs,
       same(net_closeness(dir_nf, mode = "out", normalized = FALSE),
            igraph::closeness(dir, mode = "out", normalized = FALSE)))

t <- both(function(k) net_harmonic(dir_nf, mode = "out", threads = k))
theirs <- timing(function() igraph::harmonic_centrality(dir, mode = "out"))
record(lab_dir, "harmonic centrality", t[1], t[2], "igraph", theirs,
       same(net_harmonic(dir_nf, mode = "out"),
            igraph::harmonic_centrality(dir, mode = "out")))

n_w <- round(5000 * scale)
weighted_case <- function(directed, integer) {
  g <- igraph::sample_gnm(n_w, 5L * n_w, directed = directed)
  m <- igraph::ecount(g)
  w <- if (integer) sample.int(5L, m, replace = TRUE) else round(stats::runif(m, 0.1, 3), 6)
  igraph::E(g)$weight <- w
  e <- igraph::as_edgelist(g, names = FALSE)
  nf <- net_frame(data.frame(source = e[, 1], target = e[, 2], weight = w),
                 data.frame(name = seq_len(n_w)), directed = directed, weight = "weight", weight_type = "distance")
  invisible(net_degree(nf))  # build the cache
  label <- sprintf("%s%s weights, %d / %d", if (directed) "directed, " else "",
                   if (integer) "integer (Dial)" else "real (Dijkstra)", n_w, m)
  list(g = g, nf = nf, label = label)
}
for (case in list(weighted_case(FALSE, TRUE), weighted_case(FALSE, FALSE),
                  weighted_case(TRUE, TRUE), weighted_case(TRUE, FALSE))) {
  g <- case$g
  nf <- case$nf
  w <- igraph::E(g)$weight
  t <- both(function(k) net_betweenness(nf, weighted = TRUE, threads = k))
  theirs <- timing(function() igraph::betweenness(g, weights = w))
  record(case$label, "betweenness", t[1], t[2], "igraph", theirs,
         same(net_betweenness(nf, weighted = TRUE), igraph::betweenness(g, weights = w)))
  t <- both(function(k) net_closeness(nf, mode = "out", weighted = TRUE, normalized = FALSE,
                                  threads = k))
  theirs <- timing(function() igraph::closeness(g, mode = "out", weights = w))
  record(case$label, "closeness", t[1], t[2], "igraph", theirs,
         same(net_closeness(nf, mode = "out", weighted = TRUE, normalized = FALSE),
              igraph::closeness(g, mode = "out", weights = w)))
  t <- both(function(k) net_harmonic(nf, mode = "out", weighted = TRUE, threads = k))
  theirs <- timing(function() igraph::harmonic_centrality(g, mode = "out", weights = w))
  record(case$label, "harmonic centrality", t[1], t[2], "igraph", theirs,
         same(net_harmonic(nf, mode = "out", weighted = TRUE),
              igraph::harmonic_centrality(g, mode = "out", weights = w)))
}

t <- both(function(k) net_eigen(pa_nf, threads = k))
theirs <- timing(function() igraph::eigen_centrality(pa)$vector)
record(lab_pa, "eigenvector centrality", t[1], t[2], "igraph", theirs,
       same(net_eigen(pa_nf), igraph::eigen_centrality(pa)$vector, 1e-6))

t <- both(function(k) net_constraint(und_nf, threads = k))
theirs <- timing(function() igraph::constraint(und))
record(lab_und, "constraint", t[1], t[2], "igraph", theirs,
       same(net_constraint(und_nf), igraph::constraint(und)))

# ---- Cohesion and subgroups --------------------------------------------
t <- both(function(k) net_coreness(pa_nf), threaded = FALSE)
theirs <- timing(function() igraph::coreness(pa))
record(lab_pa, "k-cores", t[1], t[2], "igraph", theirs,
       same(net_coreness(pa_nf), igraph::coreness(pa)))

t <- both(function(k) net_components(pa_nf)$membership, threaded = FALSE)
theirs <- timing(function() igraph::components(pa)$membership)
record(lab_pa, "components", t[1], t[2], "igraph", theirs,
       same_partition(net_components(pa_nf)$membership, igraph::components(pa)$membership))

t <- both(function(k) net_max_cliques(pa_nf, threads = k))
theirs <- timing(function() igraph::max_cliques(pa))
record(lab_pa, "maximal cliques", t[1], t[2], "igraph", theirs,
       length(net_max_cliques(pa_nf)) == length(igraph::max_cliques(pa)))

t <- both(function(k) net_louvain(pa_nf, seed = 1), threaded = FALSE)
theirs <- timing(function() igraph::cluster_louvain(pa))
q_ours <- attr(net_louvain(pa_nf, seed = 1), "modularity")
q_theirs <- igraph::modularity(igraph::cluster_louvain(pa))
# Louvain is a randomized heuristic, so the tools find different
# partitions. Agreement means comparable quality: this package's
# modularity is at least 98% of igraph's (finding a higher modularity is
# fine).
record(lab_pa, "Louvain communities", t[1], t[2], "igraph", theirs,
       q_ours >= 0.98 * q_theirs,
       sprintf("modularity %.4f (knappnet) vs %.4f (igraph)", q_ours, q_theirs))

t <- both(function(k) net_triad_census(dir_nf, threads = k))
theirs <- timing(function() igraph::triad_census(dir))
record(lab_dir, "triad census", t[1], t[2], "igraph", theirs,
       same(net_triad_census(dir_nf), igraph::triad_census(dir)))

# ---- Against sna (adjacency-matrix input) -------------------------------
t <- both(function(k) net_betweenness(small_nf, threads = k))
theirs <- timing(function() sna::betweenness(small_m))
record(lab_small, "betweenness", t[1], t[2], "sna", theirs,
       same(net_betweenness(small_nf), sna::betweenness(small_m)))

t <- both(function(k) net_triad_census(small_nf, threads = k))
theirs <- timing(function() sna::triad.census(small_m))
record(lab_small, "triad census", t[1], t[2], "sna", theirs,
       same(net_triad_census(small_nf), sna::triad.census(small_m)))

t <- both(function(k) net_struct_equiv(small_nf, threads = k))
theirs <- timing(function() sna::sedist(small_m, method = "euclidean"))
record(lab_small, "structural equivalence", t[1], t[2], "sna", theirs,
       same(net_struct_equiv(small_nf), sna::sedist(small_m, method = "euclidean")))

groups <- rep_len(1:4, igraph::vcount(small))
t <- both(function(k) net_brokerage(small_nf, groups, threads = k))
theirs <- timing(function() sna::brokerage(small_m, groups))
br <- sna::brokerage(small_m, groups)$raw.nli
record(lab_small, "brokerage", t[1], t[2], "sna", theirs,
       same(as.matrix(net_brokerage(small_nf, groups)[, -1]),
            br[, c("w_I", "w_O", "b_IO", "b_OI", "b_O", "t")]))

# ---- QAP ----------------------------------------------------------------
n_qap <- if (quick) 30 else 100
reps <- if (quick) 200 else 1000
y <- matrix(rbinom(n_qap^2, 1, 0.1), n_qap)
x1 <- matrix(rbinom(n_qap^2, 1, 0.1), n_qap)
x2 <- matrix(rnorm(n_qap^2), n_qap)
lab_qap <- sprintf("%d nodes, %d draws", n_qap, reps)

t <- both(function(k) net_qap_cor(y, x1, reps = reps, seed = 1, threads = k))
theirs <- timing(function() sna::qaptest(list(y, x1), sna::gcor, g1 = 1, g2 = 2, reps = reps))
record(lab_qap, "QAP correlation", t[1], t[2], "sna", theirs,
       same(net_qap_cor(y, x1, reps = 0)$statistic, sna::gcor(y, x1)))

t <- both(function(k) net_mrqap(y, list(x1, x2), reps = reps, seed = 1, threads = k))
theirs <- timing(function() sna::netlm(y, list(x1, x2), nullhyp = "qapspp", reps = reps))
record(lab_qap, "MRQAP (DSP)", t[1], t[2], "sna", theirs,
       same(net_mrqap(y, list(x1, x2), reps = 0)$coefficients,
            sna::netlm(y, list(x1, x2), nullhyp = "classical")$coefficients))

# ---- Files --------------------------------------------------------------
path <- tempfile(fileext = ".graphml")
net_write(und_nf, path)
ours <- timing(function() net_read(path))
theirs <- timing(function() igraph::read_graph(path, format = "graphml"))
record(lab_und, "read GraphML", ours, ours, "igraph", theirs,
       nrow(net_read(path)$edges) == igraph::ecount(igraph::read_graph(path, format = "graphml")))

# ---- Output -------------------------------------------------------------
res <- do.call(rbind, rows)
dir.create("bench/results", showWarnings = FALSE, recursive = TRUE)
platform <- tolower(Sys.info()[["sysname"]])
out <- file.path("bench/results", sprintf("r-%s%s.csv", platform, if (quick) "-quick" else ""))
utils::write.csv(res, out, row.names = FALSE)

fmt_s <- function(s) ifelse(s >= 1, sprintf("%.2f s", s),
                     ifelse(s >= 1e-3, sprintf("%.1f ms", s * 1e3),
                            sprintf("%.0f µs", s * 1e6)))
fmt_x <- function(x) ifelse(x >= 10, sprintf("%.0fx", x), sprintf("%.1fx", x))
cat("\n| Network | Measure | Compared with | knappnet, 1 thread | knappnet, all cores |",
    "Other tool | Speedup, 1 thread | Speedup, all cores |\n", sep = "")
cat("|---|---|---|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  cat(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |\n", r$graph, r$measure, r$tool,
              fmt_s(r$knappnet_1_thread), fmt_s(r$knappnet), fmt_s(r$tool_seconds),
              fmt_x(r$speedup_1_thread), fmt_x(r$speedup)))
}

bad <- res[!res$agree, c("graph", "measure", "tool")]
if (nrow(bad)) {
  cat("\nResults differ for:\n")
  print(bad, row.names = FALSE)
  for (k in paste(bad$graph, bad$measure, bad$tool)) {
    if (!is.null(notes[k]) && !is.na(notes[k])) cat(k, ": ", notes[[k]], "\n", sep = "")
  }
  quit(status = 1L)
}
cat("\nAll results agree.\n")
