#!/usr/bin/env Rscript
# Validation at scale against igraph (Csardi and Nepusz 2006) and sna
# (Butts 2008).
#
#   Rscript validation/scale/check_scale.R
#
# Run from the repository root. Needs Python 3 with networkx and PyYAML
# (to generate the networks and to call the engine), and the R packages
# igraph, sna, and yaml. Set PYTHON to choose the Python interpreter.
#
# Steps:
# 1. `check_scale.py --prepare-r` generates the scale networks
#    (validation/scale/generate.py) and, for each one of up to 5,000 nodes,
#    writes a directory laid out like validation/: that network as the only
#    fixture, and the manifests restricted to the measures igraph and sna
#    run on it within their cost budgets (sna's matrix routines only on
#    the smaller networks).
# 2. Every R oracle script in validation/oracles/ runs unchanged in each
#    directory (they read fixtures and manifests relative to the working
#    directory) and writes its golden files there.
# 3. `check_scale.py --compare-r` compares the engine with those golden
#    files, using the manifest's tolerances and documented differences,
#    exactly as validation/compare.py does.
#
# Exits with status 1 if an oracle script fails or any comparison fails.

root <- normalizePath(".")
check <- file.path(root, "validation", "scale", "check_scale.py")
if (!file.exists(check)) {
  stop("Run this script from the repository root.", call. = FALSE)
}
for (pkg in c("igraph", "sna", "yaml")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " is required.", call. = FALSE)
}
python <- Sys.getenv("PYTHON")
if (!nzchar(python)) {
  python <- if (nzchar(Sys.which("python3"))) "python3" else "python"
}
rscript <- file.path(R.home("bin"), "Rscript")
work <- file.path(root, "validation", "scale", "work", "r")

started <- proc.time()[["elapsed"]]
if (system2(python, c(shQuote(check), "--prepare-r", shQuote(work))) != 0) {
  stop("Preparing the scale networks failed.", call. = FALSE)
}

scripts <- sort(Sys.glob(file.path(root, "validation", "oracles", "*.R")))
trees <- sort(list.dirs(work, recursive = FALSE))
trees <- trees[dir.exists(file.path(trees, "validation", "fixtures"))]
failed <- character()
for (tree in trees) {
  t0 <- proc.time()[["elapsed"]]
  for (script in scripts) {
    old <- setwd(tree)
    status <- system2(rscript, shQuote(script))
    setwd(old)
    if (status != 0) failed <- c(failed, sprintf("%s on %s", basename(script), basename(tree)))
  }
  cat(sprintf("== %s: R oracles took %.0f s\n", basename(tree), proc.time()[["elapsed"]] - t0))
}

status <- system2(python, c(shQuote(check), "--compare-r", shQuote(work)))
cat(sprintf("Scale check with igraph and sna took %.0f s.\n",
            proc.time()[["elapsed"]] - started))
if (length(failed)) {
  cat("Oracle scripts failed:", paste(failed, collapse = "; "), "\n")
}
quit(status = if (length(failed) || status != 0) 1L else 0L)
