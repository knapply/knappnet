#!/usr/bin/env Rscript
# Runs every R code block in README.md and checks that the output shown in
# the README (lines starting with "#>") is exactly what R prints. Blocks
# preceded by the line <!-- norun --> are skipped (for example, install
# instructions).
#
#   Rscript tools/check_readme.R            # check; exits 1 on any difference
#   Rscript tools/check_readme.R --update   # rewrite README.md with fresh output
#
# Requires the package (and igraph) to be installed. Run from the repository
# root.

update <- "--update" %in% commandArgs(trailingOnly = TRUE)
path <- "README.md"
lines <- readLines(path, encoding = "UTF-8", warn = FALSE)

options(width = 80, digits = 7, useFancyQuotes = FALSE)
env <- new.env(parent = globalenv())

render <- function(block) {
  code <- block[!startsWith(block, "#>")]
  exprs <- parse(text = code, keep.source = TRUE)
  refs <- attr(exprs, "srcref")
  out <- character()
  last <- 0L
  for (i in seq_along(exprs)) {
    end <- refs[[i]][3]
    out <- c(out, code[seq.int(last + 1L, length.out = end - last)])
    last <- end
    printed <- utils::capture.output({
      v <- withVisible(eval(exprs[[i]], env))
      if (v$visible) print(v$value)
    })
    out <- c(out, ifelse(nzchar(printed), paste("#>", printed), "#>"))
  }
  if (last < length(code)) out <- c(out, code[(last + 1L):length(code)])
  out
}

result <- character()
i <- 1L
problems <- 0L
while (i <= length(lines)) {
  line <- lines[[i]]
  runnable <- line == "```r" && !(i > 1L && lines[[i - 1L]] == "<!-- norun -->")
  if (!runnable) {
    result <- c(result, line)
    i <- i + 1L
    next
  }
  close <- i + which(lines[(i + 1L):length(lines)] == "```")[1]
  block <- lines[(i + 1L):(close - 1L)]
  fresh <- render(block)
  if (!identical(fresh, block)) {
    problems <- problems + 1L
    cat(sprintf("README.md block starting at line %d differs from R's output:\n", i))
    width <- max(length(block), length(fresh))
    for (k in seq_len(width)) {
      a <- if (k <= length(block)) block[[k]] else ""
      b <- if (k <= length(fresh)) fresh[[k]] else ""
      if (!identical(a, b)) cat(sprintf("  README: %s\n  R:      %s\n", a, b))
    }
  }
  result <- c(result, line, fresh, "```")
  i <- close + 1L
}

if (update) {
  writeLines(result, path, useBytes = TRUE)
  cat(sprintf("Updated %s (%d block(s) changed).\n", path, problems))
} else if (problems > 0L) {
  cat("Run `Rscript tools/check_readme.R --update` and review the changes.\n")
  quit(status = 1L)
} else {
  cat("README.md output matches R.\n")
}
