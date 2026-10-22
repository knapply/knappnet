# Makes the repository README the home page of the website.
#
#   Rscript r/pkgdown/home.R        # from the repository root
#   Rscript pkgdown/home.R          # or from r/
#
# The README lives at the repository root, outside the R package, so it is
# copied to r/pkgdown/index.md, which pkgdown uses for the home page in
# preference to any other file. Links and images relative to the
# repository root would break on the website, so they are pointed at
# GitHub instead. The copy is ignored by git and by R CMD build; edit the
# README, never the copy.

repo <- "https://github.com/knapply/knappnet"
branch <- "main"

root <- if (file.exists(file.path("r", "DESCRIPTION"))) {
  "."
} else if (file.exists("DESCRIPTION") && file.exists(file.path("..", "README.md"))) {
  ".."
} else {
  stop("Run this script from the repository root or from r/.", call. = FALSE)
}

lines <- readLines(file.path(root, "README.md"), encoding = "UTF-8", warn = FALSE)

# A link target is relative when it has no scheme, is not an anchor on the
# same page, and is not site-absolute.
relative <- "(?!(?:[a-zA-Z][a-zA-Z0-9+.-]*:|#|/))"

# Images first, served raw so they display; then ordinary links, which
# open the file's page on GitHub.
lines <- gsub(paste0("(!\\[[^]]*\\]\\()", relative, "([^)\\s]+)"),
              sprintf("\\1https://raw.githubusercontent.com/knapply/knappnet/%s/\\2", branch),
              lines, perl = TRUE)
lines <- gsub(paste0("((?<!!)\\[[^]]*\\]\\()", relative, "([^)\\s]+)"),
              sprintf("\\1%s/blob/%s/\\2", repo, branch),
              lines, perl = TRUE)

out <- file.path(root, "r", "pkgdown", "index.md")
if (root == "..") out <- file.path("pkgdown", "index.md")
writeLines(lines, out, useBytes = TRUE)
message("Wrote ", out)
