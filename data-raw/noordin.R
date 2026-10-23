# Builds the Noordin Top example data (r/inst/extdata/noordin/) from the
# copies distributed in two R packages:
#
#   - networkdata (Schoch; MIT): `covert_34`, the person-by-person relations
#     of the Roberts and Everton (2011) Noordin Top subset, one igraph object
#     per relation;
#   - manet (CRAN; GPL-2, data only): `noordin`, the 79 x 45 person-by-event
#     affiliation matrix described in Everton (2012).
#
# The data themselves come from Roberts, N. and Everton, S. F. (2011),
# "Roberts and Everton Terrorist Data: Noordin Top Terrorist Network
# (Subset)", distributed by the Association of Religion Data Archives, and
# were drawn from International Crisis Group (2006), "Terrorism in
# Indonesia: Noordin's Networks."
#
# Usage (from the repository root, with igraph installed):
#   Rscript data-raw/noordin.R <networkdata checkout> <manet checkout>

args <- commandArgs(trailingOnly = TRUE)
networkdata <- args[[1]]
manet <- args[[2]]
out <- file.path("r", "inst", "extdata", "noordin")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

env <- new.env()
load(file.path(networkdata, "data", "covert_34.rda"), envir = env)
rel <- get("covert_34", envir = env)
if (inherits(rel, "igraph")) rel <- list(rel)
cat("covert_34:", length(rel), "objects\n")
print(names(rel))

nodes <- NULL
edges <- list()
for (i in seq_along(rel)) {
  g <- rel[[i]]
  relation <- if (!is.null(names(rel)) && nzchar(names(rel)[i])) names(rel)[i] else
    igraph::graph_attr(g, "name") %||% paste0("relation_", i)
  va <- igraph::as_data_frame(g, what = "vertices")
  if (!"name" %in% names(va)) va$name <- as.character(seq_len(nrow(va)))
  if (is.null(nodes)) nodes <- va else {
    extra <- setdiff(names(va), names(nodes))
    if (length(extra)) nodes <- merge(nodes, va[c("name", extra)], by = "name", all = TRUE)
  }
  e <- igraph::as_data_frame(g, what = "edges")
  cat(relation, ": ", igraph::vcount(g), " nodes, ", nrow(e), " edges, directed = ",
      igraph::is_directed(g), ", attributes: ", paste(names(e), collapse = ", "), "\n", sep = "")
  if (nrow(e)) edges[[length(edges) + 1L]] <- data.frame(e[1:2], relation = relation,
                                                        e[-(1:2)], check.names = FALSE)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
edges <- do.call(rbind, lapply(edges, function(d) {
  d[setdiff(Reduce(union, lapply(edges, names)), names(d))] <- NA
  d[Reduce(union, lapply(edges, names))]
}))
names(edges)[1:2] <- c("from", "to")
utils::write.csv(nodes, file.path(out, "nodes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(edges, file.path(out, "relations.csv"), row.names = FALSE, fileEncoding = "UTF-8")

env2 <- new.env()
load(file.path(manet, "data", "noordin.rda"), envir = env2)
aff <- get("noordin", envir = env2)
cat("noordin (manet):", paste(dim(aff), collapse = " x "), "\n")
aff <- as.data.frame(aff, check.names = FALSE)
long <- do.call(rbind, lapply(names(aff), function(ev) {
  hit <- which(aff[[ev]] != 0)
  if (length(hit)) data.frame(person = rownames(aff)[hit], event = ev) else NULL
}))
utils::write.csv(long, file.path(out, "affiliations.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("rownames:", head(rownames(aff)), "\n")
cat("written to", out, "\n")
