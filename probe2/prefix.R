j <- function(url) { tf <- tempfile(); download.file(url, tf, quiet = TRUE); jsonlite::fromJSON(tf, simplifyVector = FALSE) }
s <- j("https://r-universe.dev/api/search?q=exports:net_betweenness&limit=5")
r <- s$results[[1]]; cat("net_betweenness is exported by:", r$Package, "from", r[["_user"]], "\n")
all <- list()
for (u in c("cran", "bioc")) {
  skip <- 0
  repeat {
    b <- j(sprintf("https://%s.r-universe.dev/api/packages?limit=1000&skip=%d&fields=Package,_exports", u, skip))
    if (!length(b)) break
    all <- c(all, b); skip <- skip + 1000
  }
}
cat("packages read:", length(all), "\n")
ex <- do.call(rbind, lapply(all, function(p) {
  e <- unlist(p[["_exports"]]); if (!length(e)) return(NULL)
  data.frame(pkg = p$Package, fn = e)
}))
ex <- unique(ex)
cat("exports read:", nrow(ex), "\n")
for (pre in c("net_", "nw_", "nx_", "kn_", "nt_", "nk_", "ntwk_", "sn_", "kn_")) {
  hit <- ex[startsWith(ex$fn, pre), ]
  cat(sprintf("\n== prefix %-6s %4d functions in %3d packages", pre, nrow(hit), length(unique(hit$pkg))))
  if (nrow(hit)) {
    tab <- sort(table(hit$pkg), decreasing = TRUE)
    cat(":", paste(sprintf("%s (%d)", names(tab), tab), collapse = ", "))
  }
}
ours <- readLines("probe2/exports.txt")
cat("\n\n== net_<our name> already exported elsewhere:\n")
hit <- ex[ex$fn %in% paste0("net_", ours), ]; print(hit, row.names = FALSE)
cat("\n== every net_ export on CRAN/Bioconductor:\n")
h <- ex[startsWith(ex$fn, "net_"), ]; print(h[order(h$pkg, h$fn), ], row.names = FALSE)
common <- c("igraph", "sna", "network", "tidygraph", "ggraph", "dplyr", "stats", "base", "ergm",
            "RSiena", "sf", "sfnetworks", "networkDynamic", "tsna", "keyplayer", "netrankr",
            "centiserve", "qgraph", "purrr", "tidyr", "Matrix", "utils", "graphics")
x <- ex[ex$fn %in% ours & ex$pkg %in% common, ]
cat("\n== current unprefixed names that collide with commonly loaded packages:", nrow(x), "\n")
print(tapply(x$fn, x$pkg, function(v) paste(sort(v), collapse = ", ")))
cat("\n== how many of our current names exist anywhere on CRAN/Bioc:", sum(ours %in% ex$fn), "of", length(ours), "\n")
