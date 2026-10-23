# Which exported functions on CRAN and Bioconductor start with candidate
# prefixes, from r-universe's package metadata (field _exports).
get <- function(universe) {
  url <- sprintf("https://%s.r-universe.dev/api/packages?limit=100000&fields=Package,_exports", universe)
  tf <- tempfile(); download.file(url, tf, quiet = TRUE)
  jsonlite::fromJSON(tf, simplifyVector = FALSE)
}
pk <- c(get("cran"), get("bioc"))
cat("packages read:", length(pk), "\n")
ex <- do.call(rbind, lapply(pk, function(p) {
  e <- unlist(p[["_exports"]]); if (!length(e)) return(NULL)
  data.frame(pkg = p$Package, fn = e)
}))
cat("exports read:", nrow(ex), "\n")
for (pre in c("net_", "nw_", "nx_", "kn_", "knet_", "knap_", "nt_", "sna_", "ntw_", "nk_", "kp_")) {
  hit <- ex[startsWith(ex$fn, pre), ]
  cat(sprintf("\n== prefix %-6s %4d functions in %3d packages", pre, nrow(hit), length(unique(hit$pkg))))
  if (nrow(hit)) {
    tab <- sort(table(hit$pkg), decreasing = TRUE)
    cat(":", paste(sprintf("%s (%d)", names(tab), tab), collapse = ", "))
  }
}
ours <- readLines("probe2/exports.txt")
cat("\n\n== exact collisions of net_<our name> with any export:\n")
cand <- paste0("net_", ours)
hit <- ex[ex$fn %in% cand, ]
print(hit, row.names = FALSE)
cat("\n== exact collisions of our current (unprefixed) names with packages commonly loaded together:\n")
common <- c("igraph", "sna", "network", "statnet", "tidygraph", "ggraph", "dplyr", "stats", "ergm",
            "RSiena", "sf", "sfnetworks", "spdep", "networkDynamic", "tsna", "keyplayer", "influenceR",
            "netrankr", "centiserve", "qgraph", "visNetwork", "stringr", "purrr", "tidyr")
print(table(ex$pkg[ex$fn %in% ours & ex$pkg %in% common]))
x <- ex[ex$fn %in% ours & ex$pkg %in% common, ]; print(x[order(x$fn), ], row.names = FALSE)
