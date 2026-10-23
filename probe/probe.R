j <- function(url) { tf <- tempfile(); ok <- tryCatch({download.file(url, tf, quiet = TRUE); TRUE}, error = function(e) {cat("ERR", url, conditionMessage(e), "\n"); FALSE}); if (ok) jsonlite::fromJSON(tf, simplifyVector = FALSE) }
ls <- j("https://cran.r-universe.dev/api/ls"); cat("cran ls:", length(ls), "\n")
a <- j("https://cran.r-universe.dev/api/packages?limit=50000&fields=Package"); cat("packages limit 50000:", length(a), "\n")
b <- j("https://cran.r-universe.dev/api/packages?limit=1000&skip=1000&fields=Package"); cat("skip works:", length(b), if (length(b)) b[[1]]$Package, "\n")
s <- j("https://r-universe.dev/api/search?q=exports:net_betweenness&limit=5"); str(s, max.level = 2)
s2 <- j("https://r-universe.dev/api/search?q=net_&limit=3"); cat(names(s2), "\n")
x <- j("https://cran.r-universe.dev/api/packages/igraph"); cat("igraph exports:", length(x[["_exports"]]), "\n")
