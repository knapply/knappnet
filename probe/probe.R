suppressPackageStartupMessages(library(RSiena))
toggle <- function(a, mult, add, mod) {
  n <- nrow(a)
  for (i in 0:(n - 1)) for (j in 0:(n - 1))
    if (i != j && (mult * (i * n + j) + add) %% mod == 0) a[i + 1, j + 1] <- 1 - a[i + 1, j + 1]
  a
}
w1 <- matrix(0, 2, 2); w1[1, 2] <- 1
w2 <- toggle(w1, 37, 11, 7); w3 <- toggle(w2, 53, 5, 9)
print(list(w1, w2, w3))
net <- sienaDependent(array(c(w1, w2, w3), dim = c(2, 2, 3)))
dat <- sienaDataCreate(net)
eff <- getEffects(dat)
print(eff[eff$include, c("shortName", "type")])
r <- tryCatch({
  eff2 <- includeEffects(eff, "density", character = TRUE, verbose = FALSE)
  alg <- sienaAlgorithmCreate(nsub = 0, n3 = 5, cond = FALSE, simOnly = TRUE, seed = 1, projname = NULL)
  ans <- siena07(alg, data = dat, effects = eff2, batch = TRUE, silent = TRUE)
  print(ans$requestedEffects[, c("shortName", "type")]); print(ans$targets2)
}, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
