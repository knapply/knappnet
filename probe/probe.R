suppressPackageStartupMessages(library(RSiena))
cat("RSiena", as.character(packageVersion("RSiena")), "\n")
n <- 5
nets <- list(A = list(c(1, 2)), B = list(c(1, 2), c(3, 2)), C = list(c(1, 2), c(2, 3)),
             D = list(c(1, 2), c(1, 3)), E = list(c(1, 2), c(2, 1)),
             F = list(c(1, 2), c(1, 3), c(1, 4)), G = list(c(2, 1), c(3, 1), c(4, 1)),
             H = list(c(1, 2), c(2, 3), c(3, 1)), I = list(c(1, 2), c(3, 2), c(4, 2), c(2, 5)))
names_eff <- c("outPop", "outPopSqrt", "inAct", "inActSqrt", "cycle3", "inPopSqrt", "outActSqrt")
out <- NULL
for (k in names(nets)) {
  w1 <- matrix(0, n, n); w1[4, 5] <- 1
  w2 <- matrix(0, n, n)
  for (e in nets[[k]]) w2[e[1], e[2]] <- 1
  net <- sienaDependent(array(c(w1, w2), dim = c(n, n, 2)), allowOnly = FALSE)
  dat <- sienaDataCreate(net)
  eff <- getEffects(dat)
  for (nm in names_eff) eff <- includeEffects(eff, nm, character = TRUE, verbose = FALSE)
  alg <- sienaAlgorithmCreate(nsub = 0, n3 = 5, cond = FALSE, simOnly = TRUE, seed = 1,
                              projname = NULL)
  ans <- suppressWarnings(siena07(alg, data = dat, effects = eff, batch = TRUE, silent = TRUE))
  re <- ans$requestedEffects
  v <- setNames(ans$targets2[, 1], re$shortName)
  out <- rbind(out, c(net = k, round(v[c("density", names_eff)], 4)))
}
print(as.data.frame(out), right = FALSE)
eff <- getEffects(dat)
print(eff[eff$shortName %in% names_eff & eff$type == "eval",
          c("shortName", "effectName", "parm", "functionName")])
