suppressPackageStartupMessages(library(RSiena))
cat("RSiena", as.character(packageVersion("RSiena")), "\n")
# Wave 1: 0->1, 1->2, 2->0 (a 3-cycle) and 0->3; wave 2: add 3->0, 1->0.
n <- 4
w1 <- matrix(0, n, n); w1[cbind(c(1, 2, 3, 1), c(2, 3, 1, 4))] <- 1
w2 <- w1; w2[cbind(c(4, 2), c(1, 1))] <- 1
cat("wave 2 outdegrees", rowSums(w2), " indegrees", colSums(w2), "\n")
w3 <- w2; w3[2, 4] <- 1
net <- sienaDependent(array(c(w1, w2, w3), dim = c(n, n, 3)))
# Three columns given for three waves, as the validation oracle does.
ch <- varCovar(cbind(c(0, 1, 2, 3), c(1, 1, 3, 3), c(3, 0, 0, 0)))
dat <- sienaDataCreate(net, ch)
eff <- getEffects(dat)
for (nm in c("outPop", "outPopSqrt", "inAct", "inActSqrt", "cycle3", "inPop", "inPopSqrt",
             "outAct", "outActSqrt", "transTrip")) {
  eff <- includeEffects(eff, nm, character = TRUE, verbose = FALSE)
}
eff <- includeEffects(eff, egoX, altX, simX, interaction1 = "ch", verbose = FALSE)
alg <- sienaAlgorithmCreate(nsub = 0, n3 = 5, cond = FALSE, simOnly = TRUE, seed = 1,
                            projname = NULL)
ans <- siena07(alg, data = dat, effects = eff, batch = TRUE, silent = TRUE)
re <- ans$requestedEffects
print(data.frame(short = re$shortName, int = re$interaction1, parm = re$parm,
                 name = re$effectName, t1 = ans$targets2[, 1], t2 = ans$targets2[, 2]),
      right = FALSE)
cat("\ncovariate centering: mean attribute\n")
print(attributes(dat$vCovars$ch))
