suppressPackageStartupMessages(library(RSiena))
n <- 4
w1 <- matrix(0, n, n); w1[cbind(c(1, 2, 3, 1), c(2, 3, 1, 4))] <- 1
w2 <- w1; w2[cbind(c(4, 2), c(1, 1))] <- 1
w3 <- w2; w3[2, 4] <- 1
net <- sienaDependent(array(c(w1, w2, w3), dim = c(n, n, 3)))
show <- function(label, m) {
  ch <- varCovar(m)
  dat <- sienaDataCreate(net, ch)
  a <- attributes(dat$vCovars$ch)
  cat(label, ": mean", a$mean, " range", a$range, " simMeans", a$simMeans, "\n")
}
show("A same columns 0..3", cbind(0:3, 0:3))
show("B col2 1,1,3,3", cbind(0:3, c(1, 1, 3, 3)))
show("C col2 all 0", cbind(0:3, c(0, 0, 0, 0)))
show("D col1 1,1,3,3 col2 0..3", cbind(c(1, 1, 3, 3), 0:3))
show("E three cols", cbind(0:3, c(1, 1, 3, 3), c(3, 0, 0, 0)))
co <- coCovar(0:3)
dat <- sienaDataCreate(net, co)
a <- attributes(dat$cCovars$co)
cat("constant 0..3: mean", a$mean, " range", a$range, " simMean", a$simMean, "\n")
