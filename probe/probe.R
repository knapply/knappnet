remotes::install_github("knapply/knappnet", subdir = "r", upgrade = "never", quiet = TRUE)
suppressPackageStartupMessages({library(ergm); library(knappnet)})
cat("ergm", as.character(packageVersion("ergm")), " cores", parallel::detectCores(), "\n")
make <- function(n, deg, seed) {
  set.seed(seed)
  m <- n * deg / 2
  e <- unique(t(apply(matrix(sample(n, 3 * m, TRUE), ncol = 2), 1, sort)))
  e <- e[e[, 1] != e[, 2], , drop = FALSE][seq_len(m), ]
  g <- rep(c("a", "b"), length.out = n)
  list(e = e, g = g, n = n)
}
med <- function(f, k = 3) median(replicate(k, system.time(f())[["elapsed"]]))
row <- function(task, n, a, b) cat(sprintf("%-28s n=%-5d statnet %8.3fs  knappnet %8.3fs  speedup %6.1fx\n",
                                           task, n, a, b, a / b))
for (n in c(100, 300, 1000)) {
  d <- make(n, 6, n)
  net <- network::network(d$e, directed = FALSE, matrix.type = "edgelist")
  network::set.vertex.attribute(net, "g", d$g)
  nf <- netframe(data.frame(from = d$e[, 1], to = d$e[, 2]),
                 nodes = data.frame(id = seq_len(n), g = d$g), directed = FALSE)
  row("statistics", n, med(function() summary(net ~ edges + nodematch("g") + gwesp(0.5, fixed = TRUE))),
      med(function() ergm_stats(nf ~ edges + nodematch("g") + gwesp(0.5, fixed = TRUE))))
  a <- med(function() ergm(net ~ edges + nodematch("g") + gwesp(0.5, fixed = TRUE), estimate = "MPLE"))
  b <- med(function() ergm_fit(nf ~ edges + nodematch("g") + gwesp(0.5, fixed = TRUE), estimate = "MPLE"))
  row("MPLE", n, a, b)
  fm <- ergm(net ~ edges + nodematch("g") + gwesp(0.5, fixed = TRUE), estimate = "MPLE")
  th <- coef(fm)
  a <- system.time(simulate(net ~ edges + nodematch("g") + gwesp(0.5, fixed = TRUE), coef = th,
                            nsim = 100, output = "stats",
                            control = control.simulate.formula(MCMC.burnin = 1e4, MCMC.interval = 1e3)))[["elapsed"]]
  b <- system.time(ergm_simulate(nf ~ edges + nodematch("g") + gwesp(0.5, fixed = TRUE), coef = th,
                                 nsim = 100, burnin = 1e4, interval = 1e3, threads = 1,
                                 output = "stats"))[["elapsed"]]
  row("simulate 100 (1 thread)", n, a, b)
  ta <- system.time(fa <- ergm(net ~ edges + nodematch("g") + gwesp(0.5, fixed = TRUE),
                               control = control.ergm(MCMC.burnin = 1e4, MCMC.interval = 100,
                                                      MCMC.samplesize = 1000, seed = 1)))[["elapsed"]]
  tb1 <- system.time(fb <- ergm_fit(nf ~ edges + nodematch("g") + gwesp(0.5, fixed = TRUE),
                                    control = ergm_control(seed = 1, threads = 1, burnin = 1e4,
                                                           interval = 100, nsim = 1000,
                                                           loglik = FALSE)))[["elapsed"]]
  tb <- system.time(fb2 <- ergm_fit(nf ~ edges + nodematch("g") + gwesp(0.5, fixed = TRUE),
                                    control = ergm_control(seed = 1, threads = 0, burnin = 1e4,
                                                           interval = 100, nsim = 1000,
                                                           loglik = FALSE)))[["elapsed"]]
  row("MCMLE (1 thread)", n, ta, tb1)
  row("MCMLE (all cores)", n, ta, tb)
  cat("  estimates statnet :", round(coef(fa), 3), "\n")
  cat("  estimates knappnet:", round(coef(fb), 3), "\n")
  cat("  SE statnet :", round(sqrt(diag(vcov(fa))), 3), "\n")
  cat("  SE knappnet:", round(sqrt(diag(vcov(fb))), 3), "\n")
}
