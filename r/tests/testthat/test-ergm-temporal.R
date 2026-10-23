# Temporal ERGMs by pooled pseudolikelihood: the pooled fit against glm()
# on the stacked change statistics, memory terms against explicit dyadic
# covariates, bootstrap reproducibility, and btergm when installed.

waves <- function(n = 12, seed = 1) {
  set.seed(seed)
  wave <- function(prev, k) {
    e <- prev[-sample(nrow(prev), k), ]
    add <- data.frame(from = sample(n, k, TRUE), to = sample(n, k, TRUE))
    add <- add[add$from != add$to, ]
    unique(rbind(e, add))
  }
  e1 <- data.frame(from = sample(n, 35, TRUE), to = sample(n, 35, TRUE))
  e1 <- unique(e1[e1$from != e1$to, ])
  e2 <- wave(e1, 6)
  e3 <- wave(e2, 6)
  nodes <- data.frame(id = 1:n, group = rep(c("a", "b"), n / 2))
  lapply(list(e1, e2, e3), netframe, nodes = nodes, directed = TRUE)
}

test_that("the pooled pseudolikelihood matches glm() on stacked change statistics", {
  nets <- waves()
  fit <- tergm_fit(nets ~ edges + mutual + memory("stability") + delrecip + nodematch("group"),
                   boot = 20, seed = 1, threads = 2)
  expect_named(coef(fit), c("edges", "mutual", "edgecov.memory.stability", "edgecov.delrecip",
                            "nodematch.group"))
  xs <- list()
  ys <- list()
  for (t in 2:3) {
    prev <- as.matrix(nets[[t - 1]])
    prev <- (unname(prev) != 0) * 1
    mem <- 2 * prev - 1
    rec <- t(prev)
    f <- stats::as.formula("y ~ edges + mutual + edgecov(mem) + edgecov(rec) + nodematch('group')")
    y <- nets[[t]]
    model <- knappnet:::ergm_model(f, NULL, x = y, env = environment())
    cs <- knappnet:::ergm_changestats(model)
    xs[[t - 1]] <- cs$x
    ys[[t - 1]] <- cs$y
  }
  X <- do.call(rbind, xs)
  Y <- unlist(ys)
  ref <- stats::glm(Y ~ X - 1, family = stats::binomial())
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-6)
  expect_equal(fit$transitions, 2L)
  expect_equal(fit$dyads, c(132, 132))
  expect_equal(dim(fit$boot), c(20, 5))
  expect_equal(dim(fit$ci), c(5, 2))
  expect_equal(confint(fit), fit$ci)
  expect_equal(nrow(confint(fit, parm = "edges")), 1)
  expect_output(print(summary(fit)), "percentile")
  # The same seed gives the same replicates whatever the thread count.
  again <- tergm_fit(nets ~ edges + mutual + memory("stability") + delrecip + nodematch("group"),
                     boot = 20, seed = 1, threads = 1)
  expect_equal(fit$boot, again$boot)
})

test_that("memory types are the documented covariates", {
  nets <- waves(seed = 2)
  prev <- (unname(as.matrix(nets[[1]])) != 0) * 1
  for (type in c("stability", "autoregression", "innovation", "loss")) {
    cov <- switch(type, stability = 2 * prev - 1, autoregression = prev, innovation = 1 - prev,
                  loss = -prev)
    a <- tergm_fit(nets[1:2] ~ edges + memory(type), boot = 0, seed = 1, threads = 2)
    y <- nets[[2]]
    b <- ergm_fit(y ~ edges + edgecov(cov), estimate = "MPLE")
    expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-8, info = type)
  }
  # Memory terms are temporal only, and delrecip needs direction.
  expect_error(ergm_stats(nets[[1]] ~ memory("stability")), "tergm_fit")
  u <- lapply(nets, function(g) netframe(g$edges, g$nodes, directed = FALSE))
  expect_error(tergm_fit(u ~ edges + delrecip, boot = 0), "directed")
  expect_error(tergm_fit(nets[[1]] ~ edges, boot = 0), "list of at least two")
  expect_error(tergm_fit(nets ~ edges, boot = -1), "whole number")
})

test_that("pooled estimates match btergm", {
  skip_if_not_installed("btergm")
  skip_if_not_installed("network")
  nets <- waves(seed = 3)
  to_net <- function(g) {
    net <- network::network(as.matrix(g), directed = TRUE)
    network::set.vertex.attribute(net, "group", g$nodes$group)
    net
  }
  bnets <- lapply(nets, to_net)
  ours <- tergm_fit(nets ~ edges + mutual + memory("stability") + delrecip + nodematch("group"),
                    boot = 0, seed = 1, threads = 2)
  theirs <- suppressMessages(btergm::btergm(
    bnets ~ edges + mutual + memory(type = "stability") + delrecip + nodematch("group"),
    R = 10, verbose = FALSE))
  expect_equal(unname(coef(ours)), unname(theirs@coef), tolerance = 1e-5)
})
