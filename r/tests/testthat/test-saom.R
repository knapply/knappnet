# Stochastic actor-oriented models: data construction, effect parsing,
# hand-checked statistics, simulation, estimation (reproducibility,
# recovery of generating parameters), goodness of fit, and a comparison
# with RSiena when it is installed. Panels are tiny and phases short so the
# tests run in seconds; the full comparison lives in validation/.

# A directed netframe on a given node table.
net_on <- function(nodes, from, to) {
  netframe(data.frame(from = from, to = to, stringsAsFactors = FALSE), nodes = nodes)
}

# The netframe of a logical adjacency matrix.
net_of <- function(a, nodes) {
  idx <- which(a, arr.ind = TRUE)
  net_on(nodes, nodes[[1]][idx[, 1]], nodes[[1]][idx[, 2]])
}

three <- function() {
  # 0 -> 1, 1 -> 2, 0 -> 2, 2 -> 0 and an isolate: one transitive triplet,
  # one 3-cycle, one reciprocated pair.
  nodes <- data.frame(id = c("a", "b", "c", "d"), v = c(1, 2, 6, 2), g = c(0, 1, 1, 0),
                      stringsAsFactors = FALSE)
  w1 <- netframe(data.frame(from = "a", to = "b"), nodes = nodes)
  w2 <- netframe(data.frame(from = c("a", "b", "a", "c"), to = c("b", "c", "c", "a")),
                 nodes = nodes)
  saom_data(list(w1, w2), covariates = c("v", "g"))
}

# A first wave and later waves that each toggle about 12% of the ties.
drifting <- function(n, waves, seed, density = 0.2) {
  set.seed(seed)
  nodes <- data.frame(id = seq_len(n), x = round(rnorm(n), 1))
  edges <- function(a) {
    idx <- which(a, arr.ind = TRUE)
    netframe(data.frame(from = idx[, 1], to = idx[, 2]), nodes = nodes)
  }
  a <- matrix(runif(n * n) < density, n, n)
  diag(a) <- FALSE
  nets <- list(edges(a))
  for (w in seq_len(waves - 1)) {
    flip <- matrix(runif(n * n) < 0.12, n, n)
    diag(flip) <- FALSE
    a <- xor(a, flip)
    nets[[w + 1]] <- edges(a)
  }
  saom_data(nets, covariates = "x")
}

test_that("saom_data reads panels and covariates", {
  d <- three()
  expect_s3_class(d, "saom_data")
  expect_equal(d$n, 4L)
  expect_equal(sum(d$waves[[2]]), 4)
  expect_equal(names(d$covariates), c("v", "g"))
  expect_false(is.matrix(d$covariates$v))
  expect_output(print(d), "2 waves")
  # Changing covariates come from node columns that differ between waves.
  nodes2 <- data.frame(id = c("a", "b", "c", "d"), v = c(1, 2, 6, 3))
  w3 <- netframe(data.frame(from = "a", to = "c"), nodes = nodes2)
  d2 <- saom_data(list(net_on(three()$nodes, "a", "b"), w3), covariates = "v")
  expect_true(is.matrix(d2$covariates$v))
  expect_equal(d2$covariates$v[4, ], c(2, 3))
  # Nodes are matched by identifier, whatever their order.
  shuffled <- netframe(data.frame(from = "c", to = "a"),
                       nodes = three()$nodes[c(3, 1, 4, 2), ])
  d3 <- saom_data(list(net_on(three()$nodes, "a", "b"), shuffled))
  expect_true(d3$waves[[2]]["c", "a"])
  expect_error(saom_data(list(netframe(data.frame(from = 1, to = 2), directed = FALSE),
                              netframe(data.frame(from = 1, to = 2), directed = FALSE))),
               "directed")
  expect_error(saom_data(list(netframe(data.frame(from = 1, to = 2)))), "at least two")
  expect_error(saom_data(list(netframe(data.frame(from = 1, to = 2)),
                              netframe(data.frame(from = 1, to = 3)))), "different nodes")
  expect_error(saom_data(list(net_on(three()$nodes, "a", "b"), w3),
                         covariates = list(v = c(1, NA, 3, 4))), "missing")
  dy <- matrix(1, 4, 4)
  d4 <- saom_data(list(w3, w3), dyadic = list(w = dy))
  expect_equal(names(d4$dyadic), "w")
})

test_that("effects are parsed from formulas, names, and tables", {
  eff <- saom_effects(~ density + recip + gwespFF(parameter = 120) + simX("v") + X(w))
  expect_equal(eff$effect, c("density", "recip", "gwespFF", "simX", "X"))
  expect_equal(eff$covariate, c(NA, NA, NA, "v", "w"))
  expect_equal(eff$parameter, c(NA, NA, 120, NA, NA))
  expect_equal(eff$name, c("density", "recip", "gwespFF.120", "simX.v", "X.w"))
  expect_equal(saom_effects(c("density", "recip"))$effect, c("density", "recip"))
  expect_equal(saom_effects(eff)$name, eff$name)
  expect_equal(saom_effects(~ gwespFF)$parameter, 69)
  expect_true(all(saom_effect_names() %in% c(eff$effect, "transTrip", "transTies", "transRecTrip",
                                             "cycle3", "gwespFB", "gwespBF", "gwespBB", "inPop",
                                             "inPopSqrt", "outPop", "outPopSqrt", "inAct",
                                             "inActSqrt", "outAct", "outActSqrt", "egoX", "altX",
                                             "sameX", "egoXaltX")))
  expect_error(saom_effects(~ density + nope), "Unknown effect")
  expect_error(saom_effects(~ simX), "needs a covariate")
  expect_error(saom_effects(~ density("v")), "does not take a covariate")
  expect_error(saom_effects(~ recip(parameter = 2)), "does not take a parameter")
  expect_error(saom_effects(~ density + density), "included twice")
  expect_error(saom_effects(~ simX("nope"), three()), "do not have")
})

test_that("target statistics match hand-computed values", {
  d <- three()
  t <- knappnet:::saom_targets(d, ~ density + recip + transTrip + transTies + cycle3 + inPop +
                                 outAct + gwespFF + egoX("v") + altX("v") + simX("v") +
                                 sameX("g"))
  # Covariate v = (1, 2, 6, 2): mean 2.75; centered (-1.75, -0.75, 3.25, -0.75).
  # Ties at wave 2: a->b, b->c, a->c, c->a. egoX: a has out-degree 2, b 1, c 1.
  v <- c(1, 2, 6, 2) - 2.75
  expect_equal(unname(t[1, c(1:10, 12)]), c(
    4, 2, 1, 1, 3, 6, 6, 1,
    2 * v[1] + v[2] + v[3],
    v[2] + v[3] + v[3] + v[1],
    1))
  # Similarity 1 - |v_i - v_j| / 5, centered by its mean over ordered pairs.
  vv <- c(1, 2, 6, 2)
  sim <- 1 - abs(outer(vv, vv, "-")) / 5
  mean_sim <- mean(sim[row(sim) != col(sim)])
  expect_equal(unname(t[1, "simX.v"]),
               (sim[1, 2] + sim[2, 3] + sim[1, 3] + sim[3, 1]) - 4 * mean_sim)
  # Change statistics are the differences of the statistics.
  c1 <- knappnet:::saom_contributions(d, ~ density + recip + transTrip + inPop, actor = 1,
                                      wave = 2)
  # a -> b: with a -> c and b -> c it forms the triplet a -> b -> c, a -> c
  # (as its i -> h leg); b's in-degree with the tie is 1.
  expect_equal(unname(c1["b", ]), c(1, 0, 1, 1))
  # a -> c closes a -> b -> c, is reciprocated by c -> a, and c's in-degree
  # with the tie is 2.
  expect_equal(unname(c1["c", ]), c(1, 1, 1, 2))
  expect_equal(unname(c1["a", ]), c(0, 0, 0, 0))
})

test_that("simulation stops at the observed distance and is reproducible", {
  d <- drifting(12, 3, 1)
  sims <- saom_simulate(d, ~ density + recip, theta = c(-1.5, 1), n = 2, seed = 5)
  expect_length(sims, 2)
  expect_length(sims[[1]], 2)
  expect_s3_class(sims[[1]][[1]], "netframe")
  observed <- sapply(1:2, function(m) sum(d$waves[[m]] != d$waves[[m + 1]]))
  for (m in 1:2) {
    end <- as.matrix(sims[[1]][[m]]) > 0
    expect_equal(sum(end != d$waves[[m]]), observed[m])
  }
  expect_equal(dim(attr(sims[[1]], "statistics")), c(2L, 2L))
  again <- saom_simulate(d, ~ density + recip, theta = c(-1.5, 1), n = 2, seed = 5)
  expect_identical(again[[2]][[1]]$edges, sims[[2]][[1]]$edges)
  expect_false(identical(sims[[1]][[1]]$edges, sims[[2]][[1]]$edges))
  timed <- saom_simulate(d, ~ density, theta = -1, rates = 3, seed = 1)
  expect_equal(attr(timed[[1]], "times"), c(1, 1))
  expect_error(saom_simulate(d, ~ density, theta = c(1, 2)), "one value per effect")
  expect_error(saom_simulate(d, ~ density, theta = 1, rates = -1), "positive")
})

test_that("estimation is reproducible, thread-independent, and recovers parameters", {
  # A panel generated from a known model.
  d0 <- drifting(25, 2, 2, density = 0.1)
  start <- net_of(d0$waves[[1]], d0$nodes)
  scaffold <- saom_data(list(start, start))
  truth <- c(-2, 1.5)
  sim <- saom_simulate(scaffold, ~ density + recip, theta = truth, rates = 5, seed = 11)[[1]]
  data <- saom_data(c(list(start), sim))
  fit <- saom_fit(data, ~ density + recip,
                  control = saom_control(n3 = 300, nsub = 3, seed = 4, threads = 1))
  expect_s3_class(fit, "saom_fit")
  expect_named(coef(fit), c("density", "recip"))
  expect_true(all(is.finite(fit$se)) && all(fit$se > 0))
  expect_equal(dim(vcov(fit)), c(2L, 2L))
  expect_equal(nrow(fit$rates), 1L)
  expect_true(fit$rates$estimate > 0)
  expect_lt(fit$tconv_max, 0.6)
  z <- (coef(fit) - truth) / fit$se
  expect_true(all(abs(z) < 3))
  same <- saom_fit(data, ~ density + recip,
                   control = saom_control(n3 = 300, nsub = 3, seed = 4, threads = 2))
  expect_equal(coef(same), coef(fit), tolerance = 1e-10)
  expect_equal(same$deviations, fit$deviations)
  expect_output(print(fit), "conditional estimation")
  expect_output(print(summary(fit)), "Covariance matrix")
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(fit))
  expect_invisible(plot(fit, which = "deviations"))
})

test_that("unconditional estimation and initial values work", {
  d <- drifting(14, 3, 3)
  fit <- saom_fit(d, ~ density + recip,
                  control = saom_control(conditional = FALSE, n3 = 150, nsub = 2, seed = 6,
                                         threads = 2))
  expect_named(coef(fit), c("rate.1", "rate.2", "density", "recip"))
  expect_equal(fit$rates$estimate, unname(coef(fit)[1:2]))
  expect_equal(unname(fit$targets[1:2]),
               sapply(1:2, function(m) sum(d$waves[[m]] != d$waves[[m + 1]])))
  from <- saom_fit(d, ~ density + recip, initial = c(NA, 0.5),
                   control = saom_control(n3 = 100, nsub = 1, seed = 7, threads = 2))
  expect_equal(from$trajectory$subphase[1], 1L)
  expect_error(saom_fit(d, ~ density, initial = c(1, 2)), "length 1")
  expect_error(saom_fit(d, ~ density, control = list()), "saom_control")
  expect_error(saom_control(n3 = 1), "n3")
  expect_error(saom_control(diagonalize = 2), "diagonalize")
  no_change <- saom_data(list(net_on(d$nodes, 1, 2), net_on(d$nodes, 1, 2)))
  expect_error(saom_fit(no_change, ~ density, control = saom_control(n3 = 10, threads = 1)),
               "no tie changed")
})

test_that("goodness of fit compares observed with simulated statistics", {
  d <- drifting(12, 3, 8)
  fit <- saom_fit(d, ~ density + recip,
                  control = saom_control(n3 = 120, nsub = 2, seed = 9, threads = 2))
  g <- saom_gof(fit, "indegree")
  expect_s3_class(g, "saom_gof")
  expect_true(g$p >= 0 && g$p <= 1)
  expect_equal(dim(g$simulated), c(120L, 9L))
  # Observed cumulative in-degree counts over waves 2 and 3.
  indeg <- c(colSums(d$waves[[2]]), colSums(d$waves[[3]]))
  expect_equal(unname(g$observed), sapply(0:8, function(k) sum(indeg <= k)))
  raw <- saom_gof(fit, "outdegree", levels = 0:3, cumulative = FALSE)
  outdeg <- c(rowSums(d$waves[[2]]), rowSums(d$waves[[3]]))
  expect_equal(unname(raw$observed), sapply(0:3, function(k) sum(outdeg == k)))
  tri <- saom_gof(fit, "triad_census")
  expect_equal(sum(tri$observed), 2 * choose(12, 3))
  expect_equal(unname(tri$observed),
               unname(triad_census(net_of(d$waves[[2]], d$nodes), threads = 2) +
                        triad_census(net_of(d$waves[[3]], d$nodes), threads = 2)))
  geo <- saom_gof(fit, "geodesic")
  expect_equal(unname(geo$observed[length(geo$observed)]), 2 * 12 * 11)
  expect_output(print(g), "Mahalanobis")
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(g))
  expect_invisible(plot(tri, violin = FALSE))
  bare <- saom_fit(d, ~ density, control = saom_control(n3 = 20, nsub = 0, seed = 1, threads = 1,
                                                         keep_simulations = FALSE))
  expect_error(saom_gof(bare), "kept no simulations")
})

test_that("statistics and estimates agree with RSiena", {
  skip_if_not_installed("RSiena")
  skip_on_cran()
  d <- drifting(30, 3, 21, density = 0.1)
  n <- d$n
  effects <- ~ density + recip + transTrip + inPopSqrt + simX("x")
  targets <- knappnet:::saom_targets(d, effects)

  net <- RSiena::sienaDependent(array(as.numeric(unlist(d$waves)), dim = c(n, n, 3)))
  x <- RSiena::coCovar(d$covariates$x)
  dat <- RSiena::sienaDataCreate(net, x)
  eff <- RSiena::getEffects(dat)
  eff <- RSiena::includeEffects(eff, transTrip, inPopSqrt, verbose = FALSE)
  eff <- RSiena::includeEffects(eff, simX, interaction1 = "x", verbose = FALSE)
  # Target statistics: deterministic, so they must agree exactly.
  alg0 <- RSiena::sienaAlgorithmCreate(nsub = 0, n3 = 5, cond = FALSE, simOnly = TRUE,
                                       seed = 1, projname = NULL)
  ans0 <- RSiena::siena07(alg0, data = dat, effects = eff, batch = TRUE, silent = TRUE)
  re <- ans0$requestedEffects
  pick <- function(short) which(re$shortName == short & re$type == "eval")
  for (k in seq_len(ncol(targets))) {
    short <- saom_effects(effects)$effect[k]
    expect_equal(unname(targets[, k]), as.numeric(ans0$targets2[pick(short), ]),
                 tolerance = 1e-8, label = short)
  }
  # Estimates: stochastic, so they must agree within the standard errors.
  fit <- saom_fit(d, effects, control = saom_control(n3 = 1000, seed = 3, threads = 2))
  alg <- RSiena::sienaAlgorithmCreate(seed = 2, projname = NULL)
  ans <- RSiena::siena07(alg, data = dat, effects = eff, batch = TRUE, silent = TRUE)
  theirs <- ans$theta
  their_se <- ans$se
  names(theirs) <- names(their_se) <- re$shortName[re$type == "eval"]
  for (k in seq_len(ncol(targets))) {
    short <- saom_effects(effects)$effect[k]
    tol <- 2 * sqrt(fit$se[k]^2 + their_se[[short]]^2)
    expect_lt(abs(coef(fit)[[k]] - theirs[[short]]), tol, label = short)
  }
  if (!is.null(ans$rate)) {
    for (m in 1:2) {
      expect_lt(abs(fit$rates$estimate[m] - ans$rate[m]), 0.5 * ans$rate[m])
    }
  }
})
