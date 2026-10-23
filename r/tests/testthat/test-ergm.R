# Exponential random graph models: statistics against hand counts and
# ergm::summary(); the pseudolikelihood against glm() and ergm(estimate =
# "MPLE"); the Monte Carlo MLE and path-sampled log-likelihood against exact
# enumeration on a tiny network; simulation, prediction, goodness of fit,
# and diagnostics.

florentine_wealth <- function() {
  flo <- example_network("florentine")
  wealth <- c(Acciaiuoli = 10, Albizzi = 36, Barbadori = 55, Bischeri = 44, Castellani = 20,
              Ginori = 32, Guadagni = 8, Lamberteschi = 42, Medici = 103, Pazzi = 48,
              Peruzzi = 49, Ridolfi = 27, Salviati = 10, Strozzi = 146, Tornabuoni = 48)
  flo$nodes$wealth <- wealth[flo$nodes$name]
  flo$nodes$side <- ifelse(flo$nodes$wealth > 40, "rich", "poor")
  flo
}

# Two triangles sharing node 1 plus a pendant: hand-countable statistics.
bowtie <- function() {
  netframe(data.frame(from = c(1, 2, 1, 1, 4, 1), to = c(2, 3, 3, 4, 5, 5)),
           nodes = data.frame(id = 1:6, x = c(1, 2, 3, 4, 5, 6),
                              g = c("a", "a", "b", "b", "a", "b")),
           directed = FALSE)
}

directed_toy <- function() {
  set.seed(7)
  e <- data.frame(from = sample(10, 40, TRUE), to = sample(10, 40, TRUE))
  e <- unique(e[e$from != e$to, ])
  netframe(e, nodes = data.frame(id = 1:10, x = (1:10) %% 4, g = rep(c("u", "v", "w"), length.out = 10)),
           directed = TRUE)
}

gw <- function(a, k) ifelse(k > 0, exp(a) * (1 - (1 - exp(-a))^k), 0)

test_that("statistics match hand counts on a small network", {
  g <- bowtie()
  s <- ergm_stats(g ~ edges + triangle + kstar(2) + kstar(3) + degree(1) + degree(2) + isolates +
                    nodecov("x") + absdiff("x") + absdiff("x", pow = 2) + nodematch("g") +
                    nodematch("g", diff = TRUE) + nodefactor("g") + gwesp(0.5, fixed = TRUE) +
                    gwdegree(0.5, fixed = TRUE) + gwdsp(0.5, fixed = TRUE))
  deg <- c(4, 2, 2, 2, 2, 0)
  expect_equal(unname(s["edges"]), 6)
  expect_equal(unname(s["triangle"]), 2)
  expect_equal(unname(s["kstar2"]), sum(choose(deg, 2)))
  expect_equal(unname(s["kstar3"]), sum(choose(deg, 3)))
  expect_equal(unname(s["degree1"]), 0)
  expect_equal(unname(s["degree2"]), 4)
  expect_equal(unname(s["isolates"]), 1)
  expect_equal(unname(s["nodecov.x"]), sum(deg * 1:6))
  expect_equal(unname(s["absdiff.x"]), 1 + 1 + 2 + 3 + 1 + 4)
  expect_equal(unname(s["absdiff2.x"]), 1 + 1 + 4 + 9 + 1 + 16)
  # Ties 1-2 (a, a) and 1-5 (a, a) match.
  expect_equal(unname(s["nodematch.g"]), 2)
  expect_equal(unname(s[c("nodematch.g.a", "nodematch.g.b")]), c(2, 0))
  # nodefactor drops the first level: endpoints at level b are 3, 4, 6.
  expect_equal(unname(s["nodefactor.g.b"]), deg[3] + deg[4] + deg[6])
  expect_named(s, c("edges", "triangle", "kstar2", "kstar3", "degree1", "degree2", "isolates",
                    "nodecov.x", "absdiff.x", "absdiff2.x", "nodematch.g", "nodematch.g.a",
                    "nodematch.g.b", "nodefactor.g.b", "gwesp.fixed.0.5", "gwdeg.fixed.0.5",
                    "gwdsp.fixed.0.5"))
  # Every tie of a triangle has one shared partner.
  expect_equal(unname(s["gwesp.fixed.0.5"]), 6 * gw(0.5, 1))
  expect_equal(unname(s["gwdeg.fixed.0.5"]), sum(gw(0.5, deg)))
  # Dyadwise: the six ties have 1 partner each; dyads 2-4, 2-5, 3-4, 3-5
  # share node 1; the rest none.
  expect_equal(unname(s["gwdsp.fixed.0.5"]), 10 * gw(0.5, 1))
})

test_that("directed statistics match hand counts", {
  g <- netframe(data.frame(from = c(1, 2, 1, 3, 2), to = c(2, 1, 3, 2, 3)),
                nodes = data.frame(id = 1:4), directed = TRUE)
  s <- ergm_stats(g ~ edges + mutual + ttriple + ctriple + triangle + istar(2) + ostar(2) +
                    idegree(0) + odegree(2) + isolates + gwesp(0.3, fixed = TRUE))
  expect_equal(unname(s["edges"]), 5)
  expect_equal(unname(s["mutual"]), 2)
  # Transitive triples: (1->2, 2->3, 1->3), (1->3, 3->2, 1->2), and
  # (2->1, 1->3, 2->3); the one cycle is 1->3->2->1.
  expect_equal(unname(s["ttriple"]), 3)
  expect_equal(unname(s["ctriple"]), 1)
  expect_equal(unname(s["triangle"]), 4)
  expect_equal(unname(s["istar2"]), choose(2, 2) + choose(2, 2))
  expect_equal(unname(s["ostar2"]), choose(2, 2) + choose(2, 2))
  expect_equal(unname(s["idegree0"]), 1)
  expect_equal(unname(s["odegree2"]), 2)
  expect_equal(unname(s["isolates"]), 1)
  # Ties 1->2 (via 3), 1->3 (via 2), and 2->3 (via 1) each have one
  # outgoing two-path partner.
  expect_equal(unname(s["gwesp.OTP.fixed.0.3"]), 3 * gw(0.3, 1))
})

test_that("statistics match summary() in ergm", {
  skip_if_not_installed("ergm")
  skip_if_not_installed("network")
  flo <- florentine_wealth()
  net <- network::network(as.matrix(flo), directed = FALSE)
  network::set.vertex.attribute(net, "wealth", flo$nodes$wealth)
  network::set.vertex.attribute(net, "side", flo$nodes$side)
  ours <- ergm_stats(flo ~ edges + triangle + kstar(2) + kstar(3) + degree(1) + degree(2) +
                       isolates + nodecov("wealth") + absdiff("wealth") + nodematch("side") +
                       nodematch("side", diff = TRUE) + nodefactor("side") +
                       gwesp(0.5, fixed = TRUE) + gwdsp(0.25, fixed = TRUE) +
                       gwnsp(0.25, fixed = TRUE) + gwdegree(0.5, fixed = TRUE))
  theirs <- summary(net ~ edges + triangle + kstar(2) + kstar(3) + degree(1) + degree(2) +
                      isolates + nodecov("wealth") + absdiff("wealth") + nodematch("side") +
                      nodematch("side", diff = TRUE) + nodefactor("side") +
                      gwesp(0.5, fixed = TRUE) + gwdsp(0.25, fixed = TRUE) +
                      gwnsp(0.25, fixed = TRUE) + gwdegree(0.5, fixed = TRUE))
  expect_equal(unname(ours), unname(as.numeric(theirs)), tolerance = 1e-10)

  d <- directed_toy()
  net <- network::network(as.matrix(d), directed = TRUE)
  network::set.vertex.attribute(net, "x", d$nodes$x)
  network::set.vertex.attribute(net, "g", d$nodes$g)
  w <- matrix(seq(0, 1, length.out = 100), 10, 10)
  ours <- ergm_stats(d ~ edges + mutual + ttriple + ctriple + triangle + istar(2) + ostar(2) +
                       idegree(1) + odegree(1) + isolates + nodeicov("x") + nodeocov("x") +
                       nodeifactor("g") + nodeofactor("g") + edgecov(w) +
                       gwesp(0.5, fixed = TRUE) + gwdsp(0.5, fixed = TRUE) +
                       gwidegree(0.5, fixed = TRUE) + gwodegree(0.5, fixed = TRUE))
  theirs <- summary(net ~ edges + mutual + ttriple + ctriple + triangle + istar(2) + ostar(2) +
                      idegree(1) + odegree(1) + isolates + nodeicov("x") + nodeocov("x") +
                      nodeifactor("g") + nodeofactor("g") + edgecov(w) +
                      gwesp(0.5, fixed = TRUE) + gwdsp(0.5, fixed = TRUE) +
                      gwidegree(0.5, fixed = TRUE) + gwodegree(0.5, fixed = TRUE))
  expect_equal(unname(ours), unname(as.numeric(theirs)), tolerance = 1e-10)
})

test_that("the pseudolikelihood matches glm() on the change statistics", {
  flo <- florentine_wealth()
  fit <- ergm_fit(flo ~ edges + triangle + nodecov("wealth") + nodematch("side"), estimate = "MPLE")
  cs <- knappnet:::ergm_changestats(fit$model)
  # glm()'s default convergence rule leaves standard errors accurate to
  # about 1e-5, so it is tightened here.
  ref <- stats::glm(cs$y ~ cs$x - 1, family = stats::binomial(),
                    control = stats::glm.control(epsilon = 1e-14, maxit = 100))
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(fit)))), unname(sqrt(diag(vcov(ref)))), tolerance = 1e-6)
  expect_equal(fit$mple$deviance, ref$deviance, tolerance = 1e-8)
  expect_equal(fit$dyads, choose(15, 2))
  expect_false(fit$model$dyad_independent)
  expect_null(fit$loglik)
  expect_error(logLik(fit), "not available")
  s <- summary(fit)
  expect_true(all(c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in% colnames(s$coefficients)))
  expect_output(print(s), "not valid")

  # Dyad-independent: the MPLE is the MLE and the log-likelihood is available.
  fit2 <- ergm_fit(flo ~ edges + nodecov("wealth"), estimate = "MCMLE")
  expect_equal(fit2$estimate, "MPLE")
  expect_true(fit2$model$dyad_independent)
  expect_equal(as.numeric(logLik(fit2)), -fit2$mple$deviance / 2)
  expect_equal(attr(logLik(fit2), "df"), 2)
  expect_output(print(summary(fit2)), "AIC")
})

test_that("the pseudolikelihood matches ergm(estimate = 'MPLE')", {
  skip_if_not_installed("ergm")
  skip_if_not_installed("network")
  flo <- florentine_wealth()
  net <- network::network(as.matrix(flo), directed = FALSE)
  network::set.vertex.attribute(net, "wealth", flo$nodes$wealth)
  ours <- ergm_fit(flo ~ edges + gwesp(0.5, fixed = TRUE) + nodecov("wealth"), estimate = "MPLE")
  theirs <- suppressMessages(ergm::ergm(net ~ edges + gwesp(0.5, fixed = TRUE) + nodecov("wealth"),
                                        estimate = "MPLE"))
  expect_equal(unname(coef(ours)), unname(coef(theirs)), tolerance = 1e-5)
  d <- directed_toy()
  net <- network::network(as.matrix(d), directed = TRUE)
  network::set.vertex.attribute(net, "g", d$nodes$g)
  ours <- ergm_fit(d ~ edges + mutual + nodematch("g"), estimate = "MPLE")
  theirs <- suppressMessages(ergm::ergm(net ~ edges + mutual + nodematch("g"), estimate = "MPLE"))
  expect_equal(unname(coef(ours)), unname(coef(theirs)), tolerance = 1e-5)
})

test_that("exact estimation agrees with closed forms", {
  g <- bowtie()
  fit <- ergm_fit(g ~ edges, estimate = "exact")
  p <- 6 / 15
  expect_equal(unname(coef(fit)), log(p / (1 - p)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(fit)), 6 * log(p) + 9 * log(1 - p), tolerance = 1e-8)
  expect_equal(unname(vcov(fit)[1, 1]), 1 / (15 * p * (1 - p)), tolerance = 1e-6)
  expect_equal(fit$null_loglik, -15 * log(2))
  # Statistics on the boundary have no MLE.
  full <- netframe(data.frame(from = c(1, 1, 1, 2, 2, 3), to = c(2, 3, 4, 3, 4, 4)), directed = FALSE)
  expect_error(ergm_fit(full ~ edges, estimate = "exact"), "does not exist")
  big <- example_network("florentine")
  expect_error(ergm_fit(big ~ edges, estimate = "exact"), "22 free dyads")
})

test_that("the Monte Carlo MLE recovers the exact MLE within Monte Carlo error", {
  g <- bowtie()
  exact <- ergm_fit(g ~ edges + gwesp(0.5, fixed = TRUE), estimate = "exact")
  ctl <- ergm_control(seed = 3, threads = 2, burnin = 2000, interval = 20, nsim = 1500,
                      loglik_steps = 8)
  fit <- ergm_fit(g ~ edges + gwesp(0.5, fixed = TRUE), control = ctl)
  expect_equal(fit$estimate, "MCMLE")
  expect_true(fit$converged)
  se <- sqrt(diag(vcov(fit)))
  # Both the estimate and its exact counterpart are within a fraction of a
  # standard error; the bound is loose enough never to fail by chance.
  expect_lt(max(abs(coef(fit) - coef(exact)) / se), 0.5)
  expect_true(all(diag(fit$mcmc_vcov) >= 0))
  expect_equal(dim(fit$sample), c(2 * 1500, 2))
  # Path-sampled log-likelihood versus the exact one at the same coefficients.
  ptr <- knappnet:::ergm_handle(fit$model)
  ex <- .Call(knappnet:::C_ergm_exact, ptr, unname(coef(fit)))
  expect_lt(abs(as.numeric(logLik(fit)) - ex[[1]]), 4 * fit$loglik_se + 0.1)
  s <- summary(fit)
  expect_true("MCMC %" %in% colnames(s$coefficients))
  expect_output(print(s), "AIC")
  expect_output(print(fit), "Monte Carlo maximum likelihood")
})

test_that("a separated pseudolikelihood starts the Monte Carlo fit from zero", {
  # Two disjoint triangles: every tie has a shared partner and no non-tie
  # does, so the MPLE of edges + triangle is infinite while the MLE exists.
  g <- netframe(data.frame(from = c(1, 2, 3, 4, 5, 6), to = c(2, 3, 1, 5, 6, 4)),
                nodes = data.frame(id = 1:7), directed = FALSE)
  exact <- ergm_fit(g ~ edges + triangle, estimate = "exact")
  ctl <- ergm_control(seed = 5, threads = 2, burnin = 2000, interval = 20, nsim = 1500,
                      loglik = FALSE)
  expect_warning(fit <- ergm_fit(g ~ edges + triangle, control = ctl), "starts from zero")
  expect_lt(max(abs(coef(fit) - coef(exact))), 0.25)
})

test_that("simulation draws Bernoulli networks under an edges model", {
  flo <- example_network("florentine")
  theta <- -1
  nets <- ergm_simulate(flo ~ edges, coef = theta, nsim = 200, seed = 11, burnin = 1000,
                        interval = 20, threads = 2)
  expect_length(nets, 200)
  expect_true(all(vapply(nets, inherits, TRUE, "netframe")))
  counts <- attr(nets, "stats")[, "edges"]
  expect_equal(counts, vapply(nets, function(n) nrow(n$edges), 1))
  # 105 dyads with tie probability logistic(-1): mean 28.2, sd 4.4; the
  # mean of 200 draws (autocorrelated, interval 20) is within 2.
  expect_lt(abs(mean(counts) - 105 * plogis(theta)), 2)
  expect_identical(nets[[1]]$nodes, flo$nodes)
  # Reproducible by seed and independent of the thread count.
  again <- ergm_simulate(flo ~ edges, coef = theta, nsim = 5, seed = 11, burnin = 1000,
                         interval = 20, threads = 1, output = "stats")
  expect_equal(unname(again[, 1]), unname(counts[1:5]))
  fit <- ergm_fit(flo ~ edges, estimate = "MPLE")
  s <- simulate(fit, nsim = 3, seed = 2, output = "stats", threads = 2)
  expect_equal(dim(s), c(3, 1))
  expect_error(ergm_simulate(flo ~ edges, coef = c(1, 2)), "one value per statistic")
})

test_that("fixed dyads keep their state", {
  g <- bowtie()
  fixed <- matrix(FALSE, 6, 6)
  fixed[1, 2] <- fixed[2, 1] <- TRUE # tie, kept
  fixed[5, 6] <- fixed[6, 5] <- TRUE # non-tie, kept
  fit <- ergm_fit(g ~ edges, estimate = "MPLE", constraints = fixed)
  expect_equal(fit$dyads, 13)
  nets <- simulate(fit, nsim = 20, seed = 4, coef = 1, burnin = 500, interval = 10, threads = 2)
  for (net in nets) {
    m <- as.matrix(net)
    expect_equal(m["1", "2"], 1)
    expect_equal(m["5", "6"], 0)
  }
})

test_that("predict() gives conditional tie probabilities", {
  flo <- florentine_wealth()
  fit <- ergm_fit(flo ~ edges + nodecov("wealth"), estimate = "MPLE")
  p <- predict(fit)
  expect_equal(nrow(p), 105)
  expect_true(all(p$p > 0 & p$p < 1))
  # At the MLE of a dyad-independent model, expected statistics equal the
  # observed ones.
  expect_equal(sum(p$p), 20, tolerance = 1e-6)
  expect_equal(sum(p$tie), 20)
  link <- predict(fit, type = "link")
  expect_equal(plogis(link$link), p$p)
  expect_true(all(p$from %in% flo$nodes$name))
})

test_that("goodness of fit and diagnostics run and are well formed", {
  flo <- example_network("florentine")
  fit <- ergm_fit(flo ~ edges, estimate = "MPLE")
  g <- ergm_gof(fit, nsim = 30, seed = 1, burnin = 500, interval = 10, threads = 2)
  expect_s3_class(g, "ergm_gof")
  expect_named(g, c("degree", "esp", "dsp", "distance", "model", "nsim", "seed"))
  expect_equal(dim(g$degree$simulated), c(30, 15))
  expect_equal(unname(g$degree$observed), unname(tabulate(degree(flo) + 1, 15)))
  expect_equal(sum(g$esp$observed), 20)
  expect_equal(sum(g$dsp$observed), 105)
  expect_equal(sum(g$distance$observed), 105)
  expect_equal(unname(g$model$observed), 20)
  expect_true(all(g$model$summary$p >= 0 & g$model$summary$p <= 1))
  expect_output(print(g), "Degree")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(g))
  expect_error(plot(g, which = "triadcensus"), "Unknown statistic")

  d <- directed_toy()
  fitd <- ergm_fit(d ~ edges, estimate = "MPLE")
  gd <- ergm_gof(fitd, nsim = 10, seed = 2, burnin = 200, interval = 5, threads = 2)
  expect_true(all(c("idegree", "odegree", "triadcensus") %in% names(gd)))
  expect_equal(unname(gd$triadcensus$observed), unname(triad_census(d)))

  ctl <- ergm_control(seed = 1, threads = 2, burnin = 500, interval = 10, nsim = 300,
                      loglik = FALSE)
  fit2 <- ergm_fit(flo ~ edges + gwesp(0.5, fixed = TRUE), control = ctl)
  diag <- mcmc_diagnostics(fit2)
  expect_s3_class(diag, "ergm_mcmc")
  expect_equal(rownames(diag$summary), names(coef(fit2)))
  expect_true(all(diag$summary$ess > 0 & diag$summary$ess <= 600))
  expect_true(all(diag$summary$geweke_p >= 0 & diag$summary$geweke_p <= 1))
  expect_equal(dim(diag$autocorrelation), c(5, 2))
  expect_output(print(diag), "Geweke|geweke")
  expect_invisible(plot(diag))
  expect_error(mcmc_diagnostics(fit), "MCMLE")
})

test_that("formulas and arguments are checked", {
  g <- bowtie()
  expect_error(ergm_stats(g ~ edges + foo), "Unknown model term")
  expect_error(ergm_stats(g ~ nodecov("missing")), "not a column")
  expect_error(ergm_stats(g ~ gwesp(0.5, fixed = FALSE)), "fixed = TRUE")
  expect_error(ergm_stats(g ~ gwesp()), "decay")
  expect_error(ergm_stats(g ~ mutual), "directed networks")
  expect_error(ergm_stats(g ~ istar(2)), "directed networks")
  expect_error(ergm_stats(directed_toy() ~ kstar(2)), "undirected networks")
  expect_error(ergm_stats(g ~ edges(2)), "takes no arguments")
  expect_error(ergm_stats(g ~ nodefactor("g", levels = "z")), "not found")
  expect_error(ergm_stats(g ~ memory("stability")), "tergm_fit")
  expect_error(ergm_stats(~ edges), "two-sided")
  expect_error(ergm_stats(1 ~ edges), "netframe")
  loops <- netframe(data.frame(from = c(1, 2), to = c(1, 3)), directed = FALSE)
  expect_error(ergm_stats(loops ~ edges), "self-loops")
  expect_error(ergm_fit(g ~ edges, constraints = matrix(TRUE, 2, 2)), "6-by-6")
  expect_error(ergm_control(burnin = -1), "whole number")
  expect_error(ergm_control(margin = 1), "margin")
  expect_error(ergm_fit(g ~ edges, start = 1:2), "length 1")
  # Level selection by position and by name.
  expect_named(ergm_stats(g ~ nodefactor("g", levels = NULL)), c("nodefactor.g.a", "nodefactor.g.b"))
  expect_named(ergm_stats(g ~ nodefactor("g", levels = "a")), "nodefactor.g.a")
  expect_named(ergm_stats(g ~ kstar(2:3)), c("kstar2", "kstar3"))
  # edgecov from a netframe or a matrix.
  w <- as.matrix(g)
  expect_equal(unname(ergm_stats(g ~ edgecov(w))), 6)
  expect_equal(unname(ergm_stats(g ~ edgecov(g))), 6)
})
