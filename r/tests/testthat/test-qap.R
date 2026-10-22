# QAP correlation and MRQAP. Coefficients and observed statistics must
# match lm() and sna exactly; Monte Carlo p-values must agree with sna
# within Monte Carlo error. Exhaustive enumeration is checked against an
# independent implementation in validation/qap_check.py.

qap_data <- function(n = 10, seed = 3) {
  set.seed(seed)
  x1 <- matrix(rnorm(n * n), n, n)
  x2 <- matrix(rbinom(n * n, 1, 0.4), n, n)
  y <- 0.6 * x1 - 0.4 * x2 + matrix(rnorm(n * n), n, n)
  list(y = y, x1 = x1, x2 = x2)
}

# Monte Carlo p-values from two generators: with 2000 sna draws and 20000
# here, the binomial standard error of the difference is at most 0.012, so
# 0.05 is a band of more than four standard errors.
expect_close <- function(a, b, tol = 0.05) {
  expect_lt(max(abs(as.numeric(a) - as.numeric(b))), tol)
}

# Vectorizes as sna does: column-major, lower triangle for "graph".
dyads <- function(m, mode = "digraph", diag = FALSE) {
  keep <- if (mode == "graph") lower.tri(m, diag = diag) else row(m) != col(m) | diag
  m[keep]
}

test_that("mrqap estimates match lm() on the vectorized dyads", {
  d <- qap_data()
  for (mode in c("digraph", "graph")) {
    for (diag in c(FALSE, TRUE)) {
      fit <- mrqap(d$y, list(a = d$x1, b = d$x2), reps = 0, mode = mode, diag = diag)
      ref <- summary(stats::lm(dyads(d$y, mode, diag) ~ dyads(d$x1, mode, diag) +
                                 dyads(d$x2, mode, diag)))
      expect_equal(unname(fit$coefficients), unname(ref$coefficients[, 1]), tolerance = 1e-10)
      expect_equal(unname(fit$std.error), unname(ref$coefficients[, 2]), tolerance = 1e-10)
      expect_equal(unname(fit$t), unname(ref$coefficients[, 3]), tolerance = 1e-10)
      expect_equal(fit$r.squared, ref$r.squared, tolerance = 1e-10)
      expect_equal(fit$adj.r.squared, ref$adj.r.squared, tolerance = 1e-10)
      expect_equal(fit$sigma, ref$sigma, tolerance = 1e-10)
      expect_equal(unname(fit$fstatistic), unname(ref$fstatistic), tolerance = 1e-10)
      expect_named(fit$coefficients, c("(Intercept)", "a", "b"))
    }
  }
})

test_that("mrqap drops incomplete dyads as lm() does", {
  d <- qap_data()
  d$y[c(3, 17, 40)] <- NA
  d$x1[c(5, 61)] <- NA
  fit <- mrqap(d$y, list(d$x1, d$x2), reps = 0)
  ref <- stats::lm(dyads(d$y) ~ dyads(d$x1) + dyads(d$x2))
  expect_equal(unname(fit$coefficients), unname(stats::coef(ref)), tolerance = 1e-10)
  expect_equal(fit$dyads, stats::nobs(ref))
})

test_that("mrqap matches sna::netlm()", {
  skip_if_not_installed("sna")
  d <- qap_data()
  x <- list(d$x1, d$x2)
  for (mode in c("digraph", "graph")) {
    yy <- d$y
    xx <- x
    if (mode == "graph") {
      sym <- function(m) { m[upper.tri(m)] <- t(m)[upper.tri(m)]; m }
      yy <- sym(yy)
      xx <- lapply(xx, sym)
    }
    for (nullhyp in c("qap", "qapy")) {
      set.seed(11)
      ref <- sna::netlm(yy, xx, mode = mode, nullhyp = nullhyp, reps = 2000)
      fit <- mrqap(yy, xx, mode = mode, reps = 20000, seed = 7,
                   method = if (nullhyp == "qap") "dsp" else "y")
      expect_equal(unname(fit$coefficients), unname(ref$coefficients), tolerance = 1e-10)
      expect_equal(unname(fit$t), unname(ref$tstat), tolerance = 1e-10)
      expect_close(unname(fit$p.value[, 1]), unname(ref$pleeq))
      expect_close(unname(fit$p.value[, 2]), unname(ref$pgreq))
      expect_close(unname(fit$p.value[, 3]), unname(ref$pgreqabs))
    }
  }
  # Beta as the test statistic.
  set.seed(12)
  ref <- sna::netlm(d$y, x, nullhyp = "qapy", test.statistic = "beta", reps = 2000)
  fit <- mrqap(d$y, x, method = "y", statistic = "beta", reps = 20000, seed = 8)
  expect_close(unname(fit$p.value[, 3]), unname(ref$pgreqabs))
})

test_that("qap_cor matches sna::gcor() and sna::qaptest()", {
  skip_if_not_installed("sna")
  d <- qap_data()
  for (mode in c("digraph", "graph")) {
    for (diag in c(FALSE, TRUE)) {
      res <- qap_cor(d$y, d$x1, mode = mode, diag = diag, reps = 0)
      expect_equal(res$statistic,
                   as.numeric(sna::gcor(d$y, d$x1, mode = mode, diag = diag)),
                   tolerance = 1e-12)
    }
  }
  y <- d$y
  y[c(4, 25)] <- NA
  expect_equal(qap_cor(y, d$x1, reps = 0)$statistic, as.numeric(sna::gcor(y, d$x1)),
               tolerance = 1e-12)

  set.seed(5)
  ref <- sna::qaptest(list(d$x2, d$x1), sna::gcor, g1 = 1, g2 = 2, reps = 2000)
  res <- qap_cor(d$x2, d$x1, reps = 20000, seed = 9)
  expect_equal(res$statistic, as.numeric(ref$testval), tolerance = 1e-12)
  expect_close(unname(res$p.value[["Pr(>=r)"]]), ref$pgreq)
  expect_close(unname(res$p.value[["Pr(<=r)"]]), ref$pleeq)
})

test_that("exact enumeration covers every permutation", {
  d <- qap_data(n = 5)
  res <- qap_cor(d$y, d$y, exact = TRUE)
  expect_equal(res$draws, 120)
  expect_equal(res$statistic, 1)
  expect_gte(res$counts[["ge"]], 1)
  expect_length(res$dist, 120)
  fit <- mrqap(d$y, list(d$x1, d$x2), exact = TRUE)
  expect_equal(fit$draws, 120)
  expect_equal(unname(fit$counts[, "valid"]), rep(120, 3))
  expect_error(qap_cor(qap_data(11)$y, qap_data(11)$x1, exact = TRUE), "at most 10")
})

test_that("set.seed() and seed make results reproducible", {
  d <- qap_data()
  set.seed(1)
  a <- mrqap(d$y, list(d$x1, d$x2), reps = 200)
  set.seed(1)
  b <- mrqap(d$y, list(d$x1, d$x2), reps = 200)
  expect_identical(a$dist, b$dist)
  expect_identical(a$seed, b$seed)
  # CRAN allows at most two cores during checks.
  c1 <- mrqap(d$y, list(d$x1, d$x2), reps = 200, seed = 99, threads = 1)
  c2 <- mrqap(d$y, list(d$x1, d$x2), reps = 200, seed = 99, threads = 2)
  expect_identical(c1$dist, c2$dist)
  expect_identical(qap_cor(d$y, d$x1, reps = 100, seed = 3, threads = 1)$dist,
                   qap_cor(d$y, d$x1, reps = 100, seed = 3, threads = 2)$dist)
  expect_error(qap_cor(d$y, d$x1, seed = -1), "seed")
})

test_that("netframes are accepted and aligned by node name", {
  e <- data.frame(from = c("a", "b", "c", "d", "a", "e"),
                  to = c("b", "c", "d", "e", "c", "a"), w = c(1, 2, 3, 1, 2, 5))
  g <- netframe(e, weight = "w")
  h <- netframe(e[c(2, 4, 6), ], nodes = data.frame(id = c("e", "d", "c", "b", "a")))
  m <- as.matrix(g, attr = "w")
  k <- as.matrix(h)[rownames(m), colnames(m)]
  expect_equal(qap_cor(g, h, reps = 0)$statistic, qap_cor(m, k, reps = 0)$statistic)
  expect_equal(qap_cor(g, h, reps = 0, attr = NULL)$statistic,
               stats::cor(dyads(m), dyads(k)))
  fit <- mrqap(g, list(tie = h), reps = 0)
  expect_equal(unname(fit$coefficients),
               unname(stats::coef(stats::lm(dyads(m) ~ dyads(k)))), tolerance = 1e-10)
})

test_that("one-column models use Y-permutation and printing works", {
  d <- qap_data()
  fit <- mrqap(d$y, d$x1, reps = 50, seed = 1, intercept = FALSE)
  expect_equal(fit$method, "y")
  expect_output(print(fit), "Multiple regression QAP")
  expect_output(print(summary(mrqap(d$y, list(d$x1, d$x2), reps = 50, seed = 1))),
                "Pr\\(>=\\|b\\|\\)")
  expect_output(print(summary(qap_cor(d$y, d$x1, reps = 50, seed = 1))), "Pr\\(>=r\\)")
  expect_error(mrqap(d$y, list(d$x1, d$x1), reps = 0), "rank deficient")
})
