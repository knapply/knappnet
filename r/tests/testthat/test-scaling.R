# Classical multidimensional scaling, stress, and the distance helpers,
# checked against stats::cmdscale() and hand values.

kite <- function() {
  e <- data.frame(
    from = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6, 7, 8, 9),
    to   = c(2, 3, 4, 6, 4, 5, 7, 4, 6, 5, 6, 7, 7, 8, 8, 9, 10)
  )
  netframe(e, nodes = data.frame(id = 1:10), directed = FALSE)
}

# Compares two configurations up to the sign of each dimension.
expect_same_configuration <- function(a, b, tolerance = 1e-8) {
  expect_equal(dim(a), dim(b))
  for (j in seq_len(ncol(a))) {
    x <- unname(a[, j])
    y <- unname(b[, j])
    expect_true(isTRUE(all.equal(x, y, tolerance = tolerance)) ||
                  isTRUE(all.equal(x, -y, tolerance = tolerance)))
  }
}

test_that("classical scaling recovers a Euclidean configuration", {
  set.seed(11)
  pts <- matrix(rnorm(40), 20, 2)
  d <- as.matrix(dist(pts))
  fit <- classical_mds(d, k = 2)
  expect_equal(dim(fit$points), c(20L, 2L))
  expect_equal(colnames(fit$points), c("Dim1", "Dim2"))
  expect_equal(as.matrix(dist(fit$points)), d, ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(fit$stress, 0, tolerance = 1e-10)
  expect_equal(fit$gof, c(1, 1), tolerance = 1e-10)
  expect_equal(length(fit$eigenvalues), 20L)
  expect_true(all(diff(fit$eigenvalues) <= 1e-12))
  expect_true(all(abs(fit$eigenvalues[-(1:2)]) < 1e-9))
  # A "dist" object is accepted and rownames are kept.
  rownames(pts) <- letters[1:20]
  fit2 <- classical_mds(dist(pts), k = 2)
  expect_equal(rownames(fit2$points), letters[1:20])
  expect_same_configuration(fit2$points, fit$points)
  # The first non-negligible coordinate of each dimension is positive.
  for (j in 1:2) {
    v <- fit$points[, j]
    expect_true(v[which(abs(v) > 1e-8 * max(abs(v)))[1]] > 0)
  }
})

test_that("classical scaling matches stats::cmdscale()", {
  g <- kite()
  d <- mds_distances(g)
  fit <- classical_mds(d, k = 2)
  cm <- stats::cmdscale(d, k = 2, eig = TRUE)
  expect_same_configuration(fit$points, cm$points)
  expect_equal(fit$eigenvalues, cm$eig, tolerance = 1e-10)
  expect_equal(fit$gof, cm$GOF, tolerance = 1e-10)
  # Non-Euclidean input: some eigenvalues are negative, only positive ones
  # get coordinates.
  se <- mds_distances(g, "euclidean")
  fit3 <- classical_mds(se, k = 3)
  cm3 <- suppressWarnings(stats::cmdscale(se, k = 3, eig = TRUE))
  expect_same_configuration(fit3$points, cm3$points[, seq_len(ncol(fit3$points)), drop = FALSE])
  expect_equal(fit3$eigenvalues, cm3$eig, tolerance = 1e-10)
  expect_equal(fit3$gof, cm3$GOF, tolerance = 1e-10)
  expect_equal(mds_stress(d, fit$points), fit$stress)
})

test_that("stress follows Kruskal's formula 1", {
  d <- matrix(c(0, 3, 4, 3, 0, 5, 4, 5, 0), 3)
  pts <- rbind(c(0, 0), c(3, 0), c(0, 4))
  expect_equal(mds_stress(d, pts), 0)
  pts2 <- rbind(c(0, 0), c(3, 0), c(0, 3))
  # Distances 3, 3, 3 sqrt(2) against 3, 4, 5.
  dd <- c(3, 3, 3 * sqrt(2))
  expect_equal(mds_stress(d, pts2), sqrt(sum((dd - c(3, 4, 5))^2) / sum(dd^2)))
  expect_error(mds_stress(d, pts[1:2, ]), "rows")
  expect_error(mds_stress(d, "x"), "numeric matrix")
})

test_that("invalid dissimilarities are refused", {
  d <- matrix(c(0, 1, 1, 0), 2)
  expect_error(classical_mds(d, k = 2), "at most 1")
  expect_error(classical_mds(d, k = 0), "`k`")
  expect_error(classical_mds(matrix(0, 1, 1), k = 1), "two points")
  expect_error(classical_mds(matrix(c(0, 1, 2, 0), 2), k = 1), "symmetric")
  expect_error(classical_mds(matrix(c(1, 1, 1, 0), 2), k = 1), "diagonal")
  expect_error(classical_mds(matrix(c(0, NA, NA, 0), 2), k = 1), "missing")
  expect_error(classical_mds(matrix(c(0, Inf, Inf, 0), 2), k = 1), "finite")
  expect_error(classical_mds(matrix(1:6, 2), k = 1), "square")
  expect_error(classical_mds("d", k = 1), "square numeric matrix")
})

test_that("mds_distances() gives geodesic and structural-equivalence distances", {
  g <- kite()
  d <- mds_distances(g)
  expect_equal(d, distances(g))
  expect_equal(dimnames(d), list(as.character(1:10), as.character(1:10)))
  # Unreachable pairs.
  u <- netframe(data.frame(from = c(1, 2, 4), to = c(2, 3, 5)), directed = FALSE)
  du <- mds_distances(u)
  expect_equal(max(du), 3)
  expect_equal(du["1", "4"], 3)
  expect_equal(mds_distances(u, unreachable = "diameter")["1", "4"], 2)
  expect_equal(mds_distances(u, unreachable = 10)["1", "4"], 10)
  expect_error(mds_distances(u, unreachable = "far"), "`unreachable`")
  # Direction is ignored.
  dd <- netframe(data.frame(from = c(1, 2), to = c(2, 3)))
  expect_equal(unname(mds_distances(dd)), matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), 3))
  expect_equal(mds_distances(g, "euclidean"), structural_equivalence(g, "euclidean"))
  expect_equal(mds_distances(g, "hamming"), structural_equivalence(g, "hamming"))
  cr <- mds_distances(g, "correlation")
  expect_equal(unname(diag(cr)), rep(0, 10))
  expect_equal(cr[1, 2], 1 - structural_equivalence(g, "correlation")[1, 2])
  expect_equal(classical_mds(cr, k = 2)$stress, mds_stress(cr, classical_mds(cr, k = 2)$points))
})
