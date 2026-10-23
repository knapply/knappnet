# Geospatial networks: edge lengths, proximity networks, distance decay,
# and group geometry. Hand-checked against the haversine
# formula and plain Euclidean distances.

cities <- function() {
  nodes <- data.frame(id = c("London", "Paris", "Berlin", "Madrid"),
                      lon = c(-0.1276, 2.3522, 13.4050, -3.7038),
                      lat = c(51.5072, 48.8566, 52.5200, 40.4168),
                      region = c("north", "north", "north", "south"),
                      stringsAsFactors = FALSE)
  edges <- data.frame(source = c("London", "Paris", "Paris", "Madrid"),
                      target = c("Paris", "Berlin", "Madrid", "Madrid"), stringsAsFactors = FALSE)
  net_spatial(net_frame(edges, nodes, directed = FALSE), coords = c("lon", "lat"),
                  radius = 6371.0088)
}

haversine_km <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180
  h <- sin((lat2 - lat1) * r / 2)^2 +
    cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371.0088 * asin(sqrt(h))
}

test_that("net_spatial() records coordinates and checks them", {
  g <- cities()
  expect_equal(c(g$roles$x, g$roles$y), c("lon", "lat"))
  expect_equal(g$roles$crs, "lonlat")
  expect_equal(g$roles$radius, 6371.0088)
  expect_s3_class(g, "net_frame")
  plain <- net_frame(g$edges, g$nodes, directed = FALSE)
  expect_error(net_edge_lengths(plain), "spatial_network")
  expect_error(net_spatial(plain, coords = c("lon", "height")), "coords\\[2\\]")
  expect_error(net_spatial(plain, coords = "lon"), "two columns")
  expect_error(net_spatial(plain, coords = c("id", "lat")), "numeric")
  bad <- plain
  bad$nodes$lat[2] <- 95
  expect_error(net_spatial(bad, coords = c("lon", "lat")), "outside")
  bad$nodes$lat[2] <- NA
  expect_error(net_spatial(bad, coords = c("lon", "lat")), "row 2")
  expect_error(net_spatial(plain, coords = c("lon", "lat"), radius = 0), "`radius`")
  # Planar coordinates have no latitude limit.
  expect_s3_class(net_spatial(bad, coords = c("lon", "lon"), crs = "planar"), "net_frame")
})

test_that("edge lengths are great-circle or Euclidean distances", {
  g <- cities()
  len <- net_edge_lengths(g)
  expect_equal(len[1], haversine_km(-0.1276, 51.5072, 2.3522, 48.8566))
  expect_equal(len[1], 343.5303, tolerance = 1e-6)
  expect_equal(len[2], 877.4645, tolerance = 1e-6)
  expect_equal(len[3], 1052.8937, tolerance = 1e-6)
  expect_equal(len[4], 0)  # a self-loop
  # Metres by default.
  gm <- net_spatial(net_frame(g$edges, g$nodes, directed = FALSE), coords = c("lon", "lat"))
  expect_equal(net_edge_lengths(gm)[1], 1000 * len[1] * 6371008.7714 / 6371008.8,
               tolerance = 1e-8)
  # Planar: a 3-4-5 triangle.
  p <- net_frame(data.frame(source = c(1, 2, 1), target = c(2, 3, 3)),
                nodes = data.frame(id = 1:3, x = c(0, 3, 0), y = c(0, 0, 4)), directed = FALSE)
  p <- net_spatial(p, coords = c("x", "y"), crs = "planar")
  expect_equal(net_edge_lengths(p), c(3, 5, 4))
  # Antipodal points are half the circumference apart.
  a <- net_spatial(net_frame(data.frame(source = 1, target = 2),
                                nodes = data.frame(id = 1:2, lon = c(0, 180), lat = c(0, 0)),
                                directed = FALSE), coords = c("lon", "lat"), radius = 1)
  expect_equal(net_edge_lengths(a), pi)
})

test_that("distance decay falls from 1 at distance 0", {
  d <- c(0, 1, 2, 10)
  expect_equal(net_distance_decay(d, scale = 2), exp(-d / 2))
  expect_equal(net_distance_decay(d, scale = 2, form = "gaussian"), exp(-(d / 2)^2))
  expect_equal(net_distance_decay(d, scale = 2, form = "power", exponent = 3), (1 + d / 2)^-3)
  expect_equal(net_distance_decay(2, scale = 2, form = "power", exponent = 1), 0.5)
  expect_error(net_distance_decay(d, scale = 0), "`scale`")
  expect_error(net_distance_decay(-1, scale = 1), "negative")
})

test_that("proximity networks join nearby points", {
  pts <- data.frame(id = letters[1:5], x = c(0, 1, 2, 10, 10.5), y = c(0, 0, 0, 0, 0))
  g <- net_proximity(pts, coords = c("x", "y"), crs = "planar", within = 1)
  expect_false(g$directed)
  expect_equal(g$edges$source, c("a", "b", "d"))
  expect_equal(g$edges$target, c("b", "c", "e"))
  expect_equal(g$edges$distance, c(1, 1, 0.5))
  expect_equal(g$roles$weight, "distance")
  expect_equal(c(g$roles$x, g$roles$y), c("x", "y"))
  expect_equal(net_edge_lengths(g), g$edges$distance)
  # Coincident points are joined; nothing else at radius 0.
  dup <- net_proximity(rbind(pts, data.frame(id = "f", x = 0, y = 0)),
                           coords = c("x", "y"), crs = "planar", within = 0)
  expect_equal(nrow(dup$edges), 1L)
  expect_equal(dup$edges$source, "a")
  expect_equal(dup$edges$target, "f")
  # Nearest neighbours: directed, nearest first, ties by position.
  k1 <- net_proximity(pts, coords = c("x", "y"), crs = "planar", k = 1)
  expect_true(k1$directed)
  expect_equal(k1$edges$source, letters[1:5])
  expect_equal(k1$edges$target, c("b", "a", "b", "e", "d"))
  k2 <- net_proximity(pts, coords = c("x", "y"), crs = "planar", k = 2)
  expect_equal(k2$edges$target[k2$edges$source == "b"], c("a", "c"))
  expect_equal(nrow(k2$edges), 10L)
  # "either" and "mutual" symmetrize.
  either <- net_proximity(pts, coords = c("x", "y"), crs = "planar", k = 1, knn = "either")
  expect_false(either$directed)
  expect_equal(nrow(either$edges), 3L)
  mutual <- net_proximity(pts, coords = c("x", "y"), crs = "planar", k = 1, knn = "mutual")
  expect_equal(nrow(mutual$edges), 2L)
  expect_equal(paste(mutual$edges$source, mutual$edges$target), c("a b", "d e"))
  # More neighbours than other points.
  k9 <- net_proximity(pts, coords = c("x", "y"), crs = "planar", k = 9)
  expect_equal(nrow(k9$edges), 20L)
  # A spatial net_frame can be given directly.
  g2 <- net_proximity(g, within = 1)
  expect_equal(g2$edges, g$edges)
  # Longitude and latitude: the same answer as the exact distances.
  set.seed(2)
  n <- 300
  world <- data.frame(id = seq_len(n), lon = runif(n, -180, 180),
                      lat = asin(runif(n, -1, 1)) * 180 / pi)
  w <- net_proximity(world, coords = c("lon", "lat"), within = 1.5e6, threads = 2)
  d <- outer(seq_len(n), seq_len(n), function(i, j) {
    haversine_km(world$lon[i], world$lat[i], world$lon[j], world$lat[j]) * 1000 *
      6371008.7714 / 6371008.8
  })
  pairs <- which(d <= 1.5e6 & row(d) < col(d), arr.ind = TRUE)
  expect_equal(nrow(w$edges), nrow(pairs))
  expect_equal(w$edges$source, pairs[order(pairs[, 1], pairs[, 2]), 1])
  expect_equal(w$edges$target, pairs[order(pairs[, 1], pairs[, 2]), 2])
  wk <- net_proximity(world, coords = c("lon", "lat"), k = 3)
  expect_equal(nrow(wk$edges), 3L * n)
  diag(d) <- Inf
  # apply() returns one column per node, so as.vector() lists each node's three in turn.
  nearest <- as.vector(apply(d, 1, function(r) order(r)[1:3]))
  expect_equal(wk$edges$target, nearest)
  expect_error(net_proximity(pts, coords = c("x", "y"), crs = "planar"), "exactly one")
  expect_error(net_proximity(pts, coords = c("x", "y"), crs = "planar", within = 1, k = 1),
               "exactly one")
  expect_error(net_proximity(pts, coords = c("x", "y"), crs = "planar", k = 0), "`k`")
  expect_error(net_proximity(pts, coords = c("x", "y"), crs = "planar", within = -1),
               "`within`")
  expect_error(net_proximity(1:3, coords = c("x", "y"), within = 1), "data frame")
})

test_that("group geometry gives centroids and dispersion", {
  p <- net_frame(data.frame(source = integer(0), target = integer(0)),
                nodes = data.frame(id = 1:5, x = c(0, 2, 2, 0, 10), y = c(0, 0, 2, 2, 10),
                                   grp = c("a", "a", "a", "a", "b")), directed = FALSE)
  p <- net_spatial(p, coords = c("x", "y"), crs = "planar")
  gg <- net_group_geometry(p, "grp")
  expect_equal(names(gg), c("group", "size", "x", "y", "mean_distance", "standard_distance"))
  expect_equal(gg$group, c("a", "b"))
  expect_equal(gg$size, c(4L, 1L))
  expect_equal(gg$x, c(1, 10))
  expect_equal(gg$y, c(1, 10))
  expect_equal(gg$mean_distance, c(sqrt(2), 0))
  expect_equal(gg$standard_distance, c(sqrt(2), 0))
  expect_equal(net_group_geometry(p, p$nodes$grp), gg)
  # Longitude and latitude: the midpoint of two points on the equator.
  g <- cities()
  g$nodes$lon <- c(0, 90, 0, 0)
  g$nodes$lat <- c(0, 0, 0, 0)
  eq <- net_group_geometry(g, c("x", "x", "y", "y"))
  expect_equal(eq$lon, c(45, 0))
  expect_equal(eq$lat, c(0, 0))
  expect_equal(eq$mean_distance[1], haversine_km(0, 0, 45, 0))
  expect_equal(eq$mean_distance[2], 0)
  # Antipodal points have no midpoint.
  g$nodes$lon <- c(0, 180, 0, 0)
  anti <- net_group_geometry(g, c("x", "x", "y", "y"))
  expect_true(is.nan(anti$lon[1]))
  expect_error(net_group_geometry(g, 1:2), "one value per node")
})
