# Geospatial networks: edge lengths, proximity networks, distance decay,
# and group geometry. Hand-checked against the haversine
# formula and plain Euclidean distances.

cities <- function() {
  nodes <- data.frame(name = c("London", "Paris", "Berlin", "Madrid"),
                      lon = c(-0.1276, 2.3522, 13.4050, -3.7038),
                      lat = c(51.5072, 48.8566, 52.5200, 40.4168),
                      region = c("north", "north", "north", "south"),
                      stringsAsFactors = FALSE)
  edges <- data.frame(from = c("London", "Paris", "Paris", "Madrid"),
                      to = c("Paris", "Berlin", "Madrid", "Madrid"), stringsAsFactors = FALSE)
  spatial_network(netframe(edges, nodes, directed = FALSE), coords = c("lon", "lat"),
                  radius = 6371.0088)
}

haversine_km <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180
  h <- sin((lat2 - lat1) * r / 2)^2 +
    cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371.0088 * asin(sqrt(h))
}

test_that("spatial_network() records coordinates and checks them", {
  g <- cities()
  expect_equal(c(g$roles$x, g$roles$y), c("lon", "lat"))
  expect_equal(g$roles$crs, "lonlat")
  expect_equal(g$roles$radius, 6371.0088)
  expect_s3_class(g, "netframe")
  plain <- netframe(g$edges, g$nodes, directed = FALSE)
  expect_error(edge_lengths(plain), "spatial_network")
  expect_error(spatial_network(plain, coords = c("lon", "height")), "coords\\[2\\]")
  expect_error(spatial_network(plain, coords = "lon"), "two columns")
  expect_error(spatial_network(plain, coords = c("name", "lat")), "numeric")
  bad <- plain
  bad$nodes$lat[2] <- 95
  expect_error(spatial_network(bad, coords = c("lon", "lat")), "outside")
  bad$nodes$lat[2] <- NA
  expect_error(spatial_network(bad, coords = c("lon", "lat")), "row 2")
  expect_error(spatial_network(plain, coords = c("lon", "lat"), radius = 0), "`radius`")
  # Planar coordinates have no latitude limit.
  expect_s3_class(spatial_network(bad, coords = c("lon", "lon"), crs = "planar"), "netframe")
})

test_that("edge lengths are great-circle or Euclidean distances", {
  g <- cities()
  len <- edge_lengths(g)
  expect_equal(len[1], haversine_km(-0.1276, 51.5072, 2.3522, 48.8566))
  expect_equal(len[1], 343.5303, tolerance = 1e-6)
  expect_equal(len[2], 877.4645, tolerance = 1e-6)
  expect_equal(len[3], 1052.8937, tolerance = 1e-6)
  expect_equal(len[4], 0)  # a self-loop
  # Metres by default.
  gm <- spatial_network(netframe(g$edges, g$nodes, directed = FALSE), coords = c("lon", "lat"))
  expect_equal(edge_lengths(gm)[1], 1000 * len[1] * 6371008.7714 / 6371008.8,
               tolerance = 1e-8)
  # Planar: a 3-4-5 triangle.
  p <- netframe(data.frame(from = c(1, 2, 1), to = c(2, 3, 3)),
                nodes = data.frame(id = 1:3, x = c(0, 3, 0), y = c(0, 0, 4)), directed = FALSE)
  p <- spatial_network(p, coords = c("x", "y"), crs = "planar")
  expect_equal(edge_lengths(p), c(3, 5, 4))
  # Antipodal points are half the circumference apart.
  a <- spatial_network(netframe(data.frame(from = 1, to = 2),
                                nodes = data.frame(id = 1:2, lon = c(0, 180), lat = c(0, 0)),
                                directed = FALSE), coords = c("lon", "lat"), radius = 1)
  expect_equal(edge_lengths(a), pi)
})

test_that("distance decay falls from 1 at distance 0", {
  d <- c(0, 1, 2, 10)
  expect_equal(distance_decay(d, scale = 2), exp(-d / 2))
  expect_equal(distance_decay(d, scale = 2, form = "gaussian"), exp(-(d / 2)^2))
  expect_equal(distance_decay(d, scale = 2, form = "power", exponent = 3), (1 + d / 2)^-3)
  expect_equal(distance_decay(2, scale = 2, form = "power", exponent = 1), 0.5)
  expect_error(distance_decay(d, scale = 0), "`scale`")
  expect_error(distance_decay(-1, scale = 1), "negative")
})

test_that("proximity networks join nearby points", {
  pts <- data.frame(id = letters[1:5], x = c(0, 1, 2, 10, 10.5), y = c(0, 0, 0, 0, 0))
  g <- proximity_network(pts, coords = c("x", "y"), crs = "planar", within = 1)
  expect_false(g$directed)
  expect_equal(g$edges$from, c("a", "b", "d"))
  expect_equal(g$edges$to, c("b", "c", "e"))
  expect_equal(g$edges$distance, c(1, 1, 0.5))
  expect_equal(g$roles$weight, "distance")
  expect_equal(c(g$roles$x, g$roles$y), c("x", "y"))
  expect_equal(edge_lengths(g), g$edges$distance)
  # Coincident points are joined; nothing else at radius 0.
  dup <- proximity_network(rbind(pts, data.frame(id = "f", x = 0, y = 0)),
                           coords = c("x", "y"), crs = "planar", within = 0)
  expect_equal(nrow(dup$edges), 1L)
  expect_equal(dup$edges$from, "a")
  expect_equal(dup$edges$to, "f")
  # Nearest neighbours: directed, nearest first, ties by position.
  k1 <- proximity_network(pts, coords = c("x", "y"), crs = "planar", k = 1)
  expect_true(k1$directed)
  expect_equal(k1$edges$from, letters[1:5])
  expect_equal(k1$edges$to, c("b", "a", "b", "e", "d"))
  k2 <- proximity_network(pts, coords = c("x", "y"), crs = "planar", k = 2)
  expect_equal(k2$edges$to[k2$edges$from == "b"], c("a", "c"))
  expect_equal(nrow(k2$edges), 10L)
  # "either" and "mutual" symmetrize.
  either <- proximity_network(pts, coords = c("x", "y"), crs = "planar", k = 1, knn = "either")
  expect_false(either$directed)
  expect_equal(nrow(either$edges), 3L)
  mutual <- proximity_network(pts, coords = c("x", "y"), crs = "planar", k = 1, knn = "mutual")
  expect_equal(nrow(mutual$edges), 2L)
  expect_equal(paste(mutual$edges$from, mutual$edges$to), c("a b", "d e"))
  # More neighbours than other points.
  k9 <- proximity_network(pts, coords = c("x", "y"), crs = "planar", k = 9)
  expect_equal(nrow(k9$edges), 20L)
  # A spatial netframe can be given directly.
  g2 <- proximity_network(g, within = 1)
  expect_equal(g2$edges, g$edges)
  # Longitude and latitude: the same answer as the exact distances.
  set.seed(2)
  n <- 300
  world <- data.frame(id = seq_len(n), lon = runif(n, -180, 180),
                      lat = asin(runif(n, -1, 1)) * 180 / pi)
  w <- proximity_network(world, coords = c("lon", "lat"), within = 1.5e6, threads = 2)
  d <- outer(seq_len(n), seq_len(n), function(i, j) {
    haversine_km(world$lon[i], world$lat[i], world$lon[j], world$lat[j]) * 1000 *
      6371008.7714 / 6371008.8
  })
  pairs <- which(d <= 1.5e6 & row(d) < col(d), arr.ind = TRUE)
  expect_equal(nrow(w$edges), nrow(pairs))
  expect_equal(w$edges$from, pairs[order(pairs[, 1], pairs[, 2]), 1])
  expect_equal(w$edges$to, pairs[order(pairs[, 1], pairs[, 2]), 2])
  wk <- proximity_network(world, coords = c("lon", "lat"), k = 3)
  expect_equal(nrow(wk$edges), 3L * n)
  diag(d) <- Inf
  nearest <- as.vector(t(apply(d, 1, function(r) order(r)[1:3])))
  expect_equal(wk$edges$to, nearest)
  expect_error(proximity_network(pts, coords = c("x", "y"), crs = "planar"), "exactly one")
  expect_error(proximity_network(pts, coords = c("x", "y"), crs = "planar", within = 1, k = 1),
               "exactly one")
  expect_error(proximity_network(pts, coords = c("x", "y"), crs = "planar", k = 0), "`k`")
  expect_error(proximity_network(pts, coords = c("x", "y"), crs = "planar", within = -1),
               "`within`")
  expect_error(proximity_network(1:3, coords = c("x", "y"), within = 1), "data frame")
})

test_that("group geometry gives centroids and dispersion", {
  p <- netframe(data.frame(from = integer(0), to = integer(0)),
                nodes = data.frame(id = 1:5, x = c(0, 2, 2, 0, 10), y = c(0, 0, 2, 2, 10),
                                   grp = c("a", "a", "a", "a", "b")), directed = FALSE)
  p <- spatial_network(p, coords = c("x", "y"), crs = "planar")
  gg <- group_geometry(p, "grp")
  expect_equal(names(gg), c("group", "size", "x", "y", "mean_distance", "standard_distance"))
  expect_equal(gg$group, c("a", "b"))
  expect_equal(gg$size, c(4L, 1L))
  expect_equal(gg$x, c(1, 10))
  expect_equal(gg$y, c(1, 10))
  expect_equal(gg$mean_distance, c(sqrt(2), 0))
  expect_equal(gg$standard_distance, c(sqrt(2), 0))
  expect_equal(group_geometry(p, p$nodes$grp), gg)
  # Longitude and latitude: the midpoint of two points on the equator.
  g <- cities()
  g$nodes$lon <- c(0, 90, 0, 0)
  g$nodes$lat <- c(0, 0, 0, 0)
  eq <- group_geometry(g, c("x", "x", "y", "y"))
  expect_equal(eq$lon, c(45, 0))
  expect_equal(eq$lat, c(0, 0))
  expect_equal(eq$mean_distance[1], haversine_km(0, 0, 45, 0))
  expect_equal(eq$mean_distance[2], 0)
  # Antipodal points have no midpoint.
  g$nodes$lon <- c(0, 180, 0, 0)
  anti <- group_geometry(g, c("x", "x", "y", "y"))
  expect_true(is.nan(anti$lon[1]))
  expect_error(group_geometry(g, 1:2), "one value per node")
})
