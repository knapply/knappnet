# Spatial networks built on sf: geometry roles, representative points,
# geometry-derived edge lengths (checked against sf::st_length and
# sf::st_distance directly), builders from lines and polygons (checked
# against spdep::poly2nb and sfnetworks where installed), spatial
# operations that keep the network consistent, geographic plots,
# GeoPackage files, and the sfnetworks round trip. Everything here needs
# sf; the package itself must work without it.

# Four points in UTM zone 33N (metres): a unit right triangle plus a point
# 10 m away.
utm_points <- function() {
  sf::st_as_sf(data.frame(name = c("a", "b", "c", "d"), x = c(0, 3, 0, 10), y = c(0, 0, 4, 0),
                          kind = c("u", "u", "v", "v"), stringsAsFactors = FALSE),
               coords = c("x", "y"), crs = 32633, remove = FALSE)
}

utm_network <- function(directed = FALSE) {
  edges <- data.frame(from = c("a", "b", "c", "a"), to = c("b", "c", "a", "d"),
                      stringsAsFactors = FALSE)
  netframe(edges, utm_points(), directed = directed)
}

city_points <- function() {
  sf::st_as_sf(data.frame(name = c("London", "Paris", "Berlin", "Madrid"),
                          lon = c(-0.1276, 2.3522, 13.4050, -3.7038),
                          lat = c(51.5072, 48.8566, 52.5200, 40.4168),
                          stringsAsFactors = FALSE),
               coords = c("lon", "lat"), crs = 4326, remove = FALSE)
}

square <- function(x, y, side = 1) {
  sf::st_polygon(list(rbind(c(x, y), c(x + side, y), c(x + side, y + side), c(x, y + side),
                            c(x, y))))
}

# A 2-by-2 grid of unit squares: w (0,0), x (1,0), y (0,1), z (1,1).
grid_regions <- function(crs = 32633) {
  sf::st_sf(region = c("w", "x", "y", "z"), pop = c(10, 20, 30, 40),
            geometry = sf::st_sfc(square(0, 0), square(1, 0), square(0, 1), square(1, 1),
                                  crs = crs), stringsAsFactors = FALSE)
}

segment <- function(x0, y0, x1, y1) sf::st_linestring(rbind(c(x0, y0), c(x1, y1)))

# Every unit segment of a 3-by-3 lattice: 9 junctions, 12 streets.
lattice_lines <- function(crs = 32633) {
  segs <- list()
  for (i in 0:2) {
    for (j in 0:1) {
      segs <- c(segs, list(segment(j, i, j + 1, i)), list(segment(i, j, i, j + 1)))
    }
  }
  sf::st_sf(street = sprintf("s%02d", seq_along(segs)), geometry = sf::st_sfc(segs, crs = crs),
            stringsAsFactors = FALSE)
}

test_that("sf node tables record geometry roles and keep every column", {
  skip_if_not_installed("sf")
  g <- utm_network()
  expect_s3_class(g$nodes, "sf")
  expect_equal(g$roles$node_geometry, "geometry")
  expect_null(g$roles$edge_geometry)
  expect_s3_class(g$roles$crs, "crs")
  expect_equal(g$roles$crs$epsg, 32633L)
  expect_null(g$roles$x)
  expect_equal(names(g$nodes), c("name", "x", "y", "kind", "geometry"))
  expect_output(print(g), "spatial: nodes POINT; CRS")
  # Measures run as usual; the engine never sees the geometry.
  expect_equal(unname(degree(g)), c(3L, 2L, 2L, 1L))
  # as_netframe() on sf works both ways round.
  expect_equal(as_netframe(utm_points(), edges = g$edges, directed = FALSE)$edges, g$edges)
  edge_sf <- to_sf(g, "edges")
  h <- as_netframe(edge_sf, nodes = utm_points(), directed = FALSE)
  expect_equal(h$roles$edge_geometry, "geometry")
  expect_equal(h$roles$node_geometry, "geometry")
  # Nodes derived from an sf edge table have no geometry of their own.
  k <- netframe(edge_sf, directed = FALSE)
  expect_null(k$roles$node_geometry)
  expect_equal(k$roles$edge_geometry, "geometry")
  expect_error(to_sf(k, "nodes"), "no geometry")
})

test_that("nodes and edges must share a CRS", {
  skip_if_not_installed("sf")
  g <- utm_network()
  edges <- sf::st_transform(to_sf(g, "edges"), 4326)
  expect_error(netframe(edges, utm_points(), directed = FALSE), "different coordinate reference")
  expect_error(netframe(edges, utm_points(), directed = FALSE), "st_transform")
  no_crs <- sf::st_set_crs(to_sf(g, "edges"), NA)
  expect_error(netframe(no_crs, utm_points(), directed = FALSE), "st_set_crs")
  # The same CRS, or none on both, is fine.
  expect_s3_class(netframe(to_sf(g, "edges"), utm_points(), directed = FALSE), "netframe")
  expect_s3_class(netframe(no_crs, sf::st_set_crs(utm_points(), NA), directed = FALSE),
                  "netframe")
})

test_that("representative points stand in for any node geometry", {
  skip_if_not_installed("sf")
  regions <- grid_regions()
  g <- from_polygons(regions, id = "region")
  p <- knappnet:::node_coordinates(g)
  expect_false(p$lonlat)
  expect_equal(colnames(p$coords), c("x", "y"))
  expect_equal(unname(p$coords), cbind(c(0.5, 1.5, 0.5, 1.5), c(0.5, 0.5, 1.5, 1.5)))
  # Points on surface lie on the geometry, centroids need not: the centroid
  # of a U shape, (1.5, 9.5 / 7), falls in the gap between its arms.
  u_shape <- sf::st_polygon(list(rbind(c(0, 0), c(3, 0), c(3, 3), c(2, 3), c(2, 1), c(1, 1),
                                       c(1, 3), c(0, 3), c(0, 0))))
  shapes <- sf::st_sf(id = c("c", "s"), geometry = sf::st_sfc(u_shape, square(10, 10)),
                      crs = 32633)
  h <- netframe(data.frame(from = "c", to = "s"), shapes, directed = FALSE)
  pos <- knappnet:::node_points(h)
  expect_true(all(sf::st_within(pos, sf::st_geometry(shapes), sparse = FALSE)[cbind(1:2, 1:2)]))
  cen <- knappnet:::node_points(spatial_network(h, representative = "centroid"))
  expect_equal(unname(sf::st_coordinates(cen)[2, ]), c(10.5, 10.5))
  expect_equal(unname(sf::st_coordinates(cen)[1, ]), c(1.5, 9.5 / 7))
  expect_false(sf::st_within(cen[1], sf::st_geometry(shapes)[1], sparse = FALSE)[1, 1])
  expect_equal(spatial_network(h, radius = 1)$roles$radius, 1)
  expect_error(spatial_network(h, coords = c("x", "y")), "already carries geometry")
  # Point-based functions work on polygons through the representative points.
  gg <- group_geometry(g, c("low", "low", "high", "high"))
  expect_equal(gg$x, c(1, 1))
  expect_equal(gg$y, c(0.5, 1.5))
  # Multipolygons and empty geometries.
  two <- sf::st_multipolygon(list(square(0, 0), square(5, 0)))
  mp <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(two, square(0, 5)), crs = 32633)
  m <- netframe(data.frame(from = 1, to = 2), mp, directed = FALSE)
  expect_true(sf::st_within(knappnet:::node_points(m)[1], sf::st_geometry(mp)[1],
                            sparse = FALSE)[1, 1])
  empty <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(sf::st_polygon(), square(0, 5)),
                     crs = 32633)
  e <- netframe(data.frame(from = 1, to = 2), empty, directed = FALSE)
  expect_error(edge_lengths(e), "row 1")
})

test_that("edge lengths come from sf when edges have geometry", {
  skip_if_not_installed("sf")
  g <- utm_network()
  # No edge geometry: Euclidean distances between the points, in metres.
  len <- edge_lengths(g)
  expect_equal(unname(len), c(3, 5, 4, 10))
  expect_null(attr(len, "unit"))
  # Explicit straight edges agree with sf::st_length and record the unit.
  h <- netframe(to_sf(g, "edges"), utm_points(), directed = FALSE)
  hl <- edge_lengths(h)
  expect_equal(attr(hl, "unit"), "m")
  expect_equal(as.vector(hl), as.numeric(sf::st_length(sf::st_geometry(h$edges))))
  expect_equal(as.vector(hl), c(3, 5, 4, 10))
  # A bent edge is longer than the straight one.
  edges <- to_sf(g, "edges")
  bent <- sf::st_linestring(rbind(c(0, 0), c(0, 4), c(3, 4), c(3, 0)))
  sf::st_geometry(edges)[1] <- sf::st_sfc(bent, crs = 32633)
  b <- netframe(edges, utm_points(), directed = FALSE)
  expect_equal(as.vector(edge_lengths(b)), c(11, 5, 4, 10))
  # Mixed: edges without geometry get the distance between representative points.
  sf::st_geometry(edges)[2] <- sf::st_sfc(sf::st_linestring(), crs = 32633)
  mixed <- netframe(edges, utm_points(), directed = FALSE)
  expect_output(print(mixed), "edges LINESTRING, 1 without geometry")
  ml <- edge_lengths(mixed)
  expect_equal(as.vector(ml), c(11, 5, 4, 10))
  pts <- sf::st_geometry(utm_points())
  expect_equal(ml[2], as.numeric(sf::st_distance(pts[2], pts[3])), ignore_attr = TRUE)
  # Longitude and latitude: geodesic lengths from sf, and the sphere
  # otherwise, which differ by well under one percent.
  cities <- city_points()
  ce <- data.frame(from = c("London", "Paris", "Paris"), to = c("Paris", "Berlin", "Madrid"),
                   stringsAsFactors = FALSE)
  cg <- netframe(ce, cities, directed = FALSE)
  expect_true(knappnet:::node_coordinates(cg)$lonlat)
  sphere <- edge_lengths(cg)
  expect_equal(unname(sphere[1]), 343530.3, tolerance = 1e-6)
  cl <- edge_lengths(netframe(to_sf(cg, "edges"), cities, directed = FALSE))
  expect_equal(attr(cl, "unit"), "m")
  expect_equal(as.vector(cl), as.numeric(sf::st_length(sf::st_geometry(to_sf(cg, "edges")))))
  expect_true(all(abs(as.vector(cl) / sphere - 1) < 0.006))
  # Projected lengths in other units carry their unit.
  ft <- transform_network(h, "+proj=utm +zone=33 +datum=WGS84 +units=us-ft")
  expect_match(attr(edge_lengths(ft), "unit"), "foot|ft")
  expect_equal(as.vector(edge_lengths(ft)), c(3, 5, 4, 10) * 3937 / 1200, tolerance = 1e-6)
})

test_that("to_sf() returns sf tables, with implicit edges when asked", {
  skip_if_not_installed("sf")
  g <- utm_network()
  nodes <- to_sf(g, "nodes")
  expect_s3_class(nodes, "sf")
  expect_equal(nodes$name, c("a", "b", "c", "d"))
  edges <- to_sf(g, "edges")
  expect_s3_class(edges, "sf")
  expect_equal(names(edges), c("from", "to", "geometry"))
  expect_equal(as.character(sf::st_geometry_type(edges, by_geometry = TRUE)),
               rep("LINESTRING", 4))
  expect_equal(sf::st_crs(edges), sf::st_crs(32633))
  expect_equal(unname(sf::st_coordinates(edges)[1:2, 1:2]), rbind(c(0, 0), c(3, 0)))
  expect_error(to_sf(g, "edges", implicit = FALSE), "no geometry")
  expect_equal(sf::st_as_sf(g), nodes)
  expect_equal(sf::st_as_sf(g, "edges"), edges)
  # Coordinate columns become points too: WGS 84 for lon/lat, no CRS for planar.
  plain <- spatial_network(netframe(g$edges, sf::st_drop_geometry(g$nodes), directed = FALSE),
                           coords = c("x", "y"), crs = "planar")
  expect_output(print(plain), "spatial: coordinates in `x`, `y` \\(planar\\)")
  ps <- to_sf(plain, "nodes")
  expect_s3_class(ps, "sf")
  expect_true(is.na(sf::st_crs(ps)))
  expect_equal(ps$x, c(0, 3, 0, 10))
  cities <- spatial_network(netframe(data.frame(from = "London", to = "Paris"),
                                     sf::st_drop_geometry(city_points()), directed = FALSE))
  expect_equal(sf::st_crs(to_sf(cities, "nodes")), sf::st_crs(4326))
  expect_equal(sf::st_crs(to_sf(cities, "edges")), sf::st_crs(4326))
  # A netframe without any location cannot become sf.
  expect_error(to_sf(netframe(g$edges, directed = FALSE)), "no geometry")
})

test_that("check_geometry() reports edges that leave their nodes", {
  skip_if_not_installed("sf")
  g <- utm_network()
  edges <- to_sf(g, "edges")
  expect_true(all(check_geometry(netframe(edges, utm_points(), directed = FALSE))))
  expect_equal(check_geometry(g), rep(NA, 4))
  # Reversed lines pass undirected, fail directed; an offset end fails
  # unless the tolerance covers it; an empty geometry is NA.
  sf::st_geometry(edges)[1] <- sf::st_reverse(sf::st_geometry(edges)[1])
  sf::st_geometry(edges)[2] <- sf::st_sfc(segment(3, 0, 0.5, 4), crs = 32633)
  sf::st_geometry(edges)[3] <- sf::st_sfc(sf::st_linestring(), crs = 32633)
  und <- check_geometry(netframe(edges, utm_points(), directed = FALSE))
  expect_equal(und, c(TRUE, FALSE, NA, TRUE))
  dir <- check_geometry(netframe(edges, utm_points(), directed = TRUE))
  expect_equal(dir, c(FALSE, FALSE, NA, TRUE))
  expect_equal(check_geometry(netframe(edges, utm_points(), directed = TRUE), tolerance = 0.5),
               c(FALSE, TRUE, NA, TRUE))
  # A line ending inside a polygon node touches it.
  regions <- grid_regions()
  rg <- from_polygons(regions, id = "region")
  rl <- to_sf(rg, "edges")
  expect_true(all(check_geometry(netframe(rl, regions, directed = FALSE, id = "region"))))
  expect_error(check_geometry(g, tolerance = -1), "tolerance")
})

test_that("from_lines() snaps endpoints into shared nodes", {
  skip_if_not_installed("sf")
  lines <- lattice_lines()
  g <- from_lines(lines)
  expect_false(g$directed)
  expect_equal(nrow(g$nodes), 9L)
  expect_equal(nrow(g$edges), 12L)
  expect_s3_class(g$nodes, "sf")
  expect_s3_class(g$edges, "sf")
  expect_equal(g$nodes$id, 1:9)
  expect_equal(names(g$edges), c("from", "to", "street", "length", "geometry"))
  expect_equal(g$roles$weight, "length")
  expect_equal(g$roles$weight_type, "distance")
  expect_equal(g$edges$length, rep(1, 12))
  expect_equal(as.vector(edge_lengths(g)), rep(1, 12))
  expect_equal(sort(unname(degree(g))), c(2, 2, 2, 2, 3, 3, 3, 3, 4))
  expect_true(all(check_geometry(g)))
  # Nodes sit at the first endpoint that formed them, numbered by appearance.
  expect_equal(unname(sf::st_coordinates(g$nodes)[1:2, ]), rbind(c(0, 0), c(1, 0)))
  expect_equal(sf::st_crs(g$nodes), sf::st_crs(32633))
  # Shortest paths run along the streets.
  expect_equal(unname(distances(g)[1, 9]), 4)
  # Directed: from first vertex to last.
  d <- from_lines(lines, directed = TRUE)
  expect_true(d$directed)
  expect_equal(d$edges$from[1], 1L)
  expect_equal(d$edges$to[1], 2L)
  # Reserved column names are kept as attributes; a ring is a self-loop.
  ring <- sf::st_sf(from = "x", length = 9,
                    geometry = sf::st_sfc(sf::st_linestring(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 0))),
                                          crs = 32633))
  r <- from_lines(ring)
  expect_equal(nrow(r$nodes), 1L)
  expect_equal(names(r$edges), c("from", "to", "from_attr", "length_attr", "length", "geometry"))
  expect_equal(r$edges$length, 2 + sqrt(2))
  expect_output(print(r), "1 self-loop")
  # Tolerance: nearly touching ends become one node, and not without it.
  gap <- sf::st_sf(geometry = sf::st_sfc(segment(0, 0, 1, 0), segment(1.05, 0, 2, 0),
                                         segment(2, 0, 3, 0), crs = 32633))
  expect_equal(nrow(from_lines(gap)$nodes), 5L)
  s <- from_lines(gap, tolerance = 0.1)
  expect_equal(nrow(s$nodes), 4L)
  expect_equal(s$edges$from, c(1L, 2L, 3L))
  expect_equal(s$edges$to, c(2L, 3L, 4L))
  expect_true(all(check_geometry(s, tolerance = 0.1)))
  expect_false(all(check_geometry(s)))
  # An sfc works; other geometries do not.
  expect_equal(nrow(from_lines(sf::st_geometry(lines))$edges), 12L)
  expect_error(from_lines(sf::st_cast(lines, "MULTILINESTRING")), "LINESTRING")
  expect_error(from_lines(grid_regions()), "row")
  expect_error(from_lines(data.frame(a = 1)), "sf")
  # Longitude and latitude: lengths are geodesic metres.
  ll <- from_lines(sf::st_sf(geometry = sf::st_sfc(segment(0, 0, 1, 0), crs = 4326)))
  expect_equal(ll$edges$length, as.numeric(sf::st_length(sf::st_geometry(ll$edges))))
  expect_true(abs(ll$edges$length - 111195) < 200)
})

test_that("from_polygons() gives queen and rook contiguity", {
  skip_if_not_installed("sf")
  regions <- grid_regions()
  queen <- from_polygons(regions, id = "region")
  expect_false(queen$directed)
  expect_s3_class(queen$nodes, "sf")
  expect_false(inherits(queen$edges, "sf"))
  expect_equal(queen$nodes$region, c("w", "x", "y", "z"))
  expect_equal(queen$nodes$pop, c(10, 20, 30, 40))
  expect_equal(paste(queen$edges$from, queen$edges$to),
               c("w x", "w y", "w z", "x y", "x z", "y z"))
  rook <- from_polygons(regions, contiguity = "rook", id = "region")
  expect_equal(paste(rook$edges$from, rook$edges$to), c("w x", "w y", "x z", "y z"))
  # Default identifiers number the regions; an existing `id` is kept aside.
  auto <- from_polygons(regions)
  expect_equal(auto$nodes$id, 1:4)
  expect_equal(names(auto$nodes), c("id", "region", "pop", "geometry"))
  withid <- regions
  withid$id <- c("p", "q", "r", "s")
  expect_equal(names(from_polygons(withid)$nodes), c("id", "region", "pop", "id_attr", "geometry"))
  # Isolates stay; overlapping polygons are not neighbours; multipolygons work.
  far <- rbind(regions, sf::st_sf(region = "far", pop = 1,
                                  geometry = sf::st_sfc(square(10, 10), crs = 32633)))
  expect_equal(nrow(from_polygons(far, id = "region")$edges), 6L)
  expect_equal(nrow(from_polygons(far, id = "region")$nodes), 5L)
  overlap <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(square(0, 0), square(0.5, 0), crs = 32633))
  expect_equal(nrow(from_polygons(overlap, id = "id")$edges), 0L)
  two <- sf::st_multipolygon(list(square(0, 0), square(5, 5)))
  mp <- sf::st_sf(id = 1:3, geometry = sf::st_sfc(two, square(1, 0), square(6, 5)), crs = 32633)
  expect_equal(nrow(from_polygons(mp, id = "id")$edges), 2L)
  # Edges have no geometry, so lengths run between representative points.
  expect_equal(unname(edge_lengths(rook)), rep(1, 4))
  expect_error(from_polygons(lattice_lines()), "POLYGON")
  # Longitude and latitude are evaluated as planar, quietly.
  expect_silent(ll <- from_polygons(grid_regions(crs = 4326), id = "region"))
  expect_equal(nrow(ll$edges), 6L)
})

test_that("contiguity agrees with spdep::poly2nb", {
  skip_if_not_installed("sf")
  skip_if_not_installed("spdep")
  set.seed(4)
  # An irregular map with exact coordinates: a 6-by-6 grid of unit squares
  # with a third of the cells removed, so regions meet along edges, at
  # corners only, or not at all. Shared boundaries have shared vertices, the
  # case in which spdep's vertex rule and the DE-9IM rule agree.
  cells <- expand.grid(x = 0:5, y = 0:5)
  cells <- cells[sort(sample(nrow(cells), 24)), ]
  geom <- sf::st_sfc(lapply(seq_len(nrow(cells)), function(i) square(cells$x[i], cells$y[i])),
                     crs = 32633)
  tiles <- sf::st_sf(id = seq_along(geom), geometry = geom)
  for (contiguity in c("queen", "rook")) {
    g <- from_polygons(tiles, contiguity = contiguity, id = "id")
    # spdep warns that the map has more than one connected piece, as intended.
    nb <- suppressWarnings(spdep::poly2nb(tiles, queen = contiguity == "queen"))
    ours <- lapply(seq_len(nrow(tiles)), function(i) {
      sort(c(g$edges$to[g$edges$from == i], g$edges$from[g$edges$to == i]))
    })
    theirs <- lapply(nb, function(v) sort(as.integer(v[v > 0])))
    expect_equal(ours, theirs, info = contiguity)
  }
})

test_that("proximity_network() takes sf nodes of any geometry", {
  skip_if_not_installed("sf")
  pts <- utm_points()
  g <- proximity_network(pts, within = 5)
  expect_s3_class(g$nodes, "sf")
  expect_equal(g$roles$node_geometry, "geometry")
  expect_equal(g$roles$crs$epsg, 32633L)
  expect_equal(g$edges$from, c("a", "a", "b"))
  expect_equal(g$edges$to, c("b", "c", "c"))
  expect_equal(g$edges$distance, c(3, 4, 5))
  # The same from a netframe with sf nodes, and with exact distances.
  expect_equal(proximity_network(utm_network(), within = 5)$edges, g$edges)
  expect_equal(proximity_network(pts, within = 5, exact = TRUE)$edges, g$edges)
  k <- proximity_network(pts, k = 1)
  expect_true(k$directed)
  expect_equal(k$edges$to, c("b", "a", "a", "b"))
  expect_equal(proximity_network(pts, k = 1, exact = TRUE)$edges, k$edges)
  expect_equal(proximity_network(pts, k = 1, knn = "mutual")$edges$from, "a")
  # Without an identifier column the features are numbered.
  anon <- sf::st_sf(geometry = sf::st_geometry(pts))
  expect_equal(proximity_network(anon, within = 5)$edges$from, c(1L, 1L, 2L))
  expect_equal(proximity_network(anon, within = 5, id = NULL)$nodes$id, 1:4)
  # Longitude and latitude: metres on the sphere, matching edge_lengths().
  cities <- city_points()
  w <- proximity_network(cities, within = 9e5)
  expect_equal(w$edges$from, c("London", "Paris"))
  expect_equal(w$edges$to, c("Paris", "Berlin"))
  expect_equal(w$edges$distance, unname(edge_lengths(w)))
  wx <- proximity_network(cities, within = 9e5, exact = TRUE)
  expect_equal(wx$edges$to, w$edges$to)
  expect_true(all(abs(wx$edges$distance / w$edges$distance - 1) < 0.006))
  # Polygons: representative points versus the geometries themselves.
  regions <- grid_regions()
  rp <- proximity_network(regions, within = 1, id = "region")
  expect_equal(nrow(rp$edges), 4L)  # side neighbours, 1 apart at the centres
  rx <- proximity_network(regions, within = 0, id = "region", exact = TRUE)
  expect_equal(nrow(rx$edges), 6L)  # every pair touches
  expect_equal(rx$edges$distance, rep(0, 6))
  expect_error(proximity_network(data.frame(x = 1, y = 1), coords = c("x", "y"), crs = "planar",
                                 within = 1, exact = TRUE), "sf nodes")
})

test_that("spatial filters, transforms, and joins keep the network consistent", {
  skip_if_not_installed("sf")
  g <- utm_network()
  g$edges$kind <- c("p", "q", "r", "s")
  box <- sf::st_bbox(c(xmin = -1, ymin = -1, xmax = 5, ymax = 5), crs = sf::st_crs(32633))
  f <- spatial_filter(g, box)
  expect_equal(f$nodes$name, c("a", "b", "c"))
  expect_equal(f$edges$kind, c("p", "q", "r"))
  expect_equal(f$roles$node_geometry, "geometry")
  expect_equal(length(f$index$src), 3L)
  expect_equal(unname(degree(f)), c(2L, 2L, 2L))
  circle <- sf::st_sf(geometry = sf::st_sfc(sf::st_buffer(sf::st_point(c(10, 0)), 1), crs = 32633))
  expect_equal(spatial_filter(g, circle)$nodes$name, "d")
  expect_equal(nrow(spatial_filter(g, circle)$edges), 0L)
  near <- spatial_filter(g, circle, predicate = sf::st_is_within_distance, dist = 7)
  expect_equal(near$nodes$name, c("b", "d"))
  expect_equal(spatial_filter(g, g)$nodes$name, g$nodes$name)
  expect_error(spatial_filter(g, 1:3), "sf")
  expect_error(spatial_filter(netframe(g$edges, directed = FALSE), box), "node geometries")

  t <- transform_network(netframe(to_sf(g, "edges"), utm_points(), directed = FALSE), 4326)
  expect_equal(t$roles$crs$epsg, 4326L)
  expect_equal(sf::st_crs(t$nodes)$epsg, 4326L)
  expect_equal(sf::st_crs(t$edges)$epsg, 4326L)
  expect_true(knappnet:::node_coordinates(t)$lonlat)
  expect_equal(as.vector(edge_lengths(t)), c(3, 5, 4, 10), tolerance = 0.01)
  expect_equal(transform_network(t, 32633)$roles$crs$epsg, 32633L)
  expect_error(transform_network(netframe(g$edges, directed = FALSE), 4326), "no geometry")

  rect <- function(x0, y0, x1, y1) {
    sf::st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
  }
  zones <- sf::st_sf(zone = c("west", "east"), kind = c("z1", "z2"),
                     geometry = sf::st_sfc(rect(-1, -1, 2, 5), rect(2, -1, 5, 5)), crs = 32633)
  j <- spatial_join(g, zones)
  expect_equal(j$nodes$zone, c("west", "east", "west", NA))
  expect_setequal(names(j$nodes), c("name", "x", "y", "kind", "kind_y", "zone", "geometry"))
  expect_equal(j$nodes$kind, g$nodes$kind)
  expect_equal(j$nodes$kind_y, c("z1", "z2", "z1", NA))
  expect_equal(nrow(j$nodes), 4L)
  expect_equal(j$roles$node_geometry, "geometry")
  expect_equal(unname(degree(j)), unname(degree(g)))
  # A node matching two zones keeps the first, with a warning.
  both <- sf::st_sf(zone = c("w", "e"), geometry = sf::st_sfc(rect(-1, -1, 2, 5), rect(-1, -1, 2, 5)),
                    crs = 32633)
  expect_warning(jj <- spatial_join(g, both), "more than one")
  expect_equal(jj$nodes$zone, c("w", NA, "w", NA))
  expect_error(spatial_join(g, 1), "sf")
})

test_that("plot_spatial() draws sf geometries", {
  skip_if_not_installed("sf")
  on_null_device({
    g <- utm_network()
    expect_silent(out <- plot_spatial(g))
    expect_identical(out, g)
    expect_silent(plot_spatial(utm_network(directed = TRUE), node_color = "kind"))
    # Explicit and mixed edge geometry, directed and not.
    edges <- to_sf(g, "edges")
    sf::st_geometry(edges)[2] <- sf::st_sfc(sf::st_linestring(), crs = 32633)
    expect_silent(plot_spatial(netframe(edges, utm_points(), directed = TRUE), edge_width = 2))
    expect_silent(plot_spatial(netframe(edges, utm_points(), directed = FALSE)))
    # Polygons, lines, multipolygons, and geometry collections as nodes.
    regions <- from_polygons(grid_regions(), id = "region")
    expect_silent(plot_spatial(regions, node_color = "pop", main = "Regions"))
    expect_silent(plot_spatial(from_lines(lattice_lines()), edge_color = "length"))
    two <- sf::st_multipolygon(list(square(0, 0), square(5, 0)))
    coll <- sf::st_geometrycollection(list(square(0, 5), sf::st_point(c(3, 7))))
    line <- sf::st_linestring(rbind(c(8, 8), c(9, 9)))
    mixed <- sf::st_sf(id = 1:4, geometry = sf::st_sfc(two, coll, line, sf::st_point(c(2, 2))),
                       crs = 32633)
    m <- netframe(data.frame(from = c(1, 2, 3), to = c(2, 3, 4)), mixed, directed = TRUE)
    expect_silent(plot_spatial(m))
    # Longitude and latitude get a latitude-dependent aspect ratio.
    cities <- netframe(data.frame(from = "London", to = "Paris"), city_points(), directed = FALSE)
    expect_silent(plot_spatial(cities))
    expect_silent(plot_edge_lengths(netframe(to_sf(cities, "edges"), city_points(),
                                             directed = FALSE), scale = 5e5))
    expect_equal(graphics::par("mar"), c(5.1, 4.1, 4.1, 2.1))
  })
})

test_that("GeoPackage files round-trip spatial networks", {
  skip_if_not_installed("sf")
  g <- utm_network(directed = TRUE)
  g$edges$kind <- c("p", "q", "r", "s")
  g$edges$w <- c(1, 2.5, 3, 4)
  g <- netframe(g$edges, g$nodes, directed = TRUE, weight = "w", weight_type = "strength",
                graph = list(title = "utm"))
  path <- tempfile(fileext = ".gpkg")
  on.exit(unlink(path), add = TRUE)
  expect_equal(write_network(g, path), path)
  expect_true(all(c("nodes", "edges", "network") %in% sf::st_layers(path)$name))
  back <- read_network(path)
  expect_true(back$directed)
  expect_s3_class(back$nodes, "sf")
  expect_false(inherits(back$edges, "sf"))
  expect_equal(back$nodes$name, g$nodes$name)
  expect_equal(back$nodes$kind, g$nodes$kind)
  expect_equal(back$edges$from, g$edges$from)
  expect_equal(back$edges$to, g$edges$to)
  expect_equal(back$edges$kind, g$edges$kind)
  expect_equal(back$roles$weight, "w")
  expect_equal(back$roles$weight_type, "strength")
  expect_equal(back$graph$title, "utm")
  expect_equal(sf::st_crs(back$nodes)$epsg, 32633L)
  expect_equal(unname(sf::st_coordinates(back$nodes)), unname(sf::st_coordinates(g$nodes)))
  expect_equal(unname(edge_lengths(back)), unname(edge_lengths(g)))
  # Both tables spatial, undirected, and lengths as weights.
  roads <- from_lines(lattice_lines())
  write_network(roads, path)
  rb <- read_network(path)
  expect_false(rb$directed)
  expect_s3_class(rb$edges, "sf")
  expect_equal(rb$roles$weight, "length")
  expect_equal(rb$edges$street, roads$edges$street)
  expect_equal(as.vector(edge_lengths(rb)), as.vector(edge_lengths(roads)))
  expect_equal(unname(distances(rb)), unname(distances(roads)))
  # Overwriting, explicit format, and arguments that override the file.
  write_network(roads, path, format = "gpkg")
  expect_true(read_network(path, directed = TRUE)$directed)
  expect_equal(read_network(path, weight = "length", weight_type = "strength")$roles$weight_type,
               "strength")
  expect_error(read_network(path, weight = "nope"), "nope")
  # A file without the network layer still reads.
  bare <- tempfile(fileext = ".gpkg")
  on.exit(unlink(bare), add = TRUE)
  sf::st_write(roads$nodes, bare, layer = "nodes", quiet = TRUE)
  sf::st_write(roads$edges, bare, layer = "edges", quiet = TRUE)
  bb <- read_network(bare)
  expect_true(bb$directed)
  expect_equal(nrow(bb$edges), 12L)
  expect_equal(bb$roles$id, "id")
  sf::st_write(roads$nodes, bare, layer = "points", quiet = TRUE, delete_dsn = TRUE)
  expect_error(read_network(bare), "layer")
  # Writing a plain network to a GeoPackage gives aspatial tables.
  plain <- netframe(data.frame(from = c("a", "b"), to = c("b", "c")), directed = FALSE)
  write_network(plain, path)
  pb <- read_network(path)
  expect_false(pb$directed)
  expect_equal(pb$nodes$name, c("a", "b", "c"))
  expect_null(pb$roles$node_geometry)
})

test_that("sfnetworks objects convert both ways", {
  skip_if_not_installed("sf")
  skip_if_not_installed("sfnetworks")
  g <- utm_network(directed = TRUE)
  g$edges$kind <- c("p", "q", "r", "s")
  net <- to_sfnetwork(g, edges_as_lines = TRUE)
  expect_s3_class(net, "sfnetwork")
  expect_equal(igraph::vcount(net), 4)
  expect_equal(igraph::ecount(net), 4)
  expect_true(igraph::is_directed(net))
  expect_equal(igraph::vertex_attr(net, "name"), c("a", "b", "c", "d"))
  expect_equal(igraph::vertex_attr(net, "kind"), c("u", "u", "v", "v"))
  expect_equal(igraph::edge_attr(net, "kind"), c("p", "q", "r", "s"))
  back <- as_netframe(net)
  expect_true(back$directed)
  expect_equal(back$nodes$name, g$nodes$name)
  expect_equal(back$nodes$kind, g$nodes$kind)
  expect_equal(back$edges$from, g$edges$from)
  expect_equal(back$edges$to, g$edges$to)
  expect_equal(back$edges$kind, g$edges$kind)
  expect_equal(back$roles$node_geometry, "geometry")
  expect_s3_class(back$edges, "sf")
  expect_equal(as.vector(edge_lengths(back)), unname(edge_lengths(g)))
  expect_equal(betweenness(back), betweenness(g))
  # Implicit edges stay implicit; polygons become points with a message.
  imp <- to_sfnetwork(g)
  expect_false(inherits(as_netframe(imp)$edges, "sf"))
  expect_equal(as_netframe(imp)$edges$to, g$edges$to)
  regions <- from_polygons(grid_regions(), id = "region")
  expect_message(rn <- to_sfnetwork(regions), "representative points")
  expect_equal(as.character(sf::st_geometry_type(sf::st_as_sf(rn, active = "nodes"),
                                                 by_geometry = TRUE)), rep("POINT", 4))
  expect_equal(as_netframe(rn)$nodes$name, c("w", "x", "y", "z"))
  expect_s3_class(sfnetworks::as_sfnetwork(g), "sfnetwork")
  # Columns with units come back as plain numbers.
  roads <- from_lines(lattice_lines())
  roads$edges$metres <- sf::st_length(roads$edges)
  rb <- as_netframe(to_sfnetwork(roads))
  expect_false(inherits(rb$edges$metres, "units"))
  expect_type(rb$edges$metres, "double")
  expect_equal(rb$edges$metres, rep(1, 12))
  expect_equal(rb$roles$weight, "length")
})

test_that("from_lines() agrees with sfnetworks on a lattice", {
  skip_if_not_installed("sf")
  skip_if_not_installed("sfnetworks")
  lines <- lattice_lines()
  ours <- from_lines(lines)
  theirs <- sfnetworks::as_sfnetwork(lines, directed = FALSE)
  expect_equal(nrow(ours$nodes), igraph::vcount(theirs))
  expect_equal(nrow(ours$edges), igraph::ecount(theirs))
  expect_equal(sort(unname(degree(ours))), sort(unname(igraph::degree(theirs))))
  expect_equal(sort(as.vector(edge_lengths(ours))),
               sort(as.numeric(sf::st_length(sf::st_as_sf(theirs, active = "edges")))))
  # The same junctions: every node of theirs is a node of ours.
  tn <- sf::st_as_sf(theirs, active = "nodes")
  hits <- sf::st_equals(sf::st_geometry(tn), sf::st_geometry(ours$nodes))
  expect_true(all(lengths(hits) == 1L))
})
