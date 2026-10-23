# The configurations of the data model (docs/DATA_MODEL.md): classes,
# card lines, and the helpers for multilayer, ego, missing-data, signed,
# and node-spell networks, with the panel plots that show them.

path4 <- function() {
  net_frame(data.frame(source = c("a", "b", "c"), target = c("b", "c", "d")), directed = FALSE)
}

two_layers <- function() {
  edges <- data.frame(source = c("ann", "ann", "ann"), target = c("bob", "bob", "ann"),
                      from_layer = c("work", "online", "work"),
                      to_layer = c("work", "online", "online"))
  layers <- data.frame(layer = c("work", "online"), medium = c("office", "chat"))
  net_frame(edges, directed = FALSE, source_layer = "from_layer", target_layer = "to_layer",
            layers = layers)
}

test_that("classes compose in a fixed order", {
  e <- data.frame(source = c("a", "b"), target = c("b", "c"), t = c(1, 2),
                  la = c("x", "x"), lb = c("x", "y"))
  nodes <- data.frame(id = c("a", "b", "c"), x = c(0, 1, 2), y = c(0, 1, 0))
  expect_identical(class(net_frame(e)), "net_frame")
  expect_identical(class(net_frame(e, start = "t")), c("temporal_net_frame", "net_frame"))
  expect_identical(class(net_frame(e, nodes, coords = c("x", "y"), crs = "planar")),
                   c("spatial_net_frame", "net_frame"))
  expect_identical(class(net_frame(e, nodes, start = "t", coords = c("x", "y"), crs = "planar")),
                   c("spatiotemporal_net_frame", "spatial_net_frame", "temporal_net_frame",
                     "net_frame"))
  expect_identical(class(net_frame(e, source_layer = "la", target_layer = "lb", start = "t")),
                   c("multilayer_net_frame", "temporal_net_frame", "net_frame"))
})

test_that("the card names modes, times, and coordinates", {
  nodes <- data.frame(id = c("p", "q", "o", "e"), m = c("person", "person", "org", "event"))
  edges <- data.frame(source = c("p", "q", "p", "o"), target = c("o", "o", "e", "e"))
  k <- net_frame(edges, nodes, directed = FALSE, mode = "m", between_modes_only = TRUE)
  expect_output(print(k), "mode: 3-mode (2 person, 1 org, 1 event; between modes only)",
                fixed = TRUE)
  expect_error(net_frame(rbind(edges, data.frame(source = "p", target = "q")), nodes,
                         directed = FALSE, mode = "m", between_modes_only = TRUE), "same mode")
  calls <- data.frame(source = c("a", "b", "a"), target = c("b", "c", "c"), t = c(1, 2, 4))
  ct <- net_frame(calls, directed = FALSE, start = "t")
  expect_output(print(ct), "time: +contacts \\(3 events")
  expect_output(print(ct), "clock: +numbers")
  spells <- transform(calls, u = t + 1)
  expect_output(print(net_frame(spells, start = "t", end = "u")), "time: +edge spells, 1 to 5")
  panel <- net_frame(transform(calls, w = c(1L, 1L, 2L)), wave = "w")
  expect_output(print(panel), "2 waves")
  expect_output(print(panel), "wave k is [k, k + 1)", fixed = TRUE)
  nodes <- data.frame(id = c("a", "b", "c"), lon = c(0, 1, 2), lat = c(50, 51, 52))
  sp <- net_frame(calls, nodes, directed = FALSE, coords = c("lon", "lat"))
  expect_output(print(sp), "node coordinates `lon`, `lat` (longitude and latitude)",
                fixed = TRUE)
  expect_output(print(sp), "distances: great-circle")
})

test_that("multilayer networks classify edges and split into layers", {
  g <- two_layers()
  expect_true(net_is_multilayer(g))
  expect_false(net_is_multilayer(path4()))
  expect_false(net_is_multilayer(data.frame()))
  el <- net_edge_layers(g)
  expect_identical(el$class, c("intralayer", "intralayer", "coupling"))
  expect_identical(el$source_layer, c("work", "online", "work"))
  expect_identical(el$target_layer, c("work", "online", "online"))
  # A coupling is not a loop, and edges in different layers are not parallel.
  expect_true(net_is_simple(g))
  expect_output(print(g), paste0("layers:    work, online (2); 2 intralayer, 1 coupling, ",
                                 "0 interlayer (layer table: 2 x 2)"), fixed = TRUE)
  parts <- net_split_layers(g)
  expect_named(parts, c("work", "online"))
  expect_identical(nrow(parts$work$edges), 1L)
  expect_identical(nrow(parts$work$nodes), 2L)
  expect_false(net_is_multilayer(parts$online))
  expect_error(net_edge_layers(path4()), "multilayer")
  expect_error(net_frame(data.frame(source = "a", target = "b", la = "x"), source_layer = "la"),
               "layer")
})

test_that("net_as_multilayer() lifts a multiplex network", {
  mp <- net_frame(data.frame(source = c("a", "a", "b"), target = c("b", "b", "c"),
                             rel = c("x", "y", "y")), directed = FALSE, relation = "rel")
  ml <- net_as_multilayer(mp)
  expect_s3_class(ml, "multilayer_net_frame")
  expect_true(ml$declared$couplings_only)
  expect_identical(ml$layers$layer, c("x", "y"))
  cls <- net_edge_layers(ml)$class
  expect_identical(sum(cls == "intralayer"), 3L)
  expect_identical(sum(cls == "coupling"), 3L)  # one pair of layers, three actors
  expect_identical(nrow(net_as_multilayer(mp, couplings = "none")$edges), 3L)
  d <- net_frame(mp$edges, relation = "rel")
  expect_identical(sum(net_edge_layers(net_as_multilayer(d))$class == "coupling"), 6L)
  expect_error(net_as_multilayer(ml), "already")
  expect_error(net_as_multilayer(path4()), "relation")
})

test_that("missing ties and non-respondents are kept, counted, and resolved explicitly", {
  e <- data.frame(source = c("ann", "ann", "bob", "cy"), target = c("bob", "cy", "cy", "dee"),
                  unknown = c(FALSE, TRUE, FALSE, FALSE))
  g <- net_frame(e, directed = FALSE, edge_missing = "unknown", quiet = TRUE)
  expect_identical(nrow(g$edges), 4L)
  expect_output(print(g), "missing:   0 non-respondents, 1 unobserved pair", fixed = TRUE)
  expect_error(net_resolve_missing(g), "na_action")
  expect_error(net_resolve_missing(g, "maybe"), "na_action")
  a <- net_resolve_missing(g, "absent")
  expect_identical(nrow(a$edges), 3L)
  expect_null(a$roles$edge_missing)
  p <- net_resolve_missing(g, "present")
  expect_identical(nrow(p$edges), 4L)
  expect_equal(net_density(p), 4 / 6)
  expect_error(net_frame(transform(e, unknown = c(FALSE, NA, FALSE, FALSE)),
                         edge_missing = "unknown"), "unknown")
  people <- data.frame(id = c("ann", "bob", "cy", "dee"), no_reply = c(FALSE, FALSE, TRUE, FALSE))
  h <- net_frame(e[-2, ], people, node_missing = "no_reply", quiet = TRUE)
  expect_output(print(h), "1 non-respondent, 0 unobserved pairs", fixed = TRUE)
  expect_error(net_resolve_missing(h, "present"), "non-respondent")
  expect_null(net_resolve_missing(h, "absent")$roles$node_missing)
})

test_that("signed values split into sign and magnitude", {
  r <- data.frame(source = c("a", "a", "b"), target = c("b", "c", "c"), v = c(3, -2, 1))
  s <- net_as_signed(r, "v")
  expect_identical(s$sign, c(1L, -1L, 1L))
  expect_identical(s$weight, c(3, 2, 1))
  g <- net_frame(s, directed = FALSE, sign = "sign", weight = "weight", weight_type = "strength")
  expect_output(print(g), "signs:     sign (2 positive, 1 negative)", fixed = TRUE)
  expect_error(net_as_signed(transform(r, v = c(0, 1, 2)), "v"), "zero")
  expect_error(net_as_signed(transform(r, v = "x"), "v"), "numeric")
  expect_error(net_frame(transform(s, weight = -weight), sign = "sign", weight = "weight",
                         weight_type = "strength"), "negative")
})

test_that("ego networks are held in one net_frame", {
  e <- net_ego_frame(path4(), egos = c("a", "c"), threads = 1)
  expect_identical(nrow(e$nodes), 5L)
  expect_identical(nrow(e$edges), 3L)
  expect_identical(e$nodes$id, c("a", "a:b", "c:b", "c", "c:d"))
  expect_identical(e$nodes$node, c("a", "b", "b", "c", "d"))
  expect_identical(e$nodes$ego, c("a", "a", "c", "c", "c"))
  expect_identical(e$edges$ego, c("a", "c", "c"))
  expect_true(e$declared$egos_are_nodes)
  expect_output(print(e), "egos: +2 \\(egos are nodes")
  parts <- net_split_egos(e)
  expect_named(parts, c("a", "c"))
  expect_identical(nrow(parts$c$nodes), 3L)
  expect_null(parts$c$roles$ego)
  alters <- net_ego_frame(path4(), egos = "c", include_ego = FALSE, threads = 1)
  expect_identical(alters$nodes$id, c("c:b", "c:d"))
  expect_identical(nrow(alters$edges), 0L)
  expect_error(net_ego_frame(e), "already")
  expect_error(net_split_egos(path4()), "ego role")
  # Edges across egos are errors.
  nodes <- data.frame(id = c("a", "a1", "b", "b1"), ego = c("a", "a", "b", "b"))
  bad <- data.frame(source = c("a", "a1"), target = c("a1", "b1"))
  expect_error(net_frame(bad, nodes, ego = "ego"), "ego")
})

test_that("node spells and node changes", {
  nodes <- data.frame(id = c("ann", "bob", "cy"), status = "guest")
  spells <- data.frame(id = c("ann", "ann", "bob"), start = c(0, 8, 0), end = c(5, 20, 20))
  changes <- data.frame(id = c("ann", "bob"), time = c(8, 12), status = c("member", "staff"))
  edges <- data.frame(source = c("ann", "ann", "bob"), target = c("bob", "bob", "cy"),
                      from = c(1, 12, 3), to = c(4, 15, 3))
  g <- net_frame(edges, nodes, directed = FALSE, start = "from", end = "to",
                 node_spells = spells, node_changes = changes)
  expect_s3_class(g, "temporal_net_frame")
  expect_output(print(g), "node spells (table, 3 spells); 2 changes", fixed = TRUE)
  sp <- net_node_spells(g)
  ann <- sp[sp$id == "ann", ]
  expect_equal(ann$start, c(0, 8))
  expect_equal(ann$end, c(5, 20))
  expect_identical(net_nodes_at(g, 10)$status, c("member", "guest", "guest"))
  expect_identical(net_nodes_at(g, 13)$status, c("member", "staff", "guest"))
  expect_identical(net_nodes_at(g, 1)$status, c("guest", "guest", "guest"))
  edges$to[1] <- 10
  expect_error(net_frame(edges, nodes, directed = FALSE, start = "from", end = "to",
                         node_spells = spells), "row 1")
  overlap <- rbind(spells, data.frame(id = "bob", start = 10, end = 30))
  expect_error(net_frame(edges[-1, ], nodes, directed = FALSE, start = "from", end = "to",
                         node_spells = overlap), "overlap")
})

test_that("simplify combines other columns on request", {
  e <- data.frame(source = c("a", "a", "b"), target = c("b", "b", "c"), w = c(1, 2, 1),
                  note = c("met", "called", "met"))
  g <- net_frame(e, directed = FALSE, weight = "w", weight_type = "strength")
  expect_identical(net_simplify(g)$edges$note, c("met", "met"))
  expect_identical(net_simplify(g, attributes = "concat")$edges$note, c("met; called", "met"))
  n <- net_simplify(g, attributes = function(v) length(v))
  expect_equal(n$edges$note, c(2, 1))
  expect_error(net_simplify(g, attributes = "glue"), "attributes")
})

test_that("panel plots draw layers, snapshots, and egos", {
  on_null_device({
    expect_restores_par(expect_silent(net_plot_layers(two_layers(), seed = 1)))
    calls <- data.frame(source = c("a", "b", "a", "c"), target = c("b", "c", "c", "d"),
                        t = c(1, 2, 4, 6))
    ct <- net_frame(calls, directed = FALSE, start = "t")
    expect_restores_par(expect_silent(net_plot_snapshots(ct, seed = 1)))
    expect_silent(net_plot_snapshots(ct, width = 3, ncol = 1, seed = 1))
    panel <- net_frame(transform(calls, w = c(1L, 1L, 2L, 3L)), wave = "w")
    expect_silent(net_plot_snapshots(panel, seed = 1))
    nets <- net_plot_egos(path4(), egos = c("a", "c"), seed = 1)
    expect_named(nets, c("a", "c"))
    expect_silent(net_plot_egos(net_ego_frame(path4(), threads = 1), max_panels = 2))
    expect_error(net_plot_layers(path4()), "multilayer")
    expect_error(net_plot_snapshots(path4()), "temporal")
  })
})

test_that("plots show line types, arcs, loops, and unobserved pairs", {
  on_null_device({
    m <- net_frame(data.frame(source = c("a", "a", "a", "b"), target = c("b", "b", "a", "c"),
                              kind = c("x", "y", "x", "y")))
    expect_silent(plot(m, edge_linetype = "kind", seed = 1))
    expect_silent(plot(m, curved = TRUE, seed = 1))
    expect_silent(plot(m, curved = FALSE, edge_linetype = 2, seed = 1))
    expect_error(plot(m, curved = "yes"), "curved")
    s <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c"), s = c(1, -1)),
                   sign = "s")
    expect_silent(plot(s, seed = 1))
    u <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c"), na = c(FALSE, TRUE)),
                   edge_missing = "na", quiet = TRUE)
    expect_silent(plot(u, seed = 1))
    expect_silent(plot(two_layers(), seed = 1))
  })
})
