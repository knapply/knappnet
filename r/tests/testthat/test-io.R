rich_net <- function(directed = TRUE) {
  nodes <- data.frame(id = c("ana", "bo \"b\"", "céline", "d&d", "isolate"),
                      age = c(31L, 45L, NA, 22L, 60L),
                      score = c(0.5, NA, 1.25, -3, 1e-7),
                      member = c(TRUE, FALSE, NA, TRUE, FALSE),
                      team = c("red", "blue", "<blue>", NA, "green & gold"),
                      stringsAsFactors = FALSE)
  edges <- data.frame(source = c("ana", "bo \"b\"", "céline", "ana"),
                      target = c("bo \"b\"", "céline", "céline", "d&d"),
                      weight = c(1, 2.5, 0.75, 3),
                      stringsAsFactors = FALSE)
  net_frame(edges, nodes, directed = directed, weight = "weight", weight_type = "distance",
           graph = list(title = "rich"))
}

round_trip <- function(g, ext, ...) {
  path <- tempfile(fileext = ext)
  on.exit(unlink(path))
  net_write(g, path, ...)
  net_read(path, weight_type = g$roles$weight_type)
}

test_that("attribute-rich formats round-trip exactly", {
  for (ext in c(".graphml", ".gexf", ".json", ".gml")) {
    for (directed in c(TRUE, FALSE)) {
      g <- rich_net(directed)
      back <- round_trip(g, ext)
      expect_equal(back$directed, directed, info = ext)
      expect_equal(back$nodes$id, g$nodes$id, info = ext)
      expect_equal(back$edges$source, g$edges$source, info = ext)
      expect_equal(back$edges$weight, g$edges$weight, info = ext)
      expect_equal(back$nodes$score, g$nodes$score, info = ext)
      expect_equal(back$nodes$team, g$nodes$team, info = ext)
      if (ext != ".gml") {
        # GML has no boolean type, so logical columns come back as 0/1.
        expect_equal(back$nodes$member, g$nodes$member, info = ext)
        expect_identical(back$nodes$age, g$nodes$age, info = ext)
      }
    }
  }
})

test_that("structure-only formats keep structure and weights", {
  g <- net_frame(data.frame(source = c("a", "b", "c"), target = c("b", "c", "c"), w = c(2, 0.5, 1)),
                nodes = data.frame(id = c("a", "b", "c", "lonely")), weight = "w", weight_type = "distance")
  for (ext in c(".net", ".dl")) {
    back <- round_trip(g, ext)
    expect_equal(back$nodes$id, g$nodes$id, info = ext)
    expect_equal(back$edges$target, g$edges$target, info = ext)
    expect_equal(back$edges$weight, g$edges$w, info = ext)
  }
})

test_that("edge lists read and write with optional node files", {
  g <- rich_net()
  path <- tempfile(fileext = ".csv")
  nodes_path <- tempfile(fileext = ".csv")
  net_write(g, path, nodes_file = nodes_path)
  back <- net_read(path, nodes = nodes_path, weight_type = "distance")
  expect_equal(nrow(back$nodes), 5)
  expect_equal(back$edges$source, g$edges$source)
  expect_equal(back$roles$weight, "weight")
})

test_that("integer identifiers come back as integers", {
  g <- net_frame(data.frame(source = c(1L, 2L), target = c(2L, 3L)))
  back <- round_trip(g, ".graphml")
  expect_identical(back$nodes$id, 1:3)
  expect_identical(back$edges$target, c(2L, 3L))
})

test_that("DOT is write-only and formats are validated", {
  g <- rich_net()
  path <- tempfile(fileext = ".dot")
  net_write(g, path)
  expect_match(readLines(path)[1], "^digraph")
  expect_error(net_read(path), "written but not read")
  expect_error(net_write(g, tempfile(fileext = ".xyz")), "Cannot tell the format")
  expect_error(net_read("does-not-exist.graphml"), "does not exist")
})

test_that("files written here are read correctly by igraph", {
  skip_if_not_installed("igraph")
  g <- rich_net()
  for (fmt in c("graphml", "gml", "pajek")) {
    path <- tempfile()
    net_write(g, path, format = fmt)
    # igraph keeps GML's character entities (here for an accented label)
    # verbatim and warns; the counts and weights checked below are unaffected.
    ig <- tryCatch(
      if (fmt == "gml") suppressWarnings(igraph::read_graph(path, format = fmt)) else
        igraph::read_graph(path, format = fmt),
      error = function(e) {
        # Some igraph builds are compiled without GraphML support.
        if (grepl("GraphML", conditionMessage(e), ignore.case = TRUE)) NULL else stop(e)
      })
    if (is.null(ig)) next
    expect_equal(igraph::vcount(ig), nrow(g$nodes), info = fmt)
    expect_equal(igraph::ecount(ig), nrow(g$edges), info = fmt)
    expect_equal(sort(igraph::E(ig)$weight), sort(g$edges$weight), info = fmt)
  }
  path <- tempfile(fileext = ".graphml")
  net_write(g, path)
  ig <- tryCatch(igraph::read_graph(path, format = "graphml"), error = function(e) NULL)
  skip_if(is.null(ig), "this igraph build cannot read GraphML")
  # Some libxml2 versions hand igraph `&amp;` in attributes as the
  # character reference `&#38;` without decoding it; the file itself holds
  # `&amp;`, which every conforming parser reads as `&`.
  expect_equal(gsub("&#38;", "&", igraph::V(ig)$id, fixed = TRUE), g$nodes$id)
  expect_equal(igraph::V(ig)$age[-3], g$nodes$age[-3])
})

test_that("files written by igraph are read correctly here", {
  skip_if_not_installed("igraph")
  ig <- igraph::make_graph("Zachary")
  igraph::E(ig)$weight <- seq_len(igraph::ecount(ig)) / 10
  igraph::V(ig)$club <- rep(c("A", "B"), length.out = igraph::vcount(ig))
  for (fmt in c("graphml", "gml", "pajek")) {
    path <- tempfile()
    if (fmt == "gml") {
      # Some igraph versions need GML node ids given explicitly.
      igraph::write_graph(ig, path, format = fmt, id = seq_len(igraph::vcount(ig)) - 1)
    } else {
      igraph::write_graph(ig, path, format = fmt)
    }
    g <- net_read(path, format = fmt, weight_type = "distance")
    expect_equal(nrow(g$nodes), igraph::vcount(ig), info = fmt)
    expect_equal(nrow(g$edges), igraph::ecount(ig), info = fmt)
    expect_equal(g$edges$weight, igraph::E(ig)$weight, info = fmt)
    expect_equal(net_betweenness(g, weighted = FALSE),
                 unname(igraph::betweenness(ig, weights = NA)), info = fmt)
    expect_equal(net_betweenness(g), unname(igraph::betweenness(ig)), info = fmt)
  }
})

test_that("UCINET DL files are read the same way igraph reads them", {
  skip_if_not_installed("igraph")
  path <- tempfile(fileext = ".dl")
  # igraph reads full matrices as binary only, so values are compared on
  # an edge list.
  writeLines(c("dl n=4", "format = fullmatrix", "labels:", "a,b,c,d", "data:",
               "0 1 0 0", "0 0 1 0", "1 0 0 0", "0 0 0 0"), path)
  g <- net_read(path, weight_type = "strength")
  ig <- igraph::read_graph(path, format = "dl")
  expect_equal(nrow(g$edges), igraph::ecount(ig))
  expect_equal(g$nodes$id, igraph::V(ig)$name)
  expect_equal(paste(g$edges$source, g$edges$target),
               apply(igraph::as_edgelist(ig), 1, paste, collapse = " "))

  writeLines(c("dl n=4", "format = edgelist1", "labels:", "a,b,c,d", "data:",
               "1 2 1", "2 3 2.5", "3 1 1"), path)
  g <- net_read(path, weight_type = "strength")
  ig <- igraph::read_graph(path, format = "dl")
  expect_equal(nrow(g$edges), igraph::ecount(ig))
  expect_equal(g$nodes$id, igraph::V(ig)$name)
  expect_equal(g$edges$weight, igraph::E(ig)$weight)
})

multimode_net <- function() {
  nodes <- data.frame(id = c("ana", "bo", "acme", "kabul", "herat"),
                      class = c("Agent", "Agent", "Organization", "Location", "Location"),
                      title = c("Ana A.", NA, "ACME & Sons", NA, "<Herat>"),
                      age = c(31, NA, NA, NA, NA),
                      stringsAsFactors = FALSE)
  edges <- data.frame(source = c("ana", "bo", "ana", "bo", "acme"),
                      target = c("bo", "ana", "acme", "acme", "kabul"),
                      network = c("talks", "talks", "member", "member", "based"),
                      weight = c(2.5, 1, 1, 0.5, 1),
                      stringsAsFactors = FALSE)
  net_frame(edges, nodes, directed = TRUE, weight = "weight", weight_type = "distance")
}

test_that("DyNetML round-trips multimode networks", {
  g <- multimode_net()
  path <- tempfile(fileext = ".xml")
  on.exit(unlink(path))
  net_write(g, path, format = "dynetml")
  expect_match(paste(readLines(path), collapse = "\n"),
               "<nodeclass type=\"Organization\" id=\"Organization\">", fixed = TRUE)
  # A .xml file is recognized as DyNetML from its root element.
  back <- net_read(path, weight_type = "distance")
  expect_true(back$directed)
  expect_equal(back$nodes$id, g$nodes$id)
  expect_equal(back$nodes$class, g$nodes$class)
  expect_equal(back$nodes$type, g$nodes$class)
  expect_equal(back$nodes$title, g$nodes$title)
  expect_equal(back$nodes$age, g$nodes$age)
  expect_equal(back$edges$source, g$edges$source)
  expect_equal(back$edges$target, g$edges$target)
  expect_equal(back$edges$network, g$edges$network)
  expect_equal(back$edges$weight, g$edges$weight)
  expect_equal(back$roles$weight, "weight")
  ora <- net_read(path, format = "ora", weight_type = "distance")
  expect_identical(ora$edges$weight, back$edges$weight)
})

test_that("DyNetML files in ORA's layout are read", {
  path <- tempfile(fileext = ".xml")
  on.exit(unlink(path))
  writeLines(c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<DynamicNetwork>',
    '  <MetaNetwork id="embassy">',
    '    <nodes>',
    '      <nodeclass type="Agent" id="Agent">',
    '        <propertyIdentities>',
    '          <propertyIdentity id="age" type="number" singleValued="true"/>',
    '        </propertyIdentities>',
    '        <node id="A1" title="Ann"><properties><property id="age" value="34"/></properties></node>',
    '        <node id="Kabul"/>',
    '      </nodeclass>',
    '      <nodeclass type="Location" id="Location"><node id="Kabul"/></nodeclass>',
    '    </nodes>',
    '    <networks>',
    '      <network sourceType="Agent" source="Agent" targetType="Agent" target="Agent"',
    '               id="agent x agent" isDirected="true">',
    '        <link source="A1" target="Kabul" value="2"/>',
    '      </network>',
    '      <network sourceType="Agent" source="Agent" targetType="Location" target="Location"',
    '               id="agent x location" isDirected="false">',
    '        <link source="A1" target="Kabul" value="1"/>',
    '      </network>',
    '    </networks>',
    '  </MetaNetwork>',
    '</DynamicNetwork>'), path)
  g <- net_read(path, weight_type = "strength")
  # "Kabul" is both an agent and a location, so identifiers carry the class.
  expect_equal(g$nodes$id, c("Agent::A1", "Agent::Kabul", "Location::Kabul"))
  expect_equal(g$nodes$node_id, c("A1", "Kabul", "Kabul"))
  expect_equal(g$nodes$age, c(34, NA, NA))
  expect_equal(g$graph$name, "embassy")
  # The undirected network's link is stored in both directions.
  expect_true(g$directed)
  expect_equal(g$edges$source, c("Agent::A1", "Agent::A1", "Location::Kabul"))
  expect_equal(g$edges$directed, c(TRUE, FALSE, FALSE))
  expect_equal(g$edges$weight, c(2, 1, 1))
  bad <- tempfile(fileext = ".xml")
  on.exit(unlink(bad), add = TRUE)
  writeLines(c('<MetaNetwork>', '<nodes><nodeclass id="Agent" type="Agent">',
               '<node id="a"><port name="p"/></node></nodeclass></nodes></MetaNetwork>'), bad)
  expect_error(net_read(bad), "line 3.*ports")
})

test_that("writing .xml stays GraphML and reading .xml still reads GraphML", {
  g <- rich_net()
  path <- tempfile(fileext = ".xml")
  on.exit(unlink(path))
  net_write(g, path)
  expect_match(paste(readLines(path), collapse = "\n"), "<graphml")
  back <- net_read(path, weight_type = "distance")
  expect_equal(back$nodes$id, g$nodes$id)
  expect_null(back$nodes$class)
})

test_that("non-ASCII identifiers round-trip through every reader and writer", {
  # Latin, Greek, and Chinese text, written with escapes so this file is
  # ASCII; the strings are UTF-8 on every platform.
  ids <- c("José", "Zoë", "北京", "Ελλάδα",
           "São Paulo")
  g <- net_frame(data.frame(source = ids[c(1, 2, 3, 4)], target = ids[c(2, 3, 4, 1)],
                           weight = c(1, 2, 3, 4), stringsAsFactors = FALSE),
                nodes = data.frame(id = ids, place = rev(ids), stringsAsFactors = FALSE),
                weight = "weight", weight_type = "distance")
  for (ext in c(".graphml", ".gexf", ".gml", ".json", ".net", ".dl")) {
    back <- round_trip(g, ext)
    expect_equal(back$nodes$id, ids, info = ext)
    expect_equal(back$edges$source, g$edges$source, info = ext)
    expect_equal(back$edges$target, g$edges$target, info = ext)
    expect_equal(net_betweenness(back), net_betweenness(g), info = ext)
    if (ext %in% c(".graphml", ".gexf", ".gml", ".json")) {
      expect_equal(back$nodes$place, rev(ids), info = ext)
    }
  }
  path <- tempfile(fileext = ".xml")
  on.exit(unlink(path))
  net_write(g, path, format = "dynetml")
  back <- net_read(path, weight_type = "distance")
  expect_equal(back$nodes$id, ids)
  expect_equal(back$edges$target, g$edges$target)

  # Edge lists go through write.table(), which before R 4.2 on Windows (and
  # in other non-UTF-8 locales) can only write characters of the native
  # encoding.
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]), "not a UTF-8 locale")
  csv <- tempfile(fileext = ".csv")
  nodes_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(c(csv, nodes_csv)), add = TRUE)
  net_write(g, csv, nodes_file = nodes_csv)
  back <- net_read(csv, nodes = nodes_csv, weight_type = "distance")
  expect_equal(back$nodes$id, ids)
  expect_equal(back$edges$source, g$edges$source)
})

test_that("file paths may contain non-ASCII characters", {
  skip_if_not(isTRUE(l10n_info()[["UTF-8"]]) || .Platform$OS.type == "windows",
              "the locale cannot represent the path")
  g <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c")))
  dir <- file.path(tempdir(), "réseau")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE))
  for (ext in c(".graphml", ".csv")) {
    path <- file.path(dir, paste0("café", ext))
    net_write(g, path)
    expect_true(file.exists(path), info = ext)
    expect_equal(net_read(path)$edges$target, c("b", "c"), info = ext)
  }
})
