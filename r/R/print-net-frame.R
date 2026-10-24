# The print card of a net_frame (docs/DATA_MODEL.md section 18): a
# labelled card filled from the engine's summary, one line per
# configuration that is present (size, direction and modes, weights,
# signs, relations, layers, time, space, egos, missing data, loops and
# parallel edges, components, isolates, density), followed by previews of
# the node and edge tables, side by side when the console is wide enough.
# Base R only: format(), formatC(), and strtrim(), with no characters
# outside ASCII, so the output is the same on every console. The
# subclasses add a line each saying how times and distances are read.

# See man/net_frame.Rd.
print.net_frame <- function(x, ..., n = 6L, width = getOption("width")) {
  print_card(x, n, width)
}

# See man/net_temporal.Rd.
print.temporal_net_frame <- function(x, ..., n = 6L, width = getOption("width")) {
  print_card(x, n, width)
}

# See man/net_spatial.Rd.
print.spatial_net_frame <- function(x, ..., n = 6L, width = getOption("width")) {
  print_card(x, n, width)
}

# See man/net_temporal.Rd.
print.spatiotemporal_net_frame <- function(x, ..., n = 6L, width = getOption("width")) {
  print_card(x, n, width)
}

# See man/net_split_layers.Rd.
print.multilayer_net_frame <- function(x, ..., n = 6L, width = getOption("width")) {
  print_card(x, n, width)
}

# Prints the card, the subclass lines, and the table previews.
print_card <- function(x, n, width) {
  check_net_frame(x)
  n <- check_count(n, "n")
  width <- check_count(width, "width", min = 20L)
  lines <- c(card_lines(x), "", preview_lines(x, n, width - 2L))
  cat(paste0("#", ifelse(nzchar(lines), " ", ""), lines), sep = "\n")
  invisible(x)
}

# ---- The card ---------------------------------------------------------------------

# One labelled line of the card: the label padded so values line up.
card_item <- function(label, value) sprintf("  %s %s", format(paste0(label, ":"), width = 10L), value)

# The lines of the card, without the "# " prefix.
card_lines <- function(x) {
  s <- frame_summary_cached(x)
  direction <- card_item("directed", if (s$directed) "yes" else "no")
  temporal <- s$time != "static" || s$node_spells != "none" || s$change_rows > 0
  spatial <- s$node_geometry || s$edge_geometry || !is.null(x$roles$x)
  c(
    sprintf("%s: %s, %s", class(x)[[1L]], plural(as.integer(s$nodes), "node"),
            plural(as.integer(s$edges), "edge")),
    paste0(format(direction, width = 22L), sprintf(" mode: %s", mode_text(s))),
    card_item("weights", weight_text(s)),
    if (s$signed) card_item("signs", sign_text(x, s)),
    if (length(s$relations)) card_item("relations", relation_text(s)),
    if (length(s$layers)) card_item("layers", layer_text(x, s)),
    if (temporal || spatial) {
      c(card_item("time", time_text(x, s)), card_item("space", space_text(x, s)))
    } else {
      paste0(format(card_item("time", "static"), width = 22L), " space: none")
    },
    if (inherits(x, "temporal_net_frame")) card_item("clock", clock_text(x)),
    if (inherits(x, "spatial_net_frame")) card_item("distances", distance_text(x)),
    if (length(s$egos)) card_item("egos", ego_text(s)),
    if (s$missing_nodes > 0 || s$missing_edges > 0) card_item("missing", missing_text(s)),
    card_item("simple", simple_text(x, s)),
    sprintf("  components: %s    isolates: %s    density: %s", format(s$components),
            format(s$isolates), number_text3(s$density))
  )
}

number_text3 <- function(v) {
  if (is.na(v)) "NaN" else formatC(v, format = "f", digits = 3)
}

# Up to `k` labels with a count after each, then how many more.
label_counts <- function(labels, counts, k = 6L, count_first = FALSE) {
  counts <- format(as.integer(counts), trim = TRUE)
  items <- if (count_first) paste(counts, labels) else sprintf("%s (%s)", labels, counts)
  if (length(items) > k) {
    items <- c(items[seq_len(k)], sprintf("and %d more", length(items) - k))
  }
  paste(items, collapse = ", ")
}

# A numeric range as "a to b".
range_text <- function(lo, hi) {
  if (is.na(lo) || is.na(hi)) return("no values")
  if (lo == hi) return(format(lo, digits = 4))
  sprintf("%s to %s", format(lo, digits = 4), format(hi, digits = 4))
}

mode_text <- function(s) {
  k <- length(s$modes)
  if (k == 0L) return(sprintf("one-mode (%s)", plural(as.integer(s$nodes), "actor")))
  counts <- label_counts(s$modes, s$mode_counts, count_first = TRUE)
  rule <- if (isTRUE(s$between_modes_only)) {
    "; between modes only"
  } else if (isFALSE(s$between_modes_only)) {
    "; within-mode ties allowed"
  } else {
    ""
  }
  kind <- if (isTRUE(s$between_modes_only)) {
    if (k == 2L) "two-mode" else sprintf("%d-mode", k)
  } else {
    plural(k, "mode")
  }
  sprintf("%s (%s%s)", kind, counts, rule)
}

weight_text <- function(s) {
  if (!s$weighted) return("unweighted")
  sprintf("%s (%s, %s)", s$weight, s$weight_type, range_text(s$weight_min, s$weight_max))
}

sign_text <- function(x, s) {
  neg <- as.integer(s$negative)
  pos <- as.integer(s$edges) - neg
  sprintf("%s (%d positive, %d negative)", x$roles$sign, pos, neg)
}

relation_text <- function(s) {
  sprintf("%s (multiplex)", label_counts(s$relations, s$relation_counts))
}

layer_text <- function(x, s) {
  labels <- s$layers
  shown <- if (length(labels) > 6L) {
    paste0(paste(labels[1:6], collapse = ", "), sprintf(", and %d more", length(labels) - 6L))
  } else {
    paste(labels, collapse = ", ")
  }
  table <- if (!is.null(x$layers)) {
    sprintf(" (layer table: %d x %d)", nrow(x$layers), ncol(x$layers))
  } else {
    ""
  }
  sprintf("%s (%d); %s intralayer, %s coupling, %s interlayer%s", shown, length(labels),
          format(s$intralayer), format(s$coupling), format(s$interlayer), table)
}

time_text <- function(x, s) {
  range <- if (is.na(s$time_start)) "no edges" else {
    sprintf("%s to %s", format(from_seconds(s$time_start, x)),
            format(from_seconds(s$time_end, x)))
  }
  edges <- switch(s$time,
    static = "static edges",
    contacts = sprintf("contacts (%s, %s)", plural(as.integer(s$edges), "event"), range),
    spells = sprintf("edge spells, %s", range),
    waves = sprintf("%s (%s)", plural(length(unique(x$edges[[x$roles$wave]])), "wave"), range)
  )
  nodes <- switch(s$node_spells,
    none = NULL,
    roles = "node spells (one per node)",
    table = sprintf("node spells (table, %s)", plural(as.integer(s$spell_rows), "spell"))
  )
  changes <- if (s$change_rows > 0) plural(as.integer(s$change_rows), "change")
  paste(c(edges, nodes, changes), collapse = "; ")
}

space_text <- function(x, s) {
  r <- x$roles
  if (s$node_geometry || s$edge_geometry) {
    kinds <- c(
      if (s$node_geometry) sprintf("node %s", geometry_plural(x$nodes[[r$node_geometry]])),
      if (s$edge_geometry) sprintf("edge %s", geometry_plural(x$edges[[r$edge_geometry]]))
    )
    crs <- tryCatch(crs_label(r$crs), error = function(e) "a CRS")
    return(sprintf("%s (%s)", paste(kinds, collapse = ", "), crs))
  }
  if (!is.null(r$x) && !is.null(r$y)) {
    return(sprintf("node coordinates `%s`, `%s` (%s)", r$x, r$y,
                   if (identical(s$coordinates, "lonlat")) "longitude and latitude" else "planar"))
  }
  "none"
}

# The geometry type of an sfc column from its class, so sf is not needed
# to print.
geometry_kind <- function(geom) {
  cls <- class(geom)[[1L]]
  kind <- sub("^sfc_", "", cls)
  if (!startsWith(cls, "sfc_") || kind == "GEOMETRY") "mixed" else kind
}

geometry_plural <- function(geom) {
  switch(geometry_kind(geom),
    POINT = "points", MULTIPOINT = "multipoints", LINESTRING = "lines",
    MULTILINESTRING = "multilines", POLYGON = "polygons", MULTIPOLYGON = "multipolygons",
    "geometries"
  )
}

ego_text <- function(s) {
  k <- length(s$egos)
  all_nodes <- s$egos_as_nodes == k
  alters <- if (all_nodes) s$ego_nodes - 1 else s$ego_nodes
  each <- sprintf("%s alters each", range_text(min(alters), max(alters)))
  kind <- if (all_nodes) {
    "egos are nodes"
  } else if (s$egos_as_nodes == 0) {
    "alters only"
  } else {
    sprintf("%d egos are nodes", as.integer(s$egos_as_nodes))
  }
  sprintf("%d (%s; %s)", k, kind, each)
}

missing_text <- function(s) {
  sprintf("%s, %s", plural(as.integer(s$missing_nodes), "non-respondent"),
          plural(as.integer(s$missing_edges), "unobserved pair"))
}

simple_text <- function(x, s) {
  declared <- function(flag) if (isFALSE(flag)) " (declared)" else ""
  loops <- if (s$loops == 0) {
    paste0("no loops", declared(x$declared$loops))
  } else {
    plural(as.integer(s$loops), "loop")
  }
  parallel <- if (s$parallel == 0) {
    paste0("no parallel edges", declared(x$declared$multiple))
  } else {
    plural(as.integer(s$parallel), "parallel edge")
  }
  paste(loops, parallel, sep = ", ")
}

# The line a temporal network adds: how times are stored.
clock_text <- function(x) {
  r <- x$roles
  switch(if (is.null(r$time_class)) "numeric" else r$time_class,
    numeric = if (!is.null(r$wave)) "wave numbers; wave k is [k, k + 1)" else "numbers",
    Date = "dates",
    POSIXct = sprintf("date-times (%s)", if (is.null(r$tz) || !nzchar(r$tz)) "local time" else r$tz)
  )
}

# The line a spatial network adds: how distances are measured.
distance_text <- function(x) {
  r <- x$roles
  kind <- tryCatch(crs_kind(x), error = function(e) "planar")
  distances <- if (identical(kind, "lonlat")) {
    sprintf("great-circle, sphere radius %s m",
            format(if (is.null(r$radius)) 6371008.7714 else r$radius, digits = 8))
  } else {
    "planar, in the coordinate units"
  }
  points <- if (!is.null(r$node_geometry)) {
    sprintf("; representative points by %s",
            if (identical(r$representative, "centroid")) "centroid" else "point on surface")
  } else {
    ""
  }
  paste0(distances, points)
}

# ---- Table previews ----------------------------------------------------------------

# The preview of the two tables: side by side when both fit in `width`,
# else one after the other.
preview_lines <- function(x, n, width) {
  nodes <- table_block(x$nodes, "Nodes", n, width)
  edges <- table_block(x$edges, "Edges", n, width)
  if (max(nchar(nodes)) + 3L + max(nchar(edges)) > width) {
    half <- (width - 3L) %/% 2L
    nodes <- table_block(x$nodes, "Nodes", n, half)
    edges <- table_block(x$edges, "Edges", n, width - half - 3L)
    if (max(nchar(nodes)) + 3L + max(nchar(edges)) > width) {
      return(c(nodes, "", edges))
    }
  }
  k <- max(length(nodes), length(edges))
  left <- format(c(nodes, rep("", k - length(nodes))), width = max(nchar(nodes)))
  right <- c(edges, rep("", k - length(edges)))
  sub("\\s+$", "", paste(left, right, sep = "   "))
}

# The lines of one table's preview: a title with the dimensions, the
# column names, their types, the first `n` rows, and how many more rows
# and columns there are. Columns that do not fit in `width` are left out.
table_block <- function(df, title, n, width) {
  nr <- nrow(df)
  nc <- ncol(df)
  k <- min(nr, n)
  head <- sprintf("%s (%s x %s)", title, format(nr), format(nc))
  if (nc == 0L) {
    return(c(head, if (nr > k) sprintf("  ... %s more rows", format(nr - k))))
  }
  rows <- seq_len(k)
  num <- format(c("", "", as.character(rows)), width = nchar(format(k)))
  cols <- lapply(names(df), function(col) {
    v <- df[[col]]
    cells <- c(strtrim(col, 20L), column_type(v), format_cells(v, rows))
    format(cells, width = max(nchar(cells)))
  })
  widths <- vapply(cols, function(c) nchar(c[[1L]]), 0L)
  used <- nchar(num[[1L]])
  shown <- 0L
  for (w in widths) {
    if (used + 1L + w > width && shown > 0L) break
    used <- used + 1L + w
    shown <- shown + 1L
  }
  body <- num
  for (c in cols[seq_len(shown)]) body <- paste(body, c)
  body <- sub("\\s+$", "", body)
  more <- c(
    if (nr > k) sprintf("%s more rows", format(nr - k)),
    if (nc > shown) sprintf("%s more columns", format(nc - shown))
  )
  c(head, body, if (length(more)) sprintf("  ... %s", paste(more, collapse = ", ")))
}

# An abbreviation of a column's type, in the style of tibble.
column_type <- function(v) {
  if (inherits(v, "sfc")) return(sprintf("<%s>", geometry_kind(v)))
  if (is.factor(v)) return("<fct>")
  if (inherits(v, "Date")) return("<date>")
  if (inherits(v, "POSIXt")) return("<dttm>")
  if (inherits(v, "difftime")) return("<drtn>")
  if (inherits(v, "integer64")) return("<i64>")
  if (is.logical(v)) return("<lgl>")
  if (is.integer(v)) return("<int>")
  if (is.double(v)) return("<dbl>")
  if (is.character(v)) return("<chr>")
  if (is.list(v)) return("<list>")
  sprintf("<%s>", strtrim(class(v)[[1L]], 6L))
}

# The first `rows` values of a column as text of one width.
format_cells <- function(v, rows) {
  if (!length(rows)) return(character())
  if (inherits(v, "sfc")) return(rep(geometry_kind(v), length(rows)))
  v <- v[rows]
  if (is.list(v)) {
    return(vapply(v, function(e) sprintf("<%s[%d]>", class(e)[[1L]], length(e)), ""))
  }
  if (is.factor(v)) v <- as.character(v)
  out <- if (is.double(v) && !is.object(v)) {
    format(v, digits = 4L)
  } else {
    format(v)
  }
  out <- as.character(out)
  out[is.na(v)] <- "NA"
  strtrim(out, 20L)
}
