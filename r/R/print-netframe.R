# The print card of a netframe: a labelled card filled from the engine's
# summary (size, direction, modes, weights, relations, time, space,
# loops and parallel edges, components, isolates, density), followed by
# previews of the node and edge tables, side by side when the console is
# wide enough. Base R only: format(), formatC(), and strtrim(), with no
# characters outside ASCII, so the output is the same on every console.
# The subclasses add their own lines through print_card().

# See man/net_frame.Rd.
print.netframe <- function(x, ..., n = 6L, width = getOption("width")) {
  print_card(x, character(), n, width)
}

# See man/net_temporal.Rd.
print.temporal_netframe <- function(x, ..., n = 6L, width = getOption("width")) {
  print_card(x, temporal_card_lines(x), n, width)
}

# See man/net_spatial.Rd.
print.spatial_netframe <- function(x, ..., n = 6L, width = getOption("width")) {
  print_card(x, spatial_card_lines(x), n, width)
}

# See man/net_temporal.Rd.
print.spatiotemporal_netframe <- function(x, ..., n = 6L, width = getOption("width")) {
  print_card(x, c(temporal_card_lines(x), spatial_card_lines(x)), n, width)
}

# Prints the card, the subclass lines, and the table previews.
print_card <- function(x, extra, n, width) {
  check_netframe(x)
  n <- check_count(n, "n")
  width <- check_count(width, "width", min = 20L)
  lines <- c(card_lines(x), extra, "", preview_lines(x, n, width - 2L))
  cat(paste0("#", ifelse(nzchar(lines), " ", ""), lines), sep = "\n")
  invisible(x)
}

# ---- The card ---------------------------------------------------------------------

# The lines of the card, without the "# " prefix.
card_lines <- function(x) {
  s <- netframe_summary(x)
  pair <- function(left, right) {
    pad <- max(28L, nchar(left) + 4L)
    paste0(format(left, width = pad), right)
  }
  c(
    sprintf("%s: %s, %s", class(x)[[1L]], plural(as.integer(s$nodes), "node"),
            plural(as.integer(s$edges), "edge")),
    pair(sprintf("  directed: %s", if (s$directed) "yes" else "no"),
         sprintf("mode: %s", mode_text(x, s))),
    sprintf("  weights:  %s", weight_text(s)),
    sprintf("  relations: %s", relation_text(s)),
    if (s$signed) sprintf("  sign:     %s (%s of %s negative)", x$roles$sign,
                          format(s$negative), format(s$edges)),
    pair(sprintf("  time:     %s", time_text(x, s)), sprintf("space: %s", space_text(x, s))),
    sprintf("  simple:   %s", simple_text(x, s)),
    if (s$missing_nodes > 0 || s$missing_edges > 0) {
      sprintf("  missing:  %s and %s marked unobserved", plural(as.integer(s$missing_nodes), "node"),
              plural(as.integer(s$missing_edges), "edge"))
    },
    sprintf("  components: %s   isolates: %s   density: %s", format(s$components),
            format(s$isolates), number_text3(s$density))
  )
}

number_text3 <- function(v) {
  if (is.na(v)) "NaN" else formatC(v, format = "f", digits = 3)
}

label_counts <- function(labels, counts) {
  paste(sprintf("%s: %s", labels, format(as.integer(counts), trim = TRUE)), collapse = ", ")
}

mode_text <- function(x, s) {
  k <- length(s$modes)
  if (k == 0L) return("one-mode")
  counts <- label_counts(s$modes, s$mode_counts)
  if (isTRUE(s$between_modes_only)) {
    if (k == 2L) return(sprintf("two-mode (%s)", counts))
    return(sprintf("%d-mode (%s)", k, counts))
  }
  sprintf("multimode, ties within modes allowed (%s)", counts)
}

weight_text <- function(s) {
  if (!s$weighted) return("none (unweighted)")
  range <- if (is.na(s$weight_min)) "no edges" else {
    sprintf("%s to %s", format(s$weight_min, digits = 4), format(s$weight_max, digits = 4))
  }
  sprintf("%s (%s, %s)", s$weight, s$weight_type, range)
}

relation_text <- function(s) {
  if (!length(s$relations)) return("none")
  sprintf("%s; multiplex", label_counts(s$relations, s$relation_counts))
}

time_text <- function(x, s) {
  if (s$time == "static") return("static")
  r <- x$roles
  range <- if (is.na(s$time_start)) "no edges" else {
    sprintf("%s to %s", format(from_seconds(s$time_start, x)),
            format(from_seconds(s$time_end, x)))
  }
  kind <- switch(s$time,
    contacts = sprintf("contacts in `%s`", r$start),
    spells = sprintf("edge spells `%s` to `%s`", r$start,
                     if (is.null(r$end)) r$duration else r$end),
    waves = sprintf("waves in `%s`", r$wave)
  )
  sprintf("%s (%s)%s", kind, range, if (s$node_spells) ", node spells" else "")
}

space_text <- function(x, s) {
  r <- x$roles
  if (s$node_geometry || s$edge_geometry) {
    kinds <- c(
      if (s$node_geometry) sprintf("node %s", geometry_kind(x$nodes[[r$node_geometry]])),
      if (s$edge_geometry) sprintf("edge %s", geometry_kind(x$edges[[r$edge_geometry]]))
    )
    crs <- tryCatch(crs_label(r$crs), error = function(e) "a CRS")
    return(sprintf("%s geometry, %s", paste(kinds, collapse = " and "), crs))
  }
  if (!is.null(r$x) && !is.null(r$y)) {
    return(sprintf("coordinates `%s`, `%s` (%s)", r$x, r$y,
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

simple_text <- function(x, s) {
  declared <- function(flag) if (isFALSE(flag)) " (declared)" else ""
  loops <- if (s$loops == 0) {
    paste0("no loops", declared(x$declared$loops))
  } else {
    plural(as.integer(s$loops), "self-loop")
  }
  parallel <- if (s$parallel == 0) {
    paste0("no parallel edges", declared(x$declared$multiple))
  } else {
    plural(as.integer(s$parallel), "parallel edge")
  }
  paste(loops, parallel, sep = ", ")
}

# The lines a temporal network adds: how times are stored and node spells.
temporal_card_lines <- function(x) {
  r <- x$roles
  kind <- switch(if (is.null(r$time_class)) "numeric" else r$time_class,
    numeric = if (!is.null(r$wave)) "whole-numbered waves" else "numbers",
    Date = "dates",
    POSIXct = sprintf("date-times (%s)", if (is.null(r$tz) || !nzchar(r$tz)) "local time" else r$tz)
  )
  spells <- if (!is.null(r$node_start) || !is.null(r$node_end)) {
    sprintf("; node spells in %s", paste0("`", c(r$node_start, r$node_end), "`", collapse = " to "))
  } else {
    ""
  }
  sprintf("  temporal: times are %s%s", kind, spells)
}

# The lines a spatial network adds: how distances are measured.
spatial_card_lines <- function(x) {
  r <- x$roles
  kind <- tryCatch(crs_kind(x), error = function(e) "planar")
  distances <- if (identical(kind, "lonlat")) {
    sprintf("great-circle distances, sphere radius %s m",
            format(if (is.null(r$radius)) 6371008.7714 else r$radius, digits = 8))
  } else {
    "planar distances in the coordinate units"
  }
  points <- if (!is.null(r$node_geometry)) {
    sprintf("; representative points by %s",
            if (identical(r$representative, "centroid")) "centroid" else "point on surface")
  } else {
    ""
  }
  sprintf("  spatial:  %s%s", distances, points)
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
