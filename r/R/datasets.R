# Example networks shipped in inst/extdata (see inst/extdata/DATASETS.md),
# read on demand.

example_names <- c("noordin", "noordin_affiliations", "jemaah_islamiyah", "karate",
                   "florentine", "southern_women", "krackhardt", "sampson",
                   "highland_tribes", "les_miserables", "miles", "country_borders",
                   "workplace_contacts")

# Example networks (man/net_example.Rd).
net_example <- function(name = example_names, relation = NULL, weighted = FALSE) {
  name <- match.arg(name)
  weighted <- check_flag(weighted, "weighted")
  if (!is.null(relation) && !name %in% c("noordin", "krackhardt", "sampson")) {
    stop(sprintf("`relation` applies to the multiplex examples, not to \"%s\".", name),
         call. = FALSE)
  }
  w <- function(col) if (weighted) col
  wt <- function(type = "strength") if (weighted) type
  switch(name,
    noordin = net_frame(pick_relations(example_csv("noordin", "relations.csv"), relation),
                        example_csv("noordin", "nodes.csv"), directed = FALSE,
                        relation = "relation", quiet = TRUE),
    noordin_affiliations = net_twomode(example_csv("noordin", "affiliations.csv"),
                                       example_csv("noordin", "affiliation_nodes.csv"),
                                       source = "person", target = "event", quiet = TRUE),
    jemaah_islamiyah = net_frame(example_csv("jemaah_islamiyah", "edges.csv"),
                                 example_csv("jemaah_islamiyah", "nodes.csv"), directed = FALSE,
                                 weight = w("strength"), weight_type = wt(), quiet = TRUE),
    karate = net_frame(example_csv("karate", "edges.csv"), example_csv("karate", "nodes.csv"),
                       directed = FALSE, weight = w("contexts"), weight_type = wt(),
                       quiet = TRUE),
    florentine = net_frame(example_csv("florentine", "edges.csv"), directed = FALSE,
                           quiet = TRUE),
    southern_women = net_twomode(example_csv("southern_women", "edges.csv"), source = "woman",
                                 target = "event", quiet = TRUE),
    krackhardt = net_frame(pick_relations(example_csv("krackhardt", "edges.csv"), relation),
                           example_csv("krackhardt", "nodes.csv"), source = "from",
                           target = "to", id = "name", relation = "relation", quiet = TRUE),
    sampson = net_frame(pick_relations(example_csv("sampson", "edges.csv"), relation),
                        example_csv("sampson", "nodes.csv"), source = "from", target = "to",
                        id = "name", relation = "relation", sign = "sign", wave = "time",
                        weight = w("rank"), weight_type = wt(), quiet = TRUE),
    highland_tribes = net_frame(example_csv("highland_tribes", "edges.csv"),
                                example_csv("highland_tribes", "nodes.csv"), directed = FALSE,
                                source = "from", target = "to", id = "name", sign = "sign",
                                quiet = TRUE),
    les_miserables = net_frame(example_csv("les_miserables", "edges.csv"),
                               example_csv("les_miserables", "nodes.csv"), directed = FALSE,
                               source = "from", target = "to", id = "name",
                               weight = w("encounters"), weight_type = wt(), quiet = TRUE),
    miles = net_frame(example_csv("miles", "edges.csv"), example_csv("miles", "nodes.csv"),
                      directed = FALSE, source = "from", target = "to", id = "name",
                      weight = w("miles"), weight_type = wt("distance"),
                      coords = c("longitude", "latitude"), crs = "lonlat", quiet = TRUE),
    country_borders = net_frame(example_csv("country_borders", "edges.csv"),
                                example_csv("country_borders", "nodes.csv"), directed = FALSE,
                                source = "from", target = "to", id = "code",
                                weight = w("border_km"), weight_type = wt(),
                                coords = c("longitude", "latitude"), crs = "lonlat",
                                quiet = TRUE),
    workplace_contacts = net_frame(example_csv("workplace_contacts", "contacts.csv"),
                                   example_csv("workplace_contacts", "nodes.csv"),
                                   directed = FALSE, source = "from", target = "to",
                                   id = "name", start = "time", quiet = TRUE)
  )
}

# A CSV table of an example data set.
example_csv <- function(dir, file) {
  path <- system.file("extdata", dir, file, package = "knappnet")
  if (!nzchar(path)) {
    stop(sprintf(paste0("Example data file `%s/%s` is not in this installation of knappnet ",
                        "(see inst/extdata/DATASETS.md)."), dir, file), call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
}

# Whether an example data set's files are installed.
example_available <- function(name) {
  file <- switch(name,
    noordin = , noordin_affiliations = c("noordin", "relations.csv"),
    southern_women = c("southern_women", "edges.csv"),
    workplace_contacts = c("workplace_contacts", "contacts.csv"),
    c(name, "edges.csv")
  )
  nzchar(system.file("extdata", file[[1L]], file[[2L]], package = "knappnet"))
}

# The rows of the chosen relations.
pick_relations <- function(edges, relation) {
  if (is.null(relation)) return(edges)
  bad <- setdiff(relation, unique(edges$relation))
  if (length(bad)) {
    stop(sprintf("Unknown `relation` %s; choose from %s.",
                 paste0("\"", bad, "\"", collapse = ", "),
                 paste0("\"", sort(unique(edges$relation)), "\"", collapse = ", ")),
         call. = FALSE)
  }
  edges <- edges[edges$relation %in% relation, , drop = FALSE]
  rownames(edges) <- NULL
  edges
}

# See man/net_example.Rd.
net_examples <- function() {
  out <- data.frame(
    name = example_names,
    nodes = c(79L, 124L, 17L, 34L, 15L, 32L, 21L, 18L, 16L, 80L, 128L, 177L, NA),
    edges = c(1014L, 275L, 63L, 78L, 20L, 89L, 312L, 510L, 58L, 254L, 8128L, 313L, NA),
    type = c("one-mode, undirected, multiplex (10 relations)", "two-mode",
             "one-mode, undirected, valued", "one-mode, undirected, valued",
             "one-mode, undirected", "two-mode",
             "one-mode, directed, multiplex (3 relations)",
             "one-mode, directed, signed, multiplex (8 relations), 3 waves",
             "one-mode, undirected, signed", "one-mode, undirected, valued, isolates",
             "one-mode, undirected, distances, spatial (points)",
             "one-mode, undirected, valued, spatial (points), isolates",
             "one-mode, undirected, temporal (contacts)"),
    source = c("Roberts & Everton (2011)", "Everton (2012), via manet", "Koschade (2006)",
               "Zachary (1977)", "Padgett & Ansell (1993)", "Davis, Gardner & Gardner (1941)",
               "Krackhardt (1987)", "Sampson (1968)", "Read (1954)", "Knuth (1993)",
               "Knuth (1993)", "Natural Earth (2022)", "Genois et al. (2015)"),
    stringsAsFactors = FALSE
  )
  out$available <- vapply(out$name, example_available, TRUE, USE.NAMES = FALSE)
  out
}
