# Example networks shipped in inst/extdata, read on demand.

# Example networks (man/example_network.Rd).
example_network <- function(name = c("noordin", "noordin_affiliations", "jemaah_islamiyah",
                                     "karate", "florentine", "southern_women"),
                            relation = NULL, weighted = FALSE) {
  name <- match.arg(name)
  weighted <- check_flag(weighted, "weighted")
  read <- function(dir, file) {
    path <- system.file("extdata", dir, file, package = "knappnet")
    if (!nzchar(path)) stop(sprintf("Example data file `%s/%s` is missing.", dir, file), call. = FALSE)
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8")
  }
  switch(name,
    noordin = {
      edges <- read("noordin", "relations.csv")
      if (!is.null(relation)) {
        bad <- setdiff(relation, unique(edges$relation))
        if (length(bad)) {
          stop(sprintf("Unknown `relation` %s; choose from %s.",
                       paste0("\"", bad, "\"", collapse = ", "),
                       paste0("\"", sort(unique(edges$relation)), "\"", collapse = ", ")),
               call. = FALSE)
        }
        edges <- edges[edges$relation %in% relation, , drop = FALSE]
        rownames(edges) <- NULL
      }
      net_frame(edges, read("noordin", "nodes.csv"), directed = FALSE, relation = "relation",
                quiet = TRUE)
    },
    noordin_affiliations = net_twomode(read("noordin", "affiliations.csv"),
                                       read("noordin", "affiliation_nodes.csv"),
                                       source = "person", target = "event", quiet = TRUE),
    jemaah_islamiyah = net_frame(read("jemaah_islamiyah", "edges.csv"),
                                 read("jemaah_islamiyah", "nodes.csv"), directed = FALSE,
                                 weight = if (weighted) "strength",
                                 weight_type = if (weighted) "strength", quiet = TRUE),
    karate = net_frame(read("karate", "edges.csv"), read("karate", "nodes.csv"),
                       directed = FALSE, weight = if (weighted) "contexts",
                       weight_type = if (weighted) "strength", quiet = TRUE),
    florentine = net_frame(read("florentine", "edges.csv"), directed = FALSE, quiet = TRUE),
    southern_women = net_twomode(read("southern_women", "edges.csv"), source = "woman",
                                 target = "event", quiet = TRUE)
  )
}

# See man/example_network.Rd.
example_networks <- function() {
  data.frame(
    name = c("noordin", "noordin_affiliations", "jemaah_islamiyah", "karate", "florentine",
             "southern_women"),
    nodes = c(79L, 124L, 17L, 34L, 15L, 32L),
    edges = c(1014L, 275L, 63L, 78L, 20L, 89L),
    type = c("one-mode, multiplex (10 relations)", "two-mode", "one-mode, valued",
             "one-mode, valued", "one-mode", "two-mode"),
    source = c("Roberts & Everton (2011)", "Everton (2012), via manet", "Koschade (2006)",
               "Zachary (1977)", "Padgett & Ansell (1993)", "Davis, Gardner & Gardner (1941)"),
    stringsAsFactors = FALSE
  )
}
