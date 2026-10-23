#' Example networks
#'
#' Classic networks for teaching, testing, and reproducing published
#' analyses, returned as netframes. `example_networks()` lists them.
#'
#' @param name Which network; see Details.
#' @param relation For `"noordin"`, the relation(s) to keep (see Details);
#'   `NULL` keeps all, with a `relation` column on the edges.
#' @param weighted For networks with a tie-strength column, whether to use
#'   it as the network's weights (as strengths). Default `FALSE`, so
#'   measures treat the network as binary unless asked.
#'
#' @details
#' \describe{
#'   \item{`"noordin"`}{The Noordin Top terrorist network (Roberts and
#'     Everton 2011): 79 people and ten undirected relations among them,
#'     coded from International Crisis Group (2006) by Naval Postgraduate
#'     School students under Nancy Roberts and Sean Everton and cleaned by
#'     Daniel Cunningham. Relations (`relation` column): `business_finance`,
#'     `classmates`, `communications`, `friendship`, `kinship`,
#'     `logistics`, `meetings`, `operations`, `soulmates`, `training`. This
#'     is the network analyzed throughout Everton (2012) and Cunningham,
#'     Everton, and Murphy (2016). The copy here comes from the
#'     \href{https://github.com/schochastics/networkdata}{networkdata}
#'     package (`covert_34`); see `data-raw/noordin.R` in the source
#'     repository.}
#'   \item{`"noordin_affiliations"`}{A two-mode network of the same 79
#'     people and 45 events (eight organizations, five operations, eleven
#'     training events, two financial meetings, seven logistic meetings, and
#'     twelve general meetings), as distributed in the
#'     \href{https://cran.r-project.org/package=manet}{manet} package. That
#'     copy does not record names, so people and events are labelled
#'     `person_01` ... `person_79` and `event_01` ... `event_45` in its
#'     order (the node column `event` is `TRUE` for events). The last five
#'     people, the "lone wolves", attended none of the events.}
#'   \item{`"jemaah_islamiyah"`}{The Jemaah Islamiyah cell responsible for
#'     the 2002 Bali bombings (Koschade 2006): 17 people, their roles
#'     (`role` node column), and undirected ties from Koschade's Table 1,
#'     with the interaction strength from his Table 2 (1 = weakest, 5 =
#'     strongest) in the `strength` column. One cell of Table 2 is blank and
#'     is filled from its symmetric counterpart. As used in the
#'     \href{https://github.com/knapply/homophily}{homophily} package.}
#'   \item{`"karate"`}{Zachary's (1977) karate club: 34 members, 78 ties,
#'     the faction each member joined after the split (`faction`), and the
#'     number of contexts in which each pair interacted (`contexts`).}
#'   \item{`"florentine"`}{Marriage ties among 15 Florentine families
#'     (Padgett and Ansell 1993; Breiger and Pattison 1986), as distributed
#'     with NetworkX (the isolated Pucci family is not included).}
#'   \item{`"southern_women"`}{Davis, Gardner, and Gardner's (1941)
#'     Southern Women: 18 women and the 14 social events they attended, as
#'     a two-mode network.}
#' }
#'
#' @references
#' Roberts, N., and Everton, S. F. (2011). Roberts and Everton Terrorist
#' Data: Noordin Top Terrorist Network (Subset). Association of Religion
#' Data Archives.
#'
#' International Crisis Group (2006). Terrorism in Indonesia: Noordin's
#' Networks. Asia Report 114.
#'
#' Everton, S. F. (2012). *Disrupting Dark Networks*. Cambridge University
#' Press. \doi{10.1017/CBO9781139136877}
#'
#' Cunningham, D., Everton, S. F., and Murphy, P. (2016). *Understanding
#' Dark Networks: A Strategic Framework for the Use of Social Network
#' Analysis*. Rowman & Littlefield.
#'
#' Koschade, S. (2006). A social network analysis of Jemaah Islamiyah: The
#' applications to counterterrorism and intelligence. *Studies in Conflict
#' & Terrorism*, 29(6), 559-575.
#'
#' Zachary, W. W. (1977). An information flow model for conflict and fission
#' in small groups. *Journal of Anthropological Research*, 33(4), 452-473.
#' \doi{10.1086/jar.33.4.3629752}
#'
#' Padgett, J. F., and Ansell, C. K. (1993). Robust action and the rise of
#' the Medici, 1400-1434. *American Journal of Sociology*, 98(6),
#' 1259-1319. \doi{10.1086/230190}
#'
#' Breiger, R. L., and Pattison, P. E. (1986). Cumulated social roles: The
#' duality of persons and their algebras. *Social Networks*, 8(3), 215-256.
#' \doi{10.1016/0378-8733(86)90006-7}
#'
#' Davis, A., Gardner, B. B., and Gardner, M. R. (1941). *Deep South: A
#' Social Anthropological Study of Caste and Class*. University of Chicago
#' Press.
#'
#' @return `example_network()`: a netframe (two-mode for
#'   `"noordin_affiliations"` and `"southern_women"`).
#'   `example_networks()`: a data frame describing each network.
#' @export
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
      netframe(edges, read("noordin", "nodes.csv"), directed = FALSE)
    },
    noordin_affiliations = twomode(read("noordin", "affiliations.csv"),
                                   read("noordin", "affiliation_nodes.csv"), type = "event"),
    jemaah_islamiyah = netframe(read("jemaah_islamiyah", "edges.csv"),
                                read("jemaah_islamiyah", "nodes.csv"), directed = FALSE,
                                weight = if (weighted) "strength",
                                weight_type = "strength"),
    karate = netframe(read("karate", "edges.csv"), read("karate", "nodes.csv"),
                      directed = FALSE, weight = if (weighted) "contexts",
                      weight_type = "strength"),
    florentine = netframe(read("florentine", "edges.csv"), directed = FALSE),
    southern_women = twomode(read("southern_women", "edges.csv"))
  )
}

#' @rdname example_network
#' @export
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
