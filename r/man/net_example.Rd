\name{net_example}
\alias{net_example}
\alias{net_examples}
\title{Example networks}
\description{
Classic networks for teaching, testing, and reproducing published analyses,
returned as net_frames. \code{net_examples()} lists them.
}
\usage{
net_example(name = c("noordin", "noordin_affiliations", "jemaah_islamiyah",
                     "karate", "florentine", "southern_women", "krackhardt",
                     "sampson", "highland_tribes", "les_miserables", "miles",
                     "country_borders", "workplace_contacts"),
            relation = NULL, weighted = FALSE)

net_examples()
}
\arguments{
\item{name}{Which network; see Details.}

\item{relation}{For the multiplex examples (\code{"noordin"},
\code{"krackhardt"}, \code{"sampson"}), the relation(s) to keep (see
Details); \code{NULL} keeps all, with the \code{relation} role on the
edges.}

\item{weighted}{For networks with a tie-strength column, whether to use it as
the network's weights, with the meaning the data give them (strengths,
except the distances of \code{"miles"}). Default \code{FALSE}, so measures
treat the network as binary unless asked.}
}
\details{
\describe{
  \item{\code{"noordin"}}{The Noordin Top terrorist network (Roberts and
    Everton 2011): 79 people and ten undirected relations among them, coded
    from International Crisis Group (2006) by Naval Postgraduate School
    students under Nancy Roberts and Sean Everton and cleaned by Daniel
    Cunningham. Relations (\code{relation} column): \code{business_finance},
    \code{classmates}, \code{communications}, \code{friendship},
    \code{kinship}, \code{logistics}, \code{meetings}, \code{operations},
    \code{soulmates}, \code{training}. This is the network analyzed throughout
    Everton (2012) and Cunningham, Everton, and Murphy (2016). The copy here
    comes from the \href{https://github.com/schochastics/networkdata}{networkdata}
    package (\code{covert_34}); see \code{data-raw/noordin.R} in the source
    repository.}
  \item{\code{"noordin_affiliations"}}{A two-mode network of the same 79 people
    and 45 events (eight organizations, five operations, eleven training
    events, two financial meetings, seven logistic meetings, and twelve
    general meetings), as distributed in the
    \href{https://cran.r-project.org/package=manet}{manet} package. That copy
    does not record names, so people and events are labelled
    \code{person_01} \dots \code{person_79} and \code{event_01} \dots
    \code{event_45} in its order (the node column \code{mode} is
    \code{"person"} or \code{"event"}). The last five people, the "lone wolves", attended none of
    the events.}
  \item{\code{"jemaah_islamiyah"}}{The Jemaah Islamiyah cell responsible for
    the 2002 Bali bombings (Koschade 2006): 17 people, their roles
    (\code{role} node column), and undirected ties from Koschade's Table 1,
    with the interaction strength from his Table 2 (1 = weakest, 5 =
    strongest) in the \code{strength} column. One cell of Table 2 is blank and
    is filled from its symmetric counterpart. As used in the
    \href{https://github.com/knapply/homophily}{homophily} package.}
  \item{\code{"karate"}}{Zachary's (1977) karate club: 34 members, 78 ties, the
    faction each member joined after the split (\code{faction}), and the
    number of contexts in which each pair interacted (\code{contexts}).}
  \item{\code{"florentine"}}{Marriage ties among 15 Florentine families
    (Padgett and Ansell 1993; Breiger and Pattison 1986), as distributed with
    \href{https://networkx.org/}{NetworkX} (the isolated Pucci family is not
    included).}
  \item{\code{"southern_women"}}{Davis, Gardner, and Gardner's (1941) Southern
    Women: 18 women and the 14 social events they attended, as a two-mode
    network.}
  \item{\code{"krackhardt"}}{Krackhardt's (1987) 21 high-tech managers:
    three directed relations (\code{relation}: \code{advice},
    \code{friendship}, \code{reports_to}) and the node attributes
    \code{age}, \code{tenure}, \code{level}, and \code{department}, as in
    Wasserman and Faust (1994, Appendix B). Identifiers are the manager
    numbers (column \code{name}).}
  \item{\code{"sampson"}}{Sampson's (1968) 18 novices: eight directed,
    ranked relations (\code{relation}), each positive or negative
    (\code{sign}, the \code{sign} role), with liking gathered at three
    times and the other relations at the third (\code{time}, the
    \code{wave} role, so the network is temporal); \code{rank} is 3 for a
    first choice and 1 for a third (the weight with \code{weighted =
    TRUE}). Analysed by Breiger, Boorman, and Arabie (1975).}
  \item{\code{"highland_tribes"}}{Read's (1954) 16 Gahuku-Gama sub-tribes
    with their alliances (\code{sign} 1) and antagonisms (\code{sign}
    -1), the signed graph of Hage and Harary (1983).}
  \item{\code{"les_miserables"}}{The 80 characters of \emph{Les
    Miserables} and the number of times each pair meets
    (\code{encounters}, a strength), from Knuth's (1993) Stanford
    GraphBase; three characters meet no one.}
  \item{\code{"miles"}}{128 North American cities from the Stanford
    GraphBase (Knuth 1993) with the highway mileage between every pair
    (\code{miles}, a distance) and their coordinates (\code{longitude},
    \code{latitude}, recorded as the coordinate roles), so the network is
    spatial. The table is complete; keep short edges for a sparse network.}
  \item{\code{"country_borders"}}{The 177 country units of Natural
    Earth's 1:110m map and the pairs that share a land border, with the
    border length (\code{border_km}, a strength), the label point of each
    unit as coordinates, and the border itself as WKT text in the
    \code{geometry} edge column (convert with \code{sf::st_as_sfc()}).
    Identifiers are Natural Earth's \code{code}s; 21 units are isolates.}
  \item{\code{"workplace_contacts"}}{Face-to-face contacts among employees
    of an office building recorded by wearable sensors (Genois et al.
    2015): one row per 20-second contact (\code{time}, seconds, the
    \code{start} role, so the network is a contact sequence), with each
    person's \code{department}. Its files are built by the package's data
    workflow and may be absent from a development installation;
    \code{net_examples()} says whether they are installed.}
}
The newer data sets' directories under \code{system.file("extdata", package =
"knappnet")} have a README with the source, licence, and meaning of
every column.
}
\value{
\code{net_example()}: a net_frame (two-mode for
\code{"noordin_affiliations"} and \code{"southern_women"}).
\code{net_examples()}: a data frame describing each network, with
\code{available} saying whether its files are installed.
}
\references{
Roberts, N., & Everton, S. F. (2011). Roberts and Everton Terrorist Data:
Noordin Top Terrorist Network (Subset). Association of Religion Data
Archives.

International Crisis Group (2006). Terrorism in Indonesia: Noordin's Networks.
Asia Report 114.

Everton, S. F. (2012). \emph{Disrupting Dark Networks}. Cambridge University
Press. \doi{10.1017/CBO9781139136877}

Cunningham, D., Everton, S. F., & Murphy, P. (2016). \emph{Understanding
Dark Networks: A Strategic Framework for the Use of Social Network Analysis}.
Rowman & Littlefield.

Koschade, S. (2006). A social network analysis of Jemaah Islamiyah: The
applications to counterterrorism and intelligence. \emph{Studies in Conflict &
Terrorism}, 29(6), 559-575. \doi{10.1080/10576100600798418}

Zachary, W. W. (1977). An information flow model for conflict and fission in
small groups. \emph{Journal of Anthropological Research}, 33(4), 452-473.
\doi{10.1086/jar.33.4.3629752}

Padgett, J. F., & Ansell, C. K. (1993). Robust action and the rise of the
Medici, 1400-1434. \emph{American Journal of Sociology}, 98(6), 1259-1319.
\doi{10.1086/230190}

Breiger, R. L., & Pattison, P. E. (1986). Cumulated social roles: The
duality of persons and their algebras. \emph{Social Networks}, 8(3), 215-256.
\doi{10.1016/0378-8733(86)90006-7}

Davis, A., Gardner, B. B., & Gardner, M. R. (1941). \emph{Deep South: A
Social Anthropological Study of Caste and Class}. University of Chicago Press.

Krackhardt, D. (1987). Cognitive social structures. \emph{Social Networks},
9(2), 109-134. \doi{10.1016/0378-8733(87)90009-8}

Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods
and Applications}. Cambridge University Press.
\doi{10.1017/CBO9780511815478}

Sampson, S. F. (1968). \emph{A Novitiate in a Period of Change: An
Experimental and Case Study of Social Relationships}. PhD dissertation,
Cornell University.

Breiger, R. L., Boorman, S. A., & Arabie, P. (1975). An algorithm for
clustering relational data with applications to social network analysis
and comparison with multidimensional scaling. \emph{Journal of
Mathematical Psychology}, 12(3), 328-383.
\doi{10.1016/0022-2496(75)90028-0}

Read, K. E. (1954). Cultures of the Central Highlands, New Guinea.
\emph{Southwestern Journal of Anthropology}, 10(1), 1-43.
\doi{10.1086/soutjanth.10.1.3629074}

Hage, P., & Harary, F. (1983). \emph{Structural Models in Anthropology}.
Cambridge University Press. \doi{10.1017/CBO9780511659843}

Knuth, D. E. (1993). \emph{The Stanford GraphBase: A Platform for
Combinatorial Computing}. Addison-Wesley.

Natural Earth (2022). Admin 0 - Countries, 1:110m cultural vectors,
version 5.2.0. \url{https://www.naturalearthdata.com/}

Genois, M., Vestergaard, C. L., Fournet, J., Panisson, A., Bonmarin, I., &
Barrat, A. (2015). Data on face-to-face contacts in an office building
suggest a low-cost vaccination strategy based on community linkers.
\emph{Network Science}, 3(3), 326-347. \doi{10.1017/nws.2015.10}
}
\examples{
net_examples()

noordin <- net_example("noordin")
noordin
table(noordin$edges$relation)

# One relation, as in the chapter examples of Everton (2012).
trust <- net_example("noordin", relation = c("classmates", "friendship",
                                             "kinship", "soulmates"))
head(sort(net_betweenness(trust), decreasing = TRUE))

# Three relations among Krackhardt's managers, coloured by relation.
krack <- net_example("krackhardt")
plot(krack, seed = 1)

# Alliances (solid) and antagonisms (dashed) among the highland tribes.
tribes <- net_example("highland_tribes")
plot(tribes, seed = 1)

# Cities on their coordinates, joined when less than 150 miles apart.
miles <- net_example("miles", weighted = TRUE)
near <- net_frame(miles$edges[miles$edges$miles < 150, ], miles$nodes,
                  directed = FALSE, source = "from", target = "to", id = "name",
                  weight = "miles", weight_type = "distance",
                  coords = c("longitude", "latitude"))
plot(near, node_label = FALSE)
}
