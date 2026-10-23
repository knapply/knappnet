\name{net_example}
\alias{net_example}
\alias{net_examples}
\title{Example networks}
\description{
Classic networks for teaching, testing, and reproducing published analyses,
returned as netframes. \code{net_examples()} lists them.
}
\usage{
net_example(name = c("noordin", "noordin_affiliations", "jemaah_islamiyah",
                         "karate", "florentine", "southern_women"),
                relation = NULL, weighted = FALSE)

net_examples()
}
\arguments{
\item{name}{Which network; see Details.}

\item{relation}{For \code{"noordin"}, the relation(s) to keep (see Details);
\code{NULL} keeps all, with a \code{relation} column on the edges.}

\item{weighted}{For networks with a tie-strength column, whether to use it as
the network's weights (as strengths). Default \code{FALSE}, so measures treat
the network as binary unless asked.}
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
    \code{event_45} in its order (the node column \code{event} is \code{TRUE}
    for events). The last five people, the "lone wolves", attended none of
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
}
}
\value{
\code{net_example()}: a netframe (two-mode for
\code{"noordin_affiliations"} and \code{"southern_women"}).
\code{net_examples()}: a data frame describing each network.
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

ji <- net_example("jemaah_islamiyah", weighted = TRUE)
ji$nodes$role
}
