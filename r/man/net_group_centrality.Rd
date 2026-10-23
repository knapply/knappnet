\name{net_group_centrality}
\alias{net_group_centrality}
\title{Group centrality}
\description{
Everett and Borgatti's (1999) extension of degree, closeness, and
betweenness from single nodes to groups: how many outsiders a group is tied
to, how close it is to them, and how much of the traffic between outsiders
it can intercept. Applied to the members of a cell, a faction, or a set of
key players (\code{\link{net_key_players}()}), these measure the group's
position as a whole.
}
\usage{
net_group_centrality(x, group, measure = c("degree", "closeness", "betweenness"),
                 mode = c("out", "in", "all"), normalized = TRUE,
                 weighted = NULL, weights_as = NULL, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{group}{Identifiers of the group's members (values of the identifier
column of \code{x$nodes}), at least one.}
\item{measure}{\code{"degree"}, \code{"closeness"}, or
\code{"betweenness"}.}
\item{mode}{For directed networks. Degree: \code{"out"} counts
non-members the group sends ties to, \code{"in"} those sending ties to the
group, \code{"all"} either. Closeness: \code{"out"} uses distances from the
group, \code{"in"} distances to the group (NetworkX's definition), and
\code{"all"} ignores direction. Ignored by betweenness and for undirected
networks.}
\item{normalized}{Whether to normalize; see Details.}
\item{weighted}{Whether closeness and betweenness use edge weights as path
lengths. \code{NULL} (the default) uses them when the network has a weight
column. Degree is defined on binary ties and refuses \code{TRUE}.}
\item{weights_as}{Override the network's \code{weight_type}: \code{"distance"}
or \code{"strength"}.}
\item{threads}{Number of threads for betweenness, a whole number; \code{0}
uses all available cores, or 2 during \verb{R CMD check} (see \sQuote{Threads}
in \code{\link{net_frame}}).}
}
\details{
For a group \eqn{C} of \eqn{k} nodes in a network of \eqn{n}:
\describe{
  \item{\code{"degree"}}{The number of non-members tied to at least one
    member; normalized, divided by \eqn{n - k}.}
  \item{\code{"closeness"}}{From the distance \eqn{d(C, v)} between the
    group and each non-member \eqn{v}, the length of a shortest path from
    its nearest member: \eqn{1 / \sum_v d(C, v)}, or normalized
    \eqn{(n - k) / \sum_v d(C, v)}. Non-members the group cannot reach are
    left out of the sum, and the result is 0 when it reaches none; these
    are NetworkX's rules.}
  \item{\code{"betweenness"}}{\eqn{\sum_{s, t \notin C} \sigma_{st}(C) / \sigma_{st}},
    where \eqn{\sigma_{st}} is the number of shortest paths from \eqn{s} to
    \eqn{t} and \eqn{\sigma_{st}(C)} the number passing through at least one
    member. Pairs are ordered in directed networks and unordered in
    undirected ones. Normalized, divided by the number of such pairs,
    \eqn{(n-k)(n-k-1)} or half that; \code{NaN} with fewer than two
    non-members.}
}

Conventions: self-loops are ignored. Parallel edges count as distinct
shortest paths in betweenness, as in \code{\link{net_betweenness}()}, and do
not change degree or closeness. Paths use weights as for
\code{\link{net_betweenness}()}: distances, or strengths converted to lengths
\eqn{1/w}, which must be positive; in the general weighted case two path
lengths are equal when they agree to a relative tolerance of
\eqn{10^{-10}}. Betweenness needs one shortest-path search from every
non-member, run in parallel, counting all shortest paths and those that
avoid the group (Brandes 2008).

Validation: all three measures match NetworkX's
\code{group_degree_centrality()} (and its in- and out-degree versions),
\code{group_closeness_centrality()}, and
\code{group_betweenness_centrality()}, and a brute-force evaluation of the
definitions with exact path counts. NetworkX counts parallel edges once in
betweenness, and with non-integer weights it can lose shortest paths whose
floating-point lengths differ in the last bits; the exact evaluation
agrees with this function in both cases.
}
\value{
A number.
}
\references{
Everett, M. G., & Borgatti, S. P. (1999). The centrality of groups and
classes. \emph{Journal of Mathematical Sociology}, 23(3), 181-201.
\doi{10.1080/0022250X.1999.9990219}

Brandes, U. (2008). On variants of shortest-path betweenness centrality
and their generic computation. \emph{Social Networks}, 30(2), 136-145.
\doi{10.1016/j.socnet.2007.11.001}

Everton, S. F. (2012). \emph{Disrupting Dark Networks}. Cambridge
University Press. \doi{10.1017/CBO9781139136877}
}
\seealso{
\code{\link{net_key_players}}, \code{\link{net_flow_betweenness}}
}
\examples{
# The Jemaah Islamiyah cell behind the 2002 Bali bombings (Koschade 2006):
# how central is each role group as a whole?
ji <- net_example("jemaah_islamiyah")
roles <- split(ji$nodes$name, ji$nodes$role)
for (measure in c("degree", "closeness", "betweenness")) {
  print(sapply(roles, function(members) net_group_centrality(ji, members, measure, threads = 2)))
}

# In a directed network, closeness can be measured from the group or to
# it. A hypothetical chain of command: orders flow from the leaders down to
# two operational cells.
orders <- data.frame(
  source = c("boss", "boss", "lt1", "lt1", "lt2", "lt2", "op3", "lt2"),
  target   = c("lt1", "lt2", "op1", "op2", "op3", "op4", "op5", "lt1")
)
g <- net_frame(orders)
leaders <- c("boss", "lt1")
# Closeness from the group, along the direction orders travel.
net_group_centrality(g, leaders, "closeness", mode = "out")
# Closeness to the group, as NetworkX's group_closeness_centrality()
# measures it. Only lt2 can send anything back up the chain; the nodes that
# cannot reach the leaders are left out of the sum, so the score is high.
net_group_centrality(g, leaders, "closeness", mode = "in")
}
