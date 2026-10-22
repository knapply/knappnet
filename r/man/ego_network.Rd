\name{ego_network}
\alias{ego_network}
\alias{ego_networks}
\title{Ego networks}
\usage{
ego_network(x, node, order = 1L, mode = c("all", "out", "in"),
  include_ego = TRUE)

ego_networks(x, nodes = NULL, order = 1L, mode = c("all", "out", "in"),
  include_ego = TRUE, threads = 0L)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{node}{The identifier of the ego, as in the identifier column of
\code{x$nodes}.}
\item{nodes}{Identifiers of the egos; \code{NULL} (the default) uses every
node.}
\item{order}{Alters are the nodes within this many steps of ego. \code{1}
gives the first-order neighborhood.}
\item{mode}{For directed networks, \code{"all"}, \code{"out"}, or
\code{"in"}; see Details.}
\item{include_ego}{Whether ego itself is part of the result (igraph's
\code{mindist = 0}) or only its alters (\code{mindist = 1}).}
\item{threads}{Number of threads; \code{0} uses all available cores.}
}
\value{
\code{ego_network()} returns a \code{\link{netframe}}. \code{ego_networks()}
returns a list of netframes named by the egos' identifiers.
}
\description{
\code{ego_network()} extracts the ego network of one node: the node (ego),
every node within \code{order} steps of it (its alters), and all ties among
them. \code{ego_networks()} does the same for several nodes, or for every
node, at once.
}
\details{
The result is the subgraph induced by ego and its alters, as igraph's
\code{make_ego_graph()} and NetworkX's \code{ego_graph()} return: every edge
whose two endpoints are both in the ego network is kept, with all its
columns, and so is every column of the node table. Nodes and edges keep
their order in \code{x}. The weight column and its meaning carry over, and
so does the node type of a two-mode network (see \code{\link{twomode}}).

In directed networks \code{mode} chooses the alters: \code{"out"} takes the
nodes ego reaches along outgoing ties (whom ego names), \code{"in"} the
nodes that reach ego (who names ego), and \code{"all"} ignores direction.
Ties among the chosen nodes are kept whatever their direction. Undirected
networks ignore \code{mode}. Steps are counted without weights.
}
\references{
Wasserman, S., & Faust, K. (1994). \emph{Social Network Analysis: Methods
and Applications}. Cambridge University Press.
\doi{10.1017/CBO9780511815478}

Perry, B. L., Pescosolido, B. A., & Borgatti, S. P. (2018).
\emph{Egocentric Network Analysis: Foundations, Methods, and Models}.
Cambridge University Press. \doi{10.1017/9781316443255}
}
\seealso{
\code{\link{ego_summary}} for measures of every ego network at once.
}
\examples{
g <- netframe(data.frame(from = c("a", "a", "b", "c"), to = c("b", "c", "c", "d")),
              directed = FALSE)
ego_network(g, "c")
ego_network(g, "c", include_ego = FALSE)$nodes
lengths(lapply(ego_networks(g), function(e) e$nodes$name))
}
