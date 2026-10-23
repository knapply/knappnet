\name{subnetwork}
\alias{subnetwork}
\alias{filter_nodes}
\alias{remove_isolates}
\alias{remove_pendants}
\alias{largest_component}
\title{Subnetworks}
\description{
Induced subnetworks keep a set of nodes and every edge between them:
\describe{
  \item{\code{subnetwork()}}{Keeps the nodes named in \code{nodes} (UCINET:
    Data > Extract).}
  \item{\code{filter_nodes()}}{Keeps the nodes for which a condition on the
    columns of \code{x$nodes} holds, written as in \code{\link{subset}}: for
    example \code{filter_nodes(x, role == "leader")}.}
  \item{\code{remove_isolates()}}{Drops nodes with no ties to other nodes
    (UCINET: Data > Remove isolates).}
  \item{\code{remove_pendants()}}{Drops nodes with exactly one neighbor,
    once or repeatedly (UCINET: Data > Remove pendants).}
  \item{\code{largest_component()}}{Keeps the largest weak or strong
    component (UCINET: Data > Extract main component).}
}
The menu paths follow Hanneman and Riddle (2005, chapter 6).
}
\usage{
subnetwork(x, nodes)

filter_nodes(x, condition)

remove_isolates(x)

remove_pendants(x, recursive = FALSE)

largest_component(x, mode = c("weak", "strong"))
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{nodes}{Node identifiers (values of \code{x$nodes[[id]]}) to keep, or a
logical vector with one entry per row of \code{x$nodes}.}
\item{condition}{A logical expression in the columns of \code{x$nodes};
missing values count as \code{FALSE}.}
\item{recursive}{Whether to keep removing pendants until none is left.}
\item{mode}{\code{"weak"} or \code{"strong"} components; see
\code{\link{components}}.}
}
\details{
Node and edge attributes, weights, direction, and network attributes are
kept; node order is preserved. Isolates and pendants are judged on the
undirected simple graph underlying the network: a node is an isolate when
it has no neighbor other than itself (a self-loop does not save it) and a
pendant when it has exactly one such neighbor, whatever the direction and
number of edges between them. With \code{recursive = TRUE},
\code{remove_pendants()} repeats until no pendant is left, which prunes
trees down to their cycles; otherwise only the nodes that are pendants in
the original network are dropped. When several components are largest, the
one containing the earliest node in \code{x$nodes} is kept.
}
\value{A netframe.}
\references{
Hanneman, R. A., & Riddle, M. (2005). \emph{Introduction to Social Network
Methods}. University of California, Riverside.
}
\seealso{\code{\link{ego_network}}, \code{\link{components}},
\code{\link{combine_relations}}}
\examples{
noordin <- example_network("noordin")
names(noordin$nodes)
# The 1-mode kinship network, without the people it does not touch.
kin <- remove_isolates(combine_relations(noordin, method = "select", which = "kinship"))
kin

# A star 1 - {2, 3, 4} with a path 4 - 5 - 6 and an isolate 7.
g <- netframe(data.frame(from = c(1, 1, 1, 4, 5), to = c(2, 3, 4, 5, 6)),
              nodes = data.frame(id = 1:7), directed = FALSE)
subnetwork(g, c(1, 2, 3))$nodes$id
filter_nodes(g, id > 3)$nodes$id
remove_isolates(g)$nodes$id
remove_pendants(g)$nodes$id                    # 2, 3, and 6 go
remove_pendants(g, recursive = TRUE)$nodes$id  # then 1 and 5; 4 and the isolate stay
largest_component(g)$nodes$id
}
