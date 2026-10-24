\name{net_ego_frame}
\alias{net_ego_frame}
\alias{net_split_egos}
\alias{net_to_egor}
\alias{as_net_frame.egor}
\title{Ego networks in one net_frame}
\description{
The ego network of a node at order \eqn{k} is the node (ego), every node
within \eqn{k} steps of it (its alters), and the edges among them
(Marsden 2002). An egocentric dataset is a collection of ego networks,
one per sampled ego, with no ties across egos (Burt 1984; Perry et al.
2018). One net_frame holds them all: the \code{ego} role of
\code{\link{net_frame}} says which ego each node and edge belongs to.
\code{net_ego_frame()} extracts the ego networks of a whole network,
\code{net_split_egos()} returns one network per ego, and
\code{net_to_egor()} and \code{as_net_frame()} convert to and from the
three tables of the egor package.
}
\usage{
net_ego_frame(x, egos = NULL, order = 1L, mode = c("all", "out", "in"),
              include_ego = TRUE, sep = ":", threads = 0L)

net_split_egos(x)

net_to_egor(x)

\method{as_net_frame}{egor}(x, directed = FALSE, sep = ":", quiet = FALSE, ...)
}
\arguments{
\item{x}{For \code{net_ego_frame()}, a whole network (a net_frame without an
ego role); for \code{net_split_egos()} and \code{net_to_egor()}, a
net_frame with an ego role; for \code{as_net_frame()}, an \code{egor}
object.}
\item{egos}{Identifiers of the egos; \code{NULL} uses every node.}
\item{order}{How many steps from ego an alter may be (at least 1).}
\item{mode}{In a directed network, which paths count: \code{"out"}
(nodes ego reaches), \code{"in"} (nodes that reach ego), or \code{"all"}
(ignoring direction; the default).}
\item{include_ego}{Whether each ego network contains its ego and ego's
ties (\code{TRUE}, the default) or only the alters and the ties among
them (the "order 1.5" network without ego, as in structural-hole
measures; \code{FALSE}).}
\item{sep}{The separator used to build alter identifiers from the ego's
and the alter's, as in \code{"ann:bob"}, since an alter named by two
egos is a node row in each of their networks.}
\item{threads}{Number of threads; \code{0} uses all available cores.}
\item{directed}{Whether the ties of the egor object are directed.}
\item{quiet}{Whether to suppress the advisory notes.}
\item{...}{Ignored.}
}
\details{
\code{net_ego_frame()} copies each ego's alters into new node rows with
identifiers \code{ego:alter} (an ego keeps its own identifier), keeps the
original identifier in a column \code{node} and the ego in a column
\code{ego} (renamed with a trailing underscore when those names are
taken), and copies the edges among them with their attributes. The ego
column becomes the \code{ego} role of both tables and
\code{egos_are_nodes} is declared from \code{include_ego}. The networks
of different egos share no nodes, so every measure computed on the result
can be read per ego (for example by splitting its result by
\code{x$nodes$ego}).

The egor package stores egocentric data as three tables: egos, alters
with their ego, and alter-alter ties with their ego. The conversion maps
egos to node rows whose ego is themselves, alters to node rows with their
ego (identifier \code{ego:alter}), ego-alter ties to edges marked
\code{tie = "ego-alter"}, and alter-alter ties to edges marked
\code{tie = "alter-alter"}; \code{net_to_egor()} reverses it and needs
every ego as a node row. egor must be installed for \code{net_to_egor()}.
\code{\link{net_plot_egos}} draws ego networks side by side.
}
\value{
\code{net_ego_frame()} and \code{as_net_frame()}: a net_frame with an ego
role. \code{net_split_egos()}: a named list of net_frames, one per ego, in
order of first appearance. \code{net_to_egor()}: an \code{egor} object.
}
\references{
Burt, R. S. (1984). Network items and the General Social Survey.
\emph{Social Networks}, 6(4), 293-339.
\doi{10.1016/0378-8733(84)90007-8}

Marsden, P. V. (2002). Egocentric and sociocentric measures of network
centrality. \emph{Social Networks}, 24(4), 407-422.
\doi{10.1016/S0378-8733(02)00016-3}

Perry, B. L., Pescosolido, B. A., & Borgatti, S. P. (2018).
\emph{Egocentric Network Analysis: Foundations, Methods, and Models}.
Cambridge University Press. \doi{10.1017/9781316443255}
}
\seealso{
\code{\link{net_frame}} for the \code{ego} role, \code{\link{net_ego}} for
ego network measures, and \code{\link{net_plot_egos}}.
}
\examples{
# Egocentric data entered directly: two egos with their alters.
nodes <- data.frame(id = c("ann", "a1", "a2", "bob", "b1"),
                    ego = c("ann", "ann", "ann", "bob", "bob"))
edges <- data.frame(source = c("ann", "ann", "a1", "bob"),
                    target = c("a1", "a2", "a2", "b1"),
                    ego = c("ann", "ann", "ann", "bob"))
e <- net_frame(edges, nodes, directed = FALSE, ego = "ego", edge_ego = "ego",
               egos_are_nodes = TRUE)
e
names(net_split_egos(e))

# Ego networks derived from Krackhardt's advice network.
advice <- net_example("krackhardt", relation = "advice")
egos <- net_ego_frame(advice, egos = 1:3, threads = 1)
egos
sapply(net_split_egos(egos), function(g) nrow(g$nodes))
net_plot_egos(advice, egos = 1:3, node_label = FALSE, seed = 1)

if (requireNamespace("egor", quietly = TRUE)) {
  three <- net_to_egor(e)
  as_net_frame(three)
}
}
