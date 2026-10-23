\name{plot_robustness}
\alias{plot_robustness}
\title{Plot the robustness of a network to node removal}
\description{
Removes nodes one at a time in a given order and tracks how the network
falls apart: the share of nodes in the largest component, or the
fragmentation of Borgatti (2006), after each removal. Curves for several
removal strategies (highest degree first, highest betweenness first,
random) can be compared with the removal of a chosen set such as the
\code{\link{key_players}}, as in the attack-tolerance analyses of Albert,
Jeong, and Barabási (2000) and Holme et al. (2002).
}
\usage{
plot_robustness(x, nodes = NULL, strategy = c("degree", "betweenness", "random"),
                measure = c("largest_component", "fragmentation"), steps = 25,
                mode = c("weak", "strong"), seed = NULL, threads = 0L,
                label = "selected nodes", ...)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{nodes}{Node identifiers to remove, in order, for example
\code{key_players(x, k)$set}. \code{NULL} draws only the strategies.}
\item{strategy}{Removal orders to draw, any of \code{"degree"},
\code{"betweenness"}, and \code{"random"}; \code{NULL} for none.}
\item{measure}{\code{"largest_component"} (its share of the original
nodes) or \code{"fragmentation"} (\code{\link{fragmentation}}).}
\item{steps}{How many removal counts to evaluate along each strategy
curve, spread evenly from none to all but one node (the curve for
\code{nodes} is evaluated after every removal). The network is rebuilt at
each point, so the cost is \code{steps} times that of
\code{\link{components}}.}
\item{mode}{\code{"weak"} or \code{"strong"} components, as in
\code{\link{components}}.}
\item{seed}{Seed for the random order (see \code{\link{set.seed}}).}
\item{threads}{Number of threads for betweenness, a whole number; \code{0}
uses all available cores, or 2 during \verb{R CMD check} (see \sQuote{Threads}
in \code{\link{netframe}}).}
\item{label}{Legend label for the curve of \code{nodes}.}
\item{...}{Further arguments passed to \code{\link[graphics]{plot}}, such
as \code{main}.}
}
\details{
Strategy orders are fixed at the start from the initial degrees or
betweenness scores (Holme et al.'s "initial" attacks); they are not
recalculated after each removal. Degree counts every edge end and
betweenness is the unnormalized score of \code{\link{betweenness}} with
its default weighting. Ties are broken by node order. Removing a node
removes its edges; the remaining network keeps every other node, so the
largest-component share is measured against the original number of nodes.
}
\value{A data frame with one row per evaluated point, invisibly:
\code{curve} (the strategy or \code{label}), \code{removed} (nodes
removed), \code{fraction} (of the original nodes), and \code{value} (the
measure).}
\references{
Albert, R., Jeong, H., & Barabási, A.-L. (2000). Error and attack
tolerance of complex networks. \emph{Nature}, 406(6794), 378-382.
\doi{10.1038/35019019}

Holme, P., Kim, B. J., Yoon, C. N., & Han, S. K. (2002). Attack
vulnerability of complex networks. \emph{Physical Review E}, 65(5),
056109. \doi{10.1103/PhysRevE.65.056109}

Borgatti, S. P. (2006). Identifying sets of key players in a social
network. \emph{Computational and Mathematical Organization Theory}, 12(1),
21-34. \doi{10.1007/s10588-006-7084-x}
}
\seealso{\code{\link{key_players}}, \code{\link{components}},
\code{\link{fragmentation}}.}
\examples{
g <- example_network("jemaah_islamiyah")
kp <- key_players(g, 3, "fragmentation", seed = 1, threads = 2)
plot_robustness(g, nodes = kp$set, label = "key players", seed = 1, threads = 2)
plot_robustness(g, strategy = "degree", measure = "fragmentation", threads = 2)
}
