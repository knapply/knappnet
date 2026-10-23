\name{net_panel}
\alias{net_panel}
\alias{net_waves}
\title{Panel waves as a temporal network}
\description{
Longitudinal network data observed at discrete waves (the panel design of
stochastic actor-oriented models, Snijders, van de Bunt, and Steglich
2010) are stored as one temporal network in which wave \eqn{k} is the
spell \eqn{[k, k + 1)}, so windows, snapshots, rolling measures, and the
other temporal functions apply unchanged. \code{net_panel()} builds that
network from a list of static networks, one per wave, and
\code{net_waves()} splits it back.
}
\usage{
net_panel(networks, wave = "wave", node_start = "wave_first",
          node_end = "wave_last", quiet = FALSE)

net_waves(x, nodes = c("present", "all"))
}
\arguments{
\item{networks}{A list of \code{\link[=net_frame]{netframe}}s, one per
wave, in order, all with the same direction and the same identifier and
endpoint column names and the same edge columns. Names that are whole
numbers (\code{"1"}, \code{"2"}, ...) are used as the wave numbers;
otherwise the waves are numbered 1, 2, ....}
\item{wave}{The name of the edge column that will hold the wave number.}
\item{node_start, node_end}{The names of the node columns that will hold
each node's first wave and the end of its last wave (its presence spell).
They are added only when the composition changes, that is when some node
is absent from some wave.}
\item{quiet}{Whether to suppress the advisory notes of validation.}
\item{x}{A temporal netframe with a wave column, as built by
\code{net_panel()} or \code{net_temporal(wave = )}.}
\item{nodes}{Which nodes each wave gets back: those present in that wave
according to the node spells (the default), or all of them.}
}
\details{
The node table of the panel is the union of the waves' node tables, with
each node's attributes taken from the first wave it appears in. A node
present in waves 2 and 3 of 4 gets the spell \eqn{[2, 4)}: present from
the start of wave 2 to the end of wave 3. Presence is assumed contiguous;
a node that leaves and returns is recorded as present throughout. The edge
table stacks the waves' edge tables with the wave number appended, so an
edge observed in two waves is two rows (two spells of one tie, which are
not parallel edges; see \code{\link{net_simplify}}).

Wave boundaries follow the half-open convention of \code{\link{net_temporal}}:
\code{net_at(x, 2)} and \code{net_window(x, 2, 3)} both give wave 2, and
\code{net_window(x, 1, 3)} waves 1 and 2 together. The list of networks
remains the shape for exchange with RSiena and ORA.
}
\value{
\code{net_panel()} returns a netframe of class \code{temporal_netframe}
with the wave role (and node spells when the composition changes).
\code{net_waves()} returns a named list of static netframes, one per wave
in increasing order, each keeping every column and role but the wave and
the node spells.
}
\references{
Snijders, T. A. B., van de Bunt, G. G., & Steglich, C. E. G. (2010).
Introduction to stochastic actor-based models for network dynamics.
\emph{Social Networks}, 32(1), 44-60. \doi{10.1016/j.socnet.2009.02.004}

Holme, P., & Saramäki, J. (2012). Temporal networks. \emph{Physics
Reports}, 519(3), 97-125. \doi{10.1016/j.physrep.2012.03.001}
}
\seealso{
\code{\link{net_temporal}}, \code{\link{net_at}}, \code{\link{net_snapshots}},
\code{\link{net_rolling}}.
}
\examples{
w1 <- net_frame(data.frame(source = c("a", "b"), target = c("b", "c")), directed = FALSE)
w2 <- net_frame(data.frame(source = c("a", "b", "c"), target = c("b", "c", "d")),
                directed = FALSE)
w3 <- net_frame(data.frame(source = c("b", "c"), target = c("c", "d")), directed = FALSE)
panel <- net_panel(list(w1, w2, w3))
panel
panel$nodes
net_at(panel, 2)$edges
net_rolling(panel, net_density, breaks = 1:4)
waves <- net_waves(panel)
waves[["3"]]$nodes$id
}
