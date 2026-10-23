\name{net_plot_panels}
\alias{net_plot_panels}
\alias{net_plot_layers}
\alias{net_plot_snapshots}
\alias{net_plot_egos}
\title{Small multiples: layers, snapshots over time, and ego networks}
\description{
Draws one panel per layer of a multilayer network
(\code{net_plot_layers()}), per time window or wave of a temporal network
(\code{net_plot_snapshots()}), or per ego of an egocentric dataset
(\code{net_plot_egos()}). Panels of one network share its layout, so that
a node keeps its place from panel to panel and only the ties change.
}
\usage{
net_plot_layers(x, layout = NULL, ..., ncol = NULL, seed = NULL)

net_plot_snapshots(x, width = NULL, breaks = NULL, layout = NULL, ...,
                   ncol = NULL, seed = NULL)

net_plot_egos(x, egos = NULL, order = 1L, ..., ncol = NULL, max_panels = 12L,
              seed = NULL)
}
\arguments{
\item{x}{A multilayer net_frame (\code{net_plot_layers()}), a temporal
net_frame (\code{net_plot_snapshots()}), or a net_frame with or without an
ego role (\code{net_plot_egos()}).}
\item{layout}{The shared layout, as in
\code{\link[=plot.net_frame]{plot()}}; \code{NULL} lays out the network
with every edge (all layers, or all times) once.}
\item{...}{Further arguments of \code{\link[=plot.net_frame]{plot()}} for
every panel, such as \code{node_label} or \code{node_size}.}
\item{ncol}{Columns of the grid of panels; \code{NULL} chooses a near
square grid.}
\item{seed}{Seed for the layout.}
\item{width, breaks}{The time windows, as in \code{\link{net_snapshots}}:
windows of length \code{width}, or between successive \code{breaks}. When
both are \code{NULL}, a panel network gets one panel per wave and other
temporal networks four equal windows spanning their times.}
\item{egos}{Identifiers of the egos to draw; \code{NULL} draws the first
\code{max_panels} egos (or nodes, for a whole network).}
\item{order}{For a whole network, the order of the ego networks extracted
with \code{\link{net_ego_frame}}.}
\item{max_panels}{The most panels to draw.}
}
\details{
\code{net_plot_layers()} draws the intralayer edges of each layer
(\code{\link{net_split_layers}}) on every actor. \code{net_plot_snapshots()}
draws the edges active in each window (\code{\link{net_snapshots}},
keeping all nodes); when nodes have presence spells, nodes absent from a
window are drawn light gray. \code{net_plot_egos()} draws each ego network
with its own layout, ego in the accent colour; a whole network is first
split with \code{\link{net_ego_frame}}. Legends are left out, since the
panels share their encodings.
}
\value{
The shared layout, invisibly (\code{net_plot_egos()}: the list of ego
networks drawn).
}
\seealso{
\code{\link[=plot.net_frame]{plot()}}, \code{\link{net_plot_timeline}}.
}
\examples{
# Krackhardt's relations as layers.
ml <- net_as_multilayer(net_example("krackhardt"), couplings = "none")
net_plot_layers(ml, node_label = FALSE, seed = 1)

# Sampson's liking ties at three waves.
liking <- net_example("sampson", relation = "liking")
net_plot_snapshots(liking, node_label = FALSE, seed = 1)

# Contacts in windows of a chosen width.
calls <- data.frame(source = c("a", "b", "c", "a", "d"), target = c("b", "c", "d", "d", "b"),
                    time = c(1, 2, 4, 6, 7))
g <- net_frame(calls, directed = FALSE, start = "time")
net_plot_snapshots(g, width = 3, seed = 1)

# First-order ego networks of three Florentine families.
net_plot_egos(net_example("florentine"), egos = c("Medici", "Strozzi", "Albizzi"),
              seed = 1)
}
