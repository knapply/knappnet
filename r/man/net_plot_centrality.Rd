\name{net_plot_centrality}
\alias{net_plot_centrality}
\title{Plot a score on the network}
\description{
Draws the network with \code{\link{plot.net_frame}}, sizing and coloring
each node by a score such as a centrality, with a legend named after the
score.
}
\usage{
net_plot_centrality(x, scores, title = "score", color = TRUE, ...)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{scores}{A numeric vector aligned with the rows of \code{x$nodes}, as
returned by \code{\link{net_degree}}, \code{\link{net_betweenness}}, and the other
centrality functions, or the name of a numeric column of \code{x$nodes}.}
\item{title}{The name of the score, used for the legend.}
\item{color}{Whether to color the nodes by the score as well as sizing
them. \code{FALSE} keeps the default node color.}
\item{...}{Further arguments passed to \code{\link{plot.net_frame}}, such as
\code{layout}, \code{seed}, \code{node_label}, or \code{main}.}
}
\details{
Node area grows linearly with the score and the color follows the
sequential palette of \code{\link{plot.net_frame}}; undefined scores
(\code{NA}, \code{NaN}) are drawn at the smallest size in light gray.
}
\value{The layout used, invisibly, as \code{\link{plot.net_frame}} returns.}
\seealso{\code{\link{net_plot_scores}} for a ranked dot plot of the same
scores.}
\examples{
g <- net_example("florentine")
net_plot_centrality(g, net_betweenness(g), title = "betweenness", seed = 1)
net_plot_centrality(g, net_degree(g), title = "degree", color = FALSE, layout = "circle")
}
