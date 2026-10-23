\name{plot_centrality}
\alias{plot_centrality}
\title{Draw a network with nodes sized and colored by a score}
\description{
Draws the network with \code{\link{plot.netframe}}, sizing and coloring
each node by a score such as a centrality, with a legend named after the
score.
}
\usage{
plot_centrality(x, scores, title = "score", color = TRUE, ...)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{scores}{A numeric vector aligned with the rows of \code{x$nodes}, as
returned by \code{\link{degree}}, \code{\link{betweenness}}, and the other
centrality functions, or the name of a numeric column of \code{x$nodes}.}
\item{title}{The name of the score, used for the legend.}
\item{color}{Whether to color the nodes by the score as well as sizing
them. \code{FALSE} keeps the default node color.}
\item{...}{Further arguments passed to \code{\link{plot.netframe}}, such as
\code{layout}, \code{seed}, \code{node_label}, or \code{main}.}
}
\details{
Node area grows linearly with the score and the color follows the
sequential palette of \code{\link{plot.netframe}}; undefined scores
(\code{NA}, \code{NaN}) are drawn at the smallest size in light gray.
}
\value{The layout used, invisibly, as \code{\link{plot.netframe}} returns.}
\seealso{\code{\link{plot_scores}} for a ranked dot plot of the same
scores.}
\examples{
g <- example_network("florentine")
plot_centrality(g, betweenness(g), title = "betweenness", seed = 1)
plot_centrality(g, degree(g), title = "degree", color = FALSE, layout = "circle")
}
