\name{plot_scores}
\alias{plot_scores}
\title{Plot ranked scores}
\description{
Draws scores such as centralities as a dot plot sorted from the highest
down, labelled with node identifiers: the quickest way to read off who
leads on a measure and by how much.
}
\usage{
plot_scores(scores, labels = NULL, top = 25, decreasing = TRUE, ...)
}
\arguments{
\item{scores}{A numeric vector of scores, one per node, as returned by
\code{\link{degree}}, \code{\link{betweenness}}, and the other centrality
functions.}
\item{labels}{Labels for the scores: a \code{\link{netframe}} whose nodes
the scores are aligned with (its identifiers are used), a character vector
of the same length, or \code{NULL} for \code{names(scores)} or, failing
that, positions.}
\item{top}{How many of the highest scores to show; \code{Inf} shows all.}
\item{decreasing}{Whether the highest score is at the top. \code{FALSE}
ranks from the lowest instead.}
\item{...}{Further arguments passed to \code{\link[graphics]{dotchart}},
such as \code{main} or \code{xlab}.}
}
\details{
Missing and \code{NaN} scores (undefined values) are left out and counted
in the subtitle. Ties keep the order of \code{scores}.
}
\value{A data frame of the scores drawn, in the order drawn from the top,
invisibly: \code{label} and \code{score}.}
\seealso{\code{\link{plot_centrality}} to draw the scores on the network.}
\examples{
g <- example_network("florentine")
plot_scores(betweenness(g), g, main = "Betweenness")
plot_scores(closeness(g), g, top = 5, xlab = "closeness")
}
