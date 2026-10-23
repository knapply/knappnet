\name{plot_rolling}
\alias{plot_rolling}
\title{Plot a measure over time}
\description{
Draws the result of \code{\link{rolling}}, a measure computed on every
window of a temporal network, as lines over time: one line for a
network-level measure, or one per node for a node-level measure.
}
\usage{
plot_rolling(x, nodes = NULL, top = 8, ...)
}
\arguments{
\item{x}{A data frame from \code{\link{rolling}}, with columns
\code{start}, \code{end}, \code{node}, and \code{value}.}
\item{nodes}{For node-level results, the identifiers of the nodes to draw.
\code{NULL} draws the \code{top} nodes by mean value.}
\item{top}{How many nodes to draw when \code{nodes} is \code{NULL}.}
\item{...}{Further arguments passed to \code{\link[graphics]{plot}}, such
as \code{main} or \code{ylab}.}
}
\details{
Each value is drawn at the start of its window. Windows where a node's
value is undefined (\code{NaN}) leave a gap in its line. The time axis
uses the type of the window times (numbers, dates, or date-times).
}
\value{\code{x}, invisibly.}
\seealso{\code{\link{rolling}}, \code{\link{snapshots}},
\code{\link{plot_timeline}}.}
\examples{
logins <- data.frame(user = c("u1", "u2", "u1", "u3", "u2", "u3", "u1"),
                     host = c("h1", "h1", "h2", "h2", "h3", "h1", "h1"),
                     time = c(0, 1, 5, 6, 12, 13, 14))
g <- temporal_network(logins)
plot_rolling(rolling(g, degree, width = 5, step = 2), main = "Degree over time")
plot_rolling(rolling(g, function(w) nrow(w$edges), width = 5, step = 2),
             ylab = "contacts in window")
}
