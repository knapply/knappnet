\name{degree}
\alias{degree}
\title{Degree and strength}
\usage{
degree(x, mode = c("all", "out", "in"), loops = TRUE, weighted = FALSE)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{mode}{For directed networks, \code{"all"}, \code{"out"}, or
\code{"in"}. Ignored for undirected networks.}
\item{loops}{Whether self-loops count.}
\item{weighted}{Whether to sum edge weights.}
}
\value{A numeric vector aligned with the rows of \code{x$nodes}.}
\description{
Counts edge endpoints at each node. With \code{weighted = TRUE}, sums edge
weights instead (often called strength).
}
\details{
Parallel edges each count. In undirected networks a self-loop adds 2 to its
node; in directed networks it adds 1 to out-degree and 1 to in-degree. These
match igraph and NetworkX. Set \code{loops = FALSE} to ignore self-loops,
which matches the sna default (\code{diag = FALSE}).
}
\examples{
g <- netframe(data.frame(from = c(1, 1, 2), to = c(2, 3, 3)))
degree(g, mode = "out")
}
