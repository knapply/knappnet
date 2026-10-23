\name{plot_component_sizes}
\alias{plot_component_sizes}
\title{Plot the distribution of component sizes}
\description{
Draws how many components of each size the network has, with the largest
component, which most measures of connectivity are about, marked.
}
\usage{
plot_component_sizes(x, mode = c("weak", "strong"), log = FALSE, ...)
}
\arguments{
\item{x}{A \code{\link{netframe}}.}
\item{mode}{\code{"weak"} or \code{"strong"}, as in
\code{\link{components}}.}
\item{log}{Whether to use logarithmic axes (sizes and counts).}
\item{...}{Further arguments passed to \code{\link[graphics]{barplot}} (or
\code{\link[graphics]{plot}} with \code{log = TRUE}), such as
\code{main}.}
}
\value{A data frame with one row per distinct size, invisibly:
\code{size} and \code{count} (the number of components of that size).}
\seealso{\code{\link{components}}, \code{\link{fragmentation}}.}
\examples{
g <- example_network("noordin", relation = "friendship")
components(g)$no
plot_component_sizes(g)
}
