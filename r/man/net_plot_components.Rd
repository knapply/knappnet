\name{net_plot_components}
\alias{net_plot_components}
\title{Plot the distribution of component sizes}
\description{
Draws how many components of each size the network has, with the largest
component, which most measures of connectivity are about, marked.
}
\usage{
net_plot_components(x, mode = c("weak", "strong"), log = FALSE, ...)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{mode}{\code{"weak"} or \code{"strong"}, as in
\code{\link{net_components}}.}
\item{log}{Whether to use logarithmic axes (sizes and counts).}
\item{...}{Further arguments passed to \code{\link[graphics]{barplot}} (or
\code{\link[graphics]{plot}} with \code{log = TRUE}), such as
\code{main}.}
}
\value{A data frame with one row per distinct size, invisibly:
\code{size} and \code{count} (the number of components of that size).}
\seealso{\code{\link{net_components}}, \code{\link{net_fragmentation}}.}
\examples{
g <- net_example("noordin", relation = "friendship")
net_components(g)$no
net_plot_components(g)
}
