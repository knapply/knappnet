\name{net_dichotomize}
\alias{net_dichotomize}
\title{Dichotomize a weighted network}
\description{
Keeps the edges whose weight passes a test against a threshold and drops
the weight role, giving a binary network (UCINET: Transform > Dichotomize,
with the same cut-off operators; Hanneman and Riddle 2005, chapter 6). Edges
that fail the test are removed; the weight column stays on the surviving
edges as an ordinary attribute, and every other column is kept.
}
\usage{
net_dichotomize(x, threshold = 0, op = c(">", ">=", "<", "<=", "==", "!="), column = NULL)
}
\arguments{
\item{x}{A \code{\link{net_frame}}.}
\item{threshold}{The cut-off value.}
\item{op}{The comparison applied as \code{weight op threshold}: \code{">"}
(the default, UCINET's "greater than"), \code{">="}, \code{"<"},
\code{"<="}, \code{"=="}, or \code{"!="}.}
\item{column}{The numeric edge column to test, by name or position;
\code{NULL} uses the network's weight column.}
}
\value{An unweighted netframe.}
\references{
Hanneman, R. A., & Riddle, M. (2005). \emph{Introduction to Social Network
Methods}. University of California, Riverside.
}
\seealso{\code{\link{net_symmetrize}}, \code{\link{net_combine_relations}}}
\examples{
ji <- net_example("jemaah_islamiyah", weighted = TRUE)
ji
# Keep only the stronger ties (interaction strength above 2).
strong <- net_dichotomize(ji, threshold = 2)
strong
nrow(strong$edges)
nrow(net_dichotomize(ji, threshold = 2, op = ">=")$edges)
}
