\name{plot.qap_cor}
\alias{plot.qap_cor}
\title{Plot a QAP correlation test}
\description{
Draws the permutation (reference) distribution of the correlation as a
histogram, marks the observed correlation, and reports the p-values, so
that the test can be read the way Krackhardt (1987) presents it: where the
observed value falls among the correlations the relabelled networks
produce.
}
\usage{
\method{plot}{qap_cor}(x, breaks = "Sturges", ...)
}
\arguments{
\item{x}{A \code{\link{qap_cor}} result.}
\item{breaks}{Passed to \code{\link[graphics]{hist}}.}
\item{...}{Further arguments passed to \code{\link[graphics]{hist}}, such
as \code{main}.}
}
\details{
The histogram shows the finite correlations among the draws
(\code{x$dist}); draws whose correlation is undefined are left out and
counted in the subtitle. The vertical line is the observed correlation
\code{x$statistic} and the subtitle gives \code{Pr(>=|r|)}, the two-sided
p-value, together with the number of permutations.
}
\value{\code{x}, invisibly.}
\references{
Krackhardt, D. (1987). QAP partialling as a test of spuriousness.
\emph{Social Networks}, 9(2), 171-186. \doi{10.1016/0378-8733(87)90012-8}
}
\seealso{\code{\link{qap_cor}}, \code{\link{plot.mrqap}}.}
\examples{
set.seed(1)
a <- matrix(rbinom(100, 1, 0.3), 10)
b <- (a + matrix(rbinom(100, 1, 0.1), 10)) > 0
q <- qap_cor(a, b, reps = 500, seed = 1, threads = 2)
plot(q)
}
